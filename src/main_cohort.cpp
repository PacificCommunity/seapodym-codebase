#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <numeric>	  // std::accumulate
#include <algorithm>	  // std::copy
#include <map>
#include <utility>
#include <vector>
#include <mpi.h>
#include "VarParamCoupled.h"
#include "SeapodymCohortDependencyAnalyzer.h"
#include "TaskStepWorker.h"
#include "DistDataCollector.h"
#include "TaskStepManager.h"
#include "SeapodymCohort.h"
#include <CmdLineArgParser.h>
#include "ctrace.h"
#include "DataProvider.h"
#include "Tags.h"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

// Tags for the A+ ping-pong protocol between a feeder cohort and the
// dedicated A+ worker rank (must not clash with Tags.h: 0-3).
#define APLUS_NOTIFY_TAG    10
#define APLUS_ACK_TAG       11
#define APLUS_CHECKSUM_TAG  12

double time_ic_comm = 0.0, time_ic_flush = 0.0, time_ic_copy = 0.0, time_spawning = 0.0, time_getdata = 0.0, time_xreset = 0.0, time_init_cohort_spawning = 0.0, time_init_cohort_restart = 0.0, time_io_forcing = 0.0;
long   n_ic = 0;   // count of spawning-path inits, for per-init averages
double time_mpi_put = 0.0, time_send = 0.0, time_idle = 0.0;

SeapodymCohort* seapodym_cohort(const char* parfile, const int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs);
void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
void buffers_set(long int &mv, long int &mc, long int &mg);


double time_worker_init = 0.0, time_cohort_init = 0.0, time_calc = 0.0, time_mpi = 0.0, time_step = 0.0, time_overhead = 0.0;

SeapodymCohort xinit_prerun_wrapper(const char* parfile, bool useAPlus) {

	double tik = MPI_Wtime();

	SeapodymCohort cohort((char*)parfile, 0, useAPlus);

	//initialize variables of optimization
	const int nvar = cohort.nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1,nvar);

	cohort.xinit(x, x_names);

	//prepare cohort run
	cohort.prerun_model();

	double tak = MPI_Wtime();

	time_worker_init += tak - tik;

	return cohort;
}


// Task function for ordinary (non-A+) cohorts. A+ is no longer part of the
// dependency graph at all (see main()); instead, the cohort that reaches the
// oldest normal age class one time step before the last (the "feeder") hands
// its data directly to the dedicated A+ worker rank via a blocking ping-pong,
// implemented right here at the point where that step's data would otherwise
// just go to dataCollector.
void taskFunction(int task_id, int stepBeg, int stepEnd, MPI_Comm comm,
		const std::shared_ptr<spdlog::logger>& logger,
		DistDataCollector* dataCollector,
		SeapodymCohort* cohort,
		int numAgeGroups, int numTimeSteps, int aPlusWorkerRank){

	static double last_task_end = -1.0;        // per-worker process, persists across calls
	double t_in = MPI_Wtime();
	if (last_task_end >= 0.0)
		time_idle += t_in - last_task_end;     // <-- time spent in worker.run() waiting for dispatch

	double tik = MPI_Wtime();

	logger->info("> task id {} for steps {} to {}", task_id, stepBeg, stepEnd);

	logger->info("    >> initialization of task id {}", task_id);
	//initialize variables of optimization
	const int nvar = cohort->nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1,nvar);
	cohort->xinit(x, x_names);

	int cohort_id = task_id;
	cohort->restart(cohort_id);
	//initialize cohort either from restart or from spawning
	cohort->init_cohort(x,*dataCollector,numTimeSteps);
	logger->info("    << initialization of task id {}", task_id);

	// advance the cohort
	double tak = MPI_Wtime();
	time_cohort_init += tak - tik;

	for (auto step = stepBeg; step < stepEnd; ++step) {

		double tik_step = MPI_Wtime();
		logger->info("        >>> step {} of task id {}", step, task_id);
		cohort->stepForward(false);
		logger->info("End of step {} of task id {}. Checksum = {}", step, task_id, cohort->Checksum());
		logger->info("        <<< step {} of task id {}", step, task_id);
		time_step += MPI_Wtime() - tik_step;

		// Send the data to the shared collector.
		logger->info("        >>> send data for step {} of task id {}", step, task_id);
		std::vector<double> localData = cohort->GetCohortDensity();
		int chunk_id = cohort->getChunkId(step);

		double t_p = MPI_Wtime();
		dataCollector->put(chunk_id, localData.data());
		time_mpi_put += MPI_Wtime() - t_p;
		logger->info("        <<< send data for step {} of task id {}", step, task_id);

		// ------------------------------------------------------------------
		// A+ ping-pong: this cohort feeds the dedicated A+ worker once it
		// reaches the oldest normal age class (one of the numTimeSteps-1
		// "feeder" cohorts). The payload is embedded directly in the NOTIFY
		// message - there is no shared RMA slot, so a second feeder's
		// message can never clobber this one. Only after the ACK do we tell
		// the manager this step is complete, so nothing downstream can be
		// dispatched before the A+ worker has folded this contribution in -
		// that ordering is what makes an explicit A+ dependency edge in the
		// task graph unnecessary.
		//
		// The message also carries the A+ calendar step this feeds
		// (task_id+1) as its first element. Feeders can complete - and thus
		// NOTIFY - in any order (nothing in the dependency graph forces task
		// task_id+1's feeding event to happen after task_id's: a new cohort
		// is typically born numAgeGroups-2 rows before its predecessor even
		// reaches its own feeding step), but A+'s dynamics are a genuine
		// sequential recurrence, so the A+ worker must process contributions
		// in calendar order regardless of arrival order - see runAPlusWorker().
		// ------------------------------------------------------------------
		if (aPlusWorkerRank >= 0 && step == numAgeGroups - 1 && task_id <= numTimeSteps - 2) {
			logger->info("        >>> A+ notify for task id {}", task_id);
			double t_n = MPI_Wtime();
			std::vector<double> msg(localData.size() + 1);
			msg[0] = double(task_id + 1); // A+ calendar step this feeds
			std::copy(localData.begin(), localData.end(), msg.begin() + 1);
			MPI_Send(msg.data(), (int)msg.size(), MPI_DOUBLE,
					 aPlusWorkerRank, APLUS_NOTIFY_TAG, MPI_COMM_WORLD);

			int ack;
			MPI_Recv(&ack, 1, MPI_INT, aPlusWorkerRank, APLUS_ACK_TAG,
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			time_mpi += MPI_Wtime() - t_n;
			logger->info("        <<< A+ notify for task id {}", task_id);
		}

		int success = task_id;
		int output[3] = {task_id, step, success};
		logger->info("        >>> notify manager after step {} of task id {}", step, task_id);
		double tik_mpi = MPI_Wtime();
		MPI_Send(output, 3, MPI_INT, 0, END_TASK_TAG, comm);
		time_mpi += MPI_Wtime() - tik_mpi;
		logger->info("        <<< notify manager after step {} of task id {}", step, task_id);
	}

	time_calc += MPI_Wtime() - tak;

	last_task_end = MPI_Wtime();

	logger->info("< task id {} for steps {} to {}", task_id, stepBeg, stepEnd);
}

// Chunk id under which the A+ worker publishes its density at calendar step
// t, in the A+ range appended after all the normal (task,step) chunks - see
// main()'s numChunks. Newly-spawned cohorts (SeapodymCohort::InitializeCohort,
// spawning branch) read this same chunk (nb_age_class*numTimeSteps+(t-1),
// with nb_age_class==numAgeGroups when A+ is enabled) to fold A+'s
// contribution into spawning biomass.
int inline aplusChunkId(int t, int numAgeGroups, int numTimeSteps) {
	return numAgeGroups * numTimeSteps + t;
}

///////////////////////////////////////////////////////////////////////////////

// Dedicated A+ (plus group) worker loop. Runs on its own MPI rank, entirely
// outside the TaskStepManager/TaskStepWorker dependency graph - ordering
// relative to the manager is enforced by the ping-pong in taskFunction()
// above (a feeder blocks on the ACK before telling the manager it is done),
// not by an explicit dependency edge.
//
// `cohort` is a single persistent SeapodymCohort reused for every A+ time
// step: right after stepForward(), its own dvarCohortDensity *is* "last
// step's A+ pool", so unlike the dependency-graph design this replaced, no
// round trip through the data collector is needed to recover it - only the
// newly-graduating cohort's data ever crosses the wire, and it does so
// directly in the NOTIFY message rather than through a shared RMA slot.
// restartAPlus()/init_cohort_aplus() fully re-pin every field stepForward()
// reads before each call, so reusing one object across many calls is safe -
// the same pattern the ordinary-cohort path already relies on for reusing
// one SeapodymCohort across many task_ids per worker.
//
// IMPORTANT: A+'s dynamics are a genuine sequential recurrence (step t needs
// the state as of step t-1, plus the correct calendar date and forcing data
// for step t) - not a commutative accumulation. Feeders can complete, and
// therefore NOTIFY this rank, in any order (see taskFunction()'s comment),
// so this loop cannot simply process the n-th message received as step n.
// Instead each NOTIFY carries the calendar step it feeds (prepended to the
// payload by the sender), and a small reorder buffer holds early arrivals
// until their turn comes, processing steps strictly in ascending order. The
// ACK to a feeder is sent only once its step has actually been processed -
// an earlier ACK would let that feeder tell the manager it is done before
// A+ has actually caught up to it, defeating the ordering guarantee this
// whole ping-pong exists to provide.
void runAPlusWorker(SeapodymCohort& cohort, DistDataCollector* dataCollector,
		int numAgeGroups, int numData, int numTimeSteps,
		const std::shared_ptr<spdlog::logger>& logger) {

	const int nvar = cohort.nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1, nvar);
	cohort.xinit(x, x_names);

	// Runs one A+ calendar step and publishes its result. `graduating` is
	// null for t=0 (seeded from file; no upstream feeder), otherwise it is
	// the feeder's density for this step, merged with A+'s own current
	// density (its state right after the previous call - see header comment).
	auto processStep = [&](int t, const std::vector<double>* graduating) {
		if (graduating == nullptr) {
			cohort.restartAPlus(t);
			cohort.init_cohort_aplus(x, std::vector<double>(), /*seedFromFile=*/true);
		} else {
			std::vector<double> prevAPlus = cohort.GetCohortDensity();
			std::vector<double> merged(prevAPlus.size());
			for (std::size_t k = 0; k < merged.size(); ++k)
				merged[k] = prevAPlus[k] + (*graduating)[k];
			cohort.restartAPlus(t);
			cohort.init_cohort_aplus(x, merged, /*seedFromFile=*/false);
		}
		cohort.stepForward(false);
		logger->info("A+ t={} done. Checksum = {}", t, cohort.Checksum());

		// Publish this step's density so a cohort spawned right after this
		// step can fold it into its spawning-biomass sum (see
		// SeapodymCohort::InitializeCohort's spawning branch).
		std::vector<double> out = cohort.GetCohortDensity();
		int chunk_id = aplusChunkId(t, numAgeGroups, numTimeSteps);
		dataCollector->put(chunk_id, out.data());
	};

	processStep(0, nullptr);

	// Chunk 0 (A+'s file-seeded opening balance) has no feeder and therefore
	// no graph-based dependency gating it - unlike every later row, nothing
	// in the task dependency graph stops the manager from dispatching the
	// very first spawning-based cohort (which needs this exact chunk) before
	// this rank has even finished its own startup and reached this point.
	// This used to be closed with an extra MPI_Barrier(MPI_COMM_WORLD) here
	// (matched by one in main()'s manager/worker branches below), but that
	// turned out to be redundant: DistDataCollector's own post-construction
	// barrier already guarantees the window is safe to put()/get() the
	// moment construction returns, and the actual chunk-0 read is already
	// gated by the ordinary dependency chain (a newborn cohort can't be
	// dispatched until every initial cohort - including the very first
	// feeder - has completed, and a feeder can't complete without an ACK,
	// which A+'s strict-order reorder buffer only sends once row 0 has
	// already been processed).

	// Reorder buffer: row -> (source rank to ACK, its density payload).
	std::map<int, std::pair<int, std::vector<double>>> pending;
	int nextRow = 1;
	std::vector<double> buf(numData + 1);

	while (nextRow < numTimeSteps) {

		MPI_Status status;
		logger->info(">>> A+ waiting for a feeder (next row = {})", nextRow);
		MPI_Recv(buf.data(), numData + 1, MPI_DOUBLE,
				 MPI_ANY_SOURCE, APLUS_NOTIFY_TAG, MPI_COMM_WORLD, &status);
		int row = (int)std::llround(buf[0]);
		logger->info("<<< A+ received feeder data for row {}", row);
		pending.emplace(row, std::make_pair(status.MPI_SOURCE,
				std::vector<double>(buf.begin() + 1, buf.end())));

		// Drain the buffer while the next expected row is already available -
		// a single arrival can unblock several buffered rows at once.
		while (pending.count(nextRow)) {
			auto it = pending.find(nextRow);
			const int srcRank = it->second.first;
			processStep(nextRow, &it->second.second);
			pending.erase(it);

			// Only now - after nextRow has actually been processed - can the
			// feeder that fed it be told to notify the manager.
			int ack = 1;
			MPI_Send(&ack, 1, MPI_INT, srcRank, APLUS_ACK_TAG, MPI_COMM_WORLD);

			++nextRow;
		}
	}

	// Signal the manager that every A+ step has been published to
	// dataCollector, so it is safe to read the A+ chunk range for the
	// final checksum.
	int done = 1;
	MPI_Send(&done, 1, MPI_INT, 0, APLUS_CHECKSUM_TAG, MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// MPI initialization
	MPI_Init(&argc, &argv);

	time_worker_init = 0;
	time_calc = 0;

	int size, workerId;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &workerId);
	if (size < 2) {
		std::cerr << "ERROR: must have at least 2 ranks\n";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	CmdLineArgParser cmdLine;
	cmdLine.set("-s", std::string("initparfile.xml"), "Input parameter file");
	cmdLine.set("-no-aplus", false, "Disable the A+ (plus group) accumulator and reproduce "
		"pre-A+ behavior/checksum, for regression comparison.");

	// Parse the command line arguments
	bool success = cmdLine.parse(argc, argv);
	bool help = cmdLine.get<bool>("-help") || cmdLine.get<bool>("-h");
	bool useAPlus = !cmdLine.get<bool>("-no-aplus");
	if (!success) {
		std::cerr << "Error parsing command line arguments." << std::endl;
		cmdLine.help();
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (help) {
		cmdLine.help();
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// A+ needs its own dedicated worker rank, on top of the manager and at
	// least one ordinary cohort worker.
	if (useAPlus && size < 3) {
		if (workerId == 0)
			std::cerr << "ERROR: A+ requires at least 3 ranks (manager + "
				"1+ cohort worker + 1 dedicated A+ worker). Use -no-aplus "
				"to run with 2 ranks.\n";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// logger
	// Use true to let logs be overwritten, otherwise the logs will be appended
	std::string sworkerId = std::to_string(workerId);
	auto logger = spdlog::basic_logger_mt("log" + sworkerId, "log_taskfunc" + sworkerId + ".txt", true);
	logger->set_level(spdlog::level::debug);

	std::string parfile = cmdLine.get<std::string>("-s");

	// Read parfile and map
	VarParamCoupled param;
	PMap map;
	param.init_param();
	param.read(parfile);

	// Get number of time steps and number of cohorts from param.
	// numAgeGroups excludes the A+ (plus group) bin: the diagonal cohort-task
	// scheme ages "normal" cohorts through numAgeGroups steps, and the A+ bin
	// (age index sp_nb_cohorts[0]-1) is handled by a dedicated worker rank,
	// fed by a ping-pong with the oldest normal cohort each step (see
	// runAPlusWorker()/taskFunction() above) rather than being part of the
	// SeapodymCohortDependencyAnalyzer graph. With -no-aplus, numAgeGroups
	// reverts to sp_nb_cohorts[0] and A+ is just the last ordinary aging
	// cohort, reproducing pre-A+ behavior.
	int numAgeGroups = param.sp_nb_cohorts[0] - (useAPlus ? 1 : 0);
	int Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps;
	Date::init_time_variables(param, Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps, 0,0);
	//int numTasks = numAgeGroups + numTimeSteps - 1;

	// Size of map (useful to access to a specific position adress of the 4D array pointer storing density)
	map.lit_map(param);
	int numData = map.get_state_array_size();

	//Set-up the size for the shared arrays for forcing data
	std::vector<std::pair<std::string, std::size_t>> nameSizePairs = param.getDpNameSizePairs(numTimeSteps, map.get_array_size());

	// Set up the data collector. Normal (task,step) chunks come first; when
	// A+ is enabled, one extra chunk per time step is appended for the A+
	// worker to publish its density into (see aplusChunkId()) - each chunk
	// is uniquely owned by one calendar step, so unlike the shared RMA slot
	// this design replaced, there is nothing for two writers to race on.
	int numChunks = numAgeGroups * numTimeSteps + (useAPlus ? numTimeSteps : 0);

	// A+ is never part of the dependency graph now - it's handled entirely
	// by the ping-pong between the oldest normal cohort and the dedicated A+
	// worker rank, so aPlusCohort is always false here (see header comment
	// on numAgeGroups above).
	SeapodymCohortDependencyAnalyzer taskDeps(numAgeGroups, numTimeSteps, param.age_mature[0], /*aPlusCohort=*/false);
	int numCohorts = taskDeps.getNumberOfCohorts();
	std::map<int, int> stepBegMap = taskDeps.getStepBegMap();
	std::map<int, int> stepEndMap = taskDeps.getStepEndMap();
	std::map<int, std::set<std::array<int, 2>>> dependencyMap = taskDeps.getDependencyMap();

	// Dedicated A+ worker rank = last rank, only when A+ is enabled.
	// comm_farm excludes it: TaskStepManager/TaskStepWorker/dataCollect only
	// ever run across the manager + ordinary cohort workers. With -no-aplus
	// there is no dedicated rank at all, and comm_farm is just MPI_COMM_WORLD
	// (not owned - never split, never freed), reproducing the pre-A+ layout.
	const int aPlusWorkerRank = useAPlus ? size - 1 : -1;
	const bool isAPlusWorker  = useAPlus && (workerId == aPlusWorkerRank);

	MPI_Comm comm_farm;
	if (useAPlus) {
		MPI_Comm_split(MPI_COMM_WORLD, isAPlusWorker ? MPI_UNDEFINED : 0, workerId, &comm_farm);
	} else {
		comm_farm = MPI_COMM_WORLD;
	}

	int color = (workerId == 0) ? 0 : 1;
	MPI_Comm workerComm;
	MPI_Comm_split(MPI_COMM_WORLD, color, workerId, &workerComm);

	if (workerId == 0) {
		printf("[%d] Amount of data to be sent from workers to manager numData = %d numAgeGroups = %d numTimeSteps = %d numChunks = %d aPlus = %s\n", \
			workerId, numData, numAgeGroups, numTimeSteps, numChunks, useAPlus ? "on" : "off");
	}

	// dataCollect spans MPI_COMM_WORLD, not comm_farm: the dedicated A+
	// worker needs to publish into it too (see aplusChunkId()), and it sits
	// outside comm_farm entirely. TaskStepManager/TaskStepWorker still run
	// over comm_farm - only this data-sharing window includes everyone.
	DistDataCollector dataCollect(MPI_COMM_WORLD, numChunks, numData);

	if (!isAPlusWorker) {
		//
		// Farm ranks: manager (workerId 0) + ordinary cohort workers,
		// running over comm_farm (== MPI_COMM_WORLD when A+ is disabled).
		//
		if (workerId == 0) {
			//
			// Manager
			//
			double tik = MPI_Wtime();

			TaskStepManager manager(comm_farm, numCohorts, stepBegMap, stepEndMap, dependencyMap);

			// Sync the manager with the farm workers before starting to distribute the tasks
			MPI_Barrier(comm_farm);
			auto results = manager.run();

			double time_manager = MPI_Wtime() - tik;

			// Make sure the farm's own data are ready for the final checksum.
			MPI_Barrier(comm_farm);

			// The A+ worker publishes its density into dataCollect's A+ chunk
			// range as it goes (see runAPlusWorker()), but it is outside
			// comm_farm, so the barrier above doesn't cover it. Wait for its
			// explicit "done" signal instead - sent only after its very last
			// publish - before reading that range.
			if (useAPlus) {
				int done;
				MPI_Recv(&done, 1, MPI_INT, aPlusWorkerRank,
						 APLUS_CHECKSUM_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			double* data = dataCollect.getCollectedDataPtr();
			int numNormalChunks = numAgeGroups * numTimeSteps;
			double checksumNormal = std::accumulate(data, data + numNormalChunks * numData, 0.0);
			double checksumAPlus  = std::accumulate(data + numNormalChunks * numData, data + numChunks * numData, 0.0);
			double checksum = checksumNormal + checksumAPlus;
			printf("[%d] Checksum = %15.5lf (normal = %15.5lf, A+ = %15.8le) time manager = %10.5f sec\n",
				workerId, checksum, checksumNormal, checksumAPlus, time_manager);

			//free manager's singleton 'color'
			MPI_Comm_free(&workerComm);

		} else {
			//
			// Ordinary cohort worker
			//

			// Create SeapodymCohort object that will be shared among each worker
			int cmp_regime = 0;
			bool reset_buffers = false;
			//-----Memory stack sizes for dvariables and derivatives storage------
			gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
			long int gradstack_buffer, cmpdif_buffer, gs_var_buffer;
			bool grad_calc = false;
			if (cmp_regime==-1 || cmp_regime==2 || cmp_regime==4) grad_calc = true;
			buffers_init(gs_var_buffer, gradstack_buffer, cmpdif_buffer, grad_calc);
			if (reset_buffers)
				buffers_set(gs_var_buffer, gradstack_buffer, cmpdif_buffer);

			gradient_structure::set_GRADSTACK_BUFFER_SIZE(gradstack_buffer);
			gradient_structure::set_CMPDIF_BUFFER_SIZE(cmpdif_buffer);
			gradient_structure::set_NO_DERIVATIVES();

			// Does every worker need a gradiant structure object? Or does every cohort object need
			// its own gradient structure?
			gradient_structure gs(gs_var_buffer);

			{
				DataProvider dp(workerComm, nameSizePairs);

				SeapodymCohort cohort = xinit_prerun_wrapper(parfile.c_str(), useAPlus);

				cohort.setDataProvider(&dp);

				double t_shm = MPI_Wtime();
				cohort.setShmForcing();          // every node-local worker reads its timestep slice
				MPI_Barrier(dp.getShmComm());    // per-node publish-sync: all slabs visible before any read

				time_io_forcing += MPI_Wtime()-t_shm;
				// Bind the task function with the necessary parameters
				auto taskFunc = std::bind(taskFunction,
					std::placeholders::_1, // task_id
					std::placeholders::_2, // stepBeg
					std::placeholders::_3, // stepEnd
					std::placeholders::_4, // MPI communicator so we can send messages to the manager at the end of each step
					logger,
					&dataCollect,
					&cohort,
					numAgeGroups, numTimeSteps, aPlusWorkerRank);

				TaskStepWorker worker(comm_farm, taskFunc, stepBegMap, stepEndMap);

				// Sync the manager with the farm workers before starting to distribute the tasks
				MPI_Barrier(comm_farm);
				worker.run();
				MPI_Barrier(comm_farm);

				time_overhead = cohort.time_overhead;
			}
			//DataProvider clean-up
			MPI_Comm_free(&workerComm);
		}

	} else {
		//
		// Dedicated A+ worker (rank size-1, only reached when useAPlus).
		// Same forcing/shared-memory setup as an ordinary cohort worker
		// below (kept in sync with that block by inspection - both must
		// bring up an identical SeapodymCohort/DataProvider pair since A+
		// runs the exact same stepForward() dynamics), but it runs its own
		// independent loop (runAPlusWorker()) instead of TaskStepWorker.
		//
		int cmp_regime = 0;
		bool reset_buffers = false;
		gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
		long int gradstack_buffer, cmpdif_buffer, gs_var_buffer;
		bool grad_calc = false;
		if (cmp_regime==-1 || cmp_regime==2 || cmp_regime==4) grad_calc = true;
		buffers_init(gs_var_buffer, gradstack_buffer, cmpdif_buffer, grad_calc);
		if (reset_buffers)
			buffers_set(gs_var_buffer, gradstack_buffer, cmpdif_buffer);

		gradient_structure::set_GRADSTACK_BUFFER_SIZE(gradstack_buffer);
		gradient_structure::set_CMPDIF_BUFFER_SIZE(cmpdif_buffer);
		gradient_structure::set_NO_DERIVATIVES();

		gradient_structure gs(gs_var_buffer);

		{
			DataProvider dp(workerComm, nameSizePairs);

			SeapodymCohort cohort = xinit_prerun_wrapper(parfile.c_str(), useAPlus);

			cohort.setDataProvider(&dp);

			double t_shm = MPI_Wtime();
			cohort.setShmForcing();
			MPI_Barrier(dp.getShmComm());

			time_io_forcing += MPI_Wtime()-t_shm;

			runAPlusWorker(cohort, &dataCollect, numAgeGroups, numData, numTimeSteps, logger);

			time_overhead = cohort.time_overhead;
		}
		MPI_Comm_free(&workerComm);
	}

	if (workerId > 0) {
		printf("[%d] Timings calc/overhead/worker init/cohort init/comm: %10.3lf/%10.3lf/%10.3lf/%10.3lf/%10.3lf\n", workerId,
		time_calc, time_overhead, time_worker_init, time_cohort_init, time_mpi);

printf("[%d] InitCohort comm/flush/copy/init_restart,init_spawning/spawning/getdata/xreset (tot, per-init ms): "
         "%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f s ; %.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f ms over %ld inits\n", workerId,
         time_ic_comm, time_ic_flush, time_ic_copy, time_init_cohort_restart, time_init_cohort_spawning, time_spawning, time_getdata, time_xreset,
         n_ic? 1e3*time_ic_comm/n_ic:0, n_ic? 1e3*time_ic_flush/n_ic:0,
         n_ic? 1e3*time_ic_copy/n_ic:0, n_ic? 1e3*time_init_cohort_restart/numAgeGroups:0,
         n_ic? 1e3*time_init_cohort_spawning/n_ic:0, n_ic? 1e3*time_spawning/n_ic:0,
	 n_ic? 1e3*time_getdata/n_ic:0, n_ic? 1e3*time_xreset/n_ic:0,n_ic);

printf("[%d] Time IO/Put/Send/Idle: %.3f/%.3f/%.3f/%.3f ms\n",
         workerId, 1e3*time_io_forcing, 1e3*time_mpi_put, 1e3*time_send, 1e3*time_idle);
	}



	// Finalization of MPI
	////////////////////////////////////////////////////////////////////////
	dataCollect.free(); // collective on MPI_COMM_WORLD - every rank participates

	if (useAPlus && comm_farm != MPI_COMM_NULL)
		MPI_Comm_free(&comm_farm);

	MPI_Finalize();

	return 0;
}


