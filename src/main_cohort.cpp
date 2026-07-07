#include <iostream>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <numeric>	  // std::accumulate
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
	cohort->init_cohort(x,*dataCollector);
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
		// ------------------------------------------------------------------
		if (aPlusWorkerRank >= 0 && step == numAgeGroups - 1 && task_id <= numTimeSteps - 2) {
			logger->info("        >>> A+ notify for task id {}", task_id);
			double t_n = MPI_Wtime();
			MPI_Send(localData.data(), (int)localData.size(), MPI_DOUBLE,
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

///////////////////////////////////////////////////////////////////////////////

// Dedicated A+ (plus group) worker loop. Runs on its own MPI rank, entirely
// outside the TaskStepManager/TaskStepWorker dependency graph - ordering is
// enforced by the ping-pong in taskFunction() above (a feeder blocks on the
// ACK before telling the manager it is done), not by an explicit dependency
// edge.
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
void runAPlusWorker(SeapodymCohort& cohort, int numData, int numTimeSteps,
		const std::shared_ptr<spdlog::logger>& logger) {

	const int nvar = cohort.nvarcalc();
	independent_variables x(1, nvar);
	adstring_array x_names(1, nvar);
	cohort.xinit(x, x_names);

	double checksumAPlus = 0.0;

	// t=0: no upstream dependency - the opening balance is just the initial
	// condition for the A+ age bin, read from file like any other initial
	// cohort.
	cohort.restartAPlus(0);
	cohort.init_cohort_aplus(x, std::vector<double>(), /*seedFromFile=*/true);
	cohort.stepForward(false);
	checksumAPlus += cohort.Checksum();
	logger->info("A+ t=0 done. Checksum = {}", cohort.Checksum());

	std::vector<double> graduating(numData);
	for (int t = 1; t < numTimeSteps; ++t) {

		MPI_Status status;
		logger->info(">>> A+ waiting for feeder, t={}", t);
		MPI_Recv(graduating.data(), numData, MPI_DOUBLE,
				 MPI_ANY_SOURCE, APLUS_NOTIFY_TAG, MPI_COMM_WORLD, &status);
		logger->info("<<< A+ received feeder data, t={}", t);

		// Merge this step's graduating cohort with A+'s own current density.
		// No RMA fetch needed for "the previous A+ pool" - this object's own
		// state, right after the previous stepForward(), already holds it.
		std::vector<double> prevAPlus = cohort.GetCohortDensity();
		std::vector<double> merged(prevAPlus.size());
		for (std::size_t k = 0; k < merged.size(); ++k)
			merged[k] = prevAPlus[k] + graduating[k];

		cohort.restartAPlus(t);
		cohort.init_cohort_aplus(x, merged, /*seedFromFile=*/false);
		cohort.stepForward(false);
		checksumAPlus += cohort.Checksum();
		logger->info("A+ t={} done. Checksum = {}", t, cohort.Checksum());

		// Acknowledge the feeder. Only now can it notify the manager, so any
		// cohort scheduled after that point is guaranteed to see this step's
		// contribution already folded into A+.
		int ack = 1;
		MPI_Send(&ack, 1, MPI_INT, status.MPI_SOURCE, APLUS_ACK_TAG, MPI_COMM_WORLD);
	}

	// Hand the final checksum to the manager so it can print the unified
	// normal/A+ report.
	MPI_Send(&checksumAPlus, 1, MPI_DOUBLE, 0, APLUS_CHECKSUM_TAG, MPI_COMM_WORLD);
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

	// Set up the data collector for normal cohorts only - A+ no longer
	// shares this buffer (see runAPlusWorker()/taskFunction()), so its size
	// no longer needs an extra per-time-step chunk range.
	int numChunks = numAgeGroups * numTimeSteps;

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

	if (!isAPlusWorker) {
		//
		// Farm ranks: manager (workerId 0) + ordinary cohort workers,
		// running over comm_farm (== MPI_COMM_WORLD when A+ is disabled).
		//
		DistDataCollector dataCollect(comm_farm, numChunks, numData);

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

			// Make sure the data are ready for the final checksum
			MPI_Barrier(comm_farm);
			double* data = dataCollect.getCollectedDataPtr();
			double checksumNormal = std::accumulate(data, data + numChunks * numData, 0.0);

			// The A+ worker keeps its density in its own local memory (it is
			// never written through dataCollect - see runAPlusWorker()), so
			// it computes its own checksum and hands it over directly, only
			// once its last step is done.
			double checksumAPlus = 0.0;
			if (useAPlus) {
				MPI_Recv(&checksumAPlus, 1, MPI_DOUBLE, aPlusWorkerRank,
						 APLUS_CHECKSUM_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
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

		dataCollect.free(); // collective on comm_farm

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

			runAPlusWorker(cohort, numData, numTimeSteps, logger);

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
	if (useAPlus && comm_farm != MPI_COMM_NULL)
		MPI_Comm_free(&comm_farm);

	MPI_Finalize();

	return 0;
}


