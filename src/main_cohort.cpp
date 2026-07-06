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

double time_ic_comm = 0.0, time_ic_flush = 0.0, time_ic_copy = 0.0, time_spawning = 0.0, time_getdata = 0.0, time_xreset = 0.0, time_init_cohort_spawning = 0.0, time_init_cohort_restart = 0.0, time_io_forcing = 0.0;
long   n_ic = 0;   // count of spawning-path inits, for per-init averages
double time_mpi_put = 0.0, time_send = 0.0, time_idle = 0.0;
	
SeapodymCohort* seapodym_cohort(const char* parfile, const int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs);
void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
void buffers_set(long int &mv, long int &mc, long int &mg);


double time_worker_init = 0.0, time_cohort_init = 0.0, time_calc = 0.0, time_mpi = 0.0, time_step = 0.0, time_overhead = 0.0;

SeapodymCohort xinit_prerun_wrapper(const char* parfile) {

	double tik = MPI_Wtime();

	SeapodymCohort cohort((char*)parfile, 0);

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


void taskFunction(int task_id, int stepBeg, int stepEnd, MPI_Comm comm,
		const std::shared_ptr<spdlog::logger>& logger,
		DistDataCollector* dataCollector,
		SeapodymCohort* cohort,
		int firstAPlusId, int numAgeGroups, int numTimeSteps){

	static double last_task_end = -1.0;        // per-worker process, persists across calls
	double t_in = MPI_Wtime();
	if (last_task_end >= 0.0)
		time_idle += t_in - last_task_end;     // <-- time spent in worker.run() waiting for dispatch

	double tik = MPI_Wtime();

	if (task_id >= firstAPlusId) {
		// A+ (plus group) accumulator task. This is not a live biological
		// cohort: it is the spatial sum of (a) the individuals that just
		// graduated into the top age class at the previous time step and
		// (b) the existing A+ pool carried over from the previous time step.
		// stepBeg/stepEnd are always 0/1 for these tasks (see
		// SeapodymCohortDependencyAnalyzer), so there is nothing to loop over.
		int t = task_id - firstAPlusId;
		logger->info("> A+ task id {} (t={})", task_id, t);

		std::vector<double> buf(dataCollector->getNumSize());
		if (task_id == firstAPlusId) {
			// t=0: no upstream dependency; seed from the actual initial
			// condition for the A+ age bin (index sp_nb_cohorts[0]-1 == numAgeGroups).
			buf = cohort->GetInitDensity(numAgeGroups);
		} else {
			int prevAPlusChunk  = SeapodymCohort::computeAPlusChunkId(task_id - 1, firstAPlusId, numAgeGroups, numTimeSteps);
			int graduatingChunk = SeapodymCohort::computeChunkId(t - 1, numAgeGroups - 1, numAgeGroups);
			std::vector<double> prevAPlus(dataCollector->getNumSize());
			std::vector<double> graduating(dataCollector->getNumSize());
			dataCollector->get(prevAPlusChunk, prevAPlus.data());
			dataCollector->get(graduatingChunk, graduating.data());
			for (std::size_t k = 0; k < buf.size(); ++k)
				buf[k] = prevAPlus[k] + graduating[k];
		}

		int myChunk = SeapodymCohort::computeAPlusChunkId(task_id, firstAPlusId, numAgeGroups, numTimeSteps);
		dataCollector->put(myChunk, buf.data());

		int success = task_id;
		int output[3] = {task_id, stepBeg, success};
		MPI_Send(output, 3, MPI_INT, 0, END_TASK_TAG, comm);

		time_calc += MPI_Wtime() - tik;
		last_task_end = MPI_Wtime();
		logger->info("< A+ task id {} (t={})", task_id, t);
		return;
	}

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

	// Parse the command line arguments
	bool success = cmdLine.parse(argc, argv);
	bool help = cmdLine.get<bool>("-help") || cmdLine.get<bool>("-h");
	if (!success) {
		std::cerr << "Error parsing command line arguments." << std::endl;
		cmdLine.help();
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (help) {
		cmdLine.help();
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
	// (age index sp_nb_cohorts[0]-1) is modelled separately as its own chain
	// of one-step accumulator tasks (see SeapodymCohortDependencyAnalyzer).
	int numAgeGroups = param.sp_nb_cohorts[0] - 1;
	int Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps;
	Date::init_time_variables(param, Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps, 0,0);
	//int numTasks = numAgeGroups + numTimeSteps - 1;

	// Size of map (useful to access to a specific position adress of the 4D array pointer storing density)
	map.lit_map(param);
	int numData = map.get_state_array_size();

	//Set-up the size for the shared arrays for forcing data
	std::vector<std::pair<std::string, std::size_t>> nameSizePairs = param.getDpNameSizePairs(numTimeSteps, map.get_array_size());

	// set up the data collector. One extra chunk per time step is reserved
	// for the A+ (plus group) accumulator series, appended after the normal
	// (task, step) chunk range - see SeapodymCohort::computeAPlusChunkId().
	int numChunks = numAgeGroups * numTimeSteps + numTimeSteps;

	int color = (workerId == 0) ? 0 : 1;
	MPI_Comm workerComm;
	MPI_Comm_split(MPI_COMM_WORLD, color, workerId, &workerComm);

	if (workerId == 0) {
		printf("[%d] Amount of data to be sent from workers to manager numData = %d numAgeGroups = %d numTimeSteps = %d numChunks = %d\n", \
			workerId, numData, numAgeGroups, numTimeSteps, numChunks);
	}

	DistDataCollector dataCollect(MPI_COMM_WORLD, numChunks, numData);

	// analyze the cohort Id task dependencies (aPlusCohort=true adds the A+ chain)
	SeapodymCohortDependencyAnalyzer taskDeps(numAgeGroups, numTimeSteps, param.age_mature[0], /*aPlusCohort=*/true);
	int firstAPlusId = taskDeps.getFirstAPlusCohortId();
	int numCohorts = taskDeps.getNumberOfCohorts();
	std::map<int, int> stepBegMap = taskDeps.getStepBegMap();
	std::map<int, int> stepEndMap = taskDeps.getStepEndMap();
	std::map<int, std::set<std::array<int, 2>>> dependencyMap = taskDeps.getDependencyMap();

	if (workerId == 0) {
		//
		// Manager
		//
		double tik = MPI_Wtime();

		TaskStepManager manager(MPI_COMM_WORLD, numCohorts, stepBegMap, stepEndMap, dependencyMap);
		// Sync the manager with the workers before starting to distribute the tasks
		MPI_Barrier(MPI_COMM_WORLD);
		auto results = manager.run();

		double time_manager = MPI_Wtime() - tik;

		// Make sure the data are ready for the final checksum
		MPI_Barrier(MPI_COMM_WORLD);
		double* data = dataCollect.getCollectedDataPtr();
		// Split the checksum into the "normal" cohort chunk range and the A+
		// chunk range (see SeapodymCohort::computeAPlusChunkId): normal chunks
		// come first, A+ chunks are appended after. Keeping them separate lets
		// us confirm a future change to A+'s dynamics doesn't perturb normal
		// cohort results, and vice versa, instead of relying on one aggregate
		// number that could mask a regression in either half.
		int numNormalChunks = numAgeGroups * numTimeSteps;
		double checksumNormal = std::accumulate(data, data + numNormalChunks * numData, 0.0);
		double checksumAPlus  = std::accumulate(data + numNormalChunks * numData, data + numChunks * numData, 0.0);
		double checksum = checksumNormal + checksumAPlus;
		printf("[%d] Checksum = %15.5lf (normal = %15.5lf, A+ = %15.5lf) time manager = %10.5f sec\n",
			workerId, checksum, checksumNormal, checksumAPlus, time_manager);

		//free manager's singleton 'color'
		MPI_Comm_free(&workerComm);
	} else {
		//
		// Worker
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

			SeapodymCohort cohort= xinit_prerun_wrapper(parfile.c_str());

			cohort.setDataProvider(&dp);

			//MPI_Win_fence(0, dp.win());
			//if (dp.isShmRoot())
			//	cohort.setShmForcing();//needs to be in cohort where the reading is done, but will be done once
			//MPI_Win_fence(0, dp.win());
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
				firstAPlusId, numAgeGroups, numTimeSteps);

			TaskStepWorker worker(MPI_COMM_WORLD, taskFunc, stepBegMap, stepEndMap);

			// Sync the manager with the workers before starting to distribute the tasks
			MPI_Barrier(MPI_COMM_WORLD);
			worker.run();
			MPI_Barrier(MPI_COMM_WORLD);

			time_overhead = cohort.time_overhead;
		}
		//DataProvider clean-up
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
	dataCollect.free();
	MPI_Finalize();

	return 0;
}



