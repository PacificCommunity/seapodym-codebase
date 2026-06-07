#include <iostream>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <numeric>      // std::accumulate
#include <mpi.h>
#include <DataProvider.h>
#include "VarParamCoupled.h"
#include "SeapodymCohortDependencyAnalyzer.h"
#include "TaskStepWorker.h"
#include "DistDataCollector.h"
#include "TaskStepManager.h"
#include "SeapodymCohort.h"
#include <CmdLineArgParser.h>
#include "ctrace.h"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

SeapodymCohort* seapodym_cohort(const char* parfile, const int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs);
void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
void buffers_set(long int &mv, long int &mc, long int &mg);


double time_worker_init = 0.0, time_cohort_init = 0.0, time_calc = 0.0, time_mpi = 0.0, time_step = 0.0, time_overhead = 0.0;

SeapodymCohort xinit_prerun_wrapper(const char* parfile, DataProvider* dp = nullptr) {

    double tik = MPI_Wtime();

    SeapodymCohort cohort((char*)parfile, 0);

    //initialize variables of optimization
    const int nvar = cohort.nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1,nvar);

    cohort.xinit(x, x_names);

    // Attach DataProvider before prerun_model so it is available during ReadAll (Step 2+)
    cohort.setDataProvider(dp);

    //prepare cohort run
    cohort.prerun_model();

    double tak = MPI_Wtime();

    time_worker_init += tak - tik;

    return cohort;
}


void 
taskFunction(int task_id, int stepBeg, int stepEnd, MPI_Comm comm,
    const std::shared_ptr<spdlog::logger>& logger,
    DistDataCollector* dataCollector,
    SeapodymCohort* cohort)
{
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

        // Send the data to the manager.
        logger->info("        >>> send data for step {} of task id {}", step, task_id);
        std::vector<double> localData = cohort->GetCohortDensity();
        int chunk_id = cohort->getChunkId(step);

        double tik_mpi = MPI_Wtime();
        dataCollector->put(chunk_id, localData.data());
        time_mpi += MPI_Wtime() - tik_mpi;
        logger->info("        <<< send data for step {} of task id {}", step, task_id);

        int success = task_id;
        // send message to the manager that the step is complete
        int output[3] = {task_id, step, success};
        const int endTaskTag = 1;

        logger->info("        >>> notify manager after step {} of task id {}", step, task_id);
        tik_mpi = MPI_Wtime();
        MPI_Send(output, 3, MPI_INT, 0, endTaskTag, comm);
        time_mpi += MPI_Wtime() - tik_mpi;
        logger->info("        <<< notify manager after step {} of task id {}", step, task_id);
    }

    time_calc += MPI_Wtime() - tak;
    logger->info("< task id {} for steps {} to {}", task_id, stepBeg, stepEnd);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

    // MPI initialization
    MPI_Init(&argc, &argv);

    time_worker_init = 0;	
    time_calc = 0;	

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int workerId;
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

    // Get number of time steps and number of cohorts from param
    int numAgeGroups = param.sp_nb_cohorts[0];
    int Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps;
    Date::init_time_variables(param, Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps, 0,0);
    //int numTasks = numAgeGroups + numTimeSteps - 1;

    // Size of map (useful to access to a specific position adress of the 4D array pointer storing density)
    map.lit_map(param);
    int numData = map.get_state_array_size();

    // ------------------------------------------------------------------ //
    // Build a worker-only communicator so that shmRank-0 is the lowest
    // worker rank, not the manager (rank 0 in MPI_COMM_WORLD).
    // MPI_Comm_split with MPI_UNDEFINED excludes rank 0 from workerComm.
    // ------------------------------------------------------------------ //
    MPI_Comm workerComm = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD,
                   workerId == 0 ? MPI_UNDEFINED : 0,
                   workerId,
                   &workerComm);

    // Enclose DataProvider (and DistDataCollector) in a scope so their
    // destructors run before MPI_Finalize().
    {
    // set up the data collector (uses MPI_COMM_WORLD — manager + workers)
    int numChunks = numAgeGroups * numTimeSteps;

    if (workerId == 0) {
        printf("[%d] Amount of data to be sent from workers to manager numData = %d numAgeGroups = %d numTimeSteps = %d numChunks = %d\n", \
		    workerId, numData, numAgeGroups, numTimeSteps, numChunks);
    }

    DistDataCollector dataCollect(MPI_COMM_WORLD, numChunks, numData);

    // analyze the cohort Id task dependencies
    SeapodymCohortDependencyAnalyzer taskDeps(numAgeGroups, numTimeSteps);
    int numCohorts = taskDeps.getNumberOfCohorts();
    std::map<int, int> stepBegMap = taskDeps.getStepBegMap();
    std::map<int, int> stepEndMap = taskDeps.getStepEndMap();
    std::map<int, std::set<std::array<int, 2>>> dependencyMap = taskDeps.getDependencyMap();

    if (workerId == 0) {
        //
        // Manager: no DataProvider needed (manager never calls ReadAll)
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
        // print check sum
        double checksum = std::accumulate(data, data + numChunks * numData, 0.0);
        printf("[%d] Checksum = %15.5lf time manager = %10.5f sec\n", workerId, checksum, time_manager);
    } else {
        //
        // Worker: construct DataProvider on the worker-only communicator so
        // shmRank-0 is the lowest-ranked worker, not the manager.
        //
        const std::size_t mapCells = static_cast<std::size_t>(map.get_array_size());
        const std::size_t dpT  = static_cast<std::size_t>(numTimeSteps);
        const std::size_t dpNL = static_cast<std::size_t>(param.nb_layer);
        const std::size_t dpNF = static_cast<std::size_t>(param.get_nbforage());
        const std::size_t sz3  = dpT * mapCells;
        const std::size_t sz4  = dpT * dpNL * mapCells;
        const std::size_t szF  = dpT * dpNF * mapCells;

        DataProvider dataProvider(workerComm, {
            {"np1",    sz3},
            {"sst",    sz3},
            {"vld",    sz3},
            {"un",     sz4},
            {"vn",     sz4},
            {"tempn",  sz4},
            {"oxygen", sz4},
            {"forage", szF},
        });

        if (workerId == 1) {  // rank 1 = shmRank 0 = shmRoot among workers
            printf("[%d] DataProvider: mapCells=%zu T=%zu NL=%zu NF=%zu  "
                   "(np1/sst/vld=%zu  un/vn/tempn/oxygen=%zu  forage=%zu doubles)\n",
                   workerId, mapCells, dpT, dpNL, dpNF, sz3, sz4, szF);
        }

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

        SeapodymCohort cohort = xinit_prerun_wrapper(parfile.c_str(), &dataProvider);

        // Bind the task function with the necessary parameters
        auto taskFunc = std::bind(taskFunction,
            std::placeholders::_1, // task_id
            std::placeholders::_2, // stepBeg
            std::placeholders::_3, // stepEnd
            std::placeholders::_4, // MPI communicator so we can send messages to the manager at the end of each step
            logger,
            &dataCollect,
            &cohort);

        TaskStepWorker worker(MPI_COMM_WORLD, taskFunc, stepBegMap, stepEndMap);

        // Sync the manager with the workers before starting to distribute the tasks
        MPI_Barrier(MPI_COMM_WORLD);
        worker.run();
        MPI_Barrier(MPI_COMM_WORLD);

	    time_overhead = cohort.time_overhead;
    }

    if (workerId > 0) {
        printf("[%d] Timings calc/overhead/worker init/cohort init/comm: %10.3lf/%10.3lf/%10.3lf/%10.3lf/%10.3lf\n", workerId, 
        time_calc, time_overhead, time_worker_init, time_cohort_init, time_mpi);
    }

    // Finalization of MPI
    ////////////////////////////////////////////////////////////////////////
    dataCollect.free();
    } // DataProvider (worker-only) and DistDataCollector destroyed here, before MPI_Finalize

    if (workerComm != MPI_COMM_NULL)
        MPI_Comm_free(&workerComm);

    MPI_Finalize();

    return 0;
}



