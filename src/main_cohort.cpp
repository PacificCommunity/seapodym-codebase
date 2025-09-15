#include <iostream>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <numeric>      // std::accumulate
#include <mpi.h>
#include "VarParamCoupled.h"
#include "SeapodymCohortDependencyAnalyzer.h"
#include "TaskStepWorker.h"
#include "DistDataCollector.h"
#include "TaskStepManager.h"
#include "SeapodymCohort.h"
#include <CmdLineArgParser.h>
#include "ctrace.h"

using std::cout;
SeapodymCohort* seapodym_cohort(const char* parfile, const int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs);
void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
void buffers_set(long int &mv, long int &mc, long int &mg);

double tik, tak, time_init, time_calc;

SeapodymCohort xinit_prerun_wrapper(const char* parfile){

    tik = MPI_Wtime();
    
    SeapodymCohort cohort((char*)parfile, 0);

    //initialize variables of optimization
    const int nvar = cohort.nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1,nvar);

    cohort.xinit(x, x_names);
    //cout << "Total number of variables: " << nvar << '\n'<<'\n';

    //prepare cohort run
    cohort.prerun_model();

    tak = MPI_Wtime();

    time_init += tak - tik;

    return cohort;
}


void 
taskFunction(int task_id, int stepBeg, int stepEnd, MPI_Comm comm, 
    const char* parfile, int numData, 
    DistDataCollector* dataCollector,
    SeapodymCohort* cohort)
{
    //initialize variables of optimization
    const int nvar = cohort->nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1,nvar);

    int cohort_id = task_id;
    cohort->restart(cohort_id);
    //initialize cohort either from restart or from spawning
    tik = MPI_Wtime();
    cohort->init_cohort(x,*dataCollector);

    std::vector<double> localData(numData);
    
    // advance the cohort 
    for (auto step = stepBeg; step < stepEnd; ++step) {

        cohort->stepForward(false);

        // Send the data to the manager.
        std::vector<double> localData = cohort->GetCohortDensity();
        int chunk_id = cohort->getChunkId(step);
        dataCollector->put(chunk_id, localData.data());

        int success = task_id;
        // send message to the manager that the step is complete
        int output[3] = {task_id, step, success};
        const int endTaskTag = 1;
        MPI_Send(output, 3, MPI_INT, 0, endTaskTag, comm);
    }

    tak = MPI_Wtime();
    time_calc += tak - tik;

}

int main(int argc, char** argv) {

    time_init = 0;	
    time_calc = 0;	
    // MPI initialization
    MPI_Init(&argc, &argv);
    int numWorkers, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    numWorkers = size - 1;
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

    // set up the data collector
    int numChunks = numAgeGroups * numTimeSteps;
    DistDataCollector dataCollect(MPI_COMM_WORLD, numChunks, numData);
    
    // analyze the cohort Id task dependencies
    SeapodymCohortDependencyAnalyzer taskDeps(numAgeGroups, numTimeSteps);
    int numCohorts = taskDeps.getNumberOfCohorts();
    int numCohortSteps = taskDeps.getNumberOfCohortSteps();
    std::map<int, int> stepBegMap = taskDeps.getStepBegMap();
    std::map<int, int> stepEndMap = taskDeps.getStepEndMap();
    std::map<int, std::set<std::array<int, 2>>> dependencyMap = taskDeps.getDependencyMap();

/*    // print the dependencies for debugging
    if (workerId == 0) {
        for (const auto& [task_id, stepBeg] : stepBegMap) {
            int globalTimeIndex = std::max(0, task_id - numAgeGroups + 1);
            std::cout << "At time " << globalTimeIndex << " Task " << task_id << " has steps " << stepBeg << "..." << stepEndMap.at(task_id) - 1 << " and depends on: ";
            for (const auto& [task_id2, step] : dependencyMap.at(task_id)) {
                std::cout << task_id2 << ":" << step << ", ";
            }
            std::cout << std::endl;
        }
    }
*/


    if (workerId == 0) {
        // Manager
        TaskStepManager manager(MPI_COMM_WORLD, numCohorts, stepBegMap, stepEndMap, dependencyMap);
        auto results = manager.run();
        // for (const auto& [task_id, step, res] : results) {
        //     std::cout << "Task ID " << task_id << " and step " << step << ": res = " << res << std::endl;	
        // }
        // std::cout << std::endl;
//        dataCollect.displaySumChunk(5);
        double* data = dataCollect.getCollectedDataPtr();
        // print check sum
        double checksum = std::accumulate(data, data + numChunks * numData, 0.0);
        std::cout << "Checksum = " << checksum << std::endl;
    } else {
        // Worker

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
        // Des every worker need a gradiant structure object? Or does every cohort object need
        // its own gradient structure?
        gradient_structure gs(gs_var_buffer);        
        
        SeapodymCohort cohort= xinit_prerun_wrapper(parfile.c_str());

        // Bind the task function with the necessary parameters
        auto taskFunc = std::bind(taskFunction,
            std::placeholders::_1, // task_id
            std::placeholders::_2, // stepBeg
            std::placeholders::_3, // stepEnd
            std::placeholders::_4, // MPI communicator so we can send messages to the manager at the end of each step
            parfile.c_str(),
            numData,
            &dataCollect,
            &cohort);

        TaskStepWorker worker(MPI_COMM_WORLD, taskFunc, stepBegMap, stepEndMap);
        worker.run();
    }

    // Finalization of MPI
    ////////////////////////////////////////////////////////////////////////
    dataCollect.free();
    MPI_Finalize();

    TTTRACE(time_init,time_calc,time_calc/(time_calc+time_init))

    return 0;
}



