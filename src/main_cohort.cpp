#include <iostream>
#include <cstring>
#include <cstdlib>
#include <functional>
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

int getChunkId(int task_id, int step, int numAgeGroups) {
    int row = std::max(0, task_id - numAgeGroups + 1) + step;
    int col = task_id % numAgeGroups;
    return row * numAgeGroups + col;
}
void 
taskFunction(int task_id, int stepBeg, int stepEnd, MPI_Comm comm, 
    const char* parfile, int numData, 
    DistDataCollector* dataCollector,
    std::map<int, std::set<std::array<int, 2>>>* dependencyMap)
{

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

    int cohort_id = task_id; // In our case task_id is the cohort Id

    SeapodymCohort cohort((char*)parfile, cohort_id);

    //initialize variables of optimization
    const int nvar = cohort.nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1,nvar);

    cohort.xinit(x, x_names);
    //cout << "Total number of variables: " << nvar << '\n'<<'\n';

    //prepare cohort run
    cohort.prerun_model();

    //initialize cohort either from restart or from spawning
    cohort.init_cohort(x);

    int numAgeGroups = cohort.param->sp_nb_cohorts[0];

    std::vector<double> localData(numData);
    
    // advance the cohort 
    for (auto step = stepBeg; step < stepEnd; ++step) {
        /*// Fetch the data needed to create this cohort from the manager
        // and sum them up
        std::vector<double> initData(numData, 0.0);
        for (const auto& [task_id2, step] : (*dependencyMap)[task_id]) {
            int chunk_id = getChunkId(task_id2, step, numAgeGroups);
            std::vector<double> data = dataCollector->get(chunk_id);
            // check that the data are valid
            if (!data.empty() && data.back() == dataCollector->BAD_VALUE) {
                // The data have not been previously populated. This could indicate that
                // the worker has not yet produced any output for this cohort or the manager
                // has not yet received the data.
                MPI_Abort(comm, 1);
            }
            // sum up the cohort data at the previous time step
            std::transform(data.begin(), data.end(), initData.begin(), initData.begin(), std::plus<double>());
        }*/

        cohort.stepForward(false);

        // Send the data to the manager.
        std::vector<double> localData = cohort.GetCohortDensity();
        //std::fill(localData.begin(), localData.end(), 0.1);
        int chunk_id = getChunkId(task_id, step, numAgeGroups);
        dataCollector->put(chunk_id, localData.data());
        TTTRACE(cohort_id, step, chunk_id)

        int success = task_id;
        // send message to the manager that the step is complete
        int output[3] = {task_id, step, success};
        const int endTaskTag = 1;
        MPI_Send(output, 3, MPI_INT, 0, endTaskTag, comm);
    }

}

int main(int argc, char** argv) {

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
  
    // Read parfile
    VarParamCoupled param;
	param.init_param();
	param.read(parfile);

    // Get number of time steps and number of cohorts from param
    int numAgeGroups = param.sp_nb_cohorts[0];
    int Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps;
    Date::init_time_variables(param, Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps, 0,0);
    //int numTasks = numAgeGroups + numTimeSteps - 1;


    // Size of map (useful to access to a specific position adress of the 4D array pointer storing density)
    PMap map;
    map.lit_map(param);
    int numData = 0;
    const int imin = map.imin;
    const int imax = map.imax;
    for (int i = imin; i <= imax; i++){
        const int jmin = map.jinf[i];
        const int jmax = map.jsup[i];
        for (int j = jmin ; j <= jmax; j++){
            numData++;
        }
    }

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

    // print the dependencies for debugging
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

    // Bind the task function with the necessary parameters
    auto taskFunc = std::bind(taskFunction,
        std::placeholders::_1, // task_id
        std::placeholders::_2, // stepBeg
        std::placeholders::_3, // stepEnd
        std::placeholders::_4, // MPI communicator so we can send messages to the manager at the end of each step
        parfile.c_str(),
        numData,
        &dataCollect,
        &dependencyMap);

    if (workerId == 0) {
        // Manager
        TaskStepManager manager(MPI_COMM_WORLD, numCohorts, stepBegMap, stepEndMap, dependencyMap);
        auto results = manager.run();
        for (const auto& [task_id, step, res] : results) {
            std::cout << "Task ID " << task_id << " and step " << step << ": res = " << res << std::endl;	
        }
        std::cout << std::endl;
        dataCollect.displaySumChunk(12);
    } else {
        // Worker
        TaskStepWorker worker(MPI_COMM_WORLD, taskFunc, stepBegMap, stepEndMap);
        worker.run();
    }

    // Finalization of MPI
    ////////////////////////////////////////////////////////////////////////
    dataCollect.free();

    MPI_Finalize();

    return 0;
}



