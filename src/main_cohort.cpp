#include <iostream>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <mpi.h>
//#include "SeapodymCohort.h"
#include "VarParamCoupled.h"
#include "SeapodymCohortDependencyAnalyzer.h"
#include "TaskStepWorker.h"
#include "TaskStepManager.h"
#include "SeapodymCohort.h"
#include <CmdLineArgParser.h>
#include "ctrace.h"

using std::cout;
SeapodymCohort* seapodym_cohort(const char* parfile, const int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs);
void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
void buffers_set(long int &mv, long int &mc, long int &mg);

int taskFunction(int task_id, int step, int stepBeg, int stepEnd, const char* parfile) {

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

    int cohort_id = task_id; // For the time being

    static SeapodymCohort *cohort = nullptr;
    if (step==stepBeg){
        cohort = new SeapodymCohort((char*)parfile, cohort_id);

        //initialize variables of optimization
        const int nvar = cohort->nvarcalc();
        independent_variables x(1, nvar);
        adstring_array x_names(1,nvar);

        cohort->xinit(x, x_names);
        //cout << "Total number of variables: " << nvar << '\n'<<'\n';

        //initialization of simulation
        cohort->prerun_model(x);
    }

    // advance the cohort by one step
    // NOT SURE IF ALL THE cohorts share the same param? Should param be passed as an argument to the taskFunction? Or should it be computed 
    cohort->stepForward(false);

    if (step == stepEnd) {
        // remove the cohort
        delete cohort;
    }

    // Could return an error code instead
    return task_id;
}

int main(int argc, char** argv) {

    // Initialization of MPI
    int err;
    err = MPI_Init(&argc, &argv);
    int workerId = 0;
    err = MPI_Comm_rank(MPI_COMM_WORLD, &workerId);
    int size = 1;
    err = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2) {
        std::cerr << "ERROR: must have at least 2 ranks\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    CmdLineArgParser cmdLine;
    cmdLine.set("-s", std::string("initparfile.xml"), "Input parameter file");
    cmdLine.set("-ns", 5, "Number of steps for each task");

    // NO NEED TO HAVE -na, its in the parfile. However, we need that quantity before we 
    // intantiate the cohort objects.
    cmdLine.set("-nT", 1, "Number of tasks");
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
    auto taskFunc = std::bind(taskFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, parfile.c_str());
    
    // Read parfile
    VarParamCoupled* param;
    param = new VarParamCoupled();
	param->init_param();
	param->read(parfile);

    // Get number of time steps and number of cohorts from param
    int numAgeGroups = param->sp_nb_cohorts[0];
    int Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps;
    Date::init_time_variables(*param, Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, numTimeSteps, 0,0);
    int numTasks = numAgeGroups + numTimeSteps - 1;

    // analyze the conhort Id task dependencies
    SeapodymCohortDependencyAnalyzer taskDeps(numAgeGroups, numTimeSteps);
    int numCohorts = taskDeps.getNumberOfCohorts();
    int numCohortSteps = taskDeps.getNumberOfCohortSteps();
    std::map<int, int> stepBegMap = taskDeps.getStepBegMap();
    std::map<int, int> stepEndMap = taskDeps.getStepEndMap();
    std::map<int, std::set<std::array<int, 2>>> dependencyMap = taskDeps.getDependencyMap();

    // print the dependencies for debugging
    if (workerId == 0) {
        cout << "Number of age groups: " << numAgeGroups << endl;
        cout << "Number of tasks: " << numTasks << "; simulation time: " << numTimeSteps << endl << endl;
        for (const auto& [task_id, stepBeg] : stepBegMap) {
            int globalTimeIndex = std::max(0, task_id - numAgeGroups + 1);
            std::cout << "At time " << globalTimeIndex << " Task " << task_id << " has steps " << stepBeg << "..." << stepEndMap.at(task_id) - 1 << " and depends on: ";
            for (const auto& [task_id2, step] : dependencyMap.at(task_id)) {
                std::cout << task_id2 << ":" << step << ", ";
            }
            std::cout << std::endl;
        }
    }   

    if (workerId == 0) {
        // Manager
        TaskStepManager manager(MPI_COMM_WORLD, numTasks, stepBegMap, stepEndMap, dependencyMap);
        auto results = manager.run();
        for (const auto& [task_id, step, res] : results) {
            std::cout << "Task ID " << task_id << " and step " << step << ": res = " << res << std::endl;	
        }
        std::cout << std::endl;

    } else {
        // Worker

        TaskStepWorker worker(MPI_COMM_WORLD, taskFunc, stepBegMap, stepEndMap);
        worker.run();
    }

    // Finalization of MPI
    ////////////////////////////////////////////////////////////////////////
    err = MPI_Finalize();

    return 0;
}



