#include <iostream>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <mpi.h>
//#include "SeapodymCohort.h"
#include "TaskWorker.h"
#include "TaskManager.h"
#include "SeapodymCohort.h"
#include <CmdLineArgParser.h>

using std::cout;
SeapodymCohort* seapodym_cohort(const char* parfile, const int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs);
void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
void buffers_set(long int &mv, long int &mc, long int &mg);

int taskFunction(int cohortId, int step,  int numSteps, VarParamCoupled* param, int init_from_inputfile, DVAR4_ARRAY* array_ptr) {
    //array_ptr would point toward a 4D array [time; cohort_id, lon, lat]

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

    SeapodymCohort *cohort = null;
    if (step == 0) {
        cohort = new SeapodymCohort(param, cohortId);

        //initialize variables of optimization
        const int nvar = cohort.nvarcalc();
        independent_variables x(1, nvar);
        adstring_array x_names(1,nvar);

        cohort.xinit(x, x_names);
        //cout << "Total number of variables: " << nvar << '\n'<<'\n';

        //initialization of simulation
        cohort.prerun_model(x, init_from_inputfile, array_ptr);
    }

    // advance the cohort by one step
    // NOT SURE IF ALL THE cohorts share the same param? Should param be passed as an argument to the taskFunction? Or should it be computed 
    cohort.stepForward(array_ptr);

    if (step == numSteps - 1) {
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
    int numTasks = cmdLine.get<int>("-nT");

    // Read parfile once, then pass it as an argument of taskFunc
    VarParamCoupled* param;
    param = new VarParamCoupled();
	param->init_param();
	EditRunCoupled(parfile);

    // Create a 4D array to store cohort density across time: [time, age, lon, lat] (Actually 5D: [time, sp, age, lon, lat] but only one species for now)
    int nb_species = 1;
    SeapodymCohort bogus_sc = new SeapodymCohort(parfile, 0);// bogus cohort to get map object
    bogus_sc.prerun_model();
    DVAR4_ARRAY cohorts_densities;
    cohorts_densities.allocate(0, param->nbt_total - 1);
    cohorts_densities.initiliaze();
    for (int t=0; t<param->nbt_total, t++){
        cohorts_densities[t].allocate(0, nb_species - 1);
        cohorts_densities[t].initiliaze();
        for (int sp=0; sp<nb_species, sp++){
            cohorts_densities[t][sp].allocate(0, bogus_sc.nb_age_class - 1);
            for (int a=0; a<bogus_sc.nb_age_class; a++){
                cohorts_densities[t][sp][a].allocate(bogus_sc.map.imin1, bogus_sc.map.imax1, bogus_sc.map.jinf1, bogus_sc.map.jsup1);
                cohorts_densities[t][sp][a].initialize();
            }
        }
    }

    if (workerId == 0) {
        // Manager
        TaskManager manager(MPI_COMM_WORLD, numTasks);
        std::map<int, int> results = manager.run();
        for (auto [task_id, result] : results) {
            std::cout << task_id << ": " << result << ", ";	
        }
        std::cout << std::endl;

    } else {
        // Worker
        TaskStepWorker worker(MPI_COMM_WORLD, taskFunction);// How to pass arguments to taskFunction ?
        worker.run();
    }

    // Finalization of MPI
    ////////////////////////////////////////////////////////////////////////
    err = MPI_Finalize();

    return 0;
}



