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

int taskFunction(int task_id, const char* parfile) {

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
	SeapodymCohort* scp = seapodym_cohort(parfile, cmp_regime, reset_buffers, cohort_id, gs);

	scp->prerun_model();
	scp->OnRunFirstStep();

	delete scp;

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
	auto taskFunc = std::bind(taskFunction, std::placeholders::_1, parfile.c_str());
	int numTasks = cmdLine.get<int>("-nT");

	if (workerId == 0) {
		// Manager
		TaskManager manager(MPI_COMM_WORLD, numTasks);
		std::vector<int> task_ids = manager.run();
		for (auto task_id : task_ids) {
			std::cout << task_id << ", ";	
		}
		std::cout << std::endl;

	} else {
		// Worker

		TaskWorker worker(MPI_COMM_WORLD, taskFunc);
		worker.run();
	}

	// Finalization of MPI
	////////////////////////////////////////////////////////////////////////
	err = MPI_Finalize();

	return 0;
}



