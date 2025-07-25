#include <iostream>
#include <cstring>
#include <cstdlib>
#include <mpi.h>
#include "SeapodymCohort.h"
#include <CmdLineArgParser.h>

using std::cout;
SeapodymCohort* seapodym_cohort(const char* parfile, const int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs);
void buffers_init(long int &mv, long int &mc, long int &mg, const bool grad_calc);
void buffers_set(long int &mv, long int &mc, long int &mg);

int main(int argc, char** argv) {

	// Initialization of MPI
	int err;
	err = MPI_Init(&argc, &argv);
	int workerId = 0;
	err = MPI_Comm_rank(MPI_COMM_WORLD, &workerId);
	int num_workers = 1;
	err = MPI_Comm_size(MPI_COMM_WORLD, &num_workers);
	
	
	int cmp_regime = 0;
	bool reset_buffers = false;

	CmdLineArgParser cmdLine;
	cmdLine.set("-s", std::string("initparfile.xml"), "Input parameter file");
	cmdLine.set("-na", 1, "Number of age groups");
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
	gradient_structure gs(gs_var_buffer);


	// Initially
	int cohort_id = 0; //workerId;
	cout << "Worker ID: " << cohort_id << " argv[argc-1]=" << argv[argc-1] << " cmp_regime=" << cmp_regime << " reset_buffers=" << reset_buffers << std::endl;
	const char* parfile = cmdLine.get<std::string>("-s").c_str();
	SeapodymCohort* scp = seapodym_cohort(parfile, cmp_regime, reset_buffers, cohort_id, gs);
	delete scp;

	// Finalization of MPI
	////////////////////////////////////////////////////////////////////////
	err = MPI_Finalize();

	return 0;
}



