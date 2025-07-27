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

	int num_age_groups = cmdLine.get<int>("-na");
	std::vector<int> cohort_ids;
	for (int ia = 0; ia < num_age_groups; ++ia) {
		if (ia % num_workers == workerId) {
			cohort_ids.push_back(ia);
		}
	}

	// cohorts handled by this worker
	std::vector< SeapodymCohort* > cohorts;
	const char* parfile = cmdLine.get<std::string>("-s").c_str();

	// Initialize the cohorts for each age group assigned to this worker
	int out_hessian = 0;
	gradient_structure::set_USE_FOR_HESSIAN(out_hessian);

	// //initialize variables of optimization
	// VarParamCoupled var;
	// var.read(parfile);
	// const int nvar = var.nvarcalc();
	// independent_variables x(1, nvar);
	// adstring_array x_names(1, nvar);

	for (auto cohort_id : cohort_ids) {

		SeapodymCohort* scp = seapodym_cohort(parfile, cmp_regime, reset_buffers, cohort_id, gs);

		// //read parfile
		// SeapodymCohort* scp = new SeapodymCohort(parfile, cohort_id);
		// SeapodymCohort& sc = *scp;

		// sc.xinit(x, x_names);
		// cout << "Total number of variables: " << nvar << '\n'<<'\n';

		// //initialization of simulation
		// sc.prerun_model();

		cohorts.push_back(scp);
	}

	// TO DO, run a single step for each cohort habdled by this worker
	for (auto scp : cohorts) {
		//scp->prerun_model();
		scp->OnRunFirstStep();
		// scp->stepForward(false); // false means no output files written
	}

	// Clean up 
	for (auto scp : cohorts) {
		delete scp;
	}

	// Finalization of MPI
	////////////////////////////////////////////////////////////////////////
	err = MPI_Finalize();

	return 0;
}



