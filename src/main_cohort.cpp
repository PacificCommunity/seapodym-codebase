#include <iostream>
#include <cstring>
#include <cstdlib>
#include <mpi.h>
#include "SeapodymCohort.h"

using std::cout;
void help(char* argv0);
int OptionToCode(char* Option);
SeapodymCohort* seapodym_cohort(const char* parfile, const int cmp_regime, const bool reset_buffers, int cohort_id, gradient_structure& gs);
bool read_memory_options(int argc, char** argv, const bool grad_calc);
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
	SeapodymCohort* scp = seapodym_cohort(argv[argc-1],cmp_regime,reset_buffers, cohort_id, gs);
	delete scp;

	// Finalization of MPI
	////////////////////////////////////////////////////////////////////////
	err = MPI_Finalize();

	return 0;
}


int OptionToCode(char* op) {

	const int N = 12;
	const char *cmdop[N] = {"-s","-p","-H","-t","-h","-v","--simulation","--likelihood-projection","--hessian","--taylor-test","--help","--version"};
	int cmpCode[N] = {0,1,2,4,-3,-2,0,1,2,4,-3,-2};
	for (int i=0; i<N; i++)
		if (strcmp(op,cmdop[i])==0){
			if (cmpCode[i]==-2) {
				cout << "SEAPODYM 4.0 without fishing for parameter estimation using population density \n";
				cout << "Copyright (C) 2022, SPC, CLS, University of Hawaii.\n";
				exit(0);
			}
			return cmpCode[i];
		}
	return -1;//by default - optimization
}

void help(char* argv0) {

	cout << "Usage:" << argv0 << " [option] parfile \n";
	cout << "      IMPORTANT!!! If [option] is omitted, then application will start optimization run! \n"; 
	cout << "Options: \n";
	cout << "  -h, --help \t\t\t Print this message and exit.\n";
	cout << "  -H, --hessian \t\t Compute Hessian matrix.\n";
	cout << "  -t, --taylor-test   \t\t Perform Taylor derivative test with central differencing.\n";
	cout << "  -s, --simulation \t\t Run simulation without optimization.\n";
	cout << "  -v, --version \t\t Print version number and exit.\n";
	exit(0);
}



