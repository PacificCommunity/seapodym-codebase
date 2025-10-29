#include <iostream>
#include <cstring>
using std::cout;
void help(char* argv0);
void CheckInfo(char* Option, char* argv0);
int OptionToCode(char* Option,int& sub_option);


void CheckInfo(char* cmdop, char* argv0) {

	if (strcmp(cmdop,"-h")==0 || strcmp(cmdop,"--help")==0) 
		help(argv0);
	
	if (strcmp(cmdop,"-v")==0 || strcmp(cmdop,"--version")==0){
		cout << "SEAPODYM habitats with parameter estimation 5.X \n";
		cout << "Copyright (C) 2025, SPC, CLS, University of Hawaii.\n";
		exit(0);
	}
}

int OptionToCode(char* op, int& sub_option) {

	const int N = 19;
	const char *cmdop[N] = {"-ph","-s","-p","-H","-sa","-t","-mv","-mc","-mg","--phases","--simulation", "--likelihood-projection","--hessian","--local-sensitivity","--taylor-test","-sa=0","-sa=1","-sa=2","-sa=3"};
	int cmpCode[N] = {-11,0,1,2,3,4,-1,-1,-1,-11,0,1,2,3,4,3,3,3,3};
	for (int i=0; i<N; i++)
		if (strcmp(op,cmdop[i])==0){
			if (i>=N-3) sub_option = 1;
			if (i==N-2) sub_option = 2;
			if (i==N-1) sub_option = 3;			
			return cmpCode[i];
		}
	
	cout << "\nWRONG ENTRY: no such option - will run a simulation...\n";
	return 0;
}


void help(char* argv0) {

	cout << "Usage:" << argv0 << " [option] parfile \n";
	cout << "      IMPORTANT!!! If [option] is omitted, then application will start optimization run! \n"; 
	cout << "Options: \n";
	cout << "  -h, --help \t\t\t Print this message and exit.\n";
	cout << "  -v, --version \t\t Print version number and exit.\n";
	cout << "  -H, --hessian \t\t Compute Hessian matrix.\n";
	cout << "  -p, --projection \t\t Compute 2D-projection of the likelihood on a grid specified in parfile.\n";
	cout << "  -t, --taylor-test   \t\t Perform Taylor derivative test with central differencing.\n";
	cout << "  -s, --simulation \t\t Run simulation without optimization.\n";
	cout << "  -sa[=0], --local-sensitivity\t By default[FLAG=0] computes local sensitivities.\n";
	cout << "  -sa=1 \t\t\t Computes likelihood changes within parameter boundaries.\n";
	cout << "  -sa=2 \t\t\t Runs ONE-AT-A-TIME sensitivity simulations.\n";
	cout << "  -sa=3 \t\t\t Performs one ALL-AT-A-TIME simulation for GSA.\n";	
	cout << "\nIf ADDITIONAL MEMORY is needed, add after the main option (after binary name for optimization run): \n"; 
	cout << "  -mv <integer_number> \t\t Size of gs_var_buffer - the buffer for model variables. \n";	
	cout << "  -mg <integer_number> \t\t Size of gradstack_buffer - the buffer for automatic differentiation.\n";	
	cout << "  -mc <integer_number> \t\t Size of cmpdif_buffer - the buffer for differentiation with adjoint code.\n";	

	exit(0);
}

void help_sim(char* argv0) {

	cout << "Usage:" << argv0 << " [option] parfile \n";
	cout << "      IMPORTANT!!! This application does not include parameter estimation! \n"; 
	cout << "Options: \n";
	cout << "  -h, --help \t\t\t Print this message and exit.\n";
	cout << "  -v, --version \t\t Print version number and exit.\n";
	cout << "  -p, --projection \t\t Compute 2D-projection of the likelihood on a grid specified in parfile.\n";
	cout << "  -s, --simulation \t\t Run a simulation with fixed parameters.\n";
	cout << "  -sa[=1], --edge-sensitivity\t By default[FLAG=1] computes likelihood changes within parameter boundaries.\n";
	cout << "  -sa=2 \t\t\t Runs ONE-AT-A-TIME sensitivity simulations.\n";
	cout << "  -sa=3 \t\t\t Performs one ALL-AT-A-TIME simulation for GSA.\n";	
	cout << "\nIf ADDITIONAL MEMORY is needed, add after the main option: \n"; 
	cout << "  -mv <integer_number> \t\t Size of gs_var_buffer - the buffer for model variables. \n";	
	exit(0);
}

