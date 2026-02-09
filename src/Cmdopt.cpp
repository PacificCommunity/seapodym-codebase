#include <iostream>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include "Cmdopt.h"
#include <sys/stat.h>

bool read_memory_options(int argc, char** argv, const bool grad_calc);

Cmdopt::Cmdopt(int argc, char** argv){
	
	if (argc==1) help(argv[0]);

	char *cmdop = argv[1];
	if (argc==2) CheckInfo(cmdop,argv[0]);

	struct stat buffer;  
	if (stat(argv[argc-1], &buffer)==0)
		file_exists = true;

	if (!file_exists){
		cout << "The last argument must be existing parfile...\n";
		help(argv[0]);
	}
	if (argc>2 && file_exists) 
		OptionToCode(cmdop);

	if (((cmp_regime==-1 && argc>2) || (cmp_regime>=0 && argc>3)) && file_exists){
		bool grad_calc = false;
		if (cmp_regime == -1 || cmp_regime == 2)
			grad_calc = true;
		reset_buffers = read_memory_options(argc, argv, grad_calc);	
	}

	if (file_exists){
		read_sa_options(argc, argv);
	}

}

void Cmdopt::CheckInfo(char* cmdop, char* argv0) {

	if (strcmp(cmdop,"-h")==0 || strcmp(cmdop,"--help")==0) 
		help(argv0);
	
	if (strcmp(cmdop,"-v")==0 || strcmp(cmdop,"--version")==0){
		cout << "SEAPODYM habitats with parameter estimation 5.X \n";
		cout << "Copyright (C) 2025, SPC, CLS, University of Hawaii.\n";
		exit(0);
	}
}

void Cmdopt::OptionToCode(char* op) {

	const int N = 17;
	const char *cmdop[N] = {"-ph","-s","-p","-H","-sa","-na","-t","-mv","-mc","-mg","--phases","--simulation", "--likelihood-projection","--hessian","--local-sensitivity","--taylor-test","--number-aat"};
	int cmpCode[N] = {-11,0,1,2,3,-1,4,-1,-1,-1,-11,0,1,2,3,4,-1};
	for (int i=0; i<N; i++)
		if (strcmp(op,cmdop[i])==0){	
			cmp_regime = cmpCode[i];
		}
	
	cout << "\nWRONG ENTRY: no such option - will run a simulation...\n";
}


void Cmdopt::help(char* argv0) {

	cout << "Usage:" << argv0 << " [option] parfile \n";
	cout << "      IMPORTANT!!! If [option] is omitted, then application will start optimization run! \n"; 
	cout << "Options: \n";
	cout << "  -h, --help \t\t\t Print this message and exit.\n";
	cout << "  -v, --version \t\t Print version number and exit.\n";
	cout << "  -H, --hessian \t\t Compute Hessian matrix.\n";
	cout << "  -p, --projection \t\t Compute 2D-projection of the likelihood on a grid specified in parfile.\n";
	cout << "  -t, --taylor-test   \t\t Perform Taylor derivative test with central differencing.\n";
	cout << "  -s, --simulation \t\t Run simulation without optimization.\n";
	cout << "  -sa <int> \t\t\t 0 - Computes local sensititivy analysis (Default)\n\t\t\t\t 1 - Computes likelihood changes within parameter boundaries\n\t\t\t\t 2 - Runs ONE-AT-A-TIME sensitivity simulations\n\t\t\t\t 3 - Performs ALL-AT-A-TIME simulations for GSA.\n";
	cout << "  --local-sensitivity\t\t Computes local sensitivities. Equivalent to -sa 0.\n";
	cout << "  -na, --number-aat <int> \t Number of ALL-AT-A-TIME simulations for GSA (Default: 1). \n";	
	cout << "\nIf ADDITIONAL MEMORY is needed, add after the main option (after binary name for optimization run): \n"; 
	cout << "  -mv <integer_number> \t\t Size of gs_var_buffer - the buffer for model variables. \n";	
	cout << "  -mg <integer_number> \t\t Size of gradstack_buffer - the buffer for automatic differentiation.\n";	
	cout << "  -mc <integer_number> \t\t Size of cmpdif_buffer - the buffer for differentiation with adjoint code.\n";	

	exit(0);
}

void Cmdopt::help_sim(char* argv0) {

	cout << "Usage:" << argv0 << " [option] parfile \n";
	cout << "      IMPORTANT!!! This application does not include parameter estimation! \n"; 
	cout << "Options: \n";
	cout << "  -h, --help \t\t\t Print this message and exit.\n";
	cout << "  -v, --version \t\t Print version number and exit.\n";
	cout << "  -p, --projection \t\t Compute 2D-projection of the likelihood on a grid specified in parfile.\n";
	cout << "  -s, --simulation \t\t Run a simulation with fixed parameters.\n";
	cout << "  -sa <int> \t\t\t 0 - Computes local sensititivy analysis (Default)\n\t\t\t\t 1 - Computes likelihood changes within parameter boundaries\n\t\t\t\t 2 - Runs ONE-AT-A-TIME sensitivity simulations\n\t\t\t\t 3 - Performs ALL-AT-A-TIME simulations for GSA.\n";
	cout << "  --local-sensitivity\t\t Computes local sensitivities. Equivalent to -sa 0.\n";
	cout << "  -na, --number-aat <int> \t Number of ALL-AT-A-TIME simulations for GSA (Default: 1). \n";	
	cout << "\nIf ADDITIONAL MEMORY is needed, add after the main option: \n"; 
	cout << "  -mv <integer_number> \t\t Size of gs_var_buffer - the buffer for model variables. \n";	
	exit(0);
}

char* Cmdopt::getCmdOption(char ** begin, char ** end, const std::string& option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end){
        return *itr;
    }
    return 0;
}

bool Cmdopt::cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

void Cmdopt::read_sa_options(int argc, char** argv)
{
	if (cmp_regime==3){

		sftype = atol(getCmdOption(argv, argv + argc -1, "-sa"));

		// Read number of AAT
		if (sftype==3){
			if(cmdOptionExists(argv, argv+argc, "-na")){
				nb_aat = atol(getCmdOption(argv, argv + argc -1, "-na"));
			}else if(cmdOptionExists(argv, argv+argc, "--number-aat")){
				nb_aat = atol(getCmdOption(argv, argv + argc -1, "--number-aat"));
			}
		}
	}
}
