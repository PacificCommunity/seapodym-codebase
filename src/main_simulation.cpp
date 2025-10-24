#include <iostream>
#include <cstring>
#include <cstdlib>
using std::cout;
void help(char* argv0);
int OptionToCode(char* argv0, char* Option,int &sub_option);
int seapodym_coupled(const char* parfile, const int cmp_regime, const int sub_option, const bool reset_buffers);
bool read_memory_options(int argc, char** argv, const bool grad_calc);

int main(int argc, char** argv) {

	if (argc==1) help(argv[0]);
	int cmp_regime = -1;
	int sub_option = 1;
	bool reset_buffers = false;
	
	int k=1;
	char *cmdLineOption = argv[k];
	cmp_regime = OptionToCode(argv[0],cmdLineOption, sub_option);
	if (cmp_regime==-1){
		k=0;
		cout << "\n!!!SEAPODYM without parameter estimation. Running in simulation mode!!!\n"; 
		cmp_regime = 0;
	}
	if (cmp_regime==-3) help(argv[0]);
	if (argc < k+2) {
		cout << "Parfile can't be omited... \n"; 
		help(argv[0]);
	}
	if ((k==0 && argc>2) || (cmp_regime>=0 && argc>3)){
		bool grad_calc = false;
		if (cmp_regime == -1 || cmp_regime == 2)
			grad_calc = true;
		reset_buffers = read_memory_options(argc, argv, grad_calc);	
	}

	return seapodym_coupled(argv[argc-1],cmp_regime,sub_option,reset_buffers);
}


int OptionToCode(char* argv0, char* op, int &sub_option) {

	const int N = 13;
	const char *cmdop[N] = {"-s","-p","-sa","-h","-v","--simulation","--likelihood-projection","--edge-sensitivity","--help","--version","-sa=1","-sa=2","-sa=3"};
	int cmpCode[N] = {0,1,3,-3,-2,0,1,3,-3,-2,3,3};
	for (int i=0; i<N; i++)
		if (strcmp(op,cmdop[i])==0){
			if (i>=N-2) sub_option = 2;
			if (i>=N-1) sub_option = 3;
			if (cmpCode[i]==-2) {
				cout << "SEAPODYM without parameter estimation 4.0 \n";
				cout << "Copyright (C) 2022, SPC, CLS, University of Hawaii.\n";
				exit(0);
			}
			return cmpCode[i];
		}
		
	cout << "\nWrong entry: no such option.\n\n";
	help(argv0);

	return -1;//default here is simulation
}

void help(char* argv0) {

	cout << "Usage:" << argv0 << " [option] parfile \n";
	cout << "      IMPORTANT!!! This application does not include parameter estimation! \n"; 
	cout << "Options: \n";
	cout << "  -h, --help \t\t\t Print this message and exit.\n";
	cout << "  -p, --projection \t\t Compute 2D-projection of the likelihood on a grid specified in parfile.\n";
	cout << "  -s, --simulation \t\t Run a simulation with fixed parameters.\n";
	cout << "  -sa[=1], --edge-sensitivity\t By default[FLAG=1] computes likelihood changes within parameter boundaries.\n";
	cout << "  -sa=2 \t\t\t Runs ONE-AT-A-TIME sensitivity simulations.\n";
	cout << "  -sa=3 \t\t\t Performs one ALL-AT-A-TIME simulation for GSA.\n";	
	cout << "  -v, --version \t\t Print version number and exit.\n";
	exit(0);
}

