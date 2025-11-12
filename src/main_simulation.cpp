#include <iostream>
#include <cstring>
#include <sys/stat.h>
using std::cout;

void help_sim(char* argv0);
void CheckInfo(char* Option, char* argv0);
int OptionToCode(char* Option,int &sub_option);
int seapodym_coupled(const char* parfile, const int cmp_regime, const int sub_option, const bool reset_buffers);
bool read_memory_options(int argc, char** argv, const bool grad_calc);

int main(int argc, char** argv) {

	if (argc==1) help_sim(argv[0]);

	char *cmdop = argv[1];
	if (argc==2) CheckInfo(cmdop,argv[0]);

	int cmp_regime = 0;
	int sub_option = 0;
	bool reset_buffers = false;
	bool file_exists = false;

	struct stat buffer;  
	if (stat(argv[argc-1], &buffer)==0)
		file_exists = true;

	if (!file_exists){
		cout << "The last argument must be existing parfile...\n";
		help_sim(argv[0]);
	}
	if (argc>2 && file_exists) 
		cmp_regime = OptionToCode(cmdop,sub_option);
	if (cmp_regime==-1 || cmp_regime==-11 || cmp_regime==4){
		cout << "\n!!!SEAPODYM without parameter estimation. Running in simulation mode!!!\n"; 
		cmp_regime = 0;
	}
	if (cmp_regime==3 && sub_option==0)
		sub_option = 1;

	if (((cmp_regime==0 && argc>2) || (cmp_regime>=0 && argc>3)) && file_exists){
		bool grad_calc = false;
		if (cmp_regime == -1 || cmp_regime == 2)
			grad_calc = true;
		reset_buffers = read_memory_options(argc, argv, grad_calc);	
	}

	return seapodym_coupled(argv[argc-1],cmp_regime,sub_option,reset_buffers);
}


