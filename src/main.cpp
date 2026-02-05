#include <iostream>
#include <cstring>
#include <sys/stat.h>
#include "ad_options.h"
using std::cout;

void help(char* argv0);
void CheckInfo(char* Option, char* argv0);
int OptionToCode(char* Option,int &sub_option);
int seapodym_coupled(const char* parfile, const int cmp_regime, const int sub_option, const int sub_option2, const bool reset_buffers);
int seapodym_phases(const char* parfile);

int main(int argc, char** argv) {

	if (argc==1) help(argv[0]);

	char *cmdop = argv[1];
	if (argc==2) CheckInfo(cmdop,argv[0]);
	
	int cmp_regime = -1;
	int sub_option = 0;
	int sub_option2 = 0;
	bool reset_buffers = false;
	bool file_exists = false;

	struct stat buffer;  
	if (stat(argv[argc-1], &buffer)==0)
		file_exists = true;

	if (!file_exists){
		cout << "The last argument must be existing parfile...\n";
		help(argv[0]);
	}
	if (argc>2 && file_exists) 
		cmp_regime = OptionToCode(cmdop,sub_option);

	if (((cmp_regime==-1 && argc>2) || (cmp_regime>=0 && argc>3)) && file_exists){
		bool grad_calc = false;
		if (cmp_regime == -1 || cmp_regime == 2)
			grad_calc = true;
		reset_buffers = read_memory_options(argc, argv, grad_calc);	
	}

	if (cmp_regime==3 && file_exists){
		sub_option2 = read_sub_option2(argc, argv, cmp_regime);
	}

	if (cmp_regime!=-11)
		return seapodym_coupled(argv[argc-1],cmp_regime,sub_option,sub_option2,reset_buffers);
	else //not supported for now, will not work
		return seapodym_phases(argv[argc-1]);
}


