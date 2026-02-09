#include "Cmdopt.h"

int seapodym_coupled(const char* parfile, Cmdopt* cmdopt);
int seapodym_phases(const char* parfile);

int main(int argc, char** argv) {

	Cmdopt cmdopt(argc, argv);

	if (cmdopt.cmp_regime!=-11)
		return seapodym_coupled(argv[argc-1], &cmdopt);
	else //not supported for now, will not work
		return seapodym_phases(argv[argc-1]);
}


