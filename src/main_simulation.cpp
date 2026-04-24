#include "Cmdopt.h"

int seapodym_coupled(const char* parfile, Cmdopt* cmdopt);

int main(int argc, char** argv) {

	Cmdopt cmdopt(argc, argv);

	return seapodym_coupled(argv[argc-1], &cmdopt);
}


