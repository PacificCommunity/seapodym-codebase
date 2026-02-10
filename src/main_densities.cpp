#include "Cmdopt.h"

int seapodym_densities(const char* parfile, Cmdopt* cmdopt);

int main(int argc, char** argv) {

	Cmdopt cmdopt(argc, argv);

	return seapodym_densities(argv[argc-1], &cmdopt);
}

