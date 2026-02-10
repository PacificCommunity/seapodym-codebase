#include "Cmdopt.h"

int seapodym_habitats(const char* parfile, Cmdopt* cmdopt);

int main(int argc, char** argv) {

	Cmdopt cmdopt(argc, argv);

	return seapodym_habitats(argv[argc-1], &cmdopt);
}


