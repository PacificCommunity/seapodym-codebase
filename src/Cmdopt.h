#ifndef __Cmdopt_h__
#define __Cmdopt_h__

#include <string>

// Command-line option helpers, also used by ad_buffers.cpp for the memory options
// (one definition, in Cmdopt.cpp)
char* getCmdOption(char** begin, char** end, const std::string& option);
bool  cmdOptionExists(char** begin, char** end, const std::string& option);

// SEAPODYM executables. Each main_*.cpp passes its own, which sets the help text,
// the version line, the accepted options and the run mode when no option is given.
enum SeapodymBinary {
	BIN_OPTIMIZATION,	// seapodym_clte (main.cpp): with parameter estimation
	BIN_SIMULATION,		// seapodym, seapodym_fluxes (main_simulation.cpp): no parameter estimation
	BIN_HABITATS,		// seapodym_habitats (main_habitats.cpp)
	BIN_DENSITIES		// seapodym_densities (main_densities.cpp)
};

// Parsed command-line options for a SEAPODYM run.
// A plain data bundle: the constructor reads argv into these fields, and the
// worker functions (seapodym_coupled/densities/habitats and
// Sensitivity_analysis) read them back.
struct Cmdopt
{
	int  cmp_regime    = -1;    // run mode: -1 optimization, 0 simulation, 1 likelihood projection,
	                            // 2 Hessian, 3 sensitivity analysis, 4 Taylor test, -11 phases
	int  sftype        = 0;     // sensitivity-analysis subtype: 0 grad, 1 edge, 2 oat, 3 aat
	int  nb_aat        = 1;     // number of all-at-a-time GSA simulations
	bool reset_buffers = false; // whether -mv/-mg/-mc buffer sizes were given
	bool file_exists   = false; // last argument is an existing parfile

	Cmdopt(int argc, char** argv, const SeapodymBinary binary = BIN_OPTIMIZATION);
};

#endif
