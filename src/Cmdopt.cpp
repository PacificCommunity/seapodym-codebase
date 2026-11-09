#include <iostream>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include "Cmdopt.h"
#include <sys/stat.h>

using namespace std;

bool read_memory_options(int argc, char** argv, const bool grad_calc);

// ---- shared command-line helpers (declared in Cmdopt.h) --------------------

char* getCmdOption(char** begin, char** end, const std::string& option)
{
	char** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
		return *itr;
	return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
	return std::find(begin, end, option) != end;
}

// ---- file-local parsing helpers (no external callers) ----------------------
namespace {

// What each executable can do
struct BinaryInfo {
	const char* version;	// version line
	bool estimation;	// parameter estimation (adjoint code): optimization, Hessian, Taylor test, -sa grad
	bool phases;		// -ph, --phases (seapodym_clte only, not functional for now, so not in the help)
};

BinaryInfo binary_info(const SeapodymBinary bin)
{
	BinaryInfo b;
	b.estimation = true; b.phases = false;
	switch (bin) {
		case BIN_SIMULATION: b.version = "SEAPODYM without parameter estimation 5.X"; b.estimation = false; break;
		case BIN_HABITATS:   b.version = "SEAPODYM habitats with parameter estimation 5.X"; break;
		case BIN_DENSITIES:  b.version = "SEAPODYM densities with parameter estimation 5.X"; break;
		default:             b.version = "SEAPODYM with parameter estimation 5.X"; b.phases = true;
	}
	return b;
}

void help(char* argv0, const SeapodymBinary bin, const int status = 0)
{
	const BinaryInfo b = binary_info(bin);
	cout << "Usage:" << argv0 << " [option] parfile \n";
	if (b.estimation)
		cout << "      IMPORTANT!!! If [option] is omitted, then application will start optimization run! \n";
	else
		cout << "      IMPORTANT!!! This application does not include parameter estimation! If [option] is omitted, it runs a simulation. \n";
	cout << "Options: \n";
	cout << "  -h, --help \t\t\t Print this message and exit.\n";
	cout << "  -v, --version \t\t Print version number and exit.\n";
	if (b.estimation)
		cout << "  -H, --hessian \t\t Compute Hessian matrix.\n";
	cout << "  -p, --projection \t\t Compute 2D-projection of the likelihood on a grid specified in parfile.\n";
	if (b.estimation) {
		cout << "  -t, --taylor-test   \t\t Perform Taylor derivative test with central differencing.\n";
		cout << "  -s, --simulation \t\t Run simulation without optimization.\n";
		cout << "  -sa <string> \t\t\t grad - Computes local sensititivy analysis (Default)\n\t\t\t\t edge - Computes likelihood changes within parameter boundaries\n\t\t\t\t oat - Runs ONE-AT-A-TIME sensitivity simulations\n\t\t\t\t aat - Performs ALL-AT-A-TIME simulations for GSA.\n";
		cout << "  --local-sensitivity\t\t Computes local sensitivities. Equivalent to -sa grad.\n";
	} else {
		cout << "  -s, --simulation \t\t Run a simulation with fixed parameters (Default).\n";
		cout << "  -sa <string> \t\t\t edge - Computes likelihood changes within parameter boundaries (Default)\n\t\t\t\t oat - Runs ONE-AT-A-TIME sensitivity simulations\n\t\t\t\t aat - Performs ALL-AT-A-TIME simulations for GSA.\n";
	}
	cout << "  -na, --number-aat <int> \t Number of ALL-AT-A-TIME simulations for GSA (Default: 1). \n";
	if (b.estimation) {
		cout << "\nIf ADDITIONAL MEMORY is needed, add after the main option (after binary name for optimization run): \n";
		cout << "  -mv <integer_number> \t\t Size of gs_var_buffer - the buffer for model variables. \n";
		cout << "  -mg <integer_number> \t\t Size of gradstack_buffer - the buffer for automatic differentiation.\n";
		cout << "  -mc <integer_number> \t\t Size of cmpdif_buffer - the buffer for differentiation with adjoint code.\n";
	} else {
		cout << "\nIf ADDITIONAL MEMORY is needed, add after the main option: \n";
		cout << "  -mv <integer_number> \t\t Size of gs_var_buffer - the buffer for model variables. \n";
	}
	exit(status);
}

void CheckInfo(char* cmdop, char* argv0, const SeapodymBinary bin)
{
	if (strcmp(cmdop, "-h") == 0 || strcmp(cmdop, "--help") == 0)
		help(argv0, bin);

	if (strcmp(cmdop, "-v") == 0 || strcmp(cmdop, "--version") == 0) {
		cout << binary_info(bin).version << " \n";
		cout << "Copyright (C) 2025, SPC, CLS, University of Hawaii.\n";
		exit(0);
	}
}

// Decode the leading option into opt.cmp_regime; options not available in this
// executable, and unknown options, print the help and stop.
// Memory options and -na in the leading position keep the default run mode.
void option_to_code(Cmdopt& opt, char* op, char* argv0, const SeapodymBinary bin)
{
	const BinaryInfo b = binary_info(bin);
	const int KEEP = -100;
	const int N = 17;
	const char* cmdop[N] = {"-ph","-s","-p","-H","-sa","-na","-t","-mv","-mc","-mg",
	                        "--phases","--simulation","--likelihood-projection","--hessian",
	                        "--local-sensitivity","--taylor-test","--number-aat"};
	int cmpCode[N] = {-11,0,1,2,3,KEEP,4,KEEP,KEEP,KEEP,-11,0,1,2,3,4,KEEP};
	for (int i = 0; i < N; i++)
		if (strcmp(op, cmdop[i]) == 0) {
			const int code = cmpCode[i];
			const bool available = !((code == -11 && !b.phases) ||
			                         ((code == 2 || code == 4) && !b.estimation) ||
			                         (strcmp(op, "--local-sensitivity") == 0 && !b.estimation));
			if (!available) {
				cout << "\nWRONG ENTRY: option " << op << " is not available in " << argv0 << "\n\n";
				help(argv0, bin, 1);
			}
			if (code != KEEP) opt.cmp_regime = code;
			return;
		}

	cout << "\nWRONG ENTRY: no such option: " << op << "\n\n";
	help(argv0, bin, 1);
}

// Read the -sa subtype and, for aat, the number of runs, into opt.
void read_sa_options(Cmdopt& opt, int argc, char** argv, char* argv0, const SeapodymBinary bin)
{
	if (opt.cmp_regime == 3) {
		const bool estimation = binary_info(bin).estimation;
		opt.sftype = estimation ? 0 : 1;	// default: grad, or edge without the adjoint code
		char* op = getCmdOption(argv, argv + argc - 1, "-sa");
		const int N = 4;
		const char* sfop[N] = {"grad","edge","oat","aat"};
		int sftypes[N] = {0,1,2,3};
		if (op)	// '-sa' directly followed by the parfile gives no subtype: keep the default
			for (int i = 0; i < N; i++)
				if (strcmp(op, sfop[i]) == 0)
					opt.sftype = sftypes[i];
		if (opt.sftype == 0 && !estimation) {
			cout << "\nWRONG ENTRY: -sa grad needs the gradient and is not available in " << argv0 << "\n\n";
			help(argv0, bin, 1);
		}

		if (opt.sftype == 3) {
			if (cmdOptionExists(argv, argv + argc, "-na"))
				opt.nb_aat = atol(getCmdOption(argv, argv + argc - 1, "-na"));
			else if (cmdOptionExists(argv, argv + argc, "--number-aat"))
				opt.nb_aat = atol(getCmdOption(argv, argv + argc - 1, "--number-aat"));
		}
	}
}

} // anonymous namespace

// ---- the one public entry point --------------------------------------------

Cmdopt::Cmdopt(int argc, char** argv, const SeapodymBinary binary)
{
	const bool estimation = binary_info(binary).estimation;
	cmp_regime = estimation ? -1 : 0;	// run mode when no option is given

	if (argc == 1) help(argv[0], binary);

	char* cmdop = argv[1];
	if (argc == 2) CheckInfo(cmdop, argv[0], binary);

	struct stat buffer;
	if (stat(argv[argc - 1], &buffer) == 0)
		file_exists = true;

	if (!file_exists) {
		cout << "The last argument must be existing parfile...\n";
		help(argv[0], binary, 1);
	}

	if (argc > 2)
		option_to_code(*this, cmdop, argv[0], binary);

	if (!estimation && cmp_regime == 0)
		cout << "\n!!!SEAPODYM without parameter estimation. Running in simulation mode!!!\n";

	if ((cmp_regime == -1 && argc > 2) || (cmp_regime >= 0 && argc > 3)) {
		bool grad_calc = (cmp_regime == -1 || cmp_regime == 2);
		reset_buffers = read_memory_options(argc, argv, grad_calc);
	}

	read_sa_options(*this, argc, argv, argv[0], binary);
}
