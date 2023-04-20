#include <cstdlib>
#include <cstring>
#include <iostream>
using std::cout;
void help(char* argv0);
int OptionToCode(char* Option, int& sub_option);
int seapodym_flux(
    const char* parfile, const int cmp_regime, const int sub_option,
    const bool reset_buffers);
bool read_memory_options(int argc, char** argv, const bool grad_calc);

int main(int argc, char** argv) {
    if (argc == 1) help(argv[0]);
    int cmp_regime = -1;
    int sub_option = 0;
    bool reset_buffers = false;

    int k = 1;
    char* cmdLineOption = argv[k];
    cmp_regime = OptionToCode(cmdLineOption, sub_option);
    if (cmp_regime == -1) {
        k = 0;
        cout << "\n!!!SEAPODYM without parameter estimation. Running in "
                "simulation mode!!!\n";
        cmp_regime = 0;
    }
    if (cmp_regime == -3) help(argv[0]);
    if (argc < k + 2) {
        cout << "Parfile can't be omited... \n";
        help(argv[0]);
    }
    if ((k == 0 && argc > 2) || (cmp_regime >= 0 && argc > 3)) {
        bool grad_calc = false;
        if (cmp_regime == -1 || cmp_regime == 2) grad_calc = true;
        reset_buffers = read_memory_options(argc, argv, grad_calc);
    }

    return seapodym_flux(argv[argc - 1], cmp_regime, sub_option, reset_buffers);
}

int OptionToCode(char* op, int& sub_option) {
    const int N = 12;
    const char* cmdop[N] = {
        "-s",
        "-p",
        "-sa",
        "-h",
        "-v",
        "--simulation",
        "--likelihood-projection",
        "--sensitivity-analysis",
        "--help",
        "--version",
        "-sa=2",
        "-sa=3"};
    int cmpCode[N] = {0, 1, 3, -3, -2, 0, 1, 3, -3, -2, 3, 3};
    for (int i = 0; i < N; i++)
        if (strcmp(op, cmdop[i]) == 0) {
            if (i == 2 || i == 7) sub_option = 3;  // by default sa=3
            if (cmpCode[i] == -2) {
                cout << "SEAPODYM without parameter estimation 4.0 \n";
                cout << "Copyright (C) 2022, SPC, CLS, University of Hawaii.\n";
                exit(0);
            }
            return cmpCode[i];
        }

    return -1;  // no option - simulation
}

void help(char* argv0) {
    cout << "Usage:" << argv0 << " [option] parfile \n";
    cout << "      IMPORTANT!!! If [option] is omitted, then application will "
            "start optimization run! \n";
    cout << "Options: \n";
    cout << "  -h, --help \t\t\t Print this message and exit.\n";
    cout << "  -p, --projection \t\t Compute 2D-projection of the likelihood "
            "on a grid specified in parfile.\n";
    cout << "  -s, --simulation \t\t Run a simulation with fixed parameters.\n";
    cout << "  -sa, --sensitivity-analysis \t Perform global sensitivity "
            "analysis. By default performs one AAT simulations.\n";
    cout << "  -sa=2 \t\t\t ONE-AT-A-TIME sensitivity analysis.\n";
    cout << "  -sa=3 \t\t\t ALL-AT-A-TIME sensitivity analysis.\n";
    cout << "  -v, --version \t\t Print version number and exit.\n";
    exit(0);
}
