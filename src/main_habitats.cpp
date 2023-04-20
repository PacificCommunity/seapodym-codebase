#include <cstdlib>
#include <cstring>
#include <iostream>
using std::cout;
void help(char* argv0);
int OptionToCode(char* Option);
int seapodym_habitats(
    const char* parfile, const int cmp_regime, const bool reset_buffers);
bool read_memory_options(int argc, char** argv, const bool grad_calc);

int main(int argc, char** argv) {
    int cmp_regime = -1;
    bool reset_buffers = false;
    int k = 1;
    char* cmdLineOption = argv[k];
    cmp_regime = OptionToCode(cmdLineOption);
    if (cmp_regime == -1) k = 0;
    if (cmp_regime == -3) help(argv[0]);
    if (argc < k + 2) {
        cout << "Parfile can't be omited... \n";
        help(argv[0]);
    }
    if ((cmp_regime == -1 && argc > 2) || (cmp_regime >= 0 && argc > 3)) {
        bool grad_calc = false;
        if (cmp_regime == -1 || cmp_regime == 2) grad_calc = true;
        reset_buffers = read_memory_options(argc, argv, grad_calc);
    }
    return seapodym_habitats(argv[argc - 1], cmp_regime, reset_buffers);
}

int OptionToCode(char* op) {
    const int N = 8;
    const char* cmdop[N] = {"-s",           "-H",        "-h",     "-v",
                            "--simulation", "--hessian", "--help", "--version"};
    int cmpCode[N] = {0, 2, -3, -2, 0, 2, -3, -2};
    for (int i = 0; i < N; i++)
        if (strcmp(op, cmdop[i]) == 0) {
            if (cmpCode[i] == -2) {
                cout << "SEAPODYM habitats with parameter estimation 4.0 \n";
                cout << "Copyright (C) 2022, SPC, CLS, University of Hawaii.\n";
                exit(0);
            }
            return cmpCode[i];
        }

    return -1;  // by default - optimization
}

void help(char* argv0) {
    cout << "Usage:" << argv0 << " [option] parfile \n";
    cout << "      IMPORTANT!!! If [option] is omitted, then application will "
            "start optimization run! \n";
    cout << "Options: \n";
    cout << "  -h, --help \t\t\t Print this message and exit.\n";
    cout << "  -H, --hessian \t\t Compute Hessian matrix.\n";
    cout << "  -s, --simulation \t\t Run simulation without optimization.\n";
    cout << "  -v, --version \t\t Print version number and exit.\n";
    exit(0);
}
