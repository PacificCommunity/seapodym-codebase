#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>

char* getCmdOption(char** begin, char** end, const std::string& option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);
void buffers_get(const long int mv, const long int mg, const long int mc);
void buffers_init(
    long int& mv, long int& mc, long int& mg, const bool grad_calc);

// Default settings are for:
// PO-2deg bigeye config, 22-years simulation (1.6Gb)
// or 6-years optimization with CL likelihood (15Gb).
long int gs_var_size_set = 100000000L;
long int gradstack_size_set = 5000000L;
long int cmpdiff_size_set = 10000000L;

bool read_memory_options(int argc, char** argv, const bool grad_calc) {
    bool reset_buffer = false;
    long int mv, mc, mg;
    buffers_init(mv, mg, mc, grad_calc);

    if (cmdOptionExists(argv, argv + argc, "-mv")) {
        mv = atol(getCmdOption(argv, argv + argc - 1, "-mv"));
        reset_buffer = true;
    }

    if (cmdOptionExists(argv, argv + argc, "-mg")) {
        mg = atol(getCmdOption(argv, argv + argc - 1, "-mg"));
        reset_buffer = true;
    }

    if (cmdOptionExists(argv, argv + argc, "-mc")) {
        mc = atol(getCmdOption(argv, argv + argc - 1, "-mc"));
        reset_buffer = true;
    }

    buffers_get(mv, mg, mc);

    return reset_buffer;
}

void buffers_init(
    long int& mv, long int& mg, long int& mc, const bool grad_calc) {
    //-----SIMULATION MODE----------
    mv = gs_var_size_set;
    mg = gradstack_size_set;
    mc = cmpdiff_size_set;
    //------------------------------

    if (grad_calc) {
        //-------GRADCALC MODE----------
        mv *= 2;
        mg *= 25;    // 50;
        mc *= 1000;  // 2000;
        //------------------------------
    }
}

void buffers_get(const long int mv, const long int mg, const long int mc) {
    gs_var_size_set = mv;
    gradstack_size_set = mg;
    cmpdiff_size_set = mc;
}

void buffers_set(long int& mv, long int& mg, long int& mc) {
    mv = gs_var_size_set;
    mg = gradstack_size_set;
    mc = cmpdiff_size_set;
}

char* getCmdOption(char** begin, char** end, const std::string& option) {
    char** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}
