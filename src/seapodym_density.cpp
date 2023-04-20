#include <fvar.hpp>

#include "SeapodymCoupled.h"

string get_path(const char* full_path);
void Hessian_comp(const char* parfile);
void buffers_init(
    long int& mv, long int& mc, long int& mg, const bool grad_calc);
void buffers_set(long int& mv, long int& mc, long int& mg);

/*!
\brief The first function to be executed.
\details
This is the main routine that calls upper-level functions such as
1. Parfile reading and initialization of model parameters and optimization
variables;
2. Running the application in different regimes:
   a) (default) running the model in forward (simulation) and backward (gradient
computation) mode; b) simulation only with offline coupling with forage
sub-model; c) running coupled simulation for tuna-forage model; d) computing
Hessian; e) sensitivity analysis; f) computing 2d projection of likelihood
function the pair of parameters (should be specified in parfile).
*/

int seapodym_densities(
    const char* parfile, int cmp_regime, const bool reset_buffers) {
    time_t time_sec;
    time(&time_sec);
    const time_t time0 = time_sec;

    //-----Memory stack sizes for dvariables and derivatives storage------
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    long int gradstack_buffer, cmpdif_buffer, gs_var_buffer;
    bool grad_calc = false;
    if (cmp_regime == -1 || cmp_regime == 2) grad_calc = true;
    buffers_init(gs_var_buffer, gradstack_buffer, cmpdif_buffer, grad_calc);
    if (reset_buffers)
        buffers_set(gs_var_buffer, gradstack_buffer, cmpdif_buffer);

    gradient_structure::set_GRADSTACK_BUFFER_SIZE(gradstack_buffer);
    gradient_structure::set_CMPDIF_BUFFER_SIZE(cmpdif_buffer);
    gradient_structure gs(gs_var_buffer);
    //--------------------------------------------------------------------

    int out_hessian = 0;
    gradient_structure::set_USE_FOR_HESSIAN(out_hessian);

    cout << "\nstarting time: " << ctime(&time_sec) << endl;

    // if mode 2, redirecting to Hessian routine and exit.
    if (cmp_regime == 2) {
        Hessian_comp(parfile);
        return 0;
    }

    // read parfile
    SeapodymCoupled sc(parfile);

    // iniitalize variables of optimization
    const int nvar = sc.nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1, nvar);

    sc.xinit(x, x_names);
    cout << "Total number of variables: " << nvar << '\n' << '\n';

    // writing temporal output to parfile folder

    // function minimizer class
    fmm fmc(nvar);

    // flags for function minimizer
    fmc.iprint = 1;
    fmc.crit = sc.get_crit();  // 0.1;
    fmc.imax = 30;
    fmc.scroll_flag = 1;
    fmc.ifn = 0;
    fmc.maxfn = sc.get_maxfn();  // 2000;
    if (fmc.maxfn <= 0) fmc.ireturn = -1;

    // if this flag is 0 then the gradient will not be computed
    int compute_gradient = 1;
    if (cmp_regime == 0) {
        sc.param->set_gradcalc(false);
        compute_gradient = 0;
    }

    double likelihood = 0;
    double elapsed_time = 0;

    // gradient vector allocation and initialization
    dvector g(1, nvar);
    g.initialize();

    int idx = 0;
    int itr = 0;
    ios::sync_with_stdio();

    // initialization of simulation parameters
    sc.OnRunFirstStep();
    sc.ReadDensity();
    // the function is invoked in the coupled simulation only
    string tempparfile = "tempparfile.xml";
    string newparfile = "newparfile.xml";

    // clock_t time1 = clock();
    time(&time_sec);
    time_t time1 = time_sec;
    //'run_coupled' runs the model in forward (simulation) mode
    //'gradcalc' runs backward (adjoint) mode
    //'save_statistics' stores current information on function minimization
    if (compute_gradient) {
        string dirout = get_path(parfile);
        tempparfile = dirout + "/tempparfile.xml";
        newparfile = dirout + "/newparfile.xml";
        cout << "\nentering minimization loop" << endl;
        while (fmc.ireturn >= 0) {
            // call function minimizer
            fmc.fmin(likelihood, x, g);

            // update the statistics.out file if the solution has been improved
            int itn = fmc.itn;
            //			if ((itn==itr+1 && (idx>=1))){
            time(&time_sec);
            time_t time2 = time_sec;
            elapsed_time += (double)(time2 - time1) / 60.0;
            time1 = time2;
            sc.save_statistics(
                dirout, x_names, likelihood, g, elapsed_time, idx - 1, itr,
                nvar);
            sc.write(tempparfile.c_str());
            itr = itn;
            //			}
            // reset control parameters and run the model and its adjoint
            if (fmc.ireturn > 0) {
                likelihood = sc.run_density((dvar_vector)x);
                gradcalc(nvar, g);
                cout << "function evaluation " << idx++ << endl;
            }
        }
    }

    // after minimization is finished one simulation will
    // be run with estimated parameters; outputs will be saved
    gradient_structure::set_NO_DERIVATIVES();
    sc.run_density((dvar_vector)x, true);
    sc.write(newparfile.c_str());

    remove(tempparfile.c_str());

    // writes new parameters on the screen
    sc.param->outp_param(x_names, nvar);

    time(&time_sec);
    cout << "\nfinished time: " << ctime(&time_sec) << endl;

    time_t time2 = time_sec;
    double total_elapsed_time = (double)(time2 - time0) / 60.0;
    cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;
    // exit(1);

    return 0;
}

double run_model(
    SeapodymCoupled& sc, dvar_vector x, dvector& g, const int nvar) {
    double like = 0.0;
    like = sc.run_density((dvar_vector)x);
    gradcalc(nvar, g);
    return like;
}

void verify_identifier_string2(char* str1)  // ASSUME str1 is not null
{
    // Back up the stream and read the number of bytes written in the
    // ``write function'' corresponding to this ``read function''
    long int num_bytes = strlen(str1);
    char* str = new char[num_bytes + 1];
    str[num_bytes] = '\0';
    gradient_structure::get_fp()->fread(str, num_bytes);
    if (strcmp(str1, str)) {
        cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: \"" << str
             << "\" != \"" << str1 << "\"\n";
        ad_exit(1);
    }

    if (str) delete[] str;
}

string get_path(const char* parfile) {
    size_t pos;
    string full_path = string(parfile);
    pos = full_path.find_last_of("/");
    if (pos == string::npos) return ".";
    string dir = full_path.substr(0, pos);
    return dir;
}

int save_identifier_string2(char* str) {
    int length = strlen(str);
    gradient_structure::get_fp()->fwrite(str, length);
    return 0;
}

void save_long_int_value(unsigned long int x) {
    int num_bytes = sizeof(unsigned long int);
    void* y = (void*)&x;
    gradient_structure::get_fp()->fwrite(y, num_bytes);
}

unsigned long int restore_long_int_value(void) {
    void* tmpout;
    int num_bytes = sizeof(unsigned long int);
    gradient_structure::get_fp()->fread(&tmpout, num_bytes);
    return (unsigned long int)tmpout;
}
