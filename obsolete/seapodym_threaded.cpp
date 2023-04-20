#include <fvar.hpp>

#include "SeapodymCoupled.h"
#include "adthread.h"
#include "sys/stat.h"

string get_path(const char* full_path);
void Hyperspace_projection(SeapodymCoupled& sc, dvar_vector x);
void Hessian_comp(SeapodymCoupled& sc);
void Sensitivity_analysis(SeapodymCoupled& sc, const int sftype);

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
void seapodym_threaded(void* threadid) {
    long tid;
    tid = (long)threadid;
    cout << tid << endl;

    // string dirout = get_path(parfile);
    char* parfile = "./initparfile.xml";
    //   char* parfile =
    //   "/homelocal/isenina/SEAPODYM/run/YFT/ECCO/LikeComplete/I_tests/threads/initparfile.xml";

    //   std::ostringstream ostr;
    //   ostr << "s" << tid;

    // string dir = get_path(parfile);
    // chdir(dir.c_str());

    // string simdir = ostr.str();
    // mkdir(simdir.c_str(),0777);
    // chdir(simdir.c_str());

    new_thread_data* tptr = (new_thread_data*)threadid;

    time_t time_sec;
    time(&time_sec);
    const time_t time0 = time_sec;

    // initialize stack sizes for derivative information storage
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    //-----MODE OPTIMIZATION--------------
    long int gradstack_buffer = 90000000L;
    long int cmpdif_buffer = 10000000000L;
    long int gs_var_buffer = 140000000L;
    gradient_structure::set_GRADSTACK_BUFFER_SIZE(gradstack_buffer);
    gradient_structure::set_CMPDIF_BUFFER_SIZE(cmpdif_buffer);
    gradient_structure gs(gs_var_buffer);

    int cmp_regime = 3;

    int FLAG = 0;

    cout << "\nstarting time: " << ctime(&time_sec) << endl;

    ad_comm::pthread_manager->set_slave_number(tptr->thread_no);

    ad_comm::pthread_manager->cread_lock_buffer(0);
    // read the independent variables IN THE SAME ORDER AS THEY ARE SENT
    int ipar = ad_comm::pthread_manager->get_int(0);
    // release the constant buffer
    ad_comm::pthread_manager->cread_unlock_buffer(0);

    ofstream wtxt;
    char* outfile = "./sa_res.txt";

    SeapodymCoupled sc(parfile);

    // values to solve for
    const int nvar = sc.nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1, nvar);
    sc.xinit(x, x_names);
    cout << "Total number of variables: " << nvar << '\n' << '\n';

    sc.OnRunFirstStep();

    clock_t time1 = clock();

    cout << "\nstarting computing sensitivities using likelihood and gradient "
            "comp"
         << endl;
    // sc.param->like_types[0] = 7;

    dvector g(1, nvar);
    g.initialize();

    double x_init = x(ipar);
    for (int n = 0; n < 4; n++) {
        if (n == 0) x(ipar) = sc.param->par_init_step_left(ipar);
        if (n == 1) x(ipar) = sc.param->par_init_step_left(ipar);
        if (n == 2) {
            x(ipar) = x_init;
            x(ipar) = sc.param->par_init_step_right(ipar);
        }
        if (n == 3) x(ipar) = sc.param->par_init_step_right(ipar);

        double func_predict = sc.run_coupled((dvar_vector)x);
        gradcalc(nvar, g);

        wtxt.open(outfile, ios::app);
        if (wtxt) {
            wtxt << ipar << "\t" << sc.param->parfile_names[ipar] << "\t"
                 << sc.param->get_parval(ipar) << "\t";
            for (int i = 1; i <= nvar; i++) wtxt << g[i] << "\t";
            wtxt << func_predict << endl;
        }
        wtxt.close();
    }

    //	ad_comm::pthread_manager->write_lock_buffer(0);
    //	ad_comm::pthread_manager->send_dvector(g(1,nvar),0);
    //	ad_comm::pthread_manager->write_unlock_buffer(0);

    time_t time2 = clock();
    double total_elapsed_time =
        (double)((time2 - time1) / CLOCKS_PER_SEC) / 60.0;
    cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;

    pthread_exit(threadid);
}

/*
1. Mp mean exp(0) = 0.39325
2. Ms mean max(0) = 0.00422477
3. Ms mean slope(0) = 0.545754
4. M mean range(0) = 2.41356
5. a sst spawning(0) = 0.867213
6. b sst spawning(0) = 22.0033
7. opt sst larvae(0) = 29.0033
8. a sst habitat(0) = 0.250068
9. b sst habitat(0) = 11.4989
10. T age size slope(0) = 0.500026
11. a oxy habitat(0) = 1.00006e-07
12. b oxy habitat(0) = 7.53171e-05
13. eF habitat epi(0) = 3.37113e-08
14. eF habitat meso(0) = 1.21456
15. eF habitat mmeso(0) = 4.07346e-06
16. eF habitat bathy(0) = 10.601009
17. eF habitat mbathy(0) = 0.715158
18. eF habitat hmbathy(0) = 0.000234151
19. sigma species(0. ) = Mp mean exp(0) = 0.0231488
20. MSS species(0) = 0.363988
21. MSS size slope(0) = 0.100037
22. c diff fish(0) = 1.10794e-06
23. nb recruitment(0) = 0.0356171
24. spawning season peak(0) = 91.3668
25. spawning season start(0) = 1.4

 */

/// 1. Option for computing likelihood projection in 2D parametric space.
void Hyperspace_projection(SeapodymCoupled& sc, dvar_vector x) {
    //	sc.param->set_gradcalc(false);
    const int Npars = sc.param->nb_varproj - 1;
    ivector ix(0, Npars);
    ix.initialize();
    int n1 = sc.param->varproj_nsteps[0];
    int n2 = sc.param->varproj_nsteps[1];
    int nmax = max(n1, n2) - 1;
    dmatrix xvalues(0, Npars, 0, nmax);
    xvalues.initialize();
    dmatrix pars(0, Npars, 0, nmax);
    pars.initialize();
    dmatrix Lproj(0, n1 - 1, 0, n2 - 1);
    Lproj.initialize();

    sc.param->get_param_index(ix, xvalues, pars);

    sc.OnRunFirstStep();
    clock_t time1 = clock();
    cout << "\nstarting hyperspace projection computation for ";
    for (int n = 0; n < sc.param->nb_varproj; n++)
        cout << sc.param->varproj[n] << " ";
    cout << endl;

    ofstream ofs;
    const char* filename = "hyperproj.out";
    ofs.open(filename, ios::out);
    for (int n = 0; n < sc.param->nb_varproj; n++)
        ofs << sc.param->varproj[n] << " ";
    ofs << "\n" << n1 << " " << n2 << "\n";
    for (int n = 0; n <= Npars; n++) {
        for (int i = 0; i < sc.param->varproj_nsteps(n); i++)
            ofs << pars(n, i) << " ";
        ofs << "\n";
    }
    ofs.close();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            cout << j + i * n2 + 1 << ": ";
            for (int n = 0; n <= Npars; n++) {
                int k;
                if (n == 0) k = i;
                if (n == 1) k = j;
                x[ix(n)] = xvalues(n, k);
                cout << pars(n, k) << " ";
            }
            Lproj(i, j) = sc.run_coupled(x);
        }
        ofs.open(filename, ios::app);
        for (int j = 0; j < n2; j++) ofs << Lproj(i, j) << " ";
        ofs << "\n";
        ofs.close();
    }

    time_t time2 = clock();
    double total_elapsed_time =
        (double)((time2 - time1) / CLOCKS_PER_SEC) / 60.0;
    cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;
    // cleanup_temporary_files();
}

/// 2. Option for computing Hessian matrix
void Hessian_comp(SeapodymCoupled& sc) {
    // values to solve for
    const int nvar = sc.nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1, nvar);
    sc.xinit(x, x_names);
    cout << "Total number of variables: " << nvar << '\n' << '\n';

    double likelihood = 0;
    double delta = 1e-6;
    double epsilon = 0.1;
    dvector g1(1, nvar);
    g1.initialize();
    dvector g2(1, nvar);
    g2.initialize();
    dvector H1(1, nvar);
    H1.initialize();
    dvector H2(1, nvar);
    H2.initialize();
    dmatrix H(1, nvar, 1, nvar);
    H.initialize();

    sc.OnRunFirstStep();

    clock_t time1 = clock();
    cout << "\nstarting computing Hessian" << endl;

    likelihood = sc.run_coupled((dvar_vector)x);
    gradcalc(nvar, g1);
    cout << "Gradient for estimated vector: \n" << g1 << endl;
    // two-point finite-difference approximation of the Hessian
    for (int ix = 1; ix <= nvar; ix++) {
        double xs = x(ix);

        x(ix) = xs + delta;
        likelihood = sc.run_coupled((dvar_vector)x);
        gradcalc(nvar, g2);
        H1 = (g2 - g1) / delta;  // 1. step correction
        // H1 = g2; //2. central differences

        g2.initialize();
        x(ix) = xs + epsilon * delta;  // 1. step correction
        // x(ix) = xs - delta; //2. central differences
        likelihood = sc.run_coupled((dvar_vector)x);
        gradcalc(nvar, g2);
        H2 = (g2 - g1) / (epsilon * delta);  // 1. step correction
        // H2 = g2; //2. central differences

        H(ix) = (H2 - epsilon * H1) / (1 - epsilon);  // 1. step correction
        // H(ix) = (H1-H2)/(2*delta); // 2. central differences

        cout << ix << ".\t" << H(ix) << endl;

        x(ix) = xs;
        g2.initialize();
    }
    // making Hessian symmetric as it's subject to truncation error
    for (int ix = 1; ix < nvar; ix++)
        for (int kx = ix + 1; kx <= nvar; kx++) {
            double av = 0.5 * (H(ix, kx) + H(kx, ix));
            H(ix, kx) = av;
            H(kx, ix) = av;
        }

    dmatrix Cov = inv(H);
    double determ = det(H);
    dvector evalues = eigenvalues(H);
    ofstream ofs;
    const char* filename = "Hessian.out";

    ofs.open(filename, ios::out);
    ofs << nvar << "\n";

    dvector pars = sc.param->get_parvals();
    ofs << "Parameter\t"
        << "Est. value\t"
        << "Gradient"
        << "\n";
    for (int i = 1; i <= nvar; i++)
        ofs << x_names[i] << "\t" << pars(i) << "\t" << g1(i) << "\n";
    ofs << "\n";

    ofs << determ << "\n";
    ofs << evalues << "\n\n";

    for (int i = 1; i <= nvar; i++) {
        for (int j = 1; j <= nvar; j++) ofs << H(i, j) << " ";
        ofs << "\n";
    }
    ofs << "\n";
    for (int i = 1; i <= nvar; i++) {
        for (int j = 1; j <= nvar; j++) ofs << Cov(i, j) << " ";
        ofs << "\n";
    }
    ofs << "\n";

    time_t time2 = clock();
    double total_elapsed_time =
        (double)((time2 - time1) / CLOCKS_PER_SEC) / 60.0;
    cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;
}

/// 3. Option for parametric sensitivity analysis
void Sensitivity_analysis(SeapodymCoupled& sc, const int sftype) {
    // values to solve for
    const int nvar = sc.nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1, nvar);
    sc.xinit(x, x_names);
    cout << "Total number of variables: " << nvar << '\n' << '\n';

    sc.OnRunFirstStep();

    dvector s(1, nvar);
    s.initialize();

    clock_t time1 = clock();

    if (sftype == 0) {
        // sc.param->like_types[0] = 7;
        // cout << "\nstarting computing sensitivities using model predictions
        // only" << endl;

        cout << "\nstarting computing sensitivities using likelihood and "
                "gradient comp"
             << endl;

        dvector g(1, nvar);
        g.initialize();

        double func_predict = sc.run_coupled((dvar_vector)x);
        gradcalc(nvar, g);
        s = g;  // func_predict;
    } else {
        cout << "\nstarting computing sensitivities using likelihood" << endl;

        const double eps = 1e-3;
        gradient_structure::set_NO_DERIVATIVES();
        cout << "At current parameters..." << endl;
        double like_opt = sc.run_coupled((dvar_vector)x);
        cout << "At the boundaries..." << endl;
        for (int i = 1; i <= nvar; i++) {
            double xs = x(i);
            x(i) = sc.param->par_init_lo(i, eps);
            double like = sc.run_coupled((dvar_vector)x);
            double s_lo = like / like_opt - 1.0;

            x(i) = sc.param->par_init_up(i, eps);
            like = sc.run_coupled((dvar_vector)x);
            double s_up = like / like_opt - 1.0;

            x(i) = xs;
            s(i) = max(abs(s_lo), abs(s_up));
            cout << i << " " << x_names[i] << " " << s(i) << endl;
        }
    }

    cout << endl
         << "N \t"
         << "parameter \t"
         << "\trelative sensitivity" << endl;
    cout << "-------------------------------------------------------" << endl;

    for (int i = 1; i <= nvar; i++) {
        int l = length(x_names[i]);
        string tab = "\t";
        if (l < 16) tab += "\t";
        if (l < 7) tab += "\t";
        cout << i << " \t" << x_names[i] << tab << s[i] << endl;
    }
    cout << "-------------------------------------------------------" << endl;

    time_t time2 = clock();
    double total_elapsed_time =
        (double)((time2 - time1) / CLOCKS_PER_SEC) / 60.0;
    cout << "\ntotal time: " << total_elapsed_time << " minutes" << endl;
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
