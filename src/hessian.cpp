#include <fvar.hpp>

#include "SeapodymCoupled.h"

// void Hessian_comp(const char* parfile);
double run_model(
    SeapodymCoupled& sc, dvar_vector x, dvector& g, const int nvar);

/// 2. Option for computing Hessian matrix
void Hessian_comp(const char* parfile) {
    SeapodymCoupled sc(parfile);

    // values to solve for
    const int nvar = sc.nvarcalc();
    independent_variables x(1, nvar);
    adstring_array x_names(1, nvar);
    sc.xinit(x, x_names);
    cout << "Total number of variables: " << nvar << '\n' << '\n';

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

    double likelihood = 0.0;

    clock_t time1 = clock();
    cout << "\nstarting computing Hessian" << endl;
    likelihood = run_model(sc, x, g1, nvar);

    cout << "Likelihood and Gradient for estimated vector: \n"
         << likelihood << "; " << g1 << endl;
    // two-point finite-difference approximation of the Hessian
    for (int ix = 1; ix <= nvar; ix++) {
        double xs = x(ix);

        x(ix) = xs + delta;
        likelihood = run_model(sc, x, g2, nvar);

        H1 = (g2 - g1) / delta;

        g2.initialize();

        x(ix) = xs + epsilon * delta;  // 1. step correction
        likelihood = run_model(sc, x, g2, nvar);

        H2 = (g2 - g1) / (epsilon * delta);  // 1. step correction

        H(ix) = (H2 - epsilon * H1) / (1 - epsilon);  // 1. step correction

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
