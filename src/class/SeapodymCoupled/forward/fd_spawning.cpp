#include <fvar.hpp>

#include "SeapodymCoupled.h"

/// Forward main function called in simulation mode only for:
/// computing the new larval biomass in two regimes of the simulation:
/// during climatological spinup and after the population has been built.
/// See spawning.cpp

void SeapodymCoupled::Spawning(
    dvar_matrix& J, dvar_matrix& Hs, dvar_matrix& Nmature, const int jday,
    const int sp, const int pop_built, const int t_count) {
    J.initialize();

    dmatrix J_c = value(J);
    dmatrix Hs_c = value(Hs);
    dmatrix N_mat = value(Nmature);
    dvariable nb_recruitment = param->dvarsNb_recruitment[sp];

    dvar_matrix Nbr(map.imin, map.imax, map.jinf, map.jsup);
    Nbr = nb_recruitment;

    if (pop_built) {
        dvariable a_adults_spawning = param->dvarsA_adults_spawning[sp];
        dvar_matrix A_sp(map.imin, map.imax, map.jinf, map.jsup);
        A_sp = a_adults_spawning;

        spawning_built_comp(
            J_c, Hs_c, N_mat, value(nb_recruitment), value(a_adults_spawning));

        J = nograd_assign(J_c);

    } else {
        dvariable mu, sigma;
        if (!param->uncouple_sst_larvae[sp]) {
            mu = param->dvarsB_sst_spawning[sp];
            sigma = param->dvarsA_sst_spawning[sp];
        } else {
            mu = param->dvarsB_sst_larvae[sp];
            sigma = param->dvarsA_sst_larvae[sp];
        }
        dvar_matrix Sigma(map.imin, map.imax, map.jinf, map.jsup);
        dvar_matrix Mu(map.imin, map.imax, map.jinf, map.jsup);
        Sigma = sigma;
        Mu = mu;

        dvariable Tmin = mu - 2.0 * sigma;
        spawning_spinup_comp(
            J_c, Hs_c, mat.tempn[t_count][0], value(nb_recruitment),
            value(Tmin));

        J = nograd_assign(J_c);
    }
}
