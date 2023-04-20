#include "calpop.h"

/// Forward main function called in simulation mode only for:
/// predicting catch by fishery based on fishing effort.
/// See predicted_catch.cpp

void CCalpop::Predicted_Catch_Fishery(
    const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw,
    const int sp, const int f, const int k, const int year, const int month,
    const int t_count, const int step_count) {
    mat.dvarCatch_est[sp][k].initialize();
    mat.catch_est[sp][k].initialize();

    const int a0 = param.sp_a0_adult[sp];
    const int nb_ages = param.sp_nb_cohorts[sp];
    const int nb_region = param.nb_region_sp_B[sp];

    // Q:
    dvar_matrix Q, Sasympt;
    Q.allocate(map.imin, map.imax, map.jinf, map.jsup);
    Q = param.dvarsQ_sp_fishery[sp][k];
    // Sslope:
    mat.dvarsU.initialize();
    mat.dvarsV.initialize();
    mat.dvarsV = param.dvarsSslope_sp_fishery[sp][k];
    Sasympt.allocate(map.imin, map.imax, map.jinf, map.jsup);
    Sasympt.initialize();
    if (param.s_func_type[f] > 1) {
        // Slength:
        mat.dvarsU = param.dvarsSlength_sp_fishery[sp][k];
        if (param.s_func_type[f] == 3)
            Sasympt = param.dvarsSasympt_sp_fishery[sp][k];
    }
    for (int age = a0; age < nb_ages; age++) {
        predicted_catch_fishery_comp(
            map, param, mat, f, k, sp, age, value(mat.dvarDensity(sp, age)),
            step_count);
        for (int r = 0; r < nb_region; r++)
            mat.dvarLF_est[sp][age][k].elem_value(r) =
                mat.C_N_sp_age_fishery[sp][age][k][r];
    }
}
