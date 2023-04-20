#include "calpop.h"

/// Forward main function called in simulation mode only for:
/// computing the total (sum over fisheries without effort) observed catch at
/// age. See total_obs_catch_age.cpp

void CCalpop::Total_obs_catch_age_comp(
    const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw,
    const int age, const int sp, const int year, const int month,
    const int t_count) {
    mat.dvarCtot_age_obs(sp, age).initialize();
    const int nb_fishery = param.get_nbfishery();
    const int nb_region = param.nb_region;
    const double area = 1.852 * param.deltaX * 1.852 * param.deltaY;
    const double W_age_mt = param.weight[sp][age] * 0.001;

    dvar_matrix Sslope, Slen, Sasympt;
    Sslope.allocate(map.imin, map.imax, map.jinf, map.jsup);
    Sslope.initialize();
    Slen.allocate(map.imin, map.imax, map.jinf, map.jsup);
    Slen.initialize();
    Sasympt.allocate(map.imin, map.imax, map.jinf, map.jsup);

    int k = 0;
    int fne = 0;  // fisheries no effort counter
    for (int f = 0; f < nb_fishery; f++) {
        if (param.mask_fishery_sp[sp][f]) {
            if (param.mask_fishery_sp_no_effort[sp][f]) {
                // catch_units=0 if catch is in nb; =1 if in metric tons
                const int catch_units = param.fishery_catch_units[f];
                const double C2Dunits =
                    1.0 / (pow(W_age_mt, catch_units) * area);

                Sslope = param.dvarsSslope_sp_fishery[sp][k];
                if (param.s_func_type[f] > 1) {
                    // Slength:
                    // mat.dvarsU = param.dvarsSlength_sp_fishery[sp][k];
                    Slen = param.dvarsSlength_sp_fishery[sp][k];
                    if (param.s_func_type[f] == 3) {
                        // Sasympt:
                        Sasympt = param.dvarsSasympt_sp_fishery[sp][k];
                    }
                }

                total_obs_catch_age_comp(
                    map, param, mat, value(mat.dvarDensity(sp, age)),
                    mat.catch_obs[sp][k], mat.dvarCtot_age_obs(sp, age), f, fne,
                    k, age, sp, C2Dunits);

                for (int r = 0; r < nb_region; r++)
                    mat.dvarLF_est[sp][age][k].elem_value(r) =
                        mat.C_N_sp_age_fishery[sp][age][k][r];

                fne++;
            }
            k++;
        }
    }
}
