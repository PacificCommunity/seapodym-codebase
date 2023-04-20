#include "calpop.h"

/// Forward main function called in simulation mode only for:
/// computing the sum of natural and fishing mortalities
/// See total_mortality_comp.cpp

void CCalpop::Precalrec_total_mortality_comp(
    const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw,
    dvar_matrix& mortality, const int age, const int sp, const int t_count,
    const int year, const int month, const int step_count) {
    int k = 0;
    const int nb_fishery = param.get_nbfishery();
    dvar_matrix Q, Sslope;
    Q.allocate(map.imin, map.imax, map.jinf, map.jsup);
    Q.initialize();
    Sslope.allocate(map.imin, map.imax, map.jinf, map.jsup);
    Sslope.initialize();

    dmatrix effort;
    effort.allocate(map.imin, map.imax, map.jinf, map.jsup);
    effort.initialize();

    for (int f = 0; f < nb_fishery; f++) {
        if (param.mask_fishery_sp[sp][f]) {
            // 2014: no use of effort for fisheries with no or bad effort data:
            if (param.mask_fishery_sp_no_effort[sp][f]) {
                k++;
                continue;
            }

            double catchability = param.q_sp_fishery[sp][k] *
                                  (1.0 + param.q_dyn_fishery[f] * step_count);
            if (catchability < 0)
                catchability = 0;  // in case if negative trend has been chosen
            const double selectivity = Selectivity(sp, f, age);
            double sq = selectivity * catchability;

            Q = param.dvarsQ_sp_fishery[sp][k];
            Sslope = param.dvarsSslope_sp_fishery[sp][k];
            mat.dvarsU.initialize();
            mat.dvarsV.initialize();
            if (param.s_func_type[f] > 1) {
                // Slength:
                mat.dvarsU = param.dvarsSlength_sp_fishery[sp][k];
                if (param.s_func_type[f] == 3)
                    // Sasympt:
                    mat.dvarsV = param.dvarsSasympt_sp_fishery[sp][k];
            }
            if (!param.fdata_rm) {
                // effort was redistributed here using param.afcoef routine:
                rw.get_effort_rm(param, effort, f, year, month);
            } else {
                if (param.mpa_simulation)
                    effort = mat.effort[f];
                else
                    // to use GMB effort on model resolution resolution (old
                    // configuration)
                    rw.get_effort(param, effort, f, year, month);
            }

            precalrec_total_mortality_comp(
                map.carte, effort, mortality, sq, mat.lat_correction);

            k++;
        }
    }
}
