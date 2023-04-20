#include "SeapodymCoupled.h"

/// Forward main functions called in simulation mode only for:
/// complex non-linear function that depends on the density of adult biomass
/// to derive the food requirement and compute the index of depletion in case
/// if the available micronekton does not support the daily ration.
/// See food_requirement_index.cpp

void SeapodymCoupled::Food_Requirement_Index(
    dvar_matrix& IFR, dvar_matrix FR_pop, dvar_matrix ISR_denom, const int sp,
    const int age, const int t_count, const int jday) {
    IFR.initialize();

    IFR_age_comp(IFR, value(FR_pop), value(ISR_denom), age, sp, t_count);
}

// Total food requirement function, identical to Spawning_Biomass_comp, except
// for multiplier
void SeapodymCoupled::FR_pop_comp(dvar_matrix& FR_pop, const int sp) {
    FR_pop.initialize();
    const double R = param->forage_ration[sp];

    for (int a = param->sp_a0_adult[sp]; a < param->sp_nb_cohorts[sp]; a++) {
        const double W = param->weight[sp][a] * 0.001;  // in mt
        const double wrt = W * R * param->deltaT;

        for (int i = map.imin; i <= map.imax; i++) {
            const int jmin = map.jinf[i];
            const int jmax = map.jsup[i];
            for (int j = jmin; j <= jmax; j++) {
                if (map.carte(i, j)) {
                    double FR_pr = value(FR_pop(i, j));
                    FR_pop.elem_value(i, j) =
                        FR_pr + wrt * value(mat.dvarDensity(sp, a, i, j));
                }
            }
        }
    }
}

void SeapodymCoupled::ISR_denom_comp(
    dvar_matrix& ISR_denom, const int sp, const int t_count) {
    dvar_matrix dvarF_access_sum_age(map.imin, map.imax, map.jinf, map.jsup);

    ISR_denom = 1.0;
    for (int n = 0; n < nb_forage; n++) {
        mat.F_access_sum_age(sp, n, t_count).initialize();
        dvarF_access_sum_age.initialize();
        for (int a = param->sp_a0_adult[sp]; a < param->sp_nb_cohorts[sp];
             a++) {
            for (int i = map.imin; i <= map.imax; i++) {
                const int jmin = map.jinf[i];
                const int jmax = map.jsup[i];
                for (int j = jmin; j <= jmax; j++) {
                    if (map.carte(i, j)) {
                        mat.F_access_sum_age(sp, n, t_count, i, j) +=
                            value(mat.dvarF_access(n, a, i, j));
                        dvarF_access_sum_age.elem_value(i, j) =
                            mat.F_access_sum_age(sp, n, t_count, i, j);
                    }
                }
            }
        }

        for (int i = map.imin; i <= map.imax; i++) {
            const int jmin = map.jinf[i];
            const int jmax = map.jsup[i];
            for (int j = jmin; j <= jmax; j++) {
                if (map.carte(i, j)) {
                    double ISR_denom_pr = value(ISR_denom(i, j));
                    ISR_denom.elem_value(i, j) =
                        ISR_denom_pr *
                        mat.F_access_sum_age(sp, n, t_count, i, j);
                }
            }
        }
    }
}
