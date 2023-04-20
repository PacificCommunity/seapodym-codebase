#include "calpop.h"

/// Forward functions for:
/// computing the total (sum over fisheries without effort) observed catch at
/// age.

void CCalpop::total_obs_catch_age_comp(
    const PMap& map, const CParam& param, CMatrices& mat, const dmatrix& uu,
    const dmatrix& Cobs, dvar_matrix& Ctot_age_obs, const int f, const int fne,
    const int k, const int age, const int sp, const double C2Dunits) {
    const int nb_region = param.nb_region;
    const int nb_region_sp = param.nb_region_sp_B[sp];

    for (int j = b.rowmin(); j <= b.rowmax(); j++) {
        for (int i = b[j].indexmin(); i <= b[j].indexmax(); i++) {
            double EB_tot = dvarsSNsum(fne).elem_value(i, j);
            if (map.carte(i, j) && EB_tot) {
                double exploited_age_class = Selectivity(sp, f, age) * uu(i, j);
                double LF_est = exploited_age_class / EB_tot;

                // 1. compute total observed catch by age
                Ctot_age_obs.elem_value(i, j) +=
                    LF_est * Cobs(i, j) * C2Dunits * mat.lat_correction(j);

                // 2. compute predicted LF by species, fleet and region
                for (int a = 0; a < nb_region_sp; a++) {
                    // attention: area_sp_B should not contain zero IDs!
                    int r = param.area_sp_B[sp][a] - 1;
                    if (nb_region) {
                        if ((i >= map.regimin[a]) && (i < map.regimax[a]) &&
                            (j >= map.regjmin[a]) && (j < map.regjmax[a]))
                            // mat.C_N_sp_age_fishery[sp][age][k][r] += LF_est;
                            mat.C_N_sp_age_fishery[sp][age][k][r] +=
                                exploited_age_class * Cobs(i, j);
                    } else
                        // mat.C_N_sp_age_fishery[sp][age][k][0] += LF_est;
                        mat.C_N_sp_age_fishery[sp][age][k][0] +=
                            exploited_age_class * Cobs(i, j);
                }
            }
        }
    }
}

void CCalpop::Recomp_total_obs_catch_age(
    const PMap& map, CParam& param, CMatrices& mat, CReadWrite& rw,
    dmatrix& Ctot_age_obs, const int age, const int sp, const int year,
    const int month, const int t_count) {
    const int nb_fishery = param.get_nbfishery();
    const double area = 1.852 * param.deltaX * 1.852 * param.deltaY;
    const double W_age_mt = param.weight[sp][age] * 0.001;

    dmatrix EB(map.imin, map.imax, map.jinf, map.jsup);
    EB.initialize();

    dmatrix catch_obs;
    catch_obs.allocate(map.imin, map.imax, map.jinf, map.jsup);
    catch_obs.initialize();
    Ctot_age_obs.initialize();

    dmatrix uu(map.imin1, map.imax1, map.jinf1, map.jsup1);
    uu = mat.density_before(sp, t_count, age);

    int k = 0;
    for (int f = 0; f < nb_fishery; f++) {
        if (param.mask_fishery_sp[sp][f]) {
            if (param.mask_fishery_sp_no_effort[sp][f]) {
                rw.get_catch(param, catch_obs, f, year, month, sp);
                Recomp_total_exploited_biomass(
                    map, param, mat, EB, catch_obs, Selectivity(sp, f), f, sp,
                    t_count);
                // catch_units=0 if catch is in nb; =1 if in metric tons
                const int catch_units = param.fishery_catch_units[f];
                const double C2Dunits =
                    1.0 / (pow(W_age_mt, catch_units) * area);

                for (int i = map.imin; i <= map.imax; i++) {
                    const int jmin = map.jinf[i];
                    const int jmax = map.jsup[i];
                    for (int j = jmin; j <= jmax; j++) {
                        if (map.carte(i, j) && EB(i, j)) {
                            double LF_est =
                                Selectivity(sp, f, age) * uu(i, j) / EB(i, j);

                            Ctot_age_obs(i, j) += LF_est * catch_obs(i, j) *
                                                  C2Dunits *
                                                  mat.lat_correction(j);
                        }
                    }
                }
            }
            k++;
        }
    }
}
