#include "calpop.h"

/// Forward function for:
/// predicting catch by fishery based on fishing effort.

const double sq_degree =
    12347.65;  // sq.km in sq.degree (without lat_correction)

void CCalpop::predicted_catch_fishery_comp(
    const PMap& map, CParam& param, VarMatrices& mat,
    /* dmatrix total_mort, */ const int f, const int k, const int sp,
    const int age, const dmatrix& uu, const int step_count) {
    double catchability =
        param.q_sp_fishery[sp][k] * (1.0 + param.q_dyn_fishery[f] * step_count);
    if (catchability < 0)
        catchability = 0;  // in case if negative trend has been chosen
    // const double selectivity = param.selectivity_comp(sp,age,f,k);
    // const double selectivity = param.selectivity_comp(sp,age,f,k,step_count);
    const double selectivity = Selectivity(sp, f, age);

    // catch_units=0 if catch is in nb; =1 if in metric tons
    const int catch_units = param.fishery_catch_units[f];
    // Note, variable weight is either =real weight in mt or =1 if predicted
    // catch is in numbers
    const double weight = pow(1e-3 * param.weight[sp][age], catch_units);

    // compute catch in tonnes per cell Nb/km^2 * mt =
    // W(mt)*dx*dy*(1.852)^2/cell to translate length frequencies into number of
    // fish caugth in cell: Nb/km^2 = Nb*(1.852)^2*dx*dy/cell

    const int nb_region = param.nb_region;
    const int nb_region_sp = param.nb_region_sp_B[sp];
    const double s_c = selectivity * catchability;
    const double w_area = weight * sq_degree;
    dmatrix efflon(map.imin, map.imax, map.jinf, map.jsup);
    efflon = mat.efflon(f);
    dmatrix efflat(map.imin, map.imax, map.jinf, map.jsup);
    efflat = mat.efflat(f);
    const double cell_area_deg = param.deltaX * param.deltaY / (60 * 60);
    const int reso = param.fishery_reso(f);
    // Catch equation: we will compute and fit predicted catch at the resolution
    // of the fishing data assume we have E(m,n) at coarser resolution than
    // ours: rx>=deltaX, ry>=deltaY C(m,n) = q*sum_ij(E(i,j)*B(i,j)) = [if
    // E(i,j)=const=E(m,n)/(rx*ry)] = q*E/(rx*ry)*sum_ij(B(i,j)) where rx*ry is
    // the size of fisheries data cell in deltaX and deltaY correspondingly
    //(e.g. if deltaX=deltaY=2deg and fishery_reso=5deg => rx*ry = (5/2)*(5/2)
    //= 6.25)

    // note: need to make it global variable
    int m = reso * 60.0 / param.deltaX + 2;
    int n = reso * 60.0 / param.deltaY + 2;
    dmatrix af(0, m - 1, 0, n - 1);

    for (int i = map.imin; i <= map.imax; i++) {
        const int jmin = map.jinf[i];
        const int jmax = map.jsup[i];
        for (int j = jmin; j <= jmax; j++) {
            if (map.carte(i, j)) {
                // const double effort = mat.effort[f][i][j];
                double effort = mat.effort[f][i][j];
                if (effort) {
                    double fish = 0;

                    int ki = 0;
                    int kj = 0;
                    // returns fishing area in the cell to be used in biomass
                    // computation
                    param.afcoef(efflon(i, j), efflat(i, j), af, ki, kj, f);

                    // the factor 1/(rx*rw) = sum(af)/cell_area_deg
                    double afr = 0.0;
                    for (int ii = 0; ii < m; ii++)
                        for (int jj = 0; jj < n; jj++)
                            if (af(ii, jj) && map.carte(ii + ki, jj + kj))
                                afr += af(ii, jj) / cell_area_deg;

                    // Nb. of fish in a given (lon,lat)+/-freso/2 =
                    // cell_area*sum_ij(af*N(i,j)):
                    for (int ii = 0; ii < m; ii++)
                        for (int jj = 0; jj < n; jj++)
                            if (af(ii, jj) && map.carte(ii + ki, jj + kj)) {
                                fish += af(ii, jj) * uu(ii + ki, jj + kj) /
                                        (afr * mat.lat_correction[jj + kj]);
                            }

                    const double F =
                        afr * param.func_limit_one(effort * s_c / afr);

                    // total catch in weight (tonnes)
                    mat.dvarCatch_est(sp, k).elem_value(i, j) +=
                        F * fish * w_area;
                    mat.catch_est[sp][k][i][j] += F * fish * w_area;

                    // Predicted Catch-at-age by species, fleet and region
                    for (int a = 0; a < nb_region_sp; a++) {
                        // attention: area_sp_B should not contain zero IDs!
                        int r = param.area_sp_B[sp][a] - 1;
                        if (nb_region) {
                            if ((i >= map.regimin[a]) && (i < map.regimax[a]) &&
                                (j >= map.regjmin[a]) && (j < map.regjmax[a]))
                                mat.C_N_sp_age_fishery[sp][age][k][r] +=
                                    F * fish * sq_degree;
                        } else
                            mat.C_N_sp_age_fishery[sp][age][k][0] +=
                                F * fish * sq_degree;
                    }
                }
            }
        }
    }
}
