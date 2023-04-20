#include "calpop.h"

/// Main function with memory control and adjoint functions for:
/// predicting catch by fishery based on fishing effort.
/// Forward function is in predicted_catch.cpp

void dv_predicted_catch_fishery();
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);
const double sq_degree =
    12347.65;  // sq.km in sq.degree (without lat_correction)

void CCalpop::Predicted_Catch_Fishery(
    const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw,
    const int sp, const int f, const int k, const int year, const int month,
    const int t_count, const int step_count) {
    mat.dvarCatch_est[sp][k].initialize();
    mat.catch_est[sp][k].initialize();

    const int a0 = param.sp_a0_adult[sp];
    const int nb_ages = param.sp_nb_cohorts[sp];
    const int nb_region = param.nb_region_sp_B[sp];

    // Acoef (matrix of coef to compute densities on fishing area):
    // dimension of matrix
    //	const int reso = param.fishery_reso(f);
    //	m = reso*60.0/param.deltaX+2;
    //	n = reso*60.0/param.deltaY+2;
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
        mat.dvarCatch_est[sp][k].save_dvar_matrix_position();
        save_identifier_string2((char*)"before_predicted_catch_comp");

        predicted_catch_fishery_comp(
            map, param, mat, f, k, sp, age, value(mat.dvarDensity(sp, age)),
            step_count);
        for (int r = 0; r < nb_region; r++)
            mat.dvarLF_est[sp][age][k].elem_value(r) =
                mat.C_N_sp_age_fishery[sp][age][k][r];

        save_identifier_string2((char*)"Predicted_Catch_begin");
        mat.dvarCatch_est(sp, k).save_dvar_matrix_position();
        mat.dvarLF_est(sp, age, k).save_dvar_vector_position();
        mat.dvarDensity(sp, age).save_dvar_matrix_position();
        Q.save_dvar_matrix_position();
        mat.dvarsV.save_dvar_matrix_position();
        mat.dvarsU.save_dvar_matrix_position();
        Sasympt.save_dvar_matrix_position();
        save_int_value(year);
        save_int_value(month);
        save_int_value(t_count);
        save_int_value(step_count);
        save_int_value(sp);
        save_int_value(age);
        save_int_value(f);
        save_int_value(k);
        save_int_value(nb_region);
        unsigned long int pop = (unsigned long int)this;
        save_long_int_value(pop);
        unsigned long int cmat = (unsigned long int)&mat;
        save_long_int_value(cmat);
        unsigned long int pmap = (unsigned long int)&map;
        save_long_int_value(pmap);
        unsigned long int cparam = (unsigned long int)&param;
        save_long_int_value(cparam);
        unsigned long int crw = (unsigned long int)&rw;
        save_long_int_value(crw);
        save_identifier_string2((char*)"Predicted_Catch_end");

        gradient_structure::GRAD_STACK1->set_gradient_stack(
            dv_predicted_catch_fishery);
    }
}

void dv_predicted_catch_fishery() {
    verify_identifier_string2((char*)"Predicted_Catch_end");
    unsigned long int pos_rw = restore_long_int_value();
    unsigned long int pos_param = restore_long_int_value();
    unsigned long int pos_map = restore_long_int_value();
    unsigned long int pos_mat = restore_long_int_value();
    unsigned long int pos_pop = restore_long_int_value();
    unsigned nb_region = restore_int_value();
    unsigned k = restore_int_value();
    unsigned f = restore_int_value();
    unsigned age = restore_int_value();
    unsigned sp = restore_int_value();
    unsigned nstep = restore_int_value();
    unsigned t_count = restore_int_value();
    unsigned month = restore_int_value();
    unsigned year = restore_int_value();
    const dvar_matrix_position sasympt_pos = restore_dvar_matrix_position();
    const dvar_matrix_position slength_pos = restore_dvar_matrix_position();
    const dvar_matrix_position sslope_pos = restore_dvar_matrix_position();
    const dvar_matrix_position q_pos = restore_dvar_matrix_position();
    const dvar_matrix_position uu_pos = restore_dvar_matrix_position();
    const dvar_vector_position lf_pos = restore_dvar_vector_position();
    const dvar_matrix_position catch_pos = restore_dvar_matrix_position();
    verify_identifier_string2((char*)"Predicted_Catch_begin");

    verify_identifier_string2((char*)"before_predicted_catch_comp");
    const dvar_matrix_position c_pr_pos = restore_dvar_matrix_position();

    CCalpop* pop = (CCalpop*)pos_pop;
    CMatrices* mat = (CMatrices*)pos_mat;
    PMap* map = (PMap*)pos_map;
    CParam* param = (CParam*)pos_param;
    CReadWrite* rw = (CReadWrite*)pos_rw;

    const int s_func_type = param->s_func_type[f];

    const int imax = map->imax;
    const int imin = map->imin;
    const int jmax = map->jmax;
    const int jmin = map->jmin;

    dvector lat_correction(jmin, jmax);
    lat_correction.initialize();

    dmatrix effort(imin, imax, map->jinf, map->jsup);
    dmatrix efflon(imin, imax, map->jinf, map->jsup);
    dmatrix efflat(imin, imax, map->jinf, map->jsup);
    // effort.initialize();
    rw->get_effort_lonlat(*param, effort, efflon, efflat, f, year, month);

    const int reso = param->fishery_reso(f);
    int m = reso * 60.0 / param->deltaX + 2;
    int n = reso * 60.0 / param->deltaY + 2;
    dmatrix af(0, m - 1, 0, n - 1);
    af.initialize();

    dmatrix uu(map->imin1, map->imax1, map->jinf1, map->jsup1);
    uu = mat->density_before(sp, t_count, age);

    for (int j = jmin; j <= jmax; j++) {
        double lastlat = param->lastlat(j);
        lat_correction[j] = param->correction_lat(lastlat);
    }

    dmatrix dfuu = restore_dvar_matrix_derivatives(uu_pos);
    dmatrix dfC_pred = restore_dvar_matrix_derivatives(catch_pos);
    dmatrix dfC_pr_age = restore_dvar_matrix_derivatives(c_pr_pos);
    dmatrix dfQ = restore_dvar_matrix_derivatives(q_pos);
    dmatrix dfSslope = restore_dvar_matrix_derivatives(sslope_pos);
    dmatrix dfSlength = restore_dvar_matrix_derivatives(slength_pos);
    dmatrix dfSasympt = restore_dvar_matrix_derivatives(sasympt_pos);
    dvector dfLF = restore_dvar_vector_derivatives(lf_pos);

    const int catch_units = param->fishery_catch_units[f];
    const double weight = pow(1e-3 * param->weight[sp][age], catch_units);
    const double w_area = weight * sq_degree;

    const double q_dyn = 1 + param->q_dyn_fishery[f] * nstep;
    const double catchability = param->q_sp_fishery[sp][k] * q_dyn;
    // const double selectivity = param->selectivity_comp(sp,age,f,k);
    // const double selectivity = param->selectivity_comp(sp,age,f,k,nstep);
    const double selectivity = pop->Selectivity(sp, f, age);
    const double s_c = selectivity * catchability;

    double dflength, dfslope, dfasympt;
    param->dfselectivity(dfslope, dflength, dfasympt, sp, age, f, k);
    const double cell_area_deg = param->deltaX * param->deltaY / (60.0 * 60.0);

    for (int i = imax; i >= imin; i--) {
        const int jmin = map->jinf[i];
        const int jmax = map->jsup[i];
        for (int j = jmax; j >= jmin; j--) {
            if (map->carte(i, j)) {
                // const double eff = effort(i,j);
                double eff = effort(i, j);
                if (eff) {
                    // precompute biomass of fish aggregated to the resolution
                    // of fishing data
                    double fish = 0;
                    double dffish = 0;
                    int ki = 0;
                    int kj = 0;
                    param->afcoef(efflon(i, j), efflat(i, j), af, ki, kj, f);

                    double afr = 0.0;
                    for (int ii = 0; ii < m; ii++)
                        for (int jj = 0; jj < n; jj++)
                            if (af(ii, jj) && map->carte(ii + ki, jj + kj))
                                afr += af(ii, jj) / cell_area_deg;

                    for (int ii = 0; ii < m; ii++)
                        for (int jj = 0; jj < n; jj++)
                            if (af(ii, jj) && map->carte(ii + ki, jj + kj))
                                fish += af(ii, jj) * uu(ii + ki, jj + kj) /
                                        (afr * lat_correction[jj + kj]);

                    // const double F_pr = s_c*eff/(reso*reso);
                    // const double F = reso*reso*param->func_limit_one(F_pr);
                    const double F_pr = s_c * eff / afr;
                    const double F = afr * param->func_limit_one(F_pr);
                    // const double F =
                    // reso*reso*param->func_limit_one(s_c*eff/(reso*reso)); end
                    // of recomputation section

                    double dfF = 0;
                    double dfLF_pr_reg = 0;
                    for (int a = nb_region - 1; a >= 0; a--) {
                        int r = param->area_sp_B[sp][a] - 1;
                        if (param->nb_region) {
                            if ((i >= map->regimin[a]) &&
                                (i < map->regimax[a]) &&
                                (j >= map->regjmin[a]) &&
                                (j < map->regjmax[a])) {
                                // mat.C_N_sp_age_fishery[sp][age][k][r] =
                                // LF_pr_reg + F * fish * sq_degree;
                                dffish += F * sq_degree * dfLF(r);
                                dfF += fish * sq_degree * dfLF(r);
                                dfLF_pr_reg += dfLF(r);
                                dfLF(r) = 0;
                            }
                        } else {
                            // mat.C_N_sp_age_fishery[sp][age][k][0] += F * fish
                            // * sq_degree;
                            dffish += F * sq_degree * dfLF(0);
                            dfF += fish * sq_degree * dfLF(r);
                            dfLF_pr_reg += dfLF(0);
                            dfLF(0) = 0;
                        }
                        // LF_pr_reg = LF(r);
                        dfLF(r) += dfLF_pr_reg;
                        dfLF_pr_reg = 0;
                    }

                    // mat.dvarCatch_est[sp][k][i][j] = C_pr_age +
                    // F*fish*w_area;
                    dffish += F * w_area * dfC_pred(i, j);
                    dfF += fish * w_area * dfC_pred(i, j);
                    dfC_pr_age(i, j) += dfC_pred(i, j);
                    dfC_pred(i, j) = 0.0;

                    // const double F = afr*param.func_limit_one(s_c*eff/afr);
                    // double dfs_c =
                    // eff*param->dffunc_limit_one(F_pr,reso*reso*dfF)/(reso*reso);
                    double dfs_c =
                        eff * param->dffunc_limit_one(F_pr, afr * dfF) / afr;
                    dfF = 0.0;

                    /*
                    //const double F =
                    reso*reso*param.func_limit_one(s_c*eff/(reso*reso)); double
                    dfs_c = eff*param->dffunc_limit_one(F_pr)*dfF; dfF = 0.0;
                    */

                    // const double s_c = selectivity*catchability;
                    dfQ(i, j) += q_dyn * selectivity * dfs_c;
                    dfSslope(i, j) += dfslope * catchability * dfs_c;
                    if (s_func_type > 1)
                        dfSlength(i, j) += dflength * catchability * dfs_c;
                    if (s_func_type == 3)
                        dfSasympt(i, j) += dfasympt * catchability * dfs_c;
                    dfs_c = 0.0;

                    for (int jj = n - 1; jj >= 0; jj--)
                        for (int ii = m - 1; ii >= 0; ii--)
                            if (af(ii, jj) && map->carte(ii + ki, jj + kj))
                                // fish +=
                                // af(ii,jj)*uu(ii+ki,jj+kj)/lat_correction[jj+kj];
                                dfuu(ii + ki, jj + kj) +=
                                    af(ii, jj) * dffish /
                                    (afr * lat_correction(jj + kj));

                    dffish = 0;
                }
            }
        }
    }

    dfuu.save_dmatrix_derivatives(uu_pos);
    dfC_pred.save_dmatrix_derivatives(catch_pos);
    dfC_pr_age.save_dmatrix_derivatives(c_pr_pos);
    dfQ.save_dmatrix_derivatives(q_pos);
    dfSslope.save_dmatrix_derivatives(sslope_pos);
    dfSlength.save_dmatrix_derivatives(slength_pos);
    dfSasympt.save_dmatrix_derivatives(sasympt_pos);
    dfLF.save_dvector_derivatives(lf_pos);
}
