#include "calpop.h"

/// Main function with memory control and adjoint functions for:
/// predicting catch by fishery without using the effort data.
/// Forward functions are in predicted_catch_without_effort.cpp

void dv_predicted_catch_without_effort_fishery();
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);

void CCalpop::Predicted_Catch_Fishery_no_effort(
    const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw,
    const int sp, const int year, const int month) {
    const int nb_fishery = param.get_nbfishery();
    const int a0 = param.sp_a0_adult[sp];
    const int nb_ages = param.sp_nb_cohorts[sp];

    int k = 0;
    for (int f = 0; f < nb_fishery; f++) {
        if (param.mask_fishery_sp[sp][f]) {
            if (param.mask_fishery_sp_no_effort[sp][f]) {
                mat.dvarCatch_est[sp][k].initialize();
                mat.catch_est[sp][k].initialize();

                for (int age = a0; age < nb_ages; age++) {
                    predicted_catch_fishery_no_effort_comp(
                        map, param, mat, f, k, sp, age);

                    save_identifier_string2(
                        (char*)"predicted_catch_no_effort_begin");
                    mat.dvarCtot_age_est(sp, age).save_dvar_matrix_position();
                    mat.dvarCatch_est(sp, k).save_dvar_matrix_position();
                    save_int_value(year);
                    save_int_value(month);
                    save_int_value(sp);
                    save_int_value(age);
                    save_int_value(f);
                    save_int_value(k);
                    unsigned long int crw = (unsigned long int)&rw;
                    save_long_int_value(crw);
                    unsigned long int pop = (unsigned long int)this;
                    save_long_int_value(pop);
                    unsigned long int pmap = (unsigned long int)&map;
                    save_long_int_value(pmap);
                    unsigned long int cparam = (unsigned long int)&param;
                    save_long_int_value(cparam);
                    save_identifier_string2(
                        (char*)"predicted_catch_no_effort_end");

                    gradient_structure::GRAD_STACK1->set_gradient_stack(
                        dv_predicted_catch_without_effort_fishery);
                }
            }
            k++;
        }
    }
}

void dv_predicted_catch_without_effort_fishery() {
    verify_identifier_string2((char*)"predicted_catch_no_effort_end");
    unsigned long int pos_param = restore_long_int_value();
    unsigned long int pos_map = restore_long_int_value();
    unsigned long int pos_pop = restore_long_int_value();
    unsigned long int pos_rw = restore_long_int_value();
    unsigned k = restore_int_value();
    unsigned f = restore_int_value();
    unsigned age = restore_int_value();
    unsigned sp = restore_int_value();
    unsigned month = restore_int_value();
    unsigned year = restore_int_value();
    const dvar_matrix_position catch_pos = restore_dvar_matrix_position();
    const dvar_matrix_position c_age_pos = restore_dvar_matrix_position();
    verify_identifier_string2((char*)"predicted_catch_no_effort_begin");

    dmatrix dfC_age_est = restore_dvar_matrix_derivatives(c_age_pos);
    dmatrix dfC_est = restore_dvar_matrix_derivatives(catch_pos);

    CReadWrite* rw = (CReadWrite*)pos_rw;
    CCalpop* pop = (CCalpop*)pos_pop;
    PMap* map = (PMap*)pos_map;
    CParam* param = (CParam*)pos_param;

    const int imax = map->imax;
    const int imin = map->imin;
    const int jmax = map->jmax;
    const int jmin = map->jmin;

    dvector lat_correction(jmin, jmax);
    lat_correction.initialize();
    for (int j = jmin; j <= jmax; j++) {
        double lastlat = param->lastlat(j);
        lat_correction[j] = param->correction_lat(lastlat);
    }

    dmatrix Ctot_proportion_fishery;
    Ctot_proportion_fishery.allocate(imin, imax, map->jinf, map->jsup);
    pop->Recomp_C_fishery_proportion_in_Ctot(
        *map, *param, *rw, Ctot_proportion_fishery, year, month, sp, k);

    const int C_units = param->fishery_catch_units[f];
    const double weight = pow(1e-3 * param->weight[sp][age], C_units);
    const double area = 1.852 * param->deltaX * 1.852 * param->deltaY;  // sq.km

    for (int i = imax; i >= imin; i--) {
        const int jmin = map->jinf[i];
        const int jmax = map->jsup[i];
        for (int j = jmax; j >= jmin; j--) {
            const double C_proportion_f = Ctot_proportion_fishery(i, j);
            if (map->carte(i, j) && C_proportion_f) {
                // mat.dvarCatch_est(sp,k).elem_value(i,j) += C_f_age_est_cell *
                // weight;
                double dfC_f_age_est_cell = weight * dfC_est(i, j);
                // dfC_pr_est(i,j)    += dfC_est(i,j);
                // dfC_est(i,j)        = 0.0;

                // double C_f_age_est_cell =
                // value(mat.dvarCtot_age_est(sp,age,i,j)) *
                //	* mat.Ctot_proportion_fishery(k,i,j) * area /
                //mat.lat_correction(j);
                dfC_age_est(i, j) += dfC_f_age_est_cell * C_proportion_f *
                                     area / lat_correction(j);
                dfC_f_age_est_cell = 0.0;
            }
        }
    }
    dfC_age_est.save_dmatrix_derivatives(c_age_pos);
    dfC_est.save_dmatrix_derivatives(catch_pos);
}
