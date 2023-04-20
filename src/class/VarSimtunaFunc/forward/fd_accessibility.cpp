#include "VarSimtunaFunc.h"

/// Forward main functions called in simulation mode only for:
/// 1) accessibility to forage components (f_accessibility) or to
/// their respective layers (f_accessibility_layer).
/// 2) average currents given the accessibility to the layer.
/// See accessibility.cpp

void VarSimtunaFunc::Faccessibility(
    VarParamCoupled& param, VarMatrices& mat, const PMap& map, const int sp,
    const int jday, const int t_count, const int pop_built, const int tags_only,
    const ivector tags_age_solve) {
    const int a0 = param.sp_a0_adult[sp];
    const int nb_ages = param.sp_nb_cohorts[sp];

    dvariable temp_max = param.dvarsB_sst_spawning[sp];
    dvariable temp_min = param.dvarsB_sst_habitat[sp];
    dvariable oxy_teta = param.dvarsA_oxy_habitat[sp];
    dvariable oxy_cr = param.dvarsB_oxy_habitat[sp];
    dvariable temp_age_slope = param.dvarsT_age_size_slope[sp];

    if (param.gaussian_thermal_function[sp]) {
        dvariable sigma_ha = param.dvarsA_sst_habitat[sp];
        dvariable sigma_hs = param.dvarsA_sst_spawning[sp];
    }
    if (!param.gaussian_thermal_function[sp]) {
        dvariable delta3 = param.dvarsThermal_func_delta[2][sp];
    }

    for (int age = a0; age < nb_ages; age++) {
        if (param.age_compute_habitat[sp][age] !=
            param.age_compute_habitat[sp][age - 1]) {
            if (!tags_only || tags_age_solve(age)) {
                Faccessibility_comp(
                    param, mat, map, value(temp_max), value(oxy_teta),
                    value(oxy_cr), sp, age, jday, t_count);
            }
        }
    }
}

void VarSimtunaFunc::Average_currents(
    VarParamCoupled& param, VarMatrices& mat, const PMap& map, int age,
    const int t_count, const int pop_built) {
    mat.dvarsU.initialize();
    mat.dvarsV.initialize();

    Average_currents_comp(param, mat, map, age, t_count);
}
