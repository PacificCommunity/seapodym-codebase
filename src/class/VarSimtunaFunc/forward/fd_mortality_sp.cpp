#include "VarSimtunaFunc.h"

/// Forward main function called in simulation mode only for:
/// mortality rates at age. See mortality_sp.cpp

void VarSimtunaFunc::Mortality_Sp(
    VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& M,
    dvar_matrix& H, const int sp, const double mean_age_in_dtau, const int age,
    const int t_count) {
    M.initialize();

    dvariable Mp_max = param.dvarsMp_mean_max[sp];
    dvariable Mp_exp = param.dvarsMp_mean_exp[sp];
    dvariable Ms_slope = param.dvarsMs_mean_slope[sp];
    dvariable Ms_max = param.dvarsMs_mean_max[sp];
    dvariable range = param.dvarsM_mean_range[sp];

    // Note, the index '0' here assumes that all age classes except A+ have the
    // same size.
    int dtau = param.sp_unit_cohort[sp][0] * param.deltaT;

    M_sp_comp(
        map, M, value(H), value(Mp_max), value(Ms_max), value(Mp_exp),
        value(Ms_slope), value(range), mean_age_in_dtau, dtau);

    if (!param.gcalc()) {
        // pH option works only in simulation mode
        if ((age == 0) && param.use_ph1)
            M_PH_juv_comp(
                param, map, mat, M, mat.ph1[t_count], mean_age_in_dtau);

        mat.MeanVarMortality(
            map, value(M), value(Mp_max), value(Ms_max), value(Mp_exp),
            value(Ms_slope), mean_age_in_dtau, sp, age);
    }
}
