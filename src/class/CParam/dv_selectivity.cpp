#include "Param.h"

/// Adjoint code for selectivity functions (Param class)
/// Forward functions are in Param.cpp

void CParam::dfselectivity(
    double& dfslope, double& dflength, double& dfasympt, const int sp,
    const int age, const int f, const int k) {
    const double s_slope = s_slope_sp_fishery(sp, k);
    const double l = length(sp, age);

    dfslope = 0;
    dflength = 0;
    dfasympt = 0;

    double s_length, s_asympt, dff, sigma_sq, g;

    switch (s_func_type[f]) {
        case 1:  // logistic
            // selectivity = length(sp,age)/(s_slope_sp_fishery(sp,k) +
            // length(sp,age));
            dfslope = -l / pow(s_slope + l, 2);
            break;

        case 2:  // sigmoid
            s_length = s_length_sp_fishery(sp, k);
            // selectivity =
            // 1/(1+exp(s_slope_sp_fishery(sp,k)*(s_length_sp_fishery(sp,k)-length(sp,age))));
            dff = exp(s_slope * (s_length - l)) /
                  pow(1 + exp(s_slope * (s_length - l)), 2);
            dflength = -s_slope * dff;
            dfslope = (l - s_length) * dff;
            break;

        case 3:  // asymmetric Gaussian
            s_length = s_length_sp_fishery(sp, k);
            sigma_sq = pow(s_slope, 2);
            g = exp(-pow(l - s_length, 2) / sigma_sq);
            if (length(sp, age) <= s_length_sp_fishery(sp, k)) {
                // selectivity =
                // exp(-pow(length(sp,age)-s_length_sp_fishery(sp,k),2)/sigma_sq);
                dfslope = (2 * pow(l - s_length, 2) / pow(s_slope, 3)) * g;
                dflength = 2 * ((l - s_length) / sigma_sq) * g;
            } else {
                s_asympt = s_asympt_sp_fishery(sp, k);
                // selectivity =
                // (1-s_asympt_sp_fishery(sp,k))*exp(-pow(length(sp,age)-s_length_sp_fishery(sp,k),2)/sigma_sq)+s_asympt_sp_fishery(sp,k);
                dfslope = (1 - s_asympt) *
                          (2 * pow(l - s_length, 2) / pow(s_slope, 3)) * g;
                dflength = 2 * ((l - s_length) / sigma_sq) * g * (1 - s_asympt);
                dfasympt = -g + 1;
            }
            break;
        case 4:  // two-modal Gaussian: not done yet!
            const double sl1 = s_length_sp_fishery(sp, k);
            const double sl2 = s_asympt_sp_fishery(sp, k);  // temporal
            sigma_sq = pow(s_slope, 2);
            const double gm1 = exp(-pow(sl2 - sl1, 2) / sigma_sq);
            const double gm2 = exp(-pow(sl1 - sl2, 2) / sigma_sq);
            double dfgm1 = -1.0 / pow(1 + gm1 + gm2, 2);
            double dfgm2 = -1.0 / pow(1 + gm1 + gm2, 2);

            // const double gm2 = exp(-pow(sl1-sl2,2)/sigma_sq);
            dfslope += (2 * pow(sl1 - sl2, 2) / pow(s_slope, 3)) * gm2 * dfgm2;
            // const double gm1 = exp(-pow(sl2-sl1,2)/sigma_sq);
            dfslope += (2 * pow(sl2 - sl1, 2) / pow(s_slope, 3)) * gm1 * dfgm1;

            // const double g2 = exp(-pow(l-sl2,2)/sigma_sq);
            // not finished!!!
            break;
    }
}
