#include "calpop.h"

void CCalpop::tridag_GO(
    const DVECTOR& at, const DVECTOR& bt, const DVECTOR& ct, const DVECTOR& rhs,
    DVECTOR& uvec, DVECTOR& gam, const int inf, const int sup) {
    double bet = bt[inf];

    uvec[inf] = rhs[inf] / bet;
    for (int j = inf + 1; j <= sup; j++) {
        // double g=ct[j-1]/bet;
        gam[j] = ct[j - 1] / bet;
        // bet=bt[j]-at[j]*g;
        bet = bt[j] - at[j] * gam[j];
        uvec[j] = (rhs[j] - at[j] * uvec[j - 1]) / bet;
        // gam[j]=g;
    }

    for (int j = sup - 1; j >= inf; j--) uvec[j] -= gam[j + 1] * uvec[j + 1];

    return;
}

void CCalpop::tridag(
    const DVECTOR& at, const DVECTOR& bet, const DVECTOR& ct,
    const DVECTOR& rhs, DVECTOR& uvec, DVECTOR& gam, int inf, int sup) {
    gam[inf] = rhs[inf];
    for (int j = inf + 1; j <= sup; j++)
        gam[j] = rhs[j] - gam[j - 1] * at[j] * bet[j - 1];

    uvec[sup] = gam[sup] * bet[sup];
    for (int j = sup - 1; j >= inf; j--)
        uvec[j] = (gam[j] - ct[j] * uvec[j + 1]) * bet[j];

    return;
}

void CCalpop::tridag1_0(
    const dvector& at, dvar_vector& bet, const dvector& ct, dvar_vector& rhs,
    dvar_vector& uvec, dvar_vector& gam, int inf, int sup) {
    gam[inf] = rhs[inf];
    for (int j = inf + 1; j <= sup; j++)
        gam[j] = rhs[j] - gam[j - 1] * at[j] * bet[j - 1];

    uvec[sup] = gam[sup] * bet[sup];
    for (int j = sup - 1; j >= inf; j--)
        uvec[j] = (gam[j] - ct[j] * uvec[j + 1]) * bet[j];

    return;
}

void CCalpop::tridag1_1(
    const dvector& at, const dvector& bet, const dvector& ct, dvar_vector& rhs,
    dvar_vector& uvec, dvar_vector& gam, int inf, int sup) {
    gam[inf] = rhs[inf];
    for (int j = inf + 1; j <= sup; j++)
        gam[j] = rhs[j] - gam[j - 1] * at[j] * bet[j - 1];

    uvec[sup] = gam[sup] * bet[sup];
    for (int j = sup - 1; j >= inf; j--)
        uvec[j] = (gam[j] - ct[j] * uvec[j + 1]) * bet[j];

    return;
}

void CCalpop::tridag2(
    dvar_vector& at, dvar_vector& bet, dvar_vector& ct, dvar_vector& rhs,
    dvar_vector& uvec, dvar_vector& gam, int inf, int sup) {
    gam[inf] = rhs[inf];
    for (int j = inf + 1; j <= sup; j++)
        gam[j] = rhs[j] - gam[j - 1] * at[j] * bet[j - 1];

    uvec[sup] = gam[sup] * bet[sup];
    for (int j = sup - 1; j >= inf; j--)
        uvec[j] = (gam[j] - ct[j] * uvec[j + 1]) * bet[j];

    return;
}
