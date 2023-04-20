#include "calpop.h"

/// Forward main function called in simulation mode only for:
/// precaldia and caldia functions. See caldia.cpp

void CCalpop::Precaldia_Caldia(
    const PMap& map, VarParamCoupled& param, VarMatrices& mat,
    dvar_matrix& habitat, dvar_matrix& total_pop, const int sp, const int age,
    const int t_count, const int jday) {
    dvariable mss_species = param.dvarsMSS_species[sp];
    dvariable mss_size_slope = param.dvarsMSS_size_slope[sp];
    dvariable c_diff_fish = param.dvarsC_diff_fish[sp];
    dvariable sigma_species = param.dvarsSigma_species[sp];

    dvar_matrix Mss_species(map.imin, map.imax, map.jinf, map.jsup);
    Mss_species = mss_species;
    dvar_matrix Mss_size_slope(map.imin, map.imax, map.jinf, map.jsup);
    Mss_size_slope = mss_size_slope;
    dvar_matrix C_diff_fish(map.imin, map.imax, map.jinf, map.jsup);
    C_diff_fish = c_diff_fish;
    dvar_matrix Sigma_species(map.imin, map.imax, map.jinf, map.jsup);
    Sigma_species = sigma_species;

    dvar_matrix W(0, maxn - 1, 0, maxn - 1);
    W.initialize();

    precaldia_comp(
        map, param, mat, value(habitat), value(total_pop), value(mss_species),
        value(mss_size_slope), value(sigma_species), value(c_diff_fish), sp,
        age, jday);

    mat.dvarsAdvection_x = nograd_assign(mat.advection_x);
    mat.dvarsAdvection_y = nograd_assign(mat.advection_y);
    mat.dvarsDiffusion_x = nograd_assign(mat.diffusion_x);
    mat.dvarsDiffusion_y = nograd_assign(mat.diffusion_y);

    caldia(
        map, param, value(mat.dvarsDiffusion_x), value(mat.dvarsAdvection_x),
        value(mat.dvarsDiffusion_y), value(mat.dvarsAdvection_y));

    dvarsA = nograd_assign(a);
    dvarsB = nograd_assign(b);
    dvarsC = nograd_assign(c);
    dvarsD = nograd_assign(d);
    dvarsE = nograd_assign(e);
    dvarsF = nograd_assign(f);
    Ybet = nograd_assign(ybet);
    /*if (!param.gcalc()){
            // only in simulation mode: compute mean speed in BL/sec and mean
    diffusion rate in nmi^2/day
            mat.MeanVarMovement(map,value(mat.dvarsAdvection_x),value(mat.dvarsAdvection_y),
                                value(mat.dvarsDiffusion_y),value(mss_species),value(sigma_species),
                                param.length(sp,age),param.length(sp,param.sp_nb_cohorts[sp]-1),param.deltaT,sp,age);
    }*/
}
