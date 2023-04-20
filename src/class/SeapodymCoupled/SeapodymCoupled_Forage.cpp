#include "SeapodymCoupled.h"

void SeapodymCoupled::PredationMortality(int t, dmatrix total_pop) {
    d4_array consumption;
    consumption.allocate(0, nb_species - 1);
    for (int sp = 0; sp < nb_species; sp++) {
        consumption[sp].allocate(a0_adult[sp], aN_adult[sp] - 1);
        for (int a = a0_adult[sp]; a < aN_adult[sp]; a++) {
            consumption[sp][a].allocate(map.imin, map.imax, map.jinf, map.jsup);
            consumption[sp][a].initialize();
        }
    }

    for (int n = 0; n < nb_forage; n++) {
        sumF_area_pred[n] = 0;
        sumF_required_by_sp[n] = 0;  // total de la biomasse en forage requise
                                     // (predation par especes decrites)
        mean_omega_sp[n] = 0;        // mortalite due aux especes decrites

        for (int sp = 0; sp < nb_species; sp++) {
            for (int i = 1; i < nbi - 1; i++) {
                for (int j = 1; j < nbj - 1; j++) {
                    if (map.carte[i][j]) {
                        //----------------------------------------//
                        // TOTAL FORAGE BIOMASS IN PREDATOR AREA  //
                        //----------------------------------------//
                        if (total_pop[i][j] > 0) {
                            // units of F: g/m^2 = tones/km^2
                            sumF_area_pred[n] += mat.forage[t][n][i][j];
                            for (int age = a0_adult[sp]; age < aN_adult[sp];
                                 age++) {
                                // 2011/10: deleted nF_ratio (not used in the
                                // moment, to reduce number of variables
                                // allocated), to restore later if needed
                                double a = param->forage_ration
                                               [sp];  //* mat.nF_ratio[n];
                                consumption(sp, age, i, j) = a;
                                // units of tuna = thous. Nb per sq km translate
                                // to tones/km^2
                                sumF_required_by_sp[n] +=
                                    a * value(mat.dvarDensity[sp][age][i][j]) *
                                    param->weight[sp][age] * 0.001;
                            }
                        }
                    }
                }
            }
            mean_omega_sp[n] = sumF_required_by_sp[n] / sumF_area_pred[n];
        }
    }

    for (int n = 0; n < nb_forage; n++) {
        for (int i = 1; i < nbi - 1; i++) {
            for (int j = 1; j < nbj - 1; j++) {
                if (map.carte[i][j]) {
                    double F = mat.forage[t][n][i][j];
                    double FR = 0;
                    // double a_bar = 0;
                    for (int sp = 0; sp < nb_species; sp++) {
                        for (int age = a0_adult[sp]; age < aN_adult[sp];
                             age++) {
                            // units of tuna here are g/m^2 as it goes to forage
                            // eqns: Nb/km^2 * W(mt) =
                            // (Nb)*1,000,000(g)/1,000,000*m^2 = g/m^2
                            double a = consumption(sp, age, i, j);
                            FR += a * value(mat.dvarDensity[sp][age][i][j]) *
                                  param->weight[sp][age] * 0.001;
                        }
                        // a_bar /= param->sp_nb_age_class_ad[sp];
                    }
                    double mort_by_sp = 0;
                    // double h = 0.025;
                    // mort_by_sp = a_bar* FR/(1+a_bar*h*F); //if (m>1) m=1;
                    // mort_by_sp = FR; //Lotka-Volterra trophic function
                    if (F > 0)
                        mort_by_sp = FR / F;  // Patrick's trophic function

                    mat.mortality[n][i][j] =
                        func.function_lambda(*param, mat, n, i, j) +
                        mort_by_sp - mean_omega_sp[n];
                }
            }
        }
    }
}

void SeapodymCoupled::SolveADRE(d3_array F, int n) {
    pop.precaldia(*param, map, mat);
    pop.caldia(
        map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y,
        mat.advection_y);

    pop.precalrec(map, mat.mortality[n]);
    pop.calrec(map, F[n], mat.mortality[n]);
}

void SeapodymCoupled::SolveADE(d4_array S, int n, int ntimes) {
    pop.precaldia(*param, map, mat);

    pop.caldia(
        map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y,
        mat.advection_y);

    pop.precalrec(map, mat.mortality[n]);

    for (int nbt = 0; nbt < ntimes; nbt++) {
        // call "calrec"  for transport of source matrices from 0 to Tr_max
        pop.calrec(map, S[n][nbt], mat.mortality[n]);
    }
}
