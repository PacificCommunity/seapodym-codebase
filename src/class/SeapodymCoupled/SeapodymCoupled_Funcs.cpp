#include "SeapodymCoupled.h"

dmatrix SeapodymCoupled::DtoBcell(const dmatrix var) {
    // dmatrix Bcell(0,nbi-1,0,nbj-1); Bcell.initialize();
    dmatrix Bcell(map.imin, map.imax, map.jinf, map.jsup);
    Bcell.initialize();
    Bcell = var;

    const int imin = map.imin;
    const int imax = map.imax;
    for (int i = imin; i <= imax; i++) {
        const int jmin = map.jinf[i];
        const int jmax = map.jsup[i];
        for (int j = jmin; j <= jmax; j++) {
            if (map.carte[i][j]) {
                Bcell(i, j) *= cell_area / mat.lat_correction[j];
            }
        }
    }
    return (Bcell);
}

void SeapodymCoupled::getDate(int& jday) {
    int newyear;  //(Inna 12/10/11) no need for the moment
    Date::update_time_variables(
        t_count, param->deltaT, param->date_mode, jday_spinup, jday, day, month,
        year, newyear);
    qtr = (int)(((month - 1) / 3) + 1);
    date_str = Utilities::MakeDate(year, month, day);
}

void SeapodymCoupled::SaveIntermediate(const int sp, const int age) {
    // save distribution to the file
    dmatrix mat2d(0, nbi - 1, 0, nbj - 1);
    mat2d.initialize();
    for (int i = map.imin; i <= map.imax; i++) {
        for (int j = map.jinf[i]; j <= map.jsup[i]; j++) {
            if (map.carte[i][j]) {
                mat2d(i, j) = value(mat.dvarDensity(sp, age, i, j));
            }
        }
    }

    rw.wbin_mat2d("interm.tmp", mat2d, nbi, nbj, true);
}

void SeapodymCoupled::SaveDistributions(const int year, const int month) {
    for (int sp = 0; sp < nb_species; sp++) {
        // write down the state vector to be used as initial condition
        string dirname = param->strdir_output;
        string date;
        std::stringstream ss;
        ss << year << "_" << month;
        ss >> date;

        // Format of IC file: single regular DYM2 file
        string fileCohorts =
            dirname + param->sp_name[sp] + "_cohorts" + date + ".dym";
        cout << "Saving distributions to " << fileCohorts << endl;

        /// int nb_cohorts = param->sp_nb_age_class_ad[sp]+3;
        int nb_cohorts = param->sp_nb_cohorts[sp];
        double minval = .0;
        double maxval = max(value(mat.dvarDensity(sp, 0)));
        dvector zlevel(0, nb_cohorts - 1);
        for (int n = 0; n < nb_cohorts; n++) zlevel[n] = n + 1;
        rw.wbin_header(
            fileCohorts, param->idformat, param->idfunc, minval, maxval,
            param->nlong, param->nlat, nb_cohorts, zlevel[0],
            zlevel[nb_cohorts - 1], mat.xlon, mat.ylat, zlevel, mat.mask);

        SaveCohorts(fileCohorts, sp, true);
        /// SaveJuvCohorts(fileCohorts, sp, true);
        /// SaveAdultCohorts(fileCohorts, sp, true);
    }
}

void SeapodymCoupled::SaveCohorts(string fileout, int sp, bool FileMode) {
    double mult = 1.0;
    const int nb_ages = param->sp_nb_cohorts[sp];
    for (int a = 0; a < nb_ages; a++) {
        dmatrix mat2d(0, nbi - 1, 0, nbj - 1);
        mat2d.initialize();

        for (int i = map.imin; i <= map.imax; i++) {
            for (int j = map.jinf[i]; j <= map.jsup[i]; j++) {
                if (map.carte[i][j]) {
                    mat2d(i - 1, j - 1) =
                        mult * value(mat.dvarDensity(sp, a, i, j));
                }
            }
        }
        rw.wbin_transpomat2d(fileout, mat2d, nbi - 2, nbj - 2, FileMode);

        FileMode = true;
    }
}

void SeapodymCoupled::InitializeAll() {
    t_count = nbt_building + 1;
    mat.mats.initialize();
    mat.density_before.initialize();
    mat.density_after.initialize();
    mat.dvarDensity.initialize();
    mat.dvarCatch_est.initialize();

    for (int sp = 0; sp < nb_species; sp++) {
        double sfactor = 1.0;
        const int nb_ages = param->sp_nb_cohorts[sp];
        for (int a = 0; a < nb_ages; a++) {
            mat.dvarDensity(sp, a) = sfactor * mat.init_density_species(sp, a);
        }
    }
}

void SeapodymCoupled::AverageCurrents(int t, int n) {
    // return the (u,v) current which transports the forage component
    // non migrant : current in the layer
    // migrant : average current relative to time spent in day and night layers
    if (param->day_layer[n] == param->night_layer[n]) {  // non migrant forage
        int dn_lay = param->day_layer[n];
        for (int i = 1; i < nbi - 1; i++) {
            for (int j = 1; j < nbj - 1; j++) {
                if (map.carte[i][j]) {
                    mat.u[i][j] = mat.un[t][dn_lay][i][j];
                    mat.v[i][j] = mat.vn[t][dn_lay][i][j];
                }
            }
        }
    } else {  // migrant forage
        int d_lay = param->day_layer[n];
        int n_lay = param->night_layer[n];
        for (int i = 1; i < nbi - 1; i++) {
            for (int j = 1; j < nbj - 1; j++) {
                if (map.carte[i][j]) {
                    double DL = mat.daylength[i][j];
                    mat.u[i][j] = ((DL * mat.un[t][d_lay][i][j]) +
                                   ((24 - DL) * mat.un[t][n_lay][i][j])) /
                                  24;
                    mat.v[i][j] = ((DL * mat.vn[t][d_lay][i][j]) +
                                   ((24 - DL) * mat.vn[t][n_lay][i][j])) /
                                  24;
                }
            }
        }
    }
}

void SeapodymCoupled::CalcMeanTemp(
    const int t_count,
    const int
        tcur) {  // only in simulation mode: compute mean temperature at age
    if (!param->gcalc()) {
        if (t_count == 1) {
            for (int sp = 0; sp < nb_species; sp++) {
                for (int age = 0; age < aN_adult(sp); age++)
                    if (age < a0_adult(sp))
                        mat.mean_temperature(sp, age) =
                            param->b_sst_spawning[sp];
                    else
                        mat.mean_temperature(sp, age) =
                            param->temp_age[sp][age];
            }
        } else {
            for (int sp = 0; sp < nb_species; sp++) {
                // For larvae and juveniles:
                mat.MeanVarTemperature(
                    map, sp, param->sp_nb_cohort_lv[sp], a0_adult[sp], tcur);

                // Adult cohorts
                for (int age = a0_adult(sp); age < aN_adult(sp); age++)
                    mat.mean_temperature(sp, age) =
                        func.Tmean_comp(*param, mat, map, sp, age, tcur);
            }
        }
    }
}

void SeapodymCoupled::CalcSums() {
    mat.larvae.initialize();
    mat.juvenile.initialize();
    mat.young.initialize();
    mat.recruit.initialize();
    mat.adult.initialize();
    mat.total_pop.initialize();
    mat.total_pred_catch.initialize();
    if (!param->tags_only) {
        mat.total_obs_catch.initialize();
    }
    mat.sum_B_larvae.initialize();
    mat.sum_B_juv.initialize();
    mat.sum_B_young.initialize();
    mat.sum_B_recruit.initialize();
    mat.sum_B_adult.initialize();
    mat.sum_total_pop.initialize();

    // total catch per time step
    double sum_total_catch = 0;

    //------------------------------------------------------------------
    // Note, All population stages will be stored in the units of density:
    // Nb/km2 for stages Juv - Rercruits, tonnes/km2 for young - adults
    //------------------------------------------------------------------
    for (int sp = 0; sp < nb_species; sp++) {
        dvector W_mt;
        W_mt.allocate(0, aN_adult[sp] - 1);
        for (int age = 0; age < aN_adult(sp); age++) {
            W_mt[age] = param->weight[sp][age] * 0.001;
        }

        const int age_recruits = param->age_recruit[sp];
        const int nb_lv = param->sp_nb_cohort_lv[sp];
        for (int i = map.imin; i <= map.imax; i++) {
            for (int j = map.jinf[i]; j <= map.jsup[i]; j++) {
                if (map.carte[i][j]) {
                    double lat_corrected_area =
                        cell_area / mat.lat_correction[j];
                    // test:
                    // double lat_corrected_area = param->cell_surface_area(j);
                    //-----------------------------------------------------------
                    //  Larvae = age 0
                    //  units = Nb/km^2
                    //  UNITS = (Nb/km^2) * W(mt) * cell_area(km^2) = TONNES
                    //  over all domain
                    //------------------------------------------------------------
                    for (int age = 0; age < nb_lv; age++) {
                        mat.larvae[sp][i][j] =
                            value(mat.dvarDensity[sp][age][i][j]);
                        mat.sum_B_larvae[sp] +=
                            value(mat.dvarDensity[sp][age][i][j]) * W_mt(age) *
                            lat_corrected_area;
                        // mat.sum_B_larvae[sp] +=
                        // value(mat.dvarDensity[sp][age][i][j])
                        // *lat_corrected_area;
                    }
                    //-----------------------------------------------------------
                    // Juveniles  = age  1 to age_autonomous (1st quarter)
                    // biomass = sum n*weight by age class
                    // UNITS = (Nb/km^2) * W(mt) = tonnes/km^2, UNITS of sum =
                    // TONNES over entire domain
                    //------------------------------------------------------------
                    for (int age = nb_lv; age < a0_adult[sp]; age++) {
                        // for (int age=nb_lv; age< age_recruits; age++){
                        mat.juvenile[sp][i][j] +=
                            value(mat.dvarDensity[sp][age][i]
                                                 [j]);  // in numbers as larvae
                        mat.sum_B_juv[sp] +=
                            value(mat.dvarDensity[sp][age][i][j]) * W_mt(age) *
                            lat_corrected_area;
                        // mat.sum_B_juv[sp] +=
                        // value(mat.dvarDensity[sp][age][i][j])*lat_corrected_area;
                    }
                    // for (int age= 1; age < param->sp_nb_age_class_jv[sp];
                    // age++){ mat.juvenile[sp][i][j] +=
                    // value(mat.dvarJuv_species[sp][age][i][j]) *
                    // param->juv_weight[sp][age];
                    // mat.sum_B_juv[sp] +=
                    // mat.juvenile[sp][i][j]*area/mat.lat_correction[j];
                    //-----------------------------------------------------------
                    //  recruit = age  at recruitment
                    //  biomass = sum n*weight by age class
                    //  numbers per sq.km per month
                    //  UNITS = (Nb/km^2) * W(mt) = tonnes/km^2, UNITS of sum =
                    //  TONNES over entire domain
                    //------------------------------------------------------------
                    if (!param->tag_like[sp]) {
                        for (int age = age_recruits - 1;
                             age <= age_recruits + 1;
                             age++) {  // MFCL qtr class
                            mat.recruit[sp][i][j] +=
                                value(mat.dvarDensity[sp][age][i]
                                                     [j]);  // in numbers/sq.km
                            // for SumDym write weight
                            mat.sum_B_recruit[sp] +=
                                value(mat.dvarDensity[sp][age][i][j]) *
                                W_mt(age) * lat_corrected_area;  // total number
                        }
                    }
                    if (param->tag_like[sp]) {
                        for (int aa = a0_adult[sp]; aa < aN_adult[sp]; aa++) {
                            double total_tags = 0.0;
                            for (int pop = 1; pop <= param->nb_tag_files; pop++)
                                total_tags +=
                                    value(mat.dvarDensity(pop, aa, i, j));

                            mat.recruit[0][i][j] += total_tags;
                            if (value(mat.dvarDensity(1, aa, i, j)) < 0)
                                cerr << "NEGATIVE BIOMASS for " << aa << " "
                                     << value(mat.dvarDensity(1, aa, i, j))
                                     << endl;
                        }
                    }
                    // mat.sum_B_recruit[sp] +=
                    // value(mat.dvarDensity[sp][age_recruits][i][j])*W_mt(age_recruits)*lat_corrected_area;
                    //-----------------------------------------------------------
                    //  (Young) = age  autonomous to age mature
                    //  ADULT TUNA = tuna biomass from age of maturity
                    //  biomass = sum n*weight by age class
                    //  UNITS = (Nb/km^2) * W(mt) = tonnes/km^2, UNITS of sum =
                    //  TONNES over entire domain
                    //------------------------------------------------------------
                    for (int age = a0_adult[sp]; age < aN_adult[sp]; age++) {
                        if (param->maturity_age[sp][age] < 0.5)
                            mat.young[sp][i][j] +=
                                value(mat.dvarDensity[sp][age][i][j]) *
                                W_mt[age];  // * swa(sp,age);
                        else
                            mat.adult[sp][i][j] +=
                                value(mat.dvarDensity[sp][age][i][j]) *
                                W_mt[age];  // * swa(sp,age);

                        //	mat.young[sp][i][j] +=
                        //(1-param->maturity_age[sp][age])*value(mat.dvarDensity[sp][age][i][j])
                        //* W_mt[age];// * swa(sp,age); 	mat.adult[sp][i][j] +=
                        //param->maturity_age[sp][age]*value(mat.dvarDensity[sp][age][i][j])
                        //* W_mt[age];// * swa(sp,age);
                    }
                    mat.sum_B_young[sp] +=
                        mat.young[sp][i][j] * lat_corrected_area;
                    mat.sum_B_adult[sp] +=
                        mat.adult[sp][i][j] * lat_corrected_area;
                    //----------------------------------------------------------
                    // TOTAL TUNA BIOMASS feeding on the forage
                    // sum n * weight by age class
                    // UNITS = tonnes/km^2, UNITS of sum = TONNES over entire
                    // domain
                    //----------------------------------------------------------
                    mat.total_pop[sp][i][j] =
                        mat.young[sp][i][j] + mat.adult[sp][i][j];
                    mat.sum_total_pop[sp] +=
                        mat.total_pop[sp][i][j] * lat_corrected_area;

                    //----------------------------------------------------------
                    // TOTAL OBSERVED AND PREDICTED CATCH FOR SPECIES
                    // UNITS = tonnes/cell
                    //----------------------------------------------------------
                    if (!param->tags_only) {  // in this simulations will write
                                              // obs and pred tag numbers,
                        // so this loop should be omited
                        int k = 0;
                        int nbe = 0;
                        for (int f = 0; f < nb_fishery; f++) {
                            if (param->mask_fishery_sp[0][f]) {
                                if (mat.effort[f][i][j]) {
                                    mat.total_obs_catch[sp][i][j] +=
                                        mat.catch_obs[sp][k][i][j];

                                    mat.total_pred_catch[sp][i][j] +=
                                        mat.catch_est[sp][k][i][j];

                                    nbe++;
                                }
                                k++;
                            }
                        }
                    }

                    // to distinguish between 0 catch and no fishing events.
                    //					if (nbe==0) {
                    //						mat.total_obs_catch[sp][i][j]  =
                    //-1; 						mat.total_pred_catch[sp][i][j] = -1;
                    //					}
                    sum_total_catch += mat.total_obs_catch[sp][i][j];
                }
            }
        }
    }
    SUM_CATCH += sum_total_catch;
}

void SeapodymCoupled::ConsoleOutput(int flag_simulation, double like) {
    if (flag_simulation) {
        if (!param->tags_only) {
            cout << setw(4) << left << t_count << "| " << setw(14) << date_str
                 << "| " << setw(14) << mat.sum_B_larvae[0] << " | " << setw(14)
                 << mat.sum_B_juv[0] << " | " << setw(14) << mat.sum_B_young[0]
                 << " | " << setw(14) << mat.sum_B_adult[0] << " | " << setw(14)
                 << mat.sum_total_pop[0] << " | " << like << endl;
        } else {
            cout << setw(4) << left << t_count << "| " << setw(14) << date_str
                 << "| " << setw(14) << sum(mat.total_obs_catch(0)) << " | "
                 << setw(14) << sum(mat.total_pred_catch(0)) << " | " << like
                 << endl;
        }
    } else {
        if (!param->tags_only) {
            // CalcSums();
            //  Comment: Larvae biomass is the result of spawning function
            //  S(t-1) and ADRE at time t we have Larvae(t=1)=0 because S(t-1)
            //  is computed after CalcSums()
            cout << endl << "Model temporal dynamics:" << endl;
            cout << setw(4) << left << "step"
                 << "| " << setw(14) << "Date"
                 << "| " << setw(14) << "Larvae"
                 << " | " << setw(14) << "Juveniles"
                 << " | " << setw(14) << "Young"
                 << " | " << setw(14) << "Adult"
                 << " | " << setw(14) << "Total_pop"
                 << " | " << setw(14) << "Likelihood" << endl;

            cout << "----------------------------------------------------------"
                    "----------------------------------------------------------"
                 << endl;
        } else {
            cout << endl << "Model temporal dynamics:" << endl;
            cout << setw(4) << left << "step"
                 << "| " << setw(14) << "Date"
                 << "| " << setw(14) << "Observed tags"
                 << " | " << setw(14) << "Modelled tags"
                 << " | " << setw(14) << "Likelihood" << endl;

            cout << "----------------------------------------------------------"
                    "--------"
                 << endl;
        }
    }
}

void SeapodymCoupled::Set_density_region_age_zero(
    const int sp) {  // Function to be written. For the moment the old temporal
                     // code commented here

    // connectivity simulations:
    // Bismarck Sea, mask file po_soda_v2_opsat_mask_eez_PNG-AW_HSP.txt, ID=-92
    // Bismarck Sea, mask file mask_ecco_PNG_AW.txt ID=-2
    // Mariana Islands, mask file mask_Guam_NMI.txt, ID=-8
    /*
    int reg = 7;
    for (int i = map.imin; i <= map.imax; i++){
            const int jmin = map.jinf[i];
            const int jmax = map.jsup[i];
            for (int j = jmin; j <= jmax; j++){
                    if (map.carte(i,j)){
    //if (t_count==1){

                            //1. No recruitment in EEZ
                            //if (map.maskEEZ[i][j] == -92){
                            if (map.maskEEZ[i][j] == -2){
                                  for (int a=a0_adult(sp); a<aN_adult(sp); a++)
                                         if (a>=param->age_mature[sp])
                                                mat.dvarDensity(0,a,i,j) = 0.0;

                             //      for (int a=0; a<=param->age_recruit[sp];
    a++)
                             //	       mat.dvarDensity(0,a,i,j) = 0.0;
                            }

    //}

                            //2. No recruitment outsize of EEZ
                            if (map.maskEEZ[i][j] != -8){
                                  if (i>=map.regimin[reg] && i<map.regimax[reg]
    && j>=map.regjmin[reg] && j<map.regjmax[reg]){ for (int a=0;
    a<=param->age_recruit[sp]; a++) mat.dvarDensity(0,0,i,j) = 0.0;
                                  }

                                  if (i>=map.regimin[5] && i<map.regimax[5] &&
                                      j>=map.regjmin[5] && j<map.regjmax[5]){
                                          for (int a=0;
    a<=param->age_recruit[sp]; a++) mat.dvarDensity(0,0,i,j) = 0.0;
                                  }
                            }


                    }
            }
    }*/
}
