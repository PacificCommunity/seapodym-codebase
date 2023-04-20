#include "SeapodymCoupled.h"

void SeapodymCoupled::ReadTimeSeriesData(int t, int t_series) {
    int nbytetoskip = (9 + (3 * nlat * nlon) + param->nlevel +
                       ((nlat * nlon) * (t_series - 1))) *
                      4;
    rw.rbin_input2d(param->strfile_pp, map, mat.np1[t], nbi, nbj, nbytetoskip);

    for (int k = 0; k < nb_layer; k++) {
        rw.rbin_input2d(
            param->strfile_u[k], map, mat.un[t][k], nbi, nbj, nbytetoskip);
        rw.rbin_input2d(
            param->strfile_v[k], map, mat.vn[t][k], nbi, nbj, nbytetoskip);
        rw.rbin_input2d(
            param->strfile_t[k], map, mat.tempn[t][k], nbi, nbj, nbytetoskip);
        if (!param->type_oxy)
            rw.rbin_input2d(
                param->strfile_oxy[k], map, mat.oxygen[t][k], nbi, nbj,
                nbytetoskip);
    }
    UnitConversions(t);
    if (param->use_sst) {
        rw.rbin_input2d(
            param->strfile_sst, map, mat.sst[t], nbi, nbj, nbytetoskip);
        const int imin = map.imin;
        const int imax = map.imax;
        for (int i = imin; i <= imax; i++) {
            const int jmin = map.jinf[i];
            const int jmax = map.jsup[i];
            for (int j = jmin; j <= jmax; j++) {
                if (map.carte[i][j]) {
                    // BUFFER ZONE to avoid biomass accumulation effect in case
                    // of artificially closed boundary Eventually move this
                    // buffer zone declaration into preparation of model forcing
                    // or topographic index
                    if (param->longitudeMin > 80 && param->deltaX >= 60.0 &&
                        param->deltaY >= 60.0)  // indicate Pacific ocean domain
                                                // and coarse resolution
                    // The bloc below is executed only in case of the Pacific
                    // ocean domain within the IndoPacific area.
                    {
                        // to make the buffer zone in IO part of the Pacific
                        // ocean domain. This allows avoiding the problems of
                        // high densities in shallow waters and complex current
                        // system of Indonesian region which is not resolved on
                        // coarse resolutions

                        if (i < param->lontoi(118) && j <= param->lattoj(-3)) {
                            if (mat.sst(t, i, j) > 22) mat.sst(t, i, j) = 22.0;
                            if (mat.tempn(t, 0, i, j) > 20.0)
                                mat.tempn(t, 0, i, j) = 20.0;
                            if (mat.tempn(t, 1, i, j) > 15.0)
                                mat.tempn(t, 1, i, j) = 15.0;
                        }
                        if (i < param->lontoi(130) && j > param->lattoj(3)) {
                            if (mat.sst(t, i, j) > 24) mat.sst(t, i, j) = 24.0;
                        }

                        if (i < param->lontoi(140) && j >= param->lattoj(-3)) {
                            if (mat.sst(t, i, j) > 24) mat.sst(t, i, j) = 24.0;
                            if (mat.tempn(t, 0, i, j) > 20.0)
                                mat.tempn(t, 0, i, j) = 20.0;
                            if (mat.tempn(t, 1, i, j) > 15.0)
                                mat.tempn(t, 1, i, j) = 15.0;
                        }
                        if (i < param->lontoi(146) && j >= param->lattoj(-8)) {
                            if (mat.sst(t, i, j) > 22) mat.sst(t, i, j) = 22.0;
                            if (mat.tempn(t, 0, i, j) > 18.0)
                                mat.tempn(t, 0, i, j) = 18.0;
                            if (mat.tempn(t, 1, i, j) > 15.0)
                                mat.tempn(t, 1, i, j) = 15.0;
                        }
                    }
                }
            }
        }
    } else
        mat.sst[t] = mat.tempn[t][0];

    if (param->use_vld) {
        rw.rbin_input2d(
            param->strfile_vld, map, mat.vld[t], nbi, nbj, nbytetoskip);
        mat.vld[t] =
            mat.vld[t] /
            1000.0;  // in km (mtl is in mt/km2). DO NOT change the units here
                     // without changing the use of vld in Ha computation
    } else
        mat.vld[t] = 1.0;  // wont be used

    // correct the T_epi temperature by the vertical gradieng magnitude
    if (param->sp_name[0].find("skj") == 0 && param->use_vld &&
        param->use_sst) {
        const int imin = map.imin;
        const int imax = map.imax;
        for (int i = imin; i <= imax; i++) {
            const int jmin = map.jinf[i];
            const int jmax = map.jsup[i];
            for (int j = jmin; j <= jmax; j++) {
                if (map.carte[i][j]) {
                    double dTdz = 2.0 *
                                  (mat.sst(t, i, j) - mat.tempn(t, 0, i, j)) /
                                  (1000.0 * mat.vld(t, i, j));
                    if (dTdz < 0.0) dTdz = 0.0;
                    if (dTdz > 0.2) dTdz = 0.2;
                    mat.tempn(t, 0, i, j) =
                        mat.tempn(t, 0, i, j) +
                        4.0 * dTdz * (mat.sst(t, i, j) - mat.tempn(t, 0, i, j));

                    mat.Hj(i, j) = mat.tempn(t, 0, i, j);
                }
            }
        }
    }
    if (!param->flag_coupling) {
        for (int n = 0; n < nb_forage; n++) {
            rw.rbin_input2d(
                param->strfile_F[n], map, mat.forage[t][n], nbi, nbj,
                nbytetoskip);

            mat.forage(t, n) =
                mat.forage(t, n) +
                0.000001;  // to avoid zero habitat index, expecially in fine
                           // resolution simulations F can be zero!!!
        }
    } else {
        for (int n = 0; n < nb_forage; n++)
            rw.rbin_input2d(
                param->strfile_S[n], map, mat.mats[n], nbi, nbj, nbytetoskip);
    }
}

void SeapodymCoupled::ReadClimatologyData(int t, int month) {
    int nbytetoskip =
        (9 + (3 * nlat * nlon) + 12 + ((nlat * nlon) * (month - 1))) * 4;
    rw.rbin_input2d(
        param->strfile_ppmc, map, mat.np1[t], nbi, nbj, nbytetoskip);
    for (int k = 0; k < nb_layer; k++) {
        rw.rbin_input2d(
            param->strfile_umc[k], map, mat.un[t][k], nbi, nbj, nbytetoskip);
        rw.rbin_input2d(
            param->strfile_vmc[k], map, mat.vn[t][k], nbi, nbj, nbytetoskip);
        rw.rbin_input2d(
            param->strfile_tmc[k], map, mat.tempn[t][k], nbi, nbj, nbytetoskip);
        if (!param->type_oxy)
            rw.rbin_input2d(
                param->strfile_oxymc[k], map, mat.oxygen[t][k], nbi, nbj,
                nbytetoskip);
    }
    UnitConversions(t);
    if (param->use_sst)
        rw.rbin_input2d(
            param->strfile_sstmc, map, mat.sst[t], nbi, nbj, nbytetoskip);
    else
        mat.sst[t] = mat.tempn[t][0];
    if (param->use_vld) {
        rw.rbin_input2d(
            param->strfile_vldmc, map, mat.vld[t], nbi, nbj, nbytetoskip);
        mat.vld[t] = mat.vld[t] / 1000.0;  // in km
    } else
        mat.vld[t] = 1.0;

    if (!param->flag_coupling) {
        for (int n = 0; n < nb_forage; n++)
            rw.rbin_input2d(
                param->strfile_Fmc[n], map, mat.forage[t][n], nbi, nbj,
                nbytetoskip);
    } else {  // only simulation mode (no optimization yet)
        for (int n = 0; n < nb_forage; n++)
            rw.rbin_input2d(
                param->strfile_Smc[n], map, mat.mats[n], nbi, nbj, nbytetoskip);
    }
}

void SeapodymCoupled::UnitConversions(
    int t) {  // Convert units of U and V from m/s to nmi/dt
    mat.un(t) = mat.un(t) * 3600 * 24 * deltaT / 1852;
    mat.vn(t) = mat.vn(t) * (-1) * 3600 * 24 * deltaT / 1852;
    const int imin = map.imin;
    const int imax = map.imax;
    for (int i = imin; i <= imax; i++) {
        const int jmin = map.jinf[i];
        const int jmax = map.jsup[i];
        for (int j = jmin; j <= jmax; j++) {
            if (map.carte[i][j]) {
                int nlayer = map.carte(i, j);
                for (int k = 0; k < nb_layer; k++) {
                    if (k < nlayer)
                        mat.un[t][k][i][j] =
                            mat.un[t][k][i][j] * mat.lat_correction[j];
                    else {
                        mat.un[t][k][i][j] = 0.0;
                        mat.vn[t][k][i][j] = 0.0;
                    }
                }
            }
        }
    }
}

void SeapodymCoupled::ReadClimatologyOxy(int t, int t_clm) {
    int nlevel_oxy = 12;
    if (param->type_oxy == 2) nlevel_oxy = 4;
    int nbytetoskip =
        (9 + (3 * nlat * nlon) + nlevel_oxy + ((nlat * nlon) * (t_clm - 1))) *
        4;
    for (int k = 0; k < nb_layer; k++)
        rw.rbin_input2d(
            param->strfile_oxy[k], map, mat.oxygen[t][k], nbi, nbj,
            nbytetoskip);
}

void SeapodymCoupled::ReadAll() {
    cout << "Reading all forcing variables at once... " << endl;
    int jday = 0;
    int t_count_init = t_count;

    // need to read oxygen in case if month==past_month
    //(otherwise we may not have it for the first time steps)
    if (param->type_oxy == 1 && month == past_month)
        ReadClimatologyOxy(1, month);
    // need to read oxygen in case if qtr==past_qtr
    if (param->type_oxy == 2 && qtr == past_qtr) ReadClimatologyOxy(1, qtr);

    for (; t_count <= nbt_total; t_count++) {
        getDate(jday);
        //----------------------------------------------//
        //	DATA READING SECTION: U,V,T,O2,PP	//
        //----------------------------------------------//
        if (t_count > nbt_building) {
            // TIME SERIES
            t_series = t_count - nbt_building + nbt_start_series;
            ReadTimeSeriesData(t_count, t_series);
        } else if ((t_count <= nbt_building) && (month != past_month)) {
            // AVERAGED CLIMATOLOGY DATA
            ReadClimatologyData(t_count, month);
        }
        if (param->type_oxy == 1 && month != past_month) {
            // MONTHLY O2
            ReadClimatologyOxy(t_count, month);
        }
        if (param->type_oxy == 2 && qtr != past_qtr) {
            // QUARTERLY O2
            ReadClimatologyOxy(t_count, qtr);
        }
    }
    t_count = t_count_init;
}

void SeapodymCoupled::RestoreDistributions(ivector& nb_age_built) {
    // to read initial distributions
    t_count += nbt_spinup_tuna;
    for (int sp = 0; sp < nb_species; sp++) {
        const int nb_ages = param->sp_nb_cohorts[sp];
        string fileCohorts =
            param->str_dir_init + param->sp_name[sp] + "_cohorts.dym";
        int nlevel = 0;
        int nlat = param->nlat;
        int nlon = param->nlong;
        rw.rbin_headpar(fileCohorts, param->nlong, param->nlat, nlevel);
        if (nlevel != nb_ages)
            cout << "WARNING: in file " << fileCohorts
                 << " number of cohorts is " << nlevel << " != " << nb_ages
                 << " in the current parfile!" << endl;
        for (int a = 0; a < nb_ages; a++) {
            int nbytetoskip =
                (9 + (3 * nlat * nlon) + nlevel + ((nlat * nlon) * a)) * 4;
            rw.rbin_input2d(
                fileCohorts, map, mat.init_density_species[sp][a], nbi, nbj,
                nbytetoskip);
        }
        if (param->nlat != nlat) {
            cerr << "WARNING: different nlat in file " << fileCohorts << " "
                 << param->nlat << endl;
            param->nlat = nlat;
        }
        if (param->nlong != nlon) {
            cerr << "WARNING: different nlon in file " << fileCohorts << " "
                 << param->nlong << endl;
            param->nlong = nlon;
        }
        nb_age_built[sp] = nb_ages - 1;
    }
}
