#include "SeapodymCoupled.h"

void SeapodymCoupled::WriteOutput(int t, bool fishing)
{
	//Let's compute PP in three areas:
	dvector sumPP(0,3); sumPP.initialize();
	const double area = 1852*deltaX*1852*deltaY;
	for (int j=map.jmin;j<=map.jmax;j++){
	//for (int j=1;j<nbj-1;j++){
		double lat_corrected_area = area / mat.lat_correction[j];
		if (param->jtolat(j)>10.0 && param->jtolat(j)<=45.0){
			for (int i=1;i<nbi-1;i++)
				if (map.carte(i,j))
					sumPP[0] += mat.np1(t,i,j)*lat_corrected_area;
		}
		if (param->jtolat(j)>=-10.0 && param->jtolat(j)<=10.0){
			for (int i=1;i<nbi-1;i++)
				if (map.carte(i,j))
					sumPP[1] += mat.np1(t,i,j)*lat_corrected_area;
		}
		if (param->jtolat(j)>=-35.0 && param->jtolat(j)<-10.0){
			for (int i=1;i<nbi-1;i++)
				if (map.carte(i,j))
					sumPP[2] += mat.np1(t,i,j)*lat_corrected_area;
		}
	}
        sumPP[3] = sum(DtoBcell(1e6*mat.np1(t)));
//cout << sumPP << endl;
        //sumP = sum(DtoBcell(1e6*mat.np1));
        for (int n=0; n<nb_forage; n++){
                sumF[n] = sum(DtoBcell(mat.forage[t][n]));
        }

	//CalcSums();

	if (t_count > nbt_building){	
		/*to rewrite as the old code with binary fishing data files is now desactivated
		for (int sp=0; sp< nb_species;sp++){
			if (param->wbin_flag){
				// Write files with predicted fisheries data
				rw.write_fishery_data(*param, map, mat, sp, year, month, true);
				if (month==3 || month==6 || month==9 || month==12)
					rw.write_frq_data(*param, sp, year, qtr, value(mat.dvarLF_est(sp)), true);
			}
		}*/
		//save binary files with population distributions
		rw.SaveSepodymFileDym(*param, map, mat);


		if (qtr != past_qtr)
			mat.C_N_sp_age_fishery_qtr.initialize();
		//mat.C_N_sp_age_fishery_qtr = mat.C_N_sp_age_fishery;
		if (month==3 || month==6 || month==9 || month==12){
			for (int sp=0; sp< nb_species;sp++)
				for (int age=a0_adult[sp];age<aN_adult[sp];age++)
					for (int k=0; k<param->nb_fishery_by_sp[sp];k++)
						for (int r=0; r<param->nb_region_sp_B[sp];r++)
							mat.C_N_sp_age_fishery_qtr[sp][age][r][k] = mat.C_N_sp_age_fishery[sp][age][k][r];
		}
	}
	//save ascii file with sums
	rw.SaveSepodymFileTxt(*param, mat, map, sumPP, sumF, sumFprime, sumF_area_pred, sumF_required_by_sp, 
				mean_omega_sp, day, month, year, t_count,past_qtr, qtr, nbi, nbj);
        if (param->write_all_cohorts_dym){
                for (int sp=0; sp< nb_species;sp++){
                                dvector zlevel; zlevel.allocate(0, nbt_total - 1);
                                SaveCohortsDym(sp, false, zlevel);
                }
        }

	//For Joe:
	//dvector zlevel; zlevel.allocate(0, nbt_total - 1);
	//SaveOneCohortDym(0, false, zlevel);
}

void SeapodymCoupled::WriteFluxes()
{
	//Note, to put frq_like=0 and stop if regions are not defined and activated in Fluxes comp mode, then remove this condition
	if (!param->frq_like[0] && param->nb_region>0){
		if (fluxes_dt_qtr && (month==3 || month==6 || month==9 || month==12))
			rw.SaveFluxesCohortsFileTxt(*param, mat, map, day, month, year);
	}	
}

void SeapodymCoupled::WriteFileHeaders()
{
/*
	int nbt_yr = (int) (365.25/deltaT);
	int n_tt   = (int) (nbt_yr*(param->save_last_yr - param->save_first_yr + param->nb_yr_forecast) + 1.01);
	float dec  = param->save_first_yr - (int) param->save_first_yr ;
	dec = ((int)(dec/t_yrdd))* t_yrdd ;
	float firstdate_rec = ((int) param->save_first_yr) + t_yrdd/2 +  dec;

	dvector zlevel;
	zlevel.allocate(0, n_tt - 1);
	for (int n=0; n< n_tt; n++)
		zlevel[n]= firstdate_rec + (float)(n* t_yrdd);
*/
	dvector zlevel;
        zlevel.allocate(0, nbt_total - 1);
        zlevel.initialize();
        Date::zlevel_run(*param,mat.zlevel,nbt_total,zlevel,nbt_start_series);

	// Create and initialize (dym) files for saving spatial variables
	//rewrite dym mask by the mask used in the run, i.e. map.carte:
	for (int i=0; i<nbi-2; i++){
	        for (int j=0; j<nbj-2; j++){
			mat.mask[j][i] = map.carte[i+1][j+1];
		}
	}	
	rw.InitSepodymFileDym(*param, mat, nbt_total, zlevel, mat.mask);

        if (param->write_all_cohorts_dym)
                for (int sp=0; sp< nb_species;sp++)
                        SaveCohortsDym(sp, true, zlevel);

	//For Joe:
	//SaveOneCohortDym(0, true, zlevel);


	/*to rewrite as the old code with binary fishing data files is now desactivated
	if (param->wbin_flag){
		// Write header into files with predicted fisheries data
		rw.write_fishery_data(*param, map, mat, 0, 0, 0, false);
		rw.write_frq_data(*param, 0, 0, 0, value(mat.dvarLF_est(0)), false);
	}
	*/

	// Create and initialize (txt) files for saving aggregated variables
	rw.InitSepodymFileTxt(*param);

}

void SeapodymCoupled::InitFileFluxes()
{
	//Create and initialize (txt) files with cohort fluxes through regional boundaries
	//WARNING: single matrix for all species! To be modified in case of multi-species simulations
	if (!param->frq_like[0] && param->nb_region>0){
		fluxes_dt_qtr = 1;
		const int nb_reg = param->nb_region_sp_B[0];//attn sp =0!
		const int nb_ages= param->sp_nb_cohorts[0];//attn sp =0!
		Density_region.allocate(0,nb_reg-1);
		for (int r=0; r<nb_reg; r++){
			Density_region(r).allocate(0,nb_ages-1);
        		for (int a=0; a<nb_ages; a++){
				Density_region(r,a).allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
				Density_region(r,a).initialize();
			}
		}

		mat.createMatFluxes(param->nb_region_sp_B[0], param->sp_nb_cohorts[0]);
		rw.InitFluxesCohortsFileTxt(*param);
	}	
}

void SeapodymCoupled::SaveCohortsDym(int sp, bool WriteHeader, dvector zlevel)
{
        //To remove estimates in Indian ocean (boundary effect)
        //int ilim = param->lontoi(125);
        //int jlim = param->lattoj(-8);

        const int nb_ages= param->sp_nb_cohorts[sp];
        for (int a=0; a<nb_ages; a++){

		double minval = min(value(mat.dvarDensity(sp,a)));
                double maxval = max(value(mat.dvarDensity(sp,a)));
		
                std::ostringstream ostr;
                ostr << a+1;
                string fileout = param->strdir_output + param->sp_name[sp] + "_age" + ostr.str() + ".dym";
                if (WriteHeader){
                        //Write file headers during the first time step
                        rw.wbin_header(fileout, param->idformat, param->idfunc, minval, maxval,
                                        param->nlong, param->nlat, nbt_total,
                                        zlevel[0], zlevel[nbt_total-1],
                                        mat.xlon, mat.ylat, zlevel, mat.mask);

                }
                else {
                        //Append data for the current date
                        dmatrix mat2d(0, nbi - 1, 0, nbj - 1);
                        mat2d.initialize();
                        for (int i=map.imin; i <= map.imax; i++){
                                for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
                                        if (map.carte[i][j]){
                                                //Units: Nb/sq.km
                                                mat2d(i-1,j-1) = value(mat.dvarDensity(sp,a,i,j));
                                                //if (i<ilim && j>jlim) mat2d[i-1][j-1] = 0.0;
                                        }
                                }
                        }
                        rw.wbin_transpomat2d(fileout, mat2d, nbi-2, nbj-2, true);
			//update min-max values in header
			rw.rwbin_minmax(fileout, minval, maxval);
                }
        }
}

void SeapodymCoupled::SaveOneCohortDym(int sp, bool WriteHeader, dvector zlevel)
{//First created for Joe (SimpoDym) at 20160229
	double minval = min(value(mat.dvarDensity(sp,3)));
        double maxval = max(value(mat.dvarDensity(sp,3)));
	
        string fileout = param->strdir_output + param->sp_name[sp] + "_cohort_density.dym";
        if (WriteHeader){
        	//Write file headers during the first time step
		rw.wbin_header(fileout, param->idformat, param->idfunc, minval, maxval,
                                        param->nlong, param->nlat, nbt_total,
                                        zlevel[0], zlevel[nbt_total-1],
                                        mat.xlon, mat.ylat, zlevel, mat.mask);

        }
        else {
        	const int nb_ages= param->sp_nb_cohorts[sp];
		int age_extract = 3;
	        for (int age=0; age<nb_ages; age++){

			if (age == t_count+2 && age < param->sp_nb_cohorts[sp]-1) {
				cout << "age to extract " << age << endl; 
				age_extract = age;
			}
			if (t_count+2 >= param->sp_nb_cohorts[sp]-1 && age == param->sp_nb_cohorts[sp]-1){
			       	cout << "age to extract " << age << endl;
				age_extract = age;
			}
		}
	
	        //Append data for the current date
	        dmatrix mat2d(0, nbi - 1, 0, nbj - 1);
	        mat2d.initialize();
	        for (int i=map.imin; i <= map.imax; i++){
	               	for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
	                       	if (map.carte[i][j]){
	                               	//Units: Nb/sq.km
	                                mat2d(i-1,j-1) = value(mat.dvarDensity(sp,age_extract,i,j));

	                        }
	                }
		}
	        rw.wbin_transpomat2d(fileout, mat2d, nbi-2, nbj-2, true);
		//update min-max values in header
		rw.rwbin_minmax(fileout, minval, maxval);
        }
}


