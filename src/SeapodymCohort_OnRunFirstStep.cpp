#include "SeapodymCoupled.h"
#include "SeapodymCohort.h"

//Prepare cohort run: initialize control variables, set flags, allocate memory for model and data variables, read forcing and fisheries data.
void SeapodymCohort::OnRunFirstStep()
{
	sumFprime.allocate(0, nb_forage - 1); 		sumFprime.initialize();
	sumF.allocate(0, nb_forage - 1);		sumF.initialize();
	sumF_area_pred.allocate(0, nb_forage - 1); 	sumF_area_pred.initialize();
	sumF_required_by_sp.allocate(0, nb_forage - 1);	sumF_required_by_sp.initialize();
	mean_omega_sp.allocate(0, nb_forage - 1);	mean_omega_sp.initialize();

	int Tr_step = 0;
        Date::init_time_variables(*param, Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, nbt_total,0,0);

	month = 0;
	Date::idatymd(param->ndatini, year, month, day);
	param->set_nbt(nbt_total);
	nbt_building = nbt_spinup_tuna;
	tf_cohort = tstart_cohort + nb_age_class;
	if (tf_cohort > nbt_total){
		tf_cohort = nbt_total;
	}
	int nbt_cohort = tf_cohort - tstart_cohort + 1;

	//Create time-dependent forcing matrices here:
	/*int t0  = 1;
	int offset_tstart_cohort = tstart_cohort;
	//to enable reading of forcing starting from t-1:
	if (cohort_id >= nb_age_class){
		t0 = 0;
		offset_tstart_cohort = tstart_cohort-1;
		nbt_building = -1;
	}
	int nbt = nbt_cohort;*/
	nbt_building = nbt_spinup_tuna;
	t_count = nbt_building+1;
	//Create time-dependent forcing matrices here:
	int t0  = t_count;
	int nbt = nbt_total;
	//nbt_building = -1;

	mat.createMatOcean(map, t0, nbt, nbi, nbj, nb_layer, deltaT);
	mat.createMatForage(map, nb_forage, t0, nbt, nbi, nbj);
	if (!param->larvae_input_aggregated_flag[0])
		mat.createMatLarvae(map, 1, nbt, nbi, nbj, deltaT);

	int nb_pops = nb_species*(1+param->nb_tag_files);
	mat.CreateMatSpecies(map,1, 1, nbi, nbj, nb_pops, a0_adult, param->sp_nb_cohorts);
	mat.CreateMatScalingFactorLarvae(map, nb_species);
	mat.CreateMatHabitat(map,nb_species,nb_forage,nb_layer,max(param->sp_nb_cohorts),1,1,nbi,nbj,a0_adult,param->sp_nb_cohorts,param->age_compute_habitat);
	past_month=0;
	past_qtr=0;
	sumP = 0; 

	for (int j=map.jmin;j<=map.jmax;j++){
		double lat = param->lastlat(j);
		mat.lastlat[j] = lat;
		mat.lat_correction[j] = param->correction_lat(lat);
		for (int jd=1; jd<=366; jd++){
			//if (param->sp_name[0].find("skj")==0)
       			//	mat.daylength[jd][j] = func.daylength(lat,jd); 
			//else
				//function using twilight hours increases day length
       				mat.daylength[jd][j] = func.daylength_twilight(lat,jd,18); 
		}
	}
	//Reading all forcing data for the cohort lifetime window
	//ReadAll(t0, nbt, offset_tstart_cohort);
	ReadAll(t_count, nbt_total, 0);

	Habitat.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	Mortality.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Spawning_Habitat.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Total_pop.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Habitat.initialize();
	Mortality.initialize();
	Spawning_Habitat.initialize();
	Total_pop.initialize();
	mat.mats.initialize();


	rw.init_writing(*param);

	param->fdata_rm = 0;
	if (!param->flag_no_fishing){
		//Initialization of fishery files variables
		//1. Read catch and effort data
		rw.rtxt_fishery_data(*param,map,nbt_total,jday_spinup);
		//redistribute all fishing effort to the model resolution
		rw.set_effort_rm(*param,map,nbt_total,jday_spinup);
		//put all data on the model resolution in case of MPA simulations or EEZ extractions
		
		if ((min(param->fishery_reso) < param->catch_reso) 
		     && !param->mpa_simulation && !param->nb_EEZ){
			rw.degrade_fishery_reso(*param, map,nbt_total,jday_spinup);
		}
		//read LF data file if provided
		if (param->file_frq_data[0]!=""){
			bool writeobs = (!param->gcalc() && !param->scalc());
			for (int sp=0; sp<nb_species; sp++)
				rw.read_frq_data(*param, map, param->save_first_yr, param->save_last_yr, sp, writeobs);
		}
	}
	else {
		/*cout << "----------------------------------------------------" << endl;
		cout << "            MODEL RUN WITHOUT FISHING" << endl;
		cout << "----------------------------------------------------" << endl;*/
	}
	func.allocate_dvmatr(map.imin,map.imax,map.jinf,map.jsup);

	//vector of zeroes, for accessibility function
	tags_age_habitat.allocate(0,aN_adult(0));
	tags_age_habitat.initialize();
	
	
	pop.time_reading_init();
	func.time_reading_init();
	param->time_reading_init();


	//Temporarily reading the tau of the first cohort
	//Need to be just a single number for all cohorts
	dtau = param->sp_unit_cohort[0][1];
	nbt_before_first_recruitment = Date::get_nbt_before_first_recruitment(
			param->first_recruitment_date,
			param->ndatini,param->deltaT,param->date_mode); 	
	//counter of number of time steps between recruitments (survival equations)
	//once nt_dtau = dtau, recruitment occurs and nt_dtau=0
	//its initial value is dtau-nbt_before_first_recruitment_date
	nt_dtau = dtau-nbt_before_first_recruitment; 

	//routine-specific variables
	fishing = false;
	migration_flag = 0;
	step_count= 0;
	step_fishery_count= 0;
	jday = 0; 
	nbstoskip = param->nbsteptoskip; // nb of time step to skip before computing likelihood

	likelihood = 0.0;

	pop_built = 1;

	past_month=month;
	past_qtr=qtr;
	
}
