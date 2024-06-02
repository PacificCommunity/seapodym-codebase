#include "SeapodymCoupled.h"

void SeapodymCoupled::prerun_model()
{
	OnRunFirstStep();
	if (param->larvae_like[0]){
		ReadLarvae();
	}
}

///This is the main loop function. It includes the following calls:
///1- Initialising population density 
///2- Reading all forcing data (once in optimization mode, at every time step in simulation mode) 
///3- Reading fisheries data
///4- Reading tagging data if tag_like is activated
///5- Age/lifestage loop calling the ADRE solvers and ageing
///6- Predicting observed variables (catch, LF and density of tags)
///7- Likelihood computation
///8- Writing outputs (in simulation mode only)


//////////////////////////////////////////////////////////////////
//--------------------------------------------------------------//
//		     FORAGE-TUNA SIMULATION			//
//--------------------------------------------------------------//
//////////////////////////////////////////////////////////////////
/*!
\brief The tuna population main loop is in this function.
*/
double SeapodymCoupled::OnRunCoupled(dvar_vector x, const bool writeoutputfiles)
{
	InitializeAll();
	SUM_CATCH = 0.0;

	past_month=month;
	past_qtr=qtr;

	//Temporarily reading the tau of the first cohort
	//Need to be just a single number for all cohorts
	int dtau = param->sp_unit_cohort[0][1];
	int nbt_before_first_recruitment = Date::get_nbt_before_first_recruitment(
			param->first_recruitment_date,
			param->ndatini,param->deltaT,param->date_mode); 	
	//counter of number of time steps between recruitments (survival equations)
	//once nt_dtau = dtau, recruitment occurs and nt_dtau=0
	//its initial value is dtau-nbt_before_first_recruitment_date
	int nt_dtau = dtau-nbt_before_first_recruitment; 

	//routine-specific variables
	int tcur = t_count; //will be used for forcing variable time control
	int nbt_no_forecast = t_count + nbt_spinup_tuna + nbt_total - 1;
	bool fishing = false;
	bool connectivity_comp = param->connectivity_comp;
	int migration_flag = 0;
	int step_count= 0;
	int step_fishery_count= 0;
	int jday = 0; 
	int nbstoskip = param->nbsteptoskip; // nb of time step to skip before computing likelihood

	if (!param->gcalc()){
		//need to read oxygen in case if month==past_month
		//(otherwise we may not have it for the first time steps)
		if (param->type_oxy==1 && month==past_month)
			ReadClimatologyOxy(1, month);
		//need to read oxygen in case if qtr==past_qtr 
		if (param->type_oxy==2 && qtr==past_qtr)
			ReadClimatologyOxy(1, qtr);
	}

	//----------------------------------------------//
	// 	LIKELIHOOD INITIALISATION SECTION       //
	//----------------------------------------------//	
	lflike = 0.0;
	double taglike = 0;
	double stocklike = 0.0;
	double larvaelike = 0.0;
	dvariable likelihood = 0.0;
	dvariable total_stock = 0.0;
	//Reset model parameters:
	reset(x);

	//----------------------------------------------//
	// 	LOCAL MATRICES ALLOCATION SECTION       //
	//----------------------------------------------//	
	dvar_matrix Spawning_Habitat;
	dvar_matrix Total_pop;
	dvar_matrix Habitat; 
	dvar_matrix IFR; 
	dvar_matrix ISR_denom; 
	dvar_matrix FR_pop;
	dvar_matrix Mortality; 

	Habitat.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	Mortality.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Spawning_Habitat.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Total_pop.allocate(map.imin, map.imax, map.jinf, map.jsup);

	if (param->food_requirement_in_mortality(0)){ 
		//temporal, need to check memory use first 
		IFR.allocate(map.imin, map.imax, map.jinf, map.jsup);
		ISR_denom.allocate(map.imin, map.imax, map.jinf, map.jsup);
		FR_pop.allocate(map.imin, map.imax, map.jinf, map.jsup);
		IFR.initialize();
		ISR_denom.initialize();
		FR_pop.initialize();
	}
	Spawning_Habitat.initialize();
	Habitat.initialize();
	Mortality.initialize();

	// For aggregated larvae likelihood
	dvar_matrix Larvae_density_pred_at_obs;
	IVECTOR kinf;
	IVECTOR ksup;
	int nb_larvae_input_agg_groups = param->nb_larvae_input_agg_groups;
	int ntime_agg[nb_larvae_input_agg_groups];
	if (param->larvae_like[0]){
		for (int i = 0; i < nb_larvae_input_agg_groups; ++i) {
			ntime_agg[i] = 0;
		}	
		if (param->larvae_input_aggregated_flag[0]==1 && param->larvae_like[0]){
			// Aggregated larvae density over the entire period, only at obs. locations
			kinf.allocate(0, nb_larvae_input_agg_groups-1);
			ksup.allocate(0, nb_larvae_input_agg_groups-1);
			for (int k = 0; k < nb_larvae_input_agg_groups; k++){
				kinf[k] = 0;
				ksup[k] = mat.aggregated_larvae_input_vectors[k].size();
			}
			Larvae_density_pred_at_obs.allocate(0, nb_larvae_input_agg_groups-1, kinf, ksup);
			Larvae_density_pred_at_obs.initialize();	
		}
	}

	//precompute thermal habitat parameters
	for (int sp=0; sp < nb_species; sp++)
		func.Vars_at_age_precomp(*param,sp);

	//precompute seasonal switch function
	for (int sp=0; sp < nb_species; sp++){
		if (param->seasonal_migrations[sp]){
			func.Seasonal_switch_year_precomp(*param,mat,map,
						value(param->dvarsSpawning_season_peak[sp]),
						value(param->dvarsSpawning_season_start[sp]),sp);
		}
	}

	//TAG vars initialization section
	tags_age_habitat.initialize();
	if (nb_tagpops>0){
		tagpop_age_solve.initialize();
		rec_pred.initialize();
	}

	if (writeoutputfiles){
		WriteFileHeaders();
		if (!param->gcalc())
			ConsoleOutput(0,0);
	}
	//Simulation FLAGS
	bool tags_only = param->tags_only;
	int spop0 = 0;
	int mortality_off = 0;
	if (tags_only) {
		spop0 = 1;
		nb_fishery = 0;
		fishing = false;
		param->stock_like.initialize();
		param->frq_like.initialize();
		param->larvae_like.initialize();
	}
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	////------------------------------------------------------------/////
	////								/////
	////		||| START OF SIMULATION CYCLE |||		/////
	////								/////
	////------------------------------------------------------------/////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	//Add penalty function to the likelihood for sum(eF_habitat)<eF_sum
	//Not used currently
	double eFlike = 0.0;
/*	if (param->tag_like[0]){
		dvariable dvarEF_sum = sum(param->dvarsEF_habitat);
		likelihood -= 1e1*log(eF_sum-dvarEF_sum);
		eFlike -= 1e1*log(eF_sum - value(dvarEF_sum));
	}
*/
	for (;t_count <= nbt_total; t_count++)
	{
		//----------------------------------------------//
		//              INITIALISATION                  //
		//----------------------------------------------//
		sumP=0;
		sumFprime.initialize();
		sumF.initialize();
		mat.dvarCatch_est.initialize();
		if (param->tag_like[0]){
			get_tag_releases(mat.dvarDensity,tagpop_age_solve,nb_rel);
			tags_age_habitat.initialize();
			for (int spop=1; spop<nb_tagpops+1; spop++){
				for (int a=a0_adult(0); a<aN_adult(0); a++){
					//first initialize with the ages, which will be resolved
					tags_age_habitat[a] += tagpop_age_solve(spop-1,t_count-1,a);
				}
			}
			//second correct by the age_compute_habitat indices
			for (int a=a0_adult(0); a<aN_adult(0); a++){
				if (tags_age_habitat[a]){
					int aa = a;
					for (; aa>=a0_adult(0); aa--){
						if (param->age_compute_habitat(0,aa)!=param->age_compute_habitat(0,aa-1))
							break;
					}
					if (!tags_age_habitat[aa]) tags_age_habitat[aa] = -1;
				}
			}
		}
		//----------------------------------------------//
		//		       DATE			//
		//----------------------------------------------//
		getDate(jday);
		for (int sp=0; sp < nb_species; sp++){
			func.Seasonal_switch(*param,mat,map,jday,sp);	
		}
		tcur = t_count;
		if (!param->gcalc()){
			tcur = 1; 
			//----------------------------------------------//
			//	DATA READING SECTION: U,V,T,O2,PP	//
			//----------------------------------------------//
			//if ((t_count > nbt_building)) { //CC run with average effort forecast
			if ((t_count > nbt_building) && (t_count <= nbt_no_forecast)) {
				//TIME SERIES 
				t_series = t_count - nbt_building + nbt_start_series;
				ReadTimeSeriesData(tcur,t_series);	
	
			}
			else if (((t_count <= nbt_building) && (month != past_month)) || (t_count > nbt_no_forecast)) {
	
				//AVERAGED CLIMATOLOGY DATA
				ReadClimatologyData(tcur, month);
			}
			if (param->type_oxy==1 && month != past_month) {
				//MONTHLY O2
				ReadClimatologyOxy(tcur, month);
			}
			if (param->type_oxy==2 && qtr != past_qtr) {
				//QUARTERLY O2
				ReadClimatologyOxy(tcur, qtr);
			}
		}
		
		//----------------------------------------------//
		//	READING CATCH DATA: C_obs, effort	//
		//----------------------------------------------//
		if (!tags_only && (nb_fishery > 0) && sum(param->mask_fishery_sp) && ( t_count > nbt_building) 
				&& (year>=(int)param->save_first_yr) && (t_count <= nbt_no_forecast)){
			mat.catch_est.initialize();
			rw.get_fishery_data(*param, mat.effort, mat.catch_obs, mat.efflon, mat.efflat, year, month);
			fishing = true;
			if (qtr != past_qtr){
				mat.C_N_sp_age_fishery.initialize();
				mat.dvarLF_est.initialize();
				if (t_count > nbt_building+nbstoskip){
					rw.get_LF_qtr_data(*param, mat.LF_qtr_obs, year, qtr-1);
				}
			}
		}
		//MPAs are no longer supported in the code. Prepare the fishery data as designed and use as input data
                if (param->mpa_simulation){
                        rw.get_fishery_data_mpa(map,*param, mat.effort, mat.catch_obs, mat.efflon, mat.efflat, year, month);
			
                }

		if (t_count > nbt_no_forecast) {
			rw.get_average_effort(*param, mat.effort, mat.efflon, mat.efflat, 10, month);
			//to get actual fishing data if it exists
			//rw.get_fishery_data(*param, mat.effort, mat.catch_obs, year, month);
			mat.catch_obs.initialize();
			mat.catch_est.initialize();
			mat.C_N_sp_age_fishery.initialize();
			mat.dvarLF_est.initialize();
		}
		//----------------------------------------------//
		//	 SPIN-UP PARAMETERS (To be removed)	//
		//----------------------------------------------//
		int pop_built = 0; 
		if (t_count > nbt_building) 
			pop_built = 1;
		//------------------------------------------------------------------------------//
		//	TRANSPORT OF TUNA AGE CLASSES AND PREDICTED CATCH COMPUTATION		//
		//------------------------------------------------------------------------------//
		//------------------------------------------------------------------------------//
		mat.u = mat.un[tcur][0]; mat.v = mat.vn[tcur][0]; 
		//Precompute diagonal coefficients for larvae and juvenile ADREs
		pop.precaldia(*param, map, mat);
		pop.caldia(map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y, mat.advection_y);
		for (int sp=0; sp < nb_species; sp++){
			//store fish density before transport
			for (int a=a0_adult(sp); a<aN_adult(sp); a++)
				mat.density_before(sp,tcur,a) = value(mat.dvarDensity(sp,a));

			//1. Precompute some variables outside of age loop

			//1.1 Accessibility by adults (all cohorts)
 			func.Faccessibility(*param, mat, map, sp, jday, tcur, pop_built, tags_only, tags_age_habitat);//checked

			//1.2 Food Requirement by population
			if (param->food_requirement_in_mortality(sp)){
				FR_pop_comp(FR_pop, sp);
			//	ISR_denom_comp(ISR_denom, sp, tcur);
			}
			//1.3 precompute local proportions of catch by fishery
			if (fishing && param->fisheries_no_effort_exist[sp]){
				pop.Ctot_proportion_fishery_comp(map,*param,mat,rw,year,month,sp);
			}
			//1.4 precompute selectivity functions
			if (fishing)
				pop.Selectivity_comp(*param,nb_fishery,a0_adult(sp),aN_adult(sp),sp);

			//1.5 Connectivity simulations by setting the density in donor region to zero
			if (connectivity_comp)
				Set_density_region_age_zero(sp);

			//2. IMPLICIT AGE LOOP: increment age while moving through life stages 
			int age = 0;	

			//2.1 Spawning habitat	
			if (!tags_only)
				func.Spawning_Habitat(*param, mat, map, Spawning_Habitat, 1.0, sp, tcur, jday);
			//2.2 Transport and mortality of larvae (always one age class)	
			for (int n=0; n<param->sp_nb_cohort_lv[sp]; n++){
				double mean_age = mean_age_cohort[sp][age]; 
				if (!tags_only){
					func.Mortality_Sp(*param, mat, map, Mortality, Spawning_Habitat, sp, mean_age, age, tcur);
					pop.Precalrec_juv(map, mat, Mortality, tcur);//checked
					pop.Calrec_juv(map, mat, mat.dvarDensity[sp][age], Mortality, tcur);//checked
				}
				age++;
			}
			///if (t_count > nbt_spinup_forage + nt_jv){ SPINUP TO BE FIXED OR REMOVED!!!
			//2.3. Juvenile habitat	
			if (!tags_only){
				if (param->cannibalism[sp]){
					Total_Pop_comp(Total_pop,sp,jday,tcur); //adjoint
					func.Juvenile_Habitat_cannibalism(*param, mat, map, Habitat, Total_pop, sp, tcur);
				} else 
					func.Juvenile_Habitat(*param, mat, map, Habitat, sp, tcur);
			}

			//2.4. Transport and mortality of juvenile age classes	
			for (int n=0; n<param->sp_nb_cohort_jv[sp]; n++){			
				double mean_age = mean_age_cohort[sp][age];
				if (!tags_only){
					func.Mortality_Sp(*param, mat, map, Mortality, Habitat, sp, mean_age, age, tcur);
					pop.Precalrec_juv(map,  mat, Mortality, tcur);
					pop.Calrec_juv(map, mat, mat.dvarDensity[sp][age], Mortality, tcur);
				}
				age++;
			}
			///} 
			
			//----------------------------------------------//
			//	Fishing with effort data available	//
			//----------------------------------------------//
			//3. Compute predicted catch before transport & mortality	
			if (fishing){
				int k = 0;
				for(int f=0; f<nb_fishery; f++){
					if (param->mask_fishery_sp[sp][f]){
						//In Optimization Mode compute predicted (E-based) catch and LF only for
						//fisheries designed to inform parameters (flag mask_fishery_sp_like).  
						//In Simulation Mode - for all fisheries with available effort data.
						if (param->gcalc()){
							if (param->mask_fishery_sp_like[sp][f] && !param->mask_fishery_sp_no_effort[sp][f])
								pop.Predicted_Catch_Fishery(map,*param,mat,rw,sp,f,k,year,month,tcur,step_fishery_count);
						} else {
							if (!param->mask_fishery_sp_no_effort[sp][f])
								pop.Predicted_Catch_Fishery(map,*param,mat,rw,sp,f,k,year,month,tcur,step_fishery_count);
						}			
						k++;
					}
				}
				step_fishery_count++;
			}
			if (fishing && param->fisheries_no_effort_exist(sp))
				pop.Total_exploited_biomass_comp(map,*param,mat,sp,tcur);

			//4. Transport and mortality of adult cohort
			///if (t_count > nbt_spinup_forage + nt_yn){ TO BE FIXED!!!
			///for (int age=0; age<=nb_age_built[sp]; age++){///TO BE FIXED!!!	
			for (int n=0; n<param->sp_nb_cohort_ad[sp]; n++){

				if (fishing && param->fisheries_no_effort_exist(sp))
					pop.Total_obs_catch_age_comp(map,*param,mat,rw,age,sp,year,month,tcur);


				//NOTE: currently current averaging doesn't depend on seasonal migrations
				if (param->vert_movement[sp] && param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]){
					if (!tags_only || tags_age_habitat(age))
						func.Average_currents(*param, mat, map, age, tcur, pop_built);
				}
				
				//2014: catch at age computation for fisheries without effort data
				//pop.Ctot_no_effort_sp_age_comp(map, *param, mat, value(mat.dvarDensity(sp,age)), tcur, sp, age);				
				//this section will work only if seasonality switch is ON 
				//and for <1 maturity at age parameter
				if (param->migrations_by_maturity_flag && param->seasonal_migrations[sp]){
					cout << "In this version no smooth maturity; enter the age at first maturity. Exit now!" << endl; exit(1);
					
				} // end of section with seasonality switch and <1 maturity at age parameter 
				else {
					migration_flag = 0;
					if (age>=param->age_mature[sp] && param->seasonal_migrations[sp]) migration_flag = 1;

					if (param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]) {
						if (!tags_only || tags_age_habitat(age))
							func.Feeding_Habitat(*param,mat,map,Habitat,sp,age,jday,tcur,migration_flag);
					}
					double mean_age = mean_age_cohort[sp][age];

					if (!tags_only){
					if (!param->food_requirement_in_mortality(sp)){
						func.Mortality_Sp(*param, mat, map, Mortality, Habitat, sp, mean_age, age, tcur);//checked
					} else {
						Food_Requirement_Index(IFR, FR_pop, ISR_denom, sp, age, tcur, jday);
						func.Mortality_Sp(*param, mat, map, Mortality, IFR, sp, mean_age, age, tcur);//checked
					}}

					if (param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]){
						if (!tags_only || tags_age_habitat(age))
							pop.Precaldia_Caldia(map, *param, mat, Habitat, Total_pop, sp, age, tcur, jday);//checked	
					}
					if (!param->gcalc()){
						// only in simulation mode: compute mean speed 
						// in BL/sec and mean diffusion rate in nmi^2/day
						mat.MeanVarMovement(map,value(mat.dvarsAdvection_x),
							value(mat.dvarsAdvection_y),
				 			value(mat.dvarsDiffusion_y),
							param->MSS_species[sp],
							param->sigma_species[sp],
							param->length(sp,age),
							param->length(sp,param->sp_nb_cohorts[sp]-1),
							deltaT,sp,age);
					}
						
					//store adult habitat and update the counter
					if (param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]){
						//cout << age << " " << param->age_compute_habitat[sp][age] << endl;
						mat.adult_habitat(sp,tcur,param->age_compute_habitat[sp][age]) = value(Habitat);
					}
					for (int spop=spop0; spop<nb_tagpops+1; spop++){
						mortality_off = 0;
						if (sum(param->mask_fishery_sp) && !tags_only)
							fishing = true;

						if (spop>0) { 
							Mortality.initialize();
							mortality_off = 1;
							if (fishing) fishing = false;
						}
						if (spop>0){
							if (!tagpop_age_solve(spop-1,t_count-1,age))
								continue;
						}
						pop.Precalrec_Calrec_adult(map,mat,*param,rw,
							mat.dvarDensity[spop][age],Mortality,
							tcur,fishing,age,sp,year,month,jday,
							step_fishery_count,mortality_off);//checked 20150210
					}
					if (sum(param->mask_fishery_sp) && !tags_only) fishing = true;
				}
				age++;
			}//END OF IMPLICIT 'age' LOOP

			//----------------------------------------------//
			//    Now get predicted catch without effort  	//
			//----------------------------------------------//
			if (fishing && param->fisheries_no_effort_exist(sp))
				pop.Predicted_Catch_Fishery_no_effort(map,*param,mat,rw,sp,year,month);

			//store fish density after transport and mortality
			for (int a=0; a<aN_adult(sp); a++)
				mat.density_after(sp,a) = value(mat.dvarDensity(sp,a));
			if (writeoutputfiles){
				CalcMeanTemp(t_count,tcur);
				CalcSums();
			}
			//Compute total stock before the new recruitment (survival)
			if (param->stock_like[sp] && t_count > nbt_building+nbstoskip)
				Total_Stock_comp(total_stock, sp);

			//5. Compute spawning biomass (sum of young and adults density weighted by maturity-at-age)
			if (!tags_only)
				SpawningBiomass_comp(Total_pop, sp);
			//----------------------------------------------//
			//	    TUNA AGEING AND SPAWNING		//
			//----------------------------------------------//

			//6. Ageing and survival
			//6.1. First, update the age structure of tagged cohorts
			if (nt_dtau==dtau){
				for (int a=param->sp_nb_cohorts[sp]-1; a >= 1; a--){
					for (int spop=spop0; spop<nb_tagpops+1; spop++){
						if (spop>0){
							if (tagpop_age_solve(spop-1,t_count-1,a-1)==1){					
								tagpop_age_solve(spop-1,t_count,a-1)=0;
								tagpop_age_solve(spop-1,t_count,a)=1;
							}
						}
						//6.2. Ageing of population density
						Survival(mat.dvarDensity[spop][a], mat.dvarDensity[spop][a-1] , a, sp);
					}
				}
				nt_dtau=0;
			}
			//7. Spawning
			if (!tags_only)
				Spawning(mat.dvarDensity[sp][0],Spawning_Habitat,Total_pop,jday,sp,pop_built,tcur);//checked

			//8. Aggregate larvae density at larvae obs locations
			if (t_count > nbt_building+nbstoskip){
				if (param->larvae_input_aggregated_flag[0]==1 && param->larvae_like[0]){
					const int nb_lv = param->sp_nb_cohort_lv[sp];
					int iAgg = Utilities::iTimeOfYear(month, param->larvae_input_aggregation);

					// Calculate scaling factor between larvae density at 1st time step and larvae density at age_larvae_before_sst_mortality days, based on sst-dependent mortality (if applicable)
					if (param->larvae_mortality_sst[sp]){
						func.Scaling_factor_sstdep_larvae_mortality(*param, mat.sst[tcur], map, mat.dvarScaling_factor_sstdep_larvae_mortality[sp], sp);
					}

					// Aggregate larvae density at larvae obs locations
					for (auto k=0u; k<mat.aggregated_larvae_input_vectors[iAgg].size(); k++){
						if (param->larvae_mortality_sst[sp]){// In this case, it only considers the first age class (even if nb_lv>1)
							Larvae_density_pred_at_obs(iAgg, k) +=  mat.dvarDensity[sp][0][mat.aggregated_larvae_input_vectors_i[iAgg][k]][mat.aggregated_larvae_input_vectors_j[iAgg][k]] * mat.dvarScaling_factor_sstdep_larvae_mortality[sp][mat.aggregated_larvae_input_vectors_i[iAgg][k]][mat.aggregated_larvae_input_vectors_j[iAgg][k]];
						}else{
							for (int age=0; age<nb_lv; age++){
								Larvae_density_pred_at_obs(iAgg, k) += mat.dvarDensity[sp][0][mat.aggregated_larvae_input_vectors_i[iAgg][k]][mat.aggregated_larvae_input_vectors_j[iAgg][k]];
							}
						}
					}
					ntime_agg[iAgg] += 1;
				}
			}

		}//end of 'sp' loop

		//------------------------------------------------------//
		//		LIKELIHOOD SECTION			//
		//------------------------------------------------------//
		//I. Total abundance likelihood: augmented only if param->stock_like = true
		if (t_count == nbt_total)
			stocklike += get_stock_like(total_stock, likelihood);
	
		//II. Tag data likelihood
		if (param->tag_like[0])
			taglike += get_tag_like(likelihood, writeoutputfiles);

		//III. Fishing data likelihoods
		if ((fishing) && (t_count > nbt_building+nbstoskip) && (t_count <= nbt_no_forecast))
			get_catch_lf_like(likelihood);

		//------------------------------------------------------//
		//		COUPLAGE THON -FORAGE			//
		//------------------------------------------------------//
		if (param->flag_coupling){
			//----------------------------------------------//
			//		FORAGE RECRUITMENT		//
			//----------------------------------------------//
			for (int n=0; n< nb_forage; n++){
				mat.forage[tcur][n] += mat.mats[n];
				sumFprime[n] = sum(DtoBcell(mat.mats[n]));
			}
			
			PredationMortality(tcur, value(Total_pop));
			
			//----------------------------------------------//
			// 	TRANSPORT + PERTE DU FORAGE		//
			//----------------------------------------------//
			for (int n=0; n< nb_forage; n++){

				AverageCurrents(tcur,n);

				SolveADRE(mat.forage(tcur),n);
			}

		}// end of flag_coupling condition

		if (writeoutputfiles){
			if (!param->gcalc())	
				ConsoleOutput(1,value(likelihood));

			// let's count only those catch which will be used in the likelihood
			if (t_count <= nbt_building+nbstoskip) SUM_CATCH = 0.0;

			WriteOutput(tcur, fishing);

			if (tuna_spinup && t_count == nbt_building + 60){//temporal: the nb to skip can be set through parfile
				//cout << date_str << "Saving initial distributions for simulation skipping spinup" << endl;
				SaveDistributions(year, month);
			}
		}
		nt_dtau++;
		past_month=month;
		step_count++;
		if (qtr != past_qtr) past_qtr = qtr; 

	} // end of simulation loop

	// Compute larvae likelihood
	if (param->larvae_like[0]){
		// Compute the average larvae density over the entire period
		for (int iAgg=0; iAgg<nb_larvae_input_agg_groups; iAgg++){
			for (auto k=0u; k<mat.aggregated_larvae_input_vectors[iAgg].size(); k++){
				Larvae_density_pred_at_obs(iAgg, k) /= ntime_agg[iAgg];
			}
		}

		larvaelike += get_larvae_like(likelihood, Larvae_density_pred_at_obs);
		if (!param->gcalc()){
			cout << "Larvae likelihood: " << larvaelike << endl;
		}
	}

	if (writeoutputfiles) {SaveDistributions(year, month);
		cout << "total catch in simulation (optimization): " << SUM_CATCH << endl;
	}
	//if (param->tag_like[0]){
	//	delete_tag_releases();
	//}
	param->total_like = value(likelihood);
	double clike = value(likelihood)-lflike-taglike-stocklike-eFlike-taglike-larvaelike;
	if (!param->scalc()){ // all but sensitivity analysis
		cout << "end of forward run, likelihood: " << defaultfloat <<
		clike << " " << lflike << " " << taglike << " " << stocklike << " " << eFlike << " " << larvaelike << endl;
		if (clike+lflike && writeoutputfiles)
			OutputLikelihoodsFishery();
	}

	return value(likelihood);
}

void SeapodymCoupled::ReadLarvae()
{
	cout << "Reading input larvae file: " << endl;
	int nlevel = 0;
	string file_input;
	file_input = param->str_file_larvae;
    rw.rbin_headpar(file_input, nlon_input, nlat_input, nlevel);
	int nb_larvae_input_agg_groups = param->nb_larvae_input_agg_groups;
	if (nlevel != nb_larvae_input_agg_groups){
		cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: The number of nlevels in \"" << file_input << " does not match the number of groups from <larvae_input_aggregation_imonths> in the parameter file.\"\n";
		exit(1);

	}
	cout << file_input << endl;

	if (param->larvae_input_aggregated_flag[0]==1){
		mat.larvae_input.allocate(0,nb_larvae_input_agg_groups-1);
		for (int iAgg=0; iAgg<nb_larvae_input_agg_groups; iAgg++){
			mat.larvae_input[iAgg].allocate(1, nlon_input, 1, nlat_input);
			mat.larvae_input[iAgg].initialize();
		}

		for (int iAgg=0; iAgg<nb_larvae_input_agg_groups; iAgg++){
			ifstream litbin(file_input.c_str(), ios::binary | ios::in);
			const int sizeofDymInputType = sizeof(float);
			float buf;
			if (!litbin){
				cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_input << "\"\n";
				exit(1);
			}
			int nbytetoskip_aggregatedfile = (9 +(3* nlat * nlon) + nb_larvae_input_agg_groups + (nlat *nlon*iAgg)) * 4;
			litbin.seekg(nbytetoskip_aggregatedfile, ios::cur);
			for (int j=0;j<nlat_input;j++){
				for (int i=0;i<nlon_input;i++){
					litbin.read(( char *)&buf,sizeofDymInputType);
					mat.larvae_input[iAgg][i+1][j+1]= buf;
				}
			}
			litbin.close();
		}

		// Vector of non-NA observed density
		int ndata=0;
		for (int iAgg=0; iAgg<nb_larvae_input_agg_groups; iAgg++){
			for (int j=0;j<nlat_input;j++){
				for (int i=0;i<nlon_input;i++){
					if (map.carte[i+1][j+1]){
						if (mat.larvae_input[iAgg][i+1][j+1]>=0){
							mat.aggregated_larvae_input_vectors[iAgg].push_back(mat.larvae_input[iAgg][i+1][j+1]);
							mat.aggregated_larvae_input_vectors_i[iAgg].push_back(i+1);
							mat.aggregated_larvae_input_vectors_j[iAgg].push_back(j+1);
							ndata += 1;
						}
					}
				}
			}
		}
		cout << "Number of data points: " << ndata << endl;
	}
}
