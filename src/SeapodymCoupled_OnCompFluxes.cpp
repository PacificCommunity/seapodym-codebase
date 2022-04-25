#include "SeapodymCoupled.h"

///This is the main loop function including the calculation of biomass exchange
///between regions, based on the one time step simulations with non-zero biomass
///only in the donor region and quantification of biomass changes in all.
///See SeapodymCoupled_OnRunCoupled.cpp for the description of the main loop


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
	int migration_flag = 0;
	int step_count= 0;
	int step_fishery_count= 0;
	int jday = 0; 
	int nbstoskip = param->nbsteptoskip; // nb of time step to skip before computing likelihood
	ivector Nobs(0,nb_fishery-1); Nobs.initialize();

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
	dvariable likelihood = 0.0;
	dvariable total_stock = 0.0;
	//Reset model parameters:
	reset(x);

	//----------------------------------------------//
	// 	LOCAL MATRICES ALLOCATION SECTION       //
	//----------------------------------------------//	
	dvar_matrix Spawning_Habitat;
	dvar_matrix Total_pop;
	dvar_matrix mature_fish, immature_fish;
	dvar_matrix Habitat; 
	dvar_matrix IFR; 
	dvar_matrix ISR_denom; 
	dvar_matrix FR_pop;
	dvar_matrix Mortality; 

	Habitat.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	Mortality.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Spawning_Habitat.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Total_pop.allocate(map.imin, map.imax, map.jinf, map.jsup);
	mature_fish.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	immature_fish.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);

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
	mature_fish.initialize();
	immature_fish.initialize();

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
		InitFileFluxes();
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
	}

	//------------------------------------------------------//
	//	   COMPUTE FLUXES BETWEEN REGIONS		//
	//------------------------------------------------------//
	if (param->nb_region==0){
		cerr << "ERROR: In FLUXES mode the (MFCL) regional structure MUST be defined and activated. Revise the parfile, will stop now!" << endl;
		exit(1);
	}
	//Will consider that LF regional structure has always > 20 regions, while MFCL cannot have more than 20 regions:
	if (param->nb_region>20){
		cerr << "ERROR: The LF regional structure is ON, desactivate LF likelihood to compute fluxes, will stop now!" << endl;
		exit(1);
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
				&& (year>=(int)param->save_first_yr) && (t_count <= nbt_no_forecast) ){
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
//				ISR_denom_comp(ISR_denom, sp, tcur);
			}
			//1.3 precompute local proportions of catch by fishery
			if (fishing && param->fisheries_no_effort_exist[sp]){
				pop.Ctot_proportion_fishery_comp(map,*param,mat,rw,year,month,sp);
			}
			//1.4 precompute selectivity functions
			if (fishing)
				pop.Selectivity_comp(*param,nb_fishery,a0_adult(sp),aN_adult(sp),sp);

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
				//Compute fluxes here
				FluxesComp(mat.dvarDensity(sp,age), Habitat, Mortality, age, fishing, year, month, jday, step_fishery_count,tcur);
				
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
					//Compute fluxes here
					FluxesComp(mat.dvarDensity(sp,age), Habitat, Mortality, age, fishing, year, month, jday, step_fishery_count,tcur);
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

					if (!param->food_requirement_in_mortality(sp)){
						func.Mortality_Sp(*param, mat, map, Mortality, Habitat, sp, mean_age, age, tcur);//checked
					} else {
						Food_Requirement_Index(IFR, FR_pop, ISR_denom, sp, age, tcur, jday);
						func.Mortality_Sp(*param, mat, map, Mortality, IFR, sp, mean_age, age, tcur);//checked
					}
					//Compute fluxes here
					FluxesComp(mat.dvarDensity(sp,age), Habitat, Mortality, age, fishing, year, month, jday, step_fishery_count,tcur);

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

		}//end of 'sp' loop

		//------------------------------------------------------//
		//		LIKELIHOOD SECTION			//
		//------------------------------------------------------//
		//I. Total abundance likelihood
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

			WriteFluxes();

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


	if (writeoutputfiles) {SaveDistributions(year, month);
		cout << "total catch in simulation: " << SUM_CATCH << endl;
	}
	param->total_like = value(likelihood);
	double clike = value(likelihood)-lflike-taglike-stocklike-eFlike;
	if (!param->scalc()) // all but sensitivity analysis
		cout << "end of forward run, likelihood: " << 
		clike << " " << lflike << " " << taglike << " " << stocklike << " " << eFlike << endl;

	return value(likelihood);
}

void SeapodymCoupled::FluxesComp(dvar_matrix Density, dvar_matrix Habitat, dvar_matrix Mortality, const int age, const bool fishing, const int year, const int month, const int jday, const int step_fishery_count, const int tcur)
{
//temporal here:
int movement_fluxes_only = 1;
	const int sp = 0;
	const int nb_reg = param->nb_region_sp_B[sp];

	//dvar_matrix Density_region;
	dvar_matrix Mortality_copy;
	dvar_matrix Habitat_copy;
	//Density_region.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	Habitat_copy.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	Mortality_copy.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Habitat_copy = Habitat;
	Mortality_copy = Mortality;

	if (movement_fluxes_only) Mortality_copy.initialize();	

	if (fluxes_dt_qtr && (month==1 || month==4 || month==7 || month==10)){
		mat.fluxes_region(age).initialize();
	}

	if (nb_reg>0){
	    for (int r=0; r<nb_reg; r++){
		int reg = param->area_sp_B[sp][r]-1;

		//For MFCL need to extract quarterly metrics
		if (fluxes_dt_qtr && (month==1 || month==4 || month==7 || month==10)){	
			Density_region(reg,age).initialize();
			//only biomass in region reg is non-zero
	        	for (int i=map.imin; i <= map.imax; i++){
		     		for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
	        			if (map.carte[i][j]){
						if ((i>=map.regimin[reg]) && (i<map.regimax[reg])
						 && (j>=map.regjmin[reg]) && (j<map.regjmax[reg]))
							Density_region(reg,age,i,j) = Density(i,j);
					}
		  		}
		 	}
		}
		if (age<a0_adult[sp]){
			pop.Precalrec_juv(map, mat, Mortality_copy, tcur);
			pop.Calrec_juv(map, mat, Density_region(reg,age), Mortality_copy, tcur);
		}
		if (age>=a0_adult[sp]){
			pop.Precaldia_Caldia(map, *param, mat, Habitat_copy, Mortality_copy, sp, age, tcur,jday);
			pop.Precalrec_Calrec_adult(map,mat,*param,rw,Density_region(reg,age),Mortality_copy,tcur,fishing,age,sp,year,month,jday,step_fishery_count,0);
		}

		//integrate biomass in each region except 'reg'
		for (int rr=0; rr<nb_reg; rr++){
	       		int oth_reg = param->area_sp_B[sp][rr]-1;
			double tot_mass_region = 0.0;
			double weight = param->weight[sp][age]*0.001;
		
			for (int i=map.imin; i <= map.imax; i++){
				for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
	       				if (map.carte[i][j]){
						if ((i>=map.regimin[oth_reg]) && (i<map.regimax[oth_reg])
						 && (j>=map.regjmin[oth_reg]) && (j<map.regjmax[oth_reg]))
							tot_mass_region += value(Density_region(reg,age,i,j))*weight*param->cell_surface_area(j);		
					}
			     	}
			}
			mat.fluxes_region(age,reg,oth_reg) = tot_mass_region;
		}
	    }
	}
}

