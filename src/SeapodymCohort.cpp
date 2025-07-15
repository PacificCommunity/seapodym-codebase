#include "SeapodymCohort.h"
#include "Utilities.h"
#include "Date.h"
#include "sys/stat.h"



void SeapodymCohort::InitializeCohort(dvar_vector& x, bool writeoutputfiles) 
{
	age = 0;// 0 or older, needs to be defined by the CohortManager?

	t_count = nbt_building+1;
	mat.mats.initialize();

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
	nbt_no_forecast = t_count + nbt_spinup_tuna + nbt_total - 1;
	fishing = false;
	migration_flag = 0;
	step_count= 0;
	step_fishery_count= 0;
	jday = 0; 
	nbstoskip = param->nbsteptoskip; // nb of time step to skip before computing likelihood
	nbt_cohort = param->sp_nb_cohort_jv[0] + param->sp_nb_cohort_ad[0] - age;// simulation time for the cohort

	likelihood = 0.0;

	//----------------------------------------------//
	// 	LOCAL MATRICES ALLOCATION SECTION       //
	//----------------------------------------------//	
	Habitat.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	Mortality.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Spawning_Habitat.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Total_pop.allocate(map.imin, map.imax, map.jinf, map.jsup);
	dvarCohortDensity.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	for (int sp=0; sp<nb_species; sp++){
		////////////////////////////////////////////////////////////////////////
		// HERE dvarCohortDensity will be initialized from an external array ///
		dvarCohortDensity = mat.init_density_species(sp,age);
		////////////////////////////////////////////////////////////////////////
	}

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

	//Vector of zeroes, to avoid duplicating the F_accessibility function 
	tags_age_habitat;
	tags_age_habitat.allocate(0,aN_adult(0));
	tags_age_habitat.initialize();
	
	pop_built = 1;

	past_month=month;
	past_qtr=qtr;

	//if (!param->gcalc()){
		//need to read oxygen in case if month==past_month
		//(otherwise we may not have it for the first time steps)
		if (param->type_oxy==1 && month==past_month)
			ReadClimatologyOxy(1, month);
		//need to read oxygen in case if qtr==past_qtr 
		if (param->type_oxy==2 && qtr==past_qtr)
			ReadClimatologyOxy(1, qtr);
	//}

	//----------------------------------------------//
	// 	LIKELIHOOD INITIALISATION SECTION       //
	//----------------------------------------------//	
	//Reset model parameters:
	reset(x);

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

	//Output DYM file name
	string fileout;
	fileout = param->strdir_output + param->sp_name[0] + "_cohort.dym";//will need the date of birth stamp
	if (writeoutputfiles){
		WriteFileHeaders_submodel(fileout,true);	
		if (!param->gcalc())
			ConsoleOutput(0,0);
	}
}

void SeapodymCohort::stepForward(bool writeoutputfiles)
{
	//----------------------------------------------//
	//              INITIALISATION                  //
	//----------------------------------------------//
	sumP=0;
	sumFprime.initialize();
	sumF.initialize();
	mat.dvarCatch_est.initialize();
	//----------------------------------------------//
	//		       DATE			//
	//----------------------------------------------//
	getDate(jday);
	for (int sp=0; sp < nb_species; sp++){
		func.Seasonal_switch(*param,mat,map,jday,sp);	
	}

	//----------------------------------------------//
	//	DATA READING SECTION: U,V,T,O2,PP	//
	//----------------------------------------------//
	tcur = 1; 
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

	//------------------------------------------------------------------------------//
	//	TRANSPORT OF TUNA AGE CLASSES AND PREDICTED CATCH COMPUTATION		//
	//------------------------------------------------------------------------------//
	//------------------------------------------------------------------------------//

	for (int sp=0; sp < nb_species; sp++){
			
		int elarvae_model = param->elarvae_model[sp];
		//if !elarvae_model, elarvae_dt = 0 as elarvae_age = 0
		elarvae_dt = param->elarvae_age[sp]/deltaT; 
		double sigma_fcte_save = param->sigma_fcte;

		if (age <=param->sp_nb_cohort_jv[sp]){
			//Precompute diagonal coefficients for larvae and juvenile ADREs
			if (!elarvae_model){
				mat.u = mat.un[tcur][0]; mat.v = mat.vn[tcur][0]; 
				pop.precaldia(*param, map, mat);
				pop.caldia(map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y, mat.advection_y);
			}
		}

		//1. Precompute some variables outside of age loop
		
		//1.0 Forage scaling
		func.Forage_Scaling(*param, mat, map, sp, tcur);

		//No need to precompute accessibility for all ages, can be done at each time step for a given age			
		//1.1 Accessibility by adults (all cohorts)
		func.Faccessibility(*param, mat, map, sp, jday, tcur, pop_built, false, tags_age_habitat);//checked

		//1.6 Precompute Mortality range at age function
		func.mortality_range_age_comp(*param,mat,sp);

		//2. IMPLICIT AGE LOOP: increment age while moving through life stages 

		if (age==0){
			//2.0 If activated, the Early Larvae Model (ELM) solved over elarvae_age
			if (elarvae_model){
				//Do the time-splitting for the first age class: 
				//1- early (just a few days after yolk stage) and 
				//2- late (up to one month of age) larvae
		
				bool time_getpred = false;
				if (year>=param->larvae_like_firstyear
						&& year<=param->larvae_like_lastyear 
						&& t_count > nbt_building+nbstoskip)
					time_getpred = true;

					elarvae_model_run(Mortality,sp,tcur,time_getpred,writeoutputfiles);
			}

			//2.1 ADRE for late larvae in ELM or monthly larval class in default model 
			//2.1.1 Spawning habitat
			func.Spawning_Habitat(*param, mat, map, Spawning_Habitat, 1.0, sp, tcur, jday);
			
			//2.2 Transport and mortality of larvae (always one age class)	
			double mean_age = mean_age_cohort[sp][age]; 

			func.Mortality_Sp(*param, mat, map, Mortality, Spawning_Habitat, sp, mean_age, age, tcur);
			pop.Precalrec_juv(map, mat, Mortality, tcur, (1-elarvae_dt));//checked
			pop.Calrec_juv(map, mat, dvarCohortDensity, Mortality, tcur, (1-elarvae_dt));//checked
			
			//2.2.0 Only in the ELM mode need to reset movement rates for juveniles
			if (elarvae_model){
				param->sigma_fcte = sigma_fcte_save;
				mat.u = mat.un[tcur][0]; mat.v = mat.vn[tcur][0]; 
				//Precompute diagonal coefficients for juvenile ADREs
				pop.precaldia(*param, map, mat);
				pop.caldia(map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y, mat.advection_y);
			}
		}

		if (age >0 && age <=param->sp_nb_cohort_jv[sp]){
			//2.3. Juvenile habitat	
			if (param->cannibalism[sp]){
				Total_Pop_comp(Total_pop,sp,jday,tcur); //adjoint
				func.Juvenile_Habitat_cannibalism(*param, mat, map, Habitat, Total_pop, sp, tcur);
			} else 
				func.Juvenile_Habitat(*param, mat, map, Habitat, sp, tcur);
				
			//2.4. Transport and mortality of juvenile age classes	
			double mean_age = mean_age_cohort[sp][age];

			func.Mortality_Sp(*param, mat, map, Mortality, Habitat, sp, mean_age, age, tcur);
			pop.Precalrec_juv(map,  mat, Mortality, tcur, 1);
			pop.Calrec_juv(map, mat, dvarCohortDensity, Mortality, tcur, 1);
		}			

		if (age > param->sp_nb_cohort_jv[sp] && age <= param->sp_nb_cohort_jv[sp]+param->sp_nb_cohort_ad[sp]){
			//4. Transport and mortality of adult cohort

			//NOTE: currently current averaging doesn't depend on seasonal migrations
			if (param->vert_movement[sp]){
				func.Average_currents(*param, mat, map, age, tcur, pop_built);
			}
				
			//the option with smooth maturity to be revised and if necessary to be used later
			if (param->migrations_by_maturity_flag && param->seasonal_migrations[sp]){
				cout << "In this version no smooth maturity; enter the age at first maturity. Exit now!" << endl; exit(1);
				
			} // end of section with seasonality switch and <1 maturity at age parameter 
			else {
				migration_flag = 0;
				if (age>=param->age_mature[sp] && param->seasonal_migrations[sp]) 
				migration_flag = 1;

				func.Feeding_Habitat(*param,mat,map,Habitat,sp,age,jday,tcur,migration_flag);

				double mean_age = mean_age_cohort[sp][age];

				if (!param->food_requirement_in_mortality(sp)){
					func.Mortality_Sp(*param, mat, map, Mortality, Habitat, sp, mean_age, age, tcur);//checked
				} else {
					Food_Requirement_Index(IFR, FR_pop, ISR_denom, sp, age, tcur, jday);
					func.Mortality_Sp(*param, mat, map, Mortality, IFR, sp, mean_age, age, tcur);//checked
				}

				pop.Precaldia_Caldia(map, *param, mat, Habitat, Total_pop, sp, age, tcur, jday);//checked	

				if (!param->gcalc()){
					// Additional outputs:
					// only in simulation mode: compute mean speed 
					// in BL/sec and mean diffusion rate in nmi^2/day
					mat.MeanVarMovement(map,mat.advection_x,
						mat.advection_y,
						value(mat.dvarsDiffusion_y),
						param->MSS_species[sp],
						param->sigma_species[sp],
						param->length(sp,age),
						param->length(sp,param->sp_nb_cohorts[sp]-1),
						deltaT,sp,age);
				}
					
				mat.adult_habitat(sp,tcur,param->age_compute_habitat[sp][age]) = value(Habitat);

				pop.Precalrec_Calrec_adult(map,mat,*param,rw,
						dvarCohortDensity,Mortality,
						tcur,fishing,age,sp,year,month,
						jday,step_fishery_count,0);//checked 20150210	
			
			}
		}
	}//end of 'sp' loop
	cerr << setprecision(8) << "t_count = " << t_count << ": sum(density) = " << sum(dvarCohortDensity) << endl;

	if (writeoutputfiles){
		//Output DYM file name
		string fileout;
		fileout = param->strdir_output + param->sp_name[0] + "_cohort.dym";//will 
		if (!param->gcalc())	
			ConsoleOutput(1,value(likelihood));

		dmatrix mat2d(0, nbi - 1, 0, nbj - 1);
		mat2d.initialize();
		for (int i=map.imin; i <= map.imax; i++){
			for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
				if (map.carte[i][j]){
					mat2d(i-1,j-1) = value(dvarCohortDensity(i,j));
				}
			}
		}
		double minval = min(value(dvarCohortDensity));
		double maxval = max(value(dvarCohortDensity));

		rw.wbin_transpomat2d(fileout, mat2d, nbi-2, nbj-2, true);
		//update min-max values in header
		rw.rwbin_minmax(fileout, minval, maxval);
	}


	///////////////////////////////////////////
	//                                       //	
	// Needs some serialization of data here //
	//                                       //	
	///////////////////////////////////////////


	nt_dtau++;
	past_month=month;
	step_count++;
	if (qtr != past_qtr) past_qtr = qtr; 

	age++;
}

