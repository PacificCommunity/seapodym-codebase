#include "SeapodymCohort.h"
#include "Utilities.h"
#include "Date.h"
#include "sys/stat.h"
#include <chrono>
#include "DistDataCollector.h"
#include "DataProvider.h"

extern double time_ic_comm, time_ic_flush, time_ic_copy, time_spawning, time_getdata, time_xreset, time_init_cohort_restart, time_init_cohort_spawning;
extern long   n_ic;

void SeapodymCohort::prerun_model()
{
	OnRunFirstStep();
	time_overhead = 0;
}

std::vector<double> SeapodymCohort::GetCohortDensity()
{
	const int imin = map.imin1;
	const int imax = map.imax1;
	std::vector<double> vec;
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf1[i];
		const int jmax = map.jsup1[i];
		for (int j = jmin ; j <= jmax; j++){
			vec.push_back(dvarCohortDensity.elem_value(i, j));
		}
	}
	return vec;
}

std::vector<double> SeapodymCohort::GetInitDensity(int age)
{
	const int imin = map.imin1;
	const int imax = map.imax1;
	std::vector<double> vec;
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf1[i];
		const int jmax = map.jsup1[i];
		for (int j = jmin ; j <= jmax; j++){
			vec.push_back(mat.init_density_species(0, age, i, j));
		}
	}
	return vec;
}

void SeapodymCohort::InitializeCohort(dvar_vector& x, DistDataCollector& dataCollector, const bool writeoutputfiles) 
{

double t_all = MPI_Wtime();
double t_rs = 0.0;
double t_reset = MPI_Wtime();
	//Reset model parameters:
	reset(x);
time_xreset += MPI_Wtime() - t_reset;	

	//----------------------------------------------//
	//	ALLOCATE AND INITIALIZE COHORT DENSITY	//
	//----------------------------------------------//	
	dvarCohortDensity.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	if (cohort_id < nb_age_class){ 
t_rs = MPI_Wtime();

		//Initialize from restart file
		RestoreDistributions(mat.nb_age_built);
		dvarCohortDensity = mat.init_density_species(0,age_start);
t_rs = MPI_Wtime()-t_rs;
	} else {
		//Initialize from spawning
		int sp = 0;
		//int tcur = 0;
		//int tcur = t_count-1;
		int tcur = 1;
		
		std::vector<double> data( dataCollector.getNumSize() );

double t_block = MPI_Wtime();      
double t_flush_acc = 0.0, t_copy_acc = 0.0;

		dataCollector.startEpoch(); // should be as early as possible
		// Get density of all age class from dataCollector
		for (int aa=param->age_mature[sp]; aa<nb_age_class; aa++){
			
			int chunk_id = (tstart_cohort-1)*nb_age_class + aa;
			dataCollector.getAsync(chunk_id, data.data());

double t_f = MPI_Wtime();
			dataCollector.flush(); // now the data are ready to be used                       
t_flush_acc += MPI_Wtime() - t_f;

double t_c = MPI_Wtime();
			int index = 0;
			for (int i = map.imin1; i <= map.imax1; i++){
				const int jmin1 = map.jinf1[i];
				const int jmax1 = map.jsup1[i];
				for (int j = jmin1 ; j <= jmax1; j++){
					//stripping out derivatives runs faster, but 
					//need to be careful later to handle this in adjoint
					mat.dvarDensity(0,aa).elem_value(i,j) = data[index];
					index++;
				}
			}
t_copy_acc += MPI_Wtime() - t_c;
		}
		dataCollector.endEpoch(); // should be as late as possible

double t_total = MPI_Wtime() - t_block;
time_ic_flush += t_flush_acc;                      
time_ic_copy  += t_copy_acc;                       
time_ic_comm  += t_total - t_copy_acc;             
++n_ic;
		
double t_gd = MPI_Wtime();
		//Compute eggs at the end of t-1!	
		getDate(jday, tstart_cohort);

		// Get data from shared memory
		getData(true);
	
time_getdata += MPI_Wtime()-t_gd;

double t_sp = MPI_Wtime();
		//1. Spawning habitat (ToDo:IF NEEDED, see spawning_in_hs)
		func.Spawning_Habitat(*param, mat, map, Spawning_Habitat, 1.0, sp, tcur, jday);

		//2: Spawning biomass: 
		//If load Density(t-1(+deltaT),amature,..,na) passed by Manager to mat.dvarDensity(sp,a,i,j), then no need to change the function
		SpawningBiomass_comp(Total_pop, sp);

		//3: Reproduction
		Spawning(dvarCohortDensity,Spawning_Habitat,Total_pop,jday,sp,tcur);//checked
												    
time_spawning += MPI_Wtime() - t_sp;
	}

	if (writeoutputfiles){
		if (!param->gcalc())
			ConsoleOutput(0,0);
	}

	//----------------------------------------------//
	// PRECOMPUTE VARIABLES USED IN COHORT MODELING	//
	//----------------------------------------------//	
	//precompute thermal habitat parameters
	for (int sp=0; sp < nb_species; sp++){
		func.Vars_at_age_precomp(*param,sp);
		func.mortality_range_age_comp(*param,mat,sp);
	
		//precompute seasonal switch function
		if (param->seasonal_migrations[sp]){
			func.Seasonal_switch_year_precomp(*param,mat,map,
						value(param->dvarsSpawning_season_peak[sp]),
						value(param->dvarsSpawning_season_start[sp]),sp);
		}
	}	

	getDate(jday, t_count);// In order to get qtr
	if (param->type_oxy==1)
		getO2clm(month);
	if (param->type_oxy==2)
		getO2clm(qtr);

	age = age_start;
	past_month = month;
	past_qtr = qtr;

time_init_cohort_restart += t_rs;
time_init_cohort_spawning += MPI_Wtime() - t_all - t_rs;
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
	//		DATE and TIME-AGE		//
	//----------------------------------------------//
	//tcur = age - age_start + 1;
	//tcur = t_count;
	tcur = 1;

	int model_time_count = tstart_cohort + age - age_start;
	getDate(jday, model_time_count+1);
	for (int sp=0; sp < nb_species; sp++){
		func.Seasonal_switch(*param,mat,map,jday,sp);	
	}

	//----------------------------------------------//
	//	GET DATA FROM SHARED MEMORY		//
	//----------------------------------------------//
	getData();
	if ((param->type_oxy==1) && (month != past_month))
		getO2clm(month);
	else if ((param->type_oxy==2) && (qtr != past_qtr))
		getO2clm(qtr);

	//----------------------------------------------//
	//	COHORT DYNAMICS WITHOUT FISHING		//
	//----------------------------------------------//
	for (int sp=0; sp < nb_species; sp++){


		int elarvae_model = param->elarvae_model[sp];
		//if !elarvae_model, elarvae_dt = 0 as elarvae_age = 0
		elarvae_dt = param->elarvae_age[sp]/deltaT; 
		double sigma_fcte_save = param->sigma_fcte;

double t0 = MPI_Wtime();
		//In the serial version the code below is executed once 
		//for all larvae and juveniles. Here it is executed for 
		//each age in larval and juvenile stages
		if (age <=param->sp_nb_cohort_jv[sp]){
			//Precompute diagonal coefficients for larvae and juvenile ADREs
			if (!elarvae_model){
				mat.u = mat.un[tcur][0]; mat.v = mat.vn[tcur][0]; 
				pop.precaldia(*param, map, mat);
				pop.caldia(map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y, mat.advection_y);
			}
		}
time_overhead += (MPI_Wtime() - t0)*param->sp_nb_cohort_jv[sp]/(param->sp_nb_cohort_jv[sp]+1);

		//1.0 Forage scaling
		//As forage is age-independent, the code below is executed 
		//before the age loop in the serial version. Here it is
		//executed at each time/age step. Need to create a 3D
		//dvarForage(t,x,y) and rewrite the function to fill it once 
		//after we read all data. For this, need to call reset(x) 
		//earlier, in OnRunFirstStep.
t0 = MPI_Wtime();
		func.Forage_Scaling(*param, mat, map, sp, tcur);
time_overhead += (MPI_Wtime() - t0)*(param->sp_nb_cohorts[sp]-1)/param->sp_nb_cohorts[sp];

		//1.1 Accessibility by age class
		if (age >= param->sp_a0_adult[sp])
			func.Faccessibility_age(*param, mat, map, sp, age, jday, tcur, pop_built, false, tags_age_habitat);//checked

		if (age==0){
			//2.0 If activated, the Early Larvae Model (ELM) solved over elarvae_age
			if (elarvae_model){
				//Do the time-splitting for the first age class: 
				//1- early (just a few days after yolk stage) and 
				//2- late (up to one month of age) larvae
		
				bool time_getpred = false;
				if (year>=param->larvae_like_firstyear
						&& year<=param->larvae_like_lastyear 
						&& tstart_cohort > nbstoskip-1)
					time_getpred = true;

				elarvae_model_run(Mortality,dvarCohortDensity,sp,tcur,time_getpred,writeoutputfiles);
			}

			//2.1 ADRE for late larvae in ELM or monthly larval class in default model 
			//2.1.1 Spawning habitat
			func.Spawning_Habitat(*param, mat, map, Spawning_Habitat, 1.0, sp, tcur, jday);
			
			//2.2 Transport and mortality of larvae (always one age class)	
			double mean_age = mean_age_cohort[sp][age]; 

			func.Mortality_Sp(*param, mat, map, Mortality, Spawning_Habitat, sp, mean_age, age, tcur);
			pop.Precalrec_juv(map, mat, Mortality, tcur, (1-elarvae_dt));//checked
			pop.Calrec_juv(map, mat, dvarCohortDensity, Mortality, tcur, (1-elarvae_dt));//checked
			param->sigma_fcte = sigma_fcte_save;
		}

		if (age >0 && age < param->sp_nb_cohort_lv[sp] + param->sp_nb_cohort_jv[sp]){
			//2.2.0 Only in the ELM mode need to reset movement rates for juveniles
			param->sigma_fcte = sigma_fcte_save;
			mat.u = mat.un[tcur][0]; mat.v = mat.vn[tcur][0]; 
			//Precompute diagonal coefficients for juvenile ADREs
			pop.precaldia(*param, map, mat);
			pop.caldia(map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y, mat.advection_y);

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

		if (age >= param->sp_nb_cohort_lv[sp] + param->sp_nb_cohort_jv[sp] && age <= param->sp_nb_cohorts[sp]){
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

				func.Mortality_Sp(*param, mat, map, Mortality, Habitat, sp, mean_age, age, tcur);//checked

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

				// The following is commented-out because not necessary for now for simulation only (necessary for adjoint)
				//mat.adult_habitat(sp,tcur,param->age_compute_habitat[sp][age]) = value(Habitat);

				pop.Precalrec_Calrec_adult(map,mat,*param,rw,
						dvarCohortDensity,Mortality,
						tcur,fishing,age,sp,year,month,
						jday,step_fishery_count,0);//checked 20150210	
			
			}
		}
	}//end of 'sp' loop
	//int year, month, day, jday, xx;		
	//Date::update_time_variables(tcur, param->deltaT, param->date_mode, jday_spinup, jday, day, month, year, xx);
	//cerr << setprecision(8) << "cohort id: " << cohort_id << ", age = " << age << ", time = " << model_time_count << ", year = "<< year << ", month = " << month << ", sum(density) = " << sum(dvarCohortDensity) << endl;



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
	t_count++;
}

double SeapodymCohort::Checksum()
{
	const int imin = map.imin1;
	const int imax = map.imax1;
	double s=0;
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf1[i];
		const int jmax = map.jsup1[i];
		for (int j = jmin ; j <= jmax; j++){
			s += dvarCohortDensity.elem_value(i, j);
		}
	}
	return s;
}

void SeapodymCohort::setShmForcing(){ 

	if (!dp_) { cerr << "Error: setDataProvider() not called\n"; exit(1); }
	
	// Flatten the forcing series held in `mat` over [g0 .. nbt_total] 
	// into the shared window, in the canonical field order
	// that MUST match param->nforcings and the future loadSlice().
	
	double* W  = dp_->getDataPtr("forcing_allT");
	const int g0 = nbt_building + 1;	// fixed time origin (==1; spinup removed)
	const int nf = param->get_nforcings();  // fields per timestep
	const size_t cells = map.get_array_size();// active ragged cells per field
	const size_t slab  = (size_t)nf * cells;  // doubles per timestep

	// Partition the timestep range [g0, nbt_total] across node-local workers.
	// Each shm-rank reads its OWN timesteps (ReadTimeSeriesData seeks by t_series,
	// random access) into its PRIVATE mat[g0], then writes its DISJOINT window
	// slabs. No write conflicts: private scratch, disjoint offsets. Caller must
	// MPI_Barrier(dp_->getShmComm()) after this returns.
	int shmRank = dp_->getShmRank();
	int shmSize = 1;
	MPI_Comm_size(dp_->getShmComm(), &shmSize);
	const int Tsteps = nbt_total - g0 + 1;                       // number of timesteps
	const int per    = (Tsteps + shmSize - 1) / shmSize;        // ceil split
	const int tlo    = g0 + shmRank * per;
	const int thi    = (tlo + per < nbt_total + 1) ? (tlo + per) : (nbt_total + 1);

	for (int t = tlo; t < thi; ++t) {
		double* base = W + (size_t)(t - g0) * slab;  // this timestep's block
		size_t  f = 0; // running field index

		auto put = [&](const dmatrix& src) {
			double* dst = base + f * cells;
			size_t  c   = 0;
			for (int i = map.imin; i <= map.imax; ++i)
				for (int j = map.jinf[i]; j <= map.jsup[i]; ++j)
					dst[c++] = src(i, j);
			++f;
		};

		//----------------------------------------------//
		//	DATA READING SECTION: U,V,T,O2,PP	//
		//----------------------------------------------//
		//TIME SERIES 
		t_series = t - nbt_building + nbt_start_series;
		ReadTimeSeriesData(g0,t_series);	

		put(mat.np1[g0]);
		for (int n = 0; n < nb_forage; ++n) put(mat.forage[g0][n]);
		if (param->use_sst) put(mat.sst[g0]);
		for (int k = 0; k < nb_layer; ++k) put(mat.tempn[g0][k]);
		for (int k = 0; k < nb_layer; ++k) put(mat.un[g0][k]);
		for (int k = 0; k < nb_layer; ++k) put(mat.vn[g0][k]);
		if (!param->type_oxy)
			for (int k = 0; k < nb_layer; ++k) put(mat.oxygen[g0][k]);
		if (param->use_vld) put(mat.vld[g0]);
		if (param->use_ph1) put(mat.ph1[g0]);
	}

	// Set O2 from climatology
	if (param->type_oxy && dp_->isShmRoot()){
		double* W_O2clm  = dp_->getDataPtr("forcing_O2clm");
		const int nf_O2clm = param->get_nforcings_O2clm();
		const size_t slab_O2clm  = (size_t)nf_O2clm * cells;  // doubles per timestep
		int nbt_O2clm = 12;
		if (param->type_oxy==2)
			nbt_O2clm = 4;

		for (int t_clm = 1; t_clm <= nbt_O2clm; ++t_clm) {
			size_t  f = 0; // running field index
			double* base = W_O2clm + (size_t)(t_clm-1) * slab_O2clm;  // this timestep's block

			auto put = [&](const dmatrix& src) {
				double* dst = base + f * cells;
				size_t  c   = 0;
				for (int i = map.imin; i <= map.imax; ++i)
					for (int j = map.jinf[i]; j <= map.jsup[i]; ++j)
						dst[c++] = src(i, j);
				++f;
			};

			ReadClimatologyOxy(1, t_clm);
			for (int k = 0; k < nb_layer; ++k) put(mat.oxygen[1][k]);
		}
	}
}

void SeapodymCohort::getData(bool spawning_habitat_only){
	if (!dp_) { cerr << "Error: setDataProvider() not called\n"; exit(1); }

	double* W  = dp_->getDataPtr("forcing_allT");
	const int g0 = nbt_building + 1;	// fixed time origin (==1; spinup removed)
	const int nf = param->get_nforcings();  // fields per timestep
	const size_t cells = map.get_array_size();// active ragged cells per field
	const size_t slab  = (size_t)nf * cells;  // doubles per timestep
	double* base = W + (size_t)(t_count - g0) * slab;
	size_t  f = 0; // running field index

	auto get = [&](dmatrix& dst) {
		double* src = base + f * cells;
		size_t  c   = 0;
		for (int i = map.imin; i <= map.imax; ++i)
			for (int j = map.jinf[i]; j <= map.jsup[i]; ++j)
				dst(i, j) = src[c++];
		++f;
	};

	get(mat.np1[g0]);
	for (int n = 0; n < nb_forage; ++n){
		if (spawning_habitat_only && param->day_layer[n] && param->night_layer[n])
			++f;
		else
			get(mat.forage[g0][n]);
	};
	if (param->use_sst) get(mat.sst[g0]);
	if (!spawning_habitat_only){
		for (int k = 0; k < nb_layer; ++k) get(mat.tempn[g0][k]);
		for (int k = 0; k < nb_layer; ++k) get(mat.un[g0][k]);
		for (int k = 0; k < nb_layer; ++k) get(mat.vn[g0][k]);
		if (!param->type_oxy)
			for (int k = 0; k < nb_layer; ++k) get(mat.oxygen[g0][k]);
		if (param->use_vld) get(mat.vld[g0]);
		if (param->use_ph1) get(mat.ph1[g0]);
	}
}

void SeapodymCohort::getO2clm(int t_clm){
	double* W  = dp_->getDataPtr("forcing_O2clm");
	const int g0 = 1;
	const int nf = param->get_nforcings_O2clm();
	const size_t cells = map.get_array_size();// active ragged cells per field
	const size_t slab  = (size_t)nf * cells; 
	double* base = W + (size_t)(t_clm-g0) * slab;
	size_t  f = 0; // running field index

	auto get = [&](dmatrix& dst) {
		double* src = base + f * cells;
		size_t  c   = 0;
		for (int i = map.imin; i <= map.imax; ++i)
			for (int j = map.jinf[i]; j <= map.jsup[i]; ++j)
				dst(i, j) = src[c++];
		++f;
	};

	for (int k = 0; k < nb_layer; ++k) get(mat.oxygen[g0][k]);
}
