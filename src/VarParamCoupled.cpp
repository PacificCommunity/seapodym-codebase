#include "VarParamCoupled.h"
#include "Utilities.h"
//#include "Dimensions.h"
#include "Date.h"
#include "sys/stat.h"

///This routine will read parfile in xml format.
///See User Manual for more details about model parameters

//some very local variables
dvector statpars_temp;
adstring_array statpar_names_temp;

bool VarParamCoupled::read(const string& parfile) 
{

	try {
		doc.read(parfile.c_str());
	} catch (const runtime_error& e) {
		cerr << e.what() << '\n';
		return false;
	}

	////////////////////////
	// SIMULATION PARAMETERS
	////////////////////////
        latitudeMin = doc.getDouble("/latitudeMin", "value");
        latitudeMax = doc.getDouble("/latitudeMax", "value");
        longitudeMin = doc.getDouble("/longitudeMin", "value");
        longitudeMax = doc.getDouble("/longitudeMax", "value");

        nb_species = doc.getInteger("/nb_species", "value");
        nb_layer = doc.getInteger("/nb_layer", "value");
        deltaX = doc.getInteger("/deltaX", "value");
        deltaY = doc.getInteger("/deltaY", "value");
        deltaT = doc.getInteger("/deltaT", "value");

//IMPORTANT: in order to get correct dimentions (nbi,nbj)=(nlon+2,nlat+2) 
//the coordinates (latitudeMax,longitudeMin) in parfile
//should give the NORTH-WEST corner of the upper-left grid cell 
//and (latitudeMin,longitueMax) should give the SOUTH-EAST corner of the
//lower-right grid cell 
	nbi = int(((longitudeMax - longitudeMin)*60.0/deltaX)) +2;
	nbj = int(((latitudeMax  - latitudeMin)*60.0/deltaY)) +2;
	//nbi = int(((longitudeMax - longitudeMin)*60.0/deltaX)) +2;
	//nbj = int(((latitudeMax  - latitudeMin)*60.0/deltaY)) +2;
	//maxn= Utilities::MyMax( nbi, nbj );

//Dimensions *dim;
//int x = dim->getInstance().get_nbi();
//cout << x << endl;

        iterationNumber = doc.getInteger("/iterationNumber", "value");
        tuna_spinup = doc.getInteger("/tuna_spinup", "value");
        if (!doc.get("/save_first_date","year").empty()){
		int year = doc.getInteger("/save_first_date","year");
		int month = doc.getInteger("/save_first_date","month");
		int day = 15;
		save_first_yr = fdate(year,month);
		if (!doc.get("/save_first_date","day").empty())
			day = doc.getInteger("/save_first_date","day");
		ndatini = year*10000+month*100+day;

		year = doc.getInteger("/save_last_date","year");
		month = doc.getInteger("/save_last_date","month");
		save_last_yr = fdate(year,month);
		if (!doc.get("/save_last_date","day").empty())
			day = doc.getInteger("/save_last_date","day");
		ndatfin = year*10000+month*100+day;

	} else {
		//old parfiles:
		save_first_yr = doc.getDouble("/save_first_yr", "value");
        	save_last_yr = doc.getDouble("/save_last_yr", "value");
		ndatini = get_year(save_first_yr)*10000+get_month(save_first_yr)*100+15;
		ndatfin = get_year(save_last_yr)*10000+get_month(save_last_yr)*100+15;
	}
        nb_yr_forecast = doc.getInteger("/nb_yr_forecast", "value");
	nbsteptoskip   = doc.getInteger("/nb_step_to_skip", "value");

	//update ndatfin using nb_yr_forecast value:
	ndatfin += nb_yr_forecast*10000;

	date_mode = 3; //default value for monthly stepping, uses 360-days calendar
	if (deltaT<30) date_mode = 1;
	//rundates will be used for reading the data (fishing etc) for the simulation time steps
	int nbstot_max = ((int)((ndatfin-ndatini)/1e4)+1)*366/deltaT;
	rundates.allocate(0,nbstot_max);

	//The date of first recruitment (by default is the first date of run)
	first_recruitment_date = ndatini;
	if (!doc.get("/first_recruitment_date","year").empty()){
		int year = doc.getInteger("/first_recruitment_date","year");
		int month = doc.getInteger("/first_recruitment_date","month");
		int day = doc.getInteger("/first_recruitment_date","day");
		first_recruitment_date = year*10000+month*100+day;
	}

	//////////////////////////////
	// INPUT FILES and DIRECTORIES
	//////////////////////////////
	str_dir = doc.get("/strdir", "value");
        str_file_mask = str_dir + doc.get("/str_file_mask", "value");
        str_file_topo = str_dir + doc.get("/str_file_topo", "value");
	str_dir_forage = str_dir;
	if (!doc.get("/strdir_forage", "value").empty())
		str_dir_forage += doc.get("/strdir_forage", "value");
	str_dir_init = str_dir;
	if (!doc.get("/strdir_init", "value").empty())
		str_dir_init += doc.get("/strdir_init", "value");
	str_dir_fisheries = str_dir;
	if (!doc.get("/strdir_fisheries", "value").empty())
		str_dir_fisheries = doc.get("/strdir_fisheries", "value");
		//str_dir_fisheries += doc.get("/strdir_fisheries", "value");
	str_dir_tags = str_dir;
	if (!doc.get("/strdir_tags", "value").empty())
		str_dir_tags = doc.get("/strdir_tags", "value");

	strfile_pp = str_dir + doc.get("/strfile_pp", "value");

	use_sst = 0;
	if (!doc.get("/strfile_sst","value").empty()){
		strfile_sst = str_dir + doc.get("/strfile_sst", "value");
		use_sst = 1;
		cout << "SST is in use in this simulation!" << endl; 
	}
	use_vld = 0;
	if (!doc.get("/strfile_vld","value").empty()){
		strfile_vld = str_dir + doc.get("/strfile_vld", "value");
		use_vld = 1;
	}
	use_ph1 = 0;
	if (!doc.get("/strfile_ph1","value").empty()){
		strfile_ph1 = str_dir + doc.get("/strfile_ph1", "value");
		use_ph1 = 1;
	}

	type_oxy = 0;
	if (!doc.get("/type_oxy","value").empty())
		type_oxy = doc.getInteger("/type_oxy", "value");
	else cout << "WARNING: the type of oxygen data is not set up, default is monthly time series!" << endl; 
	

	for (int k=0;k<nb_layer;k++){
		std::ostringstream ostr;
		ostr << "layer" << k;
		strfile_t.push_back(str_dir+doc.get("/strfile_t", ostr.str()));
		strfile_u.push_back(str_dir+doc.get("/strfile_u", ostr.str()));
		strfile_v.push_back(str_dir+doc.get("/strfile_v", ostr.str()));
		strfile_oxy.push_back(str_dir+doc.get("/strfile_oxy", ostr.str()));
	}
	strfile_ppmc = str_dir+doc.get("/strfile_pp_clm", "value");
	for (int k=0;k<nb_layer;k++){
		std::ostringstream ostr;
		ostr << "layer" << k;
		strfile_tmc.push_back(str_dir+doc.get("/strfile_t_clm", ostr.str()));
		strfile_umc.push_back(str_dir+doc.get("/strfile_u_clm", ostr.str()));
		strfile_vmc.push_back(str_dir+doc.get("/strfile_v_clm", ostr.str()));
		if (!type_oxy)
			strfile_oxymc.push_back(str_dir+doc.get("/strfile_oxy_clm", ostr.str()));
	}

	///////////////////////////////
	// OUTPUT FILES and DIRECTORIES
	///////////////////////////////
	strdir_output = doc.get("/strdir_output", "value");
	//strdir_output = str_dir + doc.get("/strdir_output", "value");
	wbin_flag = 0;
	if (!doc.get("/wbin_flag","value").empty())
		wbin_flag = doc.getInteger("/wbin_flag", "value");

        write_all_cohorts_dym = 0;
        if (!doc.get("/write_all_cohorts_dym","value").empty())
                write_all_cohorts_dym = doc.getInteger("/write_all_cohorts_dym", "value");

	write_all_fisheries_dym = 0;
	if (!doc.get("/write_all_fisheries_dym","value").empty())
		write_all_fisheries_dym = doc.getInteger("/write_all_fisheries_dym", "value");

        //If output directory doesn't exist create it
	string test = strdir_output+"/test";
        ofstream ecritbin(test.c_str(), ios::binary|ios::out);
        if (!ecritbin){
                cerr << "WARNING : Cannot find output directory '" << strdir_output << "', will create the new folder './output'" << endl;
                strdir_output = "./output/";
                mkdir(strdir_output.c_str(),0777);
        } else remove(test.c_str());

/*
	//MPA simulations
	mpa_simulation = 0;
	mpa_scenario = 0;
	if (!doc.get("/mpa_simulation").empty()){
		mpa_simulation = doc.getInteger("/mpa_simulation/isON","value");
		if (mpa_simulation){
			mpa_scenario = doc.getInteger("/mpa_simulation/scenario", "value");
			str_file_maskMPA = str_dir + doc.get("/mpa_simulation/mpa_mask", "value");		 
		}
	}		
*/
        mpa_simulation = 0;
        nb_mpa = 0;
        if (!doc.get("/mpa_simulation").empty()){
                mpa_simulation = doc.getInteger("/mpa_simulation/isON","value");
                if (mpa_simulation)
                        nb_mpa = doc.getInteger("/mpa_simulation/nb_mpa", "value");
        }

	//Habitat parameter estimation
	habitat_run_type = 0; //0 -spawning; 1- feeding (in this case will read age in months)
	nb_habitat_run_age = 1;//one for spawning habitat, can be several for feeding habitat
			       //Note, for optimization it is indespensable to provide at least three!
	if (!doc.get("/habitat_run","type").empty()){
		habitat_run_type = doc.getInteger("/habitat_run","type");
		if (habitat_run_type>0){
			nb_habitat_run_age = doc.getInteger("/habitat_run","nb_ages");
		}
		habitat_run_age.allocate(0,nb_habitat_run_age-1);
		habitat_run_age.initialize();
		if (habitat_run_type>0){
			for (int n=0; n<nb_habitat_run_age; n++){
				habitat_run_age[n] = doc.getInteger("/habitat_run_ages",n);
			}
		}
		cout << endl << "Habitat type: " << habitat_run_type << "; ages: " << habitat_run_age << endl << endl; 		
	}

	////////////////////
	// FORAGE SECTION
	// int Tr - Time (in days) before recruitment in the forage population
	// m_inv_Lambda - inverse of Forage mortality through time transfer
	// E - Ecological transfer coefficent from new primary prod to forage
	// double C Conversion parameter from unit of PP (nitrogen or carbon) to wet weight of forage
	//(Iverson 1990): ratio C/N for fish =3.6 ; ratio Dry weight/ C =2.4 ; ratio Wet weight/Dry weight = 3.3
	// 1 mmol N -> g N (1 mole N =14g) 1 mmol N = 14/1000* 3.6 * 2.4 * 3.3 = 0.4 g WW
	// 1 mmol C -> g C (1 mole C =12g) 1 mmol C = 12/1000* 2.4 * 3.3 = 0.095 g WW
	// source_frg[nb_forage]  - percentage of new prod transferred to each forage
	// day_layer[nb_forage]   - layer where each forage is during the day
	// night_layer[nb_forage] - layer where each forage is during the night
	// tstep_forage	   - time step for the forage (days)
	// sigma_fcte - diffusion coefficient for forage (nm2.d-1)
	////////////////////

        Tr_max = doc.getInteger("/Tr_max", "value");
        Tr_exp = doc.getDouble("/Tr_exp", "value");
        inv_lambda_max = doc.getInteger("/inv_lambda_max", "value");
        inv_lambda_curv = doc.getDouble("/inv_lambda_curv", "value");
        E = doc.getDouble("/E", "value");
        c_pp = doc.getDouble("/c_pp", "value");
        nb_forage = doc.getInteger("/nb_forage", "value");
	for (int n=0;n<nb_forage;n++) {
        	frg_name.push_back(doc.get("/frg_name", n));
	}
	source_frg.allocate(0, nb_forage - 1);
	for (int n=0;n<nb_forage;n++) {
		source_frg[n] = doc.getDouble("/source_frg", frg_name[n]);
	}
	day_layer.allocate(0, nb_forage - 1);
	for (int n=0;n<nb_forage;n++) {
	        day_layer[n] = doc.getInteger("/day_layer", frg_name[n]);
	}
	night_layer.allocate(0, nb_forage - 1);
	for (int n=0;n<nb_forage;n++) {
        	night_layer[n] = doc.getInteger("/night_layer", frg_name[n]);
	}

        //tstep_forage = doc.getInteger("/tstep_forage", "value");
        sigma_fcte = doc.getDouble("/sigma_fcte", "value");
	sigma_fcte = sigma_fcte * deltaT;
	pp_transform = deltaT * c_pp * E;

	build_forage = false;
	if (!doc.get("/build_forage","flag").empty()){
		int flag = doc.getInteger("/build_forage", "flag");
		if (flag) build_forage = true;
	}
	
	//////////////////
	// SPECIES SECTION
	//////////////////
		
	for (int sp=0;sp<nb_species;sp++) {
	        sp_name.push_back(doc.get("/sp_name", sp));
	}
	if (nb_species) {
		///create vectors of model parameters
		///sp_unit_age_class_jv.allocate(0, nb_species - 1);
		///sp_nb_age_class_jv.allocate(0, nb_species - 1);
		///juv_length.allocate(0, nb_species - 1);
		///juv_weight.allocate(0, nb_species - 1);
		///sp_nb_age_class_ad.allocate(0, nb_species - 1);
		///sp_unit_age_class_ad.allocate(0, nb_species - 1);
		///sp_unit_age_class.allocate(0, nb_species - 1);
		sp_nb_cohorts.allocate(0, nb_species - 1);		//NEW
		sp_nb_cohort_lv.allocate(0, nb_species - 1);		//NEW
		sp_nb_cohort_jv.allocate(0, nb_species - 1);		//renamed sp_nb_age_class_jv;
		sp_nb_cohort_ad.allocate(0, nb_species - 1);		//renamed sp_nb_age_class_ad;
		sp_a0_adult.allocate(0, nb_species - 1);		//NEW;
		sp_unit_cohort.allocate(0, nb_species - 1);		//renamed sp_unit_age_class
		length.allocate(0, nb_species - 1);
		length_bins.allocate(0, nb_species - 1);
		weight.allocate(0, nb_species - 1);
		age_mature.allocate(0, nb_species - 1);
		maturity_age.allocate(0, nb_species - 1);
		age_autonomous.allocate(0, nb_species - 1);
		age_recruit.allocate(0, nb_species - 1);
		age_compute_habitat.allocate(0, nb_species - 1);
		seasonal_migrations.allocate(0, nb_species - 1);
		spawning_season_peak.allocate(0, nb_species - 1);
		spawning_season_start.allocate(0, nb_species - 1);
		cannibalism.allocate(0, nb_species - 1);
		nb_recruitment.allocate(0, nb_species - 1);
		a_adults_spawning.allocate(0, nb_species - 1);
		alpha_hsp_prey.allocate(0, nb_species - 1);
		alpha_hsp_predator.allocate(0, nb_species - 1);
		beta_hsp_predator.allocate(0, nb_species - 1);
		Mp_mean_max.allocate(0, nb_species - 1);
		Mp_mean_exp.allocate(0, nb_species - 1);
		Ms_mean_max.allocate(0, nb_species - 1);
		Ms_mean_slope.allocate(0, nb_species - 1);
		M_mean_range.allocate(0, nb_species - 1);
		food_requirement_in_mortality.allocate(0,nb_species-1);
		residual_competition.allocate(0,nb_species-1);
		uncouple_sst_larvae.allocate(0,nb_species-1);
		a_sst_spawning.allocate(0, nb_species - 1);
		b_sst_spawning.allocate(0, nb_species - 1);
		a_sst_larvae.allocate(0, nb_species - 1);
		b_sst_larvae.allocate(0, nb_species - 1);
		gaussian_thermal_function.allocate(0,nb_species-1);
		a_sst_habitat.allocate(0, nb_species - 1);
		b_sst_habitat.allocate(0, nb_species - 1);
		T_age_size_slope.allocate(0, nb_species - 1);
		thermal_func_delta.allocate(0,nb_forage-1);
		for (int n=0; n<nb_forage; n++){
			thermal_func_delta[n].allocate(0, nb_species - 1);
		}
		a_oxy_habitat.allocate(0, nb_species - 1);
		b_oxy_habitat.allocate(0, nb_species - 1);
		eF_habitat.allocate(0,nb_forage-1);
		for (int n=0; n<nb_forage; n++){
			eF_habitat[n].allocate(0, nb_species - 1);
		}
		hp_cannibalism.allocate(0, nb_species - 1);
		forage_ration.allocate(0, nb_species - 1);
		sigma_species.allocate(0, nb_species - 1);
		MSS_species.allocate(0, nb_species - 1);
		MSS_size_slope.allocate(0, nb_species - 1);
		c_diff_fish.allocate(0, nb_species - 1);

		sigma_ha.allocate(0, nb_species - 1);
		temp_age.allocate(0, nb_species - 1);
	}
	////////////////
	//  HABITATS  //   
	////////////////

	// 1. Spawning habitat
	for (int sp=0;sp<nb_species;sp++){
		//flag that determines whether the seasonality will be taken
		//into account in migrations of adult population
		seasonal_migrations[sp] = 0;
		if (!doc.get("/seasonal_migrations",sp_name[sp]).empty())
			seasonal_migrations[sp] = doc.getInteger("/seasonal_migrations", sp_name[sp]);
		if (seasonal_migrations[sp])
			if (gcalc() && latitudeMax>65) cout<< "ATTN: Make sure the latitudes above 65N/S are masked " << 
				       				 "to avoid NANs in the gradient of Season_peak"<< endl;	
	

		//old parameter files:
		if (doc.get("/spawning_season_peak").empty()){
			cerr << "Current version require variable parameters <spawning_season_peak> and <spawning_season_start>! Modify your parameter file accordingly. The program will exit now... " << endl; exit(1);
		}
		spawning_season_peak[sp]  = doc.getDouble("/spawning_season_peak", sp_name[sp]);
		spawning_season_start[sp] = doc.getDouble("/spawning_season_start", sp_name[sp]);

		//flag to uncouple the larvae habitat and adult habitat temperatures
		uncouple_sst_larvae.initialize(); //default value
		b_sst_larvae.initialize();
		if (!doc.get("/uncouple_sst_larvae",sp_name[sp]).empty()){
			uncouple_sst_larvae[sp] = doc.getInteger("/uncouple_sst_larvae", sp_name[sp]);
			if (uncouple_sst_larvae[sp]){
				//optimal temperature in Gaussian temperature function for larvae survival
		               	a_sst_larvae[sp] = doc.getDouble("/a_sst_larvae", sp_name[sp]);
		               	b_sst_larvae[sp] = doc.getDouble("/b_sst_larvae", sp_name[sp]);
			}
		}
		if (!uncouple_sst_larvae[sp]){
			doc.set("/a_sst_larvae/variable","use","false");		
			doc.set("/b_sst_larvae/variable","use","false");		
		}

		//standard deviation in Gaussian temperature function for spawning
               	a_sst_spawning[sp] = doc.getDouble("/a_sst_spawning", sp_name[sp]);
	
		//optimal SST for larvae and juveniles
                b_sst_spawning[sp] = doc.getDouble("/b_sst_spawning", sp_name[sp]);
	
		//previously alpha_spawning in food to predator ratio in the spawning index
		//became two parameters - one in the food (prey) function and another in predator one 
               	alpha_hsp_prey[sp] = doc.getDouble("/alpha_hsp_prey", sp_name[sp]);
               	alpha_hsp_predator[sp] = doc.getDouble("/alpha_hsp_predator", sp_name[sp]);
               	beta_hsp_predator[sp] = doc.getDouble("/beta_hsp_predator", sp_name[sp]);
	
		//maximal number of larvae (by cell) at large spawning biomass of adults
               	nb_recruitment[sp] = doc.getDouble("/nb_recruitment", sp_name[sp]);
		
		//slope coefficient in Beverton-Holt function
               	a_adults_spawning[sp] = doc.getDouble("/a_adults_spawning", sp_name[sp]);
	}

	// 2. Juvenile habitat
	for (int sp=0;sp<nb_species;sp++){	
		cannibalism[sp] = 0;
		hp_cannibalism[sp] = 0;
		if (!doc.get("/cannibalism_juv",sp_name[sp]).empty())
			cannibalism[sp] = doc.getInteger("/cannibalism_juv",sp_name[sp]);
		if (!cannibalism[sp])
			doc.set("/hp_cannibalism/variable","use","false");

		// slope coefficient of the habitat cannibalism function
	     	hp_cannibalism[sp] = doc.getDouble("/hp_cannibalism", sp_name[sp]);
		
	}

	// 3. Feeding habitat
	for (int sp=0;sp<nb_species;sp++){

		gaussian_thermal_function[sp] = 1; //default value
		if (!doc.get("/gaussian_thermal_function",sp_name[sp]).empty()){
			gaussian_thermal_function[sp] = doc.getInteger("/gaussian_thermal_function", sp_name[sp]);
		}
			
		//standard deviation in Gaussian temperature function
               	a_sst_habitat[sp] = doc.getDouble("/a_sst_habitat", sp_name[sp]);

		//optimal temperature for oldest cohort
               	b_sst_habitat[sp] = doc.getDouble("/b_sst_habitat", sp_name[sp]);
	
		//slope of the function of optimal temperature by age as a function of fish size  
               	T_age_size_slope[sp] = doc.getDouble("/T_age_size_slope", sp_name[sp]);

		if (!gaussian_thermal_function[sp]){
			doc.set("/a_sst_habitat/variable","use","false");		
			//three parameters of alternative (non-Gaussian) thermal function in feeding habitat index
			for (int n=0; n<3; n++){
				std::ostringstream ostr;
				ostr << "delta" << n+1;
		               	thermal_func_delta[n][sp] = doc.getDouble(string("/thermal_func/") + ostr.str(),sp_name[sp]);
			}
		}
		if (gaussian_thermal_function[sp]){
			for (int n=0; n<3; n++){
				std::ostringstream ostr;
				ostr << "delta" << n+1;
				string s = "/thermal_func/"+ostr.str()+"/variable";
				doc.set(s,"use","false");
			}
		}		
		//curvature coefficient of the habitat 2ml O2 function
               	a_oxy_habitat[sp] = doc.getDouble("/a_oxy_habitat", sp_name[sp]);
	
		//threshold coefficent of the habitat 2ml O2 function
               	b_oxy_habitat[sp] = doc.getDouble("/b_oxy_habitat", sp_name[sp]);

		//six eF parameters for forage groups in feeding habitat function
		for (int n=0; n<nb_forage; n++){
	               	eF_habitat[n][sp] = doc.getDouble(string("/eF_habitat/") + frg_name[n], sp_name[sp]);
		}
	}

	///////////////
	// MOVEMENT  //   
	/////////////// 

	vert_movement.allocate(0,nb_species-1);
	for (int sp=0;sp<nb_species;sp++){

		//flag to aggregate currents through vertical layers accessible to fish
		vert_movement[sp] = 0;
		if (!doc.get("/vertical_movement",sp_name[sp]).empty())
			vert_movement[sp] = doc.getInteger("/vertical_movement", sp_name[sp]);

		//multiplier to maximal diffusion coefficient which is linked to fish size
               	sigma_species[sp] = doc.getDouble("/sigma_species", sp_name[sp]);

		//curvature coefficient in the function of habitat to compute local diffusion rate
                c_diff_fish[sp] = doc.getDouble("/c_diff_fish", sp_name[sp]);

		//maximal sustainable speed in the units 'body length'/sec
                MSS_species[sp] = doc.getDouble("/MSS_species", sp_name[sp]);

		//scaling exponent in power low giving the species sustainable speed
		//in m/sec, i.e. V_mss = MSS_species * pow(L,MSS_size_slope)
                MSS_size_slope[sp] = doc.getDouble("/MSS_size_slope", sp_name[sp]);
	}

	///////////////////////////////////////////
	// Fixed parameters of population structure
	///////////////////////////////////////////
	sp_nb_cohort_lv.initialize();
	sp_nb_cohort_jv.initialize();
	sp_nb_cohort_ad.initialize();
	for (int sp=0;sp<nb_species;sp++) {

		int nb_life_stages = doc.getInteger("/nb_life_stages",sp_name[sp]);
		//life_stage = Utilities::create1d(life_stage, nb_life_stages);
		for (int n=0; n<nb_life_stages; n++)
			life_stage.push_back(doc.get("/life_stage/"+sp_name[sp],n));

		ivector nb_cohort_life_stage;
		nb_cohort_life_stage.allocate(0,nb_life_stages-1);
		for (int n=0; n<nb_life_stages; n++){
			nb_cohort_life_stage[n] = doc.getInteger("/nb_cohort_life_stage/"+sp_name[sp],n);
			if (life_stage[n].find("larvae")==0) sp_nb_cohort_lv[sp] = nb_cohort_life_stage[n];
			if (life_stage[n].find("juvenile")==0) sp_nb_cohort_jv[sp] = nb_cohort_life_stage[n];
			if (life_stage[n].find("young")==0) sp_nb_cohort_ad[sp] = nb_cohort_life_stage[n];
			if (life_stage[n].find("adult")==0) sp_nb_cohort_ad[sp] += nb_cohort_life_stage[n];
		}
		sp_nb_cohorts[sp] = sum(nb_cohort_life_stage);
		sp_a0_adult[sp] = sp_nb_cohort_lv[sp] + sp_nb_cohort_jv[sp];

		sp_unit_cohort[sp].allocate(0,sp_nb_cohorts[sp]-1);

		sigma_ha[sp].allocate(sp_a0_adult[sp],sp_nb_cohorts[sp]-1);
		temp_age[sp].allocate(sp_a0_adult[sp],sp_nb_cohorts[sp]-1);
		sigma_ha.initialize();
		temp_age.initialize();

		for (int a=0; a<sp_nb_cohorts[sp]; a++){
                	sp_unit_cohort[sp][a] = doc.getInteger("/sp_unit_cohort/"+ sp_name[sp], a);
			if (sp_unit_cohort[sp][a]<deltaT){
				cerr<<"Sorry, have to stop here, the model requires that sp_unit_cohort (size of the cohort) >= deltaT!"<< endl;
				exit(1);
			}	
			sp_unit_cohort[sp][a] = sp_unit_cohort[sp][a]/deltaT;
		}

		//length and weigth (all cohorts)
		length[sp].allocate(0, sp_nb_cohorts[sp] - 1);
		weight[sp].allocate(0, sp_nb_cohorts[sp] - 1);
		//length_bins[sp].allocate(sp_a0_adult[sp], sp_nb_cohorts[sp]);
		length_bins[sp].allocate(0, sp_nb_cohorts[sp]);
		length_bins[sp].initialize();

		for (int a = 0; a < sp_nb_cohorts[sp]; a++) {
			length[sp][a] = doc.getDouble(string("/length/") + sp_name[sp], a);
			weight[sp][a] = doc.getDouble(string("/weight/") + sp_name[sp], a);
		}
		// these length bins bounds will be used to aggregate LF data within Seapodym cohorts. 
		if (!doc.get("/length/bounds").empty()){
			for (int a = sp_a0_adult[sp]; a <= sp_nb_cohort_ad[sp]; a++) 
				length_bins[sp][a] = doc.getDouble(string("/length/") + "bounds", a);
			
		} else {  //approximate bounds linearly if they are not in the parfile
			double L_pr = length(sp,sp_a0_adult[sp]-1);
			for (int a=sp_a0_adult[sp]; a<sp_nb_cohorts[sp]-1; a++){
				length_bins[sp][a]  = 0.5*(L_pr+length(sp,a));
				L_pr = length(sp,a);
			} 
			length_bins[sp][sp_nb_cohorts[sp]-1] = 1.5*L_pr-.5*length(sp,sp_nb_cohorts[sp]-3);
			length_bins[sp][sp_nb_cohorts[sp]] = length[sp][sp_nb_cohorts[sp]-1]+2.0;
		}
//cout << sp_a0_adult << " " << length(sp,sp_a0_adult[sp]-1) << endl;
//cout << length_bins << endl; exit(1);
		//Other important life-span parameters:
		//age for larvae-juvenile phase with passive transport
		age_autonomous[sp] = 0;
	        //age_autonomous[sp] = doc.getInteger("/age_autonomous", sp_name[sp]);

		//age at recruitment
        	age_recruit[sp] = doc.getInteger("/age_recruit", sp_name[sp]);
		//ages to compute adult habitat index
		age_compute_habitat[sp].allocate(0,sp_nb_cohorts[sp]-1);
		age_compute_habitat[sp].initialize();
		if (doc.get("/age_compute_habitat/" + sp_name[sp]).empty()){                   
			for (int a=1; a<sp_nb_cohorts[sp]; a++)
				age_compute_habitat[sp][a] = age_compute_habitat[sp][a-1]+1;
		} else if (!doc.get("/age_compute_habitat/" + sp_name[sp]).empty()){
			for (int a=0; a<sp_nb_cohorts[sp]; a++)
				age_compute_habitat[sp][a] = doc.getInteger("/age_compute_habitat/"+ sp_name[sp], a);
		}
		//Check if the habitat indices are valid:
		int ind_young_cohort0 = sp_nb_cohort_lv[sp]+sp_nb_cohort_jv[sp];
		if (age_compute_habitat[sp][ind_young_cohort0]!=ind_young_cohort0){
			cout << "ERROR in age_compute_habitat vector! " <<
				"The habitat of the first young cohort should "<<
				"always be computed! Will exit now..." << endl; 
			exit(1);
		}

		//Maturity:
		//age at first maturity (usually age class with proportion of mature individuals = 0.1-0.2)
	        age_mature[sp] = doc.getInteger("/age_mature", sp_name[sp]);
		//maturity at age (=1 or estimated by MFCL proportions (optional))
		maturity_age[sp].allocate(0,sp_nb_cohorts[sp]-1);
		maturity_age[sp] = 1.0;
		for (int a=0; a<sp_nb_cohorts[sp]; a++){
			if (a<age_mature[sp]) maturity_age[sp][a] = 0.0;
		}

		if (!doc.get("/maturity_age").empty()){
			for (int a=0; a<sp_nb_cohorts[sp]; a++)
                		maturity_age[sp][a] = doc.getDouble("/maturity_age/"+ sp_name[sp], a);
			for (int a=0; a<sp_nb_cohorts[sp]; a++){
				if (maturity_age[sp][a]>0.15) {
					cout << "IMPORTANT: first mature age class is " << a << endl;
					age_mature[sp] = a;
						break;
				}
			}
		}
		migrations_by_maturity_flag = 0;
	}		

	///////////////////
	//NATURAL MORTALITY 
	///////////////////

	food_requirement_in_mortality.initialize();
	residual_competition = 1.0;
	for (int sp=0;sp<nb_species;sp++){
		if (use_ph1){
	        	M_inc_ph_a[sp]  = doc.getDouble("/M_inc_ph_a",  sp_name[sp]);
        		M_inc_ph_b[sp]  = doc.getDouble("/M_inc_ph_b",  sp_name[sp]);
		}		
		//maximal mortality rate due to predation
        	Mp_mean_max[sp]  = doc.getDouble("/Mp_mean_max",  sp_name[sp]);
		//slope coefficient in predation mortality function
        	Mp_mean_exp[sp]  = doc.getDouble("/Mp_mean_exp",  sp_name[sp]);
		//maximal mortality rate due to senescence
        	Ms_mean_max[sp]  = doc.getDouble("/Ms_mean_max",  sp_name[sp]);
		//slope coefficient in senescence mortality function
        	Ms_mean_slope[sp]= doc.getDouble("/Ms_mean_slope",sp_name[sp]);
		//variability of mortality rate locally due to habitat index
        	M_mean_range[sp] = doc.getDouble("/M_mean_range", sp_name[sp]);

		//flag for optional mortality penalizing according to food requirement
		if (!doc.get("/food_requirement_in_mortality",sp_name[sp]).empty()){
			food_requirement_in_mortality[sp] = doc.getInteger("/food_requirement_in_mortality", sp_name[sp]);
			//percentage of forage available for species assuming presence of other species in the habitat
			residual_competition[sp] = doc.getDouble("/residual_competition", sp_name[sp]);	
		}
	}

	////////////////////////////
	// COUPLING-FORAGING
	// uncoupled mode by default 
	////////////////////////////
	flag_coupling = false;
	if (!doc.get("/coupling","flag").empty()){
		int flag = doc.getInteger("/coupling", "flag");
		if (flag) flag_coupling = true;
	}
	//mean daily ration (weight forage/weight fish) for each species
	for (int sp=0;sp<nb_species;sp++) {
        	forage_ration[sp] = doc.getDouble("/forage_ration", sp_name[sp]);
	}

	/////////////////////////
	// END OF SPECIES SECTION  
	////////////////////////

	/////////////////////
	// FISHERIES SECTION    
	/////////////////////
	if ( (nb_species != 0)&& (sp_nb_cohort_ad[0]!=0) ) 
                nb_fishery = doc.getInteger("/nb_fishery", "value");
	if (nb_fishery) {
		catch_reso = deltaX/60.0;
		if (!doc.get("/degrade_fishery_reso_deg","value").empty()){
			catch_reso = doc.getDouble("/degrade_fishery_reso_deg", "value");
		}
		cpue_mult.allocate(0,nb_fishery-1);	//cpue multiplier in the likelihood (Inna 18/11/2010)
		cpue_mult.initialize();

		//fishery name (code with 1 letter and 1 number; e.g., L1 for long-line or S2 for purse-seine)
		list_fishery_name = Utilities::create1d(list_fishery_name, nb_fishery);
		for (int f=0;f<nb_fishery;f++) 
			list_fishery_name[f] = doc.get("/list_fishery_name", f);


		fishery_catch_units.allocate(0,nb_fishery-1);
		for (int f=0;f<nb_fishery;f++) 
			fishery_catch_units[f] = doc.getInteger("/fishery_catch_units",f);

		nb_catch_files = 0;
		for (int sp=0;sp<nb_species;sp++){
			if (doc.get("/file_catch_data").empty())
				file_catch_data.push_back("");
			if (!doc.get("/file_catch_data").empty())
				nb_catch_files = doc.getInteger(string("/file_catch_data/nb_files"),sp_name[sp]);
			for (int nf=0; nf<nb_catch_files; nf++){
				std::ostringstream ostr;
				ostr << "file" << nf+1;
				file_catch_data.push_back(str_dir_fisheries+doc.get("/file_catch_data/"+sp_name[sp],ostr.str()));
			}
		}
		for (int sp=0;sp<nb_species;sp++){
			nb_frq_files = 0;
			if (!doc.get("/file_frq_data").empty()){
				nb_frq_files = doc.getInteger(string("/file_frq_data/nb_files"),sp_name[sp]);
				if (nb_frq_files)
					for (int nf=0; nf<nb_frq_files; nf++){
						std::ostringstream ostr;
						ostr << "file" << nf+1;
						file_frq_data.push_back(str_dir_fisheries+doc.get("/file_frq_data/"+sp_name[sp],ostr.str()));
					}
			} else {
				cout << "WARNING: LF DATA ARE ABSENT" << endl; 
				file_frq_data.push_back("");
			}
		}

		//Number of fisheries by species (nb of fisheries catching a given species)
		nb_fishery_by_sp.allocate(0, nb_species - 1);
		nb_fishery_by_sp.initialize();

		mask_fishery_sp.allocate(0, nb_species - 1, 0, nb_fishery - 1);
			 
		for (int sp=0;sp<nb_species;sp++) {
			for (int f=0;f<nb_fishery;f++) {
				mask_fishery_sp[sp][f] = doc.getInteger(string("/mask_fishery_sp/") + sp_name[sp], f);
				if (mask_fishery_sp[sp][f]) 
					nb_fishery_by_sp[sp] ++;
			}
		}

		mask_fishery_sp_like.allocate(0, nb_species - 1, 0, nb_fishery - 1);
		mask_fishery_sp_like = mask_fishery_sp;
		if (!doc.get("/mask_fishery_likelihood").empty()){
			for (int sp=0;sp<nb_species;sp++) {
				if (!doc.get(string("/mask_fishery_likelihood/") + sp_name[sp]).empty()){
					for (int f=0;f<nb_fishery;f++) 
						mask_fishery_sp_like[sp][f] = doc.getInteger(string("/mask_fishery_likelihood/") + sp_name[sp], f);
				}
			}
		}
		mask_fishery_sp_no_effort.allocate(0,nb_species-1, 0,nb_fishery-1);
		mask_fishery_sp_no_effort.initialize(); //by default all fisheries have effort
		fisheries_no_effort_exist.allocate(0,nb_species-1);
		if (!doc.get("/mask_fishery_no_effort").empty()){
			for (int sp=0;sp<nb_species;sp++) {
				if (!doc.get(string("/mask_fishery_no_effort/") + sp_name[sp]).empty()){
					for (int f=0;f<nb_fishery;f++) 
						mask_fishery_sp_no_effort[sp][f] = doc.getInteger(string("/mask_fishery_no_effort/") + sp_name[sp], f);
				}
				//fisheries_no_effort_exist will be used as a flag to avoid catch removal from density of tags
cout << mask_fishery_sp_no_effort(sp)<< endl;				
				fisheries_no_effort_exist[sp] = sum(mask_fishery_sp_no_effort[sp]);
				
				//check: if no effort fisheries have different catch units in the parfile 
				//the catch removal method requires all fisheries catch to have the same units
				if (fisheries_no_effort_exist(sp)){
					int sum_units = 0;
					for (int f=0;f<nb_fishery;f++){
						if (mask_fishery_sp_no_effort(sp,f)){
							sum_units += fishery_catch_units(f);
						}
					}
					if (sum_units>0 && sum_units!=fisheries_no_effort_exist(sp)){
						cout <<"WRONG CONFIGURATION: Fisheries without effort in simulation - same catch units required for all 'no effort' fisheries! The program will exit, sorry."<< endl;
						exit(1);	
					}
				}
			}
		}
/*
		mask_mpa_fishery.allocate(0, nb_species - 1, 0, nb_fishery - 1);
		mask_mpa_fishery = 1;	
		if (mpa_simulation){
			for (int sp=0;sp<nb_species;sp++) {
				for (int f=0;f<nb_fishery;f++) 
					mask_mpa_fishery[sp][f] = doc.getInteger(string("/mpa_simulation/mpa_fisheries/") + sp_name[sp], f);
			}
		}
*/
		//to convert the units of LL from hooks to thous.hooks to avoid small value of catchability parameter
		//to be set by the user, by default = 1.0
		eff_units_converter.allocate(0,nb_fishery-1);
		eff_units_converter.initialize();
		eff_units_converter += 1.0;		

		if (!doc.get("/eff_units_converter").empty())
			for (int f=0;f<nb_fishery;f++) 
				eff_units_converter(f) = doc.getDouble("/eff_units_converter",f);
		
		like_types.allocate(0,nb_species-1);
		for (int sp=0;sp<nb_species;sp++){
			like_types[sp].allocate(0,nb_fishery_by_sp[sp]-1);
			int k = 0;
			for (int f=0;f<nb_fishery;f++){
				int type = doc.getInteger(string("/likelihood_types/")+sp_name[sp], f);
				if (mask_fishery_sp[sp][f]){
					like_types[sp][k] = type;
					k++;
				}
			}
		}

		//MPA simulations
                if (mpa_simulation){
                        mpa_scenario.allocate(0, nb_mpa-1);
                        mpa_ID.allocate(0, nb_mpa-1);
                        mpa_S1_X.allocate(0, nb_mpa-1);
                        mpa_fishery.allocate(0, nb_fishery - 1);
                        mpa_scenario.initialize();
                        mpa_ID.initialize();
                        mpa_S1_X.initialize();
                        mpa_fishery = 1;
                        for (int i=0;i<nb_mpa;i++){
                                mpa_ID[i] = doc.getInteger("/mpa_simulation/mpa_IDs", i);
                                mpa_scenario[i] = doc.getInteger("/mpa_simulation/mpa_scenario", i);
                                if (mpa_scenario[i]==1) mpa_S1_X[i] = doc.getInteger("/mpa_simulation/mpa_S1_X", i);
                        }
                        for (int f=0;f<nb_fishery;f++)
                                mpa_fishery[f] = doc.getInteger(string("/mpa_simulation/mpa_fisheries"), f);

                        str_file_maskMPA = str_dir + doc.get("/mpa_simulation/mpa_mask", "value");
                        actual_eff = 1;
                        if (!doc.get("/mpa_simulation/actual_eff", "use").empty())
                                if (doc.get("/mpa_simulation/actual_eff", "use") != "true")
                                        actual_eff = 0;
                }

	
		/////////////////////////
		//   LIKELIHOOD SECTION    
		/////////////////////////

		//Optimization control:

		maxfn = 1000; crit = 0.1;
		if (!doc.get("/max_nb_function_evaluations","value").empty())
			maxfn = doc.getInteger("/max_nb_function_evaluations","value");

		if (!doc.get("convergence_criterion","value").empty())
			maxfn = doc.getDouble("/convergence_criterion","value");
		//end of Optimization control:

		total_like = 0;
		if (!doc.get("/total_likelihood","value").empty())
			total_like = doc.getDouble("/total_likelihood","value");

		//1. CATCH/CPUE likelihood
		//the data that will be used in the likelihood: local catch or cpue
		cpue = false; //default is catch
		if (!doc.get("/like_c_cpue","value").empty()){
			int flag = doc.getInteger("/like_c_cpue", "value");
			if (!flag) {
				cout << "CPUE will be used in the likelihood" << endl;
				cpue = true;
			}
		}

		//2. LF likelihood
		frq_like.allocate(0,nb_species-1);
		frq_like.initialize();
		for (int sp=0;sp<nb_species;sp++)
			frq_like[sp] = doc.getInteger("/frq_likelihood",sp_name[sp]);

		//3. TAGs likelihood
		tag_like.allocate(0,nb_species-1);
		tag_like.initialize();
		tags_only = false;
		tag_gauss_kernel_on = 1;
		for (int sp=0;sp<nb_species;sp++){
			if (!doc.get("/tag_likelihood",sp_name[sp]).empty())
				tag_like[sp] = doc.getInteger("/tag_likelihood",sp_name[sp]);
			if (tag_like[sp]){
				if (!doc.get("/tag_likelihood_only","value").empty()){
					int flag = doc.getInteger("/tag_likelihood_only","value");
					if (flag) tags_only = true;
				}
				if (!doc.get("/tag_gauss_kernel_on","value").empty()){
					tag_gauss_kernel_on = doc.getInteger("/tag_gauss_kernel_on","value");
				}
			}
		}

		//4. STOCK likelihood
		stock_like.allocate(0,nb_species-1);
		mean_stock_obs.allocate(0,nb_species-1);
		stock_lonmin.allocate(0,nb_species-1);
		stock_lonmax.allocate(0,nb_species-1);
		stock_latmin.allocate(0,nb_species-1);
		stock_latmax.allocate(0,nb_species-1);
		stock_like.initialize();
		mean_stock_obs.initialize();
		stock_lonmin.initialize();
		stock_lonmax.initialize();
		stock_latmin.initialize();
		stock_latmax.initialize();

		for (int sp=0;sp<nb_species;sp++){
			if (!doc.get("/stock_likelihood",sp_name[sp]).empty())
				stock_like[sp] = doc.getInteger("/stock_likelihood",sp_name[sp]);
		
			if (stock_like[sp]){
				mean_stock_obs[sp] = doc.getDouble("/mean_stock_obs/"+sp_name[sp],"value"); 
				stock_lonmin[sp] = doc.getDouble("/mean_stock_obs/"+sp_name[sp],"lgmin"); 
				stock_lonmax[sp] = doc.getDouble("/mean_stock_obs/"+sp_name[sp],"lgmax"); 
				stock_latmin[sp] = doc.getDouble("/mean_stock_obs/"+sp_name[sp],"ltmin"); 
				stock_latmax[sp] = doc.getDouble("/mean_stock_obs/"+sp_name[sp],"ltmax"); 
			}
		}
		//likelihood parameters: variance, beta binomial 
		like_param.allocate(0, nb_species-1);
		for (int sp=0; sp<nb_species; sp++){
			like_param(sp).allocate(0, nb_fishery_by_sp[sp]-1);
			like_param(sp).initialize();
			int k=0;
			for (int f=0; f< nb_fishery; f++)
				if (mask_fishery_sp[sp][f]){
					if (like_types(sp,k)==2 || like_types(sp,k)==4 || like_types(sp,k)==5)
						like_param[sp][k] = doc.getDouble(string("/likelihood_parameters/") + list_fishery_name[f], sp_name[sp]);
					k++;
				}
		}

		for (int sp=0;sp<nb_species;sp++){
			dx_tags = 1; dy_tags = 1;
			lonmin_tags = longitudeMin; lonmax_tags = longitudeMax;
			latmin_tags = latitudeMin; latmax_tags = latitudeMax;
			if (sp_name[sp].compare("skj")==0){
				lonmin_tags = 110; lonmax_tags = 280;latmin_tags = -16;latmax_tags = 44;
			}
			if (sp_name[sp].compare("bet")==0){
				lonmin_tags = 130; lonmax_tags = 290;latmin_tags = -25;latmax_tags = 25;
			}
			if (sp_name[sp].compare("yft")==0){
				lonmin_tags = 110; lonmax_tags = 290;latmin_tags = -30;latmax_tags = 30;
			}
			if (sp_name[sp].compare("cjm")==0){
				lonmin_tags = 180; lonmax_tags = 290;latmin_tags = -50;latmax_tags = 3;
			}
			nb_tag_files = 0;
			if (tag_like(sp)){
				if (!doc.get("/tags_grid").empty()){
					dx_tags = doc.getInteger(string("/tags_grid/reso"),"dx");
					dy_tags = doc.getInteger(string("/tags_grid/reso"),"dy");
					if (dx_tags<deltaX/60){ 
						dx_tags=deltaX/60; 
						cout << "WARNING: changed the tagging data longitudinal resolution to deltaX = "
							<< deltaX/60 << " degree" << endl; 
					}
					if (dy_tags<deltaY/60){ 
						dy_tags=deltaY/60; 
						cout << "WARNING: changed the tagging data latitudinal resolution to deltaY = "
							<< deltaY/60 << " degree" <<endl; 
					}
					lonmin_tags = doc.getDouble(string("/tags_grid/longitude"),"east");
					lonmax_tags = doc.getDouble(string("/tags_grid/longitude"),"west");
					latmin_tags = doc.getDouble(string("/tags_grid/latitude"),"south");
					latmax_tags = doc.getDouble(string("/tags_grid/latitude"),"north");
				}
				nb_tag_files = doc.getInteger(string("/file_tag_data/nb_files"),sp_name[sp]);

				if (nb_tag_files){
					for (int nf=0; nf<nb_tag_files; nf++){
						std::ostringstream ostr;
						ostr << "file" << nf+1;
						file_tag_data.push_back(str_dir_tags+doc.get("/file_tag_data/"+sp_name[sp],ostr.str()));
					}
				} else {
					cout << "WARNING: TAG DATA ARE ABSENT" << endl; 
					file_tag_data.push_back("");
				}
			}
		}


		//likelihood parameters: probability of zero observation 
		prob_zero.allocate(0,nb_species-1);
		for (int sp=0;sp<nb_species;sp++){
			prob_zero(sp).allocate(0,nb_fishery_by_sp[sp]-1);
			int k = 0;
			for (int f=0; f<nb_fishery;f++){
			 	if (mask_fishery_sp[sp][f]){
					prob_zero[sp][k] = 0.5;
					if (!doc.get("/prob_zero").empty())
						prob_zero[sp][k] = doc.getDouble(string("/prob_zero/") + list_fishery_name[f],sp_name[sp]);
					k++;
				}	
			}
		}
	
		flag_twin = false;

		//catchability (q) by species and by fishery
		q_sp_fishery.allocate(0, nb_species - 1);

		for (int sp=0;sp<nb_species;sp++) {
			const int fmax = nb_fishery_by_sp[sp];
			q_sp_fishery[sp].allocate(0, fmax - 1);

			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[sp][f]){
					q_sp_fishery[sp][k] = doc.getDouble(string("/q_sp_fishery/") + list_fishery_name[f], sp_name[sp]);
					if (mask_fishery_sp_no_effort[sp][f] && doc.get("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true")
						doc.set("/q_sp_fishery/"+list_fishery_name[f]+"/variable","use","false");
					k++;
				}
			}
		}
		q_dyn_fishery.allocate(0, nb_fishery - 1);
		q_dyn_fishery.initialize();
		for (int f = 0; f < nb_fishery; f++){
			if (!doc.get("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "dyn").empty())
				q_dyn_fishery[f] = doc.getDouble("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "dyn");
		}

		//selectivity function parameters for species and fishery
		s_func_type.allocate(0,nb_fishery-1);
		for (int f = 0; f < nb_fishery; f++)
			s_func_type[f] = doc.getInteger(string("/s_sp_fishery/") + list_fishery_name[f]+string("/function_type"), "value");

		s_slope_sp_fishery.allocate(0, nb_species - 1);
		s_length_sp_fishery.allocate(0, nb_species - 1);
		s_asympt_sp_fishery.allocate(0, nb_species - 1);

		for (int sp=0;sp<nb_species;sp++) {
			const int fmax = nb_fishery_by_sp[sp];
			s_slope_sp_fishery[sp].allocate(0, fmax - 1);
			s_length_sp_fishery[sp].allocate(0, fmax - 1);
			s_asympt_sp_fishery[sp].allocate(0, fmax - 1);

			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[sp][f]){
					s_slope_sp_fishery[sp][k] = doc.getDouble(string("/s_sp_fishery/") + list_fishery_name[f], sp_name[sp]);
					if (s_func_type[f]>1)
						s_length_sp_fishery[sp][k] = doc.getDouble(string("/s_sp_fishery/")+list_fishery_name[f]+string("/length_threshold"),sp_name[sp]);
					if (s_func_type[f]>=3)
						s_asympt_sp_fishery[sp][k] = doc.getDouble(string("/s_sp_fishery/")+list_fishery_name[f]+string("/right_asymptote"),sp_name[sp]);
					k++;
				}
			}
		}
	}	
	////////////////////////
	// ZONES OF AGGREGATION 
	////////////////////////
        use_mask_catch = 0;
	flex_regstruc = 0;
	if (!doc.get("/flex_regstruc","value").empty()){
		flex_regstruc = doc.getInteger("/flex_regstruc","value");
		// regions for asismilation of LF data
		if (flex_regstruc){
			for (int sp=0;sp<nb_species;sp++) 
				if (frq_like[sp])
					define_regions();
				else {
					cout << "WARNING: LF likelihood is not active -> flex_regstruc is turned OFF!" << endl; 
					flex_regstruc = 0;
				}
		}
	}
	if (!flex_regstruc){
		if (sum(mask_fishery_sp)>0 && !doc.get("/use_mask_catch","value").empty()) {
			use_mask_catch = doc.getInteger("/use_mask_catch","value");
		}
		nb_region = doc.getInteger("/nb_region", "value");
		cout << "WARNING: flex_regstruc is OFF, areas in parfile will be used as aggregation zones and to compute LFs!" << endl; 
		for (int sp=0;sp<nb_species;sp++) 
			if (frq_like[sp])
				cout << "ATTENTION: Verify that LF data for species " << sp << " have the same regional structure!" << endl; 
		if (nb_region) {
			area = new region*[nb_region];
			for (int a=0; a< nb_region; a++) {
        	               	std::ostringstream ostr;
                	       	ostr << "/area" << a;
				area[a] = new region();
	                       	area[a]->area_id = doc.getInteger(ostr.str(), "area_id");
	                       	area[a]->lgmin	 = doc.getDouble(ostr.str(), "lgmin");
				if (area[a]->lgmin<longitudeMin) area[a]->lgmin = longitudeMin;
	                       	area[a]->lgmax	 = doc.getDouble(ostr.str(), "lgmax");
				if (area[a]->lgmax>longitudeMax) area[a]->lgmax = longitudeMax;
	                       	area[a]->ltmin	 = doc.getDouble(ostr.str(), "ltmin");
				if (area[a]->ltmin<latitudeMin)  area[a]->ltmin = latitudeMin;
	                       	area[a]->ltmax	 = doc.getDouble(ostr.str(), "ltmax");
				if (area[a]->ltmax>latitudeMax)  area[a]->ltmax = latitudeMax;
			}
			//nb regions par espece ou aggreger les donnees de biomasses
			nb_region_sp_B.allocate(0, nb_species - 1);
			for (int sp=0;sp<nb_species;sp++){
				nb_region_sp_B[sp] = doc.getInteger("/nb_region_sp_B", sp_name[sp]);
				if (nb_region_sp_B[sp]>nb_region) {
					nb_region_sp_B[sp] = nb_region;
					cout << "WARNING: nb_region_sp_B cannot be > nb_region, value " << nb_region_sp_B[sp] << " will be taken" << endl;
				}
			}
			//liste des regions par espece ou aggreger les biomasses
			area_sp_B.allocate(0, nb_species - 1);
			for (int sp=0;sp<nb_species;sp++) {
				const int amax = nb_region_sp_B[sp];
				area_sp_B[sp].allocate(0, amax - 1);
				for (int a = 0; a < amax; a++) {
	                               	area_sp_B[sp][a] = doc.getInteger(string("/area_sp_B/") + sp_name[sp], a);
				}
			}
		} else {
			nb_region_sp_B.allocate(0, nb_species - 1);
			area_sp_B.allocate(0, nb_species - 1);
			for (int sp=0;sp<nb_species;sp++){
				nb_region_sp_B[sp] = 0;//1;
				area_sp_B[sp].allocate(0, 0);
				area_sp_B[sp][0] = 0;//1;
			}
			cout << "WARNING: nb_region=0, no PREDICTED LFs will be computed"<< endl; 
		}
	}	

	// Nombre et numeros des ZEE
        nb_EEZ = nb_species <= 0 ? 0 : doc.getInteger("/nb_EEZ", "value");

	if (nb_species == 0) 
		nb_EEZ=0;

//str_file_maskEEZ = str_dir + doc.get("/str_file_maskEEZ", "value");
       	if (nb_EEZ) {
		EEZ_name = Utilities::create1d(EEZ_name, nb_EEZ);
		EEZ_ID.allocate(0, nb_EEZ - 1);
               	str_file_maskEEZ = str_dir + doc.get("/str_file_maskEEZ", "value");
        	
		for (int i=0;i<nb_EEZ;i++) {
			std::ostringstream ostr,namestr;
			ostr << i;
	                namestr << "/eez" << i;

			EEZ_name[i] = doc.get(string("/EEZ")+namestr.str(),"name");
			EEZ_ID[i] = doc.getInteger(string("/EEZ")+namestr.str(),"id");
		}
	}
	connectivity_comp = false;
	if (!doc.get("/connectivity_comp","flag").empty()){
		int flag = doc.getInteger("/connectivity_comp", "flag");
		if (flag) connectivity_comp = true;
	}

	//////////////////////////
	// BOUNDS FOR DVARIABLES 
	//////////////////////////
	_nvarcalc = 0;
	int nni = 0;

	statpars_temp.allocate(0,999);
	statpar_names_temp.allocate(0,999);

	par_read_bounds(Mp_mean_max,Mp_mean_max_min,Mp_mean_max_max,"/Mp_mean_max",nni);
	par_read_bounds(Mp_mean_exp,Mp_mean_exp_min,Mp_mean_exp_max,"/Mp_mean_exp",nni);
	par_read_bounds(Ms_mean_max,Ms_mean_max_min,Ms_mean_max_max,"/Ms_mean_max",nni);
	par_read_bounds(Ms_mean_slope,Ms_mean_slope_min,Ms_mean_slope_max,"/Ms_mean_slope",nni);
	par_read_bounds(M_mean_range,M_mean_range_min,M_mean_range_max,"/M_mean_range",nni);
	par_read_bounds(a_sst_spawning,a_sst_spawning_min,a_sst_spawning_max,"/a_sst_spawning",nni);
	par_read_bounds(b_sst_spawning,b_sst_spawning_min,b_sst_spawning_max,"/b_sst_spawning",nni);
	par_read_bounds(a_sst_larvae,a_sst_larvae_min,a_sst_larvae_max,"/a_sst_larvae",nni);
	par_read_bounds(b_sst_larvae,b_sst_larvae_min,b_sst_larvae_max,"/b_sst_larvae",nni);
	par_read_bounds(alpha_hsp_prey,alpha_hsp_prey_min,alpha_hsp_prey_max,"/alpha_hsp_prey",nni);
	par_read_bounds(alpha_hsp_predator,alpha_hsp_predator_min,alpha_hsp_predator_max,"/alpha_hsp_predator",nni);
	par_read_bounds(beta_hsp_predator,beta_hsp_predator_min,beta_hsp_predator_max,"/beta_hsp_predator",nni);
	par_read_bounds(a_sst_habitat,a_sst_habitat_min,a_sst_habitat_max,"/a_sst_habitat",nni);
	par_read_bounds(b_sst_habitat,b_sst_habitat_min,b_sst_habitat_max,"/b_sst_habitat",nni);
	par_read_bounds(T_age_size_slope,T_age_size_slope_min,T_age_size_slope_max,"/T_age_size_slope",nni);
	thermal_func_delta_min.allocate(0,2);
	thermal_func_delta_max.allocate(0,2);
	for (int n=0; n<3; n++){
		std::ostringstream ostr;
		ostr << "delta" << n+1;
		par2_read_bounds(thermal_func_delta[n],thermal_func_delta_min[n],thermal_func_delta_max[n],
				"/thermal_func",ostr.str(),nni);
	}		
	par_read_bounds(a_oxy_habitat,a_oxy_habitat_min,a_oxy_habitat_max,"/a_oxy_habitat",nni);
	par_read_bounds(b_oxy_habitat,b_oxy_habitat_min,b_oxy_habitat_max,"/b_oxy_habitat",nni);
	eF_habitat_min.allocate(0, nb_forage-1);
	eF_habitat_max.allocate(0, nb_forage-1);
	for (int n=0; n<nb_forage; n++){
		par2_read_bounds(eF_habitat[n],eF_habitat_min[n],eF_habitat_max[n],"/eF_habitat",frg_name[n],nni);
	}	
	par_read_bounds(hp_cannibalism,hp_cannibalism_min,hp_cannibalism_max,"/hp_cannibalism",nni);
	par_read_bounds(sigma_species,sigma_species_min,sigma_species_max,"/sigma_species",nni);
	par_read_bounds(MSS_species,MSS_species_min,MSS_species_max,"/MSS_species",nni);
	par_read_bounds(MSS_size_slope,MSS_size_slope_min,MSS_size_slope_max,"/MSS_size_slope",nni);
	par_read_bounds(c_diff_fish,c_diff_fish_min,c_diff_fish_max,"/c_diff_fish",nni);
	par_read_bounds(nb_recruitment,nb_recruitment_min,nb_recruitment_max,"/nb_recruitment",nni);
	par_read_bounds(a_adults_spawning,a_adults_spawning_min,a_adults_spawning_max,"/a_adults_spawning",nni);
	par_read_bounds(spawning_season_peak,spawning_season_peak_min,spawning_season_peak_max,"/spawning_season_peak",nni);
	par_read_bounds(spawning_season_start,spawning_season_start_min,spawning_season_start_max,"/spawning_season_start",nni);


	if (doc.get("/q_sp_fishery/variables", "use") == "true") {
		q_sp_fishery_min.allocate(0, nb_species-1);
		q_sp_fishery_max.allocate(0, nb_species-1);
		for (int sp=0;sp<nb_species;sp++) {
			const int fmax = nb_fishery_by_sp[sp];
			q_sp_fishery_min[sp].allocate(0, fmax-1);
			q_sp_fishery_min[sp].initialize();
			q_sp_fishery_max[sp].allocate(0, fmax-1);
			q_sp_fishery_max[sp].initialize();
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[sp][f]) {
					if (doc.get("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true"){
						_nvarcalc++; 
						q_sp_fishery_min[sp][k] = doc.getDouble("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "min");
						q_sp_fishery_max[sp][k] = doc.getDouble("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "max");
					} else {
						statpar_names_temp[nni] = "q sp fishery(" + str(sp) +","+ str(k) + ")"; 
						statpars_temp[nni] = q_sp_fishery[sp][k]; 
						nni++;
					}
					k++;
				}
			}
		}
	} else {
		for (int sp=0;sp<nb_species;sp++){
		int k = 0;
		for (int f = 0; f < nb_fishery; f++) 
			if (mask_fishery_sp[sp][f]) {
				statpar_names_temp[nni] = "q sp fishery(" + str(sp) +","+ str(k) + ")"; 
				statpars_temp[nni] = q_sp_fishery[sp][k]; 
				nni++;
				k++;
			}
		}	
	}

	if (doc.get("/s_sp_fishery/variables", "use") == "true") {
		s_slope_sp_fishery_min.allocate(0, nb_species-1);
		s_slope_sp_fishery_max.allocate(0, nb_species-1);

		s_asympt_sp_fishery_min.allocate(0, nb_species-1);
		s_asympt_sp_fishery_max.allocate(0, nb_species-1);

		for (int sp=0;sp<nb_species;sp++) {
			const int fmax = nb_fishery_by_sp[sp];
			s_slope_sp_fishery_min[sp].allocate(0, fmax-1);
			s_slope_sp_fishery_min[sp].initialize();
			s_slope_sp_fishery_max[sp].allocate(0, fmax-1);
			s_slope_sp_fishery_max[sp].initialize();
			s_asympt_sp_fishery_min[sp].allocate(0, fmax-1);
			s_asympt_sp_fishery_min[sp].initialize();
			s_asympt_sp_fishery_max[sp].allocate(0, fmax-1);
			s_asympt_sp_fishery_max[sp].initialize();
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[sp][f]) {
					if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true"){
						_nvarcalc++;
						s_slope_sp_fishery_min[sp][k] = doc.getDouble("/s_sp_fishery/"+list_fishery_name[f]+"/variable", "min");
						s_slope_sp_fishery_max[sp][k] = doc.getDouble("/s_sp_fishery/"+list_fishery_name[f]+"/variable", "max");
					} 
					else {
						statpar_names_temp[nni] = "s sp fishery(" + str(sp) +","+ str(k) + ")"; 
						statpars_temp[nni] = s_slope_sp_fishery[sp][k]; 
						nni++;
					}		
					if (s_func_type[f]>1){
						if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/length_threshold", "use") == "true")
							_nvarcalc++; 
						else {
							statpar_names_temp[nni] = "length threshold(" + str(sp) +","+ str(k) + ")"; 
							statpars_temp[nni] = s_length_sp_fishery[sp][k]; 
							nni++;
						}
					}
					if (s_func_type[f]>=3){
						if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote", "use") == "true"){
							_nvarcalc++;
							s_asympt_sp_fishery_min[sp][k] = doc.getDouble("/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote", "min");
							s_asympt_sp_fishery_max[sp][k] = doc.getDouble("/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote", "max");
						}
						else {
							statpar_names_temp[nni] = "right asymptote(" + str(sp) +","+ str(k) + ")"; 
							statpars_temp[nni] = s_asympt_sp_fishery[sp][k]; 
							nni++;
						}
					}
					k++;
				}
			}
		}
	} else {
		for (int sp=0;sp<nb_species;sp++){
		int k = 0;
		for (int f = 0; f < nb_fishery; f++) 
			if (mask_fishery_sp[sp][f]) {
				statpar_names_temp[nni] = "s sp fishery(" + str(sp) +","+ str(k) + ")"; 
				statpars_temp[nni] = s_slope_sp_fishery[sp][k]; 
				nni++;	

				if (s_func_type[f]>1){
					statpar_names_temp[nni] = "s length threshold(" + str(sp) +","+ str(k) + ")"; 
					statpars_temp[nni] = s_length_sp_fishery[sp][k]; 
					nni++;
				}
				if (s_func_type[f]>=3){
					statpar_names_temp[nni] = "s right asymptote(" + str(sp) +","+ str(k) + ")"; 
					statpars_temp[nni] = s_asympt_sp_fishery[sp][k]; 
					nni++;
				}
				k++;
			}
		}
	}

	if (doc.get("/likelihood_parameters/variables","use") == "true"){
		for (int sp=0; sp<nb_species; sp++){
			int k = 0;
			for (int f=0; f< nb_fishery; f++)
				if (mask_fishery_sp[sp][f]) {
					if (like_types[sp][k]==4||like_types[sp][k]==5)
						_nvarcalc++;
					else if (like_types[sp][k]==2) {
						statpar_names_temp[nni] = "likelihood parameter(" + str(sp) +","+ str(k) + ")";
						statpars_temp[nni] = like_param[sp][k];
						nni++;
					}	
					k++;
				}
		}
		
	} else {
		for (int sp=0;sp<nb_species;sp++){
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) 
				if (mask_fishery_sp[sp][f]) {
					statpar_names_temp[nni] = "likelihood parameter(" + str(sp) +","+ str(k) + ")"; 
					statpars_temp[nni] = like_param[sp][k]; 
					nni++;
					k++;
				}
		}
	}

	for (int sp=0; sp<nb_species; sp++)
		for (int k=0; k<nb_fishery_by_sp[sp]; k++)
			if (like_types[sp][k]==5 && doc.get("/likelihood_parameters/variables","use") == "true")
				_nvarcalc++;
	
	if (_nvarcalc==0) {cerr << "ERROR: at least one control variable required!!! Will exit now..."<< endl; exit(1);}

	parfile_names = Utilities::create1d(parfile_names, _nvarcalc);
	dvarpars.allocate(1,_nvarcalc);
	dvarpars.initialize();
	dvarpars_min.allocate(1,_nvarcalc);
	dvarpars_max.allocate(1,_nvarcalc);
	dvarpars_min.initialize();
	dvarpars_max.initialize();
	if (nni==0) nni++; 
	statpars.allocate(0,nni-1);
	statpar_names.allocate(0,nni-1);
	for (int i=0; i<nni; i++){
		statpars[i] = statpars_temp[i];
		statpar_names[i] = statpar_names_temp[i];
	}
	_nstatpars = nni;


	/////////////////////
	//OTHER RUNNING MODES 
	/////////////////////

	//2. Computing likelihood hyperspace projections
	nb_varproj = 0;
	if (!doc.get("/hyperspace_projection").empty()){
		nb_varproj = doc.getInteger("/hyperspace_projection/variables", "nb");
		varproj_nsteps.allocate(0,nb_varproj); varproj_nsteps.initialize();
		for (int n=0; n<nb_varproj; n++){
			std::ostringstream ostr;
			ostr << "var" << n+1;
			varproj.push_back(doc.get("/hyperspace_projection/" + ostr.str(),"name"));
			varproj_nsteps[n] = doc.getInteger("/hyperspace_projection/"+ ostr.str(), "nsteps");
		}	
	}
	return true;
}

void VarParamCoupled::re_read_varparam(){

	par_read(Mp_mean_max_min,Mp_mean_max_max,"/Mp_mean_max",0,1);
	par_read(Mp_mean_exp_min,Mp_mean_exp_max,"/Mp_mean_exp",0,10);
	par_read(Ms_mean_max_min,Ms_mean_max_max,"/Ms_mean_max",0,1);
	par_read(Ms_mean_slope_min,Ms_mean_slope_max,"/Ms_mean_slope",0,10);
	par_read(M_mean_range_min,M_mean_range_max,"/M_mean_range",0,100);
	par_read(a_sst_spawning_min,a_sst_spawning_max,"/a_sst_spawning",0,10);
	par_read(b_sst_spawning_min,b_sst_spawning_max,"/b_sst_spawning",0,34);
	par_read(a_sst_larvae_min,a_sst_larvae_max,"/a_sst_larvae",0,10);
	par_read(b_sst_larvae_min,b_sst_larvae_max,"/b_sst_larvae",0,34);
	par_read(alpha_hsp_prey_min,alpha_hsp_prey_max,"/alpha_hsp_prey",0,10000);
	par_read(alpha_hsp_predator_min,alpha_hsp_predator_max,"/alpha_hsp_predator",0,10000);
	par_read(beta_hsp_predator_min,beta_hsp_predator_max,"/beta_hsp_predator",0,10000);
	par_read(a_sst_habitat_min,a_sst_habitat_max,"/a_sst_habitat",0,10);
	par_read(b_sst_habitat_min,b_sst_habitat_max,"/b_sst_habitat",0,34);
	par_read(T_age_size_slope_min,T_age_size_slope_max,"/T_age_size_slope",0,10);
	for (int n=0; n<3; n++){
		std::ostringstream ostr;
		ostr << "delta" << n+1;	
		string s = "thermal_func" + ostr.str();
		par_read(thermal_func_delta_min[n],thermal_func_delta_max[n],s,0,0.5);
	}	
	par_read(a_oxy_habitat_min,a_oxy_habitat_max,"/a_oxy_habitat",0,1);
	par_read(b_oxy_habitat_min,b_oxy_habitat_max,"/b_oxy_habitat",0,10);
	for (int n=0; n<nb_forage; n++){
		string s = "/eF_habitat" + frg_name[n];
		par_read(eF_habitat_min[n],eF_habitat_max[n],s,0,100);
	}	
	par_read(hp_cannibalism_min,hp_cannibalism_max,"/hp_cannibalism",0,100);
	par_read(sigma_species_min,sigma_species_max,"/sigma_species",0,1);
	par_read(MSS_species_min,MSS_species_max,"/MSS_species",0,100);
	par_read(MSS_size_slope_min,MSS_size_slope_max,"/MSS_size_slope",0.1,2);
	par_read(c_diff_fish_min,c_diff_fish_max,"/c_diff_fish",0,1);
	par_read(nb_recruitment_min,nb_recruitment_max,"/nb_recruitment",0,100);
	par_read(a_adults_spawning_min,a_adults_spawning_max,"/a_adults_spawning",0,1000);
	par_read(spawning_season_peak_min,spawning_season_peak_max,"/spawning_season_peak",0,366);
	par_read(spawning_season_start_min,spawning_season_start_max,"/spawning_season_start",0.9,1.5);


	if (doc.get("/q_sp_fishery/variables", "use") == "true") {
		for (int sp=0;sp<nb_species;sp++) {
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[sp][f]) {
					if (doc.get("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true"){
						double min_val = doc.getDouble("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "min");
						if (min_val<0) {
							min_val = 0;
							cout << "restored min to 0" << endl;
						}
						q_sp_fishery_min[sp][k] = min_val;

						double max_val = doc.getDouble("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "max");
						if (max_val>1){
							max_val = 1;
							cout << "restored min to 1" << endl;
						}
						q_sp_fishery_max[sp][k] = max_val;
					}
					k++;
				}
			}
		}
	}

	if (doc.get("/s_sp_fishery/variables", "use") == "true") {
		for (int sp=0;sp<nb_species;sp++) {
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[sp][f]) {
					if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true"){
						double min_val = doc.getDouble("/s_sp_fishery/"+list_fishery_name[f]+"/variable", "min");
						if (min_val<0) {
							min_val = 0;
							cout << "restored min to 0" << endl;
						}
						
						double max_val = doc.getDouble("/s_sp_fishery/"+list_fishery_name[f]+"/variable", "max");					
						s_slope_sp_fishery_min[sp][k] = min_val;
						s_slope_sp_fishery_max[sp][k] = max_val;
					} 
					if (s_func_type[f]>=3){
						if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote", "use") == "true"){
							double min_val = doc.getDouble("/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote", "min");
							if (min_val<0) {
								min_val = 0;
								cout << "restored min to 0" << endl;
							}
							double max_val = doc.getDouble("/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote", "max");					
							if (max_val>1){
								max_val = 1;
								cout << "restored min to 1" << endl;
							}								
							s_asympt_sp_fishery_min[sp][k] = min_val;
							s_asympt_sp_fishery_max[sp][k] = max_val;
						}
					}
					k++;
				}
			}
		}
	}
	return;
}

void CParam::define_regions()
{ //create the regional structure from LF file

	for (int f=0; f<nb_frq_files; f++){
		string filename = file_frq_data[f];
		ifstream littxt(filename.c_str());
		if (littxt){
			int nb_fleets, nb_records, nb_region_file;

			littxt >> nb_region_file >> nb_fleets >> nb_records;
			nb_region += nb_region_file;
		} else {
			cout << endl << "WARNING : Cannot read LF data file with regions" << filename.c_str() << endl;
		}
	}
	area = new region*[nb_region];
	int a = 0; int r = 0;
	for (int f=0; f<nb_frq_files; f++){
		string filename = file_frq_data[f];
		ifstream littxt(filename.c_str());
		if (littxt){
			int nb_fleets, nb_records, nb_region_file;
			float lonmin, lonmax, latmin, latmax;
	
			littxt >> nb_region_file >> nb_fleets >> nb_records;
			for (int reg=0; reg< nb_region_file; reg++) {
				littxt >> r >> lonmin >> lonmax >> latmin >> latmax;  
				if (lonmin<0) lonmin += 360;
				if (lonmax<0) lonmax += 360;
				area[a] = new region();
	                       	area[a]->area_id = r;
	                       	area[a]->lgmin	 = lonmin;
				if (area[a]->lgmin<longitudeMin) area[a]->lgmin = longitudeMin;
	                       	area[a]->lgmax	 = lonmax;
				if (area[a]->lgmax>longitudeMax) area[a]->lgmax = longitudeMax;
	                       	area[a]->ltmin	 = latmin;
				if (area[a]->ltmin<latitudeMin)  area[a]->ltmin = latitudeMin;
	                       	area[a]->ltmax	 = latmax;
				if (area[a]->ltmax>latitudeMax)  area[a]->ltmax = latitudeMax;
				a++;
			}
		} 
	}

	//this part is temporal, just to keep old definitions to avoid a lot of changes further in the code
	//nb regions par espece ou aggreger les donnees de biomasses
	nb_region_sp_B.allocate(0, nb_species - 1);
	for (int sp=0;sp<nb_species;sp++)
		nb_region_sp_B[sp] = nb_region;
	//liste des regions par espece ou aggreger les biomasses
	area_sp_B.allocate(0, nb_species - 1);
	for (int sp=0;sp<nb_species;sp++) {
		const int amax = nb_region_sp_B[sp];
		area_sp_B[sp].allocate(0, amax - 1);
		for (int a = 0; a < amax; a++) {
			//area_sp_B[sp][a] = a+1; //wont work in multispecies mode!!!
			area_sp_B[sp][a] = area[a]->area_id; //wont work in multispecies mode!!!
		}
	}
}

void VarParamCoupled::par_read_bounds(dvector& var, double& var_min, double& var_max, string s, int& nni)
{
	string n = s;
	std::string whitespace = " ";
	n.erase(0,1);
	int pos = n.find("_");
	while (pos != -1) {
		n.replace(pos,1,whitespace);
		pos = n.find("_");
	}
	adstring name = n.c_str();

	string sv = s + "/variable";
	if (doc.get(sv, "use") == "true") {
		_nvarcalc++;
		var_min = doc.getDouble(sv, "min");
		var_max = doc.getDouble(sv, "max");
	} else {
		for (int sp=0; sp<nb_species; sp++){
			statpar_names_temp[nni] = name +"(" + str(sp) + ")"; 
			statpars_temp[nni] = var[sp]; 
			nni++;
		}
	} 
}

void VarParamCoupled::par_read(double& varmin, double& varmax, string s, const double fixmin, const double fixmax)
{
	string n = s;
	std::string whitespace = " ";
	n.erase(0,1);
	int pos = n.find("_");
	while (pos != -1) {
		n.replace(pos,1,whitespace);
		pos = n.find("_");
	}
	adstring name = n.c_str();

	string sv = s + "/variable";
	if (doc.get(sv, "use") == "true") {
		for (int sp=0; sp<nb_species; sp++){
			double var_min = doc.getDouble(sv, "min");
			if (var_min < fixmin) { 
				var_min = fixmin; 
				cout << "restored min to " << fixmin << endl;
			}
			double var_max = doc.getDouble(sv, "max");
			if (var_max > fixmax) {
				var_max = fixmax;
				cout << "restored max to " << fixmax << endl;
			}
			varmin = var_min;
			varmax = var_max;
		}
	}
}

void VarParamCoupled::par2_read_bounds(dvector& var, double& var_min, double& var_max, string s1, string s2, int& nni)
{
	string n = s1 + "_" + s2;
	std::string whitespace = " ";
	n.erase(0,1);
	int pos = n.find("_");
	while (pos != -1) {
		n.replace(pos,1,whitespace);
		pos = n.find("_");
	}
	adstring name = n.c_str();

	string sv = s1 + "/" + s2 + "/variable";
	if (doc.get(sv, "use") == "true") {
		_nvarcalc++;
		var_min = doc.getDouble(sv, "min");
		var_max = doc.getDouble(sv, "max");
	} else {
		for (int sp=0; sp<nb_species; sp++){
			statpar_names_temp[nni] = name +"(" + str(sp) + ")"; 
			statpars_temp[nni] = var[sp]; 
			nni++;
		}
	} 
}



