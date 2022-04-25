#include "SeapodymCoupled.h"

///This is the main loop for the model without fishing and fitting of density. 
///Similar to the default function, it includes the following calls:
///1- Initialising population density 
///2- Reading all forcing data (once in optimization mode, at every time step in simulation mode) 
///3- Age/lifestage loop calling the ADRE solvers and ageing
///4- Computing model density -> used as predictions
///5- Likelihood computation using input density as observations
///6- Writing outputs (in simulation mode only)
///Note, there is no modelling of tagged cohorts here!

//////////////////////////////////////////////////////////////////
//--------------------------------------------------------------//
//		     FORAGE-TUNA SIMULATION			//
//--------------------------------------------------------------//
//////////////////////////////////////////////////////////////////
/*!
\brief The tuna population simulation without fishing and density fitting.
*/
double SeapodymCoupled::OnRunDensity(dvar_vector x, const bool writeoutputfiles)
{
	InitializeAll();

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
	dvar_matrix Density_pred; 

	Habitat.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	Mortality.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Spawning_Habitat.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Total_pop.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Density_pred.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
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

	//Create time vector for output DYM files
	dvector zlevel;
        zlevel.allocate(0, nbt_total - 1);
        zlevel.initialize();

	//Output DYM file name
	string fileout;
	fileout = param->strdir_output + param->sp_name[0] + "_density_output.dym";
	if (writeoutputfiles){
	        Date::zlevel_run(*param,mat.zlevel,nbt_total,zlevel,nbt_start_series);

		// Create and initialize (dym) files for saving spatial variables
		//rewrite dym mask by the mask used in the run, i.e. map.carte:
		for (int i=0; i<nbi-2; i++){
		        for (int j=0; j<nbj-2; j++){
				mat.mask[j][i] = map.carte[i+1][j+1];
			}
		}
		double minval=0.0;
		double maxval=1.0;
		rw.wbin_header(fileout, param->idformat, param->idfunc, minval, maxval,
                                        param->nlong, param->nlat, nbt_total,
                                        zlevel[0], zlevel[nbt_total-1],
                                        mat.xlon, mat.ylat, zlevel, mat.mask);
	
		if (!param->gcalc())
			ConsoleOutput(0,0);
	}
	//Vector of zeroes, to avoid duplicating the F_accessibility function 
	ivector  tags_age_habitat;
	tags_age_habitat.allocate(0,aN_adult(0));
	tags_age_habitat.initialize();

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
		//Spin-up control: to be removed
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
 			func.Faccessibility(*param, mat, map, sp, jday, tcur, pop_built, false, tags_age_habitat);//checked

			//1.2 Food Requirement by population
			if (param->food_requirement_in_mortality(sp)){
				FR_pop_comp(FR_pop, sp);
//				ISR_denom_comp(ISR_denom, sp, tcur);
			}

			//2. IMPLICIT AGE LOOP: increment age while moving through life stages 
			int age = 0;	

			//2.1 Spawning habitat	
			func.Spawning_Habitat(*param, mat, map, Spawning_Habitat, 1.0, sp, tcur, jday);
			
			//2.2 Transport and mortality of larvae (always one age class)	
			for (int n=0; n<param->sp_nb_cohort_lv[sp]; n++){
				double mean_age = mean_age_cohort[sp][age]; 
				func.Mortality_Sp(*param, mat, map, Mortality, Spawning_Habitat, sp, mean_age, age, tcur);
				pop.Precalrec_juv(map, mat, Mortality, tcur);//checked
				pop.Calrec_juv(map, mat, mat.dvarDensity[sp][age], Mortality, tcur);//checked
				age++;
			}
			//2.3. Juvenile habitat	
			if (param->cannibalism[sp]){
				Total_Pop_comp(Total_pop,sp,jday,tcur); //adjoint
				func.Juvenile_Habitat_cannibalism(*param, mat, map, Habitat, Total_pop, sp, tcur);
			} else 
				func.Juvenile_Habitat(*param, mat, map, Habitat, sp, tcur);
				
			//2.4. Transport and mortality of juvenile age classes	
			for (int n=0; n<param->sp_nb_cohort_jv[sp]; n++){			
				double mean_age = mean_age_cohort[sp][age];

				func.Mortality_Sp(*param, mat, map, Mortality, Habitat, sp, mean_age, age, tcur);
				pop.Precalrec_juv(map,  mat, Mortality, tcur);
				pop.Calrec_juv(map, mat, mat.dvarDensity[sp][age], Mortality, tcur);
				age++;
			}
			
			//4. Transport and mortality of adult cohort
			///if (t_count > nbt_spinup_forage + nt_yn){ TO BE FIXED!!!
			///for (int age=0; age<=nb_age_built[sp]; age++){///TO BE FIXED!!!	
			for (int n=0; n<param->sp_nb_cohort_ad[sp]; n++){

				//NOTE: currently current averaging doesn't depend on seasonal migrations
				if (param->vert_movement[sp] && param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]){

					func.Average_currents(*param, mat, map, age, tcur, pop_built);
				}
					
				//this section will work only if seasonality switch is ON 
				//and for <1 maturity at age parameter
				if (param->migrations_by_maturity_flag && param->seasonal_migrations[sp]){
					cout << "In this version no smooth maturity; enter the age at first maturity. Exit now!" << endl; exit(1);
					
				} // end of section with seasonality switch and <1 maturity at age parameter 
				else {
					migration_flag = 0;
					if (age>=param->age_mature[sp] && param->seasonal_migrations[sp]) migration_flag = 1;

					if (param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]) {
						func.Feeding_Habitat(*param,mat,map,Habitat,sp,age,jday,tcur,migration_flag);
					}
					double mean_age = mean_age_cohort[sp][age];

					if (!param->food_requirement_in_mortality(sp)){
						func.Mortality_Sp(*param, mat, map, Mortality, Habitat, sp, mean_age, age, tcur);//checked
					} else {
						Food_Requirement_Index(IFR, FR_pop, ISR_denom, sp, age, tcur, jday);
						func.Mortality_Sp(*param, mat, map, Mortality, IFR, sp, mean_age, age, tcur);//checked
					}

					if (param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]){
						pop.Precaldia_Caldia(map, *param, mat, Habitat, Total_pop, sp, age, tcur, jday);//checked	
					}
					if (!param->gcalc()){
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
						
					//store adult habitat and update the counter
					if (param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]){
						//cout << age << " " << param->age_compute_habitat[sp][age] << endl;
						mat.adult_habitat(sp,tcur,param->age_compute_habitat[sp][age]) = value(Habitat);
					}
					pop.Precalrec_Calrec_adult(map,mat,*param,rw,
							mat.dvarDensity[sp][age],Mortality,
							tcur,fishing,age,sp,year,month,
							jday,step_fishery_count,0);//checked 20150210	
				}
				age++;
			}
			Density_pred.initialize();
			for (int age= a0_adult[sp]; age<aN_adult[sp]; age++){
				Density_pred += mat.dvarDensity[sp][age]* param->weight[sp][age] * 0.001;
			}
			
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
			SpawningBiomass_comp(Total_pop, sp);
			//----------------------------------------------//
			//	    TUNA AGEING AND SPAWNING		//
			//----------------------------------------------//

			//6. Ageing and survival
			if (nt_dtau==dtau){
				for (int a=param->sp_nb_cohorts[sp]-1; a >= 1; a--){
						Survival(mat.dvarDensity[sp][a], mat.dvarDensity[sp][a-1] , a, sp);
				}
				nt_dtau=0;
			}
			//7. Spawning
			Spawning(mat.dvarDensity[sp][0],Spawning_Habitat,Total_pop,jday,sp,pop_built,tcur);//checked
		}//end of 'sp' loop

		//------------------------------------------------------//
		//		COMPUTING LIKELIHOOD			//
		//------------------------------------------------------//
		//I. Total abundance likelihood: augmented only if param->stock_like = true
		//Can be useful to constrain the total abundance if the spatial distribution
		//fitting alone cannot provide the same abundance as in density_input.
		if (t_count == nbt_total)
			stocklike += get_stock_like(total_stock, likelihood);
		
		//II. Biomass density likelihood. Note, degrade it to the resolution of the density_input
		if (t_count > nbt_building+nbstoskip){
		int rr_x = (int)nlon/nlon_input; 
		int rr_y = (int)nlat/nlat_input; 
		TTRACE(rr_x,rr_y)
		if (rr_x<1 || rr_y<1){ 
			cerr << "Model resolution should be divisible without remainder by the resolution of input density field." << 
			         endl<< "Currenly nlon/nlon_input = " << rr_x << ", and nlat/nlat_input = " << rr_y << endl << "Will exit now...";
			exit(1);
		}
		int i0 = 0; int iN = nlon_input; //put 40-160 for 135E-110W in PO-2deg configuration
		int j0 = 0; int jN = nlat_input; //put 40-90 for 25S-25N in PO-2deg configuration
		for (int sp=0; sp < nb_species; sp++){
                        for (int i=i0; i < iN; i++){
				for (int j=j0; j<jN ; j++){
					if (mat.density_input(t_count,i,j)){
						dvariable Bsum = 0.0;
						for (int ii=0; ii<rr_x; ii++)
						for (int jj=0; jj<rr_y; jj++){
							if (map.carte(rr_x*i+ii,rr_y*j+jj)){
							       	Bsum += Density_pred[rr_x*i+ii][rr_y*j+jj];
							}
						}
							
						if (Bsum>0)	 
	                                		likelihood += (rr_x*rr_y*mat.density_input(t_count,i,j)-Bsum)*
							      	      (rr_x*rr_y*mat.density_input(t_count,i,j)-Bsum);
                                        }
                                }
                        }
		}}

		if (writeoutputfiles){
			if (!param->gcalc())	
				ConsoleOutput(1,value(likelihood));

			dmatrix mat2d(0, nbi - 1, 0, nbj - 1);
                        mat2d.initialize();
                        for (int i=map.imin; i <= map.imax; i++){
                                for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
                                        if (map.carte[i][j]){
                                                mat2d(i-1,j-1) = value(Density_pred(i,j));
                                        }
                                }
                        }
			double minval=0.0;
			double maxval=1.0;

                        rw.wbin_transpomat2d(fileout, mat2d, nbi-2, nbj-2, true);
			//update min-max values in header
			rw.rwbin_minmax(fileout, minval, maxval);
		}
		nt_dtau++;
		past_month=month;
		step_count++;
		if (qtr != past_qtr) past_qtr = qtr; 

	} // end of simulation loop

	param->total_like = value(likelihood);
	cout << "end of forward run, likelihood: " << value(likelihood)-eFlike << " " << eFlike << endl;

	return value(likelihood);
}

void SeapodymCoupled::ReadDensity()
{
	cout << "Reading input density file... " << endl;
	int jday = 0;
	int t_count_init = t_count;
			
	string file_input;
	file_input = param->strdir_output + param->sp_name[0] + "_density_input.dym";

	int nlevel = 0;
	rw.rbin_headpar(file_input, nlon_input, nlat_input, nlevel);
	cout << nlon_input << " "<< nlat_input << " " << nlevel << endl;

	mat.density_input.allocate(1,nbt_total);
	for (int t=1; t<=nbt_total; t++){
		mat.density_input[t].allocate(0, nlon_input, 0, nlat_input);
		mat.density_input[t].initialize();
	}
	cout << nbt_total << endl; 
	for (; t_count<=nbt_total; t_count++){
		getDate(jday);
		if (t_count > nbt_building) {
			//TIME SERIES 
			t_series = t_count - nbt_building + nbt_start_series;
			//cout << t_series << endl;
			//----------------------------------------------//
			//	READING DENSITY DATA			//
			//----------------------------------------------//
			int nbytetoskip = (9 +(3* nlat_input * nlon_input) + nlevel + ((nlat_input *nlon_input)* (t_series-1))) * 4;
			
			//rw.rbin_input2d(file_input, map, mat.density_input[t_count], nlon_input+2, nlat_input+2, nbytetoskip);
			ifstream litbin(file_input.c_str(), ios::binary | ios::in);
			if (!litbin)
			{
				cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_input << "\"\n";
				exit(1);
			}

			//---------------------------------------
			// Reading the 2d matrix
			//---------------------------------------
			litbin.seekg(nbytetoskip, ios::cur);


		 	const int sizeofDymInputType = sizeof(float);
			float buf;
			for (int j=0;j<nlat_input;j++)
			{
				for (int i=0;i<nlon_input;i++)
				{
					litbin.read(( char *)&buf,sizeofDymInputType);
					mat.density_input[t_count][i][j]= buf;
				}
			}
		
			litbin.close();

		}
	}
	t_count = t_count_init;
}


