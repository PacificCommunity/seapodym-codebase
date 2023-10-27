#include "SeapodymCoupled.h"
#include "NishikawaCategories.h"

///This is the main loop function for Habitat simulations. It is based on the default
///OnRunCoupled function, with modifications to compute either spawning of feeding habitats.
///See SeapodymCoupled_OnRunCoupled.cpp for the description of the main function.

void HabitatConsoleOutput(int t_count, int flag_simulation, string date_str, double norm_habitat_input, double norm_spawning_habitat, double like);



//////////////////////////////////////////////////////////////////
//--------------------------------------------------------------//
//		     HABITAT SIMULATION				//
//--------------------------------------------------------------//
//////////////////////////////////////////////////////////////////
/*!
\brief The main loop of habitat simulations.
*/
double SeapodymCoupled::OnRunHabitat(dvar_vector x, const bool writeoutputfiles)
{
	NishikawaCategories NshkwCat;

	t_count = nbt_building+1;
	mat.mats.initialize();

	past_month=month;
	past_qtr=qtr;

	//routine-specific variables
	int tcur = t_count; //will be used for forcing variable time control
	int nbt_no_forecast = t_count + nbt_spinup_tuna + nbt_total - 1;
	int migration_flag = 0;
	int step_count= 0;
	int jday = 0; 
	ivector Nobs(0,nb_fishery-1); Nobs.initialize();

	if (!param->gcalc()){
		
		if (!param->habitat_run_type)
			cout << endl << "Running Spawning Habitat model" << endl << endl; 		
		else 
			cout << endl << "Running Feeding Habitat model for ages: " << param->habitat_run_age << endl << endl; 		
		
		//need to read oxygen in case if month==past_month
		//(otherwise we may not have it for the first time steps)
		if (param->type_oxy==1 && month==past_month)
			ReadClimatologyOxy(1, month);
		//need to read oxygen in case if qtr==past_qtr 
		if (param->type_oxy==2 && qtr==past_qtr)
			ReadClimatologyOxy(1, qtr);
	}
	dvariable likelihood = 0.0;
	double likelihood_penalty = 50.0;
	double weight_Nobszero = 422.0/11123.0;
	//double weight_Nobszero = 0.0;

	// to remove following
	int test_inf_likelihood = 0;
	int nb_Npred_zero = 0;
	int nb_Npred_notzero = 0;
	//int nb_lkhd = 0;
	//int nb_lkhd125 = 0;
	//int nb_lkhd5 = 0;
	//double lkhd_tcount3 = 0.0;
	//int nb_lkhd_tcount3 = 0;

	reset(x);
	//----------------------------------------------//
	// 	LOCAL MATRICES ALLOCATION SECTION       //
	//----------------------------------------------//	
	
	dvar_matrix Habitat;

	Habitat.allocate(map.imin, map.imax, map.jinf, map.jsup);  
	Habitat.initialize();

	//precompute thermal habitat parameters
	for (int sp=0; sp < nb_species; sp++)
		func.Vars_at_age_precomp(*param,sp);

	//precompute seasonal switch function
	for (int sp=0; sp < nb_species; sp++){
		func.Seasonal_switch_year_precomp(*param,mat,map,
						value(param->dvarsSpawning_season_peak[sp]),
						value(param->dvarsSpawning_season_start[sp]),sp);
	}

	//DYM file names
	string fileout;
	if (param->habitat_run_type>0)
               	fileout = param->strdir_output + param->sp_name[0] + "_feeding_habitat_output.dym";
	else			
               	fileout = param->strdir_output + param->sp_name[0] + "_spawning_habitat_output.dym";
	
	//Write DYM headers
	if (writeoutputfiles){
		WriteFileHeaders_submodel(fileout);		
		if (!param->gcalc())
			HabitatConsoleOutput(0,0,date_str,0,0,0);
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
	//Not used currently, but can be useful to constrain eF parameters in habitat recalibration
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
		//spin-up flag: to be removed later
		int pop_built = 0; 
		if (t_count > nbt_building) 
			pop_built = 1;
		//------------------------------------------------------------------------------//
		//	TRANSPORT OF TUNA AGE CLASSES AND PREDICTED CATCH COMPUTATION		//
		//------------------------------------------------------------------------------//
		//------------------------------------------------------------------------------//

		for (int sp=0; sp < nb_species; sp++){

			if (param->habitat_run_type>0){
			//1. Precompute some variables outside of age loop
			//1.1 Accessibility by adults (all cohorts)
				ivector tags_age_habitat;
				tags_age_habitat.allocate(0,aN_adult(0));
				tags_age_habitat.initialize();
				func.Faccessibility(*param, mat, map, sp, jday, tcur, pop_built, 0, tags_age_habitat);//checked
			}

			//2. Starting habitat computation (implicit age loop) 
			int age = 0;	
			if (param->habitat_run_type==0){
				//2.1 Spawning habitat	
				func.Spawning_Habitat(*param, mat, map, Habitat, 1.0, sp, tcur, jday);
			}

			if (param->habitat_run_type>0){
				int rr_x = (int)nlon/nlon_input;
				int rr_y = (int)nlat/nlat_input;
				if (rr_x<1 || rr_y<1){ 
					cerr << "Model resolution should be divisible without remainder by the resolution of input density field." << 
					         endl<< "Currenly nlon/nlon_input = " << rr_x << ", and nlat/nlat_input = " << rr_y << endl << "Will exit now...";
					exit(1);
				}
			
				for (int ages=0; ages<param->nb_habitat_run_age; ages++){
					migration_flag = 0;
					if (age>=param->age_mature[sp] && param->seasonal_migrations[sp]) migration_flag = 1;

					int a;
					for (a=1; a<param->sp_nb_cohorts[sp]-1; a++){
						if (mean_age_cohort(sp,a)>param->habitat_run_age[ages]){
							age = mean_age_cohort(sp,a);
							break;
						}
					}
					if (a==param->sp_nb_cohorts[sp]-1 && age==0) age = param->sp_nb_cohorts[sp]-1;

					if (param->age_compute_habitat[sp][age]!=param->age_compute_habitat[sp][age-1]) {
						func.Feeding_Habitat(*param,mat,map,Habitat,sp,age,jday,tcur,migration_flag);
					}
					//------------------------------------------------------//
					//	     LIKELIHOOD FOR FEEDING HABITAT		//
					//------------------------------------------------------//
					for (int i=1; i<nlon_input; i++){
						for (int j=1; j<nlat_input; j++){
							double H_obs  = mat.habitat_input(ages,t_count,i,j);
							int ncount = 0;
							if (H_obs>0){
								dvariable Hmean_pred = 0.0;
								for (int ii=0; ii<rr_x; ii++)
								for (int jj=0; jj<rr_y; jj++){
									if (map.carte[rr_x*i+ii][rr_y*j+jj]){
										Hmean_pred += Habitat(rr_x*i+ii,rr_y*j+jj);
										ncount ++;
									}
								}
								if (ncount)
							       		Hmean_pred /= ncount;
								likelihood += (H_obs-Hmean_pred)*(H_obs-Hmean_pred);
							}
						}
					}			
				}
			}	
		}//end of 'sp' loop

		//------------------------------------------------------//
		//	     LIKELIHOOD FOR SPAWNING HABITAT		//
		//------------------------------------------------------//
		if (param->habitat_run_type==0){

	        int rr_x = (int)nlon/nlon_input;
        	int rr_y = (int)nlat/nlat_input;
			if (rr_x<1 || rr_y<1){ 
				cerr << "Model resolution should be divisible without remainder by the resolution of input density field." << 
				         endl<< "Currenly nlon/nlon_input = " << rr_x << ", and nlat/nlat_input = " << rr_y << endl << "Will exit now...";
				exit(1);
			}

			if (param->habitat_spawning_input_categorical_flag==0){
				for (int i=1; i<nlon_input; i++){
					for (int j=1; j<nlat_input; j++){
						double    H_obs  = mat.habitat_input[0][t_count][i][j];
						if (H_obs>0){
							int ncount = 0;
							dvariable Hmean_pred = 0.0;
							for (int ii=0; ii<rr_x; ii++)
							for (int jj=0; jj<rr_y; jj++){
								if (map.carte[rr_x*i+ii][rr_y*j+jj]){
									Hmean_pred += Habitat(rr_x*i+ii,rr_y*j+jj);
									ncount ++;
								}
							}
							if (ncount)
									Hmean_pred /= ncount;
							//if (Hmean_pred>0)
							likelihood += (H_obs-Hmean_pred)*(H_obs-Hmean_pred);
						}
					}
				}		
			}else{
				// If the habitat input is categorical
				if (month==3 || month==6 || month==9 || month ==12){
					for (int i=1; i<=nlon; i++){
						for (int j=1; j<=nlat; j++){
							if (map.carte[i][j]){
								int N_obs  = (int)mat.habitat_input[0][t_count][i][j];
								if (N_obs >= 0){
									dvariable H_pred = 0.0;
									H_pred = Habitat(i,j);

									// For categorical Poisson likelihood							
									if (H_pred == 0.0){
										nb_Npred_zero += 1;
										if (N_obs > 0){
											likelihood += likelihood_penalty;
										}
									}else{
										nb_Npred_notzero += 1;
										dvariable lkhd = NshkwCat.categorical_poisson_comp(N_obs, H_pred, weight_Nobszero, *param, 0);
										
										//double likelihood_before = value(likelihood);
										if (std::isinf(value(lkhd))){
											if (lkhd > 0){
												lkhd = likelihood_penalty;
											}else{
												lkhd = 0;
											}
										}
										likelihood += lkhd;
									}

									/*// For mixed Gaussian Kernel likelihood
									dvariable lkhd = NshkwCat.mixed_gaussian_comp(N_obs, H_pred, weight_Nobszero, *param, 0);
									likelihood += lkhd;*/
									
									if (std::isinf(value(likelihood)) && test_inf_likelihood==0){
										test_inf_likelihood += 1;
										std::cerr << "From (" << t_count << "," << i << "," << j << "): likelihood is inf" << std::endl;							
									}
								}
							}
						}
					}		
				}
			}
		}


		if (writeoutputfiles){
			if (!param->gcalc())
				HabitatConsoleOutput(t_count,1,date_str,norm(mat.habitat_input(0,t_count)),norm(value(Habitat)),value(likelihood));

			//write model habitat into a DYM file
 			dmatrix mat2d(0, nbi - 1, 0, nbj - 1);
                        mat2d.initialize();
                        for (int i=map.imin; i <= map.imax; i++){
                                for (int j=map.jinf[i] ; j<=map.jsup[i] ; j++){
                                        if (map.carte[i][j]){
                                                //mat2d(i-1,j-1) = mat.habitat_input(t_count,i,j);
                                                mat2d(i-1,j-1) = value(Habitat(i,j));//if multiple, always the last habitat is written!
                                        }
                                }
                        }
			double minval=0.0;
			double maxval=1.0;

            rw.wbin_transpomat2d(fileout, mat2d, nbi-2, nbj-2, true);
			//update min-max values in header
			rw.rwbin_minmax(fileout, minval, maxval);

		}
		past_month=month;
		step_count++;
		if (qtr != past_qtr) past_qtr = qtr; 

	} // end of simulation loop

	param->total_like = value(likelihood);
	cout << "end of forward run, likelihood: " << value(likelihood)-eFlike<< " " << eFlike <<endl;

	// to remove
	//std::cerr << "Number of null N_pred: " << nb_Npred_zero << std::endl;
	//std::cerr << "Number of non-null N_pred: " << nb_Npred_notzero << std::endl;
	//std::cerr << "Number of non-null lkhd: " << nb_lkhd << std::endl;
	//std::cerr << "Number of non-null lkhd = 0.125: " << nb_lkhd125 << std::endl;
	//std::cerr << "Number of non-null lkhd = 0.5: " << nb_lkhd5 << std::endl;
	//std::cerr << "Likelihood = " << likelihood << std::endl;
	//std::cerr << "Likelihood at timestep 3 = " << lkhd_tcount3 << std::endl;
	//std::cerr << "nb_cells with lkhd != 0 at timestep 3 = " << nb_lkhd_tcount3 << std::endl;

	return value(likelihood);
}

void SeapodymCoupled::ReadHabitat()
{
	cout << "Reading input habitat file... " << endl;
	int jday = 0;
	int t_count_init = t_count;
	//mat.createMatHabitat_input(map,param->nb_habitat_run_age,nbt_total);
	int nlevel = 0;
	string file_input;
	if (param->habitat_run_type==0)
		file_input = param->strdir_output + param->sp_name[0] + "_spawning_habitat_input.dym";
	else
		file_input = param->strdir_output + param->sp_name[0] + "_feeding_habitat_input_age1.dym";

        rw.rbin_headpar(file_input, nlon_input, nlat_input, nlevel);

	//equivalent to mat.createMatHabitat_input:
	mat.habitat_input.allocate(0,param->nb_habitat_run_age-1);
        for (int n=0; n<param->nb_habitat_run_age; n++){
		mat.habitat_input[n].allocate(1,nbt_total);
	        for (int t=1; t<=nbt_total; t++){
			mat.habitat_input[n][t].allocate(0, nlon_input, 0, nlat_input);
			mat.habitat_input[n][t].initialize();
		}
	}
	
	for (; t_count<=nbt_total; t_count++){
		getDate(jday);
		//----------------------------------------------//
		//	DATA READING SECTION: U,V,T,O2,PP	//
		//----------------------------------------------//
		if (t_count > nbt_building) {
			//TIME SERIES 
			t_series = t_count - nbt_building + nbt_start_series;

			//----------------------------------------------//
			//	READING HABITAT DATA			//
			//----------------------------------------------//

			//will need to prepare the input habitat that contains the data only for the selected time period
			int nbytetoskip = (9 +(3* nlat_input * nlon_input) + nlevel + ((nlat_input *nlon_input)* (t_count-1))) * 4;
			if (param->habitat_run_type==0){
				if (param->habitat_spawning_input_categorical_flag==0){
					file_input = param->strdir_output + param->sp_name[0] + "_spawning_habitat_input.dym";
				}else if (param->habitat_spawning_input_categorical_flag==1){
					file_input = param->str_file_larvae;
				}
				//rw.rbin_input2d(file_input, map, mat.habitat_input[0][t_count], nlon_input+2, nlat_input+2, nbytetoskip);
				ifstream litbin(file_input.c_str(), ios::binary | ios::in);
				if (!litbin){
					cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_input << "\"\n";
					exit(1);
				}

				if (param->habitat_spawning_input_categorical_flag==0){
					litbin.seekg(nbytetoskip, ios::cur);
				}else if (param->habitat_spawning_input_categorical_flag==1){
					int qtr = ((t_count - 1) % 12) / 3 + 1;
					int nbytetoskip_seasonalfile = (9 +(3* nlat * nlon) + 4 + ((nlat *nlon)* (qtr-1))) * 4;
					litbin.seekg(nbytetoskip_seasonalfile, ios::cur);
				}
				const int sizeofDymInputType = sizeof(float);
				float buf;
				for (int j=0;j<nlat_input;j++){
					for (int i=0;i<nlon_input;i++){
						litbin.read(( char *)&buf,sizeofDymInputType);
						mat.habitat_input[0][t_count][i+1][j+1]= buf;
					}
				}
				litbin.close();
			} else {
				for (int n=0; n<param->nb_habitat_run_age; n++){
					std::ostringstream ostr;
               				ostr << n+1;
					file_input = param->strdir_output + param->sp_name[0] + 
						"_feeding_habitat_input_age" + ostr.str() + ".dym";
					//rw.rbin_input2d(file_input, map, mat.habitat_input[n][t_count], nbi, nbj, nbytetoskip);
					ifstream litbin(file_input.c_str(), ios::binary | ios::in);
					if (!litbin){
						cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_input << "\"\n";
						exit(1);
					}
					litbin.seekg(nbytetoskip, ios::cur);
					const int sizeofDymInputType = sizeof(float);
		            float buf;
		            for (int j=0;j<nlat_input;j++){
						for (int i=0;i<nlon_input;i++){
							litbin.read(( char *)&buf,sizeofDymInputType);
							mat.habitat_input[n][t_count][i+1][j+1]= buf;
						}
					}
					litbin.close();
				}
			}
		}
	}
	t_count = t_count_init;
}


void HabitatConsoleOutput(int t_count, int flag_simulation, string date_str, double norm_habitat_input, double norm_Habitat, double like)
{
        if (flag_simulation){
                cout << setw(4)  << left << t_count<<"| "
                     << setw(14) << date_str<<"| "
                     << setw(14) << norm_habitat_input <<" | "
                     << setw(14) << norm_Habitat<<" | "
                     << like << endl;
        }
        else{
                //CalcSums();
                // Comment: Larvae biomass is the result of spawning function S(t-1) and ADRE at time t
                // we have Larvae(t=1)=0 because S(t-1) is computed after CalcSums()
                cout << endl << "Model temporal dynamics:" << endl;
                cout << setw(4)  << left << "step"<<"| "
                     << setw(14) << "Date"<<"| "
                     << setw(14) << "Input Habitat"<<" | "
                     << setw(14) << "Model_Habitat"<<" | "
                     << setw(14) << "Likelihood"<< endl;

                cout << "----------------------------------------------------------------------------------------------------------"<< endl;
	}
}

