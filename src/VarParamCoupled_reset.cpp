#include "VarParamCoupled.h"
void dv_Bard_pen();
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
string get_path(const char* parfile);

dvariable VarParamCoupled::reset(dvar_vector x)
{
	dvariable penalty = 0.0;
	int idx = 1;

	if (!scalc())
		cout << "total penalty: ";

//Order is important (must be the same as in xinit function)!
	if (doc.get("/Mp_mean_max/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsMp_mean_max[i] = boundp(x[idx], Mp_mean_max_min, Mp_mean_max_max, penalty);
			//large_pen(dvarsMp_mean_max[i],Mp_mean_max_min, Mp_mean_max_max,penalty,0.01);
			//Hint: last constant - critical distance to the boundary for linear scaling, i.e.
			//in order to get explicite distance multiply this constant by abs(xmax-xmin)
			Mp_mean_max[i] = value(dvarsMp_mean_max[i]);
			dvarpars[idx] = Mp_mean_max[i];
			++idx;
		}
	}

	if (doc.get("/Mp_mean_exp/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsMp_mean_exp[i] = boundp(x[idx], Mp_mean_exp_min, Mp_mean_exp_max, penalty);
			//large_pen(dvarsMp_mean_exp[i],Mp_mean_exp_min, Mp_mean_exp_max,penalty,0.01);
			Mp_mean_exp[i] = value(dvarsMp_mean_exp[i]);
			dvarpars[idx] = Mp_mean_exp[i];
			++idx;
		}
	}

	if (doc.get("/Ms_mean_max/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsMs_mean_max[i] = boundp(x[idx], Ms_mean_max_min, Ms_mean_max_max, penalty);
			//large_pen(dvarsMs_mean_max[i],Ms_mean_max_min, Ms_mean_max_max,penalty,0.01);

			Ms_mean_max[i] = value(dvarsMs_mean_max[i]);
			dvarpars[idx] = Ms_mean_max[i];
			++idx;
		}
	}

	if (doc.get("/Ms_mean_slope/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsMs_mean_slope[i] = boundp(x[idx], Ms_mean_slope_min, Ms_mean_slope_max, penalty);
			Ms_mean_slope[i] = value(dvarsMs_mean_slope[i]);
			dvarpars[idx] = Ms_mean_slope[i];
			++idx;
		}
	}

	if (doc.get("/M_mean_range/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsM_mean_range[i] = boundp(x[idx], M_mean_range_min, M_mean_range_max, penalty);
			//large_pen(dvarsM_mean_range[i],M_mean_range_min,M_mean_range_max,penalty,0.01);
			M_mean_range[i] = value(dvarsM_mean_range[i]);
			dvarpars[idx] = M_mean_range[i];
			++idx;
		}
	}

	if (doc.get("/a_sst_spawning/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsA_sst_spawning[i] = boundp(x[idx], a_sst_spawning_min, a_sst_spawning_max, penalty);
			//large_pen(dvarsA_sst_spawning[i],a_sst_spawning_min,a_sst_spawning_max,penalty,0.005);

			a_sst_spawning[i] = value(dvarsA_sst_spawning[i]);
			dvarpars[idx] = a_sst_spawning[i];
			++idx;
		}
	}

	if (doc.get("/b_sst_spawning/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {			
			dvarsB_sst_spawning[i] = boundp(x[idx], b_sst_spawning_min, b_sst_spawning_max, penalty);
			//large_pen(dvarsB_sst_spawning[i],b_sst_spawning_min,b_sst_spawning_max,penalty,0.05);

			b_sst_spawning[i] = value(dvarsB_sst_spawning[i]);
			dvarpars[idx] = b_sst_spawning[i];
			++idx;
		}
	}

	if (doc.get("/q_sp_larvae/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {			
			dvarsQ_sp_larvae[i] = boundp(x[idx], q_sp_larvae_min, q_sp_larvae_max, penalty);

			q_sp_larvae[i] = value(dvarsQ_sp_larvae[i]);
			dvarpars[idx] = q_sp_larvae[i];
			++idx;
		}
	}

	if (doc.get("/likelihood_spawning_sigma/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {			
			dvarsLikelihood_spawning_sigma[i] = boundp(x[idx], likelihood_spawning_sigma_min, likelihood_spawning_sigma_max, penalty);

			likelihood_spawning_sigma[i] = value(dvarsLikelihood_spawning_sigma[i]);
			dvarpars[idx] = likelihood_spawning_sigma[i];
			++idx;
		}
	}

	if (doc.get("/likelihood_spawning_beta/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {			
			dvarsLikelihood_spawning_beta[i] = boundp(x[idx], likelihood_spawning_beta_min, likelihood_spawning_beta_max, penalty);

			likelihood_spawning_beta[i] = value(dvarsLikelihood_spawning_beta[i]);
			dvarpars[idx] = likelihood_spawning_beta[i];
			++idx;
		}
	}

	if (doc.get("/likelihood_spawning_probzero/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {			
			dvarsLikelihood_spawning_probzero[i] = boundp(x[idx], likelihood_spawning_probzero_min, likelihood_spawning_probzero_max, penalty);

			likelihood_spawning_probzero[i] = value(dvarsLikelihood_spawning_probzero[i]);
			dvarpars[idx] = likelihood_spawning_probzero[i];
			++idx;
		}
	}

	if (doc.get("/a_sst_larvae/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {			
			dvarsA_sst_larvae[i] = boundp(x[idx], a_sst_larvae_min, a_sst_larvae_max, penalty);

			a_sst_larvae[i] = value(dvarsA_sst_larvae[i]);
			dvarpars[idx] = a_sst_larvae[i];
			++idx;
		}
	}

	if (doc.get("/b_sst_larvae/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {			
			dvarsB_sst_larvae[i] = boundp(x[idx], b_sst_larvae_min, b_sst_larvae_max, penalty);

			b_sst_larvae[i] = value(dvarsB_sst_larvae[i]);
			dvarpars[idx] = b_sst_larvae[i];
			++idx;
		}
	}

	if (doc.get("/alpha_hsp_prey/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsAlpha_hsp_prey[i] = boundp(x[idx], alpha_hsp_prey_min, alpha_hsp_prey_max, penalty);
			alpha_hsp_prey[i] = value(dvarsAlpha_hsp_prey[i]);
			dvarpars[idx] = alpha_hsp_prey[i];
			++idx;
		}
	}

	if (doc.get("/alpha_hsp_predator/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsAlpha_hsp_predator[i] = boundp(x[idx], alpha_hsp_predator_min, alpha_hsp_predator_max, penalty);
			alpha_hsp_predator[i] = value(dvarsAlpha_hsp_predator[i]);
			dvarpars[idx] = alpha_hsp_predator[i];
			++idx;
		}
	}

	if (doc.get("/beta_hsp_predator/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsBeta_hsp_predator[i] = boundp(x[idx], beta_hsp_predator_min, beta_hsp_predator_max, penalty);
			beta_hsp_predator[i] = value(dvarsBeta_hsp_predator[i]);
			dvarpars[idx] = beta_hsp_predator[i];
			++idx;
		}
	}	

	if (doc.get("/a_sst_habitat/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsA_sst_habitat[i] = boundp(x[idx], a_sst_habitat_min, a_sst_habitat_max, penalty);
			//large_pen(dvarsA_sst_habitat[i],a_sst_habitat_min,a_sst_habitat_max,penalty,0.01);

			a_sst_habitat[i] = value(dvarsA_sst_habitat[i]);
			dvarpars[idx] = a_sst_habitat[i];
			++idx;
//dvarsA_sst_spawning[i] = dvarsA_sst_habitat[i];
		}
	}


	if (doc.get("/b_sst_habitat/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsB_sst_habitat[i] = boundp(x[idx], b_sst_habitat_min, b_sst_habitat_max, penalty);
			//large_pen(dvarsB_sst_habitat[i],b_sst_habitat_min,b_sst_habitat_max,penalty,0.1);

			b_sst_habitat[i] = value(dvarsB_sst_habitat[i]);
			dvarpars[idx] = b_sst_habitat[i];
			++idx;
		}
	}
	if (doc.get("/T_age_size_slope/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsT_age_size_slope[i] = boundp(x[idx], T_age_size_slope_min, T_age_size_slope_max, penalty);

			T_age_size_slope[i] = value(dvarsT_age_size_slope[i]);
			dvarpars[idx] = T_age_size_slope[i];
			++idx;
		}
	}
	for (int n=0; n<3; n++){
		std::ostringstream ostr;
		ostr << "delta" << n+1;
		string s = "/thermal_func/"+ostr.str()+"/variable";
		if (doc.get(s,"use") == "true"){
			for (int i = 0; i < nb_species; i++) {					
				dvarsThermal_func_delta[n][i] = boundp(x[idx], thermal_func_delta_min[n], 
									thermal_func_delta_max[n], penalty);
				thermal_func_delta[n][i] = value(dvarsThermal_func_delta[n][i]);
				dvarpars[idx] = thermal_func_delta[n][i];
				++idx;				
			}
		}
	}
	if (doc.get("/a_oxy_habitat/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsA_oxy_habitat[i] = boundp(x[idx], a_oxy_habitat_min, a_oxy_habitat_max, penalty);
			//large_pen(dvarsA_oxy_habitat[i],a_oxy_habitat_min,a_oxy_habitat_max,penalty,0.1);

			a_oxy_habitat[i] = value(dvarsA_oxy_habitat[i]);
			dvarpars[idx] = a_oxy_habitat[i];
			++idx;
		}
	}
	if (doc.get("/b_oxy_habitat/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsB_oxy_habitat[i] = boundp(x[idx], b_oxy_habitat_min, b_oxy_habitat_max, penalty);
			//large_pen(dvarsB_oxy_habitat[i],b_oxy_habitat_min,b_oxy_habitat_max,penalty,0.01);

			b_oxy_habitat[i] = value(dvarsB_oxy_habitat[i]);
			dvarpars[idx] = b_oxy_habitat[i];
			++idx;
		}
	}   
	for (int n=0; n<nb_forage; n++){
		if (doc.get(string("/eF_habitat/")+frg_name[n]+"/variable", "use") == "true") {
			for (int i = 0; i < nb_species; i++) {
				dvarsEF_habitat[n][i] = boundp(x[idx], eF_habitat_min[n], eF_habitat_max[n], penalty);

				eF_habitat[n][i] = value(dvarsEF_habitat[n][i]);
				dvarpars[idx] = eF_habitat[n][i];
				++idx;
			}
		}
	}
	if (doc.get("/hp_cannibalism/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsHp_cannibalism[i] = boundp(x[idx], hp_cannibalism_min, hp_cannibalism_max, penalty);

			hp_cannibalism[i] = value(dvarsHp_cannibalism[i]);
			dvarpars[idx] = hp_cannibalism[i];
			++idx;
		}
	}

	if (doc.get("/sigma_species/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsSigma_species[i] = boundp(x[idx], sigma_species_min, sigma_species_max, penalty);
			//large_pen(dvarsSigma_species[i],sigma_species_min,sigma_species_max,penalty,0.1);

			sigma_species[i] = value(dvarsSigma_species[i]);
			dvarpars[idx] = sigma_species[i];
			++idx;
		}
	}

	if (doc.get("/MSS_species/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsMSS_species[i] = boundp(x[idx], MSS_species_min, MSS_species_max, penalty);
			//large_pen(dvarsMSS_species[i],MSS_species_min,MSS_species_max,penalty,0.1);

			MSS_species[i] = value(dvarsMSS_species[i]);
			dvarpars[idx] = MSS_species[i];
			++idx;
		}
	}

	if (doc.get("/MSS_size_slope/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsMSS_size_slope[i] = boundp(x[idx], MSS_size_slope_min, MSS_size_slope_max, penalty);

			MSS_size_slope[i] = value(dvarsMSS_size_slope[i]);
			dvarpars[idx] = MSS_size_slope[i];
			++idx;
		}
	}
	
	if (doc.get("/c_diff_fish/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsC_diff_fish[i] = boundp(x[idx], c_diff_fish_min, c_diff_fish_max, penalty);
			//large_pen(dvarsC_diff_fish[i],c_diff_fish_min,c_diff_fish_max,penalty,0.01);

			c_diff_fish[i] = value(dvarsC_diff_fish[i]);
			dvarpars[idx] = c_diff_fish[i];
			++idx;
		}
	} 

	if (doc.get("/nb_recruitment/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsNb_recruitment[i] = boundp(x[idx], nb_recruitment_min, nb_recruitment_max, penalty);
			//large_pen(dvarsNb_recruitment[i],nb_recruitment_min,nb_recruitment_max,penalty,0.01);

			nb_recruitment[i] = value(dvarsNb_recruitment[i]);
			dvarpars[idx] = nb_recruitment[i];
			++idx;
		}
	} 

	if (doc.get("/a_adults_spawning/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsA_adults_spawning[i] = boundp(x[idx], a_adults_spawning_min, a_adults_spawning_max, penalty);
			//large_pen(dvarsA_adults_spawning[i],a_adults_spawning_min,a_adults_spawning_max,penalty,0.01);

			a_adults_spawning[i] = value(dvarsA_adults_spawning[i]);
			dvarpars[idx] = a_adults_spawning[i];
			++idx;
		}
	} 

	if (doc.get("/spawning_season_peak/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsSpawning_season_peak[i] = boundp(x[idx], spawning_season_peak_min, spawning_season_peak_max, penalty);

			spawning_season_peak[i] = value(dvarsSpawning_season_peak[i]);
			dvarpars[idx] = spawning_season_peak[i];
			++idx;
		}
	} 
	if (doc.get("/spawning_season_start/variable", "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			dvarsSpawning_season_start[i] = boundp(x[idx], spawning_season_start_min, spawning_season_start_max, penalty);

			spawning_season_start[i] = value(dvarsSpawning_season_start[i]);
			dvarpars[idx] = spawning_season_start[i];
			++idx;
		}
	} 
	
	if (doc.get("/q_sp_fishery/variables", "use") == "true") {
		for (int i=0; i<nb_species; i++) {
			int ifx = 0; 
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					if (doc.get("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true"){
						dvarsQ_sp_fishery[i][ifx] = boundp(x[idx], q_sp_fishery_min[i][ifx], q_sp_fishery_max[i][ifx], penalty);
						//large_pen(dvarsQ_sp_fishery[i][ifx],q_sp_fishery_min[i][ifx],q_sp_fishery_max[i][ifx],penalty,0.005);

						q_sp_fishery[i][ifx] = value(dvarsQ_sp_fishery[i][ifx]);
						dvarpars[idx] = q_sp_fishery[i][ifx];
						++idx; 
					} ifx++;
				}
			}
		}
	}
	if (doc.get("/s_sp_fishery/variables", "use") == "true") {
		for (int i=0; i<nb_species; i++) {
			int ifx = 0; 
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true"){
						dvarsSslope_sp_fishery[i][ifx] = boundp(x[idx], s_slope_sp_fishery_min[i][ifx], s_slope_sp_fishery_max[i][ifx], penalty);
						//large_pen(dvarsSslope_sp_fishery[i][ifx],s_slope_sp_fishery_min[i][ifx],s_slope_sp_fishery_max[i][ifx],penalty,0.01);

						s_slope_sp_fishery[i][ifx] = value(dvarsSslope_sp_fishery[i][ifx]);
						dvarpars[idx] = s_slope_sp_fishery[i][ifx];
						++idx;
					}
					if (s_func_type[f]>1)
						if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/length_threshold", "use") == "true"){ 
							dvarsSlength_sp_fishery[i][ifx] = boundp(x[idx], dvarpars_min[idx], dvarpars_max[idx], penalty);
							s_length_sp_fishery[i][ifx] = value(dvarsSlength_sp_fishery[i][ifx]);
							dvarpars[idx] = s_length_sp_fishery[i][ifx];
							++idx; 
						} 
					if (s_func_type[f]>=3)
						if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote", "use") == "true"){ 
							dvarsSasympt_sp_fishery[i][ifx] = boundp(x[idx], s_asympt_sp_fishery_min[i][ifx], s_asympt_sp_fishery_max[i][ifx], penalty);
							s_asympt_sp_fishery[i][ifx] = value(dvarsSasympt_sp_fishery[i][ifx]);
							dvarpars[idx] = s_asympt_sp_fishery[i][ifx];
							++idx; 
						} 
					ifx++;
				}
			}
		}
	}

	if (doc.get("/likelihood_parameters/variables", "use") == "true") {
		for (int i=0; i<nb_species; i++) {
			int ifx = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					//if (like_types[i][ifx]==2 || like_types[i][ifx]==4 || like_types[i][ifx]==5){
					if (/*like_types[i][ifx]==2 || */like_types[i][ifx]==4 || like_types[i][ifx]==5){
					dvarsLike_param[i][ifx] = boundp(x[idx], dvarpars_min[idx], dvarpars_max[idx], penalty);
					like_param[i][ifx] = value(dvarsLike_param[i][ifx]);
					dvarpars[idx] = like_param[i][ifx];
					++idx;
					} ifx++;
				}
			}
		}
	}

	for (int i=0; i<nb_species; i++) {
		int ifx = 0;
		for (int f = 0; f < nb_fishery; f++) {
			if (mask_fishery_sp[i][f]){
				if (like_types[i][ifx]==5 && doc.get("/likelihood_parameters/variables", "use") == "true"){ 

					dvarsProb_zero[i][ifx] = boundp(x[idx], dvarpars_min[idx], dvarpars_max[idx], penalty);
					prob_zero[i][ifx] = value(dvarsProb_zero[i][ifx]);
					dvarpars[idx] = prob_zero[i][ifx];
					++idx; 
				} ifx++;
			}
		}
	}
//cout << dvarpars << endl;
//cout << statpars << endl;
//exit(1);
	if (!scalc())
		cout << value(penalty) << "; ";
	return penalty;
}

 //partially autodiff function
void VarParamCoupled::large_pen(dvariable x, double xmin, double xmax, dvariable& penalty, const double eps)
{
/*
	const double nbt = (365.25/deltaT) * (save_last_yr - save_first_yr);
	const double diff = sqrt(pow(xmax-xmin,2));
	const double pen  = 10.0*nbt/diff;
	dvariable xs = (x-min(xmin,xmax))/diff;

	if (xs<1e-7) xs = 1e-7;
	if (xs>0.9999999) xs = 0.9999999;

	if (xs>1-eps || xs<eps)
		penalty += -pen*(log(xs+1e-40)+log(1.0-xs+1.e-40)+log(4.0));

*/

	double x_c = value(x);
	dvariable Penalty = 0.0;

	const double nbt = (365.25/deltaT) * (save_last_yr - save_first_yr);
	const double diff = sqrt(pow(xmax-xmin,2));
	const double pen  = 4.0*nbt/diff;
	double xs = (x_c-min(xmin,xmax))/diff;

	if (xs>1-eps || xs<eps){
		cout << x << " ";
		Penalty = Bard_pen(xs,xmin,xmax,pen);
	}


	save_identifier_string2((char*)"Pen_comp_begin");
	x.save_prevariable_value();
	x.save_prevariable_position();
	Penalty.save_prevariable_position();
	save_double_value(eps);
	save_double_value(pen);
	save_double_value(xmin);
	save_double_value(xmax);
	save_identifier_string2((char*)"Pen_comp_end");
	
	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Bard_pen);

	//autodiff
	penalty += Penalty;
}

double VarParamCoupled::Bard_pen(double x, double xmin, double xmax, double pen)
{
	if (x<1e-7) x = 1e-7;
	if (x>0.9999999) x = 0.9999999;

	double penalty = -pen*(log(x+1e-40)+log(1.0-x+1.e-40)+log(4.0));

	return penalty;
}

void dv_Bard_pen()
{
	verify_identifier_string2((char*)"Pen_comp_end");	
	const double xmax = restore_double_value();
	const double xmin = restore_double_value();
	const double P    = restore_double_value();
	const double s    = restore_double_value();
	const prevariable_position pen_pos = restore_prevariable_position();
	const prevariable_position x_pos = restore_prevariable_position();
	const double x    = restore_prevariable_value();
	verify_identifier_string2((char*)"Pen_comp_begin");
		
	double dfx = restore_prevariable_derivative(x_pos);
	double dfpen = restore_prevariable_derivative(pen_pos);

	//recompute xs and P
	const double diff = sqrt(pow(xmax-xmin,2));
	double xs = (x-min(xmin,xmax))/diff;

	const double eps = 1e-40;

	if (xs>1-s || xs<s) {
		//double penalty = -pen*(log(xs+1e-40)+log(1.0-xs+1.e-40)+log(4.0));
		double dfxs = -P*(1.0/(xs+eps)-1.0/(1.0-xs+eps))*dfpen;

		if (xs<1e-7) //xs = 1e-7;
			dfxs = 0.0;
		
		if (xs>0.9999999) //xs = 0.9999999;
			dfxs = 0.0;


		//double xs = (x_c-min(xmin,xmax))/diff;
		dfx += (1.0/diff)*dfxs;
	}
	
	save_double_derivative(dfx,x_pos);
	save_double_derivative(dfpen,pen_pos);
}
/*
//adjoint
void VarParamCoupled::large_pen(dvariable x, double xmin, double xmax, dvariable& penalty, const double eps)
{
	double x_c = value(x);

	save_identifier_string2("before_Pen_comp");
	penalty.save_prevariable_position();

	const double nbt = (365.25/deltaT) * (save_last_yr - save_first_yr);
	const double diff = sqrt(pow(xmax-xmin,2));
	const double pen  = 10.0*nbt/diff;
	double xs = (x_c-min(xmin,xmax))/diff;

	if (xs>1-eps || xs<eps)
		penalty += Bard_pen(xs,xmin,xmax,pen);

	save_identifier_string2("Pen_comp_begin");
	x.save_prevariable_value();
	x.save_prevariable_position();
	penalty.save_prevariable_position();
	save_double_value(eps);
	save_double_value(pen);
	save_double_value(xmin);
	save_double_value(xmax);
	save_identifier_string2("Pen_comp_end");

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Bard_pen);

	
}

double VarParamCoupled::Bard_pen(double x, double xmin, double xmax, double pen)
{
	if (x<1e-7) x = 1e-7;
	if (x>0.9999999) x = 0.9999999;

	double penalty = -pen*(log(x+1e-40)+log(1.0-x+1.e-40)+log(4.0));

	return penalty;
}

void dv_Bard_pen()
{
	verify_identifier_string2("Pen_comp_end");	
	const double xmax = restore_double_value();
	const double xmin = restore_double_value();
	const double P    = restore_double_value();
	const double s    = restore_double_value();
	const prevariable_position pen_pos = restore_prevariable_position();
	const prevariable_position x_pos = restore_prevariable_position();
	const double x    = restore_prevariable_value();
	verify_identifier_string2("Pen_comp_begin");
		
	const prevariable_position pen_pr_pos = restore_prevariable_position();
	verify_identifier_string2("before_Pen_comp");

	double dfx = restore_prevariable_derivative(x_pos);
	double dfpen = restore_prevariable_derivative(pen_pos);
	double dfpen_pr = restore_prevariable_derivative(pen_pr_pos);

	double dfB = 0.0;

	//recompute xs and P
	const double diff = sqrt(pow(xmax-xmin,2));
	double xs = (x-min(xmin,xmax))/diff;

	const double eps = 1e-40;

	if (xs>1-s || xs<s) {
		//penalty += Bard_Penalty;
		dfpen_pr += dfpen;
		dfB += dfpen;
		dfpen = 0.0;

		//double penalty = -pen*(log(xs+1e-40)+log(1.0-xs+1.e-40)+log(4.0));
		double dfxs = -P*(1/(xs+eps)-1/(1-xs+eps))*dfB;

		if (xs<1e-7) //xs = 1e-7;
			dfxs = 0.0;
		
		if (xs>0.9999999) //xs = 0.9999999;
			dfxs = 0.0;

		//double xs = (x_c-min(xmin,xmax))/diff;
		dfx += (1.0/diff)*dfxs;
	}
	
	save_double_derivative(dfx,x_pos);
	save_double_derivative(dfpen,pen_pos);
	save_double_derivative(dfpen_pr,pen_pr_pos);
}
*/
void VarParamCoupled::getparam()
{	
	int idx = 1;
	while (idx <=_nvarcalc){
		for (int i = 0; i < nb_species; i++){ 
			doc.set(parfile_names[idx-1], sp_name[i], dvarpars[idx]);
			idx++;
		}
	}
	if (!doc.get("/total_likelihood","value").empty())
		doc.set("/total_likelihood","value", total_like);		
}

void VarParamCoupled::outp_param(adstring_array x_names, const int nvars)
{
	cout << '\n';
	cout << "Control variables in optimisation:" << endl;
	for (int idx=1; idx<=nvars; idx++)
		cout << idx << ". " << x_names[idx] << " = " << dvarpars[idx] << endl;
}

void VarParamCoupled::save_statistics(const string dirout, const adstring_array x_names, double likelihood, dvector g, double elapsed_time, int status, int iter, int nvars)
{
	ofstream ofs;
	double gmax;
	string filename = dirout + "/optim.rep";
//	const char* filename = ;

	if (status >= 1){
		ofs.open(filename.c_str(), ios::app);
	} else {
		string like_type = "";
		for (int sp=0; sp<nb_species; sp++)
			for (int k=0; k<nb_fishery_by_sp[sp]; k++)
				like_type += " " + str(like_types[sp][k]);
		ofs.open(filename.c_str(), ios::out);
		string text = "not included, ";
		if (tuna_spinup!=0) text = "included, period of population building: " + str((int)(save_first_yr-4)) + "-" + str(int(save_first_yr))+ "; ";
		ofs << "Simulation: domain " << latitudeMin << ", " << latitudeMax << ", " << longitudeMin << ", " << longitudeMax << 
			"; spinup " << text << "skipped after spinup or initial conditions: " << nbsteptoskip << 
			"; time period with likelihood optimization: " << save_first_yr+nbsteptoskip/12<< "-" << save_last_yr << endl;
		ofs << "Name";
		for (int idx=0; idx<_nstatpars; idx++)
			ofs << "\t" << statpar_names[idx];
		ofs << endl;
		ofs << "Value";
		for (int idx=0; idx<_nstatpars; idx++)
			ofs << "\t" << statpars[idx];
		ofs << endl;

		ofs << "Iteration\t" << "$N_{FE}$\t" << "time\t" << "L (types" << like_type << ")\t" << "$G_{max}$\t" << "||G||";
		int ncoltoskip=6;
		for (int idx=1; idx<=nvars; idx++)
			ofs << "\t" << x_names[idx];
		ofs << endl;
		for (int i=0;i<ncoltoskip;i++) ofs << "\t0";
		for (int idx=1; idx<=nvars; idx++)
			ofs << "\t" << dvarpars_min[idx];		
		ofs << endl;
		for (int i=0;i<ncoltoskip;i++) ofs << "\t0";
		for (int idx=1; idx<=nvars; idx++)
			ofs << "\t" << dvarpars_max[idx];		
		ofs << endl;	
	}
	gmax = max(fabs(g));
	ofs << iter << "\t" << status << "\t" << elapsed_time <<"\t"<< likelihood << "\t"<< gmax<< "\t" <<norm(g);
	for (int idx=1; idx<=nvars; idx++)
		ofs << "\t" << dvarpars[idx];

	ofs << endl;
	ofs.close();

}

