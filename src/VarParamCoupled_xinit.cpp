#include "VarParamCoupled.h"
bool nocase_compare(char c1, char c2);

void VarParamCoupled::xinit(dvector& x, adstring_array& x_names)
{
	int idx = 1;
	x.initialize();

	par_init(dvarsMp_mean_max,Mp_mean_max,Mp_mean_max_min,Mp_mean_max_max,"/Mp_mean_max",x,x_names,idx);
	par_init(dvarsMp_mean_exp,Mp_mean_exp,Mp_mean_exp_min,Mp_mean_exp_max,"/Mp_mean_exp",x,x_names,idx);
	par_init(dvarsMs_mean_max,Ms_mean_max,Ms_mean_max_min,Ms_mean_max_max,"/Ms_mean_max",x,x_names,idx);
	par_init(dvarsMs_mean_slope,Ms_mean_slope,Ms_mean_slope_min,Ms_mean_slope_max,"/Ms_mean_slope",x,x_names,idx);
	par_init(dvarsM_mean_range,M_mean_range,M_mean_range_min,M_mean_range_max,"/M_mean_range",x,x_names,idx);

	par_init(dvarsA_sst_spawning,a_sst_spawning,a_sst_spawning_min,a_sst_spawning_max,"/a_sst_spawning",x,x_names,idx);
	par_init(dvarsB_sst_spawning,b_sst_spawning,b_sst_spawning_min,b_sst_spawning_max,"/b_sst_spawning",x,x_names,idx);
	par_init(dvarsA_sst_larvae,a_sst_larvae,a_sst_larvae_min,a_sst_larvae_max,"/a_sst_larvae",x,x_names,idx);
	par_init(dvarsB_sst_larvae,b_sst_larvae,b_sst_larvae_min,b_sst_larvae_max,"/b_sst_larvae",x,x_names,idx);
	par_init(dvarsAlpha_hsp_prey,alpha_hsp_prey,alpha_hsp_prey_min,alpha_hsp_prey_max,"/alpha_hsp_prey",x,x_names,idx);
	par_init(dvarsAlpha_hsp_predator,alpha_hsp_predator,alpha_hsp_predator_min,alpha_hsp_predator_max,"/alpha_hsp_predator",x,x_names,idx);
	par_init(dvarsBeta_hsp_predator,beta_hsp_predator,beta_hsp_predator_min,beta_hsp_predator_max,"/beta_hsp_predator",x,x_names,idx);

	par_init(dvarsA_sst_habitat,a_sst_habitat,a_sst_habitat_min,a_sst_habitat_max,"/a_sst_habitat",x,x_names,idx);
	par_init(dvarsB_sst_habitat,b_sst_habitat,b_sst_habitat_min,b_sst_habitat_max,"/b_sst_habitat",x,x_names,idx);
	par_init(dvarsT_age_size_slope,T_age_size_slope,T_age_size_slope_min,T_age_size_slope_max,"/T_age_size_slope",x,x_names,idx);
	dvarsThermal_func_delta.allocate(0,2);
	for (int n=0; n<3; n++){
		std::ostringstream ostr;
		ostr << "delta" << n+1;
		par2_init(dvarsThermal_func_delta[n],thermal_func_delta[n],thermal_func_delta_min[n],thermal_func_delta_max[n],"/thermal_func",ostr.str(),x,x_names,idx);
	}
	par_init(dvarsA_oxy_habitat,a_oxy_habitat,a_oxy_habitat_min,a_oxy_habitat_max,"/a_oxy_habitat",x,x_names,idx);
	par_init(dvarsB_oxy_habitat,b_oxy_habitat,b_oxy_habitat_min,b_oxy_habitat_max,"/b_oxy_habitat",x,x_names,idx);
	dvarsEF_habitat.allocate(0, nb_forage - 1);
	for (int n=0; n<nb_forage; n++){
		par2_init(dvarsEF_habitat[n],eF_habitat[n],eF_habitat_min[n],eF_habitat_max[n],"/eF_habitat",frg_name[n],x,x_names,idx);
	}
	par_init(dvarsHp_cannibalism,hp_cannibalism,hp_cannibalism_min,hp_cannibalism_max,"/hp_cannibalism",x,x_names,idx);

	par_init(dvarsSigma_species,sigma_species,sigma_species_min,sigma_species_max,"/sigma_species",x,x_names,idx);
	par_init(dvarsMSS_species,MSS_species,MSS_species_min,MSS_species_max,"/MSS_species",x,x_names,idx);
	par_init(dvarsMSS_size_slope,MSS_size_slope,MSS_size_slope_min,MSS_size_slope_max,"/MSS_size_slope",x,x_names,idx);
	par_init(dvarsC_diff_fish,c_diff_fish,c_diff_fish_min,c_diff_fish_max,"/c_diff_fish",x,x_names,idx);

	par_init(dvarsNb_recruitment,nb_recruitment,nb_recruitment_min,nb_recruitment_max,"/nb_recruitment",x,x_names,idx);
	par_init(dvarsA_adults_spawning,a_adults_spawning,a_adults_spawning_min,a_adults_spawning_max,"/a_adults_spawning",x,x_names,idx);

	par_init(dvarsSpawning_season_peak,spawning_season_peak,spawning_season_peak_min,spawning_season_peak_max,"/spawning_season_peak",x,x_names,idx);
	par_init(dvarsSpawning_season_start,spawning_season_start,spawning_season_start_min,spawning_season_start_max,"/spawning_season_start",x,x_names,idx);

	dvarsQ_sp_fishery.allocate(0, nb_species - 1);
	dvarsQ_sp_fishery.initialize();
	if (doc.get("/q_sp_fishery/variables", "use") == "true") {
		for (int i=0; i<nb_species; i++) {
			const int fmax = nb_fishery_by_sp[i];
			dvarsQ_sp_fishery[i].allocate(0, fmax - 1);
			dvarsQ_sp_fishery[i].initialize();
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					if (doc.get("/q_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true"){
						const double value = q_sp_fishery[i][k];
						x[idx] = boundpin(value, q_sp_fishery_min[i][k], q_sp_fishery_max[i][k]);
						x_names[idx] = "q(" + str(i) + "," + str(k) + ")";
						dvarpars[idx] = value;
						dvarpars_min[idx] = q_sp_fishery_min[i][k];
						dvarpars_max[idx] = q_sp_fishery_max[i][k];
						parfile_names[idx-1] = "/q_sp_fishery/"+list_fishery_name[f];
						++idx; 
					}
					else dvarsQ_sp_fishery[i][k] = q_sp_fishery[i][k]; 
					k++;
				}
			}
		}
	}else {
		for (int i=0; i<nb_species; i++) {
			dvarsQ_sp_fishery[i].allocate(0, nb_fishery_by_sp[i] - 1);
			dvarsQ_sp_fishery[i].initialize();
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					dvarsQ_sp_fishery[i][k] = q_sp_fishery[i][k];
					k++;
				}
			}
		}
	}
	dvarsSslope_sp_fishery.allocate(0, nb_species - 1);
	dvarsSlength_sp_fishery.allocate(0, nb_species - 1);
	dvarsSasympt_sp_fishery.allocate(0, nb_species - 1);
	dvarsSslope_sp_fishery.initialize();
	dvarsSlength_sp_fishery.initialize();
	dvarsSasympt_sp_fishery.initialize();
	for (int i=0; i<nb_species; i++) {
		const int fmax = nb_fishery_by_sp[i];
		dvarsSslope_sp_fishery[i].allocate(0, fmax - 1);
		dvarsSlength_sp_fishery[i].allocate(0, fmax - 1);
		dvarsSasympt_sp_fishery[i].allocate(0, fmax - 1);
		dvarsSslope_sp_fishery[i].initialize();
		dvarsSasympt_sp_fishery[i].initialize();

		if (doc.get("/s_sp_fishery/variables", "use") == "true") {
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/variable", "use") == "true"){
						double value = s_slope_sp_fishery[i][k];
						x[idx] = boundpin(value, s_slope_sp_fishery_min[i][k], s_slope_sp_fishery_max[i][k]);
						x_names[idx] = "s sp fishery(" + str(i) + "," + str(k) + ")";
						//if (s_func_type[f]<=2)
						//	x_names[idx] = "s slope(" + str(i) + "," + str(k) + ")";
						//else	
						//	x_names[idx] = "s sigma(" + str(i) + "," + str(k) + ")";
						dvarpars[idx] = value;
						dvarpars_min[idx] = s_slope_sp_fishery_min[i][k];
						dvarpars_max[idx] = s_slope_sp_fishery_max[i][k];
						parfile_names[idx-1] = "/s_sp_fishery/"+list_fishery_name[f]; 
						++idx; 
					} 
					else dvarsSslope_sp_fishery[i][k] = s_slope_sp_fishery[i][k];

					if (s_func_type[f]>1){
						if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/length_threshold", "use") == "true"){
							double value = s_length_sp_fishery[i][k];
							dvarpars_min[idx] = length[i][sp_a0_adult[i]];
							dvarpars_max[idx] = length[i][sp_nb_cohorts[i]-1]+10.0;

							x[idx] = boundpin(value, dvarpars_min[idx], dvarpars_max[idx]);
							x_names[idx] = "length threshold(" + str(i) + "," + str(k) + ")";
							//if (s_func_type[f]==2)
							//	x_names[idx] = "s length(" + str(i) + "," + str(k) + ")";
							//else	x_names[idx] = "s mean length(" + str(i) + "," + str(k) + ")";
							dvarpars[idx] = value;
							parfile_names[idx-1] = "/s_sp_fishery/"+list_fishery_name[f]+"/length_threshold";
							++idx; 
						} 
						else dvarsSlength_sp_fishery[i][k] = s_length_sp_fishery[i][k];
					}

					if (s_func_type[f]>=3){
						if (doc.get("/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote", "use") == "true"){
							double value = s_asympt_sp_fishery[i][k];
							x[idx] = boundpin(value, s_asympt_sp_fishery_min[i][k], s_asympt_sp_fishery_max[i][k]);
							x_names[idx] = "right asymptote(" + str(i) + "," + str(k) + ")";
							//x_names[idx] = "s right asymptote(" + str(i) + "," + str(k) + ")";
							dvarpars[idx] = value;
							dvarpars_min[idx] = s_asympt_sp_fishery_min[i][k];
							dvarpars_max[idx] = s_asympt_sp_fishery_max[i][k];
							parfile_names[idx-1] = "/s_sp_fishery/"+list_fishery_name[f]+"/right_asymptote";
							++idx; 
						} 
						else dvarsSasympt_sp_fishery[i][k] = s_asympt_sp_fishery[i][k];
					}
					k++;
				}
			}
		} 
		else {
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					dvarsSslope_sp_fishery[i][k] = s_slope_sp_fishery[i][k];
					if (s_func_type[f]>1)
						dvarsSlength_sp_fishery[i][k] = s_length_sp_fishery[i][k];
					if (s_func_type[f]>=3)				
						dvarsSasympt_sp_fishery[i][k] = s_asympt_sp_fishery[i][k];
					k++;
				}
			}
		}
	}
	
	dvarsLike_param.allocate(0,nb_species-1);
	for (int i=0; i<nb_species; i++) {
		const int fmax = nb_fishery_by_sp[i];
		dvarsLike_param(i).allocate(0,fmax-1);
		dvarsLike_param(i).initialize();
		if (doc.get("/likelihood_parameters/variables","use") == "true") {
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					//if (like_types[i][k]==2 || like_types[i][k]==4 || like_types[i][k]==5){
					if (/*like_types[i][k]==2 ||*/ like_types[i][k]==4 || like_types[i][k]==5){
						double value = like_param(i,k);
						dvarpars_min[idx] = 0;
						dvarpars_max[idx] = 30.0;
						if (like_types[i][k]==2) dvarpars_max[idx] = 2.0;

						x[idx] = boundpin(value, dvarpars_min[idx], dvarpars_max[idx]);
						x_names[idx] = "likelihood parameter(" + str(i) + "," + str(k) + ")"; 
						dvarpars[idx] = value;
						parfile_names[idx-1] = "/likelihood_parameters/"+list_fishery_name[f];
						++idx;
					}
					else if (like_types[i][k]==2) dvarsLike_param[i][k] = like_param[i][k];
					k++;
				}
			}
		} 
		else {
			int k = 0;
			for (int f = 0; f < nb_fishery; f++) {
				if (mask_fishery_sp[i][f]){
					dvarsLike_param[i][k] = like_param[i][k];
					k++;
				}
			}
		}
	}
//cout << __LINE__ << endl;
	dvarsProb_zero.allocate(0,nb_species-1);
	for (int i=0; i<nb_species; i++){
		const int fmax = nb_fishery_by_sp[i];
		dvarsProb_zero(i).allocate(0,fmax-1);
		dvarsProb_zero(i).initialize();
		int k = 0;
		for (int f = 0; f < nb_fishery; f++) {
			if (mask_fishery_sp[i][f]){
				if (like_types[i][k]==5) {
					if (doc.get("/likelihood_parameters/variables","use") == "true"){
						double value = prob_zero(i,k);
						dvarpars_min[idx] = 0;
						dvarpars_max[idx] = 1.0;
	
						x[idx] = boundpin(value, dvarpars_min[idx], dvarpars_max[idx]);
						x_names[idx] = "probability of zero(" + str(i) + "," + str(k) + ")"; 
						dvarpars[idx] = value;
						parfile_names[idx-1] = "/prob_zero/"+list_fishery_name[f];
						++idx;
					}
					else dvarsProb_zero[i][k] = prob_zero[i][k];
				}
				k++;
			}
		}
	}
//cout << __LINE__ << endl;
//cout << dvarpars << endl;
//cout << statpars << endl;

	outp_param(x_names,idx-1);	
}

void VarParamCoupled::par_init(dvar_vector& dvarsCoef, dvector& coef, const double coef_min, const double coef_max, string s, dvector& x, adstring_array& x_names, int& idx)  
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

	dvarsCoef.allocate(0, nb_species - 1);
	dvarsCoef.initialize();
	string sv = s + "/variable";
	if (doc.get(sv, "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			const double value = coef[i];
			x[idx] = boundpin(value, coef_min, coef_max);
			x_names[idx] = name +"(" + str(i) + ")";
			dvarpars[idx] = value;
			dvarpars_min[idx] = coef_min;
			dvarpars_max[idx] = coef_max;
			parfile_names[idx-1] = s;
			++idx;
		}
	} 
	else { 
		for (int i = 0; i < nb_species; i++) 
			dvarsCoef(i) = coef(i);
	}
}

void VarParamCoupled::par2_init(dvar_vector& dvarsCoef, dvector& coef, const double coef_min, const double coef_max, string s1, string s2, dvector& x, adstring_array& x_names, int& idx)  
{

	string n = s1+"_"+s2;
	std::string whitespace = " ";
	n.erase(0,1);
	int pos = n.find("_");
	while (pos != -1) {
		n.replace(pos,1,whitespace);
		pos = n.find("_");
	}
	adstring name = n.c_str();
	dvarsCoef.allocate(0, nb_species - 1);
	dvarsCoef.initialize();
	string sv = s1 + "/" + s2 + "/variable";
	if (doc.get(sv, "use") == "true") {
		for (int i = 0; i < nb_species; i++) {
			const double value = coef[i];
			x[idx] = boundpin(value, coef_min, coef_max);
			x_names[idx] = name +"(" + str(i) + ")";
			dvarpars[idx] = value;
			dvarpars_min[idx] = coef_min;
			dvarpars_max[idx] = coef_max;
			parfile_names[idx-1] = s1+"/"+s2;
			++idx;
		}
	} 
	else { 
		for (int i = 0; i < nb_species; i++) 
			dvarsCoef(i) = coef(i);
	}
}


bool nocase_compare(char c1, char c2)
{
	return toupper(c1) == toupper(c2);
}

void VarParamCoupled::get_param_index(ivector& ix, dmatrix& xy, dmatrix& pars)
{
       for (int idx=1; idx<=_nvarcalc; idx++){
                for (int n=0; n<nb_varproj; n++){
                        string pname = parfile_names[idx-1];
			pname.erase(0,1); // parfile name strings contain "/" at the beginning
                        int res = strcmp(pname.c_str(),varproj[n].c_str());
                        if (!res){
                                double x0 = dvarpars_min[idx];
                                double deltax = (dvarpars_max[idx]-x0)/(varproj_nsteps[n]+1);
                                for (int i=0; i<varproj_nsteps(n); i++){
                                        pars(n,i) = x0+(i+1)*deltax;
                                        xy(n,i) = boundpin(pars(n,i), dvarpars_min[idx], dvarpars_max[idx]);
                                }
                                ix[n] = idx;
                        }
                }
        }
}

void VarParamCoupled::set_all_false(string* pnames)
{
	for (int idx=1; idx<=_nvarcalc; idx++){
		string sv = "/" + pnames[idx-1] + "/variable";
		doc.set(sv,"use","false");
        }
}

int VarParamCoupled::set_var_parameters(ivector phase_par_flags, string* pnames)
{//1. set the parameter to true; 2. extend the boundary by 1/10 of allowed range; 
	int nvar = 0; 
	bool par_modified = false;
	for (int idx=1; idx<=_nvarcalc; idx++){
		if (phase_par_flags[idx]){
			string sv = "/" + pnames[idx-1] + "/variable";
			doc.set(sv,"use","true");
			//extend the boundary if the previous estimate was stuck to it
			string s = "/" + pnames[idx-1];
			double pval = doc.getDouble(s,sp_name[0]);
			double pmin = doc.getDouble(sv,"min");
			double pmax = doc.getDouble(sv,"max");
			double step = 0.1*(pmax-pmin);
			double delta = 0.001*(pmax-pmin);
			if (abs(pval-pmin)<delta) {
				par_modified = true;
				double newmin = pmin-step;
				doc.set(sv,"min",newmin);
				cout << "Attn: modified minimum value of " << pnames[idx-1] << 
					" from " << pmin << " to " << newmin << endl; 
			}
			if (abs(pval-pmax)<delta) {
				par_modified = true;
				doc.set(sv,"max",pmax+step);
				cout << "Attn: modified maximal value of " << pnames[idx-1] << 
					" from " << pmax << " to " << (pmax + step) << endl; 
			}

			nvar++;
                }
        }
	if (par_modified) re_read_varparam();
	return nvar;
}


//initialize parameter near (given the chosen step size) the lower boundary
double VarParamCoupled::par_init_lo(int ix, double eps)
{
	double x = 0.0;
	double eps2 = eps*sqrt(pow(dvarpars_max[ix]-dvarpars_min[ix],2));
	x = boundpin(dvarpars_min[ix]+eps2, dvarpars_min[ix], dvarpars_max[ix]);
	return x;	
}
//initialize parameter near (given the chosen step size) the upper boundary
double VarParamCoupled::par_init_up(int ix, double eps)
{
	double x = 0.0;
	double eps2 = eps*sqrt(pow(dvarpars_max[ix]-dvarpars_min[ix],2));
	x = boundpin(dvarpars_max[ix]-eps2, dvarpars_min[ix], dvarpars_max[ix]);
	return x;	
}

//initialize parameter within min-max boundaries using the given step: for SA-OAT
double VarParamCoupled::par_init_step(int ix, double delta)
{
	double x = 0.0;
	double x_new = dvarpars_min[ix]+delta*(dvarpars_max[ix]-dvarpars_min[ix]);
	x = boundpin(x_new, dvarpars_min[ix], dvarpars_max[ix]);
	return x;	
}



//initialize parameter by stepping left half distance to lower boundary
double VarParamCoupled::par_init_step_left(int ix)
{
	double x = 0.0;
	double step = .5*(dvarpars[ix]-dvarpars_min[ix]);
	x = boundpin(dvarpars(ix)-step, dvarpars_min[ix], dvarpars_max[ix]);
	return x;	
}
//initialize parameter by stepping right half distance to upper boundary
double VarParamCoupled::par_init_step_right(int ix) {
	double x = 0.0;
	double step = .5*(dvarpars[ix]-dvarpars_min[ix]);
	x = boundpin(dvarpars_max[ix]+step, dvarpars_min[ix], dvarpars_max[ix]);
	return x;	
}



