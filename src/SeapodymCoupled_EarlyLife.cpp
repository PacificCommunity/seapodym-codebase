#include "SeapodymCoupled.h"

void SeapodymCoupled::ReadEarly(string what)
{
	int nlevel = 0;
	string file_input;
	int input_aggregated_flag, nb_input_agg_groups;
	D3_ARRAY* input = nullptr;
	std::vector<double> (*aggregated_input_vectors)[12] = nullptr;
	std::vector<int> (*aggregated_input_vectors_i)[12] = nullptr;
	std::vector<int> (*aggregated_input_vectors_j)[12] = nullptr;
	if (what == "larvae"){
		file_input = param->strfile_larvae;
		input_aggregated_flag = param->larvae_input_aggregated_flag[0];
		nb_input_agg_groups = param->nb_larvae_input_agg_groups;
		input = &mat.larvae_input;
		aggregated_input_vectors = &mat.aggregated_larvae_input_vectors;
		aggregated_input_vectors_i = &mat.aggregated_larvae_input_vectors_i;
		aggregated_input_vectors_j = &mat.aggregated_larvae_input_vectors_j;
	}else if (what == "spawning"){
		file_input = param->strfile_spawning;
		input_aggregated_flag = param->spawning_input_aggregated_flag[0];
		nb_input_agg_groups = param->nb_spawning_input_agg_groups;
		input = &mat.spawning_input;
		aggregated_input_vectors = &mat.aggregated_spawning_input_vectors;
		aggregated_input_vectors_i = &mat.aggregated_spawning_input_vectors_i;
		aggregated_input_vectors_j = &mat.aggregated_spawning_input_vectors_j;
	}else{
		cerr << "Error: in ReadEarly(), 'what' argument is not recognized" << endl;
		std::exit(EXIT_FAILURE);			
	}
	cout << "Reading input " << what << " file: "<< file_input << endl;
	rw.rbin_headpar(file_input, nlon_input, nlat_input, nlevel);

	if (input_aggregated_flag){
		if (nlevel != nb_input_agg_groups){
			cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: The number of nlevels in \"" << file_input << " does not match the number of groups from <" << what << "_input_aggregation_imonths> in the parameter file.\"\n";
			exit(1);
		}

		(*input).allocate(0,nb_input_agg_groups-1);
		for (int iAgg=0; iAgg<nb_input_agg_groups; iAgg++){
			(*input)[iAgg].allocate(1, nlon_input, 1, nlat_input);
			(*input)[iAgg].initialize();
		}

		for (int iAgg=0; iAgg<nb_input_agg_groups; iAgg++){
			ifstream litbin(file_input.c_str(), ios::binary | ios::in);
			const int sizeofDymInputType = sizeof(float);
			float buf;
			if (!litbin){
				cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: Unable to read file \"" << file_input << "\"\n";
				exit(1);
			}
			int nbytetoskip_aggregatedfile = (9 +(3* nlat * nlon) + nb_input_agg_groups + (nlat *nlon*iAgg)) * 4;
			litbin.seekg(nbytetoskip_aggregatedfile, ios::cur);
			for (int j=0;j<nlat_input;j++){
				for (int i=0;i<nlon_input;i++){
					litbin.read(( char *)&buf,sizeofDymInputType);
					(*input)[iAgg][i+1][j+1]= buf;
				}
			}
			litbin.close();
		}

		// Vector of non-NA observed density
		int ndata=0;
		for (int iAgg=0; iAgg<nb_input_agg_groups; iAgg++){
			for (int j=0;j<nlat_input;j++){
				for (int i=0;i<nlon_input;i++){
					if (map.carte[i+1][j+1]){
						if ((*input)[iAgg][i+1][j+1]>=0){
							(*aggregated_input_vectors)[iAgg].push_back((*input)[iAgg][i+1][j+1]);
							(*aggregated_input_vectors_i)[iAgg].push_back(i+1);
							(*aggregated_input_vectors_j)[iAgg].push_back(j+1);
							ndata += 1;
						}
					}
				}
			}
		}
		cout << "Number of data points: " << ndata << endl;
	}else{
		int needed_nlevels = nbt_total-nbt_building-param->nbsteptoskip;
		if (nlevel != needed_nlevels){
			cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: The number of nlevels in \"" << file_input << " (" << nlevel << ") does not match the necessary number of time steps (" << needed_nlevels << ").\"\n";
			exit(1);

		}
	}
}

void SeapodymCoupled::create_init_larvae_vars()
{
	int nbg = param->nb_larvae_input_agg_groups;
	ntime_agg_larvae.allocate(0,nbg-1);
	ntime_agg_larvae.initialize();

	qmld.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	qmld = 1.0;//.initialize();

	if (param->larvae_input_aggregated_flag[0]){
		// Aggregated larvae density over the entire period, only at obs. locations	
		kinf_larvae.allocate(0, nbg-1);
		ksup_larvae.allocate(0, nbg-1);
		for (int k = 0; k < nbg; k++){
			kinf_larvae[k] = 0;
			ksup_larvae[k] = mat.aggregated_larvae_input_vectors[k].size();
		}
		Agg_larvae_density_pred_at_obs.allocate(0, nbg-1, kinf_larvae, ksup_larvae);
		Agg_larvae_density_pred_at_obs.initialize();	
	}else{
		Larvae_density_pred.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
		Larvae_density_pred.initialize();
	}
}

void SeapodymCoupled::create_init_spawning_vars()
{
	int nbg = param->nb_spawning_input_agg_groups;
	ntime_agg_spawning.allocate(0,nbg-1);
	ntime_agg_spawning.initialize();

	if (param->spawning_input_aggregated_flag[0]){
		// Aggregated spawning index over the entire period, only at obs. locations	
		kinf_spawning.allocate(0, nbg-1);
		ksup_spawning.allocate(0, nbg-1);
		for (int k = 0; k < nbg; k++){
			kinf_spawning[k] = 0;
			ksup_spawning[k] = mat.aggregated_spawning_input_vectors[k].size();
		}
		Agg_SBHs_pred_at_obs.allocate(0, nbg-1, kinf_spawning, ksup_spawning);
		Agg_SBHs_pred_at_obs.initialize();	
	}else{
		SBHs_pred.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
		SBHs_pred.initialize();
	}
}

void SeapodymCoupled::extract_early(const int sp, const int tcur, string what)
{//Autodif function for the moment. Need to write adjoint!!!
	int input_aggregated_flag, nb_input_agg_groups;
	D3_ARRAY* input = nullptr;
	std::vector<double> (*aggregated_input_vectors)[12] = nullptr;
	std::vector<int> (*aggregated_input_vectors_i)[12] = nullptr;
	std::vector<int> (*aggregated_input_vectors_j)[12] = nullptr;
	ivector *ntime_agg = nullptr;

	if (what == "larvae"){
		if (param->q_mld_larvae && sp == 0){ //only once as it is species independent
			double slope = param->q_mld_slope;
			double depth = param->q_mld_depth/1000.0; //vld units is km 
			for (int i = map.imin; i <= map.imax; i++){	
				const int jmin = map.jinf[i];
				const int jmax = map.jsup[i];
				for (int j = jmin; j <= jmax; j++){

					qmld(i,j) = 1.0/(1.0+exp(slope*(mat.vld[tcur][i][j]-depth)));
				}
			}
		}

		input_aggregated_flag = param->larvae_input_aggregated_flag[0];
		nb_input_agg_groups = param->nb_larvae_input_agg_groups;
		input = &mat.larvae_input;
		aggregated_input_vectors = &mat.aggregated_larvae_input_vectors;
		aggregated_input_vectors_i = &mat.aggregated_larvae_input_vectors_i;
		aggregated_input_vectors_j = &mat.aggregated_larvae_input_vectors_j;
		ntime_agg = &ntime_agg_larvae;
	}else if (what == "spawning"){
		input_aggregated_flag = param->spawning_input_aggregated_flag[0];
		nb_input_agg_groups = param->nb_spawning_input_agg_groups;
		input = &mat.spawning_input;
		aggregated_input_vectors = &mat.aggregated_spawning_input_vectors;
		aggregated_input_vectors_i = &mat.aggregated_spawning_input_vectors_i;
		aggregated_input_vectors_j = &mat.aggregated_spawning_input_vectors_j;
		ntime_agg = &ntime_agg_spawning;
	}else{
		cerr << "Error: in extract_early(), 'what' argument is not recognized" << endl;
		std::exit(EXIT_FAILURE);			
	}

	if (input_aggregated_flag){
		int iAgg;
		if (what == "larvae"){
			iAgg = Utilities::iTimeOfYear(month, param->larvae_input_aggregation);
		}else{
			iAgg = Utilities::iTimeOfYear(month, param->spawning_input_aggregation);
		}
		
		// Aggregate larvae density at larvae obs locations
		for (auto k=0u; k<(*aggregated_input_vectors)[iAgg].size(); k++){
			int iv = (*aggregated_input_vectors_i)[iAgg][k];
			int jv = (*aggregated_input_vectors_j)[iAgg][k];
			
			if (what == "larvae"){
				Agg_larvae_density_pred_at_obs(iAgg, k) +=  qmld(iv,jv)*mat.dvarDensity[sp][0][iv][jv];
			}else{
				Agg_SBHs_pred_at_obs(iAgg, k) +=  Spawning_Habitat(iv,jv)*Total_pop(iv,jv);
			}
		}
		(*ntime_agg)[iAgg] += 1;
	}else{
		for (int i = map.imin; i <= map.imax; i++){
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++){
				if (what == "larvae"){
			    	Larvae_density_pred[i][j] = qmld[i][j] * mat.dvarDensity[sp][0][i][j];
				}else{
			    	SBHs_pred[i][j] = Spawning_Habitat[i][j] * Total_pop[i][j];
				}
			}
		}
	}
	
/*
//CODE COPIED FROM RUN-COUPLED, which is to be executed after spawning. 
//Without mortality_sst should be the same as above. 
					const int nb_lv = param->sp_nb_cohort_lv[sp];
					int iAgg = Utilities::iTimeOfYear(month, param->larvae_input_aggregation);

					// Calculate scaling factor between larvae density at 1st time step and larvae density at age_larvae_before_sst_mortality days, based on sst-dependent mortality (if applicable)
					if (param->larvae_mortality_sst[sp]){
						func.Scaling_factor_sstdep_larvae_mortality(*param, mat.sst[tcur], map, mat.dvarScaling_factor_sstdep_larvae_mortality[sp], sp);
					}

					// Aggregate larvae density at larvae obs locations
					for (auto k=0u; k<mat.aggregated_larvae_input_vectors[iAgg].size(); k++){
						if (param->larvae_mortality_sst[sp]){// In this case, it only considers the first age class (even if nb_lv>1)
							Agg_larvae_density_pred_at_obs(iAgg, k) +=  mat.dvarDensity[sp][0][mat.aggregated_larvae_input_vectors_i[iAgg][k]][mat.aggregated_larvae_input_vectors_j[iAgg][k]] * mat.dvarScaling_factor_sstdep_larvae_mortality[sp][mat.aggregated_larvae_input_vectors_i[iAgg][k]][mat.aggregated_larvae_input_vectors_j[iAgg][k]];
						}else{
							for (int age=0; age<nb_lv; age++){
								Agg_larvae_density_pred_at_obs(iAgg, k) += mat.dvarDensity[sp][0][mat.aggregated_larvae_input_vectors_i[iAgg][k]][mat.aggregated_larvae_input_vectors_j[iAgg][k]];
							}
						}
					}
					ntime_agg[iAgg] += 1;
 
 */
}


void SeapodymCoupled::elarvae_model_run(dvar_matrix& M, const int sp, const int tcur, bool time_getpred, bool writeoutputfiles)
{
	double sigma_fcte_save = param->sigma_fcte;

	//2.0.0 Prepare diagonal coefficients for early-life movement and mortality
	//The movement rates for early larvae:	
	param->sigma_fcte *= elarvae_dt;
	mat.u = elarvae_dt*mat.un[tcur][0]; 
	mat.v = elarvae_dt*mat.vn[tcur][0]; 

	//Precompute diagonal coefficients for early larvae ADRE		
	pop.precaldia(*param, map, mat);
	pop.caldia(map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y, mat.advection_y);

	func.Early_Mortality_Sp(*param, mat, map, M, sp, tcur);

	//Add mortality to central diagonal 
	pop.Precalrec_juv(map, mat, M, tcur,elarvae_dt);//checked

	//2.0.1 Early larvae movement and mortality ADRE solver
	pop.Calrec_juv(map, mat, mat.dvarDensity[sp][0], M, tcur,elarvae_dt);//checked

	
	//2.0.2 Aggregate larvae density at larvae obs locations for the likelihood
	if (time_getpred)
		extract_early(sp,tcur,"larvae");		
	
	if (writeoutputfiles)
		write_elarvae_dym(sp);
	//2.1.0 Now, the larval dynamics over remaining time of the first age class
	//Get the diagonal coefficients with movement rates per remaining time step
	param->sigma_fcte= (1-elarvae_dt)*sigma_fcte_save;
	mat.u = (1-elarvae_dt)*mat.un[tcur][0]; mat.v = (1-elarvae_dt)*mat.vn[tcur][0]; 
	//Precompute diagonal coefficients for larvae and juvenile ADREs
	pop.precaldia(*param, map, mat);
	pop.caldia(map, *param, mat.diffusion_x, mat.advection_x, mat.diffusion_y, mat.advection_y);
	
}


void SeapodymCoupled::get_larvae_at_obs()
{
	// Compute the average larvae density over the entire period
	
	for (int iAgg = 0; iAgg < param->nb_larvae_input_agg_groups; iAgg++){
		for (auto k=0u; k<mat.aggregated_larvae_input_vectors[iAgg].size(); k++){
			Agg_larvae_density_pred_at_obs(iAgg, k) /= ntime_agg_larvae[iAgg];
		}
	}
	
}

void SeapodymCoupled::get_SBHs_at_obs()
{
	// Compute the average SB x Hs over the entire period
	
	for (int iAgg = 0; iAgg < param->nb_spawning_input_agg_groups; iAgg++){
		for (auto k=0u; k<mat.aggregated_spawning_input_vectors[iAgg].size(); k++){
			Agg_SBHs_pred_at_obs(iAgg, k) /= ntime_agg_spawning[iAgg];
		}
	}
	
}


void SeapodymCoupled::write_elarvae_dym(const int sp)
{
	mat.larvae(sp).initialize();			
	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			
			mat.larvae(sp,i,j) = value(mat.dvarDensity(sp,0,i,j));
		}
	}
	WriteAVariableDym(mat.larvae(sp),param->sp_name[sp] + "_early_larvae.dym",false);
}

// Functions to compute the likelihood of a larvae density or SBHs observed on a continuous scale

dvariable gaussian_comp(double N_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp, string what){
	dvariable h;
	dvariable sigma;
	if (what == "larvae"){
		h = param.dvarsQ_sp_larvae[sp];
		sigma = param.dvarsLikelihood_larvae_sigma[sp];
	}else{
		h = param.dvarsQ_sp_spawning[sp];
		sigma = param.dvarsLikelihood_spawning_sigma[sp];
	}

	dvariable L_pred = N_pred * h;
	dvariable lkhd = 0.0;

	if (N_obs==0.0)
		lkhd = weight_Lobszero*L_pred*L_pred/(2.0*pow(sigma, 2.0)) ;
	else
		lkhd = pow(N_obs-L_pred, 2.0)/(2.0 * pow(sigma, 2.0));
	
	return 1000.0*lkhd;
}

dvariable poisson_comp(double L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp, string what){
	dvariable h;
	dvariable sigma;
	if (what == "larvae"){
		h = param.dvarsQ_sp_larvae[sp];
		sigma = param.dvarsLikelihood_larvae_sigma[sp];
	}else{
		h = param.dvarsQ_sp_spawning[sp];
		sigma = param.dvarsLikelihood_spawning_sigma[sp];
	}

	const double twopi = 2.0*3.141592654;
    dvariable L_pred = N_pred * h;
    dvariable lkhd = 0.0;
    if (L_obs==0){
        lkhd = weight_Lobszero * (pow(L_pred,2) / (2*pow(sigma, 2)) + log(sigma) + log(twopi)/2);
    }else{
        lkhd = L_pred - L_obs * log(L_pred) + gammln(L_obs+1);
    }
    return lkhd;
}

dvariable truncated_poisson_comp(double L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp, string what){
	dvariable h;
	if (what == "larvae"){
		h = param.dvarsQ_sp_larvae[sp];
	}else{
		h = param.dvarsQ_sp_spawning[sp];
	}

    dvariable L_pred = 1 + N_pred * h;
    dvariable lkhd = 0.0;
    L_obs += 1;
    lkhd = L_pred - L_obs * log(L_pred) + gammln(L_obs+1) + log(1-exp(-L_pred));
    if (L_obs==1.0){
        lkhd *= weight_Lobszero;
    }
    return lkhd;
}

dvariable zinb_comp(double L_obs, dvariable N_pred, VarParamCoupled& param, int sp, string what){
	dvariable h;
	dvariable beta;
	dvariable p;
	if (what == "larvae"){
		h = param.dvarsQ_sp_larvae[sp];
		beta = param.dvarsLikelihood_larvae_beta[sp];
		p = param.dvarsLikelihood_larvae_probzero[sp];
	}else{
		h = param.dvarsQ_sp_spawning[sp];
		beta = param.dvarsLikelihood_spawning_beta[sp];
		p = param.dvarsLikelihood_spawning_probzero[sp];
	}

    dvariable L_pred = N_pred * h;
    dvariable lkhd = 0.0;
    if (L_obs==0.0){
        dvariable pwr = beta*L_pred/(1-p);
        lkhd -= log(p+(1-p)*pow(beta/(1.0+beta),pwr));
    }else{
        dvariable mu = L_pred/(1-p);
        lkhd -= log(1-p) + gammln(beta*mu+L_obs) - gammln(beta*mu) -gammln(L_obs+1.0) + beta*mu*log(beta)-log(beta+1.0)*(beta*mu+L_obs);
    }
    return lkhd;
}

dvariable zip_comp(double L_obs, dvariable N_pred, VarParamCoupled& param, int sp, string what){
	dvariable h;
	dvariable p;
	if (what == "larvae"){
		h = param.dvarsQ_sp_larvae[sp];
		p = param.dvarsLikelihood_larvae_probzero[sp];
	}else{
		h = param.dvarsQ_sp_spawning[sp];
		p = param.dvarsLikelihood_spawning_probzero[sp];
	}

    dvariable L_pred = N_pred * h;
    dvariable lkhd = 0.0;
    if (L_obs==0.0){
        lkhd -= log(p + (1-p) * exp(-L_pred));
    }else{
        lkhd -= log(1-p) + L_obs * log(L_pred) - L_pred - gammln(L_obs+1.0);
    }
    return lkhd;
}

dvariable lognormal_comp(double L_obs, dvariable N_pred, VarParamCoupled& param, int sp, string what){
    const double twopi = 2.0*3.141592654;
    dvariable h     = param.dvarsQ_sp_larvae[sp];
    dvariable sigma = param.dvarsLikelihood_larvae_sigma[sp];
    dvariable L_pred = N_pred * h;
    if (value(L_pred) <= 0.0)
        return dvariable(50.0);
    dvariable lkhd = pow(log(L_obs) - log(L_pred), 2.0) / (2.0 * pow(sigma, 2.0))
                     + log(sigma) + 0.5*log(twopi) + log(L_obs);
    return lkhd;
}
