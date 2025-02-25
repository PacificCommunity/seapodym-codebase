#include "SeapodymCoupled.h"

void SeapodymCoupled::ReadLarvae()
{
	int nlevel = 0;
	string file_input;
	file_input = param->strfile_larvae;
	cout << "Reading input larvae file: "<< file_input << endl;
	rw.rbin_headpar(file_input, nlon_input, nlat_input, nlevel);

	if (param->larvae_input_aggregated_flag[0]){
		int nb_larvae_input_agg_groups = param->nb_larvae_input_agg_groups;
		if (nlevel != nb_larvae_input_agg_groups){
			cerr << "Error[" << __FILE__ << ':' << __LINE__ << "]: The number of nlevels in \"" << file_input << " does not match the number of groups from <larvae_input_aggregation_imonths> in the parameter file.\"\n";
			exit(1);
		}

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
		for (int iAgg=0; iAgg<param->nb_larvae_input_agg_groups; iAgg++){
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
	ntime_agg.allocate(0,nbg-1);
	ntime_agg.initialize();

	qmld.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	qmld = 1.0;//.initialize();

	if (param->larvae_input_aggregated_flag[0]){
		// Aggregated larvae density over the entire period, only at obs. locations	
		kinf.allocate(0, nbg-1);
		ksup.allocate(0, nbg-1);
		for (int k = 0; k < nbg; k++){
			kinf[k] = 0;
			ksup[k] = mat.aggregated_larvae_input_vectors[k].size();
		}
		Agg_larvae_density_pred_at_obs.allocate(0, nbg-1, kinf, ksup);
		Agg_larvae_density_pred_at_obs.initialize();	
	}else{
		Larvae_density_pred.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
		Larvae_density_pred.initialize();
	}
}

void SeapodymCoupled::extract_larvae(const int sp, const int tcur)
{//Autodif function for the moment. Need to write adjoint!!!

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

	if (param->larvae_input_aggregated_flag[0]){
		int iAgg = Utilities::iTimeOfYear(month, param->larvae_input_aggregation);
		
		// Aggregate larvae density at larvae obs locations
		for (auto k=0u; k<mat.aggregated_larvae_input_vectors[iAgg].size(); k++){
			int iv = mat.aggregated_larvae_input_vectors_i[iAgg][k];
			int jv = mat.aggregated_larvae_input_vectors_j[iAgg][k];
			
			Agg_larvae_density_pred_at_obs(iAgg, k) +=  qmld(iv,jv)*mat.dvarDensity[sp][0][iv][jv];
		}
		ntime_agg[iAgg] += 1;
	}else{
		for (int i = map.imin; i <= map.imax; i++){
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++){
			    Larvae_density_pred[i][j] = qmld[i][j] * mat.dvarDensity[sp][0][i][j];
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


	//Early larvae's mortality rate
	func.M_early_sp(*param, map, M,mat.sst[tcur],mat.np1[tcur]*param->pp_transform,sp);

	//Add mortality to central diagonal 
	pop.Precalrec_juv(map, mat, M, tcur,elarvae_dt);//checked

	//2.0.1 Early larvae movement and mortality ADRE solver
	pop.Calrec_juv(map, mat, mat.dvarDensity[sp][0], M, tcur,elarvae_dt);//checked

	
	//2.0.2 Aggregate larvae density at larvae obs locations for the likelihood
	if (time_getpred)
		extract_larvae(sp,tcur);		
	
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
			Agg_larvae_density_pred_at_obs(iAgg, k) /= ntime_agg[iAgg];
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

