#include "SeapodymCoupled.h"
#include "NishikawaLike.h"

void mat2vect(const PMap& map, const dmatrix effort, const dmatrix catch_obs, dvar_matrix catch_est, dvector& data_obs, dvar_vector& data_est, bool cpue, double c_units, const int catch_removal_fishery, const int cdata_temp_reso, const int month);
void LeastSquares(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood);
void Concentrated(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood, const int nobs);
void Normal(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood, dvariable& sigma, const int nobs);
void LogNormal(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood, dvariable& sigma, const int nobs);
void ZILogNormal(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood, dvariable& sigma, dvariable& prob, const int nobs);
dvariable Poisson(const dvector data_obs, dvar_vector& data_est, const double min_obs_catch, const int nobs);
dvariable TruncatedPoisson(const dvector data_obs, dvar_vector& data_est, const int nobs);
dvariable Exponential(const dvector data_obs, dvar_vector& data_est, const int nobs);
dvariable Weibull(const dvector data_obs, dvar_vector& data_est, const int nobs);
dvariable NegBinomial(const dvector data_obs, dvar_vector& data_est, dvariable& beta, const int nobs);
dvariable ZINegBinomial(const dvector data_obs, dvar_vector& data_est, dvariable& beta, dvariable& p, const int nobs);
dvariable LFlike_robust(const d3_array LF_qtr_obs, dvar3_array& dvarLF_est, const int a0, const int nb_ages, const int nb_regions_sp, const int k);

int number_data_like;
const double twopi = 2.0*3.141592654;

dvariable SeapodymCoupled::like(const int sp, const int k, const int f, const int nobs)
{
	dvariable likelihood = 0.0;
	dvariable beta, prob, sigma;

	dvector cdata_obs;
	dvar_vector cdata_est;
	cdata_obs.allocate(0,nobs-1);
	cdata_est.allocate(0,nobs-1);

	int cdata_temp_reso = param->catch_treso_likelihood[sp][f];

	bool cpue = param->cpue;
	double c_units = param->catch_units_converter(f);
	if (cpue)
		c_units = param->cpue_units_converter(f);
	
	int like_type = param->like_types[sp][k];
	const double fw_catch  = param->catch_like_weight(f);  
	const double fw_length = param->length_like_weight(f);  

	mat2vect(map, mat.effort(f), mat.catch_obs(sp,k), mat.dvarCatch_est(sp,k), cdata_obs, cdata_est, 
		cpue, c_units, param->mask_fishery_sp_no_effort(sp,f),cdata_temp_reso, month);
 
	if (cdata_temp_reso==1 || (cdata_temp_reso==3 && (month==3 || month==6 || month==9 || month==12))){
		switch (like_type) {

			case 1: LeastSquares(cdata_obs,cdata_est,likelihood);
				break;
			
			case 2: Concentrated(cdata_obs,cdata_est,likelihood,nobs);
				break;
		
			case 3: likelihood = Poisson(cdata_obs,cdata_est,param->poisson_like_min_catch,nobs);
				break;

			case 4: likelihood = TruncatedPoisson(cdata_obs,cdata_est,nobs);
				break;	
				
			case 5: sigma = param->dvarsLike_param(sp,k);
				Normal(cdata_obs,cdata_est,likelihood,sigma,nobs);
				break;

			case 6: sigma = param->dvarsLike_param(sp,k);
				LogNormal(cdata_obs,cdata_est,likelihood,sigma,nobs);
				break;

			case 7: sigma = param->dvarsLike_param(sp,k);
				prob = param->dvarsProb_zero(sp,k);
				ZILogNormal(cdata_obs,cdata_est,likelihood,sigma,prob,nobs);
				break;					
					
			case 8: beta = param->dvarsLike_param(sp,k);
				likelihood = NegBinomial(cdata_obs,cdata_est,beta,nobs);
				break;
			
			case 9: beta = param->dvarsLike_param(sp,k);
				prob = param->dvarsProb_zero(sp,k);
				likelihood = ZINegBinomial(cdata_obs,cdata_est,beta,prob,nobs);
				break;

			case 10: likelihood = Exponential(cdata_obs,cdata_est,nobs);
				break;

			case 11: likelihood = Weibull(cdata_obs,cdata_est,nobs);
				break;

        	}
		likelihood = fw_catch * likelihood;
		clike_fishery[f] += value(likelihood);
	}


	//length frequency likelihood
	if (param->frq_like[sp] && (month==3 || month==6 || month==9 || month==12)){
		const int a0 = a0_adult(sp);
		const int nb_ages = aN_adult(sp);
		const int nb_regions_sp = param->nb_region_sp_B(sp);

		//LF likelihood value for fishery 'f'
		dvariable lf_like = fw_length * LFlike_robust(mat.LF_qtr_obs(sp),mat.dvarLF_est(sp),
						  a0,nb_ages,nb_regions_sp,k);
		lflike_fishery[f] += value(lf_like);
		//lflike += value(lf_like);

		likelihood += lf_like;
	} 
	return likelihood;
}

double SeapodymCoupled::get_stock_like(dvariable& total_stock, dvariable& likelihood)
{//returns double value of stock likelihood.

	double stocklike  = 0.0;
	for (int sp=0; sp < nb_species; sp++){
		if (!param->stock_like[sp]) return 0;
		else {
			if (!param->scalc()) cout << "stock size: " << total_stock << endl;
			double mean_total_stock_obs = param->mean_stock_obs[sp];  
			if (total_stock > mean_total_stock_obs){
				dvariable slike = 0.5*pow(total_stock-mean_total_stock_obs,2);
				likelihood += slike;
				stocklike  += value(slike);
			}
		}
	}
	return stocklike;
}

double SeapodymCoupled::get_larvae_like(dvariable& likelihood, dvar_matrix& Agg_larvae_density_pred_at_obs)
{//returns double value of larvae likelihood.

	NishikawaCategories NshkwCat(*param);
	double likelihood_penalty = 50.0;
	double weight_Lobszero = 0.0;
	if (param->fit_null_larvae[0]==1){
		weight_Lobszero = param->weight_null_larvae[0];
	}

	double larvaelike  = 0.0;
	for (int sp=0; sp < nb_species; sp++){
		int like_type = param->larvae_likelihood_type[sp];
		for (int iAgg=0; iAgg<param->nb_larvae_input_agg_groups; iAgg++){
			for (auto k=0u; k<mat.aggregated_larvae_input_vectors[iAgg].size(); k++){
				// Compute likelihood
				dvariable N_pred = 0.0;
				N_pred = Agg_larvae_density_pred_at_obs(iAgg, k);
				dvariable lkhd;
				if (param->larvae_input_categorical_flag[sp]){
					int L_obs  = mat.aggregated_larvae_input_vectors[iAgg][k];
					lkhd = larvae_like(like_type, L_obs, N_pred, weight_Lobszero, likelihood_penalty, NshkwCat);

				}else{
					double L_obs  = mat.aggregated_larvae_input_vectors[iAgg][k];
					lkhd = larvae_like(like_type, L_obs, N_pred, weight_Lobszero, likelihood_penalty, NshkwCat);
				}
				likelihood += lkhd;
				larvaelike += value(lkhd);
			}
		}
	}
	return larvaelike;
}

double SeapodymCoupled::get_larvae_like(dvariable& likelihood, dvar_matrix& Larvae_density_pred, D3_ARRAY larvae_input, int t)
{//returns double value of larvae likelihood.

	NishikawaCategories NshkwCat(*param);
	double likelihood_penalty = 50.0;
	double weight_Lobszero = 0.0;
	if (param->fit_null_larvae[0]==1){
		weight_Lobszero = param->weight_null_larvae[0];
	}

	double larvaelike  = 0.0;
	for (int sp=0; sp < nb_species; sp++){
		int like_type = param->larvae_likelihood_type[sp];
		const int imin = map.imin;
		const int imax = map.imax;
		for (int i = imin; i <= imax; i++){
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin ; j <= jmax; j++){
				if (map.carte[i][j]){
					// Compute likelihood
					dvariable N_pred = Larvae_density_pred(i,j);
					dvariable lkhd;
					if (param->larvae_input_categorical_flag[sp]){
						int L_obs  = larvae_input(t,i,j);
						lkhd = larvae_like(like_type, L_obs, N_pred, weight_Lobszero, likelihood_penalty, NshkwCat);
					}else{
						double L_obs  = larvae_input(t,i,j);
						lkhd = larvae_like(like_type, L_obs, N_pred, weight_Lobszero, likelihood_penalty, NshkwCat);
					}
					likelihood += lkhd;
					larvaelike += value(lkhd);
				}
			}
		}
	}
	return larvaelike;
}


dvariable SeapodymCoupled::larvae_like(int like_type, int L_obs, dvariable N_pred, double weight_Lobszero, double likelihood_penalty, NishikawaCategories NshkwCat){
	dvariable lkhd = 0.0;
	switch (like_type){
		case 0: // Mixed Gaussian Kernel cost function
			lkhd = NshkwCat.mixed_gaussian_comp(L_obs, N_pred, weight_Lobszero, *param, 0);
			break;

		case 1:	// Categorical Poisson cost function
			if (N_pred == 0.0){
				if (L_obs > 0){
					lkhd = likelihood_penalty;
				}
			}else{
				lkhd = NshkwCat.categorical_poisson_comp(L_obs, N_pred, weight_Lobszero, *param, 0);
			}
			break;

		case 2: // Categorical Truncated Poisson cost function
			if (N_pred == 0.0){
				if (L_obs > 0){
					lkhd = likelihood_penalty;
				}
			}else{
				lkhd = NshkwCat.categorical_truncated_poisson_comp(L_obs, N_pred, weight_Lobszero, *param, 0);
			}
			break;

		case 3: // Zero-Inflated Negative Binomial cost function
			lkhd = NshkwCat.categorical_zinb_comp(L_obs, N_pred, *param, 0);
			break;

		case 4: // Zero-Inflated Poisson cost function
			lkhd = NshkwCat.categorical_zip_comp(L_obs, N_pred, *param, 0);
			break;
	}

	if (std::isinf(value(lkhd))){
		if (lkhd > 0){
			lkhd = likelihood_penalty;
		}else{
			lkhd = 0;
		}
	}

	return(lkhd);
}

dvariable SeapodymCoupled::larvae_like(int like_type, double L_obs, dvariable N_pred, double weight_Lobszero, double likelihood_penalty, NishikawaCategories NshkwCat){
	dvariable lkhd = 0.0;
	switch (like_type){
		case 0: // Gaussian cost function
			lkhd = gaussian_comp(L_obs, N_pred, weight_Lobszero, *param, 0);
			break;

		case 1:{// Poisson cost function
			if (N_pred == 0.0){
				if (L_obs > 0){
					lkhd = likelihood_penalty;
				}
			}else{
				lkhd = poisson_comp(L_obs, N_pred, weight_Lobszero, *param, 0);
			}
			break;
		}

		case 2: // Truncated Poisson cost function
			if (N_pred == 0.0){
				if (L_obs > 0){
					lkhd = likelihood_penalty;
				}
			}else{
				lkhd = truncated_poisson_comp(L_obs, N_pred, weight_Lobszero, *param, 0);
			}
			break;

		case 3: // Zero-Inflated Negative Binomial cost function
			lkhd = zinb_comp(L_obs, N_pred, *param, 0);
			break;

		case 4: // Zero-Inflated Poisson cost function
			lkhd = zip_comp(L_obs, N_pred, *param, 0);
			break;
	}
	return(lkhd);
}

void SeapodymCoupled::get_catch_lf_like(dvariable& likelihood)
{//returns double value of catch likelihood. Value of LF likelihood is stored in lflike.

	for (int sp=0; sp < nb_species; sp++){
		int k = 0;
		for (int f=0; f<nb_fishery; f++){
			if (param->mask_fishery_sp[sp][f]){
				if (param->mask_fishery_sp_like[sp][f]){
					int y = year-(int)param->save_first_yr;
					int nobs = rw.get_numrec(f,y,month);
					likelihood += like(sp,k,f,nobs);
				}
				k++;
			}
		}
	}
}

double SeapodymCoupled::get_tag_like(dvariable& likelihood, bool writeoutputs)
{//Returns double value of tag recapture likelihood.
 //This routine requires major revision

	double taglike = 0.0;
	if (month==1 || month==4 || month==7 || month==10){
		rec_obs_like.initialize();
		rec_pred_like.initialize();
	}
	for (int p=0; p<nb_tagpops; p++){
//if (sum(mat.dvarDensity(p+1))>0)		
//TTRACE(p+1,sum(mat.dvarDensity(p+1)))		
		//int nb_obs = 0;
		if (t_count==t_count_rec(p)){

			const int imin = map.imin; 
			const int imax = map.imax; 
			for (int i = imin; i <= imax; i++){
				const int jmin = map.jinf[i];
				const int jmax = map.jsup[i];
				for (int j = jmin ; j <= jmax; j++){
					if (map.carte[i][j]){
										
						double xx = param->itolon(i);
						double yy = param->jtolat(j);
						for (int ii=0; ii<nx_obs; ii++)
						for (int jj=0; jj<ny_obs; jj++){
							if (xx>xlon[ii] && xx<=xlon[ii+1]&& yy<=ylat[jj] && yy>ylat[jj+1]){
								for (int aa=a0_adult[0]; aa<aN_adult[0]; aa++)
									rec_pred(p,ii,jj) += mat.dvarDensity(p+1,aa,i,j)*cell_area/mat.lat_correction(j);
									//rec_pred(p,ii,jj) += 0.001*mat.dvarDensity(p+1,aa,i,j);
								}	
							}
						}
					}
				}
				
				mat.dvarDensity(p+1).initialize();
				tagpop_age_solve(p,t_count).initialize();

				//append to the aggregated predictions and observations
				rec_obs_like  += elem_prod(rec_obs(p),tlib_obs(p));
				rec_pred_like += elem_prod(rec_pred(p),tlib_obs(p));
				//rec_obs_like  += rec_obs(p);
				//rec_pred_like += rec_pred(p);
/*		
				//1. Concentrated
				taglike += value(norm2(rec_obs(p)-rec_pred(p)));
				likelihood += norm2(rec_obs(p)-rec_pred(p));
*/
/*				//2. Poisson
				double sf = 100.0;
				for (int ii=0; ii<nx_obs; ii++){
					for (int jj=0; jj<ny_obs; jj++){
						dvariable pred = sf*rec_pred(p,ii,jj);
						const double obs = sf*rec_obs(p,ii,jj);
						if (pred>0){
							likelihood += pred - obs*log(pred) + gammln(obs+1.0);
							taglike += value(pred - obs*log(pred) + gammln(obs+1.0));
						}
					}
				}
*/	
				//3. Weighted Logarithmic
				//spatial 2d
				//int nb_obs = sum(rec_obs(p));
				//const int ww = 2.0;
/*				double sf = ww*(1.0-nb_obs/(5.0+nb_obs));
				taglike += sf*value(norm2(log(rec_obs(p)+1e-1)-log(rec_pred(p)+1e-1)));
				likelihood += sf*norm2(log(rec_obs(p)+1e-1)-log(rec_pred(p)+1e-1));
*/				//1d (by lontigude and by latitude)
		
				if (writeoutputs){
					//Writing only in simulation mode
					if (!param->gcalc()){
		
						ofstream wtxt;
					std::ostringstream ostr;
					ostr << year;
					if (month>9) ostr << month <<15;
					else
						ostr << 0 << month <<15;
					string file_out = param->strout_tags + param->sp_name[0] + "_tags_pred_"  + ostr.str() + ".txt";
					wtxt.open(file_out.c_str(), ios::out);
					
					if (wtxt){
						wtxt << xlon << endl;
						wtxt << ylat << endl;
						wtxt << trans(value(rec_pred(p))) << endl;
					}
					wtxt.close();
			
					file_out = param->strout_tags + param->sp_name[0] + "_tags_obs_"  + ostr.str() + ".txt";
					//file_out = "./" + param->sp_name[0] + "_releases_obs_"  + ostr.str() + ".txt";
					wtxt.open(file_out.c_str(), ios::out);
					if (wtxt){
						wtxt << xlon << endl;
						wtxt << ylat << endl;
						wtxt << trans(rec_obs(p)) << endl;
					}
					wtxt.close();
			
				}
			}
		}
	}

	if ((t_count>t_count_rec[0]) && (month==3 || month==6 || month==9 || month==12)){
		//spatial 2d
		const double ww = 0.0001;

		taglike += ww*value(norm2(rec_obs_like-rec_pred_like));
		likelihood += ww*norm2(rec_obs_like-rec_pred_like);

		//1d (by lontigude and by latitude)	
		dvector obs_lon = rowsum(rec_obs_like);
		dvector obs_lat = colsum(rec_obs_like);
		dvar_vector pred_lon = rowsum(rec_pred_like);
		dvar_vector pred_lat = colsum(rec_pred_like);	
		taglike += 0.5*ww*(value(norm2(obs_lon-pred_lon)+norm2(obs_lat-pred_lat)));
		likelihood += 0.5*ww*(norm2(obs_lon-pred_lon)+norm2(obs_lat-pred_lat));		
	}

	return taglike;
}

dvariable LFlike_robust(const d3_array LF_qtr_obs, dvar3_array& dvarLF_est, const int a0, const int nb_ages, const int nb_regions_sp, const int k)
{//for explanation of this robustified likelihood see MFCL documentation
	dvariable likelihood = 0;
	const double PL      = 3.0;
	const double PLconst = 350.0;   // PL*PLconst should be > 1000 (= maximal sample size!)
	int I; 				// number of bins with data contributes to the likelihood weights
	double inv_I;

	dvector lf_obs(a0,nb_ages-1);
	dvar_vector lf_est(a0,nb_ages-1);
	lf_obs.initialize();
	lf_est.initialize();
	for (int r=0; r<nb_regions_sp; r++){
		double sum_lf_obs = 0;
		dvariable sum_lf_est = 0;
		int Nobs = 0;
		for (int age=a0; age<nb_ages; age++){
			lf_obs(age) = LF_qtr_obs(age,k,r);
			lf_est(age) = dvarLF_est(age,k,r);
			sum_lf_obs += lf_obs(age);
			if (lf_obs(age)) Nobs++;
			sum_lf_est += lf_est(age); 
		}
		if (Nobs && sum_lf_est>0){
			double tau = PL/sum_lf_obs;
			I = Nobs; inv_I = 1.0/I;
			likelihood += I*log(PLconst*tau);
			for (int a=a0; a<nb_ages; a++){
			//for (int a=a0; a<nb_ages-1; a++){//not taking the last A+ cohort in the likelihood to avoid bias!
				lf_obs(a) /= sum_lf_obs;
				lf_est(a) /= sum_lf_est;
				double ksi = lf_obs(a)*(1.0-lf_obs(a));

				if (lf_est(a) != 0)
					likelihood += 0.5*log(twopi*(ksi+inv_I))+
						      pow(lf_obs(a)-lf_est(a),2.0)/(2.0*tau*tau*(ksi+inv_I));
			}
		}
	}	

	return(likelihood);
}

void mat2vect(const PMap& map, const dmatrix effort, const dmatrix catch_obs, dvar_matrix catch_est, dvector& data_obs, dvar_vector& data_est, bool cpue, double c_units, const int catch_removal_fishery, const int cdata_temp_reso, const int month)
{
	if ((cdata_temp_reso==1) || 
			((cdata_temp_reso==3) && (month==1 || month==4 || month==7 || month==10))){
		data_obs.initialize();
		data_est.initialize();
	}

	int n=0;
	const int imin = map.imin; 
	const int imax = map.imax;
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				if (catch_removal_fishery){
					if (catch_obs(i,j)){
						data_obs[n] = c_units*catch_obs(i,j);
						data_est[n] = c_units*catch_est(i,j);
						n++;
					}
				} else {
					if (effort(i,j)){
						if (!cpue){
							data_obs[n] = c_units*catch_obs(i,j);
							data_est[n] = c_units*catch_est(i,j);
							n++;
						} else {
							data_obs[n] = c_units*catch_obs(i,j)/effort(i,j);
							data_est[n] = c_units*catch_est(i,j)/effort(i,j);
							n++;
						}
					}
				}
			}
		}
	}	
}

void LeastSquares(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood)
{
	likelihood = norm2(data_obs-data_est);
}


void Concentrated(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood, const int nobs)
{
	if (nobs)
		likelihood = 0.5*nobs*log(norm2(data_obs-data_est));
}

void Normal(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood, dvariable& sigma, const int nobs)
{	
	likelihood = 0.5*nobs*log(twopi*sigma*sigma) + norm2(data_obs-data_est)/(2.0*sigma*sigma);
}

void LogNormal(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood, dvariable& sigma, const int nobs)
{
	if (nobs)
		likelihood = nobs*(log(sigma) + 0.5*log(twopi)) + sum(log(data_obs)) + norm2(log(data_obs)-log(data_est)+0.5*sigma*sigma)/(2.0*sigma*sigma);
}


void ZILogNormal(const dvector& data_obs, dvar_vector& data_est, dvariable& likelihood, dvariable& sigma, dvariable& prob, const int nobs)
{

	for (int n=0; n<nobs; n++){
		const double obs = data_obs(n);
		dvariable est   = data_est(n);
		if (est>0){
			if (obs>0){
				likelihood += (1-prob)*(log(sigma) + 0.5*log(twopi) + log(obs) + pow(log(obs)-log(est)+0.5*sigma*sigma,2)/(2.0*sigma*sigma));
			} else if (obs==0)
				likelihood -= log(prob);
		}
	}
}

dvariable Poisson(const dvector data_obs, dvar_vector& data_est, const double min_obs_catch, const int nobs)
{
	dvariable likelihood = 0;
	
	for (int n=0; n<nobs; n++){
		dvariable pred = data_est(n);
		const double obs = data_obs(n);
		if (pred>0 && obs>min_obs_catch){
			likelihood += pred - obs*log(pred) + gammln(obs+1.0);
		}
	}
	return(likelihood);
}

dvariable TruncatedPoisson(const dvector data_obs, dvar_vector& data_est, const int nobs)
{
	dvariable likelihood = 0;

	for (int n=0; n<nobs; n++){
		dvariable pred = data_est(n);
		const double obs = data_obs(n);
		if (obs>1){
			if (pred>1)
				likelihood += pred - obs*log(pred) + gammln(obs+1.0)+log(1-exp(-pred));
		}
	}
	return(likelihood);
}

dvariable Exponential(const dvector data_obs, dvar_vector& data_est, const int nobs)
{
        dvariable likelihood = 0;
	double sigma = 2.0;

	for (int n=0; n<nobs; n++){
		dvariable pred = data_est(n);
		const double obs = data_obs(n);
		if (pred <= obs)
			likelihood += log(sigma) + obs/sigma - pred/sigma;
	}
        return(likelihood);
}

dvariable Weibull(const dvector data_obs, dvar_vector& data_est, const int nobs)
{
        dvariable likelihood = 0;
	double ks = 1.05;

	for (int n=0; n<nobs; n++){
		dvariable pred = data_est(n);
		const double obs = data_obs(n);
		if (obs>0)
			likelihood += -log(ks)-(ks-1)*log(obs)+ks*log(pred)+pow(obs/pred,ks);
	}
        return(likelihood);
}

dvariable NegBinomial(const dvector data_obs, dvar_vector& data_est, dvariable& beta, const int nobs)
{
	dvariable likelihood = 0;

	for (int n=0; n<nobs; n++){
		dvariable pred = data_est(n);
		const double obs = data_obs(n);
		if (pred>0)
			likelihood -= gammln(beta*pred+obs) - gammln(beta*pred) -gammln(obs+1.0) 
				+ beta*pred*log(beta)-log(beta+1.0)*(beta*pred+obs);
	}
	return(likelihood);
}	

dvariable ZINegBinomial(const dvector data_obs, dvar_vector& data_est, dvariable& beta, dvariable& p, const int nobs)
{
	dvariable likelihood = 0;

	for (int n=0; n<nobs; n++){
		dvariable pred = data_est(n);
		const double obs = data_obs(n);
		if (pred>0){
			if (obs>0){
				dvariable mu = pred/(1-p);
				likelihood -= log(1-p) + gammln(beta*mu+obs) - gammln(beta*mu) -gammln(obs+1.0) + beta*mu*log(beta)-log(beta+1.0)*(beta*mu+obs);
			} else if (obs==0){
				dvariable pwr = beta*pred/(1-p);
				likelihood -= log(p+(1-p)*pow(beta/(1.0+beta),pwr));
			}
		}
	}
	return(likelihood);
}	

