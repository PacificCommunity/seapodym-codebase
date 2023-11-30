#include "SeapodymCoupled.h"

void fillmatrices(const PMap& map, const dmatrix effort, const dmatrix catch_obs, dvar_matrix catch_est, dmatrix& data_obs, dvar_matrix& data_est, bool cpue, double cpue_mult, const double catch_scaling_factor);
dvariable SumSquare(const PMap& map, dvar_matrix& data_est);
void Concentrated(const dmatrix& data_obs, dvar_matrix& data_est, dvariable& likelihood, const int nobs);
void Normal(const dmatrix& data_obs, dvar_matrix& data_est, dvariable& likelihood, const double sigma_2);
void LogNormal(PMap& map, const dmatrix& data_obs, dvar_matrix& data_est, dvariable& likelihood, const int nobs, dvariable& sigma);
dvariable Poisson(PMap& map, const dmatrix data_obs, dvar_matrix& data_est);
dvariable TruncatedPoisson(PMap& map, const dmatrix data_obs, dvar_matrix& data_est);
dvariable Exponential(PMap& map, const dmatrix data_obs, dvar_matrix& data_est);
dvariable Weibull(PMap& map, const dmatrix data_obs, dvar_matrix& data_est);
double Poisson(PMap& map, const dmatrix data_obs, const dmatrix data_est);
dvariable NegBinomial(PMap& map, const dmatrix data_obs, dvar_matrix& data_est, dvariable& beta);
dvariable ZINegBinomial(PMap& map, const dmatrix data_obs, dvar_matrix& data_est, dvariable& beta, dvariable& p);
dvariable LFlike_robust(const d3_array LF_qtr_obs, dvar3_array& dvarLF_est, const int a0, const int nb_ages, const int nb_regions_sp, const int k, double& alike);
dvariable LFsum_sq(dvar3_array& dvarLF_est, const int nb_ages, const int nb_regions_sp, const int k, double& alike);
double log_factorial(const double x);

dvariable SeapodymCoupled::like(const int sp, const int k, const int f, const int nobs)
{
        double catch_mult = 1.0; // need to make it a vector of length nb_fishery 
				 // and assign values outside of this routine

	dvariable likelihood = 0.0;
	dmatrix data_obs;
	dvar_matrix data_est;
	dvariable beta, prob, sigma;
	data_obs.allocate(map.imin, map.imax, map.jinf, map.jsup);
	data_est.allocate(map.imin, map.imax, map.jinf, map.jsup);

	bool cpue = param->cpue;
	int like_type = param->like_types[sp][k];
	const double sigma_2 = pow(3.0,-2);  

	fillmatrices(map, mat.effort(f), mat.catch_obs(sp,k), mat.dvarCatch_est(sp,k), 
		     data_obs, data_est, cpue, param->cpue_mult(f),catch_mult);
 
	switch (like_type) {
	
		case 1: Normal(data_obs,data_est,likelihood,sigma_2);
			break;

		case 2: sigma = param->dvarsLike_param(sp,k);
			LogNormal(map,data_obs,data_est,likelihood,nobs,sigma);
			break;
		
		case 3: likelihood = Poisson(map,data_obs,data_est);
			break;
	
		case 4: beta = param->dvarsLike_param(sp,k);
			likelihood = NegBinomial(map,data_obs,data_est,beta);
			break;

		case 5: beta = param->dvarsLike_param(sp,k);
			prob = param->dvarsProb_zero(sp,k);
			likelihood = ZINegBinomial(map,data_obs,data_est,beta,prob);
			break;

		case 6: likelihood = TruncatedPoisson(map,data_obs,data_est);
			break;

		case 7: likelihood = SumSquare(map,data_est);
			break;

		case 8: likelihood = Exponential(map,data_obs,data_est);
			break;

		case 9: likelihood = Weibull(map,data_obs,data_est);
			break;

		case 10: Concentrated(data_obs,data_est,likelihood,nobs);
			 break;

        }
	clike_fishery[f] += value(likelihood);

	//length frequency likelihood
	if (param->frq_like[sp] && (month==3 || month==6 || month==9 || month==12)){
		const int a0 = a0_adult(sp);
		const int nb_ages = aN_adult(sp);
		const int nb_regions_sp = param->nb_region_sp_B(sp);
		if (like_type==7){//first type of sensitivity analysis
			likelihood += LFsum_sq(mat.dvarLF_est(sp),nb_ages,nb_regions_sp,k,lflike);
			return(likelihood);
		}

		//LF likelihood value for fishery 'f'
		dvariable lf_like = LFlike_robust(mat.LF_qtr_obs(sp),mat.dvarLF_est(sp),
						  a0,nb_ages,nb_regions_sp,k,lflike);
		lflike_fishery[f] += value(lf_like);

		likelihood += lf_like;
	} 
	return likelihood;
}

double SeapodymCoupled::get_stock_like(dvariable& total_stock, dvariable& likelihood)
{//returns double value of stock likelihood.

	double stocklike  = 0.0;
	for (int sp=0; sp < nb_species; sp++){
		if (!param->stock_like[sp] & !param->scalc()) return 0;
		if (param->stock_like[sp]){
			cout << "stock size: " << total_stock << endl;
			double mean_total_stock_obs = param->mean_stock_obs[sp];  
			likelihood += (total_stock-mean_total_stock_obs)*(total_stock-mean_total_stock_obs);
			stocklike += value((total_stock-mean_total_stock_obs)*(total_stock-mean_total_stock_obs));
		}
	}
	return stocklike;
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
	        mat.total_obs_catch.initialize();
		mat.total_pred_catch.initialize();
	}
	for (int p=0; p<nb_tagpops; p++){
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
				//temporally write recaptures to catch matrix (comment in CalcSums)
				for (int i = imin; i <= imax; i++){
					const int jmin = map.jinf[i];
					const int jmax = map.jsup[i];
					for (int j = jmin ; j <= jmax; j++){
						if (map.carte[i][j]){
							double xx = param->itolon(i);
							double yy = param->jtolat(j);
							if (writeoutputs){
								for (int ii=0; ii<nx_obs; ii++)
								for (int jj=0; jj<ny_obs; jj++){
									if (xx>xlon[ii] && xx<=xlon[ii+1]&& yy<=ylat[jj] && yy>ylat[jj+1]){	
										mat.total_obs_catch(0,i,j) = rec_obs(p,ii,jj)/(xr_tags*yr_tags);
										mat.total_pred_catch(0,i,j)= value(rec_pred(p,ii,jj))/(xr_tags*yr_tags);
									}
								}
							}
						}
					}
				}
				
				//TTTRACE(sum(rec_obs(p)),sum(value(rec_pred(p))),sum(mat.total_pred_catch(0)))
				mat.dvarDensity(p+1).initialize();
				tagpop_age_solve(p,t_count).initialize();

				//append to the aggregated predictions and observations
				//rec_obs_like  += elem_prod(rec_obs(p),tlib_obs(p));
				//rec_pred_like += elem_prod(rec_pred(p),tlib_obs(p));
				rec_obs_like  += rec_obs(p);
				rec_pred_like += rec_pred(p);
				//cout << norm(tlib_obs(p)) << " " << sum(rec_obs_like)<< " " << sum(value(rec_pred_like)) << endl;
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
/*				dvector obs_lon = rowsum(rec_obs(p));
				dvector obs_lat = colsum(rec_obs(p));
				dvar_vector pred_lon = rowsum(rec_pred(p));
				dvar_vector pred_lat = colsum(rec_pred(p));
		                taglike += 0.5*ww*(value(norm2(obs_lon-pred_lon)+norm2(obs_lat-pred_lat)));
		        	likelihood += 0.5*ww*(norm2(obs_lon-pred_lon)+norm2(obs_lat-pred_lat));
*/
		
				if (writeoutputs){
					//Writing only in simulation mode
					if (!param->gcalc()){
		
						ofstream wtxt;
					std::ostringstream ostr;
					ostr << year;
					if (month>9) ostr << month <<15;
					else
						ostr << 0 << month <<15;
					string file_out = "./" + param->sp_name[0] + "_tags_pred_"  + ostr.str() + ".txt";
					wtxt.open(file_out.c_str(), ios::out);
					
					if (wtxt){
						wtxt << xlon << endl;
						wtxt << ylat << endl;
						wtxt << trans(value(rec_pred(p))) << endl;
					}
					wtxt.close();
			
					file_out = "./" + param->sp_name[0] + "_tags_obs_"  + ostr.str() + ".txt";
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
		//* Note, the comments denoted '//*' is the code of SKJ taglike (verion J)
		//int nb_obs = sum(rec_obs_like);
		const float ww = 1.0; //use it if TL weight are off, but need to make it species-specific
		//const float ww = 5e-5;//20201215: increasing weight to 5e-4 for CLT experiments
		//double sf = 1.0;
		//double sf = ww*(1.0-nb_obs/(5.0+nb_obs));
		//taglike += sf*value(norm2(log(rec_obs_like+1e-4)-log(rec_pred_like+1e-4)));
//*		taglike += ww*value(norm2(rec_obs_like-rec_pred_like));
		//taglike += sf*value(norm2(log(rec_obs_like+1.0)-log(rec_pred_like+1.0)));
		//cout << sf << " " << log(rec_obs_like+1e-1) << " "<< log(rec_pred_like+1e-1)<<" " << taglike << endl;
		//likelihood += sf*norm2(log(rec_obs_like+1e-4)-log(rec_pred_like+1e-4));
//*		likelihood += ww*norm2(rec_obs_like-rec_pred_like);
		//likelihood += sf*norm2(log(rec_obs_like+1.0)-log(rec_pred_like+1.0));
		//1d (by lontigude and by latitude)
		/// comment this for e2
		
		dvector obs_lon = rowsum(rec_obs_like);
		dvector obs_lat = colsum(rec_obs_like);
		dvar_vector pred_lon = rowsum(rec_pred_like);
		dvar_vector pred_lat = colsum(rec_pred_like);	
		taglike += ww*(value(norm2(obs_lon-pred_lon)+norm2(obs_lat-pred_lat)));
		likelihood += ww*(norm2(obs_lon-pred_lon)+norm2(obs_lat-pred_lat));
		//taglike += (value(norm2(log(obs_lon+1)-log(pred_lon+1))+norm2(log(obs_lat+1)-log(pred_lat+1))));
		//likelihood += (norm2(log(obs_lon+1)-log(pred_lon+1))+norm2(log(obs_lat+1)-log(pred_lat+1)));
	}

	return taglike;
}

dvariable LFlike_robust(const d3_array LF_qtr_obs, dvar3_array& dvarLF_est, const int a0, const int nb_ages, const int nb_regions_sp, const int k, double& alike)
{//for explanation of this robustified likelihood see MFCL documentation
	dvariable likelihood = 0;
	const double twopi = 2.0*3.141592654;
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
			//for (int a=a0; a<nb_ages; a++){
			for (int a=a0; a<nb_ages-1; a++){//not taking the last A+ cohort in the likelihood to avoid bias!
				lf_obs(a) /= sum_lf_obs;
				lf_est(a) /= sum_lf_est;
				double ksi = lf_obs(a)*(1.0-lf_obs(a));

				if (lf_est(a) != 0)
					likelihood += 0.5*log(twopi*(ksi+inv_I))+
						      pow(lf_obs(a)-lf_est(a),2.0)/(2.0*tau*tau*(ksi+inv_I));
			}
		}
	}
	alike += value(likelihood);

	return(likelihood);
}

dvariable LFsum_sq(dvar3_array& dvarLF_est, const int nb_ages, const int nb_regions_sp, const int k, double& alike)
{
	dvariable likelihood = 0;

	dvar_vector lf_est(0,nb_ages-1);
	lf_est.initialize();

	for (int r=0; r<nb_regions_sp; r++){
		dvariable sum_lf_est = 0;
		for (int age=0; age<nb_ages; age++){
			lf_est(age) = dvarLF_est(age,k,r);
			sum_lf_est += lf_est(age);
		}

		if (sum_lf_est>0){
			for (int age=0; age<nb_ages; age++){
				lf_est(age) /= sum_lf_est;
			}
			likelihood += norm(lf_est);
		}
	}
	alike += value(likelihood);

	return(likelihood);

}

void fillmatrices(const PMap& map, const dmatrix effort, const dmatrix catch_obs, dvar_matrix catch_est, dmatrix& data_obs, dvar_matrix& data_est, bool cpue, double cpue_mult,const double catch_scaling_factor)
{
	data_obs.initialize();
	data_est.initialize();

	const int imin = map.imin; 
	const int imax = map.imax;

	double sf = catch_scaling_factor;
	if (cpue){
		double m = cpue_mult; // cpue multiplier
		for (int i = imin; i <= imax; i++){
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin ; j <= jmax; j++){
				if (effort(i,j)){
					data_obs(i,j) = m*catch_obs(i,j)/effort(i,j);
					data_est(i,j) = m*catch_est(i,j)/effort(i,j);
				}

			}
		}	
	} else {
		data_obs = sf*catch_obs;
		data_est = sf*catch_est;
	}
}

dvariable SumSquare(const PMap& map, dvar_matrix& data_est)
{
	dvariable likelihood = 0;

	const int imin = map.imin; 
	const int imax = map.imax; 
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				dvariable pred = data_est(i,j);
				if (pred>0)
					likelihood += pred*pred;
			}
		}
	}
	return(likelihood);
}

void Concentrated(const dmatrix& data_obs, dvar_matrix& data_est, dvariable& likelihood, const int nobs)
{
	if (nobs)
		likelihood = 0.5*nobs*log(norm2(data_obs-data_est));
}

void Normal(const dmatrix& data_obs, dvar_matrix& data_est, dvariable& likelihood, const double sigma_2)
{
	likelihood = 0.5*sigma_2*norm2(data_obs-data_est);
}


void LogNormal(PMap& map, const dmatrix& data_obs, dvar_matrix& data_est, dvariable& likelihood, const int nobs, dvariable& sigma)
{
	const int imin = map.imin; 
	const int imax = map.imax; 
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				const double obs = data_obs(i,j);
				dvariable pred   = data_est(i,j);
				likelihood += 0.5/(sigma*sigma)*pow(log(obs+1.0)-log(pred+1.0),2);
			}
		}
	}
	likelihood += nobs*log(sigma);
}

dvariable Poisson(PMap& map, const dmatrix data_obs, dvar_matrix& data_est)
{
	dvariable likelihood = 0;
	dvariable like_test = 0;
	double sf = 2.0;
	const int imin = map.imin; 
	const int imax = map.imax; 
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				dvariable pred = sf*data_est(i,j);
				const double obs = sf*data_obs(i,j);
				if (pred>0){
					if (obs>2) {
						likelihood += pred - obs*log(pred) + gammln(obs+1.0);
					}
				}
			}
		}
	}
	return(likelihood);
}

dvariable TruncatedPoisson(PMap& map, const dmatrix data_obs, dvar_matrix& data_est)
{
	dvariable likelihood = 0;
	dvariable like_test = 0;

	const int imin = map.imin; 
	const int imax = map.imax; 
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				dvariable pred = data_est(i,j);
				const double obs = data_obs(i,j);
				if (obs>1){
					if (pred>1)
						likelihood += pred - obs*log(pred) + gammln(obs+1.0)+log(1-exp(-pred));
				}
			}
		}
	}
	return(likelihood);
}

dvariable Exponential(PMap& map, const dmatrix data_obs, dvar_matrix& data_est)
{
        dvariable likelihood = 0;
	double sigma = 2.0;

        const int imin = map.imin;
        const int imax = map.imax;
        for (int i = imin; i <= imax; i++){
                const int jmin = map.jinf[i];
                const int jmax = map.jsup[i];
                for (int j = jmin ; j <= jmax; j++){
                        if (map.carte[i][j]){
                                dvariable pred = data_est(i,j);
                                const double obs = data_obs(i,j);
                                if (pred <= obs)
                                	likelihood += log(sigma) + obs/sigma - pred/sigma;
                        }
                }
        }
        return(likelihood);
}

dvariable Weibull(PMap& map, const dmatrix data_obs, dvar_matrix& data_est)
{
        dvariable likelihood = 0;
	double ks = 1.05;

        const int imin = map.imin;
        const int imax = map.imax;
        for (int i = imin; i <= imax; i++){
                const int jmin = map.jinf[i];
                const int jmax = map.jsup[i];
                for (int j = jmin ; j <= jmax; j++){
                        if (map.carte[i][j]){
                                dvariable pred = data_est(i,j);
                                const double obs = data_obs(i,j);
				if (obs>0)
                                	likelihood += -log(ks)-(ks-1)*log(obs)+ks*log(pred)+pow(obs/pred,ks);
                        }
                }
        }
        return(likelihood);
}

dvariable NegBinomial(PMap& map, const dmatrix data_obs, dvar_matrix& data_est, dvariable& beta)
{
	dvariable likelihood = 0;

	const int imin = map.imin; 
	const int imax = map.imax; 
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				dvariable pred = data_est(i,j);
				const double obs = data_obs(i,j);
				if (pred>0)
					likelihood -= gammln(beta*pred+obs) - gammln(beta*pred) -gammln(obs+1.0) 
							+ beta*pred*log(beta)-log(beta+1.0)*(beta*pred+obs);
			}
		}
	}
	return(likelihood);
}	

dvariable ZINegBinomial(PMap& map, const dmatrix data_obs, dvar_matrix& data_est, dvariable& beta, dvariable& p)
{
	dvariable likelihood = 0;
	const int imin = map.imin; 
	const int imax = map.imax; 
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				dvariable pred = data_est(i,j);
				const double obs = data_obs(i,j);
				if (pred>0){
					if (obs>0){
						dvariable mu = pred/(1-p);
						likelihood -= log(1-p) + gammln(beta*mu+obs) - gammln(beta*mu) -gammln(obs+1.0) + beta*mu*log(beta)-log(beta+1.0)*(beta*mu+obs);
					} else if (obs==0){
						dvariable pwr = beta*pred/(1-p);
						likelihood -= log(p+(1-p)*pow(beta/(1.0+beta),pwr));
					}
				} //else if (pred==0 && obs>0)
				//	likelihood -= log(1-p) + gammln(obs) - gammln(1e-10) - gammln(obs+1.0) - obs*log(beta+1.0);
				
			}
		}
	}
	return(likelihood);
}	

double Poisson(PMap& map, const dmatrix data_obs, const dmatrix data_est)
{
	double likelihood = 0;

	const int imin = map.imin; 
	const int imax = map.imax; 
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
				double pred = data_est(i,j);
				const double obs = data_obs(i,j);

				likelihood += pred - obs*log(pred+1e-10) + log_factorial(obs);
			}
		}
	}
	return(likelihood);
}


double log_factorial(const double x)
{//function works properly for values smaller than 170 only! fix later
	if (x < 1) 
		return 0;

	double fact = 1.0;
	int n = (int) x;

	if (n > 170) //beyond double precision
		n = 170;

	for (int i=1; i<=n; i++)
		fact *= i;
		
	return log(fact+1);
}

