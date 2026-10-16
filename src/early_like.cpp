#include "SeapodymCoupled.h"
#include "NishikawaLike.h"
using std::cerr;
using std::endl;

void flatten(dvector& data_obs, dvar_vector& data_pred, int nb_input_agg_groups, const std::vector<double> (*aggregated_input_vectors)[12], dvar_matrix& Agg_pred_at_obs);
void mat2vect(const PMap& map, const dmatrix input, dvar_matrix pred, dvector& data_obs, dvar_vector& data_pred);
dvariable gaussian_comp(const dvector& data_obs, dvar_vector& data_pred, double weight_Lobszero, VarParamCoupled& param, const double obs_max, int sp, string what);
dvariable poisson_comp(const dvector& data_obs, dvar_vector& data_pred, double weight_Lobszero, VarParamCoupled& param, const double obs_max, int sp, string what);
void dv_gaussian_comp();
void dv_poisson_comp();

double SeapodymCoupled::get_early_like(dvariable& likelihood, dvar_matrix& Agg_pred_at_obs, string what)
{
	NishikawaCategories NshkwCat(*param);
	double likelihood_penalty = 50.0;
	double weight_Lobszero = 0.0;
	double like_weight = 1.0;
	double like = 0.0;
	int input_categorical_flag, nb_input_agg_groups, like_type;
	std::vector<double> (*aggregated_input_vectors)[12] = nullptr;
	double obs_max = 1.0;
	if (what == "larvae"){
		if (param->fit_null_larvae[0]==1){
			weight_Lobszero = param->weight_null_larvae[0];
		}
		like_weight = param->elife_like_weight;
		like_type = param->larvae_likelihood_type[0];
		input_categorical_flag = param->larvae_input_categorical_flag[0];
		nb_input_agg_groups = param->nb_larvae_input_agg_groups;
		aggregated_input_vectors = &mat.aggregated_larvae_input_vectors;
        obs_max = larvae_obs_max;
	}else{
		weight_Lobszero = 1.0;
		like_weight = param->spawning_like_weight;
		like_type = param->spawning_likelihood_type[0];
		input_categorical_flag = 0; // categorical spawning input not implemented
		nb_input_agg_groups = param->nb_spawning_input_agg_groups;
		aggregated_input_vectors = &mat.aggregated_spawning_input_vectors;
	}

    dvector data_obs;
    dvar_vector data_pred;
    flatten(data_obs, data_pred, nb_input_agg_groups, aggregated_input_vectors, Agg_pred_at_obs);

    dvariable lkhd = 0.0;
    switch (like_type){

		case 0: // Gaussian cost function
            if (input_categorical_flag){
                lkhd = NshkwCat.mixed_gaussian_comp(data_obs, data_pred, weight_Lobszero, *param, 0);
            }else{
                lkhd = gaussian_comp(data_obs, data_pred, weight_Lobszero, *param, obs_max, 0, what);
            }
			break;

		case 1:{// Poisson cost function
            if (input_categorical_flag){
                lkhd = NshkwCat.categorical_poisson_comp(data_obs, data_pred, weight_Lobszero, *param, 0);
            }else{
                lkhd = poisson_comp(data_obs, data_pred, weight_Lobszero, *param, obs_max, 0, what);
            }
			break;
		}

		// case 2: // Truncated Poisson cost function
		// 	if (N_pred == 0.0){
		// 		if (L_obs > 0){
		// 			lkhd = likelihood_penalty;
		// 		}
		// 	}else{
		// 		lkhd = truncated_poisson_comp(L_obs, N_pred, weight_Lobszero, *param, obs_max, 0, what);
		// 	}
		// 	break;

		// case 3: // Zero-Inflated Negative Binomial cost function
		// 	lkhd = zinb_comp(L_obs, N_pred, *param, obs_max, 0, what);
		// 	break;

		// case 4: // Zero-Inflated Poisson cost function
		// 	lkhd = zip_comp(L_obs, N_pred, *param, obs_max, 0, what);
		// 	break;

		// case 5: // Logit-Normal likelihood
		// 	if (what == "spawning"){
		// 		cerr << "Logit-Normal likelihood should not be used for spawning like as observations are not probabilities." << endl;
		// 		std::exit(EXIT_FAILURE); 	

		// 	}
		// 	lkhd = logit_normal_comp(L_obs, N_pred, *param, obs_max, 0);
		// 	break;
	}

    if (input_categorical_flag){
        if (std::isinf(value(lkhd))){
            if (lkhd > 0){
                lkhd = likelihood_penalty;
            }else{
                lkhd = 0;
            }
        }
    }
    likelihood += like_weight*lkhd;
    like += like_weight*value(lkhd);

	return like;
}

double SeapodymCoupled::get_early_like(dvariable& likelihood, dvar_matrix& pred, D3_ARRAY input, int t, string what)
{//returns double value of larvae or spawning likelihood.

	NishikawaCategories NshkwCat(*param);
	double likelihood_penalty = 50.0;
	double weight_Lobszero = 0.0;
	double like  = 0.0;
	double like_weight = 1.0;
	int input_categorical_flag, like_type;
	double obs_max = 1.0;

	if (what == "larvae"){
		if (param->fit_null_larvae[0]==1){
			weight_Lobszero = param->weight_null_larvae[0];
		}
		like_weight = param->elife_like_weight;
		like_type = param->larvae_likelihood_type[0];
		input_categorical_flag = param->larvae_input_categorical_flag[0];
        obs_max = larvae_obs_max;
	}else{
		weight_Lobszero = 1.0;
		like_weight = param->spawning_like_weight;
		like_type = param->spawning_likelihood_type[0];
		input_categorical_flag = 0; // categorical spawning input not implemented
	}

    dvector data_obs;
    dvar_vector data_pred;
    mat2vect(map, input[t], pred, data_obs, data_pred);

    dvariable lkhd = 0.0;
    switch (like_type){

		case 0: // Gaussian cost function
            if (input_categorical_flag){
                lkhd = NshkwCat.mixed_gaussian_comp(data_obs, data_pred, weight_Lobszero, *param, 0);
            }else{
                lkhd = gaussian_comp(data_obs, data_pred, weight_Lobszero, *param, obs_max, 0, what);
            }
			break;

		case 1:{// Poisson cost function
            if (input_categorical_flag){
                lkhd = NshkwCat.categorical_poisson_comp(data_obs, data_pred, weight_Lobszero, *param, 0);
            }else{
                lkhd = poisson_comp(data_obs, data_pred, weight_Lobszero, *param, obs_max, 0, what);
            }
			break;
		}

		// case 2: // Truncated Poisson cost function
		// 	if (N_pred == 0.0){
		// 		if (L_obs > 0){
		// 			lkhd = likelihood_penalty;
		// 		}
		// 	}else{
		// 		lkhd = truncated_poisson_comp(L_obs, N_pred, weight_Lobszero, *param, obs_max, 0, what);
		// 	}
		// 	break;

		// case 3: // Zero-Inflated Negative Binomial cost function
		// 	lkhd = zinb_comp(L_obs, N_pred, *param, obs_max, 0, what);
		// 	break;

		// case 4: // Zero-Inflated Poisson cost function
		// 	lkhd = zip_comp(L_obs, N_pred, *param, obs_max, 0, what);
		// 	break;

		// case 5: // Logit-Normal likelihood
		// 	if (what == "spawning"){
		// 		cerr << "Logit-Normal likelihood should not be used for spawning like as observations are not probabilities." << endl;
		// 		std::exit(EXIT_FAILURE); 	

		// 	}
		// 	lkhd = logit_normal_comp(L_obs, N_pred, *param, obs_max, 0);
		// 	break;
	}

    if (input_categorical_flag){
        if (std::isinf(value(lkhd))){
            if (lkhd > 0){
                lkhd = likelihood_penalty;
            }else{
                lkhd = 0;
            }
        }
    }
    likelihood += like_weight*lkhd;
    like += like_weight*value(lkhd);

	return like;
}


void mat2vect(const PMap& map, const dmatrix input, dvar_matrix pred, dvector& data_obs, dvar_vector& data_pred)
{
	int nobs=0;
	const int imin = map.imin; 
	const int imax = map.imax;
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
                nobs++;                
			}
		}
	}	
    data_obs.allocate(0,nobs-1);
    data_pred.allocate(0,nobs-1);
    data_obs.initialize();
    data_pred.initialize();

	int n=0;
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
                data_obs[n] = input(i,j);
                data_pred[n] = pred(i,j);
                n++;                
			}
		}
	}	
}

void flatten(dvector& data_obs, dvar_vector& data_pred, int nb_input_agg_groups, const std::vector<double> (*aggregated_input_vectors)[12], dvar_matrix& Agg_pred_at_obs){

    int nobs = 0;
    for (int iAgg=0; iAgg<nb_input_agg_groups; iAgg++){
        nobs += (*aggregated_input_vectors)[iAgg].size();
    }
    if (nobs==0){
        cerr << "Error: Number of observations in early life dataset is zero." << endl;
        exit(1);
    }
    data_obs.allocate(0,nobs-1);
    data_pred.allocate(0,nobs-1);
    data_obs.initialize();
    data_pred.initialize();

    int k=0;
    for (int iAgg=0; iAgg<nb_input_agg_groups; iAgg++){
        for (auto j=0u; j<(*aggregated_input_vectors)[iAgg].size(); j++){
            data_obs(k) = (*aggregated_input_vectors)[iAgg][j];
            data_pred(k) = Agg_pred_at_obs(iAgg, j);
            k++;
        }
    }
}


static dvector obs_model_linear(double h, dvector y) { return h * y; }
static dvector obs_model_nonlinear(double h, dvector y) {
    dvector out(y.indexmin(), y.indexmax());
    for (int n=y.indexmin(); n<=y.indexmax(); n++){
        out(n) = y(n)*y(n) / (h*h + y(n)*y(n));
    }    
    return out;
} 

dvariable gaussian_comp(const dvector& data_obs, dvar_vector& data_pred, double weight_Lobszero, VarParamCoupled& param, const double obs_max, int sp, string what){

	dvariable h;
	dvariable sigma;
	int obsmodel_type;
	if (what == "larvae"){
		h = param.dvarsQ_sp_larvae[sp];
		sigma = param.dvarsLikelihood_larvae_sigma[sp];
		obsmodel_type = param.linear_larvae_obs_model;
	}else{
		h = param.dvarsQ_sp_spawning[sp];
		sigma = param.dvarsLikelihood_spawning_sigma[sp];
		obsmodel_type = 1;
	}
    dvector data_pred_transf;
    data_pred_transf.allocate(data_pred.indexmin(), data_pred.indexmax());
    data_pred_transf.initialize();
    
	if (obsmodel_type)
		data_pred_transf = obs_model_linear(value(h), value(data_pred));
	else{
		data_pred_transf = obs_max * obs_model_nonlinear(value(h), 1e-4*value(data_pred));
	}

	dvariable lkhd = 0.0;
	for (int n = data_pred_transf.indexmin(); n <= data_pred_transf.indexmax(); n++){
		const double obs = data_obs(n);
		double pred = data_pred_transf(n);
		if (obs==0)
			value(lkhd) += weight_Lobszero*pred*pred/(2.0*pow(value(sigma), 2.0));
		else
            value(lkhd) += pow(obs-pred, 2.0)/(2.0 * pow(value(sigma), 2.0));
	}
    value(lkhd) *= 1000;
	
    save_identifier_string2((char*)"gaussian_begin");
    save_double_value(obs_max);
    save_int_value(obsmodel_type);
    save_double_value(value(sigma));
    save_double_value(value(h));
    sigma.save_prevariable_position();
    h.save_prevariable_position();
    data_pred.save_dvar_vector_value();
    data_pred.save_dvar_vector_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Lobszero);
    data_obs.save_dvector_value();
    data_obs.save_dvector_position();
    save_identifier_string2((char*)"gaussian_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_gaussian_comp);

	return lkhd;
}

void dv_gaussian_comp(){
    verify_identifier_string2((char*)"gaussian_end");
    dvector_position data_obs_pos = restore_dvector_position();
    dvector data_obs = restore_dvector_value(data_obs_pos);
    const double weight_Lobszero = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    dvar_vector_position data_pred_pos = restore_dvar_vector_position();
    dvector data_pred = restore_dvar_vector_value(data_pred_pos);
    dvector dfdata_pred = restore_dvar_vector_derivatives(data_pred_pos);
    prevariable_position h_pos = restore_prevariable_position();
    double dfh = restore_prevariable_derivative(h_pos);
    prevariable_position sigma_pos = restore_prevariable_position();
    double dfsigma = restore_prevariable_derivative(sigma_pos);
    double h = restore_double_value();
    double sigma = restore_double_value();
    int obsmodel_type = restore_int_value();
    double obs_max = restore_double_value();
    verify_identifier_string2((char*)"gaussian_begin");

    dflkhd *= 1000;
    for (int n = data_pred.indexmin(); n <= data_pred.indexmax(); n++) {
        const double obs = data_obs(n);
        double D = data_pred(n);
        double pred = 0.0;
        double dfpred = 0.0;
        double dfD = 0.0;

        if (obsmodel_type){
            pred = h*D;
        }else{
            D *= 1e-4;
            pred = obs_max * D*D/(h*h+D*D);
        }

        if (obs == 0){
            //lkhd += weight_Lobszero*pred*pred/(2.0*pow(sigma, 2.0));
            dfpred += weight_Lobszero*dflkhd*pred / pow(sigma, 2);
            dfsigma -= dflkhd*weight_Lobszero*pred*pred / pow(sigma, 3);
        }else{
            //lkhd += pow(obs-pred, 2.0)/(2.0 * pow(sigma, 2.0));
            dfpred += dflkhd * (pred - obs) / pow(sigma, 2);
            dfsigma -= dflkhd * pow(pred - obs, 2) / pow(sigma, 3);
        }

        if (obsmodel_type){
            dfh += dfpred * D;
            dfdata_pred(n) += dfpred * h;
        }else{
            dfh -= 2.0*dfpred*obs_max*h*D*D/pow(h*h + D*D, 2);
            dfD += dfpred*2.0*obs_max*D*h*h/pow(h*h+D*D, 2);
            dfdata_pred(n) += dfD*1e-4;
        }
    }
    dflkhd = 0.0;

    save_double_derivative(dfsigma, sigma_pos);
    save_double_derivative(dfh, h_pos);
    dfdata_pred.save_dvector_derivatives(data_pred_pos);
    save_double_derivative(dflkhd, lkhd_pos);
}

dvariable poisson_comp(const dvector& data_obs, dvar_vector& data_pred, double weight_Lobszero, VarParamCoupled& param, const double obs_max, int sp, string what){

	dvariable h;
	dvariable sigma;
	int obsmodel_type;
	const double twopi = 2.0*3.141592654;
    double likelihood_penalty = 50.0;
	if (what == "larvae"){
		h = param.dvarsQ_sp_larvae[sp];
		sigma = param.dvarsLikelihood_larvae_sigma[sp];
		obsmodel_type = param.linear_larvae_obs_model;
	}else{
		h = param.dvarsQ_sp_spawning[sp];
		sigma = param.dvarsLikelihood_spawning_sigma[sp];
		obsmodel_type = 1;
	}
    dvector data_pred_transf;
    data_pred_transf.allocate(data_pred.indexmin(), data_pred.indexmax());
    data_pred_transf.initialize();
    
	if (obsmodel_type)
		data_pred_transf = obs_model_linear(value(h), value(data_pred));
	else{
		data_pred_transf = obs_max * obs_model_nonlinear(value(h), 1e-4*value(data_pred));
	}

	dvariable lkhd = 0.0;
	for (int n = data_pred_transf.indexmin(); n <= data_pred_transf.indexmax(); n++){
		const double obs = data_obs(n);
		double pred = data_pred_transf(n);
        if (pred == 0){
            if (obs>0)
                value(lkhd) += likelihood_penalty;
        }else{
            if (obs==0){
                value(lkhd) += weight_Lobszero * (pow(pred,2) / (2*pow(value(sigma), 2)) + log(value(sigma)) + log(twopi)/2);
            }else{
                value(lkhd) += pred - obs * log(pred) + gammln(obs+1);
            }
        }
	}

    save_identifier_string2((char*)"poisson_begin");
    save_double_value(obs_max);
    save_int_value(obsmodel_type);
    save_double_value(value(sigma));
    save_double_value(value(h));
    sigma.save_prevariable_position();
    h.save_prevariable_position();
    data_pred.save_dvar_vector_value();
    data_pred.save_dvar_vector_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Lobszero);
    data_obs.save_dvector_value();
    data_obs.save_dvector_position();
    save_identifier_string2((char*)"poisson_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_poisson_comp);

    return lkhd;
}

void dv_poisson_comp(){
    verify_identifier_string2((char*)"poisson_end");
    dvector_position data_obs_pos = restore_dvector_position();
    dvector data_obs = restore_dvector_value(data_obs_pos);
    const double weight_Lobszero = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    dvar_vector_position data_pred_pos = restore_dvar_vector_position();
    dvector data_pred = restore_dvar_vector_value(data_pred_pos);
    dvector dfdata_pred = restore_dvar_vector_derivatives(data_pred_pos);
    prevariable_position h_pos = restore_prevariable_position();
    double dfh = restore_prevariable_derivative(h_pos);
    prevariable_position sigma_pos = restore_prevariable_position();
    double dfsigma = restore_prevariable_derivative(sigma_pos);
    double h = restore_double_value();
    double sigma = restore_double_value();
    int obsmodel_type = restore_int_value();
    double obs_max = restore_double_value();
    verify_identifier_string2((char*)"poisson_begin");

    for (int n = data_pred.indexmin(); n <= data_pred.indexmax(); n++) {
        const double obs = data_obs(n);
        double D = data_pred(n);
        double pred = 0.0;
        double dfpred = 0.0;
        double dfD = 0.0;

        if (obsmodel_type){
            pred = h*D;
        }else{
            D *= 1e-4;
            pred = obs_max * D*D/(h*h+D*D);
        }

        if (pred > 0){
            if (obs==0){
                //lkhd += weight_Lobszero * (pow(pred,2) / (2*pow(sigma, 2)) + log(sigma) + log(twopi)/2);
                dfpred += dflkhd*pred*weight_Lobszero / pow(sigma, 2.0);
                dfsigma += dflkhd*weight_Lobszero*(1-pow(pred/sigma, 2.0)) / sigma;
            }else{
                //lkhd += pred - obs * log(pred) + gammln(obs+1);
                dfpred += dflkhd*(1 - obs/pred);
            }
        }
        
        if (obsmodel_type){
            dfh += dfpred * D;
            dfdata_pred(n) += dfpred * h;
        }else{
            dfh -= 2.0*dfpred*obs_max*h*D*D/pow(h*h + D*D, 2);
            dfD += dfpred*2.0*obs_max*D*h*h/pow(h*h+D*D, 2);
            dfdata_pred(n) += dfD*1e-4;
        }
    }
    dflkhd = 0.0;

    save_double_derivative(dfsigma, sigma_pos);
    save_double_derivative(dfh, h_pos);
    dfdata_pred.save_dvector_derivatives(data_pred_pos);
    save_double_derivative(dflkhd, lkhd_pos);
}
