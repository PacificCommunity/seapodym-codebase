#include "NishikawaLike.h"
#include <fvar.hpp>

// Class and function to compute the likelihood of a larvae density observed on an interval (raw Nishikawa data).

NishikawaCategories::NishikawaCategories(VarParamCoupled& param): dl(0.1) {
    int nb_cat = param.nb_larvae_cat[0];
    for (int i = 0; i < nb_cat; i++) {
        Lobs_cat[i] = param.larvae_density_bins[0][i];
        if (i < nb_cat-1){
            Lobs_diff[i] = param.larvae_density_bins[0][i+1] - param.larvae_density_bins[0][i];
        }else{
            Lobs_diff[i] = param.larvae_density_last_bin_width[0];
        }
    }
}


dvariable NishikawaCategories::mixed_gaussian_comp(dvector& data_obs, dvar_vector& data_pred, double weight_Lobszero, VarParamCoupled& param, int sp){

    dvariable sigma = param.dvarsLikelihood_larvae_sigma[sp];
    dvariable h = param.dvarsQ_sp_larvae[sp];

    ivector data_obs_i(data_obs.indexmin(), data_obs.indexmax());
    for (int i=data_obs.indexmin(); i<=data_obs.indexmax(); i++){
        data_obs_i(i) = (int)data_obs(i);
    }
    dvariable lkhd = 0.0;
    for (int n=data_pred.indexmin(); n<data_pred.indexmax(); n++){
        const int L_obs = data_obs_i(n);
		const double L_pred = value(h)*value(data_pred(n));
        if (L_obs==0){
            value(lkhd) += weight_Lobszero*L_pred*L_pred/(2*pow(value(sigma), 2)) ;
        }else if (L_obs>0){
            if (L_pred<Lobs_cat[L_obs-1]){
                value(lkhd) += pow(L_pred - Lobs_cat[L_obs-1], 2)/(2*pow(value(sigma), 2));
            }else if (L_obs!=nb_cat && L_pred>Lobs_cat[L_obs]){
                value(lkhd) += pow(L_pred - Lobs_cat[L_obs], 2)/(2*pow(value(sigma), 2));
            }
        }
    }

    save_identifier_string2((char*)"mixed_gaussian_begin");
    save_double_value(value(sigma));
    save_double_value(value(h));
    sigma.save_prevariable_position();
    h.save_prevariable_position();
    data_pred.save_dvar_vector_value();
    data_pred.save_dvar_vector_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Lobszero);
    data_obs_i.save_ivector_value();
    data_obs_i.save_ivector_position();
    for (int i = nb_cat - 1; i >= 0; i--) {
        save_double_value(Lobs_cat[i]);
    }    
    save_int_value(nb_cat);
    save_identifier_string2((char*)"mixed_gaussian_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_mixed_gaussian_comp);

    return lkhd;
}


void dv_mixed_gaussian_comp(){
    verify_identifier_string2((char*)"mixed_gaussian_end");
    int nb_cat = restore_int_value();
    double Lobs_cat[nb_cat];
    for (int i=0; i<nb_cat; i++){
        Lobs_cat[i] = restore_double_value();
    }
    ivector_position data_obs_pos = restore_ivector_position();
    ivector data_obs = restore_ivector_value(data_obs_pos);
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
    verify_identifier_string2((char*)"mixed_gaussian_begin");

    for (int n = data_pred.indexmin(); n <= data_pred.indexmax(); n++) {
        const int L_obs = data_obs(n);
        const double L_pred = h*data_pred(n);
        double dfL_pred = 0.0;
        if (L_obs==0){
            dfL_pred += weight_Lobszero*dflkhd*L_pred / pow(sigma, 2);
            dfsigma -= weight_Lobszero*dflkhd*pow(L_pred, 2) / pow(sigma, 3);
        }else if (L_obs>0){
            if (L_pred<Lobs_cat[L_obs-1]){
                //lkhd = ((L_pred - Lobs_cat[L_obs-1])**2)/(2*sigma**2)
                dfL_pred += dflkhd*(L_pred - Lobs_cat[L_obs-1]) / pow(sigma, 2);
                dfsigma -= dflkhd*pow(L_pred - Lobs_cat[L_obs-1], 2) / pow(sigma, 3);
            }else if (L_obs!=nb_cat && L_pred>Lobs_cat[L_obs]){
                //lkhd = ((L_pred - Lobs_cat[L_obs])**2)/(2*sigma**2)
                dfL_pred += dflkhd*(L_pred - Lobs_cat[L_obs]) / pow(sigma, 2);
                dfsigma -= dflkhd*pow(L_pred - Lobs_cat[L_obs], 2) / pow(sigma, 3);
            }
        }
        dfh += dfL_pred * data_pred(n);
        dfdata_pred(n) += dfL_pred * h;
    }
    dflkhd = 0.0;

    save_double_derivative(dfsigma, sigma_pos);
    save_double_derivative(dfh, h_pos);
    dfdata_pred.save_dvector_derivatives(data_pred_pos);
    save_double_derivative(dflkhd, lkhd_pos);
}

