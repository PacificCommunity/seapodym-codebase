#include "NishikawaCategories.h"
#include <fvar.hpp>

// Class and function to compute the likelihood of a larvae density observed on an interval (raw Nishikawa data).

NishikawaCategories::NishikawaCategories(VarParamCoupled& param): dl(0.1) {
    int nb_cat = param.nb_larvae_cat;
    for (int i = 0; i < nb_cat; i++) {
        Nobs_cat[i] = param.larvae_density_categories[i];
        Nobs_diff[i] = param.larvae_density_categories_width[i];
    }
}

dvariable NishikawaCategories::categorical_poisson_comp(int N_obs, dvariable H_pred, double weight_Nobszero, VarParamCoupled& param, int sp){
    // 1 - Compute N_pred from H_pred and Hs_to_larvae
    // 2 - Computes the numerical integral of a Poisson distribution function between the two N_obs bounds
    // 3 - applies -log

	const double twopi = 2.0*3.141592654;
    dvariable h = param.dvarsHs_to_larvae[sp];
    dvariable sigma = param.dvarsLikelihood_spawning_sigma[sp];
    dvariable N_pred = H_pred * h;

    double lkhd_before_log = 0.0;
    dvariable lkhd = 0.0;
    if (N_obs==0.0){
        lkhd = weight_Nobszero * (pow(N_pred,2) / (2*pow(sigma, 2)) + log(sigma) + log(twopi)/2);
    }else if (N_obs>0){
        // Numerical integral
        int nbl;
        nbl =  (int)(Nobs_diff[N_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Nobs_cat[N_obs-1] + dl * il - dl/2;
            lkhd += pow(N_pred, l) * exp(-N_pred) / std::tgamma(l+1);
        }
        lkhd *= dl;
        lkhd_before_log = value(lkhd);
        lkhd = -log(lkhd);
    }

    save_identifier_string2((char*)"categorical_poisson_begin");
    H_pred.save_prevariable_position();
    sigma.save_prevariable_position();
    h.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Nobszero);
    save_double_value(lkhd_before_log);
    save_double_value(value(N_pred));
    save_double_value(value(h));
    save_double_value(value(sigma));
    save_double_value(value(H_pred));
    save_int_value(N_obs);
    for (int i=nb_cat-1; i>=0; i--){
        save_double_value(Nobs_cat[i]);
        save_double_value(Nobs_diff[i]);
    }
    save_int_value(nb_cat);
    save_double_value(dl);
    save_identifier_string2((char*)"categorical_poisson_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_categorical_poisson_comp);

    return lkhd;
}

void dv_categorical_poisson_comp(){
    verify_identifier_string2((char*)"categorical_poisson_end");
    double dl = restore_double_value();
    int nb_cat = restore_int_value();
    double Nobs_cat[nb_cat];
    double Nobs_diff[nb_cat];
    for (int i=0; i<nb_cat; i++){
        Nobs_diff[i] = restore_double_value();
        Nobs_cat[i] = restore_double_value();
    }
    const int N_obs = restore_int_value();
    const double H_pred = restore_double_value();
    const double sigma = restore_double_value();
    const double h = restore_double_value();
    const double N_pred = restore_double_value();
    const double lkhd_before_log = restore_double_value();
    const double weight_Nobszero = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position sigma_pos = restore_prevariable_position();    
    const prevariable_position H_pred_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"categorical_poisson_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfsigma = restore_prevariable_derivative(sigma_pos);
    double dfH_pred = restore_prevariable_derivative(H_pred_pos);

    double dfN_pred = 0.0;
    if (N_obs==0.0){
        dfN_pred += weight_Nobszero * dflkhd * N_pred / pow(sigma, 2);
        dfsigma += weight_Nobszero * dflkhd * ((1/sigma) - pow(N_pred, 2) / pow(sigma, 3));
    }else if (N_obs>0){
        // lkhd = -log(lkhd)
        dflkhd = -dflkhd/lkhd_before_log;

        // lkhd *= dl
        dflkhd *= dl;

        int nbl;
        nbl =  (int)(Nobs_diff[N_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Nobs_cat[N_obs-1] + dl * il - dl/2;

            // w1 = pow(N_pred, l)/std::tgamma(l+1)
            double dfw1 = dflkhd * exp(-N_pred);

            // w2 = exp(-N_pred)
            double dfw2 = dflkhd * pow(N_pred, l) / std::tgamma(l+1);

            // w = w1*w2
            dfN_pred += dfw1 * l * pow(N_pred, l-1) / std::tgamma(l+1);
            dfN_pred += -dfw2 * exp(-N_pred);
        }
    }
    dfh += dfN_pred * H_pred;
    dfH_pred += dfN_pred * h;
    dflkhd = 0.0;

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfsigma, sigma_pos);
    save_double_derivative(dfH_pred, H_pred_pos);
}


dvariable NishikawaCategories::mixed_gaussian_comp(int N_obs, dvariable H_pred, double weight_Nobszero, VarParamCoupled& param, int sp){

    dvariable sigma = param.dvarsLikelihood_spawning_sigma[sp];
    dvariable h = param.dvarsHs_to_larvae[sp];
    dvariable N_pred = H_pred * h;

    dvariable lkhd = 0.0;
    if (N_obs==0.0){
        lkhd = weight_Nobszero*N_pred*N_pred/(2*pow(sigma, 2)) ;
    }else if (N_obs>0){
        if (N_pred<Nobs_cat[N_obs-1]){
            lkhd = pow(N_pred - Nobs_cat[N_obs-1], 2)/(2*pow(sigma, 2));
        }else if (N_obs!=nb_cat && N_pred>Nobs_cat[N_obs]){
            lkhd = pow(N_pred - Nobs_cat[N_obs], 2)/(2*pow(sigma, 2));
        }
    }

    save_identifier_string2((char*)"mixed_gaussian_begin");
    sigma.save_prevariable_position();
    h.save_prevariable_position();
    H_pred.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Nobszero);
    save_double_value(value(N_pred));
    save_double_value(value(h));
    save_double_value(value(H_pred));
    save_double_value(value(sigma));
    save_int_value(N_obs);
    for (int i=nb_cat-1; i>=0; i--){
        save_double_value(Nobs_cat[i]);
    }
    save_int_value(nb_cat);
    save_identifier_string2((char*)"mixed_gaussian_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_mixed_gaussian_comp);

    return lkhd;
}


void dv_mixed_gaussian_comp(){
    verify_identifier_string2((char*)"mixed_gaussian_end");
    int nb_cat = restore_int_value();
    double Nobs_cat[nb_cat];
    for (int i=0; i<nb_cat; i++){
        Nobs_cat[i] = restore_double_value();
    }
    const int N_obs = restore_int_value();
    const double sigma = restore_double_value();
    const double H_pred = restore_double_value();
    const double h = restore_double_value();
    const double N_pred = restore_double_value();
    const double weight_Nobszero = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position H_pred_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position sigma_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"mixed_gaussian_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfH_pred = restore_prevariable_derivative(H_pred_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfsigma = restore_prevariable_derivative(sigma_pos);

    double dfN_pred = 0.0;
    if (N_obs>=0){
        if (N_obs==0){
            dfN_pred += weight_Nobszero*dflkhd*N_pred / pow(sigma, 2);
            dfsigma += - weight_Nobszero*dflkhd*pow(N_pred, 2) / pow(sigma, 3);
        }else{
            if (N_pred<Nobs_cat[N_obs-1]){
                //lkhd = ((N_pred - Nobs_cat[N_obs-1])**2)/(2*sigma**2)
                dfN_pred += dflkhd*(N_pred - Nobs_cat[N_obs-1]) / pow(sigma, 2);
                dfsigma += - dflkhd*pow(N_pred - Nobs_cat[N_obs-1], 2) / pow(sigma, 3);
            }else if (N_obs!=nb_cat && N_pred>Nobs_cat[N_obs]){
                //lkhd = ((N_pred - Nobs_cat[N_obs])**2)/(2*sigma**2)
                dfN_pred += dflkhd*(N_pred - Nobs_cat[N_obs]) / pow(sigma, 2);
                dfsigma += - dflkhd*pow(N_pred - Nobs_cat[N_obs], 2) / pow(sigma, 3);
            }
        }

    }
    dfh += dfN_pred * H_pred;
    dfH_pred += dfN_pred * h;
    dflkhd = 0.0;

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfH_pred, H_pred_pos);
    save_double_derivative(dfsigma, sigma_pos);
}


dvariable NishikawaCategories::categorical_truncated_poisson_comp(int N_obs, dvariable H_pred, double weight_Nobszero, VarParamCoupled& param, int sp){
    dvariable h = param.dvarsHs_to_larvae[sp];
    dvariable N_pred = 1 + H_pred * h;

    double lkhd_before_log = 0.0;
    dvariable lkhd = 0.0;
    if (N_obs==0){
        lkhd = weight_Nobszero * (N_pred - Nobs_cat[0] * log(N_pred) + gammln(Nobs_cat[0]+1.0)  + log(1-exp(-N_pred)));
    }else if (N_obs>0){
        // Numerical integral
        int nbl;
        nbl =  (int)(Nobs_diff[N_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Nobs_cat[N_obs-1] + dl * il - dl/2;
            lkhd += pow(N_pred, l) * exp(-N_pred) / std::tgamma(l+1);
        }
        lkhd *= dl;
        lkhd_before_log = value(lkhd);
        lkhd = -log(lkhd);
    }

    save_identifier_string2((char*)"categorical_truncated_poisson_begin");
    H_pred.save_prevariable_position();
    h.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Nobszero);
    save_double_value(lkhd_before_log);
    save_double_value(value(N_pred));
    save_double_value(value(h));
    save_double_value(value(H_pred));
    save_int_value(N_obs);
    for (int i=nb_cat-1; i>=0; i--){
        save_double_value(Nobs_cat[i]);
        save_double_value(Nobs_diff[i]);
    }
    save_int_value(nb_cat);
    save_double_value(dl);
    save_identifier_string2((char*)"categorical_truncated_poisson_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_categorical_truncated_poisson_comp);

    return lkhd;
}

void dv_categorical_truncated_poisson_comp(){
    verify_identifier_string2((char*)"categorical_truncated_poisson_end");
    double dl = restore_double_value();
    int nb_cat = restore_int_value();
    double Nobs_cat[nb_cat];
    double Nobs_diff[nb_cat];
    for (int i=0; i<nb_cat; i++){
        Nobs_diff[i] = restore_double_value();
        Nobs_cat[i] = restore_double_value();
    }
    const int N_obs = restore_int_value();
    const double H_pred = restore_double_value();
    const double h = restore_double_value();
    const double N_pred = restore_double_value();
    const double lkhd_before_log = restore_double_value();
    const double weight_Nobszero = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position H_pred_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"categorical_truncated_poisson_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfH_pred = restore_prevariable_derivative(H_pred_pos);

    double dfN_pred = 0.0;
    if (N_obs==0){
        dfN_pred += weight_Nobszero * dflkhd * (1  - (Nobs_cat[0] / N_pred) + exp(-N_pred) / (1 - exp(-N_pred)));
    }else if (N_obs>0){
        // lkhd = -log(lkhd)
        dflkhd = -dflkhd/lkhd_before_log;

        // lkhd *= dl
        dflkhd *= dl;

        int nbl;
        nbl =  (int)(Nobs_diff[N_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Nobs_cat[N_obs-1] + dl * il - dl/2;

            // w1 = pow(N_pred, l)/std::tgamma(l+1)
            double dfw1 = dflkhd * exp(-N_pred);

            // w2 = exp(-N_pred)
            double dfw2 = dflkhd * pow(N_pred, l) / std::tgamma(l+1);

            // w = w1*w2
            dfN_pred += dfw1 * l * pow(N_pred, l-1) / std::tgamma(l+1);
            dfN_pred += -dfw2 * exp(-N_pred);
        }
    }
    dfh += dfN_pred * H_pred;
    dfH_pred += dfN_pred * h;
    dflkhd = 0.0;

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfH_pred, H_pred_pos);
}
