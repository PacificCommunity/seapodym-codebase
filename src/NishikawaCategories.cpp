#include "NishikawaCategories.h"
#include <fvar.hpp>

// Class and function to compute the likelihood of a larvae density observed on an interval (raw Nishikawa data).

NishikawaCategories::NishikawaCategories(VarParamCoupled& param): dl(0.1) {
    int nb_cat = param.nb_larvae_cat;
    for (int i = 0; i < nb_cat; i++) {
        Nobs_cat[i] = param.larvae_density_bins[i];
        if (i < nb_cat-1){
            Nobs_diff[i] = param.larvae_density_bins[i+1] - param.larvae_density_bins[i];
        }else{
            Nobs_diff[i] = param.larvae_density_last_bin_width;
        }
    }
}

dvariable NishikawaCategories::categorical_poisson_comp(int N_obs, dvariable H_pred, double weight_Nobszero, VarParamCoupled& param, int sp){
    // 1 - Compute N_pred from H_pred and q_sp_larvae
    // 2 - Computes the numerical integral of a Poisson distribution function between the two N_obs bounds
    // 3 - applies -log

	const double twopi = 2.0*3.141592654;
    dvariable h = param.dvarsQ_sp_larvae[sp];
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
    dvariable h = param.dvarsQ_sp_larvae[sp];
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
    dvariable h = param.dvarsQ_sp_larvae[sp];
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


dvariable NishikawaCategories::categorical_zinb_comp(int N_obs, dvariable H_pred, VarParamCoupled& param, int sp){
    // 1 - Compute N_pred from H_pred and q_sp_larvae
    // 2 - Computes the numerical integral of a Zero-Inflated Negative Binomial function between the two N_obs bounds
    // 3 - Apply log

    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable beta = param.dvarsLikelihood_spawning_beta[sp];
    dvariable p = param.dvarsLikelihood_spawning_probzero[sp];
    dvariable N_pred = H_pred * h;

    dvariable beforelog = 0.0;
    dvariable lkhd = 0.0;
    dvariable pwr = 0.0;
    if (N_obs==0.0){
        pwr = beta*N_pred/(1-p);
        beforelog = p + (1-p) * pow(beta/(1.0+beta), pwr);
        lkhd -= log(beforelog);
        /*if (std::isnan(value(lkhd))){
            TTRACE(N_pred, N_obs)
            TTRACE(beforelog, pwr)
            TTRACE(beta, p)
            std::cerr << std::endl;
        }*/
    }else if (N_obs>0){
        // Numerical integral
        dvariable mu = N_pred/(1-p);
        int nbl;
        nbl =  (int)(Nobs_diff[N_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Nobs_cat[N_obs-1] + dl * il - dl/2;
            //beforelog += (1-p) * std::tgamma(l + beta*mu) * pow(beta/(1.0+beta), beta*mu) * pow(1/(1.0+beta), l) / (std::tgamma(beta*mu) * std::tgamma(1.0+l));
            beforelog += (1-p) * exp(gammln(l + beta*mu)) * pow(beta/(1.0+beta), beta*mu) * pow(1/(1.0+beta), l) / (exp(gammln(beta*mu)) * exp(gammln(1.0+l)));
        }
        beforelog *= dl;
        lkhd = -log(beforelog);
    }

    /*// Adjoint code: difficult to use as it needs derivative of the Gamma function
    save_identifier_string2((char*)"categorical_zinb_begin");
    H_pred.save_prevariable_position();
    beta.save_prevariable_position();
    p.save_prevariable_position();
    h.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(pwr);
    save_double_value(value(beforelog));
    save_double_value(value(N_pred));
    save_double_value(value(h));
    save_double_value(value(beta));
    save_double_value(value(p));
    save_double_value(value(H_pred));
    save_int_value(N_obs);
    for (int i=nb_cat-1; i>=0; i--){
        save_double_value(Nobs_cat[i]);
        save_double_value(Nobs_diff[i]);
    }
    save_int_value(nb_cat);
    save_double_value(dl);
    save_identifier_string2((char*)"categorical_zinb_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_categorical_zinb_comp);*/

    return lkhd;
}

/*void dv_categorical_zinb_comp(){
    verify_identifier_string2((char*)"categorical_zinb_end");
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
    const double p = restore_double_value();
    const double beta = restore_double_value();
    const double h = restore_double_value();
    const double N_pred = restore_double_value();
    const double beforelog = restore_double_value();
    const double lkhd_before_log = restore_double_value();
    const double pwr = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position beta_pos = restore_prevariable_position();    
    const prevariable_position p_pos = restore_prevariable_position();    
    const prevariable_position H_pred_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"categorical_zinb_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfbeta = restore_prevariable_derivative(beta_pos);
    double dfp = restore_prevariable_derivative(p_pos);
    double dfH_pred = restore_prevariable_derivative(H_pred_pos);

    double dfN_pred = 0.0;
    if (N_obs >= 0){
        // lkhd = -log(y)
        double dfbeforelog = -dflkhd/beforelog;

        if (N_obs==0.0){
            // beforelog = p + (1-p) * pow(beta/(1.0+beta), pwr)
            dfp += dfbeforelog * (1 - pow(beta/(1.0+beta), pwr));
            dfbeta -= dfbeforelog * (1-p) * pwr * pow(beta/(1.0+beta), pwr-1) * beta / pow(1.0+beta, 2);

            // pwr = beta * N_pred / (1-p)
            double dfpwr = dfbeforelog * (1-p) * pow(beta/(1.0+beta), pwr) * log(beta/(1.0+beta));
            dfN_pred += dfpwr * beta/(1-p);
            dfbeta += dfpwr * N_pred/(1-p);
            dfp += dfpwr / pow(1-p, 2);
        }else{
            // beforelog *= dl
            dfbeforelog *= dl;

            int nbl;
            nbl =  (int)(Nobs_diff[N_obs-1] / dl);
            for (int il = 1; il <= nbl; il++){
                double l = Nobs_cat[N_obs-1] + dl * il - dl/2;

                // beforelog += w1*w2



            }
        }
        dfh += dfN_pred * H_pred;
        dfH_pred += dfN_pred * h;
        dflkhd = 0.0;

    }

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfsigma, sigma_pos);
    save_double_derivative(dfH_pred, H_pred_pos);
}*/
