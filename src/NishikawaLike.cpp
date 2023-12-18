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

dvariable NishikawaCategories::categorical_poisson_comp(int L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp){
    // 1 - Compute L_pred from N_pred and q_sp_larvae
    // 2 - Computes the numerical integral of a Poisson distribution function between the two L_obs bounds
    // 3 - applies -log

	const double twopi = 2.0*3.141592654;
    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable sigma = param.dvarsLikelihood_larvae_sigma[sp];
    dvariable L_pred = N_pred * h;

    double lkhd_before_log = 0.0;
    dvariable lkhd = 0.0;
    if (L_obs==0.0){
        lkhd = weight_Lobszero * (pow(L_pred,2) / (2*pow(sigma, 2)) + log(sigma) + log(twopi)/2);
    }else if (L_obs>0){
        // Numerical integral
        int nbl;
        nbl =  (int)(Lobs_diff[L_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Lobs_cat[L_obs-1] + dl * il - dl/2;
            lkhd += pow(L_pred, l) * exp(-L_pred) / std::tgamma(l+1);
        }
        lkhd *= dl;
        lkhd_before_log = value(lkhd);
        lkhd = -log(lkhd);
    }

    save_identifier_string2((char*)"categorical_poisson_begin");
    N_pred.save_prevariable_position();
    sigma.save_prevariable_position();
    h.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Lobszero);
    save_double_value(lkhd_before_log);
    save_double_value(value(L_pred));
    save_double_value(value(h));
    save_double_value(value(sigma));
    save_double_value(value(N_pred));
    save_int_value(L_obs);
    for (int i=nb_cat-1; i>=0; i--){
        save_double_value(Lobs_cat[i]);
        save_double_value(Lobs_diff[i]);
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
    double Lobs_cat[nb_cat];
    double Lobs_diff[nb_cat];
    for (int i=0; i<nb_cat; i++){
        Lobs_diff[i] = restore_double_value();
        Lobs_cat[i] = restore_double_value();
    }
    const int L_obs = restore_int_value();
    const double N_pred = restore_double_value();
    const double sigma = restore_double_value();
    const double h = restore_double_value();
    const double L_pred = restore_double_value();
    const double lkhd_before_log = restore_double_value();
    const double weight_Lobszero = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position sigma_pos = restore_prevariable_position();    
    const prevariable_position N_pred_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"categorical_poisson_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfsigma = restore_prevariable_derivative(sigma_pos);
    double dfN_pred = restore_prevariable_derivative(N_pred_pos);

    double dfL_pred = 0.0;
    if (L_obs==0.0){
        dfL_pred += weight_Lobszero * dflkhd * L_pred / pow(sigma, 2);
        dfsigma += weight_Lobszero * dflkhd * ((1/sigma) - pow(L_pred, 2) / pow(sigma, 3));
    }else if (L_obs>0){
        // lkhd = -log(lkhd)
        dflkhd = -dflkhd/lkhd_before_log;

        // lkhd *= dl
        dflkhd *= dl;

        int nbl;
        nbl =  (int)(Lobs_diff[L_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Lobs_cat[L_obs-1] + dl * il - dl/2;

            // w1 = pow(L_pred, l)/std::tgamma(l+1)
            double dfw1 = dflkhd * exp(-L_pred);

            // w2 = exp(-L_pred)
            double dfw2 = dflkhd * pow(L_pred, l) / std::tgamma(l+1);

            // w = w1*w2
            dfL_pred += dfw1 * l * pow(L_pred, l-1) / std::tgamma(l+1);
            dfL_pred += -dfw2 * exp(-L_pred);
        }
    }
    dfh += dfL_pred * N_pred;
    dfN_pred += dfL_pred * h;
    dflkhd = 0.0;

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfsigma, sigma_pos);
    save_double_derivative(dfN_pred, N_pred_pos);
}


dvariable NishikawaCategories::mixed_gaussian_comp(int L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp){

    dvariable sigma = param.dvarsLikelihood_larvae_sigma[sp];
    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable L_pred = N_pred * h;

    dvariable lkhd = 0.0;
    if (L_obs==0.0){
        lkhd = weight_Lobszero*L_pred*L_pred/(2*pow(sigma, 2)) ;
    }else if (L_obs>0){
        if (L_pred<Lobs_cat[L_obs-1]){
            lkhd = pow(L_pred - Lobs_cat[L_obs-1], 2)/(2*pow(sigma, 2));
        }else if (L_obs!=nb_cat && L_pred>Lobs_cat[L_obs]){
            lkhd = pow(L_pred - Lobs_cat[L_obs], 2)/(2*pow(sigma, 2));
        }
    }

    save_identifier_string2((char*)"mixed_gaussian_begin");
    sigma.save_prevariable_position();
    h.save_prevariable_position();
    N_pred.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Lobszero);
    save_double_value(value(L_pred));
    save_double_value(value(h));
    save_double_value(value(N_pred));
    save_double_value(value(sigma));
    save_int_value(L_obs);
    for (int i=nb_cat-1; i>=0; i--){
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
    const int L_obs = restore_int_value();
    const double sigma = restore_double_value();
    const double N_pred = restore_double_value();
    const double h = restore_double_value();
    const double L_pred = restore_double_value();
    const double weight_Lobszero = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position N_pred_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position sigma_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"mixed_gaussian_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfN_pred = restore_prevariable_derivative(N_pred_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfsigma = restore_prevariable_derivative(sigma_pos);

    double dfL_pred = 0.0;
    if (L_obs>=0){
        if (L_obs==0){
            dfL_pred += weight_Lobszero*dflkhd*L_pred / pow(sigma, 2);
            dfsigma += - weight_Lobszero*dflkhd*pow(L_pred, 2) / pow(sigma, 3);
        }else{
            if (L_pred<Lobs_cat[L_obs-1]){
                //lkhd = ((L_pred - Lobs_cat[L_obs-1])**2)/(2*sigma**2)
                dfL_pred += dflkhd*(L_pred - Lobs_cat[L_obs-1]) / pow(sigma, 2);
                dfsigma += - dflkhd*pow(L_pred - Lobs_cat[L_obs-1], 2) / pow(sigma, 3);
            }else if (L_obs!=nb_cat && L_pred>Lobs_cat[L_obs]){
                //lkhd = ((L_pred - Lobs_cat[L_obs])**2)/(2*sigma**2)
                dfL_pred += dflkhd*(L_pred - Lobs_cat[L_obs]) / pow(sigma, 2);
                dfsigma += - dflkhd*pow(L_pred - Lobs_cat[L_obs], 2) / pow(sigma, 3);
            }
        }

    }
    dfh += dfL_pred * N_pred;
    dfN_pred += dfL_pred * h;
    dflkhd = 0.0;

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfN_pred, N_pred_pos);
    save_double_derivative(dfsigma, sigma_pos);
}


dvariable NishikawaCategories::categorical_truncated_poisson_comp(int L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp){
    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable L_pred = 1 + N_pred * h;

    double lkhd_before_log = 0.0;
    dvariable lkhd = 0.0;
    if (L_obs==0){
        lkhd = weight_Lobszero * (L_pred - Lobs_cat[0] * log(L_pred) + gammln(Lobs_cat[0]+1.0)  + log(1-exp(-L_pred)));
    }else if (L_obs>0){
        // Numerical integral
        int nbl;
        nbl =  (int)(Lobs_diff[L_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Lobs_cat[L_obs-1] + dl * il - dl/2;
            lkhd += pow(L_pred, l) * exp(-L_pred) / (std::tgamma(l+1)*(1-exp(-L_pred)));
        }
        lkhd *= dl;
        lkhd_before_log = value(lkhd);
        lkhd = -log(lkhd);
    }

    save_identifier_string2((char*)"categorical_truncated_poisson_begin");
    N_pred.save_prevariable_position();
    h.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(weight_Lobszero);
    save_double_value(lkhd_before_log);
    save_double_value(value(L_pred));
    save_double_value(value(h));
    save_double_value(value(N_pred));
    save_int_value(L_obs);
    for (int i=nb_cat-1; i>=0; i--){
        save_double_value(Lobs_cat[i]);
        save_double_value(Lobs_diff[i]);
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
    double Lobs_cat[nb_cat];
    double Lobs_diff[nb_cat];
    for (int i=0; i<nb_cat; i++){
        Lobs_diff[i] = restore_double_value();
        Lobs_cat[i] = restore_double_value();
    }
    const int L_obs = restore_int_value();
    const double N_pred = restore_double_value();
    const double h = restore_double_value();
    const double L_pred = restore_double_value();
    const double lkhd_before_log = restore_double_value();
    const double weight_Lobszero = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position N_pred_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"categorical_truncated_poisson_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfN_pred = restore_prevariable_derivative(N_pred_pos);

    double dfL_pred = 0.0;
    if (L_obs==0){
        dfL_pred += weight_Lobszero * dflkhd * (1  - (Lobs_cat[0] / L_pred) + exp(-L_pred) / (1 - exp(-L_pred)));
    }else if (L_obs>0){
        // lkhd = -log(lkhd)
        dflkhd = -dflkhd/lkhd_before_log;

        // lkhd *= dl
        dflkhd *= dl;

        int nbl;
        nbl =  (int)(Lobs_diff[L_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Lobs_cat[L_obs-1] + dl * il - dl/2;

            // w1 = pow(L_pred, l)/std::tgamma(l+1)
            double dfw1 = dflkhd * exp(-L_pred) / (1-exp(-L_pred));

            // w2 = w3*w4
            double dfw2 = dflkhd * pow(L_pred, l) / std::tgamma(l+1.0);

            // w3 = exp(-L_pred); w4 = 1/(1-exp(-L_pred))
            double dfw3 = dfw2 / (1-exp(-L_pred));
            double dfw4 = dfw2 * exp(-L_pred);

            dfL_pred += dfw1 * l * pow(L_pred, l-1) / std::tgamma(l+1.0);
            dfL_pred -= dfw3 * exp(-L_pred);
            dfL_pred -= dfw4 * exp(-L_pred) / pow(1-exp(-L_pred), 2);

            /*// w1 = pow(L_pred, l)/std::tgamma(l+1)
            double dfw1 = dflkhd * exp(-L_pred);

            // w2 = exp(-L_pred)
            double dfw2 = dflkhd * pow(L_pred, l) / std::tgamma(l+1);

            // w = w1*w2
            dfL_pred += dfw1 * l * pow(L_pred, l-1) / std::tgamma(l+1);
            dfL_pred += -dfw2 * exp(-L_pred);*/
        }
    }
    dfh += dfL_pred * N_pred;
    dfN_pred += dfL_pred * h;
    dflkhd = 0.0;

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfN_pred, N_pred_pos);
}


dvariable NishikawaCategories::categorical_zinb_comp(int L_obs, dvariable N_pred, VarParamCoupled& param, int sp){
    // 1 - Compute L_pred from N_pred and q_sp_larvae
    // 2 - Computes the numerical integral of a Zero-Inflated Negative Binomial function between the two L_obs bounds
    // 3 - Apply log

    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable beta = param.dvarsLikelihood_larvae_beta[sp];
    dvariable p = param.dvarsLikelihood_larvae_probzero[sp];
    dvariable L_pred = N_pred * h;

    dvariable beforelog = 0.0;
    dvariable lkhd = 0.0;
    dvariable pwr = 0.0;
    if (L_obs==0.0){
        pwr = beta*L_pred/(1-p);
        beforelog = p + (1-p) * pow(beta/(1.0+beta), pwr);
        lkhd -= log(beforelog);
    }else if (L_obs>0){
        // Numerical integral
        dvariable mu = L_pred/(1-p);
        int nbl;
        nbl =  (int)(Lobs_diff[L_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Lobs_cat[L_obs-1] + dl * il - dl/2;
            beforelog += (1-p) * exp(gammln(l + beta*mu)) * pow(beta/(1.0+beta), beta*mu) * pow(1/(1.0+beta), l) / (exp(gammln(beta*mu)) * exp(gammln(1.0+l)));
        }
        beforelog *= dl;
        lkhd = -log(beforelog);
    }

    /*// Adjoint code: difficult to use as it needs derivative of the Gamma function
    save_identifier_string2((char*)"categorical_zinb_begin");
    N_pred.save_prevariable_position();
    beta.save_prevariable_position();
    p.save_prevariable_position();
    h.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(pwr);
    save_double_value(value(beforelog));
    save_double_value(value(L_pred));
    save_double_value(value(h));
    save_double_value(value(beta));
    save_double_value(value(p));
    save_double_value(value(N_pred));
    save_int_value(L_obs);
    for (int i=nb_cat-1; i>=0; i--){
        save_double_value(Lobs_cat[i]);
        save_double_value(Lobs_diff[i]);
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
    double Lobs_cat[nb_cat];
    double Lobs_diff[nb_cat];
    for (int i=0; i<nb_cat; i++){
        Lobs_diff[i] = restore_double_value();
        Lobs_cat[i] = restore_double_value();
    }
    const int L_obs = restore_int_value();
    const double N_pred = restore_double_value();
    const double p = restore_double_value();
    const double beta = restore_double_value();
    const double h = restore_double_value();
    const double L_pred = restore_double_value();
    const double beforelog = restore_double_value();
    const double lkhd_before_log = restore_double_value();
    const double pwr = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position beta_pos = restore_prevariable_position();    
    const prevariable_position p_pos = restore_prevariable_position();    
    const prevariable_position N_pred_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"categorical_zinb_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfbeta = restore_prevariable_derivative(beta_pos);
    double dfp = restore_prevariable_derivative(p_pos);
    double dfN_pred = restore_prevariable_derivative(N_pred_pos);

    double dfL_pred = 0.0;
    if (L_obs >= 0){
        // lkhd = -log(y)
        double dfbeforelog = -dflkhd/beforelog;

        if (L_obs==0.0){
            // beforelog = p + (1-p) * pow(beta/(1.0+beta), pwr)
            dfp += dfbeforelog * (1 - pow(beta/(1.0+beta), pwr));
            dfbeta -= dfbeforelog * (1-p) * pwr * pow(beta/(1.0+beta), pwr-1) * beta / pow(1.0+beta, 2);

            // pwr = beta * L_pred / (1-p)
            double dfpwr = dfbeforelog * (1-p) * pow(beta/(1.0+beta), pwr) * log(beta/(1.0+beta));
            dfL_pred += dfpwr * beta/(1-p);
            dfbeta += dfpwr * L_pred/(1-p);
            dfp += dfpwr / pow(1-p, 2);
        }else{
            // beforelog *= dl
            dfbeforelog *= dl;

            int nbl;
            nbl =  (int)(Lobs_diff[L_obs-1] / dl);
            for (int il = 1; il <= nbl; il++){
                double l = Lobs_cat[L_obs-1] + dl * il - dl/2;

                // beforelog += w1*w2



            }
        }
        dfh += dfL_pred * N_pred;
        dfN_pred += dfL_pred * h;
        dflkhd = 0.0;

    }

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfsigma, sigma_pos);
    save_double_derivative(dfN_pred, N_pred_pos);
}*/


dvariable NishikawaCategories::categorical_zip_comp(int L_obs, dvariable N_pred, VarParamCoupled& param, int sp){
    // 1 - Compute L_pred from N_pred and q_sp_larvae
    // 2 - Computes the numerical integral of a ZIPoisson distribution function between the two L_obs bounds
    // 3 - applies -log

    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable p = param.dvarsLikelihood_larvae_probzero[sp];
    dvariable L_pred = N_pred * h;

    double lkhd_before_log = 0.0;
    dvariable lkhd = 0.0;
    if (L_obs==0.0){
        lkhd = -log(p + (1-p) * exp(-L_pred));
    }else if (L_obs>0){
        // Numerical integral
        int nbl;
        nbl =  (int)(Lobs_diff[L_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Lobs_cat[L_obs-1] + dl * il - dl/2;
            lkhd += pow(L_pred, l) * exp(-L_pred) / std::tgamma(l+1);
        }
        lkhd *= (1-p) * dl;
        lkhd_before_log = value(lkhd);
        lkhd = -log(lkhd);
    }

    save_identifier_string2((char*)"categorical_zip_begin");
    N_pred.save_prevariable_position();
    p.save_prevariable_position();
    h.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(lkhd_before_log);
    save_double_value(value(L_pred));
    save_double_value(value(h));
    save_double_value(value(p));
    save_double_value(value(N_pred));
    save_int_value(L_obs);
    for (int i=nb_cat-1; i>=0; i--){
        save_double_value(Lobs_cat[i]);
        save_double_value(Lobs_diff[i]);
    }
    save_int_value(nb_cat);
    save_double_value(dl);
    save_identifier_string2((char*)"categorical_zip_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_categorical_zip_comp);

    return lkhd;
}

void dv_categorical_zip_comp(){
    verify_identifier_string2((char*)"categorical_zip_end");
    double dl = restore_double_value();
    int nb_cat = restore_int_value();
    double Lobs_cat[nb_cat];
    double Lobs_diff[nb_cat];
    for (int i=0; i<nb_cat; i++){
        Lobs_diff[i] = restore_double_value();
        Lobs_cat[i] = restore_double_value();
    }
    const int L_obs = restore_int_value();
    const double N_pred = restore_double_value();
    const double p = restore_double_value();
    const double h = restore_double_value();
    const double L_pred = restore_double_value();
    const double lkhd_before_log = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position h_pos = restore_prevariable_position();    
    const prevariable_position p_pos = restore_prevariable_position();    
    const prevariable_position N_pred_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"categorical_zip_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfh = restore_prevariable_derivative(h_pos);
    double dfp = restore_prevariable_derivative(p_pos);
    double dfN_pred = restore_prevariable_derivative(N_pred_pos);

    double dfL_pred = 0.0;
    if (L_obs==0.0){
        dfp -= dflkhd * (1-exp(-L_pred))/(p + (1-p)*exp(-L_pred));
        dfL_pred = dflkhd * (1-p) * exp(-L_pred) / (p + (1-p)*exp(-L_pred));
    }else if (L_obs>0){
        // lkhd = -log(lkhd)
        dflkhd = -dflkhd/lkhd_before_log;

        // lkhd *= (1-p)*dl
        dflkhd *= (1-p)*dl;
        dfp -= dflkhd*dl;

        int nbl;
        nbl =  (int)(Lobs_diff[L_obs-1] / dl);
        for (int il = 1; il <= nbl; il++){
            double l = Lobs_cat[L_obs-1] + dl * il - dl/2;
            dfL_pred += dflkhd * exp(-L_pred) * pow(L_pred, l-1) * (l - L_pred) / std::tgamma(l+1);
        }
    }
    dfh += dfL_pred * N_pred;
    dfN_pred += dfL_pred * h;
    dflkhd = 0.0;

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfh, h_pos);
    save_double_derivative(dfp, p_pos);
    save_double_derivative(dfN_pred, N_pred_pos);
}



// Functions to compute the likelihood of a larvae density observed on a continuous scale

dvariable gaussian_comp(double L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp){
    dvariable h = param.dvarsQ_sp_larvae[sp];
	dvariable L_pred = N_pred * h;
    dvariable sigma = param.dvarsLikelihood_larvae_sigma[sp];
    dvariable lkhd = 0.0;
    if (L_obs==0.0){
        lkhd = weight_Lobszero*L_pred*L_pred/(2*pow(sigma, 2)) ;
    }else{
        lkhd = pow(L_obs-L_pred, 2)/(2 * pow(sigma, 2));
    }
    return lkhd;
}

dvariable poisson_comp(double L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp){
	const double twopi = 2.0*3.141592654;
    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable sigma = param.dvarsLikelihood_larvae_sigma[sp];
    dvariable L_pred = N_pred * h;
    dvariable lkhd = 0.0;
    if (L_obs==0){
        lkhd = weight_Lobszero * (pow(L_pred,2) / (2*pow(sigma, 2)) + log(sigma) + log(twopi)/2);
    }else{
        lkhd = L_pred - L_obs * log(L_pred) + gammln(L_obs+1);
    }
    return lkhd;
}

dvariable truncated_poisson_comp(double L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp){
    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable L_pred = 1 + N_pred * h;
    dvariable lkhd = 0.0;
    L_obs += 1;
    lkhd = L_pred - L_obs * log(L_pred) + gammln(L_obs+1) + log(1-exp(-L_pred));
    if (L_obs==1.0){
        lkhd *= weight_Lobszero;
    }
    return lkhd;
}

dvariable zinb_comp(double L_obs, dvariable N_pred, VarParamCoupled& param, int sp){
    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable beta = param.dvarsLikelihood_larvae_beta[sp];
    dvariable p = param.dvarsLikelihood_larvae_probzero[sp];
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

dvariable zip_comp(double L_obs, dvariable N_pred, VarParamCoupled& param, int sp){
    dvariable h = param.dvarsQ_sp_larvae[sp];
    dvariable p = param.dvarsLikelihood_larvae_probzero[sp];
    dvariable L_pred = N_pred * h;
    dvariable lkhd = 0.0;
    if (L_obs==0.0){
        lkhd -= log(p + (1-p) * exp(-L_pred));
    }else{
        lkhd -= log(1-p) + L_obs * log(L_pred) - L_pred - gammln(L_obs+1.0);
    }
    return lkhd;
}

