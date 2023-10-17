#include "NishikawaCategories.h"
#include <fvar.hpp>

// Class and function to compute the likelihood of a larvae density observed on an interval (raw Nishikawa data).
// Intervals: 0, [0,1], [1,5], [5,10], >10

NishikawaCategories::NishikawaCategories(): Nobs_cat{0, 1, 5, 10}, Nobs_diff{1, 4, 5, 40}, dl(0.1) {}

dvariable NishikawaCategories::categorical_poisson_comp(int N_obs, dvariable N_pred){
    // 1 - Computes the numerical integral of a Poisson distribution function between the two N_obs bounds
    // 2 - applies -log

    double lkhd_before_log = 0.0;
    dvariable lkhd = 0.0;
    if (N_obs==0.0){
        lkhd = N_pred;
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
    N_pred.save_prevariable_position();
    lkhd.save_prevariable_position();
    save_double_value(lkhd_before_log);
    save_double_value(value(N_pred));
    save_int_value(N_obs);
    save_identifier_string2((char*)"categorical_poisson_end");
    
    gradient_structure::GRAD_STACK1->set_gradient_stack(NishikawaCategories::dv_categorical_poisson_comp_callback);

    return lkhd;
}

void NishikawaCategories::dv_categorical_poisson_comp(){
    verify_identifier_string2((char*)"categorical_poisson_end");
    const int N_obs = restore_int_value();
    const double N_pred = restore_double_value();
    const double lkhd_before_log = restore_double_value();
    const prevariable_position lkhd_pos = restore_prevariable_position();    
    const prevariable_position N_pred_pos = restore_prevariable_position();    
    verify_identifier_string2((char*)"categorical_poisson_begin");
    
    double dflkhd = restore_prevariable_derivative(lkhd_pos);
    double dfN_pred = restore_prevariable_derivative(N_pred_pos);

    if (N_obs==0.0){
        dfN_pred += dflkhd;
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
    dflkhd = 0.0;

    save_double_derivative(dflkhd, lkhd_pos);
    save_double_derivative(dfN_pred, N_pred_pos);
}

void NishikawaCategories::dv_categorical_poisson_comp_callback() {
    // Create an instance of the class
    NishikawaCategories instance;

    // Call the non-static member function on the instance
    instance.dv_categorical_poisson_comp();
}

void NishikawaCategories::print_fields(){
    std::cerr << "dl: " << dl << std::endl;
    std::cerr << "Categories bounds: ";
    int length = sizeof(Nobs_cat) / sizeof(Nobs_cat[0]);
    for (int i=0; i < length; i++){
        std::cerr << Nobs_cat[i] << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Categories widths: ";
    for (int i=0; i < length; i++){
        std::cerr << Nobs_diff[i] << " ";
    }
    std::cerr << std::endl;
}