#ifndef NISHIKAWACATEGORIES_H
#define NISHIKAWACATEGORIES_H

#include <cmath>
#include <iostream>
#include <fvar.hpp>
#include "VarSimtunaFunc.h"
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);

// Class and function to compute the likelihood of a larvae density observed on an interval (raw Nishikawa data).
// Intervals: 0, [0,1], [1,5], [5,10], >10

class NishikawaCategories {
public:
    // Constructor
    NishikawaCategories();

    // Destructor
    //~NishikawaCategories();

    dvariable categorical_poisson_comp(int N_obs, dvariable N_pred, double weight_Nobszero);// Function to compute the neg. log-likelihood for an observation on an interval
    void dv_categorical_poisson_comp();
    static void dv_categorical_poisson_comp_callback();

    dvariable mixed_gaussian_comp(int N_obs, dvariable N_pred, double weight_Nobszero, VarParamCoupled& param, int sp);// Function to compute the neg. log-likelihood for an observation on an interval
    void dv_mixed_gaussian_comp();
    static void dv_mixed_gaussian_comp_callback();
    

private:
    double Nobs_cat[4];// Bounds of intervals
    double Nobs_diff[4];// Widths of intervals
    double dl;// step to compute the numerical integrate
};
#endif // NISHIKAWACATEGORIES_H