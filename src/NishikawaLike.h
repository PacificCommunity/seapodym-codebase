#ifndef NISHIKAWALIKE_H
#define NISHIKAWALIKE_H

#include <cmath>
#include <iostream>
#include <fvar.hpp>
#include "VarSimtunaFunc.h"
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);

// Class and function to compute the likelihood of a larvae density observed on bins (raw Nishikawa data).

class NishikawaCategories {
public:
    // Constructor
    NishikawaCategories(VarParamCoupled& param);

    dvariable categorical_poisson_comp(int L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp);

    dvariable mixed_gaussian_comp(int L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp);

    dvariable categorical_truncated_poisson_comp(int L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp);

    dvariable categorical_zinb_comp(int L_obs, dvariable N_pred, VarParamCoupled& param, int sp);

    dvariable categorical_zip_comp(int L_obs, dvariable N_pred, VarParamCoupled& param, int sp);

private:
    static const int nb_cat = 4;
    double Lobs_cat[nb_cat];// Bounds of intervals
    double Lobs_diff[nb_cat];// Widths of intervals
    double dl;// step to compute the numerical integrate
};

void dv_categorical_poisson_comp();// gradient_structure::GRAD_STACK1->set_gradient_stack can't get a non-static member function as an argument. I would have wanted to make dv_categorical_poisson_comp() a function specific to the NishikawaCategories instance (so depending on its specific Lobs_cat, Lobs_diff, etc.), but it is not possible so I have to get this values using save_double_value and restore_double_value -> loss of memory but only solution I found 
void dv_mixed_gaussian_comp();
void dv_categorical_truncated_poisson_comp();
//void dv_categorical_zinb_comp();
void dv_categorical_zip_comp();


// Functions to compute the likelihood of a larvae density observed on a continuous scale
dvariable poisson_comp(double L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp);

dvariable gaussian_comp(double L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp);

dvariable truncated_poisson_comp(double L_obs, dvariable N_pred, double weight_Lobszero, VarParamCoupled& param, int sp);

dvariable zinb_comp(double L_obs, dvariable N_pred, VarParamCoupled& param, int sp);

dvariable zip_comp(double L_obs, dvariable N_pred, VarParamCoupled& param, int sp);

#endif // NISHIKAWALIKE_H
