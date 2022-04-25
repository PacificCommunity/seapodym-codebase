#include "calpop.h"

///Forward main function called in simulation mode only for: 
///predicting catch by fishery without using the effort data.
///See predicted_catch_without_effort.cpp

void CCalpop::Predicted_Catch_Fishery_no_effort(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int sp, const int year, const int month)
{ 
	const int nb_fishery = param.get_nbfishery();
	const int a0 = param.sp_a0_adult[sp];
	const int nb_ages = param.sp_nb_cohorts[sp];

	int k = 0; 
	for (int f=0; f<nb_fishery; f++){ 
		if (param.mask_fishery_sp[sp][f]){
			if (param.mask_fishery_sp_no_effort[sp][f]){
		
				mat.dvarCatch_est[sp][k].initialize();
				mat.catch_est[sp][k].initialize();

				for (int age=a0; age<nb_ages; age++){

					predicted_catch_fishery_no_effort_comp(map, param, mat, f, k, sp, age);

				}
			}
			k++;
		}
	}
} 
