#include "VarSimtunaFunc.h"

///Forward functions for: 
///mortality rates at age. These functions include fixed natural mortality
///rate and variable component, depending on habitat indices defined for the life stage


void VarSimtunaFunc::M_sp_comp(const PMap& map, dvar_matrix& M, const dmatrix& H, double Mp_max, double Ms_max, double Mp_exp, double Ms_slope, double range, double mean_age_in_dtau, const int dtau)
{
	const double Mp = Mp_max * exp(- Mp_exp * mean_age_in_dtau);
	const double Ms = Ms_max * pow(mean_age_in_dtau,Ms_slope);

	//to preserve variability of mortality for juvenile classes
	//will add nonlinear increment Rage, since for many cases (essic bet, skj, ncep/era40 alb) 
	//range=0 is estimated, which results in zero variation of mortality for juveniles
	double mean_age_in_month = mean_age_in_dtau*dtau/30.0;
	double	Rage = 1.0/(mean_age_in_month+2.5) + range; 
 	//this function gives +- 33%, 25% and 18%

	double Hval = 0.5;

	M = Mp + Ms;
	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				M.elem_value(i,j) *= pow(1.0+Rage,1-H(i,j)/Hval); //modified old function (Rage was added)

			}
		}
	}	
}

void VarSimtunaFunc::M_PH_juv_comp(VarParamCoupled& param, const PMap& map, CMatrices& mat, dvar_matrix& M, const dmatrix& PH, double mean_age_in_dtau)
{
	
	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				float ph_var = 1.0;
				if (PH(i,j)!=0){//masks are different, so zeroes are 
					double a = param.M_inc_ph_a[0]; 
					double b = param.M_inc_ph_b[0]; 
					//S1 only: c=8.29; S2,S3: c=7.5
					float  c = 7.5; 
					ph_var = 1.2*((1/pow(1+exp(-a*(c-PH(i,j))),b)-
						1/pow(1+exp(-a*(c-8.4)),b))/(1/pow(1+exp(-a*(c-6.9)),b)-
						1/pow(1+exp(-a*(c-8.4)),b)));
				}
				M.elem_value(i,j) +=  ph_var;

			}
		}
	}
}
