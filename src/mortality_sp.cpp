#include "VarSimtunaFunc.h"

///Forward functions for: 
///mortality rates at age. These functions include fixed natural mortality
///rate and variable component, depending on habitat indices defined for the life stage
double sigmoid1(const double tau, const double delta);

void VarSimtunaFunc::M_early_sp_comp(VarParamCoupled& param, const PMap& map, dvar_matrix& M,  const double a, const double b, const dmatrix& sst, const dmatrix& pp, const int sp)
{
	//parameters of empirical model: all fixed here
	double mort_inc = param.elarvae_mortality_inc[sp];
	double sst_low  = param.elarvae_sst_low[sp];
	double sst_high = param.elarvae_sst_high[sp];
	double slp_low  = param.elarvae_slope_low[sp];
	double slp_high = param.elarvae_slope_high[sp];

	//parameters to calibrate: can't be estimated from 
	//available data and highly correlated with R
	double mort_min = param.elarvae_mortality_min[sp];
	double mort_inc_larvae = param.elarvae_mortality_inc2[sp];

	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				//1. temperature-dependent early survival derived from observations
				double f_lowtemp = sigmoid1(slp_low,sst(i,j)-sst_low);
				double f_hightemp  = sigmoid1(slp_high,sst_high-sst(i,j));;
				double surv_early_stage = f_lowtemp * f_hightemp;

				//2. early larvae mortality due to thermal factor as in Hs
				//Gaussian
				//double f_sst2 = exp(-pow(sst(i,j)-b,2.0)/(2.0*a*a));
				//Logistic
				//double f_sst2 = sigmoid1(a,sst(i,j)-b);
				//Logistic thresholds, if active, use the same high-temp limit 
				//as for early stage (for eggs or hatched larvae)
				double surv_larvae = sigmoid1(a,sst(i,j)-b) * f_hightemp;

				//3. other factors influencing observed early larvae survival
				double f_prey = 1.0;//no other factors

				surv_larvae *= f_prey;

				M.elem_value(i,j) = mort_min + mort_inc*(1-surv_early_stage) + mort_inc_larvae*(1-surv_larvae);
			}
		}
	}
}


void VarSimtunaFunc::M_sp_comp(const PMap& map, dvar_matrix& M, const dmatrix& H, double Mp_max, double Ms_max, double Mp_exp, double Ms_slope, double range, const double Rage, const double Hval, const double mean_age_in_dtau)
{
	const double Mp = Mp_max * exp(- Mp_exp * mean_age_in_dtau);
	const double Ms = Ms_max * pow(mean_age_in_dtau,Ms_slope);

	M = Mp + Ms;
	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){

				M.elem_value(i,j) *= pow(1.0+Rage+range,1.0-H(i,j)/Hval);
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

void VarSimtunaFunc::Scaling_factor_sstdep_larvae_mortality_comp(const PMap& map, dvar_matrix& Scaling_factor, const dmatrix& sst, dvariable inv_M_max, dvariable inv_M_rate, dvariable age_larvae_before_sst_mortality, int deltaT){
	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				dvariable M = 1 / (inv_M_max*exp(inv_M_rate*sst(i,j)));
				Scaling_factor(i,j) = pow(1 - M, age_larvae_before_sst_mortality-deltaT);
			}
		}
	}
}
