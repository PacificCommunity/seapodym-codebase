#include "VarSimtunaFunc.h"

///Forward functions for: 
///mortality rates at age. These functions include fixed natural mortality
///rate and variable component, depending on habitat indices defined for the life stage
double sigmoid1(const double tau, const double delta);

double sigmoid1(const double tau, const double delta)
{
        double f = 1.0/(1.0+pow(tau,delta));
        return(f);
}

void VarSimtunaFunc::M_early_sp(VarParamCoupled& param, const PMap& map, dvar_matrix& M,  const dmatrix& sst, const dmatrix& pp, const int sp)
{
	double mort_min = param.elarvae_mortality_min[sp];
	double mort_inc = param.elarvae_mortality_inc[sp];
	double sst_low  = param.elarvae_sst_low[sp];
	double sst_high = param.elarvae_sst_high[sp];
	double slp_low  = param.elarvae_slope_low[sp];
	double slp_high = param.elarvae_slope_high[sp];

	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				//1. eggs survival function derived from observations
				double f_sst  = sigmoid1(slp_low,sst(i,j)-sst_low) + sigmoid1(slp_high,sst_high-sst(i,j))-1;
				//2. early larvae mortality due to thermal factor as in Hs
				//To be added with variable parameters if proven necessary (Nov2024)   
				double f_sst2 = exp(-pow(sst(i,j)-35.0,2.0)/(2.0*40.3225));

				//3. other factors influencing early larvae survival
				double f_prey = 1.0;//no other factors

				M.elem_value(i,j) = mort_min + mort_inc*(1-f_prey*f_sst*f_sst2);
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
