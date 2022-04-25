#include "VarSimtunaFunc.h"

///Forward functions for: 
///1) accessibility to forage components (f_accessibility) or to
///their respective layers (f_accessibility_layer).
///2) average currents given the accessibility to the layer.
///Note, these functions assume the fixed forage structure with 
///six functional groups. Hence the configuration of these 
///groups can be controled only through eF parameters.

double f_accessibility_layer(const double O2, const double T,double twosigsq, double temp_mean, double oxy_teta, double oxy_cr);
void f_accessibility(dvector& l_access, dvector& lf_access, const dvector forage, const dvector O2, const dvector T, double twosigsq, double temp_mean, double oxy_teta, double oxy_cr, const int nl, const int nb_forage, const ivector day_layer, const ivector night_layer, const double DL);

double thermal_func_type2(const double T, const double temp_min, const double temp_max, const double delta1, const double delta2, const double delta3);

double f_accessibility_layer(const double O2, const double T, double temp_age, double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr);
void f_accessibility(dvector& l_access, dvector& lf_access, const dvector forage, const dvector O2, const dvector T, double temp_mean, 
		double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr, const int nl, 
		const int nb_forage, const ivector day_layer, const ivector night_layer, const double DL);


void VarSimtunaFunc::Vars_at_age_precomp(CParam& param, const int sp)
{
	const int a0 = param.sp_a0_adult[sp];
	const int nb_ages = param.sp_nb_cohorts[sp];
	double temp_min = param.b_sst_habitat[sp];
	double temp_max = param.b_sst_spawning[sp];
	double temp_age_slope = param.T_age_size_slope[sp];
	double sigma_min = param.a_sst_spawning[sp];
	if (param.use_sst) 
		//since the estimate of sigma with sst is always lower 
		//due to larger extension of warm watermasses in the surface
		//will add 1.0 here while passing to integrated T
		sigma_min = param.a_sst_spawning[sp]+1.0;
	double sigma_max = param.a_sst_habitat[sp];

	const double W_max = param.weight[sp][param.sp_nb_cohorts[sp]-1];
	const double L_max = param.length[sp][param.sp_nb_cohorts[sp]-1];

	for (int age=a0; age<nb_ages; age++){
		const double W_age = param.weight[sp][age];
		const double L_age = param.length[sp][age];

		double R = pow(L_age/L_max,temp_age_slope);	
		double temp_sp_age = R * (temp_min-temp_max) + temp_max;

		double sigma_ha = (sigma_max-sigma_min)*W_age/W_max+sigma_min;
		//double sigma_ha = (sigma_max-sigma_min)*L_age/L_max+sigma_min;

		param.sigma_ha[sp][age] = sigma_ha;
		param.temp_age[sp][age] = temp_sp_age;
	}
}

void VarSimtunaFunc::Faccessibility_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, double temp_max, double oxy_teta, double oxy_cr, const int sp, const int age, const int jday, const int t)
{
	int Tfunc_Gaussian = param.gaussian_thermal_function[sp];
	const double delta1 = param.thermal_func_delta[0][sp];
	const double delta2 = param.thermal_func_delta[1][sp];
	const double delta3 = param.thermal_func_delta[2][sp];

	double sigma_ha = param.sigma_ha[sp][age];
	double temp_age = param.temp_age[sp][age];

	const int nb_layer = param.nb_layer;
	const int nb_forage = param.get_nbforage();
	ivector day_layer(0,nb_forage-1); day_layer = param.day_layer;
	ivector night_layer(0,nb_forage-1); night_layer = param.night_layer;

	const double twosigsq = 2.0*sigma_ha*sigma_ha;

	dvector lf_access(0,nb_layer-1);
	dvector l_access(0,nb_layer-1);
	dvector F(0,nb_forage-1);
	dvector O2(0,nb_layer-1); O2.initialize();
	dvector T(0,nb_layer-1);   T.initialize();
	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			int nl = map.carte(i,j);
			if (nl>0 && nl<=nb_layer){
				const double DL = mat.daylength[jday][j]/24.0;

				l_access.initialize();
				lf_access.initialize();
				O2.initialize();
				T.initialize();
				for (int n=0; n<nb_forage; n++)
					F(n) = mat.forage(t,n,i,j);


				for (int l=0; l<nb_layer; l++){
					O2(l) = mat.oxygen(t,l,i,j);
					T(l) = mat.tempn(t,l,i,j);
				}

				if (Tfunc_Gaussian){
					f_accessibility(l_access,lf_access,F,O2,T,twosigsq,temp_age,oxy_teta,oxy_cr,
							nb_layer,nb_forage,day_layer,night_layer,DL);
				} else {
					f_accessibility(l_access,lf_access,F,O2,T,temp_age,temp_max,delta1,delta2,delta3,
							oxy_teta,oxy_cr,nb_layer,nb_forage,day_layer,night_layer,DL);
				}

				for (int l=0; l<nb_layer; l++){
					mat.dvarZ_access(l,age).elem_value(i,j) = lf_access(l);
				}

				for (int n=0; n<nb_forage; n++){
					if (day_layer[n] < nl)
						mat.dvarF_access(n,age).elem_value(i,j) = (l_access[day_layer[n]]*DL + l_access[night_layer[n]]*(1.0-DL));
				}
			}
		}
	}
}

void VarSimtunaFunc::Average_currents_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const int age, const int t)
{//Computes tuna speed according to accessibilities to all layers

	const int nb_layer = param.nb_layer;
	double Z_access = 0.0;

	mat.u.initialize(); mat.v.initialize();
	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			int nl = map.carte(i,j);
			if (nl>0 && nl<=nb_layer){

				for (int l=0; l<nb_layer; l++){
					Z_access = value(mat.dvarZ_access(l,age,i,j));		
					mat.u(i,j) += mat.un(t,l,i,j) * Z_access;
					mat.v(i,j) += mat.vn(t,l,i,j) * Z_access;
				}
				mat.dvarsU.elem_value(i,j) = mat.u(i,j);
				mat.dvarsV.elem_value(i,j) = mat.v(i,j);
			}
		}
	}
}

double f_accessibility_layer(const double O2, const double T,double twosigsq, double temp_mean, double oxy_teta, double oxy_cr)
{

	double f_oxy = 1.0 / (1.0+ pow(oxy_teta,O2 - oxy_cr));

	double f_temp = exp(-pow(T-temp_mean,2.0) / twosigsq);

	double f_access = f_temp*f_oxy + 1e-4;
	//note, the constant is needed to avoid division by zero in weighted average calculation of lf_access
	//it does not have an impact on the absolute value of f_access

	return f_access;
}

void f_accessibility(dvector& l_access, dvector& lf_access, const dvector forage, const dvector O2, const dvector T, double twosigsq, 
		double temp_mean, double oxy_teta, double oxy_cr, const int nl, const int nb_forage, const ivector day_layer, 
		const ivector night_layer, const double DL)
{//computes accessibility functions to the vertical layer and weighted by density of forage

	lf_access.initialize();
	double sumL = 0.0;

	for (int l=0; l<nl; l++){
		l_access(l) = f_accessibility_layer(O2(l),T(l),twosigsq,temp_mean,oxy_teta,oxy_cr);
		
		// weighted by the Forage biomass
		for (int n=0;n<nb_forage;n++){
			if (day_layer[n]==l)   lf_access(l)+= l_access(l)* forage[n]*DL; 
			if (night_layer[n]==l) lf_access(l)+= l_access(l)* forage[n]*(1-DL);
		}

		sumL += lf_access(l);
	}

	for (int l=0;l<nl;l++)
		lf_access(l) = lf_access(l)/sumL;
}

double thermal_func_type2(const double T, const double temp_min, const double temp_max, const double delta1, const double delta2, const double delta3)
{//Note, the condition of delta1 and delta2 parameters: delta<0.95/max(temp_diff)	
	double temp_diff = temp_max - temp_min;
	double a1 = log(0.05+delta1*temp_diff);
	double a2 = log(0.05+delta2*temp_diff);
	double p  = -(1.0+delta3*temp_diff);
	double f_T = pow((1.0+exp(a1*(T-temp_min)))*(1.0+exp(a2*(temp_max-T))),p);

	//scaling between 0 and 1
	double xmax = (log(a2/a1)+a1*temp_min+a2*temp_max)/(a1+a2);
	f_T /= pow((1.0+exp(a1*(xmax-temp_min)))*(1.0+exp(a2*(temp_max-xmax))),p);
				
	return f_T;
}

double f_accessibility_layer(const double O2, const double T, double temp_age, double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr)
{

	double f_oxy = 1.0 / (1.0+ pow(oxy_teta,O2 - oxy_cr));

	double f_temp = thermal_func_type2(T,temp_age,temp_max,delta1,delta2,delta3);

	double f_access = f_temp*f_oxy + 1e-4;

	return f_access;
}

void f_accessibility(dvector& l_access, dvector& lf_access, const dvector forage, const dvector O2, const dvector T, double temp_age, 
		double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr, const int nl, 
		const int nb_forage, const ivector day_layer, const ivector night_layer, const double DL)
{//computes accessibility functions to the vertical layer and weighted by density of forage

	lf_access.initialize();
	double sumL = 0.0;

	for (int l=0; l<nl; l++){
		l_access(l) = f_accessibility_layer(O2(l),T(l),temp_age,temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);

		// weighted by the Forage biomass
		for (int n=0;n<nb_forage;n++){
			if (day_layer[n]==l)   lf_access(l)+= l_access(l)* forage[n]*DL; 
			if (night_layer[n]==l) lf_access(l)+= l_access(l)* forage[n]*(1-DL);
		}

		sumL += lf_access(l);
	}

	for (int l=0;l<nl;l++)
		lf_access(l) = lf_access(l)/sumL;
		//lf_access(l) = lf_access(l)/(sumL+1e-4);
}

