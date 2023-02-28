#include "VarSimtunaFunc.h"

///Forward functions for:
///feeding habitat for young and adult life stages, with or 
///without seasonal switch between habitats depending on the 
///migration flag. The function Seasonal_Habitat_Index is called
///only if the migration flag is activated.
///Note, these functions assume the fixed forage structure with 
///six functional groups. Hence the configuration of these 
///groups can be controled only through eF parameters.


double f_accessibility_layer(const double O, const double T,double twosigsq, double temp_mean, double oxy_teta, double oxy_cr);
double f_accessibility_layer(const double O2, const double T, double temp_age, double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr);


double sigma_hss_comp(double temp_age, double temp_max, const int age, const int jday);


void VarSimtunaFunc::Feeding_Habitat(VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Ha, int sp, int age, const int jday, const int t_count, const int migration_flag)
{
	
	Feeding_Habitat_Index(param,mat,map,Ha,sp,age,jday,t_count);

	if (migration_flag){
	
		dvar_matrix Hs;
		Hs.allocate(map.imin, map.imax, map.jinf, map.jsup);
		Hs.initialize();
		
		double sigma_sp_var = mat.sigma_season(sp,jday,age);//sigma_hss_comp(param.temp_age[sp][age], value(teta_max), age,jday);

		Spawning_Habitat(param,mat,map,Hs,sigma_sp_var,sp,t_count,jday);

		Seasonal_Habitat_Index(param,mat,map,Hs,Ha,sp,age,jday,t_count);
	}
}

void VarSimtunaFunc::Hf_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Ha, const int sp, const int age, const int jday, const int t)
{

	const int nb_layer = param.nb_layer;
	const int nb_forage = param.get_nbforage();
	ivector day_layer(0,nb_forage-1); day_layer = param.day_layer;
	ivector night_layer(0,nb_forage-1); night_layer = param.night_layer;

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			int nl = map.carte(i,j);
			if (nl>0 && nl<=nb_layer){
				double func_Hf = 0.0;
				double topo = map.itopo(i,j);
		
				double F = 0.0;
				double sF = 0.0;

				for (int n=0; n<nb_forage; n++){
					F = mat.forage(t,n,i,j);

					sF = param.eF_habitat[n][sp] * F;
					//accessible forage:
					func_Hf += value(mat.dvarF_access(n,age,i,j))*sF; 
				}

				//Habitat between 0 and 1
				func_Hf = param.func_limit_one(func_Hf);
				
				Ha.elem_value(i,j) = topo*func_Hf;

				//3. If not optimizing, will write habitat to DYM files
				if (!param.gcalc()){
					//last immature group if seasonal migrations, otherwise first mature
					//if (age < param.age_mature[sp])  
					if (age == 11)  
						mat.Hs[i][j] =  value(Ha(i,j));

					//always the last age computed here, will be overwritten by adult habitat
					//if seasonal migrations  
					mat.Ha[i][j] = value(Ha(i,j));	

				
				}
			}
		}
	}
}

//Inna 01/13: new function - sigma varying as a function of season peak
//Inna 05/13: the use of function is temporarily on hold, still using fixed values of sigma (see seasonal_switch)
double sigma_hss_comp(double temp_age, double temp_max, const int age, const int jday)
{
	return 0;
}


void VarSimtunaFunc::Ha_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const dmatrix Hs, dvar_matrix& Ha, const int sp, const int jday)
{
	const int nb_layer = param.nb_layer;
	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			int nl = map.carte(i,j);
			if (nl>0 && nl<=nb_layer){
				double topo = map.itopo(i,j);
				double sfunc = mat.season_switch(sp,jday,j);

				double Hf = value(Ha(i,j));
				Ha.elem_value(i,j) = topo*(Hs(i,j)*sfunc + Hf*(1.0-sfunc)) ;

				//3. If not optimizing, will write habitat to DYM files
				if (!param.gcalc()){
					//last mature group
					mat.Ha[i][j] = value(Ha(i,j));	
				}
			}
		}
	}
}

double VarSimtunaFunc::Tmean_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const int sp, const int age, const int t)
{

	const int nbl = param.nb_layer;
	dvector l_access, O2, T;
	l_access.allocate(0,nbl-1);
	O2.allocate(0,nbl-1);
	T.allocate(0,nbl-1);
	l_access.initialize();
	O2.initialize();
	T.initialize();

	double oxy_teta = param.a_oxy_habitat[sp];
	double oxy_cr   = param.b_oxy_habitat[sp];
	double temp_age = param.temp_age[sp][age];
	double sigma_ha = param.sigma_ha[sp][age];
	const double twosigsq = 2.0*sigma_ha*sigma_ha;
	const double temp_max = param.b_sst_spawning(sp);
	const double delta1   = param.thermal_func_delta[0][sp];
	const double delta2   = param.thermal_func_delta[1][sp];
	const double delta3   = param.thermal_func_delta[2][sp];

	const int Tfunc_Gaussian = param.gaussian_thermal_function[sp];

	double T_age_sum = 0.0;
	double B_age_sum = 0.0;

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			int nl = map.carte(i,j);
			if (nl>0 && nl<=nbl){


				l_access.initialize();
				for (int l=0; l<nbl; l++){
					O2(l) = mat.oxygen(t,l,i,j);
					T(l) = mat.tempn(t,l,i,j);
				}

				double l_access_sum = 0;
				double T_sum = 0;
				double T_age = temp_age;

				for (int l=0; l<nbl; l++){
					if (Tfunc_Gaussian)
						l_access(l)  = f_accessibility_layer(O2(l),T(l),twosigsq,temp_age,oxy_teta,oxy_cr);
					else
						l_access(l)  = f_accessibility_layer(O2(l),T(l),temp_age,temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);
                                        l_access_sum += l_access(l);
                                        T_sum += T(l)*l_access(l);

				}
				//only non-zero habitat, otherwise the mean will be biased
				if (l_access_sum>0.5){
					//T_age = mat.sst(t,i,j);//sst only
					T_age = T_sum/l_access_sum;
	        			T_age_sum += mat.density_after(sp,age,i,j)*T_age;
				        B_age_sum += mat.density_after(sp,age,i,j);		
				}
			}
		}
	}
        double T_mean_age = T_age_sum/B_age_sum;
	return T_mean_age;	
}

