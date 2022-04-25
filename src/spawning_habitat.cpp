#include "VarSimtunaFunc.h"

///Forward functions for: 
///spawning habitat functions. This function depends on the temperature
///in the surface layer, food for larvae and predators of larvae. 
///There is optional oxygen function, which is currently hard-coded and
///invalidated through the parameter value.


double pred_surface_comp(dvector forage, const double DL, const int nb_forage, ivector day_layer, ivector night_layer);
double hs_comp(double SST, double preys, double predators, const double a, const double b, const double c, const double d, const double e, const double ssv);
double normal(const double x, const double mu, const double sigma);
double lognormal(const double x, const double mu, const double sigma);

const double pi = 3.1415926536;

void VarSimtunaFunc::Hs_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, double a, double b, double c, double d, double e, const double ssv, const int jday, int t_count)
{
	const int nb_forage = param.get_nbforage();
	ivector dlayer   = param.day_layer;
	ivector nlayer = param.night_layer;
	dvector F(0,nb_forage-1); F.initialize();
	const double pp_transform = param.pp_transform;
	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				Hs.elem_value(i,j) = Hs_comp_elem(mat,F,pp_transform,a,b,c,d,e,ssv,nb_forage,dlayer,nlayer,jday,t_count,i,j);
			}
		}
	} 
}


double VarSimtunaFunc::Hs_comp_elem(CMatrices& mat, dvector F, const double pp_transform, const double a, const double b, const double c, const double d, const double e, const double ssv, const int nb_forage, ivector day_layer, ivector night_layer, const int jday, const int t, const int i, const int j)
{
	double Hs = 0.0;
	double pp = mat.np1(t,i,j); 
	double SST = mat.sst(t,i,j);
	double O2_l2  = mat.oxygen(t,1,i,j);

	double DL = mat.daylength(jday,j)/24.0;

	for (int n=0; n<nb_forage; n++)
		F(n) = mat.forage(t,n,i,j);		
	
	//larvae predators - micronekton at the surface
	double predators = pred_surface_comp(F,DL,nb_forage,day_layer,night_layer);

	//larvae preys - should be both phyto and zoo as a food of larvae. At the moment phyto only
	double preys = pp * pp_transform; 

	//oxygen function - for the moment static species-wise parameters
	double f_oxy = 1.0/(1.0+pow(0.01,O2_l2-0.1)); 
	//spawning habitat - product of sst, prey, predator and oxygen functions
	Hs = hs_comp(SST,preys,predators,a,b,c,d,e,ssv) * f_oxy;
	return Hs;		
}

double pred_surface_comp(dvector forage, const double DL, const int nb_forage, ivector day_layer, ivector night_layer)
{
	//predators at the surface during the day and two hours of twilight
	double frg_surf_day   = 0.0;
	double frg_surf_night = 0.0;
	for (int n=0; n<nb_forage; n++){
		if (day_layer(n)==0)   frg_surf_day   += forage(n);
		if (night_layer(n)==0) frg_surf_night += forage(n);
	}
	double pred = frg_surf_day*DL + frg_surf_night/12.0 + 1e-9; 

	return pred;		
}

double lognormal(const double x, const double mu, const double sigma)
{//20220330: function not used currently. Only the scaled version is in use, see hs_comp.
	double f_x = exp(-pow(log(x)-mu,2.0)/(2.0*sigma*sigma))/(x*sigma*sqrt(2.0*pi));
	return(f_x);
}

double normal(const double x, const double mu, const double sigma)
{
	double f_x = exp(-pow(x-mu,2.0)/(2.0*sigma*sigma));
	return(f_x);
}


double hs_comp(double SST, double preys, double predators, const double a, const double b, const double c, const double d, const double e, const double ssv)
{
	double Hs = 0.0;

	//SST function
	double f_sst = normal(SST,b,a*ssv);
	//Holling type II
	double f_prey = (preys*preys)/(c+preys*preys);
	//Optional - normal or log-normal: to do later
	//double f_pred = normal(predators,d,e);//YFT INTERIM
	//if (log-normal), devide by the value of function in the 
	//mode = exp(d-e*e) to scale it between 0 and 1 =>
	//f_pred = lognormal(predators,d,e)/lognormal(exp(d-e*e),d,e) which is:
	double f_pred = exp(d-0.5*e*e-pow(log(predators)-d,2.0)/(2.0*e*e))/predators;

	Hs = f_sst * f_prey * f_pred;

	return Hs;		
}
