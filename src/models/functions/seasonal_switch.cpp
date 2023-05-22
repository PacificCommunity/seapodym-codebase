#include "VarSimtunaFunc.h"

///Forward functions for:
///computing the seasonal switch function used to switch between habitats.


double daylength_comp(double lat, double jday, double pi);
double daylength_twilight(double lat, double jday, double pi);
double tetafunc(double arg, const double teta);


void VarSimtunaFunc::Seasonal_switch_year_precomp(CParam& param, CMatrices& mat, const PMap& map, double season_peak, double season_start, const int sp)
{

	const double pi = 4.0*(atan(0.5)+atan(1.0/3.0));
	const int jday_maxDL = 171; //jday of maximal DL
	const int ndays = 366;
	const double tau = 50.0;


	for (int i=1; i<=ndays; i++){
		double Jday = i-1  + jday_maxDL - season_peak;
		double lat = param.lastlat(map.jmin);
		for (int j=map.jmin;j<=map.jmax;j++){

			lat = param.lastlat(j);
			double DL  = daylength_comp(lat, Jday, pi)/24.0;
			double dn_ratio = DL/(1.0-DL);

			//when the value becomes 1 switches to spawning habitat
			mat.season_switch(sp,i,j) = tetafunc(dn_ratio-season_start,tau);
		}
		for (int a=param.sp_a0_adult[sp];a<param.sp_nb_cohorts[sp];a++){
			double lat_north = 58.5;
			double lat_south = -58.5;//TEMPORAL: need to find a solution for only one hemisphere
			double peak_est = 117.0;//125;//put here the estimated peak to syncronize sigma with seasonal dynamics
			double DL_north  = daylength_comp(lat_north, i-1 + jday_maxDL-peak_est, pi)/24.0;
			double DL_south  = daylength_comp(lat_south, i-1 + jday_maxDL-peak_est, pi)/24.0;
			double dn_ratio_north = DL_north/(1.0-DL_north);
			double dn_ratio_south = DL_south/(1.0-DL_south);
			double dn_ratio_max = 3.6;
			if (dn_ratio_north > dn_ratio_max || dn_ratio_south > dn_ratio_max){// && param.seasonal_migrations[sp]){
				cout << "WARNING: dn_ratio = " << dn_ratio_north << " "<< dn_ratio_south << " > dn_ratio_max" << endl;
			    	cout << "pour latitude" << lat << endl;
				cout << "Seasonal switch will not work correctly in high latitudes..." << endl;
				cout << "exiting..." << endl;
				exit(1);
			}
			
			double T_age = Topt_at_age_comp(param,17,26,sp,a);//put hear estimated parameters for a species
			//TEST ALB PAPER: constant sigma: comment this line, set sigma_season to a constant value
			mat.sigma_season(sp,i,a) = (26-T_age);//put here spawning temperature estimate
			double sigma_season_start = mat.sigma_season(sp,i,a);


			//TEST ALB PAPER: constant sigma : comment this block of sigma_season computation
 			double season_start_fix = 1.0525;
			if (dn_ratio_north>=season_start_fix)
				mat.sigma_season(sp,i,a) = 1.0*(dn_ratio_north-(season_start_fix))/(dn_ratio_max-(season_start_fix))+
					 sigma_season_start*(1-(dn_ratio_north-(season_start_fix))/(dn_ratio_max-(season_start_fix)));
			if (dn_ratio_south>=season_start_fix)
				mat.sigma_season(sp,i,a) = 1.0*(dn_ratio_south-(season_start_fix))/(dn_ratio_max-(season_start_fix))+
					 sigma_season_start*(1-(dn_ratio_south-(season_start_fix))/(dn_ratio_max-(season_start_fix)));
					 
			
		}
	}
}

//temporal function used in the TEST: using fixed values in order to avoid dependence on control parameters
double VarSimtunaFunc::Topt_at_age_comp(CParam& param, const double teta_min, const double teta_max, const int sp, const int age)
{
	//const int a0 = param.sp_a0_adult[sp];
	//const int nb_ages = param.sp_nb_cohorts[sp];

	//const double W_max = param.weight[sp][param.sp_nb_cohorts[sp]-1];
	//const double W_age = param.weight[sp][age];
	const double L_max = param.length[sp][param.sp_nb_cohorts[sp]-1];
	const double L_age = param.length[sp][age];

	double R, teta_sp_age;
	
	double pcoef = 1.0;
	R = pow(L_age/L_max,pcoef);
	
	teta_sp_age = R * (teta_min-teta_max) + teta_max;

	return teta_sp_age;
}


void VarSimtunaFunc::Seasonal_switch_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, double season_peak, double season_start, const int jday, const int sp)
{

	//const double pi = 4.0*(atan(0.5)+atan(1.0/3.0));
	//const int jday_maxDL = 171; //jday of maximal DL
	//double Jday = jday-1 + jday_maxDL - season_peak;
	
	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
			
				mat.dvarSeasonSwitch(sp).elem_value(i,j) = mat.season_switch(sp,jday,j);
			}
		}
	}
}

double daylength_twilight(double lat, double jday, double pi)
{  // The CBM model of Forsythe et al, Ecological Modelling 80 (1995) 87-95

	//trvolution angle for the day of the year
	double theta = 0.2163108 + 2*(atan(0.9671396 * tan(0.00860*(jday-186))));

	//sun's declination angle, or the angular distance at solar noon between the 
	//Sun and the equator, from the Eartch orbit revolution angle
	double phi = asin(0.39795 * cos (theta));

	//angle between the sun position and the horizon, in degrees
	//6  - civil twilight
	//12 - nautical twilight
	//18 - astronomical twilight
	double p = 6; 

	//daylength computed according to 'p'
	double arg = (sin(pi*p/180)+sin(lat*pi/180)*sin(phi))/(cos(lat*pi/180)*cos(phi));
	if (arg>1.0)  arg = 1.0;
	if (arg<-1.0) arg = -1.0;
	double DL = 24.0-(24.0/pi)*acos(arg);

	return DL;
}


double daylength_comp(double lat, double jday, double pi)
{//see original function in SimtunaFunc class
     
	double DL = 0.0;

	double rum  = 0.0;
	double delta= 0.0;
	double codel= 0.0;
	double phi  = 0.0;
	double argu = 0.0;
	
	rum   = (jday-80.0)/365.25;
	delta = sin(rum*pi*2.0)*sin(pi*23.5/180.0);
	codel = asin(delta);
	phi   = lat*pi/180.0;
    	argu  = tan(codel)*tan(phi);
	
	if (lat>66.5 || lat<-66.5){
	//function is not differentiable for these latitudes
		argu  = min(1.,argu);
		argu  = max(-1.,argu);
	}

    	DL = 24.0-2.0*acos(argu)*180.0/(pi*15.0);

	if (lat>66.5 || lat<-66.5){
	    	DL = max(DL,0.0);
	}

	return DL;
}


double tetafunc(double arg, const double teta)
{
	return 1.0/(1.0+exp(-teta*(arg)));
}

