#include "SeapodymCoupled.h"

void SeapodymCoupled::OnRunFirstStep()
{
	sumFprime.allocate(0, nb_forage - 1); 		sumFprime.initialize();
	sumF.allocate(0, nb_forage - 1);		sumF.initialize();
	sumF_area_pred.allocate(0, nb_forage - 1); 	sumF_area_pred.initialize();
	sumF_required_by_sp.allocate(0, nb_forage - 1);	sumF_required_by_sp.initialize();
	mean_omega_sp.allocate(0, nb_forage - 1);	mean_omega_sp.initialize();

	int Tr_step = 0;
        Date::init_time_variables(*param, Tr_step, nbt_spinup_tuna, jday_run, jday_spinup, nbt_total,0,0);

	month = 0;
	Date::idatymd(param->ndatini, year, month, day);
	param->set_nbt(nbt_total);
        nbt_building = nbt_spinup_tuna;
        t_count = nbt_building+1;
	//Create time-dependent forcing matrices here:
	int t0  = t_count;
	int nbt = nbt_total;
	if (!param->gcalc()) nbt = t_count;
	mat.createMatOcean(map, t0, nbt, nbi, nbj, nb_layer, deltaT);
	mat.createMatForage(map, nb_forage, t0, nbt, nbi, nbj);

	mat.CreateMatHabitat(map,nb_species,nb_forage,nb_layer,max(param->sp_nb_cohorts),t0, nbt,nbi,nbj,a0_adult,param->sp_nb_cohorts,param->age_compute_habitat);
	past_month=0;
	past_qtr=0;
	sumP = 0; 

	for (int j=map.jmin;j<=map.jmax;j++){
		double lat = param->lastlat(j);
		mat.lastlat[j] = lat;
		mat.lat_correction[j] = param->correction_lat(lat);
		for (int jd=1; jd<=366; jd++){
			//if (param->sp_name[0].find("skj")==0)
       			//	mat.daylength[jd][j] = func.daylength(lat,jd); 
			//else
				//function using twilight hours increases day length
       				mat.daylength[jd][j] = func.daylength_twilight(lat,jd,18); 
		}
	}
	
	//In Optimization Mode will read all data at once!
	if (param->gcalc())
		ReadAll();

	func.allocate_dvmatr(map.imin,map.imax,map.jinf,map.jsup);

	//to constrain the eF parameters: their sum should hold constant
	eF_sum = 6.0;//sum(param->eF_habitat);
        //cout << "Just a WARNING: the sum of eF parameters = " << eF_sum << " will be used in the optimization" << endl;

	pop.time_reading_init();
	func.time_reading_init();
	param->time_reading_init();
}
