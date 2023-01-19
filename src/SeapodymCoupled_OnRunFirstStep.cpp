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

	int nb_pops = nb_species*(1+param->nb_tag_files);
	mat.CreateMatSpecies(map,t0, nbt, nbi, nbj, nb_pops, a0_adult, param->sp_nb_cohorts);
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
	if (!tuna_spinup) 
		RestoreDistributions(mat.nb_age_built);

	rw.init_writing(*param);

	param->fdata_rm = 0;
	if (!param->flag_no_fishing){
		//Initialization of fishery files variables
		//1. Read catch and effort data
		rw.rtxt_fishery_data(*param,map,nbt_total,jday_spinup);
		//redistribute all fishing effort to the model resolution
		rw.set_effort_rm(*param,map,nbt_total,jday_spinup);
		//put all data on the model resolution in case of MPA simulations or EEZ extractions
		
		if ((min(param->fishery_reso) < param->catch_reso) 
		     && !param->mpa_simulation && !param->nb_EEZ){
			rw.degrade_fishery_reso(*param, map,nbt_total,jday_spinup);
		}
		//read LF data file if provided
		if (param->file_frq_data[0]!=""){
			for (int sp=0; sp<nb_species; sp++)
				rw.read_frq_data(*param, map, param->save_first_yr, param->save_last_yr, sp);
		}
	}
	func.allocate_dvmatr(map.imin,map.imax,map.jinf,map.jsup);
	//TAG data reading and allocation section
	nb_tagpops = param->nb_tag_files;

	if (param->tag_like[0]){
		nb_rel.allocate(0,nb_tagpops-1);
		for (int n=0; n<nb_tagpops; n++){
			nb_rel(n).allocate(0,nbt_total-1);
			nb_rel(n).initialize();
		}
		create_tag_recaptures();
		ReadTaggingData(nb_rel, t_count_rec);
	}
	
	tags_age_habitat.allocate(0,aN_adult(0));
	tags_age_habitat.initialize();
	
	if (nb_tagpops>0){
		tagpop_age_solve.allocate(0,nb_tagpops-1);
		for (int p=0; p<nb_tagpops; p++){
			tagpop_age_solve(p).allocate(0,nbt_total);
			for (int n=0; n<nbt_total+1; n++){
				tagpop_age_solve(p,n).allocate(0,aN_adult(0));
				tagpop_age_solve(p,n).initialize();
			}
		}
	}
	//END of TAG data reading and allocation section
	

	//to constrain the eF parameters: their sum should hold constant
	eF_sum = 6.0;//sum(param->eF_habitat);
        //cout << "Just a WARNING: the sum of eF parameters = " << eF_sum << " will be used in the optimization" << endl;

	pop.time_reading_init();
	func.time_reading_init();
	param->time_reading_init();
}
