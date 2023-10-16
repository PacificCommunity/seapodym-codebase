//#include "StdAfx.h"
#include "Matrices.h"
#include "Utilities.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//void CMatrices::create

void CMatrices::createMatHeader(const CParam& param)
//void CMatrices::createMatHeader(int nlong, int nlat, int nlevel)
{
	int nlevel = param.nlevel;
	int nlon = param.nlong;
	int nlat = param.nlat;

	zlevel = create1d(nlevel);
	xlon = create2d(nlat,nlon);
	ylat = create2d(nlat,nlon);
	mask.allocate(0, nlat - 1, 0, nlon - 1);
	mask.initialize();
//	mask = Utilities::create2d_int(nlat,nlon);
//cout << __LINE__ << "done"<<endl;

//	zlevel.allocate(0, nlevel - 1);
	
//	xlon.allocate(0, nlat - 1, 0, nlong - 1);
//	ylat.allocate(0, nlat - 1, 0, nlong - 1);

//	zlevel.initialize();
//	mask.initialize();
//	xlon.initialize();
//	ylat.initialize();

}

void CMatrices::createMatFluxes(const int nb_region, const int nb_cohort)
//void CMatrices::createMatTransport(const CParam& param, const PMap& map, int nbi, int nbj)
{
	fluxes_region.allocate(0,nb_cohort-1);
	fluxes_region.initialize();
	for (int n=0; n<nb_cohort; n++){
		fluxes_region(n).allocate(0,nb_region-1,0,nb_region-1);
		fluxes_region(n).initialize();
	}
}

void CMatrices::createMatTransport(const PMap& map)
//void CMatrices::createMatTransport(const CParam& param, const PMap& map, int nbi, int nbj)
{
	advection_x = create2d(map.imin1, map.imax1, map.jinf1, map.jsup1);
	advection_y = create2d(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x = create2d(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_y = create2d(map.imin1, map.imax1, map.jinf1, map.jsup1);
/*
	advection_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	advection_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
	diffusion_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);

	diffusion_x.initialize();
	diffusion_y.initialize();
	advection_x.initialize();
	advection_y.initialize();
*/
}

void CMatrices::createMatOcean(const PMap& map, int t0, int nbt, int nbi, int nbj, int nb_layer, int dt)
{
	//1D
	lastlat.allocate(map.jmin, map.jmax);				// Latitude of cell j
	lat_correction.allocate(map.jmin, map.jmax);			// Correction for cell area at latitude j
	//maxGD_lat.allocate(map.jmin, map.jmax);				// Correction for cell area at latitude j
	//dDL.allocate(map.jmin, map.jmax);				// Correction for cell area at latitude j

	lastlat.initialize();
	lat_correction.initialize();
	//maxGD_lat.initialize();
	//dDL.initialize();


	//2D
        int nbdays_max = 365+dt;
	daylength.allocate(1, nbdays_max, 0, nbj-1);			// Length of day based on latitude and date
	//grad_daylength.allocate(1, nbdays_max, 0, nbj-1);		// gradient of daylength;

	u.allocate(map.imin, map.imax, map.jinf, map.jsup);		// u(i,j) = advection est-ouest suivant les i
	v.allocate(map.imin, map.imax, map.jinf, map.jsup);		// v(i,j) = advection nord-sud suivant les j
	speed.allocate(map.imin, map.imax, map.jinf, map.jsup);		// magnitude of species velocity (without passive compnt);

	daylength.initialize();
	//grad_daylength.initialize();
	u.initialize();
	v.initialize();
	speed.initialize();

	
	//3D
	np1.allocate(t0, nbt);						// production primaire 1
	vld.allocate(t0, nbt);						// variable defining vertical layer depths (can be either MLD or ZEU)
	ph1.allocate(t0, nbt);
	sst.allocate(t0, nbt);						// SST
	for (int t=t0; t<=nbt; t++){
		np1(t).allocate(map.imin, map.imax, map.jinf, map.jsup);			
		vld(t).allocate(map.imin, map.imax, map.jinf, map.jsup);			
		ph1(t).allocate(map.imin, map.imax, map.jinf, map.jsup);			
		sst(t).allocate(map.imin, map.imax, map.jinf, map.jsup);			
	}
	np1.initialize();
	vld.initialize();
	ph1.initialize();
	sst.initialize();

	//4D
	un.allocate(t0, nbt);					
	vn.allocate(t0, nbt);					
	tempn.allocate(t0, nbt);					
	oxygen.allocate(t0, nbt);				
	for (int t=t0; t<=nbt; t++){
		un(t).allocate(0, nb_layer - 1);			// zonal current in layer n
		vn(t).allocate(0, nb_layer - 1);			// meridional current in layer n
		tempn(t).allocate(0, nb_layer - 1);			// temperature in layer n
		oxygen(t).allocate(0, nb_layer - 1);			// oxygen in layer n
		for (int n=0; n<nb_layer; n++){
			un(t,n).allocate(map.imin, map.imax, map.jinf, map.jsup);			
			vn(t,n).allocate(map.imin, map.imax, map.jinf, map.jsup);			
			tempn(t,n).allocate(map.imin, map.imax, map.jinf, map.jsup);			
			oxygen(t,n).allocate(map.imin, map.imax, map.jinf, map.jsup);			
		}
	}
	un.initialize();
	vn.initialize();
	tempn.initialize();
	oxygen.initialize();
}

void CMatrices::createMatSource(int nforage, int ntr, int nbi, int nbj)
{
	// ntr est le nombre de pas de temps entre le temps 0 et Tr (=Tr/timestep) 
//	source.allocate(0, nforage - 1, 0, ntr - 1, 0, nbi - 1, 0, nbj - 1);
//	source.initialize();

//	Tr.allocate(0, nforage - 1, 0, nbi - 1, 0, nbj - 1);
	mats.allocate(0, nforage - 1, 0, nbi - 1, 0, nbj - 1);

//	Tr.initialize();
	mats.initialize();
}

void CMatrices::createMatForage(const PMap& map, int nforage, int t0, int nbt, int nbi, int nbj)
{
	forage.allocate(t0, nbt);
	for (int t=t0; t<=nbt; t++){
		forage(t).allocate(0, nforage - 1);
		for (int n=0; n<nforage; n++){
			forage(t,n).allocate(map.imin, map.imax, map.jinf, map.jsup);
			forage(t,n).initialize();
		}
	}

//	nF_ratio.allocate(0, nforage - 1);
//	nF_ratio.initialize();
}

void CMatrices::createMatMortality(int nforage, int nbi, int nbj)
{
	mortality.allocate(0, nforage - 1, 0, nbi - 1, 0, nbj - 1);
	mortality.initialize();
}

void CMatrices::createMatNoBorder(int nbi, int nbj)
{
	mat2d_NoBorder.allocate(0, nbi - 1, 0, nbj - 1);
	mat2d_NoBorder.initialize();
}

void CMatrices::createMatHabitat(const PMap& map, const int nb_forage, const int nb_species,int t0, int nbt, const ivector sp_adult_age0, const ivector sp_nb_age_class, const imatrix age_compute_habitat)
{
	season_switch.allocate(0, nb_species);
	sigma_season.allocate(0, nb_species);
	//sigma_ha_season.allocate(0, nb_species);
	F_access_sum_age.allocate(0, nb_species);
	for (int sp=0; sp < nb_species; sp++){
		season_switch(sp).allocate(1, 366, map.jmin, map.jmax);
		const int a0_adult = sp_adult_age0(sp);
		const int nb_ages  = sp_nb_age_class(sp);
		sigma_season(sp).allocate(1, 366, a0_adult, nb_ages-1);//to be checked!!!
		//season_switch(sp).allocate(map.jmin, map.jmax);
		season_switch(sp).initialize();
		sigma_season(sp).initialize();
	//	sigma_ha_season(sp).initialize();

		F_access_sum_age(sp).allocate(0,nb_forage-1);
		for (int n=0; n < nb_forage; n++){
			F_access_sum_age(sp,n).allocate(t0, nbt);
			for (int t=t0; t<=nbt; t++)
				F_access_sum_age(sp,n,t).allocate(map.imin, map.imax, map.jinf, map.jsup);
		}
	}


	Hs.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Hj.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Ha.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Hs.initialize();
	Hj.initialize();
	Ha.initialize();

	//store the habitat indices for selected ages: age_compute_habitat (see Param class)
	adult_habitat.allocate(0,nb_species-1);
	for (int sp = 0; sp < nb_species; sp++){
		adult_habitat(sp).allocate(t0, nbt);	
		for (int t=t0; t<=nbt; t++){
			const int a0_adult = sp_adult_age0(sp);
			const int nb_ages  = age_compute_habitat(sp,sp_nb_age_class(sp)-1);
			adult_habitat(sp,t).allocate(a0_adult, nb_ages);
			//cout << "allocated adult_habitat for indices "<< a0_adult << " "<<  nb_ages << endl; 
			for (int a=a0_adult; a<nb_ages; a++){
				adult_habitat(sp,t,a).allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
			}
		}
	}
	adult_habitat.initialize();
}

void CMatrices::createMatHabitat_input(const PMap& map, const int nb_ages, const int nbt_total){

	habitat_input.allocate(0,nb_ages-1);
	for (int n=0; n<nb_ages; n++){
		habitat_input[n].allocate(0,nbt_total);
		for (int t=0; t<=nbt_total; t++){
			habitat_input(n,t).allocate(map.imin, map.imax, map.jinf, map.jsup);
			habitat_input(n,t).initialize();
		}
	}
}

void CMatrices::createMatSpecies(const PMap& map, int t0, int nbt, int nbi, int nbj, int nb_species, const ivector sp_adult_age0, const ivector sp_nb_age_class)
{
	//mean variables for computing model outputs
        mean_mortality.allocate(0,nb_species-1);
        mean_speed.allocate(0,nb_species-1);
        mean_diffusion.allocate(0,nb_species-1);
        mean_temperature.allocate(0,nb_species-1);
        for (int sp = 0; sp < nb_species; sp++){
                mean_mortality(sp).allocate(0,sp_nb_age_class(sp)-1);
                mean_speed(sp).allocate(0,sp_nb_age_class(sp)-1);
                mean_diffusion(sp).allocate(0,sp_nb_age_class(sp)-1);
                mean_temperature(sp).allocate(0,sp_nb_age_class(sp)-1);
        }
        mean_mortality.initialize();
        mean_speed.initialize();
        mean_diffusion.initialize();
        mean_temperature.initialize();
	//end of mean variables allocation&initializatconst PMap& mapion

	density_before.allocate(0,nb_species-1);
	for (int sp = 0; sp < nb_species; sp++){
		density_before(sp).allocate(t0, nbt);		//fish density before solution of ADRE
		for (int t=t0; t<=nbt; t++){
			const int a0_adult = sp_adult_age0(sp);
			const int nb_ages  = sp_nb_age_class(sp);
			density_before(sp,t).allocate(a0_adult, nb_ages - 1);	
			for (int a=0; a<nb_ages; a++){
				if (a>=a0_adult)
					density_before(sp,t,a).allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
			}
		}
	}
	density_before.initialize();

	density_after.allocate(0,nb_species-1);
	for (int sp = 0; sp < nb_species; sp++){
		const int nb_ages  = sp_nb_age_class(sp);
		density_after(sp).allocate(0, nb_ages - 1);
		for (int a=0; a<nb_ages; a++){
			density_after(sp,a).allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
		}
	}
	density_after.initialize();


/*
	pop_species.allocate(0, nb_species - 1);
	for (int sp = 0; sp < nb_species; sp++){
		const int agemax = sp_nb_age_class(sp);
		pop_species(sp).allocate(0, agemax - 1);
		pop_species(sp).initialize();

		for (int age = 0; age < agemax; age++){
			pop_species(sp,age).allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
			pop_species(sp, age).initialize();
		}
	}

	// juveniles - 3 first months
	juv_species.allocate(0, nb_species - 1);
	for (int sp = 0; sp < nb_species; sp++){
		const int agemax = 3;
		juv_species(sp).allocate(0, agemax - 1);
		juv_species(sp).initialize();

		for (int age = 0; age < agemax; age++){
			juv_species(sp,age).allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
			juv_species(sp, age).initialize();
		}
	}
*/
	nb_age_built.allocate(0, nb_species - 1); 
	nb_age_built.initialize();
	
	init_density_species.allocate(0, nb_species - 1);
	for (int sp = 0; sp < nb_species; sp++){
		const int nb_ages = sp_nb_age_class[sp];
		init_density_species[sp].allocate(0, nb_ages - 1);
		init_density_species[sp].initialize();
		for (int age = 0; age < nb_ages; age++){
			init_density_species(sp,age).allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
			init_density_species(sp,age).initialize();
		}
	}

	larvae.allocate(0, nb_species - 1);
	juvenile.allocate(0, nb_species - 1);
	young.allocate(0, nb_species - 1);
	recruit.allocate(0, nb_species - 1);
	adult.allocate(0, nb_species - 1);
	total_pop.allocate(0, nb_species - 1);
	for (int sp = 0; sp < nb_species; sp++){
		larvae(sp).allocate(map.imin, map.imax, map.jinf, map.jsup);
		juvenile(sp).allocate(map.imin, map.imax, map.jinf, map.jsup);
		young(sp).allocate(map.imin, map.imax, map.jinf, map.jsup);
		recruit(sp).allocate(map.imin, map.imax, map.jinf, map.jsup);
		adult(sp).allocate(map.imin, map.imax, map.jinf, map.jsup);
		total_pop(sp).allocate(map.imin, map.imax, map.jinf, map.jsup);
		
		larvae(sp).initialize();
		juvenile(sp).initialize();
		young(sp).initialize();
		recruit(sp).initialize();
		adult(sp).initialize();
		total_pop(sp).initialize();
	}
	// 1D
	// total de la biomasse des larves (0-1 mois)
	sum_B_larvae.allocate(0, nb_species - 1);
	// total de la biomasse des juveniles (1-4 mois)
	sum_B_juv.allocate(0, nb_species - 1);
	// total biomasse des jeunes recrues (5-9 mois)
	sum_B_young.allocate(0, nb_species - 1);
	// total de la biomasse en thon
	sum_B_recruit.allocate(0, nb_species - 1);
	// total biomasse du stock
	sum_B_adult.allocate(0, nb_species - 1);
	// total de la biomasse en thon
	sum_total_pop.allocate(0, nb_species - 1);

	sum_B_larvae.initialize();
	sum_B_juv.initialize();
	sum_B_young.initialize();
	sum_B_recruit.initialize();
	sum_B_adult.initialize();
	sum_total_pop.initialize();
}

void CMatrices::createMatEffort(const PMap& map, int nbi, int nbj, int nb_fleet)
{
	effort.allocate(0, nb_fleet - 1);
	efflon.allocate(0, nb_fleet - 1);
	efflat.allocate(0, nb_fleet - 1);
	for (int f=0; f<nb_fleet; f++){
		effort(f).allocate(map.imin,map.imax,map.jinf,map.jsup);
		efflon(f).allocate(map.imin,map.imax,map.jinf,map.jsup);
		efflat(f).allocate(map.imin,map.imax,map.jinf,map.jsup);
	}
	effort.initialize();
	efflon.initialize();
	efflat.initialize();
}

void CMatrices::createMatTotCatch(const PMap& map, int nbi, int nbj, int nb_species)
{
	total_obs_catch.allocate(0, nb_species - 1);
	total_pred_catch.allocate(0, nb_species - 1);
	for (int sp = 0; sp < nb_species; sp++){
		total_obs_catch(sp).allocate(map.imin, map.imax, map.jinf, map.jsup);
		total_pred_catch(sp).allocate(map.imin, map.imax, map.jinf, map.jsup);
		total_obs_catch(sp).initialize();
		total_pred_catch(sp).initialize();
	}
}

void CMatrices::createMatCatch(const PMap& map, int nbi, int nbj, int nb_species, const IVECTOR& nb_fleet, const ivector a0_adult,
			const IVECTOR& nb_cohorts, const IVECTOR& nb_region_sp)
{
	Ctot_proportion_fishery.allocate(0,nb_species-1);
	Ctot_proportion_fishery.initialize();
	for (int sp = 0; sp < nb_species; sp++){
		const int fleetmax = nb_fleet[sp];
		Ctot_proportion_fishery(sp).allocate(0,fleetmax-1);
		for (int f=0; f<fleetmax; f++){
			Ctot_proportion_fishery(sp,f).allocate(map.imin,map.imax,map.jinf,map.jsup);
			Ctot_proportion_fishery(sp,f).initialize();
		}
	}

	// captures observees par flotille par cellule
 	catch_obs.allocate(0, nb_species - 1);
	for (int sp = 0; sp < nb_species; sp++){

		const int fleetmax = nb_fleet[sp];
 		catch_obs[sp].allocate(0, fleetmax - 1);

		catch_obs[sp].initialize();
		for (int fleet = 0; fleet < fleetmax; fleet++){

 			catch_obs[sp][fleet].allocate(map.imin, map.imax, map.jinf, map.jsup);
			catch_obs[sp][fleet].initialize();
		}
	}
	// captures estimees par flotille par cellule
 	catch_est.allocate(0, nb_species - 1);
	catch_est.initialize();

	for (int sp = 0; sp < nb_species; sp++)
	{
		const int fleetmax = nb_fleet[sp];
 		catch_est[sp].allocate(0, fleetmax - 1);
		catch_est[sp].initialize();

		for (int fleet = 0; fleet < fleetmax; fleet++)
		{
 			catch_est[sp][fleet].allocate(map.imin, map.imax, map.jinf, map.jsup);
			catch_est[sp][fleet].initialize();
		}
	}
//cout << __FILE__ << ':' << __LINE__ << endl;

	// C_N_sp_age_fishery : somme des captures (en nombre) estimees par espece, par age
	// et par flotille a chaque pas de temps 
	// C_N_sp_age_fishery_qtr: somme des captures (en nombre) estimees par age et par flotille 
	// groupees par trimestre 
 	//C_N_sp_age_fishery=Utilities::create4d(C_N_sp_age_fishery,nb_species,nb_age_class,nb_fleet,nb_region);
	C_N_sp_age_fishery.allocate(0, nb_species - 1);
	LF_qtr_obs.allocate(0,nb_species-1);
	C_N_sp_age_fishery.initialize();
	LF_qtr_obs.initialize();	

	C_tot_no_effort_sp_age.allocate(0, nb_species - 1);
	C_tot_no_effort_sp_age.initialize();

	for (int sp = 0; sp < nb_species; sp++)
	{
		const int a0 = a0_adult(sp);
		const int agemax = nb_cohorts(sp);

		C_N_sp_age_fishery[sp].allocate(a0, agemax - 1);
		LF_qtr_obs[sp].allocate(a0, agemax - 1);
		C_N_sp_age_fishery[sp].initialize();
		LF_qtr_obs[sp].initialize();

		C_tot_no_effort_sp_age[sp].allocate(a0, agemax - 1);
		C_tot_no_effort_sp_age[sp].initialize();
		
		for (int age = a0; age < agemax; age++){

			const int fleetmax = nb_fleet(sp);
			C_N_sp_age_fishery[sp][age].allocate(0, fleetmax - 1);				
			LF_qtr_obs[sp][age].allocate(0, fleetmax - 1);				
			C_N_sp_age_fishery[sp][age].initialize();
			LF_qtr_obs[sp][age].initialize();
		
			for (int fleet = 0; fleet < fleetmax; fleet++){

				const int regionmax = nb_region_sp(sp);
				C_N_sp_age_fishery[sp][age][fleet].allocate(0, regionmax - 1);
				LF_qtr_obs[sp][age][fleet].allocate(0, regionmax - 1);
				C_N_sp_age_fishery[sp][age][fleet].initialize();
				LF_qtr_obs[sp][age][fleet].initialize();

			}
			C_tot_no_effort_sp_age[sp][age].allocate(map.imin, map.imax, map.jinf, map.jsup);
			C_tot_no_effort_sp_age[sp][age].initialize();

		}
	}
//cout << __FILE__ << ':' << __LINE__ << endl;
 	//C_N_sp_age_fishery_qtr=Utilities::create4d(C_N_sp_age_fishery_qtr,nb_species,nb_age_class,nb_region,nb_fleet);
	C_N_sp_age_fishery_qtr.allocate(0, nb_species - 1);
	C_N_sp_age_fishery_qtr.initialize();
	
	for (int sp = 0; sp < nb_species; sp++)
	{
		const int a0 = a0_adult(sp);
		const int agemax = nb_cohorts(sp);
		///const int agemax = nb_age_class(sp);
		C_N_sp_age_fishery_qtr(sp).allocate(a0, agemax - 1);
		C_N_sp_age_fishery_qtr(sp).initialize();

		for (int age = a0; age < agemax; age++)
		{
			const int regionmax = nb_region_sp(sp);
			C_N_sp_age_fishery_qtr(sp, age).allocate(0, regionmax - 1);			
			C_N_sp_age_fishery_qtr(sp, age).initialize();

			for (int region = 0; region < regionmax; region++)
			{
				const int fleetmax = nb_fleet(sp);
				C_N_sp_age_fishery_qtr(sp, age, region).allocate(0, fleetmax - 1);
				C_N_sp_age_fishery_qtr(sp, age, region).initialize();
			}
		}
	}

	// Sum_C_N_sp_age_fishery_area: sommes totales des captures (en nombre) estimees par age, espece,
	// flotille et region
	//Sum_C_N_sp_age_fishery_area=Utilities::create5d(Sum_C_N_sp_age_fishery_area,nb_species,nb_age_class,5,nb_region,nb_fleet);
//cout << __FILE__ << ':' << __LINE__ << endl;
	Sum_C_N_sp_age_fishery_area.allocate(0, nb_species - 1);
	Sum_C_N_sp_age_fishery_area.initialize();
	
	for (int sp = 0; sp < nb_species; sp++)
	{
		const int agemax = nb_cohorts(sp);
		Sum_C_N_sp_age_fishery_area(sp).allocate(0, agemax - 1);
		Sum_C_N_sp_age_fishery_area(sp).initialize();

		for (int age = 0; age < agemax; age++)
		{
			Sum_C_N_sp_age_fishery_area(sp, age).allocate(0, 5 - 1);
			Sum_C_N_sp_age_fishery_area(sp, age).initialize();

			for (int n = 0; n < 5; n++)
			{
				const int regionmax = nb_region_sp(sp);
				Sum_C_N_sp_age_fishery_area(sp, age, n).allocate(0, regionmax);
				Sum_C_N_sp_age_fishery_area(sp, age, n).initialize();

				for (int region = 0; region <= regionmax; region++)
				{
					const int fleetmax = nb_fleet(sp);
					Sum_C_N_sp_age_fishery_area(sp, age, n, region).allocate(0, fleetmax);
					Sum_C_N_sp_age_fishery_area(sp, age, n, region).initialize();
				}
			}
		}
	}
	// dim 5 = 4 quarters + the sum of quarters
//cout << __FILE__ << ':' << __LINE__ << endl;
}

/*void CMatrices::createMatLarvae(const PMap& map, int t0, int nbt)
{
	larvae_obs.allocate(t0, nbt);
	larvae_pred.allocate(t0, nbt);
	for (int t=t0; t<=nbt; t++){
		larvae_obs(t).allocate(map.imin, map.imax, map.jinf, map.jsup);	
		larvae_pred(t).allocate(map.imin, map.imax, map.jinf, map.jsup);	
	}
	larvae_obs.initialize();
	larvae_pred.initialize();
}*/

void CMatrices::createMatLarvae(const PMap& map)
{
	larvae_obs.allocate(map.imin, map.imax, map.jinf, map.jsup);	
	larvae_pred.allocate(map.imin, map.imax, map.jinf, map.jsup);	
	larvae_obs.initialize();
	larvae_pred.initialize();
}



void CMatrices::MeanVarMovement(const PMap& map, const dmatrix& Adv_x, const dmatrix& Adv_y, const dmatrix& Diff, 
		const double mss, const double sigma_species, const double length_age, const double length_age_max,
		const int dT, const int sp, const int age)
{
	const double Vvar_slope = 0.0;//0.15;
	const double Vmax_diff  = 1.9;

	double sum_density = sum(density_after(sp,age));
        if (sum_density == 0) {
        	const double lr = length_age/length_age_max;
                const double Vvar  = 1.0-Vvar_slope*lr;
                const double Vdiff = Vmax_diff- lr;
                const double V = mss*Vvar;

                mean_speed(sp,age) = sqrt(2*V*V);
                mean_diffusion(sp,age) = sigma_species*pow(Vdiff*length_age*0.01*3600*dT*24.0/1852,2)/(4*dT);
        } else {
                mean_speed(sp,age) = comp_waverage2(map, Adv_x, value(Adv_y), sp, age);
                //conversion:
                double cc = (1852.0/(length_age*0.01))/ ( dT * 3600 * 24.0);
                mean_speed(sp,age) = mean_speed(sp,age) * cc;
                mean_diffusion(sp,age) = comp_waverage(map, Diff, sp, age)/dT;
	}
}

void CMatrices::MeanVarMortality(const PMap& map, const dmatrix& M, const double Mp_max, const double Ms_max, const double Mp_exp, const double Ms_slope, const double mean_age_in_month, const int sp, const int age)
{
	double sum_density = sum(density_after(sp,age));
//	double sum_density = 0;
	if (sum_density == 0) {
        	mean_mortality(sp,age) = (Ms_max+Mp_max) * exp(- Mp_exp * mean_age_in_month) +
                                          Ms_max*pow(mean_age_in_month,Ms_slope);
        } else 
		mean_mortality(sp,age) = comp_waverage(map, M, sp, age);
		//mean_mortality(sp,age) = 0;
}

void CMatrices::MeanVarTemperature(const PMap& map, const int sp, const int sp_nb_cohort_lv, const int a0_adult, const int t_count)
{//larvae and juveniles only, see OnRunCoupled for other cohorts

	for (int a=0; a<sp_nb_cohort_lv; a++)
        	mean_temperature(sp,a) = comp_waverage(map, sst(t_count), sp, a);
	for (int a=sp_nb_cohort_lv; a<a0_adult; a++)
        	mean_temperature(sp,a) = comp_waverage(map, sst(t_count), sp, a);
        	//mean_temperature(sp,a) = comp_waverage(map, tempn(t_count,0), sp, a);
}


double CMatrices::comp_waverage(const PMap& map, const dmatrix& var, const int sp, const int age)
{
        double res = 0.0;
        double sum_density = sum(density_after(sp,age));
//      if (sum_density == 0) return(mean(var));
        for (int i = map.imin; i <= map.imax; i++){
                const int jmin = map.jinf[i];
                const int jmax = map.jsup[i];
                for (int j = jmin; j <= jmax; j++){
                        if (map.carte(i,j)){
//if (age<3)
                                res += var(i,j)*density_after(sp,age,i,j);
//else
//                                res += var(i,j)*density_before(sp,1,age,i,j);
                        }
                }
        }
        res /= sum_density;
        return(res);
}

double CMatrices::comp_waverage2(const PMap& map, const dmatrix& var1, const dmatrix& var2, const int sp, const int age)
{
        double vsum = 0.0;
        double res = 0.0;
       	double sum_density = sum(density_after(sp,age));

        for (int i = map.imin; i <= map.imax; i++){
                const int jmin = map.jinf[i];
                const int jmax = map.jsup[i];
                for (int j = jmin; j <= jmax; j++){
                        if (map.carte(i,j)){
                                vsum = sqrt(var1(i,j)*var1(i,j) + var2(i,j)*var2(i,j));
                                res += vsum*density_after(sp,age,i,j);
                        }
                }
        }
        res /= sum_density;
        return(res);
}

dvector CMatrices::create1d(int n)
{
	dvector vect;
	vect.allocate(0, n - 1);
	vect.initialize();	
	return vect;
} 

dmatrix CMatrices::create2d(int n1, int n2)
{
	dmatrix dmat;
	dmat.allocate(0,n1-1,0,n2-1);
	dmat.initialize();
	return dmat;
}

imatrix CMatrices::create2d_int(int n1, int n2)
{
	imatrix imat;
	imat.allocate(0,n1-1,0,n2-1);
	imat.initialize();
	return imat;
}

dmatrix CMatrices::create2d(int imin, int imax, ivector jinf, ivector jsup)
{
	dmatrix dmat;
	dmat.allocate(imin,imax,jinf,jsup);
	dmat.initialize();
	return dmat;
}
 
void CMatrices::create3d(d3_array& d3arr, int n1, int n2, int n3)
{
	d3arr.allocate(0,n1-1);
	for (int n=0; n<n1; n++){
		d3arr[n].allocate(0,n2-1,0,n3-1);
		d3arr[n].initialize();
	}
}

d3_array CMatrices::create3d(int n1, const int imin, const int imax, ivector& jinf, ivector& jsup)
{
	d3_array d3arr;
	d3arr.allocate(0,n1-1);
	for (int n=0; n<n1; n++){
		d3arr[n].allocate(imin,imax,jinf,jsup);
		d3arr[n].initialize();
	}
	return d3arr;
}
/*
d4_array CMatrices::create4d(int n1, int n2, int n3, int n4)
{
	d4_array d4;
	d3arr.allocate(0,n1-1);
	for (int n=0; n<n1; n++){
		d3arr[n].allocate(imin,imax,jinf,jsup);
		d3arr[n].initialize();
	}
	return d3arr;
}

*/

