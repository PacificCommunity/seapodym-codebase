#include "SeapodymCoupled.h"

int SeapodymCoupled::EditRunCoupled(const char* parfile)
{
	if (parfile == 0)
	{
		cerr << "Error[" __FILE__ << ':' << __LINE__ << "]: parameter parfile is NULL in \"int SeapodymCoupled::EditRunCoupled(const char* parfile)\"\n";
		return 0;
	}
	runtype ='C';
	
	param->str_file_param = parfile;

	param->read(param->str_file_param);
	//initialization of class variables
	nbi = param->get_nbi();    
	nbj = param->get_nbj();    
	nb_species = param->get_nbspecies();
	nb_fishery = param->get_nbfishery();
	nb_forage  = param->get_nbforage();
	nb_layer   = param->nb_layer;
	deltaX     = param->deltaX;
	deltaY     = param->deltaY;
	deltaT     = param->deltaT;
	tuna_spinup= param->tuna_spinup;
	cell_area  = 1.852*deltaX*1.852*deltaY;


	a0_adult.allocate(0,nb_species-1);
	aN_adult.allocate(0,nb_species-1);
	mean_age_cohort.allocate(0,nb_species-1);		
	for (int sp=0; sp<nb_species; sp++){
		a0_adult[sp] = param->sp_a0_adult[sp];
		aN_adult[sp] = param->sp_nb_cohorts[sp];
		mean_age_cohort[sp].allocate(0,aN_adult[sp]-1);
		mean_age_cohort[sp][0] = (double).5*param->sp_unit_cohort[sp][0];
		for (int a=1; a<aN_adult[sp]; a++)
                	mean_age_cohort[sp][a] = mean_age_cohort[sp][a-1]+
				.5*param->sp_unit_cohort[sp][a-1]+.5*param->sp_unit_cohort[sp][a];
	}
	map.lit_map(*param);
	nb_tagpops = param->nb_tag_files;

	if (param->tag_like[0] && param->use_tag_masks){
		for (int tagpop=0; tagpop<nb_tagpops; tagpop++){
			PMap tagmap;
			tagmap.lit_tagmap(*param, tagpop);
			tagmaps.push_back(tagmap);
		}
	}
	pop.InitCalPop(*param, map);
	//Notice, here param variables will be rewritten, it's better to add the test to check their identity
	//rw.rbin_headpar(param->strfile_pp, param->nlong, param->nlat, param->nlevel);
	rw.rbin_headpar(param->strfile_u[0], param->nlong, param->nlat, param->nlevel);
	if (nbi-2 != param->nlong || nbj-2!=param->nlat) {
		cerr << "INPUT files: dimensions do not match: " << 
			nbi-2 <<" != " <<param->nlong <<"; " << nbj-2 <<" != " <<param->nlat << 
			". Please correct. Will Stop now!" << endl;
		exit(1);
	}
//	rw.set_par_domain(*param);
	nlon = param->nlong;
	nlat = param->nlat;
	mat.createMatHeader(*param);

	double mval;
	rw.rbin_header(param->strfile_pp, param->idformat, param->idfunc, mval, mval, nlon, nlat, param->nlevel, 
			param->startdate, param->enddate, mat.xlon, mat.ylat, mat.zlevel, mat.mask);
	nbt_start_series = Date::dym_startdate_run(*param,mat.zlevel,nbt_total);	

	//Create only time-independent variables here
	mat.CreateMatTransport(map, nbi, nbj);
	mat.createMatNoBorder(nbi-2, nbj-2);
	int Tr_step = (int) (param->Tr_max / deltaT);
	mat.createMatSource(nb_forage, Tr_step+1, nbi, nbj);	
	mat.createMatMortality(nb_forage, nbi, nbj);

	//temporarily keep the total catch even in F0 runs
	//as it is yet used to write tag recaptures in tags_only runs
	mat.createMatTotCatch(map,nbi,nbj,nb_species);
	if (nb_fishery && !param->flag_no_fishing){
		
		mat.createMatEffort(map,nbi,nbj,nb_fishery);
		mat.CreateMatCatch(map,nbi,nbj,nb_species, param->nb_fishery_by_sp,a0_adult,aN_adult,param->nb_region_sp_B);
	}
	lflike_fishery.allocate(0,nb_fishery-1);
	clike_fishery.allocate(0,nb_fishery-1);
	param->init_param_dym();

	return 0;
}
