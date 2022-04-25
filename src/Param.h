// Param.h: interface for the CParam class.
//
//////////////////////////////////////////////////////////////////////
// Param.h
//////////////////////////////////////////////////////////////////////
#ifndef __Param_h__
#define __Param_h__


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using  namespace std;

#include "mytypes.h"
//#include "ReadWrite.h"

/*!
\brief Seapodym parameter class.
\details All static SEAPODYM parameters are defined and described here. However for DVAR parameters see class VarParamCoupled
*/
class CParam
{
protected:

	int	nbt_total;	// nb de pas de temps dans la serie
	int	nbi;		// nb cell dans la direction est-west (+ 2) zone 1
	int	nbj;		// nb cell dans la direction nord-sud (+ 2) zone 1
	//int	maxn;		// max de nbi et nbj
	int 	nb_species;	// nb of species
	int	nb_forage;	// nb of forage components
	int	nb_fishery;	// nb absolu de pecheries    	
	
public:

	//CReadWrite rw;

	bool	flag_coupling;
	bool	build_forage;
	bool    flag_twin;
	bool	connectivity_comp;
	int     tuna_spinup;
	int     wbin_flag;
	int     mpa_simulation;	
	int	nb_mpa;
//	int     mpa_scenario;
	int     type_oxy;
	int	use_sst;
	int	use_vld;
	int	use_ph1;	

	//Optimization control
	int maxfn;
	double crit;
	
	ivector vert_movement;
	ivector food_requirement_in_mortality;
	ivector uncouple_sst_larvae;	
	ivector gaussian_thermal_function;	
	ivector cannibalism;	
	string	idformat;	//Patrick 21Oct04
	int	idfunc ;	//Patrick 21Oct04
	imatrix like_types;
	bool	cpue;
	dmatrix like_param;
	dmatrix prob_zero;
	ivector tag_like;
	ivector stock_like;
	dvector mean_stock_obs;
	dvector stock_lonmin;
	dvector stock_lonmax;
	dvector stock_latmin;
	dvector stock_latmax;
	ivector frq_like;
	dvector eff_units_converter;
	dvector cpue_mult;
	double  total_like;
	int fdata_rm;

	int flex_regstruc;
	int use_mask_catch;

	string* parfile_names;
	
	int nb_varproj;
	ivector varproj_nsteps;
	vector<string> varproj;
	
	dvector statpars;
	int 	_nstatpars;
	adstring_array statpar_names;
	

	double	longitudeMin;	// longitude (0-359) minimale de la grille
	double	longitudeMax;	// longitude (0-359) maximale de la grille
	double	latitudeMin;	// latitude minimale (S negatif)  de la grille
	double	latitudeMax;	// latitude maximale  de la grille
	double	deltaX;		// dx (constant)
	double	deltaY;		// dy (constant)
	int	deltaT;		// dt (constant)
	int	nlevel;
	double	startdate, enddate;
	int	ndatini, ndatfin;
	int     date_mode;
	ivector rundates;
	int	nbytetoskip;
	double	save_first_yr;	// first year from which predictions are recorded
	double	save_last_yr;
	int first_recruitment_date;
	int	nb_yr_forecast;
	int 	nbsteptoskip;
	int	nlong, nlat;
	int	iterationNumber;// nbre d'iteration dans calcul derivees
	int	nb_layer;	// nb of vertical layers for the forage components
	DVECTOR source_frg;	// percent of source prod transferred to each forage [nb_forage]
	IVECTOR day_layer;	// layer where the forage is during the day [nb_forage]
	IVECTOR night_layer;	// layer where the forage is during the night [nb_forage]

	double	lambda;		// Forage mortality
	double	E;		// Energy transfer coefficient from Primary prod to Forage pop
	double	c_pp;		// conversion factor from PP (N or C) to wet weight of forage
	double  pp_transform;   // conversion from PP units to wet weight of zooplankton
//	int 	tstep_forage;	// time step for the forage (days)
	double 	sigma_fcte;	// diffusion coefficient for forage (nm2.d-1)
	int	inv_lambda_max;	// Forage mortality max through time transfer
	double	inv_lambda_curv;// curvature coeff of fonction Forage mortality against SST
	int	Tr_max;		// Time max(in days) before recruitment in the forage population
	double	Tr_exp;

	string	str_file_mask;
	string	str_file_topo;
	string	str_file_maskEEZ;
	string	str_file_maskMPA;
	string  str_file_param;
	string  str_dir;
	string  str_dir_forage;
	string  str_dir_init;
	string  str_dir_fisheries;
	string  str_dir_tags;
	string  strfile_pp;
	string  strfile_sst;
	string  strfile_vld;
	string  strfile_ph1;	
	vector<string> frg_name;	// array of string [nb_forage]
	vector<string> sp_name;		// array of species name [nb_species]
	vector<string> strfile_F;	// array of string [nb_forage]
	vector<string> strfile_Fmc;	// array of string [nb_forage]
	vector<string> strfile_S;	// array of string [nb_forage]
	vector<string> strfile_Smc;	// array of string [nb_forage]
	string strfile_ppmc;
	string strfile_sstmc;
	string strfile_vldmc;
	vector<string> strfile_u,strfile_v,strfile_t,strfile_oxy;	// array of string [nb_layer]
	vector<string> strfile_umc,strfile_vmc,strfile_tmc, strfile_oxymc;		// array of string [nb_layer]
	string strdir_output;
	int write_all_cohorts_dym;
	int write_all_fisheries_dym;
	
///	IVECTOR sp_nb_age_class_ad;	// number of age classes for each species [sp]
///	IVECTOR sp_unit_age_class_ad;	// time step used for the population of the species [sp] (0= pas de calcul de pop; 1=month;2=quarter )
///	IMATRIX sp_unit_age_class;	// time step (in days) used for the population of the species [sp] and cohort [a]

///	IVECTOR sp_nb_age_class_jv;	// number of age classes for each species [sp]
///	IVECTOR sp_unit_age_class_jv;	// time step used for the population of the species [sp] (0= pas de calcul de pop; 1=month;2=quarter )
///	int	max_age_class;		// max number of age classes over all species

	vector<string> life_stage;
	ivector sp_nb_cohort_life_stage;
	ivector sp_nb_cohorts; 
	ivector sp_nb_cohort_lv, sp_nb_cohort_jv, sp_nb_cohort_ad;
	ivector sp_a0_adult;
	imatrix sp_unit_cohort;
	
///	DMATRIX juv_length;		// length by age for each species (cm) for the first three months of live
///	DMATRIX juv_weight;		// weight by age for each species (kg) for the first three months of live
	DMATRIX length;			// length by age for each species (cm) [sp][sp_nb_age_class[sp]]
	DMATRIX length_bins;		// length by age for each species (cm) [sp][sp_nb_age_class[sp]]
	DMATRIX weight;			// weight by age for each species (kg) [sp][sp_nb_age_class[sp]]

	DVECTOR M_inc_ph_a;
	DVECTOR M_inc_ph_b;	
	DVECTOR Mp_mean_max;		// natural mortality: max coeff of the "predation" decreasing exponential function
	DVECTOR Mp_mean_exp;		// natural mortality: expo coeff of the "predation" decreasing exponential function
	DVECTOR Ms_mean_slope;		// natural mortality: slope coeff of the "senescence" increasing sigmoid function
	DVECTOR Ms_mean_max;		// natural mortality: max coeff of the "senescence" increasing sigmoid function
	DVECTOR M_mean_range;		// range of the variability of natural mortality around M in relation with the habitat
	dvector residual_competition;	// constant (temp) parameter accounting for competition with species, which are not in the model

	int habitat_run_type;		// these two parameters are needed to pass the info to OnRunHabitat on which 
	int nb_habitat_run_age;		// nb of habitat indices to be run: 1 for spawning (default value) and several for feeding
	ivector habitat_run_age;	// ages to compute habitat in order to estimate all parameters
	int migrations_by_maturity_flag;

	IVECTOR age_mature;		// age of first maturity in time step unit used for the species [sp]
	DMATRIX maturity_age;		// 09.09: maturity at age for the species [sp]
	IVECTOR age_autonomous;		// age for larvae-juvenile phase with passive transport [sp]
	IVECTOR age_recruit;		// age at recruitment [sp]
	imatrix age_compute_habitat;	// ages to compute habitat index, by default it is computed for all adult cohorts
	DVECTOR nb_recruitment;		// nb of fish recruited by cell [sp]
	DVECTOR a_adults_spawning;	// coefficient controlling dependence of number of spawns from number of mature fish of species [sp]
	ivector seasonal_migrations;	// flag for seasonal migrations [sp]
	dvector spawning_season_peak;	// peak of the seasonal cycle in julian day (in the North Hemisphere)
	dvector spawning_season_start;  // day/night length ratio defining the beginning of spawning migrations
	DVECTOR a_sst_spawning;		// coefficient of curvature for spawning temperature function [sp]
 	DVECTOR b_sst_spawning;		// SST mean for spawning temperature function [sp]
	DVECTOR a_sst_larvae;		// a coefficent of the larvae habitat temperature in case uncouple_sst_larvae=1, [sp]
	DVECTOR b_sst_larvae;		// b coefficent of the larvae habitat temperature in case uncouple_sst_larvae=1, [sp]
	DVECTOR alpha_hsp_prey;		// previously alpha_spawning became two parameters - one for the prey function
	DVECTOR alpha_hsp_predator;	// and another two for predators function (which can be either normal or lognormal function
	DVECTOR beta_hsp_predator;	// with alpha (mean) and beta (sigma) parameters). Names are chosen generic in case the 
					// function will be changed

	DVECTOR a_sst_habitat;		// a coefficient of the habitat temperature function [sp]
	DVECTOR b_sst_habitat;		// b coefficent of the habitat temperature function [sp]
	DVECTOR T_age_size_slope;	// T_age as a function of fish size [sp]
	dmatrix thermal_func_delta;	// alternative (non-Gaussian) thermal function parameters  
	DVECTOR a_oxy_habitat;		// a coefficient of the habitat 2ml O2 depth function [sp]
	DVECTOR b_oxy_habitat;		// b coefficent of the habitat 2ml O2 depth function [sp]
	dmatrix eF_habitat;		// Forage groups scaling parameters allowing to account 
					// for the uncertainty in energy distribution parameters 
					// in MTL model and for the species preferences of groups 
					// of MTL. Thus, 1) they are species dependent;
					// 2) they need to be control variables in optimization,
					// 3) there are six of them named after forage groups
					

	DVECTOR hp_cannibalism;		// half_pred coefficent of the juvenile habitat [sp]
	
	DVECTOR forage_ration;		// mean daily ration (weight forage/weight fish) for each species [sp]
	DVECTOR sigma_species;		// max diffusion coefficient for each species in nm2/mo [sp]
	DVECTOR MSS_species;		// max ustained speed for each species in FL/mo [sp]
	DVECTOR MSS_size_slope;		// scaling exponent of the power low to compute sustainable speed MSS*L^slope
	DVECTOR c_diff_fish;		// coefficient for the diffusion-Habitat function [sp]

	dmatrix sigma_ha;		//Gaussian std in adult habitat by sp and age
	dmatrix temp_age;		//optimal temperature by sp and age

	// Fisheries
	string	*list_fishery_name;	// liste complete des pecheries (nombre absolu) [nb_fishery]
	dvector fishery_reso;		// original resolution of fishing data 
	float catch_reso;		// resolution of catch data in degrees (can be degraded from original)
	
	ivector fishery_catch_units;	// 1- catch in metric tons, 0 - catch in number of individuals
	IVECTOR nb_fishery_by_sp;	// nb de pecheries concernees par espece [sp]
	IMATRIX mask_fishery_sp;
	IMATRIX mask_fishery_sp_no_effort;
	IMATRIX mask_fishery_sp_like;
	//IMATRIX mask_mpa_fishery;
	int nb_fishery_type;		// nb de types differents de pecheries.
	ivector fisheries_no_effort_exist;//vector of flags [sp]


        //MPA scenarios
        int actual_eff;                 // if 1: actual effort distribution will be used, if 0 - homogeneous
        ivector mpa_scenario;           // list of MPA scenarios
        ivector mpa_ID;                 // MPA IDs to be read from MPA mask
        ivector mpa_S1_X;               // percentage of E increase in case of mpa_scenario=1 
        ivector mpa_fishery;            // Fsheries to be considered in MPA scenarios	
					
	IVECTOR type_each_fishery;	// code (entier): 1 code par pecherie. [nb_fishery]
	IVECTOR list_fishery_type; 	// liste des codes des pecheries (elimine les repetitions)[nb_fishery_type]
	
	IVECTOR nb_fishery_type_sp; 	// nb de types de pecheries concernees par espece [sp]
	IMATRIX list_fishery_type_sp; 	// liste des codes des pecheries par espece [sp][nb_fishery_type_sp[sp]]
	
	DMATRIX q_sp_fishery;		// catchability by species and fishery [sp][sp_unit_age_class[sp]]
	dvector q_dyn_fishery;
	ivector s_func_type;
	DMATRIX s_slope_sp_fishery; 	// selectivity slope coefficient for species and fishery [sp][nb_fishery]
	DMATRIX s_length_sp_fishery; 	// selectivity threshold coefficient for species and fishery [sp][nb_fishery]
	DMATRIX s_asympt_sp_fishery; 	// selectivity asymmetric Gaussian second sigma [sp][nb_fishery]
	D3_ARRAY selectivity_sp_fishery_age;	 // selectivity by species, fishery and age
	vector< vector<string> > name_sp_by_file;// liste des noms (codes) d' especes pour chaque fichier de donnees de pecheries
	vector<string>  file_catch_data;
	vector<string>  file_frq_data;
	vector<string>  file_tag_data;
	int nb_catch_files, nb_frq_files, nb_tag_files;
	int tag_gauss_kernel_on;
	int dx_tags, dy_tags; 		// setup of the grid to aggregate tagging data
	float lonmin_tags,lonmax_tags, latmin_tags, latmax_tags;
	bool tags_only;			// flag to deactivate all likelihood terms except tag_like

	string  m_file_in_str;		// fichier lecture
	string  m_file_out_str; 	// fichier ecriture
	double	m_f;			// multiplicateur d'effort de peche
	int	nb_region;		// nb total de regions ou aggreger les resultats
	IVECTOR nb_region_sp_B; 	// nb de regions par espece ou aggreger les biomasses
	IVECTOR nb_region_fishery; 	// nb de regions par fishery 
	IMATRIX area_sp_B; 		// liste des regions par espece ou aggreger les biomasses [sp][area_id]
	
	struct	region			// definition de la structure region
	{
		int area_id;
		double lgmin;
		double lgmax;
		double ltmin;
		double ltmax;
	};
	
	region **area;			// structure contenant les limites des regions [nb_region](area_id, lgmin, lgmax,ltmin,ltmax)
	int	nb_EEZ;			// nb de ZEE
	IVECTOR EEZ_ID;			// code des ZEE [nb_EEZ]
	string	*EEZ_name;		// name of EEZ
public:
	CParam();
	virtual ~CParam();
	
	void init_param() ;
	void init_param_dym();
	void delete_param(bool flag) ; 
	void read_param(bool &file_found){};
	void write_param(char runtype){};
	void rbin_input2d(string file_in, const imatrix& carte, DMATRIX& mat2d, int nbi, int nbj, int nbytetoskip);
	void rbin_input2d(string file_in, DMATRIX& mat2d, int nbi, int nbj, int nbytetoskip);
	void rbin_mat2d(string file_in, const imatrix& carte, DMATRIX& mat2d, int nlat, int nlong, int nbytetoskip);
	void rbin_mat2d(string file_in, DMATRIX& mat2d, int nlat, int nlong, int nbytetoskip);
	double correction_lat(double lat);
	double lastlat(int j);
	double cell_surface_area(int j);
	double jtolat(int j);
	double itolon(int i);
	int lattoj(double lat);
	int lontoi(double lon);
	double func_limit_one(const double m);
	double dffunc_limit_one(const double x, const double dfy);
	//double dffunc_limit_one(const double m);
	void afcoef(const double lon, const double lat, dmatrix& a, int& ki, int& kj, const int reso);

	double selectivity_comp(const int sp, const int age, const int f, const int k);
	//double selectivity_comp(const int sp, const int age, const int f, const int k, const int step_count);
	void dfselectivity(double& dfslope, double& dflength, double& dfasympt, const int sp, const int age, const int f, const int k);
	void define_regions();
	float fdate(float year, float month) {return (year+(month-0.5)/12.0);} //temporal: works with monthly steps only
	int get_month(double fdate) {return ((fdate-(int)fdate)*12.0+0.505);} //temporal: works with monthly steps only
	int get_year(double fdate) {return ((int)fdate);} 
	inline int get_nbi() const {return nbi;}
	inline int get_nbj() const {return nbj;}
	inline void set_nbt(int nbt) {nbt_total = nbt;}
	inline int get_nbt() {return nbt_total;}
//	inline int get_maxn() const {return maxn;}
//	inline int get_max_age_class() const {return max_age_class;}
	inline int get_nbspecies() {return nb_species;}
	inline int get_nbfishery() const {return nb_fishery;}
	inline int get_nbforage() const {return nb_forage;}

	void time_reading_init(){elapsed_time_reading = 0;}
	double elapsed_time_reading;
};
#endif 

