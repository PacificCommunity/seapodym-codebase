// Matrices.h: interface for the CMatrices class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __Matrices_h__
#define __Matrices_h__


#include "Map.h"
#include "Param.h"

/*!
\brief Seapodym matrices class.
*/

class CMatrices
{
public:
	CMatrices() {/*DoNothing*/};
	virtual ~CMatrices() {/*DoNothing*/};

	friend class dim;

	DMATRIX xlon;	//vecteur contenant les valeurs des longitudes
	DMATRIX ylat;	//vecteur contenant les valeurs des latitudes
	DVECTOR zlevel;  //vecteur contenant les valeurs des latitudes
	IMATRIX mask;		//mask (sans les bords)
	//environmental data
	DMATRIX u;		// courant zonal moyen selon temps passe dans les differentes couches, utilise dans calpop
	DMATRIX v;
	
	DMATRIX diffusion_x;
	DMATRIX advection_x;
	DMATRIX diffusion_y;
	DMATRIX advection_y;
	DMATRIX speed; 		// magnitude of species velocity (without passive component);

	DVECTOR lastlat;	// Latitude of cell j
	DVECTOR lat_correction;	// Correction for cell area at latitude j 
//	DVECTOR maxGD_lat;	// Maximal gradient of the daylength at a given latitude
	dmatrix daylength;	// Length of day based on latitude and date
//	dmatrix grad_daylength; // gradient of daylength;
//	dvector dDL; 		// gradient of daylength;
	D3_ARRAY np1;
	D3_ARRAY sst;
	D3_ARRAY ph1;		//Inna 16/01/2017: adding PH (upper layer only) variable that will impact juvenile mortality
	D3_ARRAY vld;
	D4_ARRAY un;		// courant zonal dans la couche n
	D4_ARRAY vn;		// courant meridien dans la couche n
	D4_ARRAY tempn;		// temperature dans la couche n
	D4_ARRAY oxygen;	// oxygen dans la couche n
	D4_ARRAY forage;
	D3_ARRAY season_switch;
	D3_ARRAY sigma_season;
	//D3_ARRAY sigma_ha_season;
	//D3_ARRAY mort_by_sp;	// mort_by_sp: mortality on forage due to all species described in the model
	//DVECTOR nF_access;	// coeff d'accessibilite (absolue) a chacune des composantes forage [nb_forage]
	//DVECTOR nF_ratio;	

	d3_array fluxes_region;
	D3_ARRAY mats;
	DMATRIX Hs;
	DMATRIX Hj;
	DMATRIX Ha;
	DMATRIX mat2d_NoBorder;
	D3_ARRAY mortality;
	
	ivector nb_age_built;	// counter for the cohort 'built' during spinup
        dmatrix mean_speed;     // weigthed (by cohort distribution) average of cohort speed
        dmatrix mean_diffusion; // weigthed (by cohort distribution) average of cohort diffusion rate
        dmatrix mean_mortality; // weigthed (by cohort distribution) average of mortality-at-age
        dmatrix mean_temperature;// weigthed average of cohorts ambient temperature
	D3_ARRAY larvae;	// larvae: biomasse des larves (age0 ->1 mois)
	D3_ARRAY juvenile;	// juvenile: biomasse des juveniles (age1 ->age_ autonomous)
	D3_ARRAY young;		// young: biomasse des jeunes (age autonomous-> age mature)
	D3_ARRAY recruit;	// recruit: biomasse des recrues (age recruit)
	D3_ARRAY adult;		// adult: biomasse des adultes (age mature-> max age class)
	D3_ARRAY total_pop;	// total_pop: biomasse totale = somme de toutes les classes d'ages
	D3_ARRAY PEB;		// Population Exploitable Biomass: sum of age classes B x average selectivity function 

	d4_array habitat_input;
	d3_array density_input;	
	d4_array larvae_input;

	std::vector<double> seasonal_larvae_input_vectors[4];// Vector of non-NA observed larvae densities
	std::vector<int> seasonal_larvae_input_vectors_i[4];// Corresponding i indices
	std::vector<int> seasonal_larvae_input_vectors_j[4];// Corresponding j indices

	D3_ARRAY total_obs_catch;	
	D3_ARRAY total_pred_catch;	
	d4_array Ctot_proportion_fishery; //proportion of C(f) in total catch, defined for all (i,j)

	D4_ARRAY init_density_species;  //init distributions by species and by age
	D5_ARRAY F_access_sum_age; 	//sum of F_access over ages, used in adjoint	
	D5_ARRAY density_before; 	//fish density before solving ADRE, used in adjoint	
	D5_ARRAY adult_habitat; 	//store adult habitat index to be used in adjoint
	D4_ARRAY density_after; 	//fish density after solving ADRE, used for computation of weigthed mean variables

	DVECTOR sum_B_larvae ;	// total biomass of larvae by species [nb_species]
	DVECTOR sum_B_juv ;	// total biomass of juvenile by species [nb_species]
	DVECTOR sum_B_young ;	// total biomass of young by species [nb_species]
	DVECTOR sum_B_recruit ;	// total biomass of recruit by species [nb_species]
	DVECTOR sum_B_adult ;	// total biomass of adults by species [nb_species]
	DVECTOR sum_total_pop ;	// total biomass of pop by species [nb_species]

	D3_ARRAY effort;	// observed fishing effort by fishery and space [nb_fishery][i][j]
	D3_ARRAY efflon,efflat;	// lon-lat coordinates of observed fishing effort by fishery and space [nb_fishery][i][j]
	D4_ARRAY catch_obs;	// catch observed by [nb_species][nb_fishery_by_sp[sp]][i][j]
	D4_ARRAY catch_est;	// catch predicted by [nb_species][nb_fishery_by_sp[sp]][i][j]
	D4_ARRAY C_N_sp_age_fishery;	// somme des captures estimees [nb_species][nb_age_class[sp]][nb_fishery_by_sp[sp]]
	D4_ARRAY C_tot_no_effort_sp_age;// total catch at age for fisheries without effort data, i.e. flagged with '1' in mask_fishery_no_effort 
	D4_ARRAY LF_qtr_obs;	// observed length frequencies at age in the quarter, for given fleet and region
	D4_ARRAY C_N_sp_age_fishery_qtr; // somme par trimestre des captures estimees [nb_species][nb_age_class[sp]][nb_fishery_by_sp[sp]]
	D5_ARRAY Sum_C_N_sp_age_fishery_area; // somme totale par trimestre des captures estimees

	void createMatHeader(const CParam& param);
	//void createMatHeader(const CParam& int nlong, int nlat, int nlevel);
	void createMatOcean(const PMap& map, int t0, int nbt, int nbi, int nbj, int nb_layer, int dt);
	void createMatTransport(const PMap& map);//, int nbi, int nbj);

	void createMatFluxes(const int nb_region, const int nb_cohort);
	
	//void createMatSource(const CParam& param);
	void createMatSource(int nforage, int ntr, int nbi, int nbj);
	void createMatNoBorder(int nbi, int nbj);
	//void createMatNoBorder(const CParam& param);
	//void createMatForage(const CParam& param);
	void createMatForage(const PMap& map, int nforage, int t0, int nbt, int nbi, int nbj);
	void createMatHabitat(const PMap& map, const int nb_forage, const int nb_species,int t0, int nbt, const ivector sp_adult_age0, const ivector sp_nb_age_class, const imatrix age_compute_habitat);
	void createMatHabitat_input(const PMap& map, const int nb_ages, const int nbt_total);
	void createMatSpecies(const PMap& map, int t0, int nbt, int nbi, int nbj, int nb_species, const ivector a0_adult, const ivector sp_nb_age_class);
	void createMatEffort(const PMap& map, int nbi, int nbj, int nb_fleet);
	void createMatTotCatch(const PMap& map, int nbi, int nbj, int nb_species);
	void createMatCatch(const PMap& map,int nbi,int nbj,int nb_species,const IVECTOR& nb_fleet,const ivector a0_adult, const IVECTOR& nb_cohorts, const IVECTOR& nb_region);
	void createMatMortality(int nforage, int nbi, int nbj);
	void MeanVarMovement(const PMap& map, const dmatrix& Adv_x, const dmatrix& Adv_y, const dmatrix& Diff, 
		const double mss, const double sigma_species, const double length_age, const double length_age_max,
		const int dT, const int sp, const int age);
	void MeanVarMortality(const PMap& map, const dmatrix& M, const double Mp_max, const double Ms_max, 
		const double Mp_exp, const double Ms_slope, const double mean_age_in_month, const int sp, const int age);
	void MeanVarTemperature(const PMap& map, const int sp, const int sp_nb_cohort_lv, const int a0_adult, const int t_count);

	double comp_waverage(const PMap& map, const dmatrix& var, const int sp, const int age);
	double comp_waverage2(const PMap& map, const dmatrix& var1, const dmatrix& var2, const int sp, const int age);



private:
	dvector create1d(int n);
	dmatrix create2d(int n1, int n2);
	imatrix create2d_int(int n1, int n2);
	dmatrix create2d(int imin, int imax, ivector jinf, ivector jsup);
	//d3_array create3d(int n1, int n2, int n3);
	void create3d(d3_array& d3arr, int n1, int n2, int n3);
	d3_array create3d(int n1, const int imin, const int imax, ivector& jinf, ivector& jsup);
//	d4_array create3d(int n1, int n2, int n3, int n4);


};
/*
class dim
{
  public:
	dim() {};
	virtual ~dim() {};

  private:
    	unsigned nbt;		// nb de pas de temps dans la serie
	unsigned nbi;		// nb cell dans la direction est-west (+ 2) zone 1
	unsigned nbj;		// nb cell dans la direction nord-sud (+ 2) zone 1
	unsigned nbz;		// previously nb_layer, it is vertical dimension
	unsigned nb_species;	// nb of species
	unsigned nb_forage;	// nb of forage components
	unsigned nb_region;	//     	
	unsigned nb_fishery;	// nb absolu de pecheries    	
	ivector nb_cohorts;	// nb of cohorts by species (should be a vector!)
	ivector nbc_lv;
	ivector nbc_jv;
	ivector nbc_yn;
	ivector nbc_ad;
	ivector i0_ad;

  public:
	void set_dim(const int ni, const int nj, const int nz, const int nt,
		    const int nbf, const int nbs, const int nbc, const int nc_lv,
		    const int nc_jv, const int nc_yn, const int nc_ad, const int a0_ad,
		    const int nb_region, const int nb_fisheries){

		//spacial dimensions
		nbi = ni; nbj = nj; nbz = nz; 
		//temporal 
		nbt = nt;
		//populations
		nb_forage = nbf; nb_species = nbs; nb_cohorts = nbc;
		//cohorts
		nbc_lv = nc_lv; nbc_jv = nc_jv; nbc_yn = nc_yn; nbc_ad = nc_ad; i0_ad = a0_ad;

	}

	unsigned get_nbi(){  return nbi;}
	unsigned get_nbj(){  return nbj;}
	unsigned get_nbz(){  return nbz;}
	unsigned get_nbt(){  return nbt;}
	unsigned get_nb_species(){  return nb_species;}
	unsigned get_nb_forage(){  return nb_forage;}
	unsigned get_nb_region(){  return nb_region;}
	unsigned get_nb_fishery(){  return nb_fishery;}
	ivector get_nb_cohorts(){  return nb_cohorts;}
	ivector get_nbc_lv(){  return nbc_lv;}
	ivector get_nbc_jv(){  return nbc_jv;}
	ivector get_nbc_yn(){  return nbc_yn;}
	ivector get_nbc_ad(){  return nbc_ad;}
	ivector get_i0_ad(){  return i0_ad;}
*/
/*
	static inline dvector create1d(int n)
	{
		dvector vect;
		vect.allocate(0, n - 1);
		vect.initialize();	
		return vect;
	} 

	static inline dmatrix create2d(int n1, int n2)
	{
		dmatrix dmat;
		dmat.allocate(0,n1-1,0,n2-1);
		dmat.initialize();
		return dmat;
	}

	static inline imatrix create2d_int(int n1, int n2)
	{
		imatrix imat;
		imat.allocate(0,n1-1,0,n2-1);
		imat.initialize();
		return imat;
	}
	
	static inline dmatrix create2d(int imin, int imax, ivector jinf, ivector jsup)
	{
		dmatrix dmat;
		dmat.allocate(imin,imax,jinf,jsup);
		dmat.initialize();
		return dmat;
	}
	 
	static inline d3_array create3d(int n1, int n2, int n3)
	{
		d3_array d3arr;
		d3arr.allocate(0,n1-1);
		for (int n=0; n<n1; n++){
			d3arr[n].allocate(0,n2-1,0,n3-1);
			d3arr[n].initialize();
		}
		return d3arr;
	}

	static inline d3_array create3d(int n1, const int imin, const int imax, ivector& jinf, ivector& jsup)
	{
		d3_array d3arr;
		d3arr.allocate(0,n1-1);
		for (int n=0; n<n1; n++){
			d3arr[n].allocate(imin,imax,jinf,jsup);
			d3arr[n].initialize();
		}
		return d3arr;
	}
	*/
//};


#endif 

