//Inna 21/10/10
#ifndef DIMENSIONS_H_
#define DIMENSIONS_H_

#include "XMLDocument2.h"

class Dimensions {

private: 
  	Dimensions(); 
	~Dimensions() {} 
	Dimensions(const Dimensions &);             
	Dimensions & operator=(const Dimensions &); 

	XMLDocument2 doc;
public:
	static Dimensions& getInstance();
/*	static void setDimensions(const int ni, const int nj, const int nz, const int nt,
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
*/
   	unsigned nbt;		// nb de pas de temps dans la serie
	unsigned nbi;		// nb cell dans la direction est-west (+ 2) zone 1
	unsigned nbj;		// nb cell dans la direction nord-sud (+ 2) zone 1
	unsigned nbz;		// previously nb_layer, it is vertical dimension
	unsigned nb_species;	// nb of species
	unsigned nb_forage;	// nb of forage components
	unsigned nb_fishery;	// nb absolu de pecheries    	
	unsigned nb_region;	//     	
	unsigned nb_cohorts;	// nb of cohorts by species (should be a vector!)
	unsigned nbc_lv;
	unsigned nbc_jv;
	unsigned nbc_yn;
	unsigned nbc_ad;
	unsigned i0_ad;

	unsigned get_nbi(){  return nbi;}
	unsigned get_nbj(){  return nbj;}
	unsigned get_nbz(){  return nbz;}
	unsigned get_nbt(){  return nbt;}
	unsigned get_nb_species(){  return nb_species;}
	unsigned get_nb_forage(){  return nb_forage;}
	unsigned get_nb_cohorts(){  return nb_cohorts;}
	unsigned get_nbc_lv(){  return nbc_lv;}
	unsigned get_nbc_jv(){  return nbc_jv;}
	unsigned get_nbc_yn(){  return nbc_yn;}
	unsigned get_nbc_ad(){  return nbc_ad;}
	unsigned get_i0_ad(){  return i0_ad;}
	unsigned get_nb_region(){  return nb_region;}
	unsigned get_nb_fleet(){  return nb_fishery;}
};
#endif /* DIMENSIONS_H_ */

