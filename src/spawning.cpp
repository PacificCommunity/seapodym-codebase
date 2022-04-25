#include "SeapodymCoupled.h"
#include <fvar.hpp>

///Forward functions for: 
///computing the new larval biomass in two regimes of the simulation:
///during climatological spinup and after the population has been built.

void SeapodymCoupled::spawning_spinup_comp(dmatrix& J, dmatrix& Hs, dmatrix& SST, double nb_recruitment, double Tmin)
{

	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){
				//nb_recruitment is given in thousands of individuals per cell
				//double unit_transfer = area / mat.lat_correction[j];
				//units of J = thousand of individuals per km^2
				J(i,j) =  nb_recruitment * Hs(i,j) * tetafunc(2.0,SST(i,j)-Tmin);
			}							
		}
	}
}

double SeapodymCoupled::tetafunc(const double teta, const double arg)
{
	return 1.0/(1.0+exp(-teta*(arg)));
}


void SeapodymCoupled::spawning_built_comp(dmatrix& J, dmatrix& Hs, const dmatrix Nmature, double R, double b)
{
	//////////////////////////////////////////////////////////
	//nb_recruitment, here R: thousand of larvae per km^2 being survived
	//Units of adults, i.e. Nmature(i,j): Nb. of ind. per km^2
	//Units of J(i,j): Nb. of ind. per km^2
	for (int i = map.imin; i <= map.imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i]; 
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){ 

					double Nm = Nmature(i,j);
					double f_adults = 1000.0*R*Nm/(1.0+b*Nm);
					J(i,j) = f_adults * Hs(i,j);
			}
		}		
	}
}
