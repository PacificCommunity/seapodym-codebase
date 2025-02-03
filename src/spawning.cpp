#include "SeapodymCoupled.h"
#include <fvar.hpp>

///Forward functions for: 
///computing the new larval biomass in two regimes of the simulation:
///during climatological spinup and after the population has been built.


double SeapodymCoupled::tetafunc(const double teta, const double arg)
{
	return 1.0/(1.0+exp(-teta*(arg)));
}

void SeapodymCoupled::spawning_adult_func_comp(dmatrix& J, const dmatrix Nmature, double R, double b)
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
                                        J(i,j) = f_adults;
                        }
                }
        }
}


void SeapodymCoupled::spawning_in_hs_comp(dmatrix& J, dmatrix& Hs, const dmatrix Nmature, double R, double b)
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
