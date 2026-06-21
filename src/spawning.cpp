#include "SeapodymCoupled.h"
#include <fvar.hpp>

///Forward functions for: 
///computing the new larval biomass in two regimes of the simulation:
///during climatological spinup and after the population has been built.


double SeapodymCoupled::tetafunc(const double teta, const double arg)
{
	return 1.0/(1.0+exp(-teta*(arg)));
}

void SeapodymCoupled::spawning_adult_func_comp(dmatrix& J, const dmatrix Nmature, double R, double b, const double a)
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

					double Nm = pow(Nmature(i,j),1.0+a);
                                        J(i,j) = 1000.0*R*Nm/(1.0+b*Nm);
                        }
                }
        }
}

void SeapodymCoupled::spawning_in_hs_comp(dmatrix& J, dmatrix& Hs, const dmatrix Nmature, double R, double b, const double a)
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

					double Nm = pow(Nmature(i,j),1.0+a);
                                        double f_adults = 1000.0*R*Nm/(1.0+b*Nm);
					J(i,j) = f_adults * Hs(i,j);
			}
		}		
	}
}

// --- optional Beverton-Holt half-saturation variant: J = 1000*R*Nm/(1000*b+Nm) ---
//     R=nb_recruitment is the asymptotic recruitment; b=a_adults_spawning is the
//     half-saturation adult abundance in THOUSANDS of fish (1000*b in the denom).
//     The classical BH can be reparameterized to it as R_new = R/b; b_new = 1/b
void SeapodymCoupled::spawning_adult_func_BHsat_comp(dmatrix& J, const dmatrix Nmature, double R, double b, const double a)
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

					double Nm = pow(Nmature(i,j),1.0+a);
                                        J(i,j) = 1000.0*R*Nm/(b+Nm);
                        }
                }
        }
}

void SeapodymCoupled::spawning_in_hs_BHsat_comp(dmatrix& J, dmatrix& Hs, const dmatrix Nmature, double R, double b, const double a)
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

					double Nm = pow(Nmature(i,j),1.0+a);
                                        double f_adults = 1000.0*R*Nm/(b+Nm);
					J(i,j) = f_adults * Hs(i,j);
			}
		}		
	}
}
