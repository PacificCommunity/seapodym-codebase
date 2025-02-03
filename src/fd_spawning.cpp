#include "SeapodymCoupled.h"
#include <fvar.hpp>

///Forward main function called in simulation mode only for: 
///computing the new larval biomass in two regimes of the simulation:
///during climatological spinup and after the population has been built.
///See spawning.cpp


void SeapodymCoupled::Spawning(dvar_matrix& J, dvar_matrix& Hs, dvar_matrix& Nmature, const int jday, const int sp, const int t_count)
{

	J.initialize();

	dmatrix J_c = value(J);
	dmatrix Hs_c = value(Hs);
	dmatrix N_mat = value(Nmature);
	dvariable nb_recruitment = param->dvarsNb_recruitment[sp];

	dvar_matrix Nbr(map.imin,map.imax,map.jinf,map.jsup);
	Nbr = nb_recruitment;


	dvariable a_adults_spawning = param->dvarsA_adults_spawning[sp];
	dvar_matrix A_sp(map.imin,map.imax,map.jinf,map.jsup);
	A_sp = a_adults_spawning;

	if (param->elarvae_model[sp] | param->spawning_adult_func_only[sp])
		spawning_adult_func_comp(J_c,N_mat,value(nb_recruitment),value(a_adults_spawning));
	else 
		spawning_in_hs_comp(J_c,Hs_c,N_mat,value(nb_recruitment),value(a_adults_spawning));

	J = nograd_assign(J_c);
}

