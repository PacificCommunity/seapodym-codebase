#include "VarSimtunaFunc.h"

///Forward main function called in simulation mode only for: 
///spawning habitat functions. See spawning_habitat.cpp

void VarSimtunaFunc::Spawning_Habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, const double sigma_sp_var, int sp, const int t_count, const int jday)
{
	Hs.initialize();

	Hs_comp(param, mat, map, Hs, sigma_sp_var, sp, jday, t_count);	
}

