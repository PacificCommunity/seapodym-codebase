#include "VarSimtunaFunc.h"

///Forward main function called in simulation mode only for: 
///spawning habitat functions. See spawning_habitat.cpp

void VarSimtunaFunc::Spawning_Habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, const double sigma_sp_var, int sp, const int t_count, const int jday)
{
	Hs.initialize();

	dvariable a,b;
	if (!param.uncouple_sst_larvae[sp]){
		a = param.dvarsA_sst_spawning[sp];
		b = param.dvarsB_sst_spawning[sp];
	}
	else {
		a = param.dvarsA_sst_larvae[sp];
		b = param.dvarsB_sst_larvae[sp];
	}

	dvariable c = param.dvarsAlpha_hsp_prey[sp];
	dvariable d = param.dvarsAlpha_hsp_predator[sp];
	dvariable e = param.dvarsBeta_hsp_predator[sp];

	Hs_comp(param, mat, map, Hs, value(a), value(b), value(c), value(d), value(e), sigma_sp_var, jday, t_count);	
}

