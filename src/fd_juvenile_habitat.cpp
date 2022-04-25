#include "VarSimtunaFunc.h"

///Forward main function called in simulation mode only for: 
///juvenile habitat functions. See juvenile_habitat.cpp

void VarSimtunaFunc::Juvenile_Habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hj, int sp, const int t_count)
{
	Hj.initialize();

	dvariable a,b;

	//These two parameters will impact only the mortality range of juveniles
	//Extending the tolerance interval for juveniles compared to larve
	if (!param.uncouple_sst_larvae[sp]){
		a = param.dvarsA_sst_spawning[sp]+1.0;
		b = param.dvarsB_sst_spawning[sp];
	}
	else {
		a = param.dvarsA_sst_larvae[sp]+1.0;
		//SKJ condition temporarilly to enable the use of REF-2018 (CJFAS) solution with this version
		if (param.sp_name[0].find("skj")!=0)
			b = param.dvarsB_sst_larvae[sp];
		if (param.sp_name[0].find("skj")==0)
			b = param.dvarsB_sst_spawning[sp];
	}

	Hj_comp(param, mat, map, Hj, value(a), value(b),t_count);
	
}

void VarSimtunaFunc::Juvenile_Habitat_cannibalism(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hj, dvar_matrix& total_pop, int sp, const int t_count)
{
	Hj.initialize();

	dvariable a,b;

	//These two parameters will impact only the mortality range of juveniles
	//Extending the tolerance interval for juveniles compared to larve
	if (!param.uncouple_sst_larvae[sp]){
		a = param.dvarsA_sst_spawning[sp]+1.0;
		b = param.dvarsB_sst_spawning[sp];
	}
	else {
		a = param.dvarsA_sst_larvae[sp]+1.0;
		b = param.dvarsB_sst_larvae[sp];
		
	}
	dvariable c = param.dvarsHp_cannibalism[sp];

	Hj_cannibalism_comp(param, mat, map, Hj, value(total_pop), value(a), value(b), value(c), t_count);
	
}


