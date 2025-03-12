#include "VarSimtunaFunc.h"

///Forward main function called in simulation mode only for: 
///mortality rates at age. See mortality_sp.cpp


void VarSimtunaFunc::Mortality_Sp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& M, dvar_matrix& H, const int sp, const double mean_age_in_dtau, const int age, const int t_count)
{
	M.initialize();

	double Rage = mat.mortality_range_age[sp][age];
	double Hval = param.Hval;

	dvariable Mp_max   = param.dvarsMp_mean_max[sp];
	dvariable Mp_exp   = param.dvarsMp_mean_exp[sp];
	dvariable Ms_slope = param.dvarsMs_mean_slope[sp];
	dvariable Ms_max   = param.dvarsMs_mean_max[sp];
	dvariable range;
	if (age==0)
		range = param.dvarsM_larvae_range[sp];
	else
		range = param.dvarsM_mean_range[sp];
	
	M_sp_comp(map, M, value(H), value(Mp_max), value(Ms_max), value(Mp_exp), value(Ms_slope), value(range), Rage, Hval, mean_age_in_dtau);

    if (!param.gcalc()){
                // pH option works only in simulation mode
		if ((age==0) && param.use_ph1)
			M_PH_juv_comp(param,map,mat,M,mat.ph1[t_count],mean_age_in_dtau);	
	
		mat.MeanVarMortality(map, value(M),value(Mp_max), value(Ms_max), value(Mp_exp), value(Ms_slope), mean_age_in_dtau, sp, age);
	}
}

void VarSimtunaFunc::Scaling_factor_sstdep_larvae_mortality(VarParamCoupled& param, const dmatrix& sst, const PMap& map, dvar_matrix& S, const int sp)
{
	int deltaT = param.deltaT;
	dvariable inv_M_max = param.dvarsInv_M_max[sp];
	dvariable inv_M_rate = param.dvarsInv_M_rate[sp];
	dvariable age_larvae_before_sst_mortality = param.dvarsAge_larvae_before_sst_mortality[sp];

	Scaling_factor_sstdep_larvae_mortality_comp(map, S, sst, inv_M_max, inv_M_rate, age_larvae_before_sst_mortality, deltaT);
}
