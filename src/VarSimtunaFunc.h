#ifndef __VarSimtunaFunc_h__
#define __VarSimtunaFunc_h__

#include "SimtunaFunc.h"
#include "VarParamCoupled.h"
#include "VarMatrices.h"

/*!
\brief All SEAPODYM functions including DVAR parameters
*/

class VarSimtunaFunc : public CSimtunaFunc
{
public:
	VarSimtunaFunc() {/*DoesNothing*/}
	virtual ~VarSimtunaFunc() {/*DoesNothing*/}

public:

	//admodel functions
//	void spawning_habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& H, int sp);
//	DOUBLE Juvenile_habitat(VarParamCoupled& param, CMatrices& mat, int i, int j, int sp);
//	void mortality_sp(VarParamCoupled& param, const PMap& map, dvar_matrix& M, dvar_matrix H, const int sp, const double mean_age_in_month);
//	dvariable f_accessibility(VarParamCoupled& param, CMatrices& mat, dvariable teta, int i, int j, int sp, int age, int nf);
//	void feeding_habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Ha, int sp, int age);

	//void nF_ratio_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, const int sp, const int age, int i, int j);
	//void layer_access_comp(VarParamCoupled& param, CMatrices& mat, PMap& map, const int sp, const int age);

	//adjoint code
	void Spawning_Habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, const double sigma_sp_var, int sp, const int t_count, const int jday);
	void Hs_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, double a, double b, double c, double d, double e, const double sigma_sp_var, const int jday, int t_count);
	double Hs_comp_elem(CMatrices& mat, dvector F, const double pp_transform, const double a, const double b, const double c, const double d, const double e, const double sigma_sp_var, const int nb_forage, ivector day_layer, ivector night_layer, const int jday, const int t, const int i, const int j);
	//void Hs_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, double a, double b, double c, const int jday, int t_count);
	//void Hs_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, double a, double b, double c, double d, double e, const int jday, int t_count);
	//double Hs_comp_elem(CMatrices& mat, dvector F, const double pp_transform, const double a, const double b, const double c, const int nb_forage, ivector day_layer, ivector night_layer, const int jday, const int t, const int i, const int j);
	//double Hs_comp_elem(CMatrices& mat, dvector F, const double pp_transform, const double a, const double b, const double c, const double d, const double e, const int nb_forage, ivector day_layer, ivector night_layer, const int jday, const int t, const int i, const int j);
	//void SST_Habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, int sp, const int fpos, const int pop_built);
	//void Hsst_comp(CMatrices& mat, const PMap& map, dmatrix& Hs, double a, double b, const int t_count);
	//double Hsst_comp_elem(CMatrices& mat, const double a, const double b, const int t, const int i, const int j);

	void Juvenile_Habitat(VarParamCoupled& param,CMatrices& mat,const PMap& map,dvar_matrix& Hs, int sp, const int t_count);
	void Juvenile_Habitat_cannibalism(VarParamCoupled& param,CMatrices& mat,const PMap& map,dvar_matrix& Hs, dvar_matrix& total_pop, int sp, const int t_count);
	void Hj_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hj, double a, double b, const int t);
	void Hj_cannibalism_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hj, const dmatrix& total_pop, double a, double b, double c, const int t);
	void Faccessibility(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const int sp, const int jday, const int t_count, const int pop_built, const int tags_only, const ivector tags_age_solve);
	void Vars_at_age_precomp(CParam& param, const int sp);
	double Topt_at_age_comp(CParam& param, const double teta_min, const double teta_max, const int sp, const int age);
	void Faccessibility_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, double teta_max, double oxy_teta, double oxy_cr, const int sp, const int age, const int jday, const int t);
	void Average_currents(VarParamCoupled& param, VarMatrices& mat, const PMap& map, int age, const int t_count, const int pop_built);
	void Average_currents_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const int age, const int t);
	double Tmean_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const int sp, const int age, const int t);
	void Feeding_Habitat(VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Ha, int sp, int age, const int jday, const int t_count, const int migration_flag);
	void Hf_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Hf, const int sp, const int age, const int jday, const int t);
	void Feeding_Habitat_Index(VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Ha, int sp, int age, const int jday, const int t_count);
	void Seasonal_Habitat_Index(VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Hs, dvar_matrix& Ha, int sp, int age, const int jday, const int t_count);
	void Ha_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const dmatrix Hs, dvar_matrix& Ha, const int sp, const int jday);	
	void Seasonal_switch(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const int jday, int sp);
	void Seasonal_switch_comp(VarParamCoupled& param, VarMatrices& mat, const PMap& map, double season_peak, double season_start, const int jday, const int sp);
	void Seasonal_switch_year_precomp(CParam& param, CMatrices& mat, const PMap& map, double season_peak, double season_start,const int sp);

	void Mortality_Sp( VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& M, dvar_matrix& H, int sp, double mean_age_in_dtau, const int age, const int t_count);
	void M_sp_comp(const PMap& map, dvar_matrix& M, const dmatrix& H, double, double, double, double, double, double, const int dtau);
	void M_PH_juv_comp(VarParamCoupled& param, const PMap& map, CMatrices& mat, dvar_matrix& M, const dmatrix& PH, double mean_age_in_dtau);

	//void Diff_func(const PMap& map, dvar_matrix& F, dvar_matrix& Ha, dvariable& c_adv_fish, dvariable& chi_species);
	//void diff_func_comp(const PMap& map, dmatrix& F, dmatrix& H, double c, double move_par);

	void allocate_dvmatr(const int imin, const int imax, const ivector jinf, const ivector jsup){
		dvmatr1.allocate(imin,imax,jinf,jsup); dvmatr1.initialize();
		dvmatr2.allocate(imin,imax,jinf,jsup); dvmatr2.initialize();
		dvmatr3.allocate(imin,imax,jinf,jsup); dvmatr3.initialize();
		dvmatr4.allocate(imin,imax,jinf,jsup); dvmatr4.initialize();
		dvmatr5.allocate(imin,imax,jinf,jsup); dvmatr5.initialize();
		dvmatr6.allocate(imin,imax,jinf,jsup); dvmatr6.initialize();
		dvmatr7.allocate(imin,imax,jinf,jsup); dvmatr7.initialize();
		dvmatr8.allocate(imin,imax,jinf,jsup); dvmatr8.initialize();
	}
	dvariable adv_diff(const double H, dvariable& c)
	{
		dvariable function = 1 - (H /(c + H ));
		return function;
	}

	double elapsed_time_reading;
	void time_reading_init() {elapsed_time_reading = 0;}
	//double get_time_reading() {return elapsed_time_reading;}
private:
	dvar_matrix dvmatr1, dvmatr2, dvmatr3, dvmatr4, dvmatr5, dvmatr6, dvmatr7, dvmatr8;

};
#endif
