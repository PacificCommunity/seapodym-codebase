#ifndef __VarParamCoupled_h__
#define __VarParamCoupled_h__

#include <cstdlib>
#include "XMLDocument.h"
#include "Param.h"
#include "ReadWrite.h"


/*!
\brief Seapodym DVAR parameter class.
\details In this class we read the XML parameter file, initialize and reset variable parameters.
*/
class VarParamCoupled : public CParam
{
public:
	VarParamCoupled() {
		_nvarcalc = 0;
		_gradcalc = true;
		_saruns   = false;
	};
	virtual ~VarParamCoupled() {/*DoNothing*/};
private:
	int _nvarcalc;
	bool _gradcalc;
	bool _saruns;
	XMLDocument doc;
	CReadWrite rw;

	dvector dvarpars;
	dvector dvarpars_min, dvarpars_max;

	void large_pen(dvariable x, double xmin, double xmax, dvariable& penalty, const double eps);
	double Bard_pen(double x, double xmin, double xmax, double diff);

	void par_init(dvar_vector& dvarsCoef, dvector& coef, const double coef_min, const double coef_max, 
			string s, dvector& x, adstring_array& x_name, int& idx);
	void par2_init(dvar_vector& dvarsCoef, dvector& coef, const double coef_min, const double coef_max, 
			string s1, string s2, dvector& x, adstring_array& x_name, int& idx);
	void par_read(double& varmin, double& varmax, string s, const double fixmin, const double fixmax);
	void par_read_bounds(dvector& var, double& var_min, double& var_max, string s, int& nni);
	void par2_read_bounds(dvector& var, double& var_min, double& var_max, string s1, string s2, int& nni);


public:
	int nvarcalc() const { return _nvarcalc; }
	bool gcalc(){ return _gradcalc; }
	void set_gradcalc(bool flag){ _gradcalc = flag; }	
	bool scalc(){ return _saruns; }
	void set_scalc(bool flag){ _saruns = flag; }	
	void xinit(dvector& x, adstring_array& x_names);
	dvariable reset(dvar_vector x);

	void getparam(void);
	double get_parval(int idx){ return dvarpars[idx]; }
	dvector get_parvals(void){ return dvarpars; }

	void outp_param(adstring_array x_names, const int nvars);
	void get_param_index(ivector& ix, dmatrix& xy, dmatrix& pars);
	double par_init_lo(int ix, double eps);
	double par_init_up(int ix, double eps);
	double par_init_step(int ix, double delta);
	double par_init_step_left(int ix);
	double par_init_step_right(int ix);
	void set_all_false(string* pnames);
	int set_var_parameters(ivector phase_par_flags, string* pnames);

	bool read(const string& parfile);
	void re_read_varparam();

	void write(const char* parfile) {
		getparam();
		doc.write(parfile);
	}
//	void rbin_input2d(string file_in, const imatrix& carte, DMATRIX& mat2d, int nbi, int nbj, int nbytetoskip){
//		rw.rbin_input2d(file_in, carte, mat2d, nbi, nbj, nbytetoskip);
//	}
//	void rbin_input2d(string file_in, DMATRIX& mat2d, int nbi, int nbj, int nbytetoskip);
//	void rbin_mat2d(string file_in, const imatrix& carte, DMATRIX& mat2d, int nlat, int nlong, int nbytetoskip);
//	void rbin_mat2d(string file_in, DMATRIX& mat2d, int nlat, int nlong, int nbytetoskip);

/*
	void allocate_dvmatr(const int imin, const int imax, const ivector jinf, const ivector jsup){
		dvarsU.allocate(imin, imax, jinf, jsup);dvarsU.initialize();	
		dvarsV.allocate(imin, imax, jinf, jsup);dvarsV.initialize();
	}
*/
	void save_statistics(const string dirout, const adstring_array x_names, double likelihood, dvector g, double elapsed_time, int status, int iter, int nvars);

public:
//	dvar_matrix dvarsU;
//	dvar_matrix dvarsV;

//1. dv_mortality_sp.cpp:
	double Mp_mean_max_min;
	double Mp_mean_max_max;
	dvar_vector dvarsMp_mean_max;

//2. dv_mortality_sp.cpp:
	double Mp_mean_exp_min;
	double Mp_mean_exp_max;
	dvar_vector dvarsMp_mean_exp;

//3. dv_mortality_sp.cpp:
	double Ms_mean_max_min;
	double Ms_mean_max_max;
	dvar_vector dvarsMs_mean_max;

//4. dv_mortality_sp.cpp:
	double Ms_mean_slope_min;
	double Ms_mean_slope_max;
	dvar_vector dvarsMs_mean_slope;

//5. dv_mortality_sp.cpp:
	double M_mean_range_min;
	double M_mean_range_max;
	dvar_vector dvarsM_mean_range;

//6. dv_spawning_habitat.cpp, dv_feeding_habitat.cpp, dv_juvenile_habitat.cpp:
        double a_sst_spawning_min;
        double a_sst_spawning_max;
        dvar_vector dvarsA_sst_spawning;

//7. dv_spawning_habitat.cpp, dv_feeding_habitat.cpp, dv_juvenile_habitat.cpp:
        double b_sst_spawning_min;
        double b_sst_spawning_max;
        dvar_vector dvarsB_sst_spawning;

//8.1-8.2 dv_spawning_habitat.cpp, dv_juvenile_habitat.cpp:
        double a_sst_larvae_min;
        double a_sst_larvae_max;
        dvar_vector dvarsA_sst_larvae;

	double b_sst_larvae_min;
        double b_sst_larvae_max;
        dvar_vector dvarsB_sst_larvae;

//9-11. dv_spawning_habitat.cpp, dv_feeding_habitat.cpp:
        double alpha_hsp_prey_min;
        double alpha_hsp_prey_max;
        dvar_vector dvarsAlpha_hsp_prey;

	double alpha_hsp_predator_min;
        double alpha_hsp_predator_max;
        dvar_vector dvarsAlpha_hsp_predator;

	double beta_hsp_predator_min;
        double beta_hsp_predator_max;
        dvar_vector dvarsBeta_hsp_predator;

//12. dv_accessibility.cpp:
        double a_sst_habitat_min;
        double a_sst_habitat_max;
        dvar_vector dvarsA_sst_habitat;

//13. dv_accessibility.cpp:
        double b_sst_habitat_min;
        double b_sst_habitat_max;
        dvar_vector dvarsB_sst_habitat;

//14. dv_accessibility.cpp:
        double T_age_size_slope_min;
        double T_age_size_slope_max;
        dvar_vector dvarsT_age_size_slope;

//15-17.dv_accessibility.cpp	
        dvector thermal_func_delta_min;
        dvector thermal_func_delta_max;
        dvar_matrix dvarsThermal_func_delta;

//18. dv_feeding_habitat.cpp:
        double a_oxy_habitat_min;
        double a_oxy_habitat_max;
        dvar_vector dvarsA_oxy_habitat;

//19. dv_feeding_habitat.cpp:
        double b_oxy_habitat_min;
        double b_oxy_habitat_max;
        dvar_vector dvarsB_oxy_habitat;

//20-25. dv_feeding habitat.cpp:
        dvector eF_habitat_min;
        dvector eF_habitat_max;
        dvar_matrix dvarsEF_habitat;

//26. dv_juvenile_habitat.cpp:
        double hp_cannibalism_min;
        double hp_cannibalism_max;
        dvar_vector dvarsHp_cannibalism;

//27. dv_caldia.cpp:
        double sigma_species_min;
        double sigma_species_max;
        dvar_vector dvarsSigma_species;
                                                                               
//28. dv_caldia.cpp:
        double MSS_species_min;
        double MSS_species_max;
        dvar_vector dvarsMSS_species;

//29. dv_caldia.cpp:
        double MSS_size_slope_min;
        double MSS_size_slope_max;
        dvar_vector dvarsMSS_size_slope;

//30. dv_caldia.cpp:
        double c_diff_fish_min;
        double c_diff_fish_max;
        dvar_vector dvarsC_diff_fish;

//31. dv_spawning.cpp:
        double nb_recruitment_min;
        double nb_recruitment_max;
        dvar_vector dvarsNb_recruitment;

//32. dv_spawning.cpp
        double a_adults_spawning_min;
        double a_adults_spawning_max;
        dvar_vector dvarsA_adults_spawning;

//33. dv_feeding_habitat.cpp
        double spawning_season_peak_min;
        double spawning_season_peak_max;
        dvar_vector dvarsSpawning_season_peak;

//34. dv_feeding_habitat.cpp
        double spawning_season_start_min;
        double spawning_season_start_max;
        dvar_vector dvarsSpawning_season_start;

//35-. dv_predicted_catch.cpp, dv_calrec_precalrec.cpp:
	dmatrix q_sp_fishery_min;
	dmatrix q_sp_fishery_max;
	dvar_matrix dvarsQ_sp_fishery;

//33-. dv_predicted_catch.cpp, dv_calrec_precalrec.cpp:
	dmatrix s_slope_sp_fishery_min;
	dmatrix s_slope_sp_fishery_max;
	dmatrix s_asympt_sp_fishery_min; 
	dmatrix s_asympt_sp_fishery_max;
	dvar_matrix dvarsSslope_sp_fishery;
	dvar_matrix dvarsSlength_sp_fishery;
	dvar_matrix dvarsSasympt_sp_fishery;

// like.cpp
	dvar_matrix dvarsLike_param;
	dvar_matrix dvarsProb_zero;
};
#endif
