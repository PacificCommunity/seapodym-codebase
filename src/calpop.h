#ifndef __CCalpop_h__
#define __CCalpop_h__

//#include "mpi.h"
//#include "Map.h"
//#include "VarParamCoupled.h"
#include "VarMatrices.h"
#include "VarSimtunaFunc.h"
#include "ReadWrite.h"

/*!
\brief This is main computational class: all functions and variables to solve ADR equations are here.
*/

class CCalpop
{
public:

	CCalpop(){/*DoesNothing*/};
	virtual ~CCalpop() {/*DoesNothing*/};

//public:
	void InitCalPop(CParam& param, const PMap& map);

	void precaldia(const CParam& param, const PMap& map, CMatrices& mat);
	void precaldia_comp(const PMap& map, CParam& param, CMatrices& mat, const dmatrix& habitat, const dmatrix& total_pop,double MSS, double MSS_size_slope, double sigma_species, double c_diff_fish, const int sp, const int age, const int jday);

	void Precaldia_Caldia(const PMap& map, VarParamCoupled& param, VarMatrices& mat, dvar_matrix& habitat, dvar_matrix& total_pop, const int sp, const int age, const int t_count, const int jday);
	void caldia(const PMap& map, const CParam& param, const DMATRIX& diffusion_x,  const DMATRIX& advection_x, const DMATRIX& diffusion_y, const DMATRIX& advection_y);
	void caldia_GO(const PMap& map, const CParam& param, const DMATRIX& diffusion_x,  const DMATRIX& advection_x, const DMATRIX& diffusion_y, const DMATRIX& advection_y);

	void starvation_penalty(const PMap& map,VarParamCoupled& param, VarMatrices& mat, dvar_matrix& mortality, dvar_matrix& total_pop,dvar3_array& nF_ratio, dvar_matrix& uu, const int sp, const int age);
	void precalrec(PMap& map, const dmatrix& mortality);
	void Precalrec_juv(const PMap& map, CMatrices& mat, dvar_matrix& mortality, const int t_count);
	void precalrec_juv_comp(const PMap& map, dmatrix& bm, const dmatrix& mortality);
	//void Precalrec_total_mortality_comp(const PMap& map, VarParamCoupled& param, const CMatrices& mat, CReadWrite& rw, dvar_matrix& mortality, const int age, const int sp, const int year, const int month, const int step_count);
	void Precalrec_total_mortality_comp(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, dvar_matrix& mortality, const int age, const int sp, const int t_count, const int year, const int month, const int step_count);
	void Recomp_total_mortality_comp(const PMap& map, CParam& param, CMatrices& mat, CReadWrite& rw, dmatrix& mortality, const int age, const int sp, const int year, const int month, const int step_count);
	//void precalrec_total_mortality_comp(const PMap& map, VarParamCoupled& param, const CMatrices& mat, dvar_matrix& mortality, const int age, const int sp);
	void precalrec_total_mortality_comp(const imatrix carte, const dmatrix effort, dvar_matrix& mortality, const double sq, const dvector lat_correction);

	void Precalrec_Calrec_adult(const PMap& map, VarMatrices& mat, VarParamCoupled& param, CReadWrite& rw, dvar_matrix& uu, dvar_matrix& mortality, const int t_count, const bool fishing, const int age, const int sp, const int year, const int month, const int jday, const int step_count, const int no_mortality);

	void calrec(const PMap& map, dmatrix& uu, const dmatrix& mortality);
	void calrec1(const PMap& map, dvar_matrix& uu, const dmatrix& mortality);

	void calrec_with_catch(const PMap& map, CParam& param, dvar_matrix& uu, const dmatrix& C_obs, dvar_matrix& C_est);

	void calrec_GO(const PMap& map, dvar_matrix& uu);
	void calrec_GO_with_catch(const PMap& map, CParam& param, dvar_matrix& uu, const dmatrix& C_obs, dvar_matrix& C_est);
	void Calrec_juv(const PMap& map, CMatrices& mat, dvar_matrix& uu, dvar_matrix& mortality, const int t_count);
	void Calrec_adult(const PMap& map, dvar_matrix& uu, dvar_matrix& mortality);


	void Recomp_abc_coef(const PMap& map, CMatrices& mat, const int t_count, const dmatrix& mortality, dmatrix& aa, dmatrix& bbm, dmatrix& cc);
	void Recomp_DEF_coef(const PMap& map, CParam& param, CMatrices& mat, const int t_count, const int jday, const dmatrix& habitat, dmatrix& dd, dmatrix& ee, dmatrix& ff, dmatrix& advection_x, dmatrix& advection_y, const int sp, const int age, const double MSS, const double c_diff_fish, const double sigma_species);
	void Recomp_DEF_UV_coef(const PMap& map, CParam& param, CMatrices& mat, dmatrix& u, dmatrix& v, const dmatrix& habitat, dmatrix& dd, dmatrix& ee, dmatrix& ff, dmatrix& advection_x, dmatrix& advection_y, const int sp, const int age, const double MSS, const double c_diff_fish, const double sigma_species, const int jday);
	void RecompDiagCoef_juv(const PMap& map, CMatrices& mat, const int t_count, const dmatrix mortality, dmatrix& a, dmatrix& bm, dmatrix& c, dmatrix& d, dmatrix& e, dmatrix& f);
	void RecompDiagCoef_adult(const PMap& map, CParam& param, CMatrices& mat, const int t_count, const int jday, const dmatrix& mortality, const dmatrix& habitat, dmatrix& aa, dmatrix& bbm, dmatrix& cc, dmatrix& dd, dmatrix& ee, dmatrix& ff, const int sp, const int age, const double MSS, const double c_diff_fish, const double sigma_species);
	void RecompDiagCoef_UV_adult(const PMap& map, CParam& param, CMatrices& mat, const int t_count, const int jday, const dmatrix& mortality, const dmatrix& habitat, dmatrix& aa, dmatrix& bbm, dmatrix& cc, dmatrix& dd, dmatrix& ee, dmatrix& ff, const int sp, const int age, const double MSS, const double c_diff_fish, const double sigma_species);
	void RecompM_sp(const PMap& map, const CParam& param, dmatrix& M, const dmatrix& H, const double age, const int sp);

	void Predicted_Catch_Fishery(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int sp, const int f, const int k, const int year, const int month, const int t_count, const int step_count);
	void predicted_catch_fishery_comp(const PMap& map, CParam& param, VarMatrices& mat, const int f, const int k, const int sp, const int age, const dmatrix& uu, const int step_count);

	void Total_obs_catch_age_comp(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int age, const int sp, const int year, const int month, const int t_count);
	void Ctot_proportion_fishery_comp(const PMap& map, CParam& param, CMatrices& mat, CReadWrite& rw, const int year, const int month, const int sp);
	void Recomp_C_fishery_proportion_in_Ctot(const PMap& map, CParam& param, CReadWrite& rw, dmatrix& Ctot_proportion_fishery, const int year, const int month, const int sp, const int k);
	void Total_exploited_biomass_comp(const PMap& map, VarParamCoupled& param, VarMatrices& mat, const int sp, const int t_count);
	void Selectivity_comp(CParam& param, const int nb_fishery, const int a0, const int nb_ages, const int sp);
	void Predicted_Catch_Fishery_no_effort(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int sp, const int year, const int month);
	void predicted_catch_fishery_no_effort_comp(const PMap& map, CParam& param, VarMatrices& mat, const int f, const int k, const int sp, const int age);
	void total_exploited_biomass_comp(const imatrix carte, const dmatrix& uu, const dmatrix& Cobs, const int f, const int fne, const int age, const int sp);
	void Recomp_total_exploited_biomass(const PMap& map, CParam& param, CMatrices& mat, dmatrix& EB, const dmatrix& Cobs, const dvector& selectivity, const int f, const int sp, const int t_count);
	void total_obs_catch_age_comp(const PMap& map, const CParam& param, CMatrices& mat, const dmatrix& uu, const dmatrix& Cobs, dvar_matrix& Ctot_age_obs, const int f, const int fne, const int k, const int age, const int sp, const double C2Dunits);
	void Recomp_total_obs_catch_age(const PMap& map, CParam& param, CMatrices& mat, CReadWrite& rw, dmatrix& Ctot_age_obs, const int age, const int sp, const int year, const int month, const int t_count);

	//void predicted_catch(const CParam& param, CMatrices& mat, const int sp, const int age, const DMATRIX& uu);
	//void Predicted_Catch(const PMap& map, const VarParamCoupled& param, VarMatrices& mat, const int sp, const int age, dvar_matrix& uu);
	//void Predicted_Catch(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int sp, DVAR3_ARRAY& uu, const int year, const int month/*, bool fishing*/);
//	void Predicted_Catch(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int sp, const int age, dvar_matrix& uu, const int year, const int month/*, bool fishing*/);
//	void predicted_catch_comp(const PMap& map, const CParam& param, CMatrices& mat, const int sp, const int age, const DMATRIX& uu);
	//Predicted_LF(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int sp, const int f, const int k, const int year, const int month, const int step_count);
	//void predicted_lf_comp(const PMap& map, const CParam& param, CMatrices& mat, const int f, const int k, double catchability, double s_slope, double s_length, const int sp, const int age, const int r, const int area, const dmatrix& uu);
	//void precaldia1(VarParamCoupled& param, const PMap& map, dvar_matrix& Ha, VarSimtunaFunc& function, VarMatrices& mat, const int sp, const int age);
	//void caldia(const PMap& map, const CParam& param, dvar_matrix& diffusion_x,  dmatrix& advection_x, dvar_matrix& CHI_x, dvar_matrix& diffusion_y, dmatrix& advection_y, dvar_matrix& CHI_y, dvar_matrix& habitat);
	//void caldia(const PMap& map,const CParam& param, const DMATRIX& diffusion_x,  const DMATRIX& advection_x, const DMATRIX& CHI_x, const DMATRIX& diffusion_y, const DMATRIX& advection_y, const DMATRIX& CHI_y, const DMATRIX& habitat);
	//void caldia(const PMap& map, const VarParamCoupled& param, dvar_matrix& diffusion_x,  dvar_matrix& advection_x, dvar_matrix& diffusion_y, dvar_matrix& advection_y);
	//void precalrec(VarParamCoupled& param, PMap& map, CMatrices&, dvar_matrix& mortality, bool flag_fishing, const int sp, const int age);
	//void precalrec(PMap& map, dvar_matrix& mortality);
	//void Precalrec_adult(const PMap& map, VarParamCoupled& param, const CMatrices& mat, dvar_matrix& mortality, const int fpos, const bool pop_built, const bool fishing,const int age, const int sp);
		//void Precalrec_adult1(const PMap& map, VarParamCoupled& param, const CMatrices& mat, dvar_matrix& mortality, const dmatrix& habitat, const int fpos, const bool pop_built, const bool fishing,const int age, const int sp);

	//void Precalrec_Calrec_adult(const PMap& map, VarMatrices& mat, VarParamCoupled& param, CReadWrite& rw, dvar_matrix& uu, dvar_matrix& mortality, dvar_matrix& habitat, dvar_matrix& total_pop, dvar3_array& nF_ratio, const int fpos, const int pop_built, const bool fishing, const int age, const int sp, const int year, const int month);
	//void calrecJuv(const PMap& map, dvar_matrix& uu, dvar_matrix& mortality);
	//void calrec(const PMap& map, dvar_matrix& uu, dvar_matrix& mortality);
	//void Calrec_adult1(const PMap& map, CParam& param, dvar_matrix& uu, dvar_matrix& mortality, dvar_matrix& habitat, double c_diff, double mss, double sigma, const int fpos, const int pop_built, const bool fishing, const int age, const int sp);


	int get_iterationN(){return iterationNumber;}
	int get_maxn(){return maxn;}
	int get_Vinf(){return Vinf;}
	//void set_gradcalc(bool flag){_gradcalc = flag;}

private:
	int maxn;
	int iterationNumber;
	double sigma_fcte, deltax, deltay;
	double Vinf;
	//bool _gradcalc;

	
	// Coefficients diagonaux
	// utilisés par la fonction tridag
	DMATRIX a;
	DMATRIX b;
	DMATRIX bm;
	DMATRIX c;
	DMATRIX d;
	DMATRIX e;
	DMATRIX f;

	DMATRIX xbet;
	DMATRIX ybet;

public:
	dvar_matrix dvarsA;
	dvar_matrix dvarsB;
	dvar_matrix dvarsBM;
	dvar_matrix dvarsC;
	dvar_matrix dvarsD;
	dvar_matrix dvarsE;
	dvar_matrix dvarsF;

	dvar_matrix Xbet;
	dvar_matrix Ybet;

	//movement parameters and their derivatives
	dvar_matrix Mss_species, Mss_size_slope, C_diff_fish, Sigma_species;	

	dvar3_array dvarsSNsum;
	d3_array Selectivity;

	void Xbet_comp1(const PMap& map, int dt);
	void xbet_comp(const PMap& map, dmatrix& xbet, dmatrix& a, dmatrix& bm, dmatrix& c, int dt);
	void ybet_comp(const PMap& map, dmatrix& ybet, dmatrix& d, dmatrix& e, dmatrix& f, int dt);
	//void Xbet_comp(const PMap& map, int dt);
	//void Ybet_comp(const PMap& map, int dt);

	DMATRIX uuint; 
	
	double elapsed_time_reading;
	void time_reading_init(){elapsed_time_reading = 0; }

private:
	void xbet_comp(const PMap& map);
	void ybet_comp(const PMap& map);

	void Tridag(dvar_vector& at, dvar_vector&  bet, dvar_vector& ct, dvar_vector& rhs, dvar_vector& uvec, dvar_vector& gam, int inf, int sup);
	void tridag1_0(const dvector& at, dvar_vector&  bet, const dvector& ct, dvar_vector& rhs, dvar_vector& uvec, dvar_vector& gam, int inf, int sup);
	void tridag1_1(const dvector& at, const dvector&  bet, const dvector& ct, dvar_vector& rhs, dvar_vector& uvec, dvar_vector& gam, int inf, int sup);
	void tridag2(dvar_vector& at, dvar_vector&  bet, dvar_vector& ct, dvar_vector& rhs, dvar_vector& uvec, dvar_vector& gam, int inf, int sup);
	void tridag(const DVECTOR& at, const DVECTOR& bt, const DVECTOR& ct, const DVECTOR& r, DVECTOR& uvec, DVECTOR& gam, const int inf, const int sup);
	void tridag_GO(const DVECTOR& at, const DVECTOR& bt, const DVECTOR& ct, const DVECTOR& r, DVECTOR& uvec, DVECTOR& gam, const int inf, const int sup);
	double d1(const char pos, const double sigmam, const double sigma, const double uinf, const double twodd, const double dd, const double d);
	double d2(const char pos, const double sigmam, const double sigma, const double sigmap, const double u, const double twodd, const double dd, const double d, const int dt);
	double d3(const char pos, const double sigma, const double sigmap, const double usup, const double twodd, const double dd, const double d);


};
/////////////////////////////////////////////////////////////////////////////
#endif
