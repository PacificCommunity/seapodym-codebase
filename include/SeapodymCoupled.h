#ifndef __SeapodymCoupled_h__
#define __SeapodymCoupled_h__

#include "SeapodymDocConsole.h"
#include "Date.h"

/*!
\brief The main simulation class
*/

class SeapodymCoupled : public SeapodymDocConsole
{
public:
	SeapodymCoupled(){/*DoesNothing*/};
	SeapodymCoupled(const char* parfile) 
	{
		param = new VarParamCoupled();
		param->init_param();
		EditRunCoupled(parfile);
		elapsed_time_reading = 0;
	}
	virtual ~SeapodymCoupled(){/*DoesNothing*/};

friend class tag_release;

	int nvarcalc() const { return param->nvarcalc();};
	void xinit(dvector& x, adstring_array& names) { param->xinit(x, names); }
	double run_coupled(dvar_vector x, const bool writeoutputfiles = false) { return OnRunCoupled(x, writeoutputfiles); }
	double run_habitat(dvar_vector x, const bool writeoutputfiles = false) { return OnRunHabitat(x, writeoutputfiles); }
	double run_density(dvar_vector x, const bool writeoutputfiles = false) { return OnRunDensity(x, writeoutputfiles); }
	double run_flux(dvar_vector x, const bool writeoutputfiles = false) { return OnRunFlux(x, writeoutputfiles); }			
	dvariable reset(dvar_vector x) { return param->reset(x); }

	void write(const char* parfile) {
		param->write(parfile);
	}
	void save_statistics(const string dirout, const adstring_array x_names, double likelihood, dvector g, double elapsed_time, int status, int iter, int nvars){
		param->save_statistics(dirout, x_names, likelihood,g,elapsed_time,status,iter,nvars);
	}

	int EditRunCoupled(const char* parfile);
	double OnRunCoupled(dvar_vector x, const bool writeoutputfiles = false);
	double OnRunFlux(dvar_vector x, const bool writeoutputfiles = false);
	void OnSimulationEnd();

	double OnRunHabitat(dvar_vector x, const bool writeoutputfiles = false);
	void ReadHabitat();

	double OnRunDensity(dvar_vector x, const bool writeoutputfiles = false);
	void ReadDensity();	

	void OnRunFirstStep();
	void OnBuildForage();

	double get_total_time_reading(){
		double total_time = elapsed_time_reading
				  + pop.elapsed_time_reading
				  + func.elapsed_time_reading
				  + param->elapsed_time_reading; 
		return total_time;
	}
	//optimization control:
	int get_maxfn(){return param->maxfn;}
	double get_crit(){return param->crit;}

private:
	double dnum1, t_yrdd;
	double cell_area;
	int day, month, qtr, year, past_month, past_qtr, jday_run, jday_spinup; //date variables
	int nbt_building, nbt_spinup_tuna, nbt_forecast, nbt_total, nbt_start_series;
	//dvar_matrix F_required;
	//dmatrix F_available;
	ivector a0_adult, aN_adult; //indices of adult cohorts
	dmatrix mean_age_cohort;
	dmatrix swa; //average selectivities

	double eF_sum;
	//Fluxes comp
	int fluxes_dt_qtr;
	int fluxes_between_polygons;
	dvar4_array Density_region;	

	//tags
	d3_array rec_obs;
	d3_array tlib_obs;
	dvar3_array rec_pred;
	dmatrix rec_obs_like;
	dvar_matrix rec_pred_like;
	int nx_obs, ny_obs;  //dimension of tagging data grid
	int xr_tags, yr_tags;//ratio of tagging data to model grid resolutions
	int nb_tagpops;
	dvector xlon, ylat;
	ivector  tags_age_habitat;
	imatrix nb_rel;
	ivector t_count_rec;
	i3_array tagpop_age_solve;
	

	double catch_scaling_factor;
	double lflike; // double value of lf_like

	void get_catch_lf_like(dvariable& likelihood);
	double get_stock_like(dvariable& total_stock, dvariable& likelihood);
	double get_tag_like(dvariable& likelihood, bool writeoutputs);
	void getDate(int& jday);
	void SaveDistributions(const int year, const int month);
	void UpdateTimeVars(int& nbt_total, int& nbt_start_series);
	///void SaveJuvCohorts(string fileJuv, int sp, bool FileMode);
	///void SaveAdultCohorts(string fileAdu, int sp, bool FileMode);
	void SaveCohortsDym(int sp, bool WriteHeader, dvector zlevel);
	void SaveOneCohortDym(int sp, bool WriteHeader, dvector zlevel);
	void SaveCohorts(string fileout, int sp, bool FileMode);
//	void SaveAdultsTuna(string fileAdu, int sp, bool FileMode);
	void SaveIntermediate(const int sp, const int age);
	///void RestoreDistributions(ivector nb_juv_built, ivector nb_age_built);
	void RestoreDistributions(ivector& nb_age_built);
	void WriteFileHeaders();
	void WriteFileHeaders_submodel(const string fileout);
	void WriteOutput(int t, bool fishing);
	void InitFileFluxes();
	void WriteFluxes();
	void Set_density_region_age_zero(const int sp);
	void AverageCurrents(int t, int n);
	dmatrix DtoBcell(const dmatrix var);
	void SolveADE(d4_array var, int n, int ntimes);
	void SolveADRE(d3_array var, int n);
	void CalcSums();
	void CalcMeanTemp(const int t_count, const int tcur);
	void ConsoleOutput(int flag,double like);
	void InitializeAll();
	void ReadTimeSeriesData(int t, int t_series);
	void ReadClimatologyData(int t, int month);
	void ReadClimatologyOxy(int t, int t_clm);
	void ReadTaggingData(imatrix& nb_rel, ivector& t_count_rec);
	void get_tag_releases(dvar4_array& density, i3_array& tagpop_age_solve,imatrix& nb_rel);
	void create_tag_recaptures();
	void delete_tag_releases();
	void gaussian_kernel(dmatrix& gauss_kernel, dvector x, dvector y, double lon, double lat, double rx, double ry);
	double lon_distance(const double lon_rel, const double lon_rec, const double lat_rel, const double lat_rec);
	void ReadAll();
	void UnitConversions(int t);
	void Food_Requirement_Index(dvar_matrix& IFR, dvar_matrix FR_pop, dvar_matrix ISR_denom, const int sp, const int age, const int t_count, const int jday);
	void IFR_age_comp(dvar_matrix& IFR, dmatrix FR_pop, dmatrix ISR_denom, const int age, const int sp, const int t);
	void FR_pop_comp(dvar_matrix& FR_pop, const int sp);
	void ISR_denom_comp(dvar_matrix& ISR_denom, const int sp, const int t_count);

	void FluxesComp(dvar_matrix Density, dvar_matrix Habitat, dvar_matrix Mortality, const int age, const bool fishing, const int year, const int month, const int jday, const int step_fishery_count, const int tcur);

	void FluxesComp_polygons(dvar_matrix Density, dvar_matrix Habitat, dvar_matrix Mortality, const int age, const bool fishing, const int year, const int month, const int jday, const int step_fishery_count, const int tcur);
	
	//autodif function
	//void starvation_penalty(dvar_matrix& mortality, dvar_matrix& total_pop,dvar3_array& nF_ratio,const int sp, const int age);
	void PredationMortality(int t, dmatrix total_pop);//for forage simulation

	dvariable like(const int sp, const int k, const int f, const int nobs);
	
	void Total_Pop_comp(dvar_matrix& total_pop, const int sp, const int jday, const int t_count);
	void Total_Stock_comp(dvariable& total_stock, const int sp);
	void SpawningBiomass_comp(dvar_matrix& total_pop, const int sp);

	///void Age_class_survivals(dvar_matrix& N_a, const dmatrix N_a_1, const dmatrix Mort_a, const dmatrix Mort_a_1, const int nt_a, const int nt_a_1);
	//void Age_class_survivals(dvar_matrix& N_a, const dmatrix N_a_1, const int nt_a, const int nt_a_1);
	///void Survival(dvar_matrix& N_a, dvar_matrix& N_a_1, dvar_matrix& Mort_a, dvar_matrix& Mort_a_1, const int a, const int sp, const bool adult);
	void Survival(dvar_matrix& N_a, dvar_matrix& N_a_1, const int a, const int sp);
	void Ageing(dvar_matrix& N_a, dvar_matrix& N_a_1);
        void AgePlus(dvar_matrix& N_a, dvar_matrix& N_a_1);

	void Spawning(dvar_matrix& J, dvar_matrix& Hs, dvar_matrix& N_a, const int jday, const int sp, const int pop_built, const int t_count);
	void spawning_spinup_comp(dmatrix& J, dmatrix& Hs, dmatrix& SST, double nb_recruitment, double Tmin);
	void spawning_built_comp(dmatrix& J, dmatrix& Hs, const dmatrix  N, double nb_recruitment, double a_adults_spawning);
	double tetafunc(const double teta, const double arg);
	int get_nt(const int a, const int sp, const bool adult);
	//dvariable tetafunc(dvariable arg);

	double elapsed_time_reading;
};

/*!
\brief Class handling tag releases.
\details It stores and returns the information about release position and fish age at release.
*/
class tag_release
{
  public:
	tag_release() {/*DoNothing*/};
	virtual ~tag_release() {/*DoNothing*/};
  private:
    int i,j,a;

  public:
	int get_i(){return i;}
	int get_j(){return j;}
	int get_age(void){return a;};

	void set_release(int ii, int jj, int aa){i = ii; j = jj; a = aa;}
};



#endif
