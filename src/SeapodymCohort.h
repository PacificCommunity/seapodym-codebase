#ifndef __SeapodymCohort_h__
#define __SeapodymCohort_h__

#include <cstdlib>
//#include <fvar.hpp>
#include "XMLDocument2.h"
#include "Param.h"
#include "ReadWrite.h"
#include "SeapodymCoupled.h"


class SeapodymCohort : public SeapodymCoupled
{
public:
	SeapodymCohort(){/*DoesNothing*/};
	SeapodymCohort(const char* parfile, int cohortId) : SeapodymCoupled(parfile) {
		cohort_id = cohortId;
		nb_age_class = param->sp_nb_cohorts[0];

		// Get starting age_class and start time from cohort_id
		if (cohort_id >= nb_age_class){
			age_start = 0;
			tstart_cohort = cohort_id-nb_age_class+1;
		}else{
			age_start = nb_age_class - cohort_id - 1;
			tstart_cohort = 0;
		}
	};
	virtual ~SeapodymCohort() {/*DoNothing*/};

	//double run_cohort(dvar_vector x, const bool writeoutputfiles = false) { return OnRunCohort(x, writeoutputfiles); }		
	void init_cohort(dvar_vector x, const bool writeoutputfiles = false) { return InitializeCohort(x, writeoutputfiles); }		
	void prerun_model();
	void OnRunFirstStep();
	int nb_age_class;

private:
	int dtau;
	int nbt_before_first_recruitment; 	
	int nt_dtau; 
	int tcur; 
	int nbt_no_forecast;
	bool fishing;
	int migration_flag;
	int step_count;
	int step_fishery_count;
	int jday; 
	int nbstoskip; 
	int age;
	int tf_cohort;
	int cohort_id;
	int age_start;
	int tstart_cohort;

	dvariable likelihood;

	//DMATRIX init_state;//not needed, initialized from mat.init_density_species

	dvar_matrix Spawning_Habitat;
	dvar_matrix Total_pop;
	dvar_matrix Habitat; 
	dvar_matrix IFR; 
	dvar_matrix ISR_denom; 
	dvar_matrix FR_pop;
	dvar_matrix Mortality; 
	dvar_matrix dvarCohortDensity; //think if we want to preserve multi-species

	ivector  tags_age_habitat;
	
	int pop_built;

	void InitializeCohort(dvar_vector& x, const bool writeoutputfiles = false);

public:
	void stepForward(const bool writeoutputfiles = false);
	// Remaining to implement
	void setStateFromArray(const std::vector<double>& array);
	std::vector<double> getArrayFromState();
	void save(const std::string& restartFile);	
};
#endif
