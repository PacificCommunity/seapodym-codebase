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
	SeapodymCohort(VarParamCoupled* prm, int cohort_id) {
		// copy param
		param = prm;

		// Get starting age_class and start time from cohort_id
		int nb_age_class = param->sp_nb_cohort_jv[0] + param->sp_nb_cohort_ad[0];
		int quotient = cohort_id / nb_age_class;
		int remainder = cohort_id % nb_age_class;
		if (cohort_id >= nb_age_class){
			age_start = 0;
			t_start = 2 + (quotient - 1)*nb_age_class + remainder;
		}else{
			age_start = remainder;
			t_start = 1;
		}

		// Size of map (useful to access to a specific position adress of the 4D array pointer storing density)
		size_map = 0;
		const int imin = map.imin;
		const int imax = map.imax;
		for (int i = imin; i <= imax; i++){
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin ; j <= jmax; j++){
				size_map++;
			}
		}
	};
	virtual ~SeapodymCohort() {/*DoNothing*/};

	double run_cohort(dvar_vector x, const bool writeoutputfiles = false) { return OnRunCohort(x, writeoutputfiles); }		
	void prerun_model(dvar_vector x, int init_from_inputfile, DVAR4_ARRAY* array_ptr);
	void prerun_model(dvar_vector x);
	void OnRunFirstStep();
	void ReadAll();

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
	int nbt_cohort;
	int cohort_id;
	int age_start;
	int t_start;
	int nb_age_class;
	int size_map;

	dvariable likelihood;

	DMATRIX init_state;

	dvar_matrix Spawning_Habitat;
	dvar_matrix Total_pop;
	dvar_matrix Habitat; 
	dvar_matrix IFR; 
	dvar_matrix ISR_denom; 
	dvar_matrix FR_pop;
	dvar_matrix Mortality; 
	dvar_matrix dvarCohortDensity; 

	ivector  tags_age_habitat;
	
	int pop_built;

	double OnRunCohort(dvar_vector x, const bool writeoutputfiles);
	void InitializeCohort(dvar_vector& x, const bool writeoutputfiles);

public:
	void stepForward(DVAR4_ARRAY* array_ptr);
};
#endif
