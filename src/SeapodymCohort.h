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
	SeapodymCohort(const char* parfile) : SeapodymCoupled(parfile) {};
	virtual ~SeapodymCohort() {/*DoNothing*/};

	double run_cohort(dvar_vector x, const bool writeoutputfiles = false) { return OnRunCohort(x, writeoutputfiles); }		
	//void prerun_model(int age_start, int t_start, DMATRIX state_start);
	void prerun_model(int age_start, int t_start);

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

	void stepForward(bool writeoutputfiles);
	double OnRunCohort(dvar_vector x, const bool writeoutputfiles);
	void InitializeCohort(dvar_vector& x, const bool writeoutputfiles);

	// Remaining to implement
	void setStateFromArray(const std::vector<double>& array);
	std::vector<double> getArrayFromState();
	void save(const std::string& restartFile);	
};
#endif
