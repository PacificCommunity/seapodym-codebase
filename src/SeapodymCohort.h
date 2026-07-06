#ifndef __SeapodymCohort_h__
#define __SeapodymCohort_h__

#include <cstdlib>
//#include <fvar.hpp>
#include "XMLDocument2.h"
#include "Param.h"
#include "ReadWrite.h"
#include "SeapodymCoupled.h"
#include "DistDataCollector.h"

class DataProvider; //just a declaration - only a pointer will be stored


class SeapodymCohort : public SeapodymCoupled
{
public:
	SeapodymCohort(){/*DoesNothing*/};
	// aPlusOn selects whether the A+ (plus group) bin is carved out of the
	// normal diagonal cohort scheme (default, matches main_cohort.cpp's
	// default) or folded back in as an ordinary aging cohort, reproducing
	// pre-A+ behavior. Exposed so main_cohort.cpp's -no-aplus flag can
	// request the old behavior for regression comparison.
	SeapodymCohort(const char* parfile, int cohortId, bool aPlusOn = true) : SeapodymCoupled(parfile) {
		cohort_id = cohortId;
		aPlusEnabled = aPlusOn;
		// nb_age_class excludes the A+ (plus group) bin when aPlusEnabled: the
		// diagonal cohort-task scheme only ages "normal" cohorts through
		// nb_age_class steps; the A+ bin (age index sp_nb_cohorts[0]-1) is
		// handled separately as its own chain of accumulator tasks (see
		// main_cohort.cpp / SeapodymCohortDependencyAnalyzer). When disabled,
		// nb_age_class reverts to sp_nb_cohorts[0], i.e. A+ is just the last
		// ordinary aging cohort, as before this feature existed.
		nb_age_class = param->sp_nb_cohorts[0] - (aPlusEnabled ? 1 : 0);

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

	double time_overhead;	
	//double run_cohort(dvar_vector x, const bool writeoutputfiles = false) { return OnRunCohort(x, writeoutputfiles); }		
	void init_cohort(dvar_vector x, DistDataCollector& dataCollector, const bool writeoutputfiles = false) { return InitializeCohort(x, dataCollector, writeoutputfiles); }		
	void prerun_model();
	void OnRunFirstStep();
	double Checksum();
	std::vector<double> GetCohortDensity();

	// Chunk-id formula for a normal (non-A+) cohort task, factored out as a
	// static so callers (e.g. main_cohort.cpp) can look up the chunk of *any*
	// task/step pair, not just "this" cohort's own.
	static int computeChunkId(int taskId, int step, int numAgeGroups) {
		int row = taskId - numAgeGroups + 1 + step;
		int col = step;
		if (taskId<numAgeGroups && row==0){
			col = numAgeGroups - taskId - 1;
		}
		return row * numAgeGroups + col;
	}
	int getChunkId(int step) {
		return computeChunkId(cohort_id, step, nb_age_class);
	}

	// Chunk-id formula for an A+ (plus group) accumulator task. A+ tasks are
	// appended after all the normal (taskId, step) chunks, one per time step:
	// slot = numAgeGroups*numTimeSteps + (taskId - firstAPlusId).
	static int computeAPlusChunkId(int taskId, int firstAPlusId, int numAgeGroups, int numTimeSteps) {
		return numAgeGroups * numTimeSteps + (taskId - firstAPlusId);
	}

	void restart(int cohortId){
		cohort_id = cohortId;
		nb_age_class = param->sp_nb_cohorts[0] - (aPlusEnabled ? 1 : 0); // see constructor
		// Get starting age_class and start time from cohort_id
		if (cohort_id >= nb_age_class){
			age_start = 0;
			tstart_cohort = cohort_id-nb_age_class+1;
		}else{
			age_start = nb_age_class - cohort_id - 1;
			tstart_cohort = 0;
		}
		t_count = tstart_cohort+1;
	}

	// Sets this object up to represent the A+ (plus group) accumulator task
	// for absolute time index t (t=0 is the initial condition, before any
	// time step has elapsed). age/age_start are pinned at nb_age_class (one
	// past the last normal age, i.e. the A+ bin) rather than advancing:
	// unlike restart(), each A+(t) is an independent, one-shot task/dispatch
	// (see SeapodymCohortDependencyAnalyzer), not a persisting trajectory, so
	// there is no "next age" to advance into. tstart_cohort/t_count follow
	// the same convention restart() uses for a cohort "born" at time t, which
	// is what makes stepForward()'s existing calendar-date and forcing-data
	// lookups (both keyed off tstart_cohort/t_count) come out correct
	// without any further changes to stepForward() itself.
	void restartAPlus(int t){
		age_start = nb_age_class;
		age = nb_age_class;
		tstart_cohort = t;
		t_count = t + 1;
	}

	// Initializes the A+ task's density either from the actual initial
	// condition (seedFromFile=true, used only for t=0, which has no upstream
	// task dependency) or from mergedDensity - the previous A+ pool plus the
	// individuals that just graduated into the top age class, already summed
	// by the caller, flattened in the same map.imin1..imax1/jinf1..jsup1
	// order as GetCohortDensity(). Either way, stepForward() then runs the
	// same adult dynamics (mortality, movement, feeding habitat) on it that
	// any other cohort's step would get.
	void init_cohort_aplus(dvar_vector x, const std::vector<double>& mergedDensity, bool seedFromFile) {
		return InitializeAPlus(x, mergedDensity, seedFromFile);
	}

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
	int nb_age_class;
	bool aPlusEnabled;

	DataProvider* dp_ = nullptr;

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

	void InitializeCohort(dvar_vector& x, DistDataCollector& dataCollector, const bool writeoutputfiles = false);
	void InitializeAPlus(dvar_vector& x, const std::vector<double>& mergedDensity, bool seedFromFile);
	// Setup shared by InitializeCohort() and InitializeAPlus(): precomputed
	// per-age habitat/mortality parameters, calendar date, O2 climatology,
	// and the age/past_month/past_qtr bookkeeping stepForward() relies on.
	// Factored out so both initialization paths stay in lockstep instead of
	// risking drift between two copies of the same ~15 lines.
	void FinishInitialize(bool writeoutputfiles);

public:
	void stepForward(const bool writeoutputfiles = false);
	// Remaining to implement
	void setStateFromArray(const std::vector<double>& array);
	std::vector<double> getArrayFromState();
	void save(const std::string& restartFile);

	void setDataProvider(DataProvider* dp) {dp_ = dp;}

	void setShmForcing();
	void getData(bool spawning_habitat_only = false);
	void getO2clm(int t_clm);
};
#endif
