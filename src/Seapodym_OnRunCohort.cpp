#include "SeapodymCohort.h"

void update_density_like(dvar_matrix& Density_pred, const dmatrix density_input, const imatrix map_carte, const int nlon, const int nlat, const int nlon_input, const int nlat_input, dvariable& likelihood);
//void SeapodymCohort::prerun_model(int age_start, int t_start, DMATRIX state_start)
void SeapodymCohort::prerun_model(dvar_vector x)
{
	OnRunFirstStep();
	InitializeCohort(x, false);

}

///This is the main loop for the model without fishing and fitting of density. 
///Similar to the default function, it includes the following calls:
///1- Initialising population density 
///2- Reading all forcing data (once in optimization mode, at every time step in simulation mode) 
///3- Age/lifestage loop calling the ADRE solvers and ageing
///4- Computing model density -> used as predictions
///5- Likelihood computation using input density as observations
///6- Writing outputs (in simulation mode only)
///Note, there is no modelling of tagged cohorts here!

//////////////////////////////////////////////////////////////////
//--------------------------------------------------------------//
//		     FORAGE-TUNA SIMULATION			//
//--------------------------------------------------------------//
//////////////////////////////////////////////////////////////////
/*!
\brief The tuna population simulation without fishing and density fitting.
*/
double SeapodymCohort::OnRunCohort(dvar_vector x, const bool writeoutputfiles)
{

	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	////------------------------------------------------------------/////
	////								/////
	////		||| START OF SIMULATION CYCLE |||		/////
	////								/////
	////------------------------------------------------------------/////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	for (;t_count <= nbt_cohort; t_count++)
	{
		stepForward(writeoutputfiles);
	} // end of simulation loop

	return 0;
}



