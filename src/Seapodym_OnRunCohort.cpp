#include "SeapodymCohort.h"

void update_density_like(dvar_matrix& Density_pred, const dmatrix density_input, const imatrix map_carte, const int nlon, const int nlat, const int nlon_input, const int nlat_input, dvariable& likelihood);
void SeapodymCoupled::prerun_model()
{
	OnRunFirstStep();
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
	InitializeCohort(*param);
	past_month=month;
	past_qtr=qtr;

	if (!param->gcalc()){
		//need to read oxygen in case if month==past_month
		//(otherwise we may not have it for the first time steps)
		if (param->type_oxy==1 && month==past_month)
			ReadClimatologyOxy(1, month);
		//need to read oxygen in case if qtr==past_qtr 
		if (param->type_oxy==2 && qtr==past_qtr)
			ReadClimatologyOxy(1, qtr);
	}

	//----------------------------------------------//
	// 	LIKELIHOOD INITIALISATION SECTION       //
	//----------------------------------------------//	
	//Reset model parameters:
	reset(x);

	//precompute thermal habitat parameters
	for (int sp=0; sp < nb_species; sp++)
		func.Vars_at_age_precomp(*param,sp);

	//precompute seasonal switch function
	for (int sp=0; sp < nb_species; sp++){
		if (param->seasonal_migrations[sp]){
			func.Seasonal_switch_year_precomp(*param,mat,map,
						value(param->dvarsSpawning_season_peak[sp]),
						value(param->dvarsSpawning_season_start[sp]),sp);
		}
	}

	//Output DYM file name
	string fileout;
	fileout = param->strdir_output + param->sp_name[0] + "_cohort.dym";//will need the date of birth stamp
	if (writeoutputfiles){
		WriteFileHeaders_submodel(fileout,true);	
		if (!param->gcalc())
			ConsoleOutput(0,0);
	}

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



