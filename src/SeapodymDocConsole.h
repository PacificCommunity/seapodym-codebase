#ifndef __SeapodymDocConsole_h__
#define __SeapodymDocConsole_h__

#include "SaveTimeArea.h"
#include "Numfunc.h"
#include "calpop.h"
//#include "Dimensions.h"


/*!
\brief This class derives all necessary classes for the main simulation class
*/

class SeapodymDocConsole
{
protected:
	SeapodymDocConsole() { };

//	void UpdateDisplay();
//	CCalpop pop;

public:
	virtual ~SeapodymDocConsole() {delete param;};

	CReadWrite rw;		// objet rw derive de la classe CReadWrite
	VarParamCoupled* param;	// objet param derive de la classe CParam 
	VarMatrices mat;	// objet mat derive de la classe CMatrices
	PMap map;		// objet map derive de la classe CMap
	vector<PMap> tagmaps;		// tableau de maps for tag populations
	VarSimtunaFunc func;	// objet func derive de la classe CSimtunafunc
	CNumfunc nfunc; 	// objet nfunc derive de la classe CNumfunc
	CSaveTimeArea save;	// objet save derive de la classe CSaveTimeArea

	//Dimensions* dim;

	int nbi, nbj/*, nbt*/, nlon, nlat, deltaT;
	int nlon_input, nlat_input; // different resolution of input fields, used in seapodym_habitats and seapodym_densities only
	double deltaX, deltaY;
	double SUM_CATCH;
	int nb_fishery, nb_species, nb_forage, nb_layer, tuna_spinup;
	string date_str;
	char runtype;

	int t_count;	// nombre total de pas de temps 
	int t_series;	// nombre de pas de temps dans la serie (remis a zero chaque boucle du spinup)

	double	sumP;
	DVECTOR sumF;
	DVECTOR sumFprime;
	DVECTOR sumF_area_pred;
	DVECTOR sumF_required_by_sp;
	DVECTOR mean_omega_sp;

protected:
	//void InitDisplay() {/*DoesNothing*/};
	//void OnSimulationEnd() { param->write_param(runtype); }
	void UpdateDisplay();
	CCalpop pop;


};

#endif
