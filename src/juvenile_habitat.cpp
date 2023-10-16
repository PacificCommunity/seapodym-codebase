#include "VarSimtunaFunc.h"

///Forward functions for: 
///juvenile habitat functions. This function depends on the temperature
///in the epipelagic layer, with or without cannibalism effect 
///(note, the latter has very low sensitivity usually)

const double c_zoo = 100.0; //constant in the first tests ECCO, need to make it global or parameter as soon as early life stage data will be available


//function of temperature only
void VarSimtunaFunc::Hj_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hj, double a, double b, const int t)
{
	for (int i = map.imin; i <= map.imax; i++) {
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++) {
			if (map.carte(i,j)) {

				double SST = mat.tempn[t][0][i][j];
				//double SST = mat.sst(t,i,j);

				double zoo  = mat.np1[t][i][j]+1e-10;

				double f_food = zoo*zoo/(c_zoo+zoo*zoo);

				double f_sst  = exp(-pow((SST- b),2) / (2*a*a));

				Hj.elem_value(i,j) = f_sst * f_food; 
//mat.Hj(i,j) = value(Hj(i,j));				
			}
		}
	}
}



//function of temperature and cannibalism effect
void VarSimtunaFunc::Hj_cannibalism_comp(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hj, const dmatrix& total_pop, double a, double b, double half_predators, const int t)
{
	for (int i = map.imin; i <= map.imax; i++) {
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++) {
			if (map.carte(i,j)) {

				double SST = mat.tempn[t][0][i][j];
				//double SST = mat.sst(t,i,j);

				double predators = total_pop[i][j];
		
				double f_cannibalism = 1 - predators/(half_predators+predators);
			
				double f_sst = exp(-pow((SST- b),2) / (2*a*a));

				double zoo  = mat.np1[t][i][j]+1e-10;

				double f_food = zoo*zoo/(c_zoo+zoo*zoo);

				Hj.elem_value(i,j) = f_sst * f_food * f_cannibalism; 
			
			}
		}
	}
}
