#include "SeapodymCoupled.h"

///Forward function for:
///complex non-linear function that depends on the density of adult biomass
///to derive the food requirement and compute the index of depletion in case 
///if the available micronekton does not support the daily ration. 
///This function might be useful in multi-species models.


const double slope = 0.01;

void SeapodymCoupled::IFR_age_comp(dvar_matrix& IFR, dmatrix FR_pop, dmatrix ISR_denom, const int age, const int sp, const int t)
{//IFR_age - index of food requirement by age (accessible food / required food)
 //IFR_pop - index of population food requirement (total food resource / total requirement)
 //ISR_age - index of sharing the food resource, magnifies the intensity of competition or its absence
 //residual_competition - percentage of food resource available for species given the presence of other species in the habitat (considered constant in mono-species model) 

	const double pi = 4.0*(atan(0.5)+atan(1.0/3.0));
	const double W = param->weight[sp][age]*0.001;//in mt
	const double R = param->forage_ration[sp];
	const double wrt = W * R * param->deltaT;
	const double residual_competition = param->residual_competition[sp]; 

	dvector theta(0,nb_forage-1); theta.initialize();
	const int imin = map.imin; 
	const int imax = map.imax; 
	for (int i = imin; i <= imax; i++){
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin ; j <= jmax; j++){
			if (map.carte[i][j]){
			
				double sumF = 0.0;
				double sumF_access = 0.0;
				double F_accessible = 0.0;
				double F_required   = 0.0;
				double IFR_age_pop = 0.0;

				for (int n=0; n<nb_forage; n++){
					double f_access = value(mat.dvarF_access(n,age,i,j));
					sumF_access += f_access;	
					sumF += mat.forage(t,n,i,j);
				}

				double IFR_pop = residual_competition*sumF/FR_pop(i,j);
				//scale to vary between 0 and 1
				IFR_pop = (10/pi)*atan(IFR_pop*pi/10);

				for (int n=0; n<nb_forage; n++)
					theta(n) = value(mat.dvarF_access(n,age,i,j))/sumF_access;

				//units of F: g/m^2 = tones/km^2
				for (int n=0; n<nb_forage; n++){
					F_accessible += mat.forage[t][n][i][j] * theta(n);
				}
		
				//units of tuna = Nb per sq km; 
				//units of food requirement = mt*Nb/km^2 = tones/km^2
				F_required = value(mat.dvarDensity[sp][age][i][j]) * wrt;

				//1. compute IFR of cohort
				double IFR_age = F_accessible/(F_required+1e-4);

				IFR_age_pop = IFR_age * IFR_pop;	

				//4. Finally scaled IFR, which will be used in Mortality function
				IFR.elem_value(i,j) = 1.0/(1.0+pow(slope,IFR_age_pop-1.0));				
			}
		}
	}
}

