#include "calpop.h"

void CCalpop::precalrec(PMap& map, const DMATRIX& mortality)
{
	bm = b;

	for (int j = a.rowmin(); j <= a.rowmax(); j++)
	{
		for (int i = a[j].indexmin(); i <= a[j].indexmax(); i++)
		{
			bm[j][i] += mortality[i][j];
		}
	}
	xbet_comp(map);
} 
/*
void CCalpop::precalrec(VarParamCoupled& param, PMap& map, CMatrices& mat, dvar_matrix& mortality, bool flag_fishing, const int sp, const int age)
{
	dvarsBM = dvarsB;

	int nb_fishery = param.nb_fishery;
  
	for (int j = b.rowmin(); j <= b.rowmax(); j++){
		for (int i = b[j].indexmin(); i <= b[j].indexmax(); i++){
			if (map.carte(i,j)){
	
				if ( flag_fishing ){
					// 1- FISHING MORTALITY COEFFICIENT
					dvariable fishing_mortality = 0;
					int k = 0;
					for (int f=0; f<nb_fishery; f++){ 
						if ( param.mask_fishery_sp[sp][f]){
							dvariable catchability = param.dvarsQ_sp_fishery[sp][k];
							double selectivity  = param.selectivity_sp_fishery_age[sp][k][age];
							fishing_mortality += mat.effort[f][i][j]*catchability*selectivity*mat.lat_correction(j);
							k++;					
						}
					}			
					// 2- TOTAL MORTALITY COEFFICIENT				
					mortality[i][j] += fishing_mortality;
				}
				dvarsBM[j][i] += mortality[i][j];
			}
		}
	}

	Xbet_comp(map,2*iterationNumber);

} // END of precalrec()


void CCalpop::precalrec(PMap& map, dvar_matrix& mortality)
{
	dvarsBM = b;
	for (int j = a.rowmin(); j <= a.rowmax(); j++) {
		for (int i = a[j].indexmin(); i <= a[j].indexmax(); i++) {
			dvarsBM[j][i] += mortality[i][j];
		}
	}
	Xbet_comp1(map,2*iterationNumber);
} // END of precalrec()
*/

