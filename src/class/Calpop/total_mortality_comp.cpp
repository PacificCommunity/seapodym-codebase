#include "calpop.h"

///Forward functions for:
///computing the sum of natural and fishing mortalities

void CCalpop::precalrec_total_mortality_comp(const imatrix carte, const dmatrix effort, dvar_matrix& mortality, const double sq, const dvector lat_correction)
{
	for (int j = b.rowmin(); j <= b.rowmax(); j++) {
		for (int i = b[j].indexmin(); i <= b[j].indexmax(); i++) {
			if (carte(i,j) && effort(i,j)) {

				// 1- FISHING MORTALITY COEFFICIENT
				double fishing_mortality = effort(i,j) * sq * lat_correction[j];

				// 2- TOTAL MORTALITY COEFFICIENT				
				mortality.elem_value(i,j) += fishing_mortality;
			}
		}
	}
}


void CCalpop::Recomp_total_mortality_comp(const PMap& map, CParam& param, CMatrices& mat, CReadWrite& rw, dmatrix& mortality, const int age, const int sp, const int year, const int month, const int step_count)
{
	int k = 0;
	const int nb_fishery = param.get_nbfishery();

	dmatrix effort;
	effort.allocate(map.imin, map.imax, map.jinf, map.jsup);
	effort.initialize();

	for (int f=0; f<nb_fishery; f++){ 
		if (param.mask_fishery_sp[sp][f]){
			//2014: no use of effort for fisheries with no or bad effort data:
			if (param.mask_fishery_sp_no_effort[sp][f]){ k++; continue;}

			double catchability = param.q_sp_fishery[sp][k]*(1.0+param.q_dyn_fishery[f]*step_count);
			const double selectivity = Selectivity(sp,f,age);
			double sq  = selectivity * catchability;

			rw.get_effort_rm(param, effort, f, year, month);
			
			for (int j = b.rowmin(); j <= b.rowmax(); j++) {
				for (int i = b[j].indexmin(); i <= b[j].indexmax(); i++) {
					if (map.carte(i,j) && effort(i,j)) {
	
						// 1- FISHING MORTALITY COEFFICIENT
						double fishing_mortality = effort(i,j) * sq * mat.lat_correction[j];
	
						// 2- TOTAL MORTALITY COEFFICIENT			
						mortality(i,j) += fishing_mortality;
					}
				}
			}		
			k++;
		}
	}
}
