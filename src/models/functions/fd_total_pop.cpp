#include "SeapodymCoupled.h"
#include <fvar.hpp>

/*
//simple total_pop function
void SeapodymCoupled::Total_Pop_comp(dvar_matrix& total_pop, const int sp, const int jday, const int t_count)
{ 	
	total_pop.initialize();

	for (int a=param->sp_a0_adult[sp]; a < param->sp_nb_cohorts[sp]; a++){
		for (int i = map.imin; i <= map.imax; i++){	
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++){
				if (map.carte(i,j)){

					double total_pr = value(total_pop(i,j));				
					total_pop.elem_value(i,j) = total_pr + value(mat.dvarDensity(sp,a,i,j));
				}				 
			}
		}
	}
} 
*/

void SeapodymCoupled::Total_Pop_comp(dvar_matrix& total_pop, const int sp, const int jday, const int t_count)
{ 	
	total_pop.initialize();
	for (int a=param->sp_a0_adult[sp]+1; a < param->sp_nb_cohorts[sp]; a++){

        	int age_adult_habitat = a;
	        for (int aa=a; aa>1; aa--){
                if (param->age_compute_habitat(sp,aa-1)!=param->age_compute_habitat(sp,aa)){
                        age_adult_habitat = aa;
                        break;
                	}
        	}
		for (int i = map.imin; i <= map.imax; i++){	
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++){
				if (map.carte(i,j)){
					
					double access0 = value(mat.dvarZ_access(0,age_adult_habitat,i,j));
					double total_pr = value(total_pop(i,j));				
					total_pop.elem_value(i,j) = total_pr + access0*value(mat.dvarDensity(sp,a,i,j));
				}				 
			}
		}
	}
} 


void SeapodymCoupled::SpawningBiomass_comp(dvar_matrix& total_pop, const int sp)
{ 	
	total_pop.initialize();

	dvector maturity_age = param->maturity_age[sp];
	for (int a=param->age_mature[sp]; a < param->sp_nb_cohorts[sp]; a++){
		for (int i = map.imin; i <= map.imax; i++){	
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++){
				if (map.carte(i,j)){
					double total_pr = value(total_pop(i,j));				
					total_pop.elem_value(i,j) = total_pr + maturity_age[a] * value(mat.dvarDensity(sp,a,i,j));
				}				 
			}
		}
	}
} 

void SeapodymCoupled::Total_Stock_comp(dvariable& total_stock, const int sp)
{ 	
	const double area = 1.852*deltaX*1.852*deltaY;
	const double lonmin_stock_area = param->stock_lonmin[sp];
	const double lonmax_stock_area = param->stock_lonmax[sp];
	const double latmin_stock_area = param->stock_latmin[sp];
	const double latmax_stock_area = param->stock_latmax[sp];

	int imin = param->lontoi(lonmin_stock_area);
	int imax = param->lontoi(lonmax_stock_area)-1;//-1 to match regional computations (see SumQArea)
	if (imin<map.imin) imin = map.imin;
	if (imax>map.imax) imax = map.imax;
	int jmin_stock_area = param->lattoj(latmax_stock_area);
	int jmax_stock_area = param->lattoj(latmin_stock_area)-1;
	if (jmin_stock_area<map.jmin) jmin_stock_area = map.jmin;
	if (jmax_stock_area>map.jmax) jmax_stock_area = map.jmax;


	int nbstoskip = param->nbsteptoskip;
	int nbt = nbt_total-nbstoskip;

	for (int a=param->sp_a0_adult[sp]; a < param->sp_nb_cohorts[sp]; a++){

		double Total_Stock = value(total_stock);

		const double W_mt = param->weight[sp][a] * 0.001;
		for (int i = imin; i <= imax; i++){	
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++){
				if (map.carte(i,j) && j>jmin_stock_area && j<jmax_stock_area){
					double lat_corrected_area = area/ mat.lat_correction[j];

					//Units of Total Stock = thous.mt
					Total_Stock += 1e-3*value(mat.dvarDensity(sp,a,i,j)) * W_mt *
						       lat_corrected_area/nbt;
				}				 
			}
		}
		value(total_stock) = Total_Stock;
	}
} 
