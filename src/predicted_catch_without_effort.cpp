#include "calpop.h"

///Forward functions for:
///predicting catch by fishery without using the effort data.
///They compute local proportions of catch by fishery in the total catch over 'no effort' fisheries. 
///Note, the catches being used need to have the same units.

void CCalpop::Ctot_proportion_fishery_comp(const PMap& map, CParam& param, CMatrices& mat, CReadWrite& rw, const int year, const int month, const int sp)
{	
	const int nb_fishery = param.get_nbfishery();
	dmatrix catch_obs_total;
	catch_obs_total.allocate(map.imin,map.imax,map.jinf,map.jsup); 
	catch_obs_total.initialize();

	mat.Ctot_proportion_fishery(sp).initialize();

	//1.get total catch first
	int k = 0;
	for (int f=0; f<nb_fishery; f++){ 
		if (param.mask_fishery_sp[sp][f]){
			if (param.mask_fishery_sp_no_effort[sp][f]){
				for (int i = map.imin; i <= map.imax; i++){	
					const int jmin = map.jinf[i];
					const int jmax = map.jsup[i];
					for (int j = jmin; j <= jmax; j++){
						if (map.carte(i,j) && mat.catch_obs(sp,k,i,j)){
							catch_obs_total(i,j) += mat.catch_obs(sp,k,i,j);
							
						}
					}
				}
			}
			k++;
		}
	}
	//2.now calculate the proportions
	k = 0;
	for (int f=0; f<nb_fishery; f++){ 
		if (param.mask_fishery_sp[sp][f]){
			if (param.mask_fishery_sp_no_effort[sp][f]){
				for (int i = map.imin; i <= map.imax; i++){	
					const int jmin = map.jinf[i];
					const int jmax = map.jsup[i];
					for (int j = jmin; j <= jmax; j++){
						if (map.carte(i,j) && catch_obs_total(i,j)){
							mat.Ctot_proportion_fishery(sp,k,i,j) = mat.catch_obs(sp,k,i,j)/catch_obs_total(i,j);
						}
					}
				}
			}
			k++;
		}
	}
}

//call in run_coupled
//Recomp_C_fishery_proportion_in_Ctot(map, *param, rw, Ctot_proportion_fishery(sp), year, month, sp, k);
void CCalpop::Recomp_C_fishery_proportion_in_Ctot(const PMap& map, CParam& param, CReadWrite& rw, dmatrix& Ctot_proportion_fishery, const int year, const int month, const int sp, const int k)
{
	const int nb_fishery_sp = param.nb_fishery_by_sp[sp];
	dmatrix catch_obs_total;
	catch_obs_total.allocate(map.imin,map.imax,map.jinf,map.jsup); 
	catch_obs_total.initialize();

	d3_array catch_obs;
	catch_obs.allocate(0,nb_fishery_sp-1);
	for (int n=0; n<nb_fishery_sp; n++){
		catch_obs(n).allocate(map.imin,map.imax,map.jinf,map.jsup); 
		catch_obs(n).initialize();
	}

	Ctot_proportion_fishery.initialize();
	int k_loc = 0;
	for (int f=0; f<nb_fishery_sp; f++){ 
		if (param.mask_fishery_sp[sp][f]){
			if (param.mask_fishery_sp_no_effort[sp][f]){
				rw.get_catch(param, catch_obs(k_loc), f, year, month, sp);
				for (int i = map.imin; i <= map.imax; i++){	
					const int jmin = map.jinf[i];
					const int jmax = map.jsup[i];
					for (int j = jmin; j <= jmax; j++){
						if (map.carte(i,j) && catch_obs(k_loc,i,j)){
							catch_obs_total(i,j) += catch_obs(k_loc,i,j);
						}
					}
				}
			}
			k_loc++;
		}
	}

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j) && catch_obs_total(i,j)){
				Ctot_proportion_fishery(i,j) = catch_obs(k,i,j)/catch_obs_total(i,j);	
			}
		}
	}
}

void CCalpop::predicted_catch_fishery_no_effort_comp(const PMap& map, CParam& param, VarMatrices& mat, const int f, const int k, const int sp, const int age)
{
	//catch_units=0 if catch is in nb; =1 if in metric tons
	const int    C_units = param.fishery_catch_units[f];
	//Note, variable weight is either =real weight in mt or =1 if predicted catch is in numbers
	const double weight  = pow(1e-3*param.weight[sp][age],C_units);
	const double area = 1.852*param.deltaX*1.852*param.deltaY; //sq.km

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			const double C_est_age = value(mat.dvarCtot_age_est(sp,age,i,j));
			const double C_proportion_f = mat.Ctot_proportion_fishery(sp,k,i,j);
			if (map.carte(i,j) && C_proportion_f){
				// total catch by fishery f in nb of fish in a grid cell 
				// (dvarCtot_age_est has units of density, Nb.fish/sq.km)
				double C_f_age_est_cell = C_est_age * C_proportion_f * area/mat.lat_correction(j);

				mat.dvarCatch_est(sp,k).elem_value(i,j) += C_f_age_est_cell * weight;
				mat.catch_est[sp][k][i][j] += C_f_age_est_cell * weight;
				
			}
		}
	}
}
