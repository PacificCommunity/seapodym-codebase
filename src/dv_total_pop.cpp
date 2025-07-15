#include "SeapodymCoupled.h"
#include <fvar.hpp>

void dv_total_pop_comp();
void dv_spawning_biomass_comp();
void dv_total_stock_comp();
double f_accessibility_layer(const double O2, const double T,double twosigsq, double temp_mean, double oxy_teta, double oxy_cr);
double f_accessibility_layer(const double O2, const double T, double temp_age, double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr);
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);

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
			
		save_identifier_string2("total_pop_begin");
		mat.dvarDensity(sp,a).save_dvar_matrix_position();
		total_pop.save_dvar_matrix_position();
		unsigned long int pmap = (unsigned long int)&map;
		save_long_int_value(pmap);
		save_identifier_string2("total_pop_end");
		
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_total_pop_comp);
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
		save_identifier_string2((char*)"total_pop_begin");
		mat.dvarZ_access(0,a).save_dvar_matrix_position();
		mat.dvarDensity(sp,a).save_dvar_matrix_position();
		total_pop.save_dvar_matrix_position();
		save_int_value(t_count);
		save_int_value(jday);
		save_int_value(a);
		save_int_value(sp);
		unsigned long int cparam = (unsigned long int)*&param;
		save_long_int_value(cparam);
		unsigned long int pmap = (unsigned long int)&map;
		save_long_int_value(pmap);
		unsigned long int cmat = (unsigned long int)&mat;
		save_long_int_value(cmat);
		save_identifier_string2((char*)"total_pop_end");
		
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_total_pop_comp);
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
		
		save_identifier_string2((char*)"spawning_biomass_begin");
		save_double_value(maturity_age[a]);
		mat.dvarDensity(sp,a).save_dvar_matrix_position();
		total_pop.save_dvar_matrix_position();
		unsigned long int pmap = (unsigned long int)&map;
		save_long_int_value(pmap);
		save_identifier_string2((char*)"spawning_biomass_end");
		
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_spawning_biomass_comp);
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

		total_stock.save_prevariable_position();
		save_identifier_string2((char*)"before_total_stock_comp");

		const double W_mt = param->weight[sp][a] * 0.001;
		//for (int i = map.imin; i <= map.imax; i++){	
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
		//total_stock = nograd_assign(Total_Stock);//attn: nograd_assign doesn't work!
		value(total_stock) = Total_Stock;//works!
			
		save_identifier_string2((char*)"total_stock_begin");
		mat.dvarDensity(sp,a).save_dvar_matrix_position();
		total_stock.save_prevariable_position();
		save_int_value(a);
		save_int_value(sp);
		save_int_value(nbt);
		unsigned long int pmap = (unsigned long int)&map;
		save_long_int_value(pmap);
		unsigned long int cparam = (unsigned long int)*&param;
		save_long_int_value(cparam);
		save_identifier_string2((char*)"total_stock_end");
		
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_total_stock_comp);
	}
} 

//adjoint for Total_Stock_comp 
void dv_total_stock_comp()
{
        verify_identifier_string2((char*)"total_stock_end");
        unsigned long int pos_param = restore_long_int_value();
        unsigned long int pos_map   = restore_long_int_value();
	const int nbt_total	    = restore_int_value();
	const int sp		    = restore_int_value();
	const int age		    = restore_int_value();
	const prevariable_position ts_pos = restore_prevariable_position();
        const dvar_matrix_position Na_pos = restore_dvar_matrix_position();
        verify_identifier_string2((char*)"total_stock_begin");

        verify_identifier_string2((char*)"before_total_stock_comp");
	const prevariable_position ts_pr_pos = restore_prevariable_position();

	double dftotal_stock = restore_prevariable_derivative(ts_pos);
	double dftotal_stock_pr = restore_prevariable_derivative(ts_pr_pos);
	dmatrix dfNa = restore_dvar_matrix_derivatives(Na_pos);

	CParam* param = (CParam*) pos_param;
        PMap* map = (PMap*) pos_map;

	const double lonmin_stock_area = param->stock_lonmin[sp];
	const double lonmax_stock_area = param->stock_lonmax[sp];
	const double latmin_stock_area = param->stock_latmin[sp];
	const double latmax_stock_area = param->stock_latmax[sp];

	int imin = param->lontoi(lonmin_stock_area);
	int imax = param->lontoi(lonmax_stock_area)-1;//-1 to match regional computations (see SumQArea)
	if (imin<map->imin) imin = map->imin;
	if (imax>map->imax) imax = map->imax;

	int jmin_stock_area = param->lattoj(latmax_stock_area);
	int jmax_stock_area = param->lattoj(latmin_stock_area)-1;
	if (jmin_stock_area<map->jmin) jmin_stock_area = map->jmin;
	if (jmax_stock_area>map->jmax) jmax_stock_area = map->jmax;
	

	const double W_mt = param->weight(sp,age) * 0.001;
	const double area = 1.852*param->deltaX*1.852*param->deltaY;

	dvector lat_corrected_area(map->jmin,map->jmax);
	lat_corrected_area.initialize();
	for (int j=map->jmin;j<=map->jmax;j++){
                double lastlat = param->lastlat(j);
                lat_corrected_area[j] = area/param->correction_lat(lastlat);
        }

	double dfsum_ij = 0.0;
	
	//total_stock = total_stock_pr + sum_ij;
	dftotal_stock_pr += dftotal_stock;
	dfsum_ij	 += dftotal_stock;
	dftotal_stock     = 0.0;

        for (int i = imax; i >= imin; i--){
                const int jmin = map->jinf[i];
                const int jmax = map->jsup[i];
                for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j) && j>jmin_stock_area && j<jmax_stock_area){
			//if (map->carte(i,j)){

				//sum_ij += 1e-3*value(mat.dvarDensity(sp,a,i,j)) * W_mt * lat_corrected_area/nbt_total;
                                dfNa(i,j)  += 1e-3*dfsum_ij*W_mt*lat_corrected_area(j)/nbt_total;
			}
                }
        }
	dftotal_stock   += dftotal_stock_pr;
	dftotal_stock_pr = 0.0;

        dfNa.save_dmatrix_derivatives(Na_pos);
        save_double_derivative(dftotal_stock,ts_pos);
        save_double_derivative(dftotal_stock_pr,ts_pr_pos);
}


/*
//adjoint for simple total_pop (without accessibility to 0 layer)
void dv_total_pop_comp()
{
        verify_identifier_string2("total_pop_end");
        unsigned long int pos_map   = restore_long_int_value();
        const dvar_matrix_position total_pop_pos  = restore_dvar_matrix_position();
        const dvar_matrix_position Na_pos         = restore_dvar_matrix_position();
        verify_identifier_string2("total_pop_begin");

        dmatrix dftotal_pop = restore_dvar_matrix_derivatives(total_pop_pos);
        dmatrix dfNa        = restore_dvar_matrix_derivatives(Na_pos);
        PMap* map = (PMap*) pos_map;

        const int imax = map->imax;
        const int imin = map->imin;

        for (int i = imax; i >= imin; i--){
                const int jmin = map->jinf[i];
                const int jmax = map->jsup[i];
                for (int j = jmax; j >= jmin; j--){
                        if (map->carte(i,j)){
                                double dftotal_pr = 0.0;

                                //total_pop(i,j) = total_pr(i,j) + mat.dvarPop_species(sp,a,i,j);
                                dftotal_pr += dftotal_pop(i,j);
                                dfNa(i,j)  += dftotal_pop(i,j);
                                dftotal_pop(i,j) = 0.0;

                                //double total_pr = total_pop(i,j);
                                dftotal_pop(i,j) += dftotal_pr;
                        }
                }
        }

        dftotal_pop.save_dmatrix_derivatives(total_pop_pos);
        dfNa.save_dmatrix_derivatives(Na_pos);
}
*/

void dv_total_pop_comp()
{
	verify_identifier_string2((char*)"total_pop_end");
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	const int sp  = restore_int_value();
	const int age = restore_int_value();
	const int jday = restore_int_value();
	const int t_count = restore_int_value();
	const dvar_matrix_position total_pop_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Na_pos  	  = restore_dvar_matrix_position();
	const dvar_matrix_position Za_pos  	  = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"total_pop_begin");

	dmatrix dftotal_pop = restore_dvar_matrix_derivatives(total_pop_pos);
	dmatrix dfNa 	    = restore_dvar_matrix_derivatives(Na_pos);
	dmatrix dfZa 	    = restore_dvar_matrix_derivatives(Za_pos);

	CParam* param = (CParam*) pos_param;
	PMap* map = (PMap*) pos_map;
	CMatrices* mat = (CMatrices*) pos_mat;

	const int imax = map->imax;
	const int imin = map->imin;

	//We will need bunch of variables to recompute access0 and N(a)
	const int nbf = param->get_nbforage();
	const int nbl = param->nb_layer;

	ivector day_layer(0,nbf-1); day_layer = param->day_layer;
	ivector night_layer(0,nbf-1); night_layer = param->night_layer;

	dmatrix uu(map->imin1,map->imax1,map->jinf1,map->jsup1);
	uu = mat->density_before(sp,t_count,age);

        int age_adult_habitat = age;
        for (int aa=age; aa>1; aa--){
                if (param->age_compute_habitat(sp,aa-1)!=param->age_compute_habitat(sp,aa)){
                        age_adult_habitat = aa;
                        break;
                }
        }

	double oxy_teta = param->a_oxy_habitat[sp];
	double oxy_cr   = param->b_oxy_habitat[sp];
	double temp_age = param->temp_age[sp][age_adult_habitat];
	double sigma_ha = param->sigma_ha[sp][age_adult_habitat];
	const double twosigsq = 2.0*sigma_ha*sigma_ha;
	const double temp_max = param->b_sst_spawning(sp);
	const double delta1   = param->thermal_func_delta[0][sp];
	const double delta2   = param->thermal_func_delta[1][sp];
	const double delta3   = param->thermal_func_delta[2][sp];

	const int Tfunc_Gaussian = param->gaussian_thermal_function[sp];

	d3_array tempn,oxygen,forage;
	tempn.allocate(0,nbl-1);
	oxygen.allocate(0,nbl-1);
	for (int k=0; k<nbl;k++){
		tempn(k).allocate(imin, imax, map->jinf, map->jsup);
		oxygen(k).allocate(imin, imax, map->jinf, map->jsup);
		tempn(k) = mat->tempn(t_count,k);
		oxygen(k) = mat->oxygen(t_count,k);
	}

	forage.allocate(0,nbf-1);
	for (int f=0; f<nbf; f++){
		forage(f).allocate(imin, imax, map->jinf, map->jsup);
		forage(f) = mat->forage(t_count,f);
	}
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nlayer = map->carte(i,j);
			if (nlayer){

				//recompute access0
				const double DL = mat->daylength(jday,j)/24.0;
				dvector lf_access(0,nbl-1);
				lf_access.initialize();
				double sumL = 0.0;
				double l_access = 0.0;
				//for (int l=0; l<nlayer; l++){
				for (int l=0; l<nbl; l++){
					if (Tfunc_Gaussian){
						l_access  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),twosigsq,temp_age,oxy_teta,oxy_cr);
					} else {
						l_access  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),temp_age,temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);
					}
					// weighted by the Forage biomass
					for (int n=0;n<nbf;n++){
						if (day_layer[n]==l)   lf_access(l)+= l_access* forage(n,i,j)*DL; 
						if (night_layer[n]==l) lf_access(l)+= l_access* forage(n,i,j)*(1-DL);
					}
					sumL += lf_access(l);			
				}
				double access0 = lf_access(0)/sumL;
				//end of recomputation section	

				double dftotal_pr = 0.0;
		
				//total_pop(i,j) = total_pr(i,j) + access0*mat.dvarPop_species(sp,a,i,j);
				dftotal_pr += dftotal_pop(i,j);
				dfNa(i,j)  += access0*dftotal_pop(i,j);
				dfZa(i,j)  += uu(i,j)*dftotal_pop(i,j);
				dftotal_pop(i,j) = 0.0;

				//double total_pr = total_pop(i,j);
				dftotal_pop(i,j) += dftotal_pr;
			}
		}
	}
	dftotal_pop.save_dmatrix_derivatives(total_pop_pos); 
	dfNa.save_dmatrix_derivatives(Na_pos);
	dfZa.save_dmatrix_derivatives(Za_pos);
}

void dv_spawning_biomass_comp()
{
	verify_identifier_string2((char*)"spawning_biomass_end");
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position total_pop_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Na_pos  	  = restore_dvar_matrix_position();
	const double maturity_age = restore_double_value();
	verify_identifier_string2((char*)"spawning_biomass_begin");

	dmatrix dftotal_pop = restore_dvar_matrix_derivatives(total_pop_pos);
	dmatrix dfNa 	    = restore_dvar_matrix_derivatives(Na_pos);
	PMap* map = (PMap*) pos_map;

	const int imax = map->imax;
	const int imin = map->imin;
	
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){
				double dftotal_pr = 0.0;
		
				//total_pop(i,j) = total_pr(i,j) + maturity_age[a]* mat.dvarPop_species(sp,a,i,j);
				dftotal_pr += dftotal_pop(i,j);
				dfNa(i,j)  += maturity_age * dftotal_pop(i,j);
				dftotal_pop(i,j) = 0.0;

				//double total_pr = total_pop(i,j);
				dftotal_pop(i,j) += dftotal_pr;
			}
		}
	}

	dftotal_pop.save_dmatrix_derivatives(total_pop_pos); 
	dfNa.save_dmatrix_derivatives(Na_pos);
}


