#include "SeapodymCoupled.h"

///Main function with memory control and adjoint functions for: 
///complex non-linear function that depends on the density of adult biomass
///to derive the food requirement and compute the index of depletion in case 
///if the available micronekton does not support the daily ration. 
///This function might be useful in multi-species models.
///Forward function is in food_requirement_index.cpp


void dv_FR_pop_comp();
void dv_F_access_sum_age_comp();
void dv_ISR_denom_comp();
void dv_food_requirement_index_comp();
double f_accessibility_layer(const double O2, const double T,double twosigsq, double temp_mean, double oxy_teta, double oxy_cr);
double f_accessibility_layer(const double O2, const double T, double temp_age, double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr);


int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);

void SeapodymCoupled::Food_Requirement_Index(dvar_matrix& IFR, dvar_matrix FR_pop, dvar_matrix ISR_denom, const int sp, const int age, const int t_count, const int jday)
{
	IFR.initialize();

	IFR_age_comp(IFR, value(FR_pop), value(ISR_denom), age, sp, t_count);

	save_identifier_string2((char*)"food_requirement_comp_begin");	
	mat.dvarDensity(sp,age).save_dvar_matrix_position();
	for (int n=0; n<nb_forage; n++)
		mat.dvarF_access[n][age].save_dvar_matrix_position();
	FR_pop.save_dvar_matrix_position();
	IFR.save_dvar_matrix_position();
	save_int_value(sp);
	save_int_value(age);
	save_int_value(t_count);
	save_int_value(jday);
	unsigned long int pmap = (unsigned long int)&map;
	save_long_int_value(pmap);
	unsigned long int cparam = (unsigned long int)*&param;
	save_long_int_value(cparam);
	unsigned long int cmat = (unsigned long int)&mat;
	save_long_int_value(cmat);
	save_identifier_string2((char*)"food_requirement_comp_end");

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_food_requirement_index_comp);
}

//Total food requirement function, identical to Spawning_Biomass_comp, except for multiplier
void SeapodymCoupled::FR_pop_comp(dvar_matrix& FR_pop, const int sp)
{ 	
	FR_pop.initialize();
	const double R = param->forage_ration[sp];

	for (int a=param->sp_a0_adult[sp]; a<param->sp_nb_cohorts[sp]; a++){

		const double W = param->weight[sp][a] * 0.001; //in mt
		const double wrt = W * R * param->deltaT;

		for (int i = map.imin; i <= map.imax; i++){	
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++){
				if (map.carte(i,j)){
					double FR_pr = value(FR_pop(i,j));				
					FR_pop.elem_value(i,j) = FR_pr + wrt * value(mat.dvarDensity(sp,a,i,j));
				}				 
			}
		}
		
		save_identifier_string2((char*)"FR_pop_comp_begin");
		save_double_value(wrt);
		mat.dvarDensity(sp,a).save_dvar_matrix_position();
		FR_pop.save_dvar_matrix_position();
		unsigned long int pmap = (unsigned long int)&map;
		save_long_int_value(pmap);
		save_identifier_string2((char*)"FR_pop_comp_end");

		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_FR_pop_comp);
	}
} 


void SeapodymCoupled::ISR_denom_comp(dvar_matrix& ISR_denom, const int sp, const int t_count)
{ 	
	dvar_matrix dvarF_access_sum_age(map.imin, map.imax, map.jinf, map.jsup);

	ISR_denom = 1.0;
        for (int n=0; n<nb_forage; n++){

		mat.F_access_sum_age(sp,n,t_count).initialize();
		dvarF_access_sum_age.initialize();
 		for (int a=param->sp_a0_adult[sp]; a < param->sp_nb_cohorts[sp]; a++){

			for (int i = map.imin; i <= map.imax; i++){	
				const int jmin = map.jinf[i];
				const int jmax = map.jsup[i];
				for (int j = jmin; j <= jmax; j++){
					if (map.carte(i,j)){
						mat.F_access_sum_age(sp,n,t_count,i,j) += value(mat.dvarF_access(n,a,i,j));
						dvarF_access_sum_age.elem_value(i,j) = mat.F_access_sum_age(sp,n,t_count,i,j);
					}				 
				}
			}

	                save_identifier_string2((char*)"F_access_sum_age_comp_begin");
	                dvarF_access_sum_age.save_dvar_matrix_position();
	                mat.dvarF_access(n,a).save_dvar_matrix_position();
        	        unsigned long int pmap = (unsigned long int)&map;
	                save_long_int_value(pmap);
	                save_identifier_string2((char*)"F_access_sum_age_comp_end");

                	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_F_access_sum_age_comp);
		}

		for (int i = map.imin; i <= map.imax; i++){	
			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++){
				if (map.carte(i,j)){
					double ISR_denom_pr = value(ISR_denom(i,j));
					ISR_denom.elem_value(i,j) = ISR_denom_pr * mat.F_access_sum_age(sp,n,t_count,i,j);
				}				 
			}
		}
                save_identifier_string2((char*)"ISR_denom_comp_begin");
	        dvarF_access_sum_age.save_dvar_matrix_position();
                ISR_denom.save_dvar_matrix_position();
                save_int_value(sp);
                save_int_value(n);
                save_int_value(t_count);
                unsigned long int cmat = (unsigned long int)&mat;
                save_long_int_value(cmat);
       	        unsigned long int pmap = (unsigned long int)&map;
                save_long_int_value(pmap);
                save_identifier_string2((char*)"ISR_denom_comp_end");

               	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_ISR_denom_comp);

        }
} 

void dv_FR_pop_comp()
{
	verify_identifier_string2((char*)"FR_pop_comp_end");
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position FR_pop_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Na_pos      = restore_dvar_matrix_position();
	const double wrt = restore_double_value();
	verify_identifier_string2((char*)"FR_pop_comp_begin");

	dmatrix dfFR_pop = restore_dvar_matrix_derivatives(FR_pop_pos);
	dmatrix dfNa 	 = restore_dvar_matrix_derivatives(Na_pos);

	PMap* map = (PMap*) pos_map;

	const int imax = map->imax;
	const int imin = map->imin;
	
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){
				double dfFR_pr = 0.0;
		
				//FR_pop(i,j) = FR_pr(i,j) + W * R * mat.dvarPop_species(sp,a,i,j);
				dfFR_pr   += dfFR_pop(i,j);
				dfNa(i,j) += wrt * dfFR_pop(i,j);
				dfFR_pop(i,j) = 0.0;

				//double FR_pr = FR_pop(i,j);
				dfFR_pop(i,j) += dfFR_pr;
			}
		}
	}

	dfFR_pop.save_dmatrix_derivatives(FR_pop_pos); 
	dfNa.save_dmatrix_derivatives(Na_pos);
}


void dv_F_access_sum_age_comp()
{
	verify_identifier_string2((char*)"F_access_sum_age_comp_end");
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position Fa_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position Fasum_pos = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"F_access_sum_age_comp_begin");

	dmatrix dfF_access     = restore_dvar_matrix_derivatives(Fa_pos);
	dmatrix dfF_access_sum = restore_dvar_matrix_derivatives(Fasum_pos);

	PMap* map = (PMap*) pos_map;

	const int imax = map->imax;
	const int imin = map->imin;
	
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){
				double dfF_access_sum_pr = 0.0;
		
				//F_access_sum(i,j) = F_access_sum_pr(i,j) + mat.dvarF_access(n,a,i,j);
				dfF_access_sum_pr   += dfF_access_sum(i,j);
				dfF_access(i,j)     += dfF_access_sum(i,j);
				dfF_access_sum(i,j)  = 0.0;

				//double F_access_sum_pr = F_access_sum(i,j);
				dfF_access_sum(i,j) += dfF_access_sum_pr;
			}
		}
	}

	dfF_access.save_dmatrix_derivatives(Fa_pos); 
	dfF_access_sum.save_dmatrix_derivatives(Fasum_pos);
}

void dv_ISR_denom_comp()
{
	verify_identifier_string2((char*)"ISR_denom_comp_end");
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	const int t_count  = restore_int_value();
	const int n        = restore_int_value();
	const int sp       = restore_int_value();
	const dvar_matrix_position ISR_pos = restore_dvar_matrix_position();
	const dvar_matrix_position Fasum_pos  = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"ISR_denom_comp_begin");

	dmatrix dfISR_denom = restore_dvar_matrix_derivatives(ISR_pos);
	dmatrix dfF_access_sum = restore_dvar_matrix_derivatives(Fasum_pos);

	CMatrices* mat = (CMatrices*) pos_mat;
	PMap* map = (PMap*) pos_map;

	dvector F_access_sum(0,n);
	F_access_sum.initialize();

	const int imax = map->imax;
	const int imin = map->imin;
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){
				double dfISR_denom_pr = 0.0;
	
				//recompute ISR_denom_pr
				double ISR_denom_pr = 1.0;
				for (int f=0; f<n; f++){
					F_access_sum(f) = mat->F_access_sum_age(sp,f,t_count,i,j);
					ISR_denom_pr *= F_access_sum(f);
				}
				F_access_sum(n) = mat->F_access_sum_age(sp,n,t_count,i,j);
		
				//ISR_denom.elem_value(i,j) = ISR_denom_pr * dvarF_access_sum_age(i,j);
				dfISR_denom_pr      += F_access_sum(n) * dfISR_denom(i,j);
				dfF_access_sum(i,j) += ISR_denom_pr * dfISR_denom(i,j);
				dfISR_denom(i,j)     = 0.0;

				//double ISR_denom_pr = value(ISR_denom(i,j));
				dfISR_denom(i,j) += dfISR_denom_pr;

			}
		}
	}

	dfISR_denom.save_dmatrix_derivatives(ISR_pos); 
	dfF_access_sum.save_dmatrix_derivatives(Fasum_pos);
}

void dv_food_requirement_index_comp()
{
	verify_identifier_string2((char*)"food_requirement_comp_end");
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	const int jday    = restore_int_value();
	const int t_count  = restore_int_value();
	const int age      = restore_int_value();
	const int sp       = restore_int_value();
	const dvar_matrix_position IFR_pos       = restore_dvar_matrix_position();	
	const dvar_matrix_position FR_pop_pos    = restore_dvar_matrix_position();	
	const dvar_matrix_position Fpos5  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos4  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos3  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos2  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos1  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos0  = restore_dvar_matrix_position();
	const dvar_matrix_position Na_pos = restore_dvar_matrix_position();	
	verify_identifier_string2((char*)"food_requirement_comp_begin");

	CParam* param = (CParam*) pos_param;
	PMap* map = (PMap*) pos_map;
	CMatrices* mat = (CMatrices*) pos_mat;

	const int imax = map->imax;
	const int imin = map->imin;
	const int nb_forage = param->get_nbforage();
	d3_array dfF_access;
	dfF_access.allocate(0,nb_forage-1);
	for (int n=0; n<nb_forage; n++){
		dfF_access[n].allocate(imin, imax, map->jinf, map->jsup);
		dfF_access[n].initialize();
	}
	dfF_access[0]  = restore_dvar_matrix_derivatives(Fpos0);
	dfF_access[1]  = restore_dvar_matrix_derivatives(Fpos1);
	dfF_access[2]  = restore_dvar_matrix_derivatives(Fpos2);
	dfF_access[3]  = restore_dvar_matrix_derivatives(Fpos3);
	dfF_access[4]  = restore_dvar_matrix_derivatives(Fpos4);
	dfF_access[5]  = restore_dvar_matrix_derivatives(Fpos5);
	dmatrix dfFR_pop    = restore_dvar_matrix_derivatives(FR_pop_pos);
	dmatrix dfIFR       = restore_dvar_matrix_derivatives(IFR_pos);
	dmatrix dfNa        = restore_dvar_matrix_derivatives(Na_pos);

	const double pi = 4.0*(atan(0.5)+atan(1.0/3.0));
	double oxy_teta = param->a_oxy_habitat[sp];
	double oxy_cr = param->b_oxy_habitat[sp];
        int age_adult_habitat = age;
        for (int aa=age; aa>1; aa--){
                if (param->age_compute_habitat(sp,aa-1)!=param->age_compute_habitat(sp,aa)){
                        age_adult_habitat = aa;
                        break;
                }
        }	
	double temp_age = param->temp_age[sp][age_adult_habitat];
	double sigma_ha = param->sigma_ha[sp][age_adult_habitat];
	const double twosigsq = 2.0*sigma_ha*sigma_ha;
	const double W = param->weight[sp][age] * 0.001;
	const double R = param->forage_ration[sp];
	const double wrt = W * R * param->deltaT;
	const double residual_competition = param->residual_competition[sp];
	const double slope = 0.01;
	const double lslope = log(slope);

	const double temp_max = param->b_sst_spawning(sp);
	const double delta1   = param->thermal_func_delta[0][sp];
	const double delta2   = param->thermal_func_delta[1][sp];
	const double delta3   = param->thermal_func_delta[2][sp];

	const int Tfunc_Gaussian = param->gaussian_thermal_function[sp];

	d3_array tempn,oxygen,forage;
	const int nbl = param->nb_layer;
	tempn.allocate(0,nbl-1);
	oxygen.allocate(0,nbl-1);
	for (int k=0; k<nbl;k++){
		tempn(k).allocate(imin, imax, map->jinf, map->jsup);
		oxygen(k).allocate(imin, imax, map->jinf, map->jsup);
		tempn(k) = mat->tempn(t_count,k);
		oxygen(k) = mat->oxygen(t_count,k);
	}

	forage.allocate(0,nb_forage-1);
	for (int f=0; f<nb_forage; f++){
		forage(f).allocate(imin, imax, map->jinf, map->jsup);
		forage(f) = mat->forage(t_count,f);
	}

	ivector day_layer(0,nb_forage-1); day_layer = param->day_layer;
	ivector night_layer(0,nb_forage-1); night_layer = param->night_layer;

	d3_array N;
	const int a0_adult = param->sp_a0_adult[sp]; 
	const int nb_cohorts = param->sp_nb_cohorts[sp];

	dvector wrt_a;
	wrt_a.allocate(a0_adult,nb_cohorts-1); 
	N.allocate(a0_adult,nb_cohorts-1);
	for (int a=a0_adult; a<nb_cohorts; a++){
		wrt_a[a] = param->weight[sp][a] * 0.001 * R * param->deltaT;
		N(a).allocate(map->imin1,map->imax1,map->jinf1,map->jsup1);
		N(a) = mat->density_before(sp,t_count,a);
	}

	dvector theta(0,nb_forage-1);
	dvector dftheta(0,nb_forage-1);
	dvector F_access(0,nb_forage-1);

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nlayer = map->carte(i,j);
			if (nlayer>0){	

				//variables we will need to compute derivatives:
				double F_accessible = 0.0;
				double F_required = 0.0;
				double IFR_age = 0.0;
				double IFR_age_pop = 0.0;
				double FR_pop = 0.0;
				double IFR_pop = 0.0;
				double sumF = 0.0;
				double sumF_access = 0.0;
				theta.initialize();

				//recompute F_access, F_accessible and F_required, sumF and sumF_access
				//then FR_pop, IFR_pop, IFR_age, IFR_age_pop, ISR_denom, ISR_age and ISR_age_pr
				//1. F_access
				F_access.initialize();
				const double DL = mat->daylength(jday,j)/24.0;
				//dvector l_access(0,nlayer-1);
				dvector l_access(0,nbl-1);
				l_access.initialize();
				//for (int l=0; l<nlayer; l++)
				for (int l=0; l<nbl; l++){
					if (Tfunc_Gaussian){
						l_access(l)  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),twosigsq,
									temp_age,oxy_teta,oxy_cr);		
					} else {
						l_access(l)  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),temp_age,
									temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);		
					}
				}
				for (int n=0; n<nb_forage; n++)
					if (day_layer[n] < nlayer)
						F_access(n) = l_access[day_layer[n]]*DL + l_access[night_layer[n]]*(1-DL);

				//2.sum_F_access, sumF, ISR_age_pr
				for (int n=0; n<nb_forage; n++){
					sumF_access += F_access(n);	
					sumF += residual_competition * forage(n,i,j);
//					ISR_age_pr *= F_access(n);
				}
/*
				//3. ISR_denom
				for (int n=0; n<nb_forage; n++){
					ISR_denom *= mat->F_access_sum_age(sp,n,t_count,i,j);
				}
*/
				//4. FR_pop and IFR_pop
				for (int a=a0_adult; a <nb_cohorts; a++){
					FR_pop += wrt_a(a) * N(a,i,j);
				}
				//IFR_pop = sumF/FR_pop;
				IFR_pop = (10/pi)*atan((sumF/FR_pop)*pi/10);

				for (int n=0; n<nb_forage; n++)
					theta(n) = F_access(n)/(sumF_access);
				
				//5. F_accessible and F_required
				F_required = N(age,i,j) * wrt;
				for (int n=0; n<nb_forage; n++){
					F_accessible += forage(n,i,j) * theta(n);
				}
				//6. IRS_age
				//ISR_age = ISR_age_pr/(ISR_denom+1e-4);

				//7. IFR_age
				IFR_age = F_accessible/(F_required+1e-4);

				//8. IFR_age_pop
				//IFR_age_pop = IFR_age * pow(IFR_pop,1.0-ISR_age);
				IFR_age_pop = IFR_age * IFR_pop;
				//////////end of recomputation section//////////


				//beginning of adjoint code
				double dfIFR_pop = 0.0;
				double dfIFR_age = 0.0;
				double dfIFR_age_pop = 0.0;
				double dfF_accessible = 0.0;
				double dfF_required = 0.0;
				double dfsumF_access = 0.0;
				dftheta.initialize();

				//4. Finally scaled IFR, which will be used in Mortality function
				//IFR.elem_value(i,j) = 1.0/(1.0+pow(slope,IFR_age_pop-1.0));
				const double expr01 = pow(slope,IFR_age_pop-1.0);
				const double expr02 = (1.0+expr01)*(1.0+expr01);
				dfIFR_age_pop -=  (expr01*lslope/expr02) * dfIFR(i,j);		
				dfIFR(i,j)     = 0.0;
		
				//3. IFR of cohort taking account competition for resource with others
/*				//IFR_age_pop = IFR_age * pow(IFR_pop,1.0-ISR_age);
				double expr = pow(IFR_pop,1.0-ISR_age);
				dfIFR_age  += expr * dfIFR_age_pop;	
				dfIFR_pop  += IFR_age * (1.0-ISR_age) * pow(IFR_pop,-ISR_age) * dfIFR_age_pop;
				dfISR_age  -= IFR_age * expr * log(IFR_pop) * dfIFR_age_pop; 
				dfIFR_age_pop  = 0.0;
*/	
				//IFR_age_pop = IFR_age * IFR_pop;
				double expr = IFR_pop;
				dfIFR_age  += expr * dfIFR_age_pop;	
				dfIFR_pop  += IFR_age * dfIFR_age_pop;
				dfIFR_age_pop  = 0.0;
/*
				//2. compute IRS
				//ISR_age = ISR_age_pr/(ISR_denom(i,j)+1e-4);
				expr = ISR_denom+1e-4;
				dfISR_denom(i,j) -= dfISR_age * ISR_age_pr/(expr*expr);
				dfISR_age_pr += dfISR_age/expr;
				dfISR_age = 0.0;
*/

				//1. compute IFR of cohort
				//double IFR_age = F_accessible/(F_required+1e-4);
				expr = F_required+1e-4;
				dfF_accessible += dfIFR_age/expr;
				dfF_required   -= dfIFR_age*F_accessible/(expr*expr);
				dfIFR_age = 0.0;

				//F_required = value(mat.dvarDensity[sp][age][i][j]) * W * R;
				dfNa(i,j) +=  wrt * dfF_required;
				dfF_required = 0.0;

				for (int n=nb_forage-1; n>=0; n--){
				//	F_accessible += mat.forage[t][n][i][j] * theta(n);
					dftheta(n) += forage(n,i,j)*dfF_accessible;
				}

				expr = sumF_access;
				for (int n=nb_forage-1; n>=0; n--){
					//theta(n) = value(dvarF_access(n,age,i,j))/(sumF_access);
					dfF_access(n,i,j) += dftheta(n)/expr;
					dfsumF_access     -= dftheta(n)*F_access(n)/(expr*expr);
					dftheta(n) = 0.0;
				}

				//double IFR_pop = sum_F/FR_pop;
				//double IFR_pop = (10/pi)*atan(pi*sum_F/(FR_pop*10);
				double x = sumF/FR_pop;
				double dfx = 1.0/(1.0+pi*pi*x*x/(100));
				dfFR_pop(i,j) -= dfIFR_pop * dfx * sumF/(FR_pop*FR_pop);			
				//dfFR_pop(i,j) -= dfIFR_pop * sumF/(FR_pop*FR_pop);
				dfIFR_pop = 0.0;

				for (int n=nb_forage-1; n>=0; n--){
//					double ISR_age_pr_pr = 1.0;
//					for (int nn=0; nn<n; nn++)
//						ISR_age_pr_pr *= F_access[n];

					//ISR_age_pr *= value(dvarF_access(n,age,i,j));
//					dfF_access(n,i,j) += ISR_age_pr_pr*dfISR_age_pr;

					//sumF_access += value(dvarF_access(n,age,i,j));	
					dfF_access(n,i,j) += dfsumF_access; 
				}
			}
		}
	}
	dfF_access[0].save_dmatrix_derivatives(Fpos0);
	dfF_access[1].save_dmatrix_derivatives(Fpos1);
	dfF_access[2].save_dmatrix_derivatives(Fpos2);
	dfF_access[3].save_dmatrix_derivatives(Fpos3);
	dfF_access[4].save_dmatrix_derivatives(Fpos4);
	dfF_access[5].save_dmatrix_derivatives(Fpos5);
	dfFR_pop.save_dmatrix_derivatives(FR_pop_pos);
	dfIFR.save_dmatrix_derivatives(IFR_pos);
	dfNa.save_dmatrix_derivatives(Na_pos);			

}

