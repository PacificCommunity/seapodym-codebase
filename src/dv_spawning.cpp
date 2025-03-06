#include "SeapodymCoupled.h"
#include <fvar.hpp>

///Main function with memory control and adjoint functions for: 
///computing the new larval biomass in two regimes of the simulation:
///during climatological spinup and after the population has been built.
///Forward functions are in spawning.cpp


double pred_surface_comp(dvector forage, const double DL, const int nb_forage, ivector day_layer, ivector night_layer);
double hs_comp(double SST, double preys, double predators, const double a, const double b, const double c, const double d, const double e, const double ssv);
void dv_spawning_in_hs_comp();
void dv_spawning_adult_func_comp();
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);

void SeapodymCoupled::Spawning(dvar_matrix& J, dvar_matrix& Hs, dvar_matrix& Nmature, const int jday, const int sp, const int t_count)
{
	J.initialize();

	dmatrix J_c = value(J);
	dmatrix Hs_c = value(Hs);
	dmatrix N_mat = value(Nmature);
	dvariable nb_recruitment = param->dvarsNb_recruitment[sp];

	dvar_matrix Nbr(map.imin,map.imax,map.jinf,map.jsup); 
	Nbr = nb_recruitment;
	
	dvariable a_adults_spawning = param->dvarsA_adults_spawning[sp];
	dvar_matrix A_sp(map.imin,map.imax,map.jinf,map.jsup); 
	A_sp = a_adults_spawning;

	double a_sst, b_sst;
	if (!param->uncouple_sst_larvae[sp]){
		a_sst = param->a_sst_spawning[sp];
		b_sst = param->b_sst_spawning[sp];
	} else {
		a_sst = param->a_sst_larvae[sp];
		b_sst = param->b_sst_larvae[sp];
	}

	if (param->elarvae_model[sp] || param->spawning_adult_func_only[sp]){
		spawning_adult_func_comp(J_c,N_mat,value(nb_recruitment),value(a_adults_spawning));
		
		J = nograd_assign(J_c);

		save_identifier_string2((char*)"spawning_adult_func_begin");
		nb_recruitment.save_prevariable_value();
		Nbr.save_dvar_matrix_position();
		a_adults_spawning.save_prevariable_value();
		A_sp.save_dvar_matrix_position();
		Nmature.save_dvar_matrix_value();//TODO: recompute using density_after
		Nmature.save_dvar_matrix_position();
		J.save_dvar_matrix_position();
		unsigned long int pmap   = (unsigned long int)&map;
		save_long_int_value(pmap);
		save_identifier_string2((char*)"spawning_adult_func_end");

		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_spawning_adult_func_comp);

	} else {

		spawning_in_hs_comp(J_c,Hs_c,N_mat,value(nb_recruitment),value(a_adults_spawning));	
		J = nograd_assign(J_c);
		
		save_identifier_string2((char*)"spawning_in_hs_begin");
		nb_recruitment.save_prevariable_value();
		Nbr.save_dvar_matrix_position();
		a_adults_spawning.save_prevariable_value();
		A_sp.save_dvar_matrix_position();
		Hs.save_dvar_matrix_position();
		Nmature.save_dvar_matrix_value();//TODO: recompute using density_after
		Nmature.save_dvar_matrix_position();
		J.save_dvar_matrix_position();
		save_double_value(a_sst);
		save_double_value(b_sst);
		save_double_value(param->alpha_hsp_prey(sp));
		save_double_value(param->alpha_hsp_predator(sp));
		save_double_value(param->beta_hsp_predator(sp));
		save_int_value(jday);
		save_int_value(t_count);	
                unsigned long int cparam = (unsigned long int)*&param;
                save_long_int_value(cparam);
		unsigned long int pmap   = (unsigned long int)&map;
		save_long_int_value(pmap);
		unsigned long int pmat   = (unsigned long int)&mat;
		save_long_int_value(pmat);
		save_identifier_string2((char*)"spawning_in_hs_end");

		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_spawning_in_hs_comp);
	}
} 

void dv_spawning_adult_func_comp()
{
	verify_identifier_string2((char*)"spawning_adult_func_end");
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position J_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Nm_pos = restore_dvar_matrix_position();
	dmatrix Nm = restore_dvar_matrix_value(Nm_pos);
	const dvar_matrix_position Bsp_pos = restore_dvar_matrix_position();
	const double b = restore_prevariable_value();
	const dvar_matrix_position Nbr_pos = restore_dvar_matrix_position();
	const double R = restore_prevariable_value();
	verify_identifier_string2((char*)"spawning_adult_func_begin");

	PMap* map = (PMap*) pos_map;

	dmatrix dfNbr = restore_dvar_matrix_derivatives(Nbr_pos);
	dmatrix dfNm = restore_dvar_matrix_derivatives(Nm_pos);
	dmatrix dfBsp = restore_dvar_matrix_derivatives(Bsp_pos);
	dmatrix dfJ   = restore_dvar_matrix_derivatives(J_pos);

	const int imax = map->imax;
	const int imin = map->imin;

	const double R_nb = 1000.0*R; //units of R are thous. nb.
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){	

				double adult = Nm(i,j);
				double expr0 = 1.0+b*adult;		
				double expr1 = pow(expr0,2.0);		

				//double f_adults = 1000.0*R*adult/(1.0+b*adult);
				dfNbr(i,j) += 1000.0*(adult/expr0) * dfJ(i,j);
				dfNm(i,j)  += (R_nb/expr1) * dfJ(i,j);
				dfBsp(i,j) -= (R_nb*adult*adult/expr1) * dfJ(i,j);
				dfJ(i,j) = 0.0;
			}
		}
	}
	dfNbr.save_dmatrix_derivatives(Nbr_pos); 
	dfNm.save_dmatrix_derivatives(Nm_pos); 
	dfBsp.save_dmatrix_derivatives(Bsp_pos); 
	dfJ.save_dmatrix_derivatives(J_pos);
}

void dv_spawning_in_hs_comp()
{
	verify_identifier_string2((char*)"spawning_in_hs_end");
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned t_count = restore_int_value();
	unsigned jday    = restore_int_value();
	double e_hs = restore_double_value();
	double d_hs = restore_double_value();
	double c_hs = restore_double_value();
	double b_hs = restore_double_value();
	double a_hs = restore_double_value();
	const dvar_matrix_position J_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Nm_pos = restore_dvar_matrix_position();
	dmatrix Nm = restore_dvar_matrix_value(Nm_pos);
	const dvar_matrix_position Hs_pos = restore_dvar_matrix_position();
	const dvar_matrix_position Bsp_pos = restore_dvar_matrix_position();
	const double b = restore_prevariable_value();
	const dvar_matrix_position Nbr_pos = restore_dvar_matrix_position();
	const double R = restore_prevariable_value();
	verify_identifier_string2((char*)"spawning_in_hs_begin");

	CParam* param  = (CParam*) pos_param;
	CMatrices* mat = (CMatrices*) pos_mat;
	PMap* map = (PMap*) pos_map;

	dmatrix dfNbr = restore_dvar_matrix_derivatives(Nbr_pos);
	dmatrix dfNm = restore_dvar_matrix_derivatives(Nm_pos);
	dmatrix dfBsp = restore_dvar_matrix_derivatives(Bsp_pos);
	dmatrix dfHs  = restore_dvar_matrix_derivatives(Hs_pos);
	dmatrix dfJ   = restore_dvar_matrix_derivatives(J_pos);

	double pp_transform = param->pp_transform;
	const unsigned int nbf = param->get_nbforage();
	ivector dlayer(0,nbf-1); dlayer = param->day_layer;
	ivector nlayer(0,nbf-1); nlayer = param->night_layer;

	const int imax = map->imax;
	const int imin = map->imin;

	dmatrix SST(imin, imax, map->jinf, map->jsup);
	SST = mat->sst(t_count);

	dmatrix O2_l2(imin, imax, map->jinf, map->jsup);
        O2_l2 = mat->oxygen(t_count,1);

	dmatrix np1(imin, imax, map->jinf, map->jsup);
	np1 = mat->np1(t_count);
	d3_array forage;
	forage.allocate(0,nbf-1);
	for (unsigned int f=0; f<nbf; f++){
		forage(f).allocate(imin, imax, map->jinf, map->jsup);
		forage(f) = mat->forage(t_count,f);
	}
	dvector F(0,nbf-1);

	const double R_nb = 1000.0*R; //units of R are thous. nb.
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){	

				//recompute Hs
				const double DL = mat->daylength(jday,j)/24.0;
				double preys = np1(i,j) * pp_transform;
				for (unsigned int n=0; n<nbf; n++)
					F(n) = forage(n,i,j);		
				double predators = pred_surface_comp(F,DL,nbf,dlayer,nlayer);	
				double f_oxy = 1.0/(1.0+pow(0.01,O2_l2(i,j)-0.1));
				double Hs = hs_comp(SST(i,j),preys,predators,a_hs,b_hs,c_hs,d_hs,e_hs,1.0) * f_oxy;

				double adult = Nm(i,j);
				double expr0 = 1.0+b*adult;		
				double expr1 = pow(expr0,2.0);		

				//recompute
				double f_adults = R_nb*adult/expr0;
					
				//derivative in case of simple product
				//J(i,j)  = f_adults * Hs(i,j);
				double dffa = Hs * dfJ(i,j);
				dfHs(i,j)  += f_adults * dfJ(i,j);
				dfJ(i,j)    = 0.0;
			
				//double f_adults = 1000.0*R*adult/(1.0+b*adult);
				dfNbr(i,j) += 1000.0*(adult/expr0) * dffa;
				dfNm(i,j)  += (R_nb/expr1) * dffa;
				dfBsp(i,j) -= (R_nb*adult*adult/expr1) * dffa;
			}
		}
	}
	dfNbr.save_dmatrix_derivatives(Nbr_pos); 
	dfNm.save_dmatrix_derivatives(Nm_pos); 
	dfBsp.save_dmatrix_derivatives(Bsp_pos); 
	dfHs.save_dmatrix_derivatives(Hs_pos);
	dfJ.save_dmatrix_derivatives(J_pos);
}

