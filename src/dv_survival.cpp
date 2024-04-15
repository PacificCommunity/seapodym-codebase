#include "SeapodymCoupled.h"
#include <fvar.hpp>

///Main function with memory control, forward and adjoint functions for: 
///ageing term discretisation either for all 0..n_a-1 or for A+ age class.


void dv_age_plus_comp();
void dv_ageing_comp();
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);


void SeapodymCoupled::Survival(dvar_matrix& N_a, dvar_matrix& N_a_1, const int a, const int sp)
{
	int nb_cohorts = param->sp_nb_cohorts[sp];
	PMap* map_ptr=&map;

	if (a<nb_cohorts-1)
		Ageing(N_a, N_a_1, *map_ptr);
	else 
		AgePlus(N_a, N_a_1, *map_ptr);
} 

void SeapodymCoupled::Survival_tagpop(dvar_matrix& N_a, dvar_matrix& N_a_1, const int a, const int sp, const int tagpop)
{
	int nb_cohorts = param->sp_nb_cohorts[sp];
	PMap* map_ptr = &tagmaps[tagpop];

	if (a<nb_cohorts-1)
		Ageing(N_a, N_a_1, *map_ptr);
	else 
		AgePlus(N_a, N_a_1, *map_ptr);
} 


void SeapodymCoupled::Ageing(dvar_matrix& N_a, dvar_matrix& N_a_1, PMap& map)
{
	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){

	 			N_a.elem_value(i,j) = value(N_a_1(i,j));
			}	
		}
	}
	save_identifier_string2((char*)"ageing_begin");
	N_a_1.save_dvar_matrix_position();
	N_a.save_dvar_matrix_position();
	unsigned long int pmap = (unsigned long int)&map;
	save_long_int_value(pmap);
	save_identifier_string2((char*)"ageing_end");
	
	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_ageing_comp);

}

void SeapodymCoupled::AgePlus(dvar_matrix& N_a, dvar_matrix& N_a_1, PMap& map)
{
	N_a.save_dvar_matrix_position();
	save_identifier_string2((char*)"before_ageplus");

	for (int i = map.imin; i <= map.imax; i++){	
		const int jmin = map.jinf[i];
		const int jmax = map.jsup[i];
		for (int j = jmin; j <= jmax; j++){
			if (map.carte(i,j)){

	 			N_a.elem_value(i,j) = value(N_a(i,j)) + value(N_a_1(i,j));
			}	
		}
	}
	save_identifier_string2((char*)"ageplus_begin");
	N_a_1.save_dvar_matrix_position();
	N_a.save_dvar_matrix_position();
	unsigned long int pmap = (unsigned long int)&map;
	save_long_int_value(pmap);
	save_identifier_string2((char*)"ageplus_end");
	
	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_age_plus_comp);

}

void dv_ageing_comp()
{

	verify_identifier_string2((char*)"ageing_end");
	unsigned long int pos_map = restore_long_int_value();
	const dvar_matrix_position n_a_pos = restore_dvar_matrix_position();
	const dvar_matrix_position n_a_1_pos = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"ageing_begin");
	PMap* map = (PMap*) pos_map;

	dmatrix dfNa   = restore_dvar_matrix_derivatives(n_a_pos);
	dmatrix dfNa_1 = restore_dvar_matrix_derivatives(n_a_1_pos);

	const int imax = map->imax;
	const int imin = map->imin;
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){ 

				//N_a(t+1).elem_value(i,j) = N_a_1(t,i,j);
				dfNa_1(i,j) += dfNa(i,j);
				dfNa(i,j)    = 0.0;

			}
		}
	}

	dfNa.save_dmatrix_derivatives(n_a_pos); 
	dfNa_1.save_dmatrix_derivatives(n_a_1_pos);
}

void dv_age_plus_comp()
{
	verify_identifier_string2((char*)"ageplus_end");
	unsigned long int pos_map = restore_long_int_value();
	const dvar_matrix_position n_a_pos   = restore_dvar_matrix_position();
	const dvar_matrix_position n_a_1_pos = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"ageplus_begin");

	verify_identifier_string2((char*)"before_ageplus");
	const dvar_matrix_position n_a_t_pos = restore_dvar_matrix_position();

	PMap* map = (PMap*) pos_map;

	dmatrix dfNa   = restore_dvar_matrix_derivatives(n_a_pos);
	dmatrix dfNa_1 = restore_dvar_matrix_derivatives(n_a_1_pos);
	dmatrix dfNa_t = restore_dvar_matrix_derivatives(n_a_t_pos);

	const int imax = map->imax;
	const int imin = map->imin;
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){ 

				//N(t+1,i,j) = q_a_1 * N_a_1(t,i,j) + (1-q_a) * N_a(t,i,j);	
				//N_a(t+1).elem_value(i,j) = N_a(t,i,j) + N_a_1(t,i,j);
				dfNa_1(i,j) += dfNa(i,j);
				dfNa_t(i,j) += dfNa(i,j);
				dfNa(i,j)    = 0.0;

			}
		}
	}

	dfNa.save_dmatrix_derivatives(n_a_pos); 
	dfNa_t.save_dmatrix_derivatives(n_a_t_pos); 
	dfNa_1.save_dmatrix_derivatives(n_a_1_pos);
}


