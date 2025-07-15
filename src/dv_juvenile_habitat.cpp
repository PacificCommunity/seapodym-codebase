#include "VarSimtunaFunc.h"

///Main function with memory control and adjoint functions for: 
///juvenile habitat functions. This function depends on the temperature
///in the epipelagic layer, with or without cannibalism effect 
///(note, the latter has very low sensitivity usually)
///forward functions are in juvenile_habitat.cpp

void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);
void dv_Hj_comp(void);
void dv_Hj_cannibalism_comp(void);

const double c_zoo = 100.0; //constant in the first tests ECCO, need to make it global or parameter as soon as early life stage data will be available

void VarSimtunaFunc::Juvenile_Habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hj, int sp, const int t_count)
{
	Hj.initialize();

	dvariable a,b;

	//These two parameters will impact only the mortality range of juveniles
	//Extending the tolerance interval for juveniles compared to larve
	if (!param.uncouple_sst_larvae[sp]){
		a = param.dvarsA_sst_spawning[sp]+1.0;
		b = param.dvarsB_sst_spawning[sp];
	}
	else {
		a = param.dvarsA_sst_larvae[sp]+1.0;
		//SKJ condition temporarilly to enable the use of REF-2018 (CJFAS) solution with this version
		if (param.sp_name[0].find("skj")!=0)
			b = param.dvarsB_sst_larvae[sp];
		if (param.sp_name[0].find("skj")==0)
			b = param.dvarsB_sst_spawning[sp];
	}

	dvmatr1 = a;
	dvmatr2 = b;

	Hj_comp(param, mat, map, Hj, value(a), value(b),t_count);
	
	save_identifier_string((char*)"Hj_comp_begin");
	a.save_prevariable_value();
	dvmatr1.save_dvar_matrix_position();
	b.save_prevariable_value();
	dvmatr2.save_dvar_matrix_position();
	Hj.save_dvar_matrix_position();
	unsigned long int pmap   = (unsigned long int)&map;
	save_long_int_value(pmap);
	unsigned long int cmat   = (unsigned long int)&mat;
	save_long_int_value(cmat);
	save_int_value(t_count);
	save_identifier_string((char*)"Hj_comp_end");

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Hj_comp);

}


void VarSimtunaFunc::Juvenile_Habitat_cannibalism(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hj, dvar_matrix& total_pop, int sp, const int t_count)
{
	Hj.initialize();

	dvariable a,b;

	//These two parameters will impact only the mortality range of juveniles
	//Extending the tolerance interval for juveniles compared to larve
	if (!param.uncouple_sst_larvae[sp]){
		a = param.dvarsA_sst_spawning[sp]+1.0;
		b = param.dvarsB_sst_spawning[sp];
	}
	else {
		a = param.dvarsA_sst_larvae[sp]+1.0;
		b = param.dvarsB_sst_larvae[sp];
		
	}
	dvariable c = param.dvarsHp_cannibalism[sp];

	dvmatr1 = a;
	dvmatr2 = b;
	dvmatr3 = c;

	Hj_cannibalism_comp(param, mat, map, Hj, value(total_pop), value(a), value(b), value(c), t_count);
	
	save_identifier_string((char*)"Hj_cannibalism_comp_begin");
	a.save_prevariable_value();
	dvmatr1.save_dvar_matrix_position();
	b.save_prevariable_value();
	dvmatr2.save_dvar_matrix_position();
	c.save_prevariable_value();
	dvmatr3.save_dvar_matrix_position();
	total_pop.save_dvar_matrix_value();
	total_pop.save_dvar_matrix_position();
	Hj.save_dvar_matrix_position();
	unsigned long int pmap   = (unsigned long int)&map;
	save_long_int_value(pmap);
	unsigned long int cmat   = (unsigned long int)&mat;
	save_long_int_value(cmat);
	save_int_value(t_count);
	save_identifier_string((char*)"Hj_cannibalism_comp_end");

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Hj_cannibalism_comp);

}

void dv_Hj_comp(void)
{

	verify_identifier_string((char*)"Hj_comp_end");
	unsigned t_count    = restore_int_value();
	unsigned long int pos_mat = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position Hpos = restore_dvar_matrix_position();
	const dvar_matrix_position bpos = restore_dvar_matrix_position();
	double b = restore_prevariable_value();
	const dvar_matrix_position apos = restore_dvar_matrix_position();
	double a = restore_prevariable_value();
	
	verify_identifier_string((char*)"Hj_comp_begin");

	dmatrix dfH = restore_dvar_matrix_derivatives(Hpos);
	dmatrix dfb = restore_dvar_matrix_derivatives(bpos);
	dmatrix dfa = restore_dvar_matrix_derivatives(apos);
	
	CMatrices* mat = (CMatrices*) pos_mat;
	PMap* map = (PMap*) pos_map;

	const int imax = map->imax;
	const int imin = map->imin;

	dmatrix SST(imin,imax,map->jinf,map->jsup);
	SST = mat->tempn(t_count,0);
	//SST = mat->sst(t_count);

	dmatrix np1(imin,imax,map->jinf,map->jsup);
	np1 = mat->np1(t_count);

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){		
			if (map->carte(i,j)){
				double zoo  = np1(i,j)+1e-10;
				double f_food = zoo*zoo/(c_zoo+zoo*zoo);

				//Hj(i,j) = f_food * exp(-pow((SST(i,j)- b),2) / (2*a*a)); 
				dfa(i,j) += f_food * (pow(SST(i,j)-b,2) / (pow(a,3)))*exp(-pow(SST(i,j)-b,2)/(2*pow(a,2)))*dfH(i,j);
				dfb(i,j) += f_food * ((SST(i,j)-b) / (pow(a,2)))*exp(-pow(SST(i,j)-b,2)/(2*pow(a,2))) * dfH(i,j);
				dfH(i,j)  = 0.0;

			}
		}
	}
	dfa.save_dmatrix_derivatives(apos); 
	dfb.save_dmatrix_derivatives(bpos);
	dfH.save_dmatrix_derivatives(Hpos);
}

void dv_Hj_cannibalism_comp(void)
{

	verify_identifier_string((char*)"Hj_cannibalism_comp_end");
	unsigned t_count    = restore_int_value();
	unsigned long int pos_mat = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position Hpos = restore_dvar_matrix_position();
	const dvar_matrix_position Npos = restore_dvar_matrix_position();
	dmatrix total_pop = restore_dvar_matrix_value(Npos);
	const dvar_matrix_position cpos = restore_dvar_matrix_position();
	double c = restore_prevariable_value();
	const dvar_matrix_position bpos = restore_dvar_matrix_position();
	double b = restore_prevariable_value();
	const dvar_matrix_position apos = restore_dvar_matrix_position();
	double a = restore_prevariable_value();	
	verify_identifier_string((char*)"Hj_cannibalism_comp_begin");

	dmatrix dfH = restore_dvar_matrix_derivatives(Hpos);
	dmatrix dfN = restore_dvar_matrix_derivatives(Npos);
	dmatrix dfc = restore_dvar_matrix_derivatives(cpos);
	dmatrix dfb = restore_dvar_matrix_derivatives(bpos);
	dmatrix dfa = restore_dvar_matrix_derivatives(apos);
	
	CMatrices* mat = (CMatrices*) pos_mat;
	PMap* map = (PMap*) pos_map;

	const int imax = map->imax;
	const int imin = map->imin;

	dmatrix SST(imin,imax,map->jinf,map->jsup);
	SST = mat->tempn(t_count,0);
	//SST = mat->sst(t_count);

	dmatrix np1(imin,imax,map->jinf,map->jsup);
	np1 = mat->np1(t_count);

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){		
			if (map->carte(i,j)){

				double N = total_pop(i,j);
				double f_cannibalism = (1-N/(c+N));
				double zoo  = np1(i,j)+1e-10;

				double f_food = zoo*zoo/(c_zoo+zoo*zoo);

				//Hj(i,j) = f_food * exp(-pow((SST(i,j)- b),2) / (2*a*a)) * (1-total_pop(i,j)/(c+total_pop(i,j))); 
				dfa(i,j) += (pow(SST(i,j)-b,2) / (pow(a,3)))*exp(-pow(SST(i,j)-b,2)/(2*pow(a,2)))*f_food*f_cannibalism*dfH(i,j);
				dfb(i,j) += ((SST(i,j)-b) / (pow(a,2)))*exp(-pow(SST(i,j)-b,2)/(2*pow(a,2))) * f_food*f_cannibalism * dfH(i,j);
				dfN(i,j) -= (exp(-pow((SST(i,j)- b),2) / (2*a*a))) * f_food * c/pow(c+N,2) * dfH(i,j);
				dfc(i,j) += (exp(-pow((SST(i,j)- b),2) / (2*a*a))) * f_food * N/pow(c+N,2) * dfH(i,j);
				dfH(i,j)  = 0.0;

			}
		}
	}
	dfa.save_dmatrix_derivatives(apos); 
	dfb.save_dmatrix_derivatives(bpos);
	dfc.save_dmatrix_derivatives(cpos);
	dfN.save_dmatrix_derivatives(Npos);
	dfH.save_dmatrix_derivatives(Hpos);
}

