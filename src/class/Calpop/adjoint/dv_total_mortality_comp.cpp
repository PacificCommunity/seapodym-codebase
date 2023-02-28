#include "calpop.h"

//Main function with memory control and adjoint functions for: 
///computing the sum of natural and fishing mortalities
///Forward functions are in total_mortality_comp.cpp


void dv_total_mort_comp(void);
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);

void CCalpop::Precalrec_total_mortality_comp(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, dvar_matrix& mortality, const int age, const int sp, const int t_count, const int year, const int month, const int step_count)
{
	int k = 0;
	const int nb_fishery = param.get_nbfishery();
	dvar_matrix Q, Sslope;
	Q.allocate(map.imin, map.imax, map.jinf, map.jsup); Q.initialize();
	Sslope.allocate(map.imin, map.imax, map.jinf, map.jsup); Sslope.initialize();

	dmatrix effort;
	effort.allocate(map.imin, map.imax, map.jinf, map.jsup);
	effort.initialize();

	for (int f=0; f<nb_fishery; f++){ 
		if (param.mask_fishery_sp[sp][f]){
			//2014: no use of effort for fisheries with no or bad effort data:
			if (param.mask_fishery_sp_no_effort[sp][f]){ k++; continue;}

			double catchability = param.q_sp_fishery[sp][k]*(1.0+param.q_dyn_fishery[f]*step_count);
			if (catchability<0) catchability = 0; // in case if negative trend has been chosen
			//double selectivity = param.selectivity_comp(sp,age,f,k);
			//double selectivity = param.selectivity_comp(sp,age,f,k,step_count);
			const double selectivity = Selectivity(sp,f,age);
			double sq  = selectivity * catchability;

			Q = param.dvarsQ_sp_fishery[sp][k];
			Sslope = param.dvarsSslope_sp_fishery[sp][k];
			mat.dvarsU.initialize();
			mat.dvarsV.initialize();
			if (param.s_func_type[f]>1){
				//Slength:
				mat.dvarsU = param.dvarsSlength_sp_fishery[sp][k];
				if (param.s_func_type[f]==3)
					//Sasympt:
					mat.dvarsV = param.dvarsSasympt_sp_fishery[sp][k];
			}
			if (!param.fdata_rm){			
				//effort was redistributed here using param.afcoef routine:
				rw.get_effort_rm(param, effort, f, year, month);
			} else {
				if (param.mpa_simulation)
					effort = mat.effort[f];
				else 
					//to use GMB effort on model resolution resolution (old configuration)
					rw.get_effort(param, effort, f, year, month);
			}
			
			mortality.save_dvar_matrix_position();
			save_identifier_string2((char*)"before_total_mort_comp");
			
			
			//precalrec_total_mortality_comp(map, mat, mortality, sq, f);
			precalrec_total_mortality_comp(map.carte, effort, mortality, sq, mat.lat_correction);

			save_identifier_string2((char*)"total_mortality_comp_begin");
			Q.save_dvar_matrix_position();
			Sslope.save_dvar_matrix_position();
			mat.dvarsU.save_dvar_matrix_position();
			mat.dvarsV.save_dvar_matrix_position();
			mortality.save_dvar_matrix_position();
			//save_int_value(t_count);
			save_int_value(step_count);
			save_int_value(year);
			save_int_value(month);
			save_int_value(sp);
			save_int_value(age);
			save_int_value(f);
			save_int_value(k);
			unsigned long int pop    = (unsigned long int)this;
			save_long_int_value(pop);
			//unsigned long int cmat = (unsigned long int)&mat;
			//save_long_int_value(cmat);
			unsigned long int pmap = (unsigned long int)&map;
			save_long_int_value(pmap);
			unsigned long int cparam = (unsigned long int)&param;
			save_long_int_value(cparam);
			unsigned long int crw = (unsigned long int)&rw;
			save_long_int_value(crw);
			save_identifier_string2((char*)"total_mortality_comp_end");
			
			k++;
			gradient_structure::GRAD_STACK1->set_gradient_stack(dv_total_mort_comp);	
		}
	}
}

void dv_total_mort_comp(void)
{
	verify_identifier_string2((char*)"total_mortality_comp_end");
	unsigned long int pos_rw    = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();	
	unsigned long int pos_map   = restore_long_int_value();
	//unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_pop   = restore_long_int_value();
	unsigned k     = restore_int_value();
	unsigned f     = restore_int_value();
	unsigned age   = restore_int_value();
	unsigned sp    = restore_int_value();
	unsigned month = restore_int_value();
	unsigned year  = restore_int_value();
	unsigned step  = restore_int_value();
	//unsigned t_count  = restore_int_value();
	const dvar_matrix_position mort_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position s_asy_pos = restore_dvar_matrix_position();
	const dvar_matrix_position s_len_pos = restore_dvar_matrix_position();
	const dvar_matrix_position sslope_pos= restore_dvar_matrix_position();
	const dvar_matrix_position q_pos     = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"total_mortality_comp_begin");

	verify_identifier_string2((char*)"before_total_mort_comp");
	const dvar_matrix_position m_pr_pos = restore_dvar_matrix_position();

	CCalpop* pop = (CCalpop*) pos_pop;
	//CMatrices* mat = (CMatrices*) pos_mat;
	PMap* map = (PMap*) pos_map;
	CParam* param = (CParam*) pos_param;
	CReadWrite* rw= (CReadWrite*) pos_rw;

	const int s_func_type = param->s_func_type[f];

	dmatrix dfss = restore_dvar_matrix_derivatives(sslope_pos);
	dmatrix dfsl = restore_dvar_matrix_derivatives(s_len_pos);
	dmatrix dfsa = restore_dvar_matrix_derivatives(s_asy_pos);
	dmatrix dfq  = restore_dvar_matrix_derivatives(q_pos);
	dmatrix dfM  = restore_dvar_matrix_derivatives(mort_pos);
	dmatrix dfM_pr  = restore_dvar_matrix_derivatives(m_pr_pos);

	const int imax = map->imax;
	const int imin = map->imin;
	const int jmax = map->jmax;
	const int jmin = map->jmin;

	dmatrix effort(imin,imax,map->jinf,map->jsup); effort.initialize();
	if (!param->fdata_rm)
		//effort was redistributed here using param.afcoef routine:
		rw->get_effort_rm(*param, effort, f, year, month);
	else 
		//to use GMB effort on model resolution resolution (old configuration)
		rw->get_effort(*param, effort, f, year, month);

//	if (param->use_efr_rm)
//		rw->get_effort_rm(*param, effort, f, year, month);
//	else
//rw->get_effort(*param, effort, f, year, month);
//		rw->get_effort_rm(*param, effort, f, year, month);
	const double q_dyn = (1.0+param->q_dyn_fishery[f]*step);
	double catchability = param->q_sp_fishery[sp][k]*q_dyn;

	//const double selectivity = param->selectivity_comp(sp,age,f,k);
	//const double selectivity = param->selectivity_comp(sp,age,f,k,step);
	const double selectivity = pop->Selectivity(sp,f,age);
	double dflength, dfslope, dfasympt;
	param->dfselectivity(dfslope,dflength,dfasympt,sp,age,f,k);

	dvector lat_correction(jmin,jmax); 
	lat_correction.initialize();

	for (int j=jmin;j<=jmax;j++){
		double lastlat = param->lastlat(j);
		lat_correction[j] = param->correction_lat(lastlat);
	}
//int use_zeu = param->use_vld;
	for (int j = jmax; j >= jmin; j--){
		const int imin = map->iinf[j];
		const int imax = map->isup[j];
		const double sl = selectivity*lat_correction(j);
		const double ql = catchability*lat_correction(j);
		for (int i = imax; i >= imin; i--){
			if (map->carte(i,j) && effort(i,j)){
/*
				//test ZEU
				if (use_zeu)
					effort(i,j) *= 1.0-1000.0*pow(mat->vld(t_count,i,j),4);
*/
				
				//mortality(i,j) += mat.effort(f,i,j) * selectivity * catchability * lat_correction;
				dfq(i,j)   += effort(i,j)*q_dyn*sl*dfM(i,j);
				dfss(i,j)  += effort(i,j)*dfslope*ql*dfM(i,j);
				if (s_func_type>1)
					dfsl(i,j)  += effort(i,j)*dflength*ql*dfM(i,j);
				if (s_func_type==3)
					dfsa(i,j)  += effort(i,j)*dfasympt*ql*dfM(i,j);
				dfM_pr(i,j)+= dfM(i,j);
				dfM(i,j)    = 0.0;		
			}
		}
	}

	dfM.save_dmatrix_derivatives(mort_pos);
	dfM_pr.save_dmatrix_derivatives(m_pr_pos);
	dfq.save_dmatrix_derivatives(q_pos);
	dfss.save_dmatrix_derivatives(sslope_pos);
	dfsl.save_dmatrix_derivatives(s_len_pos);
	dfsa.save_dmatrix_derivatives(s_asy_pos);
}

