#include "calpop.h"

///Main function with memory control and adjoint functions for: 
///predicting catch by fishery based on fishing effort. 
///Forward function is in predicted_catch.cpp

void dv_predicted_catch_fishery();
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);
const double sq_degree   = 12347.65; // sq.km in sq.degree (without lat_correction)


void CCalpop::Predicted_Catch_Fishery(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int sp, const int f, const int k, const int year, const int month, const int t_count, const int step_count)
{ 
	mat.dvarCatch_est[sp][k].initialize();
	mat.catch_est[sp][k].initialize();

	const int a0 = param.sp_a0_adult[sp];
	const int nb_ages = param.sp_nb_cohorts[sp];
	const int nb_region = param.nb_region_sp_B[sp];

	//Acoef (matrix of coef to compute densities on fishing area):
	//dimension of matrix
//	const int reso = param.fishery_reso(f);
//	m = reso*60.0/param.deltaX+2; 
//	n = reso*60.0/param.deltaY+2;  
	//catchability and selectivity parameters enter the adjoint as scalars: their derivatives are
	//accumulated in dv_predicted_catch_fishery (no full-grid matrices filled with their values)
	//fishing areas of the effort cells on the grid, once for all age classes
	fishing_area_comp(map, param, mat.effort(f), mat.efflon(f), mat.efflat(f), f, fa);
	for (int age=a0; age<nb_ages; age++){

		mat.dvarCatch_est[sp][k].save_dvar_matrix_position();
		save_identifier_string2((char*)"before_predicted_catch_comp");
		
		predicted_catch_fishery_comp(map, param, mat, f, k, sp, age, value(mat.dvarDensity(sp,age)),step_count);
		for (int r=0; r<nb_region; r++)
			mat.dvarLF_est[sp][age][k].elem_value(r) = mat.C_N_sp_age_fishery[sp][age][k][r];

		save_identifier_string2((char*)"Predicted_Catch_begin");
		mat.dvarCatch_est(sp,k).save_dvar_matrix_position();
		mat.dvarLF_est(sp,age,k).save_dvar_vector_position();
		mat.dvarDensity(sp,age).save_dvar_matrix_position();
		param.dvarsQslope_fishery[k].save_prevariable_position();
		param.dvarsQ_sp_fishery[sp][k].save_prevariable_position();
		param.dvarsSslope_sp_fishery[sp][k].save_prevariable_position();
		param.dvarsSlength_sp_fishery[sp][k].save_prevariable_position();
		param.dvarsSasympt_sp_fishery[sp][k].save_prevariable_position();
		save_int_value(year);
		save_int_value(month);
		save_int_value(t_count);
		save_int_value(step_count);
		save_int_value(sp);
		save_int_value(age);
		save_int_value(f);
		save_int_value(k);
		save_int_value(nb_region);
		unsigned long int pop    = (unsigned long int)this;
		save_long_int_value(pop);
		unsigned long int cmat = (unsigned long int)&mat;
		save_long_int_value(cmat);
		unsigned long int pmap = (unsigned long int)&map;
		save_long_int_value(pmap);
		unsigned long int cparam = (unsigned long int)&param;
		save_long_int_value(cparam);
		unsigned long int crw = (unsigned long int)&rw;
		save_long_int_value(crw);
		save_identifier_string2((char*)"Predicted_Catch_end");

		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_predicted_catch_fishery);
	}
} 

void dv_predicted_catch_fishery()
{
	verify_identifier_string2((char*)"Predicted_Catch_end");
	unsigned long int pos_rw    = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_pop   = restore_long_int_value();
	unsigned nb_region = restore_int_value();
	unsigned k     = restore_int_value();
	unsigned f     = restore_int_value();
	unsigned age   = restore_int_value();
	unsigned sp    = restore_int_value();
	unsigned nstep = restore_int_value();
	unsigned t_count = restore_int_value();
	unsigned month = restore_int_value();
	unsigned year  = restore_int_value();
	const prevariable_position sasympt_pos = restore_prevariable_position();
	const prevariable_position slength_pos = restore_prevariable_position();
	const prevariable_position sslope_pos  = restore_prevariable_position();
	const prevariable_position q_pos       = restore_prevariable_position();
	const prevariable_position qslope_pos  = restore_prevariable_position();
	const dvar_matrix_position uu_pos    = restore_dvar_matrix_position();
	const dvar_vector_position lf_pos    = restore_dvar_vector_position();
	const dvar_matrix_position catch_pos = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"Predicted_Catch_begin");

	verify_identifier_string2((char*)"before_predicted_catch_comp");
	const dvar_matrix_position c_pr_pos  = restore_dvar_matrix_position();

	CCalpop* pop = (CCalpop*) pos_pop;
	CMatrices* mat = (CMatrices*) pos_mat;
	PMap* map = (PMap*) pos_map;
	CParam* param = (CParam*) pos_param;
	CReadWrite* rw= (CReadWrite*) pos_rw;

	const int s_func_type = param->s_func_type[f];

	const int imax = map->imax;
	const int imin = map->imin;
	const int jmax = map->jmax;
	const int jmin = map->jmin;

	dvector lat_correction(jmin,jmax); 
	lat_correction.initialize();

	//fishing areas of the effort cells on the grid: they depend only on the fishery and the month,
	//so they are computed once for the age classes of a fishery-step, which the reverse sweep
	//visits in a row (same values as in the forward code, see fishing_area_comp)
	static std::vector<fishing_area> fa;
	static int fa_f = -1, fa_year = -1, fa_month = -1;
	if ((int)f != fa_f || (int)year != fa_year || (int)month != fa_month){
		dmatrix effort(imin,imax,map->jinf,map->jsup);
		dmatrix efflon(imin,imax,map->jinf,map->jsup);
		dmatrix efflat(imin,imax,map->jinf,map->jsup);
		rw->get_effort_lonlat(*param, effort, efflon, efflat, f, year, month);
		pop->fishing_area_comp(*map, *param, effort, efflon, efflat, f, fa);
		fa_f = f; fa_year = year; fa_month = month;
	}

	dmatrix uu(map->imin1,map->imax1,map->jinf1,map->jsup1);
	uu = mat->density_before(sp,t_count,age);

	for (int j=jmin;j<=jmax;j++){
		double lastlat = param->lastlat(j);
		lat_correction[j] = param->correction_lat(lastlat);
	}


	dmatrix dfuu	   = restore_dvar_matrix_derivatives(uu_pos);
	dmatrix dfC_pred   = restore_dvar_matrix_derivatives(catch_pos);
	dmatrix dfC_pr_age = restore_dvar_matrix_derivatives(c_pr_pos);
	//derivatives of the catchability and selectivity parameters, accumulated over the effort cells
	double dfQ       = 0.0;
	double dfQslope  = 0.0;
	double dfSslope  = 0.0;
	double dfSlength = 0.0;
	double dfSasympt = 0.0;
	dvector dfLF	   = restore_dvar_vector_derivatives(lf_pos);

	const int catch_units = param->fishery_catch_units[f];
	const double weight = pow(1e-3*param->weight[sp][age],catch_units);	
	const double w_area  = weight*sq_degree;

	double catchability = param->q_sp_fishery[sp][k] + param->q_slope_fishery[k] * nstep;
	if (catchability<0) catchability = 0; // in case if negative trend has been chosen

	//const double selectivity = param->selectivity_comp(sp,age,f,k);
	//const double selectivity = param->selectivity_comp(sp,age,f,k,nstep);
	const double selectivity = pop->Selectivity(sp,f,age);
	const double s_c = selectivity*catchability;


	double dflength,dfslope,dfasympt;
	param->dfselectivity(dfslope,dflength,dfasympt,sp,age,f,k);

	//effort cells in the reverse of the forward order (i, then j, descending)
	for (int nc = (int)fa.size()-1; nc >= 0; nc--){
		const fishing_area& e = fa[nc];
		const int i = e.i; const int j = e.j;
		const double eff = e.effort;
		const double afr = e.afr;
		//precompute biomass of fish aggregated to the resolution of fishing data
		double fish  = 0;
		double dffish = 0;
		for (size_t p = 0; p < e.af.size(); p++)
			fish += e.af[p]*uu(e.ii[p],e.jj[p])/(afr*lat_correction[e.jj[p]]);

		//const double F_pr = s_c*eff/(reso*reso);
		//const double F = reso*reso*param->smin1(F_pr);
		const double F_pr = s_c*eff/afr;
		const double F = afr*param->smin1_c(F_pr);
		//end of recomputation section

		double dfF = 0;
		if (param->nb_region){
			//regions containing the cell, in the reverse of the forward order
			const std::vector<int>& regs = map->regions_of_cell(i,j);
			for (int n=(int)regs.size()-1; n>=0; n--){
				const int a = regs[n];
				if (a >= (int)nb_region) continue;
				int r = param->area_sp_B[sp][a]-1;
				//mat.C_N_sp_age_fishery[sp][age][k][r] += F * fish * sq_degree;
				dffish       += F*sq_degree*dfLF(r);
				dfF          += fish*sq_degree*dfLF(r);
			}
		} else if (nb_region){
			//mat.C_N_sp_age_fishery[sp][age][k][0] += F * fish * sq_degree;
			dffish       += F*sq_degree*dfLF(0);
			dfF          += fish*sq_degree*dfLF(0);
		}

		//mat.dvarCatch_est[sp][k][i][j] = C_pr_age + F*fish*w_area;
		dffish    	+= F*w_area*dfC_pred(i,j);
		dfF             += fish*w_area*dfC_pred(i,j);
		dfC_pr_age(i,j)	+= dfC_pred(i,j);
		dfC_pred(i,j)    = 0.0;	

		//double dfs_c = eff*param->dfsmin1(F_pr,reso*reso*dfF)/(reso*reso);
		double dfs_c = eff*param->dfsmin1_c(F_pr, afr*dfF)/afr;
		dfF = 0.0;

		/*
		dfF = 0.0;
		*/

		//const double s_c = selectivity*catchability;
		dfQslope  += nstep*selectivity*dfs_c;
		dfQ       += selectivity*dfs_c;
		dfSslope  += dfslope*catchability*dfs_c;
		if (s_func_type>1)
			dfSlength += dflength*catchability*dfs_c;
		if (s_func_type==3)
			dfSasympt += dfasympt*catchability*dfs_c;
		dfs_c = 0.0;

		for (int p = (int)e.af.size()-1; p >= 0; p--)
			//fish += af(ii,jj)*uu(ii+ki,jj+kj)/lat_correction[jj+kj];
			dfuu(e.ii[p],e.jj[p]) += e.af[p]*dffish/(afr*lat_correction(e.jj[p]));

		dffish = 0;
	}
	if (catchability == 0){
		dfQslope = 0.0;
		dfQ = 0.0;
	}

	dfuu.save_dmatrix_derivatives(uu_pos); 
	dfC_pred.save_dmatrix_derivatives(catch_pos);
	dfC_pr_age.save_dmatrix_derivatives(c_pr_pos);
	save_double_derivative(dfQ, q_pos);
	save_double_derivative(dfQslope, qslope_pos);
	save_double_derivative(dfSslope, sslope_pos);
	save_double_derivative(dfSlength, slength_pos);
	save_double_derivative(dfSasympt, sasympt_pos);
	dfLF.save_dvector_derivatives(lf_pos);
}


