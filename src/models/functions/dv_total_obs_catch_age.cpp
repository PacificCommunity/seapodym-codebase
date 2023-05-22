#include "calpop.h"

///Main function with memory control and adjoint functions for: 
///computing the total (sum over fisheries without effort) observed catch at age. 
///Forward functions are in total_obs_catch_age.cpp


void dv_total_obs_catch_age_comp(void);
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);


void CCalpop::Total_obs_catch_age_comp(const PMap& map, VarParamCoupled& param, VarMatrices& mat, CReadWrite& rw, const int age, const int sp, const int year, const int month, const int t_count)
{
	mat.dvarCtot_age_obs(sp,age).initialize();
	const int nb_fishery = param.get_nbfishery();
	const int nb_region  = param.nb_region;
	const double area = 1.852*param.deltaX*1.852*param.deltaY;
	const double W_age_mt = param.weight[sp][age] * 0.001; 

	dvar_matrix Sslope, Slen, Sasympt;
	Sslope.allocate(map.imin, map.imax, map.jinf, map.jsup); Sslope.initialize();
	Slen.allocate(map.imin, map.imax, map.jinf, map.jsup); Slen.initialize();
	Sasympt.allocate(map.imin, map.imax, map.jinf, map.jsup);

	int k = 0; 
	int fne = 0; //fisheries no effort counter
	for (int f=0; f<nb_fishery; f++){ 
		if (param.mask_fishery_sp[sp][f]){
			if (param.mask_fishery_sp_no_effort[sp][f]){

				//catch_units=0 if catch is in nb; =1 if in metric tons
				const int catch_units = param.fishery_catch_units[f];
				const double C2Dunits = 1.0/(pow(W_age_mt,catch_units) * area);
				
				Sslope = param.dvarsSslope_sp_fishery[sp][k];
				if (param.s_func_type[f]>1){
					//Slength:
					Slen = param.dvarsSlength_sp_fishery[sp][k];
					if (param.s_func_type[f]==3){
						//Sasympt:
						Sasympt = param.dvarsSasympt_sp_fishery[sp][k];
					}
				}
					
				total_obs_catch_age_comp(map, param, mat, value(mat.dvarDensity(sp,age)),
					mat.catch_obs[sp][k],mat.dvarCtot_age_obs(sp,age), f, fne, k, age, sp, C2Dunits);

				for (int r=0; r<nb_region; r++)
					mat.dvarLF_est[sp][age][k].elem_value(r) = mat.C_N_sp_age_fishery[sp][age][k][r];

				save_identifier_string2((char*)"total_obs_catch_comp_begin");
				Sslope.save_dvar_matrix_position();
				Slen.save_dvar_matrix_position();
				Sasympt.save_dvar_matrix_position();
				dvarsSNsum(fne).save_dvar_matrix_position();
				mat.dvarCtot_age_obs(sp,age).save_dvar_matrix_position();
				mat.dvarLF_est(sp,age,k).save_dvar_vector_position();
				mat.dvarDensity(sp,age).save_dvar_matrix_position();
				save_int_value(t_count);
				save_int_value(year);
				save_int_value(month);
				save_int_value(sp);
				save_int_value(age);
				save_int_value(f);
				save_int_value(k);
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
				save_identifier_string2((char*)"total_obs_catch_comp_end");
			
				gradient_structure::GRAD_STACK1->set_gradient_stack(dv_total_obs_catch_age_comp);
				fne++;
			}		
			k++;
		}
	}
}


void dv_total_obs_catch_age_comp(void)
{

	verify_identifier_string2((char*)"total_obs_catch_comp_end");
	unsigned long int pos_rw    = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();	
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_pop   = restore_long_int_value();
	unsigned k     = restore_int_value();
	unsigned f     = restore_int_value();
	unsigned age   = restore_int_value();
	unsigned sp    = restore_int_value();
	unsigned month = restore_int_value();
	unsigned year  = restore_int_value();
	unsigned t_count  = restore_int_value();
	const dvar_matrix_position uu_pos    = restore_dvar_matrix_position();
	const dvar_vector_position lf_pos    = restore_dvar_vector_position();
	const dvar_matrix_position Ctot_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position SNsum_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position s_asy_pos = restore_dvar_matrix_position();
	const dvar_matrix_position s_len_pos = restore_dvar_matrix_position();
	const dvar_matrix_position sslope_pos= restore_dvar_matrix_position();
	verify_identifier_string2((char*)"total_obs_catch_comp_begin");
	
	CCalpop* pop = (CCalpop*) pos_pop;
	PMap* map = (PMap*) pos_map;
	CMatrices* mat = (CMatrices*) pos_mat;
	CParam* param = (CParam*) pos_param;
	CReadWrite* rw= (CReadWrite*) pos_rw;

	const int s_func_type = param->s_func_type[f];
	const int nb_region   = param->nb_region;
	const int nb_region_sp= param->nb_region_sp_B[sp];

	dmatrix dfuu = restore_dvar_matrix_derivatives(uu_pos);
	dmatrix dfss = restore_dvar_matrix_derivatives(sslope_pos);
	dmatrix dfsl = restore_dvar_matrix_derivatives(s_len_pos);
	dmatrix dfsa = restore_dvar_matrix_derivatives(s_asy_pos);
	dmatrix dfSNsum = restore_dvar_matrix_derivatives(SNsum_pos);
	dvector dfLF	   = restore_dvar_vector_derivatives(lf_pos);
	dmatrix dfCtot  = restore_dvar_matrix_derivatives(Ctot_pos);

	const int imax = map->imax;
	const int imin = map->imin;
	const int jmax = map->jmax;
	const int jmin = map->jmin;

	dmatrix catch_obs;
	catch_obs.allocate(imin,imax,map->jinf,map->jsup); 
	catch_obs.initialize();

	rw->get_catch(*param, catch_obs, f, year, month, sp);

	const double W_age_mt = param->weight[sp][age] * 0.001; 
	const double area     = 1.852*param->deltaX*1.852*param->deltaY;
	const int catch_units = param->fishery_catch_units[f];
	const double C2Dunits = 1.0/(pow(W_age_mt,catch_units) * area);

	const double selectivity_f_age = pop->Selectivity(sp,f,age);

	double dflength, dfslope, dfasympt;
	param->dfselectivity(dfslope,dflength,dfasympt,sp,age,f,k);

	dmatrix uu;
	uu.allocate(map->imin1,map->imax1,map->jinf1,map->jsup1);
	uu = mat->density_before(sp,t_count,age);

	dmatrix EB(map->imin,map->imax,map->jinf,map->jsup);
	pop->Recomp_total_exploited_biomass(*map,*param,*mat,EB,catch_obs,pop->Selectivity(sp,f),f,sp,t_count);

	for (int j = jmax; j >= jmin; j--){
		const int imin = map->iinf[j];
		const int imax = map->isup[j];
		for (int i = imax; i >= imin; i--){
			if (map->carte(i,j) && EB(i,j)){

				double dfCa_est = 0.0;
				double dfLF_est = 0.0;
				double dfLF_pr_reg = 0.0;
				for (int a=nb_region_sp-1; a>=0; a--){	
					int r = param->area_sp_B[sp][a]-1;
					if (nb_region){ 
						if ((i>=map->regimin[a]) && (i<map->regimax[a]) &&(j>=map->regjmin[a]) && (j<map->regjmax[a])){
							//mat.C_N_sp_age_fishery[sp][age][k][r] = LF_pr_reg + LF_est; ;
							//dfLF_est   += dfLF(r);
							dfCa_est   += dfLF(r);
							dfLF_pr_reg+= dfLF(r);
							dfLF(r)     = 0;
						}
					} else {
						//mat.C_N_sp_age_fishery[sp][age][k][0] += LF_est ;
							//dfLF_est   += dfLF(0);
							dfCa_est   += dfLF(0);
							dfLF_pr_reg+= dfLF(0);
							dfLF(0)     = 0;
						}
					//LF_pr_reg = LF(r);
					dfLF(r) += dfLF_pr_reg;
					dfLF_pr_reg = 0;
				}

				//Ctot_age_obs.elem_value(i,j) = Ctot_pr(i,j) + LF_est * Cobs(i,j) * C2Dunits;
				dfLF_est      += catch_obs(i,j)*C2Dunits*mat->lat_correction(j)*dfCtot(i,j);
				//dfCtot_pr(i,j)+= dfCtot(i,j);
				//dfCtot(i,j)    = 0.0;

				//double LF_est = Selectivity(sp,f,age)*uu(i,j)/EB_tot;
				dfuu(i,j)   += dfLF_est * selectivity_f_age / EB(i,j);
				double dfS   = dfLF_est * uu(i,j) / EB(i,j);
				dfSNsum(i,j)-= dfLF_est * selectivity_f_age * uu(i,j) / pow(EB(i,j),2);
				dfLF_est     = 0.0;

				//double Ca_est = Selectivity(sp,f,age)*uu(i,j)*catch(k,i,j);
				dfuu(i,j)+= dfCa_est * selectivity_f_age * catch_obs(i,j);
				dfS      += dfCa_est * uu(i,j) * catch_obs(i,j);
				dfCa_est  = 0.0;


				//selectivity = selectivity_comp(age,f,k);
				dfss(i,j)+= dfslope*dfS;
				if (s_func_type>1)
					dfsl(i,j) += dflength*dfS;
				if (s_func_type==3)
					dfsa(i,j) += dfasympt*dfS;
				dfS = 0.0;	

			}
		}
	}
	dfCtot.save_dmatrix_derivatives(Ctot_pos);
	dfSNsum.save_dmatrix_derivatives(SNsum_pos);
	dfLF.save_dvector_derivatives(lf_pos);
	dfuu.save_dmatrix_derivatives(uu_pos);
	dfss.save_dmatrix_derivatives(sslope_pos);
	dfsl.save_dmatrix_derivatives(s_len_pos);
	dfsa.save_dmatrix_derivatives(s_asy_pos);
	
}

