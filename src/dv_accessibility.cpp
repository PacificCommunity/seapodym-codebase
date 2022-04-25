#include "VarSimtunaFunc.h"

///Main function with memory control and adjoint functions for: 
///1) accessibility to forage components (f_accessibility) or to
///their respective layers (f_accessibility_layer).
///2) average currents given the accessibility to the layer.
///Note, these functions assume the fixed forage structure with 
///six functional groups. Hence the configuration of these 
///groups can be controled only through eF parameters.
///Forward functions are in accessibility.cpp

void dv_accessibility_comp(void);
void dv_accessibility_ftype2_comp(void);
void dv_average_currents_comp(void);

double f_accessibility_layer(const double O2, const double T,double twosigsq, double temp_mean, double oxy_teta, double oxy_cr);
void dff_accessibility_layer(double dff_access, double& dfmean, double& dfsigma, double& dfteta, double& dfoxy_cr, const double T, const double O, const double sigma, const double mean, const double teta, const double oxy_cr, const double sigmasq, const double sigmacb);
double thermal_func_type2(const double T, const double temp_min, const double temp_max, const double delta1, const double delta2, const double delta3);
void dfthermal_func_type2(const double T, const double temp_min, const double temp_max, const double delta1, const double delta2, const double delta3, double& dfT, double& dftemp_min, double& dftemp_max, double& dfdelta1, double& dfdelta2, double& dfdelta3);
double f_accessibility_layer(const double O2, const double T, double temp_age, double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr);
void dff_accessibility_layer(double dff_access, double& dftemp_age, double& dftemp_max, double& dfdelta1, double& dfdelta2, double& dfdelta3, double& dfteta, double& dfoxy_cr, const double T, const double O, const double temp_age, const double temp_max, const double delta1, const double delta2, const double delta3, const double teta, const double oxy_cr);

int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);


void VarSimtunaFunc::Faccessibility(VarParamCoupled& param, VarMatrices& mat, const PMap& map, const int sp, const int jday, const int t_count, const int pop_built, const int tags_only, const ivector tags_age_solve)
{
	const int nb_forage  = param.get_nbforage();
	const int nb_layer   = param.nb_layer;
	const int a0 = param.sp_a0_adult[sp];
	const int nb_ages = param.sp_nb_cohorts[sp];

	dvariable temp_max = param.dvarsB_sst_spawning[sp];
	dvariable temp_min = param.dvarsB_sst_habitat[sp];
	dvariable oxy_teta = param.dvarsA_oxy_habitat[sp];
	dvariable oxy_cr   = param.dvarsB_oxy_habitat[sp]; 
	dvariable temp_age_slope = param.dvarsT_age_size_slope[sp];

	if (param.gaussian_thermal_function[sp]){
		dvariable sigma_ha = param.dvarsA_sst_habitat[sp];
		dvmatr3 = sigma_ha;
		dvariable sigma_hs = param.dvarsA_sst_spawning[sp];
		dvmatr7 = sigma_hs;
	}
	if (!param.gaussian_thermal_function[sp]){
		dvariable delta3 = param.dvarsThermal_func_delta[2][sp];
		dvmatr3 = delta3;

		dvariable delta1 = param.dvarsThermal_func_delta[0][sp];
		dvariable delta2 = param.dvarsThermal_func_delta[1][sp];
		dvmatr7 = delta1;
		dvmatr8 = delta2;
	}
	dvmatr1 = temp_max;
	dvmatr2 = temp_min;
	dvmatr4 = oxy_teta;
	dvmatr5 = oxy_cr;
	dvmatr6 = temp_age_slope;

	for (int age=a0; age<nb_ages; age++){
	    if (param.age_compute_habitat[sp][age]!=param.age_compute_habitat[sp][age-1]){
		if (!tags_only || tags_age_solve(age)){

			Faccessibility_comp(param,mat,map,value(temp_max),value(oxy_teta),value(oxy_cr),sp,age,jday,t_count);

			save_identifier_string2((char*)"Faccessibility_comp_begin");
			temp_max.save_prevariable_value();
			dvmatr1.save_dvar_matrix_position();
			temp_min.save_prevariable_value();
			dvmatr2.save_dvar_matrix_position();
			dvmatr3.save_dvar_matrix_position();
			oxy_teta.save_prevariable_value();
			dvmatr4.save_dvar_matrix_position();
			oxy_cr.save_prevariable_value();
			dvmatr5.save_dvar_matrix_position();
			temp_age_slope.save_prevariable_value();
			dvmatr6.save_dvar_matrix_position();
			dvmatr7.save_dvar_matrix_position();
			if (!param.gaussian_thermal_function[sp]){
				dvmatr8.save_dvar_matrix_position();
			}
			for (int n=0; n<nb_forage; n++)
				mat.dvarF_access(n,age).save_dvar_matrix_position();
			for (int l=0; l<nb_layer; l++)
				mat.dvarZ_access(l,age).save_dvar_matrix_position();
			unsigned long int cmat   = (unsigned long int)&mat;
			save_long_int_value(cmat);
			unsigned long int cparam = (unsigned long int)&param;
			save_long_int_value(cparam);
			unsigned long int pmap   = (unsigned long int)&map;
			save_long_int_value(pmap);
			save_int_value(jday);
			save_int_value(t_count);
			//save_int_value(pop_built);
			save_int_value(sp);
			save_int_value(age);
			save_identifier_string2((char*)"Faccessibility_comp_end");
			if (param.gaussian_thermal_function[sp]){
				gradient_structure::GRAD_STACK1->set_gradient_stack(dv_accessibility_comp);
			} else 	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_accessibility_ftype2_comp);
		}
	    }
	}
}

void VarSimtunaFunc::Average_currents(VarParamCoupled& param, VarMatrices& mat, const PMap& map, int age, const int t_count, const int pop_built)
{
	const int nb_layer  = param.nb_layer;
	mat.dvarsU.initialize();
	mat.dvarsV.initialize();

	Average_currents_comp(param,mat,map,age,t_count);

	save_identifier_string2((char*)"Average_currents_comp_begin");
	mat.dvarsU.save_dvar_matrix_position();	
	mat.dvarsV.save_dvar_matrix_position();	
	for (int l=0; l<nb_layer; l++)
		mat.dvarZ_access[l][age].save_dvar_matrix_position();
	unsigned long int cparam = (unsigned long int)&param;
	save_long_int_value(cparam);
	unsigned long int pmap   = (unsigned long int)&map;
	save_long_int_value(pmap);
	unsigned long int cmat   = (unsigned long int)&mat;
	save_long_int_value(cmat);
	save_int_value(t_count);
	//save_int_value(pop_built);
	save_identifier_string2((char*)"Average_currents_comp_end");

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_average_currents_comp);
}


/*
//previous (with p/2 and without re-scaling) version of dfthermal_func_type2
double dfthermal_func_type2(const double T, const double temp_min, const double temp_max, const double delta1, const double delta2, const double delta3, double& dfT, double& dftemp_min, double& dftemp_max, double& dfdelta1, double& dfdelta2, double& dfdelta3)
{
	double temp_diff = temp_max - temp_min;
	if (temp_diff<0) {cout << "temp_max cannot be smaller than temp_min! exit here..." << endl; exit(1);}
	double s_low = 0.05+delta1*temp_diff;//delta1=0.025 - test pars
	double s_hi  = 0.05+delta2*temp_diff; //delta2=0.015
	double p     =-(0.5+delta3*temp_diff); //delta3=0.2
	double dfp   = 0.0;
	double dfs_low = 0.0;
	double dfs_hi  = 0.0;

	double expr1 = pow(s_low,T-temp_min);
	double expr2 = p * pow(1.0+expr1,p-1);
	double expr3 = pow(s_hi,temp_max-T);
	double expr4 = (0.5*p) * pow(1.0+expr3,0.5*p-1);
	double y1    = pow(1+expr1,p);
	double y2    = pow(1+expr3,0.5*p);

	dftemp_min -= expr2 * expr1 * log(s_low) * y2 * dfT;
	dfs_low    += expr2 * (T-temp_min) * pow(s_low,T-temp_min-1) * y2 * dfT;
	dftemp_max += expr4 * expr3 * log(s_hi) * y1 * dfT;
	dfs_hi     += expr4 * (temp_max-T) * pow(s_hi,temp_max-T-1) * y1* dfT;
	dfp        += y1 * y2 * (log(1.0+expr1) + 0.5 * log(1.0+expr3)) * dfT;
	dfT         = 0.0;


	dftemp_min += delta3 * dfp; 
	dftemp_max -= delta3 * dfp;
	dftemp_min -= delta1 * dfs_low;
	dftemp_max += delta1 * dfs_low;
	dftemp_min -= delta2 * dfs_hi;
	dftemp_max += delta2 * dfs_hi;
	dfdelta1   += temp_diff * dfs_low;
	dfdelta2   += temp_diff * dfs_hi;
	dfdelta3   -= temp_diff * dfp;
}
*/

void dfthermal_func_type2(const double T, const double temp_min, const double temp_max, const double delta1, const double delta2, const double delta3, double& dfT, double& dftemp_min, double& dftemp_max, double& dfdelta1, double& dfdelta2, double& dfdelta3)
{
	double temp_diff = temp_max - temp_min;
	if (temp_diff<0) {cout << "temp_max cannot be smaller than temp_min! exit here..." << endl; exit(1);}
	double s_low = 0.05+delta1*temp_diff;
	double s_hi  = 0.05+delta2*temp_diff;	
	double a1    = log(s_low);
	double a2    = log(s_hi); 
	double p     = -(1.0+delta3*temp_diff);

	double E1 = 1.0+exp(a1*(T-temp_min));
	double E2 = 1.0+exp(a2*(temp_max-T));
	double f_T_pr = pow(E1*E2,p);

	double expr1 = (log(a2/a1)+a1*temp_min+a2*temp_max);
	double a1pa2 = a1+a2;
	double xmax = expr1/a1pa2;

	double E1max = 1.0+exp(a1*(xmax-temp_min));
	double E2max = 1.0+exp(a2*(temp_max-xmax));
	double f_T_xmax = pow(E1max*E2max,p);

	//f_T = f_T_pr/f_T_xmax;
	double dff_T_pr   = dfT/f_T_xmax;
	double dff_T_xmax = -f_T_pr*dfT/(f_T_xmax*f_T_xmax);
	dfT = 0.0;

	//f_T_xmax = (E1max*E2max)^p;
	double dfE1max = p*pow(E1max,p-1)*pow(E2max,p)*dff_T_xmax;
	double dfE2max = p*pow(E2max,p-1)*pow(E1max,p)*dff_T_xmax;
	double dfp     = f_T_xmax * log(E1max*E2max)  *dff_T_xmax;
        //dff_T_xmax     = 0.0; //local, no need to put to 0

	//E1max = 1+exp(a1*(xmax-temp_min));
	double expr = exp(a1*(xmax-temp_min))*dfE1max;
	double dfa1 = (xmax-temp_min) * expr;
	dftemp_min -= a1*expr;	
	double dfxmax = a1*expr;	

	//E2max = 1+exp(a2*(temp_max-xmax));
	expr = exp(a2*(temp_max-xmax))*dfE2max;
	double dfa2 = (temp_max-xmax) * expr; 
	dftemp_max += a2*expr;
	dfxmax -= a2*expr;

	//double xmax = expr1/a1pa2;
	dfa1 += dfxmax*((temp_min-1.0/a1)*a1pa2-expr1)/(a1pa2*a1pa2);
	dfa2 += dfxmax*((temp_max+1.0/a2)*a1pa2-expr1)/(a1pa2*a1pa2);
	dftemp_min += dfxmax * a1/a1pa2;
	dftemp_max += dfxmax * a2/a1pa2;

	//f_T_pr = (E1*E2)^p;
	double dfE1 = p*pow(E1,p-1)*pow(E2,p)*dff_T_pr;
	double dfE2 = p*pow(E2,p-1)*pow(E1,p)*dff_T_pr;
	dfp        += f_T_pr * log(E1*E2) * dff_T_pr;

	//double E1 = 1+exp(a1*(T-temp_min));
	expr = exp(a1*(T-temp_min))*dfE1;
	dfa1       += (T-temp_min) * expr;
	dftemp_min -= a1*expr;	

	//double E2 = 1+exp(a1*(temp_max-T));
	expr = exp(a2*(temp_max-T))*dfE2;
	dfa2       += (temp_max-T) * expr; 
	dftemp_max += a2*expr;

	//a1=log(0.05+delta1*temp_diff); a2=log(0.05+delta2*temp_diff); p=-(1.0+delta3*temp_diff);
	expr = delta1 * dfa1 / (s_low) + delta2 * dfa2 / (s_hi) - delta3 * dfp; 
	dftemp_min -= expr; 
	dftemp_max += expr;
	dfdelta1   += temp_diff * dfa1 / (s_low);
	dfdelta2   += temp_diff * dfa2 / (s_hi);
	dfdelta3   -= temp_diff * dfp;
}


//this is the adjoint for f_accessibility_layer function
void dff_accessibility_layer(double dff_access, double& dftemp_age, double& dfsigma, double& dfteta, double& dfoxy_cr, const double T, const double O, const double sigma, const double temp_age, const double teta, const double oxy_cr, const double sigmasq, const double sigmacb)
{
	const double exprO1 = pow(teta,O-oxy_cr);
	const double exprO2 = (1.0+exprO1)*(1.0+exprO1);


	//first recompute f_temp and f_oxy
	double f_oxy   = 1.0 / (1.0+exprO1);
	double f_temp  = exp(-pow(T-temp_age,2.0) / (2.0*sigmasq));

	//f_access = f_temp*f_oxy + 1e-4;
	double dfT = f_oxy * dff_access;
	double dfO = f_temp * dff_access;
	dff_access = 0.0;

	//double f_temp= exp(-pow(tempn(l,i,j)-temp_age,2) / (2*sigma*sigma));
	dftemp_age+= ((T-temp_age) / sigmasq) * f_temp * dfT;
	dfsigma   += (pow(T-temp_age,2.0) / sigmacb) * f_temp * dfT;
	dfT        = 0.0;

	//double f_oxy = 1.0 / (1.0+ pow(oxy_teta,O2 - oxy_cr)); 
	dfoxy_cr += (exprO1 * log(teta) / exprO2)  * dfO;
	dfteta   += ((oxy_cr-O)*exprO1 / (teta*exprO2)) * dfO;		
	dfO      = 0.0;
}

//this is the adjoint of the f_accessibility_layer with thermal function type2 
void dff_accessibility_layer(double dff_access, double& dftemp_age, double& dftemp_max, double& dfdelta1, double& dfdelta2, double& dfdelta3, double& dfteta, double& dfoxy_cr, const double T, const double O, const double temp_age, const double temp_max, const double delta1, const double delta2, const double delta3, const double teta, const double oxy_cr)
{
	const double exprO1 = pow(teta,O-oxy_cr);
	const double exprO2 = (1.0+exprO1)*(1.0+exprO1);

	//first recompute f_temp and f_oxy
	double f_oxy   = 1.0 / (1.0+exprO1);
	double f_temp = thermal_func_type2(T,temp_age,temp_max,delta1,delta2,delta3);

	//f_access = f_temp*f_oxy + 1e-4;
	double dfT = f_oxy * dff_access;
	double dfO = f_temp * dff_access;
	dff_access  = 0.0;

	//f_temp = thermal_func_type2(T,temp_age,temp_max,delta1,delta2,delta3);
	dfthermal_func_type2(T,temp_age,temp_max,delta1,delta2,delta3,dfT,dftemp_age,dftemp_max,dfdelta1,dfdelta2,dfdelta3);

	//double f_oxy = 1.0 / (1.0+ pow(oxy_teta,O2 - oxy_cr)); 
	dfoxy_cr += (exprO1 * log(teta) / exprO2)  * dfO;
	dfteta   += ((oxy_cr-O)*exprO1 / (teta*exprO2)) * dfO;		
	dfO      = 0.0;
}


void dv_accessibility_comp(void)
{
	verify_identifier_string2((char*)"Faccessibility_comp_end");
	const int age      = restore_int_value();
	const int sp       = restore_int_value();
	//unsigned flag      = restore_int_value();
	unsigned t_count   = restore_int_value();
	unsigned jday      = restore_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	const dvar_matrix_position Zpos2  = restore_dvar_matrix_position();
	const dvar_matrix_position Zpos1  = restore_dvar_matrix_position();
	const dvar_matrix_position Zpos0  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos5  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos4  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos3  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos2  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos1  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos0  = restore_dvar_matrix_position();
	const dvar_matrix_position Sigma_hs_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Temp_age_slope_pos    = restore_dvar_matrix_position();
	double temp_age_slope				 = restore_prevariable_value();
	const dvar_matrix_position Oxy_cr_pos    = restore_dvar_matrix_position();
	double oxy_cr 				 = restore_prevariable_value();
	const dvar_matrix_position Oxy_teta_pos  = restore_dvar_matrix_position();
	double oxy_teta 			 = restore_prevariable_value();
	const dvar_matrix_position Sigma_ha_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Temp_min_pos  = restore_dvar_matrix_position();
	double temp_min 			 = restore_prevariable_value();
	const dvar_matrix_position Temp_max_pos  = restore_dvar_matrix_position();
	double temp_max 			 = restore_prevariable_value();
	verify_identifier_string2((char*)"Faccessibility_comp_begin");


	dmatrix dfOxy_cr   = restore_dvar_matrix_derivatives(Oxy_cr_pos);
	dmatrix dfTemp_min = restore_dvar_matrix_derivatives(Temp_min_pos);
	dmatrix dfTemp_max = restore_dvar_matrix_derivatives(Temp_max_pos);
	dmatrix dfSigma_ha = restore_dvar_matrix_derivatives(Sigma_ha_pos);
	dmatrix dfSigma_hs = restore_dvar_matrix_derivatives(Sigma_hs_pos);
	dmatrix dfOxy_teta = restore_dvar_matrix_derivatives(Oxy_teta_pos);
	dmatrix dfTemp_age_slope = restore_dvar_matrix_derivatives(Temp_age_slope_pos);

	CParam* param = (CParam*) pos_param;

	const int nbf = param->get_nbforage();
	const int nbl = param->nb_layer;
	ivector day_layer(0,nbf-1); day_layer = param->day_layer;
	ivector night_layer(0,nbf-1); night_layer = param->night_layer;
	
	PMap* map = (PMap*) pos_map;
	CMatrices* mat = (CMatrices*) pos_mat;

	const int imax = map->imax;
	const int imin = map->imin;
	d3_array dfF_access, dfZ_access;
	dfF_access.allocate(0,nbf-1);
	dfZ_access.allocate(0,nbl-1);
	for (int n=0; n<nbf; n++){
		dfF_access[n].allocate(imin, imax, map->jinf, map->jsup);
		dfF_access[n].initialize();
	}
	for (int l=0; l<nbl; l++){
		dfZ_access[l].allocate(imin, imax, map->jinf, map->jsup);
		dfZ_access[l].initialize();
	}

	dfF_access[0]  = restore_dvar_matrix_derivatives(Fpos0);
	dfF_access[1]  = restore_dvar_matrix_derivatives(Fpos1);
	dfF_access[2]  = restore_dvar_matrix_derivatives(Fpos2);
	dfF_access[3]  = restore_dvar_matrix_derivatives(Fpos3);
	dfF_access[4]  = restore_dvar_matrix_derivatives(Fpos4);
	dfF_access[5]  = restore_dvar_matrix_derivatives(Fpos5);
	dfZ_access[0]  = restore_dvar_matrix_derivatives(Zpos0);
	dfZ_access[1]  = restore_dvar_matrix_derivatives(Zpos1);
	dfZ_access[2]  = restore_dvar_matrix_derivatives(Zpos2);

	dmatrix vld(imin, imax, map->jinf, map->jsup);
	vld = mat->vld(t_count);
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

	const double L_max = param->length[sp][param->sp_nb_cohorts[sp]-1];
	const double L_age = param->length[sp][age];
	const double  R    = pow(L_age/L_max,temp_age_slope);

	const double W_max = param->weight[sp][param->sp_nb_cohorts[sp]-1];
	const double W_age = param->weight[sp][age];
	const double RW    = W_age/W_max;
	//const double RW    = L_age/L_max;

	const double sigma_ha = param->sigma_ha[sp][age];
	const double temp_age = param->temp_age[sp][age];
	const double sigsq = sigma_ha*sigma_ha;
	const double sigcb = sigsq*sigma_ha;
	const double twosigsq = 2.0*sigsq;
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nlayer = map->carte(i,j);
			//if (nlayer>0){	
			if (nlayer>0 && nlayer<=nbl){

				const double DL = mat->daylength(jday,j)/24.0;

				double dftemp_age = 0.0;

                                //recompute lf_access and sumL
				dvector lf_access(0,nbl-1);
				lf_access.initialize();
				double sumL = 0.0;
				dvector l_access(0,nbl-1);
				l_access.initialize();
				//for (int l=0; l<nlayer; l++){
				for (int l=0; l<nbl; l++){
					l_access(l)  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),twosigsq,temp_age,oxy_teta,oxy_cr);
					// weighted by the Forage biomass
					for (int n=0;n<nbf;n++){
						if (day_layer[n]==l)   lf_access(l) += l_access(l)*forage(n,i,j)*DL; 
						if (night_layer[n]==l) lf_access(l) += l_access(l)*forage(n,i,j)*(1.0-DL);
					}	
					sumL += lf_access(l);			
				}
				//end of recomputation section	
	
				dvector dfl_access(0,nbl-1);
				dvector dflf_access(0,nbl-1);
				dvector dflf_access_rel(0,nbl-1);
				dfl_access.initialize();
				dflf_access.initialize();
				dflf_access_rel.initialize();

				for (int n=nbf-1; n>=0; n--){		
	
					const int dl = day_layer[n];
					const int nl = night_layer[n];

					if (dl < nlayer){
						//F_access = (l_access[dl]*DL + l_access[nl]*(1-DL));
						dfl_access(dl) += DL * dfF_access(n,i,j); 
						dfl_access(nl) += (1.0-DL) * dfF_access(n,i,j); 
						dfF_access(n,i,j) = 0.0;
					}
					dfF_access(n,i,j) = 0.0;
				}

				double expr1 = sumL;
				double dfsumL = 0.0;
				for (int l=nbl-1;l>=0;l--){	

					//Z_access(l) = lf_access(l)/(sumL);
					dfsumL -= lf_access(l)*dfZ_access(l,i,j)/(expr1*expr1);
					dflf_access(l) += dfZ_access(l,i,j)/expr1;
					dfZ_access(l,i,j) = 0.0;
				}

				double dfR = 0.0;
				double dfsigma_ha = 0.0;
				for (int l=nbl-1;l>=0;l--){

					//sumL += lf_access[l];
					dflf_access(l) += dfsumL;

					for (int n=nbf-1;n>=0;n--){
						if (night_layer[n]==l){	
							//lf_access[l]+= l_access(l) * mat.forage[n][i][j]*(1-DL);
							dfl_access(l) += forage(n,i,j) * (1.0-DL) * dflf_access(l);
						}
						if (day_layer[n]==l){
							//lf_access[l]+= l_access(l) * mat.forage[n][i][j]*DL; 
							dfl_access(l) += forage(n,i,j) * DL * dflf_access(l);
						}
					}

					//l_access(l)  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),twosigsq,teta,oxy_teta,oxy_cr);
					//dff_accessibility_layer(dfl_access(l),dftemp_age,dfSigma_ha(i,j),dfOxy_teta(i,j),dfOxy_cr(i,j),
					dff_accessibility_layer(dfl_access(l),dftemp_age,dfsigma_ha,dfOxy_teta(i,j),dfOxy_cr(i,j),
								tempn(l,i,j),oxygen(l,i,j),sigma_ha,temp_age,oxy_teta,oxy_cr,sigsq,sigcb);
								
					dfl_access(l) = 0.0;

					//double sigma_ha = RW * (sigma_max-sigma_min)+sigma_min;
					dfSigma_ha(i,j) += RW * dfsigma_ha;
					dfSigma_hs(i,j) -= (RW-1) * dfsigma_ha;
					dfsigma_ha = 0.0;


					//double temp_sp_age = R * (temp_min-temp_max) + temp_max;
					dfTemp_max(i,j) -= (R-1) * dftemp_age;
					dfTemp_min(i,j) += R * dftemp_age;
					dfR 		+= (temp_min-temp_max)* dftemp_age;
					dftemp_age 	 = 0.0;

					//double R = pow(L_age/L_max,temp_age_slope);
					dfTemp_age_slope(i,j) += R * log(L_age/L_max) * dfR;
					dfR = 0.0;
				}				
			}	
		}
	}
	dfTemp_age_slope.save_dmatrix_derivatives(Temp_age_slope_pos);
	dfOxy_cr.save_dmatrix_derivatives(Oxy_cr_pos);
	dfTemp_min.save_dmatrix_derivatives(Temp_min_pos); 
	dfTemp_max.save_dmatrix_derivatives(Temp_max_pos);
	dfSigma_ha.save_dmatrix_derivatives(Sigma_ha_pos);
	dfSigma_hs.save_dmatrix_derivatives(Sigma_hs_pos);
	dfOxy_teta.save_dmatrix_derivatives(Oxy_teta_pos);
	dfF_access[0].save_dmatrix_derivatives(Fpos0);
	dfF_access[1].save_dmatrix_derivatives(Fpos1);
	dfF_access[2].save_dmatrix_derivatives(Fpos2);
	dfF_access[3].save_dmatrix_derivatives(Fpos3);
	dfF_access[4].save_dmatrix_derivatives(Fpos4);
	dfF_access[5].save_dmatrix_derivatives(Fpos5);
	dfZ_access[0].save_dmatrix_derivatives(Zpos0);
	dfZ_access[1].save_dmatrix_derivatives(Zpos1);
	dfZ_access[2].save_dmatrix_derivatives(Zpos2);
}

void dv_accessibility_ftype2_comp(void)
{
	verify_identifier_string2((char*)"Faccessibility_comp_end");
	const int age      = restore_int_value();
	const int sp       = restore_int_value();
	//unsigned flag      = restore_int_value();
	unsigned t_count   = restore_int_value();
	unsigned jday      = restore_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	const dvar_matrix_position Zpos2  = restore_dvar_matrix_position();
	const dvar_matrix_position Zpos1  = restore_dvar_matrix_position();
	const dvar_matrix_position Zpos0  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos5  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos4  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos3  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos2  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos1  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos0  = restore_dvar_matrix_position();
	const dvar_matrix_position Tfunc2_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Tfunc1_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Temp_age_slope_pos = restore_dvar_matrix_position();
	double temp_age_slope			      = restore_prevariable_value();
	const dvar_matrix_position Oxy_cr_pos    = restore_dvar_matrix_position();
	double oxy_cr 				 = restore_prevariable_value();
	const dvar_matrix_position Oxy_teta_pos  = restore_dvar_matrix_position();
	double oxy_teta 			 = restore_prevariable_value();
	const dvar_matrix_position   Tfunc3_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position Temp_min_pos  = restore_dvar_matrix_position();
	double temp_min 			 = restore_prevariable_value();
	const dvar_matrix_position Temp_max_pos  = restore_dvar_matrix_position();
	double temp_max 			 = restore_prevariable_value();
	verify_identifier_string2((char*)"Faccessibility_comp_begin");


	dmatrix dfTemp_age_slope = restore_dvar_matrix_derivatives(Temp_age_slope_pos);
	dmatrix dfTemp_min = restore_dvar_matrix_derivatives(Temp_min_pos);
	dmatrix dfTemp_max = restore_dvar_matrix_derivatives(Temp_max_pos);
	dmatrix dfDelta1   = restore_dvar_matrix_derivatives(Tfunc1_pos);
	dmatrix dfDelta2   = restore_dvar_matrix_derivatives(Tfunc2_pos);
	dmatrix dfDelta3   = restore_dvar_matrix_derivatives(Tfunc3_pos);
	dmatrix dfOxy_cr   = restore_dvar_matrix_derivatives(Oxy_cr_pos);
	dmatrix dfOxy_teta = restore_dvar_matrix_derivatives(Oxy_teta_pos);

	CParam* param = (CParam*) pos_param;

	const int nbf = param->get_nbforage();
	const int nbl = param->nb_layer;
	ivector day_layer(0,nbf-1); day_layer = param->day_layer;
	ivector night_layer(0,nbf-1); night_layer = param->night_layer;
	
	PMap* map = (PMap*) pos_map;
	CMatrices* mat = (CMatrices*) pos_mat;

	const int imax = map->imax;
	const int imin = map->imin;
	d3_array dfF_access, dfZ_access;
	dfF_access.allocate(0,nbf-1);
	dfZ_access.allocate(0,nbl-1);
	for (int n=0; n<nbf; n++){
		dfF_access[n].allocate(imin, imax, map->jinf, map->jsup);
		dfF_access[n].initialize();
	}
	for (int l=0; l<nbl; l++){
		dfZ_access[l].allocate(imin, imax, map->jinf, map->jsup);
		dfZ_access[l].initialize();
	}

	dfF_access[0]  = restore_dvar_matrix_derivatives(Fpos0);
	dfF_access[1]  = restore_dvar_matrix_derivatives(Fpos1);
	dfF_access[2]  = restore_dvar_matrix_derivatives(Fpos2);
	dfF_access[3]  = restore_dvar_matrix_derivatives(Fpos3);
	dfF_access[4]  = restore_dvar_matrix_derivatives(Fpos4);
	dfF_access[5]  = restore_dvar_matrix_derivatives(Fpos5);
	dfZ_access[0]  = restore_dvar_matrix_derivatives(Zpos0);
	dfZ_access[1]  = restore_dvar_matrix_derivatives(Zpos1);
	dfZ_access[2]  = restore_dvar_matrix_derivatives(Zpos2);


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

	const double L_max = param->length[sp][param->sp_nb_cohorts[sp]-1];
	const double L_age = param->length[sp][age];
	const double  R    = pow(L_age/L_max,temp_age_slope);

	const double temp_age = param->temp_age[sp][age];
	const double delta1   = param->thermal_func_delta[0][sp];
	const double delta2   = param->thermal_func_delta[1][sp];
	const double delta3   = param->thermal_func_delta[2][sp];
	
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nlayer = map->carte(i,j);
			//if (nlayer>0){	
			if (nlayer>0 && nlayer<=nbl){

				const double DL = mat->daylength(jday,j)/24.0;

				double dftemp_age = 0.0;

                                //recompute lf_access and sumL
				dvector lf_access(0,nbl-1);
				lf_access.initialize();
				double sumL = 0.0;
				dvector l_access(0,nbl-1);
				l_access.initialize();			
				//for (int l=0; l<nlayer; l++){
				for (int l=0; l<nbl; l++){
					l_access(l)  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),temp_age,temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);
					// weighted by the Forage biomass
					for (int n=0;n<nbf;n++){
						if (day_layer[n]==l)   lf_access(l) += l_access(l)*forage(n,i,j)*DL; 
						if (night_layer[n]==l) lf_access(l) += l_access(l)*forage(n,i,j)*(1-DL);
					}	
					sumL += lf_access(l);	
				}
				//end of recomputation section	
	
			
				dvector dfl_access(0,nbl-1);
				dvector dflf_access(0,nbl-1);
				dvector dflf_access_rel(0,nbl-1);
				dfl_access.initialize();
				dflf_access.initialize();
				dflf_access_rel.initialize();

				for (int n=nbf-1; n>=0; n--){		
	
					const int dl = day_layer[n];
					const int nl = night_layer[n];

					if (dl < nlayer){
						//F_access = (l_access[dl]*DL + l_access[nl]*(1-DL));
						dfl_access(dl) += DL * dfF_access(n,i,j); 
						dfl_access(nl) += (1.0-DL) * dfF_access(n,i,j); 
						dfF_access(n,i,j) = 0.0;
					}
					dfF_access(n,i,j) = 0.0;
				}

				double expr1 = sumL;
				double dfsumL = 0.0;
				for (int l=nbl-1;l>=0;l--){

					//Z_access(l) = lf_access(l)/sumL;
					dfsumL -= lf_access(l)*dfZ_access(l,i,j)/(expr1*expr1);
					dflf_access(l) += dfZ_access(l,i,j)/expr1;
					dfZ_access(l,i,j) = 0.0;
				}

				double dfR = 0.0;
				for (int l=nbl-1;l>=0;l--){

					//sumL += lf_access[l];
					dflf_access(l) += dfsumL;

					for (int n=nbf-1;n>=0;n--){
						if (night_layer[n]==l){	
							//lf_access[l]+= l_access(l) * mat.forage[n][i][j]*(1-DL);
							dfl_access(l) += forage(n,i,j) * (1-DL) * dflf_access(l);
						}
						if (day_layer[n]==l){
							//lf_access[l]+= l_access(l) * mat.forage[n][i][j]*DL; 
							dfl_access(l) += forage(n,i,j) * DL * dflf_access(l);
						}
					}

					//l_access(l)  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),temp_age,temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);
					dff_accessibility_layer(dfl_access(l),dftemp_age,dfTemp_max(i,j),dfDelta1(i,j),dfDelta2(i,j),
								dfDelta3(i,j),dfOxy_teta(i,j),dfOxy_cr(i,j),tempn(l,i,j),oxygen(l,i,j),
								temp_age,temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);

					dfl_access(l) = 0.0;

					//double temp_sp_age = R * (temp_min-temp_max) + temp_max;
					dfTemp_max(i,j) -= (R-1) * dftemp_age;
					dfTemp_min(i,j) += R * dftemp_age;
					dfR 		+= (temp_min-temp_max)* dftemp_age;
					dftemp_age	 = 0.0;

					//double R = pow(L_age/L_max,temp_age_slope);
					dfTemp_age_slope(i,j) += R * log(L_age/L_max) * dfR;
					dfR = 0.0;
				}				
			}	
		}
	}

	dfTemp_age_slope.save_dmatrix_derivatives(Temp_age_slope_pos);
	dfTemp_min.save_dmatrix_derivatives(Temp_min_pos); 
	dfTemp_max.save_dmatrix_derivatives(Temp_max_pos);
	dfDelta1.save_dmatrix_derivatives(Tfunc1_pos);
	dfDelta2.save_dmatrix_derivatives(Tfunc2_pos);
	dfDelta3.save_dmatrix_derivatives(Tfunc3_pos);
	dfOxy_cr.save_dmatrix_derivatives(Oxy_cr_pos);
	dfOxy_teta.save_dmatrix_derivatives(Oxy_teta_pos);
	dfF_access[0].save_dmatrix_derivatives(Fpos0);
	dfF_access[1].save_dmatrix_derivatives(Fpos1);
	dfF_access[2].save_dmatrix_derivatives(Fpos2);
	dfF_access[3].save_dmatrix_derivatives(Fpos3);
	dfF_access[4].save_dmatrix_derivatives(Fpos4);
	dfF_access[5].save_dmatrix_derivatives(Fpos5);
	dfZ_access[0].save_dmatrix_derivatives(Zpos0);
	dfZ_access[1].save_dmatrix_derivatives(Zpos1);
	dfZ_access[2].save_dmatrix_derivatives(Zpos2);
}


void dv_average_currents_comp()
{
	verify_identifier_string2((char*)"Average_currents_comp_end");
	//unsigned flag      = restore_int_value();
	unsigned t_count   = restore_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	const dvar_matrix_position Zpos2  = restore_dvar_matrix_position();
	const dvar_matrix_position Zpos1  = restore_dvar_matrix_position();
	const dvar_matrix_position Zpos0  = restore_dvar_matrix_position();
	const dvar_matrix_position V_pos  = restore_dvar_matrix_position();	
	const dvar_matrix_position U_pos  = restore_dvar_matrix_position();	
	verify_identifier_string2((char*)"Average_currents_comp_begin");

	dmatrix dfU = restore_dvar_matrix_derivatives(U_pos);
	dmatrix dfV = restore_dvar_matrix_derivatives(V_pos);

	CParam* param = (CParam*) pos_param;
	PMap* map = (PMap*) pos_map;
	CMatrices* mat = (CMatrices*) pos_mat;
	const int nbl = param->nb_layer;
	const int imax = map->imax;
	const int imin = map->imin;

	d3_array dfZ_access;
	dfZ_access.allocate(0,nbl-1);
	for (int l=0; l<nbl; l++){
		dfZ_access[l].allocate(imin, imax, map->jinf, map->jsup);
		dfZ_access[l].initialize();
	}
	dfZ_access[0]  = restore_dvar_matrix_derivatives(Zpos0);
	dfZ_access[1]  = restore_dvar_matrix_derivatives(Zpos1);
	dfZ_access[2]  = restore_dvar_matrix_derivatives(Zpos2);

	d3_array un,vn;
	un.allocate(0,nbl-1);	
	vn.allocate(0,nbl-1);	
	
	for (int l=0; l<nbl;l++){
		un(l).allocate(imin, imax, map->jinf, map->jsup); 
		vn(l).allocate(imin, imax, map->jinf, map->jsup); 
		un(l) = mat->un(t_count,l);
		vn(l) = mat->vn(t_count,l);
	}

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nlayer = map->carte(i,j);
			if (nlayer>0 && nlayer<=nbl){
				for (int l=nbl-1; l>=0; l--){
					//param.dvarsU.elem_value(i,j) += mat.un(t,l,i,j) * Z_access;
					dfZ_access(l,i,j) += un(l,i,j) * dfU(i,j);
					//param.dvarsV.elem_value(i,j) += mat.vn(t,l,i,j) * Z_access;
					dfZ_access(l,i,j) += vn(l,i,j) * dfV(i,j);
				}
				dfU(i,j) = 0.0;
				dfV(i,j) = 0.0;
			}	
		}
	}

	dfU.save_dmatrix_derivatives(U_pos);
	dfV.save_dmatrix_derivatives(V_pos);
	dfZ_access[0].save_dmatrix_derivatives(Zpos0);
	dfZ_access[1].save_dmatrix_derivatives(Zpos1);
	dfZ_access[2].save_dmatrix_derivatives(Zpos2);
}
