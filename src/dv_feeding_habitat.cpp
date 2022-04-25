#include "VarSimtunaFunc.h"

///Main function with memory control and adjoint functions for: 
///feeding habitat for young and adult life stages, with or 
///without seasonal switch between habitats depending on the 
///migration flag. The function Seasonal_Habitat_Index is called
///only if the migration flag is activated.
///Note, these functions assume the fixed forage structure with 
///six functional groups. Hence the configuration of these 
///groups can be controled only through eF parameters.
///Forward functions are in feeding_habitat.cpp


void dv_Ha_comp(void);
void dv_Hf_comp(void);
double f_accessibility_layer(const double O, const double T,double twosigsq, double temp_mean, double oxy_teta, double oxy_cr);
double f_accessibility_layer(const double O2, const double T, double temp_age, double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr);
double hs_comp(double SST, double preys, double predators, const double a, const double b, const double c, const double d, const double e, const double ssv);
double pred_surface_comp(dvector forage, const double DL, const int nb_forage, ivector day_layer, ivector night_layer);
int save_identifier_string2(char* str);
void dfspawning_habitat_seasonality(double& dfa, double& dfb, double& dfc, double& dfH, const double a, const double b, const double c, const double pf_ratio, const double SST, const double itopo);
void dfsst_habitat(double& dfa, double& dfb, double& dfH, const double a, const double b, const double SST, const double itopo);
void dfsigma_hss_comp(dvector& length, double temp_min, double temp_max, const int sp_nb_cohorts, const int age, const int jday);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);
double daylength_comp(double lat, double jday, double pi);

void VarSimtunaFunc::Feeding_Habitat_Index(VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Hf, int sp, int age, const int jday, const int t_count){
	
	Hf.initialize();
	const int nb_forage = param.get_nbforage();

	dvmatr1 = param.dvarsEF_habitat[0][sp];
	dvmatr3 = param.dvarsEF_habitat[1][sp];
	dvmatr5 = param.dvarsEF_habitat[2][sp];
	dvmatr6 = param.dvarsEF_habitat[3][sp];
	dvmatr7 = param.dvarsEF_habitat[4][sp];
	dvmatr8 = param.dvarsEF_habitat[5][sp];

	Hf_comp(param,mat,map,Hf,sp,age,jday,t_count);

	save_identifier_string2((char*)"Hf_comp_begin");
	dvmatr1.save_dvar_matrix_position();
	dvmatr3.save_dvar_matrix_position();
	dvmatr5.save_dvar_matrix_position();
	dvmatr6.save_dvar_matrix_position();
	dvmatr7.save_dvar_matrix_position();
	dvmatr8.save_dvar_matrix_position();	
	Hf.save_dvar_matrix_position();
	for (int n=0; n<nb_forage; n++){
		mat.dvarF_access[n][age].save_dvar_matrix_position();
	}
	unsigned long int cmat   = (unsigned long int)&mat;
	save_long_int_value(cmat);
	unsigned long int cparam = (unsigned long int)&param;
	save_long_int_value(cparam);
	unsigned long int pmap   = (unsigned long int)&map;
	save_long_int_value(pmap);
	save_int_value(jday);
	save_int_value(t_count);
	save_int_value(sp);
	save_int_value(age);
	save_identifier_string2((char*)"Hf_comp_end");

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Hf_comp);
}

//Generalized habitat index with seasonal migration between spawning and feeding habitats
void VarSimtunaFunc::Seasonal_Habitat_Index(VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Hs, dvar_matrix& Ha, int sp, int age, const int jday, const int t_count){

	Ha_comp(param,mat,map,value(Hs),Ha,sp,jday);

	save_identifier_string2((char*)"Ha_comp_begin");
	Hs.save_dvar_matrix_position();
	Ha.save_dvar_matrix_position();
	mat.dvarSeasonSwitch[sp].save_dvar_matrix_position();
	unsigned long int cmat   = (unsigned long int)&mat;
	save_long_int_value(cmat);
	unsigned long int cparam = (unsigned long int)&param;
	save_long_int_value(cparam);
	unsigned long int pmap   = (unsigned long int)&map;
	save_long_int_value(pmap);
        save_int_value(jday);
        save_int_value(t_count);
        save_int_value(sp);
        save_int_value(age);
	save_identifier_string2((char*)"Ha_comp_end");

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Ha_comp);
}

void dv_Hf_comp(void)
{
	verify_identifier_string2((char*)"Hf_comp_end");
	const int age      = restore_int_value();
	const int sp       = restore_int_value();
	unsigned t_count   = restore_int_value();
	unsigned jday      = restore_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	const dvar_matrix_position Fpos5  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos4  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos3  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos2  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos1  = restore_dvar_matrix_position();
	const dvar_matrix_position Fpos0  = restore_dvar_matrix_position();
	const dvar_matrix_position Ha_pos = restore_dvar_matrix_position();
	const dvar_matrix_position E5_pos = restore_dvar_matrix_position();
	const dvar_matrix_position E4_pos = restore_dvar_matrix_position();
	const dvar_matrix_position E3_pos = restore_dvar_matrix_position();
	const dvar_matrix_position E2_pos = restore_dvar_matrix_position();
	const dvar_matrix_position E1_pos = restore_dvar_matrix_position();
	const dvar_matrix_position E0_pos = restore_dvar_matrix_position();		
	verify_identifier_string2((char*)"Hf_comp_begin");

	dmatrix dfHa 	   = restore_dvar_matrix_derivatives(Ha_pos);

	CParam* param = (CParam*) pos_param;
	PMap* map = (PMap*) pos_map;
	CMatrices* mat = (CMatrices*) pos_mat;

	double oxy_teta = param->a_oxy_habitat[sp];
	double oxy_cr = param->b_oxy_habitat[sp];
	double temp_age = param->temp_age[sp][age];
	double sigma_ha = param->sigma_ha[sp][age];
	const double twosigsq = 2.0*sigma_ha*sigma_ha;
	const double temp_max = param->b_sst_spawning(sp);
	const double delta1   = param->thermal_func_delta[0][sp];
	const double delta2   = param->thermal_func_delta[1][sp];
	const double delta3   = param->thermal_func_delta[2][sp];

	const int Tfunc_Gaussian = param->gaussian_thermal_function[sp];

	const int imax = map->imax;
	const int imin = map->imin;

	const int nbf = param->get_nbforage();
	d3_array dfF_access;
	dfF_access.allocate(0,nbf-1);
	for (int n=0; n<nbf; n++){
		dfF_access[n].allocate(imin, imax, map->jinf, map->jsup);
		dfF_access[n].initialize();
	}
	dfF_access[0]  = restore_dvar_matrix_derivatives(Fpos0);
	dfF_access[1]  = restore_dvar_matrix_derivatives(Fpos1);
	dfF_access[2]  = restore_dvar_matrix_derivatives(Fpos2);
	dfF_access[3]  = restore_dvar_matrix_derivatives(Fpos3);
	dfF_access[4]  = restore_dvar_matrix_derivatives(Fpos4);
	dfF_access[5]  = restore_dvar_matrix_derivatives(Fpos5);

	d3_array dfEF;
	dfEF.allocate(0,nbf-1);
	for (int n=0; n<nbf; n++){
		dfEF[n].allocate(imin, imax, map->jinf, map->jsup);
		dfEF[n].initialize();
	}	
	dfEF[0]  = restore_dvar_matrix_derivatives(E0_pos);
	dfEF[1]  = restore_dvar_matrix_derivatives(E1_pos);
	dfEF[2]  = restore_dvar_matrix_derivatives(E2_pos);
	dfEF[3]  = restore_dvar_matrix_derivatives(E3_pos);
	dfEF[4]  = restore_dvar_matrix_derivatives(E4_pos);
	dfEF[5]  = restore_dvar_matrix_derivatives(E5_pos);

	const int nbl = param->nb_layer;
	ivector day_layer(0,nbf-1); day_layer = param->day_layer;
	ivector night_layer(0,nbf-1); night_layer = param->night_layer;
	
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

	dvector eF;
	eF.allocate(0,nbf-1);
	eF.initialize();
	for (int n=0; n<nbf; n++){
		eF[n] = param->eF_habitat[n][sp];
	}	

	dvector T(0,nbl-1);   T.initialize();
	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nlayer = map->carte(i,j);
				if (nlayer>0){	

				const double DL = mat->daylength(jday,j)/24.0;

				double topo = map->itopo(i,j);

				for (int l=0; l<nbl; l++){
					T(l) = tempn(l,i,j);
				}

				double dffunc_Hf = 0.0;
                                //recompute func_Hf
				double Hf = 0;
				//dvector l_access(0,nlayer-1);
				dvector l_access(0,nbl-1);
				l_access.initialize();
				//for (int l=0; l<nlayer; l++){
				for (int l=0; l<nbl; l++){
					if (Tfunc_Gaussian){
						l_access(l)  = f_accessibility_layer(oxygen(l,i,j),T(l),twosigsq,
								temp_age,oxy_teta,oxy_cr);
					} else {
						l_access(l)  = f_accessibility_layer(oxygen(l,i,j),T(l),temp_age,
									temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);
					}
				}	

				for (int n=0; n<nbf; n++){
					double f_access = 0;
					if (day_layer[n] < nlayer)
						f_access = l_access[day_layer[n]]*DL + l_access[night_layer[n]]*(1-DL);

					double F = forage(n,i,j);
					F = eF[n] * F;
					Hf += f_access * F; 
				}
			
				//Ha.elem_value(i,j) = topo*func_Hf;
				dffunc_Hf += topo*dfHa(i,j);
				dfHa(i,j)  = 0.0;  

				//3. rotated hyperbola 
				//func_Hf = param.func_limit_one(Hf);
				double dfHf = param->dffunc_limit_one(Hf,dffunc_Hf);
				dffunc_Hf = 0.0;

			
				for (int n=nbf-1; n>=0; n--){		
					double f_access = 0;

					if (day_layer[n] < nlayer){
						f_access = l_access[day_layer[n]]*DL + l_access[night_layer[n]]*(1.0-DL);

						double F = forage(n,i,j);

						//Hf += f_access * eF[n] * forage[n][i][j]; 
						dfF_access(n,i,j) += eF[n] * F * dfHf;
						dfEF(n,i,j)  += f_access * F * dfHf;
					}
				}
			}	
		}
	}
//cout << age << " " << norm(dfF_access(0)) << " " << norm(dfF_access(1)) << " " << norm(dfF_access(2)) << " " << norm(dfF_access(3)) << " " << norm(dfF_access(4)) << " " << norm(dfF_access(5)) << " " << norm(dfEF(0)) << " " << norm(dfEF(0)) <<" " << norm(dfEF(1)) <<" " << norm(dfEF(2)) <<" " << norm(dfEF(3)) <<" " << norm(dfEF(4)) <<" " << norm(dfEF(5)) <<" " << norm(dfHa) <<endl;
	dfHa.save_dmatrix_derivatives(Ha_pos); 
	dfF_access[0].save_dmatrix_derivatives(Fpos0);
	dfF_access[1].save_dmatrix_derivatives(Fpos1);
	dfF_access[2].save_dmatrix_derivatives(Fpos2);
	dfF_access[3].save_dmatrix_derivatives(Fpos3);
	dfF_access[4].save_dmatrix_derivatives(Fpos4);
	dfF_access[5].save_dmatrix_derivatives(Fpos5);
	dfEF[0].save_dmatrix_derivatives(E0_pos);
	dfEF[1].save_dmatrix_derivatives(E1_pos);
	dfEF[2].save_dmatrix_derivatives(E2_pos);
	dfEF[3].save_dmatrix_derivatives(E3_pos);
	dfEF[4].save_dmatrix_derivatives(E4_pos);
	dfEF[5].save_dmatrix_derivatives(E5_pos);		
	
}
//Inna 01/13: derivative code for the new function - sigma varying as a function of season peak
//Inna 05/13: the use of function is temporarily on hold, still using fixed values of sigma (see seasonal_switch)
void dfsigma_hss_comp(double& dfTemp_max, double& dfTemp_min, double& dfsigma, double temp_min, double temp_max, const double R, const int jday)
{
	const double pi = 4.0*(atan(0.5)+atan(1.0/3.0));
	const double lat_max_sp = 60.0; //species dependent
	const int jday_maxDL = 171; //the day when DL is maximal in NH
	const int jday_minDL = 354; //the day when DL is minimal in NH
	const double season_peak = 120; //day 100 in NH
	double Jday = jday + jday_maxDL - season_peak;//tmp

	double DL       = daylength_comp(lat_max_sp,Jday, pi)/24.0;
//	double DL_max   = daylength_comp(lat_max_sp, jday_maxDL, pi)/24.0;
	double DL_min   = daylength_comp(lat_max_sp, jday_minDL, pi)/24.0;
//	double DL_start_NH  = daylength_comp(lat_max_sp, jday_maxDL-60, pi)/24.0;
	double DL_start_SH  = daylength_comp(lat_max_sp, jday_minDL-90, pi)/24.0;
	
	//double expr1_NH = (DL-DL_start_NH)/(DL_max-DL_start_NH);
//	double expr2_NH = (1-(DL-DL_start_NH)/(DL_max-DL_start_NH));

	//double expr1_SH = (DL-DL_start_SH)/(DL_min-DL_start_SH);
	double expr2_SH = (1-(DL-DL_start_SH)/(DL_min-DL_start_SH));

	double dfsigma_start = 0.0;
	double dftemp_age = 0.0;

//	if (DL>=DL_start_NH){
		//sigma = 1.0*expr1_NH+sigma_season_start*expr2_NH;
//		dfsigma_start += expr2_NH * dfsigma;
//		dfsigma = 0.0;
//	}
	if (DL<=DL_start_SH){
		//sigma = 1.0*expr1_SH+sigma_season_start*expr2_SH;
		dfsigma_start += expr2_SH * dfsigma;
		dfsigma = 0.0;
	}

	//double sigma_season_start = sigma;
	dfsigma += dfsigma_start;
	dfsigma_start = 0.0;
	
	//double sigma = 0.5*(T_a0-T_age);
	dfTemp_max    += 0.5 * dfsigma;
	dftemp_age -= 0.5 * dfsigma;
	dfsigma = 0.0;	

	//double temp_age = R * (temp_min-temp_max) + temp_max;
	dfTemp_max -= (R-1) * dftemp_age;
	dfTemp_min += R * dftemp_age;
	dftemp_age = 0.0;
}

void dv_Ha_comp(void)
{
	verify_identifier_string2((char*)"Ha_comp_end");
	const int age      = restore_int_value();
	const int sp       = restore_int_value();
	unsigned t_count   = restore_int_value();
	unsigned jday      = restore_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	const dvar_matrix_position Switch_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position Ha_pos        = restore_dvar_matrix_position();
	const dvar_matrix_position Hs_pos        = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"Ha_comp_begin");

	dmatrix dfHs 	   = restore_dvar_matrix_derivatives(Hs_pos);
	dmatrix dfHa 	   = restore_dvar_matrix_derivatives(Ha_pos);
	dmatrix dfSwitch   = restore_dvar_matrix_derivatives(Switch_pos);

	CParam* param = (CParam*) pos_param;
	PMap* map = (PMap*) pos_map;
	CMatrices* mat = (CMatrices*) pos_mat;

	double oxy_teta = param->a_oxy_habitat[sp];
        double oxy_cr   = param->b_oxy_habitat[sp];
	double temp_max = param->b_sst_spawning[sp];
        double temp_age = param->temp_age[sp][age];
        double sigma_ha = param->sigma_ha[sp][age];
        const double twosigsq = 2.0*sigma_ha*sigma_ha;
        const double delta1   = param->thermal_func_delta[0][sp];
        const double delta2   = param->thermal_func_delta[1][sp];
        const double delta3   = param->thermal_func_delta[2][sp];

	double a_hs = param->a_sst_larvae[sp];
	double b_hs = param->b_sst_larvae[sp];
	if (!param->uncouple_sst_larvae[sp]){
		a_hs = param->a_sst_spawning[sp];
		b_hs = param->b_sst_spawning[sp];
	}
	double c_hs = param->alpha_hsp_prey[sp];
	double d_hs = param->alpha_hsp_predator[sp];
	double e_hs = param->beta_hsp_predator[sp];
	double sigma_season = mat->sigma_season(sp,jday,age);//sigma_hss_comp(temp_age,temp_max, age,jday);

	const int Tfunc_Gaussian = param->gaussian_thermal_function[sp];

	const int imax = map->imax;
	const int imin = map->imin;

	const int nbf = param->get_nbforage();

	const int nbl = param->nb_layer;
        ivector day_layer(0,nbf-1); day_layer = param->day_layer;
        ivector night_layer(0,nbf-1); night_layer = param->night_layer;

        dmatrix np1(imin, imax, map->jinf, map->jsup);
        dmatrix sst(imin, imax, map->jinf, map->jsup);
//        dmatrix vld(imin, imax, map->jinf, map->jsup);
        np1 = mat->np1(t_count);
        sst = mat->sst(t_count);
//        vld = mat->vld(t_count);

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

        dvector eF;
        eF.allocate(0,nbf-1);
        eF.initialize();
        for (int n=0; n<nbf; n++){
                eF[n] = param->eF_habitat[n][sp];
        }

        dvector T(0,nbl-1);   T.initialize();
        dvector F(0,nbf-1);   F.initialize();
        const double pp_transform = param->pp_transform;

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nlayer = map->carte(i,j);
			//if (nlayer>0){	
			if (nlayer>0 && nlayer<=nbl){
				//recomputation section:
				double topo = map->itopo(i,j);
				const double DL = mat->daylength(jday,j)/24.0;
				double sfunc = mat->season_switch(sp,jday,j);


				for (int l=0; l<nbl; l++){
					T(l) = tempn(l,i,j);
				}

				double preys = np1(i,j) * pp_transform;

				for (int n=0; n<nbf; n++)
					F(n) = forage(n,i,j);		
				double predators = pred_surface_comp(F,DL,nbf,day_layer,night_layer);	

				double O2_l2  = oxygen(1,i,j);
				double f_oxy = 1.0/(1.0+pow(0.01,O2_l2-0.1));
				
				double Hs = hs_comp(sst(i,j),preys,predators,a_hs,b_hs,c_hs,d_hs,e_hs,sigma_season)*f_oxy;

                                double Hf = 0;
                                //dvector l_access(0,nlayer-1);
                                dvector l_access(0,nbl-1);
                                l_access.initialize();
                                //for (int l=0; l<nlayer; l++){
                                for (int l=0; l<nbl; l++){
                                        if (Tfunc_Gaussian){
                                                //l_access(l)  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),twosigsq,
                                                l_access(l)  = f_accessibility_layer(oxygen(l,i,j),T(l),twosigsq,
                                                                        temp_age,oxy_teta,oxy_cr);
                                        } else {
                                                //l_access(l)  = f_accessibility_layer(oxygen(l,i,j),tempn(l,i,j),temp_age,
                                                l_access(l)  = f_accessibility_layer(oxygen(l,i,j),T(l),temp_age,
                                                                        temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr);
                                        }

                                }
                                for (int n=0; n<nbf; n++){
                                        double f_access = 0;
                                        if (day_layer[n] < nlayer)
                                                f_access = l_access[day_layer[n]]*DL + l_access[night_layer[n]]*(1-DL);

                                        double F = forage(n,i,j);
                                        F = eF[n] * F;

                                        Hf += f_access * F;
                                }
                                double Hf_func = param->func_limit_one(Hf);
				

				double dfHf = 0.0;
				//Ha(i,j) = map.itopo(i,j)*(func_Hs*sfunc + func_Hf*(1-sfunc));
				dfHs(i,j) += topo*sfunc*dfHa(i,j);
				dfHf      += topo*(1.0-sfunc)*dfHa(i,j);
				//dfSwitch(jday,j) += topo*(Hs-Hf)*dfHa(i,j);
				dfSwitch(i,j) += topo*(Hs-Hf_func)*dfHa(i,j);
				dfHa(i,j)  = 0.0; 

				//double func_Hf = Ha(i,j);
				dfHa(i,j) += dfHf;				
			}	
		}
	}
//if (age>20) TTRACE(age,norm(dfTemp_max))
	dfHs.save_dmatrix_derivatives(Hs_pos); 
	dfHa.save_dmatrix_derivatives(Ha_pos); 
	dfSwitch.save_dmatrix_derivatives(Switch_pos); 
}
