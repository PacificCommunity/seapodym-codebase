#include "calpop.h"

///Main function with memory control and adjoint functions for: 
///precaldia and caldia functions, which perform discrete approximation
///of advection-diffusion equation and precomputes diagonal coefficients.
///Forward functions are in caldia.cpp

void dv_caldia(void);
void dv_caldia_UV(void);
void f_accessibility(dvector& l_access, dvector& lf_access, const dvector forage, const dvector O2, const dvector T, double twosigsq, double temp_mean, double oxy_teta, double oxy_cr, const int nl, const int nb_forage, const ivector day_layer, const ivector night_layer, const double DL);
void f_accessibility(dvector& l_access, dvector& lf_access, const dvector forage, const dvector O2, const dvector T, double temp_mean, 
		double temp_max, double delta1, double delta2, double delta3, double oxy_teta, double oxy_cr, const int nl, 
		const int nb_forage, const ivector day_layer, const ivector night_layer, const double DL);

void dfybet_comp(dmatrix& dfybet, dmatrix& dfw, dmatrix& dfd, dmatrix& dfe, dmatrix& dff, const dmatrix d, const dmatrix e, const dmatrix f, unsigned long int pos_map, const int maxn, const int dt);
void dfd1_comp(const char pos, double& dfd1, double& dfsigmam, double& dfsigma, double& dfuinf, const double uinf, const double twodd, const double d);
void dfd2_comp(const char pos, double& dfd2, double& dfsigmam, double& dfsigma, double& dfsigmap, double& dfu, const double u, const double twodd, const double d);
void dfd3_comp(const char pos, double& dfd3, double& dfsigma, double& dfsigmap, double& dfusup, const double usup, const double twodd, const double d);
void dfdh_comp(const char pos, double& dfHinf, double& dfH, double& dfHsup, double dfdh, const double d);
void dfdhx_comp(const char pos, double& dfHinf, double& dfH, double& dfHsup, double dfdh, const double habitat, const double habitat_sup, const double habitat_inf, const double d);
void dfdhy_comp(const char pos, double& dfHinf, double& dfH, double& dfHsup, double dfdh, const double habitat, const double habitat_sup, const double habitat_inf, const double d);
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);

//Some constants here waiting for revision and integration to the parfile
const double Vmax_diff  = 1.25;
const double rc = 0.0005;
const double rho = 0.99;

void CCalpop::Precaldia_Caldia(const PMap& map, VarParamCoupled& param, VarMatrices& mat, dvar_matrix& habitat, dvar_matrix& total_pop, const int sp, const int age, const int t_count, const int jday)
{
	CBord bord;
	dvariable mss_species  = param.dvarsMSS_species[sp];
	dvariable mss_size_slope  = param.dvarsMSS_size_slope[sp];
	dvariable c_diff_fish  = param.dvarsC_diff_fish[sp];
	dvariable sigma_species= param.dvarsSigma_species[sp];

	dvar_matrix Mss_species(map.imin, map.imax, map.jinf, map.jsup); Mss_species = mss_species;
	dvar_matrix Mss_size_slope(map.imin, map.imax, map.jinf, map.jsup); Mss_size_slope = mss_size_slope;
	dvar_matrix C_diff_fish(map.imin, map.imax, map.jinf, map.jsup); C_diff_fish = c_diff_fish;
	dvar_matrix Sigma_species(map.imin, map.imax, map.jinf, map.jsup); Sigma_species = sigma_species;

	dvar_matrix W(0,maxn-1,0,maxn-1);
	W.initialize();

	precaldia_comp(map, param, mat, value(habitat), value(total_pop),  value(mss_species), value(mss_size_slope), value(sigma_species), value(c_diff_fish), sp, age, jday);

	mat.dvarsAdvection_x = nograd_assign(mat.advection_x);
	mat.dvarsAdvection_y = nograd_assign(mat.advection_y);
	mat.dvarsDiffusion_x = nograd_assign(mat.diffusion_x);
	mat.dvarsDiffusion_y = nograd_assign(mat.diffusion_y);
	caldia(map, param, value(mat.dvarsDiffusion_x), value(mat.dvarsAdvection_x), value(mat.dvarsDiffusion_y), value(mat.dvarsAdvection_y));

	dvarsA = nograd_assign(a);
	dvarsB = nograd_assign(b);
	dvarsC = nograd_assign(c);
	dvarsD = nograd_assign(d);
	dvarsE = nograd_assign(e);
	dvarsF = nograd_assign(f);
	Ybet   = nograd_assign(ybet);
	/*if (!param.gcalc()){
		// only in simulation mode: compute mean speed in BL/sec and mean diffusion rate in nmi^2/day
		mat.MeanVarMovement(map,value(mat.dvarsAdvection_x),value(mat.dvarsAdvection_y),
				    value(mat.dvarsDiffusion_y),value(mss_species),value(sigma_species),
				    param.length(sp,age),param.length(sp,param.sp_nb_cohorts[sp]-1),param.deltaT,sp,age);
	}*/
	save_identifier_string2((char*)"Precaldia_Caldia_begin");
	if (param.vert_movement[sp]){
		mat.dvarsU.save_dvar_matrix_position();
		mat.dvarsV.save_dvar_matrix_position();
	}
	habitat.save_dvar_matrix_position();
	dvarsA.save_dvar_matrix_position();
	dvarsB.save_dvar_matrix_position();
	dvarsC.save_dvar_matrix_position();
	dvarsD.save_dvar_matrix_position();
	dvarsE.save_dvar_matrix_position();
	dvarsF.save_dvar_matrix_position();
	Ybet.save_dvar_matrix_position();
	W.save_dvar_matrix_position();
	mat.dvarsAdvection_x.save_dvar_matrix_position();
	mat.dvarsAdvection_y.save_dvar_matrix_position();
	mat.dvarsDiffusion_x.save_dvar_matrix_position();
	mat.dvarsDiffusion_y.save_dvar_matrix_position();
	Mss_species.save_dvar_matrix_position();
	Mss_size_slope.save_dvar_matrix_position();
	Sigma_species.save_dvar_matrix_position();
	mat.dvarSeasonSwitch[sp].save_dvar_matrix_position();
	C_diff_fish.save_dvar_matrix_position();
	c_diff_fish.save_prevariable_value();
	sigma_species.save_prevariable_value();
	mss_species.save_prevariable_value();
	mss_size_slope.save_prevariable_value();
	unsigned long int pmap   = (unsigned long int)&map;
	save_long_int_value(pmap);
	unsigned long int pop    = (unsigned long int)this;
	save_long_int_value(pop);	
	unsigned long int cparam = (unsigned long int)&param;
	save_long_int_value(cparam);
	unsigned long int cbord = (unsigned long int)&bord;
	save_long_int_value(cbord);
	unsigned long int cmat   = (unsigned long int)&mat;
	save_long_int_value(cmat);
	save_int_value(t_count);
	save_int_value(jday);
	save_int_value(sp);
	save_int_value(age);
	save_identifier_string2((char*)"Precaldia_Caldia_end");

	if (!param.vert_movement[sp])
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_caldia);
	else 
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_caldia_UV);
}

void dv_caldia()
{
	verify_identifier_string2((char*)"Precaldia_Caldia_end");
	unsigned age       = restore_int_value();
	unsigned sp        = restore_int_value();
	unsigned jday      = restore_int_value();
	unsigned t_count   = restore_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_bord  = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned long int pos_pop   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	double mss_size_slope= restore_prevariable_value();	
	double mss_species   = restore_prevariable_value();	
	double sigma_species = restore_prevariable_value();	
	double c_diff_fish   = restore_prevariable_value();		
	const dvar_matrix_position cdiff_pos= restore_dvar_matrix_position();
	const dvar_matrix_position Switch_pos= restore_dvar_matrix_position();
	const dvar_matrix_position sigma_pos= restore_dvar_matrix_position();
	const dvar_matrix_position msss_pos = restore_dvar_matrix_position();
	const dvar_matrix_position mss_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position dify_pos = restore_dvar_matrix_position();
	const dvar_matrix_position difx_pos = restore_dvar_matrix_position();
	const dvar_matrix_position advy_pos = restore_dvar_matrix_position();
	const dvar_matrix_position advx_pos = restore_dvar_matrix_position();
	const dvar_matrix_position w_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position ybet_pos = restore_dvar_matrix_position();
	const dvar_matrix_position f_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position e_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position d_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position c_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position b_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position a_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position H_pos    = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"Precaldia_Caldia_begin");

	dmatrix dfH 	 = restore_dvar_matrix_derivatives(H_pos);
	dmatrix dfAdv_x	 = restore_dvar_matrix_derivatives(advx_pos);
	dmatrix dfAdv_y	 = restore_dvar_matrix_derivatives(advy_pos);
	dmatrix dfDiff_x = restore_dvar_matrix_derivatives(difx_pos);
	dmatrix dfDiff_y = restore_dvar_matrix_derivatives(dify_pos);
	dmatrix dfMSS_slope = restore_dvar_matrix_derivatives(msss_pos);
	dmatrix dfMSS 	 = restore_dvar_matrix_derivatives(mss_pos);
	dmatrix dfSwitch = restore_dvar_matrix_derivatives(Switch_pos);
	dmatrix dfSigma  = restore_dvar_matrix_derivatives(sigma_pos);
	dmatrix dfC_diff = restore_dvar_matrix_derivatives(cdiff_pos);
	dmatrix dfa 	 = restore_dvar_matrix_derivatives(a_pos);
	dmatrix dfb 	 = restore_dvar_matrix_derivatives(b_pos);
	dmatrix dfc 	 = restore_dvar_matrix_derivatives(c_pos);
	dmatrix dfd 	 = restore_dvar_matrix_derivatives(d_pos);
	dmatrix dfe 	 = restore_dvar_matrix_derivatives(e_pos);
	dmatrix dff 	 = restore_dvar_matrix_derivatives(f_pos);
	dmatrix dfybet 	 = restore_dvar_matrix_derivatives(ybet_pos);
	dmatrix dfw 	 = restore_dvar_matrix_derivatives(w_pos);

	CMatrices* mat 	 = (CMatrices*) pos_mat;
	PMap* map 	 = (PMap*) pos_map;
	CCalpop* pop 	 = (CCalpop*) pos_pop;
	CParam* param 	 = (CParam*) pos_param;
	CBord* bord 	 = (CBord*) pos_bord;

	const int jinf = map->jmin;
	const int jsup = map->jmax;
	const int iinf = map->imin;
	const int isup = map->imax;
	const int dt   = 2*pop->get_iterationN();
	const int maxn = pop->get_maxn();
	const double Vinf = pop->get_Vinf();
	double dx      = param->deltaX;
	double twodxdx = 2*dx*dx;
	double dy      = param->deltaY;
	double twodydy = 2*dy*dy;

	dmatrix adv_x, adv_y;
	adv_x.allocate(map->imin1, map->imax1, map->jinf1, map->jsup1);
	adv_y.allocate(map->imin1, map->imax1, map->jinf1, map->jsup1);
	adv_x.initialize();
	adv_y.initialize();

	dmatrix d,e,f;

	d.allocate(iinf, isup, map->jinf, map->jsup);
	e.allocate(iinf, isup, map->jinf, map->jsup);
	f.allocate(iinf, isup, map->jinf, map->jsup);
	d.initialize();
	e.initialize();
	f.initialize();

	dmatrix habitat(H_pos);
//	dmatrix habitat;
//	habitat.allocate(map->imin1, map->imax1, map->jinf1, map->jsup1);
	int ind_adult_habitat = param->age_compute_habitat(sp,age);
	//cout << age << " " << ind_adult_habitat << endl;
	habitat = mat->adult_habitat(sp,t_count,ind_adult_habitat);

	pop->Recomp_DEF_coef(*map, *param, *mat, t_count, jday, habitat, d, e, f, adv_x, adv_y, sp, age, mss_species, c_diff_fish, sigma_species);

	//ybet_comp();
	dfybet_comp(dfybet,dfw,dfd,dfe,dff,d,e,f,pos_map,maxn,dt);

	for (int i=isup; i >= iinf; i--){
		for (int j=map->jsup[i]; j>=map->jinf[i] ; j--){
			bord->b	  = map->bord_cell[i][j];
			char pos  = bord->cotey();

			//f[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
			dfd3_comp(pos, dff(i,j), dfDiff_y(i,j), dfDiff_y(i,j+1), dfAdv_y(i,j+1), adv_y(i,j+1), twodydy, dy);

			//e[i][j] = d2(pos, sigmam, sigma, sigmap, v, twodydy, dydy, dy, dt);
			dfd2_comp(pos, dfe(i,j), dfDiff_y(i,j-1), dfDiff_y(i,j), dfDiff_y(i,j+1), dfAdv_y(i,j), adv_y(i,j), twodydy, dy);

			//d[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
			dfd1_comp(pos, dfd(i,j), dfDiff_y(i,j-1), dfDiff_y(i,j), dfAdv_y(i,j-1), adv_y(i,j-1), twodydy, dy);
		}
	}

	for (int j=jsup; j >= jinf; j--){
		for (int i=map->isup[j]; i>=map->iinf[j] ; i--){                            
			bord->b	  = map->bord_cell[i][j];
			char pos  = bord->cotex();

			//c[j][i] = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
			dfd3_comp(pos, dfc(j,i), dfDiff_x(i,j), dfDiff_x(i+1,j), dfAdv_x(i+1,j), adv_x(i+1,j), twodxdx, dx);

			//b[j][i] = d2(pos, sigmam, sigma, sigmap, u, twodxdx, dxdx, dx);
			dfd2_comp(pos, dfb(j,i), dfDiff_x(i-1,j), dfDiff_x(i,j), dfDiff_x(i+1,j), dfAdv_x(i,j), adv_x(i,j), twodxdx, dx);

			//a[j][i] = d1(pos, sigmam, sigma, uinf, twodxdx, dx);
			dfd1_comp(pos, dfa(j,i), dfDiff_x(i-1,j), dfDiff_x(i,j), dfAdv_x(i-1,j), adv_x(i-1,j), twodxdx, dx);
		}
	}

	dvector lat_correction(jinf,jsup);
	lat_correction = mat->lat_correction;

	const int    deltaT = param->deltaT;
	const double length = param->length[sp][age]*0.01;
	const double lmax   = param->length[sp][param->sp_nb_cohorts[sp]-1]*0.01;
	const double unit_x = pow(length,mss_size_slope)*(3600*24.0*deltaT/1852)*dx;
	const double unit_y = pow(length,mss_size_slope)*(3600*24.0*deltaT/1852)*dy;
	const double Dspeed = Vmax_diff-0.25*length/lmax;
	const double Dinf   = pow(Dspeed*length*3600*24.0*deltaT/1852,2)/(4.0*deltaT);
	const double Dmax   = sigma_species*Dinf;

	const int imax = map->imax;
	const int imin = map->imin;
	
	const double nb_layer = param->nb_layer;

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nl = map->carte(i,j);
			if (nl>0) {
				bord->b	   = map->bord_cell[i][j];
				//bord.b   = map->nbl_bord_cell[i][j];
				char pos_x = bord->cotex();
				char pos_y = bord->cotey();
			
				//first recompute dHdx and dHdy
				double dHdx = 0.0;
				double dHdy = 0.0;

                                //double dHx_right = (habitat[i+1][j] - habitat[i][j])/dx;
                                //double dHx_left  = (habitat[i][j] - habitat[i-1][j])/dx;

                                //double dHy_right = (habitat[i][j+1] - habitat[i][j])/dy;
                                //double dHy_left  = (habitat[i][j] - habitat[i][j-1])/dy;


				if (nl >= nb_layer){
//	                    		dHdx = (habitat[i+1][j] - habitat[i][j])/dx;
//					dHdy = (habitat[i][j+1] - habitat[i][j])/dy;

					if (pos_x == SANS)
 	                   			dHdx = (habitat[i+1][j] - habitat[i-1][j])/(2*dx);
//dHdx = max(dHx_right,dHx_left); 	                   			
					else if (pos_x == D_FERME)
						dHdx = (habitat[i][j] - habitat[i-1][j])/(dx);
					else if (pos_x == G_FERME)
						dHdx = (habitat[i+1][j] - habitat[i][j])/(dx);
					
					if (pos_y == SANS)
						dHdy = (habitat[i][j+1] - habitat[i][j-1])/(2*dy);
//dHdy = min(dHy_right,dHy_left);
					else if (pos_y == D_FERME)
						dHdy = (habitat[i][j] - habitat[i][j-1])/(dy);
					else if (pos_y == G_FERME)
						dHdy = (habitat[i][j+1] - habitat[i][j])/(dy);

				}

				double dfD    = 0.0;
				double dfdhdx = 0.0;
				double dfdhdy = 0.0;
				double dfrho_x = 0.0;
				double dfrho_y = 0.0;
/*
				//advection_y(i,j) = mat.v[i][j] + MSS * unit_y * dHdy;
				dfMSS(i,j) += unit_y * dHdy * dfAdv_y(i,j);
				dfdhdy     += mss_species * unit_y * dfAdv_y(i,j);
				dfAdv_y(i,j) = 0.0;

				//const double length = param->length[sp][age]*0.01;
				//advection_x(i,j) = mat.u[i][j] + MSS * unit_x * dHdx * mat.lat_correction[j];
				dfMSS(i,j) += unit_x * dHdx  * lat_correction(j) * dfAdv_x(i,j);
				dfdhdx     += mss_species * unit_x * lat_correction(j) * dfAdv_x(i,j);
				dfAdv_x(i,j) = 0.0;
*/

				//precompute v_x and v_y
				double v_x = mss_species * unit_x * dHdx/Vinf;				
				double v_y = mss_species * unit_y * dHdy/Vinf;
			
				double dfv_x = 0.0;
				double dfv_y = 0.0;

				//advection_y(i,j) = mat.v[i][j] + v_y;
				dfv_y += dfAdv_y(i,j);
				dfAdv_y(i,j) = 0.0;

				//advection_x(i,j) = mat.u[i][j] + v_x * mat.lat_correction[j];
				dfv_x += lat_correction(j) * dfAdv_x(i,j);
				dfAdv_x(i,j) = 0.0;


				//v_y = Vinf*param.func_limit_one(v_y/Vinf);
				if (v_y>=0) dfv_y = param->dffunc_limit_one(v_y,Vinf*dfv_y)/Vinf;
				if (v_y<0)  dfv_y = param->dffunc_limit_one(-v_y,Vinf*dfv_y)/Vinf;

				//v_x = Vinf*param.func_limit_one(v_x/Vinf);
				if (v_x>=0) dfv_x = param->dffunc_limit_one(v_x,Vinf*dfv_x)/Vinf;
				if (v_x<0)  dfv_x = param->dffunc_limit_one(-v_x,Vinf*dfv_x)/Vinf;

				//double v_y = MSS * unit_y * dHdy; 
				dfMSS(i,j) += unit_y * dHdy * dfv_y;
				dfMSS_slope(i,j) += mss_species * unit_y * log(length) * dHdy * dfv_y;
				dfdhdy     += mss_species * unit_y * dfv_y;
				
				//double v_x = MSS * unit_x * dHdx; 
				dfMSS(i,j) += unit_x * dHdx  * dfv_x;
				dfMSS_slope(i,j) += mss_species * unit_x * log(length) * dHdx  * dfv_x;
				dfdhdx     += mss_species * unit_x * dfv_x;

				double rho_x = 1.0 - rho * sqrt(dHdx*dHdx) * dx;
				double rho_y = 1.0 - rho * sqrt(dHdy*dHdy) * dy;

				//double diff_habitat = 1 - habitat(i,j)/(c_diff_fish + habitat(i,j));
				double diff_habitat = 1.0 - c_diff_fish*pow(habitat(i,j),3);
				//double diff_habitat = 1.0 - D_dec_H1*pow(habitat(i,j),c_diff_fish);
				double D = Dmax * diff_habitat;
//vary diffusion with seasons


				double sfunc = mat->season_switch(sp,jday,j);
 				D = (0.9*D*sfunc + D*(1.0-sfunc));

				//diffusion_y(i,j) = rho_y * D;
				dfrho_y += D * dfDiff_y(i,j);
				dfD     += rho_y * dfDiff_y(i,j);
				dfDiff_y(i,j) = 0.0;

				//diffusion_x(i,j) = rho_x * D * lat_correction;
				dfrho_x += D * lat_correction(j) * dfDiff_x(i,j);
				dfD     += rho_x * lat_correction(j) * dfDiff_x(i,j);
				dfDiff_x(i,j) = 0.0;

 				//D = (0.9*D*sfunc + D*(1.0-sfunc));
 				double dfD_pr  = (1.0-0.1*sfunc)*dfD;
				dfSwitch(i,j) -= 0.1*D*dfD;
				dfD = 0.0;

				//differentiate if rho depends on H!!!
				//rho_x = 1.0 - 0.9 * sqrt(dHdx*dHdx) * dx; 			
				if (dHdx > 0)
					dfdhdx -= rho*dx * dfrho_x;
				else if (dHdx < 0)
					dfdhdx += rho*dx * dfrho_x;
				
				//rho_y = 1.0 - 0.9 * sqrt(dHdy*dHdy) * dy;
				if (dHdy > 0)
					dfdhdy -= rho*dy * dfrho_y;
				else if (dHdy < 0)
					dfdhdy += rho*dy * dfrho_y;

				if (nl >= nb_layer){
					//adjoint for computing dHdx and dHdy
					//dfdh_comp(G_FERME,dfH(i,j-1),dfH(i,j),dfH(i,j+1),dfdhdy,dy);
					//dfdh_comp(G_FERME,dfH(i-1,j),dfH(i,j),dfH(i+1,j),dfdhdx,dx);	
					dfdh_comp(pos_y,dfH(i,j-1),dfH(i,j),dfH(i,j+1),dfdhdy,dy);
					dfdh_comp(pos_x,dfH(i-1,j),dfH(i,j),dfH(i+1,j),dfdhdx,dx);	
					//dfdhy_comp(pos_y,dfH(i,j-1),dfH(i,j),dfH(i,j+1),dfdhdy,habitat(i,j),habitat(i,j+1),habitat(i,j-1),dy);
					//dfdhx_comp(pos_x,dfH(i-1,j),dfH(i,j),dfH(i+1,j),dfdhdx,habitat(i,j),habitat(i+1,j),habitat(i-1,j),dx);	
				}

/*
				//double D = Dmax*(1-habitat(i,j)/(c_diff_fish+habitat(i,j)));
				dfH(i,j)     -= Dmax * c_diff_fish/pow(c_diff_fish+habitat(i,j),2) * dfD;
				//dfSigma(i,j) += 1000 * (1-habitat(i,j)/(c_diff_fish+habitat(i,j))) * dfD;
				dfSigma(i,j) += Dinf * (1-habitat(i,j)/(c_diff_fish+habitat(i,j))) * dfD;
				dfC_diff(i,j)+= Dmax * habitat(i,j)/pow(c_diff_fish+habitat(i,j),2) * dfD;
				dfD = 0.0;
*/
				//double D = Dmax*(1-c_diff_fish*pow(habitat(i,j),3));
                                dfH(i,j)     -= Dmax * 3.0 * c_diff_fish*pow(habitat(i,j),2) * dfD_pr;
                                dfSigma(i,j) += Dinf * (1-c_diff_fish*pow(habitat(i,j),3)) * dfD_pr;
                                dfC_diff(i,j)-= Dmax * pow(habitat(i,j),3) * dfD_pr;
                                dfD_pr = 0.0;

/*  				//double D = Dmax*(1-D_dec_H1*pow(habitat(i,j),c_diff_fish));
                                dfH(i,j)     -= Dmax * c_diff_fish * D_dec_H1 * pow(habitat(i,j),c_diff_fish-1) * dfD_pr;
                                dfSigma(i,j) += Dinf * diff_habitat * dfD_pr;
                                dfC_diff(i,j)-= Dmax * D_dec_H1 * pow(habitat(i,j),c_diff_fish) * log(habitat(i,j)) * dfD_pr;
                                dfD_pr = 0.0;
*/

/*				//double D = Dmax*(1-c_diff_fish*pow(habitat(i,j),3));
                                dfH(i,j)     -= Dmax * 3.0 * c_diff_fish*pow(habitat(i,j),2) * dfD;
                                dfSigma(i,j) += Dinf * (1-c_diff_fish*pow(habitat(i,j),3)) * dfD;
                                dfC_diff(i,j)-= Dmax * pow(habitat(i,j),3) * dfD;
				dfD = 0.0;
*/
			}
		}
	}
	dfa.save_dmatrix_derivatives(a_pos);
	dfb.save_dmatrix_derivatives(b_pos);
	dfc.save_dmatrix_derivatives(c_pos);
	dfd.save_dmatrix_derivatives(d_pos);
	dfe.save_dmatrix_derivatives(e_pos);
	dff.save_dmatrix_derivatives(f_pos);
	dfybet.save_dmatrix_derivatives(ybet_pos); 
	dfw.save_dmatrix_derivatives(w_pos);
	dfH.save_dmatrix_derivatives(H_pos);
	dfAdv_x.save_dmatrix_derivatives(advx_pos);
	dfAdv_y.save_dmatrix_derivatives(advy_pos);
	dfDiff_x.save_dmatrix_derivatives(difx_pos);
	dfDiff_y.save_dmatrix_derivatives(dify_pos);
	dfMSS_slope.save_dmatrix_derivatives(msss_pos);
	dfMSS.save_dmatrix_derivatives(mss_pos);
	dfSwitch.save_dmatrix_derivatives(Switch_pos);
	dfSigma.save_dmatrix_derivatives(sigma_pos);
	dfC_diff.save_dmatrix_derivatives(cdiff_pos);
}

void dv_caldia_UV()
{
	verify_identifier_string2((char*)"Precaldia_Caldia_end");
	unsigned age       = restore_int_value();
	unsigned sp        = restore_int_value();
	unsigned jday      = restore_int_value();
	unsigned t_count   = restore_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_bord  = restore_long_int_value();
	unsigned long int pos_param = restore_long_int_value();
	unsigned long int pos_pop   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	double mss_size_slope   = restore_prevariable_value();	
	double mss_species   = restore_prevariable_value();	
	double sigma_species = restore_prevariable_value();	
	double c_diff_fish   = restore_prevariable_value();		
	const dvar_matrix_position cdiff_pos= restore_dvar_matrix_position();
	const dvar_matrix_position Switch_pos= restore_dvar_matrix_position();
	const dvar_matrix_position sigma_pos= restore_dvar_matrix_position();
	const dvar_matrix_position msss_pos = restore_dvar_matrix_position();
	const dvar_matrix_position mss_pos  = restore_dvar_matrix_position();
	const dvar_matrix_position dify_pos = restore_dvar_matrix_position();
	const dvar_matrix_position difx_pos = restore_dvar_matrix_position();
	const dvar_matrix_position advy_pos = restore_dvar_matrix_position();
	const dvar_matrix_position advx_pos = restore_dvar_matrix_position();
	const dvar_matrix_position w_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position ybet_pos = restore_dvar_matrix_position();
	const dvar_matrix_position f_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position e_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position d_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position c_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position b_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position a_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position H_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position V_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position U_pos    = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"Precaldia_Caldia_begin");

	dmatrix dfH 	 = restore_dvar_matrix_derivatives(H_pos);
	dmatrix dfU	 = restore_dvar_matrix_derivatives(U_pos);
	dmatrix dfV	 = restore_dvar_matrix_derivatives(V_pos);
	dmatrix dfAdv_x	 = restore_dvar_matrix_derivatives(advx_pos);
	dmatrix dfAdv_y	 = restore_dvar_matrix_derivatives(advy_pos);
	dmatrix dfDiff_x = restore_dvar_matrix_derivatives(difx_pos);
	dmatrix dfDiff_y = restore_dvar_matrix_derivatives(dify_pos);
	dmatrix dfMSS_slope = restore_dvar_matrix_derivatives(msss_pos);
	dmatrix dfMSS 	 = restore_dvar_matrix_derivatives(mss_pos);
	dmatrix dfSwitch  = restore_dvar_matrix_derivatives(Switch_pos);
	dmatrix dfSigma  = restore_dvar_matrix_derivatives(sigma_pos);
	dmatrix dfC_diff = restore_dvar_matrix_derivatives(cdiff_pos);
	dmatrix dfa 	 = restore_dvar_matrix_derivatives(a_pos);
	dmatrix dfb 	 = restore_dvar_matrix_derivatives(b_pos);
	dmatrix dfc 	 = restore_dvar_matrix_derivatives(c_pos);
	dmatrix dfd 	 = restore_dvar_matrix_derivatives(d_pos);
	dmatrix dfe 	 = restore_dvar_matrix_derivatives(e_pos);
	dmatrix dff 	 = restore_dvar_matrix_derivatives(f_pos);
	dmatrix dfybet 	 = restore_dvar_matrix_derivatives(ybet_pos);
	dmatrix dfw 	 = restore_dvar_matrix_derivatives(w_pos);

	PMap* map 	 = (PMap*) pos_map;
	CMatrices* mat 	 = (CMatrices*) pos_mat;
	CCalpop* pop 	 = (CCalpop*) pos_pop;
	CParam* param 	 = (CParam*) pos_param;
	CBord* bord 	 = (CBord*) pos_bord;

	const int jinf = map->jmin;
	const int jsup = map->jmax;
	const int iinf = map->imin;
	const int isup = map->imax;
	const int dt   = 2*pop->get_iterationN();
	const int maxn = pop->get_maxn();
	const double Vinf = pop->get_Vinf();
	double dx      = param->deltaX;
	double twodxdx = 2*dx*dx;
	double dy      = param->deltaY;
	double twodydy = 2*dy*dy;

	//Recompute u and v matrices
	dmatrix u, v;
	u.allocate(iinf, isup, map->jinf, map->jsup);
	v.allocate(iinf, isup, map->jinf, map->jsup);
	u.initialize(); v.initialize();
	//Parameters to recompute accessibility to layers 
	const double oxy_teta = param->a_oxy_habitat[sp];
	const double oxy_cr   = param->b_oxy_habitat[sp];
	const double sigma_ha = param->sigma_ha[sp][age];
	const double temp_age = param->temp_age[sp][age];
	const double twosigsq = 2.0*sigma_ha*sigma_ha;
	const double temp_max = param->b_sst_spawning(sp);
	const double delta1   = param->thermal_func_delta[0][sp];
	const double delta2   = param->thermal_func_delta[1][sp];
	const double delta3   = param->thermal_func_delta[2][sp];
	
	const int Tfunc_Gaussian = param->gaussian_thermal_function[sp];
	
	const int nb_forage = param->get_nbforage();
	ivector day_layer(0,nb_forage-1); day_layer = param->day_layer;
	ivector night_layer(0,nb_forage-1); night_layer = param->night_layer;
	
	const int nb_layer = param->nb_layer;
	dvector lf_access(0,nb_layer-1);
	dvector l_access(0,nb_layer-1);
	dvector F(0,nb_forage-1);
	dvector O2(0,nb_forage-1); 
	dvector T(0,nb_forage-1);  
	//end of accessiblity parameters section

	for (int i = iinf; i <= isup; i++){	
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmin; j <= jmax; j++){
			const int nl = map->carte(i,j);
			if (nl>0){
				l_access.initialize();
				lf_access.initialize();
				O2.initialize();
				T.initialize();
				for (int n=0; n<nb_forage; n++)
					F(n) = mat->forage(t_count,n,i,j);
				//for (int l=0; l<nl; l++){
				for (int l=0; l<nb_layer; l++){
					O2(l) = mat->oxygen(t_count,l,i,j);
					T(l) = mat->tempn(t_count,l,i,j);
				}
				const double DL = mat->daylength(jday,j)/24.0;

				if (Tfunc_Gaussian){
					f_accessibility(l_access,lf_access,F,O2,T,twosigsq,temp_age,oxy_teta,oxy_cr,
							nb_layer,nb_forage,day_layer,night_layer,DL);	
				} else {
					f_accessibility(l_access,lf_access,F,O2,T,temp_age,temp_max,delta1,delta2,delta3,oxy_teta,oxy_cr,
							nb_layer,nb_forage,day_layer,night_layer,DL);
				}
				for (int l=0; l<nl; l++){
					u(i,j) += mat->un(t_count,l,i,j) * lf_access(l);
					v(i,j) += mat->vn(t_count,l,i,j) * lf_access(l);
				}	
			}
		}
	} //end of recalculation section
//cout << "2: " << t_count << " " << age << " "<<  l_access << " " << lf_access << " " << temp_age << " " << temp_max << " " << delta1 << " " << delta2 << " " << delta3 << " " << oxy_teta << " " << oxy_cr << " " << nb_layer << endl;			

	dmatrix adv_x, adv_y;
	adv_x.allocate(map->imin1, map->imax1, map->jinf1, map->jsup1);
	adv_y.allocate(map->imin1, map->imax1, map->jinf1, map->jsup1);
	adv_x.initialize();
	adv_y.initialize();

	dmatrix d,e,f;

	d.allocate(iinf, isup, map->jinf, map->jsup);
	e.allocate(iinf, isup, map->jinf, map->jsup);
	f.allocate(iinf, isup, map->jinf, map->jsup);
	d.initialize();
	e.initialize();
	f.initialize();

	dmatrix habitat(H_pos);
	//dmatrix habitat;
	//habitat.allocate(map->imin1, map->imax1, map->jinf1, map->jsup1);
	int ind_adult_habitat = param->age_compute_habitat(sp,age);
	//cout << age << " " << ind_adult_habitat << endl;
	habitat = mat->adult_habitat(sp,t_count,ind_adult_habitat);

	pop->Recomp_DEF_UV_coef(*map, *param, *mat, u, v, habitat, d, e, f, adv_x, adv_y, sp, age, mss_species, c_diff_fish, sigma_species,jday);

	//ybet_comp();
	dfybet_comp(dfybet,dfw,dfd,dfe,dff,d,e,f,pos_map,maxn,dt);

	for (int i=isup; i >= iinf; i--){
		for (int j=map->jsup[i]; j>=map->jinf[i] ; j--){
			bord->b	  = map->bord_cell[i][j];
			char pos  = bord->cotey();

			//f[i][j] = d3(pos, sigma, sigmap, vsup, twodydy, dydy, dy);
			dfd3_comp(pos, dff(i,j), dfDiff_y(i,j), dfDiff_y(i,j+1), dfAdv_y(i,j+1), adv_y(i,j+1), twodydy, dy);

			//e[i][j] = d2(pos, sigmam, sigma, sigmap, v, twodydy, dydy, dy, dt);
			dfd2_comp(pos, dfe(i,j), dfDiff_y(i,j-1), dfDiff_y(i,j), dfDiff_y(i,j+1), dfAdv_y(i,j), adv_y(i,j), twodydy, dy);

			//d[i][j] = d1(pos, sigmam, sigma, vinf, twodydy, dydy, dy);
			dfd1_comp(pos, dfd(i,j), dfDiff_y(i,j-1), dfDiff_y(i,j), dfAdv_y(i,j-1), adv_y(i,j-1), twodydy, dy);
		}
	}

	for (int j=jsup; j >= jinf; j--){
		for (int i=map->isup[j]; i>=map->iinf[j] ; i--){                            
			bord->b	  = map->bord_cell[i][j];
			char pos  = bord->cotex();

			//c[j][i] = d3(pos, sigma, sigmap, usup, twodxdx, dxdx, dx);
			dfd3_comp(pos, dfc(j,i), dfDiff_x(i,j), dfDiff_x(i+1,j), dfAdv_x(i+1,j), adv_x(i+1,j), twodxdx, dx);

			//b[j][i] = d2(pos, sigmam, sigma, sigmap, u, twodxdx, dxdx, dx);
			dfd2_comp(pos, dfb(j,i), dfDiff_x(i-1,j), dfDiff_x(i,j), dfDiff_x(i+1,j), dfAdv_x(i,j), adv_x(i,j), twodxdx, dx);

			//a[j][i] = d1(pos, sigmam, sigma, uinf, twodxdx, dx);
			dfd1_comp(pos, dfa(j,i), dfDiff_x(i-1,j), dfDiff_x(i,j), dfAdv_x(i-1,j), adv_x(i-1,j), twodxdx, dx);
		}
	}

	dvector lat_correction(jinf,jsup);
	lat_correction = mat->lat_correction;

	const int    deltaT = param->deltaT;
	const double length = param->length[sp][age]*0.01;
	const double lmax   = param->length[sp][param->sp_nb_cohorts[sp]-1]*0.01;
	const double unit_x = pow(length,mss_size_slope)*(3600*24.0*deltaT/1852)*dx;
	const double unit_y = pow(length,mss_size_slope)*(3600*24.0*deltaT/1852)*dy;
	const double Dspeed = Vmax_diff-0.25*length/lmax;
	const double Dinf   = pow(Dspeed*length*3600*24.0*deltaT/1852,2)/(4.0*deltaT);
	const double Dmax   = sigma_species*Dinf;
	const double rmax   = param->rmax_currents;

	const int imax = map->imax;
	const int imin = map->imin;

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			const int nl = map->carte(i,j);
			if (nl>0) {

				bord->b     = map->bord_cell[i][j];
				//bord.b   = map->nbl_bord_cell[i][j];
				char pos_x = bord->cotex();
				char pos_y = bord->cotey();
				
				//first recompute dHdx and dHdy
				double dHdx = 0.0;
				double dHdy = 0.0;

                                //double dHx_right = (habitat[i+1][j] - habitat[i][j])/dx;
                                //double dHx_left  = (habitat[i][j] - habitat[i-1][j])/dx;

                                //double dHy_right = (habitat[i][j+1] - habitat[i][j])/dy;
                                //double dHy_left  = (habitat[i][j] - habitat[i][j-1])/dy;


				if (nl >= nb_layer){
//	                    		dHdx = (habitat[i+1][j] - habitat[i][j])/dx;
//					dHdy = (habitat[i][j+1] - habitat[i][j])/dy;

					if (pos_x == SANS)
	                    			dHdx = (habitat[i+1][j] - habitat[i-1][j])/(2*dx);
//dHdx = max(dHx_right,dHx_left);	                    			
					else if (pos_x == D_FERME)
						dHdx = (habitat[i][j] - habitat[i-1][j])/(dx);
					else if (pos_x == G_FERME)
						dHdx = (habitat[i+1][j] - habitat[i][j])/(dx);
					
					if (pos_y == SANS)
						dHdy = (habitat[i][j+1] - habitat[i][j-1])/(2*dy);
//dHdy = min(dHy_right,dHy_left);						
					else if (pos_y == D_FERME)
						dHdy = (habitat[i][j] - habitat[i][j-1])/(dy);
					else if (pos_y == G_FERME)
						dHdy = (habitat[i][j+1] - habitat[i][j])/(dy);

				}

				double dfD    = 0.0;
				double dfD_season    = 0.0;
				double dfdhdx = 0.0;
				double dfdhdy = 0.0;
				double dfrho_x = 0.0;
				double dfrho_y = 0.0;
		
				//recompute c coefficient 
				const double expr1 = (u(i,j)*u(i,j)+v(i,j)*v(i,j));
				const double expr2 = sqrt(expr1);
				const double r = rmax/(1+rc*expr2);
				const double c = 1.0-r*length/lmax;
				double dfc = 0.0;
				double dfr = 0.0;

/*
				//advection_y(i,j) = c*mat.v[i][j] + MSS * unit_y * dHdy;
				dfV(i,j)   += c * dfAdv_y(i,j);
				dfc 	   += v(i,j) * dfAdv_y(i,j);
				dfMSS(i,j) += unit_y * dHdy * dfAdv_y(i,j);
				dfdhdy     += mss_species * unit_y * dfAdv_y(i,j);
				dfAdv_y(i,j) = 0.0;

				//advection_x(i,j) = c*mat.u[i][j] + MSS * unit_x * dHdx * mat.lat_correction[j];
				dfU(i,j)   += c * dfAdv_x(i,j);
				dfc 	   += u(i,j) * dfAdv_x(i,j);
				dfMSS(i,j) += unit_x * dHdx  * lat_correction(j) * dfAdv_x(i,j);
				dfdhdx     += mss_species * unit_x * lat_correction(j) * dfAdv_x(i,j);
				dfAdv_x(i,j) = 0.0;
*/
				//precompute v_x and v_y
				double v_x = mss_species * unit_x * dHdx / Vinf;				
				double v_y = mss_species * unit_y * dHdy / Vinf;

				double dfv_x = 0.0;
				double dfv_y = 0.0;

				//advection_y(i,j) = c*mat.v[i][j] + v_y;
				dfV(i,j)+= c * dfAdv_y(i,j);
 				dfc 	+= v(i,j) * dfAdv_y(i,j);
				dfv_y 	+= dfAdv_y(i,j);
				dfAdv_y(i,j) = 0.0;

				//advection_x(i,j) = c*mat.u[i][j] + v_x * mat.lat_correction[j];
				dfU(i,j)+= c * dfAdv_x(i,j);
				dfc 	+= u(i,j) * dfAdv_x(i,j);
				dfv_x   += lat_correction(j) * dfAdv_x(i,j);
				dfAdv_x(i,j) = 0.0;

				//v_y = Vinf*param.func_limit_one(v_y/Vinf);
				if (v_y>=0) dfv_y = param->dffunc_limit_one(v_y,Vinf*dfv_y)/Vinf;
				if (v_y<0)  dfv_y = param->dffunc_limit_one(-v_y,Vinf*dfv_y)/Vinf;
	
				//v_x = Vinf*param.func_limit_one(v_x/Vinf);
				if (v_x>=0) dfv_x = param->dffunc_limit_one(v_x,Vinf*dfv_x)/Vinf;
				if (v_x<0)  dfv_x = param->dffunc_limit_one(-v_x,Vinf*dfv_x)/Vinf;

				//double v_y = MSS * unit_y * dHdy; 
				dfMSS(i,j) += unit_y * dHdy * dfv_y;
				dfMSS_slope(i,j) += mss_species * unit_y * log(length) * dHdy  * dfv_y;
				dfdhdy     += mss_species * unit_y * dfv_y;
				
				//double v_x = MSS * unit_x * dHdx; 
				dfMSS(i,j) += unit_x * dHdx  * dfv_x;
				dfMSS_slope(i,j) += mss_species * unit_x * log(length) * dHdx  * dfv_x;
				dfdhdx     += mss_species * unit_x * dfv_x;

				//recompute fV
				//double fV = 1.0 - expr2/(500.0+expr2);
				double fV = 1.0 - expr2/(500.0*deltaT/30.0+expr2);

				//recompute rho with correction through passive transport
				double rho_x = fV*(1.0 - rho * sqrt(dHdx*dHdx) * dx);
				double rho_y = fV*(1.0 - rho * sqrt(dHdy*dHdy) * dy);

				//finally, contribution to derivatives dfU and dfV through c and r
				//const double c = 1.0-r*length/lmax;
				dfr -= dfc * length/lmax;
				
				//const double r = a/(1+b*sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j)));
				double expr3 = 0.0;
				if (expr2>0)
					expr3 = rmax*rc/(expr2*pow(1+rc*expr2,2.0));
				dfU(i,j) -= u(i,j) * expr3 * dfr;
				dfV(i,j) -= v(i,j) * expr3 * dfr;

				//double diff_habitat = 1 - habitat(i,j)/(c_diff_fish + habitat(i,j));
				double diff_habitat = 1.0 - c_diff_fish*pow(habitat(i,j),3);
				//double diff_habitat = 1.0 - D_dec_H1*pow(habitat(i,j),c_diff_fish);
				double D = Dmax * diff_habitat;

//vary diffusion with seasons
				double sfunc = mat->season_switch(sp,jday,j);
 				double D_season = (0.9*D*sfunc + D*(1.0-sfunc));

				//diffusion_y(i,j) = rho_y * D_season;
				dfrho_y    += D_season * dfDiff_y(i,j);
				dfD_season += rho_y * dfDiff_y(i,j);
				dfDiff_y(i,j) = 0.0;

				//diffusion_x(i,j) = rho_x * D_season * lat_correction;
				dfrho_x    += D_season * lat_correction(j) * dfDiff_x(i,j);
				dfD_season += rho_x * lat_correction(j) * dfDiff_x(i,j);
				dfDiff_x(i,j) = 0.0;

 				//D_season = (0.9*D*sfunc + D*(1.0-sfunc));
				dfSwitch(i,j)-= 0.1 *D * dfD_season;
 				dfD          += (1.0-0.1*sfunc) * dfD_season;
				dfD_season    = 0.0;

	
				//double rho_x = fV*(1.0 - rho * sqrt(dHdx*dHdx) * dx);
				double dffV = (1.0 - rho * sqrt(dHdx*dHdx) * dx) * dfrho_x;

				//double rho_y = fV*(1.0 - rho * sqrt(dHdy*dHdy) * dy); 
				dffV += (1.0 - rho * sqrt(dHdy*dHdy) * dy) * dfrho_y;
				
				//contribution to derivatives dfU and dfV through fV in rho
				//double fV = 1.0 - expr2/(600.0+expr2);
				//double dfC = -500.0*dffV/pow(500.0+expr2,2);
				double dfC = -(500.0*deltaT/30.0)*dffV/pow(500.0*deltaT/30.0+expr2,2);
 				//expr2 = sqrt(U*U+V*V);
				if (expr2>0){
					dfU(i,j) += u(i,j)*dfC/expr2; 
					dfV(i,j) += v(i,j)*dfC/expr2; 
				}
				//differentiate if rho depends on H!!!
				//rho_x = 1.0 - 0.9 * sqrt(dHdx*dHdx) * dx; 			
				if (dHdx > 0)
					//dfdhdx -= rho*dx * dfrho_x;
					dfdhdx -= rho*dx*fV * dfrho_x;
				else if (dHdx < 0)
					//dfdhdx += rho*dx * dfrho_x;
					dfdhdx += rho*dx*fV * dfrho_x;
				
				//rho_y = 1.0 - 0.9 * sqrt(dHdy*dHdy) * dy;
				if (dHdy > 0)
					//dfdhdy -= rho*dy * dfrho_y;
					dfdhdy -= rho*dy*fV * dfrho_y;
				else if (dHdy < 0)
					//dfdhdy += rho*dy * dfrho_y;
					dfdhdy += rho*dy*fV * dfrho_y;

				if (nl >= nb_layer){
					//adjoint for computing dHdx and dHdy
					//dfdh_comp(G_FERME,dfH(i,j-1),dfH(i,j),dfH(i,j+1),dfdhdy,dy);
					//dfdh_comp(G_FERME,dfH(i-1,j),dfH(i,j),dfH(i+1,j),dfdhdx,dx);	
					dfdh_comp(pos_y,dfH(i,j-1),dfH(i,j),dfH(i,j+1),dfdhdy,dy);
					dfdh_comp(pos_x,dfH(i-1,j),dfH(i,j),dfH(i+1,j),dfdhdx,dx);	
					//dfdhy_comp(pos_y,dfH(i,j-1),dfH(i,j),dfH(i,j+1),dfdhdy,habitat(i,j),habitat(i,j+1),habitat(i,j-1),dy);
					//dfdhx_comp(pos_x,dfH(i-1,j),dfH(i,j),dfH(i+1,j),dfdhdx,habitat(i,j),habitat(i+1,j),habitat(i-1,j),dx);
				}
/*
				//double D = Dmax*(1-habitat(i,j)/(c_diff_fish+habitat(i,j)));
				dfH(i,j)     -= Dmax * c_diff_fish/pow(c_diff_fish+habitat(i,j),2) * dfD;
				//dfSigma(i,j) += 1000 * (1-habitat(i,j)/(c_diff_fish+habitat(i,j))) * dfD;
				dfSigma(i,j) += Dinf * (1-habitat(i,j)/(c_diff_fish+habitat(i,j))) * dfD;
				dfC_diff(i,j)+= Dmax * habitat(i,j)/pow(c_diff_fish+habitat(i,j),2) * dfD;
				dfD = 0.0;
*/
/*
                                //double D = Dmax*(1-c_diff_fish*pow(habitat(i,j),3));
                                dfH(i,j)     -= Dmax * 3.0 * c_diff_fish*pow(habitat(i,j),2) * dfD_pr;
                                dfSigma(i,j) += Dinf * (1-c_diff_fish*pow(habitat(i,j),3)) * dfD_pr;
                                dfC_diff(i,j)-= Dmax * pow(habitat(i,j),3) * dfD_pr;
                                dfD_pr = 0.0;
*/

                                //double D = Dmax*(1-c_diff_fish*pow(habitat(i,j),3));
                                dfH(i,j)     -= Dmax * 3.0 * c_diff_fish*pow(habitat(i,j),2) * dfD;
                                dfSigma(i,j) += Dinf * (1-c_diff_fish*pow(habitat(i,j),3)) * dfD;
                                dfC_diff(i,j)-= Dmax * pow(habitat(i,j),3) * dfD;
                                dfD = 0.0;

/*  				//double D = Dmax*(1-D_dec_H1*pow(habitat(i,j),c_diff_fish));
                                dfH(i,j)     -= Dmax * c_diff_fish * D_dec_H1 * pow(habitat(i,j),c_diff_fish-1) * dfD;
                                dfSigma(i,j) += Dinf * diff_habitat * dfD;
                                dfC_diff(i,j)-= Dmax * D_dec_H1 * pow(habitat(i,j),c_diff_fish) * log(habitat(i,j)) * dfD;
                                dfD = 0.0;			
*/
			}
		}
	}
//cout << norm(dfH) << " " << norm(dfU) << " " << norm(dfV) << " "<< norm(dfd)<< " " << norm(dfe) << " "<< norm(dff)<< endl;
	dfa.save_dmatrix_derivatives(a_pos);
	dfb.save_dmatrix_derivatives(b_pos);
	dfc.save_dmatrix_derivatives(c_pos);
	dfd.save_dmatrix_derivatives(d_pos);
	dfe.save_dmatrix_derivatives(e_pos);
	dff.save_dmatrix_derivatives(f_pos);
	dfybet.save_dmatrix_derivatives(ybet_pos); 
	dfw.save_dmatrix_derivatives(w_pos);
	dfH.save_dmatrix_derivatives(H_pos);
	dfU.save_dmatrix_derivatives(U_pos);
	dfV.save_dmatrix_derivatives(V_pos);
	dfAdv_x.save_dmatrix_derivatives(advx_pos);
	dfAdv_y.save_dmatrix_derivatives(advy_pos);
//cout << __FILE__ << " "<< __LINE__ << endl;
	dfDiff_x.save_dmatrix_derivatives(difx_pos);
	dfDiff_y.save_dmatrix_derivatives(dify_pos);
	dfMSS_slope.save_dmatrix_derivatives(msss_pos);
	dfMSS.save_dmatrix_derivatives(mss_pos);
	dfSigma.save_dmatrix_derivatives(sigma_pos);
	dfSwitch.save_dmatrix_derivatives(Switch_pos);
	dfC_diff.save_dmatrix_derivatives(cdiff_pos);
}

void dfybet_comp(dmatrix& dfybet, dmatrix& dfw, dmatrix& dfd, dmatrix& dfe, dmatrix& dff, const dmatrix d, const dmatrix e, const dmatrix f, unsigned long int pos_map, const int maxn, const int dt)
{
	PMap* map = (PMap*) pos_map;
	const int imin = map->imin;
	const int imax = map->imax;

	dmatrix ybet(imin,imax,0,maxn-1);
	dmatrix w(imin,imax,0,maxn-1);
	ybet.initialize();
	w.initialize();

	//recompute w and ybet
	for (int i=imin ; i<=imax ; i++){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		ybet[i][jmin] = 1/(e[i][jmin]+dt);
        	for (int j=jmin+1; j<=jmax ; j++){
			w(i,j) = e(i,j)+dt-f(i,j-1)*d(i,j)*ybet(i,j-1);
			ybet(i,j) = 1/w(i,j);
		}
	}
	for (int i=imax; i>=imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j=jmax; j>=jmin+1; j--){
			
			//ybet[i][j] = 1/w[i][j];
			dfw(i,j)    -= (1/(w(i,j)*w(i,j)))*dfybet(i,j);
			dfybet(i,j)  = 0.0;

			//w[i][j] = e[i][j]+dt-f[i][j-1]*d[i][j]*ybet[i][j-1];
			dfe(i,j)     += dfw(i,j);
			dff(i,j-1)   -= d(i,j)*ybet(i,j-1)*dfw(i,j);
			dfd(i,j)     -= f(i,j-1)*ybet(i,j-1)*dfw(i,j);
			dfybet(i,j-1)-= f(i,j-1)*d(i,j)*dfw(i,j);
			dfw(i,j)      = 0.0;

		}
		//ybet[i][jmin] = 1/(e[i][jmin]+dt);
		dfe(i,jmin)   -= (1/((e(i,jmin)+dt)*(e(i,jmin)+dt)))*dfybet(i,jmin);
		dfybet(i,jmin) = 0.0;
	}
}

void dfd1_comp(const char pos, double& dfd1, double& dfsigmam, double& dfsigma, double& dfuinf, const double uinf, const double twodd, const double d)
{
	switch (pos){
        case SANS:		
		case D_FERME:
			if (uinf > 0){	
	                    	//d1 = -((sigmam+sigma)/twodd)-uinf/d;
				dfsigmam -= (1/twodd) * dfd1;
				dfsigma  -= (1/twodd) * dfd1;
				dfuinf   -= (1/d) * dfd1;
				dfd1      = 0.0;
			} else {
				//d1 = -((sigmam+sigma)/twodd);
				dfsigmam -= (1/twodd) * dfd1;
				dfsigma  -= (1/twodd) * dfd1;
				dfd1      = 0.0;

			}
		break;
	}
}

void dfd2_comp(const char pos, double& dfd2, double& dfsigmam, double& dfsigma, double& dfsigmap, double& dfu, const double u, const double twodd, const double d)
{
	switch (pos){
	        case SANS:
			if (u > 0){
				//d2 = ((sigmam+2*sigma+sigmap)/twodd)+u/d;
				dfsigmam += (1/twodd) * dfd2;
				dfsigma  += (2/twodd) * dfd2;
				dfsigmap += (1/twodd) * dfd2;
				dfu      += (1/d) * dfd2;
				dfd2      = 0.0;			
			} else {
				//d2 = ((sigmam+2*sigma+sigmap)/twodd)-u/d;
				dfsigmam += (1/twodd) * dfd2;
				dfsigma  += (2/twodd) * dfd2;
				dfsigmap += (1/twodd) * dfd2;
				dfu      -= (1/d) * dfd2;
				dfd2      = 0.0;			
			}
		break;

		case G_FERME:
			if (u > 0){
				//d2 = ((sigma+sigmap)/twodd)+u/d;
				dfsigma  += (1/twodd) * dfd2;
				dfsigmap += (1/twodd) * dfd2;
				dfu      += (1/d) * dfd2;
				dfd2      = 0.0;
			} else {
				//d2 = ((sigma+sigmap)/twodd);
				dfsigma  += (1/twodd) * dfd2;
				dfsigmap += (1/twodd) * dfd2;
				dfd2      = 0.0;
			}				
		break;

		case D_FERME:
			if (u > 0){
				//d2 = ((sigmam+sigma)/twodd);
				dfsigmam += (1/twodd) * dfd2;
				dfsigma  += (1/twodd) * dfd2;
				dfd2      = 0.0;
			} else {			
				//d2 = ((sigmam+sigma)/twodd)-u/d;
				dfsigmam += (1/twodd) * dfd2;
				dfsigma  += (1/twodd) * dfd2;
				dfu      -= (1/d) * dfd2;
				dfd2      = 0.0;
			}
		break;	
	}
}

void dfd3_comp(const char pos, double& dfd3, double& dfsigma, double& dfsigmap, double& dfusup, const double usup, const double twodd, const double d)
{
	switch (pos){
        	case SANS:
			case G_FERME:
			if (usup > 0){	
	                    	//d3 = -((sigma+sigmap)/twodd);
				dfsigma  -= (1/twodd) * dfd3;
				dfsigmap -= (1/twodd) * dfd3;
				dfd3      = 0.0;
			} else {
	                    	//d3 = -((sigma+sigmap)/twodd)+usup/d;
				dfsigma  -= (1/twodd) * dfd3;
				dfsigmap -= (1/twodd) * dfd3;
				dfusup   += (1/d) * dfd3;
				dfd3      = 0.0;

			}
		break;
	}
}

void dfdh_comp(const char pos, double& dfHinf, double& dfH, double& dfHsup, double dfdh, const double d)
{
	if (pos == SANS) {
		//dH = (habitat_sup - habitat_inf)/(2*d);
		dfHsup += dfdh / (2*d);
		dfHinf -= dfdh / (2*d);
		dfdh    = 0.0;
	}
	else if (pos == D_FERME) {
		//dH = (habitat - habitat_inf)/d;
		dfH    += dfdh / d;
		dfHinf -= dfdh / d;
		dfdh    = 0.0;
	}
	else if (pos == G_FERME) {
		//dH = (habitat_sup - habitat)/d;
		dfHsup += dfdh / d;
		dfH    -= dfdh / d;
		dfdh    = 0.0;
	}
}

void dfdhx_comp(const char pos, double& dfHinf, double& dfH, double& dfHsup, double dfdh, const double habitat, const double habitat_sup, const double habitat_inf, const double d)
{

	if (pos == SANS) {
		if (habitat_sup-habitat>=habitat-habitat_inf){
			//dH = (habitat_sup - habitat)/(d);
			dfHsup += dfdh / d;
			dfH    -= dfdh / d;
			dfdh    = 0.0;
		}
		else {
			//dH = (habitat - habitat_inf)/(d);
			dfH    += dfdh / d;
			dfHinf -= dfdh / d;
			dfdh    = 0.0;
		}
	}
	else if (pos == D_FERME) {
		//dH = (habitat - habitat_inf)/d;
		dfH    += dfdh / d;
		dfHinf -= dfdh / d;
		dfdh    = 0.0;
	}
	else if (pos == G_FERME) {
		//dH = (habitat_sup - habitat)/d;
		dfHsup += dfdh / d;
		dfH    -= dfdh / d;
		dfdh    = 0.0;
	}

}

void dfdhy_comp(const char pos, double& dfHinf, double& dfH, double& dfHsup, double dfdh, const double habitat, const double habitat_sup, const double habitat_inf, const double d)
{
	if (pos == SANS) {
		if (habitat_sup-habitat<=habitat-habitat_inf){
			//dH = (habitat_sup - habitat)/(d);
			dfHsup += dfdh / d;
			dfH    -= dfdh / d;
			dfdh    = 0.0;
		}
		else {
			//dH = (habitat - habitat_inf)/(d);
			dfH    += dfdh / d;
			dfHinf -= dfdh / d;
			dfdh    = 0.0;
		}
	}
	else if (pos == D_FERME) {
		//dH = (habitat - habitat_inf)/d;
		dfH    += dfdh / d;
		dfHinf -= dfdh / d;
		dfdh    = 0.0;
	}
	else if (pos == G_FERME) {
		//dH = (habitat_sup - habitat)/d;
		dfHsup += dfdh / d;
		dfH    -= dfdh / d;
		dfdh    = 0.0;
	}

}



