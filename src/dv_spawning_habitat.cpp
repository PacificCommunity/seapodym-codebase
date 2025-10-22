#include "VarSimtunaFunc.h"

///Main function with memory control and adjoint functions for: 
///spawning habitat functions. This function depends on the temperature
///in the surface layer, food for larvae and predators of larvae. 
///There is optional oxygen function, which is currently hard-coded and
///invalidated through the parameter value.
///Forward functions are in spawning_habitat.cpp


int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void dv_Hs_comp(void);
void dv_Hs_fsst2_comp(void);
void dv_Hs_fsst3_comp(void);
void dfspawning_habitat_O2(double& dfa, double& dfb, double& dfc, double& dfd, double& dfe, double& dfH, const double a, const double b, const double c, const double d, const double e, const double ssv, const double SST, const double preys, const double predators, const double O2_l2);
void dfspawning_habitat_fsst2(double& dfa, double& dfb, double& dfc, double& dfd, double& dfe, double& dff, double& dfH, const double a, const double b, const double c, const double d, const double e, const double f, const double ssv, const double SST, const double preys, const double predators, const double O2_l2);
void dfspawning_habitat_fsst3(double& dfa, double& dfb, double& dfc, double& dfd, double& dfe, double& dff, double& dfg, double& dfH, const double a, const double b, const double c, const double d, const double e, const double f, const double g, const double ssv, const double SST, const double preys, const double predators, const double O2_l2);
double pred_surface_comp(dvector forage, const double DL, const int nb_forage, ivector day_layer, ivector night_layer);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);
double gaussian(const double x, const double mu, const double sigma);
double agaussian(const double x, const double mu, const double sigma_left, const double sigma_right);
double lognormal(const double x, const double mu, const double sigma);
double sigmoid1(const double tau, const double delta);
double logistic_thresholds(const double x, const double xmin, const double xmax, const double ek1, const double ek2);

void VarSimtunaFunc::Spawning_Habitat(VarParamCoupled& param, CMatrices& mat, const PMap& map, dvar_matrix& Hs, const double sigma_sp_var, int sp, const int t_count, const int jday)
{
	Hs.initialize();
	const int nb_forage = param.get_nbforage();

	const int fsst = param.fsst_type[sp];

	dvariable a,b,f,g;
	if (!param.uncouple_sst_larvae[sp]){
		a = param.dvarsA_sst_spawning[sp];
		b = param.dvarsB_sst_spawning[sp];
	}
	else {
		a = param.dvarsA_sst_larvae[sp];
		b = param.dvarsB_sst_larvae[sp];
	
		//Different fsst used in uncoupled mode only, when
		//a_sst and b_sst control fsst of Hs and not the thermal function of Ha
		if (fsst != 1){
			f = param.dvarsC_sst_larvae[sp];
			dvmatr6 = f;
		
			if (fsst == 3){
				g = param.dvarsD_sst_larvae[sp];
				dvmatr7 = g;
			}
		}
	}

	dvariable c = param.dvarsAlpha_hsp_prey[sp];
	dvariable d = param.dvarsAlpha_hsp_predator[sp];
	dvariable e = param.dvarsBeta_hsp_predator[sp];

	dvmatr1 = a;
	dvmatr2 = b;
	dvmatr3 = c;
	dvmatr4 = d;
	dvmatr5 = e;
	
	Hs_comp(param, mat, map, Hs, sigma_sp_var, sp, jday, t_count);
	
	save_identifier_string2((char*)"Hs_comp_begin");
	a.save_prevariable_value();
	dvmatr1.save_dvar_matrix_position();
	b.save_prevariable_value();
	dvmatr2.save_dvar_matrix_position();
	c.save_prevariable_value();
	dvmatr3.save_dvar_matrix_position();
	d.save_prevariable_value();
	dvmatr4.save_dvar_matrix_position();
	e.save_prevariable_value();
	dvmatr5.save_dvar_matrix_position();
	if (fsst != 1){
		f.save_prevariable_value();
		dvmatr6.save_dvar_matrix_position();
		if (fsst == 3){
			g.save_prevariable_value();
			dvmatr7.save_dvar_matrix_position();
		}
	}
	Hs.save_dvar_matrix_position();	
	unsigned long int pmap   = (unsigned long int)&map;
	save_long_int_value(pmap);
	unsigned long int cmat = (unsigned long int)&mat;
	save_long_int_value(cmat);
	param.night_layer.save_ivector_value();
	param.night_layer.save_ivector_position();
	param.day_layer.save_ivector_value();
	param.day_layer.save_ivector_position();
	save_double_value(sigma_sp_var);
	save_double_value(param.pp_transform);
	save_int_value(nb_forage);
	save_int_value(jday);
	save_int_value(t_count);
	save_identifier_string2((char*)"Hs_comp_end");

	if (fsst == 1)
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Hs_comp);
	if (fsst == 2)
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Hs_fsst2_comp);
	if (fsst == 3)
		gradient_structure::GRAD_STACK1->set_gradient_stack(dv_Hs_fsst3_comp);
}

void dv_Hs_comp(void)
{
	verify_identifier_string2((char*)"Hs_comp_end");
	unsigned t_count = restore_int_value();
	unsigned jday    = restore_int_value();
	unsigned nbf     = restore_int_value();
	double pp_transform = restore_double_value();
	double ssv = restore_double_value();
	const ivector_position dlay_pos = restore_ivector_position();
	ivector dlayer = restore_ivector_value(dlay_pos);
	const ivector_position nlay_pos = restore_ivector_position();
	ivector nlayer = restore_ivector_value(nlay_pos);
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position Hspos = restore_dvar_matrix_position();
	const dvar_matrix_position epos  = restore_dvar_matrix_position();
	double e = restore_prevariable_value();
	const dvar_matrix_position dpos  = restore_dvar_matrix_position();
	double d = restore_prevariable_value();
	const dvar_matrix_position cpos  = restore_dvar_matrix_position();
	double c = restore_prevariable_value();
	const dvar_matrix_position bpos  = restore_dvar_matrix_position();
	double b = restore_prevariable_value();
	const dvar_matrix_position apos  = restore_dvar_matrix_position();
	double a = restore_prevariable_value();

	verify_identifier_string2((char*)"Hs_comp_begin");
	dmatrix dfHs 	 = restore_dvar_matrix_derivatives(Hspos);
	dmatrix dfe 	 = restore_dvar_matrix_derivatives(epos);
	dmatrix dfd 	 = restore_dvar_matrix_derivatives(dpos);
	dmatrix dfc 	 = restore_dvar_matrix_derivatives(cpos);
	dmatrix dfb 	 = restore_dvar_matrix_derivatives(bpos);
	dmatrix dfa 	 = restore_dvar_matrix_derivatives(apos);

	CMatrices* mat	= (CMatrices*) pos_mat;
	PMap* map	= (PMap*) pos_map;

	const int imax = map->imax;
	const int imin = map->imin;
	dmatrix SST(imin, imax, map->jinf, map->jsup);
	SST = mat->sst(t_count);
	dmatrix O2_meso(imin, imax, map->jinf, map->jsup);
	O2_meso = mat->oxygen(t_count,1);
	dmatrix np1(imin, imax, map->jinf, map->jsup);
	np1 = mat->np1(t_count);
	d3_array forage;
	forage.allocate(0,nbf-1);
	for (unsigned int f=0; f<nbf; f++){
		forage(f).allocate(imin, imax, map->jinf, map->jsup);
		forage(f) = mat->forage(t_count,f);
	}
	dvector F(0,nbf-1);

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){		

				const double DL = mat->daylength(jday,j)/24.0;
				double preys = np1(i,j) * pp_transform;
				for (unsigned int n=0; n<nbf; n++)
					F(n) = forage(n,i,j);		
				double predators = pred_surface_comp(F,DL,nbf,dlayer,nlayer);

				dfspawning_habitat_O2(dfa(i,j),dfb(i,j),dfc(i,j),dfd(i,j),dfe(i,j),dfHs(i,j),a,b,c,d,e,ssv,SST(i,j),preys,predators,O2_meso(i,j));
			}
		}
	}
	dfa.save_dmatrix_derivatives(apos); 
	dfb.save_dmatrix_derivatives(bpos);
	dfc.save_dmatrix_derivatives(cpos);
	dfd.save_dmatrix_derivatives(dpos);
	dfe.save_dmatrix_derivatives(epos);
	dfHs.save_dmatrix_derivatives(Hspos);
}

void dfspawning_habitat_O2(double& dfa, double& dfb, double& dfc, double& dfd, double& dfe, double& dfH, const double a, const double b, const double c, const double d, const double e, const double ssv, const double SST, const double preys, const double predators, const double O2_l2)
{
	double f_sst  = gaussian(SST,b,a*ssv);
	double f_prey = preys*preys/(c+preys*preys);
	double f_prey_norm = f_prey * (1.0 + c);
	//double f_pred = gaussian(predators,d,e);
	//if lognormal, comment f_pred above and uncomment the following:
	double f_pred = exp(d-0.5*e*e-pow(log(predators)-d,2.0)/(2.0*e*e))/predators;//lognormal
	double f_oxy  = (1.0/(1.0+pow(0.01,O2_l2-0.1)));
	double Hs     = f_sst * f_prey_norm * f_pred * f_oxy;

	double logp = log(predators);

        //Hs = f_sst * f_prey * f_pred * f_oxy;
	dfa += (pow(SST-b,2.0) / pow(a*ssv,3.0)) * ssv * Hs * dfH;
	dfb += ((SST-b) / (a*a*ssv*ssv))* Hs * dfH;
	
	dfc += (f_prey-1.0/(c+preys*preys)) * (Hs / (1.0 + c)) * dfH;
        //dfd += ((predators-d)/(e*e)) * Hs * dfH;
        //dfe += (pow(predators-d,2.0) / pow(e,3.0)) * Hs * dfH;
	//if lognormal, comment dfd and dfe above and uncomment the following:
	dfd += (1.0+(logp-d)/(e*e)) * Hs * dfH;
	dfe += (pow(logp-d,2.0)/pow(e,3.0)-e) * Hs * dfH;
        dfH = 0.0;
}

void dv_Hs_fsst2_comp(void)
{
	verify_identifier_string2((char*)"Hs_comp_end");
	unsigned t_count = restore_int_value();
	unsigned jday    = restore_int_value();
	unsigned nbf     = restore_int_value();
	double pp_transform = restore_double_value();
	double ssv = restore_double_value();
	const ivector_position dlay_pos = restore_ivector_position();
	ivector dlayer = restore_ivector_value(dlay_pos);
	const ivector_position nlay_pos = restore_ivector_position();
	ivector nlayer = restore_ivector_value(nlay_pos);
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position Hspos = restore_dvar_matrix_position();
	const dvar_matrix_position fpos  = restore_dvar_matrix_position();
	double f = restore_prevariable_value();
	const dvar_matrix_position epos  = restore_dvar_matrix_position();
	double e = restore_prevariable_value();
	const dvar_matrix_position dpos  = restore_dvar_matrix_position();
	double d = restore_prevariable_value();
	const dvar_matrix_position cpos  = restore_dvar_matrix_position();
	double c = restore_prevariable_value();
	const dvar_matrix_position bpos  = restore_dvar_matrix_position();
	double b = restore_prevariable_value();
	const dvar_matrix_position apos  = restore_dvar_matrix_position();
	double a = restore_prevariable_value();

	verify_identifier_string2((char*)"Hs_comp_begin");
	dmatrix dfHs 	 = restore_dvar_matrix_derivatives(Hspos);
	dmatrix dff 	 = restore_dvar_matrix_derivatives(fpos);
	dmatrix dfe 	 = restore_dvar_matrix_derivatives(epos);
	dmatrix dfd 	 = restore_dvar_matrix_derivatives(dpos);
	dmatrix dfc 	 = restore_dvar_matrix_derivatives(cpos);
	dmatrix dfb 	 = restore_dvar_matrix_derivatives(bpos);
	dmatrix dfa 	 = restore_dvar_matrix_derivatives(apos);

	CMatrices* mat	= (CMatrices*) pos_mat;
	PMap* map	= (PMap*) pos_map;

	const int imax = map->imax;
	const int imin = map->imin;
	dmatrix SST(imin, imax, map->jinf, map->jsup);
	SST = mat->sst(t_count);
	dmatrix O2_meso(imin, imax, map->jinf, map->jsup);
	O2_meso = mat->oxygen(t_count,1);
	dmatrix np1(imin, imax, map->jinf, map->jsup);
	np1 = mat->np1(t_count);
	d3_array forage;
	forage.allocate(0,nbf-1);
	for (unsigned int f=0; f<nbf; f++){
		forage(f).allocate(imin, imax, map->jinf, map->jsup);
		forage(f) = mat->forage(t_count,f);
	}
	dvector F(0,nbf-1);

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){		

				const double DL = mat->daylength(jday,j)/24.0;
				double preys = np1(i,j) * pp_transform;
				for (unsigned int n=0; n<nbf; n++)
					F(n) = forage(n,i,j);		
				double predators = pred_surface_comp(F,DL,nbf,dlayer,nlayer);

				dfspawning_habitat_fsst2(dfa(i,j),dfb(i,j),dfc(i,j),dfd(i,j),dfe(i,j),dff(i,j),dfHs(i,j),a,b,c,d,e,f,ssv,SST(i,j),preys,predators,O2_meso(i,j));
			}
		}
	}
	dfa.save_dmatrix_derivatives(apos); 
	dfb.save_dmatrix_derivatives(bpos);
	dfc.save_dmatrix_derivatives(cpos);
	dfd.save_dmatrix_derivatives(dpos);
	dfe.save_dmatrix_derivatives(epos);
	dff.save_dmatrix_derivatives(fpos);
	dfHs.save_dmatrix_derivatives(Hspos);
}

void dfspawning_habitat_fsst2(double& dfa, double& dfb, double& dfc, double& dfd, double& dfe, double& dff, double& dfH, const double a, const double b, const double c, const double d, const double e, const double f, const double ssv, const double SST, const double preys, const double predators, const double O2_l2)
{	
	double f_sst  = agaussian(SST,b,a*ssv,f);

	double f_prey = preys*preys/(c+preys*preys);
	double f_prey_norm = f_prey * (1.0 + c);
	//double f_pred = gaussian(predators,d,e);
	//if lognormal, comment f_pred above and uncomment the following:
	double f_pred = exp(d-0.5*e*e-pow(log(predators)-d,2.0)/(2.0*e*e))/predators;//lognormal
	double f_oxy  = (1.0/(1.0+pow(0.01,O2_l2-0.1)));
	double Hs     = f_sst * f_prey_norm * f_pred * f_oxy;

	double logp = log(predators);

        //Hs = f_sst * f_prey * f_pred * f_oxy;
	if (SST<b){
		dfa += (pow(SST-b,2.0) / pow(a*ssv,3.0)) * ssv * Hs * dfH;
		dfb += ((SST-b) / (a*a*ssv*ssv))* Hs * dfH;
	} else {
		dff += (pow(SST-b,2.0) / pow(f,3.0)) * Hs * dfH;
		dfb += ((SST-b) / (f*f))* Hs * dfH;
	}
	
	dfc += (f_prey-1.0/(c+preys*preys)) * (Hs / (1.0 + c)) * dfH;
        //dfd += ((predators-d)/(e*e)) * Hs * dfH;
        //dfe += (pow(predators-d,2.0) / pow(e,3.0)) * Hs * dfH;
	//if lognormal, comment dfd and dfe above and uncomment the following:
	dfd += (1.0+(logp-d)/(e*e)) * Hs * dfH;
	dfe += (pow(logp-d,2.0)/pow(e,3.0)-e) * Hs * dfH;
        dfH = 0.0;
}

void dv_Hs_fsst3_comp(void)
{
	verify_identifier_string2((char*)"Hs_comp_end");
	unsigned t_count = restore_int_value();
	unsigned jday    = restore_int_value();
	unsigned nbf     = restore_int_value();
	double pp_transform = restore_double_value();
	double ssv = restore_double_value();
	const ivector_position dlay_pos = restore_ivector_position();
	ivector dlayer = restore_ivector_value(dlay_pos);
	const ivector_position nlay_pos = restore_ivector_position();
	ivector nlayer = restore_ivector_value(nlay_pos);
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position Hspos = restore_dvar_matrix_position();
	const dvar_matrix_position gpos  = restore_dvar_matrix_position();
	double g = restore_prevariable_value();
	const dvar_matrix_position fpos  = restore_dvar_matrix_position();
	double f = restore_prevariable_value();
	const dvar_matrix_position epos  = restore_dvar_matrix_position();
	double e = restore_prevariable_value();
	const dvar_matrix_position dpos  = restore_dvar_matrix_position();
	double d = restore_prevariable_value();
	const dvar_matrix_position cpos  = restore_dvar_matrix_position();
	double c = restore_prevariable_value();
	const dvar_matrix_position bpos  = restore_dvar_matrix_position();
	double b = restore_prevariable_value();
	const dvar_matrix_position apos  = restore_dvar_matrix_position();
	double a = restore_prevariable_value();

	verify_identifier_string2((char*)"Hs_comp_begin");
	dmatrix dfHs 	 = restore_dvar_matrix_derivatives(Hspos);
	dmatrix dfg 	 = restore_dvar_matrix_derivatives(gpos);
	dmatrix dff 	 = restore_dvar_matrix_derivatives(fpos);
	dmatrix dfe 	 = restore_dvar_matrix_derivatives(epos);
	dmatrix dfd 	 = restore_dvar_matrix_derivatives(dpos);
	dmatrix dfc 	 = restore_dvar_matrix_derivatives(cpos);
	dmatrix dfb 	 = restore_dvar_matrix_derivatives(bpos);
	dmatrix dfa 	 = restore_dvar_matrix_derivatives(apos);

	CMatrices* mat	= (CMatrices*) pos_mat;
	PMap* map	= (PMap*) pos_map;

	const int imax = map->imax;
	const int imin = map->imin;
	dmatrix SST(imin, imax, map->jinf, map->jsup);
	SST = mat->sst(t_count);
	dmatrix O2_meso(imin, imax, map->jinf, map->jsup);
	O2_meso = mat->oxygen(t_count,1);
	dmatrix np1(imin, imax, map->jinf, map->jsup);
	np1 = mat->np1(t_count);
	d3_array forage;
	forage.allocate(0,nbf-1);
	for (unsigned int f=0; f<nbf; f++){
		forage(f).allocate(imin, imax, map->jinf, map->jsup);
		forage(f) = mat->forage(t_count,f);
	}
	dvector F(0,nbf-1);

	for (int i = imax; i >= imin; i--){
		const int jmin = map->jinf[i];
		const int jmax = map->jsup[i];
		for (int j = jmax; j >= jmin; j--){
			if (map->carte(i,j)){		

				const double DL = mat->daylength(jday,j)/24.0;
				double preys = np1(i,j) * pp_transform;
				for (unsigned int n=0; n<nbf; n++)
					F(n) = forage(n,i,j);		
				double predators = pred_surface_comp(F,DL,nbf,dlayer,nlayer);

				dfspawning_habitat_fsst3(dfa(i,j),dfb(i,j),dfc(i,j),dfd(i,j),dfe(i,j),dff(i,j),dfg(i,j),dfHs(i,j),a,b,c,d,e,f,g,ssv,SST(i,j),preys,predators,O2_meso(i,j));
			}
		}
	}
	dfa.save_dmatrix_derivatives(apos); 
	dfb.save_dmatrix_derivatives(bpos);
	dfc.save_dmatrix_derivatives(cpos);
	dfd.save_dmatrix_derivatives(dpos);
	dfe.save_dmatrix_derivatives(epos);
	dff.save_dmatrix_derivatives(fpos);
	dfg.save_dmatrix_derivatives(gpos);
	dfHs.save_dmatrix_derivatives(Hspos);
}

void dfspawning_habitat_fsst3(double& dfa, double& dfb, double& dfc, double& dfd, double& dfe, double& dff, double& dfg, double& dfH, const double a, const double b, const double c, const double d, const double e, const double f, const double g, const double ssv, const double SST, const double preys, const double predators, const double O2_l2)
{	
	double f_sst  = logistic_thresholds(SST,a,b,f,g);

	double f_prey = preys*preys/(c+preys*preys);
	double f_prey_norm = f_prey * (1.0 + c);
	//double f_pred = normal(predators,d,e);
	//if lognormal, comment f_pred above and uncomment the following:
	double f_pred = exp(d-0.5*e*e-pow(log(predators)-d,2.0)/(2.0*e*e))/predators;//lognormal
	double f_oxy  = (1.0/(1.0+pow(0.01,O2_l2-0.1)));
	double Hs     = f_sst * f_prey_norm * f_pred * f_oxy;

	double logp = log(predators);

	double expr11 = pow(a,SST-b);
	double expr12 = (1.0+expr11)*(1.0+expr11);
	double expr21 = pow(f,g-SST);
	double expr22 = (1.0+expr21)*(1.0+expr21);
	double expr0  = f_prey_norm * f_pred * f_oxy * dfH;

        //Hs = f_sst * f_prey * f_pred * f_oxy;
	dfa -= ((SST-b)*expr11 / (a*expr12*(1.0+expr21))) * expr0; 
	dfb += (expr11 * log(a) / (expr12*(1.0+expr21))) * expr0;
	dff -= ((g-SST)*expr21 / (f*expr22*(1.0+expr11))) * expr0; 
	dfg -= (expr21 * log(f) / (expr22*(1.0+expr11))) * expr0;
	
	dfc += (f_prey-1.0/(c+preys*preys)) * (Hs / (1.0 + c)) * dfH;
        //dfd += ((predators-d)/(e*e)) * Hs * dfH;
        //dfe += (pow(predators-d,2.0) / pow(e,3.0)) * Hs * dfH;
	//if lognormal, comment dfd and dfe above and uncomment the following:
	dfd += (1.0+(logp-d)/(e*e)) * Hs * dfH;
	dfe += (pow(logp-d,2.0)/pow(e,3.0)-e) * Hs * dfH;
        dfH = 0.0;
}

