#include "calpop.h"

///Main function with memory control and adjoint functions for: 
///precalrec for larval and juvenile life stages. 
///This routine precomputes diagonal coefficient for calrec_adre 
///Forward function is in precalrec_juv.cpp


void dv_precalrec_juv_comp(void);
void dv_precalrec_adult_comp(void);
void dv_precalrec_adult_comp1(void);
void dfxbet_juv_comp(dmatrix& dfbm, dmatrix& dfxbet, dmatrix& dfw, const dmatrix a, const dmatrix bm, const dmatrix c, unsigned long int pos_map, const int maxn, const int dt);
void dfxbet_adult_comp(dmatrix& dfa, dmatrix& dfbm, dmatrix& dfc, dmatrix& dfxbet, dmatrix& dfw, const dmatrix a, const dmatrix bm, const dmatrix c, unsigned long int pos_map, const int maxn, const int dt);
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);


void CCalpop::Precalrec_juv(const PMap& map, CMatrices& mat, dvar_matrix& mortality, const int t_count)
{
	dmatrix M_c    = value(mortality);
	dmatrix bm_c   = value(dvarsBM);
	dmatrix xbet_c = value(Xbet);

	dvar_matrix W(0,maxn-1,0,maxn-1);
	W.initialize();

	bm_c = b;

	precalrec_juv_comp(map, bm_c, M_c);

	dvarsBM = nograd_assign(bm_c);

	bm = bm_c;
	
	xbet_comp(map, xbet_c, a, bm_c, c, 2*iterationNumber);
	
	Xbet = nograd_assign(xbet_c);
	
	save_identifier_string2((char*)"Precalrec_juv_begin");
	mortality.save_dvar_matrix_value();
	mortality.save_dvar_matrix_position();
	dvarsBM.save_dvar_matrix_position();
	Xbet.save_dvar_matrix_position();
	W.save_dvar_matrix_position();
	unsigned long int pmap = (unsigned long int)&map;
	save_long_int_value(pmap);
	unsigned long int pop = (unsigned long int)this;
	save_long_int_value(pop);	
	unsigned long int cmat = (unsigned long int)&mat;
	save_long_int_value(cmat);
	save_int_value(t_count);
	save_identifier_string2((char*)"Precalrec_juv_end");

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_precalrec_juv_comp);
}

void dv_precalrec_juv_comp(void)
{
	verify_identifier_string2((char*)"Precalrec_juv_end");
	unsigned t_count   = restore_int_value();
	unsigned long int pos_mat   = restore_long_int_value();
	unsigned long int pos_pop   = restore_long_int_value();
	unsigned long int pos_map   = restore_long_int_value();
	const dvar_matrix_position w_pos    = restore_dvar_matrix_position();
	const dvar_matrix_position xbet_pos = restore_dvar_matrix_position();
	const dvar_matrix_position bm_pos   = restore_dvar_matrix_position();
	const dvar_matrix_position M_pos    = restore_dvar_matrix_position();
	dmatrix mort = restore_dvar_matrix_value(M_pos);
	verify_identifier_string2((char*)"Precalrec_juv_begin");

	dmatrix dfM 	 = restore_dvar_matrix_derivatives(M_pos);
	dmatrix dfbm 	 = restore_dvar_matrix_derivatives(bm_pos);
	dmatrix dfxbet 	 = restore_dvar_matrix_derivatives(xbet_pos);
	dmatrix dfw 	 = restore_dvar_matrix_derivatives(w_pos);

	CMatrices* mat 	 = (CMatrices*) pos_mat;
	PMap* map 	 = (PMap*) pos_map;
	CCalpop* pop 	 = (CCalpop*) pos_pop;

	const int jinf = map->jmin;
	const int jsup = map->jmax;
	const int dt   = 2*pop->get_iterationN();
	const int maxn = pop->get_maxn();

	dmatrix a, bm, c;
	a.allocate(jinf,jsup,map->iinf,map->isup);
	bm.allocate(jinf,jsup,map->iinf,map->isup);
	c.allocate(jinf,jsup,map->iinf,map->isup);

	a.initialize();
	bm.initialize();
	c.initialize();

	pop->Recomp_abc_coef(*map, *mat, t_count, mort, a, bm, c);

	//xbet_comp(map,dt);
	dfxbet_juv_comp(dfbm,dfxbet,dfw,a,bm,c,pos_map,maxn,dt);

	for (int j = jsup; j >= jinf; j--){
		const int imin = map->iinf(j);
		const int imax = map->isup(j);
		for (int i = imax; i >= imin; i--){
			//bm[j][i] = b(j,i)+mort[i][j];
			dfM(i,j) += dfbm(j,i);
			dfbm(j,i) = 0.0;

		}
	}
	dfbm.save_dmatrix_derivatives(bm_pos);
	dfxbet.save_dmatrix_derivatives(xbet_pos); 
	dfw.save_dmatrix_derivatives(w_pos);
	dfM.save_dmatrix_derivatives(M_pos);

}

void dfxbet_juv_comp(dmatrix& dfbm, dmatrix& dfxbet, dmatrix& dfw, const dmatrix a, const dmatrix bm, const dmatrix c, unsigned long int pos_map, const int maxn, const int dt)
{
	PMap* map = (PMap*) pos_map;

	const int jmin = map->jmin;
	const int jmax = map->jmax;

	dmatrix xbet(jmin,jmax,0,maxn-1);
	dmatrix w(jmin,jmax,0,maxn-1);
	xbet.initialize();
	w.initialize();

	//recompute w and xbet
	for (int j=jmin; j<=jmax ; j++){
		const int imin = map->iinf[j];
		const int imax = map->isup[j];
		xbet[j][imin] = 1/(bm[j][imin]+dt);
	        for (int i=imin+1; i<=imax ; i++){
			w[j][i] = bm[j][i]+dt-c[j][i-1]*a[j][i]*xbet[j][i-1];
			xbet[j][i] = 1/w[j][i];
		}
	}

	for (int j=jmax; j>=jmin; j--){
		const int imin = map->iinf[j];
		const int imax = map->isup[j];
		for (int i=imax; i>imin; i--){
			
			//xbet[j][i] = 1/w[j][i];
			dfw(j,i)    -= (1/(w(j,i)*w(j,i)))*dfxbet(j,i);
			dfxbet(j,i)  = 0.0;

			//w[j][i] = bm[j][i]+dt-c[j][i-1]*a[j][i]*xbet[j][i-1];
			dfbm(j,i)    += dfw(j,i);
			dfxbet(j,i-1)-= c(j,i-1)*a(j,i)*dfw(j,i);
			dfw(j,i)      = 0.0;

		}
		//xbet[j][inf] = 1/(bm[j][imin]+dt);
		dfbm(j,imin)  -= (1/((bm(j,imin)+dt)*(bm(j,imin)+dt)))*dfxbet(j,imin);
		dfxbet(j,imin) = 0.0;
	}
}


