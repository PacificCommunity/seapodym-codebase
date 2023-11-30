#include "calpop.h"

///Main function with memory control and adjoint functions for: 
///calrec for larval and juvenile life stages, i.e. with passive drift only. 
///These routines solve ADR equations for larvae and juveniles, advected 
///passively and diffused with water using the same ADI method as for adults
///with inner timestep deltaT/N.  
///Forward functions are in calrec_adre.cpp


void dv_calrec_adre();
void dftridag_y(dvector& a, dvector& bet, dvector& c, dvector& dfrhs, dvector& dfuvec, dvector& dfgam, int inf, int sup);
void dftridag_x(dvector& a, dvector& bet, dvector& c,dvector& rhs,dvector& uvec, dvector& gam, dvector& dfbet, dvector& dfrhs, dvector& dfuvec, dvector& dfgam, int inf, int sup);
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);


void CCalpop::Calrec_juv(const PMap& map, CMatrices& mat, dvar_matrix& uu, dvar_matrix& mortality, const int t_count)
{
	bm = value(dvarsBM);
	xbet = value(Xbet);
	dmatrix mort_c = value(mortality);
	calrec1(map,uu,mort_c);

	save_identifier_string2((char*)"calrec_adre_begin");
	uu.save_dvar_matrix_position();
	dvarsBM.save_dvar_matrix_position();
	Xbet.save_dvar_matrix_position();
	save_int_value(t_count);
	mortality.save_dvar_matrix_value();
	mortality.save_dvar_matrix_position();
	unsigned long int pop   = (unsigned long int)this;
	save_long_int_value(pop);
	unsigned long int pmap  = (unsigned long int)&map;
	save_long_int_value(pmap);
	unsigned long int cmat  = (unsigned long int)&mat;
	save_long_int_value(cmat);
	save_identifier_string2((char*)"calrec_adre_end");
	

	gradient_structure::GRAD_STACK1->set_gradient_stack(dv_calrec_adre);
} 

void dv_calrec_adre()
{

	verify_identifier_string2((char*)"calrec_adre_end");
	unsigned long int pos_mat = restore_long_int_value();
	unsigned long int pos_map = restore_long_int_value();
	unsigned long int pos_pop = restore_long_int_value();
	const dvar_matrix_position m_pos = restore_dvar_matrix_position();
	const dmatrix mort = restore_dvar_matrix_value(m_pos);
	const int t_count = restore_int_value();
	const dvar_matrix_position xbet_pos = restore_dvar_matrix_position();
	const dvar_matrix_position bm_pos   = restore_dvar_matrix_position();
	const dvar_matrix_position uu_pos   = restore_dvar_matrix_position();
	verify_identifier_string2((char*)"calrec_adre_begin");

	dmatrix uu(uu_pos);
	dmatrix dfbm 	= restore_dvar_matrix_derivatives(bm_pos);
	dmatrix dfxbet 	= restore_dvar_matrix_derivatives(xbet_pos);
	dmatrix dfuu	= restore_dvar_matrix_derivatives(uu_pos);	

	CMatrices* mat = (CMatrices*) pos_mat;
	PMap* map = (PMap*) pos_map;
	CCalpop* pop = (CCalpop*) pos_pop;

	const int iterationNumber  = pop->get_iterationN();  
	const int maxn = pop->get_maxn();  

	const int jinf = map->jmin;
	const int jsup = map->jmax;
	const int iinf = map->imin;
	const int isup = map->imax;

	dvector uvec(0,maxn-1);
	dvector rhs(0,maxn-1);
	dvector gam(0,maxn-1);
	uvec.initialize();
	rhs.initialize();
	gam.initialize();
	
	dvector dfuvec(0,maxn-1); 
	dvector dfxgam(0,maxn-1);
	dvector dfxrhs(0,maxn-1);
	dfxgam.initialize();
	dfxrhs.initialize();
	dfuvec.initialize();

	dvector dfygam(jinf,jsup);
	dvector dfyrhs(jinf,jsup);
	dfygam.initialize();
	dfyrhs.initialize();

	//recompute tridiag coefficients
	dmatrix a,bm,c,d,e,f;
	a.allocate(jinf,jsup,map->iinf,map->isup);
	bm.allocate(jinf,jsup,map->iinf,map->isup);
	c.allocate(jinf,jsup,map->iinf,map->isup);
	d.allocate(iinf,isup,map->jinf,map->jsup);
	e.allocate(iinf,isup,map->jinf,map->jsup);
	f.allocate(iinf,isup,map->jinf,map->jsup);

	a.initialize();
	bm.initialize();
	c.initialize();
	d.initialize();
	e.initialize();
	f.initialize();
	
	pop->RecompDiagCoef_juv(*map, *mat, t_count, mort, a, bm, c, d, e, f);

	//recompute xbet and ybet
	dmatrix xbet, ybet;
	xbet.allocate(jinf,jsup,map->iinf,map->isup);
	ybet.allocate(iinf,isup,map->jinf,map->jsup);
	xbet.initialize();
	ybet.initialize();
	pop->xbet_comp(*map, xbet, a, bm, c, 2*iterationNumber);
	pop->ybet_comp(*map, ybet, d, e, f, 2*iterationNumber);

	const dmatrix_position uuint_pos(pop->uuint);

	//recompute all 2*N-1 solutions 
	//(no need to get the last one 
	//as it is saved on the next time step)
	d3_array luu, luuint;//local uu and uuint computed from uu(t-1)
	luu.allocate(0, iterationNumber-1);
	luuint.allocate(1, iterationNumber);
	for (int itr = 1; itr <= iterationNumber; itr++){ 
		luu(itr-1).allocate(map->imin1, map->imax1, map->jinf1, map->jsup1);
		luuint(itr).allocate(map->imin1, map->imax1, map->jinf1, map->jsup1);
		luu(itr-1).initialize();
		luuint(itr).initialize();
	}

	verify_identifier_string2((char*)"One_step_calrec_uu");
	luu[0] = restore_dvar_matrix_value(uu_pos);

	//recompute all intermediate solutions:
	pop->RecompADI_step_fwd(*map, luu, luuint, a, bm, c, d, e, f, xbet, ybet);

	//ADJOINT FOR CALREC FUNCTION
	for (int itr = iterationNumber; itr >= 1; itr--){

		dmatrix dfuuint(uuint_pos);
		dfuuint.initialize();

		for (int i = isup; i >= iinf; i--){

			const int jmin = map->jinf[i];
			const int jmax = map->jsup[i];
		
			for (int j=jmin; j<=jmax; j++){
				//uu(i,j) = uvec(j); 
				dfuvec(j)+= dfuu(i,j);
				dfuu(i,j) = 0.0;
			}

			//dftridag_y(d(i),ybet(i),f(i),rhs,uvec,gam,dfyrhs,dfuvec,dfygam,jmin,jmax);
			dftridag_y(d(i),ybet(i),f(i),dfyrhs,dfuvec,dfygam,jmin,jmax);

			for (int j=jmax; j>=jmin; j--) {
				//rhs[j] = -a[j][i]*uuint[i-1][j]+(2*iterationNumber-bm[j][i])*uuint[i][j] - c[j][i]*uuint[i+1][j];
				dfuuint(i-1,j)-= a(j,i)*dfyrhs(j);
				dfbm(j,i)     -= luuint(itr,i,j)*dfyrhs(j);
				dfuuint(i,j)  += (2*iterationNumber-bm[j][i])*dfyrhs(j);
				dfuuint(i+1,j)-= c(j,i)*dfyrhs(j);
				dfyrhs(j)      = 0.0;
			}
		}

		for (int j = jsup; j >= jinf; j--){
			const int imin = map->iinf[j]; 
			const int imax = map->isup[j];
	
			// recomputing rhs(i)
			for (int i = imin; i <= imax; i++) {   

					rhs[i] = -d[i][j]*luu[itr-1][i][j-1] + (2*iterationNumber-e[i][j])*luu[itr-1][i][j] - f[i][j]*luu[itr-1][i][j+1];
				
			}

			for (int i = imax; i >= imin; i--){
				//uuint[i][j] = uvec[i];
				dfuvec(i)    += dfuuint(i,j);
				dfuuint(i,j)  = 0.0;
			}
			dftridag_x(a(j),xbet(j),c(j),rhs,uvec,gam,dfxbet(j),dfxrhs,dfuvec,dfxgam,imin,imax);

			for (int i=imax; i>=imin; i--){

					//rhs[i] = -d[i][j]*uu[i][j-1] + (2*iterationNumber-e[i][j])*uu[i][j] - f[i][j]*uu[i][j+1];
					dfuu(i,j-1) -= d(i,j)*dfxrhs(i);
					dfuu(i,j)   += (2*iterationNumber-e[i][j])*dfxrhs(i);
					dfuu(i,j+1) -= f(i,j)*dfxrhs(i);
					dfxrhs(i)    = 0.0;
			}
		}
	}
	dfbm.save_dmatrix_derivatives(bm_pos); 
	dfxbet.save_dmatrix_derivatives(xbet_pos); 
	dfuu.save_dmatrix_derivatives(uu_pos); 
		
}

void dftridag_x(dvector& a, dvector& bet, dvector& c,dvector& rhs,dvector& uvec, dvector& gam, dvector& dfbet, dvector& dfrhs, dvector& dfuvec, dvector& dfgam, int inf, int sup)
{
	int j;

	// recompute gam and uvec values
	gam[inf] = rhs[inf];
	for (int j=inf+1 ; j<=sup ; j++)
		gam[j] = rhs[j]-gam[j-1]*a[j]*bet[j-1];
	
	uvec[sup] = gam[sup]*bet[sup];
	for (int j=sup-1; j>=inf ; j--)
		uvec[j] = (gam[j]-c[j]*uvec[j+1])*bet[j];
	////////////////////////////////

 
 	for (j=inf; j<sup; j++){
		//uvec[j] = (gam[j]-c[j]*uvec[j+1])*bet[j];
		dfgam[j]   += bet[j]*dfuvec[j];
		dfuvec[j+1]-= c[j]*bet[j]*dfuvec[j];
		dfbet[j]   += (gam[j]-c[j]*uvec[j+1])*dfuvec[j]; 
		dfuvec[j]   = 0.0;
	}

	//uvec[sup] = gam[sup]*bet[sup];
  	dfgam[sup] += bet[sup]*dfuvec[sup];
	dfbet[sup] += gam[sup]*dfuvec[sup];
  	dfuvec[sup] = 0.0;

	for (j=sup;j>=inf+1;j--){
		//gam[j]   = rhs[j]-gam[j-1]*a[j]*bet[j-1];
		dfrhs[j]  += dfgam[j];
		dfgam[j-1]-= a[j]*bet[j-1]*dfgam[j];
		dfbet[j-1]-= gam[j-1]*a[j]*dfgam[j];
		dfgam[j]   = 0.0; 
	}

	//gam[inf] = rhs[inf];
	dfrhs[inf] += dfgam[inf];
	dfgam[inf]  = 0.0;


} // end of dftridag

void dftridag_y(dvector& a, dvector& bet, dvector& c, dvector& dfrhs, dvector& dfuvec, dvector& dfgam, int inf, int sup)
{
	int j;
 
 	for (j=inf; j<sup; j++){
		//uvec[j] = (gam[j]-c[j]*uvec[j+1])*bet[j];
		dfgam[j]   += bet[j]*dfuvec[j];
		dfuvec[j+1]-= c[j]*bet[j]*dfuvec[j];
		dfuvec[j]   = 0.0;
	}

	//uvec[sup] = gam[sup]*bet[sup];
  	dfgam[sup] += bet[sup]*dfuvec[sup];
  	dfuvec[sup] = 0.0;

	for (j=sup;j>=inf+1;j--){
		//gam[j]   = rhs[j]-gam[j-1]*a[j]*bet[j-1];
		dfrhs[j]  += dfgam[j];
		dfgam[j-1]-= a[j]*bet[j-1]*dfgam[j];
		dfgam[j]   = 0.0; 
	}

	//gam[inf] = rhs[inf];
	dfrhs[inf] += dfgam[inf];
	dfgam[inf]  = 0.0;



} // end of dftridag

