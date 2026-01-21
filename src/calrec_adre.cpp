#include "calpop.h"

///Forward function for:
///calrec for larval and juvenile life stages, i.e. with passive drift only. 
///These routines solve ADR equations for larvae and juveniles, advected 
///passively and diffused with water using the same ADI method as for adults
///with inner timestep deltaT/N.  


int save_identifier_string2(char* str);

void CCalpop::calrec1(const PMap& map, dvar_matrix& uu, const dmatrix& mortality)
{
	if (map.global == 1) {
		calrec_GO(map,uu);
		return;
	}
	DVECTOR uvec(0, maxn - 1);
	DVECTOR rhs(0, maxn - 1);
	DVECTOR gam(0, maxn - 1);

	uu.save_dvar_matrix_value();
	save_identifier_string2((char*)"One_step_calrec_uu");
	
	for (int itr = 1; itr <= iterationNumber; itr++) 
	{
		for (int j = map.jmin; j <= map.jmax; j++)
		{
			const int imin = map.iinf[j]; 
			const int imax = map.isup[j];

			for (int i = imin; i <= imax; i++)
			{   
				if (map.jinf[i] <= j && j <= map.jsup[i]){

					rhs[i]=-d[i][j]*value(uu(i,j-1)) + (2*iterationNumber-e[i][j])*value(uu(i,j)) - f[i][j]*value(uu(i,j+1));

                } else {
					cerr << "Error: cf file " << __FILE__ << ", line " << __LINE__ << ": itr = " << itr << ", i=" << i << ",j=" << j <<  endl;
					cout << __LINE__ << endl; exit(1);
					//rhs[i]=(2*iterationNumber)*value(uu(i,j));
				}
			}
			
			tridag(a[j],xbet[j],c[j],rhs,uvec,gam,imin,imax);

			for (int i = imin; i <= imax; i++)
				uuint(i,j) = uvec[i];
		} 
		
		for (int i = map.imin; i <= map.imax; i++)
		{
			DVECTOR& dvUUINTprev = uuint[i-1];
			DVECTOR& dvUUINT     = uuint[i];
			DVECTOR& dvUUINTnext = uuint[i+1];

			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++)
			{
				if (map.iinf[j] <= i && i <= map.isup[j])
				{
					rhs[j]=(-a[j][i]*dvUUINTprev[j])+(2*iterationNumber-bm[j][i])*dvUUINT[j] -(c[j][i]*dvUUINTnext[j]);
				}
				else
				{
					cout << __LINE__ << endl; exit(1);
					//rhs[j]=((2*iterationNumber - mortality[i][j]))*dvUUINT[j];
				}
			}

			tridag(d[i],ybet[i],f[i],rhs,uvec,gam,jmin,jmax);
			for (int j = jmin; j <= jmax; j++)
				uu.elem_value(i,j) = uvec[j];
		} 
	}
}
