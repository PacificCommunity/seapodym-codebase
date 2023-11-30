#include "calpop.h"

///Forward function for: 
///precalrec and calrec for adults functions. These routines finalize 
///the computation of diagonal elements by adding total mortalily rates 
///to coefficient 'b' (precalrec part), and then solve the discretized 
///ADE equation iteratively using ADI method with inner timestep deltaT/N.  

int save_identifier_string2(char* str);

void CCalpop::calrec_with_catch(const PMap& map, CParam& param, dvar_matrix& uu, const dmatrix& C_obs, dvar_matrix& C_est)
{
	if (map.global == 1) {
		calrec_GO_with_catch(map,param,uu,C_obs,C_est);
		return;
	}
	DVECTOR uvec(0, maxn - 1);
	DVECTOR rhs(0, maxn - 1);
	DVECTOR gam(0, maxn - 1);

	uu.save_dvar_matrix_value();
	save_identifier_string2((char*)"One_step_calrec_uu");

	for (int itr = 1; itr <= iterationNumber; itr++) 
	{
		//cout << itr << " " << norm(value(uu)) << endl;
		for (int j = map.jmin; j <= map.jmax; j++)
		{
			const int imin = map.iinf[j]; 
			const int imax = map.isup[j];

			for (int i = imin; i <= imax; i++)
			{   
				if (map.jinf[i] <= j && j <= map.jsup[i]){

					rhs[i]=-d[i][j]*value(uu(i,j-1)) + (2*iterationNumber-e[i][j])*value(uu(i,j)) - f[i][j]*value(uu(i,j+1));
                         	} else {
					cout << __LINE__ << endl; exit(1);
					//rhs[i]=(2*iterationNumber)*value(uu(i,j));
				}
			}
			
			tridag(a[j],xbet[j],c[j],rhs,uvec,gam,imin,imax);

			for (int i = imin; i <= imax; i++){
				if (C_obs(i,j)==0)
					uuint(i,j) = uvec[i];
				else {
					double C_est_ij_itr = uvec[i] * param.func_limit_one(C_obs(i,j)/(uvec[i]+1e-14)) / iterationNumber;
					uuint(i,j) = uvec[i] - C_est_ij_itr;
					C_est.elem_value(i,j) += C_est_ij_itr;
				}
				if (uvec(i)>0 && uuint(i,j)<0) {
					cout << "Negative values due to catch removal in" << 
						"i=" << i << ", j="<< j << 
						", B=" << uvec(i) << ", C="<< C_obs(i,j) << endl;
				}
			}
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
				}
			}

			tridag(d[i],ybet[i],f[i],rhs,uvec,gam,jmin,jmax);
			for (int j = jmin; j <= jmax; j++)
				uu.elem_value(i,j) = uvec[j];
		} 
	}
}
