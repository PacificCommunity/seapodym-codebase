#include "calpop.h"

/////////////////////////////////////////////////////////////////
//     -----------------------------------------------------------
//     sous programme 'calrec' : resolution du systeme d'equations
//     sur chaque demi interval de temps.
//     -----------------------------------------------------------
//     itr=boucle sur les sub-iterations
//     la partie droite de l'equation est d'abord stockee dans le vecteur RHS[i] ou [j].
//     Sur la premiere moitie de l'intervalle de temps, j systemes de i equations sont resolus:
//     j fois D[i]*U[i]=RHS[i]  
//     Sur la seconde moitie de l'intervalle de temps, i systemes de j equations sont resolus:
//     i fois D[j]*U[j]=RHS[j]  
//     les i solutions UVEC[j] sont stockees dans UU[i][j] i*(UVEC[j])=>UU[i][j]
//      ---------------------------------------------------------
/////////////////////////////////////////////////////////////////

void CCalpop::calrec(const PMap& map, dvar_matrix& uu, dvar_matrix& mortality)
{	// DEBUT DU SOUS PROGRAMME CALREC
	
	dvar_vector uvec(0, maxn - 1);
	dvar_vector rhs(0, maxn - 1);
	dvar_vector gam(0, maxn - 1);

	for (int itr = 1; itr <= iterationNumber; itr++) 
	{
		//////////////////////////////////////////////////////////////////////////////////
		// 1- CALCUL DE LA PARTIE DROITE (RHS) DE L'EQUATION POUR L'INTERVALLE DE TEMPS 0 A 1/2
		//    RHS[] correspond au vecteur g[] de Sibert 
		/////////////////////////////////////////////////////////////////////////////////
		for (int j = map.jmin; j <= map.jmax; j++)
		{
			const int imin = map.iinf[j]; 
			const int imax = map.isup[j];
			for (int i = imin; i <= imax; i++)
			{   
				dvar_vector& dvUU = uu[i];
				if (map.jinf[i] <= j && j <= map.jsup[i]) {
					rhs[i] = -dvarsD[i][j]*dvUU[j-1] + (2*iterationNumber-dvarsE[i][j])*dvUU[j] - dvarsF[i][j]*dvUU[j+1];
				}
				else {
					cout << __LINE__ << endl; exit(1);
					rhs[i] = 2*iterationNumber*dvUU[j];
				}
			}

			// appel du sous programme tridag pour la premiere moitie de l'intervalle de temps
			tridag2(dvarsA[j],Xbet[j],dvarsC[j],rhs,uvec,gam,imin,imax);

			for (int i = imin; i <= imax; i++)
				Uuint[i][j]=uvec[i];
		}

		/////////////////////////////////////////////////////////////////////////
		// calcul de la partie droite de l'equation rrhs(i,j) pour 
		// l'intervalle de temps 1/2 a 1
		for (int i = map.imin; i <= map.imax; i++)
		{
			dvar_vector& dvUUINTprev = Uuint[i-1];
			dvar_vector& dvUUINT     = Uuint[i];
			dvar_vector& dvUUINTnext = Uuint[i+1];

			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++) {
				if (map.iinf[j] <= i && i <= map.isup[j]) {
					rhs[j]=-dvarsA[j][i]*dvUUINTprev[j]+(2*iterationNumber-dvarsBM[j][i])*dvUUINT[j] - dvarsC[j][i]*dvUUINTnext[j];
				}
				else {
					cout << __LINE__ << endl; exit(1);
					rhs[j]=(2*iterationNumber - mortality[i][j])*dvUUINT[j];
				}
			}

			// appel du sous programme "tridag" pour la seconde moitie
			// de l'intervalle de temps 
			tridag2(dvarsD[i],Ybet[i],dvarsF[i],rhs,uu[i],gam,jmin,jmax);
		}
	}
} // FIN DU SOUS PROGRAMME CALREC

void CCalpop::calrecJuv(const PMap& map, dvar_matrix& uu, dvar_matrix& mortality)
{	// DEBUT DU SOUS PROGRAMME CALREC
	
	dvar_vector uvec(0, maxn - 1);
	dvar_vector rhs(0, maxn - 1);
	dvar_vector gam(0, maxn - 1);

	for (int itr = 1; itr <= iterationNumber; itr++) 
	{
		//////////////////////////////////////////////////////////////////////////////////
		// 1- CALCUL DE LA PARTIE DROITE (RHS) DE L'EQUATION POUR L'INTERVALLE DE TEMPS 0 A 1/2
		//    RHS[] correspond au vecteur g[] de Sibert 
		/////////////////////////////////////////////////////////////////////////////////
		for (int j = map.jmin; j <= map.jmax; j++)
		{
			const int imin = map.iinf[j]; 
			const int imax = map.isup[j];
			for (int i = imin; i <= imax; i++)
			{   
				dvar_vector& dvUU = uu[i];
				if (map.jinf[i] <= j && j <= map.jsup[i]) {
					rhs[i] = -d[i][j]*dvUU[j-1] + (2*iterationNumber-e[i][j])*dvUU[j] - f[i][j]*dvUU[j+1];
				}
				else {
					cout << __LINE__ << endl; exit(1);

					rhs[i] = 2*iterationNumber*dvUU[j];
				}
			}

			// appel du sous programme tridag pour la premiere moitie de l'intervalle de temps
			tridag1_0(a[j],Xbet[j],c[j],rhs,uvec,gam,imin,imax);

			for (int i = imin; i <= imax; i++)
				Uuint[i][j]=uvec[i];
		}

		/////////////////////////////////////////////////////////////////////////
		// calcul de la partie droite de l'equation rrhs(i,j) pour 
		// l'intervalle de temps 1/2 a 1
		for (int i = map.imin; i <= map.imax; i++)
		{
			dvar_vector& dvUUINTprev = Uuint[i-1];
			dvar_vector& dvUUINT     = Uuint[i];
			dvar_vector& dvUUINTnext = Uuint[i+1];

			const int jmin = map.jinf[i];
			const int jmax = map.jsup[i];
			for (int j = jmin; j <= jmax; j++) {
				if (map.iinf[j] <= i && i <= map.isup[j]) {
					rhs[j]=-a[j][i]*dvUUINTprev[j]+(2*iterationNumber-dvarsBM[j][i])*dvUUINT[j] - c[j][i]*dvUUINTnext[j];
				}
				else {
					cout << __LINE__ << endl; exit(1);
					rhs[j]=(2*iterationNumber - mortality[i][j])*dvUUINT[j];
				}
			}

			// appel du sous programme "tridag" pour la seconde moitie
			// de l'intervalle de temps 
			tridag1_1(d[i],ybet[i],f[i],rhs,uu[i],gam,jmin,jmax);
		}
	}
} // FIN DU SOUS PROGRAMME CALREC

