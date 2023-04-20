#include "calpop.h"

/////////////////////////////////////////////////////////////////
//     -----------------------------------------------------------
//     sous programme 'calrec' : resolution du systeme d'equations
//     sur chaque demi interval de temps.
//     -----------------------------------------------------------
//     itr=boucle sur les sub-iterations
//     la partie droite de l'equation est d'abord stockee dans le vecteur RHS[i]
//     ou [j]. Sur la premiere moitie de l'intervalle de temps, j systemes de i
//     equations sont resolus: j fois D[i]*U[i]=RHS[i] Sur la seconde moitie de
//     l'intervalle de temps, i systemes de j equations sont resolus: i fois
//     D[j]*U[j]=RHS[j] les i solutions UVEC[j] sont stockees dans UU[i][j]
//     i*(UVEC[j])=>UU[i][j]
//      ---------------------------------------------------------
/////////////////////////////////////////////////////////////////

void CCalpop::calrec(
    const PMap& map, DMATRIX& uu,
    const DMATRIX& mortality) {  // DEBUT DU SOUS PROGRAMME CALREC

    DVECTOR uvec(0, maxn - 1);
    DVECTOR rhs(0, maxn - 1);
    DVECTOR gam(0, maxn - 1);

    ////////////////////////////////////////
    // début de la boucle des sous-itérations
    // Resolution des systemes tridiagonaux d'equations lineaires
    // suivant la methode ADI (Press et al. Numerical Recipes)
    for (int itr = 1; itr <= iterationNumber; itr++) {
        //////////////////////////////////////////////////////////////////////////////////
        // 1- CALCUL DE LA PARTIE DROITE (RHS) DE L'EQUATION POUR L'INTERVALLE
        // DE TEMPS 0 A 1/2
        //    RHS[] correspond au vecteur g[] de Sibert
        /////////////////////////////////////////////////////////////////////////////////
        for (int j = map.jmin; j <= map.jmax; j++) {
            const int imin = map.iinf[j];
            const int imax = map.isup[j];
            for (int i = imin; i <= imax; i++) {
                DVECTOR& dvUU = uu[i];
                if (map.jinf[i] <= j && j <= map.jsup[i]) {
                    rhs[i] = -d[i][j] * dvUU[j - 1] +
                             (2 * iterationNumber - e[i][j]) * dvUU[j] -
                             f[i][j] * dvUU[j + 1];
                } else {
                    cout << __LINE__ << endl;
                    exit(1);
                    rhs[i] = 2 * iterationNumber * dvUU[j];
                }
            }

            // appel du sous programme tridag pour la premiere moitie de
            // l'intervalle de temps
            tridag(a[j], xbet[j], c[j], rhs, uvec, gam, imin, imax);

            for (int i = imin; i <= imax; i++) uuint[i][j] = uvec[i];
        }

        /////////////////////////////////////////////////////////////////////////
        // calcul de la partie droite de l'equation rrhs(i,j) pour
        // l'intervalle de temps 1/2 a 1
        for (int i = map.imin; i <= map.imax; i++) {
            DVECTOR& dvUUINTprev = uuint[i - 1];
            DVECTOR& dvUUINT = uuint[i];
            DVECTOR& dvUUINTnext = uuint[i + 1];

            const int jmin = map.jinf[i];
            const int jmax = map.jsup[i];
            for (int j = jmin; j <= jmax; j++) {
                if (map.iinf[j] <= i && i <= map.isup[j]) {
                    rhs[j] = (-a[j][i] * dvUUINTprev[j]) +
                             (2 * iterationNumber - bm[j][i]) * dvUUINT[j] -
                             (c[j][i] * dvUUINTnext[j]);
                } else {
                    cout << __LINE__ << endl;
                    exit(1);
                    rhs[j] =
                        ((2 * iterationNumber - mortality[i][j])) * dvUUINT[j];
                }
            }

            // appel du sous programme "tridag" pour la seconde moitie
            // de l'intervalle de temps
            tridag(d[i], ybet[i], f[i], rhs, uu[i], gam, jmin, jmax);
        }
    }
}

void CCalpop::calrec_GO(const PMap& map, dvar_matrix& uu) {
    // cout << __LINE__ << " ADI-solver for GLOBAL domain " << endl;
    DVECTOR uvec(0, maxn - 1);
    DVECTOR rhs(0, maxn - 1);
    DVECTOR gam(0, maxn - 1);

    int ncells = 20;
    int nti = map.imax + 2;
    int ntj = map.jmax + 2;

    // cout << nti << " " << ntj << endl; //exit(1);
    // for (int j = map.jmin; j <= map.jmax; j++)
    //	cout << j << " " << map.iinf[j] << " "<< map.isup[j] << endl;
    // cout << map.imax << endl; exit(1);

    DVECTOR uvec1(0, nti + 2 * ncells - 3);
    DVECTOR rhs1(0, nti + 2 * ncells - 3);
    DVECTOR gam1(0, nti + 2 * ncells - 3);
    uvec1.initialize();
    rhs1.initialize();
    gam1.initialize();

    dmatrix a1(0, ntj - 1, 0, 2 * nti - 1);
    dmatrix xbet1(0, ntj - 1, 0, 2 * nti - 1);
    dmatrix c1(0, ntj - 1, 0, 2 * nti - 1);
    a1.initialize();
    xbet1.initialize();
    c1.initialize();

    for (int itr = 1; itr <= iterationNumber; itr++) {
        for (int j = map.jmin; j <= map.jmax; j++) {
            const int imin = map.iinf[j];
            const int imax = map.isup[j];
            for (int i = imin; i <= imax; i++) {
                rhs[i] = -d[i][j] * value(uu(i, j - 1)) +
                         (2 * iterationNumber - e[i][j]) * value(uu(i, j)) -
                         f[i][j] * value(uu(i, j + 1));
            }

            if ((map.carte[1][j]) && (map.carte[nti - 2][j])) {
                for (int k = 0; k < ncells; k++) {
                    a1[j][k] = a[j][nti - (ncells + 1) + k];
                    xbet1[j][k] =
                        bm[j][nti - (ncells + 1) + k] + 2 * iterationNumber;
                    c1[j][k] = c[j][nti - (ncells + 1) + k];
                }
                for (int i = ncells; i < nti - 2 + ncells; i++) {
                    a1[j][i] = a[j][i - ncells + 1];
                    xbet1[j][i] = bm[j][i - ncells + 1] + 2 * iterationNumber;
                    c1[j][i] = c[j][i - ncells + 1];
                }
                for (int k = nti - 2 + ncells; k < nti - 2 + 2 * ncells; k++) {
                    a1[j][k] = a[j][k - (nti - 3) - ncells];
                    xbet1[j][k] =
                        bm[j][k - (nti - 3) - ncells] + 2 * iterationNumber;
                    c1[j][k] = c[j][k - (nti - 3) - ncells];
                }

                for (int k = 0; k < ncells; k++)
                    rhs1[k] = rhs[nti - (ncells + 1) + k];
                for (int i = ncells; i < nti - 2 + ncells; i++)
                    rhs1[i] = rhs[i - ncells + 1];
                for (int k = nti - 2 + ncells; k < nti - 2 + 2 * ncells; k++)
                    rhs1[k] = rhs[k - (nti - 3) - ncells];

                tridag_GO(
                    a1[j], xbet1[j], c1[j], rhs1, uvec1, gam1, 0,
                    (nti - 2) + 2 * ncells - 1);

                for (int i = imin; i <= imax; i++)
                    uuint[i][j] = uvec1[i + ncells - 1];

                // and these values will be used in i loop below
                uuint[0][j] = uvec1[imax + ncells - 1];
                uuint[nti - 1][j] = uvec1[ncells];

            } else {
                tridag(a[j], xbet[j], c[j], rhs, uvec, gam, imin, imax);
                for (int i = imin; i <= imax; i++) uuint[i][j] = uvec[i];
            }
        }

        for (int i = map.imin; i <= map.imax; i++) {
            DVECTOR& dvUUINTprev = uuint[i - 1];
            DVECTOR& dvUUINT = uuint[i];
            DVECTOR& dvUUINTnext = uuint[i + 1];

            const int jmin = map.jinf[i];
            const int jmax = map.jsup[i];
            for (int j = jmin; j <= jmax; j++) {
                rhs[j] = (-a[j][i] * dvUUINTprev[j]) +
                         (2 * iterationNumber - bm[j][i]) * dvUUINT[j] -
                         (c[j][i] * dvUUINTnext[j]);
            }

            tridag(d[i], ybet[i], f[i], rhs, uvec, gam, jmin, jmax);
            for (int j = jmin; j <= jmax; j++) uu.elem_value(i, j) = uvec[j];
        }
    }
    return;
}

void CCalpop::calrec_GO_with_catch(
    const PMap& map, CParam& param, dvar_matrix& uu, const dmatrix& C_obs,
    dvar_matrix& C_est) {
    // cout << __LINE__ << " ADI-solver for GLOBAL domain " << endl;
    DVECTOR uvec(0, maxn - 1);
    DVECTOR rhs(0, maxn - 1);
    DVECTOR gam(0, maxn - 1);

    int ncells = 20;
    int nti = map.imax + 2;
    int ntj = map.jmax + 2;

    // cout << nti << " " << ntj << endl; //exit(1);
    // for (int j = map.jmin; j <= map.jmax; j++)
    //	cout << j << " " << map.iinf[j] << " "<< map.isup[j] << endl;
    // cout << map.imax << endl; exit(1);

    DVECTOR uvec1(0, nti + 2 * ncells - 3);
    DVECTOR rhs1(0, nti + 2 * ncells - 3);
    DVECTOR gam1(0, nti + 2 * ncells - 3);
    uvec1.initialize();
    rhs1.initialize();
    gam1.initialize();

    dmatrix a1(0, ntj - 1, 0, 2 * nti - 1);
    dmatrix xbet1(0, ntj - 1, 0, 2 * nti - 1);
    dmatrix c1(0, ntj - 1, 0, 2 * nti - 1);
    a1.initialize();
    xbet1.initialize();
    c1.initialize();

    for (int itr = 1; itr <= iterationNumber; itr++) {
        for (int j = map.jmin; j <= map.jmax; j++) {
            const int imin = map.iinf[j];
            const int imax = map.isup[j];
            for (int i = imin; i <= imax; i++) {
                rhs[i] = -d[i][j] * value(uu(i, j - 1)) +
                         (2 * iterationNumber - e[i][j]) * value(uu(i, j)) -
                         f[i][j] * value(uu(i, j + 1));
            }

            if ((map.carte[1][j]) && (map.carte[nti - 2][j])) {
                for (int k = 0; k < ncells; k++) {
                    a1[j][k] = a[j][nti - (ncells + 1) + k];
                    xbet1[j][k] =
                        bm[j][nti - (ncells + 1) + k] + 2 * iterationNumber;
                    c1[j][k] = c[j][nti - (ncells + 1) + k];
                }
                for (int i = ncells; i < nti - 2 + ncells; i++) {
                    a1[j][i] = a[j][i - ncells + 1];
                    xbet1[j][i] = bm[j][i - ncells + 1] + 2 * iterationNumber;
                    c1[j][i] = c[j][i - ncells + 1];
                }
                for (int k = nti - 2 + ncells; k < nti - 2 + 2 * ncells; k++) {
                    a1[j][k] = a[j][k - (nti - 3) - ncells];
                    xbet1[j][k] =
                        bm[j][k - (nti - 3) - ncells] + 2 * iterationNumber;
                    c1[j][k] = c[j][k - (nti - 3) - ncells];
                }

                for (int k = 0; k < ncells; k++)
                    rhs1[k] = rhs[nti - (ncells + 1) + k];
                for (int i = ncells; i < nti - 2 + ncells; i++)
                    rhs1[i] = rhs[i - ncells + 1];
                for (int k = nti - 2 + ncells; k < nti - 2 + 2 * ncells; k++)
                    rhs1[k] = rhs[k - (nti - 3) - ncells];

                tridag_GO(
                    a1[j], xbet1[j], c1[j], rhs1, uvec1, gam1, 0,
                    (nti - 2) + 2 * ncells - 1);

                for (int i = imin; i <= imax; i++)
                    uuint[i][j] = uvec1[i + ncells - 1];

                // and these values will be used in i loop below
                uuint[0][j] = uvec1[imax + ncells - 1];
                uuint[nti - 1][j] = uvec1[ncells];

            } else {
                tridag(a[j], xbet[j], c[j], rhs, uvec, gam, imin, imax);

                for (int i = imin; i <= imax; i++) {
                    //					uuint[i][j] = uvec[i];
                    if (C_obs(i, j) == 0)
                        uuint(i, j) = uvec[i];
                    else {
                        double C_est_ij_itr =
                            uvec[i] *
                            param.func_limit_one(
                                C_obs(i, j) / (uvec[i] + 1e-14)) /
                            iterationNumber;
                        uuint(i, j) = uvec[i] - C_est_ij_itr;
                        C_est.elem_value(i, j) += C_est_ij_itr;
                    }
                    if (uvec(i) > 0 && uuint(i, j) < 0) {
                        cout << "Negative values due to catch removal in"
                             << "i=" << i << ", j=" << j << ", B=" << uvec(i)
                             << ", C=" << C_obs(i, j) << endl;
                    }
                }
            }
        }

        for (int i = map.imin; i <= map.imax; i++) {
            DVECTOR& dvUUINTprev = uuint[i - 1];
            DVECTOR& dvUUINT = uuint[i];
            DVECTOR& dvUUINTnext = uuint[i + 1];

            const int jmin = map.jinf[i];
            const int jmax = map.jsup[i];
            for (int j = jmin; j <= jmax; j++) {
                rhs[j] = (-a[j][i] * dvUUINTprev[j]) +
                         (2 * iterationNumber - bm[j][i]) * dvUUINT[j] -
                         (c[j][i] * dvUUINTnext[j]);
            }

            tridag(d[i], ybet[i], f[i], rhs, uvec, gam, jmin, jmax);
            for (int j = jmin; j <= jmax; j++) uu.elem_value(i, j) = uvec[j];
        }
    }
    return;
}

// FIN DU SOUS PROGRAMME CALREC
/////////////////////////////////////////////////////////////////////////
