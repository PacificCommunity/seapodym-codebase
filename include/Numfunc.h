// Numfunc.h: interface for the CNumfunc class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __Numfunc_h__
#define __Numfunc_h__

#include "mytypes.h"

/*!
\brief Class for computing various mathematical functions.

*/

class CNumfunc {
   public:
    CNumfunc();
    virtual ~CNumfunc();

    double sx, sy, sxx, syy, sxy;
    double nn, missval;

    // Member functions - Fonctions membres

    void corcatch(
        DMATRIX &xx, DMATRIX &yy, const int imin, const int imax,
        const ivector jinf, const ivector jsup, int &nn, double &cor, double &z,
        double &prob, const IMATRIX &mask, double missval);

    void corcpue(
        DMATRIX &xx, DMATRIX &yy, dmatrix &eff, const int imin, const int imax,
        const ivector jinf, const ivector jsup, int &nn, double &cor, double &z,
        double &prob, const IMATRIX &mask, double missval);

    void corlin(
        const DMATRIX &xx, const DMATRIX &yy, const int imin, const int imax,
        const ivector jinf, const ivector jsup, int &nn, double &cor, double &z,
        double &prob, double missval);
    // fonction qui retourne la valeur de correlation lineaire entre deux
    // matrices

    void summat(
        const DMATRIX &xx, double &sx, const int imin, const int imax,
        const ivector jinf, const ivector jsup, int &nn, const double missval);
    // fonction qui retourne la somme d'une matrice 2D

    void sumdif(
        const DMATRIX &xx, const DMATRIX &yy, const int imin, const int imax,
        const ivector jinf, const ivector jsup, const int nn, double sx,
        double sy, double &sxx, double &syy, double &sxy, const double missval);
    // fonction qui retourne la somme des carres des differences a la moyenne de
    // x et y et le produit des differences a la moyenne de x et y

    double gammln(const double xx);
    // function which returns the value ln[gamma(xx)] for xx > 0
    // Internal arithmetic will be done in double precision, a nicety that you
    // can omit if five-figure accuracy is good enough	Ref.: Numerical Recipes
    // p.214

    double betacf(double a, double b, double x);
    // Function which evaluates continuated fraction for incomplete beta
    // function by modified Lentz's method (Numerical Recipes in C 5.2) used by
    // function betai Ref.: Numerical Recipes p.227

    double betai(double a, double b, double x);
    // Return the incomplete beta function Ix(a,b).  Ref.: Numerical Recipes
    // p.227 use the functions betacf and gammln

    // double erfcc(double x);
    // Return the complementary error function erfc(x) with fractional error
    // everywhere less than 1.2 x 10-7  Ref.: Numerical Recipes p.221
    // used by function pearsn

    void pearsn(
        const double sxx, const double syy, const double sxy, const int n,
        double &r, double &prob, double &z);
    // this routine computes the correlation coefficient r
    // (returned as corr), the significance level at which the null hypothesis
    // of zero correlation is disproved (prob whose small value indicates a
    // significant correlation), and Fisher's z (returned as z), whose value can
    // be used in further statistical tests as described in 'Numerical Recipes
    // in C, Press and al. 1994 (14.5 p. 638)

    double deplete(double fish, double f, double m);
};

#endif
