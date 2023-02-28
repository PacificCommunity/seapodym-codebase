// Numfunc.cpp: implementation of the CNumfunc class.
//
//////////////////////////////////////////////////////////////////////

//#include "StdAfx.h"
#include "Numfunc.h"


//#ifdef _DEBUG
//#undef THIS_FILE
//static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
//#endif

//#include "iostream.h"
//#include "fstream.h"   
//#include "nrutil.h"       
//#include <math.h>

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define TINY 1.0e-20		// Will regularize the unusual case of complete correlation


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CNumfunc::CNumfunc()
{

}

CNumfunc::~CNumfunc()
{

}


// fonction  de calcul de correlation entre les captures observees (xx) et predites (yy)
// les zero ne sont pas pris en compte dans le calcul

void CNumfunc::corcatch(DMATRIX& xx,DMATRIX& yy, const int imin, const int imax, const ivector jinf, const ivector jsup, int &nn, double &cor, double &z, double &prob, const IMATRIX& mask, double missval) 
{

	for (int i=imin; i <=imax; i++)
	  for (int j=jinf[i]; j<=jsup[i]; j++)
		if ((mask[i][j]==0) || (xx[i][j]==missval))
		{
			xx[i][j]=0;
			yy[i][j]=0;
		}

	corlin(xx,yy,imin,imax,jinf,jsup,nn,cor,z,prob,0);
}

void CNumfunc::corcpue(DMATRIX& xx,DMATRIX& yy,dmatrix& eff, const int imin, const int imax, const ivector jinf, const ivector jsup, int &nn, double &cor, double &z, double &prob, const IMATRIX& mask, double missval)
{
	dmatrix xxx, yyy;
	xxx.allocate(imin,imax,jinf,jsup); xxx.initialize();
	yyy.allocate(imin,imax,jinf,jsup); yyy.initialize();
        for (int i=imin; i <=imax; i++)
          for (int j=jinf[i]; j<=jsup[i]; j++)
                if (mask(i,j)&&(eff(i,j)!=0)){
			xxx(i,j) = xx(i,j)/eff(i,j);
			yyy(i,j) = yy(i,j)/eff(i,j);
		}

        corlin(xxx,yyy,imin,imax,jinf,jsup,nn,cor,z,prob,0);
}

void CNumfunc::corlin(const DMATRIX& xx, const DMATRIX& yy, const int imin, const int imax, const ivector jinf, const ivector jsup, int &nn, double &cor, double &z, double &prob, double missval) 
{
	
	sx=0.f;	sy=0.f;	sxx=0.f;	syy=0.f;	sxy=0.f;
	cor=0.f;	z=0.f;	prob=1.f;	nn=0 ;
	
	summat(xx,sx,imin,imax,jinf,jsup,nn,missval);
	summat(yy,sy,imin,imax,jinf,jsup,nn,missval);

	if ((sx != 0) && (sy != 0))
	{
		sumdif(xx,yy,imin,imax,jinf,jsup,nn,sx,sy,sxx,syy,sxy,missval);
		pearsn(sxx,syy,sxy,nn,cor,prob,z);
	}
}

// fonction qui retourne la somme d'une matrice 2D
void CNumfunc::summat(const DMATRIX& xx, double &sx, const int imin, const int imax, const ivector jinf, const ivector jsup, int &nn, const double missval)
{
	nn=0;
	
        for (int i=imin; i<=imax; i++)
	        for (int j=jinf[i]; j<=jsup[i]; j++)
			if (xx[i][j]!=0)
			{
				sx += xx[i][j];
				nn++;
			}
}
    
// fonction qui retourne la somme des carres des differences a la moyenne de x et y
// et le produit des differences a la moyenne de x et y 

     
void CNumfunc::sumdif(const DMATRIX& xx, const DMATRIX& yy, const int imin, const int imax, const ivector jinf, const ivector jsup, const int nn, double sx, double sy, double &sxx, double &syy, double &sxy, const double missval)
{
	sx /= nn;	// la somme est divisee par le nombre de valeurs
	sy /= nn;	// et devient la moyenne
	    
        for (int i=imin; i<=imax; i++)
        	for (int j=jinf[i]; j<=jsup[i]; j++)
			if (xx[i][j]!=0)
			{
				double xt=xx[i][j]-sx;
				double yt=yy[i][j]-sy;
				sxx += xt*xt;
				syy += yt*yt;
				sxy += xt*yt;
			}
}
          
// this routine computes the correlation coefficient r (returned as r),
// the significance level at which the null hypothesis of zero
// correlation is disproved (prob whose small value indicates a significant correlation),
// and Fisher's z (returned as z), whose value can be used in further statistical tests 
// as described in 'Numerical Recipes in C, Press and al. 1994 (14.5 p. 638) 

void CNumfunc::pearsn(const double sxx, const double syy, const double sxy, const int n, double &r, double &prob, double &z)
{
//	double betai(double a, double b, double x);
//	double erfcc(double x);

	r= (double)(sxy/(sqrt(sxx*syy)));
	z= (double) (0.5*log((1.0+(r)+TINY)/(1.0-(r)+TINY)));	//Fisher'sz transformation
	double df=n-2;
	double t=(r)*sqrt(df/((1.0-(r)+TINY)*(1.0+(r)+TINY)));	// Equation (14.5.5)
	prob= (double)  ( betai(0.5*df,0.5, df/(df+t*t)) );	// Student's t probability */

// for large n, the easier computation of prob below, using the short routine erfcc,
// would give approximatively the same value		

//		prob=erfcc(fabs((z)*sqrt(n-1.0))/1.4142136);		

}


// function which returns the value ln[gamma(xx)] for xx > 0
// Internal arithmetic will be done in double precision, a nicety that you can omit
// if five-figure accuracy is good enough	Ref.: Numerical Recipes p.214
              
double CNumfunc::gammln(const double xx)
{
	static double cof[6]={76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};

	double x=xx;
	double y=xx;
	double tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	double ser=1.000000000190015;
	for (int j=0; j<=5; j++) ser += cof[j] /++y;

	return -tmp+log(2.5066282746310005*ser/x);
}                                                                        

// Function which evaluates continuated fraction for incomplete beta function
// by modified Lentz's method (Numerical Recipes in C 5.2)
// used in function betai
// Ref.: Numerical Recipes p.227


double CNumfunc::betacf(double a, double b, double x)
{
	//void nrerror(char error_text[]);
	int m, m2;
	double aa, c, d, del, h, qab, qam, qap;
   
	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1; m <= MAXIT; m++)
	{
		m2=2*m;
		aa=m*(b-m)*x / ((qam+m2) * (a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x / ((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}

	//if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}


// Return the incomplete beta function Ix(a,b).  Ref.: Numerical Recipes p.227
// using the functions betacf and gammln

double CNumfunc::betai(double a, double b, double x)
{
//		double betacf(double a, double b, double x);
//		double gammln(double xx);
//		void nrerror(char error_text[]);
		double bt;
	
//		if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
		if (x == 0.0 || x == 1.0) bt=0.0;
		else
			bt=exp(gammln(a+b) - gammln(a) - gammln(b) + a*log(x) + b*log(1.0-x));
		if (x < (a+1.0) / (a + b + 2.0))
			return bt*betacf(a,b,x)/a;
		else
			return 1.0 - bt*betacf(b,a,1.0-x)/b;
}

// Return the complementary error function erfc(x) with fractional error 
// everywhere less than 1.2 x 10-7
/*	
double CNumfunc::erfcc(double x)
{
		double t, z, ans;
				
		z=fabs(x);
		t=1.0/(1.0+0.5*z);
		ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
			t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
			t*(-0.82215223+t*0.17087277)))))))));
		return x >= 0.0 ? ans : 2.0-ans;
		
}
*/	
