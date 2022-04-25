// SimtunaFunc.h: interface for the CSimtunaFunc class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __SimtunaFunc_h__
#define __SimtunaFunc_h__


#include "Matrices.h"

/*!
\brief The simulation function which do not use dvariables
*/

class CSimtunaFunc  
{
public:
	CSimtunaFunc() {/*DoNothing*/};
	virtual ~CSimtunaFunc() {/*DoNothing*/};

public:
 
	double function_lambda(CParam& param, CMatrices& mat, int n, int i, int j);
	double daylength(double lat, int jday);
	double daylength_twilight(double lat, int jday, const double p);
	double grad_daylength(double lat, int jday);
	double f_accessibility_comp(const double Od, const double On, const double Td, const double Tn, 
				double twosigsq, double temp_mean, double oxy_teta, double oxy_cr, const double DL);

};

#endif 
