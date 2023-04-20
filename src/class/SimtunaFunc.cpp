// SimtunaFunc.cpp: implementation of the CSimtunaFunc class.
//
//////////////////////////////////////////////////////////////////////

// #include "StdAfx.h"
#include "SimtunaFunc.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#define new DEBUG_NEW
#endif

double CSimtunaFunc::function_lambda(
    CParam& param, CMatrices& mat, int n, int i, int j) {
    // return the value of lambda (natural mortality of forage)
    // in relation to the temperature attributed to the forage component
    // non migrant : based on temperature in the layer
    // migrant : based on the average temperature relative to time spent in day
    // and night layers
    double function;
    double a = param.inv_lambda_curv;
    int max = param.inv_lambda_max;

    double temp = 0.f;

    if (param.day_layer[n] == param.night_layer[n])
        temp = mat.tempn[1][param.day_layer[n]][i][j];
    else {  // average temperature according to time passed in both layers
        double DL = mat.daylength[i][j];
        temp = ((DL * mat.tempn[1][param.day_layer[n]][i][j]) +
                ((24 - DL) * mat.tempn[1][param.night_layer[n]][i][j])) /
               24;
    }

    if (temp > 0) {
        function = (double)(max * exp(a * temp));
        // retourne lambda en fonction du pas de temps
        function = 1 / (double)(function / param.deltaT);
    } else
        function = (double)(param.deltaT / max);

    return function;
}

double CSimtunaFunc::daylength_twilight(
    double lat, int jday,
    const double p) {  // The CBM model of Forsythe et al, Ecological Modelling
                       // 80 (1995) 87-95 p - angle between the sun position and
                       // the horizon, in degrees 6  - civil twilight 12 -
                       // nautical twilight 18 - astronomical twilight

    const double pi = 3.14159265358979;

    // revolution angle for the day of the year
    double theta =
        0.2163108 + 2 * (atan(0.9671396 * tan(0.00860 * (jday - 186))));

    // sun's declination angle, or the angular distance at solar noon between
    // the Sun and the equator, from the Eartch orbit revolution angle
    double phi = asin(0.39795 * cos(theta));

    // daylength computed according to 'p'
    double arg = (sin(pi * p / 180) + sin(lat * pi / 180) * sin(phi)) /
                 (cos(lat * pi / 180) * cos(phi));
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;
    double DL = 24.0 - (24.0 / pi) * acos(arg);

    return DL;
}

double CSimtunaFunc::daylength(double lat, int jday) {
    //-------------------------------------------------------------------
    //  New function provided by Laurent Bopp as used in the PISCES model
    //
    //  PURPOSE :compute the day length depending on latitude and the day
    //  --------
    //
    //   MODIFICATIONS:
    //   --------------
    //      original  : E. Maier-Reimer (GBC 1993)
    //      additions : C. Le Quere (1999)
    //      modifications : O. Aumont (2004)
    //	Adapted to C : P. Lehodey (2005)
    //-------------------------------------------------------------------

    double DL = 0.0;
    double pi = 3.1415926536;

    double rum = 0.0;
    double delta = 0.0;
    double codel = 0.0;
    double phi = 0.0;
    double argu = 0.0;

    rum = (jday - 80) / 365.25;
    delta = sin(rum * pi * 2) * sin(pi * 23.5 / 180);
    codel = asin(delta);
    phi = lat * pi / 180;
    argu = tan(codel) * tan(phi);

    argu = min(1., argu);
    argu = max(-1., argu);
    DL = 24.0 - 2.0 * acos(argu) * 180.0 / pi / 15;
    DL = max(DL, 0.0);

    return DL;
}

double CSimtunaFunc::grad_daylength(double lat, int jday) {
    // return the true monthly gradient of length of day (DL) for a given
    // latitude (lat)
    double DL = 0.0;
    double DL_before = 0.0;
    double grad_DL = 0.0;
    int jdaybefore;

    if (jday == 1)
        jdaybefore = 365;
    else
        jdaybefore = jday - 1;

    DL = daylength(lat, jday);
    DL_before = daylength(lat, jdaybefore);
    grad_DL = DL - DL_before;

    return grad_DL;
}

double CSimtunaFunc::f_accessibility_comp(
    const double Od, const double On, const double Td, const double Tn,
    double twosigsq, double temp_mean, double oxy_teta, double oxy_cr,
    const double DL) {
    double f_oxyday = 1.0 / (1.0 + pow(oxy_teta, Od - oxy_cr));
    double f_oxynight = 1.0 / (1.0 + pow(oxy_teta, On - oxy_cr));

    double f_oxy = f_oxyday * DL + f_oxynight * (1.0 - DL);

    double f_tempday = exp(-pow(Td - temp_mean, 2.0) / twosigsq);
    double f_tempnight = exp(-pow(Tn - temp_mean, 2.0) / twosigsq);

    double f_temp = f_tempday * DL + f_tempnight * (1.0 - DL);

    double f_access = f_temp * f_oxy;

    return f_access;
}
