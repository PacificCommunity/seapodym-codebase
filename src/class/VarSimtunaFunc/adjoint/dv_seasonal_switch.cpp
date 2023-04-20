#include "VarSimtunaFunc.h"

/// Main function with memory control and adjoint functions for:
/// computing the seasonal switch function used to switch between habitats.
/// Forward functions are in seasonal_switch.cpp

double daylength_comp(double lat, double jday, double pi);
void dv_seasonal_switch_comp(void);
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);

void VarSimtunaFunc::Seasonal_switch(
    VarParamCoupled& param, VarMatrices& mat, const PMap& map, int jday,
    int sp) {
    dvariable season_peak = param.dvarsSpawning_season_peak[sp];
    dvariable season_start = param.dvarsSpawning_season_start[sp];

    dvar_matrix season_peak_mat(map.imin, map.imax, map.jmin, map.jmax);
    season_peak_mat = season_peak;
    dvar_matrix season_start_mat(map.imin, map.imax, map.jmin, map.jmax);
    season_start_mat = season_start;

    Seasonal_switch_comp(
        param, mat, map, value(season_peak), value(season_start), jday, sp);

    save_identifier_string2((char*)"season_switch_comp_begin");
    season_peak.save_prevariable_value();
    season_peak_mat.save_dvar_matrix_position();
    season_start.save_prevariable_value();
    season_start_mat.save_dvar_matrix_position();
    mat.dvarSeasonSwitch[sp].save_dvar_matrix_position();
    unsigned long int cparam = (unsigned long int)&param;
    save_long_int_value(cparam);
    unsigned long int pmap = (unsigned long int)&map;
    save_long_int_value(pmap);
    save_int_value(jday);
    save_identifier_string2((char*)"season_switch_comp_end");

    gradient_structure::GRAD_STACK1->set_gradient_stack(
        dv_seasonal_switch_comp);
}

void dffdaylength(
    double& dfDL, double& dfJday, double lat, double jday,
    double pi) {  // adjoint for daylength function(double,double,double)

    double rum = 0.0;
    double delta = 0.0;
    double codel = 0.0;
    double argu = 0.0;

    const double expr1 = sin(pi * 23.5 / 180.0);
    const double expr2 = tan(lat * pi / 180.0);

    // recomputation
    rum = (jday - 80.0) / 365.25;
    delta = sin(rum * pi * 2.0) * expr1;
    codel = asin(delta);
    argu = tan(codel) * expr2;

    double dfargu = 0.0;
    double dfcodel = 0.0;
    double dfdelta = 0.0;
    double dfrum = 0.0;

    // DL = 24.0-2.0*acos(argu)*180.0/(pi*15.0);
    dfargu += dfDL * 360.0 / (pi * 15.0 * sqrt(1.0 - argu * argu));
    dfDL = 0.0;

    // argu  = tan(codel)*tan(phi);
    dfcodel += dfargu * expr2 / (cos(codel) * cos(codel));

    // codel = asin(delta);
    dfdelta += dfcodel / sqrt(1.0 - delta * delta);

    // delta = sin(rum*pi*2)*sin(pi*23.5/180);
    dfrum += 2.0 * pi * cos(rum * pi * 2.0) * expr1 * dfdelta;

    // rum   = (jday-80)/365.25;
    dfJday += dfrum / 365.25;
}

void dv_seasonal_switch_comp() {
    verify_identifier_string2((char*)"season_switch_comp_end");
    unsigned jday = restore_int_value();
    unsigned long int pos_map = restore_long_int_value();
    unsigned long int pos_param = restore_long_int_value();
    const dvar_matrix_position Switch_pos = restore_dvar_matrix_position();
    const dvar_matrix_position start_pos = restore_dvar_matrix_position();
    double season_start = restore_prevariable_value();
    const dvar_matrix_position peak_pos = restore_dvar_matrix_position();
    double season_peak = restore_prevariable_value();
    verify_identifier_string2((char*)"season_switch_comp_begin");

    dmatrix dfPeak = restore_dvar_matrix_derivatives(peak_pos);
    dmatrix dfStart = restore_dvar_matrix_derivatives(start_pos);
    dmatrix dfSwitch = restore_dvar_matrix_derivatives(Switch_pos);

    CParam* param = (CParam*)pos_param;
    PMap* map = (PMap*)pos_map;

    const double pi = 4.0 * (atan(0.5) + atan(1.0 / 3.0));
    const int jday_maxDL = 171;  // jday of maximal DL
    const double tau = 50.0;

    // recompute
    double Jday = jday - 1 + jday_maxDL - season_peak;
    const int jmin = map->jmin;
    const int jmax = map->jmax;
    for (int j = jmax; j >= jmin; j--) {
        // recompute
        double lat = param->lastlat(j);
        double DL = daylength_comp(lat, Jday, pi) / 24.0;
        double dn_ratio = (DL) / (1.0 - DL);

        double expr1 = exp(-tau * (dn_ratio - season_start));
        double expr2 = (1.0 + expr1) * (1.0 + expr1);
        double expr3 = tau * expr1 / expr2;

        const int imin = map->iinf[j];
        const int imax = map->isup[j];
        for (int i = imax; i >= imin; i--) {
            if (map->carte(i, j)) {
                double dfDL = 0.0;
                double dfDN = 0.0;
                double dfJday = 0.0;

                // mat.SeasonSwitch(sp).elem_value(i,j) =
                // tetafunc(dn_ratio-season_start,tau2);
                dfStart(i, j) -= expr3 * dfSwitch(i, j);
                dfDN += expr3 * dfSwitch(i, j);
                dfSwitch(i, j) = 0.0;

                // double dn_ratio = (DL/24.0)/(1.0-(DL/24.0));
                // dfDL += dfDN/((24.0-DL)*(1-DL/24.0));
                dfDL += dfDN / ((1.0 - DL) * (1.0 - DL));

                // double DL = daylength(lat, Jday, pi);
                dffdaylength(dfDL, dfJday, lat, Jday, pi);
                // double Jday = i-1 + jday_maxDL - season_peak;
                dfPeak(i, j) -= dfJday / 24.0;
            }
        }
    }

    dfPeak.save_dmatrix_derivatives(peak_pos);
    dfStart.save_dmatrix_derivatives(start_pos);
    dfSwitch.save_dmatrix_derivatives(Switch_pos);
}
