#include "VarSimtunaFunc.h"

/// Forward main function called in simulation mode only for:
/// computing the seasonal switch function used to switch between habitats.
/// See seasonal_switch.cpp

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
}
