#include "SaveTimeArea.h"

// #ifdef _DEBUG
// #undef THIS_FILE
// static char THIS_FILE[]=__FILE__;
// #define new DEBUG_NEW
// #endif

void CSaveTimeArea::SumByArea(
    const PMap& map, const dmatrix& mask_catch, const dmatrix& mat2d,
    dvector& sum_area, const dvector cell_area, const int nb_reg,
    const int nbt) {
    sum_area.initialize();

    // double cr_catch = 0.01*max(mask_catch);
    double cr_catch = 0.1 * nbt;
    for (int a = 0; a < nb_reg; a++) {
        for (int i = map.regimin[a]; i < map.regimax[a]; i++) {
            for (int j = map.regjmin[a]; j < map.regjmax[a]; j++) {
                if (map.carte(i, j) && mask_catch(i, j) > cr_catch) {
                    sum_area[a] += mat2d[i][j] * cell_area[j];
                }
            }
        }
        // total all areas
        sum_area[nb_reg] += sum_area[a];
    }
}

void CSaveTimeArea::SumByEEZ(
    const CParam& param, const PMap& map, const DMATRIX& mat2d,
    DVECTOR& sum_EEZ, const dvector cell_area) {
    sum_EEZ.initialize();
    const int nlon = param.nlong;
    const int nlat = param.nlat;

    for (int i = 0; i < nlon; i++)
        for (int j = 0; j < nlat; j++)
            if (map.carte(i, j))
                for (int a = 0; a < param.nb_EEZ; a++)
                    if (map.maskEEZ[i][j] == param.EEZ_ID[a])
                        sum_EEZ[a] += mat2d[i][j] * cell_area[j];
}

double CSaveTimeArea::SumByEEZ(
    const PMap& map, const int EEZ_ID, const DMATRIX& mat2d,
    const dvector cell_area, const int nlon, const int nlat) {
    double sum_EEZ = 0.0;
    for (int i = 0; i < nlon; i++)
        for (int j = 0; j < nlat; j++)
            if (map.carte(i, j))
                if (map.maskEEZ[i][j] == EEZ_ID)
                    sum_EEZ += mat2d[i][j] * cell_area[j];

    return sum_EEZ;
}

double CSaveTimeArea::SumByEEZ(
    const PMap& map, const int EEZ_ID, const DMATRIX& C, const DMATRIX& E,
    const int nlon, const int nlat) {
    double sum_EEZ = 0.0;
    for (int i = 0; i < nlon; i++)
        for (int j = 0; j < nlat; j++)
            if (map.carte(i, j))
                if (map.maskEEZ[i][j] == EEZ_ID)
                    if (E[i][j] > 0) {
                        sum_EEZ += C[i][j] / E[i][j];
                    }

    return sum_EEZ;
}

double CSaveTimeArea::StdCPUEByEEZ(
    const PMap& map, const int EEZ_ID, const DMATRIX& C, const DMATRIX& E,
    const double mean, const int nobs, const int nlon, const int nlat) {
    double sum_EEZ = 0.0;
    for (int i = 0; i < nlon; i++)
        for (int j = 0; j < nlat; j++)
            if (map.carte(i, j))
                if (map.maskEEZ[i][j] == EEZ_ID)
                    if (E[i][j] > 0) {
                        double val = C[i][j] / E[i][j];
                        sum_EEZ += pow(val - mean, 2);
                    }
    double std_var = 0;
    if (nobs > 1) std_var = sqrt(sum_EEZ / (nobs - 1));

    return std_var;
}

int CSaveTimeArea::NobsByEEZ(
    const PMap& map, const int EEZ_ID, const DMATRIX& E, const int nlon,
    const int nlat) {
    double sum_EEZ = 0.0;
    for (int i = 0; i < nlon; i++)
        for (int j = 0; j < nlat; j++)
            if (map.carte(i, j))
                if (map.maskEEZ[i][j] == EEZ_ID)
                    if (E[i][j] > 0) {
                        sum_EEZ++;
                    }
    return sum_EEZ;
}
