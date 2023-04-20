#ifndef __SaveTimeArea_h__
#define __SaveTimeArea_h__

#include "SimtunaFunc.h"
#include "mytypes.h"

/*!
\brief The class to aggregate variables to the regional structure.
*/
class CSaveTimeArea {
   public:
    CSaveTimeArea(){/*DoNothing*/};
    virtual ~CSaveTimeArea(){/*DoNothing*/};

    //	void SumByArea(const PMap& map, const dmatrix& mask_catch, const
    //dmatrix& mat2d, dvector& sum_area, const dvector cell_area_corr, const int
    //nb_reg, const int nbt, const int sp); 	void SumByEEZ(const CParam& param,
    //const PMap& map, const DMATRIX& mat2d, DVECTOR& sum_EEZ, const dvector
    //cell_area_corr);
    void SumByArea(
        const PMap& map, const dmatrix& mask_catch, const dmatrix& mat2d,
        dvector& sum_area, const dvector cell_area, const int nb_reg,
        const int nbt);
    void SumByEEZ(
        const CParam& param, const PMap& map, const DMATRIX& mat2d,
        DVECTOR& sum_EEZ, const dvector cell_area);
    double SumByEEZ(
        const PMap& map, const int EEZ_ID, const DMATRIX& mat2d,
        const dvector cell_area, const int nlon, const int nlat);
    int NobsByEEZ(
        const PMap& map, const int EEZ_ID, const DMATRIX& mat2d, const int nlon,
        const int nlat);
    double SumByEEZ(
        const PMap& map, const int EEZ_ID, const DMATRIX& C, const DMATRIX& E,
        const int nlon, const int nlat);
    double StdCPUEByEEZ(
        const PMap& map, const int EEZ_ID, const DMATRIX& C, const DMATRIX& E,
        const double mean, const int nobs, const int nlon, const int nlat);
};

#endif
