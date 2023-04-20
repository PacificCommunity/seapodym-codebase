#include "calpop.h"

/// Forward main function called in simulation mode only for:
/// calrec for larval and juvenile life stages, i.e. with passive drift only.
/// See calrec_adre.cpp

void CCalpop::Calrec_juv(
    const PMap& map, CMatrices& mat, dvar_matrix& uu, dvar_matrix& mortality,
    const int t_count) {
    bm = value(dvarsBM);
    xbet = value(Xbet);
    dmatrix mort_c = value(mortality);
    calrec1(map, uu, mort_c);
}
