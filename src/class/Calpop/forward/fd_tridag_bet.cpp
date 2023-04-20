#include "calpop.h"

/// Forward main function called in simulation mode only for:
/// tridag_bet function. See tridag_bet.cpp

void CCalpop::Xbet_comp1(const PMap& map, int dt) {
    dmatrix xbet_c = value(Xbet);
    dmatrix bm_c = value(dvarsBM);

    xbet_comp(map, xbet_c, a, bm_c, c, dt);

    Xbet = nograd_assign(xbet_c);
}
