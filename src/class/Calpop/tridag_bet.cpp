#include "calpop.h"

/// Forward functions for:
/// tridag_bet function. This routine precomputes an operator in
/// the Gaussian solver of the tridiagonal linear system, which
/// does not change during iterations.

void CCalpop::xbet_comp(
    const PMap& map, dmatrix& xbet, dmatrix& a, dmatrix& bm, dmatrix& c,
    int dt) {
    for (int j = map.jmin; j <= map.jmax; j++) {
        int inf = map.iinf[j];
        xbet[j][inf] = 1 / (bm[j][inf] + dt);
        for (int i = inf + 1; i <= map.isup[j]; i++) {
            double beta =
                bm[j][i] + dt - c[j][i - 1] * a[j][i] * xbet[j][i - 1];
            xbet[j][i] = 1 / beta;
        }
    }
}

void CCalpop::ybet_comp(
    const PMap& map, dmatrix& ybet, dmatrix& d, dmatrix& e, dmatrix& f,
    int dt) {
    for (int i = map.imin; i <= map.imax; i++) {
        int inf = map.jinf[i];
        ybet[i][inf] = 1 / (e[i][inf] + dt);
        for (int j = inf + 1; j <= map.jsup[i]; j++) {
            double beta = e[i][j] + dt - f[i][j - 1] * d[i][j] * ybet[i][j - 1];
            ybet[i][j] = 1 / beta;
        }
    }
}
