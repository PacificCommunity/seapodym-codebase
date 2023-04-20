#include "calpop.h"

/// Main function with memory control and adjoint functions for:
/// tridag_bet function. This routine precomputes an operator in
/// the Gaussian solver of the tridiagonal linear system, which
/// does not change during iterations.
/// Forward functions are in tridag_bet.cpp

void dv_xbet_comp1();

void CCalpop::Xbet_comp1(const PMap& map, int dt) {
    dmatrix xbet_c = value(Xbet);
    dmatrix bm_c = value(dvarsBM);

    xbet_comp(map, xbet_c, a, bm_c, c, dt);

    Xbet = nograd_assign(xbet_c);

    dvar_matrix W(0, maxn - 1, 0, maxn - 1);
    W.initialize();

    save_identifier_string((char*)"xbet_comp1_begin");
    a.save_dmatrix_value();
    a.save_dmatrix_position();
    dvarsBM.save_dvar_matrix_value();
    dvarsBM.save_dvar_matrix_position();
    c.save_dmatrix_value();
    c.save_dmatrix_position();
    Xbet.save_dvar_matrix_position();
    W.save_dvar_matrix_position();
    unsigned long int pmap = (unsigned long int)&map;
    save_int_value(pmap);
    save_int_value(maxn);
    save_int_value(dt);
    save_identifier_string((char*)"xbet_comp1_end");

    gradient_structure::GRAD_STACK1->set_gradient_stack(dv_xbet_comp1);
}

void dv_xbet_comp1() {
    verify_identifier_string((char*)"xbet_comp1_end");
    const int dt = restore_int_value();
    unsigned maxn = restore_int_value();
    unsigned long int pos_map = restore_int_value();
    const dvar_matrix_position w_pos = restore_dvar_matrix_position();
    const dvar_matrix_position xbet_pos = restore_dvar_matrix_position();
    const dmatrix_position c_pos = restore_dmatrix_position();
    dmatrix c = restore_dmatrix_value(c_pos);
    const dvar_matrix_position bm_pos = restore_dvar_matrix_position();
    dmatrix bm = restore_dvar_matrix_value(bm_pos);
    const dmatrix_position a_pos = restore_dmatrix_position();
    dmatrix a = restore_dmatrix_value(a_pos);
    verify_identifier_string((char*)"xbet_comp1_begin");

    dmatrix dfbm = restore_dvar_matrix_derivatives(bm_pos);
    dmatrix dfxbet = restore_dvar_matrix_derivatives(xbet_pos);
    dmatrix dfw = restore_dvar_matrix_derivatives(w_pos);

    PMap* map = (PMap*)pos_map;

    const int jmin = map->jmin;
    const int jmax = map->jmax;

    dmatrix xbet(jmin, jmax, 0, maxn - 1);
    dmatrix w(jmin, jmax, 0, maxn - 1);

    // recompute w and xbet
    for (int j = jmin; j <= jmax; j++) {
        const int imin = map->iinf[j];
        const int imax = map->isup[j];
        xbet[j][imin] = 1 / (bm[j][imin] + dt);
        for (int i = imin + 1; i <= imax; i++) {
            w[j][i] = bm[j][i] + dt - c[j][i - 1] * a[j][i] * xbet[j][i - 1];
            xbet[j][i] = 1 / w[j][i];
        }
    }

    for (int j = jmax; j >= jmin; j--) {
        const int imin = map->iinf[j];
        const int imax = map->isup[j];
        for (int i = imax; i >= imin + 1; i--) {
            // xbet[j][i] = 1/w[j][i];
            dfw(j, i) -= (1 / (w(j, i) * w(j, i))) * dfxbet(j, i);
            dfxbet(j, i) = 0.0;

            // w[j][i] = bm[j][i]+dt-c[j][i-1]*a[j][i]*xbet[j][i-1];
            dfbm(j, i) += dfw(j, i);
            dfxbet(j, i - 1) -= c(j, i - 1) * a(j, i) * dfw(j, i);
            dfw(j, i) = 0.0;
        }
        // xbet[j][inf] = 1/(bm[j][imin]+dt);
        dfbm(j, imin) -=
            (1 / ((bm(j, imin) + dt) * (bm(j, imin) + dt))) * dfxbet(j, imin);
        dfxbet(j, imin) = 0.0;
    }

    dfbm.save_dmatrix_derivatives(bm_pos);
    dfxbet.save_dmatrix_derivatives(xbet_pos);
    dfw.save_dmatrix_derivatives(w_pos);
}
