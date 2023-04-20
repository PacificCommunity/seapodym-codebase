#include <fvar.hpp>

#include "SeapodymCoupled.h"

/// Forward main and forward functions called in simulation mode only for:
/// ageing term discretisation either for all 0..n_a-1 or for A+ age class.

void SeapodymCoupled::Survival(
    dvar_matrix& N_a, dvar_matrix& N_a_1, const int a, const int sp) {
    int nb_cohorts = param->sp_nb_cohorts[sp];

    if (a < nb_cohorts - 1)

        Ageing(N_a, N_a_1);
    else

        AgePlus(N_a, N_a_1);
}

void SeapodymCoupled::Ageing(dvar_matrix& N_a, dvar_matrix& N_a_1) {
    for (int i = map.imin; i <= map.imax; i++) {
        const int jmin = map.jinf[i];
        const int jmax = map.jsup[i];
        for (int j = jmin; j <= jmax; j++) {
            if (map.carte(i, j)) {
                N_a.elem_value(i, j) = value(N_a_1(i, j));
            }
        }
    }
}

void SeapodymCoupled::AgePlus(dvar_matrix& N_a, dvar_matrix& N_a_1) {
    for (int i = map.imin; i <= map.imax; i++) {
        const int jmin = map.jinf[i];
        const int jmax = map.jsup[i];
        for (int j = jmin; j <= jmax; j++) {
            if (map.carte(i, j)) {
                N_a.elem_value(i, j) = value(N_a(i, j)) + value(N_a_1(i, j));
            }
        }
    }
}
