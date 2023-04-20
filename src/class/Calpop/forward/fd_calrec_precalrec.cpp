#include "calpop.h"

/// Forward main function called in simulation mode only for:
/// precalrec and calrec for adults functions. See calrec_precalrec.cpp

void CCalpop::Precalrec_Calrec_adult(
    const PMap& map, VarMatrices& mat, VarParamCoupled& param, CReadWrite& rw,
    dvar_matrix& uu, dvar_matrix& mortality, const int t_count,
    const bool fishing, const int age, const int sp, const int year,
    const int month, const int jday, const int step_count,
    const int no_mortality) {
    bool call_calrec_catch = false;
    if (param.fisheries_no_effort_exist[sp] && fishing) {
        call_calrec_catch = true;
        mat.dvarCtot_age_est[sp][age].initialize();
    }
    const int nb_fishery = param.get_nbfishery();
    if (fishing && (param.fisheries_no_effort_exist[sp] < nb_fishery)) {
        Precalrec_total_mortality_comp(
            map, param, mat, rw, mortality, age, sp, t_count, year, month,
            step_count);
    }

    a = value(dvarsA);
    c = value(dvarsC);
    d = value(dvarsD);
    e = value(dvarsE);
    f = value(dvarsF);
    xbet = value(Xbet);
    ybet = value(Ybet);
    dmatrix M = value(mortality);

    dvar_matrix W(0, maxn - 1, 0, maxn - 1);
    W.initialize();

    bm = value(dvarsB);
    precalrec_juv_comp(map, bm, M);

    dvarsBM = nograd_assign(bm);

    xbet_comp(map, xbet, a, bm, c, 2 * iterationNumber);

    Xbet = nograd_assign(xbet);

    if (!call_calrec_catch)
        calrec1(map, uu, M);
    else {
        calrec_with_catch(
            map, param, uu, value(mat.dvarCtot_age_obs[sp][age]),
            mat.dvarCtot_age_est[sp][age]);
    }
}
