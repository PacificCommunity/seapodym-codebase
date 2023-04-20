#include "calpop.h"

/// Main function with memory control and adjoint functions for:
/// precalrec and calrec for adults functions. These routines finalize
/// the computation of diagonal elements by adding total mortalily rates
/// to coefficient 'b' (precalrec part), and then solve the discretized
/// ADE equation iteratively using ADI method with inner timestep deltaT/N.
/// Forward functions are in calrec_precalrec.cpp

void dv_calrec_precalrec(void);
void dv_calrec_with_catch_precalrec();
void dftridag1(
    dvector& a, dvector& bet, dvector& c, dvector& rhs, dvector& uvec,
    dvector& gam, dvector& dfa, dvector& dfbet, dvector& dfc, dvector& dfrhs,
    dvector& dfuvec, dvector& dfgam, int inf, int sup);
void dfxbet1(
    dmatrix& dfa, dmatrix& dfbm, dmatrix& dfc, dmatrix& dfxbet, dmatrix& xbet,
    const dmatrix a, const dmatrix bm, const dmatrix c,
    unsigned long int pos_map, const int maxn, const int dt);
int save_identifier_string2(char* str);
void verify_identifier_string2(char* str);
void save_long_int_value(unsigned long int x);
unsigned long int restore_long_int_value(void);

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
    // parameters for recomputing variables
    const double c_diff = value(param.dvarsC_diff_fish[sp]);
    const double mss_sp = value(param.dvarsMSS_species[sp]);
    const double sigma_sp = value(param.dvarsSigma_species[sp]);

    a = value(dvarsA);
    c = value(dvarsC);
    d = value(dvarsD);
    e = value(dvarsE);
    f = value(dvarsF);
    xbet = value(Xbet);
    ybet = value(Ybet);
    dmatrix M = value(mortality);

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
    save_identifier_string2((char*)"Precalrec_Calrec_begin");
    if (param.vert_movement(sp)) {
        // just for recomputation of diag.coef.
        mat.u.save_dmatrix_position();
        mat.v.save_dmatrix_position();
    }
    uu.save_dvar_matrix_position();
    dvarsA.save_dvar_matrix_position();
    dvarsB.save_dvar_matrix_position();
    dvarsBM.save_dvar_matrix_position();
    dvarsC.save_dvar_matrix_position();
    dvarsD.save_dvar_matrix_position();
    dvarsE.save_dvar_matrix_position();
    dvarsF.save_dvar_matrix_position();
    Xbet.save_dvar_matrix_position();
    Ybet.save_dvar_matrix_position();
    mortality.save_dvar_matrix_position();
    save_double_value(c_diff);
    save_double_value(mss_sp);
    save_double_value(sigma_sp);
    save_int_value(step_count);
    save_int_value(t_count);
    save_int_value(jday);
    save_int_value(no_mortality);
    save_int_value(age);
    save_int_value(sp);
    save_int_value(param.vert_movement(sp));
    unsigned long int cmat = (unsigned long int)&mat;
    save_long_int_value(cmat);
    unsigned long int pop = (unsigned long int)this;
    save_long_int_value(pop);
    unsigned long int pmap = (unsigned long int)&map;
    save_long_int_value(pmap);
    unsigned long int cparam = (unsigned long int)&param;
    save_long_int_value(cparam);
    unsigned long int crw = (unsigned long int)&rw;
    save_long_int_value(crw);
    save_int_value(year);
    save_int_value(month);
    if (call_calrec_catch) {
        mat.dvarCtot_age_obs(sp, age).save_dvar_matrix_position();
        mat.dvarCtot_age_est(sp, age).save_dvar_matrix_position();
    }
    save_identifier_string2((char*)"Precalrec_Calrec_end");

    if (!call_calrec_catch)
        gradient_structure::GRAD_STACK1->set_gradient_stack(
            dv_calrec_precalrec);
    else
        gradient_structure::GRAD_STACK1->set_gradient_stack(
            dv_calrec_with_catch_precalrec);
}

void dv_calrec_precalrec() {
    verify_identifier_string2((char*)"Precalrec_Calrec_end");
    unsigned month = restore_int_value();
    unsigned year = restore_int_value();
    unsigned long int pos_rw = restore_long_int_value();
    unsigned long int pos_param = restore_long_int_value();
    unsigned long int pos_map = restore_long_int_value();
    unsigned long int pos_pop = restore_long_int_value();
    unsigned long int pos_mat = restore_long_int_value();
    unsigned vert_move = restore_int_value();
    unsigned sp = restore_int_value();
    unsigned age = restore_int_value();
    unsigned no_mort = restore_int_value();
    unsigned jday = restore_int_value();
    unsigned t_count = restore_int_value();
    unsigned step_count = restore_int_value();
    const double sigma = restore_double_value();
    const double mss_sp = restore_double_value();
    const double c_diff = restore_double_value();
    const dvar_matrix_position M_pos = restore_dvar_matrix_position();
    const dvar_matrix_position ybet_pos = restore_dvar_matrix_position();
    const dvar_matrix_position xbet_pos = restore_dvar_matrix_position();
    const dvar_matrix_position f_pos = restore_dvar_matrix_position();
    const dvar_matrix_position e_pos = restore_dvar_matrix_position();
    const dvar_matrix_position d_pos = restore_dvar_matrix_position();
    const dvar_matrix_position c_pos = restore_dvar_matrix_position();
    const dvar_matrix_position bm_pos = restore_dvar_matrix_position();
    const dvar_matrix_position b_pos = restore_dvar_matrix_position();
    const dvar_matrix_position a_pos = restore_dvar_matrix_position();
    const dvar_matrix_position uu_pos = restore_dvar_matrix_position();
    dmatrix u, v;
    if (vert_move) {
        const dmatrix_position V_pos = restore_dmatrix_position();
        const dmatrix_position U_pos = restore_dmatrix_position();
    }
    verify_identifier_string2((char*)"Precalrec_Calrec_begin");

    dmatrix uu(uu_pos);
    dmatrix dfa = restore_dvar_matrix_derivatives(a_pos);
    dmatrix dfbm = restore_dvar_matrix_derivatives(bm_pos);
    dmatrix dfc = restore_dvar_matrix_derivatives(c_pos);
    dmatrix dfd = restore_dvar_matrix_derivatives(d_pos);
    dmatrix dfe = restore_dvar_matrix_derivatives(e_pos);
    dmatrix dff = restore_dvar_matrix_derivatives(f_pos);
    dmatrix dfxbet = restore_dvar_matrix_derivatives(xbet_pos);
    dmatrix dfybet = restore_dvar_matrix_derivatives(ybet_pos);
    dmatrix dfuu = restore_dvar_matrix_derivatives(uu_pos);
    dmatrix dfM = restore_dvar_matrix_derivatives(M_pos);
    dmatrix dfb = restore_dvar_matrix_derivatives(b_pos);

    CMatrices* mat = (CMatrices*)pos_mat;
    PMap* map = (PMap*)pos_map;
    CParam* param = (CParam*)pos_param;
    CCalpop* pop = (CCalpop*)pos_pop;
    CReadWrite* rw = (CReadWrite*)pos_rw;

    const int iterationNumber = pop->get_iterationN();
    const int maxn = pop->get_maxn();
    const int jinf = map->jmin;
    const int jsup = map->jmax;
    const int iinf = map->imin;
    const int isup = map->imax;
    // recompute tridiag coefficients
    dmatrix a, bm, c, d, e, f;
    a.allocate(jinf, jsup, map->iinf, map->isup);
    bm.allocate(jinf, jsup, map->iinf, map->isup);
    c.allocate(jinf, jsup, map->iinf, map->isup);
    d.allocate(iinf, isup, map->jinf, map->jsup);
    e.allocate(iinf, isup, map->jinf, map->jsup);
    f.allocate(iinf, isup, map->jinf, map->jsup);

    a.initialize();
    bm.initialize();
    c.initialize();
    d.initialize();
    e.initialize();
    f.initialize();

    dmatrix mortality(M_pos);
    mortality.initialize();
    int ind_adult_habitat = param->age_compute_habitat(sp, age);
    int age_adult_habitat = age;
    for (int aa = age; aa > 1; aa--) {
        if (param->age_compute_habitat(sp, aa - 1) !=
            param->age_compute_habitat(sp, aa)) {
            age_adult_habitat = aa;
            break;
        }
    }

    if (!no_mort) {
        // assuming all cohorts except the last one have the same age unit:
        double mean_age = age * param->sp_unit_cohort[sp][age - 1] +
                          0.5 * param->sp_unit_cohort[sp][age];
        pop->RecompM_sp(
            *map, *param, mortality,
            mat->adult_habitat(sp, t_count, ind_adult_habitat), mean_age, sp);
        pop->Recomp_total_mortality_comp(
            *map, *param, *mat, *rw, mortality, age, sp, year, month,
            step_count);
    }
    if (!vert_move)
        pop->RecompDiagCoef_adult(
            *map, *param, *mat, t_count, jday, mortality,
            mat->adult_habitat(sp, t_count, ind_adult_habitat), a, bm, c, d, e,
            f, sp, age_adult_habitat, mss_sp, c_diff, sigma);
    else
        pop->RecompDiagCoef_UV_adult(
            *map, *param, *mat, t_count, jday, mortality,
            mat->adult_habitat(sp, t_count, ind_adult_habitat), a, bm, c, d, e,
            f, sp, age_adult_habitat, mss_sp, c_diff, sigma);

    dvector uvec(0, maxn - 1);
    dvector rhs(0, maxn - 1);
    dvector gam(0, maxn - 1);
    uvec.initialize();
    rhs.initialize();
    gam.initialize();

    dvector dfuvec(0, maxn - 1);
    dvector dfgam(0, maxn - 1);
    dvector dfrhs(0, maxn - 1);
    dfuvec.initialize();
    dfgam.initialize();
    dfrhs.initialize();

    dmatrix xbet, ybet;
    xbet.allocate(jinf, jsup, map->iinf, map->isup);
    ybet.allocate(iinf, isup, map->jinf, map->jsup);

    xbet.initialize();
    ybet.initialize();

    // recompute xbet and ybet
    pop->xbet_comp(*map, xbet, a, bm, c, 2 * iterationNumber);
    pop->ybet_comp(*map, ybet, d, e, f, 2 * iterationNumber);
    const dmatrix_position uuint_pos(pop->uuint);
    dmatrix uuint(uuint_pos);

    // ADJOINT FOR CALREC FUNCTION
    for (int itr = iterationNumber; itr >= 1; itr--) {
        verify_identifier_string2((char*)"One_step_calrec_uu");
        uu = restore_dvar_matrix_value(uu_pos);

        dmatrix dfuuint(uuint_pos);
        dfuuint.initialize();

        uuint.initialize();

        // first recompute uuint from uu
        for (int j = map->jmin; j <= map->jmax; j++) {
            const int imin = map->iinf[j];
            const int imax = map->isup[j];

            for (int i = imin; i <= imax; i++)
                rhs[i] = -d[i][j] * uu(i, j - 1) +
                         (2 * iterationNumber - e[i][j]) * uu(i, j) -
                         f[i][j] * uu(i, j + 1);

            // tridag(a[j],xbet[j],c[j],rhs,uvec,gam,imin,imax);
            gam[imin] = rhs[imin];
            for (int i = imin + 1; i <= imax; i++)
                gam[i] = rhs[i] - gam[i - 1] * a[j][i] * xbet[j][i - 1];

            uvec[imax] = gam[imax] * xbet[j][imax];
            for (int i = imax - 1; i >= imin; i--)
                uvec[i] = (gam[i] - c[j][i] * uvec[i + 1]) * xbet[j][i];

            for (int i = imin; i <= imax; i++) uuint(i, j) = uvec[i];
        }
        // end of recomputation

        for (int i = isup; i >= iinf; i--) {
            const int jmin = map->jinf[i];
            const int jmax = map->jsup[i];

            // recomputing rhs(j)
            for (int j = jmin; j <= jmax; j++) {
                rhs[j] = -a[j][i] * uuint[i - 1][j] +
                         (2 * iterationNumber - bm[j][i]) * uuint[i][j] -
                         c[j][i] * uuint[i + 1][j];
            }

            for (int j = jmin; j <= jmax; j++) {
                // uu(i,j) = uvec(j);
                dfuvec(j) += dfuu(i, j);
                dfuu(i, j) = 0.0;
            }
            dftridag1(
                d(i), ybet(i), f(i), rhs, uvec, gam, dfd(i), dfybet(i), dff(i),
                dfrhs, dfuvec, dfgam, jmin, jmax);
            for (int j = jmax; j >= jmin; j--) {
                // rhs[j] =
                // -a[j][i]*uuint[i-1][j]+(2*iterationNumber-bm[j][i])*uuint[i][j]
                // - c[j][i]*uuint[i+1][j];
                dfa(j, i) -= uuint(i - 1, j) * dfrhs(j);
                dfuuint(i - 1, j) -= a(j, i) * dfrhs(j);
                dfbm(j, i) -= uuint(i, j) * dfrhs(j);
                dfuuint(i, j) += (2 * iterationNumber - bm[j][i]) * dfrhs(j);
                dfc(j, i) -= uuint(i + 1, j) * dfrhs(j);
                dfuuint(i + 1, j) -= c(j, i) * dfrhs(j);
                dfrhs(j) = 0.0;
            }
        }

        for (int j = jsup; j >= jinf; j--) {
            const int imin = map->iinf[j];
            const int imax = map->isup[j];

            // recomputing rhs(i)
            for (int i = imin; i <= imax; i++) {
                rhs[i] = -d[i][j] * uu[i][j - 1] +
                         (2 * iterationNumber - e[i][j]) * uu[i][j] -
                         f[i][j] * uu[i][j + 1];
            }

            for (int i = imax; i >= imin; i--) {
                // uuint[i][j] = uvec[i];
                dfuvec(i) += dfuuint(i, j);
                dfuuint(i, j) = 0.0;
            }
            dftridag1(
                a(j), xbet(j), c(j), rhs, uvec, gam, dfa(j), dfxbet(j), dfc(j),
                dfrhs, dfuvec, dfgam, imin, imax);

            for (int i = imax; i >= imin; i--) {
                // rhs[i] = -d[i][j]*uu[i][j-1] +
                // (2*iterationNumber-e[i][j])*uu[i][j] - f[i][j]*uu[i][j+1];
                dfd(i, j) -= uu(i, j - 1) * dfrhs(i);
                dfuu(i, j - 1) -= d(i, j) * dfrhs(i);
                dfe(i, j) -= uu(i, j) * dfrhs(i);
                dfuu(i, j) += (2 * iterationNumber - e[i][j]) * dfrhs(i);
                dff(i, j) -= uu(i, j + 1) * dfrhs(i);
                dfuu(i, j + 1) -= f(i, j) * dfrhs(i);
                dfrhs(i) = 0.0;
            }
        }
    }

    // ADJOINT FOR PRECALREC FUNCTION

    // xbet_comp(map,dt);
    dfxbet1(
        dfa, dfbm, dfc, dfxbet, xbet, a, bm, c, pos_map, maxn,
        2 * iterationNumber);

    for (int j = jsup; j >= jinf; j--) {
        const int imin = map->iinf(j);
        const int imax = map->isup(j);
        for (int i = imax; i >= imin; i--) {
            // bm[j][i] = b(j,i)+mort[i][j];
            dfM(i, j) += dfbm(j, i);
            dfb(j, i) += dfbm(j, i);
            dfbm(j, i) = 0.0;
        }
    }
    dfa.save_dmatrix_derivatives(a_pos);
    dfbm.save_dmatrix_derivatives(bm_pos);
    dfc.save_dmatrix_derivatives(c_pos);
    dfxbet.save_dmatrix_derivatives(xbet_pos);
    dfybet.save_dmatrix_derivatives(ybet_pos);
    dfd.save_dmatrix_derivatives(d_pos);
    dfe.save_dmatrix_derivatives(e_pos);
    dff.save_dmatrix_derivatives(f_pos);
    dfuu.save_dmatrix_derivatives(uu_pos);
    dfb.save_dmatrix_derivatives(b_pos);
    dfM.save_dmatrix_derivatives(M_pos);
}

void dv_calrec_with_catch_precalrec() {
    verify_identifier_string2((char*)"Precalrec_Calrec_end");
    const dvar_matrix_position Cest_pos = restore_dvar_matrix_position();
    const dvar_matrix_position Cobs_pos = restore_dvar_matrix_position();
    unsigned month = restore_int_value();
    unsigned year = restore_int_value();
    unsigned long int pos_rw = restore_long_int_value();
    unsigned long int pos_param = restore_long_int_value();
    unsigned long int pos_map = restore_long_int_value();
    unsigned long int pos_pop = restore_long_int_value();
    unsigned long int pos_mat = restore_long_int_value();
    unsigned vert_move = restore_int_value();
    unsigned sp = restore_int_value();
    unsigned age = restore_int_value();
    unsigned no_mort = restore_int_value();
    unsigned jday = restore_int_value();
    unsigned t_count = restore_int_value();
    unsigned step_count = restore_int_value();
    const double sigma = restore_double_value();
    const double mss_sp = restore_double_value();
    const double c_diff = restore_double_value();
    const dvar_matrix_position M_pos = restore_dvar_matrix_position();
    const dvar_matrix_position ybet_pos = restore_dvar_matrix_position();
    const dvar_matrix_position xbet_pos = restore_dvar_matrix_position();
    const dvar_matrix_position f_pos = restore_dvar_matrix_position();
    const dvar_matrix_position e_pos = restore_dvar_matrix_position();
    const dvar_matrix_position d_pos = restore_dvar_matrix_position();
    const dvar_matrix_position c_pos = restore_dvar_matrix_position();
    const dvar_matrix_position bm_pos = restore_dvar_matrix_position();
    const dvar_matrix_position b_pos = restore_dvar_matrix_position();
    const dvar_matrix_position a_pos = restore_dvar_matrix_position();
    const dvar_matrix_position uu_pos = restore_dvar_matrix_position();
    dmatrix u, v;
    if (vert_move) {
        const dmatrix_position V_pos = restore_dmatrix_position();
        const dmatrix_position U_pos = restore_dmatrix_position();
    }
    verify_identifier_string2((char*)"Precalrec_Calrec_begin");

    dmatrix uu(uu_pos);
    dmatrix dfa = restore_dvar_matrix_derivatives(a_pos);
    dmatrix dfbm = restore_dvar_matrix_derivatives(bm_pos);
    dmatrix dfc = restore_dvar_matrix_derivatives(c_pos);
    dmatrix dfd = restore_dvar_matrix_derivatives(d_pos);
    dmatrix dfe = restore_dvar_matrix_derivatives(e_pos);
    dmatrix dff = restore_dvar_matrix_derivatives(f_pos);
    dmatrix dfxbet = restore_dvar_matrix_derivatives(xbet_pos);
    dmatrix dfybet = restore_dvar_matrix_derivatives(ybet_pos);
    dmatrix dfuu = restore_dvar_matrix_derivatives(uu_pos);
    dmatrix dfM = restore_dvar_matrix_derivatives(M_pos);
    dmatrix dfb = restore_dvar_matrix_derivatives(b_pos);
    dmatrix dfCobs = restore_dvar_matrix_derivatives(Cobs_pos);
    dmatrix dfCest = restore_dvar_matrix_derivatives(Cest_pos);

    CMatrices* mat = (CMatrices*)pos_mat;
    PMap* map = (PMap*)pos_map;
    CParam* param = (CParam*)pos_param;
    CCalpop* pop = (CCalpop*)pos_pop;
    CReadWrite* rw = (CReadWrite*)pos_rw;

    const int iterationNumber = pop->get_iterationN();
    const int maxn = pop->get_maxn();
    const int jinf = map->jmin;
    const int jsup = map->jmax;
    const int iinf = map->imin;
    const int isup = map->imax;
    const int nti = param->get_nbi();
    const int ntj = param->get_nbj();
    dmatrix uuint_t;  // before catch removal
    uuint_t.allocate(0, nti, 0, ntj);
    uuint_t.initialize();

    // recompute tridiag coefficients
    dmatrix a, bm, c, d, e, f;
    a.allocate(jinf, jsup, map->iinf, map->isup);
    bm.allocate(jinf, jsup, map->iinf, map->isup);
    c.allocate(jinf, jsup, map->iinf, map->isup);
    d.allocate(iinf, isup, map->jinf, map->jsup);
    e.allocate(iinf, isup, map->jinf, map->jsup);
    f.allocate(iinf, isup, map->jinf, map->jsup);

    a.initialize();
    bm.initialize();
    c.initialize();
    d.initialize();
    e.initialize();
    f.initialize();

    dmatrix C;
    C.allocate(iinf, isup, map->jinf, map->jsup);
    C.initialize();
    pop->Recomp_total_obs_catch_age(
        *map, *param, *mat, *rw, C, age, sp, year, month, t_count);

    dmatrix mortality(M_pos);
    mortality.initialize();

    int ind_adult_habitat = param->age_compute_habitat(sp, age);
    int age_adult_habitat = age;
    for (int aa = age; aa > 1; aa--) {
        if (param->age_compute_habitat(sp, aa - 1) !=
            param->age_compute_habitat(sp, aa)) {
            age_adult_habitat = aa;
            break;
        }
    }
    if (!no_mort) {
        // assuming all cohorts except the last one have the same age unit:
        double mean_age = age * param->sp_unit_cohort[sp][age - 1] +
                          0.5 * param->sp_unit_cohort[sp][age];
        pop->RecompM_sp(
            *map, *param, mortality,
            mat->adult_habitat(sp, t_count, ind_adult_habitat), mean_age, sp);
        pop->Recomp_total_mortality_comp(
            *map, *param, *mat, *rw, mortality, age, sp, year, month,
            step_count);
    }
    // cout << age << " "<<
    // norm(mat->adult_habitat(sp,t_count,ind_adult_habitat)) << " " <<
    // norm(mortality) << endl;

    if (!vert_move)
        pop->RecompDiagCoef_adult(
            *map, *param, *mat, t_count, jday, mortality,
            mat->adult_habitat(sp, t_count, ind_adult_habitat), a, bm, c, d, e,
            f, sp, age_adult_habitat, mss_sp, c_diff, sigma);
    else
        pop->RecompDiagCoef_UV_adult(
            *map, *param, *mat, t_count, jday, mortality,
            mat->adult_habitat(sp, t_count, ind_adult_habitat), a, bm, c, d, e,
            f, sp, age_adult_habitat, mss_sp, c_diff, sigma);

    dvector uvec(0, maxn - 1);
    dvector rhs(0, maxn - 1);
    dvector gam(0, maxn - 1);
    uvec.initialize();
    rhs.initialize();
    gam.initialize();

    dvector dfuvec(0, maxn - 1);
    dvector dfgam(0, maxn - 1);
    dvector dfrhs(0, maxn - 1);
    dfuvec.initialize();
    dfgam.initialize();
    dfrhs.initialize();

    dmatrix xbet, ybet;
    xbet.allocate(jinf, jsup, map->iinf, map->isup);
    ybet.allocate(iinf, isup, map->jinf, map->jsup);

    xbet.initialize();
    ybet.initialize();

    // recompute xbet and ybet
    pop->xbet_comp(*map, xbet, a, bm, c, 2 * iterationNumber);
    pop->ybet_comp(*map, ybet, d, e, f, 2 * iterationNumber);
    // if (age>25) cout << age << " " << age_adult_habitat << " " << norm(a) <<
    // " " << norm(bm) << " " << norm(c) << " " << norm(d) << " " << norm(e) <<
    // " " << norm(f) << " " << norm(xbet) << " " << norm(ybet) << endl;
    const dmatrix_position uuint_pos(pop->uuint);
    dmatrix uuint(uuint_pos);

    // ADJOINT FOR CALREC FUNCTION
    for (int itr = iterationNumber; itr >= 1; itr--) {
        verify_identifier_string2((char*)"One_step_calrec_uu");
        uu = restore_dvar_matrix_value(uu_pos);
        //		uu = uun[itr];
        // cout << itr << " " << norm(uu) << endl;

        dmatrix dfuuint(uuint_pos);
        dfuuint.initialize();

        // uuint_t.initialize();

        // first recompute uuint from uu
        for (int j = map->jmin; j <= map->jmax; j++) {
            const int imin = map->iinf[j];
            const int imax = map->isup[j];

            for (int i = imin; i <= imax; i++)
                rhs[i] = -d[i][j] * uu(i, j - 1) +
                         (2 * iterationNumber - e[i][j]) * uu(i, j) -
                         f[i][j] * uu(i, j + 1);

            // tridag(a[j],xbet[j],c[j],rhs,uvec,gam,imin,imax);
            gam[imin] = rhs[imin];
            for (int i = imin + 1; i <= imax; i++)
                gam[i] = rhs[i] - gam[i - 1] * a[j][i] * xbet[j][i - 1];

            uvec[imax] = gam[imax] * xbet[j][imax];
            for (int i = imax - 1; i >= imin; i--)
                uvec[i] = (gam[i] - c[j][i] * uvec[i + 1]) * xbet[j][i];

            for (int i = imin; i <= imax; i++) {
                if (C(i, j) == 0)
                    uuint(i, j) = uvec[i];
                else
                    uuint(i, j) =
                        uvec[i] -
                        uvec[i] *
                            param->func_limit_one(C(i, j) / (uvec[i] + 1e-14)) /
                            iterationNumber;

                // uuint(i,j) = uvec[i];
                uuint_t(i, j) = uvec[i];
            }

            // for (int i = imin; i <= imax; i++){
            //	if (C(i,j)>0)
            //		uuint(i,j) = uuint(i,j) - uuint(i,j) *
            //param->func_limit_one(C(i,j)/uuint(i,j))/iterationNumber;
            // }
        }
        // end of recomputation

        for (int i = isup; i >= iinf; i--) {
            const int jmin = map->jinf[i];
            const int jmax = map->jsup[i];

            // recomputing rhs(j)
            for (int j = jmin; j <= jmax; j++) {
                rhs[j] = -a[j][i] * uuint[i - 1][j] +
                         (2 * iterationNumber - bm[j][i]) * uuint[i][j] -
                         c[j][i] * uuint[i + 1][j];
            }

            for (int j = jmin; j <= jmax; j++) {
                // uu(i,j) = uvec(j);
                dfuvec(j) += dfuu(i, j);
                dfuu(i, j) = 0.0;
            }
            dftridag1(
                d(i), ybet(i), f(i), rhs, uvec, gam, dfd(i), dfybet(i), dff(i),
                dfrhs, dfuvec, dfgam, jmin, jmax);
            for (int j = jmax; j >= jmin; j--) {
                // rhs[j] =
                // -a[j][i]*uuint[i-1][j]+(2*iterationNumber-bm[j][i])*uuint[i][j]
                // - c[j][i]*uuint[i+1][j];
                dfa(j, i) -= uuint(i - 1, j) * dfrhs(j);
                dfuuint(i - 1, j) -= a(j, i) * dfrhs(j);
                dfbm(j, i) -= uuint(i, j) * dfrhs(j);
                dfuuint(i, j) += (2 * iterationNumber - bm[j][i]) * dfrhs(j);
                dfc(j, i) -= uuint(i + 1, j) * dfrhs(j);
                dfuuint(i + 1, j) -= c(j, i) * dfrhs(j);
                dfrhs(j) = 0.0;
            }
        }

        for (int j = jsup; j >= jinf; j--) {
            const int imin = map->iinf[j];
            const int imax = map->isup[j];

            dvector dfuuint_t;
            dfuuint_t.allocate(0, nti);
            dfuuint_t.initialize();

            // recomputing rhs(i)
            for (int i = imin; i <= imax; i++) {
                rhs[i] = -d[i][j] * uu[i][j - 1] +
                         (2 * iterationNumber - e[i][j]) * uu[i][j] -
                         f[i][j] * uu[i][j + 1];
            }

            for (int i = imax; i >= imin; i--) {
                if (C(i, j) == 0) {
                    // uuint[i][j] = uvec[i];
                    dfuvec(i) += dfuuint(i, j);
                    dfuuint(i, j) = 0.0;
                } else {
                    // recompute
                    double arg = C(i, j) / (uuint_t(i, j) + 1e-14);
                    double func = param->func_limit_one(arg);

                    // C_est(i,j) += uvec(i) * func / iterationNumber;
                    dfuvec(i) += (func / iterationNumber) * dfCest(i, j);
                    double dffunc =
                        (uuint_t(i, j) / iterationNumber) * dfCest(i, j);

                    // uuint(i,j) = uvec(i) - uvec(i) * func / iterationNumber;
                    dfuvec(i) += (1.0 - func / iterationNumber) * dfuuint(i, j);
                    dffunc -= (uuint_t(i, j) / iterationNumber) * dfuuint(i, j);
                    dfuuint(i, j) = 0.0;

                    // double func =
                    // param->func_limit_one(C(i,j)/(uvec(i)+1e-14);
                    double dfarg = param->dffunc_limit_one(arg, dffunc);
                    dfuvec(i) -= (arg / (uuint_t(i, j) + 1e-14)) * dfarg;
                    dfCobs(i, j) += dfarg / (uuint_t(i, j) + 1e-14);

                    /*					//recompute
                                                            double arg =
                       C(i,j)/(uuint_t(i,j)+1e-14); double func   =
                       param->func_limit_one(arg);

                                                            //C_est.elem_value(i,j)
                       += C_est_ij_itr; dfCest_pr(i,j) += dfCest(i,j); double
                       dfCest_itr = dfCest(i,j); dfCest(i,j)     = 0.0;

                                                            //uuint(i,j) =
                       uvec[i] - C_est_ij_itr; dfuvec(i)   += dfuuint(i,j);
                                                            dfCest_itr  -=
                       dfuuint(i,j); dfuuint(i,j) = 0.0;

                                                            //double
                       C_est_ij_itr = uvec(i) * func / iterationNumber;
                                                            dfuvec(i)   += (func
                       / iterationNumber) * dfCest_itr; double
                       dffunc=(uuint_t(i,j)/iterationNumber) * dfCest_itr;

                                                            //double func =
                       param->func_limit_one(C(i,j)/(uvec(i)+1e-14); dfuvec(i)
                       -= (arg/(uuint_t(i,j)+1e-14)) *
                       param->dffunc_limit_one(arg,dffunc);
                    */
                }
            }
            dftridag1(
                a(j), xbet(j), c(j), rhs, uvec, gam, dfa(j), dfxbet(j), dfc(j),
                dfrhs, dfuvec, dfgam, imin, imax);

            for (int i = imax; i >= imin; i--) {
                // rhs[i] = -d[i][j]*uu[i][j-1] +
                // (2*iterationNumber-e[i][j])*uu[i][j] - f[i][j]*uu[i][j+1];
                dfd(i, j) -= uu(i, j - 1) * dfrhs(i);
                dfuu(i, j - 1) -= d(i, j) * dfrhs(i);
                dfe(i, j) -= uu(i, j) * dfrhs(i);
                dfuu(i, j) += (2 * iterationNumber - e[i][j]) * dfrhs(i);
                dff(i, j) -= uu(i, j + 1) * dfrhs(i);
                dfuu(i, j + 1) -= f(i, j) * dfrhs(i);
                dfrhs(i) = 0.0;
            }
        }
    }

    // ADJOINT FOR PRECALREC FUNCTION

    // xbet_comp(map,dt);
    dfxbet1(
        dfa, dfbm, dfc, dfxbet, xbet, a, bm, c, pos_map, maxn,
        2 * iterationNumber);

    for (int j = jsup; j >= jinf; j--) {
        const int imin = map->iinf(j);
        const int imax = map->isup(j);
        for (int i = imax; i >= imin; i--) {
            // bm[j][i] = b(j,i)+mort[i][j];
            dfM(i, j) += dfbm(j, i);
            dfb(j, i) += dfbm(j, i);
            dfbm(j, i) = 0.0;
        }
    }
    // cout << age << " " << norm(dfCobs) << " " << norm(dfCest) << " " <<
    // norm(dfCest_pr) << " " << norm(dfuu) << endl;
    dfa.save_dmatrix_derivatives(a_pos);
    dfbm.save_dmatrix_derivatives(bm_pos);
    dfc.save_dmatrix_derivatives(c_pos);
    dfxbet.save_dmatrix_derivatives(xbet_pos);
    dfybet.save_dmatrix_derivatives(ybet_pos);
    dfd.save_dmatrix_derivatives(d_pos);
    dfe.save_dmatrix_derivatives(e_pos);
    dff.save_dmatrix_derivatives(f_pos);
    dfuu.save_dmatrix_derivatives(uu_pos);
    dfb.save_dmatrix_derivatives(b_pos);
    dfM.save_dmatrix_derivatives(M_pos);
    dfCobs.save_dmatrix_derivatives(Cobs_pos);
    dfCest.save_dmatrix_derivatives(Cest_pos);
}

void dftridag1(
    dvector& a, dvector& bet, dvector& c, dvector& rhs, dvector& uvec,
    dvector& gam, dvector& dfa, dvector& dfbet, dvector& dfc, dvector& dfrhs,
    dvector& dfuvec, dvector& dfgam, int inf, int sup) {
    int j;

    // recompute gam and uvec values
    gam[inf] = rhs[inf];
    for (int j = inf + 1; j <= sup; j++)
        gam[j] = rhs[j] - gam[j - 1] * a[j] * bet[j - 1];

    uvec[sup] = gam[sup] * bet[sup];
    for (int j = sup - 1; j >= inf; j--)
        uvec[j] = (gam[j] - c[j] * uvec[j + 1]) * bet[j];
    ////////////////////////////////

    for (j = inf; j < sup; j++) {
        // uvec[j] = (gam[j]-c[j]*uvec[j+1])*bet[j];
        dfgam[j] += bet[j] * dfuvec[j];
        dfc[j] -= uvec[j + 1] * bet[j] * dfuvec[j];
        dfuvec[j + 1] -= c[j] * bet[j] * dfuvec[j];
        dfbet[j] += (gam[j] - c[j] * uvec[j + 1]) * dfuvec[j];
        dfuvec[j] = 0.0;
    }

    // uvec[sup] = gam[sup]*bet[sup];
    dfgam[sup] += bet[sup] * dfuvec[sup];
    dfbet[sup] += gam[sup] * dfuvec[sup];
    dfuvec[sup] = 0.0;

    for (j = sup; j >= inf + 1; j--) {
        // gam[j]   = rhs[j]-gam[j-1]*a[j]*bet[j-1];
        dfrhs[j] += dfgam[j];
        dfgam[j - 1] -= a[j] * bet[j - 1] * dfgam[j];
        dfa[j] -= gam[j - 1] * bet[j - 1] * dfgam[j];
        dfbet[j - 1] -= gam[j - 1] * a[j] * dfgam[j];
        dfgam[j] = 0.0;
    }

    // gam[inf] = rhs[inf];
    dfrhs[inf] += dfgam[inf];
    dfgam[inf] = 0.0;

}  // end of dftridag

void dfxbet1(
    dmatrix& dfa, dmatrix& dfbm, dmatrix& dfc, dmatrix& dfxbet, dmatrix& xbet,
    const dmatrix a, const dmatrix bm, const dmatrix c,
    unsigned long int pos_map, const int maxn, const int dt) {
    PMap* map = (PMap*)pos_map;

    const int jmin = map->jmin;
    const int jmax = map->jmax;

    dmatrix w(jmin, jmax, 0, maxn - 1);
    w.initialize();
    xbet.initialize();

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
            double dfw = -(1 / (w(j, i) * w(j, i))) * dfxbet(j, i);
            dfxbet(j, i) = 0.0;

            // w[j][i] = bm[j][i]+dt-c[j][i-1]*a[j][i]*xbet[j][i-1];
            dfbm(j, i) += dfw;
            dfc(j, i - 1) -= a(j, i) * xbet(j, i - 1) * dfw;
            dfa(j, i) -= c(j, i - 1) * xbet(j, i - 1) * dfw;
            dfxbet(j, i - 1) -= c(j, i - 1) * a(j, i) * dfw;
            // dfw(j,i)      = 0.0;
        }
        // xbet[j][inf] = 1/(bm[j][imin]+dt);
        dfbm(j, imin) -=
            (1 / ((bm(j, imin) + dt) * (bm(j, imin) + dt))) * dfxbet(j, imin);
        dfxbet(j, imin) = 0.0;
    }
}
