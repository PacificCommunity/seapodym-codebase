#ifndef __VarMatrices_h__
#define __VarMatrices_h__

#include "Matrices.h"

/*!
\brief Seapodym DVAR matrices class.
*/
class VarMatrices : public CMatrices {
   public:
    VarMatrices() { /*DoesNothing*/
    }
    //	VarMatrices(const VarMatrices&) {/*DoesNothing*/}
    virtual ~VarMatrices() { /*DoesNothing*/
    }

   public:
    void CreateMatHabitat(
        PMap& map, const int nb_species, const int nforage, const int nblayer,
        const int nb_ages, int t0, int nbt, const int nbi, const int nbj,
        const ivector sp_adult_age0, const ivector sp_nb_age_class,
        const imatrix age_compute_habitat) {
        CMatrices::createMatHabitat(
            map, nforage, nb_species, t0, nbt, sp_adult_age0, sp_nb_age_class,
            age_compute_habitat);

        dvarSeasonSwitch.allocate(0, nb_species);
        dvarSigmaSeason.allocate(0, nb_species);
        for (int sp = 0; sp < nb_species; sp++) {
            dvarSeasonSwitch(sp).allocate(
                map.imin, map.imax, map.jmin, map.jmax);
            dvarSigmaSeason(sp).allocate(
                map.imin, map.imax, map.jmin, map.jmax);
            // dvarSeasonSwitch(sp).allocate(map.jmin, map.jmax);
            dvarSeasonSwitch(sp).initialize();
            dvarSigmaSeason(sp).initialize();
        }

        dvarF_access.allocate(0, nforage - 1, 0, nb_ages - 1);
        dvarZ_access.allocate(0, nblayer - 1, 0, nb_ages - 1);
        for (int a = 0; a < nb_ages; a++) {
            for (int n = 0; n < nforage; n++) {
                dvarF_access(n, a).allocate(
                    map.imin, map.imax, map.jinf, map.jsup);
                dvarF_access(n, a).initialize();
            }
            for (int l = 0; l < nblayer; l++) {
                dvarZ_access(l, a).allocate(
                    map.imin, map.imax, map.jinf, map.jsup);
                dvarZ_access(l, a).initialize();
            }
        }
    }

    void CreateMatTransport(PMap& map, const int nbi, const int nbj) {
        CMatrices::createMatTransport(map);  //,(int)nbi, (int)nbj);

        dvarsU.allocate(map.imin, map.imax, map.jinf, map.jsup);
        dvarsV.allocate(map.imin, map.imax, map.jinf, map.jsup);
        dvarsU.initialize();
        dvarsV.initialize();

        dvarsAdvection_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
        dvarsAdvection_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
        dvarsDiffusion_x.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);
        dvarsDiffusion_y.allocate(map.imin1, map.imax1, map.jinf1, map.jsup1);

        dvarsDiffusion_x.initialize();
        dvarsDiffusion_y.initialize();
        dvarsAdvection_x.initialize();
        dvarsAdvection_y.initialize();
    }

    /// void CreateMatSpecies(PMap& map, int nbi, int nbj, int nb_species, const
    /// ivector& sp_nb_age_class_jv, const ivector& sp_nb_age_class) {
    void CreateMatSpecies(
        PMap& map, int t0, int nbt, int nbi, int nbj, int nb_species,
        const ivector a0_adult, const ivector& sp_nb_cohorts) {
        // CMatrices::createMatSpecies(map, t0, nbt, nbi, nbj, nb_species,
        // a0_adult, sp_nb_cohorts);
        CMatrices::createMatSpecies(
            map, t0, nbt, nbi, nbj, 1, a0_adult, sp_nb_cohorts);
        dvarDensity.allocate(0, nb_species - 1);
        dvarDensity.allocate(0, nb_species - 1);
        for (int sp = 0; sp < nb_species; sp++) {
            // const int agemax_sp = sp_nb_cohorts(sp);
            const int agemax_sp = sp_nb_cohorts(0);
            dvarDensity(sp).allocate(0, agemax_sp - 1);
            dvarDensity(sp).initialize();

            for (int age = 0; age < agemax_sp; age++) {
                dvarDensity(sp, age).allocate(
                    map.imin1, map.imax1, map.jinf1, map.jsup1);
                dvarDensity(sp, age).initialize();
            }
        }
        /*
                        const int age_max_1 = 48-1;
                        dvarDensity_age.allocate(0,age_max-1);
                        for (int a = 0; a<age_max; a++){
                                dvarDensity_age(a).allocate(map.imin1,
           map.imax1, map.jinf1, map.jsup1); dvarDensity_age(a).initialize();
                        }
        */
    }

    void CreateMatCatch(
        PMap& map, int nbi, int nbj, int nb_species, const IVECTOR& nb_fleet,
        const ivector a0_adult, const IVECTOR& nb_cohorts,
        const IVECTOR& nb_region) {
        CMatrices::createMatCatch(
            map, nbi, nbj, nb_species, nb_fleet, a0_adult, nb_cohorts,
            nb_region);

        dvarCtot_age_obs.allocate(0, nb_species - 1);
        dvarCtot_age_est.allocate(0, nb_species - 1);
        dvarCtot_age_obs.initialize();
        dvarCtot_age_est.initialize();

        for (int sp = 0; sp < nb_species; sp++) {
            const int a0 = a0_adult(sp);
            const int agemax = nb_cohorts(sp);

            dvarCtot_age_obs[sp].allocate(a0, agemax - 1);
            dvarCtot_age_est[sp].allocate(a0, agemax - 1);
            dvarCtot_age_obs[sp].initialize();
            dvarCtot_age_est[sp].initialize();

            for (int age = a0; age < agemax; age++) {
                dvarCtot_age_obs(sp, age).allocate(
                    map.imin, map.imax, map.jinf, map.jsup);
                dvarCtot_age_est(sp, age).allocate(
                    map.imin, map.imax, map.jinf, map.jsup);
                dvarCtot_age_obs(sp, age).initialize();
                dvarCtot_age_est(sp, age).initialize();
            }
        }

        dvarCatch_est.allocate(0, nb_species - 1);
        dvarCatch_est.initialize();

        for (int sp = 0; sp < nb_species; sp++) {
            const int fleetmax = nb_fleet(sp);
            dvarCatch_est(sp).allocate(0, fleetmax - 1);
            dvarCatch_est(sp).initialize();

            for (int fleet = 0; fleet < fleetmax; fleet++) {
                dvarCatch_est(sp, fleet).allocate(
                    map.imin, map.imax, map.jinf, map.jsup);
                dvarCatch_est(sp, fleet).initialize();
            }
        }

        dvarLF_est.allocate(0, nb_species - 1);
        dvarLF_est.initialize();

        for (int sp = 0; sp < nb_species; sp++) {
            const int a0 = a0_adult(sp);
            const int agemax = nb_cohorts(sp);

            dvarLF_est[sp].allocate(a0, agemax - 1);
            dvarLF_est[sp].initialize();

            for (int age = a0; age < agemax; age++) {
                const int fleetmax = nb_fleet(sp);
                dvarLF_est[sp][age].allocate(0, fleetmax - 1);
                dvarLF_est[sp][age].initialize();

                for (int fleet = 0; fleet < fleetmax; fleet++) {
                    const int regionmax = nb_region(sp);
                    dvarLF_est[sp][age][fleet].allocate(0, regionmax - 1);
                    dvarLF_est[sp][age][fleet].initialize();
                }
            }
        }
    }
    /*	void allocate_dvmatr(const int imin, const int imax, const ivector jinf,
    const ivector jsup){ dvmatr1.allocate(imin,imax,jinf,jsup);
    dvmatr1.initialize();
    //		dvmatr2.allocate(imin,imax,jinf,jsup); dvmatr2.initialize();
    //		dvmatr3.allocate(imin,imax,jinf,jsup); dvmatr3.initialize();
    //		dvmatr4.allocate(imin,imax,jinf,jsup); dvmatr4.initialize();
            }
    */

   public:
    dvar_matrix dvarsU;
    dvar_matrix dvarsV;

    //	DVAR3_ARRAY dvarDensity_age;

    DVAR4_ARRAY dvarF_access;
    DVAR4_ARRAY dvarZ_access;
    DVAR4_ARRAY dvarDensity;
    DVAR4_ARRAY dvarCatch_est;
    DVAR4_ARRAY dvarLF_est;
    dvar4_array dvarCtot_age_obs;
    dvar4_array dvarCtot_age_est;

    dvar3_array dvarSeasonSwitch;
    dvar3_array dvarSigmaSeason;
    dvar_matrix dvarsDiffusion_x;
    dvar_matrix dvarsDiffusion_y;
    dvar_matrix dvarsAdvection_x;
    dvar_matrix dvarsAdvection_y;
    // dump matrices for computing derivatives:
    //	dvar_matrix dvmatr1, dvmatr2, dvmatr3, dvmatr4;
};
#endif
