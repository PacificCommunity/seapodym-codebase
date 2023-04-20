#include "VarSimtunaFunc.h"

/// Forward main function called in simulation mode only for:
/// feeding habitat for young and adult life stages, with or
/// without seasonal switch between habitats depending on the
/// migration flag. See feeding_habitat.cpp

void VarSimtunaFunc::Feeding_Habitat_Index(
    VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Hf,
    int sp, int age, const int jday, const int t_count) {
    Hf.initialize();

    Hf_comp(param, mat, map, Hf, sp, age, jday, t_count);
}

// Generalized habitat index with seasonal migration between spawning and
// feeding habitats
void VarSimtunaFunc::Seasonal_Habitat_Index(
    VarParamCoupled& param, VarMatrices& mat, const PMap& map, dvar_matrix& Hs,
    dvar_matrix& Ha, int sp, int age, const int jday, const int t_count) {
    Ha_comp(param, mat, map, value(Hs), Ha, sp, jday);
}
