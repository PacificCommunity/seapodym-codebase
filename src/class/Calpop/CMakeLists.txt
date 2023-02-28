set(CALPOP_DV
    adjoint/dv_caldia.cpp
    adjoint/dv_calrec_adre.cpp
    adjoint/dv_calrec_precalrec.cpp
    adjoint/dv_precalrec_juv.cpp
    adjoint/dv_predicted_catch_without_effort.cpp
    adjoint/dv_predicted_catch.cpp
    adjoint/dv_total_exploited_biomass.cpp
    adjoint/dv_total_mortality_comp.cpp
    adjoint/dv_total_obs_catch_age.cpp
    adjoint/dv_tridag_bet.cpp
)

set(CALPOP_FD
    forward/fd_caldia.cpp
    forward/fd_calrec_adre.cpp
    forward/fd_calrec_precalrec.cpp
    forward/fd_precalrec_juv.cpp
    forward/fd_predicted_catch_without_effort.cpp
    forward/fd_predicted_catch.cpp
    forward/fd_total_exploited_biomass.cpp
    forward/fd_total_mortality_comp.cpp
    forward/fd_total_obs_catch_age.cpp
    forward/fd_tridag_bet.cpp
)

set(CALPOP_SRC
    caldia.cpp
    Calpop_caldia.cpp
    Calpop_calrec.cpp
    Calpop_InitCalPop.cpp
    Calpop_precaldia.cpp
    Calpop_precalrec.cpp
    Calpop_recompute_coefs.cpp
    Calpop_tridag.cpp
    calrec_adre.cpp
    calrec_precalrec.cpp
    precalrec_juv.cpp
    predicted_catch_without_effort.cpp
    predicted_catch.cpp
    total_exploited_biomass.cpp
    total_mortality_comp.cpp
    total_obs_catch_age.cpp
    tridag_bet.cpp
)

add_library(calpop_fd ${CALPOP_SRC} ${CALPOP_FD})
target_link_libraries(calpop_fd creadwrite_fd varsimtunafunc_fd)

add_library(calpop_dv ${CALPOP_SRC} ${CALPOP_DV})
target_link_libraries(calpop_dv creadwrite_dv varsimtunafunc_dv)


























