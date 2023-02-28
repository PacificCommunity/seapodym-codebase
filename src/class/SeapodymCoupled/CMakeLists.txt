set(SEAPODYMCOUPLED_FD_SRC
    forward/fd_food_requirement_index.cpp
    forward/fd_spawning.cpp
    forward/fd_survival.cpp
    forward/fd_total_pop.cpp
)

set(SEAPODYMCOUPLED_DV_SRC
    adjoint/dv_food_requirement_index.cpp
    adjoint/dv_spawning.cpp
    adjoint/dv_survival.cpp
    adjoint/dv_total_pop.cpp
)

set(SEAPODYMCOUPLED_SRC
    food_requirement_index.cpp
    like.cpp
    Seapodym_OnRunDensity.cpp
    Seapodym_OnRunHabitat.cpp
    SeapodymCoupled_EditRunCoupled.cpp
    SeapodymCoupled_Forage.cpp
    SeapodymCoupled_Funcs.cpp
    SeapodymCoupled_OnCompFluxes.cpp
    SeapodymCoupled_OnReadForcing.cpp
    SeapodymCoupled_OnRunCoupled.cpp
    SeapodymCoupled_OnRunFirstStep.cpp
    SeapodymCoupled_OnWriteOutput.cpp
    SeapodymCoupled_ReadTags.cpp
    spawning.cpp
)

add_library(seapodymcoupled_fd ${SEAPODYMCOUPLED_SRC} ${SEAPODYMCOUPLED_FD_SRC})
target_link_libraries(seapodymcoupled_fd date_fd calpop_dv savetimearea_fd numfunc)

add_library(seapodymcoupled_dv ${SEAPODYMCOUPLED_SRC} ${SEAPODYMCOUPLED_DV_SRC})
target_link_libraries(seapodymcoupled_dv date_dv calpop_dv savetimearea_dv numfunc)