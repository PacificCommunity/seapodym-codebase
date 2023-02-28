set(VARSIMTUNAFUNC_FD_SRC
    forward/fd_accessibility.cpp
    forward/fd_feeding_habitat.cpp
    forward/fd_juvenile_habitat.cpp
    forward/fd_mortality_sp.cpp
    forward/fd_seasonal_switch.cpp
    forward/fd_spawning_habitat.cpp
)

set(VARSIMTUNAFUNC_DV_SRC
    adjoint/dv_accessibility.cpp
    adjoint/dv_feeding_habitat.cpp
    adjoint/dv_juvenile_habitat.cpp
    adjoint/dv_mortality_sp.cpp
    adjoint/dv_seasonal_switch.cpp
    adjoint/dv_spawning_habitat.cpp
)

set(VARSIMTUNAFUNC_SRC
    accessibility.cpp
    feeding_habitat.cpp
    juvenile_habitat.cpp
    mortality_sp.cpp
    seasonal_switch.cpp
    spawning_habitat.cpp
)

add_library(varsimtunafunc_fd ${VARSIMTUNAFUNC_SRC} ${VARSIMTUNAFUNC_FD_SRC})
target_link_libraries(varsimtunafunc_fd simtunafunc_fd varparamcoupled_fd)

add_library(varsimtunafunc_dv ${VARSIMTUNAFUNC_SRC} ${VARSIMTUNAFUNC_DV_SRC})
target_link_libraries(varsimtunafunc_dv simtunafunc_dv varparamcoupled_dv)