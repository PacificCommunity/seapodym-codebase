add_subdirectory(class)

# SEAPODYM MAIN
set(FUNCTIONS_SRC
    ad_buffers.cpp
    hessian.cpp
    seapodym_coupled.cpp
)
add_library(functions ${FUNCTIONS_SRC})
target_link_libraries(functions seapodymcoupled_dv)

# SEAPODYM DENSITY
set(FUNCTIONSDENSITY_SRC
    ad_buffers.cpp
    hessian.cpp
    seapodym_density.cpp
)
add_library(functionsdensity ${FUNCTIONSDENSITY_SRC})
target_link_libraries(functionsdensity seapodymcoupled_dv)

# SEAPODYM HABITAT
set(FUNCTIONSHABITAT_SRC
    ad_buffers.cpp
    hessian.cpp
    seapodym_habitat.cpp
)
add_library(functionshabitat ${FUNCTIONSHABITAT_SRC})
target_link_libraries(functionshabitat seapodymcoupled_dv)

# SEAPODYM SIMULATION ONLY
set(FUNCTIONSSIMULATION_SRC
    ad_buffers.cpp
    hessian.cpp
    seapodym_coupled.cpp
)
add_library(functionssimulation ${FUNCTIONSSIMULATION_SRC})
target_link_libraries(functionssimulation seapodymcoupled_fd)

# SEAPODYM FLUX
set(FUNCTIONSFLUX_SRC
    ad_buffers.cpp
    hessian.cpp
    seapodym_flux.cpp
)
add_library(functionsflux ${FUNCTIONSFLUX_SRC})
target_link_libraries(functionsflux seapodymcoupled_fd)


# PRODUCE LOCAL LIBRARIES AND EXECUTABLES
# ---------------------------------------

# SEAPODYM MAIN
add_executable(seapodym-cmake main.cpp)
target_link_libraries(
    seapodym-cmake
    functions
    numfunc
    cparam_dv
    date_dv
    map_dv
    matrices_dv
    simtunafunc_dv
    savetimearea_dv
    creadwrite_dv
    varparamcoupled_dv
    varsimtunafunc_dv
    calpop_dv
    functions
)

# SEAPODYM HABITAT
add_executable(seapodym-cmake-habitats main_habitats.cpp)
target_link_libraries(
    seapodym-cmake-habitats
    numfunc
    cparam_dv
    date_dv
    map_dv
    matrices_dv
    simtunafunc_dv
    savetimearea_dv
    creadwrite_dv
    varparamcoupled_dv
    varsimtunafunc_dv
    calpop_dv
    functionshabitat
)

# SEAPODYM DENSITY
add_executable(seapodym-cmake-densities main_densities.cpp)
target_link_libraries(
    seapodym-cmake-densities
    numfunc
    cparam_dv
    date_dv
    map_dv
    matrices_dv
    simtunafunc_dv
    savetimearea_dv
    creadwrite_dv
    varparamcoupled_dv
    varsimtunafunc_dv
    calpop_dv
    functionsdensity
)

# SEAPODYM SIMULATION ONLY
add_executable(seapodym-cmake-simulation-only main_simulation.cpp)
target_link_libraries(
    seapodym-cmake-simulation-only
    numfunc
    cparam_fd
    date_fd
    map_fd
    matrices_fd
    simtunafunc_fd
    savetimearea_fd
    creadwrite_fd
    varparamcoupled_fd
    varsimtunafunc_fd
    calpop_fd
    functionssimulation
)

# SEAPODYM FLUX
add_executable(seapodym-cmake-fluxes main_flux.cpp)
target_link_libraries(
    seapodym-cmake-fluxes
    numfunc
    cparam_fd
    date_fd
    map_fd
    matrices_fd
    simtunafunc_fd
    savetimearea_fd
    creadwrite_fd
    varparamcoupled_fd
    varsimtunafunc_fd
    calpop_fd
    functionsflux
)