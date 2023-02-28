add_library(numfunc Numfunc.cpp)

add_subdirectory(CParam)

add_library(date_fd Date.cpp)
target_link_libraries(date_fd cparam_fd)
add_library(date_dv Date.cpp)
target_link_libraries(date_dv cparam_dv)

add_library(map_fd Map.cpp)
target_link_libraries(map_fd cparam_fd)
add_library(map_dv Map.cpp)
target_link_libraries(map_dv cparam_dv)

add_library(matrices_fd Matrices.cpp)
target_link_libraries(matrices_fd cparam_fd map_fd)
add_library(matrices_dv Matrices.cpp)
target_link_libraries(matrices_dv cparam_dv map_dv)

add_library(simtunafunc_fd SimtunaFunc.cpp)
target_link_libraries(simtunafunc_fd matrices_fd)
add_library(simtunafunc_dv SimtunaFunc.cpp)
target_link_libraries(simtunafunc_dv matrices_dv)

add_library(savetimearea_fd SaveTimeArea.cpp)
target_link_libraries(savetimearea_fd simtunafunc_fd)
add_library(savetimearea_dv SaveTimeArea.cpp)
target_link_libraries(savetimearea_dv simtunafunc_dv)

add_subdirectory(CReadWrite)

add_subdirectory(VarParamCoupled)

add_subdirectory(VarSimtunaFunc)

add_subdirectory(Calpop)

add_subdirectory(SeapodymCoupled)



