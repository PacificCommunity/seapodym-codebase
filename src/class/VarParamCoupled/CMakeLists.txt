set(VARPARAMCOUPLED_SRC
    VarParamCoupled_reset.cpp
    VarParamCoupled_xinit.cpp
    VarParamCoupled.cpp
)

add_library(varparamcoupled_fd ${VARPARAMCOUPLED_SRC})
target_link_libraries(varparamcoupled_fd cparam_fd creadwrite_fd)

add_library(varparamcoupled_dv ${VARPARAMCOUPLED_SRC})
target_link_libraries(varparamcoupled_dv cparam_dv creadwrite_dv)