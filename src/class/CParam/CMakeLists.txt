set(CPARAM_DV_SRC
    dv_selectivity.cpp
)

set(CPARAM_FD_SRC
# It is empty for now but if something have to be added in
# "forward only" mode it should be listed here.
)

set(CPARAM_SRC
    Param.cpp
)

add_library(cparam_fd ${CPARAM_SRC} ${CPARAM_FD_SRC})

add_library(cparam_dv ${CPARAM_SRC} ${CPARAM_DV_SRC})