set(CREADWRITE_SRC
    ReadWrite_DYM.cpp
    ReadWrite_fisheries.cpp
    ReadWrite_TXT.cpp
)

add_library(creadwrite_fd ${CREADWRITE_SRC})
target_link_libraries(
    creadwrite_fd
    cparam_fd
    map_fd
    matrices_fd
    date_dv
    numfunc
    simtunafunc_fd
    savetimearea_fd
)

add_library(creadwrite_dv ${CREADWRITE_SRC})
target_link_libraries(
    creadwrite_dv
    cparam_dv
    map_dv
    matrices_dv
    date_dv
    numfunc
    simtunafunc_dv
    savetimearea_dv
)