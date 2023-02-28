cmake_minimum_required(VERSION 3.0.0)
project(seapodym-cmake VERSION 0.1.0)

# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRUE=true -DFALSE=false -D __GNUDOS__ -Dlinux -O3 -DOPT_LIB -Wall -Wno-deprecated")

# Compiler standard
set(CMAKE_CXX_STANDARD 11)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

# PATH/FIND EXTERNALS LIBRARIES
# ---------------------------

# LIBXML2 
find_package(LibXml2 REQUIRED)

#ADMB
find_library(
    ADMB_LIBRARY
    admbo
    #admb
    PATH_SUFFIXES "admb/lib"
    REQUIRED
)
find_path(ADMB_INCLUDE_DIR "fvar.hpp" PATH_SUFFIXES "admb/include" REQUIRED)

# LOG
# ---

message(STATUS "Path to the libraries :")
message(STATUS "-----------------------")
message(STATUS "\tADMB : ${ADMB_LIBRARY}")
message(STATUS "\tLIBXML2 : ${LIBXML2_LIBRARIES}")
message(STATUS "Path to the headers :")
message(STATUS "---------------------")
message(STATUS "\tADMB : ${ADMB_INCLUDE_DIR}")
message(STATUS "\tLIBXML2 : ${LIBXML2_INCLUDE_DIRS}")
message(STATUS "-----------------------")

# LINKS AND INCLUDES
# ------------------

include_directories(
    "${PROJECT_SOURCE_DIR}/include"
    "${ADMB_INCLUDE_DIR}"
    "${LIBXML2_INCLUDE_DIRS}"
)

link_libraries(
    ${ADMB_LIBRARY}
    ${LIBXML2_LIBRARIES}
)

add_subdirectory(src)