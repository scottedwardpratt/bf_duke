cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(GITHOME_MSU_PRATT $ENV{GITHOME_MSU_PRATT})
project(duke_bf CXX C)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})
include_directories(include ${GITHOME_MSU_PRATT}/commonutils/software/include ${GITHOME_MSU_PRATT}/msu_eos/software/include ${GITHOME_MSU_PRATT}/msu_sampler/software/include ${GITHOME_MSU_PRATT}/msu_boltzmann/software/include)

file(GLOB dukebf_src ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(hydro2uds STATIC ${dukebf_src})
target_link_libraries(hydro2uds ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
message("Libraries: ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}")
