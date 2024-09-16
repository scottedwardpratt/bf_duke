cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(andrew CXX)
set(GITHOME_MSU_PRATT $ENV{GITHOME_MSU_PRATT})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU_PRATT}/commonutils/software)
add_custom_target(extern_msu_eos COMMAND make WORKING_DIRECTORY ${GITHOME_MSU_PRATT}/msu_eos/software)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU_PRATT}/commonutils/software/include  ${GITHOME_MSU_PRATT}/msu_eos/software/include ${GITHOME_MSU_PRATT}/msu_sampler/software/include)

link_directories(${GITHOME_MSU_PRATT}/commonutils/software/lib ${GITHOME_MSU_PRATT}/msu_eos/software/lib)

find_library(${GITHOME_MSU_PRATT}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU_PRATT}/msu_eos/software/lib msu_eos)

###########################################################
set(SOURCE_FILES ./andrew.cc)
SET_SOURCE_FILES_PROPERTIES(./andrew.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a)
SET_SOURCE_FILES_PROPERTIES(./andrew.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a)
add_executable(andrew ${SOURCE_FILES})
add_dependencies(andrew extern_commonutils extern_msu_eos)
target_link_libraries(andrew ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY}
${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a
${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a
)