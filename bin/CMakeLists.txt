cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(bf_duke CXX)
set(GITHOME_MSU_PRATT $ENV{GITHOME_MSU_PRATT})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU_PRATT}/commonutils/software)
add_custom_target(extern_msu_boltzmann COMMAND make WORKING_DIRECTORY ${GITHOME_MSU_PRATT}/msu_boltzmann/software)
add_custom_target(extern_msu_eos COMMAND make WORKING_DIRECTORY ${GITHOME_MSU_PRATT}/msu_eos/software)
add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ${GITHOME_MSU_PRATT}/msu_sampler/software)
add_custom_target(extern_hydro2uds COMMAND make WORKING_DIRECTORY ${GITHOME_MSU_PRATT}/bf_duke/software)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU_PRATT}/commonutils/software/include  ${GITHOME_MSU_PRATT}/msu_eos/software/include ${GITHOME_MSU_PRATT}/msu_sampler/software/include ${GITHOME_MSU_PRATT}/msu_boltzmann/software/include ${GITHOME_MSU_PRATT}/bf_duke/software/include)

link_directories(${GITHOME_MSU_PRATT}/commonutils/software/lib ${GITHOME_MSU_PRATT}/msu_eos/software/lib ${GITHOME_MSU_PRATT}/msu_sampler/software/lib ${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib ${GITHOME_MSU_PRATT}/bf_duke/software/lib)

find_library(${GITHOME_MSU_PRATT}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU_PRATT}/msu_eos/software/lib msu_eos)
find_library(${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib msu_boltzmann)
find_library(${GITHOME_MSU_PRATT}/bf_duke/software/lib hydro2uds)
find_library(${GITHOME_MSU_PRATT}/msu_sampler/software/lib msu_sampler)

###########################################################

set(SOURCE_FILES ./hydro2udsmain.cc)
SET_SOURCE_FILES_PROPERTIES(hydro2udsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a)
SET_SOURCE_FILES_PROPERTIES(hydro2udsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib/libmsu_boltzmann.a)
SET_SOURCE_FILES_PROPERTIES(hydro2udsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(hydro2udsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/bf_duke/software/lib/libhydro2uds.a)
SET_SOURCE_FILES_PROPERTIES(hydro2udsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_sampler/software/lib/libmsu_sampler.a)

add_executable(hydro2uds ${SOURCE_FILES})
add_dependencies(hydro2uds extern_msu_boltzmann extern_hydro2uds extern_commonutils extern_msu_sampler extern_msu_eos)
target_link_libraries(hydro2uds GSL::gsl ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY}
${GITHOME_MSU_PRATT}/bf_duke/software/lib/libhydro2uds.a
${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib/libmsu_boltzmann.a
${GITHOME_MSU_PRATT}/msu_sampler/software/lib/libmsu_sampler.a
${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a
${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a)

###########################################################
set(SOURCE_FILES ./boltzmain.cc)
SET_SOURCE_FILES_PROPERTIES(boltzmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a)
SET_SOURCE_FILES_PROPERTIES(boltzmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib/libmsu_boltzmann.a)
SET_SOURCE_FILES_PROPERTIES(boltzmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(boltzmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/bf_duke/software/lib/libhydro2uds.a)
SET_SOURCE_FILES_PROPERTIES(boltzmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_sampler/software/lib/libmsu_sampler.a)
SET_SOURCE_FILES_PROPERTIES(boltzmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a)
add_executable(boltz ${SOURCE_FILES})
add_dependencies(boltz extern_msu_boltzmann extern_hydro2uds extern_commonutils extern_msu_sampler extern_msu_eos)
target_link_libraries(boltz GSL::gsl ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY}
${GITHOME_MSU_PRATT}/bf_duke/software/lib/libhydro2uds.a
${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib/libmsu_boltzmann.a
${GITHOME_MSU_PRATT}/msu_sampler/software/lib/libmsu_sampler.a
${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a
${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a
)
###########################################################
set(SOURCE_FILES ./boltzudsmain.cc)
SET_SOURCE_FILES_PROPERTIES(boltzudsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a)
SET_SOURCE_FILES_PROPERTIES(boltzudsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib/libmsu_boltzmann.a)
SET_SOURCE_FILES_PROPERTIES(boltzudsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(boltzudsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/bf_duke/software/lib/libhydro2uds.a)
SET_SOURCE_FILES_PROPERTIES(boltzudsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_sampler/software/lib/libmsu_sampler.a)
SET_SOURCE_FILES_PROPERTIES(boltzudsmain.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a)
add_executable(boltzuds ${SOURCE_FILES})
add_dependencies(boltzuds extern_msu_boltzmann extern_hydro2uds extern_commonutils extern_msu_sampler extern_msu_eos)
target_link_libraries(boltzuds GSL::gsl ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY}
${GITHOME_MSU_PRATT}/bf_duke/software/lib/libhydro2uds.a
${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib/libmsu_boltzmann.a
${GITHOME_MSU_PRATT}/msu_sampler/software/lib/libmsu_sampler.a
${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a
${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a
)
###########################################################
set(SOURCE_FILES ./spectrav2main.cc)
SET_SOURCE_FILES_PROPERTIES(spectrav2main.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a)
SET_SOURCE_FILES_PROPERTIES(spectrav2main.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib/libmsu_boltzmann.a)
SET_SOURCE_FILES_PROPERTIES(spectrav2main.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a)
SET_SOURCE_FILES_PROPERTIES(spectrav2main.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/bf_duke/software/lib/libhydro2uds.a)
SET_SOURCE_FILES_PROPERTIES(spectrav2main.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_sampler/software/lib/libmsu_sampler.a)
SET_SOURCE_FILES_PROPERTIES(spectrav2main.cc
	PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a)
add_executable(spectrav2 ${SOURCE_FILES})
add_dependencies(spectrav2 extern_msu_boltzmann extern_hydro2uds extern_commonutils extern_msu_sampler extern_msu_eos)
target_link_libraries(spectrav2 GSL::gsl ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY}
${GITHOME_MSU_PRATT}/bf_duke/software/lib/libhydro2uds.a
${GITHOME_MSU_PRATT}/msu_boltzmann/software/lib/libmsu_boltzmann.a
${GITHOME_MSU_PRATT}/msu_sampler/software/lib/libmsu_sampler.a
${GITHOME_MSU_PRATT}/msu_eos/software/lib/libmsu_eos.a
${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a
)
###########################################################
set(SOURCE_FILES ./summer_spectrav2.cc)
SET_SOURCE_FILES_PROPERTIES(summer_spectrav2.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a)
add_executable(summer_spectrav2 ${SOURCE_FILES})
add_dependencies(summer_spectrav2 extern_msu_boltzmann extern_hydro2uds extern_commonutils extern_msu_sampler extern_msu_eos)
target_link_libraries(summer_spectrav2 
${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a
)
###########################################################
set(SOURCE_FILES ./summer_bf.cc)
SET_SOURCE_FILES_PROPERTIES(summer_bf.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a)
add_executable(summer_bf ${SOURCE_FILES})
add_dependencies(summer_bf extern_msu_boltzmann extern_hydro2uds extern_commonutils extern_msu_sampler extern_msu_eos)
target_link_libraries(summer_bf 
${GITHOME_MSU_PRATT}/commonutils/software/lib/libmsu_commonutils.a
)
###########################################################
