## project and version
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(critic2 LANGUAGES C CXX Fortran)

## verbose make and other global things
set(CMAKE_VERBOSE_MAKEFILE ON)
include(GNUInstallDirs)
## list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## some project variables
set(critic2_URL "https://aoterodelaroza.github.io/critic2/")
set(critic2_EMAIL "aoterodelaroza@gmail.com")
set(critic2_LICENSE "GNU/GPL version 3")
set(critic2_DESCRIPTION "Analysis of quantum mechanical calculations in molecules and solids")

## interface options
option(ENABLE_GUI "Enables the critic2 GUI" OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release) ## only debug and release builds

## openmp
find_package(OpenMP MODULE COMPONENTS Fortran)

## compiler flags, fortran
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
   set(FFLAGS_DEBUG "-g -O0 -C -traceback -debug all -fpp -check all")
   set(FFLAGS_RELEASE "-O3")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
   set(FFLAGS_DEBUG "-g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp")
   set(FFLAGS_RELEASE "-O3")
   set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-rdynamic")
endif()

## compiler flags, C
if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  set(CFLAGS_DEBUG "-D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -fstack-protector-all -g -Wall")
  set(CFLAGS_RELEASE "-O3")
endif()

## compiler flags, CXX
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  set(CXXFLAGS_DEBUG "-D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -fstack-protector-all -g -Wall")
  set(CXXFLAGS_RELEASE "-O3")
endif()

## set the compiler flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BT)
set(CMAKE_Fortran_FLAGS_${BT} "${FFLAGS_${BT}}")
set(CMAKE_C_FLAGS_${BT} "${CFLAGS_${BT}}")
set(CMAKE_CXX_FLAGS_${BT} "${CXXFLAGS_${BT}}")

add_subdirectory(src)
add_subdirectory(dat)
if (ENABLE_GUI)
  add_subdirectory(src-gui)
endif()

## get_cmake_property(_variableNames VARIABLES)
## list (SORT _variableNames)
## foreach (_variableName ${_variableNames})
##     message(STATUS "${_variableName}=${${_variableName}}")
## endforeach()
