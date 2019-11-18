## project and version
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(critic2 LANGUAGES C CXX Fortran)

## some project variables
set(critic2_URL "https://aoterodelaroza.github.io/critic2/")
set(critic2_EMAIL "aoterodelaroza@gmail.com")
set(critic2_LICENSE "GNU/GPL version 3")
set(critic2_DESCRIPTION "Analysis of quantum mechanical calculations in molecules and solids")

## global flags for cmake
set(CMAKE_VERBOSE_MAKEFILE ON)
include(GNUInstallDirs)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

## default build type, and type options
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)

## interface options
option(ENABLE_GUI "Enables the critic2 GUI" OFF)
if(NOT CROSSCOMPILE)
  set(CROSSCOMPILE "OFF" CACHE STRING "Choose the cross-compilation toolchain." FORCE)
endif()
set_property(CACHE CROSSCOMPILE PROPERTY STRINGS OFF Win32 Win64)
option(USE_EXTERNAL_LAPACK "Use external LAPACK/BLAS library, if available." ON)
option(USE_EXTERNAL_QHULL "Use external QHULL library, if available." ON)
option(USE_EXTERNAL_SPGLIB "Use external spglib library, if available." ON)
option(ENABLE_OPENMP "Enable OpenMP parallelization." ON)

## openmp
if (ENABLE_OPENMP)
  find_package(OpenMP MODULE COMPONENTS Fortran)
else()
  set(OpenMP_FOUND "")
endif()

## lapack
if (USE_EXTERNAL_LAPACK)
  find_package(LAPACK)
else()
  set(LAPACK_FOUND "")
endif()

## qhull
if (USE_EXTERNAL_QHULL)
  find_package(Qhull)
else()
  set(QHULL_FOUND "")
endif()

## spglib
if (USE_EXTERNAL_SPGLIB)
  find_package(Spglib)
else()
  set(SPGLIB_FOUND "")
endif()

## compiler flags, fortran
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(FFLAGS_DEBUG "-g -O0 -C -traceback -debug all -fpp -check all")
  set(FFLAGS_RELEASE "-O3")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(FFLAGS_DEBUG "-g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp")
  set(FFLAGS_RELEASE "-O3")
  if (NOT ${CMAKE_CROSSCOMPILING})
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-rdynamic")
  endif()
endif()

## compiler flags, C
if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  set(CFLAGS_DEBUG "-D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -g -Wall")
  if (NOT ${CMAKE_CROSSCOMPILING})
    set(CFLAGS_DEBUG "${CFLAGS_DEBUG} -fstack-protector-all")
  endif()
  set(CFLAGS_RELEASE "-O3")
endif()

## compiler flags, CXX
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  set(CXXFLAGS_DEBUG "-D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -g -Wall")
  if (NOT ${CMAKE_CROSSCOMPILING})
    set(CXXFLAGS_DEBUG "${CXXFLAGS_DEBUG} -fstack-protector-all")
  endif()
  set(CXXFLAGS_RELEASE "-O3")
endif()

## set the compiler flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BT)
set(CMAKE_Fortran_FLAGS_${BT} "${FFLAGS_${BT}}")
set(CMAKE_C_FLAGS_${BT} "${CFLAGS_${BT}}")
set(CMAKE_CXX_FLAGS_${BT} "${CXXFLAGS_${BT}}")

## the syntax.txt target
add_custom_target(syntax 
  COMMAND wget -qO - 'raw.githubusercontent.com/aoterodelaroza/aoterodelaroza.github.io/master/_critic2/12_syntax.md' | awk '/~~~/{a=!a\;next}a{print}' > syntax.txt
  )

## process subdirectories
add_subdirectory(src)
add_subdirectory(dat)
if (ENABLE_GUI)
  add_subdirectory(src-gui)
endif()

#### attic ####
## get_cmake_property(_variableNames VARIABLES)
## list (SORT _variableNames)
## foreach (_variableName ${_variableNames})
##     message(STATUS "${_variableName}=${${_variableName}}")
## endforeach()
