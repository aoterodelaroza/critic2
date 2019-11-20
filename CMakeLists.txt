## project and version
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(critic2 VERSION 1.0 LANGUAGES C CXX Fortran)

## some project variables
set(critic2_URL "https://aoterodelaroza.github.io/critic2/")
set(critic2_EMAIL "aoterodelaroza@gmail.com")
set(critic2_LICENSE "GNU/GPL version 3")
set(critic2_DESCRIPTION "Analysis of quantum mechanical calculations in molecules and solids")
string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%S")

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
option(USE_LIBXC "Use the libxc library, if available." ON)
option(USE_LIBCINT "Use the libcint library, if available." ON)
option(ENABLE_OPENMP "Enable OpenMP parallelization." ON)
option(BUILD_STATIC "Build critic2 statically (tested only on linux, and no GUI)." OFF)

## perform a static build
if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  find_library(LIBDL_LIBRARY NAMES dl)
  if(NOT LIBDL_LIBRARY)
    error("Dynamic linking static library not found")
  endif()
  mark_as_advanced(LIBDL_LIBRARY)
  set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} ${LIBDL_LIBRARY}")
endif()

## compiler flags, fortran, and check versions
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(FFLAGS_DEBUG "-g -O0 -C -traceback -debug all -fpp -check all")
  set(FFLAGS_RELEASE "-O3")
  if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS 19.0.0)
    set(_message "\nA version of gfortran >=6.0 is required to compile the development version of critic2. "
      "Please, visit: https://aoterodelaroza.github.io/critic2/installation/#whichcompilerswork"
      "for more information.\n")
    message(FATAL_ERROR ${_message})
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(FFLAGS_DEBUG "-g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp")
  set(FFLAGS_RELEASE "-O3")
  if (NOT ${CMAKE_CROSSCOMPILING})
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-rdynamic")
  endif()
  if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS 6.0.0)
    set(_message "\nA version of gfortran >=6.0 is required to compile the development version of critic2. "
      "Please, visit: https://aoterodelaroza.github.io/critic2/installation/#whichcompilerswork"
      "for more information.\n")
    message(FATAL_ERROR ${_message})
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

## openmp
if (ENABLE_OPENMP)
  find_package(OpenMP MODULE COMPONENTS Fortran)
endif()

## lapack
if (USE_EXTERNAL_LAPACK)
  find_package(LAPACK)
endif()

## qhull
if (USE_EXTERNAL_QHULL)
  find_package(Qhull)
endif()

## spglib
if (USE_EXTERNAL_SPGLIB)
  find_package(Spglib)
endif()

## libxc
if (USE_LIBXC)
  find_package(Libxc)
  if (LIBXC_FOUND)
    set(MINLIBXC_VERSION 4.1)
    if (${LIBXC_VERSION} VERSION_LESS ${MINLIBXC_VERSION})
      unset(LIBXC_FOUND)
      message("!! Libxc version found (${LIBXC_VERSION},${LIBXC_xc_LIBRARY}) too old (requires ${MINLIBXC_VERSION})")
    endif()
    unset(MINLIBXC_VERSION)
  endif()
endif()

## libcint
if (USE_LIBCINT)
  find_package(Libcint)
endif()

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

## summary
set(_message
"\n\n### ${CMAKE_PROJECT_NAME}, ${CMAKE_PROJECT_VERSION} ###\n"
"System: ${CMAKE_SYSTEM} (crosscompiling? ${CROSSCOMPILE})\n"
"Date: ${DATE}\n"
"Build type: ${CMAKE_BUILD_TYPE}\n"
"Fortran compiler: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID}, ${CMAKE_Fortran_COMPILER_VERSION})\n"
"Fortran flags: ${CMAKE_Fortran_FLAGS_${BT}}\n"
"C compiler: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID}, ${CMAKE_C_COMPILER_VERSION})\n"
"C flags: ${CMAKE_C_FLAGS_${BT}}\n"
"CXX compiler: ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID}, ${CMAKE_CXX_COMPILER_VERSION})\n"
"CXX flags: ${CMAKE_CXX_FLAGS_${BT}}\n"
"Install directory: ${CMAKE_INSTALL_PREFIX}\n"
"OpenMP parallelization: ${OPENMP_FOUND}\n"
)
if (LAPACK_FOUND)
  set(_message "${_message}Lapack: lib=${LAPACK_LIBRARIES}\n")
else()
  set(_message "${_message}Lapack: <internal>\n")
endif()
if (QHULL_FOUND)
  set(_message "${_message}Qhull: lib=${QHULL_LIBRARIES} | inc=${QHULL_INCLUDE_DIRS}\n")
else()
  set(_message "${_message}Qhull: <internal>\n")
endif()
if (SPGLIB_FOUND)
  set(_message "${_message}Spglib: lib=${SPGLIB_LIBRARY}\n")
else()
  set(_message "${_message}Spglib: <internal>\n")
endif()
if (LIBXC_FOUND)
  set(_message "${_message}Libxc: ${LIBXC_VERSION} | lib=${LIBXC_xcf90_LIBRARY} ${LIBXC_xc_LIBRARY} | inc=${LIBXC_INCLUDE_DIRS}\n")
else()
  set(_message "${_message}Libxc: not used\n")
endif()
if (LIBCINT_FOUND)
  set(_message "${_message}Libcint: ${LIBCINT_VERSION} | lib=${LIBCINT_LIBRARY}\n")
else()
  set(_message "${_message}Libcint: not used\n")
endif()

set(_message "${_message}GUI build: ${ENABLE_GUI}\n")

if (ENABLE_GUI)
  set(_message 
    "${_message}gl: lib=${GL_LIBRARIES} | inc=${GL_INCLUDE_DIRS}\n"
    "glm: lib=${GLM_LIBRARIES} | inc=${GLM_INCLUDE_DIRS}\n"
    "glfw: lib=${GLFW_LIBRARIES} | inc=${GLFW_INCLUDE_DIRS}\n"
    "freetype2: lib=${FREETYPE2_LIBRARIES} | inc=${FREETYPE2_INCLUDE_DIRS}\n"
    "dl: lib=${CMAKE_DL_LIBS}\n")
endif()

set(_message "${_message}################################\n")

message(STATUS ${_message})
