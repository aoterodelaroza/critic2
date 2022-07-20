## project and version
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

## some project variables
set(critic2_MAJOR_VERSION 1)
set(critic2_MINOR_VERSION 1)
set(critic2_URL "https://aoterodelaroza.github.io/critic2/")
set(critic2_EMAIL "aoterodelaroza@gmail.com")
set(critic2_LICENSE "GNU/GPL version 3")
set(critic2_DESCRIPTION "Analysis of quantum mechanical calculations in molecules and solids")
string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%S")

## project definition
project(critic2 VERSION "${critic2_MAJOR_VERSION}.${critic2_MINOR_VERSION}" LANGUAGES C Fortran)

## global flags for cmake
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(GNUInstallDirs)
include(utils)

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
option(USE_EXTERNAL_LAPACK "Use external LAPACK/BLAS library (linear algebra operations)." ON)
#option(USE_EXTERNAL_QHULL "Use external QHULL library (calculation of convex hulls)." ON)
#option(USE_EXTERNAL_SPGLIB "Use external spglib library (crystal symmetry determination)." ON)
option(USE_LIBXC "Use the libxc library for access (exchange-correlation energies and potentials)." ON)
option(USE_LIBCINT "Use the libcint library (calculation of molecular integrals with Gaussian functions)." ON)
option(USE_READLINE "Use the readline library (improved command-line interface)." ON)
option(ENABLE_OPENMP "Enable OpenMP parallelization." ON)
option(BUILD_STATIC "Build critic2 statically (tested only on linux, and no GUI)." OFF)

## testing
option(BUILD_TESTING "Enable the regression tests for the critic2 build." OFF)
include(CTest)

## the GUI needs -fPIC
if (ENABLE_GUI)
  enable_language(CXX)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

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
  if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_EQUAL 20.2.1.20201112)
    set(FFLAGS_DEBUG "-g -O0 -traceback -debug all -fpp")
  elseif ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 19.0)
    set(_message "\nA version of intel fortran >=19.5 is required to compile the development version of critic2. "
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
  if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 6.0.0)
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
  if (OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS_${BT} "${CMAKE_Fortran_FLAGS_${BT}} ${OpenMP_Fortran_FLAGS}")
  endif()
endif()

## lapack
if (USE_EXTERNAL_LAPACK)
  set(LAPACK_FIND_QUIETLY "ON")
  find_package(LAPACK)
endif()

# ## qhull -- deactivated
# if (USE_EXTERNAL_QHULL)
#   find_package(QHULL)
# endif()

# ## spglib -- deactivated because the interface keeps changing!
# if (USE_EXTERNAL_SPGLIB)
#   find_package(SPGLIB)
# endif()

## numdiff
find_package(NUMDIFF)

## libxc
if (USE_LIBXC)
  find_package(LIBXC)
  if (LIBXC_FOUND)
    set(MINLIBXC_VERSION 4.1)
    if (${LIBXC_VERSION} VERSION_LESS ${MINLIBXC_VERSION})
      unset(LIBXC_FOUND)
      message("!! Libxc version found (${LIBXC_VERSION},${LIBXC_xc_LIBRARY}) too old (requires ${MINLIBXC_VERSION})")
    endif()
    unset(MINLIBXC_VERSION)
  endif()
endif()

## readline
if (USE_READLINE)
  find_package(READLINE)
endif()

## libcint
if (USE_LIBCINT)
  find_package(LIBCINT)
endif()

## the syntax.txt target
add_custom_target(syntax 
  COMMAND wget -qO - 'raw.githubusercontent.com/aoterodelaroza/aoterodelaroza.github.io/master/_critic2/13_syntax.md' | awk '/~~~/{a=!a\;next}a{print}' > syntax.txt
  )

## process subdirectories
add_subdirectory(src)
add_subdirectory(dat)
if (ENABLE_GUI)
  add_subdirectory(src-gui)
endif()
if (BUILD_TESTING AND NUMDIFF_FOUND)
  add_subdirectory(tests)
endif()

## summary
set(_message
"\n\n### ${CMAKE_PROJECT_NAME}, ${CMAKE_PROJECT_VERSION} (git:${critic2_GIT_COMMIT})###\n"
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
# if (QHULL_FOUND)
#   set(_message "${_message}Qhull: lib=${QHULL_LIBRARIES} | inc=${QHULL_INCLUDE_DIRS}\n")
# else()
#   set(_message "${_message}Qhull: <internal>\n")
# endif()
# if (SPGLIB_FOUND)
#   set(_message "${_message}Spglib: lib=${SPGLIB_LIBRARY}\n")
# else()
#   set(_message "${_message}Spglib: <internal>\n")
# endif()
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
if (READLINE_FOUND)
  set(_message "${_message}Readline: lib=${READLINE_LIBRARY} | inc=${READLINE_INCLUDE_DIRS}\n")
else()
  set(_message "${_message}Readline: not used\n")
endif()
if (NUMDIFF_FOUND)
  set(_message "${_message}Numdiff: ${NUMDIFF_EXE}\n")
else()
  set(_message "${_message}Numdiff: not found\n")
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
