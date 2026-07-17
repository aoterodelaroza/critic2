## project and version
cmake_minimum_required(VERSION 3.16...3.29 FATAL_ERROR)

## some project variables
file(READ "VERSION" critic2_VERSION)
string(REGEX REPLACE " *\n$" "" critic2_VERSION "${critic2_VERSION}")
set(critic2_URL "https://aoterodelaroza.github.io/critic2/")
set(critic2_EMAIL "aoterodelaroza@gmail.com")
set(critic2_LICENSE "GNU/GPL version 3")
set(critic2_DESCRIPTION "Analysis of quantum mechanical calculations in molecules and solids")
string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%S")

## For Windows cross-compilation, see cmake/toolchains/x86_64-w64-mingw32.cmake.
## For native MSYS2 builds, no toolchain file is needed (use -G "MSYS Makefiles").

## project definition
project(critic2 VERSION "${critic2_VERSION}"
  DESCRIPTION "${critic2_DESCRIPTION}"
  HOMEPAGE_URL "${critic2_URL}"
  LANGUAGES C Fortran)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_Fortran_VISIBILITY_PRESET hidden)

## global flags for cmake
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(GNUInstallDirs)
include(utils)

## The install prefix in the cache
set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "Installation directory.")

## default build type, and type options
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
else()
  set(CMAKE_BUILD_TYPE CACHE STRING "Choose the type of build.")
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Profile Release)

## interface options
option(BUILD_LIBRARY "Build the critic2 shared library." OFF)
option(ENABLE_GUI "Enables the critic2 GUI" OFF)
option(USE_EXTERNAL_LAPACK "Use external LAPACK/BLAS library (linear algebra operations)." ON)
option(USE_LIBXC "Use the libxc library for access (exchange-correlation energies and potentials)." ON)
option(USE_LIBCINT "Use the libcint library (calculation of molecular integrals with Gaussian functions)." ON)
option(USE_READLINE "Use the readline library (improved command-line interface)." ON)
option(USE_NLOPT "Use the nlopt library (global optimization for variable-cell crystal structure comparison)." ON)
option(USE_TBLITE "Use the tblite library (GFN-xTB energies and forces for the interactive dynamics)." OFF)
option(ENABLE_OPENMP "Enable OpenMP parallelization." ON)
option(BUILD_STATIC "Build critic2 statically (tested only on linux, and no GUI)." OFF)

## testing
option(BUILD_TESTING "Enable the regression tests for the critic2 build." OFF)
include(CTest)

## library requires fpic
if (BUILD_LIBRARY)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

## GUI requirements
if (ENABLE_GUI)
  ## GUI options
  option(USE_HDF5 "Use the HDF5 library (reading hdf5 phonopy vibration data)." ON)
  option(USE_FREETYPE "Use the freetype library (improved font rendering in GUI)." ON)
  option(USE_GUI_THREADS "Use threads in the GUI." ON)

  ## needs -fPIC
  enable_language(CXX)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  ## C++ standard for the bundled imgui/implot/ImGuiFileDialog
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)

  ## find GUI libraries
  if (USE_FREETYPE)
    find_package(Freetype)
  endif()

  ## HDF5
  if (USE_HDF5)
    find_package(HDF5 COMPONENTS Fortran HL)
    mark_as_advanced(HDF5_C_COMPILER_EXECUTABLE)
    mark_as_advanced(HDF5_Fortran_COMPILER_EXECUTABLE)
    mark_as_advanced(HDF5_Fortran_HL_stub)
    mark_as_advanced(HDF5_Fortran_stub)
    mark_as_advanced(SZIP_INCLUDE_DIR)
    mark_as_advanced(SZIP_LIBRARY)
    mark_as_advanced(h5_conf)
  endif()

  ## glfw3
  find_package(GLFW3 REQUIRED)

  ## opengl
  set(OpenGL_GL_PREFERENCE "GLVND")
  find_package(OpenGL REQUIRED)

  ## the GUI needs a threads library (tinycthread wraps pthreads)
  if (USE_GUI_THREADS)
    find_package(Threads REQUIRED)
  endif()
endif()

## perform a static build
if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BLA_STATIC ON)
  set(BUILD_SHARED_LIBS OFF)
  find_library(LIBDL_LIBRARY NAMES dl)
  if(NOT LIBDL_LIBRARY)
    message(FATAL_ERROR "libdl not found (required for BUILD_STATIC)")
  endif()
  mark_as_advanced(LIBDL_LIBRARY)
endif()

## compiler flags, fortran, and check versions
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(FFLAGS_DEBUG "-g -O0 -C -traceback -debug all -fpp -check all")
  set(FFLAGS_PROFILE "-g -pg -O0 -C -traceback -debug all -fpp -check all")
  set(FFLAGS_RELEASE "-O3")
  if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_EQUAL 20.2.1.20201112)
    set(FFLAGS_DEBUG "-g -O0 -traceback -debug all -fpp")
    set(FFLAGS_PROFILE "-g -pg -O0 -traceback -debug all -fpp")
  elseif ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 19.0)
    set(_message "\nA version of intel fortran >=19.5 is required to compile the development version of critic2. "
      "Please, visit: https://aoterodelaroza.github.io/critic2/installation/#whichcompilerswork"
      "for more information.\n")
    message(FATAL_ERROR ${_message})
  endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  ## ifx; -check uninit is slow and fragile in current ifx releases
  set(FFLAGS_DEBUG "-g -O0 -traceback -fpp -check all,nouninit")
  set(FFLAGS_PROFILE "-g -pg -O0 -traceback -fpp -check all,nouninit")
  set(FFLAGS_RELEASE "-O3 -fpp")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if (ENABLE_GUI)
    # -ffpe-trap=invalid makes GLFW crash in some versions
    set(FFLAGS_DEBUG "-g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=zero,overflow -fbacktrace -fdump-core -cpp -ffree-line-length-0")
    set(FFLAGS_PROFILE "-g -pg -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=zero,overflow -fbacktrace -fdump-core -cpp -ffree-line-length-0")
  else()
    set(FFLAGS_DEBUG "-g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp -ffree-line-length-0")
    set(FFLAGS_PROFILE "-g -pg -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp -ffree-line-length-0")
  endif()
  set(FFLAGS_RELEASE "-O3 -ffree-line-length-0")
  if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 6.0.0)
    set(_message "\nA version of gfortran >=6.0 is required to compile the development version of critic2. "
      "Please, visit: https://aoterodelaroza.github.io/critic2/installation/#whichcompilerswork"
      "for more information.\n")
    message(FATAL_ERROR ${_message})
  endif()
else()
  message(WARNING "Untested Fortran compiler (${CMAKE_Fortran_COMPILER_ID}); using CMake default flags")
endif()

## compiler flags, C
if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  set(CFLAGS_DEBUG "-fasynchronous-unwind-tables -fexceptions -ggdb -Wall")
  set(CFLAGS_PROFILE "-fasynchronous-unwind-tables -fexceptions -ggdb -Wall -pg")
  set(CFLAGS_RELEASE "-O3")
endif()

## compiler flags, CXX
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  set(CXXFLAGS_DEBUG "-fasynchronous-unwind-tables -fexceptions -ggdb -Wall")
  set(CXXFLAGS_PROFILE "-fasynchronous-unwind-tables -fexceptions -ggdb -Wall -pg")
  set(CXXFLAGS_RELEASE "-O3")
endif()

## set the compiler flags, but only for the configurations defined above;
## compilers not handled there keep CMake's default flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BT)
foreach(_cfg DEBUG PROFILE RELEASE)
  if (DEFINED FFLAGS_${_cfg})
    set(CMAKE_Fortran_FLAGS_${_cfg} "${FFLAGS_${_cfg}}")
  endif()
  if (DEFINED CFLAGS_${_cfg})
    set(CMAKE_C_FLAGS_${_cfg} "${CFLAGS_${_cfg}}")
  endif()
  if (DEFINED CXXFLAGS_${_cfg})
    set(CMAKE_CXX_FLAGS_${_cfg} "${CXXFLAGS_${_cfg}}")
  endif()
endforeach()

## openmp
if (ENABLE_OPENMP)
  find_package(OpenMP MODULE COMPONENTS Fortran)
endif()

## lapack
if (USE_EXTERNAL_LAPACK)
  set(LAPACK_FIND_QUIETLY "ON")
  find_package(LAPACK)
endif()

## numdiff (used by the regression tests). The test harness uses
## sh/grep/sed pipelines and symlinks and numdiff is a host tool, so the
## tests are skipped when targeting Windows or cross-compiling.
if (NOT WIN32 AND NOT CMAKE_CROSSCOMPILING)
  find_package(NUMDIFF)
endif()

## libxc
if (USE_LIBXC)
  find_package(LIBXC)
  if (LIBXC_FOUND)
    set(MINLIBXC_VERSION 5.0)
    if (${LIBXC_VERSION} VERSION_LESS ${MINLIBXC_VERSION})
      set(LIBXC_FOUND FALSE)
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

## nlopt
if (USE_NLOPT)
  find_package(NLOPT)
endif()

## tblite
if (USE_TBLITE)
  find_package(TBLITE)
endif()

## the syntax.txt target
add_custom_target(syntax
  COMMAND wget -qO - 'raw.githubusercontent.com/aoterodelaroza/aoterodelaroza.github.io/master/_critic2/13_syntax.md' | awk '/~~~/{a=!a\;next}a{print}' > syntax.txt
  )

## all Fortran modules go to a single directory; gfortran/ifort also search
## the -J/-module directory, so no extra include paths are needed
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran_modules")

## process subdirectories
add_subdirectory(src)
add_subdirectory(dat)
if (BUILD_TESTING AND NUMDIFF_FOUND)
  add_subdirectory(tests)
endif()

## configuration summary

## append one "label   value" line, with the label padded to a fixed
## width so the values line up. ;-lists (e.g. LAPACK_LIBRARIES) are shown
## space-separated so they stay readable instead of running together.
macro(summary_line label value)
  string(REPLACE ";" " " _v "${value}")
  string(LENGTH "${label}" _n)
  math(EXPR _n "18 - ${_n}")
  set(_pad "")
  if (_n GREATER 0)
    string(REPEAT " " ${_n} _pad)
  endif()
  string(APPEND _summary "  ${label}${_pad}${_v}\n")
endmacro()

## append a library line: its detail when found, "not used" otherwise
macro(summary_pkg label found detail)
  if (${found})
    summary_line("${label}" "${detail}")
  else()
    summary_line("${label}" "not used")
  endif()
endmacro()

set(_summary "\n### ${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_VERSION} ###\n")
summary_line("Date"             "${DATE}")
summary_line("Build type"       "${CMAKE_BUILD_TYPE}")
summary_line("Install prefix"   "${CMAKE_INSTALL_PREFIX}")
summary_line("Fortran compiler" "${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION})")
summary_line("Fortran flags"    "${CMAKE_Fortran_FLAGS_${BT}}")
summary_line("C compiler"       "${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")
summary_line("C flags"          "${CMAKE_C_FLAGS_${BT}}")
if (ENABLE_GUI)
  summary_line("C++ compiler"   "${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
  summary_line("C++ flags"      "${CMAKE_CXX_FLAGS_${BT}}")
endif()
summary_line("Shared library"   "${BUILD_LIBRARY}")
summary_line("Static build"     "${BUILD_STATIC}")
summary_line("OpenMP"           "${OPENMP_FOUND}")

if (LAPACK_FOUND)
  summary_line("Lapack" "${LAPACK_LIBRARIES}")
else()
  summary_line("Lapack" "<internal>")
endif()
summary_pkg("Libxc"    LIBXC_FOUND    "${LIBXC_VERSION} | ${LIBXC_xcf03_LIBRARY} ${LIBXC_xc_LIBRARY} | ${LIBXC_INCLUDE_DIRS}")
summary_pkg("Libcint"  LIBCINT_FOUND  "${LIBCINT_VERSION} | ${LIBCINT_LIBRARY}")
summary_pkg("Readline" READLINE_FOUND "${READLINE_LIBRARY} | ${READLINE_INCLUDE_DIRS}")
summary_pkg("Nlopt"    NLOPT_FOUND    "${NLOPT_VERSION} | ${NLOPT_LIBRARIES}")
summary_pkg("tblite"   TBLITE_FOUND   "${TBLITE_LIBRARIES}")
summary_pkg("Numdiff"  NUMDIFF_FOUND  "${NUMDIFF_EXE}")

summary_line("GUI" "${ENABLE_GUI}")
if (ENABLE_GUI)
  summary_line("OpenGL" "${OPENGL_LIBRARIES}")
  summary_line("GLFW"   "${GLFW3_LIBRARY} | ${GLFW3_INCLUDE_DIR}")
  summary_pkg("Freetype" FREETYPE_FOUND "${FREETYPE_LIBRARIES}")
  summary_pkg("HDF5"     HDF5_FOUND     "${HDF5_LIBRARIES}")
endif()

string(APPEND _summary "###############################")
message(STATUS "${_summary}")

## packaging
include(${CMAKE_SOURCE_DIR}/cmake/CritiCPack.cmake)
