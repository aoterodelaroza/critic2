# Copyright (c) 2011-2023, The DART development contributors
# All rights reserved.
#
# The list of contributors can be found at:
#   https://github.com/dartsim/dart/blob/master/LICENSE
#
# This file is provided under the "BSD-style" License

# Find NLOPT
#
# This sets the following variables:
#   NLOPT_FOUND
#   NLOPT_INCLUDE_DIRS
#   NLOPT_LIBRARIES
#   NLOPT_DEFINITIONS
#   NLOPT_VERSION
#
# and the following targets:
#   NLOPT::nlopt

find_package(PkgConfig QUIET)

# Check to see if pkgconfig is installed.
pkg_check_modules(PC_NLOPT nlopt QUIET)

# NLOPT_DIR may point outside the cross-compilation sysroot (same convention as
# FindLIBXC/FindLIBCINT/FindREADLINE in this directory)
if (DEFINED ENV{NLOPT_DIR})
  set(NLOPT_DIR "$ENV{NLOPT_DIR}")
endif()
set(_root_both)
if (NLOPT_DIR)
  set(_root_both CMAKE_FIND_ROOT_PATH_BOTH)
endif()

# Definitions
set(NLOPT_DEFINITIONS ${PC_NLOPT_CFLAGS_OTHER})

# Include directories
find_path(NLOPT_INCLUDE_DIRS
    NAMES nlopt.h
    PATH_SUFFIXES include
    HINTS ${NLOPT_DIR} ${PC_NLOPT_INCLUDEDIR}
    ${_root_both}
    PATHS "${CMAKE_INSTALL_PREFIX}/include")

# Libraries
find_library(NLOPT_LIBRARIES
    NAMES nlopt nlopt_cxx
    PATH_SUFFIXES lib
    HINTS ${NLOPT_DIR} ${PC_NLOPT_LIBDIR}
    ${_root_both})
unset(_root_both)

# Version. PC_NLOPT_VERSION comes from the host's pkg-config, which in a
# cross build describes the host's nlopt and not the one found above, so
# prefer the pkgconfig file that sits beside the library actually found.
set(NLOPT_VERSION ${PC_NLOPT_VERSION})
c2_pkgconfig_version(NLOPT_VERSION "${NLOPT_LIBRARIES}" nlopt "${PC_NLOPT_LIBDIR}")

# Set (NAME)_FOUND if all the variables and the version are satisfied.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLOPT
    FAIL_MESSAGE  DEFAULT_MSG
    REQUIRED_VARS NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES
    VERSION_VAR   NLOPT_VERSION)

# hide library and include variables
mark_as_advanced(NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES)

## critic2 uses the F77 interface (include 'nlopt.f'), which nlopt only
## installs when built with NLOPT_FORTRAN (Homebrew's nlopt, for one, is
## not). The nlo_* wrappers are in the library regardless, so if nlopt.f is
## missing, generate it from nlopt.h the way nlopt's own build does: the
## enumerators of nlopt_algorithm and nlopt_result, minus the NLOPT_NUM_*
## counters.
if (NLOPT_FOUND)
  find_path(NLOPT_F77_INCLUDE_DIR NAMES nlopt.f HINTS ${NLOPT_INCLUDE_DIRS} NO_DEFAULT_PATH)
  mark_as_advanced(NLOPT_F77_INCLUDE_DIR)
  if (NOT NLOPT_F77_INCLUDE_DIR)
    set(_nlopt_f77_dir "${CMAKE_BINARY_DIR}/nlopt_f77")
    file(STRINGS "${NLOPT_INCLUDE_DIRS}/nlopt.h" _nlopt_lines)
    set(_nlopt_f77 "")
    set(_inenum FALSE)
    set(_ival -1)
    foreach (_line IN LISTS _nlopt_lines)
      if (_line MATCHES "^[ \t]*typedef[ \t]+enum")
        set(_inenum TRUE)
        set(_ival -1)
      elseif (_inenum AND _line MATCHES "^[ \t]*}")
        set(_inenum FALSE)
      elseif (_inenum AND _line MATCHES "^[ \t]*(NLOPT_[A-Z0-9_]+)[ \t]*(=[ \t]*(-?[0-9]+))?")
        set(_name "${CMAKE_MATCH_1}")
        if (CMAKE_MATCH_3)
          set(_ival "${CMAKE_MATCH_3}")
        else()
          math(EXPR _ival "${_ival} + 1")
        endif()
        if (NOT _name MATCHES "^NLOPT_NUM_")
          string(APPEND _nlopt_f77 "      integer ${_name}\n      parameter (${_name}=${_ival})\n")
        endif()
      endif()
    endforeach()
    file(WRITE "${_nlopt_f77_dir}/nlopt.f" "${_nlopt_f77}")
    list(APPEND NLOPT_INCLUDE_DIRS "${_nlopt_f77_dir}")
    message(STATUS "nlopt.f not installed with nlopt; generated ${_nlopt_f77_dir}/nlopt.f from nlopt.h")
    unset(_nlopt_f77_dir)
    unset(_nlopt_lines)
    unset(_nlopt_f77)
    unset(_inenum)
    unset(_ival)
    unset(_name)
  endif()

  ## check that the F77 interface compiles and links (separate result
  ## variable, see FindLIBXC)
  try_compile(NLOPT_COMPILES "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/Modules/nlopt_test.f90"
    LINK_LIBRARIES ${NLOPT_LIBRARIES}
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${NLOPT_INCLUDE_DIRS}")
  if (NOT NLOPT_COMPILES)
    message(STATUS "Found nlopt (lib=${NLOPT_LIBRARIES} | inc=${NLOPT_INCLUDE_DIRS}) but could not compile against its Fortran interface")
    set(NLOPT_FOUND FALSE)
  endif()
endif()
