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
get_filename_component(_nlopt_libdir "${NLOPT_LIBRARIES}" DIRECTORY)
if (EXISTS "${_nlopt_libdir}/pkgconfig/nlopt.pc")
  file(STRINGS "${_nlopt_libdir}/pkgconfig/nlopt.pc" _nlopt_pcver REGEX "^Version:")
  if (_nlopt_pcver)
    string(REGEX REPLACE "^Version:[ \t]*" "" NLOPT_VERSION "${_nlopt_pcver}")
  endif()
  unset(_nlopt_pcver)
elseif (CMAKE_CROSSCOMPILING)
  # nothing to read beside the library, and the host pkg-config describes a
  # different nlopt entirely -- report no version rather than a wrong one
  set(NLOPT_VERSION "")
endif()
unset(_nlopt_libdir)

# Set (NAME)_FOUND if all the variables and the version are satisfied.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLOPT
    FAIL_MESSAGE  DEFAULT_MSG
    REQUIRED_VARS NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES
    VERSION_VAR   NLOPT_VERSION)

# hide library and include variables
mark_as_advanced(NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES)
