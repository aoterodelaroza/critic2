# Find the xtb library (https://github.com/grimme-lab/xtb)
#
# critic2 uses the xtb C API for the GFN-FF backend of the interactive
# dynamics feature, so only the shared/static library is required (the C
# bindings are declared directly in the Fortran source). A pkg-config file
# (xtb.pc) is installed by xtb's meson build and is used when available.
#
# This sets the following variables:
#   XTB_FOUND
#   XTB_LIBRARIES
#   XTB_INCLUDE_DIRS
#   XTB_VERSION

include(FindPackageHandleStandardArgs)

if (DEFINED ENV{XTB_DIR})
  set(XTB_DIR "$ENV{XTB_DIR}")
endif()

find_package(PkgConfig QUIET)
pkg_check_modules(PC_XTB xtb QUIET)

## XTB_DIR may point outside the cross-compilation sysroot
set(_root_both)
if (XTB_DIR)
  set(_root_both CMAKE_FIND_ROOT_PATH_BOTH)
endif()

# the C header (optional; the C API is declared in the Fortran source, but the
# include dir is useful if present)
find_path(XTB_INCLUDE_DIRS
  NAMES xtb.h
  HINTS ${XTB_DIR} ${PC_XTB_INCLUDEDIR} ${PC_XTB_INCLUDE_DIRS}
  PATH_SUFFIXES include
  PATHS "${CMAKE_INSTALL_PREFIX}/include"
  ${_root_both})

find_library(XTB_LIBRARIES
  NAMES xtb
  HINTS ${XTB_DIR} ${PC_XTB_LIBDIR} ${PC_XTB_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib64
  PATHS "${CMAKE_INSTALL_PREFIX}/lib" "${CMAKE_INSTALL_PREFIX}/lib64"
  ${_root_both})
unset(_root_both)

# Version. PC_XTB_VERSION comes from the host's pkg-config, which in a cross
# build describes the host's xtb and not the one found above, so prefer the
# pkgconfig file that sits beside the library actually found.
set(XTB_VERSION ${PC_XTB_VERSION})
c2_pkgconfig_version(XTB_VERSION "${XTB_LIBRARIES}" xtb "${PC_XTB_LIBDIR}")

# only the library is strictly required (the C API is bound in the Fortran source)
find_package_handle_standard_args(XTB
  FAIL_MESSAGE  DEFAULT_MSG
  REQUIRED_VARS XTB_LIBRARIES
  VERSION_VAR   XTB_VERSION)

mark_as_advanced(XTB_INCLUDE_DIRS XTB_LIBRARIES)

## Check that the C API critic2 binds is actually there: an xtb built without
## the C API, or too old to have the GFN-FF entry points, links the library
## fine but leaves those symbols undefined. Use a separate result variable:
## the normal XTB_FOUND set by find_package_handle_standard_args above would
## otherwise shadow the cache value written by try_compile, so the check
## result would be silently ignored.
if (XTB_FOUND)
  try_compile(XTB_COMPILES "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/Modules/xtb_test.f90"
    LINK_LIBRARIES ${XTB_LIBRARIES})
  if (NOT XTB_COMPILES)
    message(STATUS "Found xtb (${XTB_LIBRARIES}) but could not link against it (missing C API?)")
    set(XTB_FOUND FALSE)
  endif()
endif()

if (NOT XTB_FOUND)
  set(XTB_DIR "${XTB_DIR}" CACHE STRING "Directory containing the xtb library.")
endif()
