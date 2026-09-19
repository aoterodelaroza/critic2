# Find the tblite library (https://github.com/tblite/tblite)
#
# critic2 uses the tblite C API for the interactive dynamics feature, so only
# the shared/static library is required (the C bindings are declared directly in
# the Fortran source). A pkg-config file (tblite.pc) is installed by tblite's
# meson build and is used when available.
#
# This sets the following variables:
#   TBLITE_FOUND
#   TBLITE_LIBRARIES
#   TBLITE_INCLUDE_DIRS
#   TBLITE_VERSION

include(FindPackageHandleStandardArgs)

if (DEFINED ENV{TBLITE_DIR})
  set(TBLITE_DIR "$ENV{TBLITE_DIR}")
endif()

find_package(PkgConfig QUIET)
pkg_check_modules(PC_TBLITE tblite QUIET)

## TBLITE_DIR may point outside the cross-compilation sysroot
set(_root_both)
if (TBLITE_DIR)
  set(_root_both CMAKE_FIND_ROOT_PATH_BOTH)
endif()

# the C header (optional; the C API is declared in the Fortran source, but the
# include dir is useful if present)
find_path(TBLITE_INCLUDE_DIRS
  NAMES tblite.h
  HINTS ${TBLITE_DIR} ${PC_TBLITE_INCLUDEDIR} ${PC_TBLITE_INCLUDE_DIRS}
  PATH_SUFFIXES include
  PATHS "${CMAKE_INSTALL_PREFIX}/include"
  ${_root_both})

find_library(TBLITE_LIBRARIES
  NAMES tblite
  HINTS ${TBLITE_DIR} ${PC_TBLITE_LIBDIR} ${PC_TBLITE_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib64
  PATHS "${CMAKE_INSTALL_PREFIX}/lib" "${CMAKE_INSTALL_PREFIX}/lib64"
  ${_root_both})
unset(_root_both)

# Version. PC_TBLITE_VERSION comes from the host's pkg-config, which in a cross
# build describes the host's tblite and not the one found above, so prefer the
# pkgconfig file that sits beside the library actually found.
set(TBLITE_VERSION ${PC_TBLITE_VERSION})
c2_pkgconfig_version(TBLITE_VERSION "${TBLITE_LIBRARIES}" tblite "${PC_TBLITE_LIBDIR}")

# only the library is strictly required (the C API is bound in the Fortran source)
find_package_handle_standard_args(TBLITE
  FAIL_MESSAGE  DEFAULT_MSG
  REQUIRED_VARS TBLITE_LIBRARIES
  VERSION_VAR   TBLITE_VERSION)

mark_as_advanced(TBLITE_INCLUDE_DIRS TBLITE_LIBRARIES)

## Check that the C API critic2 binds is actually there: a tblite built without
## the C API links the library fine but leaves those symbols undefined. Use a
## separate result variable: the normal TBLITE_FOUND set by
## find_package_handle_standard_args above would otherwise shadow the cache
## value written by try_compile, so the check result would be silently ignored.
if (TBLITE_FOUND)
  try_compile(TBLITE_COMPILES "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/Modules/tblite_test.f90"
    LINK_LIBRARIES ${TBLITE_LIBRARIES})
  if (NOT TBLITE_COMPILES)
    message(STATUS "Found tblite (${TBLITE_LIBRARIES}) but could not link against it (missing C API?)")
    set(TBLITE_FOUND FALSE)
  endif()
endif()

if (NOT TBLITE_FOUND)
  set(TBLITE_DIR "${TBLITE_DIR}" CACHE STRING "Directory containing the tblite library.")
endif()
