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

find_package(PkgConfig QUIET)
pkg_check_modules(PC_TBLITE tblite QUIET)

# the C header (optional; the C API is declared in the Fortran source, but the
# include dir is useful if present)
find_path(TBLITE_INCLUDE_DIRS
  NAMES tblite.h
  HINTS ${PC_TBLITE_INCLUDEDIR} ${PC_TBLITE_INCLUDE_DIRS} $ENV{TBLITE_DIR}
  PATH_SUFFIXES include
  PATHS "${CMAKE_INSTALL_PREFIX}/include")

find_library(TBLITE_LIBRARIES
  NAMES tblite
  HINTS ${PC_TBLITE_LIBDIR} ${PC_TBLITE_LIBRARY_DIRS} $ENV{TBLITE_DIR}
  PATH_SUFFIXES lib lib64)

set(TBLITE_VERSION ${PC_TBLITE_VERSION})

include(FindPackageHandleStandardArgs)
# only the library is strictly required (the C API is bound in the Fortran source)
find_package_handle_standard_args(TBLITE
  FAIL_MESSAGE  DEFAULT_MSG
  REQUIRED_VARS TBLITE_LIBRARIES)

mark_as_advanced(TBLITE_INCLUDE_DIRS TBLITE_LIBRARIES)
