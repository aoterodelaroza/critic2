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

find_package(PkgConfig QUIET)
pkg_check_modules(PC_XTB xtb QUIET)

# the C header (optional; the C API is declared in the Fortran source, but the
# include dir is useful if present)
find_path(XTB_INCLUDE_DIRS
  NAMES xtb.h
  HINTS ${PC_XTB_INCLUDEDIR} ${PC_XTB_INCLUDE_DIRS} $ENV{XTB_DIR}
  PATH_SUFFIXES include
  PATHS "${CMAKE_INSTALL_PREFIX}/include")

find_library(XTB_LIBRARIES
  NAMES xtb
  HINTS ${PC_XTB_LIBDIR} ${PC_XTB_LIBRARY_DIRS} $ENV{XTB_DIR}
  PATH_SUFFIXES lib lib64
  PATHS "${CMAKE_INSTALL_PREFIX}/lib" "${CMAKE_INSTALL_PREFIX}/lib64")

set(XTB_VERSION ${PC_XTB_VERSION})

include(FindPackageHandleStandardArgs)
# only the library is strictly required (the C API is bound in the Fortran source)
find_package_handle_standard_args(XTB
  FAIL_MESSAGE  DEFAULT_MSG
  REQUIRED_VARS XTB_LIBRARIES)

mark_as_advanced(XTB_INCLUDE_DIRS XTB_LIBRARIES)
