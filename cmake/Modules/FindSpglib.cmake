# Find the Spglib library. This module defines:
#
#  SPGLIB_FOUND - true if the spglib library was found
#  SPGLIB_LIBRARIES - the location of the spglib library
#
find_library(SPGLIB_LIBRARY
  NAMES spglib symspg
  PATH_SUFFIXES project build bin lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPGLIB
  REQUIRED_VARS SPGLIB_LIBRARY
  )

mark_as_advanced(SPGLIB_LIBRARY)
