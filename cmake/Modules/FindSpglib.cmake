# Find the Spglib library. Sets:
#  SPGLIB_FOUND - true if the spglib library was found
#  SPGLIB_LIBRARIES - the location of the spglib library
#
include(FindPackageHandleStandardArgs)

if (DEFINED ENV{SPGLIB_DIR})
  set(SPGLIB_DIR "$ENV{SPGLIB_DIR}")
endif()

find_library(SPGLIB_LIBRARY
  NAMES spglib symspg
  PATH_SUFFIXES project build bin lib
  HINTS ${SPGLIB_DIR}
  )

find_package_handle_standard_args(SPGLIB
  REQUIRED_VARS SPGLIB_LIBRARY
  )

mark_as_advanced(SPGLIB_LIBRARY)
if (NOT SPGLIB_FOUND)
  set(SPGLIB_DIR "" CACHE STRING "Directory containing the spglib library.")
endif()
