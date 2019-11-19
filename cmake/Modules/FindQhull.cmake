# Find the Spglib library. Sets:
# QHULL_FOUND - True if QHULL was found.
# QHULL_INCLUDE_DIRS - Directories containing the QHULL include files.
# QHULL_LIBRARIES - Libraries needed to use QHULL.
#
include(FindPackageHandleStandardArgs)

if (DEFINED ENV{QHULL_DIR})
  set(QHULL_DIR "$ENV{QHULL_DIR}")
endif()

find_path(QHULL_INCLUDE_DIRS
          NAMES libqhull.h
          HINTS "${QHULL_DIR}"
          PATH_SUFFIXES qhull src/libqhull libqhull include)

find_library(QHULL_LIBRARIES
             NAMES qhull
             HINTS "${QHULL_DIR}"
             PATH_SUFFIXES project build bin lib)

find_package_handle_standard_args(QHULL
  REQUIRED_VARS QHULL_LIBRARIES QHULL_INCLUDE_DIRS
  )

mark_as_advanced(QHULL_INCLUDE_DIRS QHULL_LIBRARIES)
if (NOT QHULL_FOUND)
  set(QHULL_DIR "" CACHE STRING "Directory containing the qhull library.")
endif()

