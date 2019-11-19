# Find the Libcint library. Sets:
#  LIBCINT_FOUND - true if the libcint library was found.
#  LIBCINT_INCLUDE_DIRS - the directory containing the cint.h header.
#  LIBCINT_LIBRARY - the location of the libcint library.
#  LIBCINT_VERSION - version of the libcint library.
#
include(FindPackageHandleStandardArgs)

if (DEFINED ENV{LIBCINT_DIR})
  set(LIBCINT_DIR "$ENV{LIBCINT_DIR}")
endif()

find_path(LIBCINT_INCLUDE_DIRS
          NAMES cint.h
          PATH_SUFFIXES include
          HINTS ${LIBCINT_DIR}
	  )

find_library(LIBCINT_LIBRARY
  NAMES cint
  PATH_SUFFIXES lib
  HINTS ${LIBCINT_DIR}
  )

find_package_handle_standard_args(LIBCINT
  REQUIRED_VARS LIBCINT_LIBRARY LIBCINT_INCLUDE_DIRS
  )

mark_as_advanced(LIBCINT_LIBRARY LIBCINT_INCLUDE_DIRS)
if (NOT LIBCINT_FOUND)
  set(LIBCINT_DIR "" CACHE STRING "Directory containing the libcint library.")
endif()

## Get libcint version
if(LIBCINT_INCLUDE_DIRS)
  file(READ "${LIBCINT_INCLUDE_DIRS}/cint.h" _libcint_version_header)
  string(REGEX MATCH "define[ \t]+CINT_VERSION[ \t]+([0-9\\.]+)" LIBCINT_VERSION "${_libcint_version_header}")
  string(REGEX MATCH "([0-9\\.]+)" LIBCINT_VERSION "${LIBCINT_VERSION}")
  unset(_libcint_version_header)
endif()
