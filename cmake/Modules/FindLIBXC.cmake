## Find the libxc library. Sets:
## LIBXC_FOUND - whether the following three items were found
## LIBXC_xc_LIBRARY - the libxc library.
## LIBXC_xcf90_LIBRARY - the libxcf90 library.
## LIBXC_INCLUDE_DIRS - the directory containing the mod files.
## LIBXC_VERSION - the version of the libxc found

include(FindPackageHandleStandardArgs)

if (DEFINED ENV{LIBXC_DIR})
  set(LIBXC_DIR "$ENV{LIBXC_DIR}")
endif()

find_library(LIBXC_xc_LIBRARY
  NAMES xc
  PATH_SUFFIXES lib
  HINTS ${LIBXC_DIR}
  )
find_library(LIBXC_xcf90_LIBRARY
  NAMES xcf90
  PATH_SUFFIXES lib
  HINTS ${LIBXC_DIR}
  )

find_path(LIBXC_INCLUDE_DIRS
  NAMES xc_f90_lib_m.mod
  PATH_SUFFIXES include
  HINTS ${LIBXC_DIR}
  )

find_package_handle_standard_args(LIBXC
  REQUIRED_VARS LIBXC_xc_LIBRARY LIBXC_xcf90_LIBRARY LIBXC_INCLUDE_DIRS
  )

mark_as_advanced(LIBXC_xc_LIBRARY LIBXC_xcf90_LIBRARY LIBXC_INCLUDE_DIRS)
if (NOT LIBXC_FOUND)
  set(LIBXC_DIR "${LIBXC_DIR}" CACHE STRING "Directory containing the libxc library (>=4.1).")
endif()

## check whether we can compile against it
if (LIBXC_FOUND)
  try_compile(LIBXC_FOUND "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/Modules/libxc_test.f90"
    LINK_LIBRARIES ${LIBXC_xcf90_LIBRARY} ${LIBXC_xc_LIBRARY}
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${LIBXC_INCLUDE_DIRS}")
  if (NOT LIBXC_FOUND)
    message(STATUS "Found libxc (lib=${LIBXC_xcf90_LIBRARY} ${LIBXC_xc_LIBRARY} | inc=${LIBXC_INCLUDE_DIRS}) but could not compile against it (different compiler?)")
  endif()
endif()

## Get libxc version, from scibuilder, https://github.com/scibuilder/SciBuilder
if(LIBXC_INCLUDE_DIRS)
  file(READ "${LIBXC_INCLUDE_DIRS}/xc_version.h" _libxc_version_header)
  string(REGEX MATCH "define[ \t]+XC_VERSION[ \t]+\\\"([0-9\\.]+)\\\"" LIBXC_VERSION "${_libxc_version_header}")
  string(REGEX MATCH "([0-9\\.]+)" LIBXC_VERSION "${LIBXC_VERSION}")
  unset(_libxc_version_header)
endif()
