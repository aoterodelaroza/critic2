## Find the libxc library. Sets:
## LIBXC_FOUND - whether the following three items were found
## LIBXC_xc_LIBRARY - the libxc library.
## LIBXC_xcf03_LIBRARY - the libxcf03 library.
## LIBXC_INCLUDE_DIRS - the directory containing the mod files.
## LIBXC_VERSION - the version of the libxc found

include(FindPackageHandleStandardArgs)

if (DEFINED ENV{LIBXC_DIR})
  set(LIBXC_DIR "$ENV{LIBXC_DIR}")
endif()

## LIBXC_DIR may point outside the cross-compilation sysroot
set(_root_both)
if (LIBXC_DIR)
  set(_root_both CMAKE_FIND_ROOT_PATH_BOTH)
endif()

find_library(LIBXC_xc_LIBRARY
  NAMES xc
  PATH_SUFFIXES lib
  HINTS ${LIBXC_DIR}
  ${_root_both}
  )
find_library(LIBXC_xcf03_LIBRARY
  NAMES xcf03
  PATH_SUFFIXES lib
  HINTS ${LIBXC_DIR}
  ${_root_both}
  )

find_path(LIBXC_INCLUDE_DIRS
  NAMES xc_f03_lib_m.mod
  PATH_SUFFIXES include
  HINTS ${LIBXC_DIR}
  ${_root_both}
  )
unset(_root_both)

find_package_handle_standard_args(LIBXC
  REQUIRED_VARS LIBXC_xc_LIBRARY LIBXC_xcf03_LIBRARY LIBXC_INCLUDE_DIRS
  )

mark_as_advanced(LIBXC_xc_LIBRARY LIBXC_xcf03_LIBRARY LIBXC_INCLUDE_DIRS)
if (NOT LIBXC_FOUND)
  set(LIBXC_DIR "${LIBXC_DIR}" CACHE STRING "Directory containing the libxc library (>=4.1).")
endif()

## Check whether we can actually compile against it. libxc's Fortran .mod
## files are compiler-specific, so a libxc built with a different compiler
## (e.g. the system gfortran build used with ifort/ifx) is unusable and
## must disable libxc. Use a separate result variable: the normal
## LIBXC_FOUND set by find_package_handle_standard_args above would
## otherwise shadow the cache value written by try_compile, so the check
## result would be silently ignored.
if (LIBXC_FOUND)
  try_compile(LIBXC_COMPILES "${CMAKE_BINARY_DIR}/temp" "${CMAKE_SOURCE_DIR}/cmake/Modules/libxc_test.f90"
    LINK_LIBRARIES ${LIBXC_xcf03_LIBRARY} ${LIBXC_xc_LIBRARY}
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${LIBXC_INCLUDE_DIRS}")
  if (NOT LIBXC_COMPILES)
    message(STATUS "Found libxc (lib=${LIBXC_xcf03_LIBRARY} ${LIBXC_xc_LIBRARY} | inc=${LIBXC_INCLUDE_DIRS}) but could not compile against it (different compiler?)")
    set(LIBXC_FOUND FALSE)
  endif()
endif()

## Get libxc version, from scibuilder, https://github.com/scibuilder/SciBuilder
if(LIBXC_INCLUDE_DIRS)
  file(READ "${LIBXC_INCLUDE_DIRS}/xc_version.h" _libxc_version_header)
  string(REGEX MATCH "define[ \t]+XC_VERSION[ \t]+\\\"([0-9\\.]+)\\\"" LIBXC_VERSION "${_libxc_version_header}")
  string(REGEX MATCH "([0-9\\.]+)" LIBXC_VERSION "${LIBXC_VERSION}")
  unset(_libxc_version_header)
endif()
