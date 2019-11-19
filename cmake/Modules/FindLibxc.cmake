## Sets:
## LIBXC_FOUND - whether the following three items were found
## LIBXC_xc_LIBRARY - the libxc library.
## LIBXC_xcf90_LIBRARY - the libxcf90 library.
## LIBXC_INCLUDE_DIRS - the directory containing the mod files.
## LIBXC_VERSION - the version of the libxc found

find_library(LIBXC_xc_LIBRARY
  NAMES xc
  )
find_library(LIBXC_xcf90_LIBRARY
  NAMES xcf90
  )

find_path(LIBXC_INCLUDE_DIRS
  NAMES xc_f90_lib_m.mod
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBXC
  REQUIRED_VARS LIBXC_xc_LIBRARY LIBXC_xcf90_LIBRARY LIBXC_INCLUDE_DIRS
  )

mark_as_advanced(LIBXC_xc_LIBRARY LIBXC_xcf90_LIBRARY LIBXC_INCLUDE_DIRS)

## Get libxc version
## From scibuilder, https://github.com/scibuilder/SciBuilder
if(LIBXC_INCLUDE_DIRS)
  file(READ "${LIBXC_INCLUDE_DIRS}/xc_version.h" _libxc_version_header)
  string(REGEX MATCH "define[ \t]+XC_VERSION[ \t]+\\\"([0-9\\.]+)\\\"" LIBXC_VERSION "${_libxc_version_header}")
  string(REGEX MATCH "([0-9\\.]+)" LIBXC_VERSION "${LIBXC_VERSION}")
  unset(_libxc_version_header)
endif()
