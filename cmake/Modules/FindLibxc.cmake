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
