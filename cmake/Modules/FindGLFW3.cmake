# Locate the glfw3 library
#
# This module defines the following variables:
#
# GLFW3_LIBRARY the name of the library (or an imported target);
# GLFW3_INCLUDE_DIR where to find glfw include files.
# GLFW3_FOUND true if both the GLFW3_LIBRARY and GLFW3_INCLUDE_DIR have been found.
#
# To help locate the library and include file, you can define a
# variable called GLFW3_ROOT which points to the root of the glfw library
# installation (e.g. an unpacked glfw-*.bin.WIN64 binary distribution).

include(FindPackageHandleStandardArgs)

## prefer the config package shipped by most distributions
## (apt, Homebrew, MSYS2, vcpkg)
find_package(glfw3 CONFIG QUIET)
if (TARGET glfw)
  set(GLFW3_LIBRARY glfw)
  get_target_property(GLFW3_INCLUDE_DIR glfw INTERFACE_INCLUDE_DIRECTORIES)
  if (NOT GLFW3_INCLUDE_DIR)
    set(GLFW3_INCLUDE_DIR "")
  endif()
  find_package_handle_standard_args(GLFW3 DEFAULT_MSG GLFW3_LIBRARY)
  return()
endif()

# Check environment for root search directory
if (NOT GLFW3_ROOT AND DEFINED ENV{GLFW3_ROOT})
  set(GLFW3_ROOT "$ENV{GLFW3_ROOT}")
endif()

# GLFW3_ROOT may point outside the cross-compilation sysroot
set(_root_both)
if (GLFW3_ROOT)
  set(_root_both CMAKE_FIND_ROOT_PATH_BOTH)
endif()

# Search for the header
find_path(GLFW3_INCLUDE_DIR "GLFW/glfw3.h"
  HINTS "${GLFW3_ROOT}"
  PATH_SUFFIXES include
  ${_root_both})

# Search for the library; lib-mingw-w64 and lib-vc* are the layouts of
# the glfw binary distributions for Windows
find_library(GLFW3_LIBRARY NAMES glfw3 glfw
  HINTS "${GLFW3_ROOT}"
  PATH_SUFFIXES lib lib-mingw-w64 lib-vc2022 lib-vc2019
  ${_root_both})
unset(_root_both)

find_package_handle_standard_args(GLFW3 DEFAULT_MSG GLFW3_LIBRARY GLFW3_INCLUDE_DIR)
mark_as_advanced(GLFW3_INCLUDE_DIR GLFW3_LIBRARY)
