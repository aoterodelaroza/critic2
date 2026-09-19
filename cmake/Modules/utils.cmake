## List all variables
## sakra at https://stackoverflow.com/questions/9298278/cmake-print-out-all-accessible-variables-in-a-script
macro(list_variables)
  get_cmake_property(_variableNames VARIABLES)
  list (SORT _variableNames)
  foreach (_variableName ${_variableNames})
      message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endmacro()

## Keep a bundled/vendored numerical library optimized even in Debug
## builds by compiling it with -O2 in every configuration. lang is the
## language of the target's sources (C, Fortran, ...).
function(c2_optimize_vendored tgt lang)
  if (CMAKE_${lang}_COMPILER_ID MATCHES "GNU|Clang|Intel")
    target_compile_options(${tgt} PRIVATE -O2)
  endif()
endfunction()

## Version of a library read from the pkgconfig file that sits beside it: lib
## is the full path to the library found, pc the base name of its .pc file and
## pc_libdir the library directory the host pkg-config reported (PC_<PKG>_LIBDIR).
## The host's pkg-config cannot be trusted on its own because it describes its
## own copy of the library, which in a cross build -- or whenever <PKG>_DIR
## points somewhere else -- is not the copy found. Sets outvar to the version
## read beside the library, and leaves it undefined when the library found is
## not the one pkg-config described and carries no .pc of its own, rather than
## reporting a version that belongs to a different library.
function(c2_pkgconfig_version outvar lib pc pc_libdir)
  get_filename_component(_libdir "${lib}" DIRECTORY)
  if (EXISTS "${_libdir}/pkgconfig/${pc}.pc")
    file(STRINGS "${_libdir}/pkgconfig/${pc}.pc" _ver REGEX "^Version:")
    if (_ver)
      string(REGEX REPLACE "^Version:[ \t]*" "" _ver "${_ver}")
      set(${outvar} "${_ver}" PARENT_SCOPE)
    endif()
  elseif (NOT "${_libdir}" STREQUAL "${pc_libdir}")
    unset(${outvar} PARENT_SCOPE)
  endif()
endfunction()

## Whether a library is already available on the system, used only to choose
## the default value of the corresponding USE_* option before the real
## detection runs (e.g. USE_XTB is on by default when Debian's libxtb-dev is
## installed). pkg is the pkg-config module name and lib the library name;
## <PKG>_DIR, from the environment or the command line, is honored as it is by
## the Find modules, and so is BUILD_STATIC. The pkg-config check uses the
## PC_<PKG> prefix of the corresponding Find module, so the find_package
## further down reuses this result instead of running pkg-config all over again.
function(c2_library_available outvar pkg lib)
  string(TOUPPER "${pkg}" _up)
  set(_avail OFF)

  ## pkg-config is no help in either of two cases: in a cross build it
  ## describes the host's library and not the target's, and under BUILD_STATIC
  ## it cannot say whether a static library exists at all. Both go by
  ## find_library, which honors CMAKE_FIND_ROOT_PATH and the suffix list below.
  if (NOT CMAKE_CROSSCOMPILING AND NOT BUILD_STATIC)
    find_package(PkgConfig QUIET)
    pkg_check_modules(PC_${_up} ${pkg} QUIET)
    if (PC_${_up}_FOUND)
      set(_avail ON)
    endif()
  endif()

  ## the BUILD_STATIC block further down restricts the library suffixes; match
  ## it here (in function scope) so the default agrees with what the Find
  ## module will be able to use
  if (BUILD_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  endif()

  if (NOT _avail)
    ## same <PKG>_DIR precedence as the Find modules: the environment wins over
    ## the cache/command-line variable, and either may point outside the
    ## cross-compilation sysroot
    set(_dir "${${_up}_DIR}")
    if (DEFINED ENV{${_up}_DIR})
      set(_dir "$ENV{${_up}_DIR}")
    endif()
    set(_root_both)
    if (_dir)
      set(_root_both CMAKE_FIND_ROOT_PATH_BOTH)
    endif()
    find_library(_c2_probe_lib
      NAMES ${lib}
      HINTS "${_dir}"
      PATH_SUFFIXES lib lib64
      PATHS "${CMAKE_INSTALL_PREFIX}/lib" "${CMAKE_INSTALL_PREFIX}/lib64"
      ${_root_both})
    if (_c2_probe_lib)
      set(_avail ON)
    endif()
    ## the probe's own cache entry is dropped: the Find module does the real
    ## search, with better hints, and should not be short-circuited by this
    unset(_c2_probe_lib CACHE)
  endif()

  set(${outvar} ${_avail} PARENT_SCOPE)
endfunction()
