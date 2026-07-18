## Bundle the runtime DLLs next to critic2.exe at install time (Windows only)

## Directories searched for DLLs: the C/Fortran compiler runtime
## directories (works for any gcc version), the sysroot bin/lib dirs,
## and the directory of the glfw library for GUI builds.
set(_dll_dirs)
if (CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  foreach(_dll libgcc_s_seh-1.dll libgfortran-5.dll libwinpthread-1.dll libstdc++-6.dll)
    foreach(_cc ${CMAKE_C_COMPILER} ${CMAKE_Fortran_COMPILER})
      execute_process(COMMAND ${_cc} -print-file-name=${_dll}
        OUTPUT_VARIABLE _path OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      if (_path AND EXISTS "${_path}")
        get_filename_component(_dir "${_path}" DIRECTORY)
        list(APPEND _dll_dirs "${_dir}")
      endif()
    endforeach()
  endforeach()
  ## Native MinGW (e.g. MSYS2): the runtime DLLs (libgfortran, libstdc++,
  ## glfw3, openblas, hdf5, ...) live in the compiler's bin directory, which
  ## -print-file-name above does not report (it lists the lib/ import libs).
  ## Add it explicitly so GET_RUNTIME_DEPENDENCIES can resolve them; without
  ## this, a native build packages critic2.exe with no DLLs.
  get_filename_component(_ccbin "${CMAKE_C_COMPILER}" DIRECTORY)
  if (EXISTS "${_ccbin}")
    list(APPEND _dll_dirs "${_ccbin}")
  endif()
endif()
if (CMAKE_FIND_ROOT_PATH)
  foreach(_root ${CMAKE_FIND_ROOT_PATH})
    list(APPEND _dll_dirs "${_root}/bin" "${_root}/lib")
  endforeach()
endif()
if (ENABLE_GUI AND GLFW3_LIBRARY AND EXISTS "${GLFW3_LIBRARY}")
  get_filename_component(_dir "${GLFW3_LIBRARY}" DIRECTORY)
  list(APPEND _dll_dirs "${_dir}")
endif()
## Also search bin/ and lib/ of every prefix CMake uses to find packages. This
## bundles the DLLs of libraries found outside the compiler's own tree -- e.g.
## GLFW/HDF5/OpenBLAS from an MSYS2 prefix when building with a standalone
## (WinLibs) compiler, where the compiler runtime DLLs and the library DLLs
## live in different directories. Adding search dirs is safe: only DLLs the
## executable actually imports are bundled, and system DLLs are POST-excluded.
foreach(_p ${CMAKE_PREFIX_PATH} ${CMAKE_SYSTEM_PREFIX_PATH})
  foreach(_sub bin lib)
    if (EXISTS "${_p}/${_sub}")
      list(APPEND _dll_dirs "${_p}/${_sub}")
    endif()
  endforeach()
endforeach()
list(REMOVE_DUPLICATES _dll_dirs)

install(CODE "set(_dll_dirs \"${_dll_dirs}\")" COMPONENT runtime)
if (NOT CMAKE_HOST_WIN32)
  ## cross-compiling: tell GET_RUNTIME_DEPENDENCIES to parse PE binaries
  ## with the toolchain objdump instead of the host's default (ELF)
  install(CODE "
    set(CMAKE_GET_RUNTIME_DEPENDENCIES_PLATFORM \"windows+pe\")
    set(CMAKE_GET_RUNTIME_DEPENDENCIES_TOOL \"objdump\")
    set(CMAKE_OBJDUMP \"${CMAKE_OBJDUMP}\")
  " COMPONENT runtime)
endif()
install(CODE [[
  message(STATUS "Resolving runtime DLL dependencies of critic2.exe")
  file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES "$<TARGET_FILE:critic2>"
    RESOLVED_DEPENDENCIES_VAR _deps
    UNRESOLVED_DEPENDENCIES_VAR _undeps
    DIRECTORIES ${_dll_dirs}
    PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-"
    POST_EXCLUDE_REGEXES "[Ss]ystem32" "SysWOW64")
  foreach(_f ${_deps})
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
      TYPE SHARED_LIBRARY FOLLOW_SYMLINK_CHAIN FILES "${_f}")
  endforeach()
  ## Windows system DLLs (kernel32, msvcrt, opengl32, ...) are expected
  ## to be unresolved when cross-compiling; they ship with Windows.
  if (_undeps)
    message(STATUS "DLLs not bundled (expected for system DLLs): ${_undeps}")
  endif()
]] COMPONENT runtime)
