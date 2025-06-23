## project and version
cmake_minimum_required(VERSION 3.1...3.27 FATAL_ERROR)

# policies
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

## some project variables
file(READ "VERSION" critic2_VERSION)
string(REGEX REPLACE " *\n$" "" critic2_VERSION "${critic2_VERSION}")
set(critic2_URL "https://aoterodelaroza.github.io/critic2/")
set(critic2_EMAIL "aoterodelaroza@gmail.com")
set(critic2_LICENSE "GNU/GPL version 3")
set(critic2_DESCRIPTION "Analysis of quantum mechanical calculations in molecules and solids")
string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%S")

## build recipes
if (BUILD_RECIPE STREQUAL "MSYS2")
  ## system and compiler details
  set(CMAKE_SYSTEM_NAME Windows)
  set(TOOLCHAIN_PREFIX x86_64-w64-mingw32)

  set(CMAKE_Fortran_COMPILER ${TOOLCHAIN_PREFIX}-gfortran)
  set(CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc)
  set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)

  set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX})
  set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
  set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
  set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

  ## set the cache
  set(CMAKE_INSTALL_PREFIX "D:/msys2/home/media service/critic2-install" CACHE PATH "Install directory" FORCE)
  set(BUILD_TESTING ON CACHE BOOL "" FORCE)
  set(USE_EXTERNAL_LAPACK OFF CACHE BOOL "" FORCE)
  set(ENABLE_GUI ON CACHE BOOL "" FORCE)
  set(BUILD_TYPE "Debug" CACHE STRING "" FORCE)
  set(GLFW3_ROOT "/home/alberto/dll/glfw-3.3.8.bin.WIN64" CACHE STRING "" FORCE)
  set(USE_GUI_THREADS ON CACHE BOOL "" FORCE)
elseif (BUILD_RECIPE STREQUAL "WINE")
  ## system and compiler details
  set(CMAKE_SYSTEM_NAME Windows)
  set(TOOLCHAIN_PREFIX x86_64-w64-mingw32)
  set(CMAKE_CROSSCOMPILING_EMULATOR wine)

  set(CMAKE_Fortran_COMPILER ${TOOLCHAIN_PREFIX}-gfortran)
  set(CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc)
  set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)

  set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX})
  set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
  set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
  set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

  ## set the cache
  set(CMAKE_INSTALL_PREFIX "/home/alberto/critic2-install" CACHE PATH "Install directory" FORCE)
  set(CRITIC2_ROOT "Z:\\home\\alberto\\critic2-install\\share" CACHE STRING "Critic2 ROOT directory" FORCE)
  set(BUILD_TESTING ON CACHE BOOL "" FORCE)
  set(ENABLE_GUI ON CACHE BOOL "" FORCE)
  set(BUILD_TYPE "Debug" CACHE STRING "" FORCE)
  set(GLFW3_ROOT "/home/alberto/critic2-install/glfw-3.3.8.bin.WIN64" CACHE STRING "" FORCE)
  set(USE_GUI_THREADS ON CACHE BOOL "" FORCE)
endif()

## project definition
project(critic2 VERSION "${critic2_VERSION}" LANGUAGES C Fortran)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_Fortran_VISIBILITY_PRESET hidden)

## global flags for cmake
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_STANDARD 99)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(GNUInstallDirs)
include(utils)

## The install prefix in the cache
set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "Installation directory.")

## whether we are following a build recipe
if (NOT BUILD_RECIPE)
  set(BUILD_RECIPE OFF CACHE STRING "Choose the build recipe." FORCE)
else()
  set(BUILD_RECIPE CACHE STRING "Choose the build recipe.")
endif()
set_property(CACHE BUILD_RECIPE PROPERTY STRINGS OFF MSYS2 WINE)

## default build type, and type options
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
else()
  set(CMAKE_BUILD_TYPE CACHE STRING "Choose the type of build.")
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Profile Release)

## interface options
option(BUILD_LIBRARY "Build the critic2 shared library." OFF)
option(ENABLE_GUI "Enables the critic2 GUI" OFF)
option(USE_EXTERNAL_LAPACK "Use external LAPACK/BLAS library (linear algebra operations)." ON)
#option(USE_EXTERNAL_QHULL "Use external QHULL library (calculation of convex hulls)." ON)
#option(USE_EXTERNAL_SPGLIB "Use external spglib library (crystal symmetry determination)." ON)
option(USE_LIBXC "Use the libxc library for access (exchange-correlation energies and potentials)." ON)
option(USE_LIBCINT "Use the libcint library (calculation of molecular integrals with Gaussian functions)." ON)
option(USE_READLINE "Use the readline library (improved command-line interface)." ON)
option(USE_NLOPT "Use the nlopt library (global optimization for variable-cell crystal structure comparison)." ON)
option(ENABLE_OPENMP "Enable OpenMP parallelization." ON)
option(BUILD_STATIC "Build critic2 statically (tested only on linux, and no GUI)." OFF)

## testing
option(BUILD_TESTING "Enable the regression tests for the critic2 build." OFF)
include(CTest)

## library requires fpic
if (BUILD_LIBRARY)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

## GUI requirements
if (ENABLE_GUI)
  ## GUI options
  option(USE_HDF5 "Use the HDF5 library (reading hdf5 phonopy vibration data)." ON)
  option(USE_FREETYPE "Use the freetype library (improved font rendering in GUI)." ON)
  option(USE_GUI_THREADS "Use threads in the GUI." ON)

  ## needs -fPIC
  enable_language(CXX)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  ## find GUI libraries
  if (USE_FREETYPE)
    find_package(Freetype)
  endif()

  ## HDF5
  if (USE_HDF5)
    find_package(HDF5 COMPONENTS Fortran HL)
    mark_as_advanced(HDF5_C_COMPILER_EXECUTABLE)
    mark_as_advanced(HDF5_Fortran_COMPILER_EXECUTABLE)
    mark_as_advanced(HDF5_Fortran_HL_stub)
    mark_as_advanced(HDF5_Fortran_stub)
    mark_as_advanced(SZIP_INCLUDE_DIR)
    mark_as_advanced(SZIP_LIBRARY)
    mark_as_advanced(h5_conf)
  endif()

  ## glfw3
  find_package(GLFW3 REQUIRED)

  ## opengl
  set(OpenGL_GL_PREFERENCE "GLVND")
  find_package(OpenGL REQUIRED)

  ## the GUI requires threads
  set(ENABLE_OPENMP ON)

  ## find pthreads
  if (USE_GUI_THREADS)
    find_package(Threads REQUIRED)
    if(THREADS_HAVE_PTHREAD_ARG)
      set_property(TARGET critic2 PROPERTY COMPILE_OPTIONS "-pthread")
      set_property(TARGET critic2 PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
    endif()
  endif()
endif()

## perform a static build
if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BLA_STATIC ON)
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  find_library(LIBDL_LIBRARY NAMES dl)
  if(NOT LIBDL_LIBRARY)
    error("Dynamic linking static library not found")
  endif()
  mark_as_advanced(LIBDL_LIBRARY)
  set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} ${LIBDL_LIBRARY}")
endif()

## compiler flags, fortran, and check versions
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(FFLAGS_DEBUG "-g -O0 -C -traceback -debug all -fpp -check all")
  set(FFLAGS_PROFILE "-g -pg -O0 -C -traceback -debug all -fpp -check all")
  set(FFLAGS_RELEASE "-O3")
  if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_EQUAL 20.2.1.20201112)
    set(FFLAGS_DEBUG "-g -O0 -traceback -debug all -fpp")
    set(FFLAGS_PROFILE "-g -pg -O0 -traceback -debug all -fpp")
  elseif ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 19.0)
    set(_message "\nA version of intel fortran >=19.5 is required to compile the development version of critic2. "
      "Please, visit: https://aoterodelaroza.github.io/critic2/installation/#whichcompilerswork"
      "for more information.\n")
    message(FATAL_ERROR ${_message})
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (ENABLE_GUI)
    # -ffpe-trap=invalid makes GLFW crash in some versions
    set(FFLAGS_DEBUG "-g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=zero,overflow -fbacktrace -fdump-core -cpp -ffree-line-length-0")
    set(FFLAGS_PROFILE "-g -pg -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=zero,overflow -fbacktrace -fdump-core -cpp -ffree-line-length-0")
  else()
    set(FFLAGS_DEBUG "-g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp -ffree-line-length-0")
    set(FFLAGS_PROFILE "-g -pg -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp -ffree-line-length-0")
  endif()
  set(FFLAGS_RELEASE "-O3 -ffree-line-length-0")
  if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 6.0.0)
    set(_message "\nA version of gfortran >=6.0 is required to compile the development version of critic2. "
      "Please, visit: https://aoterodelaroza.github.io/critic2/installation/#whichcompilerswork"
      "for more information.\n")
    message(FATAL_ERROR ${_message})
  endif()
endif()

## compiler flags, C
if("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
  set(CFLAGS_DEBUG "-D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -ggdb -Wall")
  set(CFLAGS_PROFILE "-D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -ggdb -Wall -pg")
  set(CFLAGS_RELEASE "-O3")
endif()

## compiler flags, CXX
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  set(CXXFLAGS_DEBUG "-D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -ggdb -Wall")
  set(CXXFLAGS_PROFILE "-D_FORTIFY_SOURCE=2 -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -ggdb -Wall -pg")
  set(CXXFLAGS_RELEASE "-O3")
endif()

## set the compiler flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BT)
set(CMAKE_Fortran_FLAGS_${BT} "${FFLAGS_${BT}}")
set(CMAKE_C_FLAGS_${BT} "${CFLAGS_${BT}}")
set(CMAKE_CXX_FLAGS_${BT} "${CXXFLAGS_${BT}}")

## build recipes flags
if (BUILD_RECIPE STREQUAL "MSYS2" OR BUILD_RECIPE STREQUAL "WINE")
  set(CMAKE_Fortran_FLAGS_${BT} "${CMAKE_Fortran_FLAGS_${BT}} -mwindows")
  set(CMAKE_C_FLAGS_${BT} "${CMAKE_C_FLAGS_${BT}} -mwindows")
  set(CMAKE_CXX_FLAGS_${BT} "${CMAKE_CXX_FLAGS_${BT}} -mwindows")
endif()

## openmp
if (ENABLE_OPENMP)
  find_package(OpenMP MODULE COMPONENTS Fortran)
  if (OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS_${BT} "${CMAKE_Fortran_FLAGS_${BT}} ${OpenMP_Fortran_FLAGS}")
  endif()
endif()

## lapack
if (USE_EXTERNAL_LAPACK)
  set(LAPACK_FIND_QUIETLY "ON")
  find_package(LAPACK)
endif()

# ## qhull -- deactivated
# if (USE_EXTERNAL_QHULL)
#   find_package(QHULL)
# endif()

# ## spglib -- deactivated because the interface keeps changing!
# if (USE_EXTERNAL_SPGLIB)
#   find_package(SPGLIB)
# endif()

## numdiff
find_package(NUMDIFF)

## libxc
if (USE_LIBXC)
  find_package(LIBXC)
  if (LIBXC_FOUND)
    set(MINLIBXC_VERSION 5.0)
    if (${LIBXC_VERSION} VERSION_LESS ${MINLIBXC_VERSION})
      set(LIBXC_FOUND FALSE)
      message("!! Libxc version found (${LIBXC_VERSION},${LIBXC_xc_LIBRARY}) too old (requires ${MINLIBXC_VERSION})")
    endif()
    unset(MINLIBXC_VERSION)
  endif()
endif()

## readline
if (USE_READLINE)
  find_package(READLINE)
endif()

## libcint
if (USE_LIBCINT)
  find_package(LIBCINT)
endif()

## nlopt
if (USE_NLOPT)
  find_package(NLOPT)
endif()

## the syntax.txt target
add_custom_target(syntax
  COMMAND wget -qO - 'raw.githubusercontent.com/aoterodelaroza/aoterodelaroza.github.io/master/_critic2/13_syntax.md' | awk '/~~~/{a=!a\;next}a{print}' > syntax.txt
  )

## process subdirectories
add_subdirectory(src)
add_subdirectory(dat)
if (BUILD_TESTING AND NUMDIFF_FOUND)
  add_subdirectory(tests)
endif()

## summary
set(_message
"\n\n### ${CMAKE_PROJECT_NAME}, ${CMAKE_PROJECT_VERSION}###\n"
"Date: ${DATE}\n"
"Build type: ${CMAKE_BUILD_TYPE}\n"
"Fortran compiler: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID}, ${CMAKE_Fortran_COMPILER_VERSION})\n"
"Fortran flags: ${CMAKE_Fortran_FLAGS_${BT}}\n"
"C compiler: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID}, ${CMAKE_C_COMPILER_VERSION})\n"
"C flags: ${CMAKE_C_FLAGS_${BT}}\n"
)

if (ENABLE_GUI)
  set(_message
    "${_message}CXX compiler: ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID}, ${CMAKE_CXX_COMPILER_VERSION})\n"
    "CXX flags: ${CMAKE_CXX_FLAGS_${BT}}\n"
  )
endif()

set(_message
"${_message}Install directory: ${CMAKE_INSTALL_PREFIX}\n"
"OpenMP parallelization: ${OPENMP_FOUND}\n"
)

if (LAPACK_FOUND)
  set(_message "${_message}Lapack: lib=${LAPACK_LIBRARIES}\n")
else()
  set(_message "${_message}Lapack: <internal>\n")
endif()
# if (QHULL_FOUND)
#   set(_message "${_message}Qhull: lib=${QHULL_LIBRARIES} | inc=${QHULL_INCLUDE_DIRS}\n")
# else()
#   set(_message "${_message}Qhull: <internal>\n")
# endif()
# if (SPGLIB_FOUND)
#   set(_message "${_message}Spglib: lib=${SPGLIB_LIBRARY}\n")
# else()
#   set(_message "${_message}Spglib: <internal>\n")
# endif()
if (LIBXC_FOUND)
  set(_message "${_message}Libxc: ${LIBXC_VERSION} | lib=${LIBXC_xcf03_LIBRARY} ${LIBXC_xc_LIBRARY} | inc=${LIBXC_INCLUDE_DIRS}\n")
else()
  set(_message "${_message}Libxc: not used\n")
endif()
if (LIBCINT_FOUND)
  set(_message "${_message}Libcint: ${LIBCINT_VERSION} | lib=${LIBCINT_LIBRARY}\n")
else()
  set(_message "${_message}Libcint: not used\n")
endif()
if (READLINE_FOUND)
  set(_message "${_message}Readline: lib=${READLINE_LIBRARY} | inc=${READLINE_INCLUDE_DIRS}\n")
else()
  set(_message "${_message}Readline: not used\n")
endif()
if (NLOPT_FOUND)
  set(_message "${_message}Nlopt: ${NLOPT_VERSION} | lib=${NLOPT_LIBRARIES}\n")
else()
  set(_message "${_message}Nlopt: not used\n")
endif()
if (NUMDIFF_FOUND)
  set(_message "${_message}Numdiff: ${NUMDIFF_EXE}\n")
else()
  set(_message "${_message}Numdiff: not found\n")
endif()

set(_message "${_message}GUI build: ${ENABLE_GUI}\n")

if (ENABLE_GUI)
  set(_message
    "${_message}gl: lib=${OPENGL_LIBRARIES} | inc=${OPENGL_INCLUDE_DIRS}\n"
    "glm: lib=${GLM_LIBRARIES} | inc=${GLM_INCLUDE_DIRS}\n"
    "glfw: lib=${GLFW3_LIBRARY} | inc=${GLFW3_INCLUDE_DIR}\n"
    "freetype: lib=${FREETYPE_LIBRARIES} | inc=${FREETYPE_INCLUDE_DIRS}\n"
    "hdf5: lib=${HDF5_LIBRARIES} | inc=${HDF5_INCLUDE_DIRS}\n"
    "dl: lib=${CMAKE_DL_LIBS}\n")
endif()

set(_message "${_message}################################\n")

message(STATUS ${_message})
