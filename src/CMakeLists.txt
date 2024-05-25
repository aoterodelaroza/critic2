## sources and critic2 core
set(SOURCES_Fortran abinit_private.f90 abinit_private@proc.f90 autocp.f90
  autocp@proc.f90 bader.f90 bader@proc.f90 bisect.f90 bisect@proc.f90
  cfftnd.f90 c_interface_module.f90 config.f90
  crystalmod.f90 crystalmod@complex.f90 crystalmod@edit.f90 crystalmod@env.f90
  crystalmod@mols.f90
  crystalmod@proc.f90 crystalmod@powderproc.F90 crystalmod@symmetry.f90
  crystalmod@write.f90 crystalseedmod.f90
  crystalseedmod@proc.f90 dftb_private.f90
  dftb_private@proc.f90 elk_private.f90 elk_private@proc.f90 fieldmod.f90
  fieldmod@proc.f90 fieldseedmod.f90 fieldseedmod@proc.f90
  flux.f90 flux@proc.f90 fragmentmod.f90
  fragmentmod@proc.f90 global.f90 global@proc.F90
  graphics.f90 graphics@proc.f90 grid1mod.f90
  grid1mod@proc.f90 grid3mod.f90 grid3mod@proc.f90
  hashmod.f90 hashmod@proc.f90
  hirshfeld.f90 hirshfeld@proc.f90 integration.f90
  integration@proc.f90 keast.f90 meshmod.f90
  meshmod@proc.f90 molcalc.f90 nci.f90 nci@proc.f90
  pi_private.f90 pi_private@proc.f90 qe_private.f90
  qe_private@proc.f90 qtree_basic.f90
  qtree_basic@proc.f90 qtree.f90 qtree_gpaths.f90
  qtree_gpaths@proc.f90 qtree@proc.f90 qtree_tetrawork.f90
  qtree_tetrawork@proc.f90 qtree_utils.f90
  qtree_utils@proc.f90 quadpack.f90 rhoplot.f90
  rhoplot@proc.f90 sigmahole.f90 sigmahole@proc.f90
  spglib.f90 spglib@proc.f90 stm.f90
  stm@proc.f90 struct_drivers.f90 struct_drivers@proc.f90
  surface.f90 surface@proc.f90 systemmod.f90
  systemmod@proc.f90 tools.f90 tools_io.f90 tools_math.f90
  tools_math@lebedev.f90 tools_math@proc.f90 tools_math@sobol.f90
  tools@proc.f90
  tricks.F90 types.f90 types@proc.f90 wfn_private.f90
  wien_private.f90 wien_private@proc.f90 xdm.f90
  xdm@proc.f90 yt.f90 yt@proc.f90
  arithmetic.F90 arithmetic@proc.F90 config@proc.F90
  molcalc@proc.F90 param.F90 param@proc.f90
  tools_io@proc.F90 wfn_private@proc.F90
  )
set(SOURCES_C doqhull.c)
set(SOURCES_Fortran_Main critic2.F90)
set(SOURCES_Fortran_Library libcritic2.f90 libcritic2@proc.f90)
set(SOURCES ${SOURCES_Fortran} ${SOURCES_C} ${SOURCES_Fortran_Main})
add_library(critic2_core OBJECT ${SOURCES_Fortran} ${SOURCES_C})

## gfortran-10: do not check fftpack for argument mismatch
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" AND ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER 10.0.0
      OR "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_EQUAL 10.0.0))
  set_source_files_properties(cfftnd.f90 PROPERTIES COMPILE_FLAGS "-fallow-argument-mismatch")
endif()

## C standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

## header info
set(PACKAGE ${CMAKE_PROJECT_NAME})
set(VERSION "${CMAKE_PROJECT_VERSION}")
set(TARGET ${CMAKE_SYSTEM})
if (CRITIC2_ROOT)
  set_property(SOURCE config@proc.F90 APPEND PROPERTY COMPILE_FLAGS "-DDATADIR='\"${CRITIC2_ROOT}\"'")
else()
  set_property(SOURCE config@proc.F90 APPEND PROPERTY COMPILE_FLAGS "-DDATADIR='\"${CMAKE_INSTALL_FULL_DATADIR}\"'")
endif()
configure_file(config@proc.F90.in "${CMAKE_CURRENT_SOURCE_DIR}/config@proc.F90")

## executable
add_executable(critic2 ${SOURCES})
add_dependencies(critic2 critic2_core)

## library
if (BUILD_LIBRARY)
  add_library(libcritic2 SHARED ${SOURCES_Fortran} ${SOURCES_C} ${SOURCES_Fortran_Library})
  set_target_properties(libcritic2 PROPERTIES OUTPUT_NAME critic2 PUBLIC_HEADER "libcritic2.h")
  add_dependencies(libcritic2 critic2_core critic2)
endif()

## cubpack
add_subdirectory(cubpack)
include_directories("${CMAKE_BINARY_DIR}/src/cubpack")
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:cubpack>)
if (BUILD_LIBRARY)
  target_sources(libcritic2 PRIVATE $<TARGET_OBJECTS:cubpack>)
endif()
add_dependencies(critic2_core cubpack)

## openmp
if (OPENMP_FOUND)
  target_link_libraries(critic2 PRIVATE ${OpenMP_Fortran_LIBRARIES})
  if (BUILD_LIBRARY)
    target_link_libraries(libcritic2 PRIVATE ${OpenMP_Fortran_LIBRARIES})
  endif()
endif()

## lapack
if (LAPACK_FOUND)
  set_property(TARGET critic2 APPEND PROPERTY LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
  target_link_libraries(critic2 PRIVATE ${LAPACK_LIBRARIES})
  if (BUILD_LIBRARY)
    set_property(TARGET libcritic2 APPEND PROPERTY LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
    target_link_libraries(libcritic2 PRIVATE ${LAPACK_LIBRARIES})
  endif()
  if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_${BT} "${CMAKE_Fortran_FLAGS_${BT}} -fexternal-blas")
  endif()
else()
  add_subdirectory(oldlibs)
  target_sources(critic2 PRIVATE $<TARGET_OBJECTS:oldlibs>)
  if (BUILD_LIBRARY)
    target_sources(libcritic2 PRIVATE $<TARGET_OBJECTS:oldlibs>)
  endif()
  add_dependencies(critic2_core oldlibs)
endif()

## qhull
add_subdirectory(qhull)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/qhull")
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:qhull>)
if (BUILD_LIBRARY)
  target_sources(libcritic2 PRIVATE $<TARGET_OBJECTS:qhull>)
endif()
add_dependencies(critic2_core qhull)

## spglib
add_subdirectory(spglib)
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:spglib>)
if (BUILD_LIBRARY)
  target_sources(libcritic2 PRIVATE $<TARGET_OBJECTS:spglib>)
endif()
add_dependencies(critic2_core spglib)

## json-fortran
add_subdirectory(json-fortran)
include_directories("${CMAKE_BINARY_DIR}/src/json-fortran")
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:json-fortran>)
if (BUILD_LIBRARY)
  target_sources(libcritic2 PRIVATE $<TARGET_OBJECTS:json-fortran>)
endif()
add_dependencies(critic2_core json-fortran)

## libxc
if (USE_LIBXC AND LIBXC_FOUND)
  include_directories("${LIBXC_INCLUDE_DIRS}")
  target_link_libraries(critic2 PRIVATE ${LIBXC_xcf90_LIBRARY} ${LIBXC_xc_LIBRARY})
  if (BUILD_LIBRARY)
    target_link_libraries(libcritic2 PRIVATE ${LIBXC_xcf90_LIBRARY} ${LIBXC_xc_LIBRARY})
  endif()
  add_definitions("-DHAVE_LIBXC")
endif()

## libcint
if (USE_LIBCINT AND LIBCINT_FOUND)
  target_link_libraries(critic2 PRIVATE ${LIBCINT_LIBRARY})
  if (BUILD_LIBRARY)
    target_link_libraries(libcritic2 PRIVATE ${LIBCINT_LIBRARY})
  endif()
  add_definitions("-DHAVE_CINT")
endif()

## readline
if (USE_READLINE AND READLINE_FOUND)
  target_link_libraries(critic2 PRIVATE ${READLINE_LIBRARY})
  if (BUILD_LIBRARY)
    target_link_libraries(libcritic2 PRIVATE ${READLINE_LIBRARY})
  endif()
  add_definitions("-DHAVE_READLINE")
endif()

## libnlopt
if (USE_NLOPT AND NLOPT_FOUND)
  target_link_libraries(critic2 PRIVATE ${NLOPT_LIBRARIES})
  if (BUILD_LIBRARY)
    target_link_libraries(libcritic2 PRIVATE ${NLOPT_LIBRARIES})
  endif()
  add_definitions("-DHAVE_NLOPT")
  include_directories("${NLOPT_INCLUDE_DIRS}")
endif()

## emd (earth mover's distance)
add_subdirectory(emd)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/emd")
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:emd>)
if (BUILD_LIBRARY)
  target_sources(libcritic2 PRIVATE $<TARGET_OBJECTS:emd>)
endif()
add_dependencies(critic2_core emd)

## graphical user interface
if (ENABLE_GUI)
  ## GLFW
  target_include_directories(critic2 PRIVATE ${GLFW3_INCLUDE_DIR})
  target_link_libraries(critic2 PRIVATE ${GLFW3_LIBRARY})

  ## FREETYPE
  if (FREETYPE_FOUND)
   target_include_directories(critic2 PRIVATE ${FREETYPE_INCLUDE_DIRS})
   target_link_libraries(critic2 PRIVATE ${FREETYPE_LIBRARIES})
  endif()

  ## OpenGL
  target_link_libraries(critic2 PRIVATE OpenGL::GL)

  ## the dynamic linking library
  target_link_libraries(critic2 PRIVATE ${CMAKE_DL_LIBS})

  ## the threads library
  if (USE_GUI_THREADS AND CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(critic2 PRIVATE ${CMAKE_THREAD_LIBS_INIT})
  endif()

  ## process subdirectory
  add_subdirectory(gui)
  target_sources(critic2 PRIVATE $<TARGET_OBJECTS:fgui>)
  include_directories("${CMAKE_BINARY_DIR}/src/gui")

  ## add dependency to the core critic2 files
  add_dependencies(fgui critic2_core)

  ## add the HAVE_GUI definition
  add_definitions("-DHAVE_GUI")
endif()

## cross-compilation for windows, add the flag
if (WIN32)
  add_definitions("-D_WIN32")
endif()

## install the binary
install(TARGETS critic2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if (BUILD_LIBRARY)
  install(TARGETS libcritic2
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

## recipe specific copies and install
if (BUILD_RECIPE STREQUAL "MSYS2")
  foreach (i IN ITEMS libgcc_s_seh-1.dll glfw3.dll libgomp-1.dll libreadline8.dll libstdc++-6.dll libgfortran-5.dll libwinpthread-1.dll libtermcap-0.dll libquadmath-0.dll)
    execute_process (
        COMMAND bash -c "${CMAKE_C_COMPILER} -print-file-name=${i}"
        OUTPUT_VARIABLE TMP
    )
    add_custom_command(
      TARGET critic2
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${TMP} ${CMAKE_BINARY_DIR}/src)
    install(FILES ${CMAKE_BINARY_DIR}/src/${i} DESTINATION ${CMAKE_INSTALL_BINDIR})
  endforeach()
elseif (BUILD_RECIPE STREQUAL "WINE")
  foreach (i IN ITEMS libgfortran-5.dll libgcc_s_seh-1.dll libgomp-1.dll libstdc++-6.dll libquadmath-0.dll)
    add_custom_command(TARGET critic2 POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy /usr/lib/gcc/x86_64-w64-mingw32/12-win32/${i} ${CMAKE_BINARY_DIR}/src)
    install(FILES ${CMAKE_BINARY_DIR}/src/${i} DESTINATION ${CMAKE_INSTALL_BINDIR})
  endforeach()
  add_custom_command(TARGET critic2 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy /usr/x86_64-w64-mingw32/lib/libwinpthread-1.dll ${CMAKE_BINARY_DIR}/src)
  install(FILES ${CMAKE_BINARY_DIR}/src/libwinpthread-1.dll DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
