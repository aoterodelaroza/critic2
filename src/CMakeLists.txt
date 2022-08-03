## sources
set(SOURCES_Fortran abinit_private.f90 abinit_private@proc.f90 autocp.f90
  autocp@proc.f90 bader.f90 bader@proc.f90 bisect.f90 bisect@proc.f90
  cfftnd.f90 c_interface_module.f90 config.f90
  crystalmod.f90 crystalmod@proc.f90 crystalseedmod.f90
  crystalseedmod@proc.f90 dftb_private.f90
  dftb_private@proc.f90 elk_private.f90 elk_private@proc.f90
  environmod.f90 environmod@proc.f90 fieldmod.f90
  fieldmod@proc.f90 fieldseedmod.f90 fieldseedmod@proc.f90
  flux.f90 flux@proc.f90 fragmentmod.f90
  fragmentmod@proc.f90 global.f90 global@proc.F90
  graphics.f90 graphics@proc.f90 grid1mod.f90
  grid1mod@proc.f90 grid3mod.f90 grid3mod@proc.f90
  gui_interface.f90 hashmod.f90 hashmod@proc.f90
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
  tricks.f90 types.f90 types@proc.f90 wfn_private.f90
  wien_private.f90 wien_private@proc.f90 xdm.f90
  xdm@proc.f90 yt.f90 yt@proc.f90
  arithmetic.F90 arithmetic@proc.F90 config@proc.F90
  molcalc@proc.F90 param.F90 param@proc.f90
  tools_io@proc.F90 wfn_private@proc.F90
  )
set(SOURCES_C doqhull.c)
set(SOURCES_Fortran_Main critic2.F90)
set(SOURCES ${SOURCES_Fortran} ${SOURCES_C} ${SOURCES_Fortran_Main})
add_library(critic2_core OBJECT ${SOURCES_Fortran} ${SOURCES_C})

## gfortran-10: do not check fftpack for argument mismatch
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" AND "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER 10.0.0)
  set_source_files_properties(cfftnd.f90 PROPERTIES COMPILE_FLAGS -fallow-argument-mismatch)
endif()

## C standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

## header info
set(PACKAGE ${CMAKE_PROJECT_NAME})
set(VERSION "${CMAKE_PROJECT_VERSION}")
set(TARGET ${CMAKE_SYSTEM})
set_property(SOURCE config@proc.F90 APPEND PROPERTY COMPILE_FLAGS "-DDATADIR='\"${CMAKE_INSTALL_FULL_DATADIR}\"'")
configure_file(config@proc.F90.in "${CMAKE_CURRENT_SOURCE_DIR}/config@proc.F90")

## executable
add_executable(critic2 ${SOURCES})

## cubpack
add_subdirectory(cubpack)
include_directories("${CMAKE_BINARY_DIR}/src/cubpack")
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:cubpack>)

## ciftbx
add_subdirectory(ciftbx)
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:ciftbx>)

## openmp
if (OPENMP_FOUND)
  target_link_libraries(critic2 PRIVATE ${OpenMP_Fortran_LIBRARIES})
endif()

## lapack
if (LAPACK_FOUND)
  set_property(TARGET critic2 APPEND PROPERTY LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
  target_link_libraries(critic2 PRIVATE ${LAPACK_LIBRARIES})
else()
  add_subdirectory(oldlibs)
  target_sources(critic2 PRIVATE $<TARGET_OBJECTS:oldlibs>)
endif()

## qhull
add_subdirectory(qhull)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/qhull")
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:qhull>)

## spglib
add_subdirectory(spglib)
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:spglib>)

## json-fortran
add_subdirectory(json-fortran)
include_directories("${CMAKE_BINARY_DIR}/src/json-fortran")
target_sources(critic2 PRIVATE $<TARGET_OBJECTS:json-fortran>)

## libxc
if (LIBXC_FOUND)
  include_directories("${LIBXC_INCLUDE_DIRS}")
  target_link_libraries(critic2 PRIVATE ${LIBXC_xcf90_LIBRARY} ${LIBXC_xc_LIBRARY})
  add_definitions("-DHAVE_LIBXC")
endif()

## libcint
if (LIBCINT_FOUND)
  target_link_libraries(critic2 PRIVATE ${LIBCINT_LIBRARY})
  add_definitions("-DHAVE_CINT")
endif()

## readline
if (READLINE_FOUND)
  target_link_libraries(critic2 PRIVATE ${READLINE_LIBRARY})
  add_definitions("-DHAVE_READLINE")
endif()

## The library for the GUI, fortran version
if (ENABLE_GUI)
  ## GUI libraries
  target_include_directories(critic2 PRIVATE
    ${GLFW_INCLUDE_DIRS}
    ${FREETYPE2_INCLUDE_DIRS}
    ${GLM_INCLUDE_DIRS}
    ${GL_INCLUDE_DIRS}
    )
  target_link_libraries(critic2 PRIVATE
    ${GLFW_LIBRARIES}
    ${FREETYPE2_LIBRARIES}
    ${GLM_LIBRARIES}
    ${GL_LIBRARIES}
    )
  ## the dynamic linking library
  target_link_libraries(critic2 PRIVATE ${CMAKE_DL_LIBS})

  add_subdirectory(gui)
  target_sources(critic2 PRIVATE $<TARGET_OBJECTS:fgui>)
  include_directories("${CMAKE_BINARY_DIR}/src/gui")
  add_dependencies(fgui critic2_core)
endif()

## install the binary
install(TARGETS critic2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
