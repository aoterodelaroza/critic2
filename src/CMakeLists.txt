## sources and critic2 core
set(SOURCES_Fortran abinit_private.f90 abinit_private@proc.f90 autocp.f90
  autocp@proc.f90 bader.f90 bader@proc.f90 bisect.f90 bisect@proc.f90
  cfftnd.f90 c_interface_module.f90 config.f90
  crystalmod.f90 crystalmod@complex.f90 crystalmod@edit.f90 crystalmod@env.f90
  crystalmod@guessspg.f90 crystalmod@vibrations.F90 crystalmod@mols.f90
  crystalmod@proc.f90 crystalmod@powderproc.F90 crystalmod@symmetry.f90
  crystalmod@write.f90 crystalseedmod.f90
  crystalseedmod@proc.f90 dynamics.f90 dynamics@proc.f90 dftb_private.f90
  dftb_private@proc.f90 elk_private.f90 elk_private@proc.f90
  energy.f90 energy@proc.F90 fieldmod.f90
  fieldmod@proc.f90 fieldseedmod.f90 fieldseedmod@proc.f90
  flux.f90 flux@proc.f90 fragmentmod.f90
  fragmentmod@proc.f90 global.f90 global@proc.F90
  graphics.f90 graphics@proc.f90 grid1mod.f90
  grid1mod@proc.f90 grid3mod.f90 grid3mod@proc.f90
  hashmod.f90 hashmod@proc.f90
  hirshfeld.f90 hirshfeld@proc.f90 integration.f90
  integration@proc.f90 keast.f90 meshmod.f90
  meshmod@proc.f90 molcalc.f90 molsymmod.f90 molsymmod@proc.f90 nci.f90 nci@proc.f90
  pi_private.f90 pi_private@proc.f90 qe_private.f90
  qe_private@proc.f90 qtree_basic.f90
  qtree_basic@proc.f90 qtree.f90 qtree_gpaths.f90
  qtree_gpaths@proc.f90 qtree@proc.f90 qtree_tetrawork.f90
  qtree_tetrawork@proc.f90 qtree_utils.f90
  qtree_utils@proc.f90 quadpack.f90 rhoplot.f90
  rhoplot@proc.f90 sigmahole.f90 sigmahole@proc.f90
  spglib.f90 spglib@proc.f90 stm.f90
  stm@proc.f90 struct_drivers.f90 struct_drivers@proc.F90
  surface.f90 surface@proc.f90 systemmod.f90
  systemmod@proc.f90 tools.f90 tools_io.f90 tools_math.f90
  tools_math@lebedev.f90 tools_math@proc.f90 tools_math@rotations.f90 tools_math@sobol.f90
  tools@proc.f90
  tricks.F90 types.f90 types@proc.f90 wfn_private.f90
  wien_private.f90 wien_private@proc.f90 xdm.f90
  xdm@proc.f90 yt.f90 yt@proc.f90
  arithmetic.F90 arithmetic@proc.F90
  molcalc@proc.F90 param.F90 param@proc.f90
  tools_io@proc.F90 wfn_private@proc.F90
  "${CMAKE_CURRENT_BINARY_DIR}/config@proc.F90"
  )
set(SOURCES_C doqhull.c)
set(SOURCES_Fortran_Main critic2.F90)
set(SOURCES_Fortran_Library libcritic2.f90 libcritic2@proc.f90)

## the core object library, compiled once and linked into both the
## executable and the shared library
add_library(critic2_core OBJECT ${SOURCES_Fortran} ${SOURCES_C})

## interface target carrying the compile definitions, include directories,
## and external libraries shared by the core, the executable, the shared
## library, and the GUI. Toolchain-level optimization/debug flags stay in
## the global CMAKE_<LANG>_FLAGS_<CONFIG> (set in the top-level file); only
## semantic flags and definitions belong here.
add_library(critic2_defs INTERFACE)
target_link_libraries(critic2_core PUBLIC critic2_defs)

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
  set(CRITIC2_DATADIR "${CRITIC2_ROOT}")
else()
  set(CRITIC2_DATADIR "${CMAKE_INSTALL_FULL_DATADIR}")
endif()
configure_file(config@proc.F90.in "${CMAKE_CURRENT_BINARY_DIR}/config@proc.F90")
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/config@proc.F90")
  message(WARNING "Stale generated file src/config@proc.F90 found in the source tree; please delete it")
endif()

## executable (console: the command-line interface)
add_executable(critic2 ${SOURCES_Fortran_Main})
target_link_libraries(critic2 PRIVATE critic2_core)

## On Windows with the GUI, build a second, windowed executable
## (critic2-gui.exe) from the same sources. critic2.exe is the console
## CLI (output stays visible); critic2-gui.exe is the GUI with no
## console. They are otherwise identical; critic2-gui.exe just
## defaults to launching the GUI (CRITIC2_GUI_DEFAULT) and carries a
## different icon and subsystem. CRITIC2_GUI_EXE names the target if
## it exists, "" if not, so the shared link/install logic below can
## extend to it uniformly.
set(CRITIC2_GUI_EXE "")
if (WIN32 AND ENABLE_GUI)
  set(CRITIC2_GUI_EXE critic2-gui)
  add_executable(critic2-gui ${SOURCES_Fortran_Main})
  target_link_libraries(critic2-gui PRIVATE critic2_core)
  target_compile_definitions(critic2-gui PRIVATE CRITIC2_GUI_DEFAULT)
  set_target_properties(critic2-gui PROPERTIES WIN32_EXECUTABLE ON)
endif()

## library
if (BUILD_LIBRARY)
  add_library(libcritic2 SHARED ${SOURCES_Fortran_Library})
  set_target_properties(libcritic2 PROPERTIES OUTPUT_NAME critic2 PUBLIC_HEADER "libcritic2.h")
  target_link_libraries(libcritic2 PRIVATE critic2_core)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Intel")
    set_source_files_properties(libcritic2@proc.f90 libcritic2.f90 PROPERTIES COMPILE_FLAGS "-fvisibility=default")
  endif()
endif()

## link a bundled object library into the final targets. Object files do
## not chain reliably through other object libraries in older CMake, so
## each one is linked directly into the executable and the shared library.
function(link_bundled_lib lib)
  target_link_libraries(critic2 PRIVATE ${lib})
  if (CRITIC2_GUI_EXE)
    target_link_libraries(${CRITIC2_GUI_EXE} PRIVATE ${lib})
  endif()
  if (BUILD_LIBRARY)
    target_link_libraries(libcritic2 PRIVATE ${lib})
  endif()
endfunction()

## cubpack
add_subdirectory(cubpack)
link_bundled_lib(cubpack)
## the core uses cubpack modules
add_dependencies(critic2_core cubpack)

## openmp
if (OPENMP_FOUND)
  target_link_libraries(critic2_defs INTERFACE OpenMP::OpenMP_Fortran)
  target_compile_definitions(critic2_defs INTERFACE HAVE_OPENMP)
endif()

## lapack
if (LAPACK_FOUND)
  target_link_options(critic2_defs INTERFACE ${LAPACK_LINKER_FLAGS})
  target_link_libraries(critic2_defs INTERFACE ${LAPACK_LIBRARIES})
  if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(critic2_defs INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:-fexternal-blas>)
  endif()
  target_compile_definitions(critic2_defs INTERFACE HAVE_LAPACK)
else()
  add_subdirectory(oldlibs)
  link_bundled_lib(oldlibs)
endif()

## qhull
add_subdirectory(qhull)
target_include_directories(critic2_defs INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/qhull")
link_bundled_lib(qhull)

## spglib
add_subdirectory(spglib)
link_bundled_lib(spglib)

## json-fortran
add_subdirectory(json-fortran)
link_bundled_lib(json-fortran)
## the core uses json-fortran modules
add_dependencies(critic2_core json-fortran)

## libxc
if (USE_LIBXC AND LIBXC_FOUND)
  target_include_directories(critic2_defs INTERFACE "${LIBXC_INCLUDE_DIRS}")
  target_link_libraries(critic2_defs INTERFACE ${LIBXC_xcf03_LIBRARY} ${LIBXC_xc_LIBRARY})
  target_compile_definitions(critic2_defs INTERFACE HAVE_LIBXC)
endif()

## libcint
if (USE_LIBCINT AND LIBCINT_FOUND)
  target_link_libraries(critic2_defs INTERFACE ${LIBCINT_LIBRARY})
  target_compile_definitions(critic2_defs INTERFACE HAVE_CINT)
endif()

## readline
if (USE_READLINE AND READLINE_FOUND)
  target_link_libraries(critic2_defs INTERFACE ${READLINE_LIBRARY})
  target_compile_definitions(critic2_defs INTERFACE HAVE_READLINE)
endif()

## libnlopt
if (USE_NLOPT AND NLOPT_FOUND)
  target_include_directories(critic2_defs INTERFACE "${NLOPT_INCLUDE_DIRS}")
  target_link_libraries(critic2_defs INTERFACE ${NLOPT_LIBRARIES})
  target_compile_definitions(critic2_defs INTERFACE HAVE_NLOPT)
endif()

## tblite
if (USE_TBLITE AND TBLITE_FOUND)
  if (TBLITE_INCLUDE_DIRS)
    target_include_directories(critic2_defs INTERFACE "${TBLITE_INCLUDE_DIRS}")
  endif()
  target_link_libraries(critic2_defs INTERFACE ${TBLITE_LIBRARIES})
  target_compile_definitions(critic2_defs INTERFACE HAVE_TBLITE)
endif()

## hdf5
if (USE_HDF5 AND HDF5_FOUND)
  target_include_directories(critic2_defs INTERFACE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(critic2_defs INTERFACE ${HDF5_LIBRARIES})
  target_compile_definitions(critic2_defs INTERFACE HAVE_HDF5)
endif()

## emd (earth mover's distance)
add_subdirectory(emd)
target_include_directories(critic2_defs INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/emd")
link_bundled_lib(emd)

## graphical user interface. The GUI libraries (GLFW, OpenGL, freetype,
## dl, threads) are attached to the fgui target in gui/CMakeLists.txt and
## propagate to the executable through the link.
if (ENABLE_GUI)
  add_subdirectory(gui)
  target_link_libraries(critic2 PRIVATE fgui)
  if (CRITIC2_GUI_EXE)
    target_link_libraries(${CRITIC2_GUI_EXE} PRIVATE fgui)
  endif()

  ## the GUI reads the core module files
  add_dependencies(fgui critic2_core)

  ## add the HAVE_GUI definition
  target_compile_definitions(critic2_defs INTERFACE HAVE_GUI)
endif()

## cross-compilation for windows, add the flag (only needed for the
## Fortran preprocessor; C/C++ compilers predefine it)
if (WIN32)
  target_compile_definitions(critic2_defs INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:_WIN32>)
  ## critic2.exe stays a console application (WIN32_EXECUTABLE OFF) so the CLI
  ## output is visible; critic2-gui.exe is the windowed one (set above).
  ##
  ## Each executable embeds, via a resource script compiled by windres:
  ##  - the application manifest (resource type 24 = RT_MANIFEST, id 1): sets the
  ##    active code page to UTF-8 so narrow char* paths -- argv and the UTF-8
  ##    paths the GLFW drag-and-drop callback delivers -- are read as UTF-8 by the
  ##    UCRT and gfortran's OPEN (accented paths otherwise fail on CP1252).
  ##    Requires Windows 10 1903+.
  ##  - the application icon (id 1): the terminal-style icon for the console CLI,
  ##    the colour icon for the GUI. Explorer and the shortcuts use it.
  ## windres (RC) ships with the MinGW binutils.
  enable_language(RC)
  ## the icons live with the other visual assets in dat/assets; the console CLI
  ## gets the terminal icon, the GUI the colour one (same as the GUI window icon)
  set(_icondir "${CMAKE_SOURCE_DIR}/dat/assets")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/critic2_cli.rc"
    "1 24 \"${CMAKE_CURRENT_SOURCE_DIR}/critic2.manifest\"\n1 ICON \"${_icondir}/critic2.ico\"\n")
  target_sources(critic2 PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/critic2_cli.rc")
  if (CRITIC2_GUI_EXE)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/critic2_gui.rc"
      "1 24 \"${CMAKE_CURRENT_SOURCE_DIR}/critic2.manifest\"\n1 ICON \"${_icondir}/critic2-gui.ico\"\n")
    target_sources(${CRITIC2_GUI_EXE} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/critic2_gui.rc")
  endif()
endif()

## static build
if (BUILD_STATIC)
  target_link_options(critic2 PRIVATE -static)
  target_link_libraries(critic2 PRIVATE ${LIBDL_LIBRARY})
  if (CRITIC2_GUI_EXE)
    target_link_options(${CRITIC2_GUI_EXE} PRIVATE -static)
    target_link_libraries(${CRITIC2_GUI_EXE} PRIVATE ${LIBDL_LIBRARY})
  endif()
endif()

## install the binary (both executables on a Windows GUI build)
install(TARGETS critic2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
if (CRITIC2_GUI_EXE)
  install(TARGETS ${CRITIC2_GUI_EXE} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
endif()
if (BUILD_LIBRARY)
  install(TARGETS libcritic2
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT devel
	  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT devel)
endif()

## bundle the runtime DLLs on Windows
if (WIN32)
  include(${CMAKE_SOURCE_DIR}/cmake/BundleDLLs.cmake)
endif()

## software-rendering GUI build. Point MESA_DIR at a directory with Mesa's
## opengl32.dll and its dependencies (libgallium_wgl.dll, libglapi.dll,
## dxil.dll, ...). A second copy of the GUI executable is installed into a
## sibling bin-mesa/ next to those DLLs, so its static opengl32 import resolves
## to Mesa (guiForceSoftwareGLIfBundled then forces the llvmpipe software
## rasterizer). The compiler runtime DLLs (libgfortran, glfw3, ...) are NOT
## duplicated here: they are shared from bin/, which the "software rendering"
## launcher reaches by setting the working directory to bin/. When MESA_DIR is
## unset, no software build is bundled.
set(CRITIC2_HAVE_MESA OFF)
if (WIN32 AND ENABLE_GUI AND MESA_DIR)
  file(GLOB _mesa_dlls "${MESA_DIR}/*.dll")
  if (_mesa_dlls)
    install(FILES ${_mesa_dlls} DESTINATION "bin-mesa" COMPONENT runtime)
    install(TARGETS ${CRITIC2_GUI_EXE} RUNTIME DESTINATION "bin-mesa" COMPONENT runtime)
    set(CRITIC2_HAVE_MESA ON)
  else()
    message(WARNING "MESA_DIR='${MESA_DIR}' has no DLLs; the software OpenGL build will not be bundled")
  endif()
endif()
## visible to the top-level packaging config (cmake/CritiCPack.cmake), which
## adds the "critic2 (software rendering)" installer shortcut when set
set(CRITIC2_HAVE_MESA ${CRITIC2_HAVE_MESA} PARENT_SCOPE)

## Portable launcher executables at the package root, so users need not dig into
## bin/ (and never run the bin-mesa copy directly, which needs bin/ for its shared
## runtime DLLs). Each is a tiny GUI-subsystem exe (no console flash) that runs
## its target relative to its own location -- so they work from any unzip path --
## and carries the matching critic2 icon. See src/gui/launcher.c.
if (WIN32)
  function(critic2_launcher tgt outname ltarget lworkdir icon console)
    add_executable(${tgt} gui/launcher.c)
    ## a separate output dir: the GUI launcher's OUTPUT_NAME "critic2" would
    ## otherwise clobber the core critic2.exe in the build tree
    set_target_properties(${tgt} PROPERTIES WIN32_EXECUTABLE ON OUTPUT_NAME "${outname}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/launchers")
    ## -municode: use the wide entry point (wWinMain) so paths are Unicode
    target_compile_options(${tgt} PRIVATE -municode)
    target_link_options(${tgt} PRIVATE -municode)
    target_compile_definitions(${tgt} PRIVATE
      "LAUNCH_TARGET=L\"${ltarget}\"" "LAUNCH_WORKDIR=L\"${lworkdir}\"")
    if (console)
      target_compile_definitions(${tgt} PRIVATE LAUNCH_CONSOLE)
    endif()
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${tgt}.rc" "1 ICON \"${icon}\"\n")
    target_sources(${tgt} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/${tgt}.rc")
    install(TARGETS ${tgt} RUNTIME DESTINATION "." COMPONENT runtime)
  endfunction()

  set(_icondir "${CMAKE_SOURCE_DIR}/dat/assets")
  if (ENABLE_GUI)
    critic2_launcher(critic2_launch_gui "critic2" "bin/critic2-gui.exe" "bin"
                     "${_icondir}/critic2-gui.ico" OFF)
    critic2_launcher(critic2_launch_cli "critic2 (console)" "bin/critic2.exe" "bin"
                     "${_icondir}/critic2.ico" ON)
    if (CRITIC2_HAVE_MESA)
      critic2_launcher(critic2_launch_sw "critic2 (software rendering)"
                       "bin-mesa/critic2-gui.exe" "bin" "${_icondir}/critic2-gui.ico" OFF)
    endif()
  else()
    critic2_launcher(critic2_launch_cli "critic2" "bin/critic2.exe" "bin"
                     "${_icondir}/critic2.ico" ON)
  endif()
endif()
