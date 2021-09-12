## source
set(SOURCES dialog.cpp keybinding.cpp main.cpp menu.cpp message.cpp
  scene.cpp settings.cpp shapes.cpp text.cpp tree.cpp view.cpp
  imgui/imgui.cpp imgui/imgui_demo.cpp imgui/imgui_dock.cpp 
  imgui/imgui_draw.cpp imgui/imgui_impl_glfw_gl3.cpp imgui/imgui_widgets.cpp
  imgui/gl3w.c
  )

## some compilation optoins
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## executable
add_executable(gcritic2 ${SOURCES})

## pkgconfig and gl-related libraries
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
pkg_search_module(FREETYPE2 REQUIRED freetype2)
pkg_search_module(GLM REQUIRED glm)
pkg_search_module(GL REQUIRED gl)
target_include_directories(gcritic2 PRIVATE
  ${GLFW_INCLUDE_DIRS} 
  ${FREETYPE2_INCLUDE_DIRS} 
  ${GLM_INCLUDE_DIRS}
  ${GL_INCLUDE_DIRS}
  )
target_link_libraries(gcritic2 PRIVATE
  ${GLFW_LIBRARIES} 
  ${FREETYPE2_LIBRARIES} 
  ${GLM_LIBRARIES} 
  ${GL_LIBRARIES}
  )

## the critic2 libraries
target_link_libraries(gcritic2 PRIVATE lcritic2 cubpack ciftbx json-fortran)

## the dynamic linking library
target_link_libraries(gcritic2 PRIVATE ${CMAKE_DL_LIBS})

## openmp
if (OPENMP_FOUND)
  target_link_libraries(gcritic2 PRIVATE OpenMP::OpenMP_Fortran)
endif()

## lapack
if (LAPACK_FOUND)
  target_link_libraries(gcritic2 PRIVATE ${LAPACK_LIBRARIES})
else()
  target_link_libraries(gcritic2 PRIVATE oldlibs)
endif()

## qhull
if (QHULL_FOUND)
  target_link_libraries(gcritic2 PRIVATE ${QHULL_LIBRARIES})
else()
  target_link_libraries(gcritic2 PRIVATE qhull)
endif()

## spglib
if (SPGLIB_FOUND)
  target_link_libraries(gcritic2 PRIVATE ${SPGLIB_LIBRARY})
else()
  target_link_libraries(gcritic2 PRIVATE spglib)
endif()

## libxc
if (LIBXC_FOUND)
  target_link_libraries(gcritic2 PRIVATE ${LIBXC_xcf90_LIBRARY} ${LIBXC_xc_LIBRARY})
endif()

## libcint
if (LIBCINT_FOUND)
  target_link_libraries(gcritic2 PRIVATE ${LIBCINT_LIBRARY})
endif()

## readline
if (READLINE_FOUND)
  target_link_libraries(gcritic2 PRIVATE ${READLINE_LIBRARY})
endif()

## install the binary
install(TARGETS gcritic2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
