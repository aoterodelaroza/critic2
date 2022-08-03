set(SOURCES
  cimgui.cpp imgui_impl_glfw.cpp imgui_impl_opengl3.cpp
  imgui/imgui.cpp imgui/imgui_demo.cpp imgui/imgui_draw.cpp
  imgui/imgui_tables.cpp imgui/imgui_widgets.cpp imgui/imgui_freetype.cpp
  imgui/gl3w.c
  gui_main.f90 gui_main@proc.f90
  gui_interfaces.f90
  )
add_library(fgui OBJECT ${SOURCES})

target_include_directories(fgui PRIVATE "${CMAKE_BINARY_DIR}/src/")
target_include_directories(fgui PRIVATE "${FREETYPE2_INCLUDE_DIRS}")
