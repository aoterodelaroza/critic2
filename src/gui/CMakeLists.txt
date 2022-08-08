set(SOURCES
  imgui/cimgui.cpp imgui/imgui_impl_glfw.cpp imgui/imgui_impl_opengl3.cpp
  imgui/imgui.cpp imgui/imgui_demo.cpp imgui/imgui_draw.cpp
  imgui/imgui_tables.cpp imgui/imgui_widgets.cpp imgui/imgui_freetype.cpp
  imgui/gl3w.c

  gui_interfaces_cimgui.f90 gui_interfaces_glfw.f90
  gui_interfaces_opengl3.f90 interfaces/constants.cpp

  gui_main.f90 gui_main@proc.f90
  gui_keybindings.f90 gui_keybindings@proc.f90
  gui_utils.f90 gui_utils@proc.f90
  gui_window.f90 gui_window@proc.f90
  )
add_library(fgui OBJECT ${SOURCES})

target_include_directories(fgui PRIVATE "${CMAKE_BINARY_DIR}/src/")
target_include_directories(fgui PRIVATE "${FREETYPE2_INCLUDE_DIRS}")
