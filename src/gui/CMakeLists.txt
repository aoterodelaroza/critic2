set(SOURCES
  imgui/cimgui.cpp imgui/imgui_impl_glfw.cpp imgui/imgui_impl_opengl3.cpp
  imgui/imgui.cpp imgui/imgui_demo.cpp imgui/imgui_draw.cpp
  imgui/imgui_tables.cpp imgui/imgui_widgets.cpp imgui/gl3w.c

  imgui/implot.cpp imgui/implot_demo.cpp imgui/implot_items.cpp imgui/cimplot.cpp

  threads/tinycthread.c threads/tinycwrapper.c

  filedialog/ImGuiFileDialog.cpp

  gui_interfaces_cimgui.f90 gui_interfaces_glfw.f90
  gui_interfaces_opengl3.f90 gui_interfaces_threads.f90
  interfaces/constants.cpp interfaces/utils.c

  my_imgui/my_imgui.cpp

  gui_main.f90 gui_main@proc.f90
  gui_keybindings.f90 gui_keybindings@proc.f90
  gui_window.f90 gui_window@proc.f90
  gui_utils.f90 gui_utils@proc.f90
  gui_templates.f90 gui_templates@proc.f90
  shapes.f90 shapes@proc.f90
  )
if (FREETYPE2_FOUND)
  set(SOURCES ${SOURCES} imgui/imgui_freetype.cpp)
  add_definitions("-DIMGUI_ENABLE_FREETYPE")
  set_property(SOURCE imgui/cimgui.cpp APPEND PROPERTY COMPILE_FLAGS "-DCIMGUI_FREETYPE")
endif()

add_library(fgui OBJECT ${SOURCES})

target_include_directories(fgui PRIVATE "${CMAKE_BINARY_DIR}/src/")
target_include_directories(fgui PRIVATE "${FREETYPE2_INCLUDE_DIRS}")
