set(SOURCES
  cimgui/cimgui.cpp cimgui/cimguiwrapper.cpp

  imgui/imgui_impl_glfw.cpp imgui/imgui_impl_opengl3.cpp
  imgui/imgui.cpp imgui/imgui_demo.cpp imgui/imgui_draw.cpp
  imgui/imgui_tables.cpp imgui/imgui_widgets.cpp
  imgui/glfwwrapper.cpp

  implot/implot.cpp implot/implot_demo.cpp implot/implot_items.cpp implot/cimplot.cpp
  implot/implotwrapper.cpp

  tinycthread/tinycthread.c tinycthread/tinycwrapper.c

  stb/stb.c

  filedialog/ImGuiFileDialog.cpp filedialog/filedialogwrapper.cpp

  interfaces_cimgui.f90 interfaces_glfw.f90
  interfaces_opengl3.f90 interfaces_threads.f90 interfaces_stb.f90

  c_utils/utils.c

  imgui_own_widgets/imgui_own_widgets.cpp

  opengl/openglwrapper.c opengl/gl3w.c

  gui_main.f90 gui_main@proc.F90
  keybindings.f90 keybindings@proc.f90
  windows.f90 windows@proc.f90 windows@tree.f90 windows@view.f90 windows@dialog.f90
  windows@ci.f90 windows@co.f90 windows@new_struct.f90
  windows@edit_geometry.f90 windows@treeplot.f90 windows@rebond.f90 windows@scfplot.f90
  windows@load_field.f90 windows@preferences.f90 windows@about.f90
  utils.f90 utils@proc.f90 utils@math.f90
  templates.f90 templates@proc.f90
  shapes.f90 shapes@proc.f90
  shaders.f90 shaders@proc.f90
  scenes.f90 scenes@proc.f90
  )

## cross-compilation for windows, add the flag
if (WIN32)
  add_definitions("-D_WIN32")
endif()
if (USE_GUI_THREADS)
  add_definitions("-D_THREADS")
endif()

## freetype
if (FREETYPE_FOUND)
  set(SOURCES ${SOURCES} imgui/imgui_freetype.cpp)
  add_definitions("-DIMGUI_ENABLE_FREETYPE")
  set_property(SOURCE cimgui/cimgui.cpp APPEND PROPERTY COMPILE_FLAGS "-DCIMGUI_FREETYPE")
endif()

add_library(fgui OBJECT ${SOURCES})

## include directories
target_include_directories(fgui PRIVATE "${GLFW3_INCLUDE_DIR}")
target_include_directories(fgui PRIVATE "${CMAKE_BINARY_DIR}/src/")
target_include_directories(fgui PRIVATE "${FREETYPE_INCLUDE_DIRS}")
