module gui_interfaces
  use iso_c_binding
  implicit none

  public

  !xx! GLFW constants
  include "interfaces/glfw_constants.inc"

  !xx! OpenGL constants
  include "interfaces/opengl3_constants.inc"

  interface
     !xx! GLFW
     include "interfaces/glfw_proc.inc"

     !xx! GL3W
     include "interfaces/gl3w_proc.inc"

     !xx! ImGui GLFW backend (imgui_impl_glfw)
     include "interfaces/imgui_glfw_proc.inc"

     !xx! ImGui OpenGL3 backend (imgui_impl_opengl3)
     include "interfaces/imgui_opengl3_proc.inc"

     !xx! cimgui
     include "interfaces/cimgui_proc.inc"

     !xx! OpenGL
     include "interfaces/opengl3_proc.inc"
  end interface

end module gui_interfaces
