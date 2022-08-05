module gui_interfaces_cimgui
  use iso_c_binding
  implicit none

  public

  !xx! cimgui constants
  include "interfaces/cimgui_constants.inc"

  !xx! cimgui user-defined types
  include "interfaces/cimgui_types.inc"

  interface
     !xx! ImGui GLFW platform backend (imgui_impl_glfw)
     include "interfaces/imgui_glfw_proc.inc"

     !xx! ImGui OpenGL3 renderer backend (imgui_impl_opengl3)
     include "interfaces/imgui_opengl3_proc.inc"

     !xx! cimgui
     include "interfaces/cimgui_proc.inc"
  end interface

end module gui_interfaces_cimgui
