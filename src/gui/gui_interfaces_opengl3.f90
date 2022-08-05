module gui_interfaces_opengl3
  use iso_c_binding
  implicit none

  public

  !xx! OpenGL constants
  include "interfaces/opengl3_constants.inc"

  interface
     !xx! GL3W
     include "interfaces/gl3w_proc.inc"

     !xx! OpenGL
     include "interfaces/opengl3_proc.inc"
  end interface

end module gui_interfaces_opengl3
