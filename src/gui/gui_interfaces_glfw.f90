module gui_interfaces_glfw
  use iso_c_binding
  implicit none

  public

  !xx! GLFW constants
  include "interfaces/glfw_constants.inc"

  interface
     !xx! GLFW
     include "interfaces/glfw_proc.inc"
  end interface

end module gui_interfaces_glfw
