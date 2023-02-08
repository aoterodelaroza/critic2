  ! -*-F90-*-

  !xx! glfw user-defined types
  type, bind(c) :: GLFWimage
     integer(c_int) :: width
     integer(c_int) :: height
     type(c_ptr) :: pixels
  end type GLFWimage
