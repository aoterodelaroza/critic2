! -*-F90-*-

!xx! GL3W procedures
! int gl3wInit(void)
function gl3wInit() bind(c,name="gl3wInit")
  import c_int
  integer(c_int) :: gl3wInit
end function gl3wInit
! int gl3wIsSupported(int major, int minor)
function gl3wIsSupported(major, minor)
  import c_int
  integer(c_int), value :: major, minor
  integer(c_int) :: gl3wIsSupported
end function gl3wIsSupported
! void *gl3wGetProcAddress(const char *proc)
function gl3wGetProcAddress(proc)
  import c_ptr
  type(c_ptr), intent(in), value :: proc
  type(c_ptr) :: gl3wGetProcAddress
end function gl3wGetProcAddress
! void glfwGetFramebufferSize(GLFWwindow *window, int *width, int *height)
subroutine glfwGetFramebufferSize(window, width, height) bind(c,name="glfwGetFramebufferSize")
  import c_int, c_ptr
  type(c_ptr), value :: window
  integer(c_int) :: width, height
end subroutine glfwGetFramebufferSize
