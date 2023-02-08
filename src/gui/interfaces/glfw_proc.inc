! -*-F90-*-

!xx! GLFW procedures
! GLFWerrorfun glfwSetErrorCallback(GLFWerrorfun callback)
function glfwSetErrorCallback(callback) bind(c,name="glfwSetErrorCallback")
  import c_funptr
  type(c_funptr) :: glfwSetErrorCallback
  type(c_funptr), value :: callback
end function glfwSetErrorCallback
! GLFWdropfun glfwSetDropCallback (GLFWwindow * window, GLFWdropfun callback)
function glfwSetDropCallback(window,callback) bind(c,name="glfwSetDropCallback")
  import c_ptr,c_funptr
  type(c_funptr) :: glfwSetDropCallback
  type(c_ptr), value :: window
  type(c_funptr), value :: callback
end function glfwSetDropCallback
! int glfwInit(void)
function glfwInit() bind(c,name="glfwInit")
  import c_int
  integer(c_int) :: glfwInit
end function glfwInit
! void glfwWindowHint(int hint, int value)
subroutine glfwWindowHint(hint,val) bind(c,name="glfwWindowHint")
  import c_int
  integer(c_int), value :: hint, val
end subroutine glfwWindowHint
! GLFWwindow *glfwCreateWindow(int width, int height, const char *title, GLFWmonitor *monitor, GLFWwindow *share)
function glfwCreateWindow(width, height, title, monitor, share) bind(c,name="glfwCreateWindow")
  import c_int, c_ptr
  integer(c_int), value :: width, height
  type(c_ptr), value, intent(in) :: title
  type(c_ptr), value :: monitor
  type(c_ptr), value :: share
  type(c_ptr) :: glfwCreateWindow
end function glfwCreateWindow
! void glfwMakeContextCurrent(GLFWwindow *window)
subroutine glfwMakeContextCurrent(window) bind(c,name="glfwMakeContextCurrent")
  import c_ptr
  type(c_ptr), value :: window
end subroutine glfwMakeContextCurrent
! void glfwSetInputMode(GLFWwindow *window, int mode, int value)
subroutine glfwSetInputMode(window, mode, val) bind(c,name="glfwSetInputMode")
  import c_ptr, c_int
  type(c_ptr), value :: window
  integer(c_int), value :: mode, val
end subroutine glfwSetInputMode
! void glfwGetFramebufferSize(GLFWwindow *window, int *width, int *height)
subroutine glfwGetFramebufferSize(window, width, height) bind(c,name="glfwGetFramebufferSize")
  import c_int, c_ptr
  type(c_ptr), value :: window
  integer(c_int) :: width, height
end subroutine glfwGetFramebufferSize
! int glfwWindowShouldClose(GLFWwindow *window)
function glfwWindowShouldClose(window) bind(c,name="glfwWindowShouldClose")
  import c_ptr, c_int
  type(c_ptr), value :: window
  integer(c_int) :: glfwWindowShouldClose
end function glfwWindowShouldClose
! void glfwSetWindowIcon(GLFWwindow *window, int count, const GLFWimage *images)
subroutine glfwSetWindowIcon(window, count, images) bind(c,name="glfwSetWindowIcon")
  import c_ptr, c_int
  type(c_ptr), value :: window
  integer(c_int), value :: count
  type(c_ptr), intent(in), value :: images
end subroutine glfwSetWindowIcon
! void glfwSetWindowShouldClose(GLFWwindow *window, int value)
subroutine glfwSetWindowShouldClose(window, val) bind(c,name="glfwSetWindowShouldClose")
  import c_ptr, c_int
  type(c_ptr), value :: window
  integer(c_int), value :: val
end subroutine glfwSetWindowShouldClose
! void glfwSwapInterval(int interval)
subroutine glfwSwapInterval(interval) bind(c,name="glfwSwapInterval")
  import c_int
  integer(c_int), value :: interval
end subroutine glfwSwapInterval
! void glfwSwapBuffers(GLFWwindow *window)
subroutine glfwSwapBuffers(window) bind(c,name="glfwSwapBuffers")
  import c_ptr
  type(c_ptr), value :: window
end subroutine glfwSwapBuffers
! void glfwPollEvents(void)
subroutine glfwPollEvents() bind(c,name="glfwPollEvents")
end subroutine glfwPollEvents
! double glfwGetTime(void)
function glfwGetTime() bind(c,name="glfwGetTime")
  import c_double
  real(c_double) :: glfwGetTime
end function glfwGetTime
! void glfwDestroyWindow(GLFWwindow * window)
subroutine glfwDestroyWindow(window) bind(c,name="glfwDestroyWindow")
  import c_ptr
  type(c_ptr), value :: window
end subroutine glfwDestroyWindow
! void glfwTerminate(void)
subroutine glfwTerminate() bind(c,name="glfwTerminate")
end subroutine glfwTerminate
