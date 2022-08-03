module gui_interfaces
  use iso_c_binding
  implicit none

  public

  ! GLFW constants
  integer(c_int), parameter :: GLFW_SAMPLES                = int(z'0002100D',c_int)
  integer(c_int), parameter :: GLFW_CONTEXT_VERSION_MAJOR  = int(z'00022002',c_int)
  integer(c_int), parameter :: GLFW_CONTEXT_VERSION_MINOR  = int(z'00022003',c_int)
  integer(c_int), parameter :: GLFW_OPENGL_PROFILE         = int(z'00022008',c_int)
  integer(c_int), parameter :: GLFW_OPENGL_FORWARD_COMPAT  = int(z'00022006',c_int)
  integer(c_int), parameter :: GLFW_OPENGL_CORE_PROFILE    = int(z'00032001',c_int)
  integer(c_int), parameter :: GLFW_STICKY_KEYS            = int(z'00033002',c_int)

  ! OpenGL constants
  integer(c_int), parameter :: GL_COLOR_BUFFER_BIT = int(z'00004000',c_int)

  interface
!!!!!!!!!!! GLFW !!!!!!!!!!!!
     ! GLFWerrorfun glfwSetErrorCallback(GLFWerrorfun callback)
     function glfwSetErrorCallback(callback) bind(c,name="glfwSetErrorCallback")
       import c_funptr
       type(c_funptr) :: glfwSetErrorCallback
       type(c_funptr), value :: callback
     end function glfwSetErrorCallback
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
     ! int glfwWindowShouldClose(GLFWwindow *window)
     function glfwWindowShouldClose(window) bind(c,name="glfwWindowShouldClose")
       import c_ptr, c_int
       type(c_ptr), value :: window
       integer(c_int) :: glfwWindowShouldClose
     end function glfwWindowShouldClose
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
     ! void glfwDestroyWindow(GLFWwindow * window)
     subroutine glfwDestroyWindow(window) bind(c,name="glfwDestroyWindow")
       import c_ptr
       type(c_ptr), value :: window
     end subroutine glfwDestroyWindow
     ! void glfwTerminate(void)
     subroutine glfwTerminate() bind(c,name="glfwTerminate")
     end subroutine glfwTerminate

!!!!!!!!!!! GL3W !!!!!!!!!!!
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

!!!!!!!!!!! ImGui GLFW backend !!!!!!!!!!!
     ! bool ImGui_ImplGlfw_InitForOpenGL(GLFWwindow* window, bool install_callbacks);
     function ImGui_ImplGlfw_InitForOpenGL(window, install_callbacks) bind(c,name="ImGui_ImplGlfw_InitForOpenGL")
       import c_bool, c_ptr
       type(c_ptr), value :: window
       logical(c_bool), value :: install_callbacks
       logical(c_bool) :: ImGui_ImplGlfw_InitForOpenGL
     end function ImGui_ImplGlfw_InitForOpenGL
     ! void ImGui_ImplGlfw_NewFrame(void);
     subroutine ImGui_ImplGlfw_NewFrame() bind(c,name="ImGui_ImplGlfw_NewFrame")
     end subroutine ImGui_ImplGlfw_NewFrame
     ! void ImGui_ImplGlfw_Shutdown(void);
     subroutine ImGui_ImplGlfw_Shutdown() bind(c,name="ImGui_ImplGlfw_Shutdown")
     end subroutine ImGui_ImplGlfw_Shutdown
     ! bool ImGui_ImplOpenGL3_Init(const char* glsl_version);
     function ImGui_ImplOpenGL3_Init(glsl_version) bind(c,name="ImGui_ImplOpenGL3_Init")
       import c_bool, c_ptr
       type(c_ptr), value, intent(in) :: glsl_version
       logical(c_bool) :: ImGui_ImplOpenGL3_Init
     end function ImGui_ImplOpenGL3_Init
     ! void ImGui_ImplOpenGL3_NewFrame(void);
     subroutine ImGui_ImplOpenGL3_NewFrame() bind(c,name="ImGui_ImplOpenGL3_NewFrame")
     end subroutine ImGui_ImplOpenGL3_NewFrame
     ! void ImGui_ImplOpenGL3_RenderDrawData(ImDrawData* draw_data);
     subroutine ImGui_ImplOpenGL3_RenderDrawData(draw_data) bind(c,name="ImGui_ImplOpenGL3_RenderDrawData")
       import c_ptr
       type(c_ptr), value :: draw_data
     end subroutine ImGui_ImplOpenGL3_RenderDrawData
     ! void ImGui_ImplOpenGL3_Shutdown(void);
     subroutine ImGui_ImplOpenGL3_Shutdown() bind(c,name="ImGui_ImplOpenGL3_Shutdown")
     end subroutine ImGui_ImplOpenGL3_Shutdown

!!!!!!!!!!! cimgui !!!!!!!!!!!
     ! ImGuiContext* igGetCurrentContext(void);
     function igGetCurrentContext() bind(c,name="igGetCurrentContext")
       import c_ptr
       type(c_ptr) :: igGetCurrentContext
     end function igGetCurrentContext
     ! ImGuiContext* igCreateContext(ImFontAtlas* shared_font_atlas);
     function igCreateContext(shared_font_atlas) bind(c,name="igCreateContext")
       import c_ptr
       type(c_ptr), value :: shared_font_atlas
       type(c_ptr) :: igCreateContext
     end function igCreateContext
     ! void igDestroyContext(ImGuiContext* ctx);
     subroutine igDestroyContext(ctx) bind(c,name="igDestroyContext")
       import c_ptr
       type(c_ptr), value :: ctx
     end subroutine igDestroyContext
     ! void igStyleColorsDark(ImGuiStyle* dst);
     subroutine igStyleColorsDark(dst) bind(c,name="igStyleColorsDark")
       import c_ptr
       type(c_ptr), value :: dst
     end subroutine igStyleColorsDark
     ! void igNewFrame(void);
     subroutine igNewFrame() bind(c,name="igNewFrame")
     end subroutine igNewFrame
     ! void igShowDemoWindow(bool* p_open);
     subroutine igShowDemoWindow(p_open) bind(c,name="igShowDemoWindow")
       import c_bool
       logical(c_bool) :: p_open
     end subroutine igShowDemoWindow
     ! void igRender(void);
     subroutine igRender() bind(c,name="igRender")
     end subroutine igRender
     ! ImDrawData* igGetDrawData(void);
     function igGetDrawData() bind(c,name="igGetDrawData")
       import c_ptr
       type(c_ptr) :: igGetDrawData
     end function igGetDrawData

!!!!!!!!!!! OpenGL !!!!!!!!!!!
     ! GLint = GLsizei = c_int
     ! GLbitfield = c_int (unsigned int)
     ! GLfloat = c_float
     ! void glViewport (GLint x, GLint y, GLsizei width, GLsizei height);
     subroutine glViewport(x, y, width, height) bind(c,name="glViewport")
       import c_int
       integer(c_int), value :: x, y, width, height
     end subroutine glViewport
     ! void glClear (GLbitfield mask);
     subroutine glClear (mask) bind(c,name="glClear")
       import c_int
       integer(c_int), value :: mask
     end subroutine glClear
     ! void glClearColor(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha)
     subroutine glClearColor(red, green, blue, alpha) bind(c,name="glClearColor")
       import c_float
       real(c_float), value :: red, green, blue, alpha
     end subroutine glClearColor
  end interface

end module gui_interfaces
