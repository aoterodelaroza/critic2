! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

submodule (gui_main) proc
  use iso_c_binding
  implicit none

  ! opengl & shader version
  integer, parameter :: opengl_version_major = 3
  integer, parameter :: opengl_version_minor = 3
  character(len=*,kind=c_char), parameter :: shader_version = "#version 330"//c_null_char

  ! gui title
  character(len=*,kind=c_char), parameter :: gui_title = "critic2 GUI"//c_null_char

  !xx! private procedures
  ! subroutine show_main_menu()

contains

  module subroutine gui_start()
    use gui_interfaces_cimgui
    use gui_interfaces_glfw
    use gui_interfaces_opengl3
    use gui_keybindings, only: is_bind_event, BIND_QUIT, set_default_keybindings
    use c_interface_module, only: f_c_string_dup, C_string_free
    use tools_io, only: ferror, faterr, string
    integer(c_int) :: idum, display_w, display_h
    type(c_funptr) :: fdum
    type(c_ptr) :: window, ptrc
    logical(c_bool) :: ldum, show_demo_window
    character(kind=c_char,len=:), allocatable, target :: strc

    ! Initialize glfw
    fdum = glfwSetErrorCallback(c_funloc(error_callback))
    if (glfwInit() == 0) &
       call ferror('gui_start','Failed to initialize GLFW',faterr)
    call glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, opengl_version_major)
    call glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, opengl_version_minor)
    call glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE)
    call glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1)
    !xx! call glfwWindowHint(GLFW_SAMPLES, 4) ! activate multisampling

    ! set up window
    strc = gui_title
    window = glfwCreateWindow(1280, 720, c_loc(strc), c_null_ptr, c_null_ptr)
    if (.not.c_associated(window)) &
       call ferror('gui_start','Failed to create window',faterr)
    call glfwMakeContextCurrent(window)
    call glfwSwapInterval(1) ! enable vsync

    ! set up ImGui context
    ptrc = igCreateContext(c_null_ptr)
    if (.not.c_associated(ptrc)) &
       call ferror('gui_start','Failed to create ImGui context',faterr)

    ! initialize gl3w
    idum = gl3wInit()
    if (idum /= 0) &
       call ferror('gui_start','Failed to initialize OpenGL (gl3w)',faterr)
    if (gl3wIsSupported(opengl_version_major,opengl_version_minor) == 0) &
       call ferror('gui_start','gl3w: OpenGL version ' // &
       string(opengl_version_major) // '.' // string(opengl_version_minor) // ' not supported',faterr)

    ! set up ImGui style
    call glfwSetInputMode(window, GLFW_STICKY_KEYS, 1)
    call igStyleColorsDark(c_null_ptr)

    ! set up backend and renderer
    ldum = ImGui_ImplGlfw_InitForOpenGL(window, .true._c_bool)
    if (.not.ldum)&
       call ferror('gui_start','Failed to initialize ImGui (GLFW for OpenGL)',faterr)
    strc = shader_version
    ldum = ImGui_ImplOpenGL3_Init(c_loc(strc))
    if (.not.ldum)&
       call ferror('gui_start','Failed to initialize ImGui (OpenGL)',faterr)

    ! get the ImGUI IO interface and enable docking
    ptrc = igGetIO()
    call c_f_pointer(ptrc,io)
    io%configflags = ior(io%configflags,ImGuiConfigFlags_DockingEnable)

    ! set default keybindings
    call set_default_keybindings()

    ! main loop
    show_demo_window = .true.
    do while (glfwWindowShouldClose(window) == 0)
       ! poll events
       call glfwPollEvents()
       time = glfwGetTime()

       ! start the ImGui frame
       call ImGui_ImplOpenGL3_NewFrame()
       call ImGui_ImplGlfw_NewFrame()
       call igNewFrame()

       ! handle quit key binding
       if (is_bind_event(BIND_QUIT)) &
          call glfwSetWindowShouldClose(window, GLFW_TRUE)

       ! show demo window
       if (show_demo_window) &
          call igShowDemoWindow(show_demo_window)

       ! rendering
       call igRender()
       call glfwGetFramebufferSize(window, display_w, display_h)
       call glViewport(0, 0, display_w, display_h)
       call glClearColor(0.45, 0.55, 0.60, 1.00)
       call glClear(GL_COLOR_BUFFER_BIT)
       call ImGui_ImplOpenGL3_RenderDrawData(igGetDrawData());

       ! swap buffers
       call glfwSwapBuffers(window)
    end do

    ! cleanup
    call ImGui_ImplOpenGL3_Shutdown()
    call ImGui_ImplGlfw_Shutdown()
    call igDestroyContext(c_null_ptr)

    ! terminate
    call glfwDestroyWindow(window)
    call glfwTerminate()

  contains
    subroutine error_callback(error,description) bind(c)
      use c_interface_module, only: c_f_string_alloc, c_strlen
      use tools_io, only: ferror, faterr, string
      integer(c_int), value :: error
      type(c_ptr), intent(in), value :: description

      character(len=:), allocatable :: msg

      call c_f_string_alloc(description,msg)
      call ferror('glfw',"GLFW error (" // string(error) // "): " // trim(msg),faterr)

    end subroutine error_callback
  end subroutine gui_start

end submodule proc


! void show_main_menu(GLFWwindow* rootwin){
!   
!   if (BeginMainMenuBar()){
!     if (BeginMenu("File")){
!       if (MenuItem("Quit",BindKeyName(BIND_QUIT).c_str()))
! 	glfwSetWindowShouldClose(rootwin, GLFW_TRUE);
!       EndMenu();
!     }
!     if (BeginMenu("Edit")){
!       if (MenuItem("Preferences..."))
! 	OpenDialog(DLG_Preferences);
!       EndMenu();
!     }
!     if (BeginMenu("View")){
!       if (MenuItem("Tree",NULL,dlgopen[DLG_Tree]))
! 	ToggleDialog(DLG_Tree);
!       if (MenuItem("Preferences",NULL,dlgopen[DLG_Preferences]))
! 	ToggleDialog(DLG_Preferences);
!       if (MenuItem("Structural Information",NULL,dlgopen[DLG_StructInfo]))
! 	ToggleDialog(DLG_StructInfo);
!       EndMenu();
!     }
!     SameLine(0, GetContentRegionAvailWidth()-180.);
!     Text("%.3f ms/frame (%.1f FPS)", 1000.0f / GetIO().Framerate, GetIO().Framerate);
!   }
!   EndMainMenuBar();
! }
