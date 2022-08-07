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

  ! Start the critic2 GUI.
  module subroutine gui_start()
    use gui_interfaces_cimgui
    use gui_interfaces_glfw
    use gui_interfaces_opengl3
    use gui_keybindings, only: is_bind_event, BIND_QUIT, set_default_keybindings
    use c_interface_module, only: f_c_string_dup, C_string_free
    use tools_io, only: ferror, faterr, string
    integer(c_int) :: idum, display_w, display_h
    type(c_funptr) :: fdum
    type(c_ptr) :: ptrc
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
    rootwin = glfwCreateWindow(1280, 720, c_loc(strc), c_null_ptr, c_null_ptr)
    if (.not.c_associated(rootwin)) &
       call ferror('gui_start','Failed to create window',faterr)
    call glfwMakeContextCurrent(rootwin)
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
    call glfwSetInputMode(rootwin, GLFW_STICKY_KEYS, 1)
    call igStyleColorsDark(c_null_ptr)

    ! set up backend and renderer
    ldum = ImGui_ImplGlfw_InitForOpenGL(rootwin, .true._c_bool)
    if (.not.ldum)&
       call ferror('gui_start','Failed to initialize ImGui (GLFW for OpenGL)',faterr)
    strc = shader_version
    ldum = ImGui_ImplOpenGL3_Init(c_loc(strc))
    if (.not.ldum)&
       call ferror('gui_start','Failed to initialize ImGui (OpenGL)',faterr)

    ! get the ImGUI IO interface and enable docking
    ptrc = igGetIO()
    call c_f_pointer(ptrc,io)
    io%configflags = ior(io%configflags,ImGuiConfigFlags_DockingEnable) ! activate docking
    io%inifilename = c_null_ptr ! do not save the ini file, for now

    ! get the ImGUI context pointer
    ptrc = igGetCurrentContext()
    call c_f_pointer(ptrc,g)

    ! set default keybindings
    call set_default_keybindings()

    ! main loop
    show_demo_window = .true.
    do while (glfwWindowShouldClose(rootwin) == 0)
       ! poll events
       call glfwPollEvents()
       time = glfwGetTime()

       ! start the ImGui frame
       call ImGui_ImplOpenGL3_NewFrame()
       call ImGui_ImplGlfw_NewFrame()
       call igNewFrame()

       ! handle quit key binding
       if (is_bind_event(BIND_QUIT)) &
          call glfwSetWindowShouldClose(rootwin, GLFW_TRUE)

       ! show main menu
       call show_main_menu()

       ! show main dockspace
       idum = igDockSpaceOverViewport(igGetMainViewport(),&
          ImGuiDockNodeFlags_PassthruCentralNode,&
          c_null_ptr)

       ! show demo window
       if (show_demo_window) &
          call igShowDemoWindow(show_demo_window)

       ! rendering
       call igRender()
       call glfwGetFramebufferSize(rootwin, display_w, display_h)
       call glViewport(0, 0, display_w, display_h)
       call glClearColor(0.45, 0.55, 0.60, 1.00)
       call glClear(GL_COLOR_BUFFER_BIT)
       call ImGui_ImplOpenGL3_RenderDrawData(igGetDrawData());

       ! swap buffers
       call glfwSwapBuffers(rootwin)
    end do

    ! cleanup
    call ImGui_ImplOpenGL3_Shutdown()
    call ImGui_ImplGlfw_Shutdown()
    call igDestroyContext(c_null_ptr)

    ! terminate
    call glfwDestroyWindow(rootwin)
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

  !xx! private procedures

  ! Show the main menu
  subroutine show_main_menu()
    use gui_interfaces_cimgui
    use gui_utils, only: igIsItemHovered_delayed
    use gui_keybindings, only: BIND_QUIT, get_bind_keyname
    use gui_interfaces_glfw, only: GLFW_TRUE, glfwSetWindowShouldClose
    use tools_io, only: string

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: v2

    if (igBeginMainMenuBar()) then
       ! File
       str1 = "File" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then

          ! File -> Quit
          str1 = "Quit" // c_null_char
          str2 = get_bind_keyname(BIND_QUIT) // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_loc(str2),.false._c_bool,.true._c_bool)) &
             call glfwSetWindowShouldClose(rootwin, GLFW_TRUE)
          if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay)) then
             str1 = "Quit the program" // c_null_char
             call igSetTooltip(c_loc(str1))
          end if

          call igEndMenu()
       end if

       !     if (BeginMenu("Edit")){
       !       if (MenuItem("Preferences..."))
       !        OpenDialog(DLG_Preferences);
       !       EndMenu();
       !     }

       !     if (BeginMenu("View")){
       !       if (MenuItem("Tree",NULL,dlgopen[DLG_Tree]))
       !        ToggleDialog(DLG_Tree);
       !       if (MenuItem("Preferences",NULL,dlgopen[DLG_Preferences]))
       !        ToggleDialog(DLG_Preferences);
       !       if (MenuItem("Structural Information",NULL,dlgopen[DLG_StructInfo]))
       !        ToggleDialog(DLG_StructInfo);
       !       EndMenu();
       !     }

       call igGetContentRegionAvail(v2)
       call igSameLine(0._c_float, v2%x - 180._c_float)
       str1 = string(1000._c_float / io%Framerate,'f',decimal=3) // " ms/frame (" // &
          string(io%Framerate,'f',decimal=1) // " FPS)" // c_null_char
       call igText(c_loc(str1))
    end if
    call igEndMainMenuBar()

  end subroutine show_main_menu
end submodule proc
