! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
  use types, only: thread_info
  implicit none

  ! opengl & shader version
  integer, parameter :: opengl_version_major = 3
  integer, parameter :: opengl_version_minor = 3
  character(len=*,kind=c_char), parameter :: shader_version = "#version 330"//c_null_char

  ! gui title
  character(len=*,kind=c_char), parameter :: gui_title = "critic2 GUI"//c_null_char

  ! the dockspace ID
  integer(c_int) :: iddock = 0

  ! whether glfw errors are fatal
  logical :: glfw_lenient = .false.

  !xx! private procedures
  ! subroutine process_arguments()
  ! subroutine show_main_menu()
  ! function initialization_thread_worker(arg)

contains

  ! Start the critic2 GUI.
  module subroutine gui_start()
    use interfaces_threads, only: deallocate_mtx
    use interfaces_cimgui
    use interfaces_glfw
    use interfaces_opengl3
    use interfaces_stb
    use systems, only: sys, sysc, nsys, launch_initialization_thread, system_shorten_names,&
       thread, thread_ti, nthread
    use shaders, only: shaders_init, shaders_end
    use shapes, only: shapes_init, shapes_end
    use windows, only: nwin, win, wintype_tree, wintype_view, wintype_console_input,&
       wintype_console_output, wintype_about, wintype_builder, iwin_tree, iwin_view,&
       iwin_console_input, iwin_console_output, iwin_about, iwin_builder,&
       stack_create_window, stack_realloc_maybe, wpurp_view_main, windows_init
    use global, only: critic_home
    use c_interface_module, only: f_c_string_dup, C_string_free
    use tools_io, only: ferror, faterr, string, falloc, fdealloc
    use param, only: dirsep
    integer(c_int) :: idum, display_w, display_h, ileft, iright
    type(c_funptr) :: fdum
    type(c_ptr) :: ptrc
    logical(c_bool) :: ldum, show_demo_window, show_implot_demo_window
    character(kind=c_char,len=:), allocatable, target :: strc, file
    integer :: i, j, ludum(10), saveinpcon
    logical :: firstpass
    integer(c_short), allocatable, target :: range(:)
    type(GLFWimage), target :: icon

    integer, parameter :: initial_nsys = 20

    ! initialize the sys arrays
    nsys = 0
    allocate(sys(initial_nsys),sysc(initial_nsys))

    ! initialize threads: reserve some un-used LUs and then
    ! deallocate them for fopen
    if (allocated(thread)) deallocate(thread)
    allocate(thread(nthread),thread_ti(nthread))
    do i = 1, size(ludum,1)
       ludum(i) = falloc()
    end do
    do i = 1, nthread
       thread_ti(i)%id = i
       do j = 1, size(thread_ti(i)%lu,1)
          thread_ti(i)%lu(j) = falloc()
          thread(i) = c_null_ptr
       end do
    end do
    do i = 1, nthread
       do j = 1, size(thread_ti(i)%lu,1)
          call fdealloc(thread_ti(i)%lu(j))
       end do
    end do
    do i = 1, size(ludum,1)
       call fdealloc(ludum(i))
    end do

    ! Initialize glfw
    fdum = glfwSetErrorCallback(c_funloc(error_callback))
    if (glfwInit() == 0) &
       call ferror('gui_start','Failed to initialize GLFW',faterr)
    glfw_lenient = .true.
    call glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, opengl_version_major)
    call glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, opengl_version_minor)
    call glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE)
    call glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE)
    glfw_lenient = .false.
    ! call glfwWindowHint(GLFW_SAMPLES, ms_samples) ! activate multisampling

    ! set up window
    strc = gui_title
    rootwin = glfwCreateWindow(1280, 720, c_loc(strc), c_null_ptr, c_null_ptr)
    if (.not.c_associated(rootwin)) &
       call ferror('gui_start','Failed to create window',faterr)
    fdum = glfwSetDropCallback(rootwin,c_funloc(drop_callback))
    call glfwMakeContextCurrent(rootwin)
    call glfwSwapInterval(1) ! enable vsync

    ! set GUI icon
    file = trim(critic_home) // dirsep // "assets" // dirsep // "critic2_icon.png" // c_null_char
    icon%pixels = stbi_load(c_loc(file), icon%width, icon%height, idum, 4)
    if (.not.c_associated(icon%pixels)) &
       call ferror('gui_start','Could not find GUI assets: have you set CRITIC_HOME?',faterr)
    glfw_lenient = .true.
    call glfwSetWindowIcon(rootwin, 1, c_loc(icon))
    glfw_lenient = .false.
    call stbi_image_free(icon%pixels)

    ! set up ImGui context
    ptrc = igCreateContext(c_null_ptr)
    if (.not.c_associated(ptrc)) &
       call ferror('gui_start','Failed to create ImGui context',faterr)
    ptrc = ipCreateContext()
    if (.not.c_associated(ptrc)) &
       call ferror('gui_start','Failed to create ImPlot context',faterr)

    ! initialize gl3w
    idum = gl3wInit()
    if (idum /= 0) &
       call ferror('gui_start','Failed to initialize OpenGL (gl3w)',faterr)
    if (gl3wIsSupported(opengl_version_major,opengl_version_minor) == 0) &
       call ferror('gui_start','gl3w: OpenGL version ' // &
       string(opengl_version_major) // '.' // string(opengl_version_minor) // ' not supported',faterr)

    ! set glfw options
    glfw_lenient = .true.
    call glfwSetInputMode(rootwin, GLFW_STICKY_KEYS, 1)
    glfw_lenient = .false.

    ! set opengl options
    call glEnable(GL_DEPTH_TEST)
    call glEnable(GL_CULL_FACE)
    call glDisable(GL_BLEND)
    call glEnable(GL_MULTISAMPLE)

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
    ptrc = io%Fonts
    call c_f_pointer(ptrc,fonts)
    io%MouseDrawCursor = .true. ! can't get anything other than the arrow otherwise
    io%configflags = ior(io%configflags,ImGuiConfigFlags_DockingEnable) ! activate docking
    io%configflags = ior(io%configflags,ImGuiConfigFlags_DpiEnableScaleFonts)
    io%configflags = ior(io%configflags,ImGuiConfigFlags_NavEnableKeyboard)
    io%inifilename = c_null_ptr ! do not save the ini file, for now

    ! default font, with Greek letters and letter-like symbols
    range = (/32_c_short,   255_c_short,& ! default (basic latin + supplement)
             880_c_short,  1023_c_short,& ! Greek and Coptic
            8488_c_short,  8527_c_short,& ! letter-like symbols
            8704_c_short,  8959_c_short,& ! mathematical operators
            9472_c_short,  9599_c_short,& ! box drawing
            9632_c_short,  9727_c_short,& ! geometric shapes
            9728_c_short,  9983_c_short,& ! miscellaneous symbols
            9984_c_short, 10175_c_short,& ! dingbats
               0_c_short/)
    font_normal = ImFontAtlas_AddFontFromMemoryCompressedBase85TTF(io%fonts,font_dejavu_base85_ptr,&
       fontbakesize,c_null_ptr,c_loc(range))
    font_large = ImFontAtlas_AddFontFromMemoryCompressedBase85TTF(io%fonts,font_dejavu_base85_ptr,&
       fontbakesize_large,c_null_ptr,c_loc(range))

    ! bold font if using freetype
#ifdef IMGUI_ENABLE_FREETYPE
    fonts%FontBuilderFlags = ImGuiFreeTypeBuilderFlags_Bold
#endif

    ! get the ImGui context pointer and the main viewport
    ptrc = igGetCurrentContext()
    call c_f_pointer(ptrc,g)
    ptrc = igGetMainViewport()
    call c_f_pointer(ptrc,mainvwp)

    ! set the initial ImGui style
    call igStyleColorsDark(c_null_ptr)
    g%Style%FrameRounding = 3._c_float
    g%Style%Colors(ImGuiCol_TabActive+1) = g%Style%Colors(ImGuiCol_TabHovered+1)

    ! Parse the command line and read as many systems as possible
    ! Initialize the first system, if available
    call process_arguments()
    call launch_initialization_thread()
    call system_shorten_names()

    ! set default UI settings
    call set_default_ui_settings()

    ! create buffers for objects and compile and link shaders
    call shapes_init()
    call shaders_init()
    call windows_init()

    ! initialize the window stack with the toggle-able windows
    iwin_tree = stack_create_window(wintype_tree,.true.,permanent=.true.)
    iwin_view = stack_create_window(wintype_view,.true.,purpose=wpurp_view_main,permanent=.true.)
    iwin_console_input = stack_create_window(wintype_console_input,.true.,permanent=.true.)
    iwin_console_output = stack_create_window(wintype_console_output,.true.,permanent=.true.)
    iwin_builder = stack_create_window(wintype_builder,.true.,permanent=.true.)
    iwin_about = stack_create_window(wintype_about,.false.,permanent=.true.)

    ! main loop
    show_demo_window = .false.
    show_implot_demo_window = .false.
    firstpass = .true.
    do while (glfwWindowShouldClose(rootwin) == 0)
       ! poll events
       call glfwPollEvents()

       ! start the ImGui frame
       call ImGui_ImplOpenGL3_NewFrame()
       call ImGui_ImplGlfw_NewFrame()
       call igNewFrame()

       ! calculate default font size
       strc = "A" // c_null_char
       call igCalcTextSize(fontsize,c_loc(strc),c_null_ptr,.false._c_bool,-1._c_float)

       ! set the transient flags to false
       sysc(1:nsys)%highlight_transient_set = .false.

       ! show main menu
       call show_main_menu()

       ! show main dockspace
       iddock = igDockSpaceOverViewport(igGetMainViewport(),&
          ior(ImGuiDockNodeFlags_PassthruCentralNode,ImGuiDockNodeFlags_AutoHideTabBar),&
          c_null_ptr)

       ! maybe reallocate the window stack
       call stack_realloc_maybe()

       ! process the window stack
       do i = 1, nwin
          call win(i)%draw()
       end do

       ! first pass: use the dock builder routines to place the windows
       ! https://github.com/ocornut/imgui/issues/2109
       if (firstpass) then
          ileft = igDockBuilderSplitNode(iddock, ImGuiDir_Left, 0.3_c_float, idum, iright)
          ! ibottom = igDockBuilderSplitNode(iright, ImGuiDir_Down, 0.3_c_float, idum, idum2)

          call igDockBuilderDockWindow(c_loc(win(iwin_tree)%name), ileft)
          call igDockBuilderDockWindow(c_loc(win(iwin_view)%name), iright)
          call igDockBuilderDockWindow(c_loc(win(iwin_console_input)%name), ileft)
          call igDockBuilderDockWindow(c_loc(win(iwin_console_output)%name), ileft)
          call igDockBuilderDockWindow(c_loc(win(iwin_builder)%name), ileft)
          call igDockBuilderFinish(ileft)
          call igDockBuilderFinish(iright)
          call igDockBuilderFinish(iddock)
          win(iwin_console_input)%isopen = .true.
          win(iwin_console_output)%isopen = .true.
          win(iwin_builder)%isopen = .true.
          call igSetWindowFocus_Str(c_loc(win(iwin_view)%name))
       end if

       ! show demo window
       if (show_demo_window) &
          call igShowDemoWindow(show_demo_window)

       ! show implot demo window
       if (show_implot_demo_window) &
          call ipShowDemoWindow(show_implot_demo_window)

       ! if there are commands to run from the input console, set up the modal
       if (force_run_commands > 0) then
          call win(iwin_console_input)%block_gui_ci(force_run_commands == 2)
       end if

       ! rendering
       call igRender()
       call glfwGetFramebufferSize(rootwin, display_w, display_h)
       call glViewport(0, 0, display_w, display_h)
       call glClearColor(0.45, 0.55, 0.60, 1.00)
       call glClear(GL_COLOR_BUFFER_BIT)
       call ImGui_ImplOpenGL3_RenderDrawData(igGetDrawData())

       ! swap buffers
       call glfwSwapBuffers(rootwin)

       ! run commands from the input console
       if (force_run_commands == 1) then
          call win(iwin_console_input)%run_commands_ci()
          force_run_commands = 0
       elseif (force_run_commands == 2) then
          saveinpcon = win(iwin_console_input)%inpcon_selected
          do i = 1, nsys
             win(iwin_console_input)%inpcon_selected = i
             call win(iwin_console_input)%run_commands_ci()
          end do
          win(iwin_console_input)%inpcon_selected = saveinpcon
          force_run_commands = 0
       end if

       ! if the transient flag is false, clear the transient highlight
       do i = 1, nsys
          if (.not.sysc(i)%highlight_transient_set) &
             call sysc(i)%highlight_clear(.true.)
       end do

       firstpass = .false.
    end do

    ! cleanup windows
    do i = 1, nwin
       call win(i)%end()
    end do

    ! cleanup
    call shapes_end()
    call shaders_end()
    call ImGui_ImplOpenGL3_Shutdown()
    call ImGui_ImplGlfw_Shutdown()
    call ipDestroyContext(c_null_ptr)
    call igDestroyContext(c_null_ptr)

    ! cleanup mutexes
#ifdef _THREADS
    do i = 1, nsys
       if (c_associated(sysc(i)%thread_lock)) call deallocate_mtx(sysc(i)%thread_lock)
    end do
#endif

    ! terminate
    glfw_lenient = .true.
    call glfwDestroyWindow(rootwin)
    call glfwTerminate()
    glfw_lenient = .false.

  contains
    ! typedef void(* GLFWerrorfun) (int, const char *)
    subroutine error_callback(error,description) bind(c)
      use c_interface_module, only: c_f_string_alloc, c_strlen
      use tools_io, only: ferror, faterr, string, warning
      integer(c_int), value :: error
      type(c_ptr), intent(in), value :: description

      character(len=:), allocatable :: msg

      call c_f_string_alloc(description,msg)
      if (glfw_lenient) then
         call ferror('glfw',"GLFW error (" // string(error) // "): " // trim(msg),warning)
      else
         call ferror('glfw',"GLFW error (" // string(error) // "): " // trim(msg),faterr)
      end if

    end subroutine error_callback
    ! void drop_callback(GLFWwindow* window, int count, const char* paths[])
    subroutine drop_callback(window,count,ipaths) bind(c)
      use systems, only: add_systems_from_name, launch_initialization_thread,&
         system_shorten_names
      use global, only: rborder_def
      use param, only: isformat_r_unknown
      use c_interface_module, only: c_f_string_alloc
      type(c_ptr), value :: window
      integer(c_int), value :: count
      type(c_ptr), intent(in) :: ipaths(count)

      integer :: i
      character(kind=c_char,len=:), allocatable :: file

      if (count < 1) return
      do i = 1, count
         call c_f_string_alloc(ipaths(i),file)
         call add_systems_from_name(file,-1,isformat_r_unknown,.false.,rborder_def,.false.)
      end do
      call launch_initialization_thread()
      call system_shorten_names()

    end subroutine drop_callback
  end subroutine gui_start

  !> Reset all user interface settings to their default values
  module subroutine set_default_ui_settings()
    use keybindings, only: set_default_keybindings
    use param, only: JMLcol

    ! interface settings
    io%FontGlobalScale = 1._c_float
    tooltip_enabled = .true.
    tooltip_delay = 0.5_c_float
    tooltip_wrap_factor = 40._c_float
    tree_select_updates_inpcon = .true.
    tree_select_updates_view = .true.

    ! key bindings
    call set_default_keybindings()

    ! tree colors
    ColorTableCellBg = ColorTableCellBg_def
    ColorHighlightScene = ColorHighlightScene_def
    ColorHighlightSelectScene = ColorHighlightSelectScene_def
    ColorMeasureSelect = ColorMeasureSelect_def
    ColorElement = real(JMLcol,c_float) / 255._c_float

  end subroutine set_default_ui_settings

  !xx! private procedures

  ! Process the command-line arguments. Skip the options and load the files.
  subroutine process_arguments()
    use systems, only: add_systems_from_name
    use global, only: rborder_def
    use param, only: isformat_r_unknown
    integer :: argc
    integer :: i
    character(len=1024) :: argv

    argc = command_argument_count()
    do i = 1, argc
       call getarg(i,argv)
       argv = adjustl(argv)
       if (argv(1:1) == "-") cycle ! skip options
       call add_systems_from_name(argv,-1,isformat_r_unknown,.false.,rborder_def,.false.)
    end do

  end subroutine process_arguments

  ! Show the main menu
  subroutine show_main_menu()
    use interfaces_cimgui
    use systems, only: sys, sysc, sys_init, ok_system, are_threads_running, duplicate_system,&
       reread_system_from_file, remove_system, kill_initialization_thread, write_system
    use windows, only: win, iwin_tree, iwin_view, iwin_console_input,&
       iwin_console_output, iwin_builder, iwin_about, stack_create_window, wintype_dialog,&
       wpurp_dialog_openfiles, wintype_new_struct, wintype_new_struct_library,&
       wintype_preferences, wintype_view, wpurp_view_alternate, wintype_load_field,&
       wintype_about, wintype_geometry, wintype_rebond, wintype_vibrations, wintype_exportimage
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_text, iw_calcwidth, iw_menuitem
    use keybindings, only: BIND_QUIT, BIND_OPEN, BIND_CLOSE, BIND_REOPEN, BIND_NEW,&
       BIND_GEOMETRY, BIND_SAVE, get_bind_keyname, is_bind_event
    use interfaces_glfw, only: GLFW_TRUE, glfwSetWindowShouldClose
    use tools_io, only: string
    use param, only: isformat_write_from_read, isformat_w_unknown

    ! enum for the dialog types that can be launched from the menu
    integer, parameter :: d_open = 1
    integer, parameter :: d_close = 2
    integer, parameter :: d_reopen = 3
    integer, parameter :: d_new = 4
    integer, parameter :: d_newlib = 5
    integer, parameter :: d_preferences = 6
    integer, parameter :: d_geometry = 7
    integer, parameter :: d_save = 8
    integer, parameter :: D_TOTAL = 8

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    integer(c_int) :: idum
    logical :: launchquit, launch(D_TOTAL), isysok, isysvok, ifieldok, ok
    integer :: isys, isysv

    logical, save :: ttshown = .false. ! tooltip flag

    ! keybindings
    !! menu key bindings
    launch(d_open) = is_bind_event(BIND_OPEN)
    launch(d_close) = is_bind_event(BIND_CLOSE)
    launch(d_reopen) = is_bind_event(BIND_REOPEN)
    launch(d_new) = is_bind_event(BIND_NEW)
    launch(d_newlib) = .false.
    launch(d_preferences) = .false.
    launch(d_geometry) = is_bind_event(BIND_GEOMETRY)
    launch(d_save) = is_bind_event(BIND_SAVE)
    launchquit = is_bind_event(BIND_QUIT)

    ! isys is the tree selected system, isysv is the view selected system
    isys = win(iwin_tree)%tree_selected
    isysok = ok_system(isys,sys_init)
    isysv = win(iwin_view)%view_selected
    isysvok = ok_system(isysv,sys_init)

    ! start the menu
    if (igBeginMainMenuBar()) then
       ! File
       str1 = "File" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          ! File -> New
          launch(d_new) = launch(d_new) .or. iw_menuitem("New...",BIND_NEW)
          call iw_tooltip("Create a new structure from scratch",ttshown)

          ! File -> New from library
          launch(d_newlib) = iw_menuitem("New from Library...")
          call iw_tooltip("Create a new structure from the critic2 library",ttshown)

          ! File -> Open
          launch(d_open) = launch(d_open) .or. iw_menuitem("Open...",BIND_OPEN)
          call iw_tooltip("Read molecular or crystal structures from external file(s)",ttshown)

         ! File -> Duplicate
          if (iw_menuitem("Duplicate",enabled=isysok)) &
             call duplicate_system(isys)
          call iw_tooltip("Create a copy of this system",ttshown)

          ! File -> Reopen from file
          launch(d_reopen) = launch(d_reopen) .or. iw_menuitem("Restore",BIND_REOPEN,enabled=isysok)
          call iw_tooltip("Restore the system to the original geometry it had when it was first opened",ttshown)

          ! File -> Close
          launch(d_close) = launch(d_close) .or. iw_menuitem("Close",BIND_CLOSE,enabled=isysok)
          call iw_tooltip("Close the current system",ttshown)

          ! File -> Separator
          call igSeparator()

          ! File -> Save
          ok = isysok .and. (isformat_write_from_read(sysc(isys)%seed%isformat) /= isformat_w_unknown)
          launch(d_save) = launch(d_save) .or. iw_menuitem("Save",BIND_SAVE,enabled=ok)
          if (.not.ok) launch(d_save) = .false.
          call iw_tooltip("Save the current system to the file from where it was read (only input files)",ttshown)

          ! File -> Separator
          call igSeparator()

          ! File -> Load Field
          if (iw_menuitem("Load Field...",enabled=isysok)) &
             idum = stack_create_window(wintype_load_field,.true.,isys=isys,orraise=-1)
          call iw_tooltip("Load a scalar field for the current system",ttshown)

          ! File -> Remove Field
          ifieldok = isysok
          if (ifieldok) ifieldok = (sys(isys)%iref /= 0)
          if (iw_menuitem("Remove Field",enabled=ifieldok)) &
             call sys(isys)%unload_field(sys(isys)%iref)
          call iw_tooltip("Remove the reference field from the current system",ttshown)

          ! File -> Separator
          call igSeparator()

          ! File -> Quit
          launchquit = launchquit .or. iw_menuitem("Quit",BIND_QUIT)
          call iw_tooltip("Quit critic2",ttshown)

          call igEndMenu()
       else
          ttshown = .false.
       end if

       ! Edit
       str1 = "Edit" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          ! Edit -> Preferences...
          if (iw_menuitem("Preferences...")) &
             idum = stack_create_window(wintype_preferences,.true.,orraise=-1)
          call iw_tooltip("Change the user interface settings and key bindings",ttshown)

          call igEndMenu()
       else
          ttshown = .false.
       end if

       ! Windows
       str1 = "Windows" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          ! Windows -> Tree
          if (iw_menuitem("Tree",selected=logical(win(iwin_tree)%isopen))) &
             win(iwin_tree)%isopen = .not.win(iwin_tree)%isopen
          call iw_tooltip("Toggle the display of the tree window",ttshown)

          ! Windows -> Main View
          if (iw_menuitem("Main View",selected=logical(win(iwin_view)%isopen))) &
             win(iwin_view)%isopen = .not.win(iwin_view)%isopen
          call iw_tooltip("Toggle the display of the main view window",ttshown)

          ! Windows -> Input Console
          if (iw_menuitem("Input Console",selected=logical(win(iwin_console_input)%isopen))) &
             win(iwin_console_input)%isopen = .not.win(iwin_console_input)%isopen
          call iw_tooltip("Toggle the display of the input console window",ttshown)

          ! Windows -> Output Console
          if (iw_menuitem("Output Console",selected=logical(win(iwin_console_output)%isopen))) &
             win(iwin_console_output)%isopen = .not.win(iwin_console_output)%isopen
          call iw_tooltip("Toggle the display of the output console window",ttshown)

          ! Windows -> Input Console
          if (iw_menuitem("Builder",selected=logical(win(iwin_builder)%isopen))) &
             win(iwin_builder)%isopen = .not.win(iwin_builder)%isopen
          call iw_tooltip("Toggle the display of the builder window",ttshown)

          ! Windows -> Separator
          call igSeparator()

          ! Windows -> Alternate view
          if (iw_menuitem("New View Window")) &
             idum = stack_create_window(wintype_view,.true.,purpose=wpurp_view_alternate)
          call iw_tooltip("Open a new view window in addition to the main view",ttshown)

          call igEndMenu()
       else
          ttshown = .false.
       end if

       ! Windows
       str1 = "Tools" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          if (iw_menuitem("Export to Image...",enabled=isysvok)) &
             idum = stack_create_window(wintype_exportimage,.true.,idparent=iwin_view,orraise=-1)
          call iw_tooltip("Export the current view to an image file (png)",ttshown)

          ! separator
          call igSeparator()

          launch(d_geometry) = launch(d_geometry) .or. &
             iw_menuitem("View/Edit Geometry...",BIND_GEOMETRY,enabled=isysvok)
          call iw_tooltip("View and edit the atomic positions, bonds, etc.",ttshown)

          if (iw_menuitem("Recalculate Bonds...",enabled=isysvok.and..not.are_threads_running())) &
             idum = stack_create_window(wintype_rebond,.true.,isys=isysv,orraise=-1)
          call iw_tooltip("Recalculate the bonds in the current system",ttshown)

          if (iw_menuitem("Vibrations...",enabled=isysvok)) &
             idum = stack_create_window(wintype_vibrations,.true.,idparent=iwin_view,orraise=-1)
          call iw_tooltip("Display an animation showing the atomic vibrations for this system",ttshown)

          call igEndMenu()
       else
          ttshown = .false.
       end if

       ! Help
       str1 = "Help" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          ! Help -> Critic2 Manual
          if (iw_menuitem("Critic2 Manual...")) then
             str2 = "https://aoterodelaroza.github.io/critic2/" // c_null_char
             call openLink(c_loc(str2))
          end if
          call iw_tooltip("Visit the critic2 website for more information about the program",ttshown)

          ! Help -> About
          if (iw_menuitem("About...")) then
             if (win(iwin_about)%isopen) then
                call igSetWindowFocus_Str(c_loc(win(iwin_about)%name))
             else
                call win(iwin_about)%init(wintype_about,.true.,iwin_about)
             end if
          end if
          call igEndMenu()
       else
          ttshown = .false.
       end if

       ! fps message
       call igSetCursorPosX(iw_calcwidth(30,0,from_end=.true.))
       call iw_text(string(1000._c_float / io%Framerate,'f',decimal=3) // " ms/frame (" // &
          string(io%Framerate,'f',decimal=1) // " FPS)")
    end if
    call igEndMainMenuBar()

    ! process launches
    if (launch(d_new)) &
       idum = stack_create_window(wintype_new_struct,.true.,orraise=-1)
    if (launch(d_newlib)) &
       idum = stack_create_window(wintype_new_struct_library,.true.,orraise=-1)
    if (launch(d_open)) &
       idum = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openfiles,orraise=-1)
    if (isysok) then
       if (launch(d_reopen)) &
          call reread_system_from_file(isys)
       if (launch(d_close)) &
          call remove_system(isys)
       if (launch(d_save)) &
          call write_system(isys)
    end if

    if (launch(d_geometry).and.isysvok) then
       idum = stack_create_window(wintype_geometry,.true.,isys=isysv,orraise=-1)
    end if
    if (launchquit) then
       if (are_threads_running()) &
          call kill_initialization_thread()
       call glfwSetWindowShouldClose(rootwin, GLFW_TRUE)
    end if

  end subroutine show_main_menu

end submodule proc
