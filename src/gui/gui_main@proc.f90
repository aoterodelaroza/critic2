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

  ! the dockspace ID
  integer(c_int) :: iddock = 0

  ! threads in execution
  integer :: nthread = 1
  type(c_ptr), target, allocatable :: thread(:)

  !xx! private procedures
  ! subroutine process_arguments()
  ! function stack_create_window(type,isopen)
  ! subroutine show_main_menu()
  ! subroutine process_global_keybindings()
  ! function initialization_thread_worker(arg)

contains

  ! Start the critic2 GUI.
  module subroutine gui_start()
    use gui_interfaces_threads
    use gui_interfaces_cimgui
    use gui_interfaces_glfw
    use gui_interfaces_opengl3
    use gui_window, only: wintype_tree, wintype_view, wintype_console
    use gui_keybindings, only: set_default_keybindings
    use c_interface_module, only: f_c_string_dup, C_string_free
    use tools_io, only: ferror, faterr, string
    use omp_lib, only: omp_get_max_threads
    integer(c_int) :: idum, idum2, display_w, display_h, ileft, iright, ibottom
    type(c_funptr) :: fdum
    type(c_ptr) :: ptrc, pdum
    logical(c_bool) :: ldum, show_demo_window
    character(kind=c_char,len=:), allocatable, target :: strc
    integer :: i
    logical :: firstpass
    integer(c_short), allocatable, target :: range(:)

    ! initialize the sys arrays
    nsys = 0
    allocate(sys(1),sysc(1))

    ! initialize threads
     nthread = omp_get_max_threads()
    if (allocated(thread)) deallocate(thread)
    allocate(thread(nthread))

    ! Parse the command line and read as many systems as possible
    ! Initialize the first system, if available
    call process_arguments()
    call launch_initialization_thread()

    ! Initialize glfw
    fdum = glfwSetErrorCallback(c_funloc(error_callback))
    if (glfwInit() == 0) &
       call ferror('gui_start','Failed to initialize GLFW',faterr)
    call glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, opengl_version_major)
    call glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, opengl_version_minor)
    call glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE)
    call glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1)
    ! call glfwWindowHint(GLFW_SAMPLES, 4) ! activate multisampling

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
    io%configflags = ior(io%configflags,ImGuiConfigFlags_DpiEnableScaleFonts)
    io%inifilename = c_null_ptr ! do not save the ini file, for now

    ! default font, 16 pt and with Greek letters and letter-like symbols
    range = (/32_c_short,  255_c_short,& ! default (basic latin + supplement)
             880_c_short, 1023_c_short,& ! Greek and Coptic
            8488_c_short, 8527_c_short,& ! letter-like symbols
            9472_c_short, 9599_c_short,& ! box drawing
            9632_c_short, 9727_c_short,& ! geometric shapes
               0_c_short/)
    pdum = ImFontAtlas_AddFontFromMemoryCompressedBase85TTF(io%fonts,myfont_ttf_compressed_data_base85_ptr,&
       16._c_float,c_null_ptr,c_loc(range))

    ! get the ImGUI context pointer
    ptrc = igGetCurrentContext()
    call c_f_pointer(ptrc,g)

    ! set default keybindings
    call set_default_keybindings()

    ! initialize the window stack with the toggle-able windows (open, for now)
    iwin_tree = stack_create_window(wintype_tree,.true.)
    iwin_view = stack_create_window(wintype_view,.true.)
    iwin_console = stack_create_window(wintype_console,.true.)

    ! main loop
    show_demo_window = .true.
    firstpass = .true.
    do while (glfwWindowShouldClose(rootwin) == 0)
       ! poll events
       call glfwPollEvents()
       time = glfwGetTime()

       ! start the ImGui frame
       call ImGui_ImplOpenGL3_NewFrame()
       call ImGui_ImplGlfw_NewFrame()
       call igNewFrame()

       ! handle quit key binding
       call process_global_keybindings()

       ! show main menu
       call show_main_menu()

       ! show main dockspace
       iddock = igDockSpaceOverViewport(igGetMainViewport(),&
          ImGuiDockNodeFlags_PassthruCentralNode,&
          c_null_ptr)

       ! process the window stack
       do i = 1, nwin
          call win(i)%draw()
       end do

       ! first pass: use the dock builder routines to place the windows
       ! https://github.com/ocornut/imgui/issues/2109
       if (firstpass) then
          ileft = igDockBuilderSplitNode(iddock, ImGuiDir_Left, 0.25_c_float, idum, iright)
          ibottom = igDockBuilderSplitNode(iright, ImGuiDir_Down, 0.3_c_float, idum, idum2)

          call igDockBuilderDockWindow(c_loc(win(iwin_tree)%name), ileft)
          call igDockBuilderDockWindow(c_loc(win(iwin_view)%name), iright)
          call igDockBuilderDockWindow(c_loc(win(iwin_console)%name), ibottom)
          call igDockBuilderFinish(iddock);
       end if

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
       firstpass = .false.
    end do

    ! cleanup
    call ImGui_ImplOpenGL3_Shutdown()
    call ImGui_ImplGlfw_Shutdown()
    call igDestroyContext(c_null_ptr)

    ! cleanup mutexes
    do i = 1, nsys
       if (c_associated(sysc(i)%thread_lock)) call deallocate_mtx(sysc(i)%thread_lock)
    end do

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

  !> Launch the initialization threads, which will go over all systems
  !> trying to initialize them.
  module subroutine launch_initialization_thread()
    use gui_interfaces_threads

    integer :: i, idum

    do i = 1, nthread
       idum = thrd_create(c_loc(thread(i)), c_funloc(initialization_thread_worker), c_null_ptr)
    end do

  end subroutine launch_initialization_thread

  !xx! private procedures

  ! Process the command-line arguments. Skip the options and load the files.
  subroutine process_arguments()
    use grid1mod, only: grid1_register_ae
    use gui_interfaces_threads, only: allocate_mtx, mtx_init, mtx_plain
    use crystalseedmod, only: read_seeds_from_file, crystalseed
    use tools_io, only: uout
    use types, only: realloc
    use param, only: dirsep
    integer :: argc
    integer :: i, j
    character(len=1024) :: argv
    integer :: nseed, iseed, iseed_, iafield, idx
    type(crystalseed), allocatable :: seed(:)
    character(len=:), allocatable :: errmsg
    type(system), allocatable :: syaux(:)
    type(sysconf), allocatable :: syscaux(:)
    integer(c_int) :: idum
    character(len=:), allocatable :: str
    logical :: collapse

    argc = command_argument_count()
    do i = 1, argc
       call getarg(i,argv)
       argv = adjustl(argv)
       if (argv(1:1) == "-") cycle ! skip options
       call read_seeds_from_file(argv,-1,nseed,seed,collapse,errmsg,iafield)

       if (nseed > 0) then
          ! increment and reallocate if necessary
          nsys = nsys + nseed
          if (nsys > size(sys,1)) then
             allocate(syaux(2*nsys))
             syaux(1:size(sys,1)) = sys
             call move_alloc(syaux,sys)

             allocate(syscaux(2*nsys))
             syscaux(1:size(sysc,1)) = sysc
             call move_alloc(syscaux,sysc)
          end if

          do iseed = 1, nseed
             ! re-order if collapse to have the final structure be first, flag them
             if (collapse) then
                iseed_ = mod(iseed,nseed)+1
             else
                iseed_ = iseed
             end if

             ! make a new system for this seed
             idx = nsys - nseed + iseed_
             sys(idx)%isinit = .false.
             sysc(idx)%id = idx
             sysc(idx)%seed = seed(iseed)
             sysc(idx)%has_field = .false.

             ! initialization status
             if (collapse.and.iseed_ == 1) then
                ! master
                sysc(idx)%collapse = -1
                sysc(idx)%status = sys_loaded_not_init
             elseif (collapse) then
                ! dependent
                sysc(idx)%collapse = nsys - nseed + 1
                sysc(idx)%status = sys_loaded_not_init_hidden
             else
                ! independent
                sysc(idx)%collapse = 0
                sysc(idx)%status = sys_loaded_not_init
             end if

             ! initialize the mutex
             if (.not.c_associated(sysc(idx)%thread_lock)) then
                sysc(idx)%thread_lock = allocate_mtx()
                idum = mtx_init(sysc(idx)%thread_lock,mtx_plain)
             end if

             ! register all all-electron densities (global - should not be done by threads)
             do j = 1, seed(iseed)%nspc
                call grid1_register_ae(seed(iseed)%spc(j)%z)
             end do
          end do
          if (iafield > 0) &
             sysc(nsys - nseed + iafield)%has_field = .true.
       else
          write (uout,'("!! Warning !! Could not read structures from: ",A)') trim(argv)
          write (uout,'("Error: ",A)') trim(errmsg)
       end if
    end do

    ! remove all common path strings up to the last dirsep
    if (nsys > 0) then
       main: do while (.true.)
          ! grab string up to the first dirsep
          idx = index(sysc(1)%seed%name,dirsep,.true.)
          if (idx == 0) exit main
          str = sysc(1)%seed%name(1:idx)

          ! check all names start with the same string
          do i = 2, nsys
             if (len_trim(sysc(i)%seed%name) < idx) exit main
             if (sysc(i)%seed%name(1:idx) /= str) exit main
          end do

          ! remove the string
          do i = 1, nsys
             sysc(i)%seed%name = sysc(i)%seed%name(idx+1:)
          end do
       end do main
    end if

  end subroutine process_arguments

  !> Create a window in the window stack with the given type
  function stack_create_window(type,isopen)
    use gui_window, only: window
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    type(window), allocatable :: aux(:)
    integer :: stack_create_window

    nwin = nwin + 1
    if (.not.allocated(win)) then
       allocate(win(nwin))
    elseif (nwin > size(win,1)) then
       allocate(aux(2*nwin))
       aux(1:size(win,1)) = win
       call move_alloc(aux,win)
    end if
    call win(nwin)%init(type,isopen)
    stack_create_window = nwin

  end function stack_create_window

  ! Show the main menu
  subroutine show_main_menu()
    use gui_interfaces_cimgui
    use gui_utils, only: igIsItemHovered_delayed
    use gui_keybindings, only: BIND_QUIT, get_bind_keyname
    use gui_interfaces_glfw, only: GLFW_TRUE, glfwSetWindowShouldClose
    use tools_io, only: string

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: v2
    logical, save :: ttshown(2) = (/.false.,.false./) ! menu-level tooltips

    if (igBeginMainMenuBar()) then
       ! File
       str1 = "File" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then

          ! File -> Quit
          str1 = "Quit" // c_null_char
          str2 = get_bind_keyname(BIND_QUIT) // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_loc(str2),.false._c_bool,.true._c_bool)) &
             call glfwSetWindowShouldClose(rootwin, GLFW_TRUE)
          if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown(1))) then
             str1 = "Quit the program" // c_null_char
             call igSetTooltip(c_loc(str1))
             ttshown(1) = .true.
          end if

          call igEndMenu()
       else
          ttshown(1) = .false.
       end if

       ! Windows
       str1 = "Windows" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          ! File -> Tree
          str1 = "Tree" // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,win(iwin_tree)%isopen,.true._c_bool)) &
             win(iwin_tree)%isopen = .not.win(iwin_tree)%isopen
          if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown(2))) then
             str1 = "Toggle the Tree window" // c_null_char
             call igSetTooltip(c_loc(str1))
             ttshown(2) = .true.
          end if

          ! File -> View
          str1 = "View" // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,win(iwin_view)%isopen,.true._c_bool)) &
             win(iwin_view)%isopen = .not.win(iwin_view)%isopen
          if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown(2))) then
             str1 = "Toggle the View window" // c_null_char
             call igSetTooltip(c_loc(str1))
             ttshown(2) = .true.
          end if

          ! File -> Console
          str1 = "Console" // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,win(iwin_console)%isopen,.true._c_bool)) &
             win(iwin_console)%isopen = .not.win(iwin_console)%isopen
          if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown(2))) then
             str1 = "Toggle the Console window" // c_null_char
             call igSetTooltip(c_loc(str1))
             ttshown(2) = .true.
          end if
          call igEndMenu()
       else
          ttshown(2) = .false.
       end if

       ! fps message
       call igGetContentRegionAvail(v2)
       call igSameLine(0._c_float, v2%x - 180._c_float)
       str1 = string(1000._c_float / io%Framerate,'f',decimal=3) // " ms/frame (" // &
          string(io%Framerate,'f',decimal=1) // " FPS)" // c_null_char
       call igText(c_loc(str1))
    end if
    call igEndMainMenuBar()

  end subroutine show_main_menu

  ! Process the global keybindings
  subroutine process_global_keybindings()
    use gui_interfaces_cimgui
    use gui_keybindings, only: is_bind_event, BIND_QUIT, BIND_CLOSE_POPUP
    use gui_interfaces_glfw, only: glfwSetWindowShouldClose,GLFW_TRUE

    ! quit
    if (is_bind_event(BIND_QUIT)) &
       call glfwSetWindowShouldClose(rootwin, GLFW_TRUE)

    ! close last popup
    if (is_bind_event(BIND_CLOSE_POPUP)) &
       call igClosePopupsExceptModals()

  end subroutine process_global_keybindings

  ! Thread worker: run over all systems and initialize the ones that are not locked
  function initialization_thread_worker(arg)
    use gui_interfaces_threads
    use tools_io, only: string, uout
    type(c_ptr), value :: arg
    integer(c_int) :: initialization_thread_worker

    integer :: i, i0, i1
    integer(c_int) :: idum
    integer :: iff
    character(len=:), allocatable :: errmsg
    integer, pointer :: idx

    i0 = 1
    i1 = nsys
    if (c_associated(arg)) then
       call c_f_pointer(arg,idx)
       i0 = idx
       i1 = idx
    end if

    do i = i0, i1
       if (c_associated(sysc(i)%thread_lock)) then
          idum = mtx_trylock(sysc(i)%thread_lock)
          if (idum == thrd_success) then
             if (sysc(i)%status == sys_loaded_not_init) then
                sysc(i)%status = sys_initializing
                ! load the seed
                call sys(i)%new_from_seed(sysc(i)%seed)

                ! load any fields
                if (sysc(i)%has_field) then
                   call sys(i)%load_field_string(sysc(i)%seed%file,.false.,iff,errmsg)
                   if (len_trim(errmsg) > 0) then
                      write (uout,'("!! Warning !! Could not read field for system: ",A)') string(i)
                   else
                      sys(i)%f(iff)%file = sysc(i)%seed%file
                      sys(i)%f(iff)%name = sysc(i)%seed%name
                   end if
                end if

                ! this system has been initialized
                sysc(i)%status = sys_init

                ! force resize and sort of table columns (no lock needed for this)
                if (iwin_tree > 0 .and. iwin_tree <= nwin) then
                   win(iwin_tree)%forceresize = .true.
                   win(iwin_tree)%forcesort = .true.
                end if
             end if

             ! unlock
             idum = mtx_unlock(sysc(i)%thread_lock)
          end if
       end if
    end do
    initialization_thread_worker = 0

  end function initialization_thread_worker

end submodule proc
