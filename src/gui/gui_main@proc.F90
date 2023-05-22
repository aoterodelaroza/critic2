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

  ! threads in execution
  integer, parameter :: nthread = 1
  type(c_ptr), target, allocatable :: thread(:)
  type(thread_info), target, allocatable :: thread_ti(:)

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
    use shaders, only: shaders_init, shaders_end
    use shapes, only: shapes_init, shapes_end
    use windows, only: nwin, win, wintype_tree, wintype_view, wintype_console_input,&
       wintype_console_output, wintype_about, iwin_tree, iwin_view, iwin_console_input,&
       iwin_console_output, iwin_about, stack_create_window, stack_realloc_maybe,&
       wpurp_view_main
    use global, only: critic_home
    use c_interface_module, only: f_c_string_dup, C_string_free
    use tools_io, only: ferror, faterr, string, falloc, fdealloc
    use param, only: dirsep
    use omp_lib, only: omp_get_max_threads
    integer(c_int) :: idum, idum2, display_w, display_h, ileft, iright, ibottom, ileft2, iright2
    type(c_funptr) :: fdum
    type(c_ptr) :: ptrc
    logical(c_bool) :: ldum, show_demo_window, show_implot_demo_window
    character(kind=c_char,len=:), allocatable, target :: strc, file
    integer :: i, j, ludum(10), saveinpcon
    logical :: firstpass
    integer(c_short), allocatable, target :: range(:)
    type(GLFWimage), target :: icon

    ! initialize the sys arrays
    nsys = 0
    allocate(sys(1),sysc(1))

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

    ! Parse the command line and read as many systems as possible
    ! Initialize the first system, if available
    call process_arguments()
    call launch_initialization_thread()
    call system_shorten_names()

    ! Initialize glfw
    fdum = glfwSetErrorCallback(c_funloc(error_callback))
    if (glfwInit() == 0) &
       call ferror('gui_start','Failed to initialize GLFW',faterr)
    call glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, opengl_version_major)
    call glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, opengl_version_minor)
    call glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE)
    call glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE)
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
    call glfwSetWindowIcon(rootwin, 1, c_loc(icon))
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
    call glfwSetInputMode(rootwin, GLFW_STICKY_KEYS, 1)

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

    ! get the ImGui context pointer and the main viewport
    ptrc = igGetCurrentContext()
    call c_f_pointer(ptrc,g)
    ptrc = igGetMainViewport()
    call c_f_pointer(ptrc,mainvwp)

    ! set the initial ImGui style
    call igStyleColorsDark(c_null_ptr)
    g%Style%FrameRounding = 3._c_float

    ! set default UI settings
    call set_default_ui_settings()

    ! create buffers for objects and compile and link shaders
    call shapes_init()
    call shaders_init()

    ! initialize the window stack with the toggle-able windows (open, for now)
    iwin_tree = stack_create_window(wintype_tree,.true.)
    iwin_view = stack_create_window(wintype_view,.true.,purpose=wpurp_view_main)
    iwin_console_input = stack_create_window(wintype_console_input,.true.)
    iwin_console_output = stack_create_window(wintype_console_output,.true.)
    iwin_about = stack_create_window(wintype_about,.false.)

    ! main loop
    show_demo_window = .false.
    show_implot_demo_window = .false.
    firstpass = .true.
    do while (glfwWindowShouldClose(rootwin) == 0)
       ! poll events
       call glfwPollEvents()
       time = glfwGetTime()

       ! start the ImGui frame
       call ImGui_ImplOpenGL3_NewFrame()
       call ImGui_ImplGlfw_NewFrame()
       call igNewFrame()

       ! calculate default font size
       strc = "A" // c_null_char
       call igCalcTextSize(fontsize,c_loc(strc),c_null_ptr,.false._c_bool,-1._c_float)

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
          ibottom = igDockBuilderSplitNode(iright, ImGuiDir_Down, 0.3_c_float, idum, idum2)
          ileft2 = igDockBuilderSplitNode(ibottom, ImGuiDir_Left, 0.4_c_float, idum, iright2)

          call igDockBuilderDockWindow(c_loc(win(iwin_tree)%name), ileft)
          call igDockBuilderDockWindow(c_loc(win(iwin_view)%name), iright)
          call igDockBuilderDockWindow(c_loc(win(iwin_console_input)%name), ileft2)
          call igDockBuilderDockWindow(c_loc(win(iwin_console_output)%name), iright2)
          call igDockBuilderFinish(iddock)
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
    call glfwDestroyWindow(rootwin)
    call glfwTerminate()

  contains
    ! typedef void(* GLFWerrorfun) (int, const char *)
    subroutine error_callback(error,description) bind(c)
      use c_interface_module, only: c_f_string_alloc, c_strlen
      use tools_io, only: ferror, faterr, string
      integer(c_int), value :: error
      type(c_ptr), intent(in), value :: description

      character(len=:), allocatable :: msg

      call c_f_string_alloc(description,msg)
      call ferror('glfw',"GLFW error (" // string(error) // "): " // trim(msg),faterr)

    end subroutine error_callback
    ! void drop_callback(GLFWwindow* window, int count, const char* paths[])
    subroutine drop_callback(window,count,ipaths) bind(c)
      use global, only: rborder_def
      use param, only: isformat_unknown
      use c_interface_module, only: c_f_string_alloc
      type(c_ptr), value :: window
      integer(c_int), value :: count
      type(c_ptr), intent(in) :: ipaths(count)

      integer :: i
      character(kind=c_char,len=:), allocatable :: file

      if (count < 1) return
      do i = 1, count
         call c_f_string_alloc(ipaths(i),file)
         call add_systems_from_name(file,-1,isformat_unknown,.false.,rborder_def,.false.)
      end do
      call launch_initialization_thread()
      call system_shorten_names()

    end subroutine drop_callback
  end subroutine gui_start

  !> Launch the initialization threads, which will go over all systems
  !> trying to initialize them.
  module subroutine launch_initialization_thread()
    use interfaces_threads, only: thrd_create
    integer :: i, idum

    force_quit_threads = .false.
#ifdef _THREADS
    do i = 1, nthread
       idum = thrd_create(c_loc(thread(i)), c_funloc(initialization_thread_worker), c_loc(thread_ti(i)))
    end do
#else
    idum = initialization_thread_worker(c_loc(thread_ti(1)))
#endif

  end subroutine launch_initialization_thread

  !> Force the initialization threads to quit and wait for them to
  !> finish before returning.
  module subroutine kill_initialization_thread()
    use interfaces_threads, only: wrap_thrd_join
    integer :: i
    integer(c_int) :: res, idum

#ifdef _THREADS
    force_quit_threads = .true.
    do i = 1, nthread
       idum = wrap_thrd_join(c_loc(thread(i)),res)
    end do
    force_quit_threads = .false.
#endif

  end subroutine kill_initialization_thread

  ! Returns .true. if there are any initialization threads running
  module function are_threads_running()
    logical :: are_threads_running
    are_threads_running = any(thread_ti(1:nthread)%active)
  end function are_threads_running

  ! Re-write the seed names from the full-path names to remove as much
  ! of the prefixes as possible
  module subroutine system_shorten_names()
    use param, only: dirsep

    integer :: idx, i, i1
    character(len=:), allocatable :: str, removed

    ! reset all names to the full-path name
    do i = 1, nsys
       if (sysc(i)%status /= sys_empty .and..not.sysc(i)%renamed) &
          sysc(i)%seed%name = sysc(i)%fullname
    end do

    ! remove as much as possible from the beginning
    if (nsys > 0) then
       i1 = 0
       do i = 1, nsys
          if (sysc(i)%status /= sys_empty.and..not.sysc(i)%renamed) then
             i1 = i
             exit
          end if
       end do
       if (i1 > 0) then
          idx = len_trim(sysc(i1)%seed%name)+1
          removed = ""
          main: do while (.true.)
             ! grab string up to the first dirsep
             idx = index(sysc(i1)%seed%name(1:idx-1),dirsep,.true.)
             if (idx == 0) exit main
             str = sysc(i1)%seed%name(1:idx)

             ! check all names start with the same string
             do i = i1+1, nsys
                if (sysc(i)%status == sys_empty.or.sysc(i)%renamed) cycle
                if (len_trim(sysc(i)%seed%name) < idx) exit
                if (sysc(i)%seed%name(1:idx) /= str) cycle main
             end do

             ! remove the string
             removed = removed // sysc(i1)%seed%name(1:idx)
             do i = 1, nsys
                if (sysc(i)%status == sys_empty.or.sysc(i)%renamed) cycle
                sysc(i)%seed%name = sysc(i)%seed%name(idx+1:)
             end do
          end do main
       end if
    end if

  end subroutine system_shorten_names

  !> Add systems by reading them from a file, passed by name. mol = 1
  !> read as crystal, 0 read as molecule, -1 autodetect. isformat,
  !> force file format if /= 0.
  module subroutine add_systems_from_name(name,mol,isformat,readlastonly,rborder,molcubic)
    use crystalseedmod, only: crystalseed, read_seeds_from_file
    use tools_io, only: uout
    character(len=*), intent(in) :: name
    integer, intent(in) :: mol
    integer, intent(in) :: isformat
    logical, intent(in) :: readlastonly
    real*8, intent(in) :: rborder
    logical, intent(in) :: molcubic

    integer :: i, nseed
    type(crystalseed), allocatable :: seed(:)
    integer :: iafield
    logical :: collapse
    character(len=:), allocatable :: errmsg

    ! read all seeds from the file
    errmsg = ""
    nseed = 0
    allocate(seed(1))
    call read_seeds_from_file(name,mol,isformat,readlastonly,nseed,seed,collapse,errmsg,iafield)
    if (len_trim(errmsg) > 0 .or. nseed == 0) goto 999

    ! set the border and molcubic
    do i = 1, nseed
       seed(i)%border = rborder
       seed(i)%cubic = molcubic
    end do

    ! add the systems
    call add_systems_from_seeds(nseed,seed,collapse,iafield)

    return
999 continue
    write (uout,'("WARNING : Could not read structures from: ",A)') trim(name)
    write (uout,'("WARNING : ",A/)') trim(errmsg)

  end subroutine add_systems_from_name

  !> Add systems from the given seeds. If collapse is present and
  !> true, reorder to make the last system be first and collapse the
  !> other seeds in the tree view. If iafield, load a field from that
  !> seed.
  module subroutine add_systems_from_seeds(nseed,seed,collapse,iafield)
    use interfaces_cimgui, only: getCurrentWorkDir
    use grid1mod, only: grid1_register_ae
    use gui_main, only: reuse_mid_empty_systems
    use windows, only: nwin, win, iwin_tree
    use interfaces_threads, only: allocate_mtx, mtx_init, mtx_plain
    use crystalseedmod, only: read_seeds_from_file, crystalseed
    use tools_io, only: uout
    use types, only: realloc
    use param, only: dirsep
    integer, intent(in) :: nseed
    type(crystalseed), allocatable, intent(in) :: seed(:)
    logical, intent(in), optional :: collapse
    integer, intent(in), optional :: iafield

    integer :: i, j, nid
    integer :: iafield_
    integer :: iseed, iseed_, idx, in
    character(len=:), allocatable :: errmsg, str
    character(kind=c_char,len=:), allocatable, target :: strc
    type(system), allocatable :: syaux(:)
    type(sysconf), allocatable :: syscaux(:)
    integer(c_int) :: idum
    logical :: collapse_
    integer, allocatable :: id(:)

    if (nseed == 0) then
       write (uout,'("WARNING : Could not read structures.")')
       return
    end if

    ! initialize
    errmsg = ""
    collapse_ = .false.
    if (present(collapse)) collapse_ = collapse
    iafield_ = 0
    if (present(iafield)) iafield_ = iafield

    ! find contiguous IDs for the new systems
    allocate(id(nseed))
    if (reuse_mid_empty_systems) then
       ! may reuse old IDs
       nid = 0
       do i = 1, nsys
          if (sysc(i)%status == sys_empty) then
             nid = nid + 1
             id(nid) = i
             if (nid == nseed) exit
          else
             nid = 0
          end if
       end do
       do i = nid+1, nseed
          id(i) = nsys + i - nid
       end do
    else
       ! append IDs at the end, discard empty systems
       nid = 0
       do i = nsys, 1, -1
          if (sysc(i)%status /= sys_empty) then
             nid = i
             exit
          end if
       end do
       do i = 1, nseed
          id(i) = nid + i
       end do
    end if

    ! increment and reallocate if necessary
    nsys = max(nsys,id(nseed))
    if (nsys > size(sys,1)) then
       allocate(syaux(2*nsys))
       syaux(1:size(sys,1)) = sys
       call move_alloc(syaux,sys)

       allocate(syscaux(2*nsys))
       syscaux(1:size(sysc,1)) = sysc
       call move_alloc(syscaux,sysc)

       ! refresh pointers in the systems after the move_alloc
       do i = 1, nsys
          do j = 1, sys(i)%nf
             sys(i)%f(j)%sptr = c_loc(sys(i))
          end do
       end do
    end if

    do iseed = 1, nseed
       ! re-order if collapse to have the final structure be first, flag them
       if (collapse_) then
          iseed_ = mod(iseed,nseed)+1
       else
          iseed_ = iseed
       end if

       ! make a new system for this seed
       idx = id(iseed_)
       sys(idx)%isinit = .false.
       sysc(idx)%id = idx
       sysc(idx)%seed = seed(iseed)
       sysc(idx)%has_field = .false.
       sysc(idx)%renamed = .false.
       sysc(idx)%showfields = .false.
       sysc(idx)%idwin_plotscf = 0

       ! write down the full name
       str = trim(adjustl(sysc(idx)%seed%name))
       if (str(1:1) == dirsep) then
          sysc(idx)%fullname = str
       else
          if (allocated(strc)) deallocate(strc)
          allocate(character(len=2049) :: strc)
          idum = getCurrentWorkDir(c_loc(strc),2048_c_size_t)
          in = index(strc,c_null_char)-1
          if (strc(in:in) == dirsep.and.in > 0) in = in - 1
          sysc(idx)%fullname = strc(1:in) // dirsep // str
       end if

       ! initialization status
       sysc(idx)%status = sys_loaded_not_init
       if (collapse_.and.iseed_ == 1) then
          ! master
          sysc(idx)%collapse = -1
          sysc(idx)%hidden = .false.
       elseif (collapse_) then
          ! dependent
          sysc(idx)%collapse = id(1)
          sysc(idx)%hidden = .true.
       else
          ! independent
          sysc(idx)%collapse = 0
          sysc(idx)%hidden = .false.
       end if

       ! initialize the mutex
#ifdef _THREADS
       if (.not.c_associated(sysc(idx)%thread_lock)) then
          sysc(idx)%thread_lock = allocate_mtx()
          idum = mtx_init(sysc(idx)%thread_lock,mtx_plain)
       end if
#endif

       ! register all all-electron densities (global - should not
       ! be done by threads because agrid and cgrid are common)
       do j = 1, seed(iseed)%nspc
          call grid1_register_ae(seed(iseed)%spc(j)%z)
       end do

       ! set the iafield
       if (iseed == iafield_) sysc(idx)%has_field = .true.

       ! set the time
       sysc(idx)%timelastchange = time
    end do
    deallocate(id)

  end subroutine add_systems_from_seeds

  ! Remove system with index idx and leave behind a sys_empty spot. If
  ! master and collapsed, kill all dependents. If master and extended,
  ! make all dependents master.
  recursive module subroutine remove_system(idx)
    use interfaces_threads, only: deallocate_mtx
    integer, intent(in) :: idx

    integer :: i

    if (idx < 1 .or. idx > nsys) return
    if (sysc(idx)%status == sys_empty) return
    call sys(idx)%end()
    call sysc(idx)%seed%end()
#ifdef _THREADS
    if (c_associated(sysc(idx)%thread_lock)) then
       call deallocate_mtx(sysc(idx)%thread_lock)
       sysc(idx)%thread_lock = c_null_ptr
    end if
#endif
    sysc(idx)%status = sys_empty
    sysc(idx)%hidden = .false.
    sysc(idx)%showfields = .false.
    sysc(idx)%timelastchange = time

    if (sysc(idx)%collapse == -1) then
       ! kill all dependents if collapsed
       do i = 1, nsys
          if (sysc(i)%status /= sys_empty .and. sysc(i)%collapse == idx) &
             call remove_system(i)
       end do
    else
       ! make all dependents their own master if extended
       do i = 1, nsys
          if (sysc(i)%status /= sys_empty .and. sysc(i)%collapse == idx) &
             sysc(i)%collapse = 0
       end do
    end if

  end subroutine remove_system

  !> Reset all user interface settings to their default values
  module subroutine set_default_ui_settings()
    use keybindings, only: set_default_keybindings

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
    ColorTableCellBg_Mol     = ImVec4(0.43,0.8 ,0.  ,0.2)
    ColorTableCellBg_MolClus = ImVec4(0.0 ,0.8 ,0.43,0.2)
    ColorTableCellBg_MolCrys = ImVec4(0.8 ,0.43,0.0 ,0.2)
    ColorTableCellBg_Crys3d  = ImVec4(0.8 ,0.  ,0.0 ,0.2)
    ColorTableCellBg_Crys2d  = ImVec4(0.8 ,0.  ,0.43,0.2)
    ColorTableCellBg_Crys1d  = ImVec4(0.8 ,0.43,0.43,0.2)

  end subroutine set_default_ui_settings

  !xx! private procedures

  ! Process the command-line arguments. Skip the options and load the files.
  subroutine process_arguments()
    use global, only: rborder_def
    use param, only: isformat_unknown
    integer :: argc
    integer :: i
    character(len=1024) :: argv

    argc = command_argument_count()
    do i = 1, argc
       call getarg(i,argv)
       argv = adjustl(argv)
       if (argv(1:1) == "-") cycle ! skip options
       call add_systems_from_name(argv,-1,isformat_unknown,.false.,rborder_def,.false.)
    end do

  end subroutine process_arguments

  ! Show the main menu
  subroutine show_main_menu()
    use interfaces_cimgui
    use windows, only: win, iwin_tree, iwin_view, iwin_console_input,&
       iwin_console_output, iwin_about, stack_create_window, wintype_dialog,&
       wpurp_dialog_openfiles, wintype_new_struct, wintype_new_struct_library,&
       wintype_preferences, update_window_id, wintype_view, wpurp_view_alternate
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_text, iw_calcwidth
    use keybindings, only: BIND_QUIT, BIND_OPEN, BIND_NEW, get_bind_keyname, is_bind_event
    use interfaces_glfw, only: GLFW_TRUE, glfwSetWindowShouldClose
    use tools_io, only: string

    ! enum for the dialog types that can be launched from the menu
    integer, parameter :: d_open = 1
    integer, parameter :: d_new = 2
    integer, parameter :: d_newlib = 3
    integer, parameter :: d_preferences = 4

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    integer(c_int) :: idum
    logical(c_bool) :: enabled(4)
    logical :: launchquit, launch(4)
    integer :: i

    logical, save :: ttshown = .false. ! tooltip flag
    integer, save :: id(4) = (/0,0,0,0/) ! the ID for the dialogs

    ! check if the open and new dialogs are still open
    do i = 1, 4
       call update_window_id(id(i))
    end do

    ! calculate enabled and launches from keybindings
    do i = 1, 4
       enabled(i) = (id(i) == 0)
    end do

    ! keybindings
    !! menu key bindings
    launch(d_open) = (enabled(d_open) .and. is_bind_event(BIND_OPEN))
    launch(d_new) = (enabled(d_new) .and. is_bind_event(BIND_NEW))
    launch(d_newlib) = .false.
    launch(d_preferences) = .false.
    launchquit = is_bind_event(BIND_QUIT)

    ! start the menu
    if (igBeginMainMenuBar()) then
       ! File
       str1 = "File" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then

          ! File -> New
          str1 = "New..." // c_null_char
          str2 = trim(get_bind_keyname(BIND_NEW)) // c_null_char
          launch(d_new) = launch(d_new) .or. igMenuItem_Bool(c_loc(str1),c_loc(str2),.false._c_bool,enabled(d_new))
          call iw_tooltip("Create a new structure",ttshown)

          ! File -> New from library
          str1 = "New from Library..." // c_null_char
          launch(d_newlib) = igMenuItem_Bool(c_loc(str1),c_null_ptr,.false._c_bool,enabled(d_newlib))
          call iw_tooltip("Create a new structure from the critic2 library",ttshown)

          ! File -> Open
          str1 = "Open..." // c_null_char
          str2 = trim(get_bind_keyname(BIND_OPEN)) // c_null_char
          launch(d_open) = launch(d_open) .or. igMenuItem_Bool(c_loc(str1),c_loc(str2),.false._c_bool,enabled(d_open))
          call iw_tooltip("Read molecule or crystal structures from file(s)",ttshown)

          ! File -> Separator
          call igSeparator()

          ! File -> Quit
          str1 = "Quit" // c_null_char
          str2 = trim(get_bind_keyname(BIND_QUIT)) // c_null_char
          launchquit = launchquit .or. igMenuItem_Bool(c_loc(str1),c_loc(str2),.false._c_bool,.true._c_bool)
          call iw_tooltip("Quit the program",ttshown)

          call igEndMenu()
       else
          ttshown = .false.
       end if

       ! Edit
       str1 = "Edit" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          ! Edit -> Preferences...
          str1 = "Preferences..." // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,.false._c_bool,.true._c_bool)) then
             if (id(d_preferences) == 0) then
                id(d_preferences) = stack_create_window(wintype_preferences,.true.)
             else
                call igSetWindowFocus_Str(c_loc(win(id(d_preferences))%name))
             end if
          end if
          call iw_tooltip("Change the global UI settings and key bindings",ttshown)

          call igEndMenu()
       else
          ttshown = .false.
       end if

       ! Windows
       str1 = "Windows" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          ! Windows -> Tree
          str1 = "Tree" // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,win(iwin_tree)%isopen,.true._c_bool)) &
             win(iwin_tree)%isopen = .not.win(iwin_tree)%isopen
          call iw_tooltip("Toggle the tree window",ttshown)

          ! Windows -> Main View
          str1 = "Main View" // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,win(iwin_view)%isopen,.true._c_bool)) &
             win(iwin_view)%isopen = .not.win(iwin_view)%isopen
          call iw_tooltip("Toggle the main view window",ttshown)

          ! Windows -> Input Console
          str1 = "Input Console" // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,win(iwin_console_input)%isopen,.true._c_bool)) &
             win(iwin_console_input)%isopen = .not.win(iwin_console_input)%isopen
          call iw_tooltip("Toggle the input console window",ttshown)

          ! Windows -> Output Console
          str1 = "Output Console" // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,win(iwin_console_output)%isopen,.true._c_bool)) &
             win(iwin_console_output)%isopen = .not.win(iwin_console_output)%isopen
          call iw_tooltip("Toggle the output console window",ttshown)

          ! Windows -> Separator
          call igSeparator()

          ! Windows -> Alternate view
          str1 = "New Alternate View" // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,.false._c_bool,.true._c_bool)) &
             idum = stack_create_window(wintype_view,.true.,purpose=wpurp_view_alternate)
          call iw_tooltip("Open a new view window in addition to the main view",ttshown)

          call igEndMenu()
       else
          ttshown = .false.
       end if

       ! Help
       str1 = "Help" // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          ! Help -> Critic2 Manual
          str1 = "Critic2 Manual..." // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,.false._c_bool,.true._c_bool)) then
             str2 = "https://aoterodelaroza.github.io/critic2/" // c_null_char
             call openLink(c_loc(str2))
          end if
          call iw_tooltip("Visit the critic2 website for more information about the program",ttshown)

          ! Help -> About
          str1 = "About..." // c_null_char
          if (igMenuItem_Bool(c_loc(str1),c_null_ptr,.false._c_bool,.true._c_bool)) then
             if (win(iwin_about)%isopen) then
                call igSetWindowFocus_Str(c_loc(win(iwin_about)%name))
             else
                win(iwin_about)%isopen = .true.
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
       id(d_new) = stack_create_window(wintype_new_struct,.true.)
    if (launch(d_newlib)) &
       id(d_newlib) = stack_create_window(wintype_new_struct_library,.true.)
    if (launch(d_open)) &
       id(d_open) = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openfiles)
    if (launchquit) then
       if (are_threads_running()) &
          call kill_initialization_thread()
       call glfwSetWindowShouldClose(rootwin, GLFW_TRUE)
    end if

  end subroutine show_main_menu

  ! Thread worker: run over all systems and initialize the ones that are not locked
  function initialization_thread_worker(arg)
    use interfaces_threads, only: thrd_success, mtx_unlock, mtx_trylock
    use tools_io, only: string, uout
    type(c_ptr), value :: arg
    integer(c_int) :: initialization_thread_worker
    type(thread_info), pointer :: ti

    integer :: i, i0, i1
    integer(c_int) :: idum
    integer :: iff
    character(len=:), allocatable :: errmsg

    ! recover the thread info pointer
    call c_f_pointer(arg,ti)
    ti%active = .true.

    ! run over systems
    i0 = 1
    i1 = nsys
    do i = i0, i1
       ! check if they want us to quit
       if (force_quit_threads) exit
       ! try to grab the lock for this system
#ifdef _THREADS
       if (c_associated(sysc(i)%thread_lock)) then
          idum = mtx_trylock(sysc(i)%thread_lock)
          if (idum == thrd_success) then
#endif
             ! see if we can load it - only uninitialized and not hidden
             if (sysc(i)%status == sys_loaded_not_init.and..not.sysc(i)%hidden) then
                sysc(i)%status = sys_initializing
                ! load the seed
                call sys(i)%new_from_seed(sysc(i)%seed,ti=ti)

                ! load any fields
                if (sysc(i)%has_field) then
                   call sys(i)%load_field_string(sysc(i)%seed%file,.false.,iff,errmsg,ti=ti)
                   if (len_trim(errmsg) > 0) then
                      write (uout,'("!! Warning !! Could not read field for system: ",A)') string(i)
                      write (uout,'("!! Warning !! Error message: ",A)') trim(errmsg)
                   end if
                end if

                ! initialize the scene and build the initial draw list
                call sysc(i)%sc%init(i)
                call sysc(i)%sc%build_lists()

                ! copy the camera if this is a dependent
                if (sysc(i)%sc%lockedcam /= 0 .and. sysc(i)%sc%lockedcam /= i) &
                   call sysc(i)%sc%copy_cam(idx=sysc(i)%sc%lockedcam)

                ! this system has been initialized
                sysc(i)%status = sys_init
                sysc(i)%timelastchange = time
             end if

#ifdef _THREADS
             ! unlock
             idum = mtx_unlock(sysc(i)%thread_lock)
          end if
       end if
#endif
    end do
    initialization_thread_worker = 0
    ti%active = .false.

  end function initialization_thread_worker

end submodule proc
