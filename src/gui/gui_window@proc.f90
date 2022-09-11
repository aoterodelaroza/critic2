! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Some utilities for building the GUI (e.g. wrappers around ImGui routines).
submodule (gui_window) proc
  use gui_interfaces_cimgui
  use types, only: vstring
  implicit none

  ! Count unique IDs for keeping track of windows and widget
  integer :: idcount = 0

  ! column ids for the table in the tree widget
  integer(c_int), parameter :: ic_closebutton = 0
  integer(c_int), parameter :: ic_expandbutton = 1
  integer(c_int), parameter :: ic_id = 2
  integer(c_int), parameter :: ic_name = 3
  integer(c_int), parameter :: ic_spg = 4
  integer(c_int), parameter :: ic_v = 5
  integer(c_int), parameter :: ic_vmol = 6
  integer(c_int), parameter :: ic_nneq = 7
  integer(c_int), parameter :: ic_ncel = 8
  integer(c_int), parameter :: ic_nmol = 9
  integer(c_int), parameter :: ic_a = 10
  integer(c_int), parameter :: ic_b = 11
  integer(c_int), parameter :: ic_c = 12
  integer(c_int), parameter :: ic_alpha = 13
  integer(c_int), parameter :: ic_beta = 14
  integer(c_int), parameter :: ic_gamma = 15
  integer(c_int), parameter :: ic_e = 16
  integer(c_int), parameter :: ic_emol = 17
  integer(c_int), parameter :: ic_p = 18

  ! the buffer for the output console
  character(kind=c_char,len=:), allocatable, target :: outputb
  integer(c_size_t) :: lob = 0
  integer(c_size_t), parameter :: maxlob = 2000000

  ! the buffer for the input console
  character(kind=c_char,len=:), allocatable, target :: inputb
  integer(c_size_t), parameter :: maxlib = 40000

  ! command input and output
  integer, parameter :: command_inout_empty = 0
  integer, parameter :: command_inout_used = 1
  type command_inout
     integer :: id ! unique command ID
     integer :: status = command_inout_empty ! status of this command
     integer(c_size_t) :: size = 0 ! size of the output
     real(c_float) :: scrolly = 0._c_float ! scrolling position in console output view
     character(len=:,kind=c_char), allocatable :: tooltipinfo ! tooltip info for command (c_null_char term)
     character(len=:,kind=c_char), allocatable :: input ! command input (c_null_char term)
     character(len=:,kind=c_char), allocatable :: output ! command output (c_null_char term)
   contains
     procedure :: end => command_end
  end type command_inout

  ! command inout stack
  integer :: ncomid = 0 ! unique command ID generator (always incremented)
  integer :: ncom = 0 ! number of commands (com(:))
  integer :: nicom = 0 ! number of ordered commands (icom(:))
  integer :: idcom = 0 ! current output shown (0 = all)
  integer, allocatable :: icom(:) ! order in which to show the commands
  type(command_inout), allocatable, target :: com(:) ! the actual commands
  integer(c_size_t), parameter :: maxcomout = 20000000 ! maximum command size

  !xx! private procedures
  ! function tree_system_tooltip_string(i)
  ! function tree_field_tooltip_string(si,fj)
  ! subroutine get_library_structure_list(libfile,nst,st,ismol)
  ! subroutine draw_spg_table(ispg)

contains

  !xx! Methods for the command_inout type

  !> Deallocate all data in a command and reset it to emtpy
  subroutine command_end(c)
    class(command_inout), intent(inout) :: c

    c%id = 0
    c%status = command_inout_empty
    c%size = 0
    c%scrolly = 0._c_float
    if (allocated(c%tooltipinfo)) deallocate(c%tooltipinfo)
    if (allocated(c%input)) deallocate(c%input)
    if (allocated(c%output)) deallocate(c%output)

  end subroutine command_end

  !xx! Module functions and subroutines

  !> Reallocate the window stack if there are only a few windows
  !> left. This must be done *outside* any window method because it
  !> will change the location of all window pointers.
  module subroutine stack_realloc_maybe()
    type(window), allocatable :: winaux(:)

    integer, parameter :: iroom = 5

    if (nwin+iroom > size(win,1)) then
       allocate(winaux(2*(nwin+iroom)))
       winaux(1:size(win,1)) = win
       call move_alloc(winaux,win)
    end if

  end subroutine stack_realloc_maybe

  !> Create a window in the window stack with the given type. Returns
  !> the window ID.
  module function stack_create_window(type,isopen,purpose,isys)
    use tools_io, only: ferror, faterr
    use gui_window, only: window, nwin, win
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in), optional :: purpose
    integer, intent(in), optional :: isys

    integer :: stack_create_window

    integer :: i, id
    integer, parameter :: maxwin = 40

    ! find the first unused window or create a new one
    id = 0
    do i = 1, nwin
       if (.not.win(i)%isinit) then
          id = i
          exit
       end if
    end do
    if (id == 0) then
       nwin = nwin + 1
       id = nwin
    end if

    ! reallocate if necessary
    if (.not.allocated(win)) allocate(win(maxwin))

    if (nwin > size(win,1)) &
       call ferror('stack_create_window','too many windows',faterr)

    ! initialize the new window
    call win(id)%init(type,isopen,purpose,isys)
    stack_create_window = id

  end function stack_create_window

  !> Check whether the window with the given id is still open. If it is
  !> not, or id points to an invalid window, set it to id = 0. If changed
  !> is present, set it to the old id if the id has been changed.
  module subroutine update_window_id(id,changed)
    integer, intent(inout) :: id
    integer, intent(out), optional :: changed

    integer :: oid

    oid = id
    id = 0
    if (present(changed)) changed = oid
    if (oid < 1 .or. oid > nwin) return
    if (.not.win(oid)%isinit .or. .not.win(oid)%isopen) return
    id = oid
    if (present(changed)) changed = 0

  end subroutine update_window_id

  !> Initialize a window of the given type. If isiopen, initialize it
  !> as open.
  module subroutine window_init(w,type,isopen,purpose,isys)
    use gui_main, only: ColorDialogDir, ColorDialogFile
    use tools_io, only: ferror, faterr
    use param, only: bohrtoa
    class(window), intent(inout) :: w
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in), optional :: purpose
    integer, intent(in), optional :: isys

    character(kind=c_char,len=:), allocatable, target :: str1

    ! initialization of the state
    w%isinit = .true.
    w%firstpass = .true.
    w%isopen = isopen
    w%type = type
    w%id = -1
    w%name = ""
    w%table_selected = 1
    w%table_sortcid = 0
    w%table_sortdir = 1
    w%forceresize = .false.
    w%forcesort = .false.
    w%forceupdate = .false.
    w%forceinit = .false.
    w%inpcon_selected = 1
    w%okfile_set = .false. ! whether the library file has been set by the user
    w%okfile_read = .false. ! whether the structure list should be re-read from the lib
    if (allocated(w%iord)) deallocate(w%iord)
    w%dialog_data%dptr = c_null_ptr
    w%dialog_data%mol = -1
    w%dialog_data%showhidden = .false._c_bool
    w%dialog_data%isformat = isformat_unknown
    w%dialog_data%readlastonly = .false._c_bool
    w%dialog_data%purpose = wpurp_unknown
    w%dialog_data%molcubic = logical(.false.,c_bool)
    w%dialog_data%rborder = real(rborder_def*bohrtoa,c_float)

    ! type-specific initialization
    if (type == wintype_dialog) then
       ! create the dialog object and set up the style
       w%dptr = IGFD_Create()
       str1 = "+" // c_null_char
       call IGFD_SetFileStyle(w%dptr,IGFD_FileStyleByTypeDir,c_null_ptr,ColorDialogDir,c_loc(str1),c_null_ptr)
       str1 = " " // c_null_char
       call IGFD_SetFileStyle(w%dptr,IGFD_FileStyleByTypeFile,c_null_ptr,ColorDialogFile,c_loc(str1),c_null_ptr)
       if (.not.present(purpose)) &
          call ferror('window_init','dialog requires a purpose',faterr)
       w%dialog_purpose = purpose
    elseif (type == wintype_load_field) then
       if (.not.present(isys)) &
          call ferror('window_init','load_field requires isys',faterr)
       w%loadfield_isys = isys
    elseif (type == wintype_scfplot) then
       if (.not.present(isys)) &
          call ferror('window_init','scfplot requires isys',faterr)
       w%scfplot_isys = isys
    end if

  end subroutine window_init

  !> End a window and deallocate the data.
  module subroutine window_end(w)
    class(window), intent(inout) :: w

    ! window-specific destruction
    if (w%isinit .and. w%type == wintype_dialog .and. c_associated(w%dptr)) &
       call IGFD_Destroy(w%dptr)

    ! deallocate the rest of the data
    w%isinit = .false.
    w%isopen = .false.
    w%id = -1
    w%name = ""
    if (allocated(w%iord)) deallocate(w%iord)
    if (allocated(w%forceremove)) deallocate(w%forceremove)

  end subroutine window_end

  !> Return true if the root of this window is focused
  module function window_focused(w)
    use gui_main, only: g
    class(window), intent(inout) :: w
    logical :: window_focused

    type(ImGuiWindow), pointer :: wptr, navwin, navwinroot

    window_focused = .false.
    if (.not.c_associated(w%ptr).or..not.c_associated(g%NavWindow)) return
    call c_f_pointer(w%ptr,wptr)
    call c_f_pointer(g%NavWindow,navwin)
    if (.not.c_associated(navwin%RootWindow)) return
    call c_f_pointer(navwin%RootWindow,navwinroot)
    window_focused = associated(wptr,navwinroot)

  end function window_focused

  !> Draw an ImGui window.
  module subroutine window_draw(w)
    use gui_main, only: fontsize
    use gui_utils, only: iw_text
    use tools_io, only: string, ferror, faterr
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: inisize

    if (.not.w%isinit) return
    if (.not.w%isopen) return

    ! First pass on creation: assign ID, name, and flags
    if (w%id < 0) then
       w%firstpass = .true.
       idcount = idcount + 1
       w%id = idcount
       if (w%type == wintype_tree) then
          w%name = "Tree" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_view) then
          w%name = "View" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_console_input) then
          w%name = "Input Console" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_console_output) then
          w%name = "Output Console" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_dialog) then
          w%dialog_data%dptr = w%dptr
          w%dialog_data%purpose = w%dialog_purpose
          if (w%dialog_data%showhidden) then
             w%flags = ImGuiFileDialogFlags_None
          else
             w%flags = ImGuiFileDialogFlags_DontShowHiddenFiles
          end if
          str2 = "" // c_null_char ! default path
          inisize%x = 90 * fontsize%x
          inisize%y = 30 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)

          str1 = "All files (*.*){*.*}" // c_null_char
          if (w%dialog_purpose == wpurp_dialog_openfiles) then
             ! open dialog
             w%name = "Open File(s)..." // c_null_char
             str1 = &
                "&
                &All files (*.*){*.*},&
                &ABINIT (DEN...){(DEN|ELF|POT|VHA\VHXC|VXC|GDEN1|GDEN2|GDEN3|LDEN|KDEN|PAWDEN|VCLMB|VPSP)},&
                &ADF (molden){.molden},&
                &cif (cif){.cif},&
                &CRYSTAL (out){.out},&
                &cube (cube|bincube){.cube,.bincube},&
                &DFTB+ (gen){.gen},&
                &DMACRYS (dmain|16|21){.dmain,.16,.21},&
                &elk (OUT){.OUT},&
                &FHIaims (in|in.next_step|out|own){.in,.next_step,.out,.own},&
                &Gaussian (log|wfn|wfx|fchk|cube){.log,.wfn,.wfx,.fchk,.cube},&
                &ORCA (molden|molden.input){.molden,.input},&
                &postg (pgout){.pgout},&
                &psi4 (molden|dat){.molden,.dat},&
                &Quantum ESPRESSO (out|in|cube|pwc) {.out,.in,.cube,.pwc},&
                &SHELX (res|ins){.res,.ins.16},&
                &SIESTA (STRUCT_IN|STRUCT_OUT) {.STRUCT_IN,.STRUCT_OUT},&
                &TINKER (frac) {.frac},&
                &VASP (POSCAR|CONTCAR|...){(CONTCAR|CHGCAR|ELFCAR|CHG|AECCAR0|AECCAR1|AECCAR2|POSCAR)},&
                &WIEN2k (struct){.struct},&
                &Xcrysden (xsf|axsf) {.xsf,.axsf},&
                &xyz (xyz){.xyz},&
                &"// c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,0_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_savelogfile) then
             w%name = "Save Log File..." // c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile) then
             w%name = "Open Library File..." // c_null_char
             w%flags = ior(w%flags,ImGuiFileDialogFlags_Modal)
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openfieldfile) then
             w%name = "Open Field File(s)..." // c_null_char
             w%flags = ior(w%flags,ImGuiFileDialogFlags_Modal)
             str1 = &
                "&
                &All files (*.*){*.*},&
                &ABINIT (DEN...){(DEN|ELF|POT|VHA\VHXC|VXC|GDEN1|GDEN2|GDEN3|LDEN|KDEN|PAWDEN|VCLMB|VPSP)},&
                &Aimpac (qub){.qub},&
                &Binary cube file (bincube){.bincube},&
                &Cube file (cube){.cube},&
                &DFTB+ (detailed.xml){.xml},&
                &elk (grid){.grid},&
                &elk STATE.OUT (OUT){.OUT},&
                &Gaussian wavefunction (wfn|wfx|fchk){.wfn,.wfx,.fchk},&
                &Molden-style file (molden){.molden},&
                &Quantum ESPRESSO pwc file (pwc){.pwc},&
                &SIESTA (RHO...) {(RHO|BADER|DRHO|LDOS|VT|VH)},&
                &VASP ELFCAR{(ELFCAR)},&
                &VASP (CHGCAR...){(CONTCAR|CHGCAR|ELFCAR|CHG|AECCAR0|AECCAR1|AECCAR2|POSCAR)},&
                &WIEN2k (clmsum...){.clmsum,.clmup,.clmdn},&
                &Xcrysden (xsf|axsf) {.xsf,.axsf},&
                &"// c_null_char
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openonefilemodal) then
             w%name = "Open File..." // c_null_char
             w%flags = ior(w%flags,ImGuiFileDialogFlags_Modal)
             call IGFD_OpenPaneDialog2(w%dptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          else
             call ferror('window_draw','unknown dialog purpose',faterr)
          end if
       elseif (w%type == wintype_new_struct) then
          w%name = "New Structure..." // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 90 * fontsize%x
          inisize%y = 40 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_new_struct_library) then
          w%name = "New Structure from Library..." // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 90 * fontsize%x
          inisize%y = 30 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       elseif (w%type == wintype_load_field) then
          w%name = "Load Field..." // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 60 * fontsize%x
          inisize%y = 15 * fontsize%y
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       end if
    elseif (w%type == wintype_scfplot) then
       w%name = "SCF Iterations##" // string(w%scfplot_isys) // c_null_char
       w%flags = ImGuiWindowFlags_None
       inisize%x = 45 * fontsize%x
       inisize%y = inisize%x
       call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
    end if

    if (w%isopen) then
       ! assign the pointer ID for the window, if not a dialog
       if (w%type == wintype_dialog) then
          w%ptr = IGFD_GetCurrentWindow(w%dptr)
          call w%draw_dialog()
       else
          if (igBegin(c_loc(w%name),w%isopen,w%flags)) then
             ! draw the window contents, depending on type
             w%ptr = igGetCurrentWindow()
             if (w%type == wintype_tree) then
                call w%draw_tree()
             elseif (w%type == wintype_view) then
                call iw_text("Hello View!")
             elseif (w%type == wintype_console_input) then
                call w%draw_ci()
             elseif (w%type == wintype_console_output) then
                call w%draw_co()
             elseif (w%type == wintype_new_struct) then
                call w%draw_new_struct()
             elseif (w%type == wintype_new_struct_library) then
                call w%draw_new_struct_from_library()
             elseif (w%type == wintype_load_field) then
                call w%draw_load_field()
             elseif (w%type == wintype_scfplot) then
                call w%draw_scfplot()
             end if
          end if
          call igEnd()
       end if
       w%firstpass = .false.
    else
       w%firstpass = .true.
    end if

  end subroutine window_draw

  !> Draw the contents of a tree window
  module subroutine draw_tree(w)
    use gui_keybindings, only: is_bind_event, BIND_TREE_REMOVE_SYSTEM_FIELD
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button,&
       iw_text, iw_setposx_fromend, iw_calcwidth, iw_calcheight
    use gui_main, only: nsys, sys, sysc, sys_empty, sys_init,&
       sys_loaded_not_init, sys_initializing, ColorTableCellBg_Mol,&
       ColorTableCellBg_MolClus, ColorTableCellBg_MolCrys, ColorTableCellBg_Crys3d,&
       ColorTableCellBg_Crys2d, ColorTableCellBg_Crys1d, launch_initialization_thread,&
       kill_initialization_thread, system_shorten_names, remove_system, tooltip_delay,&
       ColorDangerButton, ColorFieldSelected, g, tree_select_updates_inpcon
    use fieldmod, only: type_grid
    use tools_io, only: string, uout
    use types, only: realloc
    use param, only: bohrtoa, ifformat_as_lap, ifformat_as_grad, ifformat_as_pot,&
       ifformat_as_resample
    use c_interface_module
    class(window), intent(inout), target :: w

    character(kind=c_char,len=1024), target :: txtinp
    character(kind=c_char,len=:), allocatable, target :: str, strpop, strpop2, zeroc, ch
    type(ImVec2) :: szero, sz
    integer(c_int) :: flags, color, idir
    integer :: i, j, k, nshown, newsel, jsel, ll, id, iref
    logical(c_bool) :: ldum, isel
    type(c_ptr) :: ptrc
    type(ImGuiTableSortSpecs), pointer :: sortspecs
    type(ImGuiTableColumnSortSpecs), pointer :: colspecs
    logical :: hadenabledcolumn, buttonhovered_close, buttonhovered_expand, reinit, isend, ok, found
    logical :: export
    real(c_float) :: width, pos

    type(c_ptr), save :: cfilter = c_null_ptr ! filter object (allocated first pass, never destroyed)
    logical, save :: ttshown = .false. ! tooltip flag
    integer(c_int), save :: iresample(3) = (/0,0,0/) ! for the grid resampling menu option
    integer(c_int), save :: idloadfield = 0 ! ID of the window used to load a field into the sytsem

    ! initialize
    hadenabledcolumn = .false.
    zeroc = "" // c_null_char
    szero%x = 0
    szero%y = 0
    if (.not.allocated(w%iord)) then
       w%table_sortcid = ic_id
       w%table_sortdir = 1
       w%table_selected = 1
       w%forceupdate = .true.
    end if

    ! update the window ID for the load field dialog
    call update_window_id(idloadfield)

    ! text filter
    if (.not.c_associated(cfilter)) &
       cfilter = ImGuiTextFilter_ImGuiTextFilter(c_loc(zeroc))
    str = "##treefilter" // c_null_char
    ldum = ImGuiTextFilter_Draw(cfilter,c_loc(str),0._c_float)
    call iw_tooltip("Filter systems by name in the list below. Use comma-separated fields&
       & and - for excluding. Example: inc1,inc2,-exc includes all systems&
       & with inc1 or inc2 and excludes systems with exc.",ttshown)
    if (iw_button("Clear",sameline=.true.)) then
       if (c_associated(cfilter)) &
          call ImGuiTextFilter_Clear(cfilter)
    end if
    call iw_tooltip("Clear the filter",ttshown)

    ! row of buttons
    ! button: expand
    if (iw_button("Expand")) then
       do i = 1, nsys
          call expand_system(i)
       end do
    end if
    call iw_tooltip("Expand all systems in the tree",ttshown)

    ! button: collapse
    if (iw_button("Collapse",sameline=.true.)) then
       do i = 1, nsys
          call collapse_system(i)
       end do
    end if
    call iw_tooltip("Collapse all systems in the tree",ttshown)

    ! button: export
    export = iw_button("Export",sameline=.true.)
    call iw_tooltip("Write a the current table to the output console in csv-style (for copying)",ttshown)

    ! helper
    call iw_text("(?)",sameline=.true.)
    call iw_tooltip("Right-click on the table headers for more options")

    ! right-align for the rest of the contents
    call iw_setposx_fromend(14,2)

    ! button: close
    if (iw_button("Close",danger=.true.)) then
       if (allocated(w%forceremove)) deallocate(w%forceremove)
       allocate(w%forceremove(nsys))
       k = 0
       do i = 1, nsys
          if (c_associated(cfilter)) then
             str = trim(sysc(i)%seed%name) // c_null_char
             if (.not.ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) cycle
          end if
          k = k + 1
          w%forceremove(k) = i
       end do
       call realloc(w%forceremove,k)
       if (c_associated(cfilter)) &
          call ImGuiTextFilter_Clear(cfilter)
    end if
    call iw_tooltip("Close all visible systems",ttshown)

    ! button: close all
    if (iw_button("Close All",danger=.true.,sameline=.true.)) then
       if (allocated(w%forceremove)) deallocate(w%forceremove)
       allocate(w%forceremove(nsys))
       do i = 1, nsys
          w%forceremove(i) = i
       end do
    end if
    call iw_tooltip("Close all systems",ttshown)

    ! process force options
    if (allocated(w%forceremove)) then
       ! stop all initialization threads, if any of the systems to
       ! remove may be initialized in the near future
       if (any((sysc(w%forceremove)%status == sys_loaded_not_init.and..not.sysc(w%forceremove)%hidden).or.&
                sysc(w%forceremove)%status == sys_initializing)) then
          call kill_initialization_thread()
          reinit = .true.
       else
          reinit = .false.
       end if

       ! remove a system and move the table selection if the system was selected
       do k = 1, size(w%forceremove,1)
          call remove_system(w%forceremove(k))
          ! if we removed the selected system, go to the next; either, go to the previous
          if (w%forceremove(k) == w%table_selected) then
             jsel = 0
             do j = 1, size(w%iord,1)
                if (w%iord(j) == w%table_selected) then
                   jsel = j
                   exit
                end if
             end do
             if (jsel > 0) then
                newsel = 0
                do j = jsel, size(w%iord,1)
                   i = w%iord(j)
                   if (sysc(i)%status /= sys_init) cycle
                   if (c_associated(cfilter)) then
                      str = trim(sysc(i)%seed%name) // c_null_char
                      if (.not.ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) cycle
                   end if
                   newsel = i
                   exit
                end do
                if (newsel == 0) then
                   do j = jsel, 1, -1
                      i = w%iord(j)
                      if (sysc(i)%status /= sys_init) cycle
                      if (c_associated(cfilter)) then
                         str = trim(sysc(i)%seed%name) // c_null_char
                         if (.not.ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) cycle
                      end if
                      newsel = i
                      exit
                   end do
                   if (newsel == 0) newsel = 1
                end if
             else
                newsel = 1
             end if
             w%table_selected = newsel
          end if
          ! if we removed the system for the input console, update
          if (w%forceremove(k) == win(iwin_console_input)%inpcon_selected) &
             win(iwin_console_input)%inpcon_selected = w%table_selected
       end do
       deallocate(w%forceremove)
       ! restart initialization if the threads were killed
       if (reinit) w%forceinit = .true.
    end if
    if (w%forceupdate) call w%update_tree()
    if (w%forcesort) call w%sort_tree(w%table_sortcid,w%table_sortdir)
    if (w%forceinit) then
       call kill_initialization_thread()
       call launch_initialization_thread()
       w%forceinit = .false.
    end if
    nshown = size(w%iord,1)

    ! set up the table, style and flags
    sz%x = 3._c_float
    sz%y = 1._c_float
    call igPushStyleVar_Vec2(ImGuiStyleVar_FramePadding,sz)
    sz%x = 8._c_float
    sz%y = 5._c_float
    call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz)
    sz%x = 2._c_float
    sz%y = 2._c_float
    call igPushStyleVar_Vec2(ImGuiStyleVar_CellPadding,sz)

    str = "Structures##0,0" // c_null_char
    flags = ImGuiTableFlags_Borders
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_Hideable)
    flags = ior(flags,ImGuiTableFlags_Sortable)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    if (igBeginTable(c_loc(str),19,flags,szero,0._c_float)) then
       ! force resize if asked for
       if (w%forceresize) then
          call igTableSetColumnWidthAutoAll(igGetCurrentTable())
          w%forceresize = .false.
       end if

       ! set up the columns
       ! closebutton - ID - name - spg - volume - nneq - ncel - nmol - a - b - c - alpha - beta - gamma
       str = "(close button)##0closebutton" // c_null_char
       flags = ImGuiTableColumnFlags_NoResize
       flags = ior(flags,ImGuiTableColumnFlags_NoReorder)
       flags = ior(flags,ImGuiTableColumnFlags_NoHide)
       flags = ior(flags,ImGuiTableColumnFlags_NoSort)
       flags = ior(flags,ImGuiTableColumnFlags_NoHeaderLabel)
       flags = ior(flags,ImGuiTableColumnFlags_NoHeaderWidth)
       width = max(4._c_float, g%FontSize + 2._c_float)
       call igTableSetupColumn(c_loc(str),flags,width,ic_closebutton)

       str = "(expand button)##0expandbutton" // c_null_char
       call igTableSetupColumn(c_loc(str),flags,width,ic_expandbutton)

       str = "ID##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultSort
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_id)

       str = "Name##0" // c_null_char
       flags = ImGuiTableColumnFlags_WidthStretch
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_name)

       str = "spg##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_spg)

       str = "V/Å³##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_v)

       str = "(V/Z)/Å³##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_vmol)

       str = "nneq##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_nneq)

       str = "ncel##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_ncel)

       str = "nmol##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_nmol)

       str = "a/Å##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_a)

       str = "b/Å##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_b)

       str = "c/Å##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_c)

       str = "α/°##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_alpha)

       str = "β/°##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_beta)

       str = "γ/°##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_gamma)

       str = "E/Ha##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_e)

       str = "(E/Z)/Ha##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_emol)

       str = "p/GPa##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_p)

       call igTableSetupScrollFreeze(0, 1) ! top row always visible

       ! fetch the sort specs, sort the data if necessary
       ptrc = igTableGetSortSpecs()
       if (c_associated(ptrc)) then
          call c_f_pointer(ptrc,sortspecs)
          if (c_associated(sortspecs%Specs)) then
             call c_f_pointer(sortspecs%Specs,colspecs)
             w%table_sortcid = colspecs%ColumnUserID
             w%table_sortdir = colspecs%SortDirection
             if (sortspecs%SpecsDirty .and. nshown > 1) then
                w%forcesort = .true.
                sortspecs%SpecsDirty = .false.
             end if
          else
             w%table_sortcid = ic_id
             w%table_sortdir = 1
          end if
       end if

       ! draw the header
       call igTableHeadersRow()

       ! draw the rows
       do j = 1, nshown
          i = w%iord(j)
          if (sysc(i)%status == sys_empty .or. sysc(i)%hidden) cycle
          if (c_associated(cfilter)) then
             str = trim(sysc(i)%seed%name) // c_null_char
             if (.not.ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) cycle
          end if

          ! start defining the table row
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float);
          hadenabledcolumn = .false.

          ! close button
          buttonhovered_close = .false.
          if (sysc(i)%status == sys_init) then
             if (igTableSetColumnIndex(ic_closebutton)) then
                str = "##1closebutton" // string(ic_closebutton) // "," // string(i) // c_null_char
                if (my_CloseButton(c_loc(str),ColorDangerButton)) w%forceremove = (/i/)
                buttonhovered_close = igIsItemHovered(ImGuiHoveredFlags_None)
             end if
          end if

          ! expand button
          buttonhovered_expand = .false.
          if (sysc(i)%collapse < 0) then
             if (igTableSetColumnIndex(ic_expandbutton)) then
                ! expand button for multi-seed entries
                str = "##expand" // string(ic_expandbutton) // "," // string(i) // c_null_char
                if (sysc(i)%collapse == -1) then
                   idir = ImGuiDir_Right
                else
                   idir = ImGuiDir_Down
                end if
                if (igArrowButton(c_loc(str),idir)) then
                   ! expand or collapse
                   if (sysc(i)%collapse == -1) then
                      call expand_system(i)
                   else
                      call collapse_system(i)
                   end if
                end if
                buttonhovered_expand = igIsItemHovered(ImGuiHoveredFlags_None)
             end if
          end if

          ! set background color for the name cell, if not selected
          ! if (w%table_selected /= i) then
          if (sysc(i)%seed%ismolecule) then
             color = igGetColorU32_Vec4(ColorTableCellBg_Mol)
             if (sysc(i)%status == sys_init) then
                if (sys(i)%c%nmol > 1) color = igGetColorU32_Vec4(ColorTableCellBg_MolClus)
             endif
          else
             color = igGetColorU32_Vec4(ColorTableCellBg_Crys3d)
             if (sysc(i)%status == sys_init) then
                if (sys(i)%c%ismol3d .or. sys(i)%c%nlvac == 3) then
                   color = igGetColorU32_Vec4(ColorTableCellBg_MolCrys)
                elseif (sys(i)%c%nlvac == 2) then
                   color = igGetColorU32_Vec4(ColorTableCellBg_Crys1d)
                elseif (sys(i)%c%nlvac == 1) then
                   color = igGetColorU32_Vec4(ColorTableCellBg_Crys2d)
                end if
             end if
          end if
          call igTableSetBgColor(ImGuiTableBgTarget_CellBg, color, ic_name)

          ! ID column
          if (igTableSetColumnIndex(ic_id)) then
             str = string(i)
             call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
             call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
          end if

          ! name
          if (igTableSetColumnIndex(ic_name)) then
             ! selectable
             call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)

             ! expand button
             if (sysc(i)%showfields) then
                ch = "▼"
             else
                ch = "▶"
             end if
             pos = igGetCursorPosX()
             call igSetCursorPosX(pos + g%Style%FramePadding%x)
             call iw_text(ch)
             call igSameLine(0._c_float,-1._c_float)
             call igSetCursorPosX(pos)
             str = ch // "##" // string(ic_name) // "," // string(i) // c_null_char
             sz%x = iw_calcwidth(1,1)
             sz%y = iw_calcheight(1,0)
             if (igInvisibleButton(c_loc(str),sz,ImGuiButtonFlags_None)) sysc(i)%showfields = .not.sysc(i)%showfields

             ! the actual name
             str = ""
             if (sysc(i)%collapse == -2) then
                str = "╭●─"
                do k = 2, len(string(i))
                   str = str // "─"
                end do
                str = str // "──"
             elseif (sysc(i)%collapse > 0) then
                str = "├[" // string(sysc(i)%collapse) // "]─"
             end if
             str = str // trim(sysc(i)%seed%name)
             call iw_text(str,disabled=(sysc(i)%status /= sys_init),sameline_nospace=.true.,copy_to_output=export)

             ! the fields
             if (sysc(i)%showfields) then
                do k = 0, sys(i)%nf
                   if (.not.sys(i)%f(k)%isinit) cycle

                   ! selectable
                   call igSetCursorPosX(igGetCursorPosX() + iw_calcwidth(1,1))
                   isend = (k == sys(i)%nf)
                   if (.not.isend) isend = all(.not.sys(i)%f(k+1:)%isinit)
                   if (.not.isend) call iw_text("┌",noadvance=.true.)
                   if (sys(i)%iref == k) then
                      str = "└─►(" // string(k) // ",ref): " // trim(sys(i)%f(k)%name) // "##field" // &
                         string(i) // "," // string(k) // c_null_char
                   else
                      str = "└─►(" // string(k) // "): " // trim(sys(i)%f(k)%name) // "##field" // &
                         string(i) // "," // string(k) // c_null_char
                   end if
                   isel = (w%table_selected==i) .and. (sys(i)%iref == k)
                   call igPushStyleColor_Vec4(ImGuiCol_Header,ColorFieldSelected)
                   flags = ImGuiSelectableFlags_SpanAllColumns
                   if (igSelectable_Bool(c_loc(str),isel,flags,szero)) then
                      w%table_selected = i
                      if (tree_select_updates_inpcon) &
                         win(iwin_console_input)%inpcon_selected = i
                      call sys(i)%set_reference(k,.false.)
                   end if
                   call igPopStyleColor(1)

                   ! right click to open the field context menu
                   if (igBeginPopupContextItem(c_loc(str),ImGuiPopupFlags_MouseButtonRight)) then
                      ! remove option (fields)
                      if (k > 0) then
                         strpop = "Remove" // c_null_char
                         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) &
                            call sys(i)%unload_field(k)
                         call iw_tooltip("Remove this field",ttshown)
                      end if

                      ! rename option (fields)
                      strpop = "Rename" // c_null_char
                      if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
                         strpop2 = "##inputrenamefield" // c_null_char
                         txtinp = trim(adjustl(sys(i)%f(k)%name)) // c_null_char
                         call igSetKeyboardFocusHere(0_c_int)
                         flags = ImGuiInputTextFlags_EnterReturnsTrue
                         if (igInputText(c_loc(strpop2),c_loc(txtinp),1023_c_size_t,flags,c_null_ptr,c_null_ptr)) then
                            ll = index(txtinp,c_null_char)
                            sys(i)%f(k)%name = txtinp(1:ll-1)
                            call igCloseCurrentPopup()
                         end if
                         call igEndMenu()
                      end if
                      call iw_tooltip("Rename this field",ttshown)

                      !! now the load new field options !!
                      call igSeparator()

                      ! duplicate option (fields)
                      strpop = "Duplicate" // c_null_char
                      if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) then
                         id = sys(i)%getfieldnum()
                         call sys(i)%field_copy(k,id)
                         sys(i)%f(id)%id = id
                         sys(i)%f(id)%name = trim(sys(i)%f(k)%name)
                      end if
                      call iw_tooltip("Load a copy of this field as a new field",ttshown)

                      ! grid calculation options
                      if (sys(i)%f(k)%type == type_grid) then
                         strpop = "Load gradient grid" // c_null_char
                         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) then
                            id = sys(i)%getfieldnum()
                            call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,"<generated>, gradient of $" // string(k),&
                               sys(i)%f(k)%grid,ifformat_as_grad)
                         end if
                         call iw_tooltip("Load a new grid field as the gradient of this field",ttshown)

                         strpop = "Load Laplacian grid" // c_null_char
                         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) then
                            id = sys(i)%getfieldnum()
                            call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,"<generated>, Laplacian of $" // string(k),&
                               sys(i)%f(k)%grid,ifformat_as_lap)
                         end if
                         call iw_tooltip("Load a new grid field as the Laplacian of this field",ttshown)

                         strpop = "Load potential grid" // c_null_char
                         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) then
                            id = sys(i)%getfieldnum()
                            call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,"<generated>, potential of $" // string(k),&
                               sys(i)%f(k)%grid,ifformat_as_pot)
                         end if
                         call iw_tooltip("Load a new grid field as the potential of this field",ttshown)

                         strpop = "Load resampled grid" // c_null_char
                         if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
                            flags = ImGuiInputTextFlags_None
                            strpop2 = "New size##resamplefieldmenunewsize" // c_null_char
                            call igSetNextItemWidth(iw_calcwidth(4*3,2))
                            ldum = igInputInt3(c_loc(strpop2),iresample,flags)

                            strpop2 = "OK##resamplefieldmenuok" // c_null_char
                            if (igMenuItem_Bool(c_loc(strpop2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                               id = sys(i)%getfieldnum()
                               call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,"<generated>, resample of $" // string(k),&
                                  sys(i)%f(k)%grid,ifformat_as_resample,n=iresample)
                            end if
                            call igEndMenu()
                         else
                            iresample = sys(i)%f(k)%grid%n
                         end if
                         call iw_tooltip("Load a new grid field as a resampling of this field",ttshown)

                      end if


                      call igEndPopup()
                   end if

                   ! tooltip
                   call iw_tooltip(tree_field_tooltip_string(i,k),ttshown)
                end do
             end if
          end if

          ! energy
          if (igTableSetColumnIndex(ic_e)) then
             if (sysc(i)%seed%energy /= huge(1d0)) then
                str = string(sysc(i)%seed%energy,'f',decimal=8)
             else
                str = "n/a"
             end if
             call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
             call iw_text(str,copy_to_output=export)
          end if

          if (sysc(i)%status == sys_init) then
             if (igTableSetColumnIndex(ic_spg)) then ! spg
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                elseif (.not.sys(i)%c%spgavail) then
                   str = "n/a"
                else
                   str = trim(sys(i)%c%spg%international_symbol)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if

             if (igTableSetColumnIndex(ic_v)) then ! volume
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%omega*bohrtoa**3,'f',decimal=2)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if

             if (igTableSetColumnIndex(ic_vmol)) then ! volume per molecule
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%omega*bohrtoa**3/sys(i)%c%nmol,'f',decimal=2)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if

             if (igTableSetColumnIndex(ic_nneq)) then ! nneq
                str = string(sys(i)%c%nneq)
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if

             if (igTableSetColumnIndex(ic_ncel)) then ! ncel
                str = string(sys(i)%c%ncel)
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if

             if (igTableSetColumnIndex(ic_nmol)) then ! nmol
                str = string(sys(i)%c%nmol)
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if

             if (igTableSetColumnIndex(ic_a)) then ! a
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%aa(1)*bohrtoa,'f',decimal=4)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if
             if (igTableSetColumnIndex(ic_b)) then ! b
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%aa(2)*bohrtoa,'f',decimal=4)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if
             if (igTableSetColumnIndex(ic_c)) then ! c
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%aa(3)*bohrtoa,'f',decimal=4)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if
             if (igTableSetColumnIndex(ic_alpha)) then ! alpha
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%bb(1),'f',decimal=2)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if
             if (igTableSetColumnIndex(ic_beta)) then ! beta
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%bb(2),'f',decimal=2)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if
             if (igTableSetColumnIndex(ic_gamma)) then ! gamma
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%bb(3),'f',decimal=2)
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
             end if

             if (igTableSetColumnIndex(ic_emol)) then ! energy/nmol
                if (sysc(i)%seed%energy /= huge(1d0)) then
                   str = string(sysc(i)%seed%energy/sys(i)%c%nmol,'f',decimal=8)
                else
                   str = "n/a"
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,copy_to_output=export)
             end if
             if (igTableSetColumnIndex(ic_p)) then ! pressure
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                elseif (sysc(i)%seed%pressure /= huge(1d0)) then
                   str = string(sysc(i)%seed%pressure,'f',decimal=2)
                else
                   str = "n/a"
                end if
                call write_maybe_selectable(i,buttonhovered_close,buttonhovered_expand)
                call iw_text(str,copy_to_output=export)
             end if
          end if
          ! enter new line
          if (export) write (uout,*)
       end do

       ! process the keybindings
       ok = is_bind_event(BIND_TREE_REMOVE_SYSTEM_FIELD)
       ok = ok .and. igIsWindowFocused(ImGuiFocusedFlags_None)
       if (ok) then
          jsel = w%table_selected
          iref = sys(jsel)%iref
          ok = (jsel >= 1 .and. jsel <= nsys)
          if (ok) ok = sysc(jsel)%showfields
          if (ok) ok = (sysc(jsel)%status == sys_init)
          if (ok) ok = (iref >= 1 .and. iref <= sys(jsel)%nf)
          if (ok) ok = sys(jsel)%f(iref)%isinit
          if (ok) then
             found = .false.
             do k = iref+1, sys(jsel)%nf
                if (sys(jsel)%f(k)%isinit) then
                   call sys(jsel)%set_reference(k,.false.)
                   found = .true.
                   exit
                end if
             end do
             if (.not.found) then
                do k = iref-1, 1, -1
                   if (sys(jsel)%f(k)%isinit) then
                      call sys(jsel)%set_reference(k,.false.)
                      found = .true.
                      exit
                   end if
                end do
             end if
             call sys(jsel)%unload_field(iref)
          else
             w%forceremove = (/jsel/)
          end if
       end if

       call igEndTable()
    end if
    call igPopStyleVar(3_c_int)

    ! if exporting, read the export command
    if (export) &
       ldum = win(iwin_console_input)%read_output_ci(.true.,"[Table export]")

    ! clean up
    ! call ImGuiTextFilter_destroy(cfilter)

  contains

    subroutine write_maybe_selectable(isys,bclose,bexpand)
      use gui_main, only: are_threads_running
      use gui_utils, only: iw_text
      use global, only: iunit, iunit_bohr, iunit_ang
      use tools_io, only: uout
      integer, intent(in) :: isys
      logical, intent(in) :: bclose, bexpand

      integer :: k, idx
      real(c_float) :: pos
      integer(c_int) :: flags, ll, isyscollapse
      logical(c_bool) :: selected, enabled
      logical :: ok
      character(kind=c_char,len=:), allocatable, target :: strl, strpop, strpop2
      character(kind=c_char,len=1024), target :: txtinp

      if (hadenabledcolumn) return

      ! selectable that spans all columns, with zero width
      pos = igGetCursorPosX()
      flags = ImGuiSelectableFlags_SpanAllColumns
      flags = ior(flags,ImGuiSelectableFlags_AllowItemOverlap)
      flags = ior(flags,ImGuiSelectableFlags_AllowDoubleClick)
      flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
      selected = (w%table_selected==isys)
      strl = "##selectable" // string(isys) // c_null_char
      if (igSelectable_Bool(c_loc(strl),selected,flags,szero)) then
         w%table_selected = isys
         if (tree_select_updates_inpcon) &
            win(iwin_console_input)%inpcon_selected = isys
         if (igIsMouseDoubleClicked(ImGuiPopupFlags_MouseButtonLeft)) &
            sysc(isys)%showfields = .true.
      end if
      call igSameLine(0._c_float,-1._c_float)
      call igSetCursorPosX(pos)

      ! right click to open the context menu
      if (igBeginPopupContextItem(c_loc(strl),ImGuiPopupFlags_MouseButtonRight)) then
         ! describe this system in the console output
         strpop = "Describe" // c_null_char
         enabled = (sysc(isys)%status == sys_init) .and..not.are_threads_running()
         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,enabled)) then
            idx = index(sysc(isys)%seed%name,c_null_char)
            strl = "### Describe system (" // string(isys) // "): " // trim(sysc(isys)%seed%name(1:idx-1))
            write (uout,'(/A/)') trim(strl)
            if (sys(isys)%c%ismolecule) then
               iunit = iunit_ang
            else
               iunit = iunit_bohr
            end if
            call sys(isys)%report(.true.,.true.,.true.,.true.,.true.,.true.,.true.)
            iunit = iunit_bohr
            ldum = win(iwin_console_input)%read_output_ci(.true.,"[Describe system " // string(isys) // "]")
         end if
         call iw_tooltip("Print a detailed description of this system in the Output Console",ttshown)

         ! set as current system option
         strpop = "Set as Current System" // c_null_char
         enabled = (sysc(isys)%status == sys_init)
         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,enabled)) &
            win(iwin_console_input)%inpcon_selected = isys
         call iw_tooltip("Set this system as current, the system on which commands are effected",ttshown)

         ! scf energy plot
         if (sysc(isys)%collapse /= 0) then
            strpop = "Plot SCF Iterations" // c_null_char
            if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,enabled)) then
               if (sysc(isys)%collapse < 0) then
                  isyscollapse = isys
               else
                  isyscollapse = abs(sysc(isys)%collapse)
               end if

               if (sysc(isyscollapse)%status == sys_init) then
                  if (sysc(isyscollapse)%idwin_plotscf == 0) then
                     sysc(isyscollapse)%idwin_plotscf = stack_create_window(wintype_scfplot,.true.,isys=isyscollapse)
                  else
                     call igSetWindowFocus_Str(c_loc(win(sysc(isyscollapse)%idwin_plotscf)%name))
                  end if
               end if
            end if
            call iw_tooltip("Plot the energy and other properties as a function of SCF cycle iterations",ttshown)
         end if

         ! load field
         strpop = "Load Field" // c_null_char
         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,enabled)) then
            if (idloadfield > 0) then
               win(idloadfield)%loadfield_isys = isys
            else
               idloadfield = stack_create_window(wintype_load_field,.true.,isys=isys)
            end if
         end if
         call iw_tooltip("Load a scalar field for this system",ttshown)

         ! rename option (system)
         strpop = "Rename" // c_null_char
         if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
            strpop2 = "##inputrename" // c_null_char
            txtinp = trim(adjustl(sysc(isys)%seed%name)) // c_null_char
            call igSetKeyboardFocusHere(0_c_int)
            flags = ImGuiInputTextFlags_EnterReturnsTrue
            if (igInputText(c_loc(strpop2),c_loc(txtinp),1023_c_size_t,flags,c_null_ptr,c_null_ptr)) then
               ll = index(txtinp,c_null_char)
               sysc(isys)%seed%name = txtinp(1:ll-1)
               sysc(isys)%renamed = .true.
               call igCloseCurrentPopup()
            end if
            call igEndMenu()
         end if
         call iw_tooltip("Rename this system",ttshown)

         ! remove option (system)
         ok = enabled
         if (ok) ok = (sys(isys)%nf > 0)
         if (ok) ok = any(sys(isys)%f(1:sys(isys)%nf)%isinit)
         if (ok) then
            strpop = "Remove All Fields" // c_null_char
            if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) then
               do k = 1, sys(isys)%nf
                  if (.not.sys(isys)%f(k)%isinit) cycle
                  call sys(isys)%unload_field(k)
               end do
            end if
            call iw_tooltip("Rename all fields in this system",ttshown)
         end if

         ! remove option (system)
         strpop = "Remove" // c_null_char
         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,enabled)) &
            w%forceremove = (/isys/)
         call iw_tooltip("Remove this system",ttshown)

         call igEndPopup()
      end if

      ! delayed tooltip with info about the system
      if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
         if (bclose) then
            strl = "Close this system" // c_null_char
         elseif (bexpand) then
            strl = "Expand this system" // c_null_char
         else
            strl = tree_system_tooltip_string(isys)
         end if
         call igSetTooltip(c_loc(strl))
      end if
      hadenabledcolumn = .true.

    end subroutine write_maybe_selectable

    ! un-hide the dependents and set as expanded
    subroutine expand_system(i)
      integer, intent(in) :: i

      integer :: k

      if (sysc(i)%status == sys_empty) return
      if (sysc(i)%collapse /= -1) return
      do k = 1, nsys
         if (sysc(k)%collapse == i) then
            if (sysc(k)%status == sys_loaded_not_init) w%forceinit = .true.
            sysc(k)%hidden = .false.
         end if
      end do
      sysc(i)%collapse = -2
      w%forceupdate = .true.

    end subroutine expand_system

    ! hide the dependents and set as collapsed
    subroutine collapse_system(i)
      integer, intent(in) :: i

      integer :: k

      if (sysc(i)%status == sys_empty) return
      if (sysc(i)%collapse /= -2) return
      do k = 1, nsys
         if (sysc(k)%collapse == i) sysc(k)%hidden = .true.
      end do
      sysc(i)%collapse = -1
      ! selected goes to master
      if (w%table_selected >= 1 .and. w%table_selected <= nsys) then
         if (sysc(w%table_selected)%collapse == i) w%table_selected = i
      end if
      if (win(iwin_console_input)%inpcon_selected >= 1 .and. win(iwin_console_input)%inpcon_selected <= nsys) then
         if (sysc(win(iwin_console_input)%inpcon_selected)%collapse == i) &
            win(iwin_console_input)%inpcon_selected = i
      end if
      w%forceupdate = .true.

    end subroutine collapse_system

  end subroutine draw_tree

  ! Update the table rows by building a new row index array
  ! (iord). Only the systems that are not empty are pointed by
  ! iord. This is routine is used when the systems change.
  module subroutine update_tree(w)
    use gui_main, only: sysc, nsys, sys_empty
    class(window), intent(inout) :: w

    integer :: i, n

    if (allocated(w%iord)) deallocate(w%iord)
    n = count(sysc(1:nsys)%status /= sys_empty.and..not.sysc(1:nsys)%hidden)
    allocate(w%iord(max(n,1)))
    w%iord(1) = 1
    if (n > 0) then
       n = 0
       do i = 1, nsys
          if (sysc(i)%status /= sys_empty .and..not.sysc(i)%hidden) then
             n = n + 1
             w%iord(n) = i
          end if
       end do
    end if
    w%forceupdate = .false.
    w%forcesort = .true.

  end subroutine update_tree

  ! Sort the table row order by column cid and in direction dir
  ! (ascending=1, descending=2). Modifies the w%iord.
  module subroutine sort_tree(w,cid,dir)
    use gui_main, only: sys, sysc, sys_init, sys_empty
    use tools, only: mergesort
    use tools_math, only: invert_permutation
    use tools_io, only: ferror, faterr
    use types, only: vstring
    class(window), intent(inout) :: w
    integer(c_int), intent(in) :: cid, dir

    integer :: i, n, nvalid, nnovalid
    integer, allocatable :: ival(:), iperm(:), ivalid(:), inovalid(:)
    logical, allocatable :: valid(:)
    real*8, allocatable :: rval(:)
    type(vstring), allocatable :: sval(:)
    logical :: doit

    ! initialize the identity permutation
    n = size(w%iord,1)
    allocate(iperm(n),valid(n))
    do i = 1, n
       iperm(i) = i
    end do
    valid = .true.

    ! different types, different sorts
    if (cid == ic_id .or. cid == ic_nneq .or. cid == ic_ncel .or. cid == ic_nmol) then
       ! sort by integer
       allocate(ival(n))
       do i = 1, n
          if (cid == ic_id) then
             ival(i) = w%iord(i)
          elseif (sysc(w%iord(i))%status == sys_init) then
             if (cid == ic_nneq) then
                ival(i) = sys(w%iord(i))%c%nneq
             elseif (cid == ic_ncel) then
                ival(i) = sys(w%iord(i))%c%ncel
             elseif (cid == ic_nmol) then
                ival(i) = sys(w%iord(i))%c%nmol
             end if
          else
             ival(i) = huge(1)
             valid(i) = .false.
          end if
       end do
       call mergesort(ival,iperm,1,n)
       deallocate(ival)
    elseif (cid == ic_v .or. cid == ic_a .or. cid == ic_b .or. cid == ic_c .or.&
       cid == ic_alpha .or. cid == ic_beta .or. cid == ic_gamma .or. cid == ic_vmol.or.&
       cid == ic_e .or. cid == ic_emol .or. cid == ic_p) then
       ! sort by real
       allocate(rval(n))
       do i = 1, n
          doit = sysc(w%iord(i))%status == sys_init
          if (doit) doit = (.not.sys(w%iord(i))%c%ismolecule)
          if (cid == ic_v .and. doit) then
             rval(i) = sys(w%iord(i))%c%omega
          elseif (cid == ic_a .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(1)
          elseif (cid == ic_b .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(2)
          elseif (cid == ic_c .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(3)
          elseif (cid == ic_alpha .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(1)
          elseif (cid == ic_beta .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(2)
          elseif (cid == ic_gamma .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(3)
          elseif (cid == ic_vmol .and. doit) then
             rval(i) = sys(w%iord(i))%c%omega / sys(w%iord(i))%c%nmol
          elseif (cid == ic_e .and. sysc(w%iord(i))%seed%energy /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%energy
          elseif (cid == ic_emol .and. sysc(w%iord(i))%status == sys_init .and. sysc(w%iord(i))%seed%energy /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%energy / sys(w%iord(i))%c%nmol
          elseif (cid == ic_p .and. sysc(w%iord(i))%seed%pressure /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%pressure
          else
             rval(i) = huge(1d0)
             valid(i) = .false.
          end if
       end do
       call mergesort(rval,iperm,1,n)
       deallocate(rval)
    elseif (cid == ic_name .or. cid == ic_spg) then
       ! sort by string
       allocate(sval(n))
       do i = 1, n
          if (cid == ic_name .and. sysc(w%iord(i))%status /= sys_empty .and..not.sysc(w%iord(i))%hidden) then
             sval(i)%s = trim(sysc(w%iord(i))%seed%name)
          else
             doit = (cid == ic_spg) .and. (sysc(w%iord(i))%status == sys_init)
             if (doit) doit = .not.sys(w%iord(i))%c%ismolecule
             if (doit) doit = sys(w%iord(i))%c%spgavail
             if (doit) then
                sval(i)%s = trim(sys(w%iord(i))%c%spg%international_symbol)
             else
                sval(i)%s = ""
                valid(i) = .false.
             end if
          end if
       end do
       call mergesort(sval,iperm,1,n)
       deallocate(sval)
    else
       call ferror('sort_tree','column sorting not implemented',faterr)
    end if
    valid = valid(iperm)

    ! reverse the permutation, if requested; put the no-valids at the end
    allocate(ivalid(count(valid)),inovalid(n-count(valid)))
    nvalid = 0
    nnovalid = 0
    do i = 1, n
       if (valid(i)) then
          nvalid = nvalid + 1
          ivalid(nvalid) = iperm(i)
       else
          nnovalid = nnovalid + 1
          inovalid(nnovalid) = iperm(i)
       end if
    end do
    do i = nvalid, 1, -1
       if (dir == 2) then
          iperm(nvalid-i+1) = ivalid(i)
       else
          iperm(i) = ivalid(i)
       end if
    end do
    do i = 1, nnovalid
       iperm(nvalid+i) = inovalid(i)
    end do
    deallocate(ivalid,inovalid)

    ! apply the permutation
    w%iord = w%iord(iperm)
    w%forcesort = .false.

  end subroutine sort_tree

  !> Draw the open files dialog.
  module subroutine draw_dialog(w)
    use gui_keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG
    use gui_main, only: add_systems_from_name, launch_initialization_thread,&
       system_shorten_names
    use c_interface_module, only: C_F_string_alloc, c_free
    use tools_io, only: ferror, faterr, fopen_write, fclose
    use param, only: dirsep, bohrtoa
    class(window), intent(inout), target :: w

    type(ImVec2) :: minsize, maxsize
    type(IGFD_Selection_Pair), pointer :: s(:)
    type(IGFD_Selection) :: sel
    type(c_ptr) :: cstr
    integer(c_size_t) :: i
    character(len=:), allocatable :: name, path
    logical :: readlastonly
    integer :: lu

    ! permutation for the open file format list (see dialog_user_callback)
    integer, parameter :: isperm(0:30) = (/0,7,5,1,11,4,20,3,27,8,28,29,15,17,13,&
       14,26,25,16,24,9,10,22,2,18,30,21,6,23,19,12/)

    ! set initial, minimum, and maximum sizes
    minsize%x = 0._c_float
    minsize%y = 0._c_float
    maxsize%x = 1e10_c_float
    maxsize%y = 1e10_c_float

    ! process the dialog
    if (IGFD_DisplayDialog(w%dptr,c_loc(w%name),w%flags,minsize,maxsize)) then
       ! the dialog has been closed
       if (IGFD_IsOk(w%dptr)) then
          ! with an OK, gather information
          if (w%dialog_purpose == wpurp_dialog_openfiles) then
             !! open files dialog !!
             ! open all files selected and add the new systems
             sel = IGFD_GetSelection(w%dptr)
             call c_f_pointer(sel%table,s,(/sel%count/))
             cstr = IGFD_GetCurrentPath(w%dptr)
             call C_F_string_alloc(cstr,path)
             call c_free(cstr)
             do i = 1, sel%count
                call C_F_string_alloc(s(i)%fileName,name)
                name = trim(path) // dirsep // trim(name)
                readlastonly = w%dialog_data%readlastonly
                call add_systems_from_name(name,w%dialog_data%mol,isperm(w%dialog_data%isformat),&
                   readlastonly,real(w%dialog_data%rborder/bohrtoa,8),logical(w%dialog_data%molcubic))
             end do

             ! initialize
             call launch_initialization_thread()

             ! shorten the names
             call system_shorten_names()

          elseif (w%dialog_purpose == wpurp_dialog_savelogfile) then
             !! save log file dialog !!
             cstr = IGFD_GetFilePathName(w%dptr)
             call C_F_string_alloc(cstr,name)
             call c_free(cstr)

             lu = fopen_write(name,errstop=.false.)
             if (lu < 0) &
                call ferror('draw_dialog','could not open file for writing: ' // name,faterr,syntax=.true.)
             if (idcom == 0 .and. allocated(outputb)) then
                write(lu,'(A)') outputb(1:lob)
             elseif (idcom > 0) then
                write(lu,'(A)') com(icom(idcom))%output
             end if
             call fclose(lu)
          elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile .or. &
             w%dialog_purpose == wpurp_dialog_openfieldfile .or. w%dialog_purpose == wpurp_dialog_openonefilemodal) then
             !! new structure file dialog !!
             cstr = IGFD_GetFilePathName(w%dptr)
             call C_F_string_alloc(cstr,name)
             call c_free(cstr)
             w%okfile = trim(name)
             w%okfile_set = .true.
             w%okfile_read = .true.
          else
             call ferror('draw_dialog','unknown dialog purpose',faterr)
          end if
       end if

       ! close the dialog and terminate the window
       call IGFD_CloseDialog(w%dptr)
       call w%end()
    end if

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) &
       call IGFD_ForceQuit(w%dptr)

    ! exit if focused and received the OK keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) &
       call IGFD_ForceOK(w%dptr)

  end subroutine draw_dialog

  !> Draw the contents of the input console
  module subroutine draw_ci(w)
    use gui_keybindings, only: BIND_INPCON_RUN, get_bind_keyname, is_bind_event
    use gui_main, only: sys, sysc, nsys, sys_init, g, force_run_commands
    use gui_templates, only: draw_keyword_context_menu
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text
    use systemmod, only: sy
    use tools_io, only: string
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: sz, szero, szavail
    logical(c_bool) :: ldum, is_selected
    real(c_float) :: combowidth
    logical :: ok
    integer :: i, idx

    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    szero%x = 0
    szero%y = 0

    ! allocate the input buffer if not already done
    if (.not.allocated(inputb)) then
       allocate(character(len=maxlib+1) :: inputb)
       inputb(1:1) = c_null_char
    end if

    ! first line: text
    call iw_text("Input",highlight=.true.)

    ! first line: clear button
    if (iw_button("Clear",sameline=.true.)) &
       inputb(1:1) = c_null_char
    call iw_tooltip("Clear the input text",ttshown)

    ! first line: template button
    ldum = iw_button("Template",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
    if (ok) &
       call draw_keyword_context_menu(.true.)
    call iw_tooltip("Insert a template for a critic2 command",ttshown)

    ! first line: help button
    ldum = iw_button("Help",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
    if (ok) &
       call draw_keyword_context_menu(.false.)
    call iw_tooltip("Bring up the critic2 command reference",ttshown)

    ! second line: calculate size of the RUN button
    sz%x = 2 * (igGetTextLineHeight() + 2 * g%Style%FramePadding%y) + g%Style%ItemSpacing%y
    sz%y = sz%x

    ! second line: system selector
    call igBeginGroup()
    call iw_text("System")
    call igSameLine(0._c_float,-1._c_float)

    call igGetContentRegionAvail(szavail)
    combowidth = max(szavail%x - sz%x - g%Style%ItemSpacing%x,0._c_float)

    !! set the system pointer and determine the preview string
    str2 = "" // c_null_char
    sy => null()
    if (w%inpcon_selected >= 1 .and. w%inpcon_selected <= nsys) then
       if (sysc(w%inpcon_selected)%status == sys_init) then
          str2 = string(w%inpcon_selected) // ": " // trim(sysc(w%inpcon_selected)%seed%name) // c_null_char
          sy => sys(w%inpcon_selected)
       end if
    end if
    str1 = "##systemcombo" // c_null_char
    call igSetNextItemWidth(combowidth)
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nsys
          if (sysc(i)%status == sys_init) then
             is_selected = (w%inpcon_selected == i)
             str2 = string(i) // ": " // trim(sysc(i)%seed%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) &
                w%inpcon_selected = i
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Set the current system (input commands are applied to it)",ttshown)

    ! third line: field selector
    call iw_text("Field ")
    call igSameLine(0._c_float,-1._c_float)

    str1 = "##fieldcombo" // c_null_char
    str2 = "" // c_null_char
    if (associated(sy)) &
       str2 = string(sy%iref) // ": " // trim(sy%f(sy%iref)%name) // c_null_char
    call igSetNextItemWidth(combowidth)
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       if (associated(sy)) then
          do i = 0, sy%nf
             if (.not.sy%f(i)%isinit) cycle
             is_selected = (sy%iref == i)
             str2 = string(i) // ": " // trim(sy%f(i)%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) &
                call sy%set_reference(i,.false.)
             if (is_selected) &
                call igSetItemDefaultFocus()
          end do
       end if
       call igEndCombo()
    end if
    call iw_tooltip("Set the reference field (input commands are applied to it)",ttshown)
    call igEndGroup()

    ! right-hand-side of lines 2 and 3: RUN button
    ok = iw_button("RUN",danger=.true.,sameline=.true.,siz=(/sz%x,sz%y/))
    ok = ok .or. (igIsWindowFocused(ImGuiFocusedFlags_None) .and. is_bind_event(BIND_INPCON_RUN))
    if (ok) then
       idx = index(inputb,c_null_char)
       if (idx > 1) then
          if (associated(sy)) force_run_commands = .true.
       end if
    end if
    call iw_tooltip("Run the commands (" // get_bind_keyname(BIND_INPCON_RUN) // ")",ttshown)

    ! calculate sizes and draw the multiline
    call igGetContentRegionAvail(sz)
    str1 = "##inputmultiline" // c_null_char
    ldum = igInputTextMultiline(c_loc(str1),c_loc(inputb),maxlib,sz,&
       ImGuiInputTextFlags_AllowTabInput,c_null_ptr,c_null_ptr)

  end subroutine draw_ci

  !> Run the commands from the console input
  module subroutine run_commands_ci(w)
    use gui_main, only: launch_initialization_thread, kill_initialization_thread, are_threads_running
    use global, only: critic_main
    use tools_io, only: falloc, uin, fclose, ferror, faterr
    use iso_fortran_env, only: input_unit
    class(window), intent(inout), target :: w

    integer :: idx
    integer :: ios
    logical :: reinit, ldum

    ! check we have some input
    idx = index(inputb,c_null_char)
    if (idx <= 1) return

    ! if the initialization is happening, stop it
    reinit = are_threads_running()
    if (reinit) call kill_initialization_thread()

    ! connect a scratch file to uin, write the commands, rewind, and run
    uin = falloc()
    open(unit=uin,status='scratch',form='formatted',access='stream',iostat=ios)
    if (ios /= 0) &
       call ferror("run_commands","cannot open buffer for critic2 input",faterr)
    write (uin,'(A)') inputb(1:idx-1)
    rewind(uin)
    call critic_main()

    ! read the output
    ldum = w%read_output_ci(.true.)

    ! reinitialize the threads
    if (reinit) call launch_initialization_thread()

    ! clean up
    call fclose(uin)
    uin = input_unit

  end subroutine run_commands_ci

  !> Block the GUI by dimming the background and showing the current console
  !> input. Useful for preparing the screen for running the commands
  !> in the console input.
  module subroutine block_gui_ci(w)
    use gui_main, only: mainvwp, io, g, ColorWaitBg
    use gui_utils, only: iw_text
    use param, only: newline
    class(window), intent(inout), target :: w

    integer(c_int) :: flags
    type(ImVec2) :: sz, pivot
    character(kind=c_char,len=:), allocatable, target :: str1, text
    character(kind=c_char,len=:), allocatable, target :: csystem, cfield
    logical(c_bool) :: ldum

    !! blank the background
    flags = ImGuiWindowFlags_NoDecoration
    flags = ior(flags,ImGuiWindowFlags_NoMove)
    flags = ior(flags,ImGuiWindowFlags_NoSavedSettings)
    sz%x = 0._c_float
    sz%y = 0._c_float
    call igSetNextWindowPos(mainvwp%WorkPos,0,sz)
    call igSetNextWindowSize(mainvwp%WorkSize,0)
    call igPushStyleColor_Vec4(ImGuiCol_WindowBg,ColorWaitBg)
    str1 = "##blankbackground" // c_null_char
    ldum = .true.
    call igSetNextWindowFocus()
    ldum = igBegin(c_loc(str1), ldum, flags)
    call igEnd()
    call igPopStyleColor(1)

    !! overlay
    ! set window position at the center
    sz%x = io%DisplaySize%x * 0.5_c_float
    sz%y = io%DisplaySize%y * 0.5_c_float
    pivot%x = 0.5_c_float
    pivot%y = 0.5_c_float
    call igSetNextWindowPos(sz,0,pivot)

    ! get the input details
    call w%get_input_details_ci(csystem,cfield)

    ! set window size
    text = "...Running critic2 input..." // newline //&
       "System: " // csystem // newline //&
       "Field:  " // cfield // newline //&
       "Input:  " // newline // inputb
    call igCalcTextSize(sz,c_loc(text),c_null_ptr,.false._c_bool,-1._c_float)
    sz%y = sz%y + 2 * g%Style%WindowPadding%y
    sz%x = sz%x + 2 * g%Style%WindowPadding%x
    call igSetNextWindowSize(sz,0)

    ! draw the window
    flags = ImGuiWindowFlags_NoDecoration
    flags = ior(flags,ImGuiWindowFlags_NoDocking)
    flags = ior(flags,ImGuiWindowFlags_AlwaysAutoResize)
    flags = ior(flags,ImGuiWindowFlags_NoSavedSettings)
    flags = ior(flags,ImGuiWindowFlags_NoFocusOnAppearing)
    flags = ior(flags,ImGuiWindowFlags_NoNav)
    ldum = .true.
    str1 = "##popupwait" // c_null_char
    call igSetNextWindowFocus()
    if (igBegin(c_loc(str1), ldum, flags)) then
       sz%x = 0
       sz%y = 0
       call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz)
       call iw_text("...Running critic2 input...",highlight=.true.)
       call iw_text("System: ",highlight=.true.)
       call iw_text(csystem,sameline=.true.)
       call iw_text("Field:  ",highlight=.true.)
       call iw_text(cfield,sameline=.true.)
       call iw_text("Input:  ",highlight=.true.)
       call igIndent(0._c_float)
       call igText(c_loc(inputb))
       call igPopStyleVar(1_c_int)
       call igUnindent(0._c_float)
    end if
    call igEnd()

  end subroutine block_gui_ci

  !> Read new output from the scratch LU uout. If iscom, this output corresponds
  !> to a command, so create the command i/o object. Return true if
  !> output has been read.
  module function read_output_ci(w,iscom,cominfo)
    use gui_utils, only: get_time_string
    use gui_main, only: are_threads_running
    use tools_io, only: uout, getline_raw, string, ferror, faterr
    use types, only: realloc
    use param, only: newline
    class(window), intent(inout), target :: w
    logical, intent(in) :: iscom
    character*(*), intent(in), optional :: cominfo
    logical :: read_output_ci

    character(kind=c_char,len=:), allocatable, target :: csystem, cfield
    type(command_inout), allocatable :: aux(:)
    character(len=:), allocatable :: line, commonstr, showinp
    integer(c_size_t) :: pos, lshift, ll, olob, total, newsize
    integer :: idx, ithis, i
    logical :: ok

    ! allocate the output buffer if not allocated
    read_output_ci = .false.
    if (.not.allocated(outputb)) then
       allocate(character(len=maxlob+1) :: outputb)
       outputb(1:1) = c_null_char
       lob = 0
    end if

    ! allocate the command input/output stack, if not allocated
    if (.not.allocated(com)) then
       ncom = 0
       allocate(com(10))
    end if
    if (.not.allocated(icom)) then
       nicom = 0
       allocate(icom(10))
    end if

    ! do not read if initialization threads are still running (may generate output)
    if (are_threads_running()) return

    ! get size of the new output
    inquire(uout,pos=pos)
    if (pos > 1) then
       ! there is new output, rewind, read it, and rewind again

       ! I am going to need pos-1 new characters
       ll = pos-1
       rewind(uout)

       ! do not have enough room to accomodate this much output = discard new text
       do while (ll > maxlob)
          ok = getline_raw(uout,line)
          if (.not.ok) return
          ll = ll - (len(line)+1)
       end do

       ! check whether we can accomodate the new text = discard from the beginning
       if (lob + ll > maxlob) then
          lshift = ll + lob - maxlob
          outputb(1:lob-lshift) = outputb(lshift+1:lob)
          lob = lob - lshift
          ! adjust to the next newline
          idx = index(outputb,newline)
          if (idx == 0) then
             lob = 0
          else
             outputb(1:lob-idx) = outputb(idx+1:lob)
             lob = lob - idx
          end if
       end if

       ! read the new output and rewind
       olob = lob
       do while(getline_raw(uout,line))
          ll = len(line)
          if (ll > 0) then
             outputb(lob+1:lob+ll) = line(1:ll)
             lob = lob + ll
          end if
          outputb(lob+1:lob+1) = newline
          lob = lob + 1
          inquire(uout,pos=pos)
       end do
       outputb(lob+1:lob+1) = c_null_char
       rewind(uout)

       if (iscom) then
          ! try to remove commands if we exceed the size
          newsize = lob - olob + 1
          if (newsize > maxcomout) &
             call ferror('read_output_ci','exceeded output buffer for command',faterr)
          total = newsize
          do i = 1, nicom
             total = total + com(icom(i))%size
          end do
          i = 0
          do while (total > maxcomout)
             i = i + 1
             total = total - com(icom(i))%size
             call com(icom(i))%end()
          end do
          icom(1:nicom-i) = icom(i+1:nicom)
          nicom = nicom - i

          ! get an unusued command ID
          ithis = 0
          do i = 1, ncom
             if (com(i)%status == command_inout_empty) then
                ithis = i
                exit
             end if
          end do

          ! add a new command, reallocate if necessary
          if (ithis == 0) then
             ncom = ncom + 1
             if (ncom > size(com,1)) then
                allocate(aux(2*ncom))
                aux(1:size(com,1)) = com
                call move_alloc(aux,com)
             end if
             ithis = ncom
          end if

          ! fill the new command info
          ncomid = ncomid + 1
          call w%get_input_details_ci(csystem,cfield)
          com(ithis)%id = ncomid
          com(ithis)%status = command_inout_used
          com(ithis)%size = lob - olob + 1
          if (present(cominfo)) then
             showinp = trim(cominfo)
             com(ithis)%input = c_null_char
          else
             idx = index(inputb,c_null_char)
             showinp = inputb(1:idx-1)
             com(ithis)%input = inputb(1:idx)
          end if

          commonstr = "### Command: " // string(ncomid) // " (" // get_time_string() // ")" // newline
          com(ithis)%tooltipinfo = commonstr // &
             "System: " // csystem // newline //&
             "Field: " // cfield // newline //&
             "Input: " // newline // showinp // newline //&
             "#########" // newline // newline // "[Right-click for options]"
          com(ithis)%output = commonstr // &
             "## System: " // csystem // newline //&
             "## Field: " // cfield // newline //&
             "## Input: " // newline // showinp // newline //&
             "#########" // newline // newline // outputb(olob+1:lob+1)

          ! add it to the list of active commands
          nicom = nicom + 1
          if (nicom > size(icom,1)) &
             call realloc(icom,2*nicom)
          icom(nicom) = ithis
       end if

       ! we have new data
       read_output_ci = .true.
    end if

  end function read_output_ci

  !> Fill the input buffer with the given string.
  module subroutine fill_input_ci(w,str)
    class(window), intent(inout), target :: w
    character(len=*), intent(in) :: str

    integer :: idx

    idx = index(str,c_null_char)
    if (idx == 0) then
       idx = len_trim(str)
    else
       idx = idx - 1
    end if
    inputb(1:idx+1) = trim(str(1:idx)) // c_null_char

  end subroutine fill_input_ci

  !> Get the system and field strings for current input (without null char).
  module subroutine get_input_details_ci(w,csystem,cfield)
    use gui_main, only: nsys, sysc, sys_init, sys
    use tools_io, only: string
    class(window), intent(inout), target :: w
    character(len=:), allocatable, intent(inout) :: csystem, cfield

    integer :: iref

    ! system name and field name
    csystem = "<unknown>"
    cfield = "<unknown>"
    if (w%inpcon_selected >= 1 .and. w%inpcon_selected <= nsys) then
       if (sysc(w%inpcon_selected)%status == sys_init) then
          csystem = "(" // string(w%inpcon_selected) // ") " // trim(sysc(w%inpcon_selected)%seed%name)
          iref = sys(w%inpcon_selected)%iref
          cfield = "(" // string(iref) // ") " // trim(sys(w%inpcon_selected)%f(iref)%name)
       end if
    end if

  end subroutine get_input_details_ci

  !> Draw the contents of the output console
  module subroutine draw_co(w)
    use gui_main, only: g, ColorDangerButton, ColorFrameBgAlt
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text,&
       iw_setposx_fromend, iw_calcheight, iw_calcwidth
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: i, curline, ndrawn, idx
    character(kind=c_char,len=:), allocatable, target :: str1, strpop
    type(ImVec2) :: sz, szero, szavail
    logical(c_bool) :: ldum
    logical :: setscroll, skip, pushed, ok
    real(c_float) :: itemspacing, xavail, xavail1, xx, ally
    integer :: navail, navail1

    real(c_float), save :: maxallscrolly = 0._c_float ! max scroll value of the All pane
    real(c_float), save :: allscrolly = 0._c_float ! scroll value of the All pane
    logical, save :: ttshown = .false. ! tooltip flag
    integer, save :: idsavedialog = 0 ! the ID of the save window

    ! initialize
    szero%x = 0._c_float
    szero%y = 0._c_float
    setscroll = .false.

    ! check if the save dialog is still open
    call update_window_id(idsavedialog)

    ! read new output, if available
    ldum = w%read_output_ci(.false.)

    ! get the available x
    call igGetContentRegionAvail(szavail)
    xavail = szavail%x

    ! first line: text
    call iw_text("Output",highlight=.true.)

    ! first line: clear button
    if (idcom == 0) then
       if (iw_button("Clear",sameline=.true.)) then
          outputb(1:1) = c_null_char
          lob = 0
       end if
       call iw_tooltip("Clear the output log",ttshown)
    end if

    ! first line: copy button
    if (iw_button("Copy",sameline=.true.)) then
       if (idcom == 0) then
          call igSetClipboardText(c_loc(outputb))
       else
          call igSetClipboardText(c_loc(com(icom(idcom))%output))
       end if
    end if
    call iw_tooltip("Copy the shown output log to clipboard",ttshown)

    ! first line: save button
    if (iw_button("Save",disabled=(idsavedialog > 0),sameline=.true.)) &
       idsavedialog = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savelogfile)
    call iw_tooltip("Save the shown output log to a file",ttshown)

    ! first line: remove all button
    call iw_setposx_fromend(10,1)

    if (iw_button("Remove All",danger=.true.)) then
       ! remove all command i/o information
       ncomid = 0
       ncom = 0
       nicom = 0
       idcom = 0
       icom = 0
       do i = 1, ncom
          call com(i)%end()
       end do
    end if
    call iw_tooltip("Remove all commands",ttshown)

    ! second line: all button
    xavail1 = max(xavail - iw_calcwidth(3,1) + g%Style%ItemSpacing%x,0._c_float)
    if (idcom == 0) then
       call igPushStyleColor_Vec4(ImGuiCol_Button,g%Style%Colors(ImGuiCol_ButtonActive+1))
    else
       call igPushStyleColor_Vec4(ImGuiCol_Button,ColorDangerButton)
    end if
    if (iw_button("All")) then
       idcom = 0
       setscroll = .true.
    end if
    call iw_tooltip("Show all console output",ttshown)
    call igPopStyleColor(1)

    !! second line: list of command i/os
    ! make the itemspacing zero
    sz%x = 0
    sz%y = g%Style%ItemSpacing%y
    itemspacing = g%Style%ItemSpacing%x
    call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz)

    ! calculate the button size and the number of buttons that fit
    xx = iw_calcwidth(ceiling(log10(max(maxval(com(icom(1:nicom))%id),1) + 0.1)),1)
    xavail = xavail / xx
    xavail1 = xavail1 / xx
    navail = max(floor(xavail),1)
    navail1 = max(floor(xavail1),1)

    ! render the buttons
    curline = 1
    ndrawn = 0
    do i = nicom, 1,  -1
       ! skip if no more buttons can be shown
       if (curline == 1) then
          skip = (ndrawn >= navail1)
       else
          skip = (ndrawn >= navail)
       end if

       ! skip or sameline
       if (skip) then
          ndrawn = 0
          curline = curline + 1
       else
          if (i == nicom) then
             call igSameLine(0._c_float,itemspacing)
          else
             call igSameLine(0._c_float,-1._c_float)
          end if
       end if

       ! render the button in alternate colors
       pushed = .true.
       if (idcom == i) then
          call igPushStyleColor_Vec4(ImGuiCol_Button,g%Style%Colors(ImGuiCol_ButtonActive+1))
       elseif (mod(i,2) == 0) then
          call igPushStyleColor_Vec4(ImGuiCol_Button,ColorFrameBgAlt)
       else
          pushed = .false.
       end if
       if (iw_button(string(com(icom(i))%id),popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonRight,&
          siz=(/xx,iw_calcheight(1,0)/))) then
          idcom = i
          setscroll = .true.
       end if
       if (pushed) call igPopStyleColor(1)

       ! tooltip
       call iw_tooltip(com(icom(i))%tooltipinfo)

       ! context menu
       if (ok) then
          strpop = "Edit Input" // c_null_char
          if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) then
             idx = index(com(icom(i))%input,c_null_char)
             if (idx > 0) then
                call w%fill_input_ci(com(icom(i))%input(1:idx))
             end if
             call igFocusWindow(win(iwin_console_input)%ptr)
          end if

          strpop = "Remove" // c_null_char
          if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) then
             call com(icom(i))%end()
             icom(i:nicom-1) = icom(i+1:nicom)
             nicom = nicom - 1
             if (idcom == i) then
                idcom = 0
             elseif (idcom > i) then
                idcom = idcom - 1
             end if
          end if

          call igEndPopup()
       end if

       ndrawn = ndrawn + 1
    end do
    call igPopStyleVar(1_c_int)

    ! calculate sizes and draw the multiline (with dark background and border)
    call igGetContentRegionAvail(sz)
    call igPushStyleColor_Vec4(ImGuiCol_FrameBg,g%Style%Colors(ImGuiCol_WindowBg+1))
    call igPushStyleVar_Float(ImGuiStyleVar_FrameBorderSize,1._c_float)
    str1 = "##outputmultiline" // c_null_char
    if (idcom == 0) then
       ldum = igInputTextMultiline(c_loc(str1),c_loc(outputb),lob,sz,&
          ImGuiInputTextFlags_ReadOnly,c_null_ptr,c_null_ptr)
    else
       ldum = igInputTextMultiline(c_loc(str1),c_loc(com(icom(idcom))%output),&
          com(icom(idcom))%size+1,sz,ImGuiInputTextFlags_ReadOnly,&
          c_null_ptr,c_null_ptr)
    end if
    call igPopStyleVar(1)
    call igPopStyleColor(1)

    ! Set scroll to previous value if changing output, the All view scrolls to the
    ! end on new output
    ldum = igBeginChild_Str(c_loc(str1),szero,.false._c_bool,ImGuiWindowFlags_None)
    if (setscroll) then
       if (idcom > 0) then
          call igSetScrollY_Float(com(icom(idcom))%scrolly)
       else
          call igSetScrollY_Float(allscrolly)
       end if
    endif

    ! if there is new output in all, set the scroll to the beginning
    ! of the new output
    if (idcom == 0) then
       ally = igGetScrollMaxY()
       if (abs(ally - maxallscrolly) > 1d-5) then
          if (maxallscrolly > 0._c_float) then
             call igGetItemRectSize(sz)
             allscrolly = min(maxallscrolly + sz%y - 2 * igGetTextLineHeight(),ally)
             call igSetScrollY_Float(allscrolly)
          end if
          maxallscrolly = ally
       end if
    end if

    ! write down the current scroll y
    if (idcom > 0) then
       com(icom(idcom))%scrolly = igGetScrollY()
    else
       allscrolly = igGetScrollY()
    end if
    call igEndChild()

  end subroutine draw_co

  !> Draw the contents of the new structure window.
  module subroutine draw_new_struct(w)
    use gui_keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG
    use gui_main, only: g, add_systems_from_seeds,&
       launch_initialization_thread, system_shorten_names
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text, iw_calcheight,&
       iw_calcwidth, buffer_to_string_array, iw_radiobutton, iw_combo_simple
    use crystalseedmod, only: crystalseed, realloc_crystalseed
    use global, only: rborder_def
    use tools_io, only: string, fopen_scratch, fclose, stripchar, deblank
    use types, only: vstring
    use param, only: newline, bohrtoa
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, str2, stropt
    logical(c_bool) :: ldum, doquit
    logical :: ok
    type(ImVec2) :: szero, sz, szavail
    integer :: i, idx, lu
    type(crystalseed), allocatable :: seed_(:)
    integer(c_int) :: iunitat

    ! window state
    logical, save :: ttshown = .false. ! tooltip flags
    integer(c_size_t), parameter :: maxatposbuf = 100000 !! atomic positions buffer
    character(kind=c_char,len=:), allocatable, target, save :: atposbuf
    integer(c_size_t), parameter :: maxsymopbuf = 10000 !! symmetry operations buffer
    character(kind=c_char,len=:), allocatable, target, save :: symopbuf
    integer(c_size_t), parameter :: maxlatvecbuf = 5000 !! lattice vectors buffer
    character(kind=c_char,len=:), allocatable, target, save :: latvecbuf
    integer(c_size_t), parameter :: maxnamebuf = 1024 !! name buffer
    character(len=:,kind=c_char), allocatable, target, save :: namebuf
    logical, save :: ismolecule = .false. ! whether the system is a molecule or a crystal
    integer, save :: iunitat_c = 2 ! crystal atpos units (0 = bohr, 1 = angstrom, 2 = fractional)
    integer, save :: iunitat_m = 2 ! mol atpos units (0 = bohr, 1 = angstrom)
    integer, save :: symopt = 1 ! symmetry option (1 = detect, 2 = spg, 3 = manual)
    integer, save :: ispg_selected = 1 ! selected space group
    integer, save :: cellopt = 1 ! lattice option (1 = parameters, 2 = lattice-vectors)
    integer, save :: iunitcel = 0 ! units for lattice (0 = bohr, 1 = angstrom)
    real(c_float), save :: scale = 1._c_float ! scale factor for lattice vectors
    real(c_float), save :: aa(3) = 10._c_float ! cell lengths
    real(c_float), save :: bb(3) = 90._c_float ! cell angles
    real(c_float), save :: rborder = real(rborder_def*bohrtoa,c_float) ! cell border, molecule (ang)
    logical(c_bool), save :: molcubic = .false. ! cell for molecule is cubic

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! crystal or molecule
    ldum = iw_radiobutton("Crystal",bool=ismolecule,boolval=.false.)
    call iw_tooltip("The new structure will be a periodic crystal",ttshown)
    ldum = iw_radiobutton("Molecule",bool=ismolecule,boolval=.true.,sameline=.true.)
    call iw_tooltip("The new structure will be a molecule",ttshown)

    ! name
    call iw_text("Name",highlight=.true.)
    str2 = "##name"
    ldum = igInputText(c_loc(str2),c_loc(namebuf),maxnamebuf-1,ImGuiInputTextFlags_None,c_null_ptr,c_null_ptr)

    if (.not.ismolecule) then
       ! symmetry
       call iw_text("Symmetry",highlight=.true.)

       ldum = iw_radiobutton("Detect",int=symopt,intval=1_c_int)
       call iw_tooltip("Calculate the symmetry operations from the complete list of unit cell atoms",ttshown)
       ldum = iw_radiobutton("Space group",int=symopt,intval=2_c_int,sameline=.true.)
       call iw_tooltip("Choose a space group from the list",ttshown)
       ldum = iw_radiobutton("Manual",int=symopt,intval=3_c_int,sameline=.true.)
       call iw_tooltip("Enter symmetry operations in cif-file format (e.g. -1/2+x,z,-y)",ttshown)

       ! symmetry options
       if (symopt == 2) then
          call draw_spg_table(ispg_selected)
       elseif (symopt == 3) then
          ! manual input of symmetry operations
          call igGetContentRegionAvail(szavail)
          sz%x = szavail%x
          sz%y = 9 * igGetTextLineHeightWithSpacing()
          str = "##symopsmanual" // c_null_char
          ldum = igInputTextMultiline(c_loc(str),c_loc(symopbuf),maxsymopbuf,sz,&
             ImGuiInputTextFlags_AllowTabInput,c_null_ptr,c_null_ptr)
       end if

       ! lattice: header
       call iw_text("Lattice",highlight=.true.)

       ldum = iw_radiobutton("Cell parameters",int=cellopt,intval=1_c_int)
       call iw_tooltip("Lattice is calculated from the cell lengths and angles",ttshown)
       ldum = iw_radiobutton("Lattice vectors",int=cellopt,intval=2_c_int,sameline=.true.)
       call iw_tooltip("Lattice vectors are given explicitly",ttshown)

       ! units and scale factor only if lattice vectors
       if (cellopt == 2) then
          call iw_combo_simple("Units##cel","Bohr" // c_null_char // "Angstrom" // c_null_char,&
             iunitcel,sameline=.true.)
          call iw_tooltip("Units for the cell parameters/lattice vectors",ttshown)

          call igSameLine(0._c_float,-1._c_float)
          str = "Scale" // c_null_char
          stropt = "%.3f" // c_null_char
          call igPushItemWidth(iw_calcwidth(7,1))
          ldum = igInputFloat(c_loc(str),scale,0._c_float,0._c_float,&
             c_loc(stropt),ImGuiInputTextFlags_None)
          call iw_tooltip("Scale factor to multiply the lattice vectors",ttshown)
          call igPopItemWidth()
       end if

       ! lattice: body
       if (cellopt == 1) then
          ! cell lengths
          call iw_text("Cell lengths (Å): ")
          call igSameLine(0._c_float,-1._c_float)
          str = "##celllength3" // c_null_char
          stropt = "%.3f" // c_null_char
          call igPushItemWidth(iw_calcwidth(3*8,2))
          ldum = igInputFloat3(c_loc(str),aa,c_loc(stropt),ImGuiInputTextFlags_None)

          ! cell angles
          call iw_text("Cell angles (°):  ")
          call igSameLine(0._c_float,-1._c_float)
          str = "##cellangle3" // c_null_char
          stropt = "%.3f" // c_null_char
          call igPushItemWidth(iw_calcwidth(3*8,2))
          ldum = igInputFloat3(c_loc(str),bb,c_loc(stropt),ImGuiInputTextFlags_None)
       else
          ! lattice vectors
          call igGetContentRegionAvail(szavail)
          sz%x = szavail%x
          sz%y = 4 * igGetTextLineHeightWithSpacing()
          str = "##latvecmanual" // c_null_char
          ldum = igInputTextMultiline(c_loc(str),c_loc(latvecbuf),maxlatvecbuf,sz,&
             ImGuiInputTextFlags_AllowTabInput,c_null_ptr,c_null_ptr)
       end if
    end if ! .not.ismolecule

    ! atomic positions: header
    call iw_text("Atomic positions",highlight=.true.) ! molecule
    call iw_text("(?)",sameline=.true.)
    call iw_tooltip("Give the atomic positions for this system as:" // newline //&
       "  <Sy> <x> <y> <z>" // newline //&
       "where Sy is the atomic symbol and x,y,z are the atomic coordinates.")

    ! units, does not apply in the case of crystal/cell parameters
    if (ismolecule .or. cellopt == 2) then
       call igSameLine(0._c_float,-1._c_float)
       call igSetCursorPosX(igGetCursorPosX() + 2 * g%Style%ItemSpacing%x)
       if (ismolecule) then
          call iw_combo_simple("Units","Bohr" // c_null_char // "Angstrom" // c_null_char,&
             iunitat_m,sameline=.true.)
       else
          call iw_combo_simple("Units","Bohr" // c_null_char // "Angstrom" // c_null_char // "Fractional" // c_null_char,&
             iunitat_c,sameline=.true.)
       end if
       call iw_tooltip("Units for the atomic coordinates",ttshown)
    end if

    ! atomic positions: body
    call igGetContentRegionAvail(szavail)
    sz%x = szavail%x
    if (ismolecule) then
       sz%y = max(szavail%y - iw_calcheight(2,1) - g%Style%ItemSpacing%y,igGetTextLineHeightWithSpacing())
    else
       sz%y = max(szavail%y - iw_calcheight(1,0) - g%Style%ItemSpacing%y,igGetTextLineHeightWithSpacing())
    end if
    str = "##atomicpositions" // c_null_char
    ldum = igInputTextMultiline(c_loc(str),c_loc(atposbuf),maxatposbuf,sz,&
       ImGuiInputTextFlags_AllowTabInput,c_null_ptr,c_null_ptr)

    ! options line
    if (ismolecule) then
       call iw_text("Structure options",highlight=.true.)

       ! cell border
       call igIndent(0._c_float)
       str = "Cell border (Å)" // c_null_char
       stropt = "%.3f" // c_null_char
       call igPushItemWidth(iw_calcwidth(7,1))
       ldum = igInputFloat(c_loc(str),rborder,0._c_float,0._c_float,&
          c_loc(stropt),ImGuiInputTextFlags_None)
       call iw_tooltip("Periodic cell border around new molecules",ttshown)
       call igPopItemWidth()
       call igSameLine(0._c_float,-1._c_float)

       ! cubic cell
       call igSetCursorPosX(igGetCursorPosX() + 2 * g%Style%ItemSpacing%x)
       str = "Cubic cell" // c_null_char
       ldum = igCheckbox(c_loc(str),molcubic)
       call iw_tooltip("Read new molecules inside a cubic periodic cell",ttshown)
       call igUnindent(0._c_float)
    end if

    ! right-align for the rest of the contents
    call igSetCursorPosX(iw_calcwidth(8,2,from_end=.true.))

    ! final buttons: ok
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
    ok = ok .or. iw_button("OK")
    if (ok) then
       ! build the input
       lu = fopen_scratch("formatted")

       ! symmetry and cell
       if (.not.ismolecule) then
          ! symmetry
          if (symopt == 2) then
             ! pass the spg string
             write (lu,'("spg ",A)') string(ispg_selected)
          elseif (symopt == 3) then
             ! pass the symm keywords
             call buffer_to_string_array(symopbuf,lu,prefix="symm ")
          end if

          ! cell
          if (cellopt == 1) then
             ! cell parameters
             write (lu,'("cell ",6(A," "),"ang")') (string(aa(i),'f',decimal=10),i=1,3), &
                (string(bb(i),'f',decimal=10),i=1,3)
          else
             ! lattice vectors
             write (lu,'("cartesian ",A)') string(scale,'f',decimal=10)
             if (iunitcel == 1) then
                write (lu,'("ang")')
             else
                write (lu,'("bohr")')
             end if
             call buffer_to_string_array(latvecbuf,lu)
             write (lu,'("endcartesian")')
          end if
       end if

       ! atomic positions
       if (ismolecule) then
          iunitat = iunitat_m
       else
          iunitat = iunitat_c
       end if
       if (iunitat == 0) then
          call buffer_to_string_array(atposbuf,lu,suffix=" bohr")
       elseif (iunitat == 1) then
          call buffer_to_string_array(atposbuf,lu,suffix=" ang")
       else
          call buffer_to_string_array(atposbuf,lu)
       end if

       ! molecular options
       if (ismolecule) then
          if (molcubic) write (lu,'("cubic")')
          write (lu,'("border ",A)') string(rborder/bohrtoa,'f',decimal=10)
       end if

       write (lu,'("end")')

       ! generate the seed
       rewind(lu)
       if (allocated(seed_)) deallocate(seed_)
       allocate(seed_(1))
       if (ismolecule) then
          call seed_(1)%parse_molecule_env(lu,ok)
       else
          call seed_(1)%parse_crystal_env(lu,ok)
       end if
       call fclose(lu)

       ! load the system and initialize
       if (ok.and.seed_(1)%isused) then
          idx = index(namebuf,c_null_char)
          seed_(1)%name = namebuf(1:idx-1)
          call add_systems_from_seeds(1,seed_)
          call launch_initialization_thread()
          doquit = .true.
       else
          deallocate(seed_)
       end if
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) &
       doquit = .true.

    ! quit the window
    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
      ttshown = .false.
      ismolecule = .false.
      iunitat_c = 2
      iunitat_m = 1
      symopt = 1
      ispg_selected = 1
      cellopt = 1
      iunitcel = 0
      scale = 1._c_float
      aa = 10._c_float
      bb = 90._c_float
      rborder = real(rborder_def*bohrtoa,c_float)
      molcubic = .false.

      if (allocated(atposbuf)) deallocate(atposbuf)
      allocate(character(len=maxatposbuf+1) :: atposbuf)
      atposbuf(1:26) = "H 0.00000 0.00000 0.00000" // c_null_char

      if (allocated(symopbuf)) deallocate(symopbuf)
      allocate(character(len=maxsymopbuf+1) :: symopbuf)
      symopbuf(1:70) = "## Each line is a symmetry operation of the type: " // newline //&
         "## -1/2+y,z,x" // c_null_char

      if (allocated(latvecbuf)) deallocate(latvecbuf)
      allocate(character(len=maxlatvecbuf+1) :: latvecbuf)
      latvecbuf(1:50) = "## Enter the three lattice vectors line by line" // newline // c_null_char

      if (allocated(namebuf)) deallocate(namebuf)
      allocate(character(len=maxnamebuf+1) :: namebuf)
      namebuf(1:20) = "Input structure" // c_null_char

    end subroutine init_state

    ! terminate the state for this window
    subroutine end_state()
      if (allocated(atposbuf)) deallocate(atposbuf)
      if (allocated(symopbuf)) deallocate(symopbuf)
      if (allocated(latvecbuf)) deallocate(latvecbuf)
      if (allocated(namebuf)) deallocate(namebuf)
    end subroutine end_state

  end subroutine draw_new_struct

  !> Draw the contents of the new structure from library window.
  module subroutine draw_new_struct_from_library(w)
    use gui_keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG
    use gui_main, only: g, add_systems_from_seeds,&
       launch_initialization_thread, system_shorten_names
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text, iw_calcheight,&
       iw_calcwidth, iw_radiobutton
    use crystalseedmod, only: crystalseed, realloc_crystalseed
    use spglib, only: SpglibSpaceGroupType, spg_get_spacegroup_type
    use global, only: clib_file, mlib_file, rborder_def
    use tools_io, only: string, fopen_scratch, fclose, stripchar, deblank
    use types, only: realloc
    use param, only: bohrtoa
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, stropt
    logical(c_bool) :: ldum, doquit
    logical :: ok, saveismol, doubleclicked, isset
    type(ImVec2) :: szero, sz
    integer :: i, nseed, oid
    type(crystalseed) :: seed
    type(crystalseed), allocatable :: seed_(:)

    ! window state
    logical, save :: ttshown = .false. ! tooltip flag
    logical, save :: ismolecule = .false. ! molecule/crystal choice at the top
    integer(c_int), save :: idopenlibfile = 0 ! the ID for the open library file
    integer, save :: nst = 0 ! number of library structures
    type(vstring), allocatable, save :: st(:) ! library structures
    logical(c_bool), allocatable, save :: lst(:) ! selected structures (1->nst)
    integer, save :: lastselected = 0 ! last selected structure (for shift/ctrl selection)
    real(c_float), save :: rborder = real(rborder_def*bohrtoa,c_float)
    logical(c_bool), save :: molcubic = .false.

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! check if we have info from the open library file window when it
    ! closes and recover it
    call update_window_id(idopenlibfile,oid)
    if (oid /= 0) then
       w%okfile_read = win(oid)%okfile_read
       if (w%okfile_read) then
          w%okfile_set = win(oid)%okfile_set
          w%okfile = win(oid)%okfile
       end if
    end if

    ! crystal or molecule, from library
    saveismol = ismolecule
    if (iw_radiobutton("Crystal",bool=ismolecule,boolval=.false.)) then
       if (saveismol.and..not.w%okfile_set) w%okfile = trim(clib_file)
       w%okfile_read = .true.
    end if
    call iw_tooltip("The new structure will be a periodic crystal",ttshown)
    if (iw_radiobutton("Molecule",bool=ismolecule,boolval=.true.,sameline=.true.)) then
       if (.not.saveismol.and..not.w%okfile_set) w%okfile = trim(mlib_file)
       w%okfile_read = .true.
    end if
    call iw_tooltip("The new structure will be a molecule",ttshown)

    ! render the rest of the window
    ! library file
    call iw_text("Source",highlight=.true.)
    if (iw_button("Library file")) &
       idopenlibfile = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openlibraryfile)
    call iw_tooltip("Library file from where the structures are read",ttshown)
    call iw_text(w%okfile,sameline=.true.)

    ! list box
    call iw_text("Structures to load from the library file",highlight=.true.)
    str = "##listbox" // c_null_char
    call igGetContentRegionAvail(sz)
    if (ismolecule) then
       sz%y = max(sz%y - iw_calcheight(2,1) - g%Style%ItemSpacing%y,igGetTextLineHeightWithSpacing())
    else
       sz%y = max(sz%y - iw_calcheight(1,0) - g%Style%ItemSpacing%y,igGetTextLineHeightWithSpacing())
    end if

    doubleclicked = .false.
    isset = .false.
    if (igBeginListBox(c_loc(str),sz)) then
       do i = 1, nst
          str = st(i)%s // c_null_char
          if (igSelectable_Bool(c_loc(str), lst(i), ImGuiSelectableFlags_AllowDoubleClick, szero)) then

             ! implement selection range with shift and control
             if (igIsKeyDown(ImGuiKey_ModShift).and.lastselected /= 0.and.lastselected /= i) then
                ! selecte a whole range
                lst = .false.
                if (lastselected > i) then
                   lst(i:lastselected) = .true.
                else
                   lst(lastselected:i) = .true.
                end if
             elseif (igIsKeyDown(ImGuiKey_ModCtrl)) then
                ! select and individual item and accumulate
                lst(i) = .true.
             else
                ! select one item, start range, remove all others
                lst = .false.
                lst(i) = .true.
                lastselected = i
                doubleclicked = igIsMouseDoubleClicked(ImGuiPopupFlags_MouseButtonLeft)
                isset = .true.
             end if
          end if
       end do
       call igEndListBox()
    end if

    if (ismolecule) then
       ! options line
       call iw_text("Structure options",highlight=.true.)

       ! cell border
       call igIndent(0._c_float)
       str = "Cell border (Å)" // c_null_char
       stropt = "%.3f" // c_null_char
       call igPushItemWidth(iw_calcwidth(7,1))
       ldum = igInputFloat(c_loc(str),rborder,0._c_float,0._c_float,&
          c_loc(stropt),ImGuiInputTextFlags_None)
       call iw_tooltip("Periodic cell border around new molecules",ttshown)
       call igPopItemWidth()
       call igSameLine(0._c_float,-1._c_float)

       ! cubic cell
       call igSetCursorPosX(igGetCursorPosX() + 2 * g%Style%ItemSpacing%x)
       str = "Cubic cell" // c_null_char
       ldum = igCheckbox(c_loc(str),molcubic)
       call iw_tooltip("Read new molecules inside cubic periodic cell",ttshown)
    end if

    ! right-align for the rest of the contents
    call igSetCursorPosX(iw_calcwidth(8,2,from_end=.true.))

    ! final buttons: ok
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
    ok = ok .or. doubleclicked
    ok = ok .or. iw_button("OK")
    ok = ok .and. (idopenlibfile == 0)
    if (ok) then
       nseed = count(lst(1:nst))
       if (nseed > 0) then
          ! we have systems (potentially), allocate the seed array
          allocate(seed_(nseed))
          nseed = 0
          do i = 1, nst
             if (lst(i)) then
                ! read the selected seeds from the library, check OK
                call seed%read_library(st(i)%s,ok,file=w%okfile)
                if (ok .and. ismolecule.eqv.seed%ismolecule) then
                   nseed = nseed + 1
                   seed_(nseed) = seed
                   seed_(nseed)%border = rborder/bohrtoa
                   seed_(nseed)%cubic = molcubic
                end if
             end if
          end do
          if (nseed > 0) then
             ! load the seeds as new systems
             call realloc_crystalseed(seed_,nseed)
             call add_systems_from_seeds(nseed,seed_)
             call launch_initialization_thread()
             call system_shorten_names()
          end if
       end if
       doquit = .true.
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) &
       doquit = .true.

    ! read the library file
    if (w%okfile_read) then
       call get_library_structure_list(w%okfile,nst,st,ismolecule)
       if (allocated(lst)) deallocate(lst)
       allocate(lst(nst))
       lst = .false.
       lastselected = 0
       w%okfile_read = .false.
    end if

    ! quit the window
    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
      ttshown = .false.
      ismolecule = .false.
      idopenlibfile = 0
      nst = 0
      lastselected = 0
      rborder = real(rborder_def*bohrtoa,c_float)
      molcubic = .false.
      w%okfile = trim(clib_file)
      w%okfile_set = .false.
      w%okfile_read = .true.
    end subroutine init_state

    ! terminate the state for this window
    subroutine end_state()
      nst = 0
      if (allocated(st)) deallocate(st)
      if (allocated(lst)) deallocate(lst)
    end subroutine end_state

  end subroutine draw_new_struct_from_library

  !> Draw the contents of the load field window.
  module subroutine draw_load_field(w)
    use fieldseedmod, only: field_detect_format
    use param, only: ifformat_wien, ifformat_elk, ifformat_pi,&
       ifformat_cube, ifformat_bincube, ifformat_abinit, ifformat_vasp,&
       ifformat_vaspnov, ifformat_qub, ifformat_xsf, ifformat_elkgrid,&
       ifformat_siestagrid, ifformat_dftb, ifformat_pwc, ifformat_wfn,&
       ifformat_wfx, ifformat_fchk, ifformat_molden, dirsep
    use gui_keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG
    use gui_main, only: nsys, sysc, sys, sys_init, g
    use gui_utils, only: iw_text, iw_tooltip, iw_radiobutton, iw_button,&
       iw_calcwidth
    use tools_io, only: string, uout
    class(window), intent(inout), target :: w

    logical :: oksys, ok, doquit, disabled
    integer :: isys, i, j, oid, ll, idx, iff
    type(ImVec2) :: szavail, sz, szero
    real(c_float) :: combowidth
    logical(c_bool) :: is_selected, ldum
    character(len=:,kind=c_char), allocatable, target :: str1, str2, loadstr, errmsg
    logical :: isgrid

    ! window state
    logical, save :: ttshown = .false. ! tooltip flag
    integer(c_int), save :: sourceopt = 0_c_int ! 0 = file, 1 = expression
    integer(c_int), save :: idopenfile1 = 0 ! the ID for the open field file
    integer(c_int), save :: idopenfile2 = 0 ! the ID for the first auxiliary file
    integer(c_int), save :: idopenfile3 = 0 ! the ID for the second auxiliary file
    character(len=:,kind=c_char), allocatable, target, save :: file1 ! first (main) file
    character(len=:,kind=c_char), allocatable, target, save :: file1_fmtstr ! first file format string
    logical, save :: file1_read = .false.
    integer, save :: file1_format = 0
    character(len=:,kind=c_char), allocatable, target, save :: file2 ! first auxiliary file
    logical, save :: file2_set = .false.
    character(len=:,kind=c_char), allocatable, target, save :: file3 ! second auxiliaryy file
    logical, save :: file3_set = .false.
    integer, save :: iginterp = 3 ! 0 = nearest, 1 = trilinear, 2 = trispline, 3 = tricubic, 4 = smoothrho

    ! permutation for the field format list (see dialog_user_callback)
    integer, parameter :: ifperm(0:17) = (/0,6,9,5,4,13,11,2,15,16,17,18,14,12,8,7,1,10/)

    ! DFTB+ hsd name candidates
    character*15, parameter :: hsdnames(4) = (/&
       "wfc-3ob-3-1.hsd","wfc.3ob-3-1.hsd","wfc.pbc-0-3.hsd","wfc.mio-1-1.hsd"/)

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    call update_window_id(idopenfile1,oid)
    if (oid /= 0) then
       if (win(oid)%okfile_set) then
          file1_read = .true.
          file1 = win(oid)%okfile
          file1_format = ifperm(win(oid)%dialog_data%isformat)
       end if
    end if
    call update_window_id(idopenfile2,oid)
    if (oid /= 0) then
       if (win(oid)%okfile_set) file2 = win(oid)%okfile
       file2_set = .true.
    end if
    call update_window_id(idopenfile3,oid)
    if (oid /= 0) then
       if (win(oid)%okfile_set) file3 = win(oid)%okfile
       file3_set = .true.
    end if

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    !! make sure we get a system that exists
    ! check if the system still exists
    isys = w%loadfield_isys
    oksys = system_ok(isys)
    if (.not.oksys) then
       ! reset to the table selected
       isys = win(iwin_tree)%table_selected
       oksys = system_ok(isys)
    end if
    if (.not.oksys) then
       ! reset to the input console selected
       isys = win(iwin_tree)%inpcon_selected
       oksys = system_ok(isys)
    end if
    if (.not.oksys) then
       ! check all systems and try to find one that is OK
       do i = 1, nsys
          oksys = system_ok(isys)
          if (oksys) then
             isys = i
             exit
          end if
       end do
    end if
    if (.not.oksys) then
       ! this dialog does not make sense anymore, close it and exit
       call end_state()
       call w%end()
       return
    end if

    ! system combo
    call iw_text("System")
    call igSameLine(0._c_float,-1._c_float)
    call igGetContentRegionAvail(szavail)
    combowidth = max(szavail%x - sz%x - g%Style%ItemSpacing%x,0._c_float)
    str1 = "##systemcombo" // c_null_char
    call igSetNextItemWidth(combowidth)
    if (oksys) &
       str2 = string(isys) // ": " // trim(sysc(isys)%seed%name) // c_null_char
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nsys
          if (sysc(i)%status == sys_init) then
             is_selected = (isys == i)
             str2 = string(i) // ": " // trim(sysc(i)%seed%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) then
                w%loadfield_isys = i
                isys = w%loadfield_isys
             end if
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Load a field for this system",ttshown)

    ! select the source
    ldum = iw_radiobutton("From File",int=sourceopt,intval=0_c_int)
    call iw_tooltip("Load the field from an external file",ttshown)
    ! ldum = iw_radiobutton("Expression",int=sourceopt,intval=1_c_int,sameline=.true.)
    ! call iw_tooltip("Load the field from an arithmetic expression involving other fields",ttshown)

    if (sourceopt == 0) then
       ! if a new file is read, detect the format if necessary
       if (file1_read) then
          if (file1_format == 0) &
             file1_format = -field_detect_format(file=file1)
          if (abs(file1_format) == ifformat_pi) file1_format = 0 ! aiPI deactivated for now
          if (file1_format == 0) then
             write (uout,'("!! Warning !! : Unknown field file extension")')
             file1 = ""
          end if
          file2 = ""
          file2_set = .false.
          file3 = ""
          file3_set = .false.
          file1_read = .false.
       end if

       ! source file and format
       call iw_text("Source",highlight=.true.)
       if (iw_button("File",danger=(file1_format==0))) &
          idopenfile1 = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openfieldfile)
       call iw_tooltip("File from where the field is read",ttshown)
       call iw_text(file1,sameline=.true.)

       isgrid = .false.
       if (file1_format /= 0) then
          call iw_text("Format: ")
          select case (abs(file1_format))
          case(ifformat_wien)
             file1_fmtstr = "WIEN"
             call iw_text("WIEN2k clmsum-style file (LAPW)",sameline=.true.)
             ! default struct file
             if (.not.file2_set) then
                ll = len_trim(sysc(isys)%seed%file)
                if (sysc(isys)%seed%file(ll-6:ll) == ".struct") &
                   file2 = sysc(isys)%seed%file
                file2_set = .true.
             end if
          case(ifformat_elk)
             file1_fmtstr = "ELK"
             call iw_text("elk STATE.OUT file (LAPW)",sameline=.true.)
             ! default GEOMETRY.OUT file
             if (.not.file2_set) then
                ll = len_trim(sysc(isys)%seed%file)
                if (sysc(isys)%seed%file(ll-3:ll) == ".OUT") &
                   file2 = sysc(isys)%seed%file
                file2_set = .true.
             end if
          case(ifformat_pi)
             file1_fmtstr = "PI"
             call iw_text("aiPI ion files (LCAO)",sameline=.true.)
          case(ifformat_cube)
             file1_fmtstr = "CUBE"
             call iw_text("Cube file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_bincube)
             file1_fmtstr = "BINCUBE"
             call iw_text("Binary cube file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_abinit)
             file1_fmtstr = "ABINIT"
             call iw_text("Abinit DEN-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_vasp)
             file1_fmtstr = "VASP"
             call iw_text("VASP CHGCAR-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_vaspnov)
             file1_fmtstr = "VASPNOV"
             call iw_text("VASP ELFCAR-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_qub)
             file1_fmtstr = "QUB"
             call iw_text("Aimpac qub file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_xsf)
             file1_fmtstr = "XSF"
             call iw_text("Xcrysden xsf-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_elkgrid)
             file1_fmtstr = "ELKGRID"
             call iw_text("elk grid file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_siestagrid)
             file1_fmtstr = "SIESTA"
             call iw_text("SIESTA RHO-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_dftb)
             file1_fmtstr = "DFTB"
             call iw_text("DFTB+ wavefunction (LCAO)",sameline=.true.)
             ! default eigenvec.bin file
             if (.not.file2_set) then
                idx = index(sysc(isys)%fullname,dirsep,.true.)
                file2 = sysc(isys)%fullname(1:idx) // "eigenvec.bin"
                inquire(file=file2,exist=ok)
                if (.not.ok) file2 = ""
                file2_set = .true.
             end if
             ! default hsd file (from the list of candidate names)
             if (.not.file3_set) then
                idx = index(sysc(isys)%fullname,dirsep,.true.)
                ok = .false.
                do j = 1, size(hsdnames,1)
                   file3 = sysc(isys)%fullname(1:idx) // trim(hsdnames(j))
                   inquire(file=file3,exist=ok)
                   if (ok) exit
                end do
                if (.not.ok) file3 = ""
                file3_set = .true.
             end if
          case(ifformat_pwc)
             file1_fmtstr = "PWC"
             call iw_text("Quantum ESPRESSO pwc file (grid+wfn)",sameline=.true.)
             isgrid = .true.
          case(ifformat_wfn)
             file1_fmtstr = "WFN"
             call iw_text("Gaussian wfn wavefunction file (molecular wfn)",sameline=.true.)
          case(ifformat_wfx)
             file1_fmtstr = "WFX"
             call iw_text("Gaussian wfx wavefunction file (molecular wfn)",sameline=.true.)
          case(ifformat_fchk)
             file1_fmtstr = "FCHK"
             call iw_text("Gaussian fchk wavefunction file (molecular wfn+basis set)",sameline=.true.)
          case(ifformat_molden)
             file1_fmtstr = "MOLDEN"
             call iw_text("molden wavefunction file (molecular wfn+basis set)",sameline=.true.)
          end select
          if (file1_format < 0) &
             call iw_text(" [auto-detected]",sameline=.true.)
       end if

       ! format-specific options and additional files
       select case (abs(file1_format))
       case(ifformat_wien)
          if (iw_button("File (.struct)",danger=.not.file2_set)) &
             idopenfile2 = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal)
          call iw_tooltip("WIEN2k structure file (.struct) used to interpret the clmsum-style file",ttshown)
          call iw_text(file2,sameline=.true.)
       case(ifformat_elk)
          if (iw_button("File (GEOMETRY.OUT)",danger=.not.file2_set)) &
             idopenfile2 = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal)
          call iw_tooltip("GEOMETRY.OUT structure file used to interpret the STATE.OUT",ttshown)
          call iw_text(file2,sameline=.true.)
       case(ifformat_dftb)
          if (iw_button("File (.bin)",danger=.not.file2_set)) &
             idopenfile2 = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal)
          call iw_tooltip("eigenvec.bin file for reading the DFTB+ wavefunction",ttshown)
          call iw_text(file2,sameline=.true.)

          if (iw_button("File (.hsd)",danger=.not.file3_set)) &
             idopenfile3 = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal)
          call iw_tooltip("hsd file for reading the DFTB+ wavefunction",ttshown)
          call iw_text(file3,sameline=.true.)
       end select

    elseif (sourceopt == 1) then
       ! from expression
    end if

    if (sourceopt == 0 .and. isgrid) then
       call iw_text("Options",highlight=.true.)
       call iw_text("Interpolation method")
       call iw_tooltip("Choose the interpolation method for the grid",ttshown)
       ldum = iw_radiobutton("Nearest",int=iginterp,intval=0_c_int)
       call iw_tooltip("Value of the nearest grid point",ttshown)
       ldum = iw_radiobutton("Tri-linear",int=iginterp,intval=1_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional linear interpolation",ttshown)
       ldum = iw_radiobutton("Tri-spline",int=iginterp,intval=2_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional spline interpolation",ttshown)
       ldum = iw_radiobutton("Tri-cubic",int=iginterp,intval=3_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional cubic polynomial interpolation",ttshown)
       ldum = iw_radiobutton("Smoothrho",int=iginterp,intval=4_c_int,sameline=.true.)
       call iw_tooltip("Polyharmonic splines + smoothing, all-electron densities only (recommended for them)",ttshown)
    end if

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(8,2,from_end=.true.) - g%Style%ScrollbarSize)
    call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! calculated whether we have enough info to continue to ok
    disabled = (len(file1) == 0)
    if (.not.disabled) then
       iff = abs(file1_format)
       if (iff == ifformat_wien .or. iff == ifformat_elk .or. iff == ifformat_dftb) &
          disabled = disabled .or. (len(file2) == 0)
       if (iff == ifformat_dftb) &
          disabled = disabled .or. (len(file3) == 0)
    end if
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
    ok = ok .or. iw_button("OK",disabled=disabled)
    if (ok) then
       if (sourceopt == 0) then
          loadstr = file1_fmtstr // " " // trim(file1) // " " // trim(file2) // " " // trim(file3)
          loadstr = loadstr // " notestmt"
          call sys(isys)%load_field_string(loadstr,.false.,iff,errmsg)
          if (len_trim(errmsg) > 0) then
             write (uout,'("!! Warning !! Could not read field for system: ",A)') string(isys)
             write (uout,'("!! Warning !! Load string: ",A)') trim(loadstr)
             write (uout,'("!! Warning !! Error message: ",A)') trim(errmsg)
          else
             doquit = .true.
          end if
       elseif (sourceopt == 1) then
          ! todo
       end if
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) &
       doquit = .true.

    ! quit the window
    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
      sourceopt = 0_c_int
      idopenfile1 = 0_c_int
      idopenfile2 = 0_c_int
      idopenfile3 = 0_c_int
      file1 = ""
      file1_read = .false.
      file1_format = 0
      file1_fmtstr = ""
      file2 = ""
      file2_set = .false.
      file3 = ""
      file3_set = .false.
    end subroutine init_state
    ! terminate the state for this window
    subroutine end_state()
      if (allocated(file1)) deallocate(file1)
      if (allocated(file1_fmtstr)) deallocate(file1_fmtstr)
      if (allocated(file2)) deallocate(file2)
      if (allocated(file3)) deallocate(file3)
    end subroutine end_state

    ! returns OK if system isys is initialized
    function system_ok(isys)
      integer, intent(in) :: isys
      logical :: system_ok

      system_ok = (isys > 0 .and. isys <= nsys)
      if (system_ok) system_ok = (sysc(isys)%status == sys_init)

    end function system_ok
  end subroutine draw_load_field

  !> Draw the SCF plot window.
  module subroutine draw_scfplot(w)
    use gui_utils, only: iw_text
    use tools_io, only: string
    class(window), intent(inout), target :: w

    ! xxxx ! end the window if the system is not initialized
    call iw_text("System: " // string(w%scfplot_isys))
    ! xxxx ! connect escape to close the dialog
    ! xxxx ! plot properties other than the energy?
    ! xxxx ! w%scfplot_isys = 1

    ! ! xxxx
    ! real(c_double), target, save :: x(10), y(10)
    ! type(ImVec2) :: sz
    ! type(ImVec4) :: auto
    ! character(kind=c_char,len=:), allocatable, target :: str1, str2

    ! strc = "Bleh" // c_null_char
    ! if (igBegin(c_loc(strc),show_implot_demo_window,ImGuiWindowFlags_None)) then
    !    do i = 1, 10
    !       x(i) = i
    !       y(i) = i * i
    !    end do
    !    strc = "Title" // c_null_char
    !    sz%x = 500_c_float
    !    sz%y = 500_c_float
    !    if (ipBeginPlot(c_loc(strc),sz,ImPlotFlags_None)) then
    !       str1 = "x" // c_null_char
    !       str2 = "y" // c_null_char
    !       call ipSetupAxes(c_loc(str1),c_loc(str2),ImPlotAxisFlags_None,ImPlotAxisFlags_None)
    !       str1 = "bleh" // c_null_char
    !       auto%x = 0._c_float
    !       auto%y = 0._c_float
    !       auto%z = 0._c_float
    !       auto%w = -1._c_float
    !       call ipSetNextMarkerStyle(ImPlotMarker_Circle,-1._c_float,auto,-1._c_float,auto)
    !       call ipPlotLine(c_loc(str1),c_loc(x),c_loc(y),10_c_int,ImPlotLineFlags_None,0_c_int)
    !       call ipEndPlot()
    !    end if
    !    call igEnd()
    ! end if

  end subroutine draw_scfplot

  !xx! private procedures

  ! Return the string for the tooltip shown by the tree window,
  ! corresponding to system i.
  function tree_system_tooltip_string(i) result(str)
    use crystalmod, only: pointgroup_info, holo_string
    use gui_main, only: sys, sysc, nsys, sys_init
    use tools_io, only: string
    use param, only: bohrtoa, maxzat, atmass, pcamu, bohr2cm, newline
    use tools_math, only: gcd
    integer, intent(in) :: i
    character(kind=c_char,len=:), allocatable, target :: str

    integer, allocatable :: nis(:)
    integer :: k, iz
    real*8 :: maxdv, mass, dens
    integer :: nelec
    character(len=3) :: schpg
    integer :: holo, laue
    integer :: izp0

    str = ""
    if (i < 1 .or. i > nsys) return

    ! file
    str = "||" // trim(sysc(i)%seed%name) // "||" // newline
    if (sysc(i)%status == sys_init) then
       ! file and system type
       str = str // trim(sysc(i)%seed%file) // newline
       if (sys(i)%c%ismolecule) then
          if (sys(i)%c%nmol == 1) then
             str = str // "A molecule." // newline
          else
             str = str // "A molecular cluster with " // string(sys(i)%c%nmol) // " fragments." //&
                newline
          end if
       elseif (sys(i)%c%ismol3d .or. sys(i)%c%nlvac == 3) then
          str = str // "A molecular crystal with Z=" // string(sys(i)%c%nmol)
          if (sys(i)%c%spgavail) then
             izp0 = 0
             do k = 1, sys(i)%c%nmol
                if (sys(i)%c%idxmol(k) < 0) then
                   izp0 = -1
                   exit
                elseif (sys(i)%c%idxmol(k) == 0) then
                   izp0 = izp0 + 1
                end if
             end do
             if (izp0 > 0) then
                str = str // " and Z'=" // string(izp0)
             else
                str = str // " and Z'<1"
             end if
          end if
          str = str // "." // newline
       elseif (sys(i)%c%nlvac == 2) then
          str = str // "A 1D periodic (polymer) structure." //&
             newline
       elseif (sys(i)%c%nlvac == 1) then
          str = str // "A 2D periodic (layered) structure." //&
             newline
       else
          str = str // "A crystal." //&
             newline
       end if
       str = str // newline

       ! number of atoms, electrons, molar mass
       str = str // string(sys(i)%c%ncel) // " atoms, " //&
          string(sys(i)%c%nneq) // " non-eq atoms, "//&
          string(sys(i)%c%nspc) // " species," // newline
       nelec = 0
       mass = 0d0
       do k = 1, sys(i)%c%nneq
          iz = sys(i)%c%spc(sys(i)%c%at(k)%is)%z
          if (iz >= maxzat .or. iz <= 0) cycle
          nelec = nelec + iz * sys(i)%c%at(k)%mult
          mass = mass + atmass(iz) * sys(i)%c%at(k)%mult
       end do
       str = str // string(nelec) // " electrons, " //&
          string(mass,'f',decimal=3) // " amu per cell" // newline
       ! empirical formula
       allocate(nis(sys(i)%c%nspc))
       nis = 0
       do k = 1, sys(i)%c%nneq
          nis(sys(i)%c%at(k)%is) = nis(sys(i)%c%at(k)%is) + sys(i)%c%at(k)%mult
       end do
       maxdv = gcd(nis,sys(i)%c%nspc)
       str = str // "Formula: "
       do k = 1, sys(i)%c%nspc
          str = str // string(sys(i)%c%spc(k)%name) // string(nint(nis(k)/maxdv)) // " "
       end do
       str = str // newline // newline

       if (.not.sys(i)%c%ismolecule) then
          ! cell parameters, volume, density, energy, pressure
          str = str // "a/b/c (Å): " // &
             string(sys(i)%c%aa(1)*bohrtoa,'f',decimal=4) // " " //&
             string(sys(i)%c%aa(2)*bohrtoa,'f',decimal=4) // " " //&
             string(sys(i)%c%aa(3)*bohrtoa,'f',decimal=4) //&
             newline
          str = str // "α/β/γ (°): " // &
             string(sys(i)%c%bb(1),'f',decimal=2) // " " //&
             string(sys(i)%c%bb(2),'f',decimal=2) // " " //&
             string(sys(i)%c%bb(3),'f',decimal=2) // " " //&
             newline
          str = str // "V (Å³): " // &
             string(sys(i)%c%omega*bohrtoa**3,'f',decimal=2) // newline
          dens = (mass*pcamu) / (sys(i)%c%omega*bohr2cm**3)
          str = str // "Density (g/cm³): " // string(dens,'f',decimal=3) // newline
          if (sysc(i)%seed%energy /= huge(1d0)) then
             str = str // "Energy (Ha): " // string(sysc(i)%seed%energy,'f',decimal=8) // newline
          end if
          if (sysc(i)%seed%pressure /= huge(1d0)) then
             str = str // "Pressure (GPa): " // string(sysc(i)%seed%pressure,'f',decimal=2) // newline
          end if
          str = str // newline

          ! symmetry
          if (sys(i)%c%spgavail) then
             call pointgroup_info(sys(i)%c%spg%pointgroup_symbol,schpg,holo,laue)
             str = str // "Symmetry: " // &
                string(sys(i)%c%spg%international_symbol) // " (" //&
                string(sys(i)%c%spg%spacegroup_number) // "), " //&
                string(holo_string(holo)) // "," //&
                newline

             str = str // string(sys(i)%c%neqv) // " symm-ops, " //&
                string(sys(i)%c%ncv) // " cent-vecs" //&
                newline
          else
             str = str // "Symmetry info not available" // newline
          end if
          str = str // newline
       end if

       ! number of scalar fields
       str = str // string(sys(i)%nf) // " scalar fields loaded" // newline
    else
       ! not initialized
       str = str // "Not initialized" // newline
    end if
    str = str // newline // "[Right-click for options]" // newline
    str = str // c_null_char

  end function tree_system_tooltip_string

  ! Return the string for the tooltip shown by the tree window,
  ! corresponding to system si and field fj
  function tree_field_tooltip_string(si,fj) result(str)
    use gui_main, only: sys, sysc, nsys, sys_init
    use fieldmod, only: field, type_uninit, type_promol, type_grid, type_wien,&
       type_elk, type_pi, type_wfn, type_dftb, type_promol_frag, type_ghost
    use wfn_private, only: molden_type_psi4, molden_type_orca, molden_type_adf_sto,&
       wfn_rhf, wfn_uhf, wfn_frac
    use grid3mod, only: mode_nearest, mode_trilinear, mode_trispline, mode_tricubic,&
       mode_smr
    use tools_io, only: string, nameguess
    use param, only: newline, maxzat0
    integer, intent(in) :: si, fj
    character(kind=c_char,len=:), allocatable, target :: str, aux

    integer :: i, nal
    type(field), pointer :: f

    str = ""
    if (si < 1 .or. si > nsys) return
    if (sysc(si)%status /= sys_init) return
    if (fj < 0 .or. fj > sys(si)%nf) return
    if (.not. sys(si)%f(fj)%isinit) return
    f => sys(si)%f(fj)

    ! file
    str = "||" // trim(f%name) // "||" // newline
    if (sys(si)%iref == fj) &
       str = str // "## Reference field for this system ##" // newline
    call sys(si)%aliasstring(fj,nal,aux)
    str = str // "Names: " // trim(adjustl(aux)) // newline

    ! type and type-specific info
    select case (f%type)
    case (type_uninit)
       str = str // "???" // newline // newline

    case (type_promol)
       str = str // "Promolecular density" // newline // newline

    case (type_grid)
       str = str // "Grid field, with " // string(f%grid%n(1)) // "x" // string(f%grid%n(2)) // "x" //&
          string(f%grid%n(3)) // " points" // newline
       if (f%grid%isqe) then
          str = str // "Plane-wave Kohn-Sham states available: spin=" // string(f%grid%qe%nspin) //&
             ", k-points=" // string(f%grid%qe%nks) // ", bands=" // string(f%grid%qe%nbnd) // newline
       end if
       if (f%grid%iswan) then
          str = str // "Wannier function info available: " // string(f%grid%qe%nk(1)) // "x" //&
             string(f%grid%qe%nk(2)) // "x" // string(f%grid%qe%nk(3)) // newline
       end if
       str = str // newline

       str = str // "Longest Voronoi vector (bohr): " // string(f%grid%dmax,'f',decimal=5) // newline
       str = str // "Grid integral: " // &
          string(sum(f%grid%f) * f%c%omega / real(product(f%grid%n),4),'f',decimal=8) // newline
       str = str // "Minimum value: " // string(minval(f%grid%f),'e',decimal=4) // newline
       str = str // "Average: " // string(sum(f%grid%f) / real(product(f%grid%n),4),'e',decimal=8) // newline
       str = str // "Maximum value: " // string(maxval(f%grid%f),'e',decimal=4) // newline
       str = str // "Interpolation mode: "
       select case (f%grid%mode)
       case(mode_nearest)
          str = str // "nearest grid node value" // newline
       case(mode_trilinear)
          str = str // "tri-linear" // newline
       case(mode_trispline)
          str = str // "tri-spline" // newline
       case(mode_tricubic)
          str = str // "tri-cubic" // newline
       case(mode_smr)
          str = str // "smooth all-electron density with " // string(f%grid%smr_nenv) //&
             " stencil nodes and " // string(f%grid%smr_fdmax,'f',decimal=4) //&
             " smoothing dmax factor" // newline
       end select

    case (type_wien)
       str = str // "WIEN2k spherical harmonic + plane waves" // newline // newline

       str = str // "Complex: " // string(f%wien%cmpl) // newline
       str = str // "Spherical harmonics expansion LMmax: " // string(size(f%wien%lm,2)) // newline
       str = str // "Max. points in radial grid: " // string(size(f%wien%slm,1)) // newline
       str = str // "Total number of plane waves: " // string(f%wien%nwav) // newline
       str = str // "Density-style normalization: " // string(f%wien%cnorm) // newline

    case (type_elk)
       str = str // "Elk spherical harmonic + plane waves" // newline // newline

       str = str // "Number of LM pairs: " // string(size(f%elk%rhomt,2)) // newline
       str = str // "Max. points in radial grid: " // string(size(f%elk%rhomt,1)) // newline
       str = str // "Total number of plane waves: " // string(size(f%elk%rhok)) // newline

    case (type_pi)
       str = str // "aiPI atomic radial functions" // newline // newline

       str = str // "Exact calculation: " // string(f%exact) // newline
       do i = 1, f%c%nspc
          str = str // string(f%c%spc(f%c%at(i)%is)%name,length=5) //&
             ": " // string(f%pi%bas(i)%piname) // newline
       end do

    case (type_wfn)
       str = str // "Molecular wavefunction "
       if (f%wfn%issto) then
          str = str // "(Slater-type orbitals)"
       else
          str = str // "(Gaussian-type orbitals)"
       end if
       str = str // newline // newline

       if (f%wfn%molden_type == molden_type_orca) then
          str = str // "Molden file dialect: orca" // newline
       else if (f%wfn%molden_type == molden_type_psi4) then
          str = str // "Molden file dialect: psi4" // newline
       else if (f%wfn%molden_type == molden_type_adf_sto) then
          str = str // "Molden file dialect: ADF (STOs)"  // newline
       end if
       if (f%wfn%wfntyp == wfn_rhf) then
          str = str // "Wavefunction type: restricted" // newline
          str = str // "Number of MOs (total): " // string(f%wfn%nmoall) // newline
          str = str // "Number of MOs (occupied): " // string(f%wfn%nmoocc) // newline
       elseif (f%wfn%wfntyp == wfn_uhf) then
          str = str // "Wavefunction type: unrestricted" // newline
          str = str // "Number of MOs (total): " // string(f%wfn%nmoall) //&
             " (alpha=" // string(f%wfn%nalpha+f%wfn%nalpha_virt) // ",beta=" //&
             string(f%wfn%nmoall-(f%wfn%nalpha+f%wfn%nalpha_virt)) // newline
          str = str // "Number of MOs (occupied): " // string(f%wfn%nmoocc) //&
             " (alpha=" // string(f%wfn%nalpha) // ",beta=" // string(f%wfn%nmoocc-f%wfn%nalpha) // newline
       elseif (f%wfn%wfntyp == wfn_frac) then
          str = str // "Wavefunction type: fractional occupation" // newline
          str = str // "Number of MOs: " // string(f%wfn%nmoocc) // newline
          str = str // "Number of electrons: " // string(nint(sum(f%wfn%occ(1:f%wfn%nmoocc)))) // newline
       end if
       str = str // "Number of primitives: " // string(f%wfn%npri) // newline
       str = str // "Number of EDFs: " // string(f%wfn%nedf) // newline

    case (type_dftb)
       str = str // "DFTB+, linear combination of atomic orbitals" // newline // newline

       str = str // "Number of states: " // string(f%dftb%nstates) // newline
       str = str // "Number of spin channels: " // string(f%dftb%nspin) // newline
       str = str // "Number of orbitals: " // string(f%dftb%norb) // newline
       str = str // "Number of kpoints: " // string(f%dftb%nkpt) // newline
       str = str // "Real wavefunction? " // string(f%dftb%isreal) // newline
       str = str // "Exact calculation? " // string(f%exact) // newline

    case (type_promol_frag)
       str = str // "Promolecular with fragment" // newline // newline

       str = str // "Number of atoms in fragment: " // string(f%fr%nat) // newline

    case (type_ghost)
       str = str // "Ghost field" // newline // newline

       str = str // "Expression: " // string(f%expr)

    case default
       str = str // "???" // newline // newline
    end select
    str = str // "Use core densities: " // string(f%usecore) // newline
    if (any(f%zpsp > 0)) then
       str = str // "Core charges (ZPSP): "
       do i = 1, maxzat0
          if (f%zpsp(i) > 0) then
             str = str // string(nameguess(i,.true.)) // "(" // string(f%zpsp(i)) // ") "
          end if
       end do
       str = str // newline
    end if
    str = str // "Numerical derivatives: " // string(f%numerical) // newline // newline

    ! critical points
    str = str // "Non-equivalent critical points: " // string(f%ncp) // newline
    str = str // "Cell critical points: " // string(f%ncpcel) // newline

    ! last message
    str = str // newline // "[Right-click for options]" // newline
    str = str // c_null_char

  end function tree_field_tooltip_string

  ! the callback for the right-hand-side pane of the dialog
  subroutine dialog_user_callback(vFilter, vUserData, vCantContinue) bind(c)
    use gui_main, only: g
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_text, iw_radiobutton,&
       iw_combo_simple
    use gui_interfaces_cimgui
    use tools_io, only: string
    type(c_ptr), intent(in), value :: vFilter ! const char *
    type(c_ptr), value :: vUserData ! void *
    logical(c_bool) :: vCantContinue ! bool *

    character(kind=c_char,len=:), allocatable, target :: str, stropt, strex
    type(dialog_userdata), pointer :: data
    logical(c_bool) :: ldum
    type(ImVec2) :: sz
    integer(c_int) :: flags

    logical, save :: ttshown = .false. ! tooltip flag

    ! generate the data pointer
    call c_f_pointer(vUserData,data)

    !! common options !!

    ! header
    call iw_text("Options",highlight=.true.)

    ! show hidden files
    str = "Show hidden files" // c_null_char
    if (igCheckbox(c_loc(str),data%showhidden)) then
       flags = IGFD_GetFlags(data%dptr)
       if (data%showhidden) then
          flags = iand(flags,not(ImGuiFileDialogFlags_DontShowHiddenFiles))
       else
          flags = ior(flags,ImGuiFileDialogFlags_DontShowHiddenFiles)
       end if
       call IGFD_SetFlags(data%dptr,flags)
    end if
    call iw_tooltip("Show the OS hidden files and directories in this dialog",ttshown)

    !! options specific to the open files dialog !!
    if (data%purpose == wpurp_dialog_openfiles) then
       str = "Read last structure only" // c_null_char
       ldum = igCheckbox(c_loc(str),data%readlastonly)
       call iw_tooltip("Read only the last structure in the file",ttshown)
       call igNewLine()

       ! radio buttons for auto/crystal/molecule
       call iw_text("Read structures as...",highlight=.true.)
       ldum = iw_radiobutton("Auto-detect",int=data%mol,intval=-1_c_int)
       call iw_tooltip("Auto-detect whether new structures are read as crystals or molecules",ttshown)
       ldum = iw_radiobutton("Crystal",int=data%mol,intval=0_c_int)
       call iw_tooltip("Force new structures to be read as crystals",ttshown)
       ldum = iw_radiobutton("Molecule",int=data%mol,intval=2_c_int)
       call iw_tooltip("Force new structures to be read as molecules",ttshown)

       ! molecular options
       call igIndent(0._c_float)
       str = "Cell border (Å)" // c_null_char
       stropt = "%.3f" // c_null_char
       strex = string(data%rborder,'f',decimal=3) // c_null_char
       call igCalcTextSize(sz,c_loc(strex),c_null_ptr,.false._c_bool,-1._c_float)
       call igPushItemWidth(sz%x + 2 * g%Style%FramePadding%x)
       ldum = igInputFloat(c_loc(str),data%rborder,0._c_float,0._c_float,&
          c_loc(stropt),ImGuiInputTextFlags_None)
       call iw_tooltip("Periodic cell border around new molecules",ttshown)
       call igPopItemWidth()

       str = "Cubic cell" // c_null_char
       ldum = igCheckbox(c_loc(str),data%molcubic)
       call iw_tooltip("Read new molecules inside cubic periodic cell",ttshown)
       call igUnindent(0._c_float)
       call igNewLine()

       ! Input structure format (isformat)
       call iw_text("Read file format",highlight=.true.)
       stropt = "" &
          // "Auto-detect" // c_null_char &             ! isformat_unknown = 0
          // "Abinit DEN-style file" // c_null_char &   ! isformat_abinit = 7
          // "Binary cube" // c_null_char &             ! isformat_bincube = 5
          // "CIF file" // c_null_char &                ! isformat_cif = 1
          // "CRYSTAL output" // c_null_char &          ! isformat_crystal = 11
          // "Cube" // c_null_char &                    ! isformat_cube = 4
          // "DFTB+ gen file" // c_null_char &          ! isformat_gen = 20
          // "DMACRYS .21 file" // c_null_char &        ! isformat_f21 = 3
          // "DMACRYS dmain file" // c_null_char &      ! isformat_dmain = 27
          // "elk GEOMETRY.OUT" // c_null_char &        ! isformat_elk = 8
          // "FHIaims input" // c_null_char &           ! isformat_aimsin = 28
          // "FHIaims output" // c_null_char &          ! isformat_aimsout = 29
          // "Gaussian fchk" // c_null_char &           ! isformat_fchk = 15
          // "Gaussian output" // c_null_char &         ! isformat_gaussian = 17
          // "Gaussian wfn" // c_null_char &            ! isformat_wfn = 13
          // "Gaussian wfx" // c_null_char &            ! isformat_wfx = 14
          // "ORCA molden file" // c_null_char &        ! isformat_orca = 26
          // "postg output" // c_null_char &            ! isformat_pgout = 25
          // "psi4 molden file" // c_null_char &        ! isformat_molden = 16
          // "psi4 output" // c_null_char &             ! isformat_dat = 24
          // "Quantum ESPRESSO input" // c_null_char &  ! isformat_qein = 9
          // "Quantum ESPRESSO output" // c_null_char & ! isformat_qeout = 10
          // "Quantum ESPRESSO pwc" // c_null_char &    ! isformat_pwc = 22
          // "SHELX" // c_null_char &                   ! isformat_shelx = 2
          // "SIESTA IN/OUT file" // c_null_char &      ! isformat_siesta = 18
          // "TINKER frac file" // c_null_char &        ! isformat_tinkerfrac = 30
          // "VASP" // c_null_char &                    ! isformat_vasp = 21
          // "WIEN2k struct file" // c_null_char &      ! isformat_struct = 6
          // "Xcrysden axsf file" // c_null_char &      ! isformat_axsf = 23
          // "Xcrysden xsf file" // c_null_char &       ! isformat_xsf = 19
          // "xyz file" // c_null_char                  ! isformat_xyz = 12
       call iw_combo_simple("##formatcombo",stropt,data%isformat)
       call iw_tooltip("Force the new structure to be read with this file format, or auto-detect&
          &from the extension",ttshown)
       call igNewLine()
    elseif (data%purpose == wpurp_dialog_openfieldfile) then
       ! Input field format (ifformat)
       call iw_text("Read file format",highlight=.true.)
       stropt = "" &
          // "Auto-detect" // c_null_char &              ! ifformat_unknown = 0
          // "ABINIT DEN-style file" // c_null_char &    ! ifformat_abinit = 6
          // "Aimpac qub" // c_null_char &               ! ifformat_qub = 9
          // "Binary cube" // c_null_char &              ! ifformat_bincube = 5
          // "Cube file" // c_null_char &                ! ifformat_cube = 4
          // "DFTB+ detailed.xml" // c_null_char &       ! ifformat_dftb = 13
          // "elk grid" // c_null_char &                 ! ifformat_elkgrid = 11
          // "elk STATE.OUT" // c_null_char &            ! ifformat_elk = 2
          // "Gaussian wfn" // c_null_char &             ! ifformat_wfn = 15
          // "Gaussian wfx" // c_null_char &             ! ifformat_wfx = 16
          // "Gaussian fchk" // c_null_char &            ! ifformat_fchk = 17
          // "Molden file" // c_null_char &              ! ifformat_molden = 18
          // "Quantum ESPRESSO pwc" // c_null_char &     ! ifformat_pwc = 14
          // "SIESTA RHO-style file" // c_null_char &    ! ifformat_siestagrid = 12
          // "VASP ELFCAR-style file" // c_null_char &   ! ifformat_vaspnov = 8
          // "VASP CHGCAR-style file" // c_null_char &   ! ifformat_vasp = 7
          // "WIEN2k clmsum-style file" // c_null_char & ! ifformat_wien = 1
          // "Xcrysden xsf" // c_null_char               ! ifformat_xsf = 10
       call iw_combo_simple("##formatcombo",stropt,data%isformat)
       call iw_tooltip("Force the new field to be read with the given file format, or auto-detect&
          &from the extension",ttshown)
       call igNewLine()
    end if

  end subroutine dialog_user_callback

  !> Get the structure list from the library file. If ismol = .true.,
  !> read only the moelcules. If ismol = .false., read only the
  !> crystals.
  subroutine get_library_structure_list(libfile,nst,st,ismol)
    use tools_io, only: fopen_read, fclose, lgetword, getword, equal, getline
    use types, only: vstring, realloc
    character(len=:), allocatable, intent(in) :: libfile
    integer, intent(out) :: nst
    type(vstring), allocatable, intent(inout) :: st(:)
    logical, intent(in) :: ismol

    integer :: lu, lp
    character(len=:), allocatable :: word, name, line
    logical :: ok

    ! open the file
    nst = 0
    inquire(file=libfile,exist=ok)
    if (.not.ok) return
    lu = fopen_read(libfile)
    if (lu < 0) return

    ! preallocate
    if (allocated(st)) deallocate(st)
    allocate(st(10))

    ! read the structures
    main: do while (getline(lu,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,'structure')) then
          name = getword(line,lp)

          ok = getline(lu,line)
          if (.not.ok) exit main
          lp = 1
          word = lgetword(line,lp)
          ok = equal(word,'crystal').and..not.ismol
          ok = ok .or. (equal(word,'molecule').and.ismol)
          if (.not.ok) cycle main

          nst = nst + 1
          if (nst > size(st,1)) call realloc(st,2*nst)
          st(nst)%s = name
       endif
    end do main
    if (nst > 0) &
       call realloc(st,nst)
    ! clean up
    call fclose(lu)

  end subroutine get_library_structure_list

  !> Draw a space group table. Entry ispg (Hall number) is selected.
  subroutine draw_spg_table(ispg)
    use gui_utils, only: iw_text
    use spglib, only: SpglibSpaceGroupType, spg_get_spacegroup_type
    use tools_io, only: ioj_left, string, deblank, stripchar
    integer, intent(inout) :: ispg

    integer :: i
    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: flags
    type(ImVec2) :: sz, szero
    type(SpglibSpaceGroupType) :: sa

    ! initialize
    szero%x = 0
    szero%y = 0

    ! space group table
    str = "##spacegrouptable" // c_null_char
    flags = ImGuiTableFlags_Borders
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_ScrollX)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)

    call igGetContentRegionAvail(sz)
    sz%y = 8 * igGetTextLineHeightWithSpacing()
    if (igBeginTable(c_loc(str),7,flags,sz,0._c_float)) then
       ! set up columns
       str = "Hall" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "ITA" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,1)

       str = "HM short" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,2)

       str = "HM long" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,3)

       str = "choice" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,4)

       str = "system" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,5)

       str = "Hall symbol" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,6)
       call igTableSetupScrollFreeze(0, 1) ! top row always visible
       call igTableHeadersRow()

       ! table body
       do i = 1, 530
          sa = spg_get_spacegroup_type(i)

          ! Hall
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float);
          if (igTableSetColumnIndex(0)) then
             str = string(i) // c_null_char
             flags = ImGuiSelectableFlags_SpanAllColumns
             flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
             if (igSelectable_Bool(c_loc(str),logical(ispg == i,c_bool),flags,szero)) &
                ispg = i
          end if

          ! ITA
          if (igTableNextColumn()) call iw_text(string(sa%number,5,ioj_left))

          ! HM-short
          if (igTableNextColumn()) then
             str = deblank(sa%international_short)
             call iw_text(trim(stripchar(str,"_")))
          end if

          ! HM-long
          if (igTableNextColumn()) then
             str = deblank(sa%international_full)
             call iw_text(trim(stripchar(str,"_")))
          end if

          ! choice
          if (igTableNextColumn()) call iw_text(string(sa%choice,6,ioj_left))

          ! system
          if (igTableNextColumn()) then
             if (sa%number >= 1 .and. sa%number <= 2) then
                str = "triclinic"
             elseif (sa%number >= 3 .and. sa%number <= 15) then
                str = "monoclinic"
             elseif (sa%number >= 16 .and. sa%number <= 74) then
                str = "orthorhombic"
             elseif (sa%number >= 75 .and. sa%number <= 142) then
                str = "tetragonal"
             elseif (sa%number >= 143 .and. sa%number <= 167) then
                str = "trigonal"
             elseif (sa%number >= 168 .and. sa%number <= 194) then
                str = "hexagonal"
             elseif (sa%number >= 195 .and. sa%number <= 230) then
                str = "cubic"
             else
                str = "??"
             end if
             call iw_text(str)
          end if

          ! hall symbol
          if (igTableNextColumn()) call iw_text(string(sa%hall_symbol))
       end do

       call igTableSetColumnWidthAutoAll(igGetCurrentTable())
       call igEndTable()
    end if

    ! current choice
    sa = spg_get_spacegroup_type(ispg)
    str = deblank(sa%international_full)
    str = "Current choice: " // trim(stripchar(str,"_"))
    if (len_trim(sa%choice) > 0) str = str // " (" // trim(sa%choice) // ")"
    call iw_text(str)

  end subroutine draw_spg_table

end submodule proc
