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
  implicit none

  ! Count unique IDs for keeping track of windows and widget
  integer :: idcount = 0

  ! column ids for the table in the tree widget
  integer(c_int), parameter :: ic_closebutton = 0
  integer(c_int), parameter :: ic_id = 1
  integer(c_int), parameter :: ic_name = 2
  integer(c_int), parameter :: ic_spg = 3
  integer(c_int), parameter :: ic_v = 4
  integer(c_int), parameter :: ic_nneq = 5
  integer(c_int), parameter :: ic_ncel = 6
  integer(c_int), parameter :: ic_nmol = 7
  integer(c_int), parameter :: ic_a = 8
  integer(c_int), parameter :: ic_b = 9
  integer(c_int), parameter :: ic_c = 10
  integer(c_int), parameter :: ic_alpha = 11
  integer(c_int), parameter :: ic_beta = 12
  integer(c_int), parameter :: ic_gamma = 13

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
  ! function tree_tooltip_string(i)
  ! subroutine dialog_user_callback(vFilter, vUserData, vCantContinue)

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

  !> Create a window in the window stack with the given type. Returns
  !> the window ID.
  module function stack_create_window(type,isopen,purpose)
    use gui_window, only: window, nwin, win
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in), optional :: purpose

    type(window), allocatable :: aux(:)
    integer :: stack_create_window

    integer :: i, id

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
    if (.not.allocated(win)) then
       allocate(win(nwin))
    elseif (nwin > size(win,1)) then
       allocate(aux(2*nwin))
       aux(1:size(win,1)) = win
       call move_alloc(aux,win)
    end if

    ! initialize the new window
    call win(id)%init(type,isopen,purpose)
    stack_create_window = id

  end function stack_create_window

  !> Initialize a window of the given type. If isiopen, initialize it
  !> as open.
  module subroutine window_init(w,type,isopen,purpose)
    use gui_main, only: ColorDialogDir, ColorDialogFile
    use tools_io, only: ferror, faterr
    use param, only: bohrtoa
    class(window), intent(inout) :: w
    integer, intent(in) :: type
    logical, intent(in) :: isopen
    integer, intent(in), optional :: purpose

    character(kind=c_char,len=:), allocatable, target :: str1

    ! initialization
    w%isinit = .true.
    w%isopen = isopen
    w%type = type
    w%id = -1
    w%name = ""
    if (allocated(w%iord)) deallocate(w%iord)
    w%dialog_data%ptr = c_null_ptr
    w%dialog_data%mol = -1
    w%dialog_data%showhidden = .false._c_bool
    w%dialog_data%isformat = isformat_unknown
    w%dialog_data%readlastonly = .false._c_bool
    w%dialog_data%purpose = wpurp_unknown
    w%dialog_data%molcubic = logical(.false.,c_bool)
    w%dialog_data%rborder = real(rborder_def*bohrtoa,c_float)
    w%dialog_purpose = wpurp_unknown

    ! type-specific initialization
    if (type == wintype_dialog) then
       ! create the dialog object and set up the style
       w%ptr = IGFD_Create()
       str1 = "+" // c_null_char
       call IGFD_SetFileStyle(w%ptr,IGFD_FileStyleByTypeDir,c_null_ptr,ColorDialogDir,c_loc(str1),c_null_ptr)
       str1 = " " // c_null_char
       call IGFD_SetFileStyle(w%ptr,IGFD_FileStyleByTypeFile,c_null_ptr,ColorDialogFile,c_loc(str1),c_null_ptr)
       if (.not.present(purpose)) &
          call ferror('window_init','dialog requires a purpose',faterr)
       w%dialog_purpose = purpose
    end if

  end subroutine window_init

  !> End a window and deallocate the data.
  module subroutine window_end(w)
    class(window), intent(inout) :: w

    ! window-specific destruction
    if (w%isinit .and. w%type == wintype_dialog .and. c_associated(w%ptr)) &
       call IGFD_Destroy(w%ptr)

    ! deallocate the rest of the data
    w%dialog_purpose = wpurp_unknown
    w%isinit = .false.
    w%isopen = .false.
    w%id = -1
    w%name = ""
    if (allocated(w%iord)) deallocate(w%iord)

  end subroutine window_end

  !> Draw an ImGui window.
  module subroutine window_draw(w)
    use tools_io, only: string, ferror, faterr
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: inisize

    if (.not.w%isinit) return
    if (.not.w%isopen) return

    ! First pass: assign ID, name, and flags
    if (w%id < 0) then
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
          w%dialog_data%ptr = w%ptr
          w%dialog_data%purpose = w%dialog_purpose
          w%flags = ImGuiFileDialogFlags_DontShowHiddenFiles
          str2 = "" // c_null_char ! default path
          inisize%x = 800._c_float
          inisize%y = 480._c_float
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)

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
             call IGFD_OpenPaneDialog2(w%ptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,0_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_savelogfile) then
             w%name = "Save Log File..." // c_null_char
             str1 = "All files (*.*){*.*}" // c_null_char
             call IGFD_OpenPaneDialog2(w%ptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile) then
             w%name = "Open Library File..." // c_null_char
             str1 = "All files (*.*){*.*}" // c_null_char
             call IGFD_OpenPaneDialog2(w%ptr,c_loc(w%name),c_loc(w%name),c_loc(str1),c_loc(str2),&
                c_funloc(dialog_user_callback),280._c_float,1_c_int,c_loc(w%dialog_data),w%flags)
          else
             call ferror('window_draw','unknown dialog purpose',faterr)
          end if
       elseif (w%type == wintype_new) then
          w%name = "New Structure..." // c_null_char
          w%flags = ImGuiWindowFlags_None
          inisize%x = 800._c_float
          inisize%y = 680._c_float
          call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)
       end if
    end if

    if (w%isopen) then
       if (w%type == wintype_dialog) then
          call w%draw_dialog()
       else
          if (igBegin(c_loc(w%name),w%isopen,w%flags)) then
             ! assign the pointer ID for the window, if not a dialog
             w%ptr = igGetCurrentWindow()

             ! draw the window contents, depending on type
             if (w%type == wintype_tree) then
                call w%draw_tree()
             elseif (w%type == wintype_view) then
                str1 = "Hello View!"
                call igText(c_loc(str1))
             elseif (w%type == wintype_console_input) then
                call w%draw_ci()
             elseif (w%type == wintype_console_output) then
                call w%draw_co()
             elseif (w%type == wintype_new) then
                call w%draw_new()
             end if
          end if
          call igEnd()
       end if
    end if

  end subroutine window_draw

  !> Draw the contents of a tree window
  module subroutine draw_tree(w)
    use gui_keybindings, only: is_bind_event, BIND_TREE_REMOVE_SYSTEM
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button
    use gui_main, only: nsys, sys, sysc, sys_empty, sys_init,&
       sys_loaded_not_init, sys_initializing, ColorTableCellBg_Mol,&
       ColorTableCellBg_MolClus, ColorTableCellBg_MolCrys, ColorTableCellBg_Crys3d,&
       ColorTableCellBg_Crys2d, ColorTableCellBg_Crys1d, launch_initialization_thread,&
       kill_initialization_thread, system_shorten_names, remove_system, tooltip_delay,&
       ColorDangerButton, g
    use tools_io, only: string
    use types, only: realloc
    use param, only: bohrtoa
    use c_interface_module
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, zeroc
    type(ImVec2) :: szero, sz, sztext, szavail
    integer(c_int) :: flags, color, idir
    integer :: i, j, k, nshown, newsel, jsel
    logical(c_bool) :: ldum
    type(c_ptr) :: ptrc
    type(ImGuiTableSortSpecs), pointer :: sortspecs
    type(ImGuiTableColumnSortSpecs), pointer :: colspecs
    logical :: hadenabledcolumn, buttonhovered_close, buttonhovered_expand, reinit, dopop
    type(c_ptr), save :: cfilter = c_null_ptr
    logical, save :: ttshown = .false.
    real(c_float) :: rshift

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

    ! text filter
    if (.not.c_associated(cfilter)) &
       cfilter = ImGuiTextFilter_ImGuiTextFilter(c_loc(zeroc))
    str = "##treefilter" // c_null_char
    ldum = ImGuiTextFilter_Draw(cfilter,c_loc(str),0._c_float)
    call iw_tooltip("Filter systems by name in the list below. Use comma-separated fields&
       &and - for excluding. Example: inc1,inc2,-exc includes all systems&
       &with inc1 or inc2 and excludes systems with exc.",ttshown)
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
    call igSameLine(0._c_float,-1._c_float)

    ! button: collapse
    if (iw_button("Collapse")) then
       do i = 1, nsys
          call collapse_system(i)
       end do
    end if
    call iw_tooltip("Collapse all systems in the tree",ttshown)
    call igSameLine(0._c_float,-1._c_float)

    ! insert spacing for red buttons on the right
    call igGetContentRegionAvail(szavail)
    sz%x = g%Style%ItemSpacing%x
    str = "Close" // c_null_char
    call igCalcTextSize(sztext,c_loc(str),c_null_ptr,.false._c_bool,-1._c_float)
    sz%x = sz%x + sztext%x + 2 * g%Style%FramePadding%x
    str = "Close All" // c_null_char
    call igCalcTextSize(sztext,c_loc(str),c_null_ptr,.false._c_bool,-1._c_float)
    sz%x = sz%x + sztext%x + 2 * g%Style%FramePadding%x
    rshift = szavail%x - sz%x
    if (rshift > 0) &
       call igSetCursorPosX(igGetCursorPosX() + rshift)

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

    ! set up the table
    str = "Structures##0,0" // c_null_char
    flags = ImGuiTableFlags_Borders
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_Hideable)
    flags = ior(flags,ImGuiTableFlags_Sortable)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    if (igBeginTable(c_loc(str),14,flags,szero,0._c_float)) then
       ! force resize if asked for
       if (w%forceresize) then
          call igTableSetColumnWidthAutoAll(igGetCurrentTable())
          w%forceresize = .false.
       end if

       ! set up the columns
       ! closebutton - ID - name - spg - volume - nneq - ncel - nmol - a - b - c - alpha - beta - gamma
       str = "x##0closebutton" // c_null_char
       flags = ImGuiTableColumnFlags_NoResize
       flags = ior(flags,ImGuiTableColumnFlags_NoReorder)
       flags = ior(flags,ImGuiTableColumnFlags_NoHide)
       flags = ior(flags,ImGuiTableColumnFlags_NoSort)
       flags = ior(flags,ImGuiTableColumnFlags_NoHeaderLabel)
       flags = ior(flags,ImGuiTableColumnFlags_NoHeaderWidth)
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_closebutton)

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

          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float);
          hadenabledcolumn = .false.

          ! close button
          buttonhovered_close = .false.
          if (igTableSetColumnIndex(ic_closebutton)) then
             str = "##1closebutton" // string(ic_closebutton) // "," // string(i) // c_null_char
             if (my_CloseButton(c_loc(str),ColorDangerButton)) w%forceremove = (/i/)
             buttonhovered_close = igIsItemHovered(ImGuiHoveredFlags_None)
          end if

          ! set background color for the name cell, if not selected
          if (w%table_selected /= i) then
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
          end if

          ! name
          buttonhovered_expand = .false.
          if (igTableSetColumnIndex(ic_name)) then
             dopop = .false.
             if (sysc(i)%collapse < 0) then
                ! push to reduce the spacing between button and name
                dopop = .true.
                sz%x = 1._c_float
                sz%y = 0._c_float
                call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz)

                ! expand button for multi-seed entries
                str = "##expand" // string(ic_name) // "," // string(i) // c_null_char
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
                call igSameLine(0._c_float,-1._c_float)
             end if

             ! the actual name
             str = ""
             if (sysc(i)%collapse > 0) then
                str = "├[" // string(sysc(i)%collapse) // "]─"
             end if
             str = str // trim(sysc(i)%seed%name)
             call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             if (dopop) &
                call igPopStyleVar(1_c_int)
          end if

          ! ID column
          if (igTableSetColumnIndex(ic_id)) then
             str = string(i)
             call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
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
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if

             if (igTableSetColumnIndex(ic_v)) then ! volume
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%omega*bohrtoa**3,'f',decimal=2)
                end if
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if

             if (igTableSetColumnIndex(ic_nneq)) then ! nneq
                str = string(sys(i)%c%nneq)
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if

             if (igTableSetColumnIndex(ic_ncel)) then ! ncel
                str = string(sys(i)%c%ncel)
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if

             if (igTableSetColumnIndex(ic_nmol)) then ! nmol
                str = string(sys(i)%c%nmol)
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if

             if (igTableSetColumnIndex(ic_a)) then ! a
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%aa(1)*bohrtoa,'f',decimal=4)
                end if
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if
             if (igTableSetColumnIndex(ic_b)) then ! b
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%aa(2)*bohrtoa,'f',decimal=4)
                end if
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if
             if (igTableSetColumnIndex(ic_c)) then ! c
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%aa(3)*bohrtoa,'f',decimal=4)
                end if
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if
             if (igTableSetColumnIndex(ic_alpha)) then ! alpha
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%bb(1),'f',decimal=2)
                end if
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if
             if (igTableSetColumnIndex(ic_beta)) then ! beta
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%bb(2),'f',decimal=2)
                end if
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if
             if (igTableSetColumnIndex(ic_gamma)) then ! gamma
                if (sys(i)%c%ismolecule) then
                   str = "<mol>"
                else
                   str = string(sys(i)%c%bb(3),'f',decimal=2)
                end if
                call write_text_maybe_selectable(i,str,buttonhovered_close,buttonhovered_expand)
             end if

          end if
       end do

       ! process the keybindings
       if (igIsWindowFocused(ImGuiFocusedFlags_None)) then
          if (is_bind_event(BIND_TREE_REMOVE_SYSTEM)) &
             w%forceremove = (/w%table_selected/)
       end if

       call igEndTable()
    end if

    ! clean up
    ! call ImGuiTextFilter_destroy(cfilter)

  contains

    subroutine write_text_maybe_selectable(isys,str,bclose,bexpand)
      use gui_main, only: tree_select_updates_inpcon
      use global, only: iunit, iunit_bohr, iunit_ang
      use tools_io, only: uout
      integer, intent(in) :: isys
      character(kind=c_char,len=:), allocatable, target :: str
      logical, intent(in) :: bclose, bexpand

      integer(c_int) :: flags, ll
      logical(c_bool) :: selected, enabled
      logical, save :: ttshown = .false. ! delayed tooltips
      character(kind=c_char,len=:), allocatable, target :: strl, strpop, strpop2
      character(kind=c_char,len=1024), target :: txtinp

      if (.not.hadenabledcolumn) then
         ! selectable that spans all columns
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
         end if
         call igSameLine(0._c_float,-1._c_float)

         ! right click to open the context menu
         if (igBeginPopupContextItem(c_loc(strl),ImGuiPopupFlags_MouseButtonRight)) then
            ! describe this system in the console output
            strpop = "Describe (Output Console)" // c_null_char
            enabled = (sysc(isys)%status == sys_init)
            if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,enabled)) then
               write (uout,'(/"### Describe system (",A,"): ",A/)') string(isys),&
                  trim(sysc(isys)%seed%name)
               if (sys(isys)%c%ismolecule) then
                  iunit = iunit_ang
               else
                  iunit = iunit_bohr
               end if
               call sys(isys)%report(.true.,.true.,.true.,.true.,.true.,.true.,.true.)
               iunit = iunit_bohr
            end if

            ! set as current system option
            strpop = "Set as Current System" // c_null_char
            enabled = (sysc(isys)%status == sys_init)
            if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,enabled)) &
               win(iwin_console_input)%inpcon_selected = isys

            ! remove option
            strpop = "Remove" // c_null_char
            if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,.true._c_bool)) &
               w%forceremove = (/isys/)

            ! rename option
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

            call igEndPopup()
         end if

         ! delayed tooltip with info about the system
         if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
            if (bclose) then
               strl = "Close this system" // c_null_char
            elseif (bexpand) then
               strl = "Expand this system" // c_null_char
            else
               strl = tree_tooltip_string(isys)
            end if
            call igSetTooltip(c_loc(strl))
         end if
      end if
      str = str // c_null_char
      if (sysc(isys)%status == sys_init) then
         call igText(c_loc(str))
      else
         call igTextDisabled(c_loc(str))
      end if
      hadenabledcolumn = .true.

    end subroutine write_text_maybe_selectable

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
       cid == ic_alpha .or. cid == ic_beta .or. cid == ic_gamma) then
       ! sort by real
       allocate(rval(n))
       do i = 1, n
          doit = sysc(w%iord(i))%status == sys_init
          if (doit) doit = (.not.sys(w%iord(i))%c%ismolecule)
          if (doit) then
             if (cid == ic_v) then
                rval(i) = sys(w%iord(i))%c%omega
             elseif (cid == ic_a) then
                rval(i) = sys(w%iord(i))%c%aa(1)
             elseif (cid == ic_b) then
                rval(i) = sys(w%iord(i))%c%aa(2)
             elseif (cid == ic_c) then
                rval(i) = sys(w%iord(i))%c%aa(3)
             elseif (cid == ic_alpha) then
                rval(i) = sys(w%iord(i))%c%bb(1)
             elseif (cid == ic_beta) then
                rval(i) = sys(w%iord(i))%c%bb(2)
             elseif (cid == ic_gamma) then
                rval(i) = sys(w%iord(i))%c%bb(3)
             end if
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
    if (IGFD_DisplayDialog(w%ptr,c_loc(w%name),ImGuiWindowFlags_None,minsize,maxsize)) then
       ! the dialog has been closed
       if (IGFD_IsOk(w%ptr)) then
          ! with an OK, gather information
          if (w%dialog_purpose == wpurp_dialog_openfiles) then
             !! open files dialog !!
             ! open all files selected and add the new systems
             sel = IGFD_GetSelection(w%ptr)
             call c_f_pointer(sel%table,s,(/sel%count/))
             cstr = IGFD_GetCurrentPath(w%ptr)
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
             cstr = IGFD_GetFilePathName(w%ptr)
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
          elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile) then
             !! new structure file dialog !!
             cstr = IGFD_GetFilePathName(w%ptr)
             call C_F_string_alloc(cstr,name)
             call c_free(cstr)
             w%libraryfile = trim(name) // c_null_char
             w%libraryfile_set = .true.
          else
             call ferror('draw_dialog','unknown dialog purpose',faterr)
          end if
       end if

       ! close the dialog and terminate the window
       call IGFD_CloseDialog(w%ptr)
       call w%end()
    end if

  end subroutine draw_dialog

  !> Draw the contents of the input console
  module subroutine draw_ci(w)
    use gui_keybindings, only: BIND_INPCON_RUN, get_bind_keyname, is_bind_event
    use gui_main, only: ColorHighlightText, sys, sysc, nsys, sys_init, g,&
       force_run_commands
    use gui_templates, only: draw_keyword_context_menu
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button
    use systemmod, only: sy
    use tools_io, only: string
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: sz, szero, szavail
    logical(c_bool) :: ldum, is_selected
    logical, save :: ttshown = .false.
    logical :: ok
    integer :: i, idx

    ! initialize
    szero%x = 0
    szero%y = 0

    ! allocate the input buffer if not already done
    if (.not.allocated(inputb)) then
       allocate(character(len=maxlib+1) :: inputb)
       inputb(1:1) = c_null_char
    end if

    ! first line: text
    str1 = "Input" // c_null_char
    call igTextColored(ColorHighlightText,c_loc(str1))

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
    str1 = "System" // c_null_char
    call igText(c_loc(str1))
    call igSameLine(0._c_float,-1._c_float)

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
    call igGetContentRegionAvail(szavail)
    call igSetNextItemWidth(szavail%x - sz%x - g%Style%ItemSpacing%x)
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
    str1 = "Field " // c_null_char
    call igText(c_loc(str1))
    call igSameLine(0._c_float,-1._c_float)

    str1 = "##fieldcombo" // c_null_char
    str2 = "" // c_null_char
    if (associated(sy)) &
       str2 = string(sy%iref) // ": " // trim(sy%f(sy%iref)%name) // c_null_char
    call igGetContentRegionAvail(szavail)
    call igSetNextItemWidth(szavail%x - sz%x - g%Style%ItemSpacing%x)
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
    use gui_main, only: mainvwp, io, g, ColorWaitBg, ColorHighlightText
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
    text = "Running critic2 input..." // newline //&
       "System: " // csystem // newline //&
       "Field: " // cfield // newline //&
       "Input: " // newline // inputb
    call igCalcTextSize(sz,c_loc(text),c_null_ptr,.false._c_bool,-1._c_float)
    sz%y = sz%y + 2 * g%Style%WindowPadding%y
    sz%x = sz%x + 2 * g%Style%WindowPadding%x + 3 * igGetTextLineHeight()
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
       str1 = "...Running critic2 input..." // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str1))

       str1 = "System: " // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str1))
       call igSameLine(0._c_float,-1._c_float)
       str1 = csystem // c_null_char
       call igText(c_loc(str1))

       str1 = "Field: " // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str1))
       call igSameLine(0._c_float,-1._c_float)
       str1 = cfield // c_null_char
       call igText(c_loc(str1))

       str1 = "Input: " // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str1))
       call igIndent(0._c_float)
       call igText(c_loc(inputb))

       call igPopStyleVar(1_c_int)
    end if
    call igEnd()

  end subroutine block_gui_ci

  !> Read new output from the scratch LU uout. If iscom, this output corresponds
  !> to a command, so create the command i/o object. Return true if
  !> output has been read.
  module function read_output_ci(w,iscom)
    use gui_utils, only: get_time_string
    use gui_main, only: are_threads_running
    use tools_io, only: uout, getline_raw, string, ferror, faterr
    use types, only: realloc
    use param, only: newline
    class(window), intent(inout), target :: w
    logical, intent(in) :: iscom
    logical :: read_output_ci

    character(kind=c_char,len=:), allocatable, target :: csystem, cfield
    type(command_inout), allocatable :: aux(:)
    character(len=:), allocatable :: line, commonstr
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
          idx = index(inputb,c_null_char)
          com(ithis)%id = ncomid
          com(ithis)%input = inputb(1:idx)
          com(ithis)%status = command_inout_used
          com(ithis)%size = lob - olob + 1

          commonstr = "### Command: " // string(ncomid) // " (" // get_time_string() // ")" // newline
          com(ithis)%tooltipinfo = commonstr // &
             "System: " // csystem // newline //&
             "Field: " // cfield // newline //&
             "Input: " // newline // inputb(1:idx-1) // newline //&
             "#########" // newline // newline // "[Right-click for options]"
          com(ithis)%output = commonstr // &
             "## System: " // csystem // newline //&
             "## Field: " // cfield // newline //&
             "## Input: " // newline // inputb(1:idx-1) // newline //&
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
    use gui_main, only: ColorHighlightText, g, ColorDangerButton,&
       ColorFrameBgAlt
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: i, curline, ndrawn, idx
    character(kind=c_char,len=:), allocatable, target :: str1, strpop
    type(ImVec2) :: sz, szero, sztext, szavail
    logical(c_bool) :: ldum
    logical, save :: ttshown = .false.
    logical :: setscroll, skip, pushed, ok
    integer, save :: idsavedialog = 0
    real(c_float) :: itemspacing, xavail, xavail1, rshift
    real(c_float), save :: allscrolly
    integer :: navail, navail1

    ! check if the save dialog is still open
    if (idsavedialog > 0) then
       if (idsavedialog < 1 .or. idsavedialog > nwin) then
          idsavedialog = 0
       elseif (.not.win(idsavedialog)%isinit .or. .not.win(idsavedialog)%isopen) then
          idsavedialog = 0
       end if
    end if

    ! read new output, if available
    ldum = w%read_output_ci(.false.)

    ! initialize
    szero%x = 0._c_float
    szero%y = 0._c_float
    setscroll = .false.

    ! get the available x
    call igGetContentRegionAvail(szavail)
    xavail = szavail%x

    ! first line: text
    str1 = "Output" // c_null_char
    call igTextColored(ColorHighlightText,c_loc(str1))
    call igSameLine(0._c_float,-1._c_float)

    ! first line: clear button
    if (idcom == 0) then
       if (iw_button("Clear")) then
          outputb(1:1) = c_null_char
          lob = 0
       end if
       call iw_tooltip("Clear the output log",ttshown)
       call igSameLine(0._c_float,-1._c_float)
    end if

    ! first line: copy button
    if (iw_button("Copy")) then
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
    call igSameLine(0._c_float,-1._c_float)
    call igGetContentRegionAvail(szavail)
    str1 = "Remove All" // c_null_char
    call igCalcTextSize(sztext,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
    rshift = szavail%x - (sztext%x + 2 * g%Style%FramePadding%x)
    if (rshift > 0) &
       call igSetCursorPosX(igGetCursorPosX() + rshift)

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
    str1 = "All" // c_null_char
    call igCalcTextSize(sztext,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
    xavail1 = xavail - (sztext%x + 2 * g%Style%FramePadding%x + g%Style%ItemSpacing%x)
    if (idcom == 0) then
       call igPushStyleColor_Vec4(ImGuiCol_Button,g%Style%Colors(ImGuiCol_ButtonActive+1))
    else
       call igPushStyleColor_Vec4(ImGuiCol_Button,ColorDangerButton)
    end if
    if (igButton(c_loc(str1),szero)) then
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
    sz%y = igGetTextLineHeight() + 2 * g%Style%FramePadding%y
    str1 = string(maxval(com(icom(1:nicom))%id)) // c_null_char
    call igCalcTextSize(sztext,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
    sz%x = max(sz%y,sztext%x + 2 * g%Style%FramePadding%x)
    xavail = xavail / sz%x
    xavail1 = xavail1 / sz%x
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
       if (iw_button(string(com(icom(i))%id),popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonRight)) then
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
    if (idcom > 0) then
       com(icom(idcom))%scrolly = igGetScrollY()
    else
       allscrolly = igGetScrollY()
    end if
    call igEndChild()

  end subroutine draw_co

  !> Draw the contents of the new structure window.
  module subroutine draw_new(w)
    use gui_main, only: g, ColorHighlightText, add_systems_from_seeds,&
       launch_initialization_thread, system_shorten_names
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button
    use crystalseedmod, only: crystalseed, realloc_crystalseed
    use spglib, only: SpglibSpaceGroupType, spg_get_spacegroup_type
    use global, only: clib_file, mlib_file, rborder_def
    use tools_io, only: string, fopen_scratch, fclose, ioj_left, stripchar, deblank
    use types, only: vstring
    use param, only: newline, bohrtoa
    class(window), intent(inout), target :: w

    integer(c_size_t), parameter :: namebufsiz = 1024

    character(kind=c_char,len=:), allocatable, target :: str, str2, stropt, strex, left
    logical(c_bool) :: ldum, doquit
    logical :: changed, readlib, ok, exloop
    type(ImVec2) :: szero, sz, sztext, szavail
    integer :: i, nseed, idx, lu
    real(c_float) :: rshift
    type(crystalseed) :: seed
    type(crystalseed), allocatable :: seed_(:)
    integer(c_int) :: flags
    type(SpglibSpaceGroupType) :: sa

    logical(c_bool), save :: ismolecule = .false.
    logical(c_bool), save :: fromlibrary = .false.
    logical, save :: ttshown = .false.
    integer(c_int), save :: idum = 0
    logical, save :: firstpass = .true.
    integer, save :: nst = 0
    integer, save :: lastselected = 0
    real(c_float), save :: rborder = real(rborder_def*bohrtoa,c_float)
    logical(c_bool), save :: molcubic = .false.
    integer(c_int), save :: iunit = 2 ! 0 = bohr, 1 = angstrom, 2 = fractional
    integer(c_int), save :: iunitcel = 0 ! 0 = bohr, 1 = angstrom
    character(len=namebufsiz,kind=c_char), target, save :: namebuf = c_null_char
    type(vstring), allocatable, save :: st(:)
    logical(c_bool), allocatable, save :: lst(:)
    character(kind=c_char,len=:), allocatable, target, save :: atposbuf
    character(kind=c_char,len=:), allocatable, target, save :: symopbuf
    character(kind=c_char,len=:), allocatable, target, save :: latvecbuf
    integer, save :: symopt = 1 ! 1 = detect, 2 = spg, 3 = manual
    integer, save :: cellopt = 1 ! 1 = parameters, 2 = lattice-vectors
    integer, save :: ispg_selected = 1
    real(c_float), save :: aa(3) = 0, bb(3) = 0, scale = 1

    integer(c_size_t), parameter :: maxatposbuf = 100000
    integer(c_size_t), parameter :: maxsymopbuf = 10000
    integer(c_size_t), parameter :: maxlatvecbuf = 5000

    ! initialize
    szero%x = 0
    szero%y = 0
    readlib = .false.
    if (firstpass) then
       w%libraryfile = trim(clib_file) // c_null_char
       firstpass = .false.
       readlib = .true.
       if (allocated(atposbuf)) deallocate(atposbuf)
       allocate(character(len=maxatposbuf+1) :: atposbuf)
       atposbuf(1:26) = "H 0.00000 0.00000 0.00000" // c_null_char
       if (allocated(symopbuf)) deallocate(symopbuf)
       allocate(character(len=maxsymopbuf+1) :: symopbuf)
       symopbuf(1:66) = "## Each line is a symmetry operation of the type: " // newline //&
          "## -1/2+y,z,x" // c_null_char
       if (allocated(latvecbuf)) deallocate(latvecbuf)
       allocate(character(len=maxlatvecbuf+1) :: latvecbuf)
       latvecbuf(1:50) = "## Enter the three lattice vectors line by line" // newline // c_null_char
       aa = 10._c_float
       bb = 90._c_float
       scale = 1._c_float
    end if

    ! check if we have info from the open library file window when it
    ! closes
    if (idum > 0) then
       if (.not.win(idum)%isopen) then
          if (win(idum)%libraryfile_set) then
             w%libraryfile_set = .true.
             w%libraryfile = win(idum)%libraryfile
             call win(idum)%end()
             readlib = .true.
          end if
          idum = 0
       end if
    end if

    ! crystal or molecule, from library
    changed = .false.
    str = "Crystal" // c_null_char
    if (igRadioButton_Bool(c_loc(str),.not.ismolecule)) then
       ismolecule = .false.
       changed = .true.
       iunit = 2
       namebuf = "Crystal structure (manual)" // c_null_char
    end if
    call iw_tooltip("The new structure will be a periodic crystal",ttshown)
    call igSameLine(0._c_float,-1._c_float)
    str = "Molecule" // c_null_char
    if (igRadioButton_Bool(c_loc(str),ismolecule)) then
       ismolecule = .true.
       changed = .true.
       iunit = 1
       namebuf = "Molecule structure (manual)" // c_null_char
    end if
    call iw_tooltip("The new structure will be a molecule",ttshown)
    call igSameLine(0._c_float,-1._c_float)
    call igSetCursorPosX(igGetCursorPosX() + 4 * g%Style%ItemSpacing%x)
    str = "From library" // c_null_char
    ldum = igCheckbox(c_loc(str),fromlibrary)
    call iw_tooltip("Read the structure from the critic2 library",ttshown)
    if (changed.and..not.w%libraryfile_set) then
       if (ismolecule) then
          w%libraryfile = trim(mlib_file) // c_null_char
       else
          w%libraryfile = trim(clib_file) // c_null_char
       end if
       readlib = .true.
    end if

    ! render the rest of the window
    doquit = .false.
    if (fromlibrary) then
       ! library file
       str = "Source" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       if (iw_button("Library file",disabled=(idum /= 0))) &
          idum = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openlibraryfile)
       call iw_tooltip("Library file from where the structures are read",ttshown)
       call igSameLine(0._c_float,-1._c_float)
       call igText(c_loc(w%libraryfile))

       ! list box
       str = "Structures to load from the library file" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       str = "##listbox" // c_null_char
       call igGetContentRegionAvail(sz)
       if (ismolecule) then
          sz%y = sz%y - (3 * igGetTextLineHeight() + 6 * g%Style%FramePadding%y + &
             2 * g%Style%ItemSpacing%y + g%Style%WindowPadding%y)
       else
          sz%y = sz%y - (2 * igGetTextLineHeight() + g%Style%FramePadding%y + g%Style%ItemSpacing%y + g%Style%WindowPadding%y)
       end if
       ldum = igBeginListBox(c_loc(str),sz)
       do i = 1, nst
          str = st(i)%s // c_null_char
          if (igSelectable_Bool(c_loc(str), lst(i), ImGuiSelectableFlags_None, szero)) then

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
             end if
          end if
       end do
       call igEndListBox()

       if (ismolecule) then
          ! options line
          str = "Structure options" // c_null_char
          call igTextColored(ColorHighlightText,c_loc(str))

          ! cell border
          call igIndent(0._c_float)
          str = "Cell border (Å)" // c_null_char
          stropt = "%.3f" // c_null_char
          strex = string(rborder,'f',decimal=3) // c_null_char
          call igCalcTextSize(sz,c_loc(strex),c_null_ptr,.false._c_bool,-1._c_float)
          call igPushItemWidth(sz%x + 2 * g%Style%FramePadding%x)
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

       ! insert spacing for buttons on the right
       call igGetContentRegionAvail(szavail)
       sz%x = g%Style%ItemSpacing%x
       str = "OK" // c_null_char
       call igCalcTextSize(sztext,c_loc(str),c_null_ptr,.false._c_bool,-1._c_float)
       sz%x = sz%x + sztext%x + 2 * g%Style%FramePadding%x
       str = "Cancel" // c_null_char
       call igCalcTextSize(sztext,c_loc(str),c_null_ptr,.false._c_bool,-1._c_float)
       sz%x = sz%x + sztext%x + 2 * g%Style%FramePadding%x
       rshift = szavail%x - sz%x
       if (rshift > 0) &
          call igSetCursorPosX(igGetCursorPosX() + rshift)

       ! final buttons: ok
       if (iw_button("OK",disabled=(idum /= 0))) then
          nseed = count(lst(1:nst))
          if (nseed > 0) then
             ! we have systems (potentially), allocate the seed array
             allocate(seed_(nseed))
             nseed = 0
             do i = 1, nst
                if (lst(i)) then
                   ! read the selected seeds from the library, check OK
                   call seed%read_library(st(i)%s,logical(ismolecule),ok)
                   if (ok) then
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

    elseif (ismolecule) then
       ! name
       str = "Name" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       str2 = "##name"
       ldum = igInputText(c_loc(str2),c_loc(namebuf),namebufsiz-1,ImGuiInputTextFlags_None,c_null_ptr,c_null_ptr)

       ! atomic positions: header
       str = "Atomic positions" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       call igSameLine(0._c_float,-1._c_float)
       str = "(?)" // c_null_char
       call igText(c_loc(str))
       call iw_tooltip("Give the atomic positions for this system as:" // newline //&
          "  <Sy> <x> <y> <z>" // newline //&
          "where Sy is the atomic symbol and x,y,z are the atomic coordinates.")

       ! units
       call igSameLine(0._c_float,-1._c_float)
       call igSetCursorPosX(igGetCursorPosX() + 2 * g%Style%ItemSpacing%x)
       str = "Units" // c_null_char
       stropt = "Bohr" // c_null_char // "Angstrom"
       strex = "Angstrom    " // c_null_char
       call igCalcTextSize(sz,c_loc(strex),c_null_ptr,.false._c_bool,-1._c_float)
       call igSetNextItemWidth(sz%x)
       ldum = igCombo_Str(c_loc(str), iunit, c_loc(stropt), -1_c_int)
       call iw_tooltip("Units for the atomic coordinates",ttshown)

       ! atomic positions: body
       call igGetContentRegionAvail(szavail)
       sz%x = szavail%x
       sz%y = szavail%y - (3 * igGetTextLineHeight() + 6 * g%Style%FramePadding%y + &
          2 * g%Style%ItemSpacing%y + g%Style%WindowPadding%y)
       str = "##atomicpositions" // c_null_char
       ldum = igInputTextMultiline(c_loc(str),c_loc(atposbuf),maxatposbuf,sz,&
          ImGuiInputTextFlags_AllowTabInput,c_null_ptr,c_null_ptr)

       ! options line
       str = "Structure options" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))

       ! cell border
       call igIndent(0._c_float)
       str = "Cell border (Å)" // c_null_char
       stropt = "%.3f" // c_null_char
       strex = string(rborder,'f',decimal=3) // c_null_char
       call igCalcTextSize(sz,c_loc(strex),c_null_ptr,.false._c_bool,-1._c_float)
       call igPushItemWidth(sz%x + 2 * g%Style%FramePadding%x)
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
       call igUnindent(0._c_float)

       ! insert spacing for buttons on the right
       call igGetContentRegionAvail(szavail)
       sz%x = g%Style%ItemSpacing%x
       str = "OK" // c_null_char
       call igCalcTextSize(sztext,c_loc(str),c_null_ptr,.false._c_bool,-1._c_float)
       sz%x = sz%x + sztext%x + 2 * g%Style%FramePadding%x
       str = "Cancel" // c_null_char
       call igCalcTextSize(sztext,c_loc(str),c_null_ptr,.false._c_bool,-1._c_float)
       sz%x = sz%x + sztext%x + 2 * g%Style%FramePadding%x
       rshift = szavail%x - sz%x
       if (rshift > 0) &
          call igSetCursorPosX(igGetCursorPosX() + rshift)

       ! final buttons: ok
       if (iw_button("OK")) then
          ! build the input
          lu = fopen_scratch("formatted")
          left = atposbuf(1:index(atposbuf,c_null_char)-1)
          idx = 0
          exloop = .false.
          do while (.not.exloop)
             idx = index(left,newline)
             if (idx == 0) then
                idx = len_trim(left) + 1
                exloop = .true.
             end if
             if (len_trim(left(:idx-1)) > 0) then
                if (iunit == 0) then
                   write (lu,'(A," bohr")') left(:idx-1)
                else
                   write (lu,'(A," ang")') left(:idx-1)
                end if
             end if
             if(.not.exloop) left = left(idx+1:)
          end do
          if (molcubic) write (lu,'("cubic")')
          write (lu,'("border ",A)') string(rborder/bohrtoa,'f',decimal=10)
          write (lu,'("end")')

          ! generate the seed
          rewind(lu)
          if (allocated(seed_)) deallocate(seed_)
          allocate(seed_(1))
          call seed_(1)%parse_molecule_env(lu,ok)
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
       call igSameLine(0._c_float,-1._c_float)

       ! final buttons: cancel
       if (iw_button("Cancel")) doquit = .true.
    else
       ! name
       str = "Name" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       str2 = "##name"
       ldum = igInputText(c_loc(str2),c_loc(namebuf),namebufsiz-1,ImGuiInputTextFlags_None,c_null_ptr,c_null_ptr)

       ! symmetry
       str = "Symmetry" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))

       str = "Detect" // c_null_char
       if (igRadioButton_Bool(c_loc(str),logical(symopt == 1,c_bool))) symopt = 1
       call iw_tooltip("Calculate the symmetry operations from the complete list of unit cell atoms",ttshown)
       call igSameLine(0._c_float,-1._c_float)
       str = "Space group" // c_null_char
       if (igRadioButton_Bool(c_loc(str),logical(symopt == 2,c_bool))) symopt = 2
       call iw_tooltip("Choose a space group from the list",ttshown)
       call igSameLine(0._c_float,-1._c_float)
       str = "Manual" // c_null_char
       if (igRadioButton_Bool(c_loc(str),logical(symopt == 3,c_bool))) symopt = 3
       call iw_tooltip("Enter symmetry operations in cif-file format (e.g. -1/2+x,z,-y)",ttshown)

       ! symmetry options
       if (symopt == 2) then
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
                   if (igSelectable_Bool(c_loc(str),logical(ispg_selected == i,c_bool),flags,szero)) &
                      ispg_selected = i
                end if

                ! ITA
                if (igTableNextColumn()) then
                   str = string(sa%number,5,ioj_left) // c_null_char
                   call igText(c_loc(str))
                end if

                ! HM-short
                if (igTableNextColumn()) then
                   str = deblank(sa%international_short)
                   str = trim(stripchar(str,"_")) // c_null_char
                   call igText(c_loc(str))
                end if

                ! HM-long
                if (igTableNextColumn()) then
                   str = deblank(sa%international_full)
                   str = trim(stripchar(str,"_")) // c_null_char
                   call igText(c_loc(str))
                end if

                ! choice
                if (igTableNextColumn()) then
                   str = string(sa%choice,6,ioj_left) // c_null_char
                   call igText(c_loc(str))
                end if

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
                   str = str // c_null_char
                   call igText(c_loc(str))
                end if

                ! hall symbol
                if (igTableNextColumn()) then
                   str = string(sa%hall_symbol) // c_null_char
                   call igText(c_loc(str))
                end if
             end do

             call igTableSetColumnWidthAutoAll(igGetCurrentTable())
             call igEndTable()
          end if

          ! current choice
          sa = spg_get_spacegroup_type(ispg_selected)
          str = deblank(sa%international_full)
          str = "Current choice: " // trim(stripchar(str,"_"))
          if (len_trim(sa%choice) > 0) str = str // " (" // trim(sa%choice) // ")"
          str = str // c_null_char
          call igText(c_loc(str))
       elseif (symopt == 3) then
          ! manual input of symmetry operations
          call igGetContentRegionAvail(szavail)
          sz%x = szavail%x
          sz%y = 8 * igGetTextLineHeightWithSpacing()
          str = "##symopsmanual" // c_null_char
          ldum = igInputTextMultiline(c_loc(str),c_loc(symopbuf),maxsymopbuf,sz,&
             ImGuiInputTextFlags_AllowTabInput,c_null_ptr,c_null_ptr)
       end if

       ! lattice: header
       str = "Lattice" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))

       str = "Cell parameters" // c_null_char
       if (igRadioButton_Bool(c_loc(str),logical(cellopt == 1,c_bool))) then
          cellopt = 1
          iunit = 2
       end if
       call iw_tooltip("Lattice is calculated from the cell lengths and angles",ttshown)
       call igSameLine(0._c_float,-1._c_float)

       str = "Lattice vectors" // c_null_char
       if (igRadioButton_Bool(c_loc(str),logical(cellopt == 2,c_bool))) cellopt = 2
       call iw_tooltip("Lattice vectors are given explicitly",ttshown)
       call igSameLine(0._c_float,-1._c_float)

       str = "Units##cel" // c_null_char
       stropt = "Bohr" // c_null_char // "Angstrom" // c_null_char // c_null_char
       strex = "Angstrom    " // c_null_char
       call igCalcTextSize(sz,c_loc(strex),c_null_ptr,.false._c_bool,-1._c_float)
       call igSetNextItemWidth(sz%x)
       ldum = igCombo_Str(c_loc(str), iunitcel, c_loc(stropt), -1_c_int)
       call iw_tooltip("Units for the cell parameters/lattice vectors",ttshown)

       if (cellopt == 2) then
          call igSameLine(0._c_float,-1._c_float)
          str = "Scale" // c_null_char
          stropt = "%.3f" // c_null_char
          strex = string(scale,'f',decimal=3) // c_null_char
          call igCalcTextSize(sz,c_loc(strex),c_null_ptr,.false._c_bool,-1._c_float)
          call igPushItemWidth(sz%x + 2 * g%Style%FramePadding%x)
          ldum = igInputFloat(c_loc(str),scale,0._c_float,0._c_float,&
             c_loc(stropt),ImGuiInputTextFlags_None)
          call iw_tooltip("Scale factor to multiply the lattice vectors",ttshown)
          call igPopItemWidth()
       end if

       ! lattice: body
       if (cellopt == 1) then
          ! cell lengths
          str = "Cell lengths (Å): " // c_null_char
          call igText(c_loc(str))
          call igSameLine(0._c_float,-1._c_float)
          str = "##celllength3" // c_null_char
          stropt = "%.3f" // c_null_char
          ldum = igInputFloat3(c_loc(str),aa,c_loc(stropt),ImGuiInputTextFlags_None)

          ! cell angles
          str = "Cell angles (°):  " // c_null_char
          call igText(c_loc(str))
          call igSameLine(0._c_float,-1._c_float)
          str = "##cellangle3" // c_null_char
          stropt = "%.3f" // c_null_char
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

       ! atomic positions: header
       str = "Atomic positions" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       call igSameLine(0._c_float,-1._c_float)
       str = "(?)" // c_null_char
       call igText(c_loc(str))
       if (symopt == 1) then
          str = "Enter all atoms in the unit cell as:"
       else
          str = "Enter the non-equivalent (symmetry-unique) atoms in the unit cell as:"
       end if
       str = str // newline // "  <Sy> <x> <y> <z>" // newline //&
          "where Sy is the atomic symbol and x,y,z are the "
       if (cellopt == 1) then
          str = str // "fractional atomic coordinates."
       else
          str = str // "atomic coordinates."
       end if
       call iw_tooltip(str)

       if (cellopt == 2) then
          call igSameLine(0._c_float,-1._c_float)
          call igSetCursorPosX(igGetCursorPosX() + 2 * g%Style%ItemSpacing%x)
          str = "Units" // c_null_char
          stropt = "Bohr" // c_null_char // "Angstrom" // c_null_char // "Fractional" // c_null_char // c_null_char
          strex = "Fractional    " // c_null_char
          call igCalcTextSize(sz,c_loc(strex),c_null_ptr,.false._c_bool,-1._c_float)
          call igSetNextItemWidth(sz%x)
          ldum = igCombo_Str(c_loc(str), iunit, c_loc(stropt), -1_c_int)
          call iw_tooltip("Units for the atomic coordinates",ttshown)
       end if

       ! atomic positions: body
       call igGetContentRegionAvail(szavail)
       sz%x = szavail%x
       sz%y = szavail%y - (2 * igGetTextLineHeight() + 4 * g%Style%FramePadding%y + &
          2 * g%Style%ItemSpacing%y + g%Style%WindowPadding%y)
       str = "##atomicpositions" // c_null_char
       ldum = igInputTextMultiline(c_loc(str),c_loc(atposbuf),maxatposbuf,sz,&
          ImGuiInputTextFlags_AllowTabInput,c_null_ptr,c_null_ptr)

       ! insert spacing for buttons on the right
       call igGetContentRegionAvail(szavail)
       sz%x = g%Style%ItemSpacing%x
       str = "OK" // c_null_char
       call igCalcTextSize(sztext,c_loc(str),c_null_ptr,.false._c_bool,-1._c_float)
       sz%x = sz%x + sztext%x + 2 * g%Style%FramePadding%x
       str = "Cancel" // c_null_char
       call igCalcTextSize(sztext,c_loc(str),c_null_ptr,.false._c_bool,-1._c_float)
       sz%x = sz%x + sztext%x + 2 * g%Style%FramePadding%x
       rshift = szavail%x - sz%x
       if (rshift > 0) &
          call igSetCursorPosX(igGetCursorPosX() + rshift)

       ! final buttons: ok
       if (iw_button("OK")) then
          ! build the input
          lu = fopen_scratch("formatted")

          ! symmetry
          if (symopt == 2) then
             ! pass the spg string
             write (lu,'("spg ",A)') string(ispg_selected)
          elseif (symopt == 3) then
             ! pass the symm keywords
             idx = 0
             left = symopbuf(1:index(symopbuf,c_null_char)-1)
             exloop = .false.
             do while (.not.exloop)
                idx = index(left,newline)
                if (idx == 0) then
                   idx = len_trim(left) + 1
                   exloop = .true.
                end if
                if (len_trim(left(:idx-1)) > 0) &
                   write (lu,'("symm ",A)') left(:idx-1)
                if(.not.exloop) left = left(idx+1:)
             end do
          end if

          ! cell
          if (cellopt == 1) then
             ! cell parameters
             if (iunitcel == 1) then
                write (lu,'("cell ",6(A," "),"ang")') (string(aa(i),'f',decimal=10),i=1,3), &
                   (string(bb(i),'f',decimal=10),i=1,3)
             else
                write (lu,'("cell ",6(A," "),"bohr")') (string(aa(i),'f',decimal=10),i=1,3), &
                   (string(bb(i),'f',decimal=10),i=1,3)
             end if
          else
             ! lattice vectors
             write (lu,'("cartesian ",A)') string(scale,'f',decimal=10)
             if (iunitcel == 1) then
                write (lu,'("ang")')
             else
                write (lu,'("bohr")')
             end if
             idx = 0
             left = latvecbuf(1:index(latvecbuf,c_null_char)-1)
             exloop = .false.
             do while (.not.exloop)
                idx = index(left,newline)
                if (idx == 0) then
                   idx = len_trim(left) + 1
                   exloop = .true.
                end if
                if (len_trim(left(:idx-1)) > 0) &
                   write (lu,'(A)') left(:idx-1)
                if(.not.exloop) left = left(idx+1:)
             end do
             write (lu,'("endcartesian")')
          end if

          ! atoms
          left = atposbuf(1:index(atposbuf,c_null_char)-1)
          idx = 0
          exloop = .false.
          do while (.not.exloop)
             idx = index(left,newline)
             if (idx == 0) then
                idx = len_trim(left) + 1
                exloop = .true.
             end if
             if (len_trim(left(:idx-1)) > 0) then
                if (iunit == 0) then
                   write (lu,'(A," bohr")') left(:idx-1)
                elseif (iunit == 1) then
                   write (lu,'(A," ang")') left(:idx-1)
                else
                   write (lu,'(A)') left(:idx-1)
                end if
             end if
             if(.not.exloop) left = left(idx+1:)
          end do
          write (lu,'("end")')

          ! generate the seed
          rewind(lu)
          if (allocated(seed_)) deallocate(seed_)
          allocate(seed_(1))
          call seed_(1)%parse_crystal_env(lu,ok)
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
    end if

    ! read the library file
    if (readlib) then
       ! use doquit
       call get_library_structure_list(w%libraryfile,nst,st)
       if (allocated(lst)) deallocate(lst)
       allocate(lst(nst))
       lst = .false.
       lastselected = 0
    end if

    ! quit the window
    if (doquit) then
       ! reset the state and kill the window
       ismolecule = .false.
       fromlibrary = .false.
       ttshown = .false.
       idum = 0
       firstpass = .true.
       lastselected = 0
       nst = 0
       rborder = real(rborder_def*bohrtoa,c_float)
       molcubic = .false.
       iunit = 0
       iunitcel = 0
       namebuf(1:1) = c_null_char
       symopt = 1
       cellopt = 1
       ispg_selected = 1
       aa = 10._c_float
       bb = 90._c_float
       scale = 1._c_float
       if (allocated(st)) deallocate(st)
       if (allocated(lst)) deallocate(lst)
       if (allocated(atposbuf)) deallocate(atposbuf)
       call w%end()
    end if

  end subroutine draw_new

  !xx! private procedures

  ! Return the string for the tooltip shown by the tree window,
  ! corresponding to system i.
  function tree_tooltip_string(i) result(str)
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
          ! cell parameters, volume, density
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
          str = str // "Density (g/cm³): " // string(dens,'f',decimal=3) // newline &
             // newline

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

  end function tree_tooltip_string

  ! the callback for the right-hand-side pane of the dialog
  subroutine dialog_user_callback(vFilter, vUserData, vCantContinue) bind(c)
    use gui_main, only: ColorHighlightText, g
    use gui_utils, only: igIsItemHovered_delayed, iw_tooltip
    use gui_interfaces_cimgui
    use tools_io, only: string
    type(c_ptr), intent(in), value :: vFilter ! const char *
    type(c_ptr), value :: vUserData ! void *
    logical(c_bool) :: vCantContinue ! bool *

    character(kind=c_char,len=:), allocatable, target :: str, stropt, strex
    type(dialog_userdata), pointer :: data
    logical(c_bool) :: ldum
    logical, save :: ttshown = .false.
    type(ImVec2) :: sz

    ! generate the data pointer
    call c_f_pointer(vUserData,data)

    !! common options !!

    ! header
    str = "Options" // c_null_char
    call igTextColored(ColorHighlightText,c_loc(str))

    ! show hidden files
    str = "Show hidden files" // c_null_char
    if (igCheckbox(c_loc(str),data%showhidden)) then
       if (data%showhidden) then
          call IGFD_SetFlags(data%ptr,ImGuiFileDialogFlags_None)
       else
          call IGFD_SetFlags(data%ptr,ImGuiFileDialogFlags_DontShowHiddenFiles)
       end if
    end if
    call iw_tooltip("Show the OS hidden files and directories in this dialog",ttshown)

    !! options specific to the open files dialog !!
    if (data%purpose == wpurp_dialog_openfiles) then
       str = "Read last structure only" // c_null_char
       ldum = igCheckbox(c_loc(str),data%readlastonly)
       call iw_tooltip("Read only the last structure in the file",ttshown)
       call igNewLine()

       ! radio buttons for auto/crystal/molecule
       str = "Read structures as..." // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       str = "Auto-detect" // c_null_char
       ldum = igRadioButton_IntPtr(c_loc(str),data%mol,-1)
       call iw_tooltip("Auto-detect whether new structures are read as crystals or molecules",ttshown)
       str = "Crystal" // c_null_char
       ldum = igRadioButton_IntPtr(c_loc(str),data%mol,0)
       call iw_tooltip("Force new structures to be read as crystals",ttshown)
       str = "Molecule" // c_null_char
       ldum = igRadioButton_IntPtr(c_loc(str),data%mol,1)
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
       str = "Read file format" // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       str = "##formatcombo" // c_null_char
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
          // "xyz file" // c_null_char &                ! isformat_xyz = 12
          // c_null_char
       ldum = igCombo_Str(c_loc(str), data%isformat,c_loc(stropt),-1_c_int);
       call iw_tooltip("Force new structures read with a given file format, or auto-detect&
          &from the extension",ttshown)
       call igNewLine()
    end if

  end subroutine dialog_user_callback

  !> Get the structure list from the library file
  subroutine get_library_structure_list(libfile,nst,st)
    use tools_io, only: fopen_read, fclose, lgetword, getword, equal, getline
    use types, only: vstring, realloc
    character(len=:), allocatable, intent(in) :: libfile
    integer, intent(out) :: nst
    type(vstring), allocatable, intent(inout) :: st(:)

    integer :: idx, lu, lp
    character(len=:), allocatable :: file, word, line
    logical :: ok

    ! get the file name
    idx = index(libfile,c_null_char)
    if (idx > 0) then
       file = libfile(1:idx-1)
    else
       file = libfile
    end if

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
          word = getword(line,lp)
          nst = nst + 1
          if (nst > size(st,1)) call realloc(st,2*nst)
          st(nst)%s = word
       endif
    end do main
    if (nst > 0) &
       call realloc(st,nst)

    ! clean up
    call fclose(lu)

  end subroutine get_library_structure_list

end submodule proc
