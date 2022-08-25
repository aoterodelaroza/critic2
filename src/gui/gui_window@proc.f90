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

  !xx! private procedures
  ! function tree_tooltip_string(i)
  ! subroutine dialog_user_callback(vFilter, vUserData, vCantContinue)

contains

  !> Create a window in the window stack with the given type. Returns
  !> the window ID.
  function stack_create_window(type,isopen,purpose)
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
          else
             call ferror('window_draw','unknown dialog purpose',faterr)
          end if
       end if
    end if

    if (w%isopen) then
       if (w%type == wintype_tree .or. w%type == wintype_console_input .or.&
          w%type == wintype_console_output .or. w%type == wintype_view) then
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
                call w%draw_console_input()
             elseif (w%type == wintype_console_output) then
                call w%draw_console_output()
             end if
          end if
          call igEnd()
       elseif (w%type == wintype_dialog) then
          call w%draw_dialog()
       end if
    end if

  end subroutine window_draw

  !> Draw the contents of a tree window
  module subroutine draw_tree(w)
    use gui_keybindings, only: is_bind_event, BIND_TREE_REMOVE_SYSTEM
    use gui_utils, only: igIsItemHovered_delayed
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
    integer(c_int) :: flags, color
    integer :: i, j, k, nshown, newsel, jsel
    logical(c_bool) :: ldum
    type(c_ptr) :: ptrc
    type(ImGuiTableSortSpecs), pointer :: sortspecs
    type(ImGuiTableColumnSortSpecs), pointer :: colspecs
    logical :: hadenabledcolumn, buttonhovered_close, buttonhovered_expand, reinit
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
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str = &
          "Filter systems by name in the list below. Use comma-separated fields" // new_line('a') //&
          "and - for excluding. Example: inc1,inc2,-exc includes all systems   " // new_line('a') //&
          "with inc1 or inc2 and excludes systems with exc." // c_null_char
       call igSetTooltip(c_loc(str))
    end if
    call igSameLine(0._c_float,-1._c_float)
    str = "Clear" // c_null_char
    if (igButton(c_loc(str),szero)) then
       if (c_associated(cfilter)) &
          call ImGuiTextFilter_Clear(cfilter)
    end if
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str = "Clear the filter" // c_null_char
       call igSetTooltip(c_loc(str))
    end if

    ! row of buttons
    ! button: expand
    str = "Expand" // c_null_char
    if (igButton(c_loc(str),szero)) then
       do i = 1, nsys
          call expand_system(i)
       end do
    end if
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str = "Expand all systems in the tree" // c_null_char
       call igSetTooltip(c_loc(str))
    end if
    call igSameLine(0._c_float,-1._c_float)

    ! button: collapse
    str = "Collapse" // c_null_char
    if (igButton(c_loc(str),szero)) then
       do i = 1, nsys
          call collapse_system(i)
       end do
    end if
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str = "Collapse all systems in the tree" // c_null_char
       call igSetTooltip(c_loc(str))
    end if
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
    str = "Close" // c_null_char
    call igPushStyleColor_Vec4(ImGuiCol_Button,ColorDangerButton)
    if (igButton(c_loc(str),szero)) then
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
    call igPopStyleColor(1)
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str = "Close all visible systems" // c_null_char
       call igSetTooltip(c_loc(str))
    end if
    call igSameLine(0._c_float,-1._c_float)

    ! button: close all
    str = "Close All" // c_null_char
    call igPushStyleColor_Vec4(ImGuiCol_Button,ColorDangerButton)
    if (igButton(c_loc(str),szero)) then
       if (allocated(w%forceremove)) deallocate(w%forceremove)
       allocate(w%forceremove(nsys))
       do i = 1, nsys
          w%forceremove(i) = i
       end do
    end if
    call igPopStyleColor(1)
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str = "Close all systems" // c_null_char
       call igSetTooltip(c_loc(str))
    end if

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
             str = "✕##" // string(ic_closebutton) // "," // string(i) // c_null_char
             if (igSmallButton(c_loc(str))) w%forceremove = (/i/)
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
             if (sysc(i)%collapse < 0) then
                ! extend button for multi-seed entries
                if (sysc(i)%collapse == -1) then
                   str = "▶##" // string(ic_name) // "," // string(i) // c_null_char ! collapsed
                else
                   str = "▼##" // string(ic_name) // "," // string(i) // c_null_char ! extended
                end if
                if (igSmallButton(c_loc(str))) then
                   ! extend or collapse
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
            ! set as current system option
            strpop = "Set as current system" // c_null_char
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

    ! un-hide the dependents and set as extended
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
    use param, only: dirsep
    class(window), intent(inout), target :: w

    type(ImVec2) :: minsize, maxsize, inisize
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
    inisize%x = 800._c_float
    inisize%y = 480._c_float
    minsize%x = 0._c_float
    minsize%y = 0._c_float
    maxsize%x = 1e10_c_float
    maxsize%y = 1e10_c_float
    call igSetNextWindowSize(inisize,ImGuiCond_FirstUseEver)

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
                call add_systems_from_name(name,w%dialog_data%mol,isperm(w%dialog_data%isformat),readlastonly)
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
             if (allocated(outputb)) then
                lu = fopen_write(name,errstop=.false.)
                if (lu >= 0) then
                   write(lu,'(A)') outputb(1:lob)
                   call fclose(lu)
                else
                   call ferror('draw_dialog','could not open file for writing: ' // name,faterr,syntax=.true.)
                end if
             end if
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
  module subroutine draw_console_input(w)
    use gui_main, only: ColorHighlightText, tooltip_delay, sys, sysc, nsys, sys_init, g,&
       ColorDangerButton
    use gui_utils, only: igIsItemHovered_delayed
    use systemmod, only: sy
    use tools_io, only: string
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: sz, szero, szavail
    logical(c_bool) :: ldum, is_selected
    logical, save :: ttshown = .false.
    integer :: i
    ! the input buffer
    character(kind=c_char,len=:), allocatable, target, save :: inputb
    integer(c_size_t), parameter :: maxlib = 40000

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
    call igSameLine(0._c_float,-1._c_float)

    ! first line: clear button
    str1 = "Clear" // c_null_char
    if (igButton(c_loc(str1),szero)) &
       inputb(1:1) = c_null_char
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str1 = "Clear the input text" // c_null_char
       call igSetTooltip(c_loc(str1))
    end if

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
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str1 = "Set the current system (input commands are applied to it)" // c_null_char
       call igSetTooltip(c_loc(str1))
    end if

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
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str1 = "Set the reference field (input commands are applied to it)" // c_null_char
       call igSetTooltip(c_loc(str1))
    end if
    call igEndGroup()

    ! right-hand-side of lines 2 and 3: RUN button
    call igSameLine(0._c_float,-1._c_float)
    str1 = "RUN" // c_null_char
    call igPushStyleColor_Vec4(ImGuiCol_Button,ColorDangerButton)
    if (igButton(c_loc(str1),sz)) then
       write (*,*) "bleh RUN!"
    end if
    call igPopStyleColor(1)


    ! calculate sizes and draw the multiline
    call igGetContentRegionAvail(sz)
    str1 = "##inputmultiline" // c_null_char
    ldum = igInputTextMultiline(c_loc(str1),c_loc(inputb),maxlib,sz,&
       ImGuiInputTextFlags_AllowTabInput,c_null_ptr,c_null_ptr)

  end subroutine draw_console_input

  !> Draw the contents of the output console
  module subroutine draw_console_output(w)
    use gui_main, only: ColorHighlightText, tooltip_delay, g
    use gui_utils, only: igIsItemHovered_delayed
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1
    type(ImVec2) :: sz, szero
    logical(c_bool) :: ldum
    logical, save :: ttshown = .false.
    logical :: doscroll
    integer, save :: idsavedialog = 0

    ! check if the save dialog is still open
    if (idsavedialog > 0) then
       if (idsavedialog < 1 .or. idsavedialog > nwin) then
          idsavedialog = 0
       elseif (.not.win(idsavedialog)%isinit .or. .not.win(idsavedialog)%isopen) then
          idsavedialog = 0
       end if
    end if

    ! read new output, if available
    call read_output_unit()

    ! initialize
    szero%x = 0._c_float
    szero%y = 0._c_float

    ! first line: text
    str1 = "Output" // c_null_char
    call igTextColored(ColorHighlightText,c_loc(str1))
    call igSameLine(0._c_float,-1._c_float)

    ! first line: clear button
    str1 = "Clear" // c_null_char
    if (igButton(c_loc(str1),szero)) then
       outputb(1:1) = c_null_char
       lob = 0
    end if
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str1 = "Clear the output log" // c_null_char
       call igSetTooltip(c_loc(str1))
    end if
    call igSameLine(0._c_float,-1._c_float)

    ! first line: copy button
    str1 = "Copy" // c_null_char
    if (igButton(c_loc(str1),szero)) &
       call igSetClipboardText(c_loc(outputb))
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str1 = "Copy the output log to clipboard" // c_null_char
       call igSetTooltip(c_loc(str1))
    end if
    call igSameLine(0._c_float,-1._c_float)

    ! Renamed to BeginDisabled() / EndDisabled() and pushed on master.
    ! Added style.DisabledAlpha and ImGuiStyleVar_DisabledAlpha, defaulting to 0.60f.

    ! first line: save button
    str1 = "Save" // c_null_char
    call igBeginDisabled(logical(idsavedialog > 0,c_bool))
    if (igButton(c_loc(str1),szero)) &
       idsavedialog = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savelogfile)
    call igEndDisabled()
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str1 = "Save the output log to a file" // c_null_char
       call igSetTooltip(c_loc(str1))
    end if

    ! calculate sizes and draw the multiline (with dark background and border)
    call igGetContentRegionAvail(sz)
    call igPushStyleColor_Vec4(ImGuiCol_FrameBg,g%Style%Colors(ImGuiCol_WindowBg+1))
    call igPushStyleVar_Float(ImGuiStyleVar_FrameBorderSize,1._c_float)
    str1 = "##outputmultiline" // c_null_char
    ldum = igInputTextMultiline(c_loc(str1),c_loc(outputb),lob,sz,ImGuiInputTextFlags_ReadOnly,c_null_ptr,c_null_ptr)
    call igPopStyleVar(1)
    call igPopStyleColor(1)

    ! auto-scroll to the end if we have new output
    if (doscroll) then
       ldum = igBeginChild_Str(c_loc(str1),szero,.false._c_bool,ImGuiWindowFlags_None)
       call igSetScrollHereY(1._c_float)
       call igEndChild()
    end if

  contains
    ! read new output from the scratch LU uout
    subroutine read_output_unit()
      use gui_main, only: are_threads_running
      use tools_io, only: uout, getline_raw

      character(len=:), allocatable :: line
      integer(c_size_t) :: pos, lshift, ll
      integer :: idx
      logical :: ok

      ! do not scroll for now
      doscroll = .false.

      ! allocate the output buffer if not allocated
      if (.not.allocated(outputb)) then
         allocate(character(len=maxlob+1) :: outputb)
         outputb(1:1) = c_null_char
         lob = 0
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
            idx = index(outputb,new_line('a'))
            if (idx == 0) then
               lob = 0
            else
               outputb(1:lob-idx) = outputb(idx+1:lob)
               lob = lob - idx
            end if
         end if

         ! read the new output and rewind
         do while(getline_raw(uout,line))
            ll = len(line)
            if (ll > 0) then
               outputb(lob+1:lob+ll) = line(1:ll)
               lob = lob + ll
            end if
            outputb(lob+1:lob+1) = new_line('a')
            lob = lob + 1
            inquire(uout,pos=pos)
         end do
         outputb(lob+1:lob+1) = c_null_char
         rewind(uout)

         ! scroll to keep up with the new output
         doscroll = .true.
      end if

    end subroutine read_output_unit
  end subroutine draw_console_output

  !xx! private procedures

  ! Return the string for the tooltip shown by the tree window,
  ! corresponding to system i.
  function tree_tooltip_string(i) result(str)
    use crystalmod, only: pointgroup_info, holo_string
    use gui_main, only: sys, sysc, nsys, sys_init
    use tools_io, only: string
    use param, only: bohrtoa, maxzat, atmass, pcamu, bohr2cm
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
    str = "||" // trim(sysc(i)%seed%name) // "||" // new_line(str)
    if (sysc(i)%status == sys_init) then
       ! file and system type
       str = str // trim(sysc(i)%seed%file) // new_line(str)
       if (sys(i)%c%ismolecule) then
          if (sys(i)%c%nmol == 1) then
             str = str // "A molecule." // new_line(str)
          else
             str = str // "A molecular cluster with " // string(sys(i)%c%nmol) // " fragments." //&
                new_line(str)
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
          str = str // "." // new_line(str)
       elseif (sys(i)%c%nlvac == 2) then
          str = str // "A 1D periodic (polymer) structure." //&
             new_line(str)
       elseif (sys(i)%c%nlvac == 1) then
          str = str // "A 2D periodic (layered) structure." //&
             new_line(str)
       else
          str = str // "A crystal." //&
             new_line(str)
       end if
       str = str // new_line(str)

       ! number of atoms, electrons, molar mass
       str = str // string(sys(i)%c%ncel) // " atoms, " //&
          string(sys(i)%c%nneq) // " non-eq atoms, "//&
          string(sys(i)%c%nspc) // " species," // new_line(str)
       nelec = 0
       mass = 0d0
       do k = 1, sys(i)%c%nneq
          iz = sys(i)%c%spc(sys(i)%c%at(k)%is)%z
          if (iz >= maxzat .or. iz <= 0) cycle
          nelec = nelec + iz * sys(i)%c%at(k)%mult
          mass = mass + atmass(iz) * sys(i)%c%at(k)%mult
       end do
       str = str // string(nelec) // " electrons, " //&
          string(mass,'f',decimal=3) // " amu per cell" // new_line(str)
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
       str = str // new_line(str) // new_line(str)

       if (.not.sys(i)%c%ismolecule) then
          ! cell parameters, volume, density
          str = str // "a/b/c (Å): " // &
             string(sys(i)%c%aa(1)*bohrtoa,'f',decimal=4) // " " //&
             string(sys(i)%c%aa(2)*bohrtoa,'f',decimal=4) // " " //&
             string(sys(i)%c%aa(3)*bohrtoa,'f',decimal=4) //&
             new_line(str)
          str = str // "α/β/γ (°): " // &
             string(sys(i)%c%bb(1),'f',decimal=2) // " " //&
             string(sys(i)%c%bb(2),'f',decimal=2) // " " //&
             string(sys(i)%c%bb(3),'f',decimal=2) // " " //&
             new_line(str)
          str = str // "V (Å³): " // &
             string(sys(i)%c%omega*bohrtoa**3,'f',decimal=2) // new_line(str)
          dens = (mass*pcamu) / (sys(i)%c%omega*bohr2cm**3)
          str = str // "Density (g/cm³): " // string(dens,'f',decimal=3) // new_line(str) &
             // new_line(str)

          ! symmetry
          if (sys(i)%c%spgavail) then
             call pointgroup_info(sys(i)%c%spg%pointgroup_symbol,schpg,holo,laue)
             str = str // "Symmetry: " // &
                string(sys(i)%c%spg%international_symbol) // " (" //&
                string(sys(i)%c%spg%spacegroup_number) // "), " //&
                string(holo_string(holo)) // "," //&
                new_line(str)

             str = str // string(sys(i)%c%neqv) // " symm-ops, " //&
                string(sys(i)%c%ncv) // " cent-vecs" //&
                new_line(str)
          else
             str = str // "Symmetry info not available" // new_line(str)
          end if
          str = str // new_line(str)
       end if

       ! number of scalar fields
       str = str // string(sys(i)%nf) // " scalar fields loaded" // new_line(str)
    else
       ! not initialized
       str = str // "Not initialized" // new_line(str)
    end if
    str = str // c_null_char

  end function tree_tooltip_string

  ! the callback for the right-hand-side pane of the dialog
  subroutine dialog_user_callback(vFilter, vUserData, vCantContinue) bind(c)
    use gui_main, only: ColorHighlightText, tooltip_delay
    use gui_utils, only: igIsItemHovered_delayed
    use gui_interfaces_cimgui
    type(c_ptr), intent(in), value :: vFilter ! const char *
    type(c_ptr), value :: vUserData ! void *
    logical(c_bool) :: vCantContinue ! bool *

    character(kind=c_char,len=:), allocatable, target :: str, stropt
    type(dialog_userdata), pointer :: data
    logical(c_bool) :: ldum
    logical, save :: ttshown = .false.

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
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       str = "Show the OS hidden files and directories in this dialog" // c_null_char
       call igSetTooltip(c_loc(str))
    end if

    !! options specific to the open files dialog !!
    if (data%purpose == wpurp_dialog_openfiles) then
       str = "Read last structure only" // c_null_char
       ldum = igCheckbox(c_loc(str),data%readlastonly)
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
          str = "Read only the last structure in the file" // c_null_char
          call igSetTooltip(c_loc(str))
       end if
       call igNewLine()

       ! radio buttons for auto/crystal/molecule
       str = "Read structures as..." // c_null_char
       call igTextColored(ColorHighlightText,c_loc(str))
       str = "Auto-detect" // c_null_char
       ldum = igRadioButton_IntPtr(c_loc(str),data%mol,-1)
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
          str = "Auto-detect whether new structures are read as crystals or molecules" // c_null_char
          call igSetTooltip(c_loc(str))
       end if
       str = "Crystal" // c_null_char
       ldum = igRadioButton_IntPtr(c_loc(str),data%mol,0)
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
          str = "Force new structures to be read as crystals" // c_null_char
          call igSetTooltip(c_loc(str))
       end if
       str = "Molecule" // c_null_char
       ldum = igRadioButton_IntPtr(c_loc(str),data%mol,1)
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
          str = "Force new structures to be read as molecules" // c_null_char
          call igSetTooltip(c_loc(str))
       end if
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
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
          str = &
             "Force new structures read with a given file format, or auto-detect"//new_line('a')//&
             "from the extension"//c_null_char

          call igSetTooltip(c_loc(str))
       end if
       call igNewLine()
    end if

  end subroutine dialog_user_callback

end submodule proc
