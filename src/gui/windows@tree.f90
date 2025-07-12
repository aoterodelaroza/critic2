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

! Windows, tree & associated dialogs: load_field and scfplot windows.
submodule (windows) tree
  use interfaces_cimgui
  implicit none

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
  integer(c_int), parameter :: ic_NUMCOLUMNS = 19 ! keep up to date

  ! color for vibrations and fields in tree
  real(c_float), parameter :: rgba_vibrations(4) = (/0.84_c_float,0.86_c_float,0.00_c_float,1.00_c_float/)
  real(c_float), parameter :: rgba_fields(4) = (/0.00_c_float,1.00_c_float,0.50_c_float,1.00_c_float/)
  real(c_float), parameter :: rgba_reference(4) = (/0.890_c_float,0.706_c_float,0.129_c_float,1._c_float/)


  !xx! private procedures
  ! function tree_system_tooltip(i)
  ! function tree_field_tooltip_string(si,fj)

contains

  !> Draw the contents of a tree window
  module subroutine draw_tree(w)
    use interfaces_glfw, only: glfwGetTime
    use windows, only: win, iwin_view
    use keybindings, only: is_bind_event, BIND_TREE_REMOVE_SYSTEM_FIELD, BIND_TREE_MOVE_UP,&
       BIND_TREE_MOVE_DOWN
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button,&
       iw_text, iw_setposx_fromend, iw_calcwidth, iw_calcheight, iw_menuitem
    use gui_main, only: nsys, sys, sysc, sys_empty, sys_init, sys_ready,&
       sys_loaded_not_init, sys_initializing,&
       launch_initialization_thread, ColorTableCellBg,&
       kill_initialization_thread, system_shorten_names, remove_system, tooltip_delay,&
       ColorDangerButton, ColorFieldSelected, g, tree_select_updates_inpcon,&
       tree_select_updates_view, fontsize, ok_system
    use fieldmod, only: type_grid
    use tools_io, only: string, uout
    use types, only: realloc
    use param, only: bohrtoa, ifformat_as_resample, ifformat_as_ft_x, ifformat_as_ft_y,&
       ifformat_as_ft_z, ifformat_as_ft_xx, ifformat_as_ft_xy, ifformat_as_ft_xz,&
       ifformat_as_ft_yy, ifformat_as_ft_yz, ifformat_as_ft_zz, ifformat_as_ft_grad,&
       ifformat_as_ft_lap, ifformat_as_ft_pot
    use c_interface_module
    class(window), intent(inout), target :: w

    character(kind=c_char,len=1024), target :: txtinp
    character(kind=c_char,len=:), allocatable, target :: str, strpop, strpop2, zeroc, ch
    character(kind=c_char,len=:), allocatable :: tooltipstr
    type(ImVec2) :: szero, sz
    type(ImVec4) :: col4
    integer(c_int) :: flags, color, idir
    integer :: i, j, k, nshown, newsel, jsel, ll, id, iref, inext, iprev, iaux, nfreal
    logical(c_bool) :: ldum, isel
    type(c_ptr) :: ptrc
    type(ImGuiTableSortSpecs), pointer :: sortspecs
    type(ImGuiTableColumnSortSpecs), pointer :: colspecs
    logical :: hadenabledcolumn, reinit, isend, ok, found
    logical :: export, didtableselected, hadrow
    real(c_float) :: width, pos
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f

    type(c_ptr), save :: cfilter = c_null_ptr ! filter object (allocated first pass, never destroyed)
    logical, save :: ttshown = .false. ! tooltip flag
    integer(c_int), save :: iresample(3) = (/0,0,0/) ! for the grid resampling menu option
    integer(c_int), save :: shown_after_filter = 0 ! number of systems shown after the filter

    ! initialize
    hadenabledcolumn = .false.
    zeroc = "" // c_null_char
    tooltipstr = ""
    szero%x = 0
    szero%y = 0
    if (.not.allocated(w%iord)) then
       w%tree_sortcid = ic_id
       w%tree_sortdir = 1
       w%tree_selected = 1
       w%forceupdate = .true.
       w%forceresize = .true.
       w%forcesort = .true.
    end if

    ! update the tree based on time signals between dependent windows
    do i = 1, nsys
       ! if a system has changed fundamentally, the table needs an update (maybe)
       if (w%timelast_tree_update < sysc(i)%timelastchange_geometry) w%forceupdate = .true.
       ! if a system has been rebonded, the "nmol" column may have changed
       if (w%timelast_tree_resize < sysc(i)%timelastchange_rebond) w%forceresize = .true.
       if (w%timelast_tree_sort < sysc(i)%timelastchange_rebond) w%forcesort = .true.
    end do

    ! Tree options button
    export = .false.
    ldum = (iw_button("❇",popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft))
    if (ok) then
       ! info at the top
       ldum = iw_menuitem("Right-click header for more",enabled=.false.)
       call igSeparator()

       ! button: expand
       if (iw_menuitem("Expand All")) then
          do i = 1, nsys
             call expand_system(i)
          end do
       end if
       call iw_tooltip("Expand all systems in the tree (show all SCF iterations)",ttshown)

       ! button: collapse
       if (iw_menuitem("Collapse All")) then
          do i = 1, nsys
             call collapse_system(i)
          end do
       end if
       call iw_tooltip("Collapse all systems in the tree (hide all SCF iterations)",ttshown)

       ! button: export
       export = iw_menuitem("Export Tree Table")
       call iw_tooltip("Write the current tree to the output console in csv-style (for copying)",ttshown)

       ! button: plot
       if (iw_menuitem("Plot Tree Data...")) &
          iaux = stack_create_window(wintype_treeplot,.true.,idparent=w%id,orraise=-1)
       call iw_tooltip("Plot the tree data",ttshown)
       call igSeparator()

       ! close visible
       if (iw_menuitem("Close Visible")) then
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

       ! close all systems
       if (iw_menuitem("Close All")) then
          if (allocated(w%forceremove)) deallocate(w%forceremove)
          allocate(w%forceremove(nsys))
          do i = 1, nsys
             w%forceremove(i) = i
          end do
       end if
       call iw_tooltip("Close all systems",ttshown)

       call igEndPopup()
    end if
    call igSameLine(0._c_float,-1._c_float)

    ! text filter
    call igGetContentRegionAvail(sz) ! Clear ... x/x shown
    width = max(sz%x - iw_calcwidth(16+5,1) - g%Style%WindowPadding%x,1._c_float)
    if (.not.c_associated(cfilter)) &
       cfilter = ImGuiTextFilter_ImGuiTextFilter(c_loc(zeroc))
    str = "##treefilter" // c_null_char
    ldum = ImGuiTextFilter_Draw(cfilter,c_loc(str),width)
    call iw_tooltip("Filter systems by name in the list below. Use comma-separated fields&
       & and - for excluding. Example: inc1,inc2,-exc includes all systems&
       & with inc1 or inc2 and excludes systems with exc.",ttshown)
    if (iw_button("Clear",sameline=.true.)) then
       if (c_associated(cfilter)) &
          call ImGuiTextFilter_Clear(cfilter)
    end if
    call iw_tooltip("Clear the filter",ttshown)

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
          if (w%forceremove(k) == w%tree_selected) then
             jsel = 0
             do j = 1, size(w%iord,1)
                if (w%iord(j) == w%tree_selected) then
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
             call select_system(newsel,.false.)
          end if
          ! if we removed the system for the input console or the view, update
          if (w%forceremove(k) == win(iwin_console_input)%inpcon_selected) &
             win(iwin_console_input)%inpcon_selected = w%tree_selected
          if (w%forceremove(k) == win(iwin_view)%view_selected) &
             call win(iwin_view)%select_view(w%tree_selected)
       end do
       deallocate(w%forceremove)
       ! restart initialization if the threads were killed
       if (reinit) w%forceinit = .true.
       w%forceupdate = .true.
    end if
    if (w%forceupdate) &
       call w%update_tree()
    if (w%forcesort) &
       call w%sort_tree(w%tree_sortcid,w%tree_sortdir)
    if (w%forceinit) then
       call kill_initialization_thread()
       call launch_initialization_thread()
       w%forceinit = .false.
    end if
    nshown = size(w%iord,1)

    ok = (nshown > 1)
    if (.not.ok) &
       ok = (sysc(w%iord(1))%status /= sys_empty)
    ! final message in the header line
    if (ok) &
       call iw_text(" " // string(shown_after_filter) // "/" // string(nshown) // " shown",&
       sameline=.true.)

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

    ! next and previous systems
    didtableselected = .false.
    inext = 0
    iprev = 0

    str = "Structures##0,0" // c_null_char
    flags = ImGuiTableFlags_Borders
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_Hideable)
    flags = ior(flags,ImGuiTableFlags_Sortable)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    if (igBeginTable(c_loc(str),ic_NUMCOLUMNS,flags,szero,0._c_float)) then
       ! force resize if asked for
       if (w%forceresize) then
          call igTableSetColumnWidthAutoAll(igGetCurrentTable())
          w%forceresize = .false.
          w%timelast_tree_resize = glfwGetTime()
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
       width = max(4._c_float, fontsize%y + 2._c_float)
       call igTableSetupColumn(c_loc(str),flags,width,ic_closebutton)

       str = "(expand button)##0expandbutton" // c_null_char
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_expandbutton)

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

       str = "nat/ncel##0" // c_null_char
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
             w%tree_sortcid = colspecs%ColumnUserID
             w%tree_sortdir = colspecs%SortDirection
             if (sortspecs%SpecsDirty .and. nshown > 1) then
                w%forcesort = .true.
                sortspecs%SpecsDirty = .false.
             end if
          else
             w%tree_sortcid = ic_id
             w%tree_sortdir = 1
          end if
       end if

       ! draw the header
       call igTableHeadersRow()

       ! start the clipper
       clipper = ImGuiListClipper_ImGuiListClipper()
       call ImGuiListClipper_Begin(clipper,nshown,-1._c_float)

       ! draw the rows
       hadrow = .false.
       shown_after_filter = 0
       do while(ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do j = clipper_f%DisplayStart+1, clipper_f%DisplayEnd

             i = w%iord(j)
             if (sysc(i)%status == sys_empty .or. sysc(i)%hidden) cycle
             if (c_associated(cfilter)) then
                str = trim(sysc(i)%seed%name) // c_null_char
                if (.not.ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) cycle
             end if
             shown_after_filter = shown_after_filter + 1

             ! start defining the table row
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             hadenabledcolumn = .false.
             hadrow = .true.

             ! close button
             if (sysc(i)%status == sys_init) then
                if (igTableSetColumnIndex(ic_closebutton)) then
                   call igAlignTextToFramePadding()
                   str = "##1closebutton" // string(ic_closebutton) // "," // string(i) // c_null_char
                   if (my_CloseButton(c_loc(str),ColorDangerButton)) w%forceremove = (/i/)
                   if (igIsItemHovered(ImGuiHoveredFlags_None)) &
                      tooltipstr = "Close this system"
                end if
             end if

             ! expand button
             if (igTableSetColumnIndex(ic_expandbutton)) then
                if (sysc(i)%collapse < 0) then
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
                   if (igIsItemHovered(ImGuiHoveredFlags_None)) &
                      tooltipstr = "Expand this system"
                end if
             end if

             ! set background color for the name cell, if not selected
             if (sysc(i)%status >= sys_ready) then
                col4 = ImVec4(ColorTableCellBg(1,sys(i)%c%iperiod),ColorTableCellBg(2,sys(i)%c%iperiod),&
                   ColorTableCellBg(3,sys(i)%c%iperiod),ColorTableCellBg(4,sys(i)%c%iperiod))
                color = igGetColorU32_Vec4(col4)
                call igTableSetBgColor(ImGuiTableBgTarget_CellBg, color, ic_name)
             end if

             ! ID column
             if (igTableSetColumnIndex(ic_id)) then
                str = string(i)
                call write_maybe_selectable(i,tooltipstr)
                ok = (sysc(i)%status >= sys_ready)
                if (ok) ok = sys(i)%c%vib%hasvibs
                if (ok) then
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export,&
                      rgba=rgba_vibrations)
                else
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
             end if

             ! name
             if (igTableSetColumnIndex(ic_name)) then
                ! selectable
                call write_maybe_selectable(i,tooltipstr)

                ! expand button
                if (sysc(i)%showfields) then
                   ch = "▼"
                else
                   ch = "▶"
                end if
                pos = igGetCursorPosX()
                call igSetCursorPosX(pos + g%Style%FramePadding%x)

                nfreal = 0
                do k = 1, sys(i)%nf
                   if (sys(i)%f(k)%isinit) nfreal = nfreal + 1
                end do
                if (nfreal > 0) then
                   call iw_text(ch,rgba=rgba_fields)
                else
                   call iw_text(ch)
                end if
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
                      isel = (w%tree_selected==i) .and. (sys(i)%iref == k)
                      call igPushStyleColor_Vec4(ImGuiCol_Header,ColorFieldSelected)
                      flags = ImGuiSelectableFlags_SpanAllColumns
                      if (igSelectable_Bool(c_loc(str),isel,flags,szero)) then
                         call select_system(i,.false.)
                         call sys(i)%set_reference(k,.false.)
                      end if
                      call igPopStyleColor(1)

                      ! right click to open the field context menu
                      if (igBeginPopupContextItem(c_loc(str),ImGuiPopupFlags_MouseButtonRight)) then
                         ! remove option (fields)
                         if (k > 0) then
                            if (iw_menuitem("Remove")) &
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
                            if (igInputText(c_loc(strpop2),c_loc(txtinp),1023_c_size_t,flags,c_null_funptr,c_null_ptr)) then
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
                         if (iw_menuitem("Duplicate")) then
                            id = sys(i)%getfieldnum()
                            call sys(i)%field_copy(k,id)
                            sys(i)%f(id)%id = id
                            sys(i)%f(id)%name = trim(sys(i)%f(k)%name)
                         end if
                         call iw_tooltip("Load a copy of this field as a new field",ttshown)

                         ! grid calculation options
                         if (sys(i)%f(k)%type == type_grid) then
                            strpop = "Load Fourier-Transformed Grid" // c_null_char
                            if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
                               ldum = iw_menuitem("[First derivatives]",enabled=.false.)
                               if (iw_menuitem("x")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, x-derivative of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_x)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the x-component of &
                                  &this field's gradient",ttshown)
                               if (iw_menuitem("y")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, y-derivative of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_y)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the y-component of &
                                  &this field's gradient",ttshown)
                               if (iw_menuitem("z")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, z-derivative of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_z)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the z-component of &
                                  &this field's gradient",ttshown)
                               if (iw_menuitem("Gradient Norm")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, grad of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_grad)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the norm of &
                                  &this field's gradient",ttshown)
                               call igSeparator()

                               ldum = iw_menuitem("[Second derivatives]",enabled=.false.)
                               if (iw_menuitem("xx")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, xx-derivative of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_xx)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the xx component of &
                                  &this field's Hessian matrix",ttshown)
                               if (iw_menuitem("xy")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, xy-derivative of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_xy)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the xy component of &
                                  &this field's Hessian matrix",ttshown)
                               if (iw_menuitem("xz")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, xz-derivative of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_xz)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the xz component of &
                                  &this field's Hessian matrix",ttshown)
                               if (iw_menuitem("yy")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, yy-derivative of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_yy)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the yy component of &
                                  &this field's Hessian matrix",ttshown)
                               if (iw_menuitem("yz")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, yz-derivative of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_yz)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the yz component of &
                                  &this field's Hessian matrix",ttshown)
                               if (iw_menuitem("zz")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,"<generated>, zz-derivative of $" // string(k),&
                                     sys(i)%f(k)%grid,ifformat_as_ft_zz)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the zz component of &
                                  &this field's Hessian matrix",ttshown)
                               if (iw_menuitem("Laplacian")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, lap of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_lap)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the Laplacian of this field",ttshown)

                               call igSeparator()
                               if (iw_menuitem("Potential")) then
                                  id = sys(i)%getfieldnum()
                                  call sys(i)%f(id)%load_as_fftgrid(sys(i)%c,id,&
                                     "<generated>, potential of $" // string(k),sys(i)%f(k)%grid,ifformat_as_ft_pot)
                               end if
                               call iw_tooltip("Load a new grid field using FFT as the potential that generates&
                                  &this field (via Poisson's equation)",ttshown)
                               call igEndMenu()
                            end if

                            strpop = "Load Resampled Grid" // c_null_char
                            if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
                               flags = ImGuiInputTextFlags_None
                               strpop2 = "New Size##resamplefieldmenunewsize" // c_null_char
                               call igSetNextItemWidth(iw_calcwidth(4*3,2))
                               ldum = igInputInt3(c_loc(strpop2),iresample,flags)

                               if (iw_menuitem("OK##resamplefieldmenuok")) then
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
                      if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
                         if (igIsMouseHoveringRect(g%LastItemData%NavRect%min,&
                            g%LastItemData%NavRect%max,.false._c_bool)) then
                            call tree_field_tooltip_string(i,k)
                         end if
                      end if

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
                call write_maybe_selectable(i,tooltipstr)
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
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_v)) then ! volume
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%omega*bohrtoa**3,'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_vmol)) then ! volume per molecule
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%omega*bohrtoa**3/sys(i)%c%nmol,'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_nneq)) then ! nneq
                   str = string(sys(i)%c%nneq)
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_ncel)) then ! ncel
                   str = string(sys(i)%c%ncel)
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_nmol)) then ! nmol
                   str = string(sys(i)%c%nmol)
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_a)) then ! a
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%aa(1)*bohrtoa,'f',decimal=4)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_b)) then ! b
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%aa(2)*bohrtoa,'f',decimal=4)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_c)) then ! c
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%aa(3)*bohrtoa,'f',decimal=4)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_alpha)) then ! alpha
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%bb(1),'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_beta)) then ! beta
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%bb(2),'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_gamma)) then ! gamma
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%bb(3),'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_emol)) then ! energy/nmol
                   if (sysc(i)%seed%energy /= huge(1d0)) then
                      str = string(sysc(i)%seed%energy/sys(i)%c%nmol,'f',decimal=8)
                   else
                      str = "n/a"
                   end if
                   call write_maybe_selectable(i,tooltipstr)
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
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,copy_to_output=export)
                end if
             end if
             ! enter new line
             if (export) write (uout,*)
          end do ! clipper indices

          if (.not.hadrow) then
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             if (igTableSetColumnIndex(ic_name)) then
                call igAlignTextToFramePadding()
                call iw_text("No systems loaded")
             end if
             hadrow = .true.
          end if
       end do ! clipper step
       call ImGuiListClipper_End(clipper)

       ! process the keybindings
       ok = is_bind_event(BIND_TREE_REMOVE_SYSTEM_FIELD)
       ok = ok .and. igIsWindowFocused(ImGuiFocusedFlags_None)
       if (ok) then
          jsel = w%tree_selected
          iref = sys(jsel)%iref
          ok = ok_system(jsel,sys_init)
          if (ok) ok = sysc(jsel)%showfields
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

    ! process the keybindings
    !! up and down the tree
    if (is_bind_event(BIND_TREE_MOVE_UP)) then
       if (iprev > 0) &
          w%forceselect = iprev
    elseif (is_bind_event(BIND_TREE_MOVE_DOWN)) then
       if (inext > 0) &
          w%forceselect = inext
    end if

    ! if exporting, read the export command
    if (export) &
       ldum = win(iwin_console_input)%read_output_ci(.true.,"[Table export]")

    ! clean up
    ! call ImGuiTextFilter_destroy(cfilter)

  contains

    subroutine write_maybe_selectable(isys,tooltipstr)
      use gui_main, only: are_threads_running, duplicate_system, reread_system_from_file
      use keybindings, only: BIND_GEOMETRY
      use utils, only: iw_text
      use global, only: iunit, iunit_bohr, iunit_ang
      use tools_io, only: uout
      integer, intent(in) :: isys
      character(kind=c_char,len=:), allocatable, intent(in) :: tooltipstr

      integer :: k, idx, iaux
      real(c_float) :: pos
      integer(c_int) :: flags, ll, isyscollapse, idum
      logical(c_bool) :: selected
      logical :: enabled, enabled_no_threads
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
      selected = (w%tree_selected==isys)
      strl = "##selectable" // string(isys) // c_null_char
      ok = igSelectable_Bool(c_loc(strl),selected,flags,szero)
      ok = ok .or. (w%forceselect == isys)
      if (ok) then
         call select_system(isys,.false.)
         if (w%forceselect > 0) then
            w%forceselect = 0
            call igSetKeyboardFocusHere(0)
         end if
         if (igIsMouseDoubleClicked(ImGuiPopupFlags_MouseButtonLeft)) &
            sysc(isys)%showfields = .not.sysc(isys)%showfields
      end if
      call igSameLine(0._c_float,-1._c_float)
      call igSetCursorPosX(pos)

      ! update the iprev and inext
      if (inext==0.and.didtableselected) inext = isys
      if (isys == w%tree_selected) didtableselected = .true.
      if (.not.didtableselected) iprev = isys

      enabled = (sysc(isys)%status == sys_init)
      enabled_no_threads = enabled.and..not.are_threads_running()

      ! right click to open the context menu
      if (igBeginPopupContextItem(c_loc(strl),ImGuiPopupFlags_MouseButtonRight)) then
         ! scf energy plot
         if (sysc(isys)%collapse /= 0) then
            if (iw_menuitem("Plot SCF Iterations...",enabled=enabled_no_threads)) then
               if (sysc(isys)%collapse < 0) then
                  isyscollapse = isys
               else
                  isyscollapse = abs(sysc(isys)%collapse)
               end if

               if (sysc(isyscollapse)%status == sys_init) &
                  iaux = stack_create_window(wintype_scfplot,.true.,isys=isyscollapse,orraise=-1)
            end if
            call iw_tooltip("Plot the energy and other properties as a function of SCF cycle iterations",ttshown)
         end if

         ! geometry submenu (system)
         strpop = "System" // c_null_char
         if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
            ! Geometry
            if (iw_menuitem("View/Edit Geometry...",BIND_GEOMETRY,enabled=enabled)) &
               idum = stack_create_window(wintype_geometry,.true.,isys=isys,orraise=-1)
            call iw_tooltip("View and edit the atomic positions, bonds, etc.",ttshown)

            ! describe this system in the console output
            if (iw_menuitem("Describe in Output Window",enabled=enabled)) then
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
            call iw_tooltip("Print a detailed description of this system in the output console",ttshown)

            ! rebond
            if (iw_menuitem("Recalculate Bonds...",enabled=enabled_no_threads)) &
               idum = stack_create_window(wintype_rebond,.true.,isys=isys,orraise=-1)
            call iw_tooltip("Recalculate the bonds in this system",ttshown)

            call igEndMenu()
         end if

         ! fields submenu (system)
         strpop = "Fields" // c_null_char
         if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
            ! load field
            if (iw_menuitem("Load Field...",enabled=enabled)) &
               iaux = stack_create_window(wintype_load_field,.true.,isys=isys,orraise=-1)
            call iw_tooltip("Load a scalar field for this system",ttshown)

            call igEndMenu()
         end if

         ! vibrations submenu (system)
         strpop = "Vibrations" // c_null_char
         if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
            ! load vibration data
            if (iw_menuitem("Load Vibration Data...",enabled=enabled)) then
               iaux = stack_create_window(wintype_dialog,.true.,purpose=wpurp_dialog_openvibfile,&
                  isys=isys,orraise=-1)
            end if
            call iw_tooltip("Load vibration data from an external file for this system",ttshown)

            ! clear vibration data
            if (iw_menuitem("Clear Vibration Data",enabled=enabled)) &
               call sys(isys)%c%vib%end()
            call iw_tooltip("Clear the vibration data for this system",ttshown)

            call igEndMenu()
         end if

         ! rename option (system)
         strpop = "Rename" // c_null_char
         if (igBeginMenu(c_loc(strpop),.true._c_bool)) then
            strpop2 = "##inputrename" // c_null_char
            txtinp = trim(adjustl(sysc(isys)%seed%name)) // c_null_char
            call igSetKeyboardFocusHere(0_c_int)
            flags = ImGuiInputTextFlags_EnterReturnsTrue
            if (igInputText(c_loc(strpop2),c_loc(txtinp),1023_c_size_t,flags,c_null_funptr,c_null_ptr)) then
               ll = index(txtinp,c_null_char)
               sysc(isys)%seed%name = txtinp(1:ll-1)
               sysc(isys)%renamed = .true.
               call igCloseCurrentPopup()
            end if
            call igEndMenu()
         end if
         call iw_tooltip("Rename this system",ttshown)

         ! set as current system option
         if (iw_menuitem("Set as Current System",enabled=enabled_no_threads)) &
            call select_system(isys,.true.)
         call iw_tooltip("Set this system as the current system",ttshown)

         ! set as current system option
         if (iw_menuitem("Display in New View",enabled=enabled_no_threads)) then
            idum = stack_create_window(wintype_view,.true.,purpose=wpurp_view_alternate)
            win(idum)%sc = sysc(isys)%sc
            win(idum)%view_selected = isys
         end if
         call iw_tooltip("Display this system in a new view window",ttshown)

         ! duplicate system
         if (iw_menuitem("Duplicate",enabled=enabled_no_threads)) &
            call duplicate_system(isys)
         call iw_tooltip("Create a copy of this system",ttshown)

         ! reopen from file
         if (iw_menuitem("Reopen from File",enabled=enabled_no_threads)) then
            call reread_system_from_file(isys)
         end if
         call iw_tooltip("Read the file for this system and reopen it (only last structure is read)",ttshown)

         ! remove option (system)
         ok = enabled
         if (ok) ok = (sys(isys)%nf > 0)
         if (ok) ok = any(sys(isys)%f(1:sys(isys)%nf)%isinit)
         if (ok) then
            if (iw_menuitem("Remove All Fields")) then
               do k = 1, sys(isys)%nf
                  if (.not.sys(isys)%f(k)%isinit) cycle
                  call sys(isys)%unload_field(k)
               end do
            end if
            call iw_tooltip("Remove all fields in this system",ttshown)
         end if

         ! remove option (system)
         if (iw_menuitem("Close",enabled=enabled)) &
            w%forceremove = (/isys/)
         call iw_tooltip("Close this system",ttshown)

         call igEndPopup()
      end if

      ! delayed tooltip with info about the system
      if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
         if (igIsMouseHoveringRect(g%LastItemData%NavRect%min,g%LastItemData%NavRect%max,.false._c_bool)) then
            if (len(tooltipstr) > 0) then
               strl = tooltipstr // c_null_char
               call igSetTooltip(c_loc(strl))
            else
               call tree_system_tooltip(isys)
            end if
         end if
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
      if (w%tree_selected >= 1 .and. w%tree_selected <= nsys) then
         if (sysc(w%tree_selected)%collapse == i) call select_system(i,.true.)
      end if
      w%forceupdate = .true.

    end subroutine collapse_system

    ! Select table system i. If force, change inpcon and view selections
    ! even if disabled by the preferences.
    subroutine select_system(i,force)
      integer, intent(in) :: i
      logical, intent(in) :: force

      w%tree_selected = i
      if (tree_select_updates_inpcon .or. force) &
         win(iwin_console_input)%inpcon_selected = i
      if (tree_select_updates_view .or. force) &
         call win(iwin_view)%select_view(i)

    end subroutine select_system

  end subroutine draw_tree

  ! Update the table rows by building a new row index array
  ! (iord). Only the systems that are not empty are pointed by
  ! iord. This routine is used when the systems change.
  module subroutine update_tree(w)
    use interfaces_glfw, only: glfwGetTime
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
    w%timelast_tree_update = glfwGetTime()

  end subroutine update_tree

  ! Sort the table row order by column cid and in direction dir
  ! (ascending=1, descending=2). Modifies the w%iord.
  module subroutine sort_tree(w,cid,dir)
    use interfaces_glfw, only: glfwGetTime
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

    ! update the time
    w%timelast_tree_sort = glfwGetTime()

  end subroutine sort_tree

  !xx! private procedures

  ! Return the string for the tooltip shown by the tree window,
  ! corresponding to system i.
  subroutine tree_system_tooltip(i)
    use utils, only: iw_text
    use crystalmod, only: pointgroup_info, holo_string, iperiod_3d_crystal,&
       iperiod_3d_layered, iperiod_3d_chain, iperiod_3d_molecular,&
       iperiod_2d, iperiod_1d, iperiod_0d, iperiod_mol_single,&
       iperiod_mol_cluster
    use gui_main, only: sys, sysc, sys_init, ColorTableCellBg, ok_system, sys_empty
    use tools_io, only: string
    use param, only: bohrtoa, maxzat, atmass, pcamu, bohrtocm
    use tools_math, only: gcd
    integer, intent(in) :: i

    character(kind=c_char,len=:), allocatable, target :: str
    integer, allocatable :: nis(:)
    integer :: k, iz
    real*8 :: maxdv, mass, dens
    real(c_float) :: rgba(4)
    integer :: nelec
    character(len=3) :: schpg
    integer :: holo, laue
    integer :: izp0

    str = ""
    if (.not.ok_system(i,sys_empty)) return

    ! begin the tooltip
    call igBeginTooltip()

    ! file
    call iw_text("[",danger=.true.)
    call iw_text(trim(sysc(i)%seed%name),highlight=.true.,sameline_nospace=.true.)
    call iw_text("]",danger=.true.,sameline_nospace=.true.)

    if (sysc(i)%status == sys_init) then
       ! file and system type
       str = trim(adjustl(sysc(i)%seed%file))
       call iw_text("File: ",highlight=.true.)
       call iw_text(str,sameline_nospace=.true.)

       if (sys(i)%c%iperiod == iperiod_3d_crystal) then
          str = "A crystal"
       elseif (sys(i)%c%iperiod == iperiod_3d_layered) then
          str = "A layered crystal"
       elseif (sys(i)%c%iperiod == iperiod_3d_chain) then
          str = "A crystal with one-dimensional bonded chains"
       elseif (sys(i)%c%iperiod == iperiod_3d_molecular) then
          str = "A molecular crystal with Z=" // string(sys(i)%c%nmol)
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
       elseif (sys(i)%c%iperiod == iperiod_2d) then
          str = "A surface (2D periodic structure)"
       elseif (sys(i)%c%iperiod == iperiod_1d) then
          str = "A one-dimensional periodic structure"
       elseif (sys(i)%c%iperiod == iperiod_0d) then
          str = "A molecule in a big periodic box"
       elseif (sys(i)%c%iperiod == iperiod_mol_single) then
          str = "A molecule"
       elseif (sys(i)%c%iperiod == iperiod_mol_cluster) then
          str = "A molecular cluster with " // string(sys(i)%c%nmol) // " fragments"
       end if
       rgba = ColorTableCellBg(:,sys(i)%c%iperiod)
       rgba(4) = 1.0_c_float
       call iw_text(str,rgba=rgba)

       ! number of atoms
       call iw_text("")
       call iw_text("Atoms: ",highlight=.true.)
       str = string(sys(i)%c%ncel) // " ("
       if (.not.sys(i)%c%ismolecule) &
          str = str // string(sys(i)%c%nneq) // " non-eq, "
       str = str // string(sys(i)%c%nspc) // " species)"
       call iw_text(str,sameline_nospace=.true.)

       ! electrons, molar mass
       nelec = 0
       mass = 0d0
       do k = 1, sys(i)%c%nneq
          iz = sys(i)%c%spc(sys(i)%c%at(k)%is)%z
          if (iz >= maxzat .or. iz <= 0) cycle
          nelec = nelec + iz * sys(i)%c%at(k)%mult
          mass = mass + atmass(iz) * sys(i)%c%at(k)%mult
       end do
       call iw_text("Electrons: ",highlight=.true.)
       call iw_text(string(nelec),sameline_nospace=.true.)
       call iw_text(" Mass (amu): ",highlight=.true.,sameline=.true.)
       call iw_text(string(mass,'f',decimal=3),sameline_nospace=.true.)

       ! empirical formula
       allocate(nis(sys(i)%c%nspc))
       nis = 0
       do k = 1, sys(i)%c%nneq
          nis(sys(i)%c%at(k)%is) = nis(sys(i)%c%at(k)%is) + sys(i)%c%at(k)%mult
       end do
       maxdv = gcd(nis,sys(i)%c%nspc)
       call iw_text("Formula: ",highlight=.true.)
       str = ""
       do k = 1, sys(i)%c%nspc
          str = str // string(sys(i)%c%spc(k)%name) // string(nint(nis(k)/maxdv)) // " "
       end do
       if (len(str) > 30) str = str(1:27) // "..."
       call iw_text(str,sameline_nospace=.true.)

       if (.not.sys(i)%c%ismolecule) then
          ! cell parameters, volume, density, energy, pressure
          call iw_text("")
          call iw_text("a/b/c (Å): ",highlight=.true.)
          str = string(sys(i)%c%aa(1)*bohrtoa,'f',decimal=4) // " " //&
             string(sys(i)%c%aa(2)*bohrtoa,'f',decimal=4) // " " //&
             string(sys(i)%c%aa(3)*bohrtoa,'f',decimal=4)
          call iw_text(str,sameline_nospace=.true.)

          call iw_text("α/β/γ (°): ",highlight=.true.)
          str = string(sys(i)%c%bb(1),'f',decimal=2) // " " //&
             string(sys(i)%c%bb(2),'f',decimal=2) // " " //&
             string(sys(i)%c%bb(3),'f',decimal=2)
          call iw_text(str,sameline_nospace=.true.)

          call iw_text("V (Å³): ",highlight=.true.)
          str = string(sys(i)%c%omega*bohrtoa**3,'f',decimal=2)
          call iw_text(str,sameline_nospace=.true.)

          call iw_text("Density (g/cm³): ",highlight=.true.)
          dens = (mass*pcamu) / (sys(i)%c%omega*bohrtocm**3)
          call iw_text(string(dens,'f',decimal=3),sameline_nospace=.true.)

          if (sysc(i)%seed%energy /= huge(1d0)) then
             call iw_text("Energy (Ha): ",highlight=.true.)
             call iw_text(string(sysc(i)%seed%energy,'f',decimal=8),sameline_nospace=.true.)
          end if
          if (sysc(i)%seed%pressure /= huge(1d0)) then
             call iw_text("Pressure (GPa): ",highlight=.true.)
             call iw_text(string(sysc(i)%seed%pressure,'f',decimal=2),sameline_nospace=.true.)
          end if

          ! symmetry
          if (sys(i)%c%spgavail) then
             call pointgroup_info(sys(i)%c%spg%pointgroup_symbol,schpg,holo,laue)
             call iw_text("")
             call iw_text("Space group: ",highlight=.true.)
             str = string(sys(i)%c%spg%international_symbol) // " (" //&
                string(sys(i)%c%spg%spacegroup_number) // ", " //&
                string(holo_string(holo)) // ")"
             call iw_text(str,sameline_nospace=.true.)
          else
             call iw_text("")
             call iw_text("Space group: ",highlight=.true.)
             call iw_text("not available",sameline_nospace=.true.)
          end if
       end if

       ! number of scalar fields
       if (sys(i)%nf > 0) then
          call iw_text("")
          call iw_text(string(sys(i)%nf) // " scalar fields loaded",rgba=rgba_fields)
       end if

       ! vibrations
       if (sys(i)%c%vib%hasvibs) then
          call iw_text("")
          call iw_text("Vibration data available",rgba=rgba_vibrations)
       end if
    else
       ! not initialized
       call iw_text("")
       call iw_text("Not initialized",danger=.true.)
    end if
    call iw_text("")
    call iw_text("[Right-click for options]",highlight=.true.)

    call igEndTooltip()

  end subroutine tree_system_tooltip

  ! Return the string for the tooltip shown by the tree window,
  ! corresponding to system si and field fj
  subroutine tree_field_tooltip_string(si,fj)
    use utils, only: iw_text
    use gui_main, only: sys, sys_init, ok_system
    use fieldmod, only: field, type_uninit, type_promol, type_grid, type_wien,&
       type_elk, type_pi, type_wfn, type_dftb, type_promol_frag, type_ghost
    use wfn_private, only: molden_type_psi4, molden_type_orca, molden_type_adf_sto,&
       wfn_rhf, wfn_uhf, wfn_frac
    use grid3mod, only: mode_nearest, mode_trilinear, mode_trispline, mode_tricubic,&
       mode_smr
    use tools_io, only: string, nameguess
    use param, only: maxzat0
    integer, intent(in) :: si, fj

    character(kind=c_char,len=:), allocatable, target :: str, aux
    integer :: i, nal
    type(field), pointer :: f

    str = ""
    if (.not.ok_system(si,sys_init)) return
    if (fj < 0 .or. fj > sys(si)%nf) return
    if (.not. sys(si)%f(fj)%isinit) return
    f => sys(si)%f(fj)

    ! begin the tooltip
    call igBeginTooltip()

    ! file
    call iw_text("[",danger=.true.)
    call iw_text(trim(f%name),highlight=.true.,sameline_nospace=.true.)
    call iw_text("]",danger=.true.,sameline_nospace=.true.)

    ! field type
    select case (f%type)
    case (type_uninit)
       str = "???"
    case (type_promol)
       str = "Promolecular density"
    case (type_grid)
       str = "Grid field, with " // string(f%grid%n(1)) // "x" // string(f%grid%n(2)) // "x" //&
          string(f%grid%n(3)) // " points"
    case (type_wien)
       str = "WIEN2k spherical harmonic + plane waves"
    case (type_elk)
       str = "Elk spherical harmonic + plane waves"
    case (type_pi)
       str = "aiPI atomic radial functions"
    case (type_wfn)
       str = "Molecular wavefunction "
       if (f%wfn%issto) then
          str = str // "(Slater-type orbitals)"
       else
          str = str // "(Gaussian-type orbitals)"
       end if
    case (type_dftb)
       str = "DFTB+, linear combination of atomic orbitals"
    case (type_promol_frag)
       str = "Promolecular with fragment"
    case (type_ghost)
       str = "Ghost field"
    case default
       str = "???"
    end select
    call iw_text(str,rgba=rgba_fields)

    ! reference field
    if (sys(si)%iref == fj) &
       call iw_text("Reference field for this system",rgba=rgba_reference)

    ! aliases
    call sys(si)%aliasstring(fj,nal,aux)
    call iw_text("Names: ",highlight=.true.)
    call iw_text(trim(adjustl(aux)),sameline_nospace=.true.)

    ! field details
    select case (f%type)
    case (type_grid)
       if (f%grid%isqe) then
          call iw_text("Plane-wave Kohn-Sham states available: ",highlight=.true.)
          call iw_text("spin=" // string(f%grid%qe%nspin) // ", k-points=" // string(f%grid%qe%nks) //&
             ", bands=" // string(f%grid%qe%nbnd),sameline_nospace=.true.)
       end if
       if (f%grid%iswan) then
          call iw_text("Wannier function info available: ",highlight=.true.)
          call iw_text(string(f%grid%qe%nk(1)) // "x" // string(f%grid%qe%nk(2)) // "x" //&
             string(f%grid%qe%nk(3)),sameline_nospace=.true.)
       end if

       call iw_text("Longest Voronoi vector (bohr): ",highlight=.true.)
       call iw_text(string(f%grid%dmax,'f',decimal=5),sameline_nospace=.true.)
       call iw_text("Grid integral: ",highlight=.true.)
       call iw_text(string(sum(f%grid%f) * f%c%omega / real(product(f%grid%n),4),'f',decimal=8),sameline_nospace=.true.)
       call iw_text("Minimum value: ",highlight=.true.)
       call iw_text(string(minval(f%grid%f),'e',decimal=4),sameline_nospace=.true.)
       call iw_text("Average: ",highlight=.true.)
       call iw_text(string(sum(f%grid%f) / real(product(f%grid%n),4),'e',decimal=8),sameline_nospace=.true.)
       call iw_text("Maximum value: ",highlight=.true.)
       call iw_text(string(maxval(f%grid%f),'e',decimal=4),sameline_nospace=.true.)
       call iw_text("Interpolation mode: ",highlight=.true.)
       select case (f%grid%mode)
       case(mode_nearest)
          str = "nearest grid node value"
       case(mode_trilinear)
          str = "tri-linear"
       case(mode_trispline)
          str = "tri-spline"
       case(mode_tricubic)
          str = "tri-cubic"
       case(mode_smr)
          str = "smooth all-electron density with " // string(f%grid%smr_nenv) //&
             " stencil nodes and " // string(f%grid%smr_fdmax,'f',decimal=4) //&
             " smoothing dmax factor"
       end select
       call iw_text(str,sameline_nospace=.true.)

    case (type_wien)
       call iw_text("Complex: ",highlight=.true.)
       call iw_text(string(f%wien%cmpl),sameline_nospace=.true.)
       call iw_text("Spherical harmonics expansion LMmax: ",highlight=.true.)
       call iw_text(string(size(f%wien%lm,2)),sameline_nospace=.true.)
       call iw_text("Max. points in radial grid: ",highlight=.true.)
       call iw_text(string(size(f%wien%slm,1)),sameline_nospace=.true.)
       call iw_text("Total number of plane waves: ",highlight=.true.)
       call iw_text(string(f%wien%nwav),sameline_nospace=.true.)
       call iw_text("Density-style normalization: ",highlight=.true.)
       call iw_text(string(f%wien%cnorm),sameline_nospace=.true.)

    case (type_elk)
       call iw_text("Number of LM pairs: ",highlight=.true.)
       call iw_text(string(size(f%elk%rhomt,2)),sameline_nospace=.true.)
       call iw_text("Max. points in radial grid: ",highlight=.true.)
       call iw_text(string(size(f%elk%rhomt,1)),sameline_nospace=.true.)
       call iw_text("Total number of plane waves: ",highlight=.true.)
       call iw_text(string(size(f%elk%rhok)),sameline_nospace=.true.)

    case (type_pi)
       call iw_text("Exact calculation: ",highlight=.true.)
       call iw_text(string(f%exact),sameline_nospace=.true.)
       do i = 1, f%c%nspc
          call iw_text(string(f%c%spc(f%c%at(i)%is)%name,length=5) // ": ")
          call iw_text(string(f%pi%bas(i)%piname),sameline_nospace=.true.)
       end do

    case (type_wfn)
       if (f%wfn%molden_type == molden_type_orca) then
          call iw_text("Molden file dialect: ",highlight=.true.)
          call iw_text("orca",sameline_nospace=.true.)
       else if (f%wfn%molden_type == molden_type_psi4) then
          call iw_text("Molden file dialect: ",highlight=.true.)
          call iw_text("psi4",sameline_nospace=.true.)
       else if (f%wfn%molden_type == molden_type_adf_sto) then
          call iw_text("Molden file dialect: ",highlight=.true.)
          call iw_text("ADF (STOs)",sameline_nospace=.true.)
       end if
       if (f%wfn%wfntyp == wfn_rhf) then
          call iw_text("Wavefunction type: ",highlight=.true.)
          call iw_text("restricted",sameline_nospace=.true.)
          call iw_text("Number of MOs (total): ",highlight=.true.)
          call iw_text(string(f%wfn%nmoall),sameline_nospace=.true.)
          call iw_text("Number of MOs (occupied): ",highlight=.true.)
          call iw_text(string(f%wfn%nmoocc),sameline_nospace=.true.)
       elseif (f%wfn%wfntyp == wfn_uhf) then
          call iw_text("Wavefunction type: ",highlight=.true.)
          call iw_text("unrestricted",sameline_nospace=.true.)
          call iw_text("Number of MOs (total): ",highlight=.true.)
          call iw_text(string(f%wfn%nmoall) //&
             " (alpha=" // string(f%wfn%nalpha+f%wfn%nalpha_virt) // ",beta=" //&
             string(f%wfn%nmoall-(f%wfn%nalpha+f%wfn%nalpha_virt)),sameline_nospace=.true.)
          call iw_text("Number of MOs (occupied): ",highlight=.true.)
          call iw_text(string(f%wfn%nmoocc) // " (alpha=" // string(f%wfn%nalpha) //&
             ",beta=" // string(f%wfn%nmoocc-f%wfn%nalpha),sameline_nospace=.true.)
       elseif (f%wfn%wfntyp == wfn_frac) then
          call iw_text("Wavefunction type: ",highlight=.true.)
          call iw_text("fractional occupation",sameline_nospace=.true.)
          call iw_text("Number of MOs: ",highlight=.true.)
          call iw_text(string(f%wfn%nmoocc),sameline_nospace=.true.)
          call iw_text("Number of electrons: ",highlight=.true.)
          call iw_text(string(nint(sum(f%wfn%occ(1:f%wfn%nmoocc)))),sameline_nospace=.true.)
       end if
       call iw_text("Number of primitives: ",highlight=.true.)
       call iw_text(string(f%wfn%npri),sameline_nospace=.true.)
       call iw_text("Number of EDFs: ",highlight=.true.)
       call iw_text(string(f%wfn%nedf),sameline_nospace=.true.)

    case (type_dftb)
       call iw_text("Number of states: ",highlight=.true.)
       call iw_text(string(f%dftb%nstates),sameline_nospace=.true.)
       call iw_text("Number of spin channels: ",highlight=.true.)
       call iw_text(string(f%dftb%nspin),sameline_nospace=.true.)
       call iw_text("Number of orbitals: ",highlight=.true.)
       call iw_text(string(f%dftb%norb),sameline_nospace=.true.)
       call iw_text("Number of kpoints: ",highlight=.true.)
       call iw_text(string(f%dftb%nkpt),sameline_nospace=.true.)
       call iw_text("Real wavefunction? ",highlight=.true.)
       call iw_text(string(f%dftb%isreal),sameline_nospace=.true.)
       call iw_text("Exact calculation? ",highlight=.true.)
       call iw_text(string(f%exact),sameline_nospace=.true.)

    case (type_promol_frag)
       call iw_text("Number of atoms in fragment: ",highlight=.true.)
       call iw_text(string(f%fr%nat),sameline_nospace=.true.)

    case (type_ghost)
       call iw_text("Expression: ",highlight=.true.)
       call iw_text(string(f%expr),sameline_nospace=.true.)
    end select

    call iw_text("Use core densities: ",highlight=.true.)
    call iw_text(string(f%usecore),sameline_nospace=.true.)
    if (any(f%zpsp > 0)) then
       call iw_text("Core charges (ZPSP): ",highlight=.true.)
       str = ""
       do i = 1, maxzat0
          if (f%zpsp(i) > 0) then
             str = str // string(nameguess(i,.true.)) // "(" // string(f%zpsp(i)) // ") "
          end if
       end do
       call iw_text(str,sameline_nospace=.true.)
    end if
    call iw_text("Numerical derivatives: ",highlight=.true.)
    call iw_text(string(f%numerical),sameline_nospace=.true.)

    ! critical points
    call igNewLine()
    call iw_text("Non-equivalent critical points: ",highlight=.true.)
    call iw_text(string(f%ncp),sameline_nospace=.true.)
    call iw_text("Cell critical points: ",highlight=.true.)
    call iw_text(string(f%ncpcel),sameline_nospace=.true.)

    ! last message
    call igNewLine()
    call iw_text("[Right-click for options]",highlight=.true.)

    call igEndTooltip()

  end subroutine tree_field_tooltip_string

  !> Draw the contents of the load field window.
  module subroutine draw_load_field(w)
    use fieldseedmod, only: field_detect_format
    use param, only: ifformat_wien, ifformat_elk, ifformat_pi,&
       ifformat_cube, ifformat_bincube, ifformat_abinit, ifformat_vasp,&
       ifformat_vaspnov, ifformat_qub, ifformat_xsf, ifformat_elkgrid,&
       ifformat_siestagrid, ifformat_dftb, ifformat_pwc, ifformat_wfn,&
       ifformat_wfx, ifformat_fchk, ifformat_molden, dirsep
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use gui_main, only: nsys, sysc, sys, sys_init, g, ok_system
    use utils, only: iw_text, iw_tooltip, iw_radiobutton, iw_button,&
       iw_calcwidth
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: oksys, ok, doquit, disabled
    integer :: isys, i, j, ll, idx, iff, iaux
    type(ImVec2) :: szavail, szero
    real(c_float) :: combowidth
    logical(c_bool) :: is_selected, ldum
    character(len=:,kind=c_char), allocatable, target :: str1, str2, loadstr
    logical :: isgrid

    integer, parameter :: itoken_file1 = 1
    integer, parameter :: itoken_file2 = 2
    integer, parameter :: itoken_file3 = 3

    ! window state
    logical, save :: ttshown = .false. ! tooltip flag
    integer(c_int), save :: sourceopt = 0_c_int ! 0 = file, 1 = expression
    character(len=:,kind=c_char), allocatable, target, save :: file1 ! first (main) file
    character(len=:,kind=c_char), allocatable, target, save :: file1_fmtstr ! first file format string
    logical, save :: file1_set = .false.
    integer, save :: file1_format = 0
    character(len=:,kind=c_char), allocatable, target, save :: file2 ! first auxiliary file
    logical, save :: file2_set = .false.
    character(len=:,kind=c_char), allocatable, target, save :: file3 ! second auxiliaryy file
    logical, save :: file3_set = .false.
    integer, save :: iginterp = 3 ! 0 = nearest, 1 = trilinear, 2 = trispline, 3 = tricubic, 4 = smoothrho

    ! DFTB+ hsd name candidates
    character*15, parameter :: hsdnames(4) = (/&
       "wfc-3ob-3-1.hsd","wfc.3ob-3-1.hsd","wfc.pbc-0-3.hsd","wfc.mio-1-1.hsd"/)

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    !! make sure we get a system that exists
    ! check if the system still exists
    isys = w%isys
    oksys = ok_system(isys,sys_init)
    if (.not.oksys) then
       ! reset to the table selected
       isys = win(iwin_tree)%tree_selected
       oksys = ok_system(isys,sys_init)
    end if
    if (.not.oksys) then
       ! reset to the input console selected
       isys = win(iwin_tree)%inpcon_selected
       oksys = ok_system(isys,sys_init)
    end if
    if (.not.oksys) then
       ! check all systems and try to find one that is OK
       do i = 1, nsys
          oksys = ok_system(isys,sys_init)
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
    call iw_text("System",highlight=.true.)
    call igSameLine(0._c_float,-1._c_float)
    call igGetContentRegionAvail(szavail)
    combowidth = max(szavail%x - g%Style%ItemSpacing%x,0._c_float)
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
                w%isys = i
                isys = w%isys
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
       if (w%okfile_set) then
          w%errmsg = ""
          if (w%itoken == itoken_file1) then
             file1 = w%okfile
             file1_format = w%okfile_format
             if (file1_format == 0) &
                file1_format = -field_detect_format(file=w%okfile)
             if (abs(file1_format) == ifformat_pi) file1_format = 0 ! aiPI deactivated for now
             if (file1_format == 0) then
                w%errmsg = "Unknown field file extension"
                w%okfile = ""
             end if
             file1_set = .true.
             file2 = ""
             file2_set = .false.
             file3 = ""
             file3_set = .false.
          elseif (w%itoken == itoken_file2) then
             file2 = w%okfile
             file2_set = .true.
          elseif (w%itoken == itoken_file3) then
             file3 = w%okfile
             file3_set = .true.
          end if
       end if
       ! source file and format
       call iw_text("Source",highlight=.true.)
       if (iw_button("File",danger=(file1_format==0))) then
          iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openfieldfile,&
             idparent=w%id,itoken=itoken_file1)
       end if
       call iw_tooltip("File from where the field is read",ttshown)
       if (len(file1) > 0) then
          if (iw_button("Clear##clearfile1",sameline=.true.,danger=.true.)) then
             file1_format = 0
             file1 = ""
             file1_set = .false.
             file1_fmtstr = ""
             file2 = ""
             file2_set = .false.
             file3 = ""
             file3_set = .false.
          end if
       end if
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
          if (iw_button("File (.struct)",danger=.not.file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("WIEN2k structure file (.struct) used to interpret the clmsum-style file",ttshown)
          if (len(file2) > 0) then
             if (iw_button("Clear##clearfile2",sameline=.true.,danger=.true.)) &
                file2 = ""
          end if
          call iw_text(file2,sameline=.true.)
       case(ifformat_elk)
          if (iw_button("File (GEOMETRY.OUT)",danger=.not.file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("GEOMETRY.OUT structure file used to interpret the STATE.OUT",ttshown)
          if (len(file2) > 0) then
             if (iw_button("Clear##clearfile3",sameline=.true.,danger=.true.)) &
                file2 = ""
          end if
          call iw_text(file2,sameline=.true.)
       case(ifformat_dftb)
          if (iw_button("File (.bin)",danger=.not.file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("eigenvec.bin file for reading the DFTB+ wavefunction",ttshown)
          if (len(file2) > 0) then
             if (iw_button("Clear##clearfile4",sameline=.true.,danger=.true.)) &
                file2 = ""
          end if
          call iw_text(file2,sameline=.true.)

          if (iw_button("File (.hsd)",danger=.not.file3_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file3)
          end if
          call iw_tooltip("hsd file for reading the DFTB+ wavefunction",ttshown)
          if (len(file3) > 0) then
             if (iw_button("Clear##clearfile5",sameline=.true.,danger=.true.)) &
                file3 = ""
          end if
          call iw_text(file3,sameline=.true.)
       end select

       ! error message, if applicable
       if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    elseif (sourceopt == 1) then
       ! from expression
       w%errmsg = ""
    end if
    w%okfile_set = .false.

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
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
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
          call sys(isys)%load_field_string(loadstr,.false.,iff,w%errmsg)
          if (len_trim(w%errmsg) == 0) doquit = .true.
       elseif (sourceopt == 1) then
          ! todo
       end if
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) .or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) &
       doquit = .true.

    ! quit the window
    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
      w%errmsg = ""
      sourceopt = 0_c_int
      file1 = ""
      file1_set = .false.
      file1_format = 0
      file1_fmtstr = ""
      file2 = ""
      file2_set = .false.
      file3 = ""
      file3_set = .false.
    end subroutine init_state
    ! terminate the state for this window
    subroutine end_state()
      file1 = ""
      file1_set = .false.
      file1_format = 0
      file1_fmtstr = ""
      file2 = ""
      file2_set = .false.
      file3 = ""
      file3_set = .false.
      if (allocated(file1)) deallocate(file1)
      if (allocated(file1_fmtstr)) deallocate(file1_fmtstr)
      if (allocated(file2)) deallocate(file2)
      if (allocated(file3)) deallocate(file3)
    end subroutine end_state
  end subroutine draw_load_field

  !> Update tasks for the SCF window, before the window is created.
  module subroutine update_scfplot(w)
    use gui_main, only: sysc, sys_init, ok_system
    class(window), intent(inout), target :: w

    integer :: isys
    logical :: doquit

    isys = w%isys
    doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = (sysc(isys)%collapse >= 0)
    if (doquit) call w%end()

  end subroutine update_scfplot

  !> Draw the SCF plot window.
  module subroutine draw_scfplot(w)
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use gui_main, only: nsys, sysc, sys_init, ok_system
    use utils, only: iw_text
    use types, only: realloc
    use tools_io, only: string
    class(window), intent(inout), target :: w

    real(c_double) :: xmin, xmax
    integer(c_int) :: nticks
    real*8 :: dy
    logical :: doquit
    integer :: i, isys, num
    type(ImVec2) :: sz
    type(ImVec4) :: auto
    character(len=:,kind=c_char), allocatable, target :: str1, str2

    integer, parameter :: maxxticks = 10

    isys = w%isys
    doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = (sysc(isys)%collapse >= 0)
    if (w%firstpass) w%plotn = 0

    if (.not.doquit) then
       ! get the property and prepare x and y
       if (w%firstpass) then
          if (allocated(w%plotx)) deallocate(w%plotx)
          if (allocated(w%ploty)) deallocate(w%ploty)
          num = count(sysc(1:nsys)%collapse == isys) + 1
          allocate(w%plotx(num),w%ploty(num))
          w%plotn = 0
          do i = 1, nsys
             if (sysc(i)%collapse == isys .and. sysc(i)%seed%energy /= huge(1d0)) then
                w%plotn = w%plotn + 1
                w%plotx(w%plotn) = w%plotn
                w%ploty(w%plotn) = sysc(i)%seed%energy
             end if
          end do
          if (sysc(isys)%seed%energy /= huge(1d0)) then
             w%plotn = w%plotn + 1
             w%plotx(w%plotn) = w%plotn
             w%ploty(w%plotn) = sysc(isys)%seed%energy
          end if
          call realloc(w%plotx,w%plotn)
          call realloc(w%ploty,w%plotn)
          w%ymax = maxval(w%ploty)
          w%ymin = minval(w%ploty)
       end if

       ! make the plot
       if (w%plotn > 0) then
          str1 = "##scfiterations" // string(isys) // c_null_char
          call igGetContentRegionAvail(sz)
          if (ipBeginPlot(c_loc(str1),sz,ImPlotFlags_None)) then
             str1 = "SCF Iteration" // c_null_char
             str2 = "Energy (Ha)" // c_null_char
             call ipSetupAxes(c_loc(str1),c_loc(str2),ImPlotAxisFlags_AutoFit,ImPlotAxisFlags_AutoFit)

             str1 = "%d" // c_null_char
             nticks = size(w%plotx,1)
             xmin = w%plotx(1)
             xmax = w%plotx(nticks)
             if (nticks > maxxticks) then
                nticks = maxxticks + 1
                xmax = xmin + ceiling((xmax - xmin) / maxxticks) * maxxticks
             end if
             call ipSetupAxisTicks(ImAxis_X1,xmin,xmax,nticks)
             call ipSetupAxisFormat(ImAxis_X1,c_loc(str1))

             dy = (w%ymax - w%ymin)
             str1 = "%." // string(min(ceiling(max(abs(log10(dy)),0d0)) + 1,10)) // "f" // c_null_char
             call ipSetupAxisFormat(ImAxis_Y1,c_loc(str1))
             call ipGetPlotCurrentLimits(xmin,xmax,w%ymin,w%ymax) ! these need to be switched

             str1 = "Energy" // c_null_char
             auto%x = 0._c_float
             auto%y = 0._c_float
             auto%z = 0._c_float
             auto%w = -1._c_float
             call ipSetNextMarkerStyle(ImPlotMarker_Circle,-1._c_float,auto,-1._c_float,auto)

             call ipPlotLine(c_loc(str1),c_loc(w%plotx),c_loc(w%ploty),w%plotn,ImPlotLineFlags_None,0_c_int)
             call ipEndPlot()
          end if
       end if
    end if

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    if (doquit) call w%end()

  end subroutine draw_scfplot

  !> Update tasks for the rebond window, before the window is
  !> created.
  module subroutine update_rebond(w)
    use tools_io, only: string
    class(window), intent(inout), target :: w

    if (w%firstpass.or.w%tied_to_tree) then
       w%name = "Recalculate Bonds###rebond"  // string(w%id) // c_null_char
    else
       w%name = "Recalculate Bonds [detached]###rebond"  // string(w%id) // c_null_char
    end if

  end subroutine update_rebond

  !> Draw the contents of the rebond window.
  module subroutine draw_rebond(w)
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use gui_main, only: nsys, sysc, sys, sys_init, g, ok_system, lastchange_rebond
    use utils, only: iw_text, iw_tooltip, iw_calcwidth, iw_button, iw_calcheight
    use global, only: bondfactor_def
    use tools_io, only: string, nameguess
    use param, only: atmcov0, maxzat0, bohrtoa, newline
    class(window), intent(inout), target :: w

    logical :: doquit, ch
    integer :: i, iz, isys, natused
    type(ImVec2) :: szavail, szero, sz0
    integer(c_int) :: flags
    real(c_float) :: combowidth, rad, bf, bfmin, bfmax
    logical(c_bool) :: is_selected
    character(len=:,kind=c_char), allocatable, target :: str1, str2, str3
    integer, allocatable :: iat(:)
    logical :: atused(maxzat0)
    real*8 :: mrad

    integer, parameter :: ic_name = 0
    integer, parameter :: ic_z = 1
    integer, parameter :: ic_radius = 2

    ! window state
    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    ! first pass
    if (w%firstpass) then
       w%tied_to_tree = (w%isys == win(iwin_tree)%tree_selected)
    end if

    ! if tied to tree, update the isys
    if (w%tied_to_tree .and. (w%isys /= win(iwin_tree)%tree_selected)) &
       w%isys = win(iwin_tree)%tree_selected

    ! check if the system still exists
    if (.not.ok_system(w%isys,sys_init)) then
       ! this dialog does not make sense anymore, close it and exit
       call w%end()
       return
    end if
    isys = w%isys

    ! system combo
    call iw_text("System",highlight=.true.)
    call igSameLine(0._c_float,-1._c_float)
    call igGetContentRegionAvail(szavail)
    combowidth = max(szavail%x - g%Style%ItemSpacing%x,0._c_float)
    str1 = "##systemcombo" // c_null_char
    call igSetNextItemWidth(combowidth)
    str2 = string(isys) // ": " // trim(sysc(isys)%seed%name) // c_null_char
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nsys
          if (sysc(i)%status == sys_init) then
             is_selected = (isys == i)
             str2 = string(i) // ": " // trim(sysc(i)%seed%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) then
                w%isys = i
                isys = w%isys
                w%tied_to_tree = w%tied_to_tree .and. (w%isys == win(iwin_tree)%tree_selected)
             end if
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Recalculate the bonds in this system",ttshown)

    ! determine which atoms are to be used
    atused = .false.
    do i = 1, sys(isys)%c%nspc
       atused(sys(isys)%c%spc(i)%z) = .true.
    end do
    allocate(iat(count(atused)))
    natused = 0
    do i = 1, maxzat0
       if (atused(i)) then
          natused = natused + 1
          iat(natused) = i
       end if
    end do

    ! the radii table
    mrad = 0d0
    call iw_text("Atomic Radii",highlight=.true.)
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tableatomrcov" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(5,natused)+1,0,.false.)
    if (igBeginTable(c_loc(str1),3,flags,sz0,0._c_float)) then
       ! header setup
       str2 = "Atom" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_name)

       str2 = "Z" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_z)

       str2 = "Radius (Å)" // c_null_char
       flags = ImGuiTableColumnFlags_WidthStretch
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_radius)
       call igTableSetupScrollFreeze(0, 1) ! top row always visible

       ! draw the header
       call igTableHeadersRow()
       call igTableSetColumnWidthAutoAll(igGetCurrentTable())

       ! draw the rows
       do i = 1, natused
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
          iz = iat(i)

          ! name
          if (igTableSetColumnIndex(ic_name)) then
             call igAlignTextToFramePadding()
             call iw_text(string(nameguess(iz,.true.)))
          end if

          ! Z
          if (igTableSetColumnIndex(ic_z)) then
             call igAlignTextToFramePadding()
             call iw_text(string(iz))
          end if

          ! radius
          if (igTableSetColumnIndex(ic_radius)) then
             str2 = "##tableradius" // string(i) // c_null_char
             str3 = "%.3f" // c_null_char
             call igPushItemWidth(iw_calcwidth(5,1))
             mrad = max(mrad,sysc(isys)%atmcov(iz))
             rad = real(sysc(isys)%atmcov(iz) * bohrtoa,c_float)
             ch = igDragFloat(c_loc(str2),rad,0.01_c_float,0._c_float,2.65_c_float,c_loc(str3),&
                ImGuiSliderFlags_AlwaysClamp)
             if (ch) sysc(isys)%atmcov(iz) = rad / bohrtoa
             call igPopItemWidth()
          end if
       end do
       call igEndTable()
    end if

    ! bond factor
    call iw_text("Bond factor",highlight=.true.)
    str2 = "##bondfactor" // c_null_char
    str3 = "%.4f" // c_null_char
    call igPushItemWidth(iw_calcwidth(6,1))
    bf = real(sysc(isys)%bondfactor,c_float)
    call igSameLine(0._c_float,-1._c_float)
    bfmin = 1.0_c_float
    bfmax = 2.0_c_float

    ch = igDragFloat(c_loc(str2),bf,0.001_c_float,bfmin,bfmax,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Bond factor parameter for connectivity calculation",ttshown)
    if (ch) sysc(isys)%bondfactor = bf
    call igPopItemWidth()

    ! explanation message
    call iw_text("Atoms i and j are bonded if:"//newline//&
       "  d < (ri+rj)*(bond factor)"//newline//&
       "where d = atom distance; ri, rj = atom radii")

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(34,4,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! reset to the default bonds
    if (iw_button("Reset",danger=.true.)) then
       sysc(isys)%atmcov = atmcov0
       sysc(isys)%bondfactor = bondfactor_def
    end if
    call iw_tooltip("Reset the default bonding parameters and recalculate the system bonds",ttshown)

    ! apply the changes to all systems
    if (iw_button("Apply to All Systems",sameline=.true.)) then
       do i = 1, nsys
          if (sysc(i)%status >= sys_init) then
             sysc(i)%atmcov = sysc(isys)%atmcov
             sysc(i)%bondfactor = sysc(isys)%bondfactor
             call sys(i)%c%find_asterisms(sys(i)%c%nstar,sysc(i)%atmcov,sysc(i)%bondfactor)
             call sys(i)%c%fill_molecular_fragments()
             call sys(i)%c%calculate_molecular_equivalence()
             call sys(i)%c%calculate_periodicity()
             call sysc(i)%set_timelastchange(lastchange_rebond)
          end if
       end do
    end if

    ! apply the changes
    if (iw_button("Apply",sameline=.true.)) then
       ! find the atomic connectivity and the molecular fragments
       call sys(isys)%c%find_asterisms(sys(isys)%c%nstar,sysc(isys)%atmcov,sysc(isys)%bondfactor)
       call sys(isys)%c%fill_molecular_fragments()
       call sys(isys)%c%calculate_molecular_equivalence()
       call sys(isys)%c%calculate_periodicity()
       call sysc(isys)%set_timelastchange(lastchange_rebond)
    end if
    call iw_tooltip("Recalculate the system bonds with the selected parameters",ttshown)

    ! close button
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.
    doquit = doquit .or. iw_button("Close",sameline=.true.)

    ! quit the window
    if (doquit) &
       call w%end()

  end subroutine draw_rebond

  !> Draw the tree plot window
  module subroutine draw_treeplot(w)
    use interfaces_glfw, only: glfwGetTime
    use gui_main, only: sysc, sys_empty, sys_init, sys
    use windows, only: win, iwin_tree
    use utils, only: iw_calcwidth, iw_combo_simple, iw_tooltip, iw_text
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use tools_io, only: string
    use types, only: realloc
    use param, only: bohrtoa
    class(window), intent(inout), target :: w

    type(ImVec2) :: sz
    logical :: ok
    integer :: i, j, nshown
    real*8 :: valx, valy
    character(len=:,kind=c_char), allocatable, target :: str1, str2
    type(ImVec4) :: auto
    logical(c_bool) :: ch, forceupdate

    integer, save :: ic_plotx = 0, ic_ploty = 0
    logical, save :: ttshown = .false. ! tooltip flag

    integer, parameter :: ictrans(0:14) = (/ic_id,ic_e,ic_v,ic_vmol,ic_nneq,ic_ncel,&
       ic_nmol,ic_a,ic_b,ic_c,ic_alpha,ic_beta,ic_gamma,ic_emol,ic_p/)

    ! initialize
    forceupdate = .false.
    if (w%firstpass) then
       ic_plotx = 2
       ic_ploty = 1
       forceupdate = .true.
    end if

    ! update the plot based on time signals between dependent windows
    !! if the tree info has been updated, also update the plot
    if (w%timelast_plot_update < win(iwin_tree)%timelast_tree_update) forceupdate = .true.

    ! x-choose and y-choose
    str1 = "ID"//c_null_char//"Energy (Ha)"//c_null_char//"Volume (Å³)"//c_null_char//&
       "Volume/Z (Å³)"//c_null_char//&
       "Number of symmetry-unique atoms"//c_null_char//"Number of cell atoms"//c_null_char//&
       "Number of molecules"//c_null_char//"a (Å)"//c_null_char//&
       "b (Å)"//c_null_char//"c (Å)"//c_null_char//"α (°)"//c_null_char//"β (°)"//c_null_char//&
       "γ (°)"//c_null_char//"Energy/Z (Ha)"//c_null_char//"Pressure (GPa)"//c_null_char//c_null_char
    call iw_text("x: ")
    call iw_combo_simple("##xselect",str1,ic_plotx,changed=ch,sameline=.true.)
    call iw_tooltip("Property to represent on the x-axis",ttshown)
    forceupdate = forceupdate .or. ch

    call iw_text("y: ")
    call iw_combo_simple("##yselect",str1,ic_ploty,changed=ch,sameline=.true.)
    call iw_tooltip("Property to represent on the y-axis",ttshown)
    forceupdate = forceupdate .or. ch

    ! update the plot data if necessary
    if (forceupdate) then
       w%plotn = 0
       if (allocated(w%plotx)) deallocate(w%plotx)
       if (allocated(w%ploty)) deallocate(w%ploty)
       nshown = size(win(iwin_tree)%iord,1)
       allocate(w%plotx(nshown),w%ploty(nshown))
       do j = 1, nshown
          i = win(iwin_tree)%iord(j)
          if (sysc(i)%status == sys_empty .or. sysc(i)%hidden) cycle

          ok = getvalue(valx,ictrans(ic_plotx),i)
          ok = ok .and. getvalue(valy,ictrans(ic_ploty),i)
          if (ok) then
             w%plotn = w%plotn + 1
             w%plotx(w%plotn) = valx
             w%ploty(w%plotn) = valy
          end if
       end do
       if (w%plotn > 0) then
          call realloc(w%plotx,w%plotn)
          call realloc(w%ploty,w%plotn)
       end if
       w%timelast_plot_update = glfwGetTime()
    end if

    ! make the plot
    if (w%plotn > 0) then
       str1 = "##treeplot" // c_null_char
       call igGetContentRegionAvail(sz)
       if (ipBeginPlot(c_loc(str1),sz,ImPlotFlags_None)) then
          call getname(str1,ictrans(ic_plotx))
          call getname(str2,ictrans(ic_ploty))
          call ipSetupAxes(c_loc(str1),c_loc(str2),ImPlotAxisFlags_AutoFit,ImPlotAxisFlags_AutoFit)

          auto%x = 0._c_float
          auto%y = 0._c_float
          auto%z = 0._c_float
          auto%w = -1._c_float
          call ipSetNextMarkerStyle(ImPlotMarker_Circle,-1._c_float,auto,-1._c_float,auto)

          call ipPlotLine(c_loc(str2),c_loc(w%plotx),c_loc(w%ploty),w%plotn,ImPlotLineFlags_None,0_c_int)
          call ipEndPlot()
       end if
    end if

    ! close
    if ((w%focused() .and. (is_bind_event(BIND_OK_FOCUSED_DIALOG) .or. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG) .or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)))) &
       call w%end()

  contains
    ! Get value ic for system i and return it in val. The function returns ok if it
    ! is a valid value.
    function getvalue(val,ic,i)
      real*8, intent(out) :: val
      integer, intent(in) :: ic, i
      logical :: getvalue

      val = 0d0
      getvalue = .false.
      if (ic == ic_id) then
         val = real(i,8)
      elseif (ic == ic_e) then
         if (sysc(i)%seed%energy == huge(1d0)) return
         val = sysc(i)%seed%energy
      elseif (ic == ic_v) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%omega*bohrtoa**3
      elseif (ic == ic_vmol) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%omega*bohrtoa**3/sys(i)%c%nmol
      elseif (ic == ic_nneq) then
         if (sysc(i)%status /= sys_init) return
         val = sys(i)%c%nneq
      elseif (ic == ic_ncel) then
         if (sysc(i)%status /= sys_init) return
         val = sys(i)%c%ncel
      elseif (ic == ic_nmol) then
         if (sysc(i)%status /= sys_init) return
         val = sys(i)%c%nmol
      elseif (ic == ic_a) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%aa(1)*bohrtoa
      elseif (ic == ic_b) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%aa(2)*bohrtoa
      elseif (ic == ic_c) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%aa(3)*bohrtoa
      elseif (ic == ic_alpha) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%bb(1)
      elseif (ic == ic_beta) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%bb(2)
      elseif (ic == ic_gamma) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%bb(3)
      elseif (ic == ic_emol) then
         if (sysc(i)%status /= sys_init .or. sysc(i)%seed%energy == huge(1d0)) return
         val = sysc(i)%seed%energy/sys(i)%c%nmol
      elseif (ic == ic_p) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule .or.&
             sysc(i)%seed%pressure == huge(1d0)) return
         val = sysc(i)%seed%pressure
      end if
      getvalue = .true.

    end function getvalue

    ! Get value ic for system i and return it in val. The function returns ok if it
    ! is a valid value.
    subroutine getname(str,ic)
      character(len=:,kind=c_char), allocatable, target :: str
      integer, intent(in) :: ic

      if (ic == ic_id) then
         str = "ID" // c_null_char
      elseif (ic == ic_e) then
         str = "Energy (Ha)" // c_null_char
      elseif (ic == ic_v) then
         str = "Volume (Å³)" // c_null_char
      elseif (ic == ic_vmol) then
         str = "Volume/Z (Å³)" // c_null_char
      elseif (ic == ic_nneq) then
         str = "Number of symmetry-unique atoms" // c_null_char
      elseif (ic == ic_ncel) then
         str = "Number of cell atoms" // c_null_char
      elseif (ic == ic_nmol) then
         str = "Number of molecules" // c_null_char
      elseif (ic == ic_a) then
         str = "a (Å)" // c_null_char
      elseif (ic == ic_b) then
         str = "b (Å)" // c_null_char
      elseif (ic == ic_c) then
         str = "c (Å)" // c_null_char
      elseif (ic == ic_alpha) then
         str = "α (°)" // c_null_char
      elseif (ic == ic_beta) then
         str = "β (°)" // c_null_char
      elseif (ic == ic_gamma) then
         str = "γ (°)" // c_null_char
      elseif (ic == ic_emol) then
         str = "Energy/Z (Ha)" // c_null_char
      elseif (ic == ic_p) then
         str = "Pressure (GPa)" // c_null_char
      end if

    end subroutine getname

  end subroutine draw_treeplot

  !> Update tasks for the geometry window, before the window is
  !> created.
  module subroutine update_geometry(w)
    use tools_io, only: string
    class(window), intent(inout), target :: w

    if (w%firstpass.or.w%tied_to_tree) then
       w%name = "View/Edit Geometry###view_edit_geometry"  // string(w%id) // c_null_char
    else
       w%name = "View/Edit Geometry [detached]###view_edit_geometry"  // string(w%id) // c_null_char
    end if

  end subroutine update_geometry

  !> Draw the geometry window.
  module subroutine draw_geometry(w)
    use scenes, only: reptype_atoms
    use windows, only: iwin_view, iwin_tree
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS, BIND_EDITGEOM_REMOVE
    use gui_main, only: nsys, sysc, sys, sys_init, g, ok_system, ColorHighlightScene,&
       ColorHighlightSelectScene
    use utils, only: iw_text, iw_tooltip, iw_calcwidth, iw_button, iw_calcheight, iw_calcwidth,&
       iw_combo_simple, iw_highlight_selectable, iw_coloredit
    use global, only: dunit0, iunit_ang
    use tools_io, only: string, nameguess, ioj_right
    class(window), intent(inout), target :: w

    logical :: domol, dowyc, doidx, havesel, removehighlight
    logical :: doquit, clicked
    integer :: ihighlight, iclicked, nhigh
    logical(c_bool) :: is_selected, redo_highlights
    integer(c_int) :: atompreflags, flags, ntype, ncol, ndigit, ndigitm, ndigitidx
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2, suffix
    character(len=:), allocatable :: name
    integer, allocatable :: ihigh(:)
    real(c_float), allocatable :: irgba(:,:)
    type(ImVec2) :: szavail, szero, sz0
    real(c_float) :: combowidth, rgb(3)
    integer :: i, j, isys, icol, ispc, iz, iview
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    logical :: havergb, ldum, ok
    real*8 :: x0(3)
    type(ImVec4) :: col4
    integer(c_int) :: color

    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    ihighlight = 0
    iclicked = 0
    doquit = .false.
    szero%x = 0
    szero%y = 0
    redo_highlights = .false.
    removehighlight = .false.

    ! first pass
    if (w%firstpass) then
       w%geometry_atomtype = 1
       w%tied_to_tree = (w%isys == win(iwin_tree)%tree_selected)
       if (allocated(w%geometry_selected)) deallocate(w%geometry_selected)
       if (allocated(w%geometry_rgba)) deallocate(w%geometry_rgba)
       w%geometry_select_rgba = ColorHighlightSelectScene
    end if

    ! if tied to tree, update the isys
    if (w%tied_to_tree .and. (w%isys /= win(iwin_tree)%tree_selected)) &
       call change_system(win(iwin_tree)%tree_selected)

    ! check if the system still exists
    if (.not.ok_system(w%isys,sys_init)) then
       ! this dialog does not make sense anymore, close it and exit
       call w%end()
       return
    end if
    isys = w%isys

    ! system combo
    atompreflags = ImGuiTabItemFlags_None
    call iw_text("System",highlight=.true.)
    call igSameLine(0._c_float,-1._c_float)
    call igGetContentRegionAvail(szavail)
    combowidth = max(szavail%x - g%Style%ItemSpacing%x,0._c_float)
    str1 = "##systemcombo" // c_null_char
    call igSetNextItemWidth(combowidth)
    str2 = string(isys) // ": " // trim(sysc(isys)%seed%name) // c_null_char
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nsys
          if (sysc(i)%status == sys_init) then
             is_selected = (isys == i)
             str2 = string(i) // ": " // trim(sysc(i)%seed%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) then
                call change_system(i)
                isys = w%isys
                atompreflags = ImGuiTabItemFlags_SetSelected
             end if
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Recalculate the bonds in this system",ttshown)

    ! show the tabs
    str1 = "##drawgeometry_tabbar" // c_null_char
    flags = ImGuiTabBarFlags_Reorderable
    call igBeginGroup()
    if (igBeginTabBar(c_loc(str1),flags)) then
       !! atoms tab !!
       str2 = "Atoms##drawgeometry_atomstab" // c_null_char
       flags = atompreflags
       if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
          ! group atom types
          if (.not.sys(isys)%c%ismolecule) then
             call iw_combo_simple("Atom types##atomtypeselectgeom","Species"//c_null_char//&
                "Symmetry unique" //c_null_char//"Cell"//c_null_char//c_null_char,&
                w%geometry_atomtype)
          else
             call iw_combo_simple("Atom types##atomtypeselectgeom","Species"//c_null_char//"Atoms"//c_null_char//&
                c_null_char,w%geometry_atomtype)
          end if
          call iw_tooltip("Group atoms by these categories",ttshown)
          if (w%geometry_atomtype == 0) then
             ntype = sys(isys)%c%nspc
          elseif (w%geometry_atomtype == 1) then
             ntype = sys(isys)%c%nneq
          elseif (w%geometry_atomtype == 2) then
             ntype = sys(isys)%c%ncel
          end if

          ! reallocate if ntype has changed and redo highlights
          if (allocated(w%geometry_selected)) then
             if (size(w%geometry_selected,1) /= ntype) then
                deallocate(w%geometry_selected)
                if (allocated(w%geometry_rgba)) deallocate(w%geometry_rgba)
             end if
          end if
          if (.not.allocated(w%geometry_selected)) then
             allocate(w%geometry_selected(ntype))
             w%geometry_selected = .false.
             redo_highlights = .true.
          end if
          if (.not.allocated(w%geometry_rgba)) then
             allocate(w%geometry_rgba(4,ntype))
             w%geometry_rgba = 0._c_float
          end if

          ! whether to do the molecule column, wyckoff, nneq index
          domol = (w%geometry_atomtype == 2 .or. (w%geometry_atomtype == 1 .and. sys(isys)%c%ismolecule))
          dowyc = (w%geometry_atomtype == 1 .and..not.sys(isys)%c%ismolecule)
          doidx = (w%geometry_atomtype == 2 .and..not.sys(isys)%c%ismolecule)

          ! number of columns
          ncol = 3
          if (domol) ncol = ncol + 1 ! mol
          if (dowyc) ncol = ncol + 1 ! wyckoff/multiplicity
          if (doidx) ncol = ncol + 1 ! nneq idx
          if (w%geometry_atomtype > 0) ncol = ncol + 1 ! coordinates

          ! atom style table, for atoms
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_Resizable)
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          str1="##tableatomstyles" // c_null_char
          sz0%x = 0
          sz0%y = iw_calcheight(min(10,ntype)+1,0,.false.)
          if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
             icol = -1

             ! header setup
             icol = icol + 1
             str2 = "Id" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)

             icol = icol + 1
             str2 = "Atom" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)

             icol = icol + 1
             str2 = "Z " // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)

             if (domol) then
                icol = icol + 1
                str2 = "mol" // c_null_char
                flags = ImGuiTableColumnFlags_None
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
             end if

             if (dowyc) then
                icol = icol + 1
                if (sys(isys)%c%havesym > 0 .and. sys(isys)%c%spgavail) then
                   str2 = "Wyc" // c_null_char
                else
                   str2 = "Mul" // c_null_char
                end if
                flags = ImGuiTableColumnFlags_None
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
             end if

             if (doidx) then
                icol = icol + 1
                str2 = "idx" // c_null_char
                flags = ImGuiTableColumnFlags_None
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
             end if

             if (w%geometry_atomtype > 0) then
                icol = icol + 1
                if (sys(isys)%c%ismolecule) then
                   str2 = "Coordinates (Å)" // c_null_char
                else
                   str2 = "Coordinates (fractional)" // c_null_char
                end if
                flags = ImGuiTableColumnFlags_WidthStretch
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
             end if
             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! start the clipper
             clipper = ImGuiListClipper_ImGuiListClipper()
             call ImGuiListClipper_Begin(clipper,ntype,-1._c_float)

             ! calculate the number of digits for output
             ndigit = ceiling(log10(ntype+0.1d0))
             ndigitm = 0
             ndigitidx = 0
             if (domol) ndigitm = ceiling(log10(sys(isys)%c%nmol+0.1d0))
             if (doidx) ndigitidx = ceiling(log10(sys(isys)%c%nneq+0.1d0))

             ! get the current view, if available
             iview = 0
             if (isys == win(iwin_view)%view_selected) then
                iview = iwin_view
             else
                do j = 1, nwin
                   if (.not.win(j)%isinit) cycle
                   if (win(j)%type /= wintype_view.or..not.associated(win(j)%sc)) cycle
                   if (isys == win(j)%view_selected) then
                      iview = j
                      exit
                   end if
                end do
             end if

             ! draw the rows
             do while(ImGuiListClipper_Step(clipper))
                call c_f_pointer(clipper,clipper_f)
                do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                   suffix = "_" // string(i)
                   icol = -1

                   ! start the table and identify the species and Z
                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   if (w%geometry_atomtype == 0) then ! species
                      ispc = i
                      name = string(sys(isys)%c%spc(ispc)%name,2)
                   elseif (w%geometry_atomtype == 1) then ! nneq
                      ispc = sys(isys)%c%at(i)%is
                      name = trim(sys(isys)%c%at(i)%name)
                   elseif (w%geometry_atomtype == 2) then ! ncel
                      ispc = sys(isys)%c%atcel(i)%is
                      name = trim(sys(isys)%c%at(sys(isys)%c%atcel(i)%idx)%name)
                   end if
                   iz = sys(isys)%c%spc(ispc)%z

                   ! get the color from the first active atoms representation in the main view
                   havergb = .false.
                   if (iview > 0) then
                      do j = 1, win(iview)%sc%nrep
                         if (win(iview)%sc%rep(j)%type == reptype_atoms.and.win(iview)%sc%rep(j)%isinit.and.&
                            win(iview)%sc%rep(j)%shown) then
                            if (win(iview)%sc%rep(j)%atom_style%type == 0) then ! color by species
                               rgb = win(iview)%sc%rep(j)%atom_style%rgb(:,ispc)
                               havergb = .true.
                            elseif (win(iview)%sc%rep(j)%atom_style%type == w%geometry_atomtype) then ! color by nneq or ncel
                               rgb = win(iview)%sc%rep(j)%atom_style%rgb(:,i)
                               havergb = .true.
                            elseif (win(iview)%sc%rep(j)%atom_style%type == 1 .and. w%geometry_atomtype == 2) then ! color by nneq, select by ncel
                               rgb = win(iview)%sc%rep(j)%atom_style%rgb(:,sys(isys)%c%atcel(i)%idx)
                               havergb = .true.
                            end if
                            if (havergb) exit
                         end if
                      end do
                   end if

                   ! background color for the table row
                   if (w%geometry_selected(i)) then
                      col4 = ImVec4(w%geometry_rgba(1,i),w%geometry_rgba(2,i),&
                         w%geometry_rgba(3,i),w%geometry_rgba(4,i))
                      color = igGetColorU32_Vec4(col4)
                      call igTableSetBgColor(ImGuiTableBgTarget_RowBg0, color, icol)
                   end if

                   ! id
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      call igAlignTextToFramePadding()
                      call iw_text(string(i,ndigit))

                      ! the highlight selectable: hover and click
                      clicked = .false.
                      ok = iw_highlight_selectable("##selectablemoltable" // suffix,&
                         selected=w%geometry_selected(i),clicked=clicked)
                      if (ok) ihighlight = i
                      if (clicked) iclicked = i
                   end if

                   ! name
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      if (havergb) then
                         ldum = iw_coloredit("##tablecolorg" // suffix,rgb=rgb,nointeraction=.true.)
                         call iw_text(name,sameline=.true.)
                      else
                         call iw_text(name)
                      end if
                   end if

                   ! Z
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) call iw_text(string(iz))

                   ! molecule
                   if (domol) then
                      icol = icol + 1
                      ! i is a complete list index in this case
                      if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%idatcelmol(1,i),ndigitm))
                   end if

                   ! multiplicity
                   if (dowyc) then
                      icol = icol + 1
                      ! i is an nneq index in this case
                      if (igTableSetColumnIndex(icol)) then
                         if (sys(isys)%c%havesym > 0 .and. sys(isys)%c%spgavail) then
                            call iw_text(string(sys(isys)%c%at(i)%mult) // sys(isys)%c%at(i)%wyc)
                         else
                            call iw_text(string(sys(isys)%c%at(i)%mult,2))
                         end if
                      end if
                   end if

                   ! multiplicity
                   if (doidx) then
                      icol = icol + 1
                      ! i is cell index in this case
                      if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%atcel(i)%idx,ndigitidx))
                   end if

                   ! coordinates
                   if (w%geometry_atomtype > 0) then
                      icol = icol + 1
                      if (igTableSetColumnIndex(icol)) then
                         s = ""
                         if (w%geometry_atomtype > 0) then
                            if (sys(isys)%c%ismolecule) then
                               x0 = (sys(isys)%c%atcel(i)%r+sys(isys)%c%molx0) * dunit0(iunit_ang)
                            elseif (w%geometry_atomtype == 1) then
                               x0 = sys(isys)%c%at(i)%x
                            else
                               x0 = sys(isys)%c%atcel(i)%x
                            endif
                            s = string(x0(1),'f',8,4,ioj_right) //" "// string(x0(2),'f',8,4,ioj_right) //" "//&
                               string(x0(3),'f',8,4,ioj_right)
                         end if
                         call iw_text(s)
                      end if
                   end if
                end do ! clipper indices
             end do ! clipper step

             ! end the clipper and the table
             call ImGuiListClipper_End(clipper)
             call igEndTable()
          end if

          !! highlight/selection row
          ! highlight color
          call igAlignTextToFramePadding()
          call iw_text("Selection",highlight=.true.)
          call igSameLine(0._c_float,-1._c_float)
          ldum = iw_coloredit("##drawgeometryhighlightcolor",rgba=w%geometry_select_rgba)
          call iw_tooltip("Color used for highlighting atoms")

          ! Highlight buttons: all, none, toggle
          if (iw_button("All##highlightall",sameline=.true.)) then
             w%geometry_selected = .true.
             do i = 1, size(w%geometry_selected,1)
                w%geometry_rgba(:,i) = w%geometry_select_rgba
             end do
             redo_highlights = .true.
          end if
          call iw_tooltip("Select all atoms in the system",ttshown)
          if (iw_button("None##highlightnone",sameline=.true.)) then
             w%geometry_selected = .false.
             redo_highlights = .true.
          end if
          call iw_tooltip("Deselect all atoms",ttshown)
          if (iw_button("Toggle##highlighttoggle",sameline=.true.)) then
             w%geometry_selected = .not.w%geometry_selected
             do i = 1, size(w%geometry_selected,1)
                if (w%geometry_selected(i)) &
                   w%geometry_rgba(:,i) = w%geometry_select_rgba
             end do
             redo_highlights = .true.
          end if
          call iw_tooltip("Toggle atomic selection",ttshown)

          !! edit row
          ! highlight color
          call igAlignTextToFramePadding()
          call iw_text("Edit",highlight=.true.)

          ! Remove button
          havesel = any(w%geometry_selected)
          if (iw_button("Remove##removeselection",sameline=.true.,disabled=.not.havesel)) &
             removehighlight = .true.
          call iw_tooltip("Remove selected atoms (" // trim(get_bind_keyname(BIND_EDITGEOM_REMOVE)) // ")",ttshown)

          call igEndTabItem()
       end if

       !! cell tab !!
       if (.not.sys(isys)%c%ismolecule) then
          str2 = "Cell##drawgeometry_celltab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
             call iw_text("blah")
             call igEndTabItem()
          end if
       end if

       !! molecules tab !!
       str2 = "Molecules##drawgeometry_molstab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
          call iw_text("blah")
          call igEndTabItem()
       end if

       !! bonds tab !!
       str2 = "Bonds##drawgeometry_bondstab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
          call iw_text("blah")
          call igEndTabItem()
       end if

       call igEndTabBar()
    end if
    call igEndGroup()

    ! hover highlight
    if (ihighlight > 0) then
       call sysc(isys)%highlight_atoms(.true.,(/ihighlight/),w%geometry_atomtype,&
          reshape(ColorHighlightScene,(/4,1/)))
    end if

    ! process clicked
    if (iclicked > 0) then
       w%geometry_selected(iclicked) = .not.w%geometry_selected(iclicked)
       w%geometry_rgba(:,iclicked) = w%geometry_select_rgba
       redo_highlights = .true.
    end if

    ! redo highlights
    if (redo_highlights.and.allocated(w%geometry_selected).and.allocated(w%geometry_rgba)) then
       call sysc(isys)%highlight_clear(.false.)
       allocate(ihigh(count(w%geometry_selected)),irgba(4,count(w%geometry_selected)))
       nhigh = 0
       do i = 1, size(w%geometry_selected,1)
          if (w%geometry_selected(i)) then
             nhigh = nhigh + 1
             ihigh(nhigh) = i
             irgba(:,nhigh) = w%geometry_rgba(:,i)
          end if
       end do
       call sysc(isys)%highlight_atoms(.false.,ihigh,w%geometry_atomtype,irgba)
       deallocate(ihigh,irgba)
    end if

    ! remove highlighted atoms
    removehighlight = removehighlight .or. (w%focused() .and. is_bind_event(BIND_EDITGEOM_REMOVE))
    if (removehighlight) &
       call sysc(isys)%remove_highlighted_atoms()

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(5,1,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! close button
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) &
       doquit = .true.
    doquit = doquit .or. iw_button("Close")

    ! quit the window
    if (doquit) &
       call w%end()

  contains
    subroutine change_system(i)
      integer, intent(in) :: i

      ! do nothing if we are already in the same system
      if (w%isys == i) return

      ! clear the highlights for the current system
      call sysc(w%isys)%highlight_clear(.false.)

      ! remove the selecion and highlights
      if (allocated(w%geometry_selected)) deallocate(w%geometry_selected)
      if (allocated(w%geometry_rgba)) deallocate(w%geometry_rgba)

      ! change the system
      w%isys = i
      w%tied_to_tree = w%tied_to_tree .and. (w%isys == win(iwin_tree)%tree_selected)

    end subroutine change_system

  end subroutine draw_geometry

end submodule tree
