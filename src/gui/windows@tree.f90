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

  !xx! private procedures
  ! function tree_system_tooltip_string(i)
  ! function tree_field_tooltip_string(si,fj)

contains

  !> Draw the contents of a tree window
  module subroutine draw_tree(w)
    use keybindings, only: is_bind_event, BIND_TREE_REMOVE_SYSTEM_FIELD, BIND_TREE_MOVE_UP,&
       BIND_TREE_MOVE_DOWN
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button,&
       iw_text, iw_setposx_fromend, iw_calcwidth, iw_calcheight
    use gui_main, only: nsys, sys, sysc, sys_empty, sys_init,&
       sys_loaded_not_init, sys_initializing, ColorTableCellBg_Mol,&
       ColorTableCellBg_MolClus, ColorTableCellBg_MolCrys, ColorTableCellBg_Crys3d,&
       ColorTableCellBg_Crys2d, ColorTableCellBg_Crys1d, launch_initialization_thread,&
       kill_initialization_thread, system_shorten_names, remove_system, tooltip_delay,&
       ColorDangerButton, ColorFieldSelected, g, tree_select_updates_inpcon,&
       tree_select_updates_view, io
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
    integer :: i, j, k, nshown, newsel, jsel, ll, id, iref, inext, iprev
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
    call iw_tooltip("Write the current table to the output console in csv-style (for copying)",ttshown)

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
          ! if we removed the system for the input console or the view, update
          if (w%forceremove(k) == win(iwin_console_input)%inpcon_selected) &
             win(iwin_console_input)%inpcon_selected = w%table_selected
          if (w%forceremove(k) == win(iwin_view)%view_selected) then
             win(iwin_view)%view_selected = w%table_selected
             win(iwin_view)%forcerender = .true.
          end if
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

    ! next and previous systems
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
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
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
                      if (tree_select_updates_view) then
                         if (win(iwin_view)%view_selected /= i) &
                            win(iwin_view)%forcerender = .true.
                         win(iwin_view)%view_selected = i
                      end if
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

    ! process the key bindings
    if (is_bind_event(BIND_TREE_MOVE_UP)) then
       if (iprev > 0) then
          w%forceselect = iprev
          call igSetWindowFocus_Str(c_loc(w%name))
       end if
    elseif (is_bind_event(BIND_TREE_MOVE_DOWN)) then
       if (inext > 0) then
          w%forceselect = inext
          call igSetWindowFocus_Str(c_loc(w%name))
       end if
    end if

    ! if exporting, read the export command
    if (export) &
       ldum = win(iwin_console_input)%read_output_ci(.true.,"[Table export]")

    ! clean up
    ! call ImGuiTextFilter_destroy(cfilter)

  contains

    subroutine write_maybe_selectable(isys,bclose,bexpand)
      use gui_main, only: are_threads_running
      use utils, only: iw_text
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
      ok = igSelectable_Bool(c_loc(strl),selected,flags,szero)
      ok = ok .or. (w%forceselect == isys)
      if (ok) then
         w%table_selected = isys
         if (w%forceselect > 0) then
            w%forceselect = 0
            call igSetKeyboardFocusHere(0)
         end if
         if (tree_select_updates_inpcon) &
            win(iwin_console_input)%inpcon_selected = isys
         if (tree_select_updates_view) then
            if (win(iwin_view)%view_selected /= isys) &
               win(iwin_view)%forcerender = .true.
            win(iwin_view)%view_selected = isys
         end if
         if (igIsMouseDoubleClicked(ImGuiPopupFlags_MouseButtonLeft)) &
            sysc(isys)%showfields = .true.
      end if
      call igSameLine(0._c_float,-1._c_float)
      call igSetCursorPosX(pos)

      ! update the iprev and inext
      if (isys < w%table_selected) then
         iprev = isys
      elseif (isys > w%table_selected .and. inext == 0) then
         inext = isys
      end if

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
         if (igMenuItem_Bool(c_loc(strpop),c_null_ptr,.false._c_bool,enabled)) then
            win(iwin_console_input)%inpcon_selected = isys
            if (win(iwin_view)%view_selected /= isys) &
               win(iwin_view)%forcerender = .true.
            win(iwin_view)%view_selected = isys
            w%table_selected = isys
         end if
         call iw_tooltip("Set this system as current",ttshown)

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
         if (igIsMouseHoveringRect(g%LastItemData%NavRect%min,g%LastItemData%NavRect%max,.false._c_bool)) then
            if (bclose) then
               strl = "Close this system" // c_null_char
            elseif (bexpand) then
               strl = "Expand this system" // c_null_char
            else
               strl = tree_system_tooltip_string(isys)
            end if
            call igSetTooltip(c_loc(strl))
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
      if (w%table_selected >= 1 .and. w%table_selected <= nsys) then
         if (sysc(w%table_selected)%collapse == i) w%table_selected = i
      end if
      if (win(iwin_console_input)%inpcon_selected >= 1 .and. win(iwin_console_input)%inpcon_selected <= nsys) then
         if (sysc(win(iwin_console_input)%inpcon_selected)%collapse == i) &
            win(iwin_console_input)%inpcon_selected = i
      end if
      if (win(iwin_view)%view_selected >= 1 .and. win(iwin_view)%view_selected <= nsys) then
         if (sysc(win(iwin_view)%view_selected)%collapse == i) then
            if (win(iwin_view)%view_selected /= i) &
               win(iwin_view)%forcerender = .true.
            win(iwin_view)%view_selected = i
         end if
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

  !> Draw the contents of the load field window.
  module subroutine draw_load_field(w)
    use fieldseedmod, only: field_detect_format
    use param, only: ifformat_wien, ifformat_elk, ifformat_pi,&
       ifformat_cube, ifformat_bincube, ifformat_abinit, ifformat_vasp,&
       ifformat_vaspnov, ifformat_qub, ifformat_xsf, ifformat_elkgrid,&
       ifformat_siestagrid, ifformat_dftb, ifformat_pwc, ifformat_wfn,&
       ifformat_wfx, ifformat_fchk, ifformat_molden, dirsep
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG
    use gui_main, only: nsys, sysc, sys, sys_init, g
    use utils, only: iw_text, iw_tooltip, iw_radiobutton, iw_button,&
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

       ! disabled buttons?
       disabled = (idopenfile1 /= 0) .or. (idopenfile2 /= 0) .or. (idopenfile3 /= 0)

       ! source file and format
       call iw_text("Source",highlight=.true.)
       if (iw_button("File",danger=(file1_format==0),disabled=disabled)) &
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
          if (iw_button("File (.struct)",danger=.not.file2_set,disabled=disabled)) &
             idopenfile2 = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal)
          call iw_tooltip("WIEN2k structure file (.struct) used to interpret the clmsum-style file",ttshown)
          call iw_text(file2,sameline=.true.)
       case(ifformat_elk)
          if (iw_button("File (GEOMETRY.OUT)",danger=.not.file2_set,disabled=disabled)) &
             idopenfile2 = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal)
          call iw_tooltip("GEOMETRY.OUT structure file used to interpret the STATE.OUT",ttshown)
          call iw_text(file2,sameline=.true.)
       case(ifformat_dftb)
          if (iw_button("File (.bin)",danger=.not.file2_set,disabled=disabled)) &
             idopenfile2 = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal)
          call iw_tooltip("eigenvec.bin file for reading the DFTB+ wavefunction",ttshown)
          call iw_text(file2,sameline=.true.)

          if (iw_button("File (.hsd)",danger=.not.file3_set,disabled=disabled)) &
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
       if (idopenfile1 /= 0) win(idopenfile1)%forcequitdialog = .true.
       if (idopenfile2 /= 0) win(idopenfile2)%forcequitdialog = .true.
       if (idopenfile3 /= 0) win(idopenfile3)%forcequitdialog = .true.
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
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG
    use gui_main, only: nsys, sysc, sys_init
    use utils, only: iw_text
    use tools_io, only: string
    class(window), intent(inout), target :: w

    real(c_double), allocatable, target, save :: x(:), y(:)
    real(c_double) :: xmin, xmax
    real(c_double), save :: ymin, ymax
    real*8 :: dy
    logical :: doquit
    integer :: i, isys, num, n
    type(ImVec2) :: sz
    type(ImVec4) :: auto
    character(len=:,kind=c_char), allocatable, target :: str1, str2

    isys = w%scfplot_isys
    doquit = (isys < 1 .or. isys > nsys)
    if (.not.doquit) doquit = (sysc(isys)%status /= sys_init)
    if (.not.doquit) doquit = (sysc(isys)%collapse >= 0)

    if (.not.doquit) then
       ! get the property and prepare x and y
       if (.not.allocated(x) .or..not.allocated(y)) then
          if (allocated(x)) deallocate(x)
          if (allocated(y)) deallocate(y)
          num = count(sysc(1:nsys)%collapse == isys) + 1
          allocate(x(num),y(num))
          n = 0
          do i = 1, nsys
             if (sysc(i)%collapse == isys) then
                n = n + 1
                x(n) = n
                y(n) = sysc(i)%seed%energy
             end if
          end do
          x(num) = num
          y(num) = sysc(isys)%seed%energy
          ymax = maxval(y)
          ymin = minval(y)
       end if

       ! make the plot
       str1 = "##scfiterations" // string(isys) // c_null_char
       call igGetContentRegionAvail(sz)
       if (ipBeginPlot(c_loc(str1),sz,ImPlotFlags_None)) then
          str1 = "SCF Iteration" // c_null_char
          str2 = "Energy (Ha)" // c_null_char
          call ipSetupAxes(c_loc(str1),c_loc(str2),ImPlotAxisFlags_None,ImPlotAxisFlags_None)

          str1 = "%.0f" // c_null_char
          call ipSetupAxisTicks(ImAxis_X1,x(1),x(size(x,1)),size(x,1))
          call ipSetupAxisFormat(ImAxis_X1,c_loc(str1))

          dy = (ymax - ymin)
          str1 = "%." // string(min(ceiling(max(abs(log10(dy)),0d0)) + 1,10)) // "f" // c_null_char
          call ipSetupAxisFormat(ImAxis_Y1,c_loc(str1))
          call ipGetPlotCurrentLimits(xmin,xmax,ymin,ymax) ! these need to be switched

          str1 = "Energy" // c_null_char
          auto%x = 0._c_float
          auto%y = 0._c_float
          auto%z = 0._c_float
          auto%w = -1._c_float
          call ipSetNextMarkerStyle(ImPlotMarker_Circle,-1._c_float,auto,-1._c_float,auto)

          call ipPlotLine(c_loc(str1),c_loc(x),c_loc(y),size(x,1),ImPlotLineFlags_None,0_c_int)
          call ipEndPlot()
       end if
    end if

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) doquit = .true.

    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains

    ! terminate the state for this window
    subroutine end_state()
      if (allocated(x)) deallocate(x)
      if (allocated(y)) deallocate(y)
    end subroutine end_state

  end subroutine draw_scfplot

end submodule tree
