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
       sys_loaded_not_init, launch_initialization_thread, ColorTableCellBg,&
       kill_initialization_thread, system_shorten_names, tooltip_delay,&
       ColorDangerButton, ColorFieldSelected, g, fontsize, ok_system
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
    integer :: i, j, k, nshown, jsel, ll, id, iref, inext, iprev, iaux, nfreal
    logical(c_bool) :: ldum, isel
    type(c_ptr) :: ptrc
    type(ImGuiTableSortSpecs), pointer :: sortspecs
    type(ImGuiTableColumnSortSpecs), pointer :: colspecs
    logical :: hadenabledcolumn, isend, ok, found
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
       w%tree_sortcid = ic_tree_id
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
       ! if a system has been rebonded, the "nmol" column may have changed: sort and resize
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

    !! process force options !!
    ! remove systems from the tree
    if (allocated(w%forceremove)) &
       call w%remove_systems_tree(cfilter)
    ! remap the tree index (iord)
    if (w%forceupdate) &
       call w%remap_tree()
    ! sort the tree
    if (w%forcesort) &
       call w%sort_tree()
    ! re-run the initialization threads
    if (w%forceinit) then
       call kill_initialization_thread()
       call launch_initialization_thread()
       w%forceinit = .false.
    end if
    nshown = size(w%iord,1)

    ! final message in the header line
    ok = (nshown > 1)
    if (.not.ok) &
       ok = (sysc(w%iord(1))%status /= sys_empty)
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
    if (igBeginTable(c_loc(str),ic_tree_NUMCOLUMNS,flags,szero,0._c_float)) then
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
       call igTableSetupColumn(c_loc(str),flags,width,ic_tree_closebutton)

       str = "(expand button)##0expandbutton" // c_null_char
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_expandbutton)

       str = "ID##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultSort
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_id)

       str = "Name##0" // c_null_char
       flags = ImGuiTableColumnFlags_WidthStretch
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_name)

       str = "spg##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_spg)

       str = "V/Å³##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_v)

       str = "(V/Z)/Å³##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_vmol)

       str = "nneq##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_nneq)

       str = "nat/ncel##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_ncel)

       str = "nmol##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_nmol)

       str = "a/Å##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_a)

       str = "b/Å##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_b)

       str = "c/Å##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_c)

       str = "α/°##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_alpha)

       str = "β/°##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_beta)

       str = "γ/°##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_gamma)

       str = "E/Ha##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_e)

       str = "(E/Z)/Ha##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_emol)

       str = "p/GPa##0" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,ic_tree_p)

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
             w%tree_sortcid = ic_tree_id
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
                if (igTableSetColumnIndex(ic_tree_closebutton)) then
                   call igAlignTextToFramePadding()
                   str = "##1closebutton" // string(ic_tree_closebutton) // "," // string(i) // c_null_char
                   if (my_CloseButton(c_loc(str),ColorDangerButton)) w%forceremove = (/i/)
                   if (igIsItemHovered(ImGuiHoveredFlags_None)) &
                      tooltipstr = "Close this system"
                end if
             end if

             ! expand button
             if (igTableSetColumnIndex(ic_tree_expandbutton)) then
                if (sysc(i)%collapse < 0) then
                   ! expand button for multi-seed entries
                   str = "##expand" // string(ic_tree_expandbutton) // "," // string(i) // c_null_char
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
                call igTableSetBgColor(ImGuiTableBgTarget_CellBg, color, ic_tree_name)
             end if

             ! ID column
             if (igTableSetColumnIndex(ic_tree_id)) then
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
             if (igTableSetColumnIndex(ic_tree_name)) then
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
                str = ch // "##" // string(ic_tree_name) // "," // string(i) // c_null_char
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
                         call w%select_system_tree(i)
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
             if (igTableSetColumnIndex(ic_tree_e)) then
                if (sysc(i)%seed%energy /= huge(1d0)) then
                   str = string(sysc(i)%seed%energy,'f',decimal=8)
                else
                   str = "n/a"
                end if
                call write_maybe_selectable(i,tooltipstr)
                call iw_text(str,copy_to_output=export)
             end if

             if (sysc(i)%status == sys_init) then
                if (igTableSetColumnIndex(ic_tree_spg)) then ! spg
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

                if (igTableSetColumnIndex(ic_tree_v)) then ! volume
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%omega*bohrtoa**3,'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_tree_vmol)) then ! volume per molecule
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%omega*bohrtoa**3/sys(i)%c%nmol,'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_tree_nneq)) then ! nneq
                   str = string(sys(i)%c%nneq)
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_tree_ncel)) then ! ncel
                   str = string(sys(i)%c%ncel)
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_tree_nmol)) then ! nmol
                   str = string(sys(i)%c%nmol)
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_tree_a)) then ! a
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%aa(1)*bohrtoa,'f',decimal=4)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_tree_b)) then ! b
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%aa(2)*bohrtoa,'f',decimal=4)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_tree_c)) then ! c
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%aa(3)*bohrtoa,'f',decimal=4)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_tree_alpha)) then ! alpha
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%bb(1),'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_tree_beta)) then ! beta
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%bb(2),'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_tree_gamma)) then ! gamma
                   if (sys(i)%c%ismolecule) then
                      str = "<mol>"
                   else
                      str = string(sys(i)%c%bb(3),'f',decimal=2)
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,disabled=(sysc(i)%status /= sys_init),copy_to_output=export)
                end if

                if (igTableSetColumnIndex(ic_tree_emol)) then ! energy/nmol
                   if (sysc(i)%seed%energy /= huge(1d0)) then
                      str = string(sysc(i)%seed%energy/sys(i)%c%nmol,'f',decimal=8)
                   else
                      str = "n/a"
                   end if
                   call write_maybe_selectable(i,tooltipstr)
                   call iw_text(str,copy_to_output=export)
                end if
                if (igTableSetColumnIndex(ic_tree_p)) then ! pressure
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
             if (igTableSetColumnIndex(ic_tree_name)) then
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
         call w%select_system_tree(isys)
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
            call w%select_system_tree(isys)
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
         if (sysc(w%tree_selected)%collapse == i) &
            call w%select_system_tree(i)
      end if
      w%forceupdate = .true.

    end subroutine collapse_system

  end subroutine draw_tree

  ! Update the table rows by building a new row index array
  ! (iord). Only the systems that are not empty are pointed by
  ! iord. This routine is used when the systems change.
  module subroutine remap_tree(w)
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

  end subroutine remap_tree

  ! Sort the table row order by column cid and in direction dir
  ! (ascending=1, descending=2). Modifies the w%iord.
  module subroutine sort_tree(w)
    use interfaces_glfw, only: glfwGetTime
    use gui_main, only: sys, sysc, sys_init, sys_empty
    use tools, only: mergesort
    use tools_math, only: invert_permutation
    use tools_io, only: ferror, faterr
    use types, only: vstring
    class(window), intent(inout) :: w

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
    if (w%tree_sortcid == ic_tree_id.or.w%tree_sortcid == ic_tree_nneq.or.w%tree_sortcid == ic_tree_ncel.or.&
       w%tree_sortcid == ic_tree_nmol) then
       ! sort by integer
       allocate(ival(n))
       do i = 1, n
          if (w%tree_sortcid == ic_tree_id) then
             ival(i) = w%iord(i)
          elseif (sysc(w%iord(i))%status == sys_init) then
             if (w%tree_sortcid == ic_tree_nneq) then
                ival(i) = sys(w%iord(i))%c%nneq
             elseif (w%tree_sortcid == ic_tree_ncel) then
                ival(i) = sys(w%iord(i))%c%ncel
             elseif (w%tree_sortcid == ic_tree_nmol) then
                ival(i) = sys(w%iord(i))%c%nmol
             end if
          else
             ival(i) = huge(1)
             valid(i) = .false.
          end if
       end do
       call mergesort(ival,iperm,1,n)
       deallocate(ival)
    elseif (w%tree_sortcid == ic_tree_v.or.w%tree_sortcid == ic_tree_a.or.w%tree_sortcid == ic_tree_b.or.&
       w%tree_sortcid == ic_tree_c.or.w%tree_sortcid == ic_tree_alpha.or.w%tree_sortcid == ic_tree_beta.or.&
       w%tree_sortcid == ic_tree_gamma.or.w%tree_sortcid == ic_tree_vmol.or.&
       w%tree_sortcid == ic_tree_e.or.w%tree_sortcid == ic_tree_emol.or.w%tree_sortcid == ic_tree_p) then
       ! sort by real
       allocate(rval(n))
       do i = 1, n
          doit = sysc(w%iord(i))%status == sys_init
          if (doit) doit = (.not.sys(w%iord(i))%c%ismolecule)
          if (w%tree_sortcid == ic_tree_v .and. doit) then
             rval(i) = sys(w%iord(i))%c%omega
          elseif (w%tree_sortcid == ic_tree_a .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(1)
          elseif (w%tree_sortcid == ic_tree_b .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(2)
          elseif (w%tree_sortcid == ic_tree_c .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(3)
          elseif (w%tree_sortcid == ic_tree_alpha .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(1)
          elseif (w%tree_sortcid == ic_tree_beta .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(2)
          elseif (w%tree_sortcid == ic_tree_gamma .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(3)
          elseif (w%tree_sortcid == ic_tree_vmol .and. doit) then
             rval(i) = sys(w%iord(i))%c%omega / sys(w%iord(i))%c%nmol
          elseif (w%tree_sortcid == ic_tree_e .and. sysc(w%iord(i))%seed%energy /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%energy
          elseif (w%tree_sortcid == ic_tree_emol .and. sysc(w%iord(i))%status == sys_init .and.&
             sysc(w%iord(i))%seed%energy /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%energy / sys(w%iord(i))%c%nmol
          elseif (w%tree_sortcid == ic_tree_p .and. sysc(w%iord(i))%seed%pressure /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%pressure
          else
             rval(i) = huge(1d0)
             valid(i) = .false.
          end if
       end do
       call mergesort(rval,iperm,1,n)
       deallocate(rval)
    elseif (w%tree_sortcid == ic_tree_name .or. w%tree_sortcid == ic_tree_spg) then
       ! sort by string
       allocate(sval(n))
       do i = 1, n
          if (w%tree_sortcid == ic_tree_name .and. sysc(w%iord(i))%status /= sys_empty.and.&
             .not.sysc(w%iord(i))%hidden) then
             sval(i)%s = trim(sysc(w%iord(i))%seed%name)
          else
             doit = (w%tree_sortcid == ic_tree_spg) .and. (sysc(w%iord(i))%status == sys_init)
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
       if (w%tree_sortdir == 2) then
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

  !> Remove systems given by index idx from the tree.
  module subroutine remove_systems_tree(w,cfilter)
    use gui_main, only: sysc, sys_init, sys_initializing, sys_loaded_not_init,&
       kill_initialization_thread, remove_system
    class(window), intent(inout) :: w
    type(c_ptr), intent(inout) :: cfilter

    logical :: reinit
    integer :: i, j, jsel, k, newsel
    character(kind=c_char,len=:), allocatable, target :: str

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
          call w%select_system_tree(newsel)
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

  end subroutine remove_systems_tree

  !> Select system idx from the tree
  module subroutine select_system_tree(w,idx)
    use gui_main, only: tree_select_updates_inpcon, tree_select_updates_view
    class(window), intent(inout) :: w
    integer, intent(in) :: idx

    w%tree_selected = idx
    if (tree_select_updates_inpcon) &
       win(iwin_console_input)%inpcon_selected = idx
    if (tree_select_updates_view) &
       call win(iwin_view)%select_view(idx)

  end subroutine select_system_tree

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

end submodule tree
