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
  real(c_float), parameter :: rgba_partialocc(4) = (/0.94_c_float,0.60_c_float,0.00_c_float,1.00_c_float/)
  real(c_float), parameter :: rgba_fields(4) = (/0.00_c_float,1.00_c_float,0.50_c_float,1.00_c_float/)
  real(c_float), parameter :: rgba_reference(4) = (/0.890_c_float,0.706_c_float,0.129_c_float,1._c_float/)
  ! color for the tree expand/collapse control-button icon
  real(c_float), parameter :: rgba_expand(4) = (/0.72_c_float,0.72_c_float,0.78_c_float,1.00_c_float/)
  ! group headers, which are not systems and are not waiting to initialize
  real(c_float), parameter :: rgba_group(4) = (/0.95_c_float,0.80_c_float,0.25_c_float,1.00_c_float/)

  ! vertical item spacing pushed around the systems table; the height
  ! reserved under the table has to use it too, because ImGui puts one
  ! of these gaps between the table and whatever follows it
  real(c_float), parameter :: table_itemspacing_y = 5._c_float

  ! deferred multi-selection actions: they act on the systems shown in
  ! the tree, which is known only after the filter has been applied
  integer, parameter :: sel_none = 0
  integer, parameter :: sel_all = 1
  integer, parameter :: sel_invert = 2

  !xx! private procedures
  ! subroutine tree_select_none()
  ! function tree_system_tooltip(i)
  ! function tree_field_tooltip_string(si,fj)

contains

  !> Draw the contents of a tree window
  module subroutine draw_tree(w)
    use interfaces_glfw, only: glfwGetTime
    use windows, only: win
    use keybindings, only: is_bind_event, BIND_TREE_REMOVE_SYSTEM_FIELD, BIND_TREE_MOVE_UP,&
       BIND_TREE_MOVE_DOWN, BIND_TREE_SELECT_ALL
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_inputtext, iw_text,&
       iw_setposx_fromend, iw_calcwidth, iw_calcheight, iw_menuitem, iw_inputint3, iw_icon_button,&
       iw_close_button, iw_table_column, iw_beginmenu, iw_inputfloat, iw_inputint, iw_checkbox
    use systems, only: nsys, sys, sysc, sys_empty, sys_group, sys_init, sys_ready,&
       sys_loaded_not_init, launch_initialization_thread, are_threads_running,&
       kill_initialization_thread, system_shorten_names, ok_system, sys_initializing,&
       remove_systems, reload_field_with_virtuals
    use interfaces_glfw, only: glfwGetTime
    use gui_main, only: ColorTableCellBg, tooltip_delay, errmsg_linger,&
       ColorFieldSelected, ColorTableHighlightRow,&
       ColorTableSelectedBorder, g, fontsize, io
    use icons, only: icon_tex, icon_tex_fmt, rgba_icon_fmt, icon_fmt_MAX, icon_ui_group,&
       format_name, icon_prop_fields, icon_prop_vib, icon_prop_occ,&
       icon_ui_expand, icon_ui_collapse
    use fieldmod, only: type_grid, type_wien, type_pi, type_dftb
    use representations, only: reptype_isosurface, repflavor_isosurface
    use grid3mod, only: mode_nearest, mode_trilinear, mode_trispline, mode_tricubic, mode_smr
    use tools_io, only: string, uout, nameguess
    use types, only: realloc
    use param, only: mlen, bohrtoa, ifformat_as_resample, ifformat_as_ft_x, ifformat_as_ft_y,&
       ifformat_as_ft_z, ifformat_as_ft_xx, ifformat_as_ft_xy, ifformat_as_ft_xz,&
       ifformat_as_ft_yy, ifformat_as_ft_yz, ifformat_as_ft_zz, ifformat_as_ft_grad,&
       ifformat_as_ft_lap, ifformat_as_ft_pot
    use c_interface_module
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, zeroc, ch
    character(kind=c_char,len=:), allocatable :: tooltipstr
    type(ImVec2) :: szero, sz, uv0, uv1
    type(ImVec2) :: selrect_min, selrect_max, fieldrect_min, fieldrect_max
    type(ImVec2) :: tableclip_min, tableclip_max
    logical :: havesel, havefield
    type(c_ptr) :: tabledl
    type(ImVec4) :: col4, tintcol, nobord
    integer(c_int) :: flags, color, itex
    integer :: i, j, k, jsel, iref, inext, iprev, ithis, iaux, ifmt
    integer :: nshown, nshown_after_filter, maxprops, nprop, nrow, ithis_row
    logical :: hasfield, hasvib, hasocc
    logical(c_bool) :: ldum
    type(c_ptr) :: ptrc
    type(ImGuiTableSortSpecs), pointer :: sortspecs
    type(ImGuiTableColumnSortSpecs), pointer :: colspecs
    logical :: hadenabledcolumn, ok, found, reinit
    integer :: ndrawn ! icons already drawn in the current cell
    logical :: export
    real(c_float) :: width, pos, hbottom
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    integer, allocatable :: ishown(:)
    integer, allocatable :: ishown_row(:), ifield_row(:)

    integer, allocatable, save :: forceremove(:) ! enter integers to remove one or more systems
    integer, save :: forceselect = 0 ! force selection of a system
    logical, save :: forcereassign = .false. ! force a check of reassign current selected system
    logical, save :: forceremap = .false. ! force a remap of the tree
    logical, save :: forcesort = .false. ! force a sort of the tree
    logical, save :: forceresize = .false. ! force a resize of the tree columns
    logical, save :: forceinit = .false. ! run the initialization threads
    integer, save :: forceselaction = sel_none ! deferred multi-selection action
    type(c_ptr), save :: cfilter = c_null_ptr ! filter object (allocated first pass, never destroyed)
    logical, save :: ttshown = .false. ! tooltip flag
    integer, save :: lastmaxprops = -1 ! last max number of simultaneous property icons (to detect resize)
    integer(c_int), save :: iresample(3) = (/0,0,0/) ! for the grid resampling menu option
    real(c_float), save :: fieldnorm = 1._c_float ! for the field normalize menu option
    real*8, save :: fieldintg = 0d0 ! current integral shown in the normalize menu option
    logical, save :: fieldintg_open = .false. ! whether the normalize menu was open last frame
    logical, save :: fielddocore = .false. ! for the field core-augmentation menu option
    integer(c_int), allocatable, save :: fieldzpsp(:) ! for the field core-augmentation menu option

    ! initialize
    hadenabledcolumn = .false.
    zeroc = "" // c_null_char
    tooltipstr = ""
    szero%x = 0
    szero%y = 0
    uv0 = ImVec2(0._c_float,0._c_float)
    uv1 = ImVec2(1._c_float,1._c_float)
    nobord = ImVec4(0._c_float,0._c_float,0._c_float,0._c_float)
    if (w%firstpass) then
       w%sortcid = ic_tree_id
       w%sortdir = 1
       w%isys = 1
       forceremap = .true.
       forcesort = .true.
    end if

    ! update the tree based on time signals between dependent windows;
    ! also find the maximum number of property icons shown simultaneously
    ! in any row, to size the properties column
    maxprops = 0
    do i = 1, nsys
       ! if a system has changed fundamentally, the table needs an update (maybe)
       if (w%timelast_tree_update < sysc(i)%timelastchange_geometry) forceremap = .true.
       ! if the currently selected system has changed, maybe we need to reassign
       if (w%isys == i .and. w%timelast_assign < sysc(i)%timelastchange_geometry) &
          forcereassign = .true.
       ! if a system has been rebonded, the "nmol" column may have changed: sort and resize
       if (w%timelast_tree_resize < sysc(i)%timelastchange_rebond) forceresize = .true.
       if (w%timelast_sort < sysc(i)%timelastchange_rebond) forcesort = .true.
       ! count property icons for this system
       if (sysc(i)%status >= sys_ready) then
          nprop = 0
          if (any(sys(i)%f(1:sys(i)%nf)%isinit)) nprop = nprop + 1
          if (sys(i)%c%vib%hasvibs) nprop = nprop + 1
          if (sys(i)%c%haveocc) nprop = nprop + 1
          maxprops = max(maxprops,nprop)
       end if
    end do
    ! resize the columns if the properties column needs a different width
    if (maxprops /= lastmaxprops) then
       forceresize = .true.
       lastmaxprops = maxprops
    end if

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

       ! button: save multiple. It reads the whole tree, so it waits for
       ! the initialization threads, like the entry in the row context menu
       if (iw_menuitem("Save Multiple...",enabled=(nsys > 0 .and..not.are_threads_running()))) &
          iaux = stack_create_window(wintype_save_multiple,.true.,idparent=w%id,orraise=-1)
       call iw_tooltip("Write the structures of several systems in the tree to files",ttshown)
       call igSeparator()

       ! close the selection. A selected group header goes in as well: it
       ! closes the systems it heads, and remove_system returns early on a
       ! member that the header already took with it
       if (iw_menuitem("Close Selected",enabled=(tree_nselected() > 0))) then
          if (allocated(forceremove)) deallocate(forceremove)
          allocate(forceremove(nsys))
          k = 0
          do i = 1, nsys
             if (sysc(i)%status == sys_empty .or. .not.sysc(i)%tselected) cycle
             k = k + 1
             forceremove(k) = i
          end do
          call realloc(forceremove,k)
       end if
       call iw_tooltip("Close the systems selected in the tree",ttshown)

       ! close visible
       if (iw_menuitem("Close Visible")) then
          if (allocated(forceremove)) deallocate(forceremove)
          allocate(forceremove(nsys))
          k = 0
          do i = 1, nsys
             if (c_associated(cfilter)) then
                str = trim(sysc(i)%seed%name) // c_null_char
                if (.not.ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) cycle
             end if
             k = k + 1
             forceremove(k) = i
          end do
          call realloc(forceremove,k)
          if (c_associated(cfilter)) &
             call ImGuiTextFilter_Clear(cfilter)
       end if
       call iw_tooltip("Close all visible systems",ttshown)

       ! close all systems
       if (iw_menuitem("Close All")) then
          if (allocated(forceremove)) deallocate(forceremove)
          allocate(forceremove(nsys))
          do i = 1, nsys
             forceremove(i) = i
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
    if (allocated(forceremove)) then
       reinit = any((sysc(forceremove)%status == sys_loaded_not_init.and..not.sysc(forceremove)%hidden).or.&
          sysc(forceremove)%status == sys_initializing)
       ! closing a group header takes its members with it, so they decide
       ! whether the initialization thread has to stand down as well
       if (.not.reinit) then
          do k = 1, size(forceremove,1)
             if (sysc(forceremove(k))%status /= sys_group) cycle
             reinit = any((sysc(1:nsys)%collapse == forceremove(k)) .and. &
                (sysc(1:nsys)%status == sys_loaded_not_init .or. &
                sysc(1:nsys)%status == sys_initializing))
             if (reinit) exit
          end do
       end if
       if (reinit) call kill_initialization_thread()
       call remove_systems(forceremove)
       if (reinit) call launch_initialization_thread()
       forceremap = .true.
       forcereassign = .true.
       deallocate(forceremove)
    end if
    ! reassign the tree (needs to go before regenerating iord)
    if (forcereassign) then
       call w%reassign_tree(cfilter)
       forcereassign = .false.
    end if
    ! remap the tree index (regenerate iord)
    if (forceremap) then
       call w%remap_tree()
       forceremap = .false.
       forcesort = .true.
    end if
    ! sort the tree
    if (forcesort) then
       call w%sort_tree()
       forcesort = .false.
    end if
    ! re-run the initialization threads
    if (forceinit) then
       call kill_initialization_thread()
       call launch_initialization_thread()
       forceinit = .false.
    end if

    ! external request to select a system
    if (w%forceselect > 0) then
       call w%select_system_tree(w%forceselect)
       forceselect = w%forceselect
       w%forceselect = 0
    end if

    ! count systems, select which ones will be shown; next and previous systems
    nshown_after_filter = 0
    nshown = 0
    inext = 0
    iprev = 0
    ithis = 0
    if (allocated(w%iord)) then
       nshown = size(w%iord,1)
       allocate(ishown(nshown))
       if (c_associated(cfilter)) then
          do j = 1, size(w%iord,1)
             i = w%iord(j)
             if (sysc(i)%status == sys_empty .or. sysc(i)%hidden) cycle
             str = trim(sysc(i)%seed%name) // c_null_char
             if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
                nshown_after_filter = nshown_after_filter + 1
                ishown(nshown_after_filter) = i
                if (i == w%isys) ithis = nshown_after_filter
             end if
          end do
       else
          do i = 1, nshown
             ishown(i) = i
          end do
       end if
       if (ithis > 0) then
          ! neighbors for keyboard navigation: a group header is not a system
          ! to select, so stepping through the list passes over it
          do j = ithis-1, 1, -1
             if (sysc(ishown(j))%status /= sys_group) then
                iprev = ishown(j)
                exit
             end if
          end do
          do j = ithis+1, nshown_after_filter
             if (sysc(ishown(j))%status /= sys_group) then
                inext = ishown(j)
                exit
             end if
          end do
       end if
    end if

    ! Build the flattened list of table rows: one row per shown system,
    ! plus one row per initialized field of the expanded systems.
    havesel = .false.
    havefield = .false.
    nrow = 0
    ithis_row = 0
    if (nshown_after_filter > 0) then
       ! upper bound for the number of rows
       nrow = nshown_after_filter
       do j = 1, nshown_after_filter
          i = ishown(j)
          if (sysc(i)%showfields .and. sysc(i)%status >= sys_ready) &
             nrow = nrow + count(sys(i)%f(0:sys(i)%nf)%isinit)
       end do
       allocate(ishown_row(nrow),ifield_row(nrow))

       ! fill the row list
       nrow = 0
       do j = 1, nshown_after_filter
          i = ishown(j)
          nrow = nrow + 1
          ishown_row(nrow) = i
          ifield_row(nrow) = -1
          if (i == w%isys) ithis_row = nrow
          if (sysc(i)%showfields .and. sysc(i)%status >= sys_ready) then
             do k = 0, sys(i)%nf
                if (.not.sys(i)%f(k)%isinit) cycle
                nrow = nrow + 1
                ishown_row(nrow) = i
                ifield_row(nrow) = k
             end do
          end if
       end do
    end if

    ! apply the deferred multi-selection action, now that the list of
    ! systems shown in the tree is known
    if (forceselaction /= sel_none) then
       ! select all starts from an empty selection; invert does not
       if (forceselaction == sel_all) call tree_select_none()
       do j = 1, nshown_after_filter
          i = ishown(j)
          if (sysc(i)%status == sys_empty) cycle
          call set_selected(i,(forceselaction == sel_all) .or. .not.sysc(i)%tselected)
       end do
       w%tree_selanchor = 0
       forceselaction = sel_none
    end if

    ! final message in the header line
    if (nshown > 1) then
       call iw_text(" " // string(nshown_after_filter) // "/" // string(nshown) // " shown",&
          sameline=.true.)
       k = tree_nselected()
       if (k > 0) &
          call iw_text(", " // string(k) // " selected",sameline_nospace=.true.)
    end if

    ! Height to keep free under the table for the row of buttons. Each
    ! half of it has to be measured with the style that is in force
    ! where it is actually drawn, which is not the same style: the
    ! button comes after the table's style has been popped again, so it
    ! has the window's frame height, while the gap above it is the
    ! spacing the table itself pushed. Getting either from the other's
    ! style leaves the content a few pixels too tall for the window,
    ! and the window grows a scroll bar that has no business there.
    hbottom = igGetFrameHeight() + table_itemspacing_y

    ! set up the table, style and flags
    sz%x = 3._c_float
    sz%y = 1._c_float
    call igPushStyleVar_Vec2(ImGuiStyleVar_FramePadding,sz)
    sz%x = 8._c_float
    sz%y = table_itemspacing_y
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
    sz%x = 0._c_float
    sz%y = -hbottom
    if (igBeginTable(c_loc(str),ic_tree_NUMCOLUMNS,flags,sz,0._c_float)) then
       ! The table scrolls, so it is its own child window: its draw
       ! list is captured for the selection borders drawn at the end,
       ! and with it the clip rect in force here, which is the table's
       ! visible area. By the time the borders are drawn the table has
       ! ended and that clip rect is gone from the list, so a selected
       ! row scrolled out of sight would otherwise put its border over
       ! the tab bar above the table.
       tabledl = igGetWindowDrawList()
       call ImDrawList_GetClipRectMin(tableclip_min,tabledl)
       call ImDrawList_GetClipRectMax(tableclip_max,tabledl)
       ! force resize if asked for
       if (forceresize) then
          call igTableSetColumnWidthAutoAll(igGetCurrentTable())
          forceresize = .false.
          w%timelast_tree_resize = glfwGetTime()
       end if

       ! set up the columns
       ! closebutton - expandbutton - ID - format - properties - name - spg - volume -
       ! nneq - ncel - nmol - zprime - a - b - c - alpha - beta - gamma - e - emol - p
       flags = ImGuiTableColumnFlags_NoResize
       flags = ior(flags,ImGuiTableColumnFlags_NoReorder)
       flags = ior(flags,ImGuiTableColumnFlags_NoHide)
       flags = ior(flags,ImGuiTableColumnFlags_NoSort)
       flags = ior(flags,ImGuiTableColumnFlags_NoHeaderLabel)
       flags = ior(flags,ImGuiTableColumnFlags_NoHeaderWidth)
       call iw_table_column("(close button)##0closebutton",id=ic_tree_closebutton,flags=flags,&
          width=max(4._c_float, fontsize%y + 2._c_float))

       call iw_table_column("(expand button)##0expandbutton",id=ic_tree_expandbutton,flags=flags)

       call iw_table_column("Id##0",id=ic_tree_id,flags=ImGuiTableColumnFlags_DefaultSort)

       flags = ImGuiTableColumnFlags_NoResize
       flags = ior(flags,ImGuiTableColumnFlags_NoSort)
       flags = ior(flags,ImGuiTableColumnFlags_NoHeaderLabel)
       flags = ior(flags,ImGuiTableColumnFlags_NoHeaderWidth)
       call iw_table_column("format##0",id=ic_tree_format,flags=flags,&
          width=max(4._c_float, fontsize%y + 4._c_float))

       call iw_table_column("Properties##0",id=ic_tree_props,flags=flags)

       call iw_table_column("Name##0",id=ic_tree_name,flags=ImGuiTableColumnFlags_WidthStretch)

       call iw_table_column("Sym##0",id=ic_tree_spg,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("V/Å³##0",id=ic_tree_v,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("(V/Z)/Å³##0",id=ic_tree_vmol,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("Nneq##0",id=ic_tree_nneq,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("Nat/ncel##0",id=ic_tree_ncel,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("Z##0",id=ic_tree_nmol,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("Z'##0",id=ic_tree_zprime,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("a/Å##0",id=ic_tree_a,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("b/Å##0",id=ic_tree_b,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("c/Å##0",id=ic_tree_c,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("α/°##0",id=ic_tree_alpha,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("β/°##0",id=ic_tree_beta,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("γ/°##0",id=ic_tree_gamma,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("E/Ha##0",id=ic_tree_e,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("(E/Z)/Ha##0",id=ic_tree_emol,flags=ImGuiTableColumnFlags_DefaultHide)

       call iw_table_column("p/GPa##0",id=ic_tree_p,flags=ImGuiTableColumnFlags_DefaultHide)

       call igTableSetupScrollFreeze(0, 1) ! top row always visible

       ! fetch the sort specs, sort the data if necessary
       ptrc = igTableGetSortSpecs()
       if (c_associated(ptrc)) then
          call c_f_pointer(ptrc,sortspecs)
          if (c_associated(sortspecs%Specs)) then
             call c_f_pointer(sortspecs%Specs,colspecs)
             w%sortcid = colspecs%ColumnUserID
             w%sortdir = colspecs%SortDirection
             if (sortspecs%SpecsDirty .and. nshown > 1) then
                forcesort = .true.
                sortspecs%SpecsDirty = .false.
             end if
          else
             w%sortcid = ic_tree_id
             w%sortdir = 1
          end if
       end if

       ! draw the header
       call igTableHeadersRow()

       ! the big table
       if (allocated(w%iord) .and. nshown_after_filter > 0) then
          ! start the clipper over the flattened row list (systems + fields
          ! of the expanded systems)
          clipper = ImGuiListClipper_ImGuiListClipper()
          call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)
          call ImGuiListClipper_ForceDisplayRangeByIndices(clipper,ithis_row-4,ithis_row+4)

          ! draw the rows
          do while(ImGuiListClipper_Step(clipper))
             call c_f_pointer(clipper,clipper_f)
             do j = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                i = ishown_row(j)

                ! field rows of an expanded system
                if (ifield_row(j) >= 0) then
                   call draw_field_row(i,ifield_row(j))
                   cycle
                end if

                ! start defining the table row
                call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                hadenabledcolumn = .false.

                ! close button
                if (sysc(i)%status == sys_init .or. sysc(i)%status == sys_group) then
                   if (igTableSetColumnIndex(ic_tree_closebutton)) then
                      str = "##1closebutton" // string(ic_tree_closebutton) // "," // string(i)
                      if (iw_close_button(str)) forceremove = (/i/)
                      if (igIsItemHovered(ImGuiHoveredFlags_None)) then
                         if (sysc(i)%status /= sys_group) then
                            tooltipstr = "Close this system"
                         elseif (sysc(i)%collapse == -1) then
                            tooltipstr = "Close this group and all the systems in it"
                         else
                            tooltipstr = "Close this group, keeping its systems"
                         end if
                      end if
                   end if
                end if

                ! expand/collapse button for multi-seed entries (SCF iterations)
                if (igTableSetColumnIndex(ic_tree_expandbutton)) then
                   if (sysc(i)%collapse < 0) then
                      str = "##expand" // string(ic_tree_expandbutton) // "," // string(i)
                      if (sysc(i)%collapse == -1) then
                         itex = icon_tex(icon_ui_expand)
                         ch = "▶"
                      else
                         itex = icon_tex(icon_ui_collapse)
                         ch = "▼"
                      end if
                      if (iw_icon_button(str,itex,rgba_expand,ch)) call toggle_group(i)
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
                   call tree_cell(str)
                end if

                ! format column
                if (igTableSetColumnIndex(ic_tree_format)) then
                   call write_maybe_selectable(i,tooltipstr)
                   ndrawn = 0
                   if (sysc(i)%status == sys_group) then
                      ! a group header has no format: it gets its own icon
                      call draw_icon_cell(.true.,icon_tex(icon_ui_group),rgba_group,&
                         "A group of systems",ndrawn)
                   else
                      ifmt = sysc(i)%seed%isformat
                      if (ifmt < 0 .or. ifmt > icon_fmt_MAX) ifmt = 0
                      call draw_icon_cell(.true.,icon_tex_fmt(ifmt),rgba_icon_fmt(:,ifmt),&
                         "Format: " // format_name(ifmt),ndrawn)
                   end if
                end if

                ! properties column (fields, vibrations, disorder)
                if (igTableSetColumnIndex(ic_tree_props)) then
                   call write_maybe_selectable(i,tooltipstr)
                   hasfield = .false.
                   hasvib = .false.
                   hasocc = .false.
                   if (sysc(i)%status >= sys_ready) then
                      hasfield = any(sys(i)%f(1:sys(i)%nf)%isinit)
                      hasvib = sys(i)%c%vib%hasvibs
                      hasocc = sys(i)%c%haveocc
                   end if
                   ! icons, then pad to maxprops
                   nprop = 0
                   call draw_icon_cell(hasfield,icon_tex(icon_prop_fields),rgba_fields,&
                      "Scalar fields loaded (click to show or hide the field list)",nprop)
                   if (hasfield) then
                      if (igIsItemHovered(ImGuiHoveredFlags_None) .and.&
                         igIsMouseClicked(ImGuiMouseButton_Left,.false._c_bool) .and.&
                         sysc(i)%status /= sys_group) &
                         sysc(i)%showfields = .not.sysc(i)%showfields
                   end if
                   call draw_icon_cell(hasvib,icon_tex(icon_prop_vib),rgba_vibrations,&
                      "Vibration data available",nprop)
                   call draw_icon_cell(hasocc,icon_tex(icon_prop_occ),rgba_partialocc,&
                      "Partial site occupancies (disorder)",nprop)
                   sz = ImVec2(fontsize%y,fontsize%y)
                   do k = nprop+1, maxprops
                      if (k > 1) call igSameLine(0._c_float,2._c_float)
                      call igDummy(sz)
                   end do
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

                   call iw_text(ch)
                   call igSameLine(0._c_float,-1._c_float)
                   call igSetCursorPosX(pos)
                   str = ch // "##" // string(ic_tree_name) // "," // string(i) // c_null_char
                   sz%x = iw_calcwidth(1,1)
                   sz%y = iw_calcheight(1,0)
                   if (igInvisibleButton(c_loc(str),sz,ImGuiButtonFlags_None) .and.&
                      sysc(i)%status /= sys_group) sysc(i)%showfields = .not.sysc(i)%showfields

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
                   if (sysc(i)%status == sys_group) then
                      ! a group header is not pending initialization: colored
                      ! instead of grayed out, so it does not read as inactive
                      call iw_text(str,rgba=rgba_group,sameline_nospace=.true.,copy_to_output=export)
                   else
                      call iw_text(str,disabled=(sysc(i)%status /= sys_init),sameline_nospace=.true.,copy_to_output=export)
                   end if
                end if

                ! energy
                if (igTableSetColumnIndex(ic_tree_e)) then
                   if (sysc(i)%seed%energy /= huge(1d0)) then
                      str = string(sysc(i)%seed%energy,'f',decimal=8)
                   else
                      str = "n/a"
                   end if
                   call tree_cell(str,disable=.false.)
                end if

                if (sysc(i)%status == sys_init) then
                   if (igTableSetColumnIndex(ic_tree_spg)) then ! spg
                      if (sys(i)%c%ismolecule) then
                         if (sys(i)%c%pg%avail) then
                            str = trim(sys(i)%c%pg%symbol)
                         else
                            str = "n/a"
                         end if
                      else
                         if (.not.sys(i)%c%spgavail) then
                            str = "n/a"
                         else
                            str = trim(sys(i)%c%spg%international_symbol)
                         end if
                      end if
                      call tree_cell(str)
                   end if

                   if (igTableSetColumnIndex(ic_tree_v)) then ! volume
                      call tree_cell(cell_or_mol(sys(i)%c%omega*bohrtoa**3,2))
                   end if

                   if (igTableSetColumnIndex(ic_tree_vmol)) then ! volume per molecule
                      ! cell_or_mol evaluates its argument even for a molecule,
                      ! where the result is discarded; guard the division
                      call tree_cell(cell_or_mol(sys(i)%c%omega*bohrtoa**3/max(sys(i)%c%nmol,1),2))
                   end if

                   if (igTableSetColumnIndex(ic_tree_nneq)) then ! nneq
                      call tree_cell(string(sys(i)%c%nneq))
                   end if

                   if (igTableSetColumnIndex(ic_tree_ncel)) then ! ncel
                      call tree_cell(string(sys(i)%c%ncel))
                   end if

                   if (igTableSetColumnIndex(ic_tree_nmol)) then ! Z (nmol)
                      call tree_cell(string(sys(i)%c%nmol))
                   end if

                   if (igTableSetColumnIndex(ic_tree_zprime)) then ! Z'
                      if (sys(i)%c%ismolecule) then
                         str = "<mol>"
                      elseif (sys(i)%c%spgavail) then
                         str = string(nint(real(sys(i)%c%nmol,8) / real(sys(i)%c%neqv,8)))
                      else
                         str = "n/a"
                      end if
                      call tree_cell(str)
                   end if

                   if (igTableSetColumnIndex(ic_tree_a)) then ! a
                      call tree_cell(cell_or_mol(sys(i)%c%aa(1)*bohrtoa,4))
                   end if
                   if (igTableSetColumnIndex(ic_tree_b)) then ! b
                      call tree_cell(cell_or_mol(sys(i)%c%aa(2)*bohrtoa,4))
                   end if
                   if (igTableSetColumnIndex(ic_tree_c)) then ! c
                      call tree_cell(cell_or_mol(sys(i)%c%aa(3)*bohrtoa,4))
                   end if
                   if (igTableSetColumnIndex(ic_tree_alpha)) then ! alpha
                      call tree_cell(cell_or_mol(sys(i)%c%bb(1),2))
                   end if
                   if (igTableSetColumnIndex(ic_tree_beta)) then ! beta
                      call tree_cell(cell_or_mol(sys(i)%c%bb(2),2))
                   end if
                   if (igTableSetColumnIndex(ic_tree_gamma)) then ! gamma
                      call tree_cell(cell_or_mol(sys(i)%c%bb(3),2))
                   end if

                   if (igTableSetColumnIndex(ic_tree_emol)) then ! energy/nmol
                      if (sysc(i)%seed%energy /= huge(1d0)) then
                         str = string(sysc(i)%seed%energy/sys(i)%c%nmol,'f',decimal=8)
                      else
                         str = "n/a"
                      end if
                      call tree_cell(str,disable=.false.)
                   end if
                   if (igTableSetColumnIndex(ic_tree_p)) then ! pressure
                      if (sys(i)%c%ismolecule) then
                         str = "<mol>"
                      elseif (sysc(i)%seed%pressure /= huge(1d0)) then
                         str = string(sysc(i)%seed%pressure,'f',decimal=2)
                      else
                         str = "n/a"
                      end if
                      call tree_cell(str,disable=.false.)
                   end if
                end if
                ! enter new line
                if (export) write (uout,*)
             end do ! row indices
          end do ! clipper step
          call ImGuiListClipper_End(clipper)
          call ImGuiListClipper_destroy(clipper)
       else
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
          if (igTableSetColumnIndex(ic_tree_name)) then
             call iw_text("No systems loaded",alignframe=.true.)
          end if
       end if
       call igEndTable()

       ! selection borders, drawn over the finished table; the system
       ! border comes last so it stays on top of the field border
       if (havefield) call draw_selection_border(fieldrect_min,fieldrect_max,ColorFieldSelected)
       if (havesel) call draw_selection_border(selrect_min,selrect_max,ColorTableSelectedBorder)
    end if
    call igPopStyleVar(3_c_int)

    ! transient error from the field-settings menu items
    if (len_trim(w%errmsg) > 0 .and. glfwGetTime() - w%timelast_errmsg < errmsg_linger) &
       call iw_text(w%errmsg,danger=.true.,wrap=.true.)

    ! multi-selection buttons, under the table. They carry no label, so
    ! that they still fit when the tree panel is narrow; the tooltips and
    ! the "n selected" counter in the header line say what they act on.
    ! The actions that need the list of shown systems are deferred to the
    ! next pass, where it is built
    if (iw_button("All")) forceselaction = sel_all
    call iw_tooltip("Select all the systems visible in the tree&
       & (control-click and shift-click select them by hand)",ttshown)
    if (iw_button("None",sameline=.true.)) then
       call tree_select_none()
       w%tree_selanchor = 0
    end if
    call iw_tooltip("Deselect all systems in the tree",ttshown)
    if (iw_button("Toggle",sameline=.true.)) forceselaction = sel_invert
    call iw_tooltip("Select the visible systems that are not selected, and vice versa",ttshown)

    ! drop any pending selection that was not consumed above
    forceselect = 0

    !! process the keybindings
    ! select all systems
    if (is_bind_event(BIND_TREE_SELECT_ALL).and.w%focused().and..not.io%WantTextInput) &
       forceselaction = sel_all

    ! remove system or field
    if (is_bind_event(BIND_TREE_REMOVE_SYSTEM_FIELD).and.w%focused()) then
       jsel = w%isys
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
          forceremove = (/jsel/)
       end if
    end if
    ! up and down the tree
    if (is_bind_event(BIND_TREE_MOVE_UP)) then
       if (iprev > 0) &
          forceselect = iprev
    elseif (is_bind_event(BIND_TREE_MOVE_DOWN)) then
       if (inext > 0) &
          forceselect = inext
    end if

    ! if exporting, read the export command
    if (export) &
       ldum = read_output_uout(.true.,"[Table export]")

  contains

    !> Write str into the tree table cell for the current row and column: the
    !> row-spanning selectable first, then the text, echoed to the output when
    !> the tree is being exported. disable = grey the text out when the system
    !> is not initialized (default true; false for the columns whose value comes
    !> from the seed and is therefore available beforehand).
    subroutine tree_cell(str,disable)
      character(len=*), intent(in) :: str
      logical, intent(in), optional :: disable

      logical :: disable_

      disable_ = .true.
      if (present(disable)) disable_ = disable

      call write_maybe_selectable(i,tooltipstr)
      call iw_text(str,disabled=(disable_ .and. sysc(i)%status /= sys_init),copy_to_output=export)

    end subroutine tree_cell

    !> Format the cell value val with dec decimals, or "<mol>" if the current
    !> system is a molecule and the quantity therefore has no meaning.
    function cell_or_mol(val,dec) result(str)
      real*8, intent(in) :: val
      integer, intent(in) :: dec
      character(len=:), allocatable :: str

      if (sys(i)%c%ismolecule) then
         str = "<mol>"
      else
         str = string(val,'f',decimal=dec)
      end if

    end function cell_or_mol

    subroutine write_maybe_selectable(isys,tooltipstr)
      use systems, only: are_threads_running, duplicate_system, reread_system_from_file,&
         sys_group, group_is_scf, group_master
      use keybindings, only: BIND_GEOMETRY
      use utils, only: iw_text, iw_inputtext, iw_beginmenu
      use global, only: iunit, iunit_bohr, iunit_ang
      use tools_io, only: uout
      use param, only: mlen
      integer, intent(in) :: isys
      character(kind=c_char,len=:), allocatable, intent(in) :: tooltipstr

      integer :: k, idx, iaux, ianchor
      real(c_float) :: pos
      integer(c_int) :: flags, isyscollapse, idum
      integer :: npop
      logical(c_bool) :: selected
      logical :: enabled, enabled_no_threads
      logical :: ok, okmouse
      character(kind=c_char,len=:), allocatable, target :: strl

      if (hadenabledcolumn) return

      ! selectable that spans all columns, with zero width
      pos = igGetCursorPosX()
      flags = ImGuiSelectableFlags_SpanAllColumns
      flags = ior(flags,ImGuiSelectableFlags_AllowItemOverlap)
      flags = ior(flags,ImGuiSelectableFlags_AllowDoubleClick)
      flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
      selected = (w%isys==isys) .or. sysc(isys)%tselected
      strl = "##selectable" // string(isys) // c_null_char
      ! any row in the multi-selection gets the highlight-row fill; the
      ! current system, when not selected, gets no fill. It also gets a
      ! border (drawn after the table, so it lands on top of everything)
      if (sysc(isys)%tselected) then
         col4 = ImVec4(ColorTableHighlightRow(1),ColorTableHighlightRow(2),&
            ColorTableHighlightRow(3),ColorTableHighlightRow(4))
         call igPushStyleColor_Vec4(ImGuiCol_Header,col4)
         col4%w = min(col4%w * 1.4_c_float,1._c_float)
         call igPushStyleColor_Vec4(ImGuiCol_HeaderHovered,col4)
         col4%w = min(col4%w * 1.4_c_float,1._c_float)
         call igPushStyleColor_Vec4(ImGuiCol_HeaderActive,col4)
         npop = 3
      elseif (w%isys == isys) then
         col4 = ImVec4(0._c_float,0._c_float,0._c_float,0._c_float)
         call igPushStyleColor_Vec4(ImGuiCol_Header,col4)
         npop = 1
      else
         npop = 0
      end if
      ok = igSelectable_Bool(c_loc(strl),selected,flags,szero)
      if (npop > 0) call igPopStyleColor(int(npop,c_int))
      if (w%isys == isys) then
         call igGetItemRectMin(selrect_min)
         call igGetItemRectMax(selrect_max)
         havesel = .true.
      end if
      okmouse = ok
      ok = ok .or. (forceselect == isys)

      ! multi-selection: control-click toggles one row, shift-click selects
      ! the range from the anchor. Neither changes the current system.
      ! the window type starts with isys = 1 and no anchor, so fall back to
      ! the current system: a shift-click must work on the first click too
      ianchor = w%tree_selanchor
      if (ianchor == 0) ianchor = w%isys
      if (okmouse .and. igIsKeyDown(ImGuiKey_ModShift) .and.&
         range_is_selectable(ianchor,isys)) then
         call select_range(ianchor,isys)
      elseif (okmouse .and. igIsKeyDown(ImGuiKey_ModCtrl)) then
         call set_selected(isys,.not.sysc(isys)%tselected)
         w%tree_selanchor = isys
      elseif (ok) then
         ! a plain activation selects this row alone and starts a new
         ! range; a selection forced by another window (a new structure
         ! being opened, say) leaves the multi-selection alone
         if (okmouse) call tree_select_none()
         w%tree_selanchor = isys
         if (sysc(isys)%status /= sys_group) then
            if (okmouse) call set_selected(isys,.true.)
            call w%select_system_tree(isys)
            if (forceselect > 0) then
               forceselect = 0
               call igSetKeyboardFocusHere(0)
            end if
            if (igIsMouseDoubleClicked(ImGuiPopupFlags_MouseButtonLeft)) &
               sysc(isys)%showfields = .not.sysc(isys)%showfields
         elseif (okmouse) then
            ! a group header is not a system to display: a click on it opens
            ! or closes the group instead
            call toggle_group(isys)
         else
            ! ... and it is never the current system, so the keyboard moves
            ! through it instead of stopping on it
            forceselect = 0
         end if
      end if
      call igSameLine(0._c_float,-1._c_float)
      call igSetCursorPosX(pos)

      enabled = (sysc(isys)%status == sys_init)
      enabled_no_threads = enabled.and..not.are_threads_running()

      ! right click to open the context menu
      if (igBeginPopupContextItem(c_loc(strl),ImGuiPopupFlags_MouseButtonRight)) then
         ! set as current system, with a mark on the current one
         if (iw_menuitem("Set as Current System",enabled=enabled_no_threads,&
            selected=(w%isys == isys))) &
            call w%select_system_tree(isys)
         call iw_tooltip("Set this system as the current system",ttshown)

         ! rename option (system)
         if (iw_beginmenu("Rename",sysc(isys)%status /= sys_group)) then
            if (iw_inputtext("##inputrename",bufsize=mlen-1,textf=sysc(isys)%seed%name,grabfocus=.true.,&
               notlive=.true.,flags=ImGuiInputTextFlags_AutoSelectAll)) then
               sysc(isys)%renamed = .true.
               call igCloseCurrentPopup()
            end if
            call igEndMenu()
         end if
         call iw_tooltip("Rename this system",ttshown)

         call igSeparator()

         ! duplicate system
         if (iw_menuitem("Duplicate",enabled=enabled_no_threads)) &
            call duplicate_system(isys)
         call iw_tooltip("Create a copy of this system",ttshown)

         ! display in a new view window
         if (iw_menuitem("Display in New View",enabled=enabled_no_threads)) then
            idum = stack_create_window(wintype_view,.true.,purpose=wpurp_view_alternate)
            win(idum)%sc = sysc(isys)%sc
            ! the value copy aliased the source scene's GL handles; detach so the
            ! new view lazily builds its own instance buffers
            call win(idum)%sc%gl%detach()
            win(idum)%isys = isys
         end if
         call iw_tooltip("Display this system in a new view window",ttshown)

         ! scf energy plot: only for a group of systems that are the steps
         ! of one calculation, which is headed by the last of those steps.
         ! A group headed by a header entry has no such sequence
         isyscollapse = 0
         if (group_is_scf(isys)) isyscollapse = group_master(isys)
         if (isyscollapse > 0) then
            if (iw_menuitem("Plot SCF Iterations...",enabled=enabled_no_threads)) then
               if (sysc(isyscollapse)%status == sys_init) &
                  iaux = stack_create_window(wintype_scfplot,.true.,isys=isyscollapse,orraise=-1)
            end if
            call iw_tooltip("Plot the energy and other properties as a function of SCF cycle iterations",ttshown)
         end if

         call igSeparator()
         ldum = iw_menuitem("[Contents]",enabled=.false.)

         ! geometry: the geometry window edits what its view shows, so
         ! point the main view at this row before opening it there
         if (iw_menuitem("View/Edit Geometry...",BIND_GEOMETRY,enabled=enabled)) then
            call win(iwin_view)%select_view(isys)
            idum = stack_create_window(wintype_geometry,.true.,idparent=iwin_view,orraise=-1)
         end if
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
            ldum = read_output_uout(.true.,"[Describe system " // string(isys) // "]")
         end if
         call iw_tooltip("Print a detailed description of this system in the output console",ttshown)

         ! fields submenu (system)
         if (iw_beginmenu("Fields")) then
            ! load field
            if (iw_menuitem("Load Field...",enabled=enabled)) &
               iaux = stack_create_window(wintype_load_field,.true.,isys=isys,orraise=-1)
            call iw_tooltip("Load a scalar field for this system",ttshown)

            ! load field from file
            if (iw_menuitem("Load Field From File...",enabled=enabled)) then
               iaux = stack_create_window(wintype_load_field,.true.,isys=isys,orraise=-1)
               win(iaux)%lf%gotofile = .true.
            end if
            call iw_tooltip("Load a scalar field for this system from a file,&
               & choosing the file first",ttshown)

            ! remove all fields, destructive
            ok = enabled
            if (ok) ok = (sys(isys)%nf > 0)
            if (ok) ok = any(sys(isys)%f(1:sys(isys)%nf)%isinit)
            if (iw_menuitem("Remove All Fields",enabled=ok,danger=.true.)) then
               do k = 1, sys(isys)%nf
                  if (.not.sys(isys)%f(k)%isinit) cycle
                  call sys(isys)%unload_field(k)
               end do
            end if
            call iw_tooltip("Remove all fields in this system",ttshown)

            call igEndMenu()
         end if

         ! vibrations submenu (system)
         if (iw_beginmenu("Vibrations")) then
            ! load vibration data
            if (iw_menuitem("Load Vibration Data...",enabled=enabled)) then
               iaux = stack_create_window(wintype_dialog,.true.,purpose=wpurp_dialog_openvibfile,&
                  isys=isys,orraise=-1)
            end if
            call iw_tooltip("Load vibration data from an external file for this system",ttshown)

            ! clear vibration data, destructive
            if (iw_menuitem("Clear Vibration Data",enabled=enabled,danger=.true.)) &
               call sys(isys)%c%vib%end()
            call iw_tooltip("Clear the vibration data for this system",ttshown)

            call igEndMenu()
         end if

         ! new isosurface of the reference field: the new object defaults
         ! to the reference field, so pointing the view at this system
         ! before adding the representation is all that is needed
         if (iw_menuitem("New Isosurface",enabled=enabled)) then
            call win(iwin_view)%select_view(isys)
            call w%select_system_tree(isys)
            call win(iwin_view)%add_rep_and_edit(reptype_isosurface,repflavor_isosurface)
         end if
         call iw_tooltip("Display an isosurface of the reference field in the view",ttshown)

         ! save several systems to files. It acts on the tree and not on
         ! this system, so it goes apart from the entries above
         call igSeparator()
         if (iw_menuitem("Save Multiple...",enabled=(nsys > 0 .and..not.are_threads_running()))) &
            idum = stack_create_window(wintype_save_multiple,.true.,idparent=w%id,orraise=-1)
         call iw_tooltip("Write the structures of several systems in the tree to files",ttshown)

         ! destructive actions, in red at the bottom
         call igSeparator()

         ! reopen from file
         if (iw_menuitem("Restore from File",enabled=enabled_no_threads,danger=.true.)) then
            call reread_system_from_file(isys)
         end if
         call iw_tooltip("Read the file for this system and reopen it, discarding&
            & the current structure (only the last structure is read)",ttshown)

         ! close option (system). A group header can be closed even though
         ! it is not a system: that closes the systems it heads
         if (iw_menuitem("Close",enabled=(enabled .or. sysc(isys)%status == sys_group),&
            danger=.true.)) &
            forceremove = (/isys/)
         if (sysc(isys)%status == sys_group) then
            call iw_tooltip("Close all the systems in this group",ttshown)
         else
            call iw_tooltip("Close this system",ttshown)
         end if

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

    !> Draw one property icon (tinted texture tex, color rgba, with the
    !> given tooltip) if show is true and the texture is available. Icons
    !> are packed with no gaps: ndrawn is the number of icons already drawn
    !> in this cell; it is incremented when this icon is drawn.
    subroutine draw_icon_cell(show,tex,rgba,tooltip,ndrawn)
      logical, intent(in) :: show
      integer(c_int), intent(in) :: tex
      real(c_float), intent(in) :: rgba(4)
      character(len=*), intent(in) :: tooltip
      integer, intent(inout) :: ndrawn

      type(ImVec2) :: sz
      type(ImVec4) :: tintcol
      character(kind=c_char,len=:), allocatable, target :: strl

      if (.not.show .or. tex == 0) return
      if (ndrawn > 0) call igSameLine(0._c_float,2._c_float)
      ndrawn = ndrawn + 1
      sz = ImVec2(fontsize%y,fontsize%y)
      tintcol = ImVec4(rgba(1),rgba(2),rgba(3),rgba(4))
      call igImage(int(tex,c_intptr_t),sz,uv0,uv1,tintcol,nobord)
      if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
         strl = tooltip // c_null_char
         call igSetTooltip(c_loc(strl))
      end if

    end subroutine draw_icon_cell

    !> Apply a SETFIELD-style option string opt to field k of system i,
    !> reporting any error to the output console.
    subroutine field_setopt(i,k,opt)
      integer, intent(in) :: i, k
      character(len=*), intent(in) :: opt

      character(len=:), allocatable :: errmsg

      call sys(i)%f(k)%set_options(opt,errmsg)
      if (len_trim(errmsg) > 0) then
         write (uout,'("!! Warning !! Error setting field option: ",A)') trim(errmsg)
         ! show the error transiently in the tree window, too
         w%errmsg = "Error setting option for field " // string(k) // ": " // trim(errmsg)
         w%timelast_errmsg = glfwGetTime()
      end if

    end subroutine field_setopt

    ! Draw a border around a row of the tree table given by its corners
    ! (pmin,pmax) in the color rgba. The thickness is 8% of the row
    ! height, so it scales with the font size and the DPI, with a
    ! 2-pixel floor so the line does not vanish into the table's row
    ! separators. The draw list is clipped to the current column, so the
    ! clip rect is overridden while the border is drawn -- but only
    ! within the table's own visible area, or a row scrolled out of
    ! sight would draw its border on the widgets above or below the
    ! table. Nothing is drawn once the row has scrolled out entirely.
    subroutine draw_selection_border(pmin,pmax,rgba)
      type(ImVec2), intent(in) :: pmin, pmax
      real(c_float), intent(in) :: rgba(4)

      type(ImVec4) :: colb
      real(c_float) :: thick
      type(ImVec2) :: cmin, cmax

      thick = max(2._c_float,0.12_c_float * (pmax%y - pmin%y))
      cmin%x = max(pmin%x-thick,tableclip_min%x)
      cmin%y = max(pmin%y-thick,tableclip_min%y)
      cmax%x = min(pmax%x+thick,tableclip_max%x)
      cmax%y = min(pmax%y+thick,tableclip_max%y)
      if (cmax%x <= cmin%x .or. cmax%y <= cmin%y) return
      colb = ImVec4(rgba(1),rgba(2),rgba(3),rgba(4))
      call ImDrawList_PushClipRect(tabledl,cmin,cmax,.false._c_bool)
      call ImDrawList_AddRect(tabledl,pmin,pmax,igGetColorU32_Vec4(colb),&
         0._c_float,0_c_int,thick)
      call ImDrawList_PopClipRect(tabledl)

    end subroutine draw_selection_border

    !> Draw the table row corresponding to field k of system i (one of
    !> the rows below an expanded system). The row contains a
    !> full-width selectable in the name column with the field name,
    !> the field context menu, and the field tooltip.
    subroutine draw_field_row(i,k)
      integer, intent(in) :: i, k

      character(kind=c_char,len=:), allocatable, target :: str
      character(len=:), allocatable :: errmsg
      logical(c_bool) :: isel, ldum
      logical :: isend, okz
      integer(c_int) :: flags, color
      type(ImVec4) :: col4
      integer :: id, is, irep

      ! start the row and color the name cell like the parent system's
      call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
      if (sysc(i)%status >= sys_ready) then
         col4 = ImVec4(ColorTableCellBg(1,sys(i)%c%iperiod),ColorTableCellBg(2,sys(i)%c%iperiod),&
            ColorTableCellBg(3,sys(i)%c%iperiod),ColorTableCellBg(4,sys(i)%c%iperiod))
         color = igGetColorU32_Vec4(col4)
         call igTableSetBgColor(ImGuiTableBgTarget_CellBg, color, ic_tree_name)
      end if

      ! red X to remove this field, in the properties column; the
      ! promolecular field (0) cannot be removed
      if (k > 0) then
         if (igTableSetColumnIndex(ic_tree_props)) then
            str = "##fieldclose" // string(i) // "," // string(k)
            if (iw_close_button(str)) then
               call sys(i)%unload_field(k)
               return
            end if
            call iw_tooltip("Remove this field",ttshown)
         end if
      end if

      if (.not.igTableSetColumnIndex(ic_tree_name)) return

      ! selectable; the indent leaves the vertical bar of the corner
      ! character aligned with the expand triangle of the system row
      call igSetCursorPosX(igGetCursorPosX() + 2._c_float * g%Style%FramePadding%x)
      isend = (k == sys(i)%nf)
      if (.not.isend) isend = all(.not.sys(i)%f(k+1:)%isinit)
      if (.not.isend) call iw_text("┌",noadvance=.true.)
      ! the reference field carries no text marker: the border marks it
      str = "└─►(" // string(k) // "): " // trim(sys(i)%f(k)%name) // "##field" // &
         string(i) // "," // string(k) // c_null_char
      isel = (w%isys==i) .and. (sys(i)%iref == k)
      col4 = ImVec4(0._c_float,0._c_float,0._c_float,0._c_float)
      call igPushStyleColor_Vec4(ImGuiCol_Header,col4)
      flags = ImGuiSelectableFlags_SpanAllColumns
      if (igSelectable_Bool(c_loc(str),isel,flags,szero)) then
         call w%select_system_tree(i)
         call sys(i)%set_reference(k,.false.)
      end if
      call igPopStyleColor(1)

      ! the reference field of the current system gets a border, drawn
      ! after the table so the system border lands on top of it
      if (isel) then
         call igGetItemRectMin(fieldrect_min)
         call igGetItemRectMax(fieldrect_max)
         havefield = .true.
      end if

      ! right click to open the field context menu
      if (igBeginPopupContextItem(c_loc(str),ImGuiPopupFlags_MouseButtonRight)) then
         ! set as reference (same as a left-click on the row)
         if (iw_menuitem("Set as Reference",selected=(sys(i)%iref == k))) then
            call w%select_system_tree(i)
            call sys(i)%set_reference(k,.false.)
         end if
         call iw_tooltip("Make this the reference field of the system",ttshown)

         ! rename option (fields)
         if (iw_beginmenu("Rename")) then
            if (iw_inputtext("##inputrenamefield",bufsize=mlen-1,textf=sys(i)%f(k)%name,grabfocus=.true.,&
               notlive=.true.,flags=ImGuiInputTextFlags_AutoSelectAll)) &
               call igCloseCurrentPopup()
            call igEndMenu()
         end if
         call iw_tooltip("Rename this field",ttshown)

         call igSeparator()

         ! duplicate option (fields)
         if (iw_menuitem("Duplicate")) then
            id = sys(i)%getfieldnum()
            call sys(i)%field_copy(k,id)
            sys(i)%f(id)%id = id
            sys(i)%f(id)%name = trim(sys(i)%f(k)%name) // " (copy)"
         end if
         call iw_tooltip("Load a copy of this field as a new field",ttshown)

         ! reload with virtual orbitals (wavefunction fields whose file
         ! has virtuals that were not read)
         if (iw_menuitem("Reload with Virtual Orbitals",enabled=sys(i)%f(k)%has_unread_virtuals())) then
            call reload_field_with_virtuals(i,k,errmsg)
            ! the occupied orbitals are unchanged: keep their cached grids
            if (len_trim(errmsg) == 0) call mo_cache_keep_after_reload(i,k)
            if (len_trim(errmsg) > 0) then
               write (uout,'("!! Warning !! Could not reload the field: ",A)') trim(errmsg)
               ! show the error transiently in the tree window, too
               w%errmsg = "Error reloading field " // string(k) // ": " // trim(errmsg)
               w%timelast_errmsg = glfwGetTime()
            end if
         end if
         call iw_tooltip("Re-read the wavefunction file for this field, this time&
            & including the virtual (unoccupied) orbitals",ttshown)

         !! field settings !!
         call igSeparator()
         ldum = iw_menuitem("[Settings]",enabled=.false.)

         ! interpolation mode and normalization (grids)
         if (sys(i)%f(k)%type == type_grid) then
            if (iw_beginmenu("Interpolation")) then
               if (iw_menuitem("Nearest",selected=(sys(i)%f(k)%grid%mode == mode_nearest))) &
                  call field_setopt(i,k,"nearest")
               call iw_tooltip("Value of the nearest grid point",ttshown)
               if (iw_menuitem("Tri-linear",selected=(sys(i)%f(k)%grid%mode == mode_trilinear))) &
                  call field_setopt(i,k,"trilinear")
               call iw_tooltip("Three-dimensional linear interpolation",ttshown)
               if (iw_menuitem("Tri-spline",selected=(sys(i)%f(k)%grid%mode == mode_trispline))) &
                  call field_setopt(i,k,"trispline")
               call iw_tooltip("Three-dimensional spline interpolation",ttshown)
               if (iw_menuitem("Tri-cubic",selected=(sys(i)%f(k)%grid%mode == mode_tricubic))) &
                  call field_setopt(i,k,"tricubic")
               call iw_tooltip("Three-dimensional cubic polynomial interpolation",ttshown)
               if (iw_menuitem("Smoothrho",selected=(sys(i)%f(k)%grid%mode == mode_smr))) &
                  call field_setopt(i,k,"smoothrho")
               call iw_tooltip("Polyharmonic splines + smoothing, all-electron densities&
                  & only (recommended for them)",ttshown)
               call igEndMenu()
            end if
            call iw_tooltip("Interpolation method used to evaluate the grid between its points",ttshown)

            if (iw_beginmenu("Normalize")) then
               ! seed the displayed integral and the value when the
               ! submenu opens (a full pass over the grid)
               if (.not.fieldintg_open) then
                  fieldintg = sys(i)%f(k)%grid%cell_integral(sys(i)%c%omega)
                  fieldnorm = real(fieldintg,c_float)
                  fieldintg_open = .true.
               end if
               call iw_text("Current integral: " // string(fieldintg,'f',decimal=6))
               ldum = iw_inputfloat("Value##fieldmenunorm",fieldnorm,width=10)
               if (iw_menuitem("OK##fieldmenunormok")) &
                  call field_setopt(i,k,"normalize " // string(fieldnorm,'f',decimal=6))
               call igEndMenu()
            else
               fieldintg_open = .false.
            end if
            call iw_tooltip("Rescale the grid values so its cell integral equals the given value",ttshown)
         end if

         ! core augmentation (any field type)
         if (iw_beginmenu("Core Augmentation")) then
            ldum = iw_checkbox("Activate core augmentation##fieldmenudocore",fielddocore)
            call iw_tooltip("Augment the field with core density contributions, built from&
               & the pseudopotential charges (ZPSP) below: the number of valence electrons&
               & of each species (leave at -1 to skip a species)",ttshown)
            if (fielddocore) then
               do is = 1, sys(i)%c%nspc
                  ldum = iw_inputint(trim(sys(i)%c%spc(is)%name) // "##fieldmenuzpsp" // string(is),&
                     fieldzpsp(is),width=5)
               end do
            end if
            okz = .not.fielddocore
            if (.not.okz) okz = any(fieldzpsp >= 0)
            if (iw_menuitem("OK##fieldmenucoreok",enabled=okz)) then
               if (fielddocore) then
                  ! the SETFIELD parser accepts species names or element symbols;
                  ! emit one pair per element (first species wins)
                  str = "zpsp"
                  do is = 1, sys(i)%c%nspc
                     if (fieldzpsp(is) < 0) cycle
                     if (any(sys(i)%c%spc(1:is-1)%z == sys(i)%c%spc(is)%z .and. fieldzpsp(1:is-1) >= 0)) cycle
                     str = str // " " // trim(nameguess(sys(i)%c%spc(is)%z,.true.)) // " " //&
                        string(fieldzpsp(is))
                  end do
                  call field_setopt(i,k,str)
               else
                  call field_setopt(i,k,"nocore")
               end if
            end if
            call igEndMenu()
         else
            ! seed the menu state from the field while the menu is closed
            fielddocore = sys(i)%f(k)%usecore
            if (allocated(fieldzpsp)) then
               if (size(fieldzpsp,1) /= sys(i)%c%nspc) deallocate(fieldzpsp)
            end if
            if (.not.allocated(fieldzpsp)) allocate(fieldzpsp(sys(i)%c%nspc))
            fieldzpsp = -1
            if (allocated(sys(i)%f(k)%zpsp)) then
               is = min(size(sys(i)%f(k)%zpsp,1),sys(i)%c%nspc)
               fieldzpsp(1:is) = sys(i)%f(k)%zpsp(1:is)
            end if
         end if
         call iw_tooltip("Augment the field with core density contributions",ttshown)

         ! derivatives
         if (iw_beginmenu("Derivatives")) then
            if (iw_menuitem("Numerical",selected=sys(i)%f(k)%numerical)) &
               call field_setopt(i,k,"numerical")
            call iw_tooltip("Calculate the field derivatives numerically",ttshown)
            if (iw_menuitem("Analytical",selected=.not.sys(i)%f(k)%numerical)) &
               call field_setopt(i,k,"analytical")
            call iw_tooltip("Calculate the field derivatives analytically",ttshown)
            call igEndMenu()
         end if
         call iw_tooltip("Method for calculating the field derivatives",ttshown)

         ! typnuc
         if (iw_beginmenu("Nuclear critical point signature")) then
            if (iw_menuitem("-3 (maximum)",selected=(sys(i)%f(k)%typnuc == -3))) &
               call field_setopt(i,k,"typnuc -3")
            if (iw_menuitem("-1",selected=(sys(i)%f(k)%typnuc == -1))) &
               call field_setopt(i,k,"typnuc -1")
            if (iw_menuitem("+1",selected=(sys(i)%f(k)%typnuc == 1))) &
               call field_setopt(i,k,"typnuc 1")
            if (iw_menuitem("+3 (minimum)",selected=(sys(i)%f(k)%typnuc == 3))) &
               call field_setopt(i,k,"typnuc 3")
            call igEndMenu()
         end if
         call iw_tooltip("Signature of the critical points at the nuclear positions&
            & (-3 for fields with maxima at the nuclei, like the density)",ttshown)

         ! exact/approximate sums (aiPI and DFTB+)
         if (sys(i)%f(k)%type == type_pi .or. sys(i)%f(k)%type == type_dftb) then
            if (iw_beginmenu("Sum Evaluation")) then
               if (iw_menuitem("Exact",selected=sys(i)%f(k)%exact)) &
                  call field_setopt(i,k,"exact")
               if (iw_menuitem("Approximate",selected=.not.sys(i)%f(k)%exact)) &
                  call field_setopt(i,k,"approximate")
               call igEndMenu()
            end if
            call iw_tooltip("Evaluate the lattice sums exactly or approximately",ttshown)
         end if

         ! wien2k normalization convention
         if (sys(i)%f(k)%type == type_wien) then
            if (allocated(sys(i)%f(k)%wien)) then
               if (iw_beginmenu("Normalization Convention")) then
                  if (iw_menuitem("Density (RHONORM)",selected=sys(i)%f(k)%wien%cnorm)) &
                     call field_setopt(i,k,"rhonorm")
                  if (iw_menuitem("Potential (VNORM)",selected=.not.sys(i)%f(k)%wien%cnorm)) &
                     call field_setopt(i,k,"vnorm")
                  call igEndMenu()
               end if
               call iw_tooltip("Normalization convention of the WIEN2k field&
                  & (density-style or potential-style)",ttshown)
            end if
         end if

         !! load new field options !!
         call igSeparator()
         ldum = iw_menuitem("[New field from this one]",enabled=.false.)

         ! grid calculation options
         if (sys(i)%f(k)%type == type_grid) then
            if (iw_beginmenu("Load Fourier-Transformed Grid")) then
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

            if (iw_beginmenu("Load Resampled Grid")) then
               ldum = iw_inputint3("New Size##resamplefieldmenunewsize",iresample,width=4*3)

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

         ! new isosurface object bound to this field
         if (iw_menuitem("New Isosurface",enabled=ok_system(i,sys_init))) then
            ! the main view must show this system: the object lives in
            ! its scene and the object editor follows the view
            call win(iwin_view)%select_view(i)
            call w%select_system_tree(i)
            call win(iwin_view)%add_rep_and_edit(reptype_isosurface,repflavor_isosurface,id=irep)
            if (irep > 0) call win(iwin_view)%sc%rep(irep)%iso%set_field(i,k)
         end if
         call iw_tooltip("Display an isosurface of this field in the view",ttshown)

         ! remove option (fields), at the bottom and in red: it destroys
         ! the field with no confirmation
         if (k > 0) then
            call igSeparator()
            if (iw_menuitem("Remove",danger=.true.)) &
               call sys(i)%unload_field(k)
            call iw_tooltip("Remove this field",ttshown)
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

    end subroutine draw_field_row

    ! Set the tree multi-selection flag of system i to val. A group
    ! header is not a system: selecting it selects all its members.
    subroutine set_selected(i,val)
      integer, intent(in) :: i
      logical, intent(in) :: val

      integer :: k

      sysc(i)%tselected = val
      if (sysc(i)%status /= sys_group) return
      do k = 1, nsys
         if (sysc(k)%status == sys_empty) cycle
         if (sysc(k)%collapse == i) sysc(k)%tselected = val
      end do

    end subroutine set_selected

    ! Whether both systems are rows currently shown in the tree, and can
    ! therefore delimit a selection range
    function range_is_selectable(ianchor,i)
      integer, intent(in) :: ianchor, i
      logical :: range_is_selectable

      integer :: j

      range_is_selectable = .false.
      if (ianchor == 0) return
      do j = 1, nshown_after_filter
         if (ishown(j) == ianchor) then
            range_is_selectable = .true.
            exit
         end if
      end do
      if (.not.range_is_selectable) return
      range_is_selectable = .false.
      do j = 1, nshown_after_filter
         if (ishown(j) == i) then
            range_is_selectable = .true.
            exit
         end if
      end do

    end function range_is_selectable

    ! Select the rows between systems ianchor and i (both included), in
    ! the order in which they are shown in the tree.
    subroutine select_range(ianchor,i)
      integer, intent(in) :: ianchor, i

      integer :: j, j1, j2

      j1 = 0
      j2 = 0
      do j = 1, nshown_after_filter
         if (ishown(j) == ianchor) j1 = j
         if (ishown(j) == i) j2 = j
      end do
      if (j1 == 0 .or. j2 == 0) return

      call tree_select_none()
      do j = min(j1,j2), max(j1,j2)
         call set_selected(ishown(j),.true.)
      end do

    end subroutine select_range

    ! open a collapsed group, or close an open one
    subroutine toggle_group(i)
      integer, intent(in) :: i

      if (sysc(i)%collapse == -1) then
         call expand_system(i)
      elseif (sysc(i)%collapse == -2) then
         call collapse_system(i)
      end if

    end subroutine toggle_group

    ! un-hide the dependents and set as expanded
    subroutine expand_system(i)
      integer, intent(in) :: i

      integer :: k

      if (sysc(i)%status == sys_empty) return
      if (sysc(i)%collapse /= -1) return
      do k = 1, nsys
         if (sysc(k)%collapse == i) then
            if (sysc(k)%status == sys_loaded_not_init) forceinit = .true.
            sysc(k)%hidden = .false.
         end if
      end do
      sysc(i)%collapse = -2
      forceremap = .true.

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
      if (w%isys >= 1 .and. w%isys <= nsys) then
         if (sysc(w%isys)%collapse == i .and. sysc(i)%status /= sys_group) &
            call w%select_system_tree(i)
      end if
      forceremap = .true.

    end subroutine collapse_system

  end subroutine draw_tree

  ! Update the table rows by building a new row index array
  ! (iord). Only the systems that are not empty are pointed by
  ! iord. This routine is used when the systems change.
  module subroutine remap_tree(w)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sysc, nsys, sys_empty
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
    w%timelast_tree_update = glfwGetTime()
    if (n == 0) deallocate(w%iord)

  end subroutine remap_tree

  ! Sort the table row order by column cid and in direction dir
  ! (ascending=1, descending=2). Modifies the w%iord.
  module subroutine sort_tree(w)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sysc, sys_init, sys_empty
    use tools, only: mergesort
    use tools_math, only: invert_permutation
    use tools_io, only: ferror, faterr, string
    use types, only: vstring
    class(window), intent(inout) :: w

    integer :: i, n, nvalid, nnovalid
    integer, allocatable :: ival(:), iperm(:), ivalid(:), inovalid(:)
    logical, allocatable :: valid(:)
    real*8, allocatable :: rval(:)
    type(vstring), allocatable :: sval(:)
    logical :: doit

    ! update the time
    w%timelast_sort = glfwGetTime()

    ! check if we have something to sort
    if (.not.allocated(w%iord)) return

    ! initialize the identity permutation
    n = size(w%iord,1)
    allocate(iperm(n),valid(n))
    do i = 1, n
       iperm(i) = i
    end do
    valid = .true.

    ! different types, different sorts
    if (w%sortcid == ic_tree_id.or.w%sortcid == ic_tree_nneq.or.w%sortcid == ic_tree_ncel.or.&
       w%sortcid == ic_tree_nmol.or.w%sortcid == ic_tree_zprime) then
       ! sort by integer
       allocate(ival(n))
       do i = 1, n
          if (w%sortcid == ic_tree_id) then
             ival(i) = w%iord(i)
          elseif (sysc(w%iord(i))%status == sys_init) then
             if (w%sortcid == ic_tree_nneq) then
                ival(i) = sys(w%iord(i))%c%nneq
             elseif (w%sortcid == ic_tree_ncel) then
                ival(i) = sys(w%iord(i))%c%ncel
             elseif (w%sortcid == ic_tree_nmol) then
                ival(i) = sys(w%iord(i))%c%nmol
             elseif (w%sortcid == ic_tree_zprime) then
                if (sys(w%iord(i))%c%spgavail .and. .not.sys(w%iord(i))%c%ismolecule) then
                   ival(i) = nint(real(sys(w%iord(i))%c%nmol,8) / real(sys(w%iord(i))%c%neqv,8))
                else
                   ival(i) = huge(1)
                   valid(i) = .false.
                end if
             end if
          else
             ival(i) = huge(1)
             valid(i) = .false.
          end if
       end do
       call mergesort(ival,iperm,1,n)
       deallocate(ival)
    elseif (w%sortcid == ic_tree_v.or.w%sortcid == ic_tree_a.or.w%sortcid == ic_tree_b.or.&
       w%sortcid == ic_tree_c.or.w%sortcid == ic_tree_alpha.or.w%sortcid == ic_tree_beta.or.&
       w%sortcid == ic_tree_gamma.or.w%sortcid == ic_tree_vmol.or.&
       w%sortcid == ic_tree_e.or.w%sortcid == ic_tree_emol.or.w%sortcid == ic_tree_p) then
       ! sort by real
       allocate(rval(n))
       do i = 1, n
          doit = sysc(w%iord(i))%status == sys_init
          if (doit) doit = (.not.sys(w%iord(i))%c%ismolecule)
          if (w%sortcid == ic_tree_v .and. doit) then
             rval(i) = sys(w%iord(i))%c%omega
          elseif (w%sortcid == ic_tree_a .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(1)
          elseif (w%sortcid == ic_tree_b .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(2)
          elseif (w%sortcid == ic_tree_c .and. doit) then
             rval(i) = sys(w%iord(i))%c%aa(3)
          elseif (w%sortcid == ic_tree_alpha .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(1)
          elseif (w%sortcid == ic_tree_beta .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(2)
          elseif (w%sortcid == ic_tree_gamma .and. doit) then
             rval(i) = sys(w%iord(i))%c%bb(3)
          elseif (w%sortcid == ic_tree_vmol .and. doit) then
             rval(i) = sys(w%iord(i))%c%omega / sys(w%iord(i))%c%nmol
          elseif (w%sortcid == ic_tree_e .and. sysc(w%iord(i))%seed%energy /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%energy
          elseif (w%sortcid == ic_tree_emol .and. sysc(w%iord(i))%status == sys_init .and.&
             sysc(w%iord(i))%seed%energy /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%energy / sys(w%iord(i))%c%nmol
          elseif (w%sortcid == ic_tree_p .and. sysc(w%iord(i))%seed%pressure /= huge(1d0)) then
             rval(i) = sysc(w%iord(i))%seed%pressure
          else
             rval(i) = huge(1d0)
             valid(i) = .false.
          end if
       end do
       call mergesort(rval,iperm,1,n)
       deallocate(rval)
    elseif (w%sortcid == ic_tree_name .or. w%sortcid == ic_tree_spg) then
       ! sort by string
       allocate(sval(n))
       do i = 1, n
          if (w%sortcid == ic_tree_name .and. sysc(w%iord(i))%status /= sys_empty.and.&
             .not.sysc(w%iord(i))%hidden) then
             sval(i)%s = trim(sysc(w%iord(i))%seed%name)
          else
             doit = (w%sortcid == ic_tree_spg) .and. (sysc(w%iord(i))%status == sys_init)
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
       call ferror('sort_tree','column sorting not implemented for column ' // string(w%sortcid),faterr)
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
       if (w%sortdir == 2) then
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

  end subroutine sort_tree

  !> Reassign the currently selected system
  module subroutine reassign_tree(w,cfilter)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sysc, sys_empty, sys_init, kill_initialization_thread, remove_system
    class(window), intent(inout) :: w
    type(c_ptr), intent(inout) :: cfilter

    integer :: idx, i, j, jsel, newsel
    character(kind=c_char,len=:), allocatable, target :: str

    ! this routine only works if the selected tree system is empty
    idx = w%isys
    if (sysc(idx)%status > sys_empty) then
       w%timelast_assign = glfwGetTime()
       return
    end if
    if (.not.allocated(w%iord)) then
       call w%select_system_tree(1)
       return
    end if

    ! find jsel, the currently selected system in the table
    jsel = 0
    do j = 1, size(w%iord,1)
       if (w%iord(j) == idx) then
          jsel = j
          exit
       end if
    end do
    if (jsel > 0) then
       ! search for the next initialized system
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
          ! not found: search backwards for an initialized system
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
          ! well, then set it to 1
          if (newsel == 0) newsel = 1
       end if
    else
       newsel = 1
    end if
    call w%select_system_tree(newsel)

  end subroutine reassign_tree

  !> Select system idx from the tree
  module subroutine select_system_tree(w,idx)
    use interfaces_glfw, only: glfwGetTime
    class(window), intent(inout) :: w
    integer, intent(in) :: idx

    w%isys = idx
    ! the anchor follows the current system, so that a shift-click works
    ! even when the selection was made by another window (or at startup)
    w%tree_selanchor = idx
    w%timelast_assign = glfwGetTime()

  end subroutine select_system_tree

  !> Number of systems currently selected in the tree (group headers,
  !> which are not systems, do not count).
  module function tree_nselected()
    use systems, only: nsys, sysc, sys_empty, sys_group
    integer :: tree_nselected

    tree_nselected = count(sysc(1:nsys)%tselected .and. sysc(1:nsys)%status /= sys_empty .and.&
       sysc(1:nsys)%status /= sys_group)

  end function tree_nselected

  !> Build the list of systems in the tree, in tree order, and return it
  !> in idx. If onlyselected, include only the systems that are selected
  !> in the tree. Group headers are not systems and are never included,
  !> but the members of a collapsed group are: they are not rows of their
  !> own in the tree (iord skips them) and they would otherwise be lost.
  module subroutine tree_system_list(idx,onlyselected)
    use systems, only: nsys, sysc, sys_empty, sys_group
    use types, only: realloc
    integer, allocatable, intent(inout) :: idx(:)
    logical, intent(in) :: onlyselected

    integer :: i, j, k, n

    if (allocated(idx)) deallocate(idx)
    allocate(idx(nsys))
    n = 0
    if (allocated(win(iwin_tree)%iord)) then
       do j = 1, size(win(iwin_tree)%iord,1)
          i = win(iwin_tree)%iord(j)
          if (i < 1 .or. i > nsys) cycle
          if (sysc(i)%status == sys_empty) cycle
          if (sysc(i)%status /= sys_group) then
             if (.not.(onlyselected .and. .not.sysc(i)%tselected)) then
                n = n + 1
                idx(n) = i
             end if
          end if
          ! the systems this row heads, if it is a collapsed master
          if (sysc(i)%collapse /= -1) cycle
          do k = 1, nsys
             if (sysc(k)%status == sys_empty .or. sysc(k)%status == sys_group) cycle
             if (sysc(k)%collapse /= i .or. .not.sysc(k)%hidden) cycle
             if (onlyselected .and. .not.sysc(k)%tselected) cycle
             n = n + 1
             idx(n) = k
          end do
       end do
    end if
    call realloc(idx,n)

  end subroutine tree_system_list

  !xx! private procedures

  ! Clear the tree multi-selection.
  subroutine tree_select_none()
    use systems, only: nsys, sysc

    sysc(1:nsys)%tselected = .false.

  end subroutine tree_select_none

  ! Return the string for the tooltip shown by the tree window,
  ! corresponding to system i.
  subroutine tree_system_tooltip(i)
    use utils, only: iw_text
    use crystalmod, only: pointgroup_info, holo_string, iperiod_3d_crystal,&
       iperiod_3d_layered, iperiod_3d_chain, iperiod_3d_molecular,&
       iperiod_2d, iperiod_1d, iperiod_0d, iperiod_mol_single,&
       iperiod_mol_cluster
    use systems, only: sys, sysc, sys_init, ok_system, sys_empty, sys_group, group_nmembers
    use gui_main, only: ColorTableCellBg
    use tools_io, only: string
    use param, only: bohrtoa, maxzat, atmass, pcamu, bohrtocm
    integer, intent(in) :: i

    character(kind=c_char,len=:), allocatable, target :: str
    real*8, allocatable :: nis(:)
    integer :: k, iz, nsysgroup
    real*8 :: mass, dens
    real(c_float) :: rgba(4)
    real*8 :: nelec
    character(len=3) :: schpg
    integer :: holo, laue
    integer :: izp0

    str = ""

    ! a group header is not a system: say what it holds and where it came from
    if (sysc(i)%status == sys_group) then
       call igBeginTooltip()
       call iw_text(trim(sysc(i)%group_label) // " from " // trim(sysc(i)%group_parent),&
          highlight=.true.)
       nsysgroup = group_nmembers(i)
       call iw_text(string(nsysgroup) // " systems in this group")
       call igEndTooltip()
       return
    end if
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
             izp0 = nint(real(sys(i)%c%nmol,8) / real(sys(i)%c%neqv,8))
             str = str // " and Z'=" // string(izp0)
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

       ! electrons, molar mass (weighted by site occupancy)
       call sys(i)%c%composition(nis)
       nelec = 0d0
       mass = 0d0
       do k = 1, sys(i)%c%nspc
          iz = sys(i)%c%spc(k)%z
          if (iz >= maxzat .or. iz <= 0) cycle
          nelec = nelec + iz * nis(k)
          mass = mass + atmass(iz) * nis(k)
       end do
       call iw_text("Electrons: ",highlight=.true.)
       if (abs(nelec - anint(nelec)) < 1d-6) then
          call iw_text(string(nint(nelec)),sameline_nospace=.true.)
       else
          call iw_text(string(nelec,'f',decimal=4),sameline_nospace=.true.)
       end if
       call iw_text(" Mass (amu): ",highlight=.true.,sameline=.true.)
       call iw_text(string(mass,'f',decimal=3),sameline_nospace=.true.)

       ! empirical formula (weighted by site occupancy)
       call iw_text("Formula: ",highlight=.true.)
       str = sys(i)%c%formula_string(.false.)
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

       ! partial occupancies
       if (sys(i)%c%haveocc) then
          call iw_text("")
          call iw_text("Partial site occupancies",rgba=rgba_partialocc)
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
    use systems, only: sys, sys_init, ok_system
    use fieldmod, only: field, type_uninit, type_promol, type_grid, type_wien,&
       type_elk, type_pi, type_wfn, type_dftb, type_promol_frag, type_ghost
    use wfn_private, only: molden_type_psi4, molden_type_orca, molden_type_adf_sto,&
       wfn_rhf, wfn_uhf, wfn_frac
    use grid3mod, only: mode_nearest, mode_trilinear, mode_trispline, mode_tricubic,&
       mode_smr
    use tools_io, only: string
    integer, intent(in) :: si, fj

    character(kind=c_char,len=:), allocatable, target :: str, aux
    integer :: i, nal, nmoa, nmob, nocca, noccb
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
       call iw_text(string(f%grid%cell_integral(f%c%omega),'f',decimal=8),sameline_nospace=.true.)
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
          nmoa = f%wfn%get_mo_nspin(1,nocca)
          nmob = f%wfn%get_mo_nspin(2,noccb)
          call iw_text("Number of MOs (total): ",highlight=.true.)
          call iw_text(string(f%wfn%nmoall) //&
             " (alpha=" // string(nmoa) // ",beta=" // string(nmob) // ")",sameline_nospace=.true.)
          call iw_text("Number of MOs (occupied): ",highlight=.true.)
          call iw_text(string(f%wfn%nmoocc) // " (alpha=" // string(nocca) //&
             ",beta=" // string(noccb) // ")",sameline_nospace=.true.)
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
    if (allocated(f%zpsp)) then
       ! f%zpsp is indexed by species (1:nspc), not by atomic number
       if (any(f%zpsp > 0)) then
          call iw_text("Core charges (ZPSP): ",highlight=.true.)
          str = ""
          do i = 1, min(size(f%zpsp,1),sys(si)%c%nspc)
             if (f%zpsp(i) > 0) then
                str = str // trim(sys(si)%c%spc(i)%name) // "(" // string(f%zpsp(i)) // ") "
             end if
          end do
          call iw_text(str,sameline_nospace=.true.)
       end if
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
