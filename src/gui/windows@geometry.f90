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

! Routines for the view/edit geometry window.
submodule (windows) geometry
  use interfaces_cimgui
  implicit none
contains

  !> Update tasks for the view/edit geometry window. This is run right
  !> before the window is created and drawn.
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
    use representations, only: reptype_atoms
    use windows, only: iwin_view, iwin_tree
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS, BIND_EDITGEOM_REMOVE
    use systems, only: nsys, sysc, sys, sys_init, ok_system, reread_system_from_file
    use gui_main, only: g, ColorHighlightScene, ColorHighlightSelectScene
    use utils, only: iw_text, iw_tooltip, iw_calcwidth, iw_button, iw_calcheight, iw_calcwidth,&
       iw_combo_simple, iw_highlight_selectable, iw_coloredit
    use global, only: dunit0, iunit_ang
    use tools_io, only: string, nameguess, ioj_center
    class(window), intent(inout), target :: w

    logical :: domol, dowyc, doidx, havesel, removehighlight
    logical :: doquit, clicked, forcesort
    integer :: ihighlight, iclicked, nhigh, dec, icolsort(0:9)
    logical(c_bool) :: is_selected, redo_highlights
    integer(c_int) :: atompreflags, flags, ntype, ncol, ndigit, ndigitm, ndigitidx, color
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2, suffix
    character(len=:), allocatable :: name
    integer, allocatable :: ihigh(:)
    real(c_float), allocatable :: irgba(:,:)
    type(ImVec2) :: szavail, szero, sz0
    real(c_float) :: combowidth, rgb(3)
    integer :: ii, i, j, isys, icol, ispc, iz, iview
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    logical :: havergb, ldum, ok
    real*8 :: x0(3)
    type(ImVec4) :: col4
    type(c_ptr) :: ptrc
    type(ImGuiTableSortSpecs), pointer :: sortspecs
    type(ImGuiTableColumnSortSpecs), pointer :: colspecs

    logical, save :: ttshown = .false. ! tooltip flag

    integer, parameter :: ic_id = 0
    integer, parameter :: ic_atom = 1
    integer, parameter :: ic_zat = 2
    integer, parameter :: ic_mol = 3
    integer, parameter :: ic_mul = 4
    integer, parameter :: ic_wyc = 5
    integer, parameter :: ic_idx = 6
    integer, parameter :: ic_x = 7
    integer, parameter :: ic_y = 8
    integer, parameter :: ic_z = 9

    ! initialize
    ihighlight = 0
    iclicked = 0
    doquit = .false.
    szero%x = 0
    szero%y = 0
    redo_highlights = .false.
    removehighlight = .false.
    forcesort = .false.

    ! first pass
    if (w%firstpass) then
       w%geometry_atomtype = 1
       w%tied_to_tree = (w%isys == win(iwin_tree)%tree_selected)
       if (allocated(w%geometry_selected)) deallocate(w%geometry_selected)
       if (allocated(w%geometry_rgba)) deallocate(w%geometry_rgba)
       w%geometry_select_rgba = ColorHighlightSelectScene
       w%geometry_sortcid = 0
       w%geometry_sortdir = 1
       if (allocated(w%iord)) deallocate(w%iord)
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

    ! force sort if the system has been rebonded or changed geometry
    if (w%timelast_geometry_sort < sysc(isys)%timelastchange_rebond) then
       if (allocated(w%iord)) deallocate(w%iord)
       forcesort = .true.
    end if

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

    !! line of global buttons
    ! restore, only if system is independent or master
    if (iw_button("Restore")) &
       call reread_system_from_file(isys)
    call iw_tooltip("Restore the system to the original geometry it had when it was first opened",ttshown)

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
             call iw_combo_simple("Types##atomtypeselectgeom","Species"//c_null_char//&
                "Symmetry unique" //c_null_char//"Cell"//c_null_char//c_null_char,&
                w%geometry_atomtype)
          else
             call iw_combo_simple("Types##atomtypeselectgeom","Species"//c_null_char//"Atoms"//c_null_char//&
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

          ! if the order array is not allocated or if its size is wrong, force a sort
          if (.not.allocated(w%iord)) then
             forcesort = .true.
          elseif (size(w%iord,1) /= ntype) then
             forcesort = .true.
             if (allocated(w%iord)) deallocate(w%iord)
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
          if (w%geometry_atomtype > 0) ncol = ncol + 3 ! coordinates

          ! atom style table, for atoms
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_Resizable)
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          flags = ior(flags,ImGuiTableFlags_Sortable)
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
             icolsort(icol) = ic_id

             icol = icol + 1
             str2 = "Atom" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
             icolsort(icol) = ic_atom

             icol = icol + 1
             str2 = "Z " // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
             icolsort(icol) = ic_zat

             if (domol) then
                icol = icol + 1
                str2 = "mol" // c_null_char
                flags = ImGuiTableColumnFlags_None
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
                icolsort(icol) = ic_mol
             end if

             if (dowyc) then
                icol = icol + 1
                if (sys(isys)%c%havesym > 0 .and. sys(isys)%c%spgavail) then
                   str2 = "Wyc" // c_null_char
                   icolsort(icol) = ic_wyc
                else
                   str2 = "Mul" // c_null_char
                   icolsort(icol) = ic_mul
                end if
                flags = ImGuiTableColumnFlags_None
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
             end if

             if (doidx) then
                icol = icol + 1
                str2 = "idx" // c_null_char
                flags = ImGuiTableColumnFlags_None
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
                icolsort(icol) = ic_idx
             end if

             if (w%geometry_atomtype > 0) then
                icol = icol + 1
                if (sys(isys)%c%ismolecule) then
                   str2 = "x (Å)" // c_null_char
                else
                   str2 = "x" // c_null_char
                end if
                flags = ImGuiTableColumnFlags_WidthStretch
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
                icolsort(icol) = ic_x

                icol = icol + 1
                if (sys(isys)%c%ismolecule) then
                   str2 = "y (Å)" // c_null_char
                else
                   str2 = "y" // c_null_char
                end if
                flags = ImGuiTableColumnFlags_WidthStretch
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
                icolsort(icol) = ic_y

                icol = icol + 1
                if (sys(isys)%c%ismolecule) then
                   str2 = "z (Å)" // c_null_char
                else
                   str2 = "z" // c_null_char
                end if
                flags = ImGuiTableColumnFlags_WidthStretch
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
                icolsort(icol) = ic_z
             end if
             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! fetch the sort specs, sort the data if necessary
             ptrc = igTableGetSortSpecs()
             if (c_associated(ptrc)) then
                call c_f_pointer(ptrc,sortspecs)
                if (c_associated(sortspecs%Specs)) then
                   call c_f_pointer(sortspecs%Specs,colspecs)
                   w%geometry_sortcid = colspecs%ColumnUserID
                   w%geometry_sortdir = colspecs%SortDirection
                   if (sortspecs%SpecsDirty .and. ntype > 1) then
                      forcesort = .true.
                      sortspecs%SpecsDirty = .false.
                   end if
                else
                   w%geometry_sortcid = 0
                   w%geometry_sortdir = 1
                end if
             end if

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! sort
             if (forcesort) call table_sort()

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
                do ii = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                   i = w%iord(ii)
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
                      ok = iw_highlight_selectable("##selectablemoltable" // suffix,clicked=clicked)
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

                   ! cell index
                   if (doidx) then
                      icol = icol + 1
                      ! i is cell index in this case
                      if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%atcel(i)%idx,ndigitidx))
                   end if

                   ! coordinates
                   if (w%geometry_atomtype > 0) then
                      dec = 6
                      if (sys(isys)%c%ismolecule) then
                         x0 = (sys(isys)%c%atcel(i)%r+sys(isys)%c%molx0) * dunit0(iunit_ang)
                         dec = 4
                      elseif (w%geometry_atomtype == 1) then
                         x0 = sys(isys)%c%at(i)%x
                      else
                         x0 = sys(isys)%c%atcel(i)%x
                      endif
                      do j = 1, 3
                         icol = icol + 1
                         if (igTableSetColumnIndex(icol)) then
                            s = string(x0(j),'f',dec+4,dec,ioj_center)
                            call iw_text(s)
                         end if
                      end do
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

       !! symmetry tab !!
       if (.not.sys(isys)%c%ismolecule) then
          str2 = "Symmetry##drawgeometry_symmetrytab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
             call iw_text("blah")
             call igEndTabItem()
          end if
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

    subroutine table_sort()
      use tools, only: mergesort
      use interfaces_glfw, only: glfwGetTime
      use types, only: vstring
      integer :: ii, i, ispc
      integer, allocatable :: ival(:), iperm(:)
      real*8, allocatable :: rval(:)
      type(vstring), allocatable :: sval(:)

      ! update the time
      w%timelast_geometry_sort = glfwGetTime()

      ! reallocate the iord
      if (ntype <= 1) then
         if (allocated(w%iord)) deallocate(w%iord)
         allocate(w%iord(1))
         w%iord = 1
         return
      else
         if (allocated(w%iord)) then
            if (size(w%iord,1) /= ntype) deallocate(w%iord)
         end if
         if (.not.allocated(w%iord)) then
            allocate(w%iord(ntype))
            do i = 1, ntype
               w%iord(i) = i
            end do
         end if
      end if

      ! carry out the sort
      allocate(iperm(ntype))
      do i = 1, ntype
         iperm(i) = i
      end do
      if (icolsort(w%geometry_sortcid)==ic_id .or. icolsort(w%geometry_sortcid)==ic_zat .or.&
         icolsort(w%geometry_sortcid)==ic_mol .or. icolsort(w%geometry_sortcid)==ic_mul .or.&
         icolsort(w%geometry_sortcid)==ic_idx) then
         ! integers
         allocate(ival(ntype))
         do ii = 1, ntype
            i = w%iord(ii)
            if (icolsort(w%geometry_sortcid) == ic_id) then
               ival(ii) = i
            elseif (icolsort(w%geometry_sortcid) == ic_zat) then
               if (w%geometry_atomtype == 0) then ! species
                  ispc = i
               elseif (w%geometry_atomtype == 1) then ! nneq
                  ispc = sys(isys)%c%at(i)%is
               elseif (w%geometry_atomtype == 2) then ! ncel
                  ispc = sys(isys)%c%atcel(i)%is
               end if
               ival(ii) = sys(isys)%c%spc(ispc)%z
            elseif (icolsort(w%geometry_sortcid) == ic_mol) then
               ival(ii) = sys(isys)%c%idatcelmol(1,i)
            elseif (icolsort(w%geometry_sortcid) == ic_mul) then
               ival(ii) = sys(isys)%c%at(i)%mult
            elseif (icolsort(w%geometry_sortcid) == ic_idx) then
               ival(ii) = sys(isys)%c%atcel(i)%idx
            end if
         end do
         call mergesort(ival,iperm,1,ntype)
         deallocate(ival)
      elseif (icolsort(w%geometry_sortcid)==ic_x.or.icolsort(w%geometry_sortcid)==ic_y.or.&
         icolsort(w%geometry_sortcid)==ic_z) then
         ! real
         allocate(rval(ntype))
         do ii = 1, ntype
            i = w%iord(ii)
            if (sys(isys)%c%ismolecule) then
               x0 = (sys(isys)%c%atcel(i)%r+sys(isys)%c%molx0) * dunit0(iunit_ang)
            elseif (w%geometry_atomtype == 1) then
               x0 = sys(isys)%c%at(i)%x
            else
               x0 = sys(isys)%c%atcel(i)%x
            endif

            if (icolsort(w%geometry_sortcid) == ic_x) then
               rval(ii) = x0(1)
            elseif (icolsort(w%geometry_sortcid) == ic_y) then
               rval(ii) = x0(2)
            elseif (icolsort(w%geometry_sortcid) == ic_z) then
               rval(ii) = x0(3)
            end if
         end do
         call mergesort(rval,iperm,1,ntype)
         deallocate(rval)
      else
         ! strings
         allocate(sval(ntype))
         do ii = 1, ntype
            i = w%iord(ii)
            if (icolsort(w%geometry_sortcid) == ic_atom) then
               if (w%geometry_atomtype == 0) then ! species
                  sval(ii)%s = string(sys(isys)%c%spc(i)%name,2)
               elseif (w%geometry_atomtype == 1) then ! nneq
                  sval(ii)%s = trim(sys(isys)%c%at(i)%name)
               elseif (w%geometry_atomtype == 2) then ! ncel
                  sval(ii)%s = trim(sys(isys)%c%at(sys(isys)%c%atcel(i)%idx)%name)
               end if
            elseif (icolsort(w%geometry_sortcid) == ic_wyc) then
               sval(ii)%s = string(sys(isys)%c%at(i)%mult) // sys(isys)%c%at(i)%wyc
            end if
         end do
         call mergesort(sval,iperm,1,ntype)
         deallocate(sval)
      end if

      ! reverse if necessary
      if (w%geometry_sortdir == 2) then
         allocate(ival(ntype))
         do i = 1, ntype
            ival(i) = iperm(ntype-i+1)
         end do
         iperm = ival
      end if

      ! apply the permutation
      w%iord = w%iord(iperm)

    end subroutine table_sort

  end subroutine draw_geometry

end submodule geometry
