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
    use representations, only: reptype_atoms, reptype_symelem, repflavor_symelem
    use crystalmod, only: symop_kind_plane, symop_kind_axis
    use windows, only: iwin_view, iwin_tree
    use interfaces_glfw, only: glfwGetTime
    use crystalmod, only: holo_string, laue_string, pointgroup_info
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS, BIND_EDITSELECT_REMOVE,&
       BIND_EDITSELECT_DESELECT
    use global, only: bondfactor_def, bonddelta_def
    use systems, only: nsys, sysc, sys, sys_init, ok_system, reread_system_from_file,&
       atlisttype_species, atlisttype_nneq, atlisttype_ncel_frac, atlisttype_ncel_bohr,&
       atlisttype_ncel_ang, atlisttype_nmol, celltransform_standard,&
       celltransform_primstd, celltransform_niggli, celltransform_delaunay
    use gui_main, only: g, ColorHighlightScene, ColorHighlightSelectScene, ColorHighlightBondScene,&
       ColorHighlightBondScene2,&
       ColorTableHighlightRow, ColorBlack, ColorWhite, ColorButtonHoverFactor,&
       ColorButtonActiveFactor, lumweights
    use utils, only: iw_text, iw_tooltip, iw_helpermark, iw_calcwidth, iw_button, iw_calcheight, iw_calcwidth,&
       iw_combo_simple, iw_highlight_selectable, iw_coloredit, iw_dragfloat_real8, iw_checkbox,&
       iw_inputtext, iw_periodictable, iw_menuitem, iw_radiobutton, iw_intstepper, iw_inputint,&
       iw_inputint3
    use types, only: realloc, molsymop_rotation, molsymop_plane, molsymop_imp_rotation
    use tools_io, only: string, nameguess, ioj_center, ioj_right, isinteger, isreal
    use param, only: newline, bohrtoa, pi, atmcov0, maxzat0
    class(window), intent(inout), target :: w

    logical :: domol, dowyc, doidx, docoord, havesel, haveexpr, doocc
    logical :: doquit, clicked, forcesort, ch, lch, deselected, chvol, iactive
    integer :: ihighlight, iclicked, iclicked_ini, iclicked_end, nhigh, dec, icolsort(0:17)
    integer :: ihlbond, ihlbtn ! bonds tab: hovered central atom and hovered neighbor button (cell ids)
    integer :: ipickhl ! cell id of the atom awaiting an add-bond pick (0 = none), highlighted
    integer :: ibrm1, ibrm2, lbrm(3) ! bonds tab: deferred bond removal (cell ids + lattice vector)
    integer :: ibord1, ibord2, lbord(3), ibordval ! bonds tab: deferred bond-order change (cell ids + lvec + order)
    integer :: table_hltype
    logical(c_bool) :: is_selected
    integer(c_int) :: atompreflags, flags, ntype, ncol, ndigit, ndigitm, ndigitidx, color
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2, suffix
    character(len=:), allocatable :: smix
    character(kind=c_char,len=:), allocatable, target :: strx, stry, strz, stropt
    character(len=:), allocatable :: name
    integer, allocatable :: ihigh(:)
    real(c_float), allocatable :: irgba(:,:)
    ! per-frame selection state for the aggregate views (species/nneq); not persisted
    integer, allocatable :: rowstate(:) ! 0 = not selected, 1 = partial, 2 = selected
    integer, allocatable :: rowntot(:)  ! number of cell atoms in each row
    real(c_float), allocatable :: rowrgba(:,:) ! highlight color of each row
    integer :: istate ! scratch selection state
    real(c_float) :: rrgba(4) ! scratch highlight color
    type(ImVec2) :: szavail, szero, sz0
    real(c_float) :: combowidth, rgb(3), lum
    integer :: ii, i, j, jj, isys, icol, ispc, iz, izout, iview, im, jm, ncon
    integer :: ord, zi, zj ! bonds tab: bond order and atomic numbers of the two bonded atoms
    integer :: natused_bonds ! bonds tab: number of distinct atomic species present
    integer, allocatable :: iat_bonds(:) ! bonds tab: Z values of distinct species
    character(len=2), allocatable :: name_bonds(:) ! bonds tab: element symbols of distinct species
    logical :: atused_bonds(maxzat0) ! bonds tab: species presence flags
    real*8 :: dbond
    character(len=:), allocatable :: bondglyph, bondword ! bonds tab: bond-type glyph and word
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    logical :: havergb, havergb_, ldum, ok, oksys, ballow
    real*8 :: x0(3), x6(6), xold(3), x6old(6), res, vol, volold, scal
    real*8 :: occval, occold
    real*8 :: stdrot(3,3), stdcom(3), stdext, stdaxlen ! transient std-orientation axes
    integer :: ieuler_drag ! which Euler angle (1/2/3) is being dragged (0 = none)
    real*8 :: rotdir(3), rotlen ! transient rotation-axis direction and half-length
    type(ImVec4) :: col4
    type(c_ptr) :: ptrc
    type(ImGuiTableSortSpecs), pointer :: sortspecs
    type(ImGuiTableColumnSortSpecs), pointer :: colspecs
    character(len=3) :: schpg
    integer :: holo, laue
    ! symmetry tab
    character(len=:), allocatable :: str_eps ! buffer for the symmetry-tolerance input
    character(len=:), allocatable :: saxc, saxx ! axis strings (crystallographic / Cartesian)
    real*8 :: raxc(3), raxx(3) ! rotation axis in crystallographic / Cartesian coordinates
    integer :: neqv, ncv ! number of operations / centering vectors
    integer :: ihl_symop ! operations table: hovered operation (0 = none)
    integer :: ioptype ! molecular operation type (molsymop_*)

    ! actions at the end of the window draw
    integer :: iaction, iaction_i1, iaction_i2
    real*8 :: iaction_x(3), iaction_x6(6), iaction_m(3,3)
    logical :: iaction_l
    character(len=:), allocatable :: iaction_str
    integer, parameter :: iaction_set_attype_name = 0
    integer, parameter :: iaction_set_atomic_number = 1
    integer, parameter :: iaction_add_species = 2
    integer, parameter :: iaction_set_attype_species = 3
    integer, parameter :: iaction_set_atom_position = 4
    integer, parameter :: iaction_add_species_change_atom = 5
    integer, parameter :: iaction_add_atom = 6
    integer, parameter :: iaction_edit_highlighted = 7
    integer, parameter :: iaction_reorder_highlighted = 8
    integer, parameter :: iaction_swap_atom_ids = 9
    integer, parameter :: iaction_change_cell = 10
    integer, parameter :: iaction_transform_cell = 11
    integer, parameter :: iaction_transform_matrix = 12
    integer, parameter :: iaction_swap_mol_ids = 13
    integer, parameter :: iaction_set_molecule_position = 14
    integer, parameter :: iaction_set_molecule_rotation = 15
    integer, parameter :: iaction_remove_molecules = 16
    integer, parameter :: iaction_reorder_molecules = 17
    integer, parameter :: iaction_restore = 18
    integer, parameter :: iaction_sym_recalc = 19
    integer, parameter :: iaction_sym_clear = 20
    integer, parameter :: iaction_sym_refine = 21
    integer, parameter :: iaction_sym_wholemols = 22
    integer, parameter :: iaction_sym_analyze = 23
    integer, parameter :: iaction_set_atom_occupancy = 24
    integer, parameter :: iaction_sym_deleteops = 25

    ! edit actions on highglighted atoms
    integer, parameter :: edit_remove = 1
    integer, parameter :: edit_merge = 2
    integer, parameter :: edit_duplicate = 3

    ! table column IDs
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
    integer, parameter :: ic_nat = 10
    integer, parameter :: ic_expr = 11
    integer, parameter :: ic_ea = 12 ! Euler angle alpha (molecules tab)
    integer, parameter :: ic_eb = 13 ! Euler angle beta
    integer, parameter :: ic_eg = 14 ! Euler angle gamma
    integer, parameter :: ic_sym = 15 ! point-group symbol (molecules tab)
    integer, parameter :: ic_occ = 16 ! site occupancy

    ! allowed atom list types in tables
    integer, parameter :: atlisttype_allowed(4) = (/atlisttype_nneq,&
       atlisttype_ncel_frac,atlisttype_ncel_ang,atlisttype_ncel_bohr/)

    ! allowed coordinate types for the molecular center of mass combo
    integer, parameter :: atlisttype_coord_allowed(3) = (/atlisttype_ncel_frac,&
       atlisttype_ncel_bohr,atlisttype_ncel_ang/)

    ! threshold for resetting atom position
    real*8, parameter :: epsmoved = 1d-8

    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    ihighlight = 0
    ihl_symop = 0
    ihlbond = 0
    ihlbtn = 0
    ibrm1 = 0
    ibord1 = 0
    ieuler_drag = 0
    iclicked = 0
    doquit = .false.
    szero%x = 0
    szero%y = 0
    forcesort = .false.
    iaction = -1
    table_hltype = w%geometry_atomtype

    ! first pass
    if (w%firstpass) then
       w%tied_to_tree = (w%isys == win(iwin_tree)%tree_selected)
       call clear_highlights_table()
       call reset_sort()
       w%geometry_select_rgba = ColorHighlightSelectScene
       w%lastselected = 0
       w%tabselected = ""
       w%geometry_expression = ""
       w%geometry_expression_ok = .false.
       w%geometry_expr_error = ""
       w%geometry_atomtype = 1
       w%geometry_moltype = atlisttype_ncel_frac
       w%geometry_forcewyc = .true.
       w%geometry_sortcid = 0
       w%geometry_sortdir = 1
       w%geometry_input_coord = 0d0
       w%geometry_input_species = 1
       w%geometry_addbond_iat = 0
       w%geometry_addbond_iview = 0
       w%geometry_addbond_time = 0d0
       w%geometry_cell_simple = .true.
       w%geometry_cell_nrep = 1_c_int
       w%geometry_cell_intmat = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
       w%geometry_cell_cen = 1
       w%geometry_cell_origin = 0d0
       w%geometry_cell_inice = 10
       if (allocated(w%geometry_cell_nice_rmax)) deallocate(w%geometry_cell_nice_rmax)
       if (allocated(w%geometry_cell_nice_mmax)) deallocate(w%geometry_cell_nice_mmax)
       call clear_sym_cache()
    end if

    ! if tied to tree, update the isys
    if (w%tied_to_tree .and. (w%isys /= win(iwin_tree)%tree_selected)) then
       call change_system(win(iwin_tree)%tree_selected)
    end if

    ! check if the system still exists
    if (.not.ok_system(w%isys,sys_init)) then
       ! this dialog does not make sense anymore, close it and exit
       call w%end()
       return
    end if
    isys = w%isys

    ! handle a pending add-bond pick commanded to a view window
    ipickhl = 0
    if (w%geometry_addbond_iat > 0) then
       iview = w%geometry_addbond_iview
       oksys = (iview >= 1 .and. iview <= nwin)
       if (oksys) oksys = win(iview)%isinit .and. win(iview)%isopen .and. win(iview)%type == wintype_view
       ok = oksys
       if (ok) ok = win(iview)%view_selected == isys .and. win(iview)%vmdata%owner == w%id .and.&
          sysc(isys)%timelastchange_geometry < w%geometry_addbond_time
       if (.not.ok) then
          ! the view is gone, shows another system, another window took over
          ! the pick, or the geometry changed (stale cell-atom ids): cancel
          if (oksys) then
             if (win(iview)%vmdata%owner == w%id) then
                win(iview)%viewmode = vm_navigate
                win(iview)%viewmode_transient = .false.
             end if
          end if
          w%geometry_addbond_iat = 0
          w%geometry_addbond_iview = 0
       elseif (win(iview)%viewmode >= 0) then
          ! the pick finished: add a single bond if a valid atom was clicked
          ! (self-bonds and duplicates are rejected by add_bond)
          if (win(iview)%vmdata%idx(1) > 0) &
             call sysc(isys)%add_bond(w%geometry_addbond_iat,win(iview)%vmdata%idx(1),&
                win(iview)%vmdata%idx(2:4),1)
          win(iview)%vmdata%idx = 0
          w%geometry_addbond_iat = 0
          w%geometry_addbond_iview = 0
       else
          ! the pick is in progress: highlight the atom receiving the bond
          ! (applied in the hover-highlight section at the end of the draw)
          ipickhl = w%geometry_addbond_iat
       end if
    end if

    ! set the initial atomtype
    if (w%firstpass) then
       if (sys(isys)%c%ismolecule) then
          w%geometry_atomtype = atlisttype_ncel_ang
       else
          w%geometry_atomtype = atlisttype_ncel_frac
       end if
    end if

    ! force sort if the system has been rebonded or changed geometry
    if (w%timelast_geometry_sort < sysc(isys)%timelastchange_rebond) then
       call reset_sort()
    end if
    if (w%timelast_geometry_clearhighlights < sysc(isys)%timelastchange_geometry) then
       ! the geometry changed: cell-atom indices are stale, so clear the
       ! system selection and force the window cache to be rebuilt
       call sysc(isys)%highlight_clear(.false.)
       call clear_highlights_table()
       ! the symmetry-vs-epsilon analysis is now stale
       call clear_sym_cache()
    end if

    ! system combo
    if (w%firstpass) then
       atompreflags = ImGuiTabItemFlags_SetSelected
    else
       atompreflags = ImGuiTabItemFlags_None
    end if
    call igAlignTextToFramePadding()
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
    call iw_tooltip("View and edit the geometry in this system",ttshown)

    !! line of global buttons
    ! restore, only if system is independent or master
    if (iw_button("Restore",danger=.true.)) iaction = iaction_restore
    call iw_tooltip("Restore the system to the original geometry it had when it was first opened",ttshown)
    ldum = iw_checkbox("Keep Bonding",w%geometry_keepbonding,sameline=.true.)
    call iw_tooltip("Preserve bonding when the system geometry changes",ttshown)

    ! show the tabs
    str1 = "##drawgeometry_tabbar" // c_null_char
    flags = ImGuiTabBarFlags_None
    call igBeginGroup()
    if (igBeginTabBar(c_loc(str1),flags)) then
       !! species tab !!
       str2 = "Species##drawgeometry_speciestab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
          ! check if the tab changed
          call check_changed_tab("species")

          w%geometry_atomtype = atlisttype_species
          table_hltype = atlisttype_species
          ntype = sys(isys)%c%nspc

          ! build the per-row selection state (aggregate view: needs one pass)
          call build_aggregate_state()

          ! if the order array is not allocated or if its size is wrong, force a sort
          if (.not.allocated(w%iord)) then
             call reset_sort()
          elseif (size(w%iord,1) /= ntype) then
             call reset_sort()
          end if

          ! number of columns
          ncol = 4

          ! atom style table, for atoms
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_Resizable)
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          flags = ior(flags,ImGuiTableFlags_Sortable)
          str1="##tablespctyles_" // string(isys) // c_null_char
          call igGetContentRegionAvail(sz0)
          sz0%x = 0
          sz0%y = sz0%y - iw_calcheight(4,0,.true.)
          if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
             icol = -1

             ! header setup
             icol = icol + 1
             str2 = "Id" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_id

             icol = icol + 1
             str2 = "Atom" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_atom

             icol = icol + 1
             str2 = "Z " // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_zat

             icol = icol + 1
             str2 = "Num. Atoms" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_nat

             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! fetch the sort specs, sort the data if necessary
             call fetch_sort_specs()

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! sort
             if (forcesort) call table_sort()

             ! calculate the number of digits for output
             ndigit = ceiling(log10(ntype+0.1d0))
             ndigitm = 0
             ndigitidx = 0

             ! get the current view, if available
             call get_current_view()

             ! draw the rows
             do ii = 1, ntype
                i = w%iord(ii)
                suffix = "_" // string(i)
                icol = -1

                ! start the table and identify the species and Z
                call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                ispc = sysc(isys)%attype_species(w%geometry_atomtype,i)
                name = sysc(isys)%attype_name(w%geometry_atomtype,i)
                iz = sys(isys)%c%spc(ispc)%z

                ! get the color from the first active atoms representation in the main view
                havergb = color_from_view(w%geometry_atomtype,i,rgb)

                ! background color for the table row (full or partial selection)
                call current_row_state(i,istate,rrgba)
                if (istate == 2) then
                   col4 = ImVec4(rrgba(1),rrgba(2),rrgba(3),rrgba(4))
                   color = igGetColorU32_Vec4(col4)
                   call igTableSetBgColor(ImGuiTableBgTarget_RowBg0, color, icol)
                elseif (istate == 1) then
                   col4 = ImVec4(rrgba(1),rrgba(2),rrgba(3),0.4_c_float*rrgba(4))
                   color = igGetColorU32_Vec4(col4)
                   call igTableSetBgColor(ImGuiTableBgTarget_RowBg0, color, icol)
                end if

                ! id
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   call igAlignTextToFramePadding()
                   call iw_text(string(i,ndigit))
                end if
                ! emit even when column 0 is clipped off-screen
                call process_selectable_clicks()

                ! name
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   if (havergb) then
                      ldum = iw_coloredit("##tablecolorg" // suffix,rgb=rgb,nointeraction=.true.)
                      call igSameLine(0._c_float,-1._c_float)
                   end if
                   if (iw_inputtext("##nametextinput" // string(i),bufsize=11,texta=name,width=max(3,len(name)))) then
                      iaction = iaction_set_attype_name
                      iaction_i1 = i
                      iaction_str = name
                   end if
                end if

                ! Z
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   ldum = iw_button(string(iz,3) // "##Z" // string(i),popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
                   if (ok) then
                      izout = iw_periodictable()
                      if (izout >= 0) then
                         iaction = iaction_set_atomic_number
                         iaction_i1 = i
                         iaction_i2 = izout
                         call igCloseCurrentPopup()
                      end if
                      call igEndPopup()
                   end if
                end if

                ! nat (reuse the count from the aggregate pass)
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   call iw_text(string(rowntot(i)))
                end if
             end do

             ! last table row (new species)
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             if (igTableSetColumnIndex(0)) then
                ldum = iw_button("New",popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
                if (ok) then
                   izout = iw_periodictable()
                   if (izout >= 0) then
                      iaction = iaction_add_species
                      iaction_i1 = izout
                      call igCloseCurrentPopup()
                   end if
                   call igEndPopup()
                end if
             end if

             ! end the table
             call igEndTable()
          end if

          ! highlight/selection row
          call draw_highlight_buttons()

          ! edit row
          call draw_edit_buttons()

          call igEndTabItem()
       end if

       !! atoms tab !!
       str2 = "Atoms##drawgeometry_atomstab" // c_null_char
       flags = atompreflags
       if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
          ! check if the tab changed
          call check_changed_tab("atoms")

          ! group atom types
          if (sysc(isys)%attype_combo_simple("Types##atomtypeselectgeom",w%geometry_atomtype,atlisttype_allowed)) then
             call reset_sort()
          end if
          call iw_tooltip("Group atoms by these categories",ttshown)
          table_hltype = w%geometry_atomtype
          ntype = sysc(isys)%attype_number(w%geometry_atomtype)

          if (w%geometry_atomtype == atlisttype_nneq) then
             ldum = iw_checkbox("Keep Symm.",w%geometry_forcewyc,sameline=.true.)
             call iw_tooltip("If checked, changes to the atomic coordinates force the atom to maintain &
                &its current site symmetry (Wyckoff position)",ttshown)
          end if

          if (.not.sys(isys)%c%ismolecule) then
             if (sys(isys)%c%spgavail) then
                call iw_text("  Spg: " // trim(sys(isys)%c%spg%international_symbol),sameline=.true.)
             else
                call iw_text("  Spg: n/a",sameline=.true.)
             end if
          end if

          ! nneq is an aggregate view (needs one pass); cell views compute the
          ! selection state inline per visible row
          if (w%geometry_atomtype == atlisttype_nneq) then
             call build_aggregate_state()
          else
             call clear_aggregate_state()
          end if

          ! if the order array is not allocated or if its size is wrong, force a sort
          if (.not.allocated(w%iord)) then
             call reset_sort()
          elseif (size(w%iord,1) /= ntype) then
             call reset_sort()
          end if

          ! whether to do the molecule column, wyckoff, nneq index
          domol = (w%geometry_atomtype == atlisttype_ncel_frac.or.w%geometry_atomtype == atlisttype_ncel_bohr.or.&
             w%geometry_atomtype == atlisttype_ncel_ang)
          dowyc = (w%geometry_atomtype == atlisttype_nneq .and..not.sys(isys)%c%ismolecule)
          doidx = ((w%geometry_atomtype == atlisttype_ncel_frac.or.w%geometry_atomtype == atlisttype_ncel_bohr.or.&
             w%geometry_atomtype == atlisttype_ncel_ang).and..not.sys(isys)%c%ismolecule)
          docoord = w%geometry_atomtype /= atlisttype_species
          doocc = (w%geometry_atomtype /= atlisttype_species)

          ! number of columns
          ncol = 3
          if (domol) ncol = ncol + 1 ! mol
          if (dowyc) ncol = ncol + 1 ! wyckoff/multiplicity
          if (doidx) ncol = ncol + 1 ! nneq idx
          if (doocc) ncol = ncol + 1 ! occupancy
          if (docoord) ncol = ncol + 3 ! coordinates

          ! whether we have an expression
          haveexpr = (w%geometry_expression_ok .and. len(w%geometry_expression) > 0 .and.&
             len(w%geometry_expr_error) == 0)
          if (haveexpr) ncol = ncol + 1

          ! atom style table, for atoms
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_Resizable)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_Sortable)
          str1="##tableatomstyles_" // string(isys) // "_" // string(w%geometry_atomtype) // c_null_char
          call igGetContentRegionAvail(sz0)
          sz0%x = 0
          sz0%y = sz0%y - iw_calcheight(5,0,.true.)
          if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
             icol = -1

             ! TableSetupColumn(const char* label, ImGuiTableColumnFlags flags = 0, float init_width_or_weight = 0.0f, ImGuiID user_id = 0);
             ! header setup
             icol = icol + 1
             str2 = "Id" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_id

             icol = icol + 1
             str2 = "Atom" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_atom

             icol = icol + 1
             str2 = "Z " // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_zat

             if (domol) then
                icol = icol + 1
                str2 = "mol" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
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
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             end if

             if (doidx) then
                icol = icol + 1
                str2 = "idx" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
                icolsort(icol) = ic_idx
             end if

             if (doocc) then
                icol = icol + 1
                str2 = "Occ." // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
                icolsort(icol) = ic_occ
             end if

             if (docoord) then
                icol = icol + 1
                if (w%geometry_atomtype == atlisttype_ncel_ang .or. w%geometry_atomtype == atlisttype_ncel_bohr) then
                   str2 = "/" // sysc(isys)%attype_coordinates_units(w%geometry_atomtype)
                   strx = "x" // str2 // c_null_char
                   stry = "y" // str2 // c_null_char
                   strz = "z" // str2 // c_null_char
                else
                   strx = "x" // c_null_char
                   stry = "y" // c_null_char
                   strz = "z" // c_null_char
                end if
                call igTableSetupColumn(c_loc(strx),ImGuiTableColumnFlags_None,0.0_c_float,icol)
                icolsort(icol) = ic_x

                icol = icol + 1
                call igTableSetupColumn(c_loc(stry),ImGuiTableColumnFlags_None,0.0_c_float,icol)
                icolsort(icol) = ic_y

                icol = icol + 1
                call igTableSetupColumn(c_loc(strz),ImGuiTableColumnFlags_None,0.0_c_float,icol)
                icolsort(icol) = ic_z
             end if

             if (haveexpr) then
                icol = icol + 1
                str2 = "Expression" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_NoSort,0.0_c_float,icol)
                icolsort(icol) = ic_expr
             end if
             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! fetch the sort specs, sort the data if necessary
             call fetch_sort_specs()

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! sort
             if (forcesort) then
                call table_sort()
             end if

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
             call get_current_view()

             ! draw the rows
             do while(ImGuiListClipper_Step(clipper))
                call c_f_pointer(clipper,clipper_f)
                do ii = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                   i = w%iord(ii)
                   suffix = "_" // string(i)
                   icol = -1

                   ! start the table and identify the species and Z
                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   ispc = sysc(isys)%attype_species(w%geometry_atomtype,i)
                   name = sysc(isys)%attype_name(w%geometry_atomtype,i)
                   iz = sys(isys)%c%spc(ispc)%z

                   ! get the color from the first active atoms representation in the main view
                   havergb = color_from_view(w%geometry_atomtype,i,rgb)

                   ! background color for the table row (full or partial selection)
                   call current_row_state(i,istate,rrgba)
                   if (istate == 2) then
                      col4 = ImVec4(rrgba(1),rrgba(2),rrgba(3),rrgba(4))
                      color = igGetColorU32_Vec4(col4)
                      call igTableSetBgColor(ImGuiTableBgTarget_RowBg0, color, icol)
                   elseif (istate == 1) then
                      col4 = ImVec4(rrgba(1),rrgba(2),rrgba(3),0.4_c_float*rrgba(4))
                      color = igGetColorU32_Vec4(col4)
                      call igTableSetBgColor(ImGuiTableBgTarget_RowBg0, color, icol)
                   end if

                   ! id
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      call igAlignTextToFramePadding()
                      str1 = string(i,ndigit)
                      if (iw_inputtext("##idinput" // suffix,bufsize=11,texta=str1,width=3,notlive=.true.)) then
                         if (isinteger(j,str1)) then
                            if (j > 0 .and. j <= ntype) then
                               iaction = iaction_swap_atom_ids
                               iaction_i1 = i
                               iaction_i2 = j
                            end if
                         end if
                      end if
                   end if
                   ! emit even when column 0 is clipped off-screen
                   call process_selectable_clicks()

                   ! atom name
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      if (havergb) then
                         ldum = iw_coloredit("##tablecolorg" // suffix,rgb=rgb,nointeraction=.true.)
                         call igSameLine(0._c_float,-1._c_float)
                      end if
                      if (iw_inputtext("##nametextinput" // suffix,bufsize=11,texta=name,width=max(3,len(name)))) then
                         iaction = iaction_set_attype_name
                         iaction_i1 = i
                         iaction_str = name
                      end if
                   end if

                   ! Z
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      ldum = iw_button(string(iz,3) // "##Z" // string(i),popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
                      if (ok) then
                         ldum = iw_menuitem("Species",enabled=.false.)
                         call igSeparator()
                         do j = 1, sys(isys)%c%nspc
                            if (iw_menuitem(string(j) // ": " // trim(sys(isys)%c%spc(j)%name))) then
                               iaction = iaction_set_attype_species
                               iaction_i1 = i
                               iaction_i2 = j
                            end if
                         end do
                         call igSeparator()
                         str1 = "New" // c_null_char
                         if (igBeginMenu(c_loc(str1),.true._c_bool)) then
                            izout = iw_periodictable()
                            if (izout >= 0) then
                               iaction = iaction_add_species_change_atom
                               iaction_i1 = izout
                               iaction_i2 = i
                               call igCloseCurrentPopup()
                            end if
                            call igEndMenu()
                         end if
                         call igEndPopup()
                      end if
                   end if

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

                   ! occupancy
                   if (doocc) then
                      icol = icol + 1
                      if (igTableSetColumnIndex(icol)) then
                         occval = sysc(isys)%attype_occupancy(w%geometry_atomtype,i)
                         occold = occval
                         if (iw_dragfloat_real8("##occ" // string(isys) // "_" // string(i),&
                            x1=occval,speed=0.001d0,min=1d-10,max=1d0,decimal=3,notlive=.true.,&
                            flags=ImGuiSliderFlags_AlwaysClamp)) then
                            if (abs(occval-occold) > 1d-10) then
                               iaction = iaction_set_atom_occupancy
                               iaction_i1 = i
                               iaction_x(1) = occval
                            end if
                         end if
                         ! mixed site: list the occupants on hover
                         smix = sysc(isys)%attype_mixed(w%geometry_atomtype,i)
                         if (len_trim(smix) > 0) &
                            call iw_tooltip("Mixed site: " // trim(smix))
                      end if
                   end if

                   ! coordinates
                   if (docoord) then
                      dec = sysc(isys)%attype_coordinates_decimals(w%geometry_atomtype)
                      x0 = sysc(isys)%attype_coordinates(w%geometry_atomtype,i)
                      xold = x0
                      ch = .false.
                      do j = 1, 3
                         icol = icol + 1
                         if (igTableSetColumnIndex(icol)) then
                            s = string(x0(j),'f',dec+4,dec,ioj_center)
                            ch = ch .or. iw_dragfloat_real8("##x" // string(isys) // "_" // string(i) // "_" // string(j),&
                               x1=x0(j),speed=0.001d0,decimal=dec,notlive=.true.)
                         end if
                      end do
                      if (ch .and. any(abs(x0-xold) > epsmoved)) then
                         iaction = iaction_set_atom_position
                         iaction_i1 = i
                         iaction_x = x0
                         iaction_l = w%geometry_forcewyc
                      end if
                   end if

                   ! expression
                   if (haveexpr) then
                      icol = icol + 1
                      if (igTableSetColumnIndex(icol)) then
                         if (w%geometry_expression_ok) then
                            if (w%geometry_atomtype == atlisttype_nneq) then
                               x0 = sys(isys)%c%at(i)%r
                            else
                               x0 = sys(isys)%c%atcel(i)%r
                            end if
                            res = sys(isys)%eval(w%geometry_expression,w%geometry_expr_error,x0)
                            if (len(w%geometry_expr_error) > 0) then
                               w%geometry_expression_ok = .false.
                            else
                               call iw_text(string(res,'f',decimal=8))
                            end if
                         end if
                      end if
                   end if
                end do ! clipper indices
             end do ! clipper step

             ! end the clipper
             call ImGuiListClipper_End(clipper)

             ! last table row (new atom)
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             if (igTableSetColumnIndex(0)) then
                ldum = iw_button("New##newatom",popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
                if (ldum) then
                   w%geometry_input_coord = 0d0
                   w%geometry_input_species = 1
                end if
                if (ok) then
                   call draw_addatom_popup()
                   call igEndPopup()
                end if
             end if

             ! end the table
             call igEndTable()
          end if

          ! highlight/selection row
          call draw_highlight_buttons()

          ! edit row
          call draw_edit_buttons()

          ! expression row
          call igAlignTextToFramePadding()
          call iw_text("Expression",highlight=.true.)

          ! filter text input
          call iw_helpermark("Examples:"//newline//&
             "- '@dnuc:2' = distance to cell atom 2"//newline//&
             "- 'log($0)' = log of the promolecular density"//newline//&
             "- 'abs(@x) < 2 && abs(@y) < 2 && abs(@z) < 2' = atoms in the (-2,2) box"//newline//&
             "Click the Help button for more info.")
          call igSameLine(0._c_float,-1._c_float)
          if (iw_inputtext("##filtertext",bufsize=1023,width=30,texta=w%geometry_expression,&
             notlive=.true.)) then
             w%geometry_expression_ok = .true.
             w%geometry_expr_error = ""
          end if
          call iw_tooltip("Show on the table the result of using this expression at the atomic positions.",ttshown)
          if (iw_button("Help##helpfilter",sameline=.true.)) then
             str2 = "https://aoterodelaroza.github.io/critic2/manual/arithmetics" // c_null_char
             call openLink(c_loc(str2))
          end if
          call iw_tooltip("Open the manual page regarding arithmetic expressions.",ttshown)
          if (iw_button("Clear",sameline=.true.)) then
             w%geometry_expression_ok = .false.
             w%geometry_expression = ""
             w%geometry_expr_error = ""
          end if
          call iw_tooltip("Clear the expression",ttshown)
          if (len_trim(w%geometry_expr_error) > 0) then
             call iw_text("Error: " // trim(w%geometry_expr_error),danger=.true.)
          end if
          call igEndTabItem()
       end if

       !! cell tab !!
       if (.not.sys(isys)%c%ismolecule) then
          str2 = "Cell##drawgeometry_celltab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
             ! check if the tab changed
             call check_changed_tab("cell")

             call iw_text("Lattice Parameters",highlight=.true.)

             ! keep symmetry and space group
             ldum = iw_checkbox("Keep Symmetry",w%geometry_forcewyc)
             call iw_tooltip("If checked, force system to maintain the current space group",ttshown)
             if (sys(isys)%c%spgavail) then
                call pointgroup_info(sys(isys)%c%spg%pointgroup_symbol,schpg,holo,laue)
                call iw_text("  Spg: " // trim(sys(isys)%c%spg%international_symbol),sameline=.true.)
                call iw_text("(" // string(holo_string(holo)) // ")",sameline=.true.)
             else
                call iw_text("  Spg: n/a",sameline=.true.)
             end if

             x6(1:3) = sys(isys)%c%aa
             x6(4:6) = sys(isys)%c%bb
             x6old = x6

             call igAlignTextToFramePadding()
             call iw_text("a/b/c (Å): ")
             ch = .false.
             ch = ch .or. iw_dragfloat_real8("##celllengthsa",x1=x6(1),speed=0.005d0,decimal=6,scale=bohrtoa,&
                notlive=.true.,sameline=.true.)
             ch = ch .or. iw_dragfloat_real8("##celllengthsb",x1=x6(2),speed=0.005d0,decimal=6,scale=bohrtoa,&
                notlive=.true.,sameline=.true.)
             ch = ch .or. iw_dragfloat_real8("##celllengthsc",x1=x6(3),speed=0.005d0,decimal=6,scale=bohrtoa,&
                notlive=.true.,sameline=.true.)

             call igAlignTextToFramePadding()
             call iw_text("α/β/γ (°): ")
             ch = ch .or. iw_dragfloat_real8("##cellangsa",x1=x6(4),speed=0.01d0,decimal=4,sameline=.true.,&
                notlive=.true.)
             ch = ch .or. iw_dragfloat_real8("##cellanbsb",x1=x6(5),speed=0.01d0,decimal=4,sameline=.true.,&
                notlive=.true.)
             ch = ch .or. iw_dragfloat_real8("##cellanbsc",x1=x6(6),speed=0.01d0,decimal=4,sameline=.true.,&
                notlive=.true.)

             ! cell volume
             call igAlignTextToFramePadding()
             call iw_text("Volume (Å³): ")
             vol = sys(isys)%c%omega * bohrtoa**3
             volold = vol
             chvol = iw_dragfloat_real8("##cellvolume",x1=vol,speed=0.5d0,decimal=4,min=1d-6,&
                notlive=.true.,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Cell volume in Å³. Changing it scales the cell isotropically, &
                &compressing or expanding the crystal at constant fractional coordinates",ttshown)

             if (ch .and. any(abs(x6-x6old) > epsmoved)) then
                iaction = iaction_change_cell
                iaction_x6 = x6
                iaction_l = w%geometry_forcewyc
             elseif (chvol .and. abs(vol-volold) > epsmoved .and. vol > 1d-6) then
                ! isotropic scale factor on the lengths to reach the target volume
                scal = (vol / volold)**(1d0/3d0)
                iaction = iaction_change_cell
                iaction_x6(1:3) = sys(isys)%c%aa * scal
                iaction_x6(4:6) = sys(isys)%c%bb
                iaction_l = w%geometry_forcewyc
             end if

             ! standardization and reduction of the cell
             call iw_text("Cell Transformations",highlight=.true.)
             if (iw_button("Conventional##celltransfstd")) then
                iaction = iaction_transform_cell
                iaction_i1 = celltransform_standard
             end if
             call iw_tooltip("Transform to the conventional cell",ttshown)
             if (iw_button("Primitive##celltransfprim",sameline=.true.)) then
                iaction = iaction_transform_cell
                iaction_i1 = celltransform_primstd
             end if
             call iw_tooltip("Transform to the standard primitive cell",ttshown)
             if (iw_button("Niggli##celltransfnig",sameline=.true.)) then
                iaction = iaction_transform_cell
                iaction_i1 = celltransform_niggli
             end if
             call iw_tooltip("Transform to the primitive Niggli cell",ttshown)
             if (iw_button("Delaunay##celltransfdel",sameline=.true.)) then
                iaction = iaction_transform_cell
                iaction_i1 = celltransform_delaunay
             end if
             call iw_tooltip("Transform to the primitive Delaunay cell",ttshown)

             ! cell transformation (simple supercell or general matrix)
             call igAlignTextToFramePadding()
             call iw_text("Supercells",highlight=.true.)
             ldum = iw_radiobutton("Simple",bool=w%geometry_cell_simple,boolval=.true.,sameline=.true.)
             call iw_tooltip("Build a supercell by repeating the cell an integer number of&
                & times along each lattice vector",ttshown)
             ldum = iw_radiobutton("Full",bool=w%geometry_cell_simple,boolval=.false.,sameline=.true.)
             call iw_tooltip("General transformation: each row is a lattice vector of the new&
                & cell, written in crystallographic coordinates of the current cell",ttshown)

             if (w%geometry_cell_simple) then
                ! simple transformation: integer multiples of the a, b, c axes
                call igAlignTextToFramePadding()
                call iw_text("na/nb/nc: ")
                do jm = 1, 3
                   ldum = iw_intstepper("cellnrep" // string(jm),w%geometry_cell_nrep(jm),&
                      minval=1_c_int,sameline=.true.,&
                      tooltip="Number of times the cell is repeated along the a, b, and c lattice vectors")
                end do

                if (iw_button("Reset##cellresettransf",sameline=.true.)) then
                   if (w%geometry_cell_simple) then
                      w%geometry_cell_nrep = 1
                   else
                      w%geometry_cell_intmat = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
                      w%geometry_cell_cen = 1
                   end if
                end if
                call iw_tooltip("Reset the supercell transformation",ttshown)
             else
                ! arbitrary transformation matrix: each row is a lattice vector of the new
                ! cell = integer triplet + (for centered cells) a centering vector
                w%geometry_cell_cen = min(max(w%geometry_cell_cen,1),max(sys(isys)%c%ncv,1))
                if (sys(isys)%c%ncv > 1) then
                   ! build the centering options once
                   stropt = ""
                   do j = 1, sys(isys)%c%ncv
                      stropt = stropt // cell_cen_label(sys(isys)%c%cen(:,j)) // c_null_char
                   end do
                end if
                ! track changes to the vectors/origin to clear a stale error message
                ch = .false.
                do im = 1, 3
                   call igAlignTextToFramePadding()
                   if (im == 1) then
                      call iw_text("a' ")
                   elseif (im == 2) then
                      call iw_text("b' ")
                   else
                      call iw_text("c' ")
                   end if
                   ch = ch .or. iw_inputint3("##cellmat" // string(im),&
                      w%geometry_cell_intmat(:,im),width=3*3,sameline=.true.,notlive=.true.)
                   if (sys(isys)%c%ncv > 1) then
                      call iw_text(" + ",sameline=.true.)
                      call iw_combo_simple("##cellcen" // string(im),stropt,w%geometry_cell_cen(im),&
                         changed=lch,sameline=.true.,startsatone=.true.)
                      ch = ch .or. lch
                      call iw_tooltip("Centering vector added to this lattice vector",ttshown)
                   end if
                end do
                call igAlignTextToFramePadding()
                call iw_text("Origin")
                do jm = 1, 3
                   ch = ch .or. iw_dragfloat_real8("##cellorigin" // string(jm),&
                      x1=w%geometry_cell_origin(jm),speed=0.01d0,decimal=4,sameline=.true.,&
                      notlive=.true.)
                end do
                ! if the user changed the matrix or the origin, clear the error message
                if (ch) w%errmsg = ""
             end if
             ok = iw_button("Apply##celltransfapply",danger=.true.)
             call iw_tooltip("Apply the supercell transformation",ttshown)
             if (ok) then
                iaction = iaction_transform_matrix
                if (w%geometry_cell_simple) then
                   iaction_m = 0d0
                   iaction_m(1,1) = real(w%geometry_cell_nrep(1),8)
                   iaction_m(2,2) = real(w%geometry_cell_nrep(2),8)
                   iaction_m(3,3) = real(w%geometry_cell_nrep(3),8)
                   iaction_x = 0d0
                else
                   do im = 1, 3
                      iaction_m(:,im) = real(w%geometry_cell_intmat(:,im),8) + &
                         sys(isys)%c%cen(:,w%geometry_cell_cen(im))
                   end do
                   iaction_x = w%geometry_cell_origin
                end if
                iaction_l = .false.
             end if
             if (len_trim(w%errmsg) > 0) &
                call iw_text(w%errmsg,danger=.true.,sameline=.true.)

             ! nice supercell search
             call iw_text("Nice supercells",highlight=.true.)
             call iw_helpermark(&
                "For a given n, search for the supercell containing n primitive cells that "//&
                "fits inside the largest possible sphere (with radius rmax). The niceness is defined as the "//&
                "ratio between the radius of the inscribed sphere and the maximum possible "//&
                "ratio, equal to that of a cubic cell. Select the maximum value of n "//&
                "and click search to find all the transformations (expensive for large maximum "//&
                "values). Click on any of the table rows to effect the transformation.")
             call iw_tooltip("Search for the most cube-like supercells up to the given size&
                & (number of times the current cell). Click a row to transform to that supercell.",ttshown)
             ldum = iw_intstepper("cellnicesize",w%geometry_cell_inice,label="Max. size",minval=1_c_int,&
                tooltip="Maximum supercell size (number of times the current cell) to consider in the search")
             if (iw_button("Search##cellnicesearch",sameline=.true.)) then
                w%geometry_cell_inice = max(w%geometry_cell_inice,1_c_int)
                call sysc(isys)%cell_nice_list(int(w%geometry_cell_inice),&
                   w%geometry_cell_nice_rmax,w%geometry_cell_nice_mmax)
             end if
             call iw_tooltip("Search for the nicest shells up to the indicated&
                & maximum size (can be expensive)",ttshown)

             if (allocated(w%geometry_cell_nice_rmax)) then
                flags = ImGuiTableFlags_RowBg
                flags = ior(flags,ImGuiTableFlags_Borders)
                flags = ior(flags,ImGuiTableFlags_ScrollY)
                flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
                str1 = "##cellnicetable" // c_null_char
                sz0%x = iw_calcwidth(2,1) + iw_calcwidth(8,1) + iw_calcwidth(8,1) + iw_calcwidth(27,1) +&
                   g%Style%ScrollbarSize + 4._c_float
                sz0%y = igGetFrameHeight() + 11._c_float * (igGetTextLineHeight() + 2._c_float*g%Style%CellPadding%y)
                if (igBeginTable(c_loc(str1),4,flags,sz0,0._c_float)) then
                   str2 = "n" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,0)
                   str2 = "rmax (Å)" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,1)
                   str2 = "niceness" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,2)
                   str2 = "transformation" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,3)
                   call igTableSetupScrollFreeze(0,1)
                   call igTableHeadersRow()

                   do i = 1, size(w%geometry_cell_nice_rmax,1)
                      if (w%geometry_cell_nice_rmax(i) <= 0d0) cycle
                      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)

                      ! clickable row spanning all columns: apply this supercell
                      if (igTableSetColumnIndex(0)) then
                         str2 = string(i) // "##cellnicerow" // string(i) // c_null_char
                         if (igSelectable_Bool(c_loc(str2),logical(.false.,c_bool),&
                            ImGuiSelectableFlags_SpanAllColumns,szero)) then
                            iaction = iaction_transform_matrix
                            iaction_m = w%geometry_cell_nice_mmax(:,:,i)
                            iaction_x = 0d0
                            iaction_l = .false.
                         end if
                      end if
                      if (igTableSetColumnIndex(1)) &
                         call iw_text(string(w%geometry_cell_nice_rmax(i)*bohrtoa,'f',7,3))
                      if (igTableSetColumnIndex(2)) &
                         call iw_text(string(8d0*w%geometry_cell_nice_rmax(i)**3/(i*sys(isys)%c%omega),'f',7,5))
                      if (igTableSetColumnIndex(3)) then
                         s = ""
                         do im = 1, 3
                            do jm = 1, 3
                               s = s // string(nint(w%geometry_cell_nice_mmax(jm,im,i)),length=3,justify=ioj_right)
                            end do
                         end do
                         call iw_text(s)
                      end if
                   end do
                   call igEndTable()
                end if
             end if

             call igEndTabItem()
          end if
       end if

       !! molecules tab !!
       str2 = "Molecules##drawgeometry_molstab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
          ! check if the tab changed
          call check_changed_tab("molecules")

          ! the selection/highlights in this tab refer to molecules
          table_hltype = atlisttype_nmol

          ! coordinate type for the center of mass
          ldum = sysc(isys)%attype_combo_simple("Center of mass##molcoordselectgeom",&
             w%geometry_moltype,atlisttype_coord_allowed)
          call iw_tooltip("Coordinate system for the molecular center of mass",ttshown)
          call iw_text("(Sym./angles for discrete molecules only)",sameline=.true.)

          ntype = sys(isys)%c%nmol

          ! molecules compute the selection state inline per visible row
          call clear_aggregate_state()

          ! if the order array is not allocated or if its size is wrong, force a sort
          if (.not.allocated(w%iord)) then
             call reset_sort()
          elseif (size(w%iord,1) /= ntype) then
             call reset_sort()
          end if

          ! whether to show the idx (symmetry-equivalence) column
          doidx = (.not.sys(isys)%c%ismolecule .and. allocated(sys(isys)%c%idxmol))

          ! number of columns
          ncol = 2 ! id, nat
          ncol = ncol + 1 ! point-group symbol
          if (doidx) ncol = ncol + 1 ! idx
          ncol = ncol + 3 ! center of mass
          ncol = ncol + 3 ! euler angles

          ! molecule table
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_Resizable)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_Sortable)
          str1="##tablemolecules_" // string(isys) // "_" // string(w%geometry_moltype) // c_null_char
          call igGetContentRegionAvail(sz0)
          sz0%x = 0
          sz0%y = sz0%y - iw_calcheight(5,0,.true.)
          if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
             icol = -1

             ! header setup
             icol = icol + 1
             str2 = "Id" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_id

             icol = icol + 1
             str2 = "nat" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_nat

             icol = icol + 1
             str2 = "Sym." // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_NoSort,0.0_c_float,icol)
             icolsort(icol) = ic_sym

             if (doidx) then
                icol = icol + 1
                str2 = "idx" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
                icolsort(icol) = ic_idx
             end if

             icol = icol + 1
             if (w%geometry_moltype == atlisttype_ncel_ang .or. w%geometry_moltype == atlisttype_ncel_bohr) then
                str2 = "/" // sysc(isys)%attype_coordinates_units(w%geometry_moltype)
                strx = "x" // str2 // c_null_char
                stry = "y" // str2 // c_null_char
                strz = "z" // str2 // c_null_char
             else
                strx = "x" // c_null_char
                stry = "y" // c_null_char
                strz = "z" // c_null_char
             end if
             call igTableSetupColumn(c_loc(strx),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_x
             icol = icol + 1
             call igTableSetupColumn(c_loc(stry),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_y
             icol = icol + 1
             call igTableSetupColumn(c_loc(strz),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_z

             ! Euler angles (ZYZ, degrees) of the standard orientation
             icol = icol + 1
             str2 = "α/°" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_ea
             icol = icol + 1
             str2 = "β/°" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_eb
             icol = icol + 1
             str2 = "γ/°" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,icol)
             icolsort(icol) = ic_eg

             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! fetch the sort specs, sort the data if necessary
             call fetch_sort_specs()

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! sort
             if (forcesort) call table_sort()

             ! start the clipper
             clipper = ImGuiListClipper_ImGuiListClipper()
             call ImGuiListClipper_Begin(clipper,ntype,-1._c_float)

             ! number of digits for the id output
             ndigit = ceiling(log10(ntype+0.1d0))

             ! decimals for the coordinate output
             if (w%geometry_moltype == atlisttype_ncel_frac) then
                dec = 6
             else
                dec = 4
             end if

             ! draw the rows
             do while(ImGuiListClipper_Step(clipper))
                call c_f_pointer(clipper,clipper_f)
                do ii = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                   i = w%iord(ii)
                   suffix = "_" // string(i)
                   icol = -1

                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)

                   ! background color for the table row (full or partial selection)
                   call current_row_state(i,istate,rrgba)
                   if (istate == 2) then
                      col4 = ImVec4(rrgba(1),rrgba(2),rrgba(3),rrgba(4))
                      color = igGetColorU32_Vec4(col4)
                      call igTableSetBgColor(ImGuiTableBgTarget_RowBg0, color, icol)
                   elseif (istate == 1) then
                      col4 = ImVec4(rrgba(1),rrgba(2),rrgba(3),0.4_c_float*rrgba(4))
                      color = igGetColorU32_Vec4(col4)
                      call igTableSetBgColor(ImGuiTableBgTarget_RowBg0, color, icol)
                   end if

                   ! id
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      call igAlignTextToFramePadding()
                      str1 = string(i,ndigit)
                      if (iw_inputtext("##idinput" // suffix,bufsize=11,texta=str1,width=3,notlive=.true.)) then
                         if (isinteger(j,str1)) then
                            if (j > 0 .and. j <= ntype) then
                               iaction = iaction_swap_mol_ids
                               iaction_i1 = i
                               iaction_i2 = j
                            end if
                         end if
                      end if
                   end if
                   ! emit even when column 0 is clipped off-screen
                   call process_selectable_clicks()

                   ! number of atoms
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%mol(i)%nat))

                   ! point-group symbol (only meaningful for discrete fragments)
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      if (sys(isys)%c%mol(i)%discrete) then
                         call iw_text(sys(isys)%c%mol(i)%pgsymbol())
                      else
                         call iw_text("--")
                      end if
                   end if

                   ! idx (symmetry equivalence)
                   if (doidx) then
                      icol = icol + 1
                      if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%idxmol(i)))
                   end if

                   ! center of mass, converted to the chosen coordinate type;
                   ! editable (translation) only for discrete fragments
                   x0 = mol_com_coords(i)
                   if (sys(isys)%c%mol(i)%discrete) then
                      xold = x0
                      ch = .false.
                      do j = 1, 3
                         icol = icol + 1
                         if (igTableSetColumnIndex(icol)) then
                            ch = ch .or. iw_dragfloat_real8("##com" // string(isys) // "_" // string(i) // "_" // string(j),&
                               x1=x0(j),speed=0.001d0,decimal=dec,notlive=.true.)
                            ! keep the molecule highlighted (and its axes shown)
                            ! while its center of mass is being dragged
                            if (logical(igIsItemActive())) ihighlight = i
                         end if
                      end do
                      if (ch .and. any(abs(x0-xold) > epsmoved)) then
                         iaction = iaction_set_molecule_position
                         iaction_i1 = i
                         iaction_x = x0
                      end if
                   else
                      do j = 1, 3
                         icol = icol + 1
                         if (igTableSetColumnIndex(icol)) &
                            call iw_text(string(x0(j),'f',dec+4,dec,ioj_right))
                      end do
                   end if

                   ! Euler angles (ZYZ, degrees) of the standard orientation;
                   ! editable (rigid rotation) for discrete molecules with >1 atom
                   x0 = mol_euler_angles(i)
                   if (sys(isys)%c%mol(i)%discrete .and. sys(isys)%c%mol(i)%nat > 1) then
                      ! while this molecule's Euler angles are being dragged, edit the unwrapped angles.
                      if (w%geometry_euler_drag_mol == i) x0 = w%geometry_euler_drag_val
                      xold = x0
                      ch = .false.
                      iactive = .false.
                      do j = 1, 3
                         icol = icol + 1
                         if (igTableSetColumnIndex(icol)) then
                            ch = ch .or. iw_dragfloat_real8("##euler" // string(isys) // "_" // string(i) // "_" // string(j),&
                               x1=x0(j),speed=0.5d0,decimal=2,notlive=.true.)
                            ! keep the molecule highlighted, axes shown while dragging
                            if (logical(igIsItemActive())) then
                               ihighlight = i
                               ieuler_drag = j
                               iactive = .true.
                            end if
                         end if
                      end do
                      if (ch .and. any(abs(x0-xold) > epsmoved)) then
                         iaction = iaction_set_molecule_rotation
                         iaction_i1 = i
                         iaction_x = x0 * (pi / 180d0)
                      end if
                      ! persist the unwrapped angles for the next frame of the drag
                      if (iactive) then
                         w%geometry_euler_drag_mol = i
                         w%geometry_euler_drag_val = x0
                      end if
                   elseif (.not.sys(isys)%c%mol(i)%discrete) then
                      ! non-discrete fragment: the standard orientation (and its
                      ! Euler angles) are not meaningful, so do not show them
                      do j = 1, 3
                         icol = icol + 1
                         if (igTableSetColumnIndex(icol)) call iw_text("--")
                      end do
                   else
                      do j = 1, 3
                         icol = icol + 1
                         if (igTableSetColumnIndex(icol)) &
                            call iw_text(string(x0(j),'f',6,2,ioj_right))
                      end do
                   end if
                end do
             end do
             call igEndTable()
          end if

          ! highlight/selection row
          call draw_highlight_buttons()

          ! edit row (remove/reorder molecules)
          call draw_mol_edit_buttons()

          call igEndTabItem()
       end if

       !! bonds tab !!
       str2 = "Bonds##drawgeometry_bondstab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
          ! check if the tab changed
          call check_changed_tab("bonds")

          ! blue header
          call iw_text("System Bonds",highlight=.true.)
          call iw_helpermark("This table shows the bonds in the current system. The bonds drawn on the view window &
             &are usually the ones in this list. However, the user can choose to draw different bonding arrangements &
             &using the Bonds tab in the atoms object (Recalculate Bonds). In that case, adding or removing bonds in this table &
             &will not affect what is shown on the view.",sameline=.true.)

          ! legend for the bond-type glyphs used in the Bonded atoms column
          call iw_text("─ single   ═ double   ≡ triple   ○ aromatic   ┄ dashed",sameline=.true.)

          ! number of cell atoms and digits for the Id column
          ntype = sys(isys)%c%ncel
          ndigit = ceiling(log10(ntype+0.1d0))
          table_hltype = atlisttype_ncel_frac

          ! determine which atomic species are present (needed for rebond radii table height)
          atused_bonds = .false.
          do i = 1, sys(isys)%c%nspc
             if (sys(isys)%c%spc(i)%z <= 0) cycle
             atused_bonds(sys(isys)%c%spc(i)%z) = .true.
          end do
          natused_bonds = count(atused_bonds)
          allocate(iat_bonds(natused_bonds),name_bonds(natused_bonds))
          natused_bonds = 0
          do i = 1, maxzat0
             if (atused_bonds(i)) then
                natused_bonds = natused_bonds + 1
                iat_bonds(natused_bonds) = i
                name_bonds(natused_bonds) = nameguess(i,.true.)
             end if
          end do

          ! locate the view holding the atom colors for the bonded-atom buttons
          call get_current_view()

          ! bonds table (height leaves room for rebond controls below)
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_Resizable)
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          flags = ior(flags,ImGuiTableFlags_ScrollX)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          str1 = "##tablebonds_" // string(isys) // c_null_char
          call igGetContentRegionAvail(sz0)
          sz0%x = 0
          sz0%y = sz0%y - iw_calcheight(1,0,.true.) &
             - iw_calcheight(3,4,.false.) &
             - iw_calcheight(min(5,natused_bonds)+1,0,.false.)
          if (igBeginTable(c_loc(str1),3,flags,sz0,0._c_float)) then
             str2 = "Id" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,0)
             str2 = "Atom" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,1)
             str2 = "Bonded atoms" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,2)
             call igTableSetupScrollFreeze(0, 1) ! top row always visible
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! draw the rows (clipped for performance)
             clipper = ImGuiListClipper_ImGuiListClipper()
             call ImGuiListClipper_Begin(clipper,ntype,-1._c_float)
             do while(ImGuiListClipper_Step(clipper))
                call c_f_pointer(clipper,clipper_f)
                do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                   suffix = "_" // string(i)
                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   icol = -1

                   ! id (with a row-spanning selectable to detect hover)
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      call igAlignTextToFramePadding()
                      call iw_text(string(i,ndigit))
                   end if
                   ! emit even when column 0 is clipped off-screen
                   if (iw_highlight_selectable("##bondselect" // suffix)) ihlbond = i

                   ! atom name
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      call igAlignTextToFramePadding()
                      call iw_text(trim(sys(isys)%c%at(sys(isys)%c%atcel(i)%idx)%name))
                   end if

                   ! bonded atoms (one colored button per neighbor)
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      ! "+" button: pick an atom in the view to add a bond to this atom
                      if (iw_button("+##addbond" // suffix,disabled=(iview == 0))) then
                         w%geometry_addbond_iat = i
                         w%geometry_addbond_iview = iview
                         w%geometry_addbond_time = glfwGetTime()
                         call win(iview)%viewmode_set_forced(vm_pick_atom,&
                            "Please pick an atom to bond to atom " // string(i) // "...",w%id)
                      end if
                      call iw_tooltip("Add a bond to this atom: click, then pick an atom in the view window",ttshown)

                      ncon = 0
                      if (allocated(sys(isys)%c%nstar)) ncon = sys(isys)%c%nstar(i)%ncon
                      do j = 1, ncon
                         jj = sys(isys)%c%nstar(i)%idcon(j)
                         name = trim(sys(isys)%c%at(sys(isys)%c%atcel(jj)%idx)%name)

                         ! bond order -> glyph that visually depicts the bond, and word
                         ord = sys(isys)%c%nstar(i)%ordcon(j)
                         select case (ord)
                         case (-1) ! aromatic (1.5)
                            bondglyph = "○"
                            bondword = "aromatic"
                         case (0) ! dashed (drawing only)
                            bondglyph = "┄"
                            bondword = "dashed"
                         case (1) ! single
                            bondglyph = "─"
                            bondword = "single"
                         case (2) ! double
                            bondglyph = "═"
                            bondword = "double"
                         case (3) ! triple
                            bondglyph = "≡"
                            bondword = "triple"
                         case default
                            bondglyph = "?"
                            bondword = "unknown"
                         end select

                         ! color the button with the bonded atom's color from the view
                         havergb_ = color_from_view(atlisttype_ncel_frac,jj,rgb)
                         if (havergb_) then
                            col4 = ImVec4(rgb(1),rgb(2),rgb(3),1._c_float)
                            call igPushStyleColor_Vec4(ImGuiCol_Button,col4)
                            col4 = ImVec4(min(rgb(1)*ColorButtonHoverFactor,1._c_float),&
                               min(rgb(2)*ColorButtonHoverFactor,1._c_float),&
                               min(rgb(3)*ColorButtonHoverFactor,1._c_float),1._c_float)
                            call igPushStyleColor_Vec4(ImGuiCol_ButtonHovered,col4)
                            col4 = ImVec4(rgb(1)*ColorButtonActiveFactor,rgb(2)*ColorButtonActiveFactor,&
                               rgb(3)*ColorButtonActiveFactor,1._c_float)
                            call igPushStyleColor_Vec4(ImGuiCol_ButtonActive,col4)
                            ! readable label: black on light atoms, white on dark ones
                            lum = lumweights(1)*rgb(1)+lumweights(2)*rgb(2)+lumweights(3)*rgb(3)
                            if (lum > 0.5_c_float) then
                               col4 = ColorBlack
                            else
                               col4 = ColorWhite
                            end if
                            call igPushStyleColor_Vec4(ImGuiCol_Text,col4)
                         end if

                         ! the button (colored by the bonded atom, prefixed by
                         ! the bond-type glyph). Left-click opens a menu,
                         ! right-click removes the bond.
                         ldum = iw_button(bondglyph // " " // string(jj) //&
                            "##bond" // suffix // "_" // string(j),sameline=.true.)
                         if (havergb_) call igPopStyleColor(4)

                         ! hovering highlights the row's atom and marks this
                         ! neighbor, and shows a tooltip with the bond details;
                         ! right-click removes the bond (deferred)
                         if (igIsItemHovered(ImGuiHoveredFlags_None)) then
                            ihlbond = i
                            ihlbtn = jj
                            ! highlight the row in the table
                            color = igGetColorU32_Vec4(ColorTableHighlightRow)
                            call igTableSetBgColor(ImGuiTableBgTarget_RowBg0, color, -1)
                            ! bond distance and sum of covalent radii (in angstrom)
                            dbond = norm2(sys(isys)%c%atcel(jj)%r + sys(isys)%c%x2c(dble(sys(isys)%c%nstar(i)%lcon(:,j))) -&
                               sys(isys)%c%atcel(i)%r)
                            zi = sys(isys)%c%spc(sys(isys)%c%atcel(i)%is)%z
                            zj = sys(isys)%c%spc(sys(isys)%c%atcel(jj)%is)%z
                            call iw_tooltip(trim(name) // " (cell id " // string(jj) // "), " // bondword // " bond" // newline //&
                               "Bond length: " // string(dbond*bohrtoa,'f',decimal=4) // " Å")
                            if (igIsMouseClicked(ImGuiMouseButton_Right,.false._c_bool)) then
                               ibrm1 = i
                               ibrm2 = jj
                               lbrm = sys(isys)%c%nstar(i)%lcon(:,j)
                            end if
                         end if

                         ! left-click opens the bond menu (Order submenu + Remove)
                         str2 = "##bondmenu" // suffix // "_" // string(j) // c_null_char
                         if (ldum) call igOpenPopup_Str(c_loc(str2),ImGuiPopupFlags_None)
                         if (igBeginPopup(c_loc(str2),ImGuiWindowFlags_None)) then
                            str1 = "Order" // c_null_char
                            if (igBeginMenu(c_loc(str1),.true._c_bool)) then
                               ! tick marks the current ordcon value for this bond
                               if (iw_menuitem("Single",selected=(sys(isys)%c%nstar(i)%ordcon(j)==1)))&
                                  call defer_setorder(1)
                               if (iw_menuitem("Double",selected=(sys(isys)%c%nstar(i)%ordcon(j)==2)))&
                                  call defer_setorder(2)
                               if (iw_menuitem("Triple",selected=(sys(isys)%c%nstar(i)%ordcon(j)==3)))&
                                  call defer_setorder(3)
                               if (iw_menuitem("Dashed",selected=(sys(isys)%c%nstar(i)%ordcon(j)==0)))&
                                  call defer_setorder(0)
                               if (iw_menuitem("Aromatic",selected=(sys(isys)%c%nstar(i)%ordcon(j)==-1)))&
                                  call defer_setorder(-1)
                               call igEndMenu()
                            end if
                            if (iw_menuitem("Remove")) then
                               ibrm1 = i
                               ibrm2 = jj
                               lbrm = sys(isys)%c%nstar(i)%lcon(:,j)
                            end if
                            call igEndPopup()
                         end if
                      end do
                   end if
                end do
             end do
             call igEndTable()
          end if

          ! perform deferred bond edits (after the table loop, so nstar is not
          ! edited while being iterated)
          if (ibrm1 > 0) call sysc(isys)%remove_bond(ibrm1,ibrm2,lbrm)
          if (ibord1 > 0) call sysc(isys)%set_bond_order(ibord1,ibord2,lbord,ibordval)

          ! separator + recalculate bonds section
          call igSeparator()
          call iw_text("Recalculate Bonds",highlight=.true.)
          call iw_text("Atoms A and B are bonded if:")
          call iw_text("  A and B non-metals: d < (r_cov(i)+r_cov(j))*f"//newline//&
             "  A or B metal:       d < d_NN + δ"//newline//&
             "r_cov = covalent radius. d_NN = nearest-neighbor distance.")

          ! atomic radii + allowed-bond matrix table
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_Resizable)
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          flags = ior(flags,ImGuiTableFlags_ScrollX)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          str1 = "##tableatomrcov_geom" // c_null_char
          sz0%x = 0
          sz0%y = iw_calcheight(min(5,natused_bonds)+1,0,.false.)
          if (igBeginTable(c_loc(str1),3+natused_bonds,flags,sz0,0._c_float)) then
             str2 = "Atom" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,0)
             str2 = "Z" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,1)
             str2 = "Radius (Å)" // c_null_char
             call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,2)
             do j = 1, natused_bonds
                str2 = trim(name_bonds(j)) // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,2+j)
             end do
             call igTableSetupScrollFreeze(0,1)
             call igTableHeadersRow()

             do i = 1, natused_bonds
                call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                iz = iat_bonds(i)
                if (igTableSetColumnIndex(0)) then
                   call igAlignTextToFramePadding()
                   call iw_text(trim(name_bonds(i)))
                end if
                if (igTableSetColumnIndex(1)) then
                   call igAlignTextToFramePadding()
                   call iw_text(string(iz))
                end if
                if (igTableSetColumnIndex(2)) then
                   ldum = iw_dragfloat_real8("##tableradius_geom" // string(i),x1=sysc(isys)%atmcov(iz),&
                      speed=0.01d0,min=0d0,max=2.65d0,scale=bohrtoa,decimal=3,&
                      flags=ImGuiSliderFlags_AlwaysClamp)
                end if
                ! allowed-bond checkboxes, one per element (symmetric mask by Z)
                do j = 1, natused_bonds
                   if (igTableSetColumnIndex(2+j)) then
                      ! the mask is allocated lazily: unallocated means all pairs
                      ! allowed, so default to checked and only allocate on a change
                      if (allocated(sysc(isys)%bondallowed)) then
                         ballow = sysc(isys)%bondallowed(iz,iat_bonds(j))
                      else
                         ballow = .true.
                      end if
                      if (iw_checkbox("##bondallowed_" // string(i) // "_" // string(j),ballow)) then
                         if (.not.allocated(sysc(isys)%bondallowed)) then
                            allocate(sysc(isys)%bondallowed(maxzat0,maxzat0))
                            sysc(isys)%bondallowed = .true.
                         end if
                         sysc(isys)%bondallowed(iz,iat_bonds(j)) = ballow
                         sysc(isys)%bondallowed(iat_bonds(j),iz) = ballow
                      end if
                      call iw_tooltip("Allow " // trim(name_bonds(i)) // "-" //&
                         trim(name_bonds(j)) // " bonds to form when recalculating",ttshown)
                   end if
                end do
             end do
             call igEndTable()
          end if
          deallocate(iat_bonds,name_bonds)

          ! bond factor drag
          call igAlignTextToFramePadding()
          ldum = iw_dragfloat_real8("Bond factor (f)##bondfactor_geom",x1=sysc(isys)%bondfactor,speed=0.001d0,&
             min=1d0,max=4d0,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Bond factor parameter (multiplicative) for non-metal bonding (see formula above)",ttshown)
          call iw_text(" ",sameline=.true.)

          ! bond delta drag
          ldum = iw_dragfloat_real8("Bond delta δ (Å)##bonddelta_geom",x1=sysc(isys)%bonddelta,speed=0.001d0,&
             min=0d0,max=2d0,scale=bohrtoa,decimal=4,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Distance tolerance (additive) for metal-non-metal bonding (see formula above)",ttshown)

          ! allow all / allow none / reset / apply buttons. The mask is freed
          ! when it returns to fully-allowed (the unallocated = all-allowed default)
          if (iw_button("Allow All##allowallbondtypes")) then
             if (allocated(sysc(isys)%bondallowed)) deallocate(sysc(isys)%bondallowed)
          end if
          call iw_tooltip("Allow all bond types to form",ttshown)
          if (iw_button("Allow None##allownonebondtypes",sameline=.true.)) then
             if (.not.allocated(sysc(isys)%bondallowed)) allocate(sysc(isys)%bondallowed(maxzat0,maxzat0))
             sysc(isys)%bondallowed = .false.
          end if
          call iw_tooltip("Forbid all bond types from forming",ttshown)
          if (iw_button("Reset",sameline=.true.)) then
             sysc(isys)%atmcov = atmcov0
             sysc(isys)%bondfactor = bondfactor_def
             sysc(isys)%bonddelta = bonddelta_def
             if (allocated(sysc(isys)%bondallowed)) deallocate(sysc(isys)%bondallowed)
          end if
          call iw_tooltip("Reset atomic radii, bond factor, bond delta, and allowed bond types to defaults",ttshown)
          if (iw_button("Apply",danger=.true.,sameline=.true.)) &
             call sysc(isys)%rebond()
          call iw_tooltip("Recalculate bonds with the parameters above",ttshown)

          call igEndTabItem()
       end if

       !! symmetry tab !!
       str2 = "Symmetry##drawgeometry_symmetrytab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str2),c_null_ptr,flags)) then
          ! check if the tab changed
          call check_changed_tab("symmetry")

          if (sys(isys)%c%ismolecule) then
             !! symmetry: molecule !!

             ! current point group
             call iw_text("Point Group",highlight=.true.)
             if (sys(isys)%c%pg%avail) then
                call iw_text("  " // trim(sys(isys)%c%pg%symbol))
             else
                call iw_text("  Not available (symmetry not computed)")
             end if

             ! symmetry operations table
             if (sys(isys)%c%pg%avail) then
                call iw_text("Operations",highlight=.true.)
                call iw_text("(" // string(sys(isys)%c%pg%nop) // ")",sameline=.true.)

                flags = ImGuiTableFlags_None
                flags = ior(flags,ImGuiTableFlags_RowBg)
                flags = ior(flags,ImGuiTableFlags_Borders)
                flags = ior(flags,ImGuiTableFlags_ScrollY)
                flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
                str1 = "##symopstablemol" // c_null_char
                sz0%x = 0
                ! size to fit all operations, but cap at the space left above the
                ! Close button (then scroll) so the table never stretches past it
                sz0%y = iw_calcheight(sys(isys)%c%pg%nop+1,0,.false.)
                call igGetContentRegionAvail(szavail)
                szavail%y = szavail%y - iw_calcheight(5,0,.true.)
                if (sz0%y > szavail%y) sz0%y = szavail%y
                call symop_ensure_sel(sys(isys)%c%pg%nop)
                if (igBeginTable(c_loc(str1),3,flags,sz0,0._c_float)) then
                   str2 = "#" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,0)
                   str2 = "Sym" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,1)
                   str2 = "Axis (Å)" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,2)
                   call igTableSetupScrollFreeze(0,1)
                   call igTableHeadersRow()

                   do i = 1, sys(isys)%c%pg%nop
                      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                      call symop_selectbg(i)
                      ! axis/normal (Cartesian unit vector); identity, inversion
                      ! and unknown operations have no meaningful axis
                      ioptype = sys(isys)%c%pg%op(i)%type
                      if (ioptype == molsymop_rotation .or. ioptype == molsymop_imp_rotation .or.&
                         ioptype == molsymop_plane) then
                         raxx = sys(isys)%c%pg%op(i)%axis
                         saxx = "[" // string(raxx(1),'f',length=6,decimal=3,justify=ioj_right) // "," //&
                            string(raxx(2),'f',length=6,decimal=3,justify=ioj_right) // "," //&
                            string(raxx(3),'f',length=6,decimal=3,justify=ioj_right) // "]"
                      else
                         saxx = ""
                      end if
                      if (igTableSetColumnIndex(0)) call iw_text(string(i))
                      ! whole-row select/hover: emit even when column 0 is clipped off-screen
                      if (iw_highlight_selectable("##symopselectmol_" // string(i),clicked=clicked)) ihl_symop = i
                      if (clicked) call symop_click(i)
                      if (igTableSetColumnIndex(1)) call iw_text(trim(sys(isys)%c%pg%op(i)%sym))
                      if (igTableSetColumnIndex(2)) call iw_text(saxx)
                   end do
                   call igEndTable()
                end if

                ! draw the selected/hovered symmetry elements and the selection row
                call symop_display_and_buttons(sys(isys)%c%pg%nop)
             end if
          else
             !! symmetry: crystal !!
             ! current space group
             call iw_text("Space Group",highlight=.true.)
             if (sys(isys)%c%spgavail) then
                call pointgroup_info(sys(isys)%c%spg%pointgroup_symbol,schpg,holo,laue)
                call iw_text("  Hermann-Mauguin: " // trim(sys(isys)%c%spg%international_symbol) //&
                   " (number " // string(sys(isys)%c%spg%spacegroup_number) // ")")
                call iw_text("  Point group: " // trim(sys(isys)%c%spg%pointgroup_symbol) //&
                   " (" // trim(schpg) // ")")
                call iw_text("  Holohedry: " // trim(holo_string(holo)) //&
                   ", Laue class: " // trim(laue_string(laue)))
             else
                call iw_text("  Not available (symmetry not computed)")
             end if

             call igAlignTextToFramePadding()
             if (iw_button("Recalculate##symrecalc")) iaction = iaction_sym_recalc
             call iw_tooltip("Recalculate crystal symmetry at the given tolerance (symeps) level",ttshown)
             call iw_text("at tolerance:",sameline=.true.)
             str_eps = string(sysc(isys)%symeps,'e',decimal=2)
             if (iw_inputtext("##symepsinput",bufsize=32,texta=str_eps,width=8,sameline=.true.,notlive=.true.)) then
                if (isreal(res,str_eps)) then
                   if (res > 0d0) sysc(isys)%symeps = res
                end if
             end if
             call iw_tooltip("Distance tolerance (bohr) used to detect the symmetry (symeps)",ttshown)

             if (iw_button("Clear##symclear",danger=.true.)) iaction = iaction_sym_clear
             call iw_tooltip("Remove crystal symmetry",ttshown)
             if (iw_button("Refine##symrefine",sameline=.true.)) iaction = iaction_sym_refine
             call iw_tooltip("Idealize the cell parameters and move the atoms to their ideal &
                &symmetry positions.",ttshown)
             if (iw_button("Whole molecules##symwholemols",sameline=.true.,&
                disabled=.not.sys(isys)%c%ismol3d)) iaction = iaction_sym_wholemols
             call iw_tooltip("Choose a symmetry subgroup such that the asymmetric unit contains whole molecules &
                &(molecular crystals only)",ttshown)
             if (len_trim(w%errmsg) > 0) &
                call iw_text(w%errmsg,danger=.true.)

             ! symmetry operations table
             neqv = sys(isys)%c%neqv
             ncv = sys(isys)%c%ncv
             call iw_text("Operations",highlight=.true.)
             call iw_text("(" // string(neqv) // ")",sameline=.true.)

             if (neqv >= 1) then
                ! the operation strings are expensive to build (axis analysis
                ! does an eigendecomposition per operation), so cache them and
                ! rebuild only when the symmetry changes (clear_sym_cache
                ! invalidates the cache on any geometry change)
                if (allocated(w%geometry_sym_ops)) then
                   if (size(w%geometry_sym_ops,1) /= neqv*max(ncv,1)) then
                      deallocate(w%geometry_sym_ops)
                      if (allocated(w%geometry_sym_hm)) deallocate(w%geometry_sym_hm)
                      if (allocated(w%geometry_sym_axes)) deallocate(w%geometry_sym_axes)
                   end if
                end if
                if (.not.allocated(w%geometry_sym_ops)) then
                   allocate(w%geometry_sym_ops(neqv*max(ncv,1)),w%geometry_sym_hm(neqv*max(ncv,1)),&
                      w%geometry_sym_axes(3,neqv*max(ncv,1)))
                   call sys(isys)%c%struct_report_symxyz(w%geometry_sym_ops,hmsym=w%geometry_sym_hm,&
                      axcr=w%geometry_sym_axes)
                end if

                flags = ImGuiTableFlags_None
                flags = ior(flags,ImGuiTableFlags_RowBg)
                flags = ior(flags,ImGuiTableFlags_Borders)
                flags = ior(flags,ImGuiTableFlags_ScrollY)
                flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
                str1 = "##symopstable" // c_null_char
                sz0%x = 0
                sz0%y = iw_calcheight(min(neqv,8)+1,0,.false.)
                call symop_ensure_sel(neqv)
                if (igBeginTable(c_loc(str1),5,flags,sz0,0._c_float)) then
                   str2 = "#" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,0)
                   str2 = "HM" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,1)
                   str2 = "Operation" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,2)
                   str2 = "Axis (cryst.)" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,3)
                   str2 = "Axis (Cartesian)" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,4)
                   call igTableSetupScrollFreeze(0,1)
                   call igTableHeadersRow()

                   ! show only the operations under the identity centering (1..neqv)
                   do i = 1, neqv
                      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                      call symop_selectbg(i)
                      ! axis in crystallographic coordinates (cached) and Cartesian
                      ! (computed on the fly); identity/inversion have a zero axis
                      raxc = w%geometry_sym_axes(:,i)
                      if (norm2(raxc) < 1d-10) then
                         saxc = ""
                         saxx = ""
                      else
                         if (all(abs(raxc - nint(raxc)) < 1d-10)) then
                            saxc = "[" // string(nint(raxc(1)),2,justify=ioj_right) // "," //&
                               string(nint(raxc(2)),2,justify=ioj_right) // "," //&
                               string(nint(raxc(3)),2,justify=ioj_right) // "]"
                         else
                            saxc = "[" // string(raxc(1),'f',length=4,decimal=1,justify=ioj_right) // "," //&
                               string(raxc(2),'f',length=4,decimal=1,justify=ioj_right) // "," //&
                               string(raxc(3),'f',length=4,decimal=1,justify=ioj_right) // "]"
                         end if
                         raxx = sys(isys)%c%x2c(raxc)
                         raxx = raxx / norm2(raxx)
                         saxx = "[" // string(raxx(1),'f',length=6,decimal=3,justify=ioj_right) // "," //&
                            string(raxx(2),'f',length=6,decimal=3,justify=ioj_right) // "," //&
                            string(raxx(3),'f',length=6,decimal=3,justify=ioj_right) // "]"
                      end if
                      if (igTableSetColumnIndex(0)) call iw_text(string(i))
                      ! whole-row select/hover: the row-spanning selectable must be
                      ! emitted even when column 0 is clipped (window dragged past
                      ! the left edge), otherwise the highlight stops working
                      if (iw_highlight_selectable("##symopselect_" // string(i),clicked=clicked)) ihl_symop = i
                      if (clicked) call symop_click(i)
                      if (igTableSetColumnIndex(1)) call iw_text(string(w%geometry_sym_hm(i),length=3,justify=ioj_right))
                      if (igTableSetColumnIndex(2)) call iw_text(trim(w%geometry_sym_ops(i)))
                      if (igTableSetColumnIndex(3)) call iw_text(saxc)
                      if (igTableSetColumnIndex(4)) call iw_text(saxx)
                   end do
                   call igEndTable()
                end if

                ! draw the selected/hovered symmetry elements and the selection row
                call symop_display_and_buttons(neqv)
             end if

             ! centering vectors table
             call iw_text("Centering Vectors",highlight=.true.)
             call iw_text("(" // string(ncv) // ")",sameline=.true.)

             flags = ImGuiTableFlags_None
             flags = ior(flags,ImGuiTableFlags_RowBg)
             flags = ior(flags,ImGuiTableFlags_Borders)
             flags = ior(flags,ImGuiTableFlags_ScrollY)
             flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
             str1 = "##symcentable" // c_null_char
             sz0%x = 0
             sz0%y = iw_calcheight(min(ncv,4)+1,0,.false.)
             if (igBeginTable(c_loc(str1),2,flags,sz0,0._c_float)) then
                str2 = "#" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,0)
                str2 = "Coordinates (fractional)" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,1)
                call igTableSetupScrollFreeze(0,1)
                call igTableHeadersRow()
                do i = 1, ncv
                   call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                   if (igTableSetColumnIndex(0)) call iw_text(string(i))
                   if (igTableSetColumnIndex(1)) call iw_text(cell_cen_label(sys(isys)%c%cen(:,i)))
                end do
                call igEndTable()
             end if

             ! space group vs tolerance analysis
             call igAlignTextToFramePadding()
             call iw_text("Space Group Analysis",highlight=.true.)
             call iw_helpermark("Calculate the space group as a function of the symmetry &
                &tolerance (symprec). Click a row to adopt that tolerance and recalculate.")
             if (iw_button("Analyze##symanalyze",sameline=.true.)) iaction = iaction_sym_analyze
             call iw_tooltip("Scan the space group over a range of symmetry tolerances",ttshown)

             if (allocated(w%geometry_sym_analyze_eps)) then
                flags = ImGuiTableFlags_None
                flags = ior(flags,ImGuiTableFlags_RowBg)
                flags = ior(flags,ImGuiTableFlags_Borders)
                flags = ior(flags,ImGuiTableFlags_ScrollY)
                flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
                str1 = "##symanaltable" // c_null_char
                sz0%x = 0
                sz0%y = iw_calcheight(size(w%geometry_sym_analyze_eps,1)+1,0,.false.)
                if (igBeginTable(c_loc(str1),3,flags,sz0,0._c_float)) then
                   str2 = "symprec" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,0)
                   str2 = "Space group" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,1)
                   str2 = "Number" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,2)
                   call igTableSetupScrollFreeze(0,1)
                   call igTableHeadersRow()
                   do i = 1, size(w%geometry_sym_analyze_eps,1)
                      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                      if (igTableSetColumnIndex(0)) then
                         str2 = string(w%geometry_sym_analyze_eps(i),'e',10,2) //&
                            "##symanalrow" // string(i) // c_null_char
                         if (igSelectable_Bool(c_loc(str2),logical(.false.,c_bool),&
                            ImGuiSelectableFlags_SpanAllColumns,szero)) then
                            ! adopt this row's tolerance and recalculate
                            sysc(isys)%symeps = w%geometry_sym_analyze_eps(i)
                            iaction = iaction_sym_recalc
                         end if
                      end if
                      if (igTableSetColumnIndex(1)) call iw_text(trim(w%geometry_sym_analyze_sym(i)))
                      if (igTableSetColumnIndex(2)) call iw_text(string(w%geometry_sym_analyze_num(i)))
                   end do
                   call igEndTable()
                end if
             end if
          end if ! .not. ismolecule

          call igEndTabItem()
       end if

       call igEndTabBar()
    end if

    call igEndGroup()

    ! hover highlight, plus the atom awaiting an add-bond pick (transient
    ! highlights clear each frame, so they must go in a single call)
    if (ihlbond > 0) then
       ! bonds tab: highlight the hovered atom (central color) and its bonded
       ! neighbors (a different color; lighter for a hovered neighbor button)
       ncon = 0
       if (allocated(sys(isys)%c%nstar)) ncon = sys(isys)%c%nstar(ihlbond)%ncon
       allocate(ihigh(2+ncon),irgba(4,2+ncon))
       ihigh(1) = ihlbond
       irgba(:,1) = ColorHighlightScene
       do j = 1, ncon
          jj = sys(isys)%c%nstar(ihlbond)%idcon(j)
          ihigh(1+j) = jj
          if (jj == ihlbtn) then
             irgba(:,1+j) = ColorHighlightBondScene2
          else
             irgba(:,1+j) = ColorHighlightBondScene
          end if
       end do
       nhigh = 1 + ncon
       if (ipickhl > 0) then
          nhigh = nhigh + 1
          ihigh(nhigh) = ipickhl
          irgba(:,nhigh) = ColorHighlightScene
       end if
       call sysc(isys)%highlight_atoms(.true.,ihigh(1:nhigh),atlisttype_ncel_frac,irgba(:,1:nhigh))
       deallocate(ihigh,irgba)
    elseif (ihighlight > 0 .and. ipickhl == 0) then
       call sysc(isys)%highlight_atoms(.true.,(/ihighlight/),table_hltype,&
          reshape(ColorHighlightScene,(/4,1/)))
    elseif (ipickhl > 0) then
       call sysc(isys)%highlight_atoms(.true.,(/ipickhl/),atlisttype_ncel_frac,&
          reshape(ColorHighlightScene,(/4,1/)))
    end if

    ! hovering a molecule row: show its standard-orientation axes at the COM
    if (ihighlight > 0 .and. table_hltype == atlisttype_nmol) then
       if (sysc(isys)%sc%isinit /= 0 .and. sys(isys)%c%mol(ihighlight)%discrete) then
          call sys(isys)%c%mol(ihighlight)%standard_axes(m_std=stdrot,xcm=stdcom)
          ! axis length scaled to the molecule's size (floor for tiny ones)
          stdext = 0d0
          do jm = 1, sys(isys)%c%mol(ihighlight)%nat
             stdext = max(stdext,norm2(sys(isys)%c%mol(ihighlight)%at(jm)%r - stdcom))
          end do
          stdaxlen = max(0.5d0 * stdext, 1.5d0)
          stdcom = stdcom + sys(isys)%c%molx0

          ! standard-orientation axes
          call sysc(isys)%sc%show_transient_axes(ihighlight,stdcom,stdrot,stdaxlen)

          if (ieuler_drag /= 0) then
             ! dragging a Euler angle: also show the rotation axis for that angle
             ! (ZYZ convention, in the scene/lab cartesian frame):
             !   1 (alpha) -> lab Z; 2 (beta) -> line of nodes; 3 (gamma) -> body Z
             if (ieuler_drag == 1) then
                rotdir = (/0d0,0d0,1d0/)
             elseif (ieuler_drag == 2) then
                rotdir = (/-stdrot(2,3),stdrot(1,3),0d0/) ! zlab x bodyZ
                if (norm2(rotdir) < 1d-6) then
                   rotdir = (/0d0,1d0,0d0/) ! gimbal pole: line of nodes undefined
                else
                   rotdir = rotdir / norm2(rotdir)
                end if
             else
                rotdir = stdrot(:,3) ! molecule body Z (already a unit vector)
             end if
             rotlen = max(stdext,1.5d0) * 1.2d0
             call sysc(isys)%sc%show_transient_rotaxis(-(ihighlight*4+ieuler_drag),stdcom,rotdir,rotlen)
          end if
       end if
    end if

    ! the Euler-angle drag has ended (no drag widget active this frame): drop the
    ! persisted unwrapped angles so the table shows the canonical decomposition again
    if (ieuler_drag == 0) w%geometry_euler_drag_mol = 0

    ! process clicked: write the selection directly to the system
    if (iclicked > 0) then
       ! single click: toggle all cell atoms of this row (a fully selected row
       ! is cleared; a partial or unselected row is fully selected)
       call current_row_state(iclicked,istate,rrgba)
       if (istate == 2) then
          call sysc(isys)%highlight_clear(.false.,(/iclicked/),table_hltype)
       else
          call sysc(isys)%highlight_atoms(.false.,(/iclicked/),table_hltype,&
             reshape(w%geometry_select_rgba,(/4,1/)))
       end if
    elseif (iclicked == -1) then
       ! range (shift) or single (ctrl) selection: select the rows
       nhigh = iclicked_end - iclicked_ini + 1
       allocate(ihigh(nhigh),irgba(4,nhigh))
       nhigh = 0
       do ii = iclicked_ini, iclicked_end
          nhigh = nhigh + 1
          ihigh(nhigh) = w%iord(ii)
          irgba(:,nhigh) = w%geometry_select_rgba
       end do
       call sysc(isys)%highlight_atoms(.false.,ihigh,table_hltype,irgba)
       deallocate(ihigh,irgba)
    end if

    ! remove/merge highlighted atoms
    if (w%focused() .and. is_bind_event(BIND_EDITSELECT_REMOVE)) then
       iaction = iaction_edit_highlighted
       iaction_i1 = edit_remove
    end if

    ! deselect all highlighted atoms
    deselected = .false.
    if (w%focused() .and. is_bind_event(BIND_EDITSELECT_DESELECT)) then
       if (allocated(sysc(isys)%highlight_rgba)) then
          if (any(sysc(isys)%highlight_rgba >= 0._c_float)) then
             call sysc(isys)%highlight_clear(.false.)
             deselected = .true.
          end if
       end if
    end if

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(5,1,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! close button
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if (.not.deselected .and. ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS))) &
       doquit = .true.
    doquit = doquit .or. iw_button("Close")

    ! quit the window
    if (doquit) &
       call w%end()

    ! process actions at the end; clear any previous error when a new action runs
    if (iaction /= -1) w%errmsg = ""
    if (iaction == iaction_restore) then
       call sysc(isys)%reread_geometry_from_file()

    elseif (iaction == iaction_set_attype_name) then
       call sysc(isys)%set_attype_name(w%geometry_atomtype,iaction_i1,iaction_str)

    elseif (iaction == iaction_set_atomic_number) then
       call sysc(isys)%set_atomic_number(w%geometry_atomtype,iaction_i1,iaction_i2,setatomnames=.true.)

    elseif (iaction == iaction_add_species) then
       call sysc(isys)%add_species(iaction_i1)

    elseif (iaction == iaction_set_attype_species) then
       call sysc(isys)%set_attype_species(w%geometry_atomtype,iaction_i1,iaction_i2,&
          copybonding=w%geometry_keepbonding)

    elseif (iaction == iaction_set_atom_position) then
       call sysc(isys)%set_atom_position(w%geometry_atomtype,iaction_i1,iaction_x,iaction_l,&
          copybonding=w%geometry_keepbonding)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_set_atom_occupancy) then
       call sysc(isys)%set_attype_occupancy(w%geometry_atomtype,iaction_i1,iaction_x(1))
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_add_species_change_atom) then
       call sysc(isys)%add_species(iaction_i1)
       call sysc(isys)%set_attype_species(w%geometry_atomtype,iaction_i2,sys(isys)%c%nspc,&
          copybonding=w%geometry_keepbonding)

    elseif (iaction == iaction_add_atom) then
       call sysc(isys)%attype_add_atom(w%geometry_atomtype,iaction_i1,iaction_x)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_edit_highlighted) then
       if (w%geometry_atomtype == atlisttype_species) then
          ! rowstate is the per-species selection state built this frame in the
          ! species tab (size nspc); fully selected species have rowstate == 2
          call sysc(isys)%edit_highlighted_species(rowstate == 2,remove=(iaction_i1==edit_remove),&
             merge=(iaction_i1==edit_merge),duplicate=(iaction_i1==edit_duplicate),errmsg=w%errmsg)
       else
          call sysc(isys)%edit_highlighted_atoms(remove=(iaction_i1==edit_remove),&
             merge=(iaction_i1==edit_merge),duplicate=(iaction_i1==edit_duplicate),errmsg=w%errmsg)
       end if
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_reorder_highlighted) then
       call sysc(isys)%attype_reorder(w%geometry_atomtype,w%iord)

    elseif (iaction == iaction_swap_atom_ids) then
       call sysc(isys)%attype_swap_atoms(w%geometry_atomtype,iaction_i1,iaction_i2)

    elseif (iaction == iaction_swap_mol_ids) then
       call sysc(isys)%swap_molecules(iaction_i1,iaction_i2)

    elseif (iaction == iaction_set_molecule_position) then
       call sysc(isys)%set_molecule_position(w%geometry_moltype,iaction_i1,iaction_x,&
          copybonding=w%geometry_keepbonding)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_set_molecule_rotation) then
       call sysc(isys)%set_molecule_rotation(iaction_i1,iaction_x,copybonding=w%geometry_keepbonding)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_remove_molecules) then
       call sysc(isys)%edit_highlighted_atoms(remove=.true.,errmsg=w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_reorder_molecules) then
       call sysc(isys)%attype_reorder(atlisttype_nmol,w%iord)

    elseif (iaction == iaction_change_cell) then
       call sysc(isys)%move_cell(iaction_x6(1:3),iaction_x6(4:6),iaction_l,&
          copybonding=w%geometry_keepbonding)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_transform_cell) then
       call sysc(isys)%transform_cell(iaction_i1,.false.,errmsg=w%errmsg)
       if (allocated(w%geometry_cell_nice_rmax)) deallocate(w%geometry_cell_nice_rmax)
       if (allocated(w%geometry_cell_nice_mmax)) deallocate(w%geometry_cell_nice_mmax)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_transform_matrix) then
       call sysc(isys)%transform_cell_matrix(iaction_m,iaction_x,iaction_l,errmsg=w%errmsg)
       if (allocated(w%geometry_cell_nice_rmax)) deallocate(w%geometry_cell_nice_rmax)
       if (allocated(w%geometry_cell_nice_mmax)) deallocate(w%geometry_cell_nice_mmax)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_sym_recalc) then
       call sysc(isys)%recalc_symmetry(w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_sym_clear) then
       call sysc(isys)%clear_symmetry()
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_sym_refine) then
       call sysc(isys)%refine_symmetry(w%errmsg)
       ! refine changes the cell metric and atomic positions, so the
       ! nice-supercell search results are stale
       if (allocated(w%geometry_cell_nice_rmax)) deallocate(w%geometry_cell_nice_rmax)
       if (allocated(w%geometry_cell_nice_mmax)) deallocate(w%geometry_cell_nice_mmax)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_sym_wholemols) then
       call sysc(isys)%wholemols_op(w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_sym_analyze) then
       call sysc(isys)%spg_analysis(w%geometry_sym_analyze_eps,w%geometry_sym_analyze_sym,&
          w%geometry_sym_analyze_num)

    elseif (iaction == iaction_sym_deleteops) then
       if (allocated(w%geometry_sym_sel)) then
          call sysc(isys)%reduce_symmetry(w%geometry_sym_sel,w%errmsg)
          ! neqv changed: drop the cached operations table + selection so they
          ! are rebuilt at the new operation count
          call clear_sym_cache()
          sysc(isys)%sc%nextbuildlists_fixcam = .true.
       end if

    end if

  contains
    ! format a centering vector as "(f1,f2,f3)" with common fractions
    function cell_cen_label(x) result(str)
      real*8, intent(in) :: x(3)
      character(len=:,kind=c_char), allocatable :: str
      str = "(" // frac_str(x(1)) // "," // frac_str(x(2)) // "," // frac_str(x(3)) // ")"
    end function cell_cen_label

    ! format a fractional coordinate in [0,1) as a small fraction when possible
    function frac_str(v) result(str)
      real*8, intent(in) :: v
      character(len=:,kind=c_char), allocatable :: str
      real*8 :: w
      real*8, parameter :: eps = 1d-4
      w = v - floor(v)
      if (abs(w) < eps .or. abs(w-1d0) < eps) then
         str = "0"
      elseif (abs(w-0.5d0) < eps) then
         str = "1/2"
      elseif (abs(w-1d0/3d0) < eps) then
         str = "1/3"
      elseif (abs(w-2d0/3d0) < eps) then
         str = "2/3"
      elseif (abs(w-0.25d0) < eps) then
         str = "1/4"
      elseif (abs(w-0.75d0) < eps) then
         str = "3/4"
      else
         str = string(w,'f',decimal=3)
      end if
    end function frac_str

    ! change the system on which the geometry window operates
    subroutine change_system(i)
      integer, intent(in) :: i

      ! do nothing if we are already in the same system
      if (w%isys == i) return

      ! reset the last-selected row and the table sort
      w%lastselected = 0

      ! remove the cached cell-transformation data and reorder the table
      if (allocated(w%geometry_cell_nice_rmax)) deallocate(w%geometry_cell_nice_rmax)
      if (allocated(w%geometry_cell_nice_mmax)) deallocate(w%geometry_cell_nice_mmax)
      call clear_sym_cache()
      call reset_sort()

      ! clear the cell transformations
      w%geometry_cell_nrep = 1_c_int
      w%geometry_cell_intmat = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
      w%geometry_cell_cen = 1
      w%geometry_cell_origin = 0d0

      ! change the system
      w%isys = i
      w%tied_to_tree = w%tied_to_tree .and. (w%isys == win(iwin_tree)%tree_selected)

    end subroutine change_system

    ! deallocate the sort array and force a sort of the table
    subroutine reset_sort()

      w%lastselected = 0
      if (allocated(w%iord)) deallocate(w%iord)
      forcesort = .true.

    end subroutine reset_sort

    ! Ensure the symmetry-operation selection array exists and is sized to nop
    ! (the number of operations in the table); reset it whenever the size changes.
    subroutine symop_ensure_sel(nop)
      integer, intent(in) :: nop

      if (allocated(w%geometry_sym_sel)) then
         if (size(w%geometry_sym_sel,1) /= nop) deallocate(w%geometry_sym_sel)
      end if
      if (.not.allocated(w%geometry_sym_sel)) then
         allocate(w%geometry_sym_sel(nop))
         w%geometry_sym_sel = .false.
      end if

    end subroutine symop_ensure_sel

    ! Set the background color of the current symmetry-operation table row when
    ! operation i is selected. Call right after igTableNextRow.
    subroutine symop_selectbg(i)
      integer, intent(in) :: i

      integer(c_int) :: color
      type(ImVec4) :: col4

      if (.not.allocated(w%geometry_sym_sel)) return
      if (i < 1 .or. i > size(w%geometry_sym_sel,1)) return
      if (.not.w%geometry_sym_sel(i)) return
      col4 = ImVec4(w%geometry_select_rgba(1),w%geometry_select_rgba(2),&
         w%geometry_select_rgba(3),w%geometry_select_rgba(4))
      color = igGetColorU32_Vec4(col4)
      call igTableSetBgColor(ImGuiTableBgTarget_RowBg0,color,-1_c_int)

    end subroutine symop_selectbg

    ! Process a click on symmetry-operation row i: plain click toggles the row,
    ! shift-click selects the range from the last-selected row.
    subroutine symop_click(i)
      integer, intent(in) :: i

      integer :: lo, hi

      if (.not.allocated(w%geometry_sym_sel)) return
      if (igIsKeyDown(ImGuiKey_ModShift).and.w%lastselected /= 0.and.w%lastselected /= i) then
         lo = min(i,w%lastselected)
         hi = max(i,w%lastselected)
         w%geometry_sym_sel(lo:hi) = .true.
      else
         w%geometry_sym_sel(i) = .not.w%geometry_sym_sel(i)
      end if
      w%lastselected = i
      w%geometry_sym_selgen = w%geometry_sym_selgen + 1

    end subroutine symop_click

    ! Push the selected (and hovered) symmetry elements to the view, then draw
    ! the Selection row (color picker + All/None/Toggle). nop is the number of
    ! operations shown in the table.
    subroutine symop_display_and_buttons(nop)
      integer, intent(in) :: nop

      integer :: i, n, tag, hovadd, kind1, order1, lioptype, idobj, iview, nsel
      integer, allocatable :: skind(:), sorder(:)
      real*8, allocatable :: sorig(:,:), sdir(:,:)
      real*8 :: orig1(3), dir1(3), lraxx(3), lraxc(3)
      character(len=1) :: lhm1, lcdig
      logical :: idnew

      ! display the selected and hovered symmetry elements in the view
      if (sysc(isys)%sc%isinit /= 0) then
         allocate(skind(nop),sorder(nop),sorig(3,nop),sdir(3,nop))
         n = 0
         do i = 1, nop
            if (.not.(w%geometry_sym_sel(i).or.i == ihl_symop)) cycle

            ! symmetry element for operation i: kind (symelem_kind_*; 0 = nothing
            ! to draw, i.e. identity/inversion/unknown), origin and direction
            ! (cartesian, bohr) and rotation order; molecule and crystal cases
            kind1 = 0
            order1 = 0
            orig1 = 0d0
            dir1 = 0d0
            if (sys(isys)%c%ismolecule) then
               lioptype = sys(isys)%c%pg%op(i)%type
               if (lioptype == molsymop_plane) then
                  kind1 = symop_kind_plane
               elseif (lioptype == molsymop_rotation .or. lioptype == molsymop_imp_rotation) then
                  kind1 = symop_kind_axis
               end if
               if (kind1 == 0) cycle
               lraxx = sys(isys)%c%pg%op(i)%axis
               if (norm2(lraxx) > 1d-10) lraxx = lraxx / norm2(lraxx)
               orig1 = sys(isys)%c%pg%xcm + sys(isys)%c%molx0
               dir1 = lraxx
               order1 = sys(isys)%c%pg%op(i)%opn
            else
               lraxc = w%geometry_sym_axes(:,i)
               if (norm2(lraxc) < 1d-10) cycle
               lhm1 = w%geometry_sym_hm(i)(1:1) ! HM symbol is stored left-aligned
               if (lhm1 >= "a" .and. lhm1 <= "z") then
                  kind1 = symop_kind_plane
               else
                  kind1 = symop_kind_axis
                  ! rotation order from the symbol: the digit, after an optional
                  ! leading "-" (rotoinversion: "-3"/"-4"/"-6")
                  if (lhm1 == "-") then
                     lcdig = w%geometry_sym_hm(i)(2:2)
                  else
                     lcdig = lhm1
                  end if
                  if (lcdig >= "0" .and. lcdig <= "9") order1 = ichar(lcdig) - ichar("0")
               end if
               lraxx = sys(isys)%c%x2c(lraxc)
               if (norm2(lraxx) > 1d-10) lraxx = lraxx / norm2(lraxx)
               ! crystal symmetry elements pass through the origin
               dir1 = lraxx
            end if

            n = n + 1
            skind(n) = kind1
            sorig(:,n) = orig1
            sdir(:,n) = dir1
            sorder(n) = order1
         end do
         ! unique tag: the selection generation combined with the hovered row
         ! (only when it adds to the selection), so the scene rebuilds the
         ! element list only when the displayed set actually changes
         hovadd = 0
         if (ihl_symop > 0 .and. ihl_symop <= nop) then
            if (.not.w%geometry_sym_sel(ihl_symop)) hovadd = ihl_symop
         end if
         tag = w%geometry_sym_selgen * (nop + 2) + (hovadd + 1)
         call sysc(isys)%sc%show_symelems(tag,n,skind(1:n),sorig(:,1:n),sdir(:,1:n),sorder(1:n))
         deallocate(skind,sorder,sorig,sdir)
      end if

      ! Selection row: all/none/toggle buttons (elements use the default colors)
      call igAlignTextToFramePadding()
      call iw_text("Selection",highlight=.true.)
      if (iw_button("All##symselall",sameline=.true.)) then
         w%geometry_sym_sel = .true.
         w%geometry_sym_selgen = w%geometry_sym_selgen + 1
      end if
      call iw_tooltip("Select all symmetry operations",ttshown)
      if (iw_button("None##symselnone",sameline=.true.)) then
         w%geometry_sym_sel = .false.
         w%geometry_sym_selgen = w%geometry_sym_selgen + 1
      end if
      call iw_tooltip("Deselect all symmetry operations",ttshown)
      if (iw_button("Toggle##symseltoggle",sameline=.true.)) then
         w%geometry_sym_sel = .not.w%geometry_sym_sel
         w%geometry_sym_selgen = w%geometry_sym_selgen + 1
      end if
      call iw_tooltip("Toggle the symmetry-operation selection",ttshown)

      ! create a persistent symmetry-elements object from the current selection,
      ! or reuse the existing one if there already is a symmetry-elements object
      if (iw_button("Create Object##symcreateobj",sameline=.true.,disabled=(sysc(isys)%sc%isinit==0))) then
         idobj = 0
         do i = 1, sysc(isys)%sc%nrep
            if (sysc(isys)%sc%rep(i)%isinit .and. sysc(isys)%sc%rep(i)%type == reptype_symelem) then
               idobj = i
               exit
            end if
         end do
         idnew = (idobj == 0)
         if (idnew) &
            call sysc(isys)%sc%add_representation(reptype_symelem,repflavor_symelem,id=idobj)
         if (idobj > 0) then
            ! ensure the operation snapshot/visibility style is initialized
            call sysc(isys)%sc%rep(idobj)%update()
            associate (rr => sysc(isys)%sc%rep(idobj))
              if (rr%symelem%style%isinit) then
                 if (size(rr%symelem%style%shown,1) == size(w%geometry_sym_sel,1)) then
                    if (idnew) then
                       ! new object: show exactly the selected operations
                       rr%symelem%style%shown = w%geometry_sym_sel
                    else
                       ! reuse: mark the selected operations as shown in it
                       where (w%geometry_sym_sel) rr%symelem%style%shown = .true.
                    end if
                    sysc(isys)%sc%forcebuildlists = .true.
                 end if
              end if
            end associate
         end if
         if (idobj > 0) then
            ! open the editor if a view window for this system is available
            iview = 0
            do i = 1, nwin
               if (win(i)%isinit .and. win(i)%isopen .and. win(i)%type == wintype_view) then
                  if (win(i)%view_selected == isys .and. associated(win(i)%sc)) then
                     iview = i
                     exit
                  end if
               end if
            end do
            if (iview > 0) &
               idobj = stack_create_window(wintype_editrep,.true.,isys=isys,irep=idobj,&
                  idparent=iview,orraise=-1)

            ! clear the table selection so the transient preview does not
            ! duplicate the elements now drawn by the persistent object
            w%geometry_sym_sel = .false.
            w%geometry_sym_selgen = w%geometry_sym_selgen + 1
         end if
      end if
      call iw_tooltip("Create a symmetry-elements object from the current selection",ttshown)

      ! Delete: reduce the crystal symmetry to the largest subgroup that excludes
      ! the selected operations (crystals only; the identity cannot be deleted)
      if (.not.sys(isys)%c%ismolecule) then
         call igAlignTextToFramePadding()
         call iw_text("Transform",highlight=.true.)

         nsel = 0
         do i = 2, min(nop,size(w%geometry_sym_sel,1))
            if (w%geometry_sym_sel(i)) nsel = nsel + 1
         end do
         if (iw_button("Delete##symseldelete",sameline=.true.,disabled=(nsel==0))) &
            iaction = iaction_sym_deleteops
         call iw_tooltip("Delete the selected operations and rebuild with a maximal "//&
            "subgroup that excludes them.",ttshown)
      end if

    end subroutine symop_display_and_buttons

    ! deallocate the cached symmetry data (operations table, selection, and the
    ! symmetry-vs-epsilon analysis arrays), so they are rebuilt on next use
    subroutine clear_sym_cache()

      if (allocated(w%geometry_sym_ops)) deallocate(w%geometry_sym_ops)
      if (allocated(w%geometry_sym_hm)) deallocate(w%geometry_sym_hm)
      if (allocated(w%geometry_sym_axes)) deallocate(w%geometry_sym_axes)
      if (allocated(w%geometry_sym_sel)) deallocate(w%geometry_sym_sel)
      w%geometry_sym_selgen = w%geometry_sym_selgen + 1
      if (allocated(w%geometry_sym_analyze_eps)) deallocate(w%geometry_sym_analyze_eps)
      if (allocated(w%geometry_sym_analyze_sym)) deallocate(w%geometry_sym_analyze_sym)
      if (allocated(w%geometry_sym_analyze_num)) deallocate(w%geometry_sym_analyze_num)

    end subroutine clear_sym_cache

    ! record the time the selection was last cleared (used to detect geometry
    ! changes that invalidate the system selection)
    subroutine clear_highlights_table()
      use interfaces_glfw, only: glfwGetTime

      w%timelast_geometry_clearhighlights = glfwGetTime()

    end subroutine clear_highlights_table

    ! center of mass of molecule imol, in the currently selected coordinate type
    function mol_com_coords(imol) result(x0)
      integer, intent(in) :: imol
      real*8 :: x0(3)

      x0 = sys(isys)%c%mol(imol)%cmass()
      if (w%geometry_moltype == atlisttype_ncel_frac) then
         x0 = sys(isys)%c%c2x(x0)
      else
         if (sys(isys)%c%ismolecule) x0 = x0 + sys(isys)%c%molx0
         if (w%geometry_moltype == atlisttype_ncel_ang) x0 = x0 * bohrtoa
      end if

    end function mol_com_coords

    ! Euler angles (ZYZ convention, in degrees) of the standard
    ! orientation of molecule imol, relative to the Cartesian axes.
    function mol_euler_angles(imol) result(eul)
      integer, intent(in) :: imol
      real*8 :: eul(3)

      call sys(isys)%c%mol(imol)%standard_axes(euler=eul)
      eul = eul * (180d0 / pi)

    end function mol_euler_angles

    ! calculate the w%iord to sort the table
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
      w%lastselected = 0

      ! reallocate the iord
      if (ntype <= 1) then
         if (allocated(w%iord)) deallocate(w%iord)
         allocate(w%iord(1))
         w%iord = 1
         return
      else
         if (allocated(w%iord)) deallocate(w%iord)
         allocate(w%iord(ntype),iperm(ntype))
         do i = 1, ntype
            w%iord(i) = i
            iperm(i) = i
         end do
      end if

      ! carry out the sort
      if (icolsort(w%geometry_sortcid)==ic_id .or. icolsort(w%geometry_sortcid)==ic_zat .or.&
         icolsort(w%geometry_sortcid)==ic_mol .or. icolsort(w%geometry_sortcid)==ic_mul .or.&
         icolsort(w%geometry_sortcid)==ic_idx .or. icolsort(w%geometry_sortcid)==ic_nat) then
         ! integers
         allocate(ival(ntype))
         do ii = 1, ntype
            i = w%iord(ii)
            if (table_hltype == atlisttype_nmol) then
               if (icolsort(w%geometry_sortcid) == ic_id) then
                  ival(ii) = i
               elseif (icolsort(w%geometry_sortcid) == ic_nat) then
                  ival(ii) = sys(isys)%c%mol(i)%nat
               elseif (icolsort(w%geometry_sortcid) == ic_idx) then
                  ival(ii) = sys(isys)%c%idxmol(i)
               end if
            else
               if (icolsort(w%geometry_sortcid) == ic_id) then
                  ival(ii) = i
               elseif (icolsort(w%geometry_sortcid) == ic_zat) then
                  ispc = sysc(isys)%attype_species(w%geometry_atomtype,i)
                  ival(ii) = sys(isys)%c%spc(ispc)%z
               elseif (icolsort(w%geometry_sortcid) == ic_mol) then
                  ival(ii) = sys(isys)%c%idatcelmol(1,i)
               elseif (icolsort(w%geometry_sortcid) == ic_mul) then
                  ival(ii) = sys(isys)%c%at(i)%mult
               elseif (icolsort(w%geometry_sortcid) == ic_idx) then
                  ival(ii) = sys(isys)%c%atcel(i)%idx
               end if
            end if
         end do
         call mergesort(ival,iperm,1,ntype)
         deallocate(ival)
      elseif (icolsort(w%geometry_sortcid)==ic_x.or.icolsort(w%geometry_sortcid)==ic_y.or.&
         icolsort(w%geometry_sortcid)==ic_z.or.icolsort(w%geometry_sortcid)==ic_ea.or.&
         icolsort(w%geometry_sortcid)==ic_eb.or.icolsort(w%geometry_sortcid)==ic_eg.or.&
         icolsort(w%geometry_sortcid)==ic_occ) then
         ! real
         allocate(rval(ntype))
         do ii = 1, ntype
            i = w%iord(ii)
            if (icolsort(w%geometry_sortcid)==ic_ea.or.icolsort(w%geometry_sortcid)==ic_eb.or.&
               icolsort(w%geometry_sortcid)==ic_eg) then
               ! Euler angles (molecules tab)
               x0 = mol_euler_angles(i)
               if (icolsort(w%geometry_sortcid) == ic_ea) then
                  rval(ii) = x0(1)
               elseif (icolsort(w%geometry_sortcid) == ic_eb) then
                  rval(ii) = x0(2)
               else
                  rval(ii) = x0(3)
               end if
            elseif (icolsort(w%geometry_sortcid) == ic_occ) then
               ! occupancy
               rval(ii) = sysc(isys)%attype_occupancy(w%geometry_atomtype,i)
            else
               if (table_hltype == atlisttype_nmol) then
                  x0 = mol_com_coords(i)
               else
                  x0 = sysc(isys)%attype_coordinates(w%geometry_atomtype,i)
               end if

               if (icolsort(w%geometry_sortcid) == ic_x) then
                  rval(ii) = x0(1)
               elseif (icolsort(w%geometry_sortcid) == ic_y) then
                  rval(ii) = x0(2)
               elseif (icolsort(w%geometry_sortcid) == ic_z) then
                  rval(ii) = x0(3)
               end if
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
               sval(ii)%s = sysc(isys)%attype_name(w%geometry_atomtype,i)
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

    ! Build the per-row selection state (rowstate = 0/1/2) and color (rowrgba)
    ! and the per-row cell-atom count (rowntot) for the current aggregate view
    ! (species or nneq), from the system-wide per-cell-atom selection
    ! (sysc%highlight_rgba). This is the only case that needs a full O(ncel)
    ! pass per frame, because species/nneq rows have no reverse map (row ->
    ! cell atoms). Cell and molecule views compute their state inline per
    ! visible row instead (see row_select_state). A row is 2 (selected) if all
    ! its cell atoms are highlighted, 1 (partial) if only some are, 0 otherwise.
    subroutine build_aggregate_state()
      integer :: i, j, nat
      integer, allocatable :: nselrow(:)

      ! (re)allocate the per-row arrays to the current number of rows
      if (allocated(rowstate)) then
         if (size(rowstate,1) /= ntype) then
            deallocate(rowstate)
            if (allocated(rowrgba)) deallocate(rowrgba)
            if (allocated(rowntot)) deallocate(rowntot)
            w%lastselected = 0
         end if
      end if
      if (.not.allocated(rowstate)) allocate(rowstate(ntype))
      if (.not.allocated(rowrgba)) allocate(rowrgba(4,ntype))
      if (.not.allocated(rowntot)) allocate(rowntot(ntype))
      rowstate = 0
      rowrgba = 0._c_float
      rowntot = 0

      ! nothing to do if there is no selection for this system
      nat = sys(isys)%c%ncel
      if (nat <= 0) return

      ! count selected vs. total cell atoms per row, and pick up the color.
      ! Each cell atom maps to its row in the current table type (O(ncel)).
      allocate(nselrow(ntype))
      nselrow = 0
      do j = 1, nat
         i = sysc(isys)%attype_celatom_to_id(table_hltype,j)
         if (i < 1 .or. i > ntype) cycle
         rowntot(i) = rowntot(i) + 1
         if (allocated(sysc(isys)%highlight_rgba)) then
            if (any(sysc(isys)%highlight_rgba(:,j) >= 0._c_float)) then
               nselrow(i) = nselrow(i) + 1
               rowrgba(:,i) = sysc(isys)%highlight_rgba(:,j)
            end if
         end if
      end do

      ! classify each row: 2 = fully selected, 1 = partially selected, 0 = unselected
      do i = 1, ntype
         if (rowntot(i) > 0 .and. nselrow(i) == rowntot(i)) then
            rowstate(i) = 2
         elseif (nselrow(i) > 0) then
            rowstate(i) = 1
         end if
      end do

    end subroutine build_aggregate_state

    ! deallocate the aggregate per-row arrays, so current_row_state falls back
    ! to the inline (cell/molecule) computation
    subroutine clear_aggregate_state()

      if (allocated(rowstate)) deallocate(rowstate)
      if (allocated(rowrgba)) deallocate(rowrgba)
      if (allocated(rowntot)) deallocate(rowntot)

    end subroutine clear_aggregate_state

    ! Compute the selection state (0/1/2) and highlight color of a single row
    ! id of the given table type, directly from sysc%highlight_rgba. Cell views
    ! are O(1); molecules iterate their fragment members (O(molecule size));
    ! species/nneq scan the cell atoms (O(ncel), only used for single-row
    ! queries such as a click, never per visible row in the draw loop).
    subroutine row_select_state(type,id,state,rgba)
      use systems, only: atlisttype_ncel_frac, atlisttype_ncel_bohr,&
         atlisttype_ncel_ang, atlisttype_nmol
      integer, intent(in) :: type, id
      integer, intent(out) :: state
      real(c_float), intent(out) :: rgba(4)

      integer :: j, k, cidx, nsel, ntot

      state = 0
      rgba = 0._c_float
      if (.not.allocated(sysc(isys)%highlight_rgba)) return

      if (type == atlisttype_ncel_frac .or. type == atlisttype_ncel_bohr .or.&
         type == atlisttype_ncel_ang) then
         ! cell views: the row is a single cell atom
         if (id >= 1 .and. id <= sys(isys)%c%ncel) then
            if (any(sysc(isys)%highlight_rgba(:,id) >= 0._c_float)) then
               state = 2
               rgba = sysc(isys)%highlight_rgba(:,id)
            end if
         end if
      elseif (type == atlisttype_nmol) then
         ! molecules: iterate the fragment members via their cell-atom index
         nsel = 0
         ntot = sys(isys)%c%mol(id)%nat
         do k = 1, ntot
            cidx = sys(isys)%c%mol(id)%at(k)%cidx
            if (cidx < 1 .or. cidx > sys(isys)%c%ncel) cycle
            if (any(sysc(isys)%highlight_rgba(:,cidx) >= 0._c_float)) then
               nsel = nsel + 1
               rgba = sysc(isys)%highlight_rgba(:,cidx)
            end if
         end do
         if (ntot > 0 .and. nsel == ntot) then
            state = 2
         elseif (nsel > 0) then
            state = 1
         end if
      else
         ! species/nneq: scan the cell atoms mapping to this row
         nsel = 0
         ntot = 0
         do j = 1, sys(isys)%c%ncel
            if (sysc(isys)%attype_celatom_to_id(type,j) /= id) cycle
            ntot = ntot + 1
            if (any(sysc(isys)%highlight_rgba(:,j) >= 0._c_float)) then
               nsel = nsel + 1
               rgba = sysc(isys)%highlight_rgba(:,j)
            end if
         end do
         if (ntot > 0 .and. nsel == ntot) then
            state = 2
         elseif (nsel > 0) then
            state = 1
         end if
      end if

    end subroutine row_select_state

    ! Selection state (0/1/2) and color of row i in the current table. Uses the
    ! per-frame aggregate arrays if present (species/nneq), otherwise computes
    ! inline for the current table type (cell/molecule).
    subroutine current_row_state(i,state,rgba)
      integer, intent(in) :: i
      integer, intent(out) :: state
      real(c_float), intent(out) :: rgba(4)

      if (allocated(rowstate)) then
         state = rowstate(i)
         rgba = rowrgba(:,i)
      else
         call row_select_state(table_hltype,i,state,rgba)
      end if

    end subroutine current_row_state

    ! get the sort specs from the imgui table
    subroutine fetch_sort_specs()

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

    end subroutine fetch_sort_specs

    ! get the view associated with the currently selected system
    subroutine get_current_view()

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

    end subroutine get_current_view

    ! color of atom iat (of atom-list type itype) from the first shown atoms
    ! representation in the current view (iview); returns .true. and fills rgbo
    ! when a color is found, .false. otherwise
    function color_from_view(itype,iat,rgbo) result(have)
      integer, intent(in) :: itype, iat
      real(c_float), intent(out) :: rgbo(3)
      logical :: have

      integer :: jrep, idd

      have = .false.
      rgbo = 0._c_float
      if (iview > 0) then
         do jrep = 1, win(iview)%sc%nrep
            if (win(iview)%sc%rep(jrep)%type == reptype_atoms.and.win(iview)%sc%rep(jrep)%isinit.and.&
               win(iview)%sc%rep(jrep)%shown) then
               idd = sysc(isys)%attype_type_id_to_id(itype,iat,win(iview)%sc%rep(jrep)%atoms%style%type)
               if (idd /= 0) then
                  have = .true.
                  rgbo = win(iview)%sc%rep(jrep)%atoms%style%rgb(:,idd)
                  exit
               end if
            end if
         end do
      end if

    end function color_from_view

    ! record a deferred bond-order change for the current bonds-tab bond
    ! (cell atom i, neighbor jj, neighbor-star entry j) to order val
    subroutine defer_setorder(val)
      integer, intent(in) :: val
      ibord1 = i
      ibord2 = jj
      lbord = sys(isys)%c%nstar(i)%lcon(:,j)
      ibordval = val
    end subroutine defer_setorder

    ! draw the row of buttons controlling the highlights
    subroutine draw_highlight_buttons()
      integer, allocatable :: tstate(:)

      ! highlight color
      call igAlignTextToFramePadding()
      call iw_text("Selection",highlight=.true.)
      call igSameLine(0._c_float,-1._c_float)
      ldum = iw_coloredit("##drawgeometryhighlightcolor",rgba=w%geometry_select_rgba)
      call iw_tooltip("Color used for highlighting atoms")

      ! Highlight buttons: all, none, toggle.
      if (iw_button("All##highlightall",sameline=.true.)) then
         allocate(ihigh(ntype),irgba(4,ntype))
         do i = 1, ntype
            ihigh(i) = i
            irgba(:,i) = w%geometry_select_rgba
         end do
         call sysc(isys)%highlight_atoms(.false.,ihigh,table_hltype,irgba)
         deallocate(ihigh,irgba)
      end if
      call iw_tooltip("Select all atoms in the system",ttshown)
      if (iw_button("None##highlightnone",sameline=.true.)) then
         call sysc(isys)%highlight_clear(.false.)
      end if
      call iw_tooltip("Deselect all atoms (" // trim(get_bind_keyname(BIND_EDITSELECT_DESELECT)) // ")",ttshown)
      if (iw_button("Toggle##highlighttoggle",sameline=.true.)) then
         ! compute the per-row selection state once
         allocate(tstate(ntype))
         do i = 1, ntype
            call current_row_state(i,tstate(i),rrgba)
         end do
         ! select the rows that are not fully selected
         nhigh = count(tstate /= 2)
         if (nhigh > 0) then
            allocate(ihigh(nhigh),irgba(4,nhigh))
            nhigh = 0
            do i = 1, ntype
               if (tstate(i) /= 2) then
                  nhigh = nhigh + 1
                  ihigh(nhigh) = i
                  irgba(:,nhigh) = w%geometry_select_rgba
               end if
            end do
            call sysc(isys)%highlight_atoms(.false.,ihigh,table_hltype,irgba)
            deallocate(ihigh,irgba)
         end if
         ! clear the rows that were fully selected
         nhigh = count(tstate == 2)
         if (nhigh > 0) then
            allocate(ihigh(nhigh))
            nhigh = 0
            do i = 1, ntype
               if (tstate(i) == 2) then
                  nhigh = nhigh + 1
                  ihigh(nhigh) = i
               end if
            end do
            call sysc(isys)%highlight_clear(.false.,ihigh,table_hltype)
            deallocate(ihigh)
         end if
         deallocate(tstate)
      end if
      call iw_tooltip("Toggle atomic selection",ttshown)

    end subroutine draw_highlight_buttons

    ! draw the contents of the add-atom popup: position and species
    ! selectors plus the Add button. Shared by the Add button in the edit
    ! row and the New row at the bottom of the atoms table.
    subroutine draw_addatom_popup()

      logical :: ldum, ok
      integer :: izout
      character(len=:), allocatable :: str

      ! the input position is interpreted in the same coordinate type as
      ! shown in the table, so the new atom's row matches the user input
      call igAlignTextToFramePadding()
      call iw_text("Position (" // sysc(isys)%attype_coordinates_units(w%geometry_atomtype) // ")")
      ldum = iw_dragfloat_real8("##xaddcoord",x1=w%geometry_input_coord(1),speed=0.001d0,decimal=6,&
         notlive=.true.,sameline=.true.)
      ldum = iw_dragfloat_real8("##yaddcoord",x1=w%geometry_input_coord(2),speed=0.001d0,decimal=6,&
         notlive=.true.,sameline=.true.)
      ldum = iw_dragfloat_real8("##zaddcoord",x1=w%geometry_input_coord(3),speed=0.001d0,decimal=6,&
         notlive=.true.,sameline=.true.)

      call igAlignTextToFramePadding()
      call iw_text("Species")
      if (w%geometry_input_species > 0) then
         str = string(w%geometry_input_species) // ": " // trim(sys(isys)%c%spc(w%geometry_input_species)%name)
      else
         str = trim(nameguess(abs(w%geometry_input_species),.true.))
      end if
      ldum = iw_button(str // "##speciesaddcoord",popupcontext=ok,&
         popupflags=ImGuiPopupFlags_MouseButtonLeft,sameline=.true.)
      if (ok) then
         ldum = iw_menuitem("Species ",enabled=.false.)
         call igSeparator()
         do j = 1, sys(isys)%c%nspc
            if (iw_menuitem(string(j) // ": " // trim(sys(isys)%c%spc(j)%name))) then
               w%geometry_input_species = j
            end if
         end do
         call igSeparator()
         str1 = "New" // c_null_char
         if (igBeginMenu(c_loc(str1),.true._c_bool)) then
            izout = iw_periodictable()
            if (izout >= 0) then
               w%geometry_input_species = -izout
               call igCloseCurrentPopup()
            end if
            call igEndMenu()
         end if
         call igEndPopup()
      end if

      if (iw_button("Add")) then
         iaction = iaction_add_atom
         iaction_i1 = w%geometry_input_species
         iaction_x = w%geometry_input_coord
         call igCloseCurrentPopup()
      end if

    end subroutine draw_addatom_popup

    ! draw the row of buttons controlling the edition of the system
    subroutine draw_edit_buttons()

      logical :: ldum, ok
      integer :: izout

      ! highlight color
      call igAlignTextToFramePadding()
      call iw_text("Edit",highlight=.true.)

      ! Add button: reset the input fields when the popup opens
      ldum = iw_button("Add##addatom",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
      if (ldum) then
         w%geometry_input_coord = 0d0
         w%geometry_input_species = 1
      end if
      if (ok) then
         if (w%geometry_atomtype == atlisttype_species) then
            izout = iw_periodictable()
            if (izout >= 0) then
               iaction = iaction_add_species
               iaction_i1 = izout
               call igCloseCurrentPopup()
            end if
         else
            call draw_addatom_popup()
         end if
         call igEndPopup()
      end if

      ! Duplicate button
      havesel = .false.
      if (allocated(sysc(isys)%highlight_rgba)) &
         havesel = any(sysc(isys)%highlight_rgba >= 0._c_float)
      if (iw_button("Duplicate##duplicateselection",sameline=.true.,disabled=.not.havesel)) then
         iaction = iaction_edit_highlighted
         iaction_i1 = edit_duplicate
      end if
      call iw_tooltip("Duplicate selected atoms",ttshown)

      ! Remove button
      if (iw_button("Remove##removeselection",sameline=.true.,disabled=.not.havesel)) then
         iaction = iaction_edit_highlighted
         iaction_i1 = edit_remove
      end if
      call iw_tooltip("Remove selected atoms (" // trim(get_bind_keyname(BIND_EDITSELECT_REMOVE)) // ")",ttshown)

      ! Merge button
      if (iw_button("Merge##mergeselection",sameline=.true.,disabled=.not.havesel)) then
         iaction = iaction_edit_highlighted
         iaction_i1 = edit_merge
      end if
      call iw_tooltip("Merge selected atoms",ttshown)

      ! Merge button
      if (iw_button("Reorder##reorderselection",sameline=.true.)) then
         iaction = iaction_reorder_highlighted
      end if
      call iw_tooltip("Relabel the atoms/species so their IDs are in the same order as shown in the table",ttshown)

    end subroutine draw_edit_buttons

    ! draw the row of buttons for editing the molecules (molecules tab)
    subroutine draw_mol_edit_buttons()

      ! Edit label
      call igAlignTextToFramePadding()
      call iw_text("Edit",highlight=.true.)

      ! is there a selection?
      havesel = .false.
      if (allocated(sysc(isys)%highlight_rgba)) &
         havesel = any(sysc(isys)%highlight_rgba >= 0._c_float)

      ! Remove button: remove the cell atoms of the selected molecules
      if (iw_button("Remove##removemol",sameline=.true.,disabled=.not.havesel)) then
         iaction = iaction_remove_molecules
      end if
      call iw_tooltip("Remove the selected molecules",ttshown)

      ! Reorder button: relabel so the molecule IDs follow the table order
      if (iw_button("Reorder##reordermol",sameline=.true.)) then
         iaction = iaction_reorder_molecules
      end if
      call iw_tooltip("Relabel the molecules so their IDs are in the same order as shown in the table",ttshown)

    end subroutine draw_mol_edit_buttons

    ! check whether the tab has changed: reset highlights and sort
    subroutine check_changed_tab(tab)

      character*(*) :: tab

      if (tab /= w%tabselected) then
         ! the selection persists across tabs (it lives on the system, indexed
         ! by cell atom); the window cache is rebuilt from it for the new tab
         call reset_sort()
         ! the shift-range anchor is per-table, so drop it when switching tabs
         w%lastselected = 0
         w%tabselected = tab
      end if

    end subroutine check_changed_tab

    ! add the selectable to a table row and process clicks
    subroutine process_selectable_clicks()

      ! the highlight selectable: hover and click
      clicked = .false.
      if (iw_highlight_selectable("##selectableintable" // suffix,clicked=clicked)) &
         ihighlight = i
      if (clicked) then
         ! implement selection range with shift and control
         if (igIsKeyDown(ImGuiKey_ModShift).and.w%lastselected /= 0.and.w%lastselected /= i) then
            ! selecte a whole range
            iclicked = -1
            iclicked_ini = min(ii,w%lastselected)
            iclicked_end = max(ii,w%lastselected)
         elseif (igIsKeyDown(ImGuiKey_ModCtrl)) then
            iclicked = -1
            iclicked_ini = ii
            iclicked_end = ii
         else
            iclicked = i
            w%lastselected = ii
         end if
      end if

    end subroutine process_selectable_clicks

  end subroutine draw_geometry

end submodule geometry
