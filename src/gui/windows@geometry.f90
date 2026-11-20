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

  !> Draw the geometry window.
  module subroutine draw_geometry(w)
    use representations, only: display_cell_window
    use display, only: rep_display
    use representations, only: reptype_symelem, repflavor_symelem, rep_shape, shapekind_arrow,&
       shapekind_cylinder, axes_length_def, axes_radius_def, axes_conelength_def,&
       axes_coneradius_def, rotaxis_radius_def, replat_primitive, shapekind_sphere, shapekind_box,&
       latvec_rad_def
    use crystalmod, only: symelem_type_mask, origin_atom, origin_cm, origin_center, origin_inversion
    use interfaces_glfw, only: glfwGetTime
    use crystalmod, only: holo_string, laue_string, pointgroup_info
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_OK_FOCUSED_DIALOG,&
       BIND_EDITSELECT_REMOVE, BIND_CANCEL
    use global, only: bondfactor_def, bonddelta_def
    use systems, only: nsys, sysc, sys, sys_init, ok_system, reread_system_from_file,&
       atlisttype_species, atlisttype_nneq, atlisttype_ncel_frac, atlisttype_ncel_bohr,&
       atlisttype_ncel_ang, atlisttype_nmol, celltransform_standard,&
       celltransform_primstd, celltransform_niggli, celltransform_delaunay
    use gui_main, only: g, ColorHighlightScene, ColorHighlightSelectScene, ColorHighlightBondScene,&
       ColorHighlightBondScene2, ColorTableHighlightRow, ColorAxes_def, ColorRotaxis_def
    use utils, only: iw_close_event, iw_table_headers_row, iw_text, iw_tooltip, iw_helpermark, iw_expression_helper, iw_exprctx_atoms, iw_calcwidth,&
       iw_button, iw_calcheight, iw_capture_text,&
       iw_atom_button, iw_combo_simple, iw_highlight_selectable, iw_coloredit, iw_dragfloat_real8,&
       iw_checkbox, iw_inputtext, iw_periodictable, iw_menuitem, iw_radiobutton, iw_intstepper,&
       iw_inputint, iw_inputint3, iw_icon_togglebutton, iw_setpos_bottomright, iw_table_column,&
       iw_beginmenu, iw_begintabitem, iw_separatortext, iw_separatortext_height
    use types, only: realloc, molsymop_rotation, molsymop_plane, molsymop_imp_rotation
    use tools_io, only: string, nameguess, ioj_center, ioj_right, isinteger, isreal
    use tools_math, only: det3, matinv, cellpar_from_metric
    use param, only: newline, bohrtoa, pi, atmcov0, maxzat0, eye
    class(window), intent(inout), target :: w

    logical :: domol, dowyc, doidx, docoord, havesel, haveexpr, doocc
    logical :: doquit, clicked, forcesort, ch, lch, deselected, chvol, iactive, syschanged
    integer :: ihighlight, iclicked, iclicked_ini, iclicked_end, nhigh, dec, icolsort(0:17)
    integer :: iwrite
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
    real(c_float) :: combowidth, rgb(3)
    integer :: ii, i, j, jj, isys, icol, isort, ispc, iz, izout, iview, im, jm, ncon
    integer :: ord, zi, zj ! bonds tab: bond order and atomic numbers of the two bonded atoms
    integer :: natused_bonds ! bonds tab: number of distinct atomic species present
    integer, allocatable :: iat_bonds(:) ! bonds tab: Z values of distinct species
    character(len=2), allocatable :: name_bonds(:) ! bonds tab: element symbols of distinct species
    logical :: atused_bonds(maxzat0) ! bonds tab: species presence flags
    real*8 :: dbond
    character(len=:), allocatable :: bondglyph, bondword ! bonds tab: bond-type glyph and word
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    logical :: havergb, havergb_, ldum, ok, ballow
    real*8 :: x0(3), x6(6), xold(3), x6old(6), res, vol, volold, scal
    real*8 :: occval, occold
    real*8 :: stdrot(3,3), stdcom(3), stdext, stdaxlen ! transient std-orientation axes
    integer :: ieuler_drag ! which Euler angle (1/2/3) is being dragged (0 = none)
    real*8 :: rotdir(3), rotlen ! transient rotation-axis direction and half-length
    type(rep_shape) :: stdshp(3) ! transient std-orientation axes (arrows)
    type(rep_shape) :: rotshp(1) ! transient rotation axis (cylinder)
    type(rep_shape) :: lvshp(7) ! transient preview of the new cell (spheres, arrows, box)
    logical :: lvshown ! what the cell tab pick needs was shown this frame (lattice points; the tab, for the origin)
    real*8 :: xb(3)
    real*8 :: pm(3,3), pt(3) ! cell tab: matrix and origin of the supercell being prepared
    logical :: pok ! cell tab: whether that transformation is valid
    character(len=:), allocatable :: pmsg ! cell tab: why it is not
    real*8 :: xcl(3), xca(3) ! cell tab: lengths and angles of the new cell
    integer :: ib, nsl, ncan ! cell tab: supercell candidate, selected size, and its number of sublattices
    ! preview of the new cell, in lattice-point sphere radii
    real*8, parameter :: pv_radfac = 1.5d0 ! the new lattice vectors' spheres
    real*8, parameter :: pv_arrowr = 0.3d0 ! arrow shaft radius
    real*8, parameter :: pv_headr = 2.5d0 ! arrowhead radius (shaft radii)
    real*8, parameter :: pv_headl = 0.15d0 ! arrowhead length (fraction of the arrow)
    real*8, parameter :: pv_edger = 0.25d0 ! outline radius
    real(c_float), parameter :: pv_edge_rgb(3) = (/1._c_float,0.6_c_float,0._c_float/) ! outline color
    real*8, parameter :: pv_molfac = 1d0 ! the spheres at the centers of the molecules
    real(c_float), parameter :: pv_mol_rgb(3) = (/0.70_c_float,0.93_c_float,0.55_c_float/) ! their color (bordered, as the lattice points)
    character(len=*), parameter :: errorigin = "The origin is not a lattice vector: pick another lattice point"
    ! preview of the plane where the cell is cut to set the vacuum
    type(rep_shape) :: cutshp(2) ! transient plane of the cut (a flat box) and outline of the new cell
    logical :: lvac ! the vacuum starts over: the axis was chosen, or the lattice changed elsewhere
    logical :: lcut ! the cut of the vacuum changed
    logical :: lwid ! a widget of the vacuum section changed
    logical :: lvacapply ! an edit of the vacuum is applied in this frame
    logical :: lvaccommit ! the value of the vacuum was entered or reset
    real*8 :: vacini ! the vacuum that Reset goes back to (bohr)
    character(len=:), allocatable :: vactip ! tooltip of the vacuum field
    real*8 :: vacwin(3) ! start of the cell in which the view draws the atoms
    real*8 :: vacx0(3), vacm(3,3), vachnew ! the new cell: origin, lattice vectors and height
    real*8 :: vacz ! the cut, from the start of the cell in which the atoms are drawn
    real*8 :: vacperp(3) ! vacuum along each axis, perpendicular to the other two
    integer :: ivj, ivk ! the two axes of the plane
    logical :: lmolcrys ! a molecular crystal: the molecules are kept whole unless told otherwise
    ! tags of the transient shapes of the vacuum: the plane and the
    ! outline, and the centers of the molecules kept whole (the other
    ! transient shapes of this window use 0, the positive atom ids and
    ! -5 or less)
    integer, parameter :: tag_vaccut = -1
    integer, parameter :: tag_vacmol = -2
    ! two vacuums closer than this (Å) show as the same value in the
    ! field (just over half its last digit, decimal=4)
    real*8, parameter :: vacsame = 0.6d-4
    real(c_float), parameter :: pv_cut_alpha = 0.3_c_float ! opacity of the plane
    integer :: istat, latpick(4) ! lattice-vector pick: status and the picked point
    integer :: atpick(4) ! origin pick: the picked atom
    real*8 :: xpick(3)
    type(ImVec4) :: col4
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
    character(len=:), allocatable :: iaction_str, experr
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
    integer, parameter :: iaction_shift_origin = 26
    integer, parameter :: iaction_set_vacuum = 27

    ! edit actions on highglighted atoms
    integer, parameter :: edit_remove = 1
    integer, parameter :: edit_merge = 2
    integer, parameter :: edit_duplicate = 3

    ! the tables of the window that can be written as text
    integer, parameter :: itable_atoms = 1
    integer, parameter :: itable_species = 2
    integer, parameter :: itable_cellnice = 3
    integer, parameter :: itable_molecules = 4
    integer, parameter :: itable_bonds = 5
    integer, parameter :: itable_symopsmol = 6
    integer, parameter :: itable_symops = 7
    integer, parameter :: itable_symcen = 8
    integer, parameter :: itable_symanal = 9
    integer, parameter :: itable_cellsub = 10

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
       call clear_highlights_table()
       call clear_empty_rows()
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
       w%sortcid = 0
       w%sortdir = 1
       w%geometry_input_coord = 0d0
       w%geometry_input_species = 1
       call w%geometry_addbond%clear()
       w%geometry_cell = cellform_state()
       ! the table of sublattices is long, and seldom needed
       w%geometry_sec(geomsec_sublat)%open = .false.
       call clear_sym_cache()
    end if

    ! This window edits the system shown in its anchor view, and commands that
    ! same view for its picks. If the view is gone, or it moved to a system that
    ! no longer exists, the window has nothing left to act on.
    if (.not.w%anchor(iview,isys,syschanged)) then
       call w%end()
       return
    end if
    if (syschanged .and. .not.w%firstpass) &
       call reset_for_new_system()
    if (.not.ok_system(isys,sys_init)) then
       call w%end()
       return
    end if

    ! handle an active bond operation commanded to the anchor view
    call bondmode_poll(w,iview,.true.)

    ! handle a pending add-bond pick commanded to a view window
    ipickhl = 0
    if (w%geometry_addbond%is_staged()) then
       ! the isys test matters for the frame after this window retargets its
       ! view: the picked index would belong to the system the view moved to
       ok = (win(iview)%isys == isys) .and. (win(iview)%viewmode_owner_is(w%id)) .and.&
          .not.w%geometry_addbond%is_stale(sysc(isys)%timelastchange_geometry)
       if (.not.ok) then
          ! the view moved to another system, another window took over the pick,
          ! or the geometry changed (stale cell-atom ids): cancel
          call win(iview)%viewmode_release_forced(w%id)
          call w%geometry_addbond%clear()
       elseif (win(iview)%viewmode >= 0) then
          ! the pick finished: add a single bond if a valid atom was clicked
          ! (self-bonds and duplicates are rejected by add_bond)
          if (win(iview)%vmdata%idx(1) > 0) &
             call sysc(isys)%add_bond(w%geometry_addbond%idx(1),win(iview)%vmdata%idx(1),&
                win(iview)%vmdata%idx(2:4),1)
          win(iview)%vmdata%idx = 0
          call w%geometry_addbond%clear()
       else
          ! the pick is in progress: highlight the atom receiving the bond
          ! (applied in the hover-highlight section at the end of the draw)
          ipickhl = w%geometry_addbond%idx(1)
       end if
    end if

    ! handle a pending lattice-vector pick (matrix cell transformation): the
    ! picked point is the integer part and the centering vector of the
    ! row; with Pick All, the pick moves on to the next rows. The origin
    ! is not a lattice vector: picking it keeps the pick waiting.
    lvshown = .false.
    if (w%geometry_lvpick_row == lvpick_origin) then
       ! the atom at the new origin
       call view_pick_result(iview,w%id,isys,w%geometry_lvpick,istat,xpick,atidx=atpick)
       if (istat == ipick_point .and. atpick(1) > 0) then
          iaction = iaction_shift_origin
          iaction_i1 = origin_atom
          iaction_i2 = atpick(1)
       end if
       if (istat /= ipick_pending) w%geometry_lvpick_row = 0
    elseif (w%geometry_lvpick_row > 0) then
       call view_pick_result(iview,w%id,isys,w%geometry_lvpick,istat,xpick,latidx=latpick)
       if (istat == ipick_point .and. latpick(4) > 0 .and. latpick(4) <= sys(isys)%c%ncv) then
          if (all(abs(latpick(1:3) + sys(isys)%c%cen(:,latpick(4))) < 1d-5)) then
             w%errmsg = errorigin
             call lvpick_start(w%geometry_lvpick_row,w%geometry_lvpick_last)
          else
             i = w%geometry_lvpick_row
             w%geometry_cell%intmat(:,i) = latpick(1:3)
             w%geometry_cell%cen(i) = latpick(4)
             w%errmsg = ""
             if (i < w%geometry_lvpick_last) then
                call lvpick_start(i+1,w%geometry_lvpick_last)
             else
                w%geometry_lvpick_row = 0
             end if
          end if
       elseif (istat /= ipick_pending) then
          w%geometry_lvpick_row = 0
          if (w%errmsg == errorigin) w%errmsg = ""
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
    if (w%timelast_sort < sysc(isys)%timelastchange_rebond) then
       call reset_sort()
    end if
    if (w%timelast_geometry_clearhighlights < sysc(isys)%timelastchange_geometry) then
       ! the geometry changed: cell-atom indices are stale, so clear the
       ! system selection and force the window cache to be rebuilt
       call sysc(isys)%highlight_clear(.false.)
       call clear_highlights_table()
       call clear_empty_rows()
       ! the symmetry-vs-epsilon analysis and the supercell search
       ! (its symmetry counts depend on the operations and the atoms)
       ! are now stale
       call clear_sym_cache()
       call clear_search_results()
    end if

    ! system combo. The atoms tab opens by default and whenever another
    ! window asks for it (geometry_seltab, cleared after the tab bar)
    if (w%firstpass .or. w%geometry_seltab == geomtab_atoms) then
       atompreflags = ImGuiTabItemFlags_SetSelected
    else
       atompreflags = ImGuiTabItemFlags_None
    end if
    call iw_text("System",highlight=.true.,alignframe=.true.)
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
                ! move the view this window is anchored to, so that the two can
                ! never disagree about which system is being edited
                call w%retarget(i)
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
    ldum = iw_checkbox("Keep bonding",w%geometry_keepbonding,labelleft=.true.,gap=.true.)
    call iw_tooltip("Preserve bonding when the system geometry changes",ttshown)

    ! show the tabs
    str1 = "##drawgeometry_tabbar" // c_null_char
    flags = ImGuiTabBarFlags_None
    call igBeginGroup()
    if (igBeginTabBar(c_loc(str1),flags)) then
       !! species tab !!
       if (iw_begintabitem("Species##drawgeometry_speciestab")) then
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
          sz0%y = sz0%y - iw_calcheight(2,0,.true.)
          if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
             icol = -1

             ! header setup
             call iw_table_column("Id",icol=icol,sortid=ic_id,icolsort=icolsort)

             call iw_table_column("Atom",icol=icol,sortid=ic_atom,icolsort=icolsort)

             call iw_table_column("Z ",icol=icol,sortid=ic_zat,icolsort=icolsort)

             call iw_table_column("Num. atoms",icol=icol,sortid=ic_nat,icolsort=icolsort)

             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! fetch the sort specs, sort the data if necessary
             call fetch_sort_specs()

             ! draw the header
             call iw_table_headers_row(autofit=.true.,writemenu=iwrite,ttshown=ttshown)
             call w%table_write_begin(itable_species,iwrite)

             ! sort
             if (forcesort) call table_sort()

             ! calculate the number of digits for output
             ndigit = ceiling(log10(ntype+0.1d0))
             ndigitm = 0
             ndigitidx = 0

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
                   call iw_text(string(i,ndigit),alignframe=.true.)
                end if
                ! emit even when column 0 is clipped off-screen
                call process_selectable_clicks()

                ! name
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   if (havergb) then
                      ldum = iw_coloredit("##tablecolorg" // suffix,rgb=rgb,nointeraction=.true.)
                   end if
                   if (iw_inputtext("##nametextinput" // string(i),bufsize=11,texta=name,width=max(3,len(name)),&
                      sameline=havergb)) then
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
                      if (sys(isys)%c%nspc > 1) then
                         ldum = iw_menuitem("Move the atoms to",enabled=.false.)
                         call igSeparator()
                         do j = 1, sys(isys)%c%nspc
                            if (j == i) cycle
                            if (iw_menuitem(string(j) // ": " // trim(sys(isys)%c%spc(j)%name))) then
                               iaction = iaction_set_attype_species
                               iaction_i1 = i
                               iaction_i2 = j
                            end if
                         end do
                         call igSeparator()
                      end if
                      if (iw_beginmenu("Change the element to")) then
                         izout = iw_periodictable()
                         if (izout >= 0) then
                            iaction = iaction_set_atomic_number
                            iaction_i1 = i
                            iaction_i2 = izout
                            call igCloseCurrentPopup()
                         end if
                         call igEndMenu()
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

             call w%table_write_end(itable_species,table_title("Species",isys))

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
       if (iw_begintabitem("Atoms##drawgeometry_atomstab",flags=atompreflags)) then
          ! check if the tab changed
          call check_changed_tab("atoms")

          ! group atom types
          call iw_text("Types",alignframe=.true.)
          if (sysc(isys)%attype_combo_simple("##atomtypeselectgeom",w%geometry_atomtype,atlisttype_allowed,&
             sameline=.true.)) then
             call reset_sort()
             call w%table_write_cancel()
          end if
          call iw_tooltip("Group atoms by these categories",ttshown)
          table_hltype = w%geometry_atomtype
          ntype = sysc(isys)%attype_number(w%geometry_atomtype)

          if (w%geometry_atomtype == atlisttype_nneq) then
             ldum = iw_checkbox("Keep symmetry",w%geometry_forcewyc,labelleft=.true.,gap=.true.)
             call iw_tooltip("If checked, changes to the atomic coordinates force the atom to maintain &
                &its current site symmetry (Wyckoff position)",ttshown)
          end if

          if (.not.sys(isys)%c%ismolecule) then
             call iw_text("Spg",alignframe=.true.,gap=.true.)
             if (sys(isys)%c%spgavail) then
                call iw_text(trim(sys(isys)%c%spg%international_symbol),sameline=.true.)
             else
                call iw_text("n/a",sameline=.true.)
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
          sz0%y = sz0%y - iw_calcheight(3,0,.true.)
          if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
             icol = -1

             ! header setup
             call iw_table_column("Id",icol=icol,sortid=ic_id,icolsort=icolsort)

             call iw_table_column("Atom",icol=icol,sortid=ic_atom,icolsort=icolsort)

             call iw_table_column("Z ",icol=icol,sortid=ic_zat,icolsort=icolsort)

             if (domol) then
                call iw_table_column("Mol",icol=icol,sortid=ic_mol,icolsort=icolsort)
             end if

             if (dowyc) then
                if (sys(isys)%c%havesym > 0 .and. sys(isys)%c%spgavail) then
                   str2 = "Wyc"
                   isort = ic_wyc
                else
                   str2 = "Mul"
                   isort = ic_mul
                end if
                call iw_table_column(str2,icol=icol,sortid=isort,icolsort=icolsort)
             end if

             if (doidx) then
                call iw_table_column("Idx",icol=icol,sortid=ic_idx,icolsort=icolsort)
             end if

             if (doocc) then
                call iw_table_column("Occ.",icol=icol,sortid=ic_occ,icolsort=icolsort)
             end if

             if (docoord) then
                if (w%geometry_atomtype == atlisttype_ncel_ang .or. w%geometry_atomtype == atlisttype_ncel_bohr) then
                   str2 = "/" // sysc(isys)%attype_coordinates_units(w%geometry_atomtype)
                   strx = "x" // str2
                   stry = "y" // str2
                   strz = "z" // str2
                else
                   strx = "x"
                   stry = "y"
                   strz = "z"
                end if
                call iw_table_column(strx,icol=icol,sortid=ic_x,icolsort=icolsort)

                call iw_table_column(stry,icol=icol,sortid=ic_y,icolsort=icolsort)

                call iw_table_column(strz,icol=icol,sortid=ic_z,icolsort=icolsort)
             end if

             if (haveexpr) then
                call iw_table_column("Expression",icol=icol,sortid=ic_expr,icolsort=icolsort,flags=ImGuiTableColumnFlags_NoSort)
             end if
             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! fetch the sort specs, sort the data if necessary
             call fetch_sort_specs()

             ! draw the header
             call iw_table_headers_row(autofit=.true.,writemenu=iwrite,ttshown=ttshown)

             ! sort
             if (forcesort) then
                call table_sort()
             end if

             ! start the clipper; if the table is being written as text,
             ! its rows are captured as they are drawn (all of them)
             clipper = ImGuiListClipper_ImGuiListClipper()
             call ImGuiListClipper_Begin(clipper,ntype,-1._c_float)
             call w%table_write_begin(itable_atoms,iwrite,clipper)

             ! calculate the number of digits for output
             ndigit = ceiling(log10(ntype+0.1d0))
             ndigitm = 0
             ndigitidx = 0
             if (domol) ndigitm = ceiling(log10(sys(isys)%c%nmol+0.1d0))
             if (doidx) ndigitidx = ceiling(log10(sys(isys)%c%nneq+0.1d0))

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
                      end if
                      if (iw_inputtext("##nametextinput" // suffix,bufsize=11,texta=name,width=max(3,len(name)),&
                         sameline=havergb)) then
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
                         if (iw_beginmenu("New")) then
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
             call ImGuiListClipper_destroy(clipper)
             call w%table_write_end(itable_atoms,table_title("Atoms",isys))

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
          call iw_text("Expression",alignframe=.true.)

          ! filter text input
          call iw_helpermark("Examples:"//newline//&
             "- '@dnuc:2' = distance to cell atom 2"//newline//&
             "- 'log($0)' = log of the promolecular density"//newline//&
             "- 'abs(@x) < 2 && abs(@y) < 2 && abs(@z) < 2' = atoms in the (-2,2) box"//newline//&
             "Use the Insert button for fields, functions, and examples.")
          if (iw_inputtext("##filtertext",bufsize=1023,width=30,texta=w%geometry_expression,&
             notlive=.true.,sameline=.true.)) then
             w%geometry_expression_ok = .true.
             w%geometry_expr_error = ""
          end if
          call iw_tooltip("Show on the table the result of using this expression at the atomic positions.",ttshown)
          if (iw_expression_helper("##helpfilter",isys,iw_exprctx_atoms,texta=w%geometry_expression,&
             ttshown=ttshown)) then
             ! an expression built by inserting is applied when it is
             ! valid; until then, show what is wrong with it
             call sys(isys)%check_expression(w%geometry_expression,experr)
             w%geometry_expr_error = experr
             w%geometry_expression_ok = (len_trim(experr) == 0)
          end if
          if (iw_button("Clear",sameline=.true.)) then
             w%geometry_expression_ok = .false.
             w%geometry_expression = ""
             w%geometry_expr_error = ""
          end if
          call iw_tooltip("Clear the expression",ttshown)
          if (len_trim(w%geometry_expr_error) > 0) then
             call iw_text("Error: " // trim(w%geometry_expr_error),danger=.true.,wrap=.true.)
          end if
          call igEndTabItem()
       end if

       !! cell tab !!
       if (.not.sys(isys)%c%ismolecule) then
          flags = ImGuiTabItemFlags_None
          if (w%geometry_seltab == geomtab_cell) flags = ImGuiTabItemFlags_SetSelected
          if (iw_begintabitem("Cell##drawgeometry_celltab",flags=flags)) then
             ! check if the tab changed
             call check_changed_tab("cell")

             ! while this tab is open, the view shows the whole cell in
             ! the vacuum directions: the cell cut around the slab looks
             ! the same whatever the vacuum and the cell length (not in
             ! the frame after this window moves its view to another system)
             if (associated(win(iview)%sc) .and. win(iview)%isys == isys) call win(iview)%sc%show_full_cell()

             call iw_separatortext("Lattice parameters",w%geometry_sec(geomsec_latpar))
             if (w%geometry_sec(geomsec_latpar)%open) then
                ! keep symmetry and space group
                ldum = iw_checkbox("Keep symmetry",w%geometry_forcewyc,labelleft=.true.)
                call iw_tooltip("If checked, force system to maintain the current space group",ttshown)
                call iw_text("Spg",alignframe=.true.,gap=.true.)
                if (sys(isys)%c%spgavail) then
                   call pointgroup_info(sys(isys)%c%spg%pointgroup_symbol,schpg,holo,laue)
                   call iw_text(trim(sys(isys)%c%spg%international_symbol) // " (" //&
                      string(holo_string(holo)) // ")",sameline=.true.)
                else
                   call iw_text("n/a",sameline=.true.)
                end if

                x6(1:3) = sys(isys)%c%aa
                x6(4:6) = sys(isys)%c%bb
                x6old = x6

                call iw_text("a/b/c (Å)",alignframe=.true.)
                ch = .false.
                lch = iw_dragfloat_real8("##celllengthsa",x1=x6(1),speed=0.005d0,decimal=6,scale=bohrtoa,&
                   notlive=.true.,sameline=.true.)
                ch = ch .or. lch
                lch = iw_dragfloat_real8("##celllengthsb",x1=x6(2),speed=0.005d0,decimal=6,scale=bohrtoa,&
                   notlive=.true.,sameline=.true.)
                ch = ch .or. lch
                lch = iw_dragfloat_real8("##celllengthsc",x1=x6(3),speed=0.005d0,decimal=6,scale=bohrtoa,&
                   notlive=.true.,sameline=.true.)
                ch = ch .or. lch

                call iw_text("α/β/γ (°)",alignframe=.true.)
                lch = iw_dragfloat_real8("##cellangsa",x1=x6(4),speed=0.01d0,decimal=4,sameline=.true.,&
                   notlive=.true.)
                ch = ch .or. lch
                lch = iw_dragfloat_real8("##cellanbsb",x1=x6(5),speed=0.01d0,decimal=4,sameline=.true.,&
                   notlive=.true.)
                ch = ch .or. lch
                lch = iw_dragfloat_real8("##cellanbsc",x1=x6(6),speed=0.01d0,decimal=4,sameline=.true.,&
                   notlive=.true.)
                ch = ch .or. lch

                ! cell volume
                call iw_text("Volume (Å³)",alignframe=.true.)
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
             end if

             ! standardization and reduction of the cell
             call iw_separatortext("Standard cells",w%geometry_sec(geomsec_celltr),hasmore=.true.)
             if (w%geometry_sec(geomsec_celltr)%open) then
                if (iw_button("Conventional##celltransfstd")) then
                   iaction = iaction_transform_cell
                   iaction_i1 = celltransform_standard
                   iaction_l = w%geometry_cell%refine
                end if
                call iw_tooltip("Transform to the conventional cell",ttshown)
                if (iw_button("Primitive##celltransfprim",sameline=.true.)) then
                   iaction = iaction_transform_cell
                   iaction_i1 = celltransform_primstd
                   iaction_l = w%geometry_cell%refine
                end if
                call iw_tooltip("Transform to the standard primitive cell",ttshown)
                if (iw_button("Niggli##celltransfnig",sameline=.true.)) then
                   iaction = iaction_transform_cell
                   iaction_i1 = celltransform_niggli
                   iaction_l = .false.
                end if
                call iw_tooltip("Transform to the primitive Niggli cell",ttshown)
                if (iw_button("Delaunay##celltransfdel",sameline=.true.)) then
                   iaction = iaction_transform_cell
                   iaction_i1 = celltransform_delaunay
                   iaction_l = .false.
                end if
                call iw_tooltip("Transform to the primitive Delaunay cell",ttshown)
                ! (an option that is on stays in view when More is closed)
                if (w%geometry_sec(geomsec_celltr)%more .or. w%geometry_cell%refine) then
                   ldum = iw_checkbox("Refine positions##cellrefine",w%geometry_cell%refine,labelleft=.true.)
                   call iw_tooltip("Move the atoms to their ideal symmetry positions when transforming&
                      & to the conventional or the primitive cell",ttshown)
                end if
             end if

             ! the supercell: a simple repetition, a general matrix, or a
             ! sublattice chosen from the supercell search
             call iw_separatortext("Supercell",w%geometry_sec(geomsec_super),&
                hasmore=(w%geometry_cell%mode == cellmode_matrix))
             if (w%geometry_sec(geomsec_super)%open) then
                call iw_text("Mode",alignframe=.true.)
                ldum = iw_radiobutton("Simple",int=w%geometry_cell%mode,intval=cellmode_simple,sameline=.true.)
                call iw_tooltip("Build a supercell by repeating the cell an integer number of&
                   & times along each lattice vector",ttshown)
                ldum = iw_radiobutton("Matrix",int=w%geometry_cell%mode,intval=cellmode_matrix,sameline=.true.)
                call iw_tooltip("General transformation: each row is a lattice vector of the new&
                   & cell, written in crystallographic coordinates of the current cell",ttshown)
                ldum = iw_radiobutton("Search",int=w%geometry_cell%mode,intval=cellmode_search,sameline=.true.)
                call iw_tooltip("Search the supercells of each size and choose one of them: the nicest&
                   & (the closest to a cube) or the one with the longest reach",ttshown)
                ldum = iw_radiobutton("Orient",int=w%geometry_cell%mode,intval=cellmode_orient,sameline=.true.)
                call iw_tooltip("Find a cell with two axes in a family of lattice planes and the third&
                   & perpendicular to them, or with one axis along a lattice direction",ttshown)

                if (w%geometry_cell%mode == cellmode_simple) then
                   ! simple transformation: integer multiples of the a, b, c axes
                   call iw_text("na/nb/nc",alignframe=.true.)
                   do jm = 1, 3
                      ldum = iw_intstepper("cellnrep" // string(jm),w%geometry_cell%nrep(jm),&
                         minval=1_c_int,sameline=.true.,&
                         tooltip="Number of times the cell is repeated along the a, b, and c lattice vectors")
                   end do
                elseif (w%geometry_cell%mode == cellmode_matrix) then
                   ! arbitrary transformation matrix: each row is a lattice vector of the new
                   ! cell = integer triplet + (for centered cells) a centering vector
                   w%geometry_cell%cen = min(max(w%geometry_cell%cen,1),max(sys(isys)%c%ncv,1))
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
                      if (im == 1) then
                         call iw_text("a' ",alignframe=.true.)
                      elseif (im == 2) then
                         call iw_text("b' ",alignframe=.true.)
                      else
                         call iw_text("c' ",alignframe=.true.)
                      end if
                      lch = iw_inputint3("##cellmat" // string(im),&
                         w%geometry_cell%intmat(:,im),width=3*3,sameline=.true.,notlive=.true.)
                      ch = ch .or. lch
                      if (sys(isys)%c%ncv > 1) then
                         call iw_text(" + ",sameline=.true.)
                         call iw_combo_simple("##cellcen" // string(im),stropt,w%geometry_cell%cen(im),&
                            changed=lch,sameline=.true.,startsatone=.true.)
                         ch = ch .or. lch
                         call iw_tooltip("Centering vector added to this lattice vector",ttshown)
                      end if
                      ! no picks with the inverse: the points in the view are
                      ! lattice vectors of the new cell, not rows of its inverse
                      if (iw_button("Pick##cellpick" // string(im),sameline=.true.,&
                         disabled=(w%geometry_lvpick_row == im .or. w%geometry_cell%inverse))) &
                         call lvpick_start(im,im)
                      call iw_tooltip("Pick this lattice vector in the view: click the lattice point &
                         &it ends at (a click anywhere else cancels). Not available with the inverse&
                         & matrix",ttshown,whendisabled=.true.)
                   end do
                   if (iw_button("Pick All##cellpickall",disabled=((w%geometry_lvpick_row > 0 .and.&
                      w%geometry_lvpick_row /= lvpick_origin) .or.&
                      w%geometry_cell%inverse))) call lvpick_start(1,3)
                   call iw_tooltip("Pick the three lattice vectors in the view, one after the other &
                      &(a click anywhere else than a lattice point cancels). Not available with the&
                      & inverse matrix",ttshown,whendisabled=.true.)
                   call iw_text("Origin",alignframe=.true.)
                   do jm = 1, 3
                      lch = iw_dragfloat_real8("##cellorigin" // string(jm),&
                         x1=w%geometry_cell%origin(jm),speed=0.01d0,decimal=4,sameline=.true.,&
                         notlive=.true.)
                      ch = ch .or. lch
                   end do
                   call iw_tooltip("Origin of the new cell, in crystallographic coordinates of the&
                      & current cell",ttshown)
                   if (w%geometry_sec(geomsec_super)%more .or. w%geometry_cell%inverse) then
                      lch = iw_checkbox("Inverse##cellinverse",w%geometry_cell%inverse,labelleft=.true.)
                      ch = ch .or. lch
                      call iw_tooltip("Use the inverse of the matrix: the rows are then the lattice vectors&
                         & of the current cell written in crystallographic coordinates of the new cell",ttshown)
                   end if
                   ! if the user changed the matrix or the origin, clear the error message
                   if (ch) w%errmsg = ""
                elseif (w%geometry_cell%mode == cellmode_orient) then
                   ! a plane or a direction, and the axis of the new cell that goes with it
                   call iw_text("Type",alignframe=.true.)
                   ch = iw_radiobutton("Plane (hkl)",bool=w%geometry_cell%orient_dir,boolval=.false.,sameline=.true.)
                   call iw_tooltip("The indices are the Miller indices of a family of lattice planes. Two&
                      & axes of the new cell are the shortest lattice vectors in the plane, and the third&
                      & is the shortest lattice vector perpendicular to it. If the lattice has none, the&
                      & third axis goes to the next plane, as close to the normal as possible",ttshown)
                   lch = iw_radiobutton("Direction [uvw]",bool=w%geometry_cell%orient_dir,boolval=.true.,&
                      sameline=.true.)
                   ch = ch .or. lch
                   call iw_tooltip("The indices are a lattice direction. One axis of the new cell is the&
                      & shortest lattice vector along it, and the other two are the shortest lattice&
                      & vectors that complete the cell",ttshown)
                   call iw_text("Indices",alignframe=.true.)
                   lch = iw_inputint3("##cellorientidx",w%geometry_cell%orient_ivec,width=3*3,sameline=.true.)
                   ch = ch .or. lch
                   call iw_tooltip("Indices of the plane or the direction, referred to the current cell",ttshown)
                   call iw_text("Axis",alignframe=.true.,gap=.true.)
                   call iw_combo_simple("##cellorientaxis","a" // c_null_char // "b" // c_null_char // "c" //&
                      c_null_char,w%geometry_cell%orient_axis,changed=lch,sameline=.true.,startsatone=.true.)
                   ch = ch .or. lch
                   call iw_tooltip("Axis of the new cell that goes out of the plane, or along the&
                      & direction",ttshown)
                   if (.not.w%geometry_cell%orient_dir) then
                      lch = iw_checkbox("Primitive cell",w%geometry_cell%orient_prim,labelleft=.true.,gap=.true.)
                      ch = ch .or. lch
                      call iw_tooltip("Make the new cell primitive: the third axis goes to the next plane,&
                         & as close to the normal as possible, even if the lattice has a vector&
                         & perpendicular to the plane",ttshown)
                   end if
                   if (ch) w%errmsg = ""
                else
                   ! supercell search: the criterion and the largest size
                   call iw_text("Criterion",alignframe=.true.)
                   call iw_combo_simple("##cellsearchcrit","Niceness" // c_null_char // "Reach" // c_null_char,&
                      w%geometry_cell%crit,changed=lch,sameline=.true.,startsatone=.true.,tooltips=&
                      "The supercell that contains the largest sphere, which is the closest to a cube" //&
                      c_null_char //&
                      "The supercell with the longest shortest lattice vector: two atoms closer than half&
                      & of it are never periodic images of each other" // c_null_char,ttshown=ttshown)
                   ! the criterion changed: the results would no longer match it
                   if (lch) call clear_search_results()
                   if (iw_checkbox("Fewest displacements##cellsearchmindisp",w%geometry_cell%mindisp,&
                      labelleft=.true.,gap=.true.)) call clear_search_results()
                   call iw_tooltip("Choose, for each size, only among the supercells that need the fewest&
                      & displaced structures in a phonon calculation. They keep the most symmetry, but&
                      & can be very elongated",ttshown)
                   call iw_text("Max. size",alignframe=.true.)
                   ldum = iw_intstepper("cellsearchsize",w%geometry_cell%nmax,minval=1_c_int,sameline=.true.,&
                      tooltip="Largest supercell size (number of current cells) considered in the search")
                   if (iw_button("Search##cellsearch",sameline=.true.)) then
                      w%geometry_cell%nmax = max(w%geometry_cell%nmax,1_c_int)
                      call clear_search_results()
                      w%errmsg = ""
                      call sysc(isys)%cell_nice_list(int(w%geometry_cell%nmax),&
                         w%geometry_cell%crit == cellcrit_reach,w%geometry_cell%mindisp,&
                         w%geometry_cell%nc,w%geometry_cell%ic0,w%geometry_cell%cand,&
                         w%geometry_cell%ibest,w%errmsg)
                   end if
                   call iw_tooltip("Search the supercells of every size up to the maximum. The cost&
                      & grows quickly with the size",ttshown)

                   ! the supercell chosen for each size: click a row to select it
                   if (allocated(w%geometry_cell%cand)) then
                      if (cell_table_begin("##cellsearchtable","n","Size of the supercell: the number of&
                         & current cells it contains",count(w%geometry_cell%ibest > 0),10)) then
                         call w%table_write_begin(itable_cellnice,iwrite)
                         do i = 1, size(w%geometry_cell%ibest,1)
                            ib = w%geometry_cell%ibest(i)
                            if (ib == 0) cycle
                            call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                            if (igTableSetColumnIndex(0)) call iw_text(string(i))
                            ldum = iw_highlight_selectable("##cellsearchrow" // string(i),clicked=clicked,&
                               selected=(w%geometry_cell%nsel == i))
                            if (clicked) call select_candidate(ib,i)
                            call cell_candidate_columns(ib,i)
                         end do
                         if (w%geometry_cell%crit == cellcrit_reach) then
                            call w%table_write_end(itable_cellnice,table_title("Longest-reach supercells",isys))
                         else
                            call w%table_write_end(itable_cellnice,table_title("Nicest supercells",isys))
                         end if
                         call igEndTable()
                      end if
                   end if
                end if
             end if

             ! the transformation the supercell section describes now
             call pending_cell(pm,pt,pok,pmsg)

             if (w%geometry_sec(geomsec_super)%open) then
                ! the cell the transformation leads to
                call iw_text("New cell:",alignframe=.true.)
                if (pok) then
                   call pending_metrics(pm,s,xcl,xca)
                   call iw_text(s,sameline=.true.)
                   call iw_text("a/b/c (Å):",alignframe=.true.)
                   call iw_text(string(xcl(1)*bohrtoa,'f',decimal=4) // "  " // string(xcl(2)*bohrtoa,'f',decimal=4) //&
                      "  " // string(xcl(3)*bohrtoa,'f',decimal=4),sameline=.true.)
                   call iw_text("α/β/γ (°):",alignframe=.true.,gap=.true.)
                   call iw_text(string(xca(1),'f',decimal=2) // "  " // string(xca(2),'f',decimal=2) //&
                      "  " // string(xca(3),'f',decimal=2),sameline=.true.)
                else
                   call iw_text(pmsg,disabled=.true.,sameline=.true.)
                end if

                ok = iw_button("Apply##celltransfapply",danger=.true.,disabled=.not.pok)
                call iw_tooltip("Transform the crystal to the new cell",ttshown)
                if (ok) then
                   iaction = iaction_transform_matrix
                   iaction_m = pm
                   iaction_x = pt
                   iaction_l = .false.
                end if
                if (iw_button("Reset##cellresettransf",sameline=.true.)) then
                   if (w%geometry_cell%mode == cellmode_simple) then
                      w%geometry_cell%nrep = 1_c_int
                   elseif (w%geometry_cell%mode == cellmode_matrix) then
                      w%geometry_cell%intmat = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
                      w%geometry_cell%cen = 1
                      w%geometry_cell%origin = 0d0
                      w%geometry_cell%inverse = .false.
                   elseif (w%geometry_cell%mode == cellmode_orient) then
                      call reset_orient_inputs()
                   else
                      w%geometry_cell%isel = 0
                      w%geometry_cell%nsel = 0
                   end if
                   w%errmsg = ""
                end if
                call iw_tooltip("Reset the supercell transformation",ttshown)
             end if

             ! every sublattice of the size selected in the search, to
             ! choose another one than the search did
             if (w%geometry_cell%mode == cellmode_search .and. allocated(w%geometry_cell%cand) .and.&
                w%geometry_cell%nsel > 0) then
                nsl = w%geometry_cell%nsel
                ncan = w%geometry_cell%nc(nsl)
                call iw_separatortext("Sublattices of size " // string(nsl) // " (" // string(ncan) // ")",&
                   w%geometry_sec(geomsec_sublat))
                if (w%geometry_sec(geomsec_sublat)%open) then
                   if (cell_table_begin("##cellsublattable","#","Distinct sublattices of this size, the best&
                      & first; the asterisk marks the one chosen by the search",ncan,6)) then
                      ! start the clipper; if the table is being written as text,
                      ! its rows are captured as they are drawn (all of them)
                      clipper = ImGuiListClipper_ImGuiListClipper()
                      call ImGuiListClipper_Begin(clipper,ncan,-1._c_float)
                      call w%table_write_begin(itable_cellsub,iwrite,clipper)
                      do while(ImGuiListClipper_Step(clipper))
                         call c_f_pointer(clipper,clipper_f)
                         do ii = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                            ib = w%geometry_cell%ic0(nsl) + ii
                            ! the displacements are counted when a row is first drawn
                            if (w%geometry_cell%cand(ib)%ndisp == -1) then
                               experr = ""
                               call sysc(isys)%cell_nice_count(w%geometry_cell%cand(ib:ib),experr)
                               if (len_trim(experr) > 0) w%errmsg = experr
                            end if
                            call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                            if (igTableSetColumnIndex(0)) then
                               if (ib == w%geometry_cell%ibest(nsl)) then
                                  call iw_text(string(ii) // "*")
                               else
                                  call iw_text(string(ii))
                               end if
                            end if
                            ldum = iw_highlight_selectable("##cellsublatrow" // string(ii),clicked=clicked,&
                               selected=(w%geometry_cell%isel == ib))
                            if (clicked) call select_candidate(ib,nsl)
                            call cell_candidate_columns(ib,nsl)
                         end do
                      end do
                      call ImGuiListClipper_End(clipper)
                      call ImGuiListClipper_destroy(clipper)
                      call w%table_write_end(itable_cellsub,table_title("Sublattices of size " // string(nsl),isys))
                      call igEndTable()
                   end if
                end if
             end if

             ! the vacuum along an axis: the cell is cut and made longer
             ! or shorter there, at constant atomic positions
             associate (gc => w%geometry_cell)
               lvac = .false.
               lcut = .false.
               lwid = .false.
               lvacapply = .false.

               ! first time with this system: the axis with the most vacuum,
               ! or c, and the molecules whole in a molecular crystal
               lmolcrys = (sys(isys)%c%ismol3d .and. sys(isys)%c%nmol < sys(isys)%c%ncel)
               if (gc%vac_axis == 0) then
                  gc%vac_keepmol = lmolcrys
                  vacperp = sys(isys)%c%vaclength / (sys(isys)%c%aa * sys(isys)%c%ar)
                  gc%vac_axis = 3
                  do im = 1, 2
                     if (vacperp(im) > vacperp(gc%vac_axis) + 1d-6) gc%vac_axis = im
                  end do
                  lvac = .true.
               end if

               call iw_separatortext("Vacuum",w%geometry_sec(geomsec_vacuum),hasmore=.true.)
               if (w%geometry_sec(geomsec_vacuum)%open) then
                  call iw_text("Axis",alignframe=.true.)
                  call iw_combo_simple("##cellvacaxis","a" // c_null_char // "b" // c_null_char // "c" // c_null_char,&
                     gc%vac_axis,changed=lch,sameline=.true.,startsatone=.true.)
                  call iw_tooltip("Axis of the cell along which the vacuum is set. The cell is cut by a&
                     & plane parallel to the other two axes and grows or shrinks there",ttshown)
                  ! another axis: its largest gap and its current vacuum
                  lwid = lch
                  if (lch) then
                     gc%vac_cutmode = vaccut_gap
                     lvac = .true.
                  end if
                  call iw_text("Slab position",alignframe=.true.,gap=.true.)
                  call iw_combo_simple("##cellvacpos","Center" // c_null_char // "Bottom" // c_null_char,&
                     gc%vac_pos,changed=lch,sameline=.true.,startsatone=.true.,tooltips=&
                     "The atoms are in the middle of the new cell, with half of the vacuum on each side" //&
                     c_null_char //&
                     "The lowest atom is at the bottom of the new cell, with all the vacuum above" // c_null_char,&
                     ttshown=ttshown)
                  lwid = lwid .or. lch
                  if (lch) lvacapply = .true.

                  ! where the cell is cut
                  call iw_text("Cut",alignframe=.true.)
                  call iw_combo_simple("##cellvaccutmode","Largest gap" // c_null_char // "Position" // c_null_char,&
                     gc%vac_cutmode,changed=lch,sameline=.true.,startsatone=.true.,tooltips=&
                     "Cut the cell in the middle of the largest empty gap along the axis: the vacuum of a&
                     & slab, or the widest spacing between atomic planes" // c_null_char //&
                     "Cut the cell at the given position along the axis: choose the position, the vacuum&
                     & and the options, then press Apply" // c_null_char,ttshown=ttshown)
                  lcut = lch
                  if (gc%vac_cutmode == vaccut_position) then
                     lch = iw_dragfloat_real8("##cellvaccut",x1=gc%vac_cut,speed=0.002d0,min=0d0,max=1d0,decimal=4,&
                        sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
                     lcut = lcut .or. lch
                  else
                     call iw_text("at " // string(gc%vac_cut,'f',decimal=4),sameline=.true.)
                  end if
                  call iw_tooltip("Fractional coordinate along the axis at which the cell is cut. The atoms&
                     & below it go to the top of the slab, and an atom at the cut stays at the bottom",ttshown)

                  ! the options (one that is on stays in view when More is
                  ! closed, and so do they in a molecular crystal, where
                  ! the molecules option starts on)
                  if (w%geometry_sec(geomsec_vacuum)%more .or. gc%vac_ortho .or. gc%vac_keepmol .or. lmolcrys) then
                     lch = iw_checkbox("Orthogonal axis##cellvacortho",gc%vac_ortho,labelleft=.true.)
                     lwid = lwid .or. lch
                     if (lch .and. gc%vac_ortho) lvacapply = .true.
                     call iw_tooltip("Make the axis perpendicular to the other two. Use it only with enough&
                        & vacuum: the periodic images of the slab along the axis slide over each other",ttshown)
                     lch = iw_checkbox("Keep molecules whole##cellvackeepmol",gc%vac_keepmol,labelleft=.true.,&
                        gap=.true.)
                     call iw_tooltip("Do not cut the molecules: the cut divides the centers of the&
                        & molecules (the spheres in the view) instead of the atoms. On at the start in a&
                        & molecular crystal",ttshown)
                     lcut = lcut .or. lch
                     if (lch) lvacapply = .true.
                  end if
                  lwid = lwid .or. lcut
               end if

               ! the current vacuum, found again when the axis, the cut,
               ! the geometry or the bonds change (not at every step of a
               ! running simulation: when it stops), with the section open
               ! or not (the outline in the view uses it). A position of
               ! the cut is not carried over to another lattice (a cell
               ! transformation from any window, undo). With the cut at
               ! the largest gap nothing is left pending: the vacuum of
               ! the field is the current one except while it is edited.
               ! Reset goes back to the height the cell has when the axis
               ! is chosen or the lattice is changed by anything but an
               ! edit of this section (which notes its own lattice): the
               ! height, not the vacuum, because the vacuum depends on how
               ! the atoms are cut (the molecules option, the gap).
               if (lvac .or. lcut .or. gc%vac_time < 0d0 .or.&
                  (sysc(isys)%timelastchange_rebond > gc%vac_time .and. .not.sysc(isys)%md_run)) then
                  if (any(abs(sys(isys)%c%m_x2c - gc%vac_m) > 1d-8)) then
                     gc%vac_cutmode = vaccut_gap
                     gc%vac_m = sys(isys)%c%m_x2c
                     lvac = .true.
                  end if
                  if (gc%vac_cutmode == vaccut_position) then
                     call sys(isys)%c%vacuum_info(gc%vac_axis,gc%vac_keepmol,gc%vac_cur,gc%vac_thick,&
                        cut=gc%vac_cut,ubot=gc%vac_ubot)
                  else
                     call sys(isys)%c%vacuum_info(gc%vac_axis,gc%vac_keepmol,gc%vac_cur,gc%vac_thick,&
                        cutpos=gc%vac_cut,ubot=gc%vac_ubot)
                  end if
                  if (lvac .or. gc%vac_time < 0d0 .or. gc%vac_cutmode == vaccut_gap) gc%vac_new = gc%vac_cur
                  if (lvac) gc%vac_hini = 1d0 / sys(isys)%c%ar(gc%vac_axis)
                  gc%vac_time = glfwGetTime()
               end if

               if (w%geometry_sec(geomsec_vacuum)%open) then
                  ! the vacuum: the outline of the new cell follows the value,
                  ! and the cell changes when the drag ends or the value is
                  ! entered (the structure is rebuilt once, and a cell that
                  ! was cut is not cut again at every step of the drag). A
                  ! value that the field shows as the current one is the
                  ! current one (Escape gives it back rounded). The vacuum
                  ! is not negative, unless it is already (molecules kept
                  ! whole that interlock): then it does not go below that.
                  call iw_text("Vacuum (Å)",alignframe=.true.)
                  lch = iw_dragfloat_real8("##cellvacnew",x1=gc%vac_new,speed=0.01d0,&
                     min=min(0d0,gc%vac_cur*bohrtoa),decimal=4,scale=bohrtoa,sameline=.true.,notlive=.true.,&
                     committed=lvaccommit,flags=ImGuiSliderFlags_AlwaysClamp)
                  lwid = lwid .or. lch
                  vactip = "Distance between the atoms on both sides of the cut, measured along the&
                     & perpendicular to the other two axes (negative if molecules kept whole interlock).&
                     & The new cell is outlined in the view, and the cell changes when the drag ends or&
                     & the value is entered"
                  if (gc%vac_cutmode == vaccut_position) vactip = vactip // " and Apply is pressed (now it is " //&
                     string(gc%vac_cur*bohrtoa,'f',decimal=4) // " Å at this cut)"
                  call iw_tooltip(vactip,ttshown)

                  ! back to the vacuum that gives the cell its height at the start
                  vacini = gc%vac_hini - gc%vac_thick
                  if (iw_button("Reset##cellvacreset",danger=.true.,sameline=.true.,&
                     disabled=(abs(gc%vac_new - vacini) * bohrtoa <= vacsame))) then
                     gc%vac_new = vacini
                     lwid = .true.
                     lvaccommit = .true.
                  end if
                  call iw_tooltip("Go back to the starting vacuum (" // string(vacini*bohrtoa,'f',decimal=4) //&
                     " Å): the cell gets the height it had along this axis when the system or the axis was&
                     & chosen here, or when the lattice was last changed by another tool or by undo",ttshown,&
                     whendisabled=.true.)

                  ! a new value is applied; one that the field shows as the
                  ! current one is the current one
                  if (lvaccommit) then
                     if (abs(gc%vac_new - gc%vac_cur) * bohrtoa > vacsame) then
                        lvacapply = .true.
                     elseif (gc%vac_cutmode == vaccut_gap) then
                        gc%vac_new = gc%vac_cur
                     end if
                  end if

                  ! a cut at a position is prepared and applied once (the
                  ! position refers to the atoms as they are); otherwise the
                  ! edits apply at once
                  if (gc%vac_cutmode == vaccut_position) then
                     if (iw_button("Apply##cellvacapply",danger=.true.)) iaction = iaction_set_vacuum
                     call iw_tooltip("Cut the cell at this position and set the vacuum, keeping the positions&
                        & of the atoms relative to each other",ttshown)
                  elseif (lvacapply) then
                     iaction = iaction_set_vacuum
                  end if

                  ! the error of the last attempt referred to the old inputs
                  if (lwid) w%errmsg = ""
               end if

               ! the plane of the cut and the outline of the new cell in the
               ! view, when the cut is chosen by hand or a new vacuum is
               ! being prepared, while the values of the section are those
               ! of this lattice. Both are drawn as the cell is (its origin,
               ! not the origin of the Display, along the other two axes).
               lch = (gc%vac_axis > 0 .and. gc%vac_time >= 0d0)
               if (lch) lch = associated(win(iview)%sc) .and. all(abs(sys(isys)%c%m_x2c - gc%vac_m) < 1d-8) .and.&
                  (gc%vac_cutmode == vaccut_position .or. abs(gc%vac_new - gc%vac_cur) > 1d-6)
               if (lch) then
                  ivj = mod(gc%vac_axis,3) + 1
                  ivk = mod(ivj,3) + 1

                  ! the plane, along the axis in the cell where the view
                  ! draws the atoms (the display origin or, with vacuum,
                  ! the reassembled slab); a cut at the edge of that cell,
                  ! the middle of the vacuum, stays where it is
                  vacwin = display_cell_window(isys,win(iview)%sc%disp,rep_display())
                  vacz = gc%vac_cut - vacwin(gc%vac_axis)
                  xb = 0d0
                  xb(gc%vac_axis) = gc%vac_cut
                  if (abs(vacz - nint(vacz)) > 1d-6) xb(gc%vac_axis) = gc%vac_cut - floor(vacz)
                  cutshp(1) = rep_shape(kind=shapekind_box,x1=sys(isys)%c%x2c(xb),&
                     rad=pv_edger*latvec_rad_def,rgb=pv_edge_rgb,alpha=pv_cut_alpha)
                  cutshp(1)%v(:,ivj) = sys(isys)%c%m_x2c(:,ivj)
                  cutshp(1)%v(:,ivk) = sys(isys)%c%m_x2c(:,ivk)

                  ! the new cell, around the atoms drawn in the view (not a
                  ! cell too flat to be made, see set_vacuum)
                  call vacuum_newcell(vacx0,vacm,vachnew)
                  cutshp(2) = rep_shape(kind=shapekind_box,x1=vacx0,v=vacm,rad=pv_edger*latvec_rad_def,&
                     rgb=pv_edge_rgb,alpha=0._c_float,shown=(vachnew >= 1d0))
                  call win(iview)%sc%show_transient_shapes(w%uid,tag_vaccut,cutshp)
               end if

               ! the centers of the molecules, while they are kept whole
               ! (not with the lattice points of the matrix mode, which
               ! they would hide where both coincide, nor during the pick
               ! of the atom at the new origin: they would hide atoms)
               if (gc%vac_keepmol .and. gc%vac_axis > 0 .and. associated(win(iview)%sc) .and.&
                  gc%mode /= cellmode_matrix .and. w%geometry_lvpick_row /= lvpick_origin) &
                  call vacuum_show_molcenters()
             end associate

             ! the origin of the cell: at a picked atom, or at a point
             ! calculated from the structure
             call iw_separatortext("Origin",w%geometry_sec(geomsec_origin))
             if (w%geometry_sec(geomsec_origin)%open) then
                if (iw_button("Pick Atom##celloriginpick",disabled=(iview == 0 .or.&
                   w%geometry_lvpick_row == lvpick_origin))) call lvpick_start(lvpick_origin,lvpick_origin)
                call iw_tooltip("Move the origin of the cell to an atom: click it in the view (a click&
                   & anywhere else cancels)",ttshown,whendisabled=.true.)
                if (iw_button("Center of Mass##cellorigincm",sameline=.true.)) then
                   iaction = iaction_shift_origin
                   iaction_i1 = origin_cm
                   iaction_i2 = 0
                end if
                call iw_tooltip("Move the origin of the cell to the center of mass of the atoms",ttshown)
                if (iw_button("Center Slab##cellorigincenter",sameline=.true.)) then
                   iaction = iaction_shift_origin
                   iaction_i1 = origin_center
                   iaction_i2 = 0
                end if
                call iw_tooltip("Move the origin so that the atoms are in the middle of the cell, along&
                   & the axes that have vacuum",ttshown)
                if (iw_button("Inversion##celloriginiinv",sameline=.true.)) then
                   iaction = iaction_shift_origin
                   iaction_i1 = origin_inversion
                   iaction_i2 = 0
                end if
                call iw_tooltip("Move the origin of the cell to the closest inversion center of the&
                   & space group",ttshown)
             end if
             ! the pick of the new origin lasts while this tab is shown
             if (w%geometry_lvpick_row == lvpick_origin) lvshown = .true.

             ! the preview of the new cell in the view: its lattice vectors
             ! (spheres in the axis colors at their ends, and arrows from
             ! the origin ending at their surface) and its outline, at its
             ! origin. In the matrix mode, also the lattice points around
             ! the displayed cells, for the picks (not during the pick of
             ! the atom at the new origin: they would hide the atoms).
             if (associated(win(iview)%sc)) then
                if (w%geometry_cell%mode == cellmode_matrix .and. w%geometry_lvpick_row /= lvpick_origin) then
                   call win(iview)%sc%show_transient_latvec(w%uid,1,replat_primitive,1)
                   lvshown = .true.
                end if
                if (pok .and. (w%geometry_cell%mode == cellmode_matrix .or. any(abs(pm - eye) > 1d-8) .or.&
                   any(abs(pt) > 1d-8))) then
                   do im = 1, 3
                      xb = sys(isys)%c%x2c(pm(:,im))
                      lvshp(im) = rep_shape(kind=shapekind_sphere,x1=xb,rad=pv_radfac*latvec_rad_def,&
                         rgb=ColorAxes_def(:,im))
                      lvshp(3+im) = rep_shape(kind=shapekind_arrow,shown=(norm2(xb) > 2d0*lvshp(im)%rad),&
                         rad=pv_arrowr*latvec_rad_def,headr=pv_headr,headl=pv_headl,rgb=ColorAxes_def(:,im))
                      if (lvshp(3+im)%shown) lvshp(3+im)%v(:,1) = xb * (1d0 - lvshp(im)%rad / norm2(xb))
                   end do
                   lvshp(7) = rep_shape(kind=shapekind_box,x1=sys(isys)%c%x2c(pt),&
                      v=matmul(sys(isys)%c%m_x2c,pm),rad=pv_edger*latvec_rad_def,rgb=pv_edge_rgb,alpha=0._c_float)
                   call win(iview)%sc%show_transient_shapes(w%uid,0,lvshp)
                end if
             end if

             ! the error of the last cell transformation or supercell search
             if (len_trim(w%errmsg) > 0) &
                call iw_text(w%errmsg,danger=.true.,wrap=.true.)

             call igEndTabItem()
          end if
       end if

       !! molecules tab !!
       if (iw_begintabitem("Molecules##drawgeometry_molstab")) then
          ! check if the tab changed
          call check_changed_tab("molecules")

          ! the selection/highlights in this tab refer to molecules
          table_hltype = atlisttype_nmol

          ! coordinate type for the center of mass
          call iw_text("Center of mass",alignframe=.true.)
          ldum = sysc(isys)%attype_combo_simple("##molcoordselectgeom",&
             w%geometry_moltype,atlisttype_coord_allowed,sameline=.true.)
          call iw_tooltip("Coordinate system for the molecular center of mass. The point group (Sym) and &
             &the Euler angles are shown for discrete molecules only",ttshown)

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
          sz0%y = sz0%y - iw_calcheight(3,0,.true.)
          if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
             icol = -1

             ! header setup
             call iw_table_column("Id",icol=icol,sortid=ic_id,icolsort=icolsort)

             call iw_table_column("Nat",icol=icol,sortid=ic_nat,icolsort=icolsort)

             call iw_table_column("Sym",icol=icol,sortid=ic_sym,icolsort=icolsort,flags=ImGuiTableColumnFlags_NoSort)

             if (doidx) then
                call iw_table_column("Idx",icol=icol,sortid=ic_idx,icolsort=icolsort)
             end if

             if (w%geometry_moltype == atlisttype_ncel_ang .or. w%geometry_moltype == atlisttype_ncel_bohr) then
                str2 = "/" // sysc(isys)%attype_coordinates_units(w%geometry_moltype)
                strx = "x" // str2
                stry = "y" // str2
                strz = "z" // str2
             else
                strx = "x"
                stry = "y"
                strz = "z"
             end if
             call iw_table_column(strx,icol=icol,sortid=ic_x,icolsort=icolsort)
             call iw_table_column(stry,icol=icol,sortid=ic_y,icolsort=icolsort)
             call iw_table_column(strz,icol=icol,sortid=ic_z,icolsort=icolsort)

             ! Euler angles (ZYZ, degrees) of the standard orientation
             call iw_table_column("α/°",icol=icol,sortid=ic_ea,icolsort=icolsort)
             call iw_table_column("β/°",icol=icol,sortid=ic_eb,icolsort=icolsort)
             call iw_table_column("γ/°",icol=icol,sortid=ic_eg,icolsort=icolsort)

             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! fetch the sort specs, sort the data if necessary
             call fetch_sort_specs()

             ! draw the header
             call iw_table_headers_row(autofit=.true.,writemenu=iwrite,ttshown=ttshown)

             ! sort
             if (forcesort) call table_sort()

             ! start the clipper; if the table is being written as text,
             ! its rows are captured as they are drawn (all of them)
             clipper = ImGuiListClipper_ImGuiListClipper()
             call ImGuiListClipper_Begin(clipper,ntype,-1._c_float)
             call w%table_write_begin(itable_molecules,iwrite,clipper)

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
             ! Step() ends the clipper on its own, but the object still has to go
             call ImGuiListClipper_destroy(clipper)
             call w%table_write_end(itable_molecules,table_title("Molecules",isys))
             call igEndTable()
          end if

          ! highlight/selection row
          call draw_highlight_buttons()

          ! edit row (remove/reorder molecules)
          call draw_mol_edit_buttons()

          call igEndTabItem()
       end if

       !! bonds tab !!
       flags = ImGuiTabItemFlags_None
       if (w%geometry_seltab == geomtab_bonds) flags = ImGuiTabItemFlags_SetSelected
       if (iw_begintabitem("Bonds##drawgeometry_bondstab",flags=flags)) then
          ! check if the tab changed
          call check_changed_tab("bonds")

          ! bond operations driven by clicking in the view window: the same
          ! three the builder offers, armed on the view showing this system,
          ! and shown as the same icon toggles (highlighted while armed)
          call iw_separatortext("Bond operations",w%geometry_sec(geomsec_bondops))
          if (w%geometry_sec(geomsec_bondops)%open) then
             call iw_text("In the view",alignframe=.true.)
             call bondmode_button("##geombondcreate",vm_builder_bond,&
                "Create bonds: bond two atoms clicked one after the other in the view"//&
                " window. Only the connectivity changes: no atom is moved or deleted")
             call bondmode_button("##geombondremove",vm_builder_bondremove,&
                "Remove bonds: remove the bond clicked in the view window")
             call bondmode_button("##geombondorder",vm_builder_bondorder,&
                "Bond order: cycle the order of the bond clicked in the view window:"//&
                " single, double, triple, aromatic, dashed")
             if (iview == 0) &
                call iw_text("(no view window open for this system)",sameline=.true.,&
                   alignframe=.true.)
          end if

          ! the bonds of the system
          call iw_separatortext("System bonds",w%geometry_sec(geomsec_sysbonds))

          ! number of cell atoms and digits for the Id column
          ntype = sys(isys)%c%ncel
          ndigit = ceiling(log10(ntype+0.1d0))
          table_hltype = atlisttype_ncel_frac

          ! the atomic species present, for the radii table of the bond
          ! calculation section (and its height, taken from the bonds table)
          natused_bonds = 0
          if (w%geometry_sec(geomsec_bondcalc)%open) then
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
          end if

          if (w%geometry_sec(geomsec_sysbonds)%open) then
             ! legend for the bond-type glyphs used in the Bonded atoms column
             call iw_text("─ single   ═ double   ≡ triple   ○ aromatic   ┄ dashed")
             call iw_helpermark("This table shows the bonds in the current system. The bonds drawn on the view window &
                &are usually the ones in this list. However, the user can choose to draw different bonding arrangements &
                &using the Bonds tab in the atoms object (Recalculate Bonds). In that case, adding or removing bonds in this table &
                &will not affect what is shown on the view.",sameline=.true.)

             ! bonds table (height leaves room for the bond calculation section below)
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
             sz0%y = sz0%y - iw_calcheight(1,0,.true.) - iw_separatortext_height()
             if (w%geometry_sec(geomsec_bondcalc)%open) &
                sz0%y = sz0%y - iw_calcheight(3,4,.false.) - iw_calcheight(min(5,natused_bonds)+1,0,.false.)
             sz0%y = max(sz0%y,iw_calcheight(3,0,.false.))
             if (igBeginTable(c_loc(str1),3,flags,sz0,0._c_float)) then
                call iw_table_column("Id",id=0)
                call iw_table_column("Atom",id=1)
                call iw_table_column("Bonded atoms",id=2)
                call iw_table_headers_row(freezetop=.true.,autofit=.true.,writemenu=iwrite,ttshown=ttshown)

                ! draw the rows (clipped for performance; all of them if the
                ! table is being written as text)
                clipper = ImGuiListClipper_ImGuiListClipper()
                call ImGuiListClipper_Begin(clipper,ntype,-1._c_float)
                call w%table_write_begin(itable_bonds,iwrite,clipper)
                do while(ImGuiListClipper_Step(clipper))
                   call c_f_pointer(clipper,clipper_f)
                   do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                      suffix = "_" // string(i)
                      call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                      icol = -1

                      ! id (with a row-spanning selectable to detect hover)
                      icol = icol + 1
                      if (igTableSetColumnIndex(icol)) then
                         call iw_text(string(i,ndigit),alignframe=.true.)
                      end if
                      ! emit even when column 0 is clipped off-screen
                      if (iw_highlight_selectable("##bondselect" // suffix)) ihlbond = i

                      ! atom name
                      icol = icol + 1
                      if (igTableSetColumnIndex(icol)) then
                         call iw_text(trim(sys(isys)%c%at(sys(isys)%c%atcel(i)%idx)%name),alignframe=.true.)
                      end if

                      ! bonded atoms (one colored button per neighbor)
                      icol = icol + 1
                      if (igTableSetColumnIndex(icol)) then
                         ! "+" button: pick an atom in the view to add a bond to this atom
                         if (iw_button("+##addbond" // suffix,disabled=(iview == 0),nocapture=.true.)) then
                            call w%geometry_addbond%stage((/i,0,0,0/))
                            w%geometry_lvpick_row = 0
                            call win(iview)%viewmode_set_forced(vm_pick_atom,&
                               "Pick an atom to bond to atom " // string(i),w%id)
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

                            ! the button (colored by the bonded atom's color in the
                            ! view, prefixed by the bond-type glyph). Left-click
                            ! opens a menu, right-click removes the bond.
                            havergb_ = color_from_view(atlisttype_ncel_frac,jj,rgb)
                            ldum = iw_atom_button(bondglyph // " " // string(jj) //&
                               "##bond" // suffix // "_" // string(j),rgb,havergb=havergb_,&
                               sameline=.true.)

                            ! hovering highlights the row's atom and marks this
                            ! neighbor, and shows a tooltip with the bond details;
                            ! right-click removes the bond (deferred)
                            if (igIsItemHovered(ImGuiHoveredFlags_None)) then
                               ihlbond = i
                               ihlbtn = jj
                               ! highlight the row in the table
                               color = igGetColorU32_Vec4(ImVec4(ColorTableHighlightRow(1),&
                                  ColorTableHighlightRow(2),ColorTableHighlightRow(3),&
                                  ColorTableHighlightRow(4)))
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
                               if (iw_beginmenu("Order")) then
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
                ! Step() ends the clipper on its own, but the object still has to go
                call ImGuiListClipper_destroy(clipper)
                call w%table_write_end(itable_bonds,table_title("Bonds",isys))
                call igEndTable()
             end if

             ! perform deferred bond edits (after the table loop, so nstar is not
             ! edited while being iterated)
             if (ibrm1 > 0) call sysc(isys)%remove_bond(ibrm1,ibrm2,lbrm)
             if (ibord1 > 0) call sysc(isys)%set_bond_order(ibord1,ibord2,lbord,ibordval)
          end if

          ! recalculate bonds section
          call iw_separatortext("Bond calculation",w%geometry_sec(geomsec_bondcalc))
          if (w%geometry_sec(geomsec_bondcalc)%open) then
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
                call iw_table_column("Atom",id=0)
                call iw_table_column("Z",id=1)
                call iw_table_column("Radius (Å)",id=2)
                do j = 1, natused_bonds
                   call iw_table_column(trim(name_bonds(j)),id=2+j)
                end do
                call iw_table_headers_row(freezetop=.true.)

                do i = 1, natused_bonds
                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   iz = iat_bonds(i)
                   if (igTableSetColumnIndex(0)) then
                      call iw_text(trim(name_bonds(i)),alignframe=.true.)
                   end if
                   if (igTableSetColumnIndex(1)) then
                      call iw_text(string(iz),alignframe=.true.)
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

             ! bond factor and bond delta drags
             call iw_text("Bond factor (f)",alignframe=.true.)
             ldum = iw_dragfloat_real8("##bondfactor_geom",x1=sysc(isys)%bondfactor,speed=0.001d0,&
                min=1d0,max=4d0,decimal=4,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Bond factor parameter (multiplicative) for non-metal bonding (see formula above)",ttshown)
             call iw_text("Bond delta δ (Å)",alignframe=.true.,gap=.true.)
             ldum = iw_dragfloat_real8("##bonddelta_geom",x1=sysc(isys)%bonddelta,speed=0.001d0,&
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
             deallocate(iat_bonds,name_bonds)
          end if

          call igEndTabItem()
       end if

       !! symmetry tab !!
       flags = ImGuiTabItemFlags_None
       if (w%geometry_seltab == geomtab_symmetry) flags = ImGuiTabItemFlags_SetSelected
       if (iw_begintabitem("Symmetry##drawgeometry_symmetrytab",flags=flags)) then
          ! check if the tab changed
          call check_changed_tab("symmetry")

          if (sys(isys)%c%ismolecule) then
             !! symmetry: molecule !!

             ! current point group
             call iw_separatortext("Point group",w%geometry_sec(geomsec_symgroup))
             if (w%geometry_sec(geomsec_symgroup)%open) then
                call iw_text("Schoenflies",alignframe=.true.)
                if (sys(isys)%c%pg%avail) then
                   call iw_text(trim(sys(isys)%c%pg%symbol),sameline=.true.)
                else
                   call iw_text("not available (symmetry not computed)",sameline=.true.)
                end if

                ! symmetrize the molecular geometry
                if (iw_button("Refine##symrefinemol")) iaction = iaction_sym_refine
                call iw_tooltip("Move atoms to occupy their symmetry positions exactly",ttshown)
             end if

             ! symmetry operations table
             if (sys(isys)%c%pg%avail) then
                call iw_separatortext("Operations (" // string(sys(isys)%c%pg%nop) // ")",w%geometry_sec(geomsec_symops))
                call symop_ensure_sel(sys(isys)%c%pg%nop)
                if (w%geometry_sec(geomsec_symops)%open) then
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
                   szavail%y = szavail%y - iw_calcheight(3,0,.true.)
                   if (sz0%y > szavail%y) sz0%y = szavail%y
                   if (igBeginTable(c_loc(str1),3,flags,sz0,0._c_float)) then
                      call iw_table_column("#",id=0)
                      call iw_table_column("Sym",id=1)
                      call iw_table_column("Axis (Å)",id=2)
                      call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)
                      call w%table_write_begin(itable_symopsmol,iwrite)

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
                      call w%table_write_end(itable_symopsmol,table_title("Symmetry operations",isys))
                      call igEndTable()
                   end if

                   ! the selection row
                   call symop_buttons(sys(isys)%c%pg%nop)
                end if

                ! draw the selected/hovered symmetry elements in the view
                call symop_display(sys(isys)%c%pg%nop)
             end if
          else
             !! symmetry: crystal !!
             ! current space group
             call iw_separatortext("Space group",w%geometry_sec(geomsec_symgroup))
             if (w%geometry_sec(geomsec_symgroup)%open) then
                if (sys(isys)%c%spgavail) then
                   call pointgroup_info(sys(isys)%c%spg%pointgroup_symbol,schpg,holo,laue)
                   call iw_text("Hermann-Mauguin",alignframe=.true.)
                   call iw_text(trim(sys(isys)%c%spg%international_symbol) //&
                      " (number " // string(sys(isys)%c%spg%spacegroup_number) // ")",sameline=.true.)
                   call iw_text("Point group",alignframe=.true.)
                   call iw_text(trim(sys(isys)%c%spg%pointgroup_symbol) // " (" // trim(schpg) // ")",&
                      sameline=.true.)
                   call iw_text("Holohedry",alignframe=.true.)
                   call iw_text(trim(holo_string(holo)),sameline=.true.)
                   call iw_text("Laue class",alignframe=.true.,gap=.true.)
                   call iw_text(trim(laue_string(laue)),sameline=.true.)
                else
                   call iw_text("Not available (symmetry not computed)")
                end if

                call iw_text("Tolerance (Å)",alignframe=.true.)
                str_eps = string(sysc(isys)%symeps*bohrtoa,'e',decimal=2)
                if (iw_inputtext("##symepsinput",bufsize=32,texta=str_eps,width=8,sameline=.true.,notlive=.true.)) then
                   if (isreal(res,str_eps)) then
                      if (res > 0d0) sysc(isys)%symeps = res / bohrtoa
                   end if
                end if
                call iw_tooltip("Distance tolerance (Å) used to detect the symmetry (symeps)",ttshown)
                if (iw_button("Recalculate##symrecalc",sameline=.true.)) iaction = iaction_sym_recalc
                call iw_tooltip("Recalculate crystal symmetry at the given tolerance (symeps) level",ttshown)

                if (iw_button("Clear##symclear",danger=.true.)) iaction = iaction_sym_clear
                call iw_tooltip("Remove crystal symmetry",ttshown)
                if (iw_button("Refine##symrefine",sameline=.true.)) iaction = iaction_sym_refine
                call iw_tooltip("Idealize the cell parameters and move the atoms to their ideal &
                   &symmetry positions.",ttshown)
                if (iw_button("Whole Molecules##symwholemols",sameline=.true.,&
                   disabled=.not.sys(isys)%c%ismol3d)) iaction = iaction_sym_wholemols
                call iw_tooltip("Choose a symmetry subgroup such that the asymmetric unit contains whole molecules &
                   &(molecular crystals only)",ttshown)
             end if

             ! symmetry operations table
             neqv = sys(isys)%c%neqv
             ncv = sys(isys)%c%ncv
             call iw_separatortext("Operations (" // string(neqv*max(ncv,1)) // ")",w%geometry_sec(geomsec_symops))

             if (neqv >= 1) then
                call symop_ensure_sel(neqv*max(ncv,1))
                if (w%geometry_sec(geomsec_symops)%open) then
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
                   sz0%y = iw_calcheight(min(neqv*max(ncv,1),8)+1,0,.false.)
                   if (igBeginTable(c_loc(str1),5,flags,sz0,0._c_float)) then
                      call iw_table_column("#",id=0)
                      call iw_table_column("HM",id=1)
                      call iw_table_column("Operation",id=2)
                      call iw_table_column("Axis (cryst.)",id=3)
                      call iw_table_column("Axis (Cartesian)",id=4)
                      call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)
                      call w%table_write_begin(itable_symops,iwrite)

                      ! every operation, centered copies included: they are
                      ! distinct isometries with their own symbol and elements
                      do i = 1, neqv*max(ncv,1)
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
                      call w%table_write_end(itable_symops,table_title("Symmetry operations",isys))
                      call igEndTable()
                   end if

                   ! the selection row
                   call symop_buttons(neqv*max(ncv,1))
                end if

                ! draw the selected/hovered symmetry elements in the view
                call symop_display(neqv*max(ncv,1))
             end if

             ! centering vectors table
             call iw_separatortext("Centering vectors (" // string(ncv) // ")",w%geometry_sec(geomsec_symcen))
             if (w%geometry_sec(geomsec_symcen)%open) then
                flags = ImGuiTableFlags_None
                flags = ior(flags,ImGuiTableFlags_RowBg)
                flags = ior(flags,ImGuiTableFlags_Borders)
                flags = ior(flags,ImGuiTableFlags_ScrollY)
                flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
                str1 = "##symcentable" // c_null_char
                sz0%x = 0
                sz0%y = iw_calcheight(min(ncv,4)+1,0,.false.)
                if (igBeginTable(c_loc(str1),2,flags,sz0,0._c_float)) then
                   call iw_table_column("#",id=0)
                   call iw_table_column("Coordinates (fractional)",id=1)
                   call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)
                   call w%table_write_begin(itable_symcen,iwrite)
                   do i = 1, ncv
                      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                      if (igTableSetColumnIndex(0)) call iw_text(string(i))
                      if (igTableSetColumnIndex(1)) call iw_text(cell_cen_label(sys(isys)%c%cen(:,i)))
                   end do
                   call w%table_write_end(itable_symcen,table_title("Centering vectors",isys))
                   call igEndTable()
                end if
             end if

             ! space group vs tolerance analysis
             call iw_separatortext("Space group analysis",w%geometry_sec(geomsec_symanal))
             if (w%geometry_sec(geomsec_symanal)%open) then
                call iw_text("Calculate the space group as a function of the symmetry &
                   &tolerance (symprec). Click a row to adopt that tolerance and recalculate.",&
                   disabled=.true.,wrap=.true.)
                if (iw_button("Analyze##symanalyze")) iaction = iaction_sym_analyze
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
                      call iw_table_column("Symprec (Å)",id=0)
                      call iw_table_column("Space group",id=1)
                      call iw_table_column("Number",id=2)
                      call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)
                      call w%table_write_begin(itable_symanal,iwrite)
                      do i = 1, size(w%geometry_sym_analyze_eps,1)
                         call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                         if (igTableSetColumnIndex(0)) then
                            str2 = string(w%geometry_sym_analyze_eps(i)*bohrtoa,'e',10,2) //&
                               "##symanalrow" // string(i) // c_null_char
                            call iw_capture_text(str2,label=.true.)
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
                      call w%table_write_end(itable_symanal,table_title("Space group analysis",isys))
                      call igEndTable()
                   end if
                end if
             end if
          end if ! .not. ismolecule

          ! the error of the last symmetry operation
          if (len_trim(w%errmsg) > 0) &
             call iw_text(w%errmsg,danger=.true.,wrap=.true.)

          call igEndTabItem()
       end if

       call igEndTabBar()
    end if
    w%geometry_seltab = geomtab_none ! the request has been served

    call igEndGroup()

    ! the lattice-vector pick ends with the lattice points it picks
    ! from, and the pick of the origin with the cell tab
    if (.not.lvshown .and. w%geometry_lvpick_row > 0) call lvpick_cancel()

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

          ! standard-orientation axes: one arrow per axis (columns of stdrot),
          ! the Axes object proportions scaled to the axis length
          do jm = 1, 3
             stdshp(jm) = rep_shape(kind=shapekind_arrow,x1=stdcom,v=0d0,&
                rad=axes_radius_def * stdaxlen / axes_length_def,&
                headr=axes_coneradius_def / axes_radius_def,&
                headl=min(axes_conelength_def,axes_length_def) / axes_length_def,&
                rgb=ColorAxes_def(:,jm))
             stdshp(jm)%v(:,1) = stdaxlen * stdrot(:,jm)
          end do
          call sysc(isys)%sc%show_transient_shapes(w%uid,ihighlight,stdshp)

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

             ! a black cylinder spanning stdcom +/- rotlen*rotdir (the
             ! tags are negative, apart from the axes above)
             rotshp(1) = rep_shape(kind=shapekind_cylinder,x1=stdcom - rotlen * rotdir,v=0d0,&
                rad=rotaxis_radius_def,rgb=ColorRotaxis_def)
             rotshp(1)%v(:,1) = 2d0 * rotlen * rotdir
             call sysc(isys)%sc%show_transient_shapes(w%uid,-(ihighlight*4+ieuler_drag),rotshp)
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
       if (empty_row(iclicked)) then
          w%geometry_spcsel(iclicked) = (istate /= 2)
       elseif (istate == 2) then
          call sysc(isys)%highlight_clear(.false.,(/iclicked/),table_hltype)
       else
          call sysc(isys)%highlight_atoms(.false.,(/iclicked/),table_hltype,&
             reshape(w%geometry_select_rgba,(/4,1/)))
       end if
    elseif (iclicked == -1) then
       ! range (shift) or single (ctrl) selection: select the rows. A row with
       ! no atoms has nothing to highlight, so it is flagged on the side
       nhigh = iclicked_end - iclicked_ini + 1
       allocate(ihigh(nhigh),irgba(4,nhigh))
       nhigh = 0
       do ii = iclicked_ini, iclicked_end
          if (empty_row(w%iord(ii))) then
             w%geometry_spcsel(w%iord(ii)) = .true.
             cycle
          end if
          nhigh = nhigh + 1
          ihigh(nhigh) = w%iord(ii)
          irgba(:,nhigh) = w%geometry_select_rgba
       end do
       if (nhigh > 0) &
          call sysc(isys)%highlight_atoms(.false.,ihigh(1:nhigh),table_hltype,irgba(:,1:nhigh))
       deallocate(ihigh,irgba)
    end if

    ! remove/merge highlighted atoms
    if (w%focused() .and. is_bind_event(BIND_EDITSELECT_REMOVE)) then
       iaction = iaction_edit_highlighted
       iaction_i1 = edit_remove
    end if

    ! deselect all highlighted atoms
    deselected = .false.
    if (w%focused() .and. is_bind_event(BIND_CANCEL)) then
       if (allocated(sysc(isys)%highlight_rgba)) then
          if (any(sysc(isys)%highlight_rgba >= 0._c_float)) then
             call sysc(isys)%highlight_clear(.false.)
             deselected = .true.
          end if
       end if
       if (any_empty_row_selected()) then
          call clear_empty_rows()
          deselected = .true.
       end if
    end if

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(5,1)

    ! close button
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if (.not.deselected .and. iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.
    doquit = doquit .or. iw_button("Close")

    ! quit the window
    if (doquit) &
       call w%end()

    ! process actions at the end; clear any previous error when a new action runs
    if (iaction /= -1) w%errmsg = ""
    if (iaction == iaction_restore) then
       call sysc(isys)%reread_geometry_from_file(errmsg=w%errmsg)

    elseif (iaction == iaction_set_attype_name) then
       call sysc(isys)%set_attype_name(w%geometry_atomtype,iaction_i1,iaction_str)

    elseif (iaction == iaction_set_atomic_number) then
       call sysc(isys)%set_atomic_number(w%geometry_atomtype,iaction_i1,iaction_i2,setatomnames=.true.,&
          copybonding=w%geometry_keepbonding)

    elseif (iaction == iaction_add_species) then
       call sysc(isys)%add_species(iaction_i1)

    elseif (iaction == iaction_set_attype_species) then
       call sysc(isys)%set_attype_species(w%geometry_atomtype,iaction_i1,iaction_i2,&
          copybonding=w%geometry_keepbonding,errmsg=w%errmsg)

    elseif (iaction == iaction_set_atom_position) then
       call sysc(isys)%set_atom_position(w%geometry_atomtype,iaction_i1,iaction_x,iaction_l,&
          copybonding=w%geometry_keepbonding,errmsg=w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_set_atom_occupancy) then
       call sysc(isys)%set_attype_occupancy(w%geometry_atomtype,iaction_i1,iaction_x(1))
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_add_species_change_atom) then
       call sysc(isys)%add_species(iaction_i1)
       call sysc(isys)%set_attype_species(w%geometry_atomtype,iaction_i2,sys(isys)%c%nspc,&
          copybonding=w%geometry_keepbonding,errmsg=w%errmsg)

    elseif (iaction == iaction_add_atom) then
       call sysc(isys)%attype_add_atom(w%geometry_atomtype,iaction_i1,iaction_x,errmsg=w%errmsg)
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
       call sysc(isys)%attype_reorder(w%geometry_atomtype,w%iord,w%errmsg)

    elseif (iaction == iaction_swap_atom_ids) then
       call sysc(isys)%attype_swap_atoms(w%geometry_atomtype,iaction_i1,iaction_i2,w%errmsg)

    elseif (iaction == iaction_swap_mol_ids) then
       call sysc(isys)%swap_molecules(iaction_i1,iaction_i2,w%errmsg)

    elseif (iaction == iaction_set_molecule_position) then
       call sysc(isys)%set_molecule_position(w%geometry_moltype,iaction_i1,iaction_x,&
          copybonding=w%geometry_keepbonding,errmsg=w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_set_molecule_rotation) then
       call sysc(isys)%set_molecule_rotation(iaction_i1,iaction_x,copybonding=w%geometry_keepbonding,errmsg=w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_remove_molecules) then
       call sysc(isys)%edit_highlighted_atoms(remove=.true.,errmsg=w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_reorder_molecules) then
       call sysc(isys)%attype_reorder(atlisttype_nmol,w%iord,w%errmsg)

    elseif (iaction == iaction_change_cell) then
       call sysc(isys)%move_cell(iaction_x6(1:3),iaction_x6(4:6),iaction_l,&
          copybonding=w%geometry_keepbonding,errmsg=w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_transform_cell) then
       call sysc(isys)%transform_cell(iaction_i1,iaction_l,errmsg=w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_transform_matrix) then
       call sysc(isys)%transform_cell_matrix(iaction_m,iaction_x,iaction_l,errmsg=w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.
       ! the inputs described a transformation of the cell that was replaced
       if (len_trim(w%errmsg) == 0) call reset_cell_inputs()

    elseif (iaction == iaction_shift_origin) then
       call sysc(isys)%shift_origin(iaction_i1,iaction_i2,w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.
       ! the position of the cut of the vacuum referred to the old origin
       if (len_trim(w%errmsg) == 0) call reset_vacuum_inputs()

    elseif (iaction == iaction_set_vacuum) then
       associate (gc => w%geometry_cell)
         ! where the new cell is now, for the atoms the view draws
         call vacuum_newcell(vacx0,vacm,vachnew)
         vacz = 0d0
         if (gc%vac_pos == vacpos_center .and. vachnew > 0d0) vacz = 0.5d0 * gc%vac_new / vachnew
         im = gc%vac_axis
         call sysc(isys)%set_vacuum(gc%vac_axis,gc%vac_new,gc%vac_pos,gc%vac_ortho,gc%vac_keepmol,&
            gc%vac_cutmode == vaccut_position,gc%vac_cut,w%geometry_keepbonding,w%errmsg,lch)
         if (len_trim(w%errmsg) > 0) then
            ! a vacuum that could not be set goes back to the current one
            if (gc%vac_cutmode == vaccut_gap) gc%vac_new = gc%vac_cur
         else
            ! the atoms stay where they are on the screen, and the cell
            ! lands on its outline: the lowest atom of the slab, at vacz
            ! in the new cell, is drawn in the cell that the view gives
            ! it now (in the main view of the system in real space: the
            ! translation is for the main scene)
            if (lch .and. win(iview)%ismain .and. .not.win(iview)%kspace .and. associated(win(iview)%sc)) then
               vacwin = display_cell_window(isys,win(iview)%sc%disp,rep_display())
               call sysc(isys)%view_shift(vacx0 + floor(vacz - vacwin(im) + 1d-6) * sys(isys)%c%m_x2c(:,im))
            end if
            ! the position of the cut referred to the atoms before the
            ! edit; the new lattice is this section's (Reset keeps its
            ! starting point)
            call reset_vacuum_inputs()
            gc%vac_m = sys(isys)%c%m_x2c
         end if
       end associate
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_sym_recalc) then
       call sysc(isys)%recalc_symmetry(w%errmsg)
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_sym_clear) then
       call sysc(isys)%clear_symmetry()
       sysc(isys)%sc%nextbuildlists_fixcam = .true.

    elseif (iaction == iaction_sym_refine) then
       call sysc(isys)%refine_symmetry(w%errmsg)
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
    ! One of the bond operations driven by clicking in the view window:
    ! an icon toggle showing the mode's own glyph, highlighted while the
    ! mode is armed on the view showing this system.
    subroutine bondmode_button(strid,mode,ttip)
      character(len=*), intent(in) :: strid, ttip
      integer, intent(in) :: mode

      logical :: lstate
      integer(c_int) :: itex
      character(len=:), allocatable :: fallback

      call viewmode_icon(mode,itex,fallback)
      lstate = (w%builder_vm == mode)
      if (iw_icon_togglebutton(strid,itex,fallback,state=lstate,disabled=(iview == 0),&
         sameline=.true.)) then
         w%errmsg = ""
         w%geometry_lvpick_row = 0
         call bondmode_toggle(w,iview,isys,mode)
      end if
      call iw_tooltip(ttip,ttshown,whendisabled=.true.)

    end subroutine bondmode_button

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
    !> Drop the state cached for the previous system, after the anchor view has
    !> moved to a different one.
    subroutine reset_for_new_system()

      ! reset the last-selected row, the empty-species selection and the sort
      w%lastselected = 0
      call clear_empty_rows()
      call w%table_write_cancel()

      ! remove the cached cell-transformation data and reorder the table
      call clear_search_results()
      call clear_sym_cache()
      call reset_sort()

      ! clear the cell transformations
      call reset_cell_inputs()
      call reset_vacuum_inputs()
      w%geometry_cell%vac_axis = 0
      w%geometry_lvpick_row = 0

    end subroutine reset_for_new_system

    ! Start a pick of the cell tab in the anchor view: the lattice
    ! vector in row irow of the matrix cell transformation (the pick
    ! goes on to the next rows up to ilast), or the atom at the new
    ! origin of the cell (irow = lvpick_origin). It replaces a pending
    ! add-bond pick and the other pick of the cell tab.
    subroutine lvpick_start(irow,ilast)
      integer, intent(in) :: irow, ilast

      character(len=3), parameter :: lbl(3) = (/"a' ","b' ","c' "/)

      call w%geometry_addbond%clear()
      w%geometry_lvpick_row = irow
      w%geometry_lvpick_last = ilast
      call w%geometry_lvpick%arm()
      if (irow == lvpick_origin) then
         if (w%errmsg == errorigin) w%errmsg = ""
         call win(iview)%viewmode_set_forced(vm_pick_atom,"Click the atom at the new origin",w%id)
      else
         call win(iview)%viewmode_set_forced(vm_pick_atom,"Click the lattice point of " // trim(lbl(irow)),&
            w%id,latonly=.true.)
      end if

    end subroutine lvpick_start

    ! Cancel a pending lattice-vector or origin pick and release the view.
    subroutine lvpick_cancel()

      call win(iview)%viewmode_release_forced(w%id)
      w%geometry_lvpick_row = 0

    end subroutine lvpick_cancel

    ! The inputs of the supercell transformation, in the simple,
    ! matrix and orient modes, back to their defaults.
    subroutine reset_cell_inputs()

      w%geometry_cell%nrep = 1_c_int
      w%geometry_cell%intmat = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
      w%geometry_cell%cen = 1
      w%geometry_cell%origin = 0d0
      w%geometry_cell%inverse = .false.
      call reset_orient_inputs()

    end subroutine reset_cell_inputs

    ! The cell that the vacuum section gives with its current values,
    ! in the frame of the cell now: x0 is its origin (Cartesian, bohr),
    ! m are its lattice vectors (columns) and hnew is its height over
    ! the other two axes. The rules are those of crystal%set_vacuum:
    ! the height is the thickness of the slab plus the vacuum, the axis
    ! keeps its direction or is perpendicular to the other two, and the
    ! lowest atom of the slab is at the bottom face or half the vacuum
    ! over it. That atom is taken in the cell where the view draws the
    ! atoms (the display origin or, along an axis with vacuum, the
    ! reassembled slab).
    subroutine vacuum_newcell(x0,m,hnew)
      real*8, intent(out) :: x0(3), m(3,3), hnew

      real*8 :: z0, xwin(3), d(3)

      associate (gc => w%geometry_cell)
        ! direction of the new axis, of unit height over the other two
        if (gc%vac_ortho) then
           d = sys(isys)%c%m_c2x(gc%vac_axis,:) / sys(isys)%c%ar(gc%vac_axis)
        else
           d = sys(isys)%c%m_x2c(:,gc%vac_axis) * sys(isys)%c%ar(gc%vac_axis)
        end if
        hnew = gc%vac_thick + gc%vac_new
        m = sys(isys)%c%m_x2c
        m(:,gc%vac_axis) = hnew * d
        z0 = 0d0
        if (gc%vac_pos == vacpos_center) z0 = 0.5d0 * gc%vac_new
        xwin = 0d0
        if (associated(win(iview)%sc)) xwin = display_cell_window(isys,win(iview)%sc%disp,rep_display())
        x0 = (gc%vac_ubot - floor(gc%vac_ubot - xwin(gc%vac_axis) + 1d-6)) *&
           sys(isys)%c%m_x2c(:,gc%vac_axis) - z0 * d
      end associate

    end subroutine vacuum_newcell

    ! Show the centers of the units that the vacuum section keeps whole
    ! as spheres in the view. They are the points that the cut divides
    ! in crystal%set_vacuum (see vacuum_layout): the average position
    ! of the atoms of each molecule, layer or chain that is not
    ! periodic along the axis; a unit of one atom has no sphere. A
    ! molecule that the view draws whole (the connected molecules of
    ! the Display) has the sphere at the center of the atoms as drawn,
    ! with the rules of the atoms in add_draw_elements: the origin of
    ! the Display moves each atom, and along an axis with vacuum the
    ! molecule is in the reassembled slab. The other units have it in
    ! the cell where the view draws the atoms.
    subroutine vacuum_show_molcenters()
      integer :: imol, k, i, n
      real*8 :: xwin(3), tshift(3), xc(3), xk(3)
      logical :: whole, asdrawn, usetshift, isvac(3)
      type(rep_shape), allocatable :: shp(:)

      if (.not.allocated(sys(isys)%c%mol) .or. .not.allocated(sys(isys)%c%lvecmolc)) return
      associate (c => sys(isys)%c, iaxis => w%geometry_cell%vac_axis, disp => win(iview)%sc%disp)
        allocate(shp(c%nmol))
        ! the cell of the view: along an axis with vacuum it is not at
        ! the origin of the Display
        xwin = display_cell_window(isys,disp,rep_display())
        tshift = disp%tshift
        usetshift = any(abs(tshift) > 1d-5)
        isvac = (abs(xwin - tshift) > 1d-6)
        n = 0
        do imol = 1, c%nmol
           if (c%mol(imol)%nat < 2) cycle
           whole = c%mol(imol)%discrete
           if (.not.whole) then
              whole = .true.
              if (allocated(c%mol(imol)%lvec)) whole = all(c%mol(imol)%lvec(iaxis,1:c%mol(imol)%nlvec) == 0)
           end if
           if (.not.whole) cycle

           asdrawn = (c%mol(imol)%discrete .and. disp%onemotif)
           xc = 0d0
           do k = 1, c%mol(imol)%nat
              i = c%mol(imol)%at(k)%cidx
              if (asdrawn) then
                 xk = c%mol(imol)%at(k)%x
                 if (usetshift) then
                    where (.not.isvac) xk = xk - floor(c%atcel(i)%x - tshift)
                 end if
              else
                 xk = c%atcel(i)%x + c%lvecmolc(:,i)
              end if
              xc = xc + xk
           end do
           xc = xc / real(c%mol(imol)%nat,8)
           where (isvac .or. .not.asdrawn) xc = xc - floor(xc - xwin + 1d-6)
           n = n + 1
           shp(n) = rep_shape(kind=shapekind_sphere,x1=c%x2c(xc),rad=pv_molfac*latvec_rad_def,&
              rgb=pv_mol_rgb,rim=.true.)
        end do
        if (n > 0) call win(iview)%sc%show_transient_shapes(w%uid,tag_vacmol,shp(1:n))
      end associate

    end subroutine vacuum_show_molcenters

    ! The vacuum and the cut of the vacuum section refer to a structure
    ! that was replaced: find them again, with the cut at the largest gap.
    subroutine reset_vacuum_inputs()

      w%geometry_cell%vac_cutmode = vaccut_gap
      w%geometry_cell%vac_time = -1d0

    end subroutine reset_vacuum_inputs

    ! The inputs of the orient mode of the supercell transformation back
    ! to their defaults.
    subroutine reset_orient_inputs()

      w%geometry_cell%orient_dir = .false.
      w%geometry_cell%orient_ivec = (/0_c_int,0_c_int,1_c_int/)
      w%geometry_cell%orient_axis = 3
      w%geometry_cell%orient_prim = .false.

    end subroutine reset_orient_inputs

    ! The supercell transformation that the supercell section of the
    ! cell tab describes in its current mode: the lattice vectors of
    ! the new cell (columns of m, crystallographic coordinates of the
    ! current cell) and its origin (t), as newcell will use them (a
    ! left-handed set of vectors is inverted). isok is false if there
    ! is none (no supercell selected in the search) or the vectors do
    ! not make a cell, and msg says why.
    subroutine pending_cell(m,t,isok,msg)
      real*8, intent(out) :: m(3,3), t(3)
      logical, intent(out) :: isok
      character(len=:), allocatable, intent(out) :: msg

      integer :: k

      m = eye
      t = 0d0
      isok = .true.
      msg = ""
      if (w%geometry_cell%mode == cellmode_simple) then
         do k = 1, 3
            m(k,k) = real(w%geometry_cell%nrep(k),8)
         end do
      elseif (w%geometry_cell%mode == cellmode_matrix) then
         ! integer part plus centering
         do k = 1, 3
            m(:,k) = real(w%geometry_cell%intmat(:,k),8) +&
               sys(isys)%c%cen(:,min(max(w%geometry_cell%cen(k),1),max(sys(isys)%c%ncv,1)))
         end do
         t = w%geometry_cell%origin
         if (w%geometry_cell%inverse) then
            isok = (abs(det3(m)) > 1d-8)
            if (isok) call matinv(m,3)
         end if
      elseif (w%geometry_cell%mode == cellmode_orient) then
         call sys(isys)%c%cell_orient(int(w%geometry_cell%orient_ivec),w%geometry_cell%orient_dir,&
            w%geometry_cell%orient_axis,w%geometry_cell%orient_prim,m,msg)
         isok = (len_trim(msg) == 0)
      else
         isok = allocated(w%geometry_cell%cand) .and. w%geometry_cell%isel > 0
         if (isok) isok = (w%geometry_cell%isel <= size(w%geometry_cell%cand,1))
         if (isok) m = real(w%geometry_cell%cand(w%geometry_cell%isel)%m,8)
         if (.not.isok) msg = "Select a supercell in the table"
      end if
      ! the same threshold newcell uses, which also flips a left-handed cell
      if (isok) isok = (abs(det3(m)) * max(sys(isys)%c%ncv,1) > 1d-2)
      if (isok) then
         if (det3(m) < 0d0) m = -m
      elseif (len_trim(msg) == 0) then
         msg = "The lattice vectors are not valid"
      end if

    end subroutine pending_cell

    ! The cell that the valid supercell transformation m (see
    ! pending_cell) leads to: its size and number of atoms (str), cell
    ! lengths (xl, bohr) and angles (xa, degrees).
    subroutine pending_metrics(m,str,xl,xa)
      use tools_math, only: gcd
      real*8, intent(in) :: m(3,3)
      character(kind=c_char,len=:), allocatable, intent(inout) :: str
      real*8, intent(out) :: xl(3), xa(3)

      real*8 :: d, r(3,3)
      integer :: num, den, ig

      ! the size of the new cell, a fraction of the current one if it
      ! is smaller than a cell of a centered lattice
      d = abs(det3(m))
      den = max(sys(isys)%c%ncv,1)
      if (abs(d * den - nint(d * den)) < 1d-6) then
         num = nint(d * den)
         ig = gcd(num,den)
         str = string(num / ig)
         if (den / ig > 1) str = str // "/" // string(den / ig)
      else
         str = string(d,'f',decimal=4)
      end if
      if (d < 1d0 + 1d-6) then
         str = str // " cell, "
      else
         str = str // " cells, "
      end if
      str = str // string(nint(d * sys(isys)%c%ncel)) // " atoms"

      r = matmul(sys(isys)%c%m_x2c,m)
      call cellpar_from_metric(matmul(transpose(r),r),xl,xa)

    end subroutine pending_metrics

    ! Select supercell candidate ic (index in the search results), of
    ! size nsz: it becomes the transformation of the search mode, and its
    ! lattice vectors are copied to the matrix mode, where they can be
    ! modified.
    subroutine select_candidate(ic,nsz)
      integer, intent(in) :: ic, nsz

      w%geometry_cell%isel = ic
      w%geometry_cell%nsel = nsz
      w%geometry_cell%intmat = w%geometry_cell%cand(ic)%m
      w%geometry_cell%cen = 1
      w%geometry_cell%origin = 0d0
      w%geometry_cell%inverse = .false.
      w%errmsg = ""

    end subroutine select_candidate

    ! Begin one of the tables of supercells of the cell tab (the search
    ! results and the sublattices of one size), with nrow rows, at most
    ! maxrow of them in view: set up the columns and draw the header
    ! row. first is the header of the first column and ttfirst its
    ! tooltip. Returns whether the table is open; if it is, the caller
    ! draws the rows and ends it.
    function cell_table_begin(strid,first,ttfirst,nrow,maxrow) result(isopen)
      character(len=*), intent(in) :: strid, first, ttfirst
      integer, intent(in) :: nrow, maxrow
      logical :: isopen

      integer(c_int) :: tflags
      type(ImVec2) :: tsz
      character(kind=c_char,len=:), allocatable, target :: tid

      tflags = ImGuiTableFlags_RowBg
      tflags = ior(tflags,ImGuiTableFlags_Borders)
      tflags = ior(tflags,ImGuiTableFlags_ScrollX)
      tflags = ior(tflags,ImGuiTableFlags_ScrollY)
      tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
      tid = strid // c_null_char
      tsz%x = 0._c_float
      tsz%y = (min(max(nrow,1),maxrow) + 1) * (igGetTextLineHeight() + 2._c_float * g%Style%CellPadding%y)
      isopen = igBeginTable(c_loc(tid),8,tflags,tsz,0._c_float)
      if (.not.isopen) return

      call iw_table_column(first,id=0)
      call iw_table_column("Rmax (Å)",id=1)
      call iw_table_column("Niceness",id=2)
      call iw_table_column("Reach (Å)",id=3)
      call iw_table_column("Nops",id=4)
      call iw_table_column("Nindep",id=5)
      call iw_table_column("Ndisp",id=6)
      call iw_table_column("Transformation",id=7)
      call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown,tooltips=&
         ttfirst // c_null_char //&
         "Radius of the largest sphere that fits in the supercell" // c_null_char //&
         "Volume of the cube around that sphere divided by the volume of the supercell: 1 for a&
         & cube, lower for a skewed or elongated cell" // c_null_char //&
         "Half the length of the shortest lattice vector of the supercell: two atoms closer than&
         & this are never periodic images of each other" // c_null_char //&
         "Symmetry operations of the crystal that the supercell keeps" // c_null_char //&
         "Independent atoms that a finite-displacement phonon calculation displaces in the&
         & supercell" // c_null_char //&
         "Displaced structures that a finite-displacement phonon calculation needs in the&
         & supercell" // c_null_char //&
         "Lattice vectors of the supercell in crystallographic coordinates of the current cell&
         & (the arguments of NEWCELL)" // c_null_char)

    end function cell_table_begin

    ! The columns after the first of supercell candidate ic (index in
    ! the search results), of size nsz, in the current row of one of the
    ! tables of supercells of the cell tab.
    subroutine cell_candidate_columns(ic,nsz)
      integer, intent(in) :: ic, nsz

      integer :: k1, k2
      character(kind=c_char,len=:), allocatable :: sm

      if (igTableSetColumnIndex(1)) &
         call iw_text(string(w%geometry_cell%cand(ic)%r*bohrtoa,'f',7,3))
      if (igTableSetColumnIndex(2)) &
         call iw_text(string(8d0*w%geometry_cell%cand(ic)%r**3/(nsz*sys(isys)%c%omega),'f',7,5))
      if (igTableSetColumnIndex(3)) then
         if (w%geometry_cell%cand(ic)%reach > 0d0) &
            call iw_text(string(w%geometry_cell%cand(ic)%reach*bohrtoa,'f',7,3))
      end if
      if (igTableSetColumnIndex(4)) &
         call iw_text(string(w%geometry_cell%cand(ic)%nops))
      if (igTableSetColumnIndex(5)) then
         if (w%geometry_cell%cand(ic)%nindep >= 0) &
            call iw_text(string(w%geometry_cell%cand(ic)%nindep))
      end if
      if (igTableSetColumnIndex(6)) then
         if (w%geometry_cell%cand(ic)%ndisp >= 0) &
            call iw_text(string(w%geometry_cell%cand(ic)%ndisp))
      end if
      if (igTableSetColumnIndex(7)) then
         sm = ""
         do k1 = 1, 3
            do k2 = 1, 3
               sm = sm // string(w%geometry_cell%cand(ic)%m(k2,k1),length=3,justify=ioj_right)
            end do
         end do
         call iw_text(sm)
      end if

    end subroutine cell_candidate_columns

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

    ! Push the selected (and hovered) symmetry elements to the view. nop is
    ! the number of operations shown in the table.
    subroutine symop_display(nop)
      integer, intent(in) :: nop

      integer :: i, n, tag, hovadd
      integer, allocatable :: sop(:)

      ! display the selected and hovered symmetry elements in the view; the
      ! element geometry (kind, direction, and position in the displayed cells)
      ! is calculated from the operation when the draw lists are built
      if (sysc(isys)%sc%isinit /= 0) then
         allocate(sop(nop))
         n = 0
         do i = 1, nop
            if (.not.(w%geometry_sym_sel(i).or.i == ihl_symop)) cycle
            n = n + 1
            sop(n) = i
         end do
         ! unique tag: the selection generation combined with the hovered row
         ! (only when it adds to the selection), so the scene rebuilds the
         ! element list only when the displayed set actually changes
         hovadd = 0
         if (ihl_symop > 0 .and. ihl_symop <= nop) then
            if (.not.w%geometry_sym_sel(ihl_symop)) hovadd = ihl_symop
         end if
         tag = w%geometry_sym_selgen * (nop + 2) + (hovadd + 1)
         call sysc(isys)%sc%show_transient_symelems(w%uid,tag,n,sop(1:n))
         deallocate(sop)
      end if

    end subroutine symop_display

    ! the selection row of the symmetry operations table (nop operations):
    ! all/none/toggle, the object creation, and the Delete transform
    subroutine symop_buttons(nop)
      integer, intent(in) :: nop

      integer :: i, idobj, nsel
      logical, allocatable :: smask(:)
      logical :: idnew

      ! Selection row: all/none/toggle buttons (elements use the default colors)
      call iw_text("Selection",alignframe=.true.)
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
      ! the object goes into the anchor view's scene, which is the one this
      ! window commands and the one the editor will be bound to; for an
      ! alternate view that is its private scene, not sysc(isys)%sc
      if (iw_button("Create Object##symcreateobj",sameline=.true.,&
         disabled=(.not.associated(win(iview)%sc)))) then
         idobj = 0
         associate (sc => win(iview)%sc)
           do i = 1, sc%nrep
              if (sc%rep(i)%isinit .and. sc%rep(i)%type == reptype_symelem) then
                 idobj = i
                 exit
              end if
           end do
           idnew = (idobj == 0)
           if (idnew) &
              call sc%add_representation(reptype_symelem,repflavor_symelem,id=idobj)
           if (idobj > 0) then
              ! ensure the operation snapshot/visibility style is initialized
              call sc%rep(idobj)%update()
              associate (rr => sc%rep(idobj))
                ! the object's visibility list is per element type, while the
                ! table selection is per symmetry operation, so translate
                if (rr%symelem%style%isinit .and. rr%symelem%style%se%ntype > 0) then
                   allocate(smask(rr%symelem%style%se%ntype))
                   call symelem_type_mask(rr%symelem%style%se,&
                      pack((/(i,i=1,nop)/),w%geometry_sym_sel),smask)
                   if (idnew) then
                      ! new object: show exactly the selected operations
                      rr%symelem%style%shown = smask
                   else
                      ! reuse: mark the selected operations as shown in it
                      where (smask) rr%symelem%style%shown = .true.
                   end if
                   deallocate(smask)
                   sc%forcebuildlists = .true.
                   call sc%undo_note("Show symmetry elements")
                end if
              end associate
           end if
         end associate
         if (idobj > 0) then
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
         call iw_text("Transform",alignframe=.true.)

         nsel = 0
         do i = 2, min(nop,size(w%geometry_sym_sel,1))
            if (w%geometry_sym_sel(i)) nsel = nsel + 1
         end do
         if (iw_button("Delete##symseldelete",sameline=.true.,disabled=(nsel==0))) &
            iaction = iaction_sym_deleteops
         call iw_tooltip("Delete the selected operations and rebuild with a maximal "//&
            "subgroup that excludes them.",ttshown)
      end if

    end subroutine symop_buttons

    ! deallocate the results of the supercell search (stale after any
    ! change of the cell or of the search criterion), and the selection
    subroutine clear_search_results()

      if (allocated(w%geometry_cell%nc)) deallocate(w%geometry_cell%nc)
      if (allocated(w%geometry_cell%ic0)) deallocate(w%geometry_cell%ic0)
      if (allocated(w%geometry_cell%ibest)) deallocate(w%geometry_cell%ibest)
      if (allocated(w%geometry_cell%cand)) deallocate(w%geometry_cell%cand)
      w%geometry_cell%isel = 0
      w%geometry_cell%nsel = 0

    end subroutine clear_search_results

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
      w%timelast_sort = glfwGetTime()
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
      if (icolsort(w%sortcid)==ic_id .or. icolsort(w%sortcid)==ic_zat .or.&
         icolsort(w%sortcid)==ic_mol .or. icolsort(w%sortcid)==ic_mul .or.&
         icolsort(w%sortcid)==ic_idx .or. icolsort(w%sortcid)==ic_nat) then
         ! integers
         allocate(ival(ntype))
         do ii = 1, ntype
            i = w%iord(ii)
            if (table_hltype == atlisttype_nmol) then
               if (icolsort(w%sortcid) == ic_id) then
                  ival(ii) = i
               elseif (icolsort(w%sortcid) == ic_nat) then
                  ival(ii) = sys(isys)%c%mol(i)%nat
               elseif (icolsort(w%sortcid) == ic_idx) then
                  ival(ii) = sys(isys)%c%idxmol(i)
               end if
            else
               if (icolsort(w%sortcid) == ic_id) then
                  ival(ii) = i
               elseif (icolsort(w%sortcid) == ic_zat) then
                  ispc = sysc(isys)%attype_species(w%geometry_atomtype,i)
                  ival(ii) = sys(isys)%c%spc(ispc)%z
               elseif (icolsort(w%sortcid) == ic_mol) then
                  ival(ii) = sys(isys)%c%idatcelmol(1,i)
               elseif (icolsort(w%sortcid) == ic_mul) then
                  ival(ii) = sys(isys)%c%at(i)%mult
               elseif (icolsort(w%sortcid) == ic_idx) then
                  ival(ii) = sys(isys)%c%atcel(i)%idx
               end if
            end if
         end do
         call mergesort(ival,iperm,1,ntype)
         deallocate(ival)
      elseif (icolsort(w%sortcid)==ic_x.or.icolsort(w%sortcid)==ic_y.or.&
         icolsort(w%sortcid)==ic_z.or.icolsort(w%sortcid)==ic_ea.or.&
         icolsort(w%sortcid)==ic_eb.or.icolsort(w%sortcid)==ic_eg.or.&
         icolsort(w%sortcid)==ic_occ) then
         ! real
         allocate(rval(ntype))
         do ii = 1, ntype
            i = w%iord(ii)
            if (icolsort(w%sortcid)==ic_ea.or.icolsort(w%sortcid)==ic_eb.or.&
               icolsort(w%sortcid)==ic_eg) then
               ! Euler angles (molecules tab)
               x0 = mol_euler_angles(i)
               if (icolsort(w%sortcid) == ic_ea) then
                  rval(ii) = x0(1)
               elseif (icolsort(w%sortcid) == ic_eb) then
                  rval(ii) = x0(2)
               else
                  rval(ii) = x0(3)
               end if
            elseif (icolsort(w%sortcid) == ic_occ) then
               ! occupancy
               rval(ii) = sysc(isys)%attype_occupancy(w%geometry_atomtype,i)
            else
               if (table_hltype == atlisttype_nmol) then
                  x0 = mol_com_coords(i)
               else
                  x0 = sysc(isys)%attype_coordinates(w%geometry_atomtype,i)
               end if

               if (icolsort(w%sortcid) == ic_x) then
                  rval(ii) = x0(1)
               elseif (icolsort(w%sortcid) == ic_y) then
                  rval(ii) = x0(2)
               elseif (icolsort(w%sortcid) == ic_z) then
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
            if (icolsort(w%sortcid) == ic_atom) then
               sval(ii)%s = sysc(isys)%attype_name(w%geometry_atomtype,i)
            elseif (icolsort(w%sortcid) == ic_wyc) then
               sval(ii)%s = string(sys(isys)%c%at(i)%mult) // sys(isys)%c%at(i)%wyc
            end if
         end do
         call mergesort(sval,iperm,1,ntype)
         deallocate(sval)
      end if

      ! reverse if necessary
      if (w%sortdir == 2) then
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
      logical :: isspc

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

      ! the species view is the only one whose rows can hold no atoms
      if (table_hltype == atlisttype_species) then
         if (allocated(w%geometry_spcsel)) then
            if (size(w%geometry_spcsel,1) /= ntype) deallocate(w%geometry_spcsel)
         end if
         if (.not.allocated(w%geometry_spcsel)) then
            allocate(w%geometry_spcsel(ntype))
            w%geometry_spcsel = .false.
         end if
      end if

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

      ! classify each row: 2 = fully selected, 1 = partially selected, 0 =
      ! unselected
      isspc = (table_hltype == atlisttype_species) .and. allocated(w%geometry_spcsel)
      do i = 1, ntype
         if (rowntot(i) > 0 .and. nselrow(i) == rowntot(i)) then
            rowstate(i) = 2
         elseif (nselrow(i) > 0) then
            rowstate(i) = 1
         elseif (isspc) then ! rowntot(i) == 0 here
            if (w%geometry_spcsel(i)) then
               rowstate(i) = 2
               rowrgba(:,i) = w%geometry_select_rgba
            end if
         end if
      end do

    end subroutine build_aggregate_state

    ! Whether row i of the current table is a species row with no atoms,
    ! whose selection the window tracks instead of the atom highlight.
    function empty_row(i)
      integer, intent(in) :: i
      logical :: empty_row

      ! the aggregate arrays outlive the tab that built them, so neither
      ! their presence nor their size can be assumed here
      empty_row = (table_hltype == atlisttype_species) .and. allocated(rowntot) .and.&
         allocated(w%geometry_spcsel)
      if (empty_row) empty_row = (i >= 1 .and. i <= size(rowntot,1) .and.&
         i <= size(w%geometry_spcsel,1))
      if (empty_row) empty_row = (rowntot(i) == 0)

    end function empty_row

    ! Whether any row with no atoms is selected. build_aggregate_state has
    ! already folded the window's flags into rowstate, so this reads the same
    ! state the table drew rather than a second copy of it.
    function any_empty_row_selected()
      logical :: any_empty_row_selected

      ! only while the species tab is the one being drawn: the Bonds and
      ! Symmetry tabs build no aggregate arrays and would read the leftovers
      any_empty_row_selected = .false.
      if (w%tabselected /= "species") return
      if (allocated(rowstate) .and. allocated(rowntot)) &
         any_empty_row_selected = any(rowstate == 2 .and. rowntot == 0)

    end function any_empty_row_selected

    ! Deselect every species row that has no atoms.
    subroutine clear_empty_rows()

      if (allocated(w%geometry_spcsel)) w%geometry_spcsel = .false.

    end subroutine clear_empty_rows

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
      use utils, only: iw_table_sort_specs

      if (iw_table_sort_specs(w%sortcid,w%sortdir,0_c_int,ntype > 1) .and. ntype > 1) &
         forcesort = .true.

    end subroutine fetch_sort_specs

    ! color of atom iat (of atom-list type itype) from the first shown atoms
    ! representation in the current view (iview); returns .true. and fills rgbo
    ! when a color is found, .false. otherwise
    function color_from_view(itype,iat,rgbo) result(have)
      integer, intent(in) :: itype, iat
      real(c_float), intent(out) :: rgbo(3)
      logical :: have

      have = atom_view_rgb(iview,isys,itype,iat,rgbo)

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
      logical, allocatable :: tatom(:) ! the row is backed by atoms

      ! highlight color
      call iw_text("Selection",alignframe=.true.)
      ldum = iw_coloredit("##drawgeometryhighlightcolor",rgba=w%geometry_select_rgba,sameline=.true.)
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
         ! only the rows with no atoms: the rest are in the atom highlight.
         ! empty_row carries the guards -- this button is shared by the tabs
         ! whose tables have no such rows and no aggregate arrays at all
         do i = 1, ntype
            if (empty_row(i)) w%geometry_spcsel(i) = .true.
         end do
      end if
      call iw_tooltip("Select all atoms in the system",ttshown)
      if (iw_button("None##highlightnone",sameline=.true.)) then
         call sysc(isys)%highlight_clear(.false.)
         call clear_empty_rows()
      end if
      call iw_tooltip("Deselect all atoms (" // trim(get_bind_keyname(BIND_CANCEL)) // ")",ttshown)
      if (iw_button("Toggle##highlighttoggle",sameline=.true.)) then
         ! compute the per-row selection state once
         ! state of every row, and the rows with no atoms flip their own flag
         ! (they are not in the atom highlight, so they are filtered out of
         ! the two calls below and counted out of them as well)
         allocate(tstate(ntype),tatom(ntype))
         do i = 1, ntype
            call current_row_state(i,tstate(i),rrgba)
            tatom(i) = .not.empty_row(i)
            if (.not.tatom(i)) w%geometry_spcsel(i) = (tstate(i) /= 2)
         end do
         ! select the rows that are not fully selected
         nhigh = count(tstate /= 2 .and. tatom)
         if (nhigh > 0) then
            allocate(ihigh(nhigh),irgba(4,nhigh))
            nhigh = 0
            do i = 1, ntype
               if (tstate(i) /= 2 .and. tatom(i)) then
                  nhigh = nhigh + 1
                  ihigh(nhigh) = i
                  irgba(:,nhigh) = w%geometry_select_rgba
               end if
            end do
            call sysc(isys)%highlight_atoms(.false.,ihigh,table_hltype,irgba)
            deallocate(ihigh,irgba)
         end if
         ! clear the rows that were fully selected
         nhigh = count(tstate == 2 .and. tatom)
         if (nhigh > 0) then
            allocate(ihigh(nhigh))
            nhigh = 0
            do i = 1, ntype
               if (tstate(i) == 2 .and. tatom(i)) then
                  nhigh = nhigh + 1
                  ihigh(nhigh) = i
               end if
            end do
            call sysc(isys)%highlight_clear(.false.,ihigh,table_hltype)
            deallocate(ihigh)
         end if
         deallocate(tstate,tatom)
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

      ! if the selected species does not exist, default to C
      if (w%geometry_input_species > sys(isys)%c%nspc) &
         w%geometry_input_species = -6

      ! the input position is interpreted in the same coordinate type as
      ! shown in the table, so the new atom's row matches the user input
      call iw_text("Position (" // sysc(isys)%attype_coordinates_units(w%geometry_atomtype) // ")",alignframe=.true.)
      ldum = iw_dragfloat_real8("##xaddcoord",x1=w%geometry_input_coord(1),speed=0.001d0,decimal=6,&
         notlive=.true.,sameline=.true.)
      ldum = iw_dragfloat_real8("##yaddcoord",x1=w%geometry_input_coord(2),speed=0.001d0,decimal=6,&
         notlive=.true.,sameline=.true.)
      ldum = iw_dragfloat_real8("##zaddcoord",x1=w%geometry_input_coord(3),speed=0.001d0,decimal=6,&
         notlive=.true.,sameline=.true.)

      call iw_text("Species",alignframe=.true.)
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
         if (iw_beginmenu("New")) then
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

      call iw_text("Edit",alignframe=.true.)

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

      ! Duplicate button. A selected species with no atoms counts: it is not
      ! in the atom highlight, but the edit buttons act on it all the same
      havesel = .false.
      if (allocated(sysc(isys)%highlight_rgba)) &
         havesel = any(sysc(isys)%highlight_rgba >= 0._c_float)
      havesel = havesel .or. any_empty_row_selected()
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

      call iw_text("Edit",alignframe=.true.)

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
