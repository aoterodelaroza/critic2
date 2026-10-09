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

! Routines for the edit representation window.
submodule (windows) editrep
  use interfaces_cimgui
  implicit none

contains

  !> Update tasks for the edit representation window, before the
  !> window is created.
  module subroutine update_editrep(w)
    use systems, only: sys_init, ok_system, sysc
    use representations, only: reptype_gpaths
    use windows, only: win
    class(window), intent(inout), target :: w

    integer :: isys, iview, i
    logical :: doquit

    ! This window edits one object living in its anchor view's scene, so unlike
    ! the other view tools it does not follow the view to another system: the
    ! object belongs to the scene it was created in. Check the anchor, the
    ! system and the object are all still there, and close otherwise.
    isys = w%isys
    iview = w%anchor_view()
    doquit = (iview == 0)
    if (.not.doquit) doquit = .not.ok_system(isys,sys_init)
    ! the scene has to be checked before the representation: the view
    ! owns the list rep points into, and tearing it down (another
    ! system in this view, a scene reset) does not clear the pointer
    if (.not.doquit) doquit = .not.associated(win(iview)%sc)
    if (.not.doquit) doquit = (win(iview)%isys /= isys)
    if (.not.doquit) doquit = .not.associated(w%rep)
    if (.not.doquit) doquit = .not.w%rep%isinit
    if (.not.doquit) doquit = (w%rep%type <= 0)

    ! if they aren't, quit the window; a gradient path highlighted from
    ! the table of this editor goes back to normal (the object may be
    ! gone, so the system's scene is swept instead of w%rep)
    if (doquit) then
       if (ok_system(isys,sys_init)) then
          do i = 1, sysc(isys)%sc%nrep
             if (.not.sysc(isys)%sc%rep(i)%isinit) cycle
             if (sysc(isys)%sc%rep(i)%type /= reptype_gpaths) cycle
             if (all(sysc(isys)%sc%rep(i)%gpaths%ihover == 0)) cycle
             sysc(isys)%sc%rep(i)%gpaths%ihover = 0
             sysc(isys)%sc%forcebuildlists = .true.
          end do
       end if
       call w%end()
    end if

  end subroutine update_editrep

  !> Draw the edit represenatation window.
  module subroutine draw_editrep(w)
    use representations, only: representation, reptype_atoms, reptype_bonds, reptype_labels,&
       reptype_polyhedra, reptype_unitcell, reptype_axes, reptype_symelem, reptype_text,&
       reptype_measure, reptype_isosurface, reptype_shapes, reptype_cps, reptype_gpaths, iso_map_color,&
       reptype_planar, reptype_legend
    use windows, only: win
    use keybindings, only: is_bind_event, BIND_OK_FOCUSED_DIALOG
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_text, iw_tooltip, iw_button, iw_calcheight, iw_checkbox,&
       iw_inputtext, iw_close_event, iw_setpos_bottomright, iw_periodicity_widget
    use grid3mod, only: hscale_num, hscale_log, hscale_asinh
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: isys, iview
    logical :: doquit, ok, doper, ch
    logical :: changed, renamed

    logical, save :: ttshown = .false. ! tooltip flag

    ! check the anchor, the system and the object are all still there (see
    ! update_editrep); the last test catches the anchor view moving to another
    ! system, which leaves this object behind in the old scene
    isys = w%isys
    iview = w%anchor_view()
    doquit = (iview == 0)
    if (.not.doquit) doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = .not.win(iview)%isopen
    ! the scene has to be checked before the representation (see update_editrep)
    if (.not.doquit) doquit = .not.associated(win(iview)%sc)
    if (.not.doquit) doquit = (win(iview)%isys /= isys)
    if (.not.doquit) doquit = .not.associated(w%rep)
    if (.not.doquit) doquit = .not.w%rep%isinit
    if (.not.doquit) doquit = (w%rep%type <= 0)

    if (.not.doquit) then
       ! whether the rep has changed
       changed = .false.

       ! system
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! name block (the representation type is fixed at creation and cannot be
       ! changed here)
       call iw_text("Name",highlight=.true.,alignframe=.true.)
       renamed = iw_inputtext("##nametextinput",bufsize=1023,texta=w%rep%name,width=30,sameline=.true.)
       call iw_tooltip("Name of this object",ttshown)

       ! shown checkbox
       changed = changed .or. iw_checkbox("Show",w%rep%shown,sameline=.true.)
       call iw_tooltip("Toggle show/hide this object",ttshown)

       ! periodicity override, for the objects drawn over the cells of the
       ! Display (atoms, unit cell, and a whole-cell isosurface) in crystals
       doper = .not.sys(isys)%c%ismolecule
       if (doper) doper = w%rep%uses_periodicity()
       if (doper) then
          ch = iw_periodicity_widget(w%rep%disp%ncell,ttshown,pertype=w%rep%disp%pertype)
          changed = changed .or. ch
       end if

       ! update the representation to respond to changes in the number of
       ! atoms and molecules, then draw the type-dependent items
       call w%rep%update()
       if (w%rep%type == reptype_atoms) then
          changed = changed .or. w%draw_editrep_atoms(ttshown)
       elseif (w%rep%type == reptype_bonds) then
          changed = changed .or. w%draw_editrep_bonds(ttshown)
       elseif (w%rep%type == reptype_labels) then
          changed = changed .or. w%draw_editrep_labels(ttshown)
       elseif (w%rep%type == reptype_polyhedra) then
          changed = changed .or. w%draw_editrep_polyhedra(ttshown)
       elseif (w%rep%type == reptype_unitcell) then
          changed = changed .or. w%draw_editrep_unitcell(ttshown)
       elseif (w%rep%type == reptype_axes) then
          changed = changed .or. w%draw_editrep_axes(ttshown)
       elseif (w%rep%type == reptype_symelem) then
          changed = changed .or. w%draw_editrep_symelem(ttshown)
       elseif (w%rep%type == reptype_text) then
          changed = changed .or. w%draw_editrep_text(ttshown)
       elseif (w%rep%type == reptype_measure) then
          changed = changed .or. w%draw_editrep_measure(ttshown)
       elseif (w%rep%type == reptype_isosurface) then
          changed = changed .or. w%draw_editrep_isosurface(ttshown)
       elseif (w%rep%type == reptype_shapes) then
          changed = changed .or. w%draw_editrep_shapes(ttshown)
       elseif (w%rep%type == reptype_planar) then
          changed = changed .or. w%draw_editrep_planar(ttshown)
       elseif (w%rep%type == reptype_cps) then
          changed = changed .or. w%draw_editrep_cps(ttshown)
       elseif (w%rep%type == reptype_gpaths) then
          changed = changed .or. w%draw_editrep_gpaths(ttshown)
       elseif (w%rep%type == reptype_legend) then
          changed = changed .or. w%draw_editrep_legend(ttshown)
       end if

       ! rebuild draw lists if necessary
       if (changed) win(iview)%sc%forcebuildlists = .true.
       if (changed) then
          call win(iview)%sc%undo_note("Edit " // trim(w%rep%name))
       elseif (renamed) then
          call win(iview)%sc%undo_note("Rename object")
       end if

       ! right-align and bottom-align for the rest of the contents
       call iw_setpos_bottomright(10,2)

       ! reset button
       if (iw_button("Reset",danger=.true.)) then
          call w%rep%set_defaults(0)
          win(iview)%sc%forcebuildlists = .true.
          call win(iview)%sc%undo_note("Reset " // trim(w%rep%name))
       end if

       ! close button
       ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
       ok = ok .or. iw_button("Close",sameline=.true.)
       if (ok) doquit = .true.
    end if

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window; a mapped-isosurface line still held commits
    ! its pending rebuild here, because the release edge never fires once
    ! the window is gone
    if (doquit) then
       if (w%editrep_isoline >= 1 .and. iview > 0 .and. associated(w%rep)) then
          if (w%rep%isinit .and. w%rep%type == reptype_isosurface .and.&
             win(iview)%isopen .and. associated(win(iview)%sc)) then
             if (w%editrep_isoline <= w%rep%iso%niso) then
                if (w%rep%iso%slot(w%editrep_isoline)%imap_mode /= iso_map_color) &
                   win(iview)%sc%forcebuildlists = .true.
             end if
          end if
       end if

       ! a gradient path highlighted from the table goes back to normal
       if (iview > 0 .and. associated(w%rep)) then
          if (w%rep%isinit .and. w%rep%type == reptype_gpaths .and. any(w%rep%gpaths%ihover /= 0)) then
             w%rep%gpaths%ihover = 0
             if (win(iview)%isopen .and. associated(win(iview)%sc)) win(iview)%sc%forcebuildlists = .true.
          end if
       end if
       call w%end()
    end if

  end subroutine draw_editrep

  !> Draw the editrep (Object) window, atoms class. Returns true if
  !> the scene needs rendering again. ttshown = the tooltip flag.
  !> Draw the editrep window, atoms class. Returns true if the
  !> representation has changed.
  module function draw_editrep_atoms(w,ttshown) result(changed)
    use systems, only: sys, sysc, atlisttype_species
    use gui_main, only: ColorHighlightScene, ColorElement
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_button,&
       iw_checkbox, iw_coloredit, iw_dragfloat_real8
    use param, only: atmcov, atmvdw, jmlcol, jmlcol2, bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: ispc, isys, iz
    logical :: ch, typechanged
    integer :: i, ihighlight, highlight_type

    ! initialize
    ihighlight = 0
    highlight_type = atlisttype_species
    changed = .false.
    isys = w%isys

    ! global options for atoms
    call iw_text("Global Options",highlight=.true.,alignframe=.true.)
    if (iw_button("Reset##resetglobalatoms",sameline=.true.,danger=.true.)) then
       call w%rep%set_defaults(1)
       changed = .true.
    end if
    call iw_tooltip("Reset to the default settings for the atom representation")
    call iw_combo_simple("Radii ##atomradiicombo","Covalent"//c_null_char//"Van der Waals"//c_null_char//&
       "Constant"//c_null_char,w%rep%atoms%radii_type,changed=ch)
    call iw_tooltip("Set atomic radii to the tabulated values of this type",ttshown)

    if (w%rep%atoms%radii_type == 2) then
       ! constant size
       ch = ch .or. iw_dragfloat_real8("Value##atomradii",x1=w%rep%atoms%radii_value,speed=0.01d0,&
          min=0d0,max=5d0,scale=bohrtoa,decimal=3,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Atomic radii (Å)",ttshown)

       if (ch) then
          w%rep%atoms%style%rad(1:w%rep%atoms%style%ntype) = w%rep%atoms%radii_value
          changed = .true.
       end if
    else
       ! variable size
       ch = ch .or. iw_dragfloat_real8("Scale##atomradiiscale",x1=w%rep%atoms%radii_scale,speed=0.01d0,&
          min=0d0,max=5d0,decimal=3,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Scale factor for the tabulated atomic radii",ttshown)

       if (ch) then
          do i = 1, w%rep%atoms%style%ntype
             ispc = sysc(isys)%attype_species(w%rep%atoms%style%type,i)
             iz = sys(isys)%c%spc(ispc)%z
             if (w%rep%atoms%radii_type == 0) then
                w%rep%atoms%style%rad(i) = atmcov(iz)
             else
                w%rep%atoms%style%rad(i) = atmvdw(iz)
             end if
             w%rep%atoms%style%rad(i) = w%rep%atoms%style%rad(i) * w%rep%atoms%radii_scale
          end do
          changed = .true.
       end if
    end if

    ! style buttons: set color
    call iw_combo_simple("Colors ##atomcolorselect","Current defaults" // c_null_char //&
       "jmol (light)" // c_null_char // "jmol2 (dark)" // c_null_char,&
       w%rep%atoms%color_type,changed=ch)
    call iw_tooltip("Set the color of all atoms to the tabulated values",ttshown)
    if (ch) then
       do i = 1, w%rep%atoms%style%ntype
          ispc = sysc(isys)%attype_species(w%rep%atoms%style%type,i)
          iz = sys(isys)%c%spc(ispc)%z
          if (w%rep%atoms%color_type == 0) then
             w%rep%atoms%style%rgb(:,i) = ColorElement(:,iz)
          elseif (w%rep%atoms%color_type == 1) then
             w%rep%atoms%style%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
          else
             w%rep%atoms%style%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
          end if
       end do
       changed = .true.
    end if

    ! border size
    changed = changed .or. iw_dragfloat_real8("Border Size (Å)",x1=w%rep%atoms%border_size,&
       speed=0.002d0,min=0d0,max=1d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Change the thickness of the atom borders",ttshown)

    ! color
    changed = changed .or. iw_coloredit("Border Color",rgb=w%rep%atoms%border_rgb,sameline=.true.)
    call iw_tooltip("Color of the border for the atoms",ttshown)

    ! the atoms at the centers of coordination polyhedra
    changed = changed .or. iw_checkbox("Shrink atoms to fit inside polyhedra",w%rep%atoms%fitpoly)
    call iw_tooltip("Draw the atoms at the centers of the coordination polyhedra small &
       &enough to fit inside them (but never smaller than about a third of their size)",ttshown)

    ! occupancy sectors
    if (sys(isys)%c%haveocc) then
       changed = changed .or. iw_checkbox("Occupancy sectors",w%rep%atoms%occ_sectors)
       call iw_tooltip("Draw partially occupied sites as spheres with a filled "//&
          "sector proportional to the occupancy",ttshown)
       if (w%rep%atoms%occ_sectors) then
          changed = changed .or. iw_coloredit("Vacancy color",rgb=w%rep%atoms%occ_empty_rgb,sameline=.true.)
          call iw_tooltip("Color of the sector corresponding to a vacancy drawn on partially occupied atoms",&
             ttshown)
       end if
    end if

    ! the atom and molecule style tables (a change of grouping
    ! remakes the style arrays, and the table is drawn next frame)
    ch = atom_table_widget(isys,w%rep%atoms%style%type,typechanged,ihighlight,highlight_type,&
       rgb=w%rep%atoms%style%rgb,rad=w%rep%atoms%style%rad)
    if (typechanged) call w%rep%atoms%style%reset(w%rep)
    changed = changed .or. ch
    if (w%rep%mols%style%isinit) then
       ch = mol_table_widget(isys,ihighlight,highlight_type,&
          tint=w%rep%mols%style%tint_rgb,scale=w%rep%mols%style%scale_rad)
       changed = changed .or. ch
    end if

    ! process transient highlighs
    if (ihighlight > 0) then
       call sysc(isys)%highlight_atoms(.true.,(/ihighlight/),highlight_type,&
          reshape(ColorHighlightScene,(/4,1/)))
    end if

  end function draw_editrep_atoms

  !> Draw the editrep window, bonds class. Returns true if the
  !> representation has changed.
  module function draw_editrep_bonds(w,ttshown) result(changed)
    use systems, only: sys
    use tools_io, only: string
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_radiobutton, iw_calcheight, iw_checkbox, iw_coloredit,&
       iw_dragfloat_real8, iw_table_column
    use param, only: atmcov0, newline, bohrtoa
    use global, only: bondfactor_def, bonddelta_def
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: isys, iz
    character(kind=c_char,len=:), allocatable, target :: str1, suffix
    logical :: ch, ldum
    integer(c_int) :: flags, nspcpair
    integer :: i, j, k, nrow
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    integer, allocatable :: indi(:), indj(:)
    type(ImVec2) :: sz

    integer(c_int), parameter :: ic_sp1 = 0
    integer(c_int), parameter :: ic_sp2 = 1
    integer(c_int), parameter :: ic_shown = 2

    ! initialize
    changed = .false.
    isys = w%isys

    call iw_text("Global Options",highlight=.true.,alignframe=.true.)
    if (iw_button("Reset##resetglobal",sameline=.true.,danger=.true.)) then
       call w%rep%set_defaults(2)
       changed = .true.
    end if
    call iw_tooltip("Reset to the covalent bonding for this system and the default settings")

    ! rest of the options (record changes)
    ch = .false.
    call iw_text("Style",alignframe=.true.)
    call iw_combo_simple("##tablebondstyleglobalselect",&
       "Single color"//c_null_char//"Two colors"//c_null_char,w%rep%bonds%color_style,sameline=.true.,changed=ch,&
       tooltips="Draw each bond in a single color, the bond color"//c_null_char//&
       "Draw each bond in the colors of its two atoms, each on its side"//c_null_char,ttshown=ttshown)

    call iw_text(" Radius (Å)",sameline=.true.)
    ch = ch .or. iw_dragfloat_real8("##radiusbondtableglobal",x1=w%rep%bonds%rad,speed=0.005d0,&
       min=0d0,max=2d0,scale=bohrtoa,decimal=3,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Radius of the bonds",ttshown)

    ! border size
    ch = ch .or. iw_dragfloat_real8("Border Size (Å)",x1=w%rep%bonds%border_size,speed=0.002d0,&
       min=0d0,max=1d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Change the thickness of the bond borders",ttshown)

    ! color
    ch = ch .or. iw_coloredit("Border Color",rgb=w%rep%bonds%border_rgb,sameline=.true.)
    call iw_tooltip("Color of the border for the bonds",ttshown)

    ! color
    call iw_text("Color",alignframe=.true.)
    ch = ch .or. iw_coloredit("##colorbondtableglobal",rgb=w%rep%bonds%rgb,sameline=.true.)
    call iw_tooltip("Color of the bonds",ttshown)

    ! order
    call iw_text(" Order",sameline=.true.)
    call iw_combo_simple("##tablebondorderselectglobal",&
       "Dashed"//c_null_char//"Single"//c_null_char//"Double"//c_null_char//"Triple"//c_null_char//&
       "Calculated"//c_null_char,&
       w%rep%bonds%order,sameline=.true.,changed=ldum)
    ch = ch .or. ldum
    call iw_tooltip("Bond order: a fixed order for all bonds (dashed, single, double, triple) or the&
       & order determined by critic2 for each bond (calculated)",ttshown)

    ! both atoms
    call iw_text(" Both Atoms",sameline=.true.)
    ch = ch .or. iw_checkbox("##bothatomstableglobal",w%rep%bonds%bothends,sameline=.true.)
    call iw_tooltip("Represent a bond if both end-atoms are in the scene (checked) or if only &
       &one end-atom is in the scene (unchecked)",ttshown)

    ! Jeffrey-Steiner hydrogen-bond strength classification
    call iw_text("Hydrogen bonds",alignframe=.true.)
    ch = ch .or. iw_checkbox("##hbondclassify",w%rep%bonds%hbond_classify,sameline=.true.)
    call iw_tooltip("Represent only hydrogen bonds and color them by their Jeffrey-Steiner strength.",ttshown)

    if (w%rep%bonds%hbond_classify) then
       ch = ch .or. iw_coloredit("Strong##hbstrongcolor",rgb=w%rep%bonds%hbond_rgb(:,1))
       ch = ch .or. iw_coloredit("Moderate##hbmodcolor",rgb=w%rep%bonds%hbond_rgb(:,2),sameline=.true.)
       ch = ch .or. iw_coloredit("Weak##hbweakcolor",rgb=w%rep%bonds%hbond_rgb(:,3),sameline=.true.)

       call iw_text("H...A Distance (Å)",alignframe=.true.)
       ch = ch .or. iw_dragfloat_real8("strong|moderate##hbdist1",x1=w%rep%bonds%hbond_dist(1),speed=0.01d0,&
          min=0d0,max=5d0,scale=bohrtoa,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("H...A distance class boundaries",ttshown)
       ch = ch .or. iw_dragfloat_real8("moderate|weak##hbdist2",x1=w%rep%bonds%hbond_dist(2),speed=0.01d0,&
          min=0d0,max=5d0,scale=bohrtoa,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("H...A distance class boundaries",ttshown)

       call iw_text("D-H...A Angle (°)",alignframe=.true.)
       ch = ch .or. iw_dragfloat_real8("weak|moderate##hbang1",x1=w%rep%bonds%hbond_ang(1),speed=0.5d0,&
          min=0d0,max=180d0,decimal=1,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("D-H...A angle class boundaries",ttshown)
       ch = ch .or. iw_dragfloat_real8("moderate|strong##hbang2",x1=w%rep%bonds%hbond_ang(2),speed=0.5d0,&
          min=0d0,max=180d0,decimal=1,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("D-H...A angle class boundaries",ttshown)
    end if

    call iw_text("Atom Pair Selection",highlight=.true.)

    nspcpair = min(5,sys(isys)%c%nspc*(sys(isys)%c%nspc+1)/2+1)
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tablespeciesbonding" // c_null_char
    sz%x = iw_calcwidth(17,3)
    sz%y = iw_calcheight(nspcpair,0,.false.)
    if (igBeginTable(c_loc(str1),3,flags,sz,0._c_float)) then
       ! header setup
       call iw_table_column("Atom 1",id=ic_sp1,flags=ImGuiTableColumnFlags_WidthFixed)

       call iw_table_column("Atom 2",id=ic_sp2,flags=ImGuiTableColumnFlags_WidthFixed)

       call iw_table_column("Show",id=ic_shown,flags=ImGuiTableColumnFlags_WidthFixed)

       ! draw the header
       call iw_table_headers_row(freezetop=.true.,autofit=.true.)

       ! start the clipper
       nrow = sys(isys)%c%nspc * (sys(isys)%c%nspc + 1) / 2
       clipper = ImGuiListClipper_ImGuiListClipper()
       call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)
       allocate(indi(nrow),indj(nrow))
       k = 0
       do i = 1, sys(isys)%c%nspc
          do j = i, sys(isys)%c%nspc
             k = k + 1
             indi(k) = i
             indj(k) = j
          end do
       end do

       ! draw the rows
       do while(ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do k = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
             i = indi(k)
             j = indj(k)

             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             suffix = "_" // string(i) // "_" // string(j)

             ! species
             if (igTableSetColumnIndex(ic_sp1)) then
                call iw_text(trim(sys(isys)%c%spc(i)%name),alignframe=.true.)
             end if
             if (igTableSetColumnIndex(ic_sp2)) &
                call iw_text(trim(sys(isys)%c%spc(j)%name))

             ! shown
             if (igTableSetColumnIndex(ic_shown)) then
                if (iw_checkbox("##bondtableshown" // suffix,w%rep%bonds%style%shown(i,j))) then
                   ch = .true.
                   w%rep%bonds%style%shown(j,i) = w%rep%bonds%style%shown(i,j)
                end if
                call iw_tooltip("Toggle display of bonds connecting these atom types",ttshown)
             end if
          end do ! clipper range
       end do ! clipper step

       ! end the clipper and the table
       deallocate(indi,indj)
       call ImGuiListClipper_End(clipper)
       call ImGuiListClipper_destroy(clipper)
       call igEndTable()
    end if ! begintable

    ! style buttons: show/hide
    if (iw_button("Show All##showallbonds")) then
       w%rep%bonds%style%shown = .true.
       ch = .true.
    end if
    call iw_tooltip("Show all bonds in the system",ttshown)
    if (iw_button("Hide All##hideallbonds",sameline=.true.)) then
       w%rep%bonds%style%shown = .false.
       ch = .true.
    end if
    call iw_tooltip("Hide all bonds in the system",ttshown)
    if (iw_button("Toggle Show/Hide##toggleallatoms",sameline=.true.)) then
       do i = 1, sys(isys)%c%nspc
          do j = i, sys(isys)%c%nspc
             w%rep%bonds%style%shown(j,i) = .not.w%rep%bonds%style%shown(j,i)
             w%rep%bonds%style%shown(i,j) = w%rep%bonds%style%shown(j,i)
          end do
       end do
       ch = .true.
    end if
    call iw_tooltip("Toggle the show/hide status for all bonds",ttshown)

    !! recalculate bonds block !!
    call iw_text("Recalculate Bonds",highlight=.true.)

    ! choose between the system's bonds and per-representation custom bonds
    if (iw_radiobutton("System bonds",bool=w%rep%bonds%style%use_sys_nstar,boolval=.true.)) then
       ! switched back to system bonds: re-read the system connectivity and redraw
       call w%rep%bonds%style%copy_neighstars_from_system(w%rep%id)
       changed = .true.
    end if
    call iw_tooltip("Draw the bonds calculated for the system. Use view/edit geometry window to modify.",ttshown)
    ldum = iw_radiobutton("Custom bonds",bool=w%rep%bonds%style%use_sys_nstar,boolval=.false.,sameline=.true.)
    call iw_tooltip("Draw bonds computed for this representation with custom distance criteria",ttshown)

    ! recalculation controls, only shown for custom bonds (same
    ! criteria as the geometry window's Recalculate Bonds section)
    if (.not.w%rep%bonds%style%use_sys_nstar) then
       ! explanation of the bonding criteria
       call iw_text("Atoms A and B are bonded if:")
       call iw_text("  A and B non-metals: d < (r_cov(i)+r_cov(j))*f"//newline//&
          "  A or B metal:       d < d_NN + δ"//newline//&
          "r_cov = covalent radius. d_NN = nearest-neighbor distance.")

       ! per-species covalent radii table
       flags = ImGuiTableFlags_None
       flags = ior(flags,ImGuiTableFlags_Resizable)
       flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       flags = ior(flags,ImGuiTableFlags_ScrollX)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       str1 = "##tableatomrcov_editrep" // c_null_char
       sz%x = 0
       sz%y = iw_calcheight(min(5,sys(isys)%c%nspc)+1,0,.false.)
       if (igBeginTable(c_loc(str1),3,flags,sz,0._c_float)) then
          call iw_table_column("Atom",id=0)
          call iw_table_column("Z",id=1)
          call iw_table_column("Radius (Å)",id=2)
          call iw_table_headers_row(freezetop=.true.)

          do i = 1, sys(isys)%c%nspc
             iz = sys(isys)%c%spc(i)%z
             if (iz <= 0) cycle
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             if (igTableSetColumnIndex(0)) then
                call iw_text(trim(sys(isys)%c%spc(i)%name),alignframe=.true.)
             end if
             if (igTableSetColumnIndex(1)) then
                call iw_text(string(iz),alignframe=.true.)
             end if
             if (igTableSetColumnIndex(2)) then
                ldum = iw_dragfloat_real8("##tableradius_editrep" // string(i),x1=w%rep%bonds%atmrad(iz),&
                   speed=0.01d0,min=0d0,max=2.65d0,scale=bohrtoa,decimal=3,&
                   flags=ImGuiSliderFlags_AlwaysClamp)
             end if
          end do
          call igEndTable()
       end if

       ! bond factor and bond delta
       call igAlignTextToFramePadding()
       ldum = iw_dragfloat_real8("Bond factor (f)##bondfactor_editrep",x1=w%rep%bonds%bfactor,&
          speed=0.001d0,min=1d0,max=4d0,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Bond factor parameter (multiplicative) for non-metal bonding (see formula above)",ttshown)
       call iw_text(" ",sameline=.true.)
       ldum = iw_dragfloat_real8("Bond delta δ (Å)##bonddelta_editrep",x1=w%rep%bonds%bdelta,&
          speed=0.001d0,min=0d0,max=2d0,scale=bohrtoa,decimal=4,sameline=.true.,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Distance tolerance (additive) for metal bonding (see formula above)",ttshown)
       call iw_text(" ",sameline=.true.)

       ! reset button (same line as the factor/delta drags)
       if (iw_button("Reset##resetcustombond",sameline=.true.)) then
          w%rep%bonds%atmrad = atmcov0
          w%rep%bonds%bfactor = bondfactor_def
          w%rep%bonds%bdelta = bonddelta_def
       end if
       call iw_tooltip("Reset covalent radii, bond factor, and bond delta to defaults",ttshown)

       ! apply button
       if (iw_button("Apply##applyglobal",danger=.true.)) then
          call w%rep%bonds%style%generate_neighstars(w%rep)
          changed = .true.
       end if
       call iw_tooltip("Recalculate and draw bonds using the criteria above",ttshown)
    end if

    ! immediately update if non-distances have changed
    if (ch) changed = .true.

  end function draw_editrep_bonds

  !> Draw the editrep window, labels class. Returns true if the
  !> representation has changed.
  module function draw_editrep_labels(w,ttshown) result(changed)
    use systems, only: sys, sysc, atlisttype_species, atlisttype_ncel_ang, atlisttype_nmol,&
       atlisttype_nneq
    use gui_main, only: ColorHighlightScene
    use tools_io, only: string
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_calcheight, iw_checkbox, iw_coloredit, iw_highlight_selectable,&
       iw_dragfloat_real8, iw_inputtext, iw_table_column, iw_field_combo
    use representations, only: cps_field, cps_abbrev
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: isys
    character(kind=c_char,len=:), allocatable, target :: str1, suffix
    logical :: ch
    integer(c_int) :: lst, flags
    integer :: i, intable, nrow, is, ncol, ihighlight, highlight_type, ifield, ncp
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    type(ImVec2) :: sz

    integer(c_int), parameter :: lsttrans(0:7) = (/0,1,2,2,2,3,4,5/)
    integer(c_int), parameter :: lsttransi(0:5) = (/0,1,2,5,6,7/)

    ! initialize
    ihighlight = 0
    highlight_type = atlisttype_species
    changed = .false.
    isys = w%isys

    call iw_text("Global Options",highlight=.true.,alignframe=.true.)
    if (iw_button("Reset##resetglobal",sameline=.true.,danger=.true.)) then
       w%rep%labels%type = 0
       call w%rep%set_defaults(3)
       call w%rep%labels%style%reset(w%rep)
       changed = .true.
    end if
    call iw_tooltip("Reset to the labels to the default settings")

    if (sys(isys)%c%ismolecule) then
       lst = lsttrans(w%rep%labels%type)
       call iw_combo_simple("Text##labelcontentselect","Atomic symbol"//c_null_char//&
          "Atom name"// c_null_char//"Atom ID"// c_null_char//&
          "Species ID"// c_null_char// "Atomic number"// c_null_char// "Molecule ID"// c_null_char,&
          lst,changed=ch)
       w%rep%labels%type = lsttransi(lst)
    else
       call iw_combo_simple("Text##labelcontentselect","Atomic symbol"//c_null_char//&
          "Atom name"//c_null_char//"Cell atom ID"//c_null_char//&
          "Cell atom ID + lattice vector"//c_null_char//"Symmetry-unique atom ID"//c_null_char//&
          "Species ID"//c_null_char//"Atomic number"//c_null_char//"Molecule ID"//c_null_char//&
          "Wyckoff position"//c_null_char,&
          w%rep%labels%type,changed=ch)
    end if
    if (ch) call w%rep%labels%style%reset(w%rep)
    call iw_tooltip("Text to display in the atom labels",ttshown)
    changed = changed .or. ch

    ! scale, constant size, color
    changed = changed .or. iw_dragfloat_real8("Scale##labelscale",x1=w%rep%labels%scale,speed=0.01d0,&
       min=0d0,max=10d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Scale factor for the atom labels",ttshown)

    changed = changed .or. iw_checkbox("Constant size##labelconstsize",&
       w%rep%labels%const_size,sameline=.true.)
    call iw_tooltip("Labels have constant size (on) or labels scale with the&
       & size of the associated atom (off)",ttshown)

    changed = changed .or. iw_coloredit("Color##labelcolor",rgb=w%rep%labels%rgb,sameline=.true.)
    call iw_tooltip("Color of the atom labels",ttshown)

    ! offset
    changed = changed .or. iw_dragfloat_real8("Offset (Å)",x3=w%rep%labels%offset,&
       speed=0.001d0,min=99.999d0,max=99.999d0,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Offset the position of the labels relative to the atom center",ttshown)

    ! field whose critical points are labeled, if any field has them
    ! and the label type has CP rows (the rows follow in update_styles)
    if (any(w%rep%labels%type == (/0,1,2,3,4,8/)) .and. cps_field(isys) >= 0) then
       call iw_text("Critical points of",alignframe=.true.)
       call igSameLine(0._c_float,-1._c_float)
       ifield = w%rep%labels%ifield
       if (iw_field_combo("##labelcpfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
          nonestr="<field not available>")) then
          w%rep%labels%ifield = ifield
          changed = .true.
       end if
       call iw_tooltip("Field whose critical points can be labeled (rows below the atoms in the table)",ttshown)

       changed = changed .or. iw_dragfloat_real8("CP scale##labelscalecp",x1=w%rep%labels%scale_cp,&
          speed=0.01d0,min=0d0,max=10d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Scale factor for the critical point labels",ttshown)
    end if

    ! table for label selection
    call iw_text("Label Selection",highlight=.true.)

    ! number of entries in the table
    select case(w%rep%labels%type)
    case (0,5,6)
       intable = atlisttype_species
       nrow = sys(isys)%c%nspc
       ncol = 5
       call iw_text("(per species)",sameline=.true.)
    case (2,3)
       intable = atlisttype_ncel_ang
       nrow = sys(isys)%c%ncel
       ncol = 5
       call iw_text("(per atom)",sameline=.true.)
    case (1,4,8)
       intable = atlisttype_nneq
       nrow = sys(isys)%c%nneq
       ncol = 5
       call iw_text("(per symmetry-unique atom)",sameline=.true.)
    case (7)
       intable = atlisttype_nmol
       nrow = sys(isys)%c%nmol
       ncol = 3
       call iw_text("(per molecule)",sameline=.true.)
    end select

    ! the table itself
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tablespecieslabels" // c_null_char
    sz%x = iw_calcwidth(30,ncol)
    ncp = w%rep%labels%style%ncp
    sz%y = iw_calcheight(min(8,nrow+ncp+1),0,.false.)
    if (igBeginTable(c_loc(str1),ncol,flags,sz,0._c_float)) then
       ncol = -1

       ! header setup
       call iw_table_column("Id",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)

       if (intable /= atlisttype_nmol) then
          call iw_table_column("Atom",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)

          call iw_table_column("Z ",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
       end if

       call iw_table_column("Show",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)

       call iw_table_column("Text",icol=ncol,flags=ImGuiTableColumnFlags_WidthStretch)

       ! draw the header
       call iw_table_headers_row(freezetop=.true.,autofit=.true.)

       ! start the clipper
       clipper = ImGuiListClipper_ImGuiListClipper()
       call ImGuiListClipper_Begin(clipper,nrow+ncp,-1._c_float)

       ! draw the rows
       do while(ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd

             ! set up the next row
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             suffix = "_" // string(i)
             ncol = -1

             ! the critical point rows come after the atoms
             if (i > nrow) then
                call draw_cp_row(i-nrow)
                cycle
             end if

             ! id
             ncol = ncol + 1
             if (igTableSetColumnIndex(ncol)) then
                call iw_text(string(i),alignframe=.true.)

                ! the highlight selectable
                if (iw_highlight_selectable("##selectablelabeltable" // suffix)) then
                   ihighlight = i
                   highlight_type = intable
                end if
             end if

             is = sysc(isys)%attype_species(intable,i)

             ! atom
             if (intable /= atlisttype_nmol) then
                ncol = ncol + 1
                if (igTableSetColumnIndex(ncol)) &
                   call iw_text(trim(sys(isys)%c%spc(is)%name))

                ! Z
                ncol = ncol + 1
                if (igTableSetColumnIndex(ncol)) &
                   call iw_text(string(sys(isys)%c%spc(is)%z))
             end if

             ! shown
             ncol = ncol + 1
             if (igTableSetColumnIndex(ncol)) then
                changed = changed .or. iw_checkbox("##labeltableshown" // suffix,w%rep%labels%style%shown(i))
                call iw_tooltip("Toggle display of labels for these atoms/molecules",ttshown)
             end if

             ! text
             ncol = ncol + 1
             if (igTableSetColumnIndex(ncol)) then
                changed = changed .or. iw_inputtext("##labeltabletext" // string(i),bufsize=32,&
                   textf=w%rep%labels%style%str(i),width=15)
                call iw_tooltip("Text for the atomic labels",ttshown)
             end if
          end do ! table rows: clipper range
       end do ! table rows: clipper step

       ! end the clipper and the table
       call ImGuiListClipper_End(clipper)
       call ImGuiListClipper_destroy(clipper)
       call igEndTable()
    end if ! begintable

    ! style buttons: show/hide, for the atoms and for the CPs (the
    ! labels padded to the same width)
    call iw_text("Atoms:",alignframe=.true.)
    call igSameLine(0._c_float,-1._c_float)
    ch = showhide_buttons(w%rep%labels%style%shown,"labels","labels",ttshown)
    changed = changed .or. ch
    if (ncp > 0) then
       call iw_text("CPs:  ",alignframe=.true.)
       call igSameLine(0._c_float,-1._c_float)
       ch = showhide_buttons(w%rep%labels%style%cpshown,"cplabels","critical point labels",ttshown)
       changed = changed .or. ch
    end if

    ! process transient highlighs
    if (ihighlight > 0) then
       call sysc(isys)%highlight_atoms(.true.,(/ihighlight/),highlight_type,&
          reshape(ColorHighlightScene,(/4,1/)))
    end if

  contains
    !> Draw row j of the critical point rows of the label table.
    subroutine draw_cp_row(j)
      integer, intent(in) :: j

      integer :: it

      ! id: the cell CP (types 2,3) or symmetry-unique CP index; none
      ! for the per-type rows
      ncol = ncol + 1
      if (igTableSetColumnIndex(ncol)) then
         select case(w%rep%labels%type)
         case (2,3)
            call iw_text(string(sys(isys)%c%ncel+j),alignframe=.true.)
         case (1,4,8)
            call iw_text(string(sys(isys)%c%nneq+j),alignframe=.true.)
         case default
            call iw_text("",alignframe=.true.)
         end select
      end if

      ! CP type in the atom column, and nothing in the Z column
      if (intable /= atlisttype_nmol) then
         ncol = ncol + 1
         it = w%rep%labels%style%cptyp(j)
         if (igTableSetColumnIndex(ncol) .and. it >= 0 .and. it <= 3) &
            call iw_text(cps_abbrev(it))
         ncol = ncol + 1
      end if

      ! shown
      ncol = ncol + 1
      if (igTableSetColumnIndex(ncol)) then
         changed = changed .or. iw_checkbox("##labeltablecpshown" // suffix,w%rep%labels%style%cpshown(j))
         call iw_tooltip("Toggle display of labels for these critical points",ttshown)
      end if

      ! text
      ncol = ncol + 1
      if (igTableSetColumnIndex(ncol)) then
         changed = changed .or. iw_inputtext("##labeltablecptext" // suffix,bufsize=32,&
            textf=w%rep%labels%style%cpstr(j),width=15)
         call iw_tooltip("Text for the critical point labels",ttshown)
      end if

    end subroutine draw_cp_row
  end function draw_editrep_labels

  !> Draw the editrep window, coordination polyhedra class. Returns
  !> true if the representation has changed.
  module function draw_editrep_polyhedra(w,ttshown) result(changed)
    use systems, only: sys, sysc, atlisttype_species, atlisttype_nneq
    use gui_main, only: ColorHighlightScene
    use tools_io, only: string
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_calcheight, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_table_column,&
       iw_helpermark
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: ispc, isys, iz
    character(kind=c_char,len=:), allocatable, target :: str1
    logical :: ch, typechanged
    integer(c_int) :: flags
    integer :: i, j, ncol, ihighlight, highlight_type
    integer :: itype_combo, newtype
    type(ImVec2) :: sz

    ! initialize
    changed = .false.
    ihighlight = 0
    highlight_type = atlisttype_species
    isys = w%isys

    ! make sure the style is initialized and matches the current system
    if (.not.w%rep%poly%style%isinit) then
       call w%rep%poly%style%reset(w%rep)
    elseif (w%rep%poly%style%ntype /= sysc(isys)%attype_number(w%rep%poly%style%type) .or.&
       size(w%rep%poly%style%corner,1) /= sys(isys)%c%nspc) then
       call w%rep%poly%style%reset(w%rep)
    end if

    ! show the corner atoms of coordination polyhedra, with the colors and
    ! radii of this object (the table below, shown only when they are drawn)
    changed = changed .or. iw_checkbox("Show atoms at polyhedra corners##polyshowcorners",&
       w%rep%poly%showcorners)
    call iw_tooltip("Also draw the atoms at the polyhedra corners that the Display leaves out, "//&
       "with the colors and radii in the table below.",ttshown)
    if (w%rep%poly%showcorners) then
       ch = atom_table_widget(isys,w%rep%atoms%style%type,typechanged,ihighlight,highlight_type,&
          rgb=w%rep%atoms%style%rgb,rad=w%rep%atoms%style%rad)
       if (typechanged) call w%rep%atoms%style%reset(w%rep)
       changed = changed .or. ch .or. typechanged
    end if

    ! center type selector
    call iw_text("Centers and Corners",highlight=.true.,alignframe=.true.)
    call iw_helpermark("In the table, each row corresponds to a polyhedron center. &
       &For each center, the columns show the atom name, whether the polyhedron is shown (Show),&
       & the distance range to the corners (Min/Max), and which atomic species are allowed as corners.",sameline=.true.)
    itype_combo = 0
    if (w%rep%poly%style%type == atlisttype_nneq) itype_combo = 1
    if (w%rep%poly%style%type == sysc(isys)%attype_celatom_type()) itype_combo = 2
    call iw_combo_simple("Centers##polycentertype","Species" // c_null_char //&
       "Non-equivalent atoms" // c_null_char // "Cell atoms" // c_null_char,itype_combo)
    call iw_tooltip("How to group the atoms that act as polyhedra centers",ttshown)
    newtype = atlisttype_species
    if (itype_combo == 1) newtype = atlisttype_nneq
    if (itype_combo == 2) newtype = sysc(isys)%attype_celatom_type()
    if (newtype /= w%rep%poly%style%type) then
       w%rep%poly%style%type = newtype
       call w%rep%poly%style%reset(w%rep)
       changed = .true.
    end if

    ! per-center table: each row is a center atom; the species columns
    ! on the right select which species are its corners
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_ScrollX)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    ncol = 5 + sys(isys)%c%nspc
    str1 = "##tablepolycenters_editrep" // c_null_char
    sz%x = 0
    sz%y = iw_calcheight(min(8,w%rep%poly%style%ntype)+1,0,.false.)
    if (igBeginTable(c_loc(str1),ncol,flags,sz,0._c_float)) then
       call iw_table_column("Id",id=0)
       call iw_table_column("Atom",id=1)
       call iw_table_column("Show",id=2)
       call iw_table_column("Min (Å)",id=3)
       call iw_table_column("Max (Å)",id=4)
       do j = 1, sys(isys)%c%nspc
          call iw_table_column(trim(sys(isys)%c%spc(j)%name),id=4+j)
       end do
       call igTableSetupScrollFreeze(1,1)
       call iw_table_headers_row()

       do i = 1, w%rep%poly%style%ntype
          ispc = sysc(isys)%attype_species(w%rep%poly%style%type,i)
          iz = sys(isys)%c%spc(ispc)%z
          if (iz <= 0) cycle
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
          if (igTableSetColumnIndex(0)) then
             call iw_text(string(i),alignframe=.true.)
          end if
          if (igTableSetColumnIndex(1)) &
             call iw_text(sysc(isys)%attype_name(w%rep%poly%style%type,i))
          if (igTableSetColumnIndex(2)) &
             changed = changed .or. iw_checkbox("##polyshown" // string(i),w%rep%poly%style%shown(i))
          if (igTableSetColumnIndex(3)) &
             changed = changed .or. iw_dragfloat_real8("##polydmin" // string(i),&
                x1=w%rep%poly%style%dmin(i),speed=0.01d0,min=0d0,max=20d0,scale=bohrtoa,&
                decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
          if (igTableSetColumnIndex(4)) &
             changed = changed .or. iw_dragfloat_real8("##polydmax" // string(i),&
                x1=w%rep%poly%style%dmax(i),speed=0.01d0,min=0d0,max=20d0,scale=bohrtoa,&
                decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
          do j = 1, sys(isys)%c%nspc
             if (igTableSetColumnIndex(4+j)) &
                changed = changed .or. iw_checkbox("##polycorner" // string(i) // "_" // string(j),&
                   w%rep%poly%style%corner(j,i))
          end do
       end do
       call igEndTable()
    end if

    if (iw_button("Reset##resetpolycenters",danger=.true.)) then
       call w%rep%poly%style%reset(w%rep)
       changed = .true.
    end if
    call iw_tooltip("Reset the centers, corners, and distances to defaults",ttshown)

    ! appearance
    call iw_text("Appearance",highlight=.true.)
    changed = changed .or. iw_dragfloat_real8("Face opacity##polyalpha",x1=w%rep%poly%alpha,&
       speed=0.005d0,min=0d0,max=1d0,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Opacity of the polyhedron faces (0 = transparent, 1 = opaque)",ttshown)

    changed = changed .or. iw_checkbox("Faces use the central atom color##polyusecen",w%rep%poly%usecentercolor)
    call iw_tooltip("Color the faces with the shade of the central (cation) atom",ttshown)
    if (.not.w%rep%poly%usecentercolor) &
       changed = changed .or. iw_coloredit("Face color##polyfacecolor",rgb=w%rep%poly%rgb,sameline=.true.)

    changed = changed .or. iw_dragfloat_real8("Edge radius (Å)##polyedgerad",x1=w%rep%poly%edge_rad,&
       speed=0.002d0,min=0d0,max=1d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Radius of the polyhedron edge cylinders",ttshown)

    changed = changed .or. iw_checkbox("Edges use the central atom color##polyusecenedge",&
       w%rep%poly%usecentercolor_edge)
    call iw_tooltip("Color the edges with the shade of the central (cation) atom",ttshown)
    if (.not.w%rep%poly%usecentercolor_edge) &
       changed = changed .or. iw_coloredit("Edge color##polyedgecolor",rgb=w%rep%poly%edge_rgb,sameline=.true.)

    ! process transient highlighs
    if (ihighlight > 0) then
       call sysc(isys)%highlight_atoms(.true.,(/ihighlight/),highlight_type,&
          reshape(ColorHighlightScene,(/4,1/)))
    end if

  end function draw_editrep_polyhedra

  !> Draw the editrep window, unit cell class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_unitcell(w,ttshown) result(changed)
    use utils, only: iw_text, iw_tooltip, iw_clamp_color3, iw_checkbox,&
       iw_coloredit, iw_dragfloat_real8
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch

    ! initialize
    changed = .false.

    !! styles
    ! each widget is evaluated into ch first: .or. is allowed to short-circuit,
    ! and a widget skipped because changed is already true is a widget not drawn
    call iw_text("Style",highlight=.true.)
    ch = iw_checkbox("Color crystallographic axes",w%rep%uc%coloraxes)
    call iw_tooltip("Represent crystallographic axes with colors (a=red,b=green,c=blue)",ttshown)
    changed = changed .or. ch
    ch = iw_checkbox("Hide axes in vacuum directions",w%rep%uc%vaccutsticks)
    call iw_tooltip("In systems with vacuum direction(s), do not show the unit cell in the vacuum region",ttshown)
    changed = changed .or. ch

    ch = iw_dragfloat_real8("Radius (Å)##outer",x1=w%rep%uc%radius,speed=0.005d0,&
       min=0d0,max=5d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Radii of the unit cell edges",ttshown)
    changed = changed .or. ch

    ch = iw_coloredit("Color",rgb=w%rep%uc%rgb,sameline=.true.)
    call iw_tooltip("Color of the unit cell edges",ttshown)
    if (ch) then
       w%rep%uc%rgb = min(w%rep%uc%rgb,1._c_float)
       w%rep%uc%rgb = max(w%rep%uc%rgb,0._c_float)
       changed = .true.
    end if

    ch = iw_dragfloat_real8("Translate Origin (fractional)##originuc",&
       x3=w%rep%uc%origin,speed=0.001d0,decimal=5)
    call iw_tooltip("Translation vector for the origin of the cell drawn (moves the cell sticks only)",ttshown)
    changed = changed .or. ch

    !! inner divisions
    call iw_text("Inner Divisions",highlight=.true.)
    ch = iw_checkbox("Display inner divisions",w%rep%uc%inner)
    call iw_tooltip("Represent the inner divisions inside a supercell",ttshown)
    changed = changed .or. ch
    if (w%rep%uc%inner) then
       ch = iw_dragfloat_real8("Radius (Å)##inner",x1=w%rep%uc%radiusinner,speed=0.005d0,&
          min=0d0,max=5d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Radii of the inner unit cell edges",ttshown)
       changed = changed .or. ch

       ch = iw_checkbox("Use dashed lines",w%rep%uc%innerstipple)
       call iw_tooltip("Use dashed lines for the inner cell divisions",ttshown)
       changed = changed .or. ch

       if (w%rep%uc%innerstipple) then
          ch = iw_dragfloat_real8("Dash length (Å)",x1=w%rep%uc%innersteplen,speed=0.1d0,&
             min=0.001d0,max=100d0,scale=bohrtoa,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Length of the dashed lines for the inner cell divisions (in Å)",ttshown)
          changed = changed .or. ch
       end if
    end if

  end function draw_editrep_unitcell

  !> Draw the editrep window, legend class. Returns true if the scene
  !> needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_legend(w,ttshown) result(changed)
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8,&
       iw_combo_simple, iw_table_column, iw_table_headers_row, iw_inputtext, iw_calcwidth,&
       iw_calcheight
    use systems, only: sys
    use shapes, only: legend_label_len
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch
    real(c_float) :: rgba(4)
    character(kind=c_char,len=:), allocatable, target :: str1
    integer(c_int) :: flags
    integer :: is, ncol
    type(ImVec2) :: sz

    ! initialize
    changed = .false.

    ! each widget is evaluated into ch first: .or. is allowed to short-circuit,
    ! and a widget skipped because changed is already true is a widget not drawn
    call iw_text("Placement",highlight=.true.)
    call iw_combo_simple("Corner##legendcorner","Top left" // c_null_char // "Top right" // c_null_char //&
       "Bottom left" // c_null_char // "Bottom right" // c_null_char,w%rep%legend%corner,changed=ch)
    call iw_tooltip("Corner of the view where the legend is drawn",ttshown)
    changed = changed .or. ch

    ch = iw_dragfloat_real8("Size##legendsize",x1=w%rep%legend%scale,speed=0.01d0,&
       min=0.2d0,max=5d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Size of the legend (1 = default)",ttshown)
    changed = changed .or. ch

    ! the species rows: whether they are shown, and their text
    call iw_text("Species",highlight=.true.)
    if (w%rep%legend%nspc == sys(w%isys)%c%nspc .and. allocated(w%rep%legend%shown)) then
       flags = ImGuiTableFlags_None
       flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
       flags = ior(flags,ImGuiTableFlags_RowBg)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       str1 = "##tablelegendspecies" // c_null_char
       sz%x = iw_calcwidth(30,3)
       sz%y = iw_calcheight(min(8,w%rep%legend%nspc+1),0,.false.)
       if (igBeginTable(c_loc(str1),3,flags,sz,0._c_float)) then
          ncol = -1
          call iw_table_column("Species",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Show",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Text",icol=ncol,flags=ImGuiTableColumnFlags_WidthStretch)
          call iw_table_headers_row(freezetop=.true.,autofit=.true.)

          do is = 1, w%rep%legend%nspc
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)

             ! species
             if (igTableSetColumnIndex(0_c_int)) &
                call iw_text(trim(sys(w%isys)%c%spc(is)%name),alignframe=.true.)

             ! shown
             if (igTableSetColumnIndex(1_c_int)) then
                ch = iw_checkbox("##legendtableshown" // string(is),w%rep%legend%shown(is))
                call iw_tooltip("Show this species in the legend",ttshown)
                changed = changed .or. ch
             end if

             ! text
             if (igTableSetColumnIndex(2_c_int)) then
                ch = iw_inputtext("##legendtabletext" // string(is),bufsize=legend_label_len,&
                   textf=w%rep%legend%label(is),width=15)
                call iw_tooltip("Text for this species in the legend",ttshown)
                changed = changed .or. ch
             end if
          end do
          call igEndTable()
       end if
    end if

    call iw_text("Colors",highlight=.true.)
    ch = iw_coloredit("Text",rgb=w%rep%legend%textrgb)
    call iw_tooltip("Color of the legend text",ttshown)
    changed = changed .or. ch

    rgba = (/w%rep%legend%bgrgb,w%rep%legend%bgalpha/)
    ch = iw_coloredit("Background",rgba=rgba,sameline=.true.)
    call iw_tooltip("Color and opacity of the background of the legend box",ttshown)
    if (ch) then
       w%rep%legend%bgrgb = rgba(1:3)
       w%rep%legend%bgalpha = rgba(4)
       changed = .true.
    end if

    ch = iw_checkbox("Border",w%rep%legend%border)
    call iw_tooltip("Draw the border of the legend box",ttshown)
    changed = changed .or. ch
    if (w%rep%legend%border) then
       ch = iw_coloredit("Border color",rgb=w%rep%legend%borderrgb,sameline=.true.)
       call iw_tooltip("Color of the border of the legend box",ttshown)
       changed = changed .or. ch
    end if

  end function draw_editrep_legend

  !> Draw the editrep window, cartesian axes class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_axes(w,ttshown) result(changed)
    use windows, only: win
    use representations, only: axplace_scene, axplace_window
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_inputtext,&
       iw_radiobutton, iw_combo_simple
    use systems, only: sys
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch
    integer :: icoord, iview
    real*8 :: zf

    ! initialize
    changed = .false.

    ! the toolbar, which places and edits the axes in the view
    if (w%anchor_view() > 0) &
       call w%editrep_toolbar(ttshown)

    !! axes kind: cartesian (x/y/z) or crystallographic (a/b/c). Only
    !! meaningful for crystals (molecules have no lattice vectors).
    if (.not.sys(w%isys)%c%ismolecule) then
       call iw_text("Type",highlight=.true.)
       icoord = w%rep%axes%kind
       call iw_combo_simple("##axeskind","Cartesian" // c_null_char // &
          "Crystallographic" // c_null_char,icoord,changed=ch)
       call iw_tooltip("Represent the cartesian (x/y/z) or crystallographic (a/b/c) axes",ttshown)
       if (ch) then
          w%rep%axes%kind = icoord
          if (icoord == 0) then
             w%rep%axes%labelstr(1) = "x"
             w%rep%axes%labelstr(2) = "y"
             w%rep%axes%labelstr(3) = "z"
          else
             w%rep%axes%labelstr(1) = "a"
             w%rep%axes%labelstr(2) = "b"
             w%rep%axes%labelstr(3) = "c"
          end if
          changed = .true.
       end if
    end if

    !! position
    call iw_text("Position",highlight=.true.)
    ch = iw_radiobutton("Fixed in Window",int=w%rep%axes%placement,intval=int(axplace_window,c_int))
    if (ch) w%rep%axes%scale_auto = .true. ! re-size the gizmo for the scene
    changed = changed .or. ch
    call iw_tooltip("Anchor the axes at a fixed position in the window",ttshown)
    changed = changed .or. iw_radiobutton("At Position",int=w%rep%axes%placement,&
       intval=int(axplace_scene,c_int),sameline=.true.)
    call iw_tooltip("Place the axes at the cartesian origin",ttshown)
    if (w%rep%axes%placement == axplace_scene) then
       if (sys(w%isys)%c%ismolecule) then
          ! molecules: only cartesian options, referred to the molecular center
          icoord = max(w%rep%axes%coordtype,1) - 1
          call iw_combo_simple("Coordinates##axescoord","Cartesian (Å)" // c_null_char // &
             "Cartesian (bohr)" // c_null_char,icoord,changed=ch)
          if (ch .or. (icoord+1) /= w%rep%axes%coordtype) changed = .true.
          w%rep%axes%coordtype = icoord + 1
       else
          icoord = w%rep%axes%coordtype
          call iw_combo_simple("Coordinates##axescoord","Crystallographic" // c_null_char // &
             "Cartesian (Å)" // c_null_char // "Cartesian (bohr)" // c_null_char,icoord,changed=ch)
          if (ch .or. icoord /= w%rep%axes%coordtype) changed = .true.
          w%rep%axes%coordtype = icoord
       end if
       call iw_tooltip("Coordinate system in which the origin is given",ttshown)
       changed = changed .or. iw_dragfloat_real8("##originaxes",x3=w%rep%axes%origin,speed=0.001d0,decimal=5)
       call iw_tooltip("Coordinates for the origin of the axes",ttshown)
    else
       changed = changed .or. iw_dragfloat_real8("from left##axeswinx",x1=w%rep%axes%winpos(1),speed=0.005d0,&
          min=0d0,max=1d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Horizontal position of the axes, as a fraction of the window width from the left",ttshown)
       changed = changed .or. iw_dragfloat_real8("from bottom##axeswiny",x1=w%rep%axes%winpos(2),speed=0.005d0,&
          min=0d0,max=1d0,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Vertical position of the axes, as a fraction of the window height from the bottom",ttshown)
       if (iw_checkbox("Scale with zoom##axesscalewithzoom",w%rep%axes%scalewithzoom)) then
          changed = .true.
          ! keep the apparent on-screen size unchanged across the toggle
          iview = w%anchor_view()
          if (iview > 0) then
             zf = real(win(iview)%sc%overlay_zoom_factor(),8)
             if (zf > 1d-10) then
                if (w%rep%axes%scalewithzoom) then
                   w%rep%axes%scale = w%rep%axes%scale * zf
                else
                   w%rep%axes%scale = w%rep%axes%scale / zf
                end if
                w%rep%axes%scale_auto = .false.
             end if
          end if
       end if
       call iw_tooltip("Let the gizmo grow and shrink as the scene is zoomed in and out (on), or keep&
          & it at a constant size on the window (off)",ttshown)
    end if

    !! global scale
    call iw_text("Scale",highlight=.true.)
    ch = iw_dragfloat_real8("Scale##axesscale",x1=w%rep%axes%scale,speed=0.01d0,&
       min=0.01d0,max=100d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
    if (ch) w%rep%axes%scale_auto = .false. ! a manual edit disables auto-sizing
    changed = changed .or. ch
    call iw_tooltip("Global scale factor applied to the whole gizmo (arrows and labels)",ttshown)

    !! geometry
    call iw_text("Arrow shaft",highlight=.true.)
    changed = changed .or. iw_dragfloat_real8("Length (Å)##arrowshaftlength",x1=w%rep%axes%length,speed=0.01d0,&
       min=0d0,max=100d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Length of each cartesian axis",ttshown)
    changed = changed .or. iw_dragfloat_real8("Radius (Å)##axesradius",x1=w%rep%axes%radius,speed=0.005d0,&
       min=0d0,max=5d0,scale=bohrtoa,decimal=3,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Radius of the axis shafts",ttshown)

    !! arrowheads
    call iw_text("Arrow head",highlight=.true.)
    changed = changed .or. iw_dragfloat_real8("Length (Å)##arrowheadlength",x1=w%rep%axes%conelength,speed=0.005d0,&
       min=0d0,max=10d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Length of the arrowhead cones (controls how pointy the arrows are)",ttshown)
    changed = changed .or. iw_dragfloat_real8("Radius (Å)##conerad",x1=w%rep%axes%coneradius,speed=0.005d0,&
       min=0d0,max=5d0,scale=bohrtoa,decimal=3,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Base radius of the arrowhead cones (controls how wide the arrows are)",ttshown)

    !! colors
    call iw_text("Colors",highlight=.true.)
    changed = changed .or. iw_coloredit("x",rgb=w%rep%axes%rgb(:,1))
    call iw_tooltip("Color of the x axis",ttshown)
    changed = changed .or. iw_coloredit("y",rgb=w%rep%axes%rgb(:,2),sameline=.true.)
    call iw_tooltip("Color of the y axis",ttshown)
    changed = changed .or. iw_coloredit("z",rgb=w%rep%axes%rgb(:,3),sameline=.true.)
    call iw_tooltip("Color of the z axis",ttshown)

    !! labels
    call iw_text("Labels",highlight=.true.)
    changed = changed .or. iw_checkbox("Show x/y/z labels",w%rep%axes%showlabels)
    call iw_tooltip("Draw the x/y/z labels at the axis tips",ttshown)
    if (w%rep%axes%showlabels) then
       changed = changed .or. iw_dragfloat_real8("Label size",x1=w%rep%axes%labelscale,speed=0.01d0,&
          min=0.1d0,max=5d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Scale of the axis labels",ttshown)
       changed = changed .or. iw_coloredit("Color##axeslabel",rgb=w%rep%axes%labelrgb,sameline=.true.)
       call iw_tooltip("Color of the axis labels",ttshown)
       if (w%rep%axes%placement == axplace_scene) then
          ! for the window-anchored gizmo the label zoom behavior follows the
          ! "Scale with zoom" option above
          changed = changed .or. iw_checkbox("Constant size##axeslabelconstsize",&
             w%rep%axes%labelconstsize,sameline=.true.)
          call iw_tooltip("Labels have constant size (on) or labels scale with the&
             & size of the arrowhead (off)",ttshown)
       end if
       changed = changed .or. iw_inputtext("x##lblx",bufsize=31,textf=w%rep%axes%labelstr(1),width=5)
       changed = changed .or. iw_inputtext("y##lbly",bufsize=31,textf=w%rep%axes%labelstr(2),width=5,sameline=.true.)
       changed = changed .or. iw_inputtext("z##lblz",bufsize=31,textf=w%rep%axes%labelstr(3),width=5,sameline=.true.)
       call iw_tooltip("Text shown for each axis label",ttshown)

       changed = changed .or. iw_dragfloat_real8("Distance to arrow head (Å)",x1=w%rep%axes%labeldistance,&
          speed=0.01d0,scale=bohrtoa,decimal=3)
       call iw_tooltip("Distance from the arrowhead to the label, along the axis (all axes)",ttshown)
       changed = changed .or. iw_dragfloat_real8("x-axis offset##axeslbloffx",x3=w%rep%axes%labeloffset(:,1),&
          speed=0.01d0,scale=bohrtoa,decimal=3)
       call iw_tooltip("Cartesian offset (Å) of the x-axis label from its position on the axis",ttshown)
       changed = changed .or. iw_dragfloat_real8("y-axis offset##axeslbloffy",x3=w%rep%axes%labeloffset(:,2),&
          speed=0.01d0,scale=bohrtoa,decimal=3)
       call iw_tooltip("Cartesian offset (Å) of the y-axis label from its position on the axis",ttshown)
       changed = changed .or. iw_dragfloat_real8("z-axis offset##axeslbloffz",x3=w%rep%axes%labeloffset(:,3),&
          speed=0.01d0,scale=bohrtoa,decimal=3)
       call iw_tooltip("Cartesian offset (Å) of the z-axis label from its position on the axis",ttshown)
    end if

  end function draw_editrep_axes

  !xx! shared table widgets (also used by the display window)

  !> Draw the table of atom groups of system isys, grouped by itype
  !> (atlisttype_*), with the columns the caller asks for: a Show
  !> checkbox column when shown is given (with Show All / Hide All /
  !> Toggle buttons under the table), and Col / Radius columns when rgb
  !> and rad are given. Returns whether anything changed. A change of
  !> grouping is reported in typechanged, and the table is still drawn
  !> over the previous grouping (which the arrays describe); the caller
  !> remakes its arrays for the new one. ihighlight/highlight_type are set to the
  !> row under the mouse (and its grouping) for the scene highlight, and
  !> left alone otherwise.
  module function atom_table_widget(isys,itype,typechanged,ihighlight,highlight_type,shown,rgb,rad,&
     shift) result(changed)
    use systems, only: sys, sysc, atlisttype_species, atlisttype_nneq, atlisttype_ncel_ang,&
       atlisttype_ncel_frac
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_calcheight, iw_checkbox, iw_button, iw_coloredit,&
       iw_highlight_selectable, iw_dragfloat_real8, iw_table_column, iw_intstepper
    use tools_io, only: string, ioj_right
    use param, only: bohrtoa
    use display, only: disp_maxshift
    integer, intent(in) :: isys
    integer, intent(inout) :: itype
    logical, intent(out) :: typechanged
    integer, intent(inout) :: ihighlight
    integer, intent(inout) :: highlight_type
    logical, intent(inout), optional :: shown(:)
    real(c_float), intent(inout), optional :: rgb(:,:)
    real*8, intent(inout), optional :: rad(:)
    integer, intent(inout), optional :: shift(:,:)
    logical :: changed

    logical :: domol, docoord, doshown, dostyle, doshift, ch
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2, suffix
    real*8 :: x0(3)
    type(ImVec2) :: sz0
    integer :: ispc, i, j, iz, ncol, icol, ntype, itab
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f

    logical, save :: ttshown = .false. ! tooltip flag

    integer, parameter :: atlisttype_allowed_crys(3) = (/atlisttype_species,atlisttype_nneq,&
       atlisttype_ncel_frac/)
    integer, parameter :: atlisttype_allowed_mol(2) = (/atlisttype_species,atlisttype_ncel_ang/)
    integer, allocatable :: atlisttype_allowed(:)

    ! initialize
    changed = .false.
    typechanged = .false.
    doshown = present(shown)
    dostyle = present(rgb) .and. present(rad)
    doshift = present(shift)
    if (doshift) doshift = .not.sys(isys)%c%ismolecule
    if (doshown) then
       call iw_text("Shown Atoms",highlight=.true.)
    elseif (dostyle) then
       call iw_text("Atom Style",highlight=.true.)
    end if
    if (sys(isys)%c%ismolecule) then
       atlisttype_allowed = atlisttype_allowed_mol
    else
       atlisttype_allowed = atlisttype_allowed_crys
    end if

    ! the grouping; a change is the caller's to apply, so this frame the
    ! table is still drawn over the grouping the arrays describe
    itab = itype
    typechanged = sysc(isys)%attype_combo_simple("Atom types##atomtypeselection",itype,&
       atlisttype_allowed,units=.false.)
    call iw_tooltip("Group atoms by these categories",ttshown)
    changed = changed .or. typechanged

    ! the arrays must describe this grouping
    ntype = sysc(isys)%attype_number(itab)
    if (doshown) then
       if (size(shown,1) /= ntype) return
    end if
    if (dostyle) then
       if (size(rgb,2) /= ntype .or. size(rad,1) /= ntype) return
    end if
    if (doshift) then
       if (size(shift,2) /= ntype) return
    end if

    ! whether to do the molecule column and the coordinates
    domol = (itab == atlisttype_ncel_ang)
    docoord = (itab == atlisttype_nneq .or. itab == atlisttype_ncel_ang .or. itab == atlisttype_ncel_frac)
    ncol = 3
    if (doshown) ncol = ncol + 1 ! show
    if (dostyle) ncol = ncol + 2 ! col, radius
    if (doshift) ncol = ncol + 1 ! lattice-vector shift
    if (domol) ncol = ncol + 1 ! mol
    if (docoord) ncol = ncol + 1 ! coordinates

    ! atom table
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tableatomstyles" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(5,ntype)+1,0,.false.)
    if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
       icol = -1

       ! header setup
       call iw_table_column("Id",icol=icol)
       call iw_table_column("Atom",icol=icol)
       call iw_table_column("Z ",icol=icol)
       if (doshown) call iw_table_column("Show",icol=icol)
       if (dostyle) then
          call iw_table_column("Col",icol=icol)
          call iw_table_column("Radius",icol=icol)
       end if
       if (doshift) call iw_table_column("Shift",icol=icol)
       if (domol) call iw_table_column("Mol",icol=icol)
       if (docoord) then
          if (itab == atlisttype_ncel_ang) then
             str2 = "Coordinates (Å)"
          else
             str2 = "Coordinates (fractional)"
          end if
          call iw_table_column(str2,icol=icol,flags=ImGuiTableColumnFlags_WidthStretch)
       end if

       ! draw the header
       call iw_table_headers_row(freezetop=.true.,autofit=.true.)

       ! start the clipper
       clipper = ImGuiListClipper_ImGuiListClipper()
       call ImGuiListClipper_Begin(clipper,ntype,-1._c_float)

       ! draw the rows
       do while(ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
             suffix = "_" // string(i)
             icol = -1

             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             ispc = sysc(isys)%attype_species(itab,i)
             iz = sys(isys)%c%spc(ispc)%z

             ! id
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) then
                call iw_text(string(i),alignframe=.true.)

                ! the highlight selectable
                if (iw_highlight_selectable("##selectableatomtable" // suffix)) then
                   ihighlight = i
                   highlight_type = itab
                end if
             end if

             ! name
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) call iw_text(sysc(isys)%attype_name(itab,i))

             ! Z
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) call iw_text(string(iz))

             ! shown
             if (doshown) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   changed = changed .or. iw_checkbox("##tableshown" // suffix ,shown(i))
                   call iw_tooltip("Toggle display of the atom/bond/label associated to this atom",ttshown)
                end if
             end if

             ! color and radius
             if (dostyle) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   ch = iw_coloredit("##tablecolor" // suffix,rgb=rgb(:,i))
                   call iw_tooltip("Atom color",ttshown)
                   changed = changed .or. ch
                end if

                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   ch = iw_dragfloat_real8("##tableradius" // string(i),x1=rad(i),speed=0.01d0,&
                      min=0d0,max=5d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
                   call iw_tooltip("Radius of the sphere representing the atom",ttshown)
                   if (ch) then
                      rad(i) = max(rad(i),0d0)
                      changed = .true.
                   end if
                end if
             end if

             ! lattice-vector shift
             if (doshift) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   do j = 1, 3
                      ch = iw_intstepper("tableshift" // string(j) // suffix,shift(j,i),&
                         minval=-disp_maxshift,maxval=disp_maxshift,ndigit=3,notlive=.true.,&
                         sameline=(j > 1),tooltip="Draw these atoms this lattice vector away &
                         &(fractional, integers)")
                      changed = changed .or. ch
                   end do
                end if
             end if

             ! molecule
             if (domol) then
                icol = icol + 1
                ! i is a complete list index in this case
                if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%idatcelmol(1,i)))
             end if

             ! coordinates
             if (docoord) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   s = ""
                   if (itab > 0) then
                      x0 = sysc(isys)%attype_coordinates(itab,i)
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
       call ImGuiListClipper_destroy(clipper)
       call igEndTable()
    end if

    ! show/hide buttons
    if (doshown) changed = changed .or. showhide_buttons(shown,"atoms","atoms/bonds/labels",ttshown)

  end function atom_table_widget

  !> Draw the table of molecules of system isys (nothing if it has one
  !> molecule or fewer), with the columns the caller asks for: a Show
  !> checkbox column when shown is given (with Show All / Hide All /
  !> Toggle buttons under the table), and Tint / Scale columns when tint
  !> and scale are given. Returns whether anything changed.
  !> ihighlight/highlight_type are set to the row under the mouse (and
  !> the molecule grouping), and left alone otherwise.
  module function mol_table_widget(isys,ihighlight,highlight_type,shown,tint,scale) result(changed)
    use systems, only: sys, atlisttype_nmol
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_calcheight, iw_checkbox, iw_button, iw_coloredit,&
       iw_highlight_selectable, iw_dragfloat_real8, iw_table_column
    use global, only: iunit_ang, dunit0
    use tools_io, only: string, ioj_right
    integer, intent(in) :: isys
    integer, intent(inout) :: ihighlight
    integer, intent(inout) :: highlight_type
    logical, intent(inout), optional :: shown(:)
    real(c_float), intent(inout), optional :: tint(:,:)
    real*8, intent(inout), optional :: scale(:)
    logical :: changed

    logical :: doshown, dostyle, ch
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2
    real*8 :: x0(3)
    type(ImVec2) :: sz0
    integer :: i, ncol, icol, nmol
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f

    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    changed = .false.
    doshown = present(shown)
    dostyle = present(tint) .and. present(scale)
    nmol = sys(isys)%c%nmol
    if (nmol <= 1) return
    if (doshown) then
       if (size(shown,1) /= nmol) return
    end if
    if (dostyle) then
       if (size(tint,2) /= nmol .or. size(scale,1) /= nmol) return
    end if
    if (doshown) then
       call iw_text("Shown Molecules",highlight=.true.)
    elseif (dostyle) then
       call iw_text("Molecule Style",highlight=.true.)
    end if

    ncol = 2
    if (doshown) ncol = ncol + 1 ! show
    if (dostyle) ncol = ncol + 2 ! tint, scale
    ncol = ncol + 1 ! center of mass

    ! molecule table
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tablemolstyles" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(5,nmol)+1,0,.false.)
    if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
       icol = -1

       ! header setup
       call iw_table_column("Id",icol=icol)
       call iw_table_column("Nat",icol=icol)
       if (doshown) call iw_table_column("Show",icol=icol)
       if (dostyle) then
          call iw_table_column("Tint",icol=icol)
          call iw_table_column("Scale",icol=icol)
       end if
       if (sys(isys)%c%ismolecule) then
          str2 = "Center of mass (Å)"
       else
          str2 = "Center of mass (fractional)"
       end if
       call iw_table_column(str2,icol=icol,flags=ImGuiTableColumnFlags_WidthStretch)

       ! draw the header
       call iw_table_headers_row(freezetop=.true.,autofit=.true.)

       ! start the clipper
       clipper = ImGuiListClipper_ImGuiListClipper()
       call ImGuiListClipper_Begin(clipper,nmol,-1._c_float)

       ! draw the rows
       do while(ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
             icol = -1
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)

             ! id
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) then
                call iw_text(string(i),alignframe=.true.)

                ! the highlight selectable
                if (iw_highlight_selectable("##selectablemoltable_" // string(i))) then
                   ihighlight = i
                   highlight_type = atlisttype_nmol
                end if
             end if

             ! nat
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%mol(i)%nat))

             ! shown
             if (doshown) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   changed = changed .or. iw_checkbox("##tablemolshown" // string(i),shown(i))
                   call iw_tooltip("Toggle display of all atoms in this molecule",ttshown)
                end if
             end if

             ! tint and scale
             if (dostyle) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   ch = iw_coloredit("##tablemolcolor" // string(i),rgb=tint(:,i))
                   call iw_tooltip("Molecule color tint",ttshown)
                   changed = changed .or. ch
                end if

                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   changed = changed .or. iw_dragfloat_real8("##tablemolradius" // string(i),&
                      x1=scale(i),speed=0.005d0,min=0d0,max=5d0,decimal=3,&
                      flags=ImGuiSliderFlags_AlwaysClamp)
                   call iw_tooltip("Scale factor for the atomic radii in this molecule",ttshown)
                end if
             end if

             ! center of mass
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) then
                x0 = sys(isys)%c%mol(i)%cmass(.false.)
                if (sys(isys)%c%ismolecule) then
                   x0 = (x0+sys(isys)%c%molx0) * dunit0(iunit_ang)
                else
                   x0 = sys(isys)%c%c2x(x0)
                endif
                s = string(x0(1),'f',8,4,ioj_right) //" "// string(x0(2),'f',8,4,ioj_right) //" "//&
                   string(x0(3),'f',8,4,ioj_right)
                call iw_text(s)
             end if
          end do ! clipper indices
       end do ! clipper step

       ! end the clipper and the table
       call ImGuiListClipper_End(clipper)
       call ImGuiListClipper_destroy(clipper)
       call igEndTable()
    end if

    ! show/hide buttons
    if (doshown) changed = changed .or. showhide_buttons(shown,"molecules","molecules",ttshown)

  end function mol_table_widget

  !xx! private procedures

  !> The Show All / Hide All / Toggle Show/Hide button row over the mask
  !> shown; idsuffix makes the button IDs unique and what names the
  !> items in the tooltips. Returns whether the mask changed.
  function showhide_buttons(shown,idsuffix,what,ttshown) result(changed)
    use utils, only: iw_button, iw_tooltip
    logical, intent(inout) :: shown(:)
    character(len=*), intent(in) :: idsuffix
    character(len=*), intent(in) :: what
    logical, intent(inout) :: ttshown
    logical :: changed

    changed = .false.
    if (iw_button("Show All##showall" // idsuffix)) then
       shown = .true.
       changed = .true.
    end if
    call iw_tooltip("Show all " // what // " in the system",ttshown)
    if (iw_button("Hide All##hideall" // idsuffix,sameline=.true.)) then
       shown = .false.
       changed = .true.
    end if
    call iw_tooltip("Hide all " // what // " in the system",ttshown)
    if (iw_button("Toggle Show/Hide##toggleall" // idsuffix,sameline=.true.)) then
       shown = .not.shown
       changed = .true.
    end if
    call iw_tooltip("Toggle the show/hide status for all " // what,ttshown)

  end function showhide_buttons

  !> The Apply to Type / Apply to All button pair under the options of
  !> the selected item of a list: the first copies the style of the
  !> selected item to the other items of the same type (typename, plural)
  !> and the second to every other item (items, plural). idsuffix makes
  !> the button IDs unique. Returns 0 (no click), 1 (apply to type), or
  !> 2 (apply to all).
  function apply_style_buttons(idsuffix,typename,items,ttshown) result(imode)
    use utils, only: iw_button, iw_tooltip
    character(len=*), intent(in) :: idsuffix
    character(len=*), intent(in) :: typename
    character(len=*), intent(in) :: items
    logical, intent(inout) :: ttshown
    integer :: imode

    imode = 0
    if (iw_button("Apply to Type##applytype" // idsuffix,danger=.true.)) imode = 1
    call iw_tooltip("Apply these style options to all the " // typename,ttshown)
    if (iw_button("Apply to All##applyall" // idsuffix,danger=.true.,sameline=.true.)) imode = 2
    call iw_tooltip("Apply these style options to all the " // items // ". Options that do &
       &not apply to an item's type are left as they are.",ttshown)

  end function apply_style_buttons

  !> Draw the editrep (Object) window, symmetry-elements class. Returns true if
  !> the scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_symelem(w,ttshown) result(changed)
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_combo_simple,&
       iw_button, iw_calcheight, iw_table_column
    use systems, only: sys
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch, hasdir
    integer :: i, icoord, ntype, ncol, icolop, iwrite
    ! the table, for table_write_* (the measurement tables are 2-4)
    integer, parameter :: itable_symelem = 1
    integer(c_int) :: flags
    type(ImVec2) :: sz0
    character(kind=c_char,len=:), allocatable, target :: str1

    changed = .false.

    ! the symmetry-element style (element types + visibility) was refreshed by
    ! the caller if the geometry changed since the last reset
    ntype = w%rep%symelem%style%se%ntype

    ! the direction column has the axis/plane-normal indices, only for crystals
    hasdir = .not.sys(w%isys)%c%ismolecule

    !! origin (molecules only; in a crystal the element positions are calculated)
    if (sys(w%isys)%c%ismolecule) then
       call iw_text("Origin",highlight=.true.)
       ! cartesian options only (coordtype 1/2), referred to the molecular
       ! center; the combo is 0-based, hence the +/-1 offset
       icoord = max(w%rep%symelem%coordtype,1) - 1
       call iw_combo_simple("Coordinates##symelemcoord","Cartesian (Å)" // c_null_char // &
          "Cartesian (bohr)" // c_null_char,icoord,changed=ch)
       if (ch .or. (icoord+1) /= w%rep%symelem%coordtype) changed = .true.
       w%rep%symelem%coordtype = icoord + 1
       call iw_tooltip("Coordinate system in which the origin is given",ttshown)
       changed = changed .or. iw_dragfloat_real8("##originsymelem",x3=w%rep%symelem%origin,speed=0.001d0,decimal=5)
       call iw_tooltip("Point all the symmetry elements pass through",ttshown)
    end if

    !! color
    call iw_text("Color",highlight=.true.)
    changed = changed .or. iw_checkbox("Custom color##symelemcustomrgb",w%rep%symelem%usecustomrgb)
    call iw_tooltip("Color all elements with a single custom color. If off, mirror planes and glide &
       &planes have colors of their own and rotation axes are colored by rotation order. A glide &
       &plane and a screw axis also carry a dashed arrow along the translation that goes with it.",&
       ttshown)
    if (w%rep%symelem%usecustomrgb) &
       changed = changed .or. iw_coloredit("##symelemrgb",rgb=w%rep%symelem%rgb,sameline=.true.)

    !! elements
    call iw_text("Elements",highlight=.true.)
    call iw_tooltip("The symmetry elements of this system. Operations lists the symmetry operations &
       &that generate each element, numbered as in the symmetry table of the View/Edit Geometry &
       &window; one operation can generate elements with different symbols.",ttshown)
    if (w%rep%symelem%style%isinit) then
       ! all / none / toggle
       if (iw_button("All##symelemall")) then
          w%rep%symelem%style%shown = .true.
          changed = .true.
       end if
       call iw_tooltip("Show all symmetry elements",ttshown)
       if (iw_button("None##symelemnone",sameline=.true.)) then
          w%rep%symelem%style%shown = .false.
          changed = .true.
       end if
       call iw_tooltip("Hide all symmetry elements",ttshown)
       if (iw_button("Toggle##symelemtoggle",sameline=.true.)) then
          w%rep%symelem%style%shown = .not.w%rep%symelem%style%shown
          changed = .true.
       end if
       call iw_tooltip("Toggle the symmetry element selection",ttshown)

       ! per-element-type table
       flags = ImGuiTableFlags_None
       flags = ior(flags,ImGuiTableFlags_RowBg)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       icolop = merge(3,2,hasdir)
       ncol = icolop + 1
       str1 = "##symelemtable" // c_null_char
       sz0%x = 0
       sz0%y = iw_calcheight(min(ntype,10)+1,0,.false.)
       if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
          call iw_table_column("Show",id=0)
          call iw_table_column("Symbol",id=1)
          if (hasdir) call iw_table_column("Direction",id=2)
          call iw_table_column("Operations",id=icolop)
          call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)
          call w%table_write_begin(itable_symelem,iwrite)

          do i = 1, ntype
             call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
             if (igTableSetColumnIndex(0)) then
                if (iw_checkbox("##symelemshow" // string(i),w%rep%symelem%style%shown(i))) changed = .true.
             end if
             if (igTableSetColumnIndex(1)) call iw_text(trim(w%rep%symelem%style%se%label(i)))
             if (hasdir) then
                if (igTableSetColumnIndex(2)) call iw_text(trim(w%rep%symelem%style%se%dirlabel(i)))
             end if
             if (igTableSetColumnIndex(icolop)) &
                call iw_text(oplist_string(w%rep%symelem%style%se,i))
          end do
          call w%table_write_end(itable_symelem,"Symmetry elements")
          call igEndTable()
       end if
    end if

  end function draw_editrep_symelem

  !> Blank-separated list of the symmetry operations that generate element
  !> type it of the element list se, for the Operations column.
  function oplist_string(se,it) result(str)
    use crystalmod, only: symelem_list
    use tools_io, only: string
    type(symelem_list), intent(in) :: se
    integer, intent(in) :: it
    character(len=:), allocatable :: str

    integer :: i

    str = ""
    do i = 1, se%nop(it)
       if (i == 1) then
          str = string(se%iop(i,it))
       else
          str = str // " " // string(se%iop(i,it))
       end if
    end do

  end function oplist_string

  !> Draw the editrep window, text annotations. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_text(w,ttshown) result(changed)
    use representations, only: textpos_screen, textpos_point, textpos_atom, textpos_bond,&
       text_delete, text_copy_style
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_combo_simple,&
       iw_button, iw_calcheight, iw_inputtext, iw_close_button, iw_highlight_selectable, iw_radiobutton,&
       iw_table_column
    use systems, only: sys, sysc
    use tools_io, only: string
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch, ok, okp, ldum, focus
    integer :: i, iview, isel, idel, ipl, iplpick, imode
    real*8 :: xc(3)
    integer(c_int) :: flags
    type(ImVec2) :: sz0
    character(kind=c_char,len=:), allocatable, target :: str1

    ! names of the placements (textpos_*), for the table and the Apply to Type tooltip
    character(len=11), parameter :: placename(0:3) = (/"on-screen  ","3D position",&
       "atom       ","bond       "/)

    ! initialize
    changed = .false.
    iview = w%anchor_view()

    ! handle a pending atom pick commanded to the parent view
    if (w%editrep_pick_item > 0) then
       ok = (w%editrep_pick_item <= w%rep%text%ntext)
       if (ok) ok = (win(iview)%vmdata%owner == w%id) .and.&
          .not.w%editrep_pick%is_stale(sysc(w%isys)%timelastchange_geometry)
       if (.not.ok) then
          ! the item was deleted, another window took over the pick, or the
          ! geometry changed (stale cell-atom ids): cancel
          call win(iview)%viewmode_release_forced(w%id)
          w%editrep_pick_item = 0
       elseif (win(iview)%viewmode >= 0) then
          ! the pick finished
          i = w%editrep_pick_item
          iplpick = w%rep%text%t(i)%placement
          ! the delivery has to be the kind of anchor the item still wants:
          ! the placement can be changed while the pick is out
          if (win(iview)%vmdata%bidx(1) > 0 .and. iplpick == textpos_bond) then
             ! a bond: its two atom images anchor the text
             w%rep%text%t(i)%idx1 = win(iview)%vmdata%bidx(1:4)
             w%rep%text%t(i)%idx2 = win(iview)%vmdata%bidx(5:8)
             changed = .true.
          elseif (win(iview)%vmdata%idx(1) > 0 .and. iplpick == textpos_atom) then
             ! an atom
             w%rep%text%t(i)%idx1 = win(iview)%vmdata%idx(1:4)
             changed = .true.
          elseif (win(iview)%vmdata%flag == 1 .and. win(iview)%vmdata%bidx(1) == 0 .and.&
             iplpick == textpos_point) then
             ! a 3D position: the clicked atom, or the empty-space click
             ! unprojected onto the plane of the scene center. A bond delivery
             ! (the placement was changed while a bond pick was out) carries no
             ! click position, so it is left to be discarded below
             call view_pick_point(iview,w%isys,xc,okp)
             if (okp) then
                if (sys(w%isys)%c%ismolecule) then
                   w%rep%text%t(i)%pos = (xc + sys(w%isys)%c%molx0) * bohrtoa
                else
                   w%rep%text%t(i)%pos = sys(w%isys)%c%c2x(xc)
                end if
                changed = .true.
             end if
          elseif (win(iview)%vmdata%flag == 1 .and. win(iview)%vmdata%bidx(1) == 0 .and.&
             iplpick == textpos_screen) then
             ! an on-screen position: the click as a fraction of the view window
             call view_texpos_to_winfrac(iview,win(iview)%vmdata%xpos,w%rep%text%t(i)%winpos)
             changed = .true.
          end if
          ! the pick is over either way; nothing delivered means it was cancelled
          w%editrep_pick_item = 0
          win(iview)%vmdata%idx = 0
          win(iview)%vmdata%bidx = 0
          win(iview)%vmdata%flag = 0
       end if
    end if

    ! the toolbar, which edits the texts in the view: select, remove,
    ! and one tool per placement
    if (iview > 0) &
       call w%editrep_toolbar(ttshown)

    ! table of text items
    call iw_text("Text objects",highlight=.true.)
    idel = 0
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    str1 = "##texttable" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(w%rep%text%ntext,5)+1,0,.false.)
    if (igBeginTable(c_loc(str1),4,flags,sz0,0._c_float)) then
       call iw_table_column("",id=0,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Show",id=1,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Text",id=2,flags=ImGuiTableColumnFlags_WidthStretch)
       call iw_table_column("Placement",id=3,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_headers_row(freezetop=.true.)

       do i = 1, w%rep%text%ntext
          call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
          ! delete button
          if (igTableSetColumnIndex(0)) then
             call igAlignTextToFramePadding()
             if (iw_close_button("##textdel" // string(i))) idel = i
             call iw_tooltip("Remove this text",ttshown)
          end if
          ! shown checkbox
          if (igTableSetColumnIndex(1)) then
             if (iw_checkbox("##textshow" // string(i),w%rep%text%t(i)%shown)) changed = .true.
             call iw_tooltip("Toggle show/hide this text",ttshown)
          end if
          ! the text (first line only; edited in the box below the table)
          if (igTableSetColumnIndex(2)) then
             call iw_text(first_line(w%rep%text%t(i)%str),alignframe=.true.)
          end if
          ! placement summary; a row-spanning selectable picks the edited item,
          ! and the row being edited is shown highlighted
          if (igTableSetColumnIndex(3)) then
             call iw_text(trim(placename(w%rep%text%t(i)%placement)),alignframe=.true.)
             ldum = iw_highlight_selectable("##textsel" // string(i),clicked=ch,&
                selected=(i == w%rep%text%isel))
             if (ch) w%rep%text%isel = i
          end if
       end do
       call igEndTable()
    end if

    ! process a deletion, dropping a mouse drag on the texts in the view
    if (idel > 0) then
       call text_delete(w%rep%text,idel)
       if (iview > 0) then
          if (win(iview)%annot%irep == w%irep .and. win(iview)%oe%op == objop_move) &
             win(iview)%oe%op = objop_none
       end if
       if (w%editrep_pick_item == idel) then
          ! cancel a pick pending on the deleted item
          call win(iview)%viewmode_release_forced(w%id)
          w%editrep_pick_item = 0
       elseif (w%editrep_pick_item > idel) then
          w%editrep_pick_item = w%editrep_pick_item - 1
       end if
       changed = .true.
    end if

    !! section for changing the selected text; a text placed or
    !! double-clicked in the view asks for the keyboard (an old request,
    !! made while this window was not drawn, is dropped)
    focus = (w%editrep_focustext >= 0 .and. igGetFrameCount() - w%editrep_focustext <= 2)
    w%editrep_focustext = -1
    if (w%rep%text%isel >= 1 .and. w%rep%text%isel <= w%rep%text%ntext) then
       isel = w%rep%text%isel

       ! the text, edited in a multiline box
       call iw_text("Text",highlight=.true.)
       if (focus) call igSetWindowFocus_Nil()
       if (iw_inputtext("##textstredit",bufsize=4095,texta=w%rep%text%t(isel)%str,nlines=3,&
          grabfocus=focus,selectall=.true.)) changed = .true.
       call iw_tooltip("Text of the selected annotation (multiple lines allowed)",ttshown)

       ! placement
       call iw_text("Placement",highlight=.true.,alignframe=.true.)
       ipl = w%rep%text%t(isel)%placement
       call iw_combo_simple("##textplacement","On-screen" // c_null_char // "3D position" // c_null_char //&
          "Atom" // c_null_char // "Bond" // c_null_char,ipl,changed=ch,sameline=.true.)
       changed = changed .or. ch
       w%rep%text%t(isel)%placement = ipl
       call iw_tooltip("Place the text: anchored to the view (on-screen), at a &
          &3D position in the system, or tied to an atom or a bond",ttshown)

       if (ipl == textpos_screen) then
          ch = iw_dragfloat_real8("Position##textwinpos",x2=w%rep%text%t(isel)%winpos,&
             speed=0.005d0,min=0d0,max=1d0,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Position of the text in the viewport, as fractions of the window &
             &size from the left and bottom borders",ttshown)
          changed = changed .or. ch
          call text_pick_button("textpickwinpos",vm_pick_atom,"Pick the on-screen position for the text",&
             "Click, then pick the position in the view window",acceptempty=.true.,sameline=.true.)
          ch = iw_radiobutton("In front##textfront",bool=w%rep%text%t(isel)%infront,&
             boolval=.true.)
          call iw_tooltip("Draw the text in front of the scene",ttshown)
          changed = changed .or. ch
          ch = iw_radiobutton("Behind##textbehind",bool=w%rep%text%t(isel)%infront,&
             boolval=.false.,sameline=.true.)
          call iw_tooltip("Draw the text behind the scene",ttshown)
          changed = changed .or. ch
       elseif (ipl == textpos_point) then
          ch = iw_dragfloat_real8("Position##textpos",x3=w%rep%text%t(isel)%pos,&
             speed=0.001d0,decimal=5)
          if (sys(w%isys)%c%ismolecule) then
             call iw_tooltip("Position of the text (Cartesian coordinates, Å)",ttshown)
          else
             call iw_tooltip("Position of the text (fractional coordinates)",ttshown)
          end if
          changed = changed .or. ch
          call text_pick_button("textpickpos",vm_pick_atom,"Pick the text position",&
             "Click, then pick the position in the view window: an atom to use its &
             &position, or empty space for the clicked point on the plane through the &
             &scene center",acceptempty=.true.,sameline=.true.)
       elseif (ipl == textpos_atom) then
          call text_pick_button("textpickatom",vm_pick_atom,"Pick the atom to anchor the text to",&
             "Click, then pick the anchor atom in the view window")
          call anchor_atom_button(w%rep%text%t(isel)%idx1,"##textanchor1","Anchor atom of the text")
       else ! textpos_bond
          call text_pick_button("textpickbond",vm_pick_bond,"Pick the bond to anchor the text to",&
             "Click, then pick the bond in the view window")
          call anchor_atom_button(w%rep%text%t(isel)%idx1,"##textanchor1",&
             "First atom of the anchor bond")
          call anchor_atom_button(w%rep%text%t(isel)%idx2,"##textanchor2",&
             "Second atom of the anchor bond")
       end if

       ! 3D placements: on-screen offset from the anchor and depth toggle
       if (ipl /= textpos_screen) then
          ch = iw_dragfloat_real8("Offset (Å)##textoffset",x2=w%rep%text%t(isel)%offset,&
             speed=0.01d0,decimal=2)
          call iw_tooltip("Offset of the text from its anchor, in on-screen (in-plane) coordinates",ttshown)
          changed = changed .or. ch
          ch = iw_checkbox("Depth##textdepth",w%rep%text%t(isel)%depth)
          call iw_tooltip("The text is hidden by objects in front of it. If unchecked, the &
             &text is drawn on top of all objects.",ttshown)
          changed = changed .or. ch
       end if

       ! style
       call iw_text("Style",highlight=.true.)
       ch = iw_dragfloat_real8("Size##textscale",x1=w%rep%text%t(isel)%scale,&
          speed=0.01d0,min=0.05d0,max=10d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Size of the text",ttshown)
       changed = changed .or. ch
       ch = iw_coloredit("Color##textcolor",rgb=w%rep%text%t(isel)%rgb,sameline=.true.)
       call iw_tooltip("Color of the text",ttshown)
       changed = changed .or. ch
       ch = iw_checkbox("Scale with zoom##textzoom",w%rep%text%t(isel)%scalewithzoom)
       call iw_tooltip("Whether the text size scales when zooming in and out, or stays &
          &at a constant on-screen size",ttshown)
       changed = changed .or. ch

       ! apply these style options to the texts with this placement, or to all
       imode = apply_style_buttons("text","texts with " // trim(placename(ipl)) // " placement",&
          "texts",ttshown)
       if (imode > 0) then
          do i = 1, w%rep%text%ntext
             if (i == isel .or. (imode == 1 .and. w%rep%text%t(i)%placement /= ipl)) cycle
             call text_copy_style(w%rep%text%t(i),w%rep%text%t(isel))
          end do
          changed = .true.
       end if
    end if

  contains
    !> Draw the Pick button for the placement of the selected text: arms a
    !> pick of mode mode (vm_pick_atom or vm_pick_bond) in the anchor view
    !> that fills the anchor or the position of the item. id is the ImGui
    !> id suffix, msg the prompt shown in the view bar, and tt the button
    !> tooltip; acceptempty makes an empty-space click deliver a position
    !> instead of cancelling. The result is received at the top of this
    !> routine. One pick can be pending at a time.
    subroutine text_pick_button(id,mode,msg,tt,acceptempty,sameline)
      character(len=*), intent(in) :: id, msg, tt
      integer, intent(in) :: mode
      logical, intent(in), optional :: acceptempty, sameline

      logical :: aempty, sline

      aempty = .false.
      if (present(acceptempty)) aempty = acceptempty
      sline = .false.
      if (present(sameline)) sline = sameline

      if (iw_button("Pick##" // id,sameline=sline,disabled=(w%editrep_pick_item > 0))) then
         w%editrep_pick_item = isel
         call w%editrep_pick%arm()
         call win(iview)%viewmode_set_forced(mode,msg,w%id,acceptempty=aempty)
      end if
      call iw_tooltip(tt,ttshown)

    end subroutine text_pick_button

    !> First line of a text (with a continuation mark if there are more lines).
    function first_line(str) result(s)
      character(len=*), intent(in) :: str
      character(len=:), allocatable :: s
      integer :: idx
      idx = index(str,new_line('a'))
      if (idx > 0) then
         s = str(1:idx-1) // " (...)"
      else
         s = str
      end if
    end function first_line

    !> Inert colored button naming the anchor atom idx, in the atom color
    !> and the species+cell-index form the rest of the GUI uses ("?" when
    !> the anchor is unset or no longer names a cell atom). idn is the
    !> ImGui id suffix, the anchors of one item being told apart by it,
    !> and tt the button tooltip.
    subroutine anchor_atom_button(idx,idn,tt)
      integer(c_int), intent(in) :: idx(4)
      character(len=*), intent(in) :: idn, tt

      logical :: ldum2

      ldum2 = draw_anchor_button(iview,w%isys,idx,idn,sameline=.true.,inert=.true.)
      call iw_tooltip(tt,ttshown)

    end subroutine anchor_atom_button

  end function draw_editrep_text

  !> Draw the editrep (Object) window, measurements class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_measure(w,ttshown) result(changed)
    use representations, only: measurement_item, measure_delete
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_button,&
       iw_calcheight, iw_close_button, iw_intstepper, iw_highlight_selectable,&
       iw_helpermark, iw_table_column, iw_begintabitem
    use keybindings, only: get_bind_keyname, BIND_NAV_MEASURE, BIND_NAV_MEASURE_TOGGLE
    use systems, only: sys, sysc
    use tools_io, only: string
    use param, only: pi, bohrtoa
    use interfaces_glfw, only: glfwGetTime
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer(c_int) :: flags
    integer :: idel, k, iview, i, itab
    logical :: ok
    type(ImVec2) :: sz0
    character(kind=c_char,len=:), allocatable, target :: str1

    changed = .false.
    idel = 0
    iview = w%anchor_view()

    ! keep the selected item in range
    if (w%rep%measure%isel > w%rep%measure%nitem) w%rep%measure%isel = w%rep%measure%nitem
    if (w%rep%measure%isel < 0) w%rep%measure%isel = 0

    ! handle a pending atom pick commanded to the parent view: when the user
    ! clicks an atom button in a table, the view enters pick mode and the picked
    ! atom replaces that measurement atom (mirrors the text-anchor pick above)
    if (w%editrep_pick_item > 0) then
       ok = (w%editrep_pick_item <= w%rep%measure%nitem)
       if (ok) ok = (win(iview)%vmdata%owner == w%id) .and.&
          .not.w%editrep_pick%is_stale(sysc(w%isys)%timelastchange_geometry)
       if (.not.ok) then
          ! the item was deleted, another window took over the pick, or the
          ! geometry changed (stale cell-atom ids): cancel
          call win(iview)%viewmode_release_forced(w%id)
          w%editrep_pick_item = 0
          w%editrep_pick_slot = 0
       elseif (win(iview)%viewmode >= 0) then
          ! the pick finished: replace the atom if one was picked (else cancelled)
          i = w%editrep_pick_item
          k = w%editrep_pick_slot
          if (win(iview)%vmdata%bidx(1) > 0) then
             ! a bond: its two atom images become the ends of the distance
             if (w%rep%measure%item(i)%n == 2) then
                call w%rep%measure%item(i)%set_anchor(1,w%isys,win(iview)%vmdata%bidx(1:4))
                call w%rep%measure%item(i)%set_anchor(2,w%isys,win(iview)%vmdata%bidx(5:8))
                changed = .true.
             end if
          elseif (win(iview)%vmdata%idx(1) > 0 .and. k >= 1 .and. k <= w%rep%measure%item(i)%n) then
             call w%rep%measure%item(i)%set_anchor(k,w%isys,win(iview)%vmdata%idx(1:4))
             changed = .true.
          end if
          w%editrep_pick_item = 0
          w%editrep_pick_slot = 0
          win(iview)%vmdata%idx = 0
          win(iview)%vmdata%bidx = 0
       end if
    end if

    ! the toolbar, which makes and edits the measurements in the view:
    ! select, remove, and one tool per kind
    if (iview > 0) &
       call w%editrep_toolbar(ttshown)

    ! usage hint
    call iw_text("Measurements",highlight=.true.)
    call iw_helpermark("Measure with the toolbar: pick the atoms in the view window, two for a "//&
       "distance (or a bond), three for an angle (the second atom is the vertex), four for a "//&
       "dihedral. In navigation, "//trim(get_bind_keyname(BIND_NAV_MEASURE))//" selects atoms and "//&
       trim(get_bind_keyname(BIND_NAV_MEASURE_TOGGLE))//" on the last atom, or on empty space, "//&
       "adds the measurement to the first measurements object; the same click again removes it.")

    ! a measurement selected elsewhere (the view) brings up its tab
    itab = 0
    if (w%rep%measure%isel /= w%lastselected) then
       if (w%rep%measure%isel >= 1) itab = w%rep%measure%item(w%rep%measure%isel)%n
       w%lastselected = w%rep%measure%isel
    end if
    ! one tab per category: a selectable table, then the options for the
    ! selected row underneath
    str1 = "##editrepmeasuretabbar" // c_null_char
    flags = ImGuiTabBarFlags_None
    if (igBeginTabBar(c_loc(str1),flags)) then
       if (iw_begintabitem("Distances##editrepmeasure_disttab",flags=tabflags(2))) then
          call cat_table(2,"dist")
          call item_options(2)
          call igEndTabItem()
       end if
       if (iw_begintabitem("Angles##editrepmeasure_angtab",flags=tabflags(3))) then
          call cat_table(3,"ang")
          call item_options(3)
          call igEndTabItem()
       end if
       if (iw_begintabitem("Dihedrals##editrepmeasure_dihtab",flags=tabflags(4))) then
          call cat_table(4,"dih")
          call item_options(4)
          call igEndTabItem()
       end if
       call igEndTabBar()
    end if

    ! process a deletion (global item index), dropping a mouse drag on the
    ! measurements in the view
    if (idel > 0) then
       call measure_delete(w%rep%measure,idel)
       w%lastselected = w%rep%measure%isel ! the same selection: the tab stays
       if (iview > 0) then
          if (win(iview)%annot%irep == w%irep .and. win(iview)%oe%op /= objop_none) &
             win(iview)%oe%op = objop_none
       end if
       ! keep a pending atom pick bound to the right item (or cancel it if that
       ! item was the one deleted), so the pick can't commit to the wrong row
       if (w%editrep_pick_item == idel) then
          call win(iview)%viewmode_release_forced(w%id)
          w%editrep_pick_item = 0
          w%editrep_pick_slot = 0
       elseif (w%editrep_pick_item > idel) then
          w%editrep_pick_item = w%editrep_pick_item - 1
       end if
       changed = .true.
    end if

  contains

    !> Flags of the tab of the measurements with n atoms: selected if the
    !> selection moved to one of them.
    function tabflags(n) result(fl)
      integer, intent(in) :: n
      integer(c_int) :: fl

      fl = ImGuiTabItemFlags_None
      if (itab == n) fl = ImGuiTabItemFlags_SetSelected

    end function tabflags

    !> Draw the table of measurement items of the given category (ncat
    !> = number of atoms: 2 distances, 3 angles, 4 dihedrals). Every
    !> atom is a colored button that commands the view into pick mode
    !> to replace it. Every row has the value and is selectable (sets
    !> that editor).
    subroutine cat_table(ncat,tabidn)
      integer, intent(in) :: ncat
      character(len=*), intent(in) :: tabidn

      integer(c_int) :: tflags
      integer :: i, k, nrow, ncol, icvalue, iwrite
      logical :: ch, ldum

      ! count items of this category (for the table height)
      nrow = 0
      do i = 1, w%rep%measure%nitem
         if (w%rep%measure%item(i)%n == ncat) nrow = nrow + 1
      end do

      ! columns: delete, show, one button per atom (ncat of them), value
      ncol = ncat + 3
      icvalue = ncat + 2

      tflags = ImGuiTableFlags_None
      tflags = ior(tflags,ImGuiTableFlags_RowBg)
      tflags = ior(tflags,ImGuiTableFlags_Borders)
      tflags = ior(tflags,ImGuiTableFlags_ScrollY)
      tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
      str1 = "##measuretable" // tabidn // c_null_char
      sz0%x = 0
      sz0%y = iw_calcheight(min(nrow,6)+1,0,.false.)
      if (igBeginTable(c_loc(str1),ncol,tflags,sz0,0._c_float)) then
         call iw_table_column("",id=0,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Show",id=1,flags=ImGuiTableColumnFlags_WidthFixed)
         do k = 1, ncat
            call iw_table_column("Atom " // string(k),id=k+1,flags=ImGuiTableColumnFlags_WidthFixed)
         end do
         call iw_table_column("Value",id=icvalue,flags=ImGuiTableColumnFlags_WidthStretch)
         call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)
         ! (the table id is the number of atoms: 2-4; 1 is symmetry elements)
         call w%table_write_begin(ncat,iwrite)

         do i = 1, w%rep%measure%nitem
            if (w%rep%measure%item(i)%n /= ncat) cycle
            call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
            ! delete button
            if (igTableSetColumnIndex(0)) then
               call igAlignTextToFramePadding()
               if (iw_close_button("##measuredel" // tabidn // string(i))) idel = i
               call iw_tooltip("Remove this measurement",ttshown)
            end if
            ! shown checkbox
            if (igTableSetColumnIndex(1)) then
               if (iw_checkbox("##measureshow" // tabidn // string(i),w%rep%measure%item(i)%shown)) &
                  changed = .true.
               call iw_tooltip("Toggle show/hide this measurement",ttshown)
            end if
            ! one colored atom button per atom
            do k = 1, ncat
               if (igTableSetColumnIndex(k+1)) &
                  call atom_button(i,k,"##a" // string(k) // tabidn // string(i))
            end do
            ! value, plus a row-spanning selectable that picks the edited item
            if (igTableSetColumnIndex(icvalue)) then
               call iw_text(item_value_str(w%rep%measure%item(i)),alignframe=.true.)
               ldum = iw_highlight_selectable("##measuresel" // tabidn // string(i),&
                  clicked=ch,selected=(i == w%rep%measure%isel))
               if (ch) w%rep%measure%isel = i
            end if
         end do
         call w%table_write_end(ncat,"Measurements (" // string(ncat) // " atoms)")
         call igEndTable()
      end if
    end subroutine cat_table

    !> Draw the per-item options for the selected measurement, underneath its
    !> table, but only when the selected item belongs to this category (ncat).
    subroutine item_options(ncat)
      integer, intent(in) :: ncat

      integer :: is, j, imode
      integer(c_int) :: idec

      character(len=9), parameter :: catname(2:4) = (/"distances","angles   ","dihedrals"/)

      is = w%rep%measure%isel
      if (is < 1 .or. is > w%rep%measure%nitem) return
      if (w%rep%measure%item(is)%n /= ncat) return

      associate (it => w%rep%measure%item(is))
         ! color (the table has no color column for any category)
         changed = changed .or. iw_coloredit("Color##measureitemcol",rgb=it%rgb)
         call iw_tooltip("Color of this measurement",ttshown)
         ! decimals (with a trailing label to match the other options)
         idec = int(it%ndec,c_int)
         ! the tooltip goes on the stepper: a plain text has no item id, so a
         ! delayed (ttshown) tooltip attached to it never fires
         if (iw_intstepper("measureitemdec",idec,minval=0_c_int,maxval=8_c_int,sameline=.true.,&
            tooltip="Decimal places shown for this value")) then
            it%ndec = idec
            changed = .true.
         end if
         call iw_text("Decimals",sameline=.true.)
         ! segment/arm/edge radius
         changed = changed .or. iw_dragfloat_real8("Segment radius (Å)##measureitemrad",&
            x1=it%rad,speed=0.001d0,min=0.005d0,max=1d0,scale=bohrtoa,decimal=3,&
            flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Radius of the segments, arms, and dihedral edges",ttshown)
         ! distances: a bond fills both ends in one pick, then the line style
         if (ncat == 2) then
            if (iw_button("Pick bond##measurepickbond",disabled=(w%editrep_pick_item > 0))) then
               w%editrep_pick_item = is
               w%editrep_pick_slot = 0 ! a bond fills both ends, not one slot
               call w%editrep_pick%arm()
               call win(iview)%viewmode_set_forced(vm_pick_bond,&
                  "Pick the bond to measure",w%id)
            end if
            call iw_tooltip("Click, then pick a bond in the view window: its two atoms become&
               & the ends of this distance",ttshown)
            changed = changed .or. iw_checkbox("Dashed line##measureitemdash",it%dashed)
            call iw_tooltip("Draw the segment as a dashed line instead of a solid one",ttshown)
            if (it%dashed) then
               changed = changed .or. iw_dragfloat_real8("Dash length (Å)##measureitemdashlen",&
                  x1=it%dashlen,speed=0.005d0,min=0.02d0,max=5d0,scale=bohrtoa,decimal=2,&
                  flags=ImGuiSliderFlags_AlwaysClamp)
               call iw_tooltip("Length of each dash (and gap) along the segment",ttshown)
            end if
         end if
         ! sector (angles and dihedrals)
         if (ncat >= 3) then
            changed = changed .or. iw_dragfloat_real8("Sector radius (Å)##measureitemsrad",&
               x1=it%sectorrad,speed=0.01d0,min=0.05d0,max=10d0,scale=bohrtoa,decimal=2,&
               flags=ImGuiSliderFlags_AlwaysClamp)
            call iw_tooltip("Radius of the angle/dihedral sector",ttshown)
            changed = changed .or. iw_dragfloat_real8("Sector opacity##measureitemsalpha",&
               x1=it%sectoralpha,speed=0.01d0,min=0d0,max=1d0,decimal=2,&
               flags=ImGuiSliderFlags_AlwaysClamp)
            call iw_tooltip("Opacity of the angle/dihedral sector",ttshown)
         end if
         ! plane opacity and extension (dihedrals)
         if (ncat == 4) then
            changed = changed .or. iw_dragfloat_real8("Plane opacity##measureitempalpha",&
               x1=it%planealpha,speed=0.01d0,min=0d0,max=1d0,decimal=2,&
               flags=ImGuiSliderFlags_AlwaysClamp)
            call iw_tooltip("Opacity of the dihedral planes",ttshown)
            changed = changed .or. iw_dragfloat_real8("Plane extension (%)##measureitempext",&
               x1=it%planeext,speed=1d0,min=0d0,max=100d0,scale=100d0,decimal=0,&
               flags=ImGuiSliderFlags_AlwaysClamp)
            call iw_tooltip("How far each plane rectangle extends beyond the atoms, "//&
               "as a percentage of its size on each side",ttshown)
         end if
         ! label orientation (distances): along the segment, or always horizontal
         if (ncat == 2) then
            changed = changed .or. iw_checkbox("Orient label along segment##measureitemorient",it%orient)
            call iw_tooltip("Run the label along the segment; if unchecked, keep it horizontal",ttshown)
         end if
         ! label scaling: constant on-screen size, or scaling with the system
         changed = changed .or. iw_checkbox("Scale label with system size##measureitemscsys",it%scalesystem)
         call iw_tooltip("Checked: the label scales with the system (grows when zooming in). "//&
            "Unchecked: it keeps a constant, legible on-screen size for any system size.",ttshown)
         ! label size
         changed = changed .or. iw_dragfloat_real8("Label size##measureitemtscale",&
            x1=it%textscale,speed=0.01d0,min=0.05d0,max=10d0,decimal=2,&
            flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Size of the numeric label",ttshown)
         ! in-plane (screen-relative) offset of the label from its anchor
         changed = changed .or. iw_dragfloat_real8("Label offset (Å)##measureitemoffset",&
            x2=it%offset,speed=0.01d0,decimal=2)
         call iw_tooltip("Offset of the label from its anchor, in the screen plane (angstrom)",ttshown)
         ! apply these style options to the measurements in this tab, or to all
         imode = apply_style_buttons("measure",trim(catname(ncat)),"measurements",ttshown)
         if (imode > 0) then
            do j = 1, w%rep%measure%nitem
               if (j == is .or. (imode == 1 .and. w%rep%measure%item(j)%n /= ncat)) cycle
               call w%rep%measure%item(j)%copy_style(it)
            end do
            changed = .true.
         end if
      end associate
    end subroutine item_options

    !> Draw a button for atom islot of measurement item iitem (cell
    !> atom id + lattice vector), tinted with its color in the parent
    !> view. Label = species name + cell id ("?" if the atom is
    !> stale). Clicking it commands the parent view into pick mode.
    subroutine atom_button(iitem,islot,idn)
      integer, intent(in) :: iitem, islot
      character(len=*), intent(in) :: idn

      logical :: clicked
      character(len=:), allocatable :: lbl

      clicked = draw_anchor_button(iview,w%isys,w%rep%measure%item(iitem)%idx(:,islot),idn,&
         disabled=(w%editrep_pick_item > 0),lbl=lbl,stamp=w%rep%measure%item(iitem)%stamp(islot))
      if (w%rep%measure%item(iitem)%idx(1,islot) < 0) then
         call iw_tooltip("Critical point " // lbl // " (click to pick an atom to replace it in the view)",ttshown)
      else
         call iw_tooltip("Atom " // lbl // " (click to pick a replacement in the view)",ttshown)
      end if

      ! command the parent view into pick mode; the poll at the top of
      ! draw_editrep_measure commits the picked atom into this slot
      if (clicked) then
         w%editrep_pick_item = iitem
         w%editrep_pick_slot = islot
         call w%editrep_pick%arm()
         call win(iview)%viewmode_set_forced(vm_pick_atom,&
            "Pick the replacement atom",w%id)
      end if
    end subroutine atom_button

    !> Current value of a measurement, formatted with units (or "(stale)" if an
    !> anchor atom no longer exists).
    function item_value_str(item) result(s)
      type(measurement_item), intent(in) :: item
      character(len=:), allocatable :: s
      real*8 :: f(3,4), val
      if (.not.item%anchors_xfrac(w%isys,f)) then
         s = "(stale)"
         return
      end if
      if (item%n == 2) then
         val = sys(w%isys)%c%distance(f(:,1),f(:,2)) * bohrtoa
         s = string(val,'f',decimal=item%ndec) // " Å"
      elseif (item%n == 3) then
         val = sys(w%isys)%c%angle(f(:,1),f(:,2),f(:,3)) * 180d0 / pi
         s = string(val,'f',decimal=item%ndec) // "°"
      else
         val = sys(w%isys)%c%dihedral(f(:,1),f(:,2),f(:,3),f(:,4)) * 180d0 / pi
         s = string(val,'f',decimal=item%ndec) // "°"
      end if
    end function item_value_str

  end function draw_editrep_measure

  !> Draw the editrep window, shapes class. Returns true if the scene needs
  !> rendering again.
  module function draw_editrep_shapes(w,ttshown) result(changed)
    use representations, only: rep_shape, shapekind_sphere, shapekind_box, shapekind_arrow,&
       shapekind_cone, shapekind_NUM, shapekind_name, shapekind_combostr,&
       shape_size_def, shape_edge_def, arrow_length_def, arrow_radius_def, shapes_delete,&
       shape_copy_style, shape_proportions
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_checkbox, iw_coloredit,&
       iw_dragfloat_real8, iw_dragfloat_realc, iw_combo_simple, iw_button, iw_calcheight,&
       iw_close_button, iw_highlight_selectable, iw_table_column
    use systems, only: sys
    use tools_io, only: string, lower
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch, ldum
    integer :: i, k, iview, isel, idel, ikind, istat, ndec, imode
    real*8 :: xc(3), xdsp(3)
    integer(c_int) :: flags
    type(ImVec2) :: sz0
    character(kind=c_char,len=:), allocatable, target :: str1
    character(len=:), allocatable :: strunit

    ! initialize
    changed = .false.
    iview = w%anchor_view()

    ! positions are shown the way the atomic coordinates are shown elsewhere:
    ! fractional in a crystal, cartesian angstrom in a molecule
    if (sys(w%isys)%c%ismolecule) then
       strunit = "Å"
       ndec = 4
    else
       strunit = "fractional"
       ndec = 6
    end if

    ! handle a pending position pick commanded to the parent view
    if (w%editrep_pick_item > 0) then
       if (w%editrep_pick_item > w%rep%shapes%nshape) then
          ! the shape was deleted under the pick
          call win(iview)%viewmode_release_forced(w%id)
          w%editrep_pick_item = 0
       else
          call view_pick_result(iview,w%id,w%isys,w%editrep_pick,istat,xc)
          if (istat == ipick_point) then
             ! the shapes are anchored in the absolute frame
             if (sys(w%isys)%c%ismolecule) xc = xc + sys(w%isys)%c%molx0
             associate (sh => w%rep%shapes%shape(w%editrep_pick_item))
               if (w%editrep_pick_slot == 0) then
                  ! the anchor: the end points stay where they are, so the
                  ! vectors to them absorb the move (a sphere has none)
                  if (sh%kind == shapekind_box) then
                     do k = 1, 3
                        sh%v(:,k) = sh%v(:,k) + sh%x1 - xc
                     end do
                  elseif (sh%kind /= shapekind_sphere) then
                     sh%v(:,1) = sh%v(:,1) + sh%x1 - xc
                  end if
                  sh%x1 = xc
               else
                  ! an end point: only that vector changes
                  sh%v(:,w%editrep_pick_slot) = xc - sh%x1
               end if
             end associate
             changed = .true.
          end if
          if (istat /= ipick_pending) w%editrep_pick_item = 0
       end if
    end if

    ! the toolbar, which edits the shapes in the view: select, remove,
    ! and one tool per kind of shape
    if (iview > 0) &
       call w%editrep_toolbar(ttshown)

    ! table of shapes
    call iw_text("3D Shapes",highlight=.true.)
    idel = 0
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    str1 = "##shapestable" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(w%rep%shapes%nshape,5)+1,0,.false.)
    if (igBeginTable(c_loc(str1),4,flags,sz0,0._c_float)) then
       call iw_table_column("",id=0,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Show",id=1,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Kind",id=2,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Position (" // strunit // ")",id=3,&
          flags=ImGuiTableColumnFlags_WidthStretch)
       call iw_table_headers_row(freezetop=.true.)

       do i = 1, w%rep%shapes%nshape
          call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)

          ! delete button
          if (igTableSetColumnIndex(0)) then
             call igAlignTextToFramePadding()
             if (iw_close_button("##shapedel" // string(i))) idel = i
             call iw_tooltip("Remove this shape",ttshown)
          end if

          ! shown checkbox
          if (igTableSetColumnIndex(1)) then
             if (iw_checkbox("##shapeshow" // string(i),w%rep%shapes%shape(i)%shown)) changed = .true.
             call iw_tooltip("Toggle show/hide this shape",ttshown)
          end if

          ! kind
          if (igTableSetColumnIndex(2)) &
             call iw_text(kind_label(w%rep%shapes%shape(i)%kind),alignframe=.true.)

          ! position, and the row-spanning selectable that picks the edited shape
          if (igTableSetColumnIndex(3)) then
             xdsp = pos_to_display(w%rep%shapes%shape(i)%x1)
             call iw_text("(" // string(xdsp(1),'f',decimal=ndec) // ", " //&
                string(xdsp(2),'f',decimal=ndec) // ", " // string(xdsp(3),'f',decimal=ndec) //&
                ")",alignframe=.true.)
             ldum = iw_highlight_selectable("##shapesel" // string(i),clicked=ch,&
                selected=(i == w%rep%shapes%isel))
             if (ch) w%rep%shapes%isel = i
          end if
       end do
       call igEndTable()
    end if

    ! process a deletion, dropping a mouse drag on the shapes in the view
    if (idel > 0) then
       call shapes_delete(w%rep%shapes,idel)
       if (win(iview)%oe%op == objop_handle .or. win(iview)%oe%op == objop_move) &
          win(iview)%oe%op = objop_none
       if (w%editrep_pick_item == idel) then
          ! cancel a pick pending on the deleted shape
          call win(iview)%viewmode_release_forced(w%id)
          w%editrep_pick_item = 0
       elseif (w%editrep_pick_item > idel) then
          w%editrep_pick_item = w%editrep_pick_item - 1
       end if
       changed = .true.
    end if
    if (w%rep%shapes%isel < 1 .or. w%rep%shapes%isel > w%rep%shapes%nshape) return

    ! options for the selected shape. Each widget is evaluated into ch first:
    ! .or. is allowed to short-circuit, and a widget skipped because changed is
    ! already true is a widget not drawn
    isel = w%rep%shapes%isel
    associate (sh => w%rep%shapes%shape(isel))
      call iw_text("Shape " // string(isel) // " (" // kind_label(sh%kind) // "), positions in " //&
         strunit,highlight=.true.)

      ! the kind; a shape that changes kind keeps its position and gets a
      ! default vector if it has none to show
      ikind = sh%kind
      call iw_combo_simple("Kind##shapekind",shapekind_combostr,ikind,changed=ch,&
         startsatone=.true.)
      if (ch) then
         sh%kind = ikind
         ! the radius means a different thing for every kind (sphere radius,
         ! box edge thickness, shaft thickness...), so it always goes back to
         ! the default; the position is kept, and so is every vector the new
         ! kind can use (seed_vectors only fills in the ones it lacks); a cone
         ! or a cylinder then takes its width from its length
         call seed_radius(sh)
         call seed_vectors(sh)
         call shape_proportions(sh)
      end if
      changed = changed .or. ch
      call iw_tooltip("Kind of 3D shape",ttshown)

      ! the anchor; moving it translates the whole shape, since the end
      ! points below are stored as vectors from it
      xdsp = pos_to_display(sh%x1)
      ch = iw_dragfloat_real8("Position##shapepos",x3=xdsp,speed=0.001d0,&
         decimal=ndec,notlive=.true.)
      if (ch) sh%x1 = pos_from_display(xdsp)
      call iw_tooltip("Position of the shape: the center of a sphere, a corner of a box, the &
         &tail of an arrow, the center of the base of a cone or the first end of a cylinder",ttshown)
      changed = changed .or. ch
      call shape_pick_button("shapepick",0,"Pick the position for the shape")

      ! the end points, one per axis for a box and one for the rest
      if (sh%kind == shapekind_box) then
         do k = 1, 3
            xdsp = pos_to_display(sh%x1 + sh%v(:,k))
            ch = iw_dragfloat_real8("Axis-" // string(k) // " end##shapeend" // string(k),&
               x3=xdsp,speed=0.001d0,decimal=ndec,notlive=.true.)
            if (ch) sh%v(:,k) = pos_from_display(xdsp) - sh%x1
            call iw_tooltip("End point of axis " // string(k) // " of the box, i.e. the corner &
               &reached from the position along that axis",ttshown)
            changed = changed .or. ch
            call shape_pick_button("shapeendpick" // string(k),k,&
               "Pick the end point of axis " // string(k))
         end do
      elseif (sh%kind /= shapekind_sphere) then
         xdsp = pos_to_display(sh%x1 + sh%v(:,1))
         ch = iw_dragfloat_real8("End point##shapeend1",x3=xdsp,&
            speed=0.001d0,decimal=ndec,notlive=.true.)
         if (ch) sh%v(:,1) = pos_from_display(xdsp) - sh%x1
         call iw_tooltip("Tip of the arrow, apex of the cone, or the other end of the cylinder",&
            ttshown)
         changed = changed .or. ch
         call shape_pick_button("shapeendpick1",1,"Pick the end point")
      end if

      ! the radius/thickness
      if (sh%kind == shapekind_sphere) then
         str1 = "Radius (Å)##shaperad"
      elseif (sh%kind == shapekind_box) then
         str1 = "Edge thickness (Å)##shaperad"
      elseif (sh%kind == shapekind_arrow) then
         str1 = "Shaft thickness (Å)##shaperad"
      elseif (sh%kind == shapekind_cone) then
         str1 = "Base width (Å)##shaperad"
      else
         str1 = "Thickness (Å)##shaperad"
      end if
      ch = iw_dragfloat_real8(str1,x1=sh%rad,speed=0.002d0,min=0d0,max=20d0,scale=bohrtoa,&
         decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
      call iw_tooltip("Radius of the sphere, thickness of the box edges or of the cylinder &
         &and the arrow shaft, or width of the base of the cone",ttshown)
      changed = changed .or. ch

      ! the arrowhead
      if (sh%kind == shapekind_arrow) then
         ch = iw_dragfloat_real8("Head Size##shapeheadr",x1=sh%headr,speed=0.02d0,min=1d0,&
            max=6d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Width of the arrowhead, in units of the shaft thickness",ttshown)
         changed = changed .or. ch

         ch = iw_dragfloat_real8("Head Length##shapeheadl",x1=sh%headl,speed=0.005d0,min=0d0,&
            max=1d0,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Length of the arrowhead, as a fraction of the whole arrow",ttshown)
         changed = changed .or. ch
      end if

      ! the color, and the opacity for the kinds that are not drawn opaque
      ch = iw_coloredit("Color##shapecolor",rgb=sh%rgb)
      call iw_tooltip("Color of the shape",ttshown)
      changed = changed .or. ch

      if (sh%kind == shapekind_sphere .or. sh%kind == shapekind_box) then
         ch = iw_dragfloat_realc("Opacity##shapealpha",x1=sh%alpha,speed=0.01_c_float,&
            min=0._c_float,max=1._c_float,decimal=2,sameline=.true.,&
            flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Opacity of the object",ttshown)
         changed = changed .or. ch
      end if

      ! apply these style options to the shapes of this kind, or to all
      imode = apply_style_buttons("shapes",lower(kind_label(sh%kind)) // " shapes","shapes",ttshown)
      if (imode > 0) then
         do i = 1, w%rep%shapes%nshape
            if (i == isel .or. (imode == 1 .and. w%rep%shapes%shape(i)%kind /= sh%kind)) cycle
            call shape_copy_style(w%rep%shapes%shape(i),sh)
         end do
         changed = .true.
      end if
    end associate

  contains
    !> Convert the absolute-frame cartesian position xc (bohr) to the units
    !> the editor shows: fractional for a crystal, angstrom for a molecule.
    function pos_to_display(xc) result(x)
      real*8, intent(in) :: xc(3)
      real*8 :: x(3)

      if (sys(w%isys)%c%ismolecule) then
         x = xc * bohrtoa
      else
         x = sys(w%isys)%c%c2x(xc)
      end if

    end function pos_to_display

    !> Inverse of pos_to_display.
    function pos_from_display(x) result(xc)
      real*8, intent(in) :: x(3)
      real*8 :: xc(3)

      if (sys(w%isys)%c%ismolecule) then
         xc = x / bohrtoa
      else
         xc = sys(w%isys)%c%x2c(x)
      end if

    end function pos_from_display

    !> Button that commands the parent view to pick a position for the shape
    !> being edited. islot is the destination: 0 the anchor, k the end point
    !> of vector k.
    subroutine shape_pick_button(id,islot,tt)
      character(len=*), intent(in) :: id
      integer, intent(in) :: islot
      character(len=*), intent(in) :: tt

      if (iw_button("Pick##" // id,sameline=.true.,disabled=(w%editrep_pick_item > 0))) then
         w%editrep_pick_item = isel
         w%editrep_pick_slot = islot
         call w%editrep_pick%arm()
         call win(iview)%viewmode_set_forced(vm_pick_atom,tt,w%id,acceptempty=.true.)
      end if
      call iw_tooltip(tt // " (click an atom, or empty space for a point on the plane through &
         &the scene center)",ttshown)

    end subroutine shape_pick_button

    !> Name of shape kind k, for the table and the per-shape header.
    function kind_label(k) result(str)
      integer, intent(in) :: k
      character(len=:), allocatable :: str

      if (k >= 1 .and. k <= shapekind_NUM) then
         str = trim(shapekind_name(k))
      else
         str = "???"
      end if

    end function kind_label

    !> Fill in the geometry vectors that the kind of the shape sh needs and
    !> does not have yet; a vector it already has is kept, so a shape keeps
    !> its size across a kind change. A zero vector would draw nothing (a box
    !> with one would be a flat or collapsed parallelepiped).
    subroutine seed_vectors(sh)
      type(rep_shape), intent(inout) :: sh

      integer :: j

      if (sh%kind == shapekind_sphere) return
      if (sh%kind == shapekind_box) then
         do j = 1, 3
            if (all(abs(sh%v(:,j)) < 1d-10)) sh%v(j,j) = shape_size_def
         end do
      elseif (all(abs(sh%v(:,1)) < 1d-10)) then
         sh%v(1,1) = arrow_length_def
      end if

    end subroutine seed_vectors

    !> Default radius for the kind of the shape sh: the sphere radius, the
    !> thickness of the box edges, or the shaft/base radius of the rest.
    subroutine seed_radius(sh)
      type(rep_shape), intent(inout) :: sh

      if (sh%kind == shapekind_sphere) then
         sh%rad = 0.5d0 * shape_size_def
      elseif (sh%kind == shapekind_box) then
         sh%rad = shape_edge_def
      else
         sh%rad = arrow_radius_def
      end if

    end subroutine seed_radius

  end function draw_editrep_shapes

  !> Draw the editrep window, planar shapes class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_planar(w,ttshown) result(changed)
    use representations, only: planar_shape, planarkind_name,&
       planarkind_ellipse, planarkind_rect, planarkind_arrow, planarkind_freehand,&
       planarkind_curve, planarheads_combostr, planardash_combostr, planarfill_combostr,&
       planarfill_hatched, planarfill_crosshatched, planar_isclosed, planar_haspoints, planar_delete,&
       planar_copy_style
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_checkbox, iw_coloredit,&
       iw_dragfloat_real8, iw_dragfloat_realc, iw_combo_simple, iw_button, iw_calcheight,&
       iw_close_button, iw_highlight_selectable, iw_table_column
    use tools_io, only: string, lower
    use param, only: pi
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch, ldum
    integer :: i, k, iview, isel, idel, iswap, ihead, idash, ifill, imode
    integer(c_int) :: flags
    real*8 :: xdsp(2), pxs, wpx, angd, dx, dy
    type(ImVec2) :: sz0
    character(kind=c_char,len=:), allocatable, target :: str1
    type(planar_shape) :: shaux

    character(len=6), parameter :: curvept(3) = (/"Start ","Middle","End   "/)

    ! initialize
    changed = .false.
    iview = w%anchor_view()
    if (iview == 0) return

    ! NDC per screen pixel in the view: the render buffer spans the
    ! longest side of the view image
    dx = win(iview)%v_rmax%x - win(iview)%v_rmin%x
    dy = win(iview)%v_rmax%y - win(iview)%v_rmin%y
    pxs = 2d0 / max(dx,dy,1d0)

    ! the toolbar, which edits the drawing in the view: select, remove,
    ! and one tool per kind of shape
    call w%editrep_toolbar(ttshown)

    ! table of shapes
    call iw_text("2D Drawing",highlight=.true.)
    idel = 0
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    str1 = "##planartable" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(w%rep%planar%nshape,5)+1,0,.false.)
    if (igBeginTable(c_loc(str1),4,flags,sz0,0._c_float)) then
       call iw_table_column("",id=0,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Show",id=1,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Kind",id=2,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Position (view)",id=3,flags=ImGuiTableColumnFlags_WidthStretch)
       call iw_table_headers_row(freezetop=.true.)

       do i = 1, w%rep%planar%nshape
          call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)

          ! delete button
          if (igTableSetColumnIndex(0)) then
             call igAlignTextToFramePadding()
             if (iw_close_button("##planardel" // string(i))) idel = i
             call iw_tooltip("Remove this shape",ttshown)
          end if

          ! shown checkbox
          if (igTableSetColumnIndex(1)) then
             if (iw_checkbox("##planarshow" // string(i),w%rep%planar%shape(i)%shown)) changed = .true.
             call iw_tooltip("Toggle show/hide this shape",ttshown)
          end if

          ! kind
          if (igTableSetColumnIndex(2)) &
             call iw_text(kind_label(w%rep%planar%shape(i)%kind),alignframe=.true.)

          ! position, and the row-spanning selectable that picks the edited shape
          if (igTableSetColumnIndex(3)) then
             xdsp = anchor_pos(w%rep%planar%shape(i))
             call iw_text("(" // string(xdsp(1),'f',decimal=3) // ", " //&
                string(xdsp(2),'f',decimal=3) // ")",alignframe=.true.)
             ldum = iw_highlight_selectable("##planarsel" // string(i),clicked=ch,&
                selected=(i == w%rep%planar%isel))
             if (ch) then
                w%rep%planar%isel = i
                win(iview)%forcerender = .true.
             end if
          end if
       end do
       call igEndTable()
    end if

    ! process a deletion, dropping a mouse drag on the shapes in the view
    if (idel > 0) then
       call planar_delete(w%rep%planar,idel)
       if (win(iview)%oe%op == objop_handle .or. win(iview)%oe%op == objop_move) &
          win(iview)%oe%op = objop_none
       changed = .true.
    end if
    if (w%rep%planar%nshape == 0) return
    if (w%rep%planar%isel < 1 .or. w%rep%planar%isel > w%rep%planar%nshape) return

    ! options for the selected shape. Each widget is evaluated into ch first:
    ! .or. is allowed to short-circuit, and a widget skipped because changed is
    ! already true is a widget not drawn
    isel = w%rep%planar%isel
    iswap = 0
    associate (sh => w%rep%planar%shape(isel))
      call iw_text("Shape " // string(isel) // " (" // kind_label(sh%kind) // ")",highlight=.true.)

      ! drawing order
      if (iw_button("Raise##planarraise",sameline=.true.,disabled=(isel == w%rep%planar%nshape))) &
         iswap = isel + 1
      call iw_tooltip("Draw this shape after (on top of) the next one",ttshown)
      if (iw_button("Lower##planarlower",sameline=.true.,disabled=(isel == 1))) &
         iswap = isel - 1
      call iw_tooltip("Draw this shape before (under) the previous one",ttshown)

      ! geometry, in view coordinates (the render buffer spans -1 to 1)
      if (sh%kind == planarkind_ellipse .or. sh%kind == planarkind_rect) then
         ch = iw_dragfloat_real8("Center##planarxc",x2=sh%xc,speed=0.002d0,min=-10d0,max=10d0,&
            decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Position of the center (view coordinates: -1 to 1 across the render &
            &buffer, which spans the longest side of the view)",ttshown)
         changed = changed .or. ch
         ch = iw_dragfloat_real8("Half-sizes##planarhs",x2=sh%hs,speed=0.002d0,min=0d0,max=4d0,&
            decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Half of the width and of the height of the shape, before the rotation &
            &(view coordinates)",ttshown)
         changed = changed .or. ch
         angd = sh%ang * 180d0 / pi
         ch = iw_dragfloat_real8("Angle (°)##planarang",x1=angd,speed=0.5d0,min=-180d0,max=180d0,&
            decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
         if (ch) sh%ang = angd * pi / 180d0
         call iw_tooltip("Rotation of the shape, counterclockwise",ttshown)
         changed = changed .or. ch
      elseif (sh%kind == planarkind_freehand) then
         call iw_text("Freehand line with " // string(sh%npt) // " points")
      elseif (planar_haspoints(sh)) then
         do k = 1, sh%npt
            if (sh%kind == planarkind_curve) then
               str1 = trim(curvept(k))
            else
               str1 = "Point " // string(k)
            end if
            ch = iw_dragfloat_real8(str1 // "##planarpt" // string(k),x2=sh%x(:,k),&
               speed=0.002d0,min=-10d0,max=10d0,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
            if (sh%kind == planarkind_curve .and. k == 2) then
               call iw_tooltip("Point midway along the curve: it sets the bend (view coordinates)",&
                  ttshown)
            else
               call iw_tooltip("Position of " // lower(str1) // " (view coordinates)",ttshown)
            end if
            changed = changed .or. ch
         end do
      end if

      ! the outline
      ch = iw_checkbox("Outline##planarstroke",sh%stroke)
      call iw_tooltip("Draw the outline of the shape",ttshown)
      changed = changed .or. ch
      wpx = sh%width / pxs
      ch = iw_dragfloat_real8("Width (px)##planarwidth",x1=wpx,speed=0.1d0,min=0.5d0,max=100d0,&
         decimal=1,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
      if (ch) sh%width = wpx * pxs
      call iw_tooltip("Width of the outline, in pixels of the view at its current size (it &
         &scales with the view and the exported image)",ttshown)
      changed = changed .or. ch
      idash = sh%dash
      call iw_combo_simple("Style##planardash",planardash_combostr,idash,changed=ch,sameline=.true.)
      if (ch) sh%dash = idash
      call iw_tooltip("Style of the outline: solid, dashed, or dotted (the dash and dot spacing &
         &scales with the width)",ttshown)
      changed = changed .or. ch
      ch = iw_coloredit("Color##planarrgb",rgb=sh%rgb)
      call iw_tooltip("Color of the outline",ttshown)
      changed = changed .or. ch
      ch = iw_dragfloat_realc("Opacity##planaralpha",x1=sh%alpha,speed=0.01_c_float,&
         min=0._c_float,max=1._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
      call iw_tooltip("Opacity of the outline",ttshown)
      changed = changed .or. ch

      ! the arrowheads
      if (sh%kind == planarkind_arrow .or. sh%kind == planarkind_curve) then
         ihead = sh%heads
         call iw_combo_simple("Heads##planarheads",planarheads_combostr,ihead,changed=ch)
         if (ch) sh%heads = ihead
         call iw_tooltip("Arrowheads: none (a line), at the end, at the start, or at both ends",&
            ttshown)
         changed = changed .or. ch
         ch = iw_dragfloat_real8("Head Length##planarheadl",x1=sh%headl,speed=0.05d0,min=1d0,&
            max=20d0,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Length of the arrowhead, in units of the outline width",ttshown)
         changed = changed .or. ch
         ch = iw_dragfloat_real8("Head Width##planarheadw",x1=sh%headw,speed=0.05d0,min=1d0,&
            max=20d0,decimal=1,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Width of the arrowhead, in units of the outline width",ttshown)
         changed = changed .or. ch
      end if

      ! the fill
      if (planar_isclosed(sh)) then
         ifill = sh%filltype
         call iw_combo_simple("Fill##planarfilltype",planarfill_combostr,ifill,changed=ch)
         if (ch) sh%filltype = ifill
         call iw_tooltip("Fill of the inside of the shape: none, solid, or lines (hatched) or &
            &crossed lines (cross-hatched) in the fill color",ttshown)
         changed = changed .or. ch
         ch = iw_coloredit("Fill Color##planarfillrgb",rgb=sh%fillrgb,sameline=.true.)
         call iw_tooltip("Color of the inside of the shape",ttshown)
         changed = changed .or. ch
         ch = iw_dragfloat_realc("Fill Opacity##planarfillalpha",x1=sh%fillalpha,speed=0.01_c_float,&
            min=0._c_float,max=1._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Opacity of the inside of the shape",ttshown)
         changed = changed .or. ch

         ! the hatch lines
         if (sh%filltype == planarfill_hatched .or. sh%filltype == planarfill_crosshatched) then
            angd = sh%hatchang * 180d0 / pi
            ch = iw_dragfloat_real8("Angle (°)##planarhatchang",x1=angd,speed=0.5d0,min=-90d0,&
               max=90d0,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
            if (ch) sh%hatchang = angd * pi / 180d0
            call iw_tooltip("Angle of the hatch lines, counterclockwise from the horizontal",ttshown)
            changed = changed .or. ch
            wpx = sh%hatchsp / pxs
            ch = iw_dragfloat_real8("Spacing (px)##planarhatchsp",x1=wpx,speed=0.1d0,min=2d0,&
               max=200d0,decimal=1,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
            if (ch) sh%hatchsp = wpx * pxs
            call iw_tooltip("Distance between the hatch lines, in pixels of the view at its &
               &current size",ttshown)
            changed = changed .or. ch
            wpx = sh%hatchw / pxs
            ch = iw_dragfloat_real8("Line Width (px)##planarhatchw",x1=wpx,speed=0.05d0,min=0.5d0,&
               max=50d0,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
            if (ch) sh%hatchw = wpx * pxs
            call iw_tooltip("Width of the hatch lines, in pixels of the view at its current size",&
               ttshown)
            changed = changed .or. ch
         end if
      end if

      ! in front of or behind the scene
      ch = iw_checkbox("In front of the scene##planarinfront",sh%infront)
      call iw_tooltip("Draw the shape on top of the scene (checked) or behind it, only over &
         &the background (unchecked)",ttshown)
      changed = changed .or. ch

      ! apply these style options to the shapes of this kind, or to all
      imode = apply_style_buttons("planar",lower(kind_label(sh%kind)) // " shapes","shapes",ttshown)
      if (imode > 0) then
         do i = 1, w%rep%planar%nshape
            if (i == isel .or. (imode == 1 .and. w%rep%planar%shape(i)%kind /= sh%kind)) cycle
            call planar_copy_style(w%rep%planar%shape(i),sh)
         end do
         changed = .true.
      end if
    end associate

    ! a change in the drawing order, keeping the shape selected
    if (iswap > 0) then
       shaux = w%rep%planar%shape(iswap)
       w%rep%planar%shape(iswap) = w%rep%planar%shape(isel)
       w%rep%planar%shape(isel) = shaux
       w%rep%planar%isel = iswap
       changed = .true.
    end if

  contains
    !> Name of planar shape kind k, for the table and the per-shape header.
    function kind_label(k) result(str)
      integer, intent(in) :: k
      character(len=:), allocatable :: str

      str = trim(planarkind_name(k))

    end function kind_label

    !> The position shown for a shape in the table: the center of an
    !> ellipse or rectangle, the first point of the rest.
    function anchor_pos(sh) result(x)
      type(planar_shape), intent(in) :: sh
      real*8 :: x(2)

      if (planar_haspoints(sh) .and. sh%npt > 0) then
         x = sh%x(:,1)
      else
         x = sh%xc
      end if

    end function anchor_pos

  end function draw_editrep_planar

  !> Draw the toolbar of the object editor w: one icon button per tool
  !> of its object type (objtool_list). A tool button arms the object
  !> editing mode of the view on this object (annot_set_tool), and
  !> clicking the armed tool again turns it off. The view owns the tool,
  !> so a tool shows as armed here when the view works on this object
  !> with it, whether it was armed here or in the annotation row of the
  !> view. The cancel bind turns the tool off while this window has the
  !> focus (the view handles it when it does).
  module subroutine editrep_toolbar(w,ttshown)
    use utils, only: iw_text, iw_tooltip, iw_icon_togglebutton, iw_push_iconrow_frame,&
       iw_pop_iconrow_frame
    use gui_main, only: tooltip_enabled
    use keybindings, only: is_bind_event, BIND_CANCEL
    use icons, only: icon_tex
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown

    integer :: iview, i, k, itool, itype
    logical :: ldum
    integer, allocatable :: tools(:), icons(:)
    character(len=2), allocatable :: falls(:)

    iview = w%anchor_view()
    if (iview == 0) return
    itype = w%rep%type

    ! the tool of the view, if it works on this object
    itool = objtool_none
    if (win(iview)%viewmode == vm_objedit .and. win(iview)%vmdata%owner == iview .and.&
       win(iview)%annot%irep == w%irep) itool = win(iview)%annot%itool
    if (itool /= objtool_none .and. w%focused()) then
       if (is_bind_event(BIND_CANCEL,norepeat=.true.)) then
          call win(iview)%annot_set_tool(0,objtool_none,0)
          itool = objtool_none
       end if
    end if

    call iw_text("Toolbar",highlight=.true.)
    call iw_push_iconrow_frame()
    call objtool_list(itype,tools,icons,falls)
    do i = 1, size(tools)
       k = tools(i)
       ldum = (itool == k)
       if (iw_icon_togglebutton("##objtool" // string(k),icon_tex(icons(i)),trim(falls(i)),&
          state=ldum,sameline=(i > 1))) then
          ! a tool takes the view from a pick commanded by this window's
          ! tables, which is then abandoned
          w%editrep_pick_item = 0
          if (ldum) then
             call win(iview)%annot_set_tool(itype,k,w%irep)
          else
             call win(iview)%annot_set_tool(0,objtool_none,0)
          end if
       end if
       if (tooltip_enabled) then
          if (igIsItemHovered(ImGuiHoveredFlags_None)) call iw_tooltip(objtool_hint(itype,k),ttshown)
       end if
    end do
    call iw_pop_iconrow_frame()

  end subroutine editrep_toolbar

  !> The tools of the object type itype (reptype_planar, _shapes, _text,
  !> _measure, _axes): select, remove (not for the axes) and one per kind
  !> of item it draws (objtool_kind0 + kind), with their icons (icon_*
  !> ids) and the texts drawn instead if they did not load.
  module subroutine objtool_list(itype,tools,icons,falls)
    use representations, only: planarkind_NUM, shapekind_NUM, textpos_bond, axplace_scene,&
       axplace_window
    use icons, only: icon_pl_select, icon_vm_remove, icon_pl_ellipse, icon_pl_rect,&
       icon_pl_polygon, icon_pl_polyline, icon_pl_arrow, icon_pl_freehand, icon_pl_curve,&
       icon_sh_sphere, icon_ui_cell, icon_sh_cone, icon_sh_cylinder, icon_tx_screen,&
       icon_tx_point, icon_tx_atom, icon_tx_bond, icon_ms_distance, icon_ms_angle,&
       icon_ms_dihedral, icon_ax_scene, icon_ax_window, icon_at_paint, icon_at_enlarge,&
       icon_at_shrink, icon_at_hide, icon_ui_polyhedra, icon_ui_labels, icon_at_shift
    integer, intent(in) :: itype
    integer, allocatable, intent(out) :: tools(:)
    integer, allocatable, intent(out) :: icons(:)
    character(len=2), allocatable, intent(out) :: falls(:)

    integer :: k

    if (itype == reptype_atoms) then
       ! the atom tools (atomtool_*), with no select or remove
       tools = (/(k, k = objtool_kind0+1, objtool_kind0+atomtool_NUM)/)
       icons = (/icon_at_paint,icon_at_enlarge,icon_at_shrink,icon_at_hide,icon_ui_polyhedra,&
          icon_ui_labels,icon_at_shift/)
       falls = (/"Pa","En","Sh","Hi","Po","La","Lv"/)
    elseif (itype == reptype_planar) then
       ! the kinds are planarkind_*
       tools = (/(k, k = objtool_select, objtool_kind0+planarkind_NUM)/)
       icons = (/icon_pl_select,icon_vm_remove,icon_pl_ellipse,icon_pl_rect,icon_pl_polygon,&
          icon_pl_polyline,icon_pl_arrow,icon_pl_curve,icon_pl_freehand/)
       falls = (/"Se","Rm","El","Re","Pg","Pl","Ar","Cu","Fh"/)
    elseif (itype == reptype_shapes) then
       ! the kinds are shapekind_*
       tools = (/(k, k = objtool_select, objtool_kind0+shapekind_NUM)/)
       icons = (/icon_pl_select,icon_vm_remove,icon_sh_sphere,icon_ui_cell,icon_pl_arrow,&
          icon_sh_cone,icon_sh_cylinder/)
       falls = (/"Se","Rm","Sp","Bx","Ar","Co","Cy"/)
    elseif (itype == reptype_text) then
       ! the tool of placement ipl is objtool_kind0 + ipl + 1
       tools = (/(k, k = objtool_select, objtool_kind0+textpos_bond+1)/)
       icons = (/icon_pl_select,icon_vm_remove,icon_tx_screen,icon_tx_point,icon_tx_atom,&
          icon_tx_bond/)
       falls = (/"Se","Rm","Sc","Pt","At","Bd"/)
    elseif (itype == reptype_measure) then
       ! the tool for n atoms is objtool_kind0 + n - 1
       tools = (/(k, k = objtool_select, objtool_kind0+3)/)
       icons = (/icon_pl_select,icon_vm_remove,icon_ms_distance,icon_ms_angle,icon_ms_dihedral/)
       falls = (/"Se","Rm","Di","An","Dh"/)
    elseif (itype == reptype_axes) then
       ! select, and the two placements (objtool_kind0 + 1 + axplace_*)
       tools = (/objtool_select,objtool_kind0+1+axplace_scene,objtool_kind0+1+axplace_window/)
       icons = (/icon_pl_select,icon_ax_scene,icon_ax_window/)
       falls = (/"Se","Sc","Wi"/)
    else
       ! select and remove, for all the object types (the annotation row)
       tools = (/objtool_select,objtool_remove/)
       icons = (/icon_pl_select,icon_vm_remove/)
       falls = (/"Se","Rm"/)
    end if

  end subroutine objtool_list

  !> Whether tool itool (objtool_kind0 + a kind) of the object type
  !> itype acts on a click alone, with no meaning for a drag: the atom
  !> tools except paint, and the tools that place texts, measure, or
  !> place the axes. A drag with them rotates the camera, as in
  !> navigation.
  module function objtool_isclick(itype,itool) result(ok)
    integer, intent(in) :: itype, itool
    logical :: ok

    ok = .false.
    if (itool <= objtool_kind0) return
    if (itype == reptype_atoms) then
       ok = (itool /= objtool_kind0 + atomtool_paint)
    else
       ok = (itype == reptype_text .or. itype == reptype_measure .or. itype == reptype_axes)
    end if

  end function objtool_isclick

  !> The tooltip for tool itool (objtool_*, or objtool_kind0 + a kind) of
  !> the object type itype in the toolbars; itype = 0 for the select and
  !> remove tools of the annotation row, which work on all the types.
  module function objtool_hint(itype,itool) result(str)
    use keybindings, only: BIND_OBJEDIT_DRAW, BIND_OBJEDIT_DELETE
    integer, intent(in) :: itype, itool
    character(len=:), allocatable :: str

    if (itype == reptype_atoms) then
       str = atom_tool_hint(itool)
    elseif (itype == reptype_planar) then
       str = planar_tool_hint(itool)
    elseif (itype == reptype_shapes) then
       str = shapes_tool_hint(itool)
    elseif (itype == reptype_text) then
       str = text_tool_hint(itool)
    elseif (itype == reptype_measure) then
       str = measure_tool_hint(itool)
    elseif (itype == reptype_axes) then
       str = axes_tool_hint(itool)
    elseif (itool == objtool_select) then
       str = "Select: click (" // kn(BIND_OBJEDIT_DRAW) // ") a 2D or 3D shape, a text, a " //&
          "measurement, or the axes to select it, and drag it or its handles to edit it; " //&
          "double-click it to open its object editor. " // kn(BIND_OBJEDIT_DELETE) //&
          " removes the selected item. " // exit_hint()
    elseif (itool == objtool_remove) then
       str = "Remove: click (" // kn(BIND_OBJEDIT_DRAW) // ") a 2D or 3D shape, a text, or a " //&
          "measurement to remove it. " // exit_hint()
    else
       str = ""
    end if

  end function objtool_hint

  !> The prompt for tool itool (objtool_*, or objtool_kind0 + a kind) of
  !> the object type itype in the view bar; itype = 0 for the select and
  !> remove tools of the annotation row, which work on all the types.
  module function objtool_prompt(itype,itool) result(str)
    integer, intent(in) :: itype, itool
    character(len=:), allocatable :: str

    if (itype == reptype_atoms) then
       str = atom_tool_prompt(itool)
    elseif (itype == reptype_planar) then
       str = planar_tool_prompt(itool)
    elseif (itype == reptype_shapes) then
       str = shapes_tool_prompt(itool)
    elseif (itype == reptype_text) then
       str = text_tool_prompt(itool)
    elseif (itype == reptype_measure) then
       str = measure_tool_prompt(itool)
    elseif (itype == reptype_axes) then
       str = axes_tool_prompt(itool)
    elseif (itool == objtool_select) then
       str = "Click an item to select it; drag to edit it"
    elseif (itool == objtool_remove) then
       str = "Click an item to remove it"
    else
       str = ""
    end if

  end function objtool_prompt

  !> Draw the editrep window, isosurface class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_isosurface(w,ttshown) result(changed)
    use systems, only: sys
    use representations, only: iso_isgridfield, iso_level_custom,&
       iso_npts_custom_min, iso_npts_custom_max, iso_level_optstr_custom, iso_defaultlevel,&
       iso_level_tipstr_custom, iso_custom_mode_tipstr,&
       iso_region_cell, iso_region_frac, iso_region_ortho, iso_region_parallel,&
       iso_region_simplebox, iso_region_cube, iso_region_bbox, iso_region_name,&
       iso_knd_cry, iso_knd_mol, iso_region_modes_cry, iso_region_desc, iso_nhist,&
       iso_hscale_name,&
       iso_hscale_y, iso_map_color, iso_map_field, iso_map_expr, iso_map_optstr, iso_explen,&
       iso_region_modes_mol, iso_region_to_box, iso_region_seed,&
       iso_region_point_from_cart, iso_level_ptsang,&
       iso_custom_mode_optstr, iso_custom_ptsang, iso_ptsang_min, iso_ptsang_max,&
       iso_ptsang_from_npts, rep_shape, shapekind_box
    use grid3mod, only: hscale_num, hscale_log, hscale_asinh
    use utils, only: iw_table_headers_row, iw_text, iw_tooltip, iw_coloredit, iw_dragfloat_real8, iw_textwidth,&
       iw_calcwidth, iw_calcheight, iw_combo_simple, iw_button, iw_intstepper, iw_checkbox,&
       iw_close_button, iw_table_column, iw_highlight_selectable,&
       iw_field_combo, iw_cmap_optstr, iw_ncmap, iw_colormap_lut, iw_arith_help,&
       iw_arith_help_button, iw_helpermark, iw_inputtext, igIsItemHovered_delayed,&
       duration_string
    use gui_main, only: fontsize, g, tooltip_enabled, tooltip_delay
    use param, only: newline
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: i, isys, iview, nstage(3), nshow(3), nrow, nmode, istat, ilevprev, ipad, ihb, iheld
    integer :: rmodes(size(iso_region_modes_mol)), iknd, nsc, isc, scmap(hscale_num), idel, iline, imode, ifield, ntick
    integer(c_int) :: ncus(3), tflags, dtflags

    real(c_float), parameter :: hist_lwidth = 4.5_c_float ! histogram curve width (px)
    real(c_float), parameter :: hist_lwdrag = 5.5_c_float ! isovalue line width (px)
    real*8 :: box(3,0:3), prev0(3), prevv(3,3), flo(3), fhi(3), xpick(3), xlo, xhi, xeps
    real*8 :: xnew, dx, dmin, alpha8, maprspeed, rlo, rhi, reps
    logical :: ch, ch2, ldum, goodf, isgrid, navail, okbox, capped, lapply, haverange
    logical :: ismapped, hascolors, plotted
    logical(c_bool) :: is_selected
    real(c_float) :: rgba(4)
    real*8 :: speed
    real(c_double) :: xdrag
    real(c_float) :: pxline, pyline
    real(c_double), target :: hx(2*iso_nhist), hy(2*iso_nhist)
    character(kind=c_char,len=:), allocatable, target :: str1, str2
    character(len=:), allocatable :: str3
    type(ImVec2) :: szero, szplot, sztab, mpos, szavail
    real(c_float) :: plotw
    real(c_double), target :: tickv(32)
    type(ImVec4) :: colline, colhist

    ! initialize
    changed = .false.
    isys = w%isys
    iview = w%anchor_view()
    szero%x = 0
    szero%y = 0

    ! handle a pending region-point pick commanded to the anchor view
    if (w%editrep_isopick >= 0) then
       if (w%rep%iso%iregion /= w%editrep_isopick_mode .or.&
          .not.sys(isys)%goodfield(w%rep%iso%ifield)) then
          call win(iview)%viewmode_release_forced(w%id)
          w%editrep_isopick = -1
       else
          call view_pick_result(iview,w%id,isys,w%editrep_pick,istat,xpick)
          if (istat == ipick_point) &
             call iso_region_point_from_cart(isys,w%rep%iso%iregion,xpick,&
                w%rep%iso%rgn_x(:,w%editrep_isopick))
          if (istat /= ipick_pending) w%editrep_isopick = -1
       end if
    end if

    ! field selector
    call iw_text("Field",highlight=.true.)
    goodf = sys(isys)%goodfield(w%rep%iso%ifield)
    call igSameLine(0._c_float,-1._c_float)
    ifield = w%rep%iso%ifield
    if (iw_field_combo("##isofieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
       nonestr="<field not available>")) then
       call w%rep%iso%set_field(isys,ifield)
       w%editrep_isoline = 0
       changed = .true.
    end if
    call iw_tooltip("Field whose isosurface is displayed",ttshown)
    if (.not.goodf) then
       call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       w%editrep_isoline = 0
       return
    end if

    ! read the field: it may have stopped being a grid (e.g. a
    ! different field reloaded into the same slot), which invalidates
    ! the native level
    isgrid = iso_isgridfield(isys,w%rep%iso%ifield)
    if (.not.isgrid .and. w%rep%iso%ilevel == 0) then
       call w%rep%iso%set_field(isys,w%rep%iso%ifield)
       w%editrep_isoline = 0
       changed = .true.
    end if

    ! region section
    call iw_text("Region",highlight=.true.)
    if (sys(isys)%c%ismolecule) then
       iknd = iso_knd_mol
       nmode = size(iso_region_modes_mol)
       rmodes(1:nmode) = iso_region_modes_mol
    else
       iknd = iso_knd_cry
       nmode = size(iso_region_modes_cry)
       rmodes(1:nmode) = iso_region_modes_cry
    end if
    ! reset a staged mode this system kind does not offer
    if (all(rmodes(1:nmode) /= w%rep%iso%iregion)) then
       w%rep%iso%iregion = rmodes(1)
       call iso_region_seed(isys,w%rep%iso%iregion,w%rep%iso%rgn_x)
    end if
    call iw_text("Mode",alignframe=.true.)
    str1 = "##isoregion" // c_null_char
    str2 = trim(iso_region_name(w%rep%iso%iregion,iknd)) // c_null_char
    call igSameLine(0._c_float,-1._c_float)
    call igSetNextItemWidth(iw_calcwidth(len(iso_region_name)+4,0))
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nmode
          is_selected = (w%rep%iso%iregion == rmodes(i))
          str2 = trim(iso_region_name(rmodes(i),iknd)) // c_null_char
          if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) then
             if (w%rep%iso%iregion /= rmodes(i)) then
                w%rep%iso%iregion = rmodes(i)
                call iso_region_seed(isys,w%rep%iso%iregion,w%rep%iso%rgn_x)
             end if
          end if
          if (is_selected) &
             call igSetItemDefaultFocus()
       end do
       call igEndCombo()
    end if
    ! the modes explained in a tooltip, one per line
    if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
       if (tooltip_enabled .and. igIsMouseHoveringRect(g%LastItemData%NavRect%min,&
          g%LastItemData%NavRect%max,.false._c_bool)) then
          call igBeginTooltip()
          call iw_text("Region over which the isosurface is calculated:")
          call iw_text("")
          do i = 1, nmode
             call iw_text(trim(iso_region_name(rmodes(i),iknd)) // " = " //&
                trim(iso_region_desc(rmodes(i))),wrap=.true.)
          end do
          call igEndTooltip()
       end if
    end if

    ! Region coordinates
    if (w%rep%iso%iregion == iso_region_bbox) then
       call iw_text("Buffer (Å)",alignframe=.true.)
       ldum = iw_dragfloat_real8("##isorgnbuf",x1=w%rep%iso%rgn_x(1,1),sameline=.true.,&
          speed=0.05d0,decimal=3,min=0d0,flags=ImGuiSliderFlags_AlwaysClamp)
       w%rep%iso%rgn_x(2:3,1) = w%rep%iso%rgn_x(1,1)
       call iw_tooltip("Distance added around the bounding box of the atoms, in&
          & every direction",ttshown)
    elseif (w%rep%iso%iregion == iso_region_simplebox .or.&
       w%rep%iso%iregion == iso_region_cube) then
       ! both labels padded to the width of the longest, so the drags line up
       ! and do not shift when the mode changes
       if (w%rep%iso%iregion == iso_region_cube) then
          str3 = "Half-length (Å)"
       else
          str3 = "Half-lengths (Å)"
       end if
       call iw_text(string("x0 (Å)",length=len("Half-lengths (Å)")),alignframe=.true.)
       ldum = iw_dragfloat_real8("##isorgn0",x3=w%rep%iso%rgn_x(:,0),sameline=.true.,&
          speed=0.05d0,decimal=3)
       call iw_tooltip(region_point_tooltip(w%rep%iso%iregion,0),ttshown)
       call region_pick_button(0)
       call iw_text(string(str3,length=len("Half-lengths (Å)")),alignframe=.true.)
       if (w%rep%iso%iregion == iso_region_cube) then
          ldum = iw_dragfloat_real8("##isorgn1",x1=w%rep%iso%rgn_x(1,1),sameline=.true.,&
             speed=0.05d0,decimal=3,min=0d0,flags=ImGuiSliderFlags_AlwaysClamp)
          w%rep%iso%rgn_x(2:3,1) = w%rep%iso%rgn_x(1,1)
       else
          ldum = iw_dragfloat_real8("##isorgn1",x3=w%rep%iso%rgn_x(:,1),sameline=.true.,&
             speed=0.05d0,decimal=3,min=0d0,flags=ImGuiSliderFlags_AlwaysClamp)
       end if
       call iw_tooltip(region_point_tooltip(w%rep%iso%iregion,1),ttshown)
    elseif (w%rep%iso%iregion /= iso_region_cell) then
       nrow = 1 ! last drag row: the far corner, or the third edge endpoint
       if (w%rep%iso%iregion == iso_region_parallel) nrow = 3
       do i = 0, nrow
          str1 = "##isorgn" // string(i)
          if (w%rep%iso%iregion == iso_region_frac) then
             call iw_text("x" // string(i),alignframe=.true.)
             ldum = iw_dragfloat_real8(str1,x3=w%rep%iso%rgn_x(:,i),sameline=.true.,&
                speed=0.01d0,decimal=4)
          else
             call iw_text("x" // string(i) // " (Å)",alignframe=.true.)
             ldum = iw_dragfloat_real8(str1,x3=w%rep%iso%rgn_x(:,i),sameline=.true.,&
                speed=0.05d0,decimal=3)
          end if
          call iw_tooltip(region_point_tooltip(w%rep%iso%iregion,i),ttshown)
          call region_pick_button(i)
       end do
    end if

    ! the region box
    call iso_region_to_box(isys,w%rep%iso%iregion,w%rep%iso%rgn_x,box,okbox)
    if (.not.okbox) &
       call iw_text("The region is degenerate (zero volume)",danger=.true.,wrap=.true.)

    ! grid section
    call iw_text("Grid",highlight=.true.)

    ! coarseness level combo ("Native grid" only for grid fields over the
    ! whole cell)
    navail = isgrid .and. (w%rep%iso%iregion == iso_region_cell)
    if (.not.navail .and. w%rep%iso%ilevel == 0) &
       w%rep%iso%ilevel = iso_defaultlevel
    if (navail) then
       str1 = "Native grid" // c_null_char // iso_level_optstr_custom
       str2 = "The grid of the field itself" // c_null_char // iso_level_tipstr_custom
    else
       str1 = iso_level_optstr_custom
       str2 = iso_level_tipstr_custom
    end if
    ilevprev = w%rep%iso%ilevel
    call iw_text("Level",alignframe=.true.)
    call iw_combo_simple("##isolevel",str1,w%rep%iso%ilevel,startsatone=.not.navail,&
       sameline=.true.,changed=ch,tooltips=str2,ttshown=ttshown)

    ! custom level: the grid given by the user, either as the number of
    ! points along each axis or as a resolution in points per angstrom
    if (w%rep%iso%ilevel == iso_level_custom) then
       ! arriving at the custom level: start from the grid that was on
       ! screen, written both ways, so either mode continues it. The
       ! resolution is the one those dimensions amount to, not the
       ! level's tabulated value: a named level is floored per axis and
       ! capped in total, so the two need not describe the same grid
       if (ch .and. ilevprev /= iso_level_custom) then
          if (ilevprev == 0) then
             w%rep%iso%nptscustom = sys(isys)%f(w%rep%iso%ifield)%grid%n
          elseif (okbox) then
             w%rep%iso%nptscustom = w%rep%iso%staged_grid_size(isys,ilevprev,box)
          end if
          ! only where the dimensions were re-seeded: a degenerate region
          ! leaves both alone, rather than deriving one from the other's
          ! stale value over a zero-volume box
          if (okbox .or. ilevprev == 0) &
             w%rep%iso%ptsangcustom = npts_to_ptsang(w%rep%iso%nptscustom)
       end if

       call iw_text("Given as",alignframe=.true.)
       call iw_combo_simple("##isocustom",iso_custom_mode_optstr,w%rep%iso%icustom,&
          sameline=.true.,changed=ch2,tooltips=iso_custom_mode_tipstr,ttshown=ttshown)
       ! switching between the two carries the grid over, so the mode
       ! change by itself does not move the isosurface
       if (ch2 .and. okbox) then
          if (w%rep%iso%icustom == iso_custom_ptsang) then
             w%rep%iso%ptsangcustom = npts_to_ptsang(w%rep%iso%nptscustom)
          else
             ! the dimensions the resolution on screen gives, which the
             ! staged options no longer decode to now that the mode changed
             w%rep%iso%nptscustom = w%rep%iso%staged_grid_size(isys,iso_level_custom,box,&
                icustom=iso_custom_ptsang)
          end if
       end if

       if (w%rep%iso%icustom == iso_custom_ptsang) then
          call iw_text("Resolution (pts/Å)",alignframe=.true.)
          ldum = iw_dragfloat_real8("##isoptsang",x1=w%rep%iso%ptsangcustom,sameline=.true.,&
             speed=0.05d0,decimal=2,min=iso_ptsang_min,max=iso_ptsang_max,&
             flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Number of grid points per angstrom along each axis of the grid",ttshown)
       else
          ncus = w%rep%iso%nptscustom
          ipad = ceiling(log10(max(maxval(ncus),1) + 0.1))
          ldum = iw_intstepper("n1##isonpts1",ncus(1),label="n1:",minval=int(iso_npts_custom_min,c_int),&
             maxval=int(iso_npts_custom_max,c_int),ndigit=ipad,notlive=.true.)
          ldum = iw_intstepper("n2##isonpts2",ncus(2),label="n2:",minval=int(iso_npts_custom_min,c_int),&
             maxval=int(iso_npts_custom_max,c_int),ndigit=ipad,notlive=.true.,sameline=.true.)
          ldum = iw_intstepper("n3##isonpts3",ncus(3),label="n3:",minval=int(iso_npts_custom_min,c_int),&
             maxval=int(iso_npts_custom_max,c_int),ndigit=ipad,notlive=.true.,sameline=.true.)
          call iw_tooltip("Number of grid points along each axis of the grid",ttshown)
          w%rep%iso%nptscustom = ncus
       end if
    end if

    ! grid dimensions
    lapply = .false.
    nstage = 0
    if (okbox) then
       nstage = w%rep%iso%staged_grid_size(isys,w%rep%iso%ilevel,box,capped=capped)
       nshow = nstage
       str2 = ""
       if (all(nstage == 0)) then
          nshow = sys(isys)%f(w%rep%iso%ifield)%grid%n
          str2 = " (native)"
       elseif (capped) then
          str2 = " (capped)"
       end if
       if (w%rep%iso%ilevel /= iso_level_custom .or.&
          w%rep%iso%icustom == iso_custom_ptsang .or. capped) then
          call iw_text("Grid: " // string(nshow(1)) // " x " // string(nshow(2)) // " x " //&
             string(nshow(3)) // str2)
          call iw_tooltip("Dimensions of the grid that will support the isosurface",ttshown)
       end if
       lapply = .not.w%rep%iso%grid_isapplied(nstage,w%rep%iso%iregion,w%rep%iso%rgn_x)
    end if
    if (iw_button("Calculate grid",danger=.true.,disabled=.not.lapply)) then
       w%errmsg = ""
       if (all(nstage == 0)) then
          ! the native grid of a grid field: nothing to sample
          call w%rep%iso%apply_grid(nstage,w%rep%iso%iregion,w%rep%iso%rgn_x)
          changed = .true.
       elseif (w%request_block()) then
          ! the field is sampled in a blocking job (run_editrep)
          w%editrep_pending_n = nstage
       else
          w%errmsg = "Another calculation is starting: press the button again"
       end if
    end if
    call iw_tooltip("Use the selected grid and region for the isosurface. The options&
       & above have no effect on the isosurface until this button is pressed",ttshown)
    if (iw_button("Estimate cost",sameline=.true.,disabled=all(nstage == 0))) &
       w%rep%iso%costest = w%rep%iso%measure_cost(isys,nstage)
    call iw_tooltip("Estimate the time needed to sample the field with the selected options",ttshown)
    if (w%rep%iso%costest >= 0d0 .and. any(nstage /= 0)) &
       call iw_text("~" // duration_string(w%rep%iso%costest * product(real(nstage,8))),&
          sameline=.true.)

    ! the state of the options, on the same line as the button that fixes it
    if (okbox .and. .not.w%rep%iso%isgenerated(isys)) then
       call iw_text("not generated",danger=.true.,sameline=.true.)
    elseif (lapply) then ! false for a degenerate region
       call iw_text("settings changed",danger=.true.,sameline=.true.)
    end if
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)
    if (w%rep%iso%outdomain) &
       call iw_text("The field could not be evaluated in part of the region",&
          danger=.true.,wrap=.true.)

    ! show the staged region in the view for as long as this editor is
    ! open
    if (okbox .and. win(iview)%sc%isinit /= 0 .and.&
       .not.(w%rep%iso%iregion == iso_region_cell .and. .not.sys(isys)%c%ismolecule)) then
       if (w%rep%iso%iregion == iso_region_cell .and. isgrid) then
          call sys(isys)%f(w%rep%iso%ifield)%grid%get_domain(prevv,prev0,flo=flo,fhi=fhi)
          prev0 = prev0 + matmul(prevv,flo)
          do i = 1, 3
             prevv(:,i) = prevv(:,i) * (fhi(i) - flo(i))
          end do
       else
          prev0 = box(:,0)
          prevv = box(:,1:3)
       end if
       if (sys(isys)%c%ismolecule) prev0 = prev0 + sys(isys)%c%molx0
       call win(iview)%sc%show_transient_shapes(w%id,1,(/rep_shape(kind=shapekind_box,&
          x1=prev0,v=prevv,rad=region_edgerad,rgb=region_rgb,alpha=region_alpha)/))
    end if

    ! isosurface
    call iw_text("Isosurface(s)",highlight=.true.)
    call iw_helpermark("One row per isosurface: whether it is drawn, where its color&
       & comes from, the value of the field on it, and what it encloses -- the fraction&
       & of the integral of |f| and the fraction of the volume, both over the region&
       & the isosurface was calculated on. Click a row&
       & to edit the colors of a mapped isosurface underneath the table.",sameline=.true.)

    ! histogram of the values backing the mesh, with one draggable line
    ! per isosurface
    plotted = .false.
    if (w%rep%iso%nhist > 0) then
       ! the width the value axis will get, for the tick spacing below
       call igGetContentRegionAvail(szavail)
       plotw = szavail%x - iw_calcwidth(8,0)
       str1 = "##isohistogram" // c_null_char
       szplot%x = -1._c_float
       szplot%y = iw_calcheight(1,5,.false.)
       if (ipBeginPlot(c_loc(str1),szplot,ior(ior(ImPlotFlags_NoTitle,ImPlotFlags_NoLegend),&
          ior(ImPlotFlags_NoInputs,ImPlotFlags_NoMouseText)))) then
          plotted = .true.
          ! no tick marks
          str1 = "Field value" // c_null_char
          str2 = "Points" // c_null_char
          call ipSetupAxes(c_loc(str1),c_loc(str2),ImPlotAxisFlags_NoTickMarks,&
             ImPlotAxisFlags_NoTickMarks)
          call ipSetupAxisScale(ImAxis_X1,implot_scale(w%rep%iso%hist_xscale))
          call ipSetupAxisScale(ImAxis_Y1,implot_scale(w%rep%iso%hist_yscale))

          ! the staircase binned to match the value axis
          hx = w%rep%iso%hist_x(:,w%rep%iso%hist_xscale)
          hy = w%rep%iso%hist_y(:,w%rep%iso%hist_xscale)

          ! explicit limits, always locked
          xlo = w%rep%iso%frange(1)
          xhi = w%rep%iso%frange(2)
          if (w%rep%iso%hist_xscale == hscale_asinh .and. xlo*xhi >= 0d0) then
             xeps = 1d-6 * max(abs(xlo),abs(xhi))
             xlo = min(xlo,-xeps)
             xhi = max(xhi,xeps)
          end if
          call ipSetupAxisLimits(ImAxis_X1,real(xlo,c_double),real(xhi,c_double),&
             ImPlotCond_Always)

          ! decade ticks, spaced so their labels do not run into each other.
          if (w%rep%iso%hist_xscale == hscale_log .and. xlo > 0d0) then
             ntick = log_ticks(xlo,xhi,plotw,tickv)
             if (ntick >= 2) &
                call ipSetupAxisTicksValues(ImAxis_X1,c_loc(tickv),int(ntick,c_int),.false._c_bool)
          end if
          if (w%rep%iso%hist_yscale == hscale_log) then
             call ipSetupAxisLimits(ImAxis_Y1,0.5_c_double,real(maxval(hy)*2d0,c_double),&
                ImPlotCond_Always)
          else
             call ipSetupAxisLimits(ImAxis_Y1,0._c_double,real(maxval(hy)*1.05d0,c_double),&
                ImPlotCond_Always)
          end if

          ! outline only
          str1 = "Distribution" // c_null_char
          colhist = ImVec4(0.78_c_float,0.78_c_float,0.78_c_float,1._c_float)
          call ipSetNextLineStyle(colhist,hist_lwidth)
          call ipPlotLine(c_loc(str1),c_loc(hx),c_loc(hy),&
             int(w%rep%iso%nhist,c_int),ImPlotLineFlags_None,0_c_int)

          ! one draggable line per isosurface, in its own colour
          iheld = w%editrep_isoline ! the line held in the previous frame (0 = none)
          iline = iheld
          if (iline < 1 .or. iline > w%rep%iso%niso) then
             iline = 1
             call igGetMousePos(mpos)
             dmin = huge(dmin)
             do i = 1, w%rep%iso%niso
                call ipPlotToPixels(real(w%rep%iso%slot(i)%isoval,c_double),1._c_double,&
                   ImAxis_X1,ImAxis_Y1,pxline,pyline)
                dx = abs(real(pxline - mpos%x,8))
                if (dx < dmin) then
                   dmin = dx
                   iline = i
                end if
             end do
          end if
          w%editrep_isoline = 0
          do i = 1, w%rep%iso%niso
             associate (s => w%rep%iso%slot(i))
               xdrag = real(s%isoval,c_double)
               colline = ImVec4(s%rgb(1),s%rgb(2),s%rgb(3),1._c_float)
               dtflags = ImPlotDragToolFlags_None
               if (i /= iline) dtflags = ImPlotDragToolFlags_NoInputs
               if (logical(ipDragLineX(int(i-1,c_int),xdrag,colline,hist_lwdrag,dtflags))) then
                  ! keep the grab on this line for as long as it is held,
                  ! wherever the cursor goes
                  w%editrep_isoline = i
                  xnew = max(min(real(xdrag,8),w%rep%iso%frange(2)),w%rep%iso%frange(1))
                  ! a held but motionless line changes nothing: rebuilding
                  ! the scene for it would be pure waste
                  if (xnew /= s%isoval) then
                     s%isoval = xnew
                     ! a flat isosurface is only re-triangulated, and follows
                     ! the line as it slides. A mapped one also pays a map
                     ! evaluation at every vertex of the new mesh, which is too
                     ! much for every frame of a drag: it waits for the release
                     if (s%imap_mode == iso_map_color) changed = .true.
                  end if
               end if
             end associate
          end do
          ! the level of a mapped isosurface lands when its line is let go
          if (iheld >= 1 .and. w%editrep_isoline == 0) then
             if (iheld <= w%rep%iso%niso) then
                if (w%rep%iso%slot(iheld)%imap_mode /= iso_map_color) changed = .true.
             end if
          end if
          call ipEndPlot()
       end if
       call iw_tooltip("Distribution of the values of the field over the grid data.&
          & Drag the vertical line to change the isovalue of the isosurface drawn in that color.",ttshown)
    end if
    if (.not.plotted .and. w%editrep_isoline /= 0) then
       ! the histogram did not draw this frame: treat the held line as
       ! released so a mapped isosurface commits its pending rebuild
       if (w%editrep_isoline >= 1 .and. w%editrep_isoline <= w%rep%iso%niso) then
          if (w%rep%iso%slot(w%editrep_isoline)%imap_mode /= iso_map_color) changed = .true.
       end if
       w%editrep_isoline = 0
    end if
    if (w%rep%iso%nhist > 0) then

       ! axis scales
       call iw_text("Scale x:",alignframe=.true.)
       nsc = 0
       do i = 1, hscale_num
          if (.not.w%rep%iso%hist_have(i)) cycle
          nsc = nsc + 1
          scmap(nsc) = i
          if (i == w%rep%iso%hist_xscale) isc = nsc
       end do
       str1 = ""
       do i = 1, nsc
          str1 = str1 // trim(iso_hscale_name(scmap(i))) // c_null_char
       end do
       call iw_combo_simple("##isohistxscale",str1,isc,sameline=.true.,changed=ch,startsatone=.true.)
       if (ch) w%rep%iso%hist_xscale = scmap(max(min(isc,nsc),1))
       call iw_tooltip("Scale of the value axis",ttshown)

       call iw_text("y:",sameline=.true.,alignframe=.true.)
       isc = 1
       str1 = ""
       do i = 1, size(iso_hscale_y)
          str1 = str1 // trim(iso_hscale_name(iso_hscale_y(i))) // c_null_char
          if (iso_hscale_y(i) == w%rep%iso%hist_yscale) isc = i
       end do
       call iw_combo_simple("##isohistyscale",str1,isc,sameline=.true.,changed=ch,startsatone=.true.)
       if (ch) w%rep%iso%hist_yscale = iso_hscale_y(max(min(isc,size(iso_hscale_y)),1))
       call iw_tooltip("Scale of the count axis",ttshown)
    else
       call iw_text("The distribution of the values of the field appears here once the&
          & grid has been calculated",wrap=.true.)
    end if

    ! the isosurfaces of this object: one row each
    idel = 0
    tflags = ImGuiTableFlags_None
    tflags = ior(tflags,ImGuiTableFlags_RowBg)
    tflags = ior(tflags,ImGuiTableFlags_Borders)
    tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
    str1 = "##isotable" // c_null_char
    ! the table grows with the list
    sztab%x = 0
    sztab%y = 0
    ! the range prefix of the level tooltips does not depend on the row
    haverange = (w%rep%iso%frange(1) <= w%rep%iso%frange(2))
    if (haverange) &
       str3 = "Value of the field on this isosurface, in atomic units (data range: " //&
       string(w%rep%iso%frange(1),'e',decimal=4) // " to " //&
       string(w%rep%iso%frange(2),'e',decimal=4) // ")"
    if (igBeginTable(c_loc(str1),5,tflags,sztab,0._c_float)) then
       ! the row-control column keeps the width of its widgets even when it is empty
       call iw_table_column("",id=0,flags=ImGuiTableColumnFlags_WidthFixed,&
          width=max(4._c_float,iw_calcwidth(1,1,ncheck=1)))
       call iw_table_column("Mode",id=1,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Color/map",id=2,flags=ImGuiTableColumnFlags_WidthFixed,&
          width=iw_calcwidth(12,1))
       ! the isovalue drag sizes itself from the integer digits of the data
       ! range, so it takes the slack; the readout beside it is a fixed width
       call iw_table_column("Isovalue",id=3,flags=ImGuiTableColumnFlags_WidthStretch)
       call iw_table_column("Encloses",id=4,flags=ImGuiTableColumnFlags_WidthFixed,&
          width=iw_calcwidth(13,0))
       call iw_table_headers_row()

       do i = 1, w%rep%iso%niso
          associate (s => w%rep%iso%slot(i))
            call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
            ! draw this isosurface, and remove it (never the last one)
            if (igTableSetColumnIndex(0)) then
               if (iw_checkbox("##isoshown" // string(i),s%shown)) changed = .true.
               call iw_tooltip("Draw this isosurface in the view",ttshown)
               if (w%rep%iso%niso > 1) then
                  call igSameLine(0._c_float,-1._c_float)
                  if (iw_close_button("##isodel" // string(i))) idel = i
                  call iw_tooltip("Remove this isosurface",ttshown)
               end if
            end if
            ! where the color of this isosurface comes from
            if (igTableSetColumnIndex(1)) then
               imode = s%imap_mode
               call iw_combo_simple("##isomapmode" // string(i),iso_map_optstr,imode,changed=ch)
               if (ch .and. imode /= s%imap_mode) then
                  s%imap_mode = imode
                  ! first time this isosurface is mapped: start from the
                  ! field it is an isosurface of, which always exists
                  if (s%imap_mode == iso_map_field .and. s%imap < 0) s%imap = w%rep%iso%ifield
                  w%rep%iso%isel = i
                  changed = .true.
               end if
               call iw_tooltip("Whether this isosurface has a single color or takes it&
                  & from the values of another field or expression",ttshown)
            end if
            ! the color itself
            if (igTableSetColumnIndex(2)) then
               if (s%imap_mode == iso_map_field) then
                  if (iw_field_combo("##isomapfield" // string(i),isys,s%imap)) then
                     s%icmap_auto = .true.
                     s%maprange_auto = .true.
                     w%rep%iso%isel = i
                     changed = .true.
                  end if
                  call iw_tooltip("Field whose values color this isosurface",ttshown)
               elseif (s%imap_mode == iso_map_expr) then
                  call igGetContentRegionAvail(szavail)
                  call igSetNextItemWidth(szavail%x)
                  if (iw_inputtext("##isomapexpr" // string(i),bufsize=iso_explen,&
                     textf=s%mapexpr,notlive=.true.)) then
                     ! a new expression is a new quantity, like a new field
                     s%maperr = ""
                     s%icmap_auto = .true.
                     s%maprange_auto = .true.
                     w%rep%iso%isel = i
                     changed = .true.
                  end if
                  call iw_tooltip("Expression whose values color this isosurface."//newline//&
                     iw_arith_help,ttshown)
               else
                  rgba(1:3) = s%rgb
                  rgba(4) = s%alpha
                  if (iw_coloredit("##isocol" // string(i),rgba=rgba)) then
                     s%rgb = rgba(1:3)
                     s%alpha = rgba(4)
                     changed = .true.
                  end if
                  call iw_tooltip("Color and opacity of this isosurface",ttshown)
               end if
            end if
            ! the level: this widget is notlive, so its re-triangulation
            ! runs on commit rather than on every frame of the drag
            if (igTableSetColumnIndex(3)) then
               speed = max(0.01d0 * abs(s%isoval),1d-4)
               ! same as the histogram line: a mapped isosurface is rebuilt on
               ! the commit (the drag let go, or the typed value entered),
               ! a flat one on every step of the drag
               if (haverange) then
                  ldum = iw_dragfloat_real8("##isoval" // string(i),x1=s%isoval,&
                     speed=speed,decimal=6,notlive=.true.,min=w%rep%iso%frange(1),&
                     max=w%rep%iso%frange(2),flags=ImGuiSliderFlags_AlwaysClamp,committed=ch)
                  if (s%imap_mode /= iso_map_color) ldum = ch
                  changed = changed .or. ldum
                  call iw_tooltip(str3,ttshown)
               else
                  ldum = iw_dragfloat_real8("##isoval" // string(i),x1=s%isoval,&
                     speed=speed,decimal=6,notlive=.true.,committed=ch)
                  if (s%imap_mode /= iso_map_color) ldum = ch
                  changed = changed .or. ldum
                  call iw_tooltip("Value of the field on this isosurface (atomic units)",ttshown)
               end if
            end if
            ! what this isosurface encloses, the number that decides a level
            if (igTableSetColumnIndex(4)) then
               ihb = hist_bin(s%isoval)
               if (ihb >= 1) then
                  call iw_text(string(100d0*w%rep%iso%hist_q(ihb),'f',decimal=1) // "%/" //&
                     string(100d0*w%rep%iso%hist_v(ihb),'f',decimal=1) // "%",alignframe=.true.)
               else
                  call iw_text("--",alignframe=.true.)
               end if
               ! a row-spanning selectable picks the isosurface whose map
               ! options are shown under the table
               ldum = iw_highlight_selectable("##isosel" // string(i),clicked=ch,&
                  selected=(i == w%rep%iso%isel))
               if (ch) w%rep%iso%isel = i
            end if
          end associate
       end do

       ! the last row adds an isosurface, with a level and a color picked
       ! from the ones already there
       call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
       if (igTableSetColumnIndex(1)) then
          if (iw_button("Add##isoadd")) then
             call w%rep%iso%add_iso()
             changed = .true.
          end if
          call iw_tooltip("Add an isosurface to this object",ttshown)
       end if
       call igEndTable()
    end if
    if (idel > 0) then
       call w%rep%iso%del_iso(idel)
       ! the slots above the deleted one moved down
       if (idel < w%rep%iso%isel) w%rep%iso%isel = w%rep%iso%isel - 1
       changed = .true.
    end if

    ! the options of the field map of the selected isosurface, underneath the
    ! table: colormap, opacity, the range the colormap spans and its legend.
    ! Only a mapped isosurface has any (a flat one carries its color and
    ! opacity in the table swatch), but the block is drawn either way, greyed
    ! out: its height must not depend on which row is selected, or picking a
    ! row would shift everything under it. Only an expression adds a line to
    ! it, for the help of what is being written
    if (w%rep%iso%niso < 1) return
    w%rep%iso%isel = max(min(w%rep%iso%isel,w%rep%iso%niso),1)
    associate (s => w%rep%iso%slot(w%rep%iso%isel))
       ! two questions, and they are not the same one. Whether this isosurface
       ! is mapped at all decides what the user may edit -- the colormap, the
       ! opacity (which in map mode lives nowhere else) and the range policy
       ! stay editable while the surface is hidden or not yet calculated.
       ! Whether its colors have actually been computed decides what may be
       ! reported: asking the state, rather than clearing it wherever it could
       ! go stale, is what keeps the legend from describing a surface that is
       ! gone -- a field change, or a grid that has not been calculated yet
       ismapped = (s%imap_mode == iso_map_field)
       if (ismapped) ismapped = sys(isys)%goodfield(s%imap)
       if (s%imap_mode == iso_map_expr) ismapped = (len_trim(s%mapexpr) > 0)
       hascolors = ismapped
       if (hascolors) hascolors = s%built .and. allocated(s%mesh%rgbv)

       ! the help and the error are messages: they stay legible
       if (s%imap_mode == iso_map_expr) then
          call iw_text("Expression",alignframe=.true.)
          call iw_helpermark(iw_arith_help)
          call iw_arith_help_button("##isoexprhelp",ttshown)
       end if
       if (len_trim(s%maperr) > 0) &
          call iw_text("Error: " // trim(s%maperr),danger=.true.,wrap=.true.)

       if (.not.ismapped) call igBeginDisabled(.true._c_bool)
       ! colormap
       call iw_text("Colormap",alignframe=.true.)
       isc = max(min(s%icmap,iw_ncmap),1)
       call iw_combo_simple("##isocmap",iw_cmap_optstr,isc,sameline=.true.,changed=ch,startsatone=.true.)
       if (ch) then
          s%icmap = isc
          s%icmap_auto = .false. ! the user has chosen: stop suggesting
          changed = .true.
       end if
       call iw_tooltip("Colors the values of the map field are drawn with",ttshown)

       ! opacity (in color mode this lives in the color picker)
       alpha8 = real(s%alpha,8)
       call iw_text("Opacity",sameline=.true.,alignframe=.true.)
       if (iw_dragfloat_real8("##isomapalpha",x1=alpha8,speed=0.01d0,min=0d0,max=1d0,&
          decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)) then
          s%alpha = real(alpha8,c_float)
          changed = .true.
       end if
       call iw_tooltip("Opacity of this isosurface (1 = opaque)",ttshown)

       ! the range the colormap spans, and whether it follows the data
       ch = iw_checkbox("Auto range##isomapauto",s%maprange_auto)
       if (ch) then
          ! the range is recomputed by the pass that evaluates the map
          ! field, so switching it back on asks for that pass again
          if (s%maprange_auto) s%imap_built = -1
          changed = .true.
       end if
       call iw_tooltip("Span the colormap over the values the map field takes on this isosurface",ttshown)
       if (.not.s%maprange_auto) then
          maprspeed = max(0.01d0 * maxval(abs(s%maprange)),1d-6)
          if (iw_dragfloat_real8("##isomaprange",x2=s%maprange,speed=maprspeed,decimal=6,&
             sameline=.true.)) changed = .true.
          call iw_tooltip("Values at the two ends of the colormap",ttshown)
       elseif (hascolors) then
          call iw_text("[" // string(s%maprange(1),'e',decimal=4) // ", " //&
             string(s%maprange(2),'e',decimal=4) // "]",sameline=.true.)
       else
          call iw_text("[--, --]",sameline=.true.)
       end if

       ! the legend
       rlo = s%maprange(1)
       rhi = s%maprange(2)
       if (rhi <= rlo) then
          reps = max(1d-6 * abs(rlo),1d-12)
          rlo = rlo - reps
          rhi = rhi + reps
       end if
       call draw_colorbar(s%icmap,rlo,rhi,hascolors)
       if (.not.ismapped) call igEndDisabled()

       ! what the colors could not cover, where there are colors
       if (hascolors .and. s%mapoutdomain) then
          if (s%imap_mode == iso_map_expr) then
             call iw_text("The expression could not be evaluated on part of this&
                & isosurface (shown in grey)",danger=.true.,wrap=.true.)
          else
             call iw_text("Some of this isosurface is outside the domain of the map&
                & field (shown in grey)",danger=.true.,wrap=.true.)
          end if
       end if
    end associate

  contains
    !> A "nice" step (1, 2 or 5 times a power of ten) that divides span
    !> into about ntarget intervals.
    function nice_step(span,ntarget) result(step)
      real*8, intent(in) :: span
      integer, intent(in) :: ntarget
      real*8 :: step

      real*8 :: raw, f
      integer :: ie

      step = 1d0
      if (span <= 0d0) return
      raw = span / real(max(ntarget,1),8)
      ie = floor(log10(raw))
      f = raw / 10d0**ie
      if (f <= 1d0) then
         step = 10d0**ie
      elseif (f <= 2d0) then
         step = 2d0 * 10d0**ie
      elseif (f <= 5d0) then
         step = 5d0 * 10d0**ie
      else
         step = 10d0**(ie+1)
      end if

    end function nice_step

    !> Decade tick values between vlo and vhi (both positive) for a
    !> logarithmic axis wpix pixels long, thinned so that the labels do
    !> not overlap. Returns how many were written into tv.
    function log_ticks(vlo,vhi,wpix,tv) result(nt)
      real*8, intent(in) :: vlo, vhi
      real(c_float), intent(in) :: wpix
      real(c_double), intent(out) :: tv(:)
      integer :: nt

      integer :: e, emin, emax, estep, nmax

      nt = 0
      if (vlo <= 0d0 .or. vhi <= vlo) return
      emin = ceiling(log10(vlo) - 1d-10)
      emax = floor(log10(vhi) + 1d-10)
      if (emax < emin) return

      ! how many labels fit: a decade label is about six characters wide,
      ! and they need a gap of their own between them
      nmax = max(int(wpix / iw_calcwidth(9,0)),2)
      estep = max((emax - emin + nmax - 1) / nmax,1)

      do e = emin, emax, estep
         if (nt >= size(tv)) exit
         nt = nt + 1
         tv(nt) = 10d0**e
      end do
      if (nt < 2) nt = 0

    end function log_ticks

    !> Draw the color legend of a field map: the colormap icmap as a
    !> horizontal bar from vlo (left) to vhi (right), with tick labels at
    !> round values underneath. With filled false the bar is drawn empty:
    !> nothing is mapped, but the legend keeps its height so that picking
    !> a row in the table does not shift what is under it. Drawn here
    !> rather than with implot's ColormapScale, which labels the ticks its
    !> locator places outside the range and then clips them.
    subroutine draw_colorbar(icmap,vlo,vhi,filled)
      integer, intent(in) :: icmap
      real*8, intent(in) :: vlo, vhi
      logical, intent(in) :: filled

      integer, parameter :: nband = 64 ! bands the gradient is drawn with

      integer :: i, ndec, nlab
      integer*8 :: i8
      logical :: usesci
      real*8 :: step, v, span
      real(c_float) :: barw, barh, pad, xa, xb, wl
      real(c_float) :: lut(3,nband)
      type(ImVec2) :: p0, pa, pb, szdum
      type(ImVec4) :: col
      type(c_ptr) :: dl
      character(kind=c_char,len=:), allocatable, target :: strl

      span = vhi - vlo

      ! a wide, short bar: the window is short of height, not of width
      barh = 1.2_c_float * fontsize%y
      pad = 0.3_c_float * fontsize%y
      call igGetContentRegionAvail(szdum)
      barw = max(szdum%x - 2._c_float * fontsize%x,iw_calcwidth(20,0))
      call igGetCursorScreenPos(p0)
      dl = igGetWindowDrawList()

      ! the gradient, in bands from the left (vlo) to the right (vhi), then
      ! labels at round values under the bar: as many as fit without touching,
      ! centered on their value but kept inside the bar, so the ends are not
      ! cut off. One format for the whole scale, so they read as one row: the
      ! values are round multiples of the step, so fixed point needs exactly
      ! the step's decimals; very large or very small ones go scientific
      ! instead of growing a tail of zeros
      if (filled .and. span > 0d0) then
         call iw_colormap_lut(icmap,lut)
         do i = 1, nband
            xa = p0%x + barw * real(i-1,c_float) / real(nband,c_float)
            xb = p0%x + barw * real(i,c_float) / real(nband,c_float)
            pa = ImVec2(xa,p0%y)
            pb = ImVec2(xb + 1._c_float,p0%y + barh)
            col = ImVec4(lut(1,i),lut(2,i),lut(3,i),1._c_float)
            call ImDrawList_AddRectFilled(dl,pa,pb,igGetColorU32_Vec4(col),0._c_float,0_c_int)
         end do

         nlab = max(int(barw / iw_calcwidth(12,0)),1)
         step = nice_step(span,nlab)
         ndec = max(-floor(log10(step)),0)
         usesci = (ndec > 8 .or. max(abs(vlo),abs(vhi)) >= 1d5)
         ! 64-bit bounds: v/step is unbounded for a small user-typed range
         ! at a large offset, even though the tick count itself is bounded
         do i8 = ceiling(vlo/step - 1d-10,8), floor(vhi/step + 1d-10,8)
            v = real(i8,8) * step
            if (v < vlo .or. v > vhi) cycle
            if (usesci) then
               strl = string(v,'e',decimal=2) // c_null_char
            elseif (ndec == 0) then
               strl = string(nint(v)) // c_null_char ! whole numbers, without a trailing dot
            else
               strl = string(v,'f',decimal=ndec) // c_null_char
            end if
            wl = iw_textwidth(strl)
            xa = p0%x + barw * real((v - vlo) / span,c_float) - 0.5_c_float * wl
            xa = min(max(xa,p0%x),p0%x + barw - wl)
            pa = ImVec2(xa,p0%y + barh + pad)
            call ImDrawList_AddText_Vec2(dl,pa,igGetColorU32_Col(ImGuiCol_Text,1._c_float),&
               c_loc(strl),c_null_ptr)
         end do
      end if

      ! the frame goes on top of the bands, and stands alone when there are none
      pa = ImVec2(p0%x,p0%y)
      pb = ImVec2(p0%x + barw,p0%y + barh)
      call ImDrawList_AddRect(dl,pa,pb,igGetColorU32_Col(ImGuiCol_Border,1._c_float),&
         0._c_float,0_c_int,1._c_float)

      ! reserve the space the bar and its labels occupy
      szdum = ImVec2(barw,barh + pad + fontsize%y)
      call igDummy(szdum)

    end subroutine draw_colorbar

    !> ImPlot axis scale for one of the hscale_* ids.
    function implot_scale(is) result(isc_)
      integer, intent(in) :: is
      integer(c_int) :: isc_

      if (is == hscale_log) then
         isc_ = ImPlotScale_Log10
      elseif (is == hscale_asinh) then
         isc_ = ImPlotScale_SymLog
      else
         isc_ = ImPlotScale_Linear
      end if

    end function implot_scale

    !> Index of the bin containing the value x in the binning the
    !> cumulative arrays use: logarithmic in |x| over hist_cumrange.
    !> Returns 0 when there is no histogram or x is zero (a zero
    !> isosurface encloses everything, which the caller reports as
    !> unavailable rather than as 100%).
    function hist_bin(x) result(ib)
      real*8, intent(in) :: x
      integer :: ib

      real*8 :: t, lo, hi

      ib = 0
      if (w%rep%iso%nhist <= 0) return
      if (x == 0d0) return
      lo = w%rep%iso%hist_cumrange(1)
      hi = w%rep%iso%hist_cumrange(2)
      if (hi <= lo) return
      t = (log10(abs(x)) - lo) / (hi - lo)
      ib = min(max(1 + int(t * iso_nhist),1),iso_nhist)

    end function hist_bin

    !> The resolution (points per angstrom) that grid dimensions nn
    !> amount to over the staged region, within the bounds a custom
    !> resolution obeys, so the two ways of writing a custom grid agree
    !> on the total number of points.
    function npts_to_ptsang(nn) result(pa)
      integer, intent(in) :: nn(3)
      real*8 :: pa

      pa = iso_ptsang_from_npts(product(real(nn,8)),&
         w%rep%iso%staged_box_lengths(isys,box))
      if (pa <= 0d0) pa = iso_level_ptsang(iso_defaultlevel)
      pa = min(max(pa,iso_ptsang_min),iso_ptsang_max)

    end function npts_to_ptsang

    !> Draw the Pick button of region coordinate row irow: arms a pick
    !> in the anchor view that fills the row with the clicked position
    !> (an atom, or the unprojected point for empty space); the result
    !> is received at the top of draw_editrep_isosurface, which also
    !> cancels the pick if the region mode stops matching the one
    !> stamped here. One pick can be pending at a time.
    subroutine region_pick_button(irow)
      integer, intent(in) :: irow

      if (iw_button("Pick##isorgnpick" // string(irow),sameline=.true.,&
         disabled=(w%editrep_isopick >= 0))) then
         w%editrep_isopick = irow
         w%editrep_isopick_mode = w%rep%iso%iregion
         call w%editrep_pick%arm()
         call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the region point",w%id,&
            acceptempty=.true.)
      end if
      call iw_tooltip("Pick this point in the view: click an atom to use its position,&
         & or empty space to use the clicked point",ttshown)

    end subroutine region_pick_button

    !> Tooltip text for region coordinate row i (0 = origin/first
    !> corner/center) of region mode iregion.
    function region_point_tooltip(iregion,i) result(str)
      integer, intent(in) :: iregion
      integer, intent(in) :: i
      character(len=:), allocatable :: str

      if (iregion == iso_region_frac) then
         if (i == 0) then
            str = "First corner of the region (fractional coordinates)"
         else
            str = "Second corner of the region (fractional coordinates)"
         end if
      elseif (iregion == iso_region_ortho) then
         if (i == 0) then
            str = "First corner of the region (Cartesian coordinates, Å)"
         else
            str = "Second corner of the region (Cartesian coordinates, Å)"
         end if
      elseif (iregion == iso_region_parallel) then
         if (i == 0) then
            str = "Origin of the region parallelepiped (Cartesian coordinates, Å)"
         else
            str = "Endpoint of edge " // string(i) // " of the region parallelepiped (Cartesian, Å)"
         end if
      elseif (iregion == iso_region_cube) then
         if (i == 0) then
            str = "Center of the region (Cartesian coordinates, Å); seeded with the molecular centroid"
         else
            str = "Half the edge length of the cube"
         end if
      else ! iso_region_simplebox
         if (i == 0) then
            str = "Center of the region (Cartesian coordinates, Å); seeded with the molecular centroid"
         else
            str = "Half the box length along each Cartesian axis"
         end if
      end if

    end function region_point_tooltip
  end function draw_editrep_isosurface


  !> Draw the editrep window, critical points class. Returns true if
  !> the representation has changed. ttshown = the tooltip flag.
  module function draw_editrep_cps(w,ttshown) result(changed)
    use systems, only: sys
    use representations, only: cps_rad_def, cps_name
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_realc,&
       iw_field_combo, iw_calcwidth
    use param, only: bohrtoa
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: isys, ifield, i, it, ncount(0:3)
    logical :: ch

    ! initialize
    changed = .false.
    isys = w%isys

    ! field selector
    call iw_text("Field",highlight=.true.)
    call igSameLine(0._c_float,-1._c_float)
    ifield = w%rep%cps%ifield
    if (iw_field_combo("##cpsfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
       nonestr="<field not available>")) then
       w%rep%cps%ifield = ifield
       changed = .true.
    end if
    call iw_tooltip("Field whose critical points are displayed",ttshown)
    if (.not.sys(isys)%goodfield(w%rep%cps%ifield)) then
       call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       return
    end if

    ! count the critical points in the cell, by type (not the nuclei)
    ncount = 0
    associate(f => sys(isys)%f(w%rep%cps%ifield))
      do i = sys(isys)%c%nneq+1, f%ncp
         it = f%cp(i)%typind
         if (it >= 0 .and. it <= 3) ncount(it) = ncount(it) + f%cp(i)%mult
      end do
    end associate
    if (all(ncount == 0)) &
       call iw_text("This field has no critical points other than the nuclei (run AUTO)",&
          disabled=.true.,wrap=.true.)

    ! global radius scale
    ch = iw_dragfloat_realc("Radius scale",x1=w%rep%cps%radscale,speed=0.01_c_float,&
       min=0.01_c_float,max=10._c_float,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Scale factor applied to the radii of all critical points",ttshown)
    changed = changed .or. ch

    ! per-type color, radius, and show
    call iw_text("Critical points",highlight=.true.)
    do it = 0, 3
       ch = iw_coloredit("##cpscolor" // string(it),rgb=w%rep%cps%rgb(:,it))
       call iw_tooltip("Color of the critical points of this type",ttshown)
       changed = changed .or. ch
       ch = iw_dragfloat_realc("##cpsrad" // string(it),x1=w%rep%cps%rad(it),speed=0.002_c_float,&
          min=0.01_c_float,max=2._c_float,scale=real(bohrtoa,c_float),decimal=3,&
          sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Radius of the critical points of this type in Å, before the scale (default " //&
          string(real(cps_rad_def,8)*bohrtoa,'f',decimal=3) // " Å)",ttshown)
       changed = changed .or. ch
       ch = iw_checkbox(trim(cps_name(it)) // "##cpsshow" // string(it),w%rep%cps%show(it),sameline=.true.)
       call iw_tooltip("Show the critical points of this type (count: number in the cell)",ttshown)
       changed = changed .or. ch
       call iw_text("(" // string(ncount(it)) // " in the cell)",disabled=.true.,sameline=.true.)
    end do

  end function draw_editrep_cps

  !> Draw the editrep window, gradient paths class. Returns true if
  !> the representation has changed.
  module function draw_editrep_gpaths(w,ttshown) result(changed)
    use systems, only: sys, sysc, atlisttype_species, atlisttype_nneq
    use gui_main, only: ColorHighlightScene
    use representations, only: gpaths_rad_def, field_has_cps
    use utils, only: iw_text, iw_tooltip, iw_coloredit, iw_dragfloat_realc,&
       iw_field_combo, iw_calcwidth, iw_calcheight, iw_combo_simple, iw_checkbox,&
       iw_table_column, iw_table_headers_row, iw_atom_button, iw_highlight_selectable, iw_button, iw_inputfloat
    use param, only: bohrtoa
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: isys, iview, ifield, ihighlight, highlight_type, i, j, k, ncol, ihover(3)
    integer :: nuniq, ncell
    integer, allocatable :: iuniq(:,:), icell(:,:)
    integer(c_int) :: istyle, itable
    logical, pointer :: pshown1(:)
    logical :: ldum
    logical :: ch, typechanged, pathsok
    real(c_float) :: xcol
    character(kind=c_char,len=:), allocatable, target :: str1, suffix
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    type(ImVec2) :: sz

    ! initialize
    changed = .false.
    isys = w%isys
    iview = w%anchor_view()
    ihighlight = 0
    highlight_type = atlisttype_species
    ihover = 0

    ! the widgets start at one column, after the longest label
    xcol = igGetCursorPosX() + iw_calcwidth(len("Style "),0)

    ! field selector
    call label("Field")
    ifield = w%rep%gpaths%ifield
    if (iw_field_combo("##gpathsfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
       nonestr="<field not available>")) then
       w%rep%gpaths%ifield = ifield
       changed = .true.
    end if
    call iw_tooltip("Field whose gradient paths are displayed",ttshown)
    if (.not.sys(isys)%goodfield(w%rep%gpaths%ifield)) then
       call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       return
    end if
    if (.not.field_has_cps(isys,w%rep%gpaths%ifield,withpaths=.true.)) &
       call iw_text("This field has no gradient paths (run AUTO, or load a checkpoint that has them)",&
          disabled=.true.,wrap=.true.)

    ! style of all paths
    call label("Style")
    istyle = int(w%rep%gpaths%style,c_int)
    call iw_combo_simple("##gpathsstyle","Continuous path" // c_null_char //&
       "String of spheres" // c_null_char,istyle,changed=ch,tooltips=&
       "Draw each path as a continuous tube" // c_null_char //&
       "Draw each path as a string of spheres, one per path point" // c_null_char,ttshown=ttshown)
    ! the combo items are in the order of the gpaths_style_* values
    if (ch) then
       w%rep%gpaths%style = int(istyle)
       changed = .true.
    end if

    ! global color and radius: set those of all paths
    pathsok = w%rep%gpaths%paths_ok(isys)
    call label("Color")
    ch = iw_coloredit("##gpathscolor",rgb=w%rep%gpaths%rgb)
    call iw_tooltip("Color of all gradient paths (sets the color of every path in the table below)",ttshown)
    if (ch) then
       call w%rep%gpaths%fill_rgb()
       changed = .true.
    end if
    ch = iw_dragfloat_realc("Radius (Å)##gpathsrad",x1=w%rep%gpaths%rad,speed=0.002_c_float,&
       min=0.005_c_float,max=1._c_float,scale=real(bohrtoa,c_float),decimal=3,&
       sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Radius of the tubes or spheres of all gradient paths in Å (sets the radius of every "//&
       "path in the table below; default " // string(real(gpaths_rad_def,8)*bohrtoa,'f',decimal=3) // " Å)",&
       ttshown)
    if (ch) then
       if (pathsok) w%rep%gpaths%prad = w%rep%gpaths%rad
       changed = .true.
    end if

    ! the tables of gradient paths (for now, the two bond paths of each
    ! BCP): one row per path of each symmetry-unique CP, and one per
    ! path of each CP in the cell, in a second tab when there are more
    nuniq = 0
    ncell = 0
    if (pathsok) then
       associate(c => sys(isys)%c, f => sys(isys)%f(w%rep%gpaths%ifield))
         ! the rows
         allocate(iuniq(2,2*f%ncp),icell(2,2*f%ncpcel))
         do i = c%nneq+1, f%ncp
            do j = 1, 2
               if (f%cpgp(j,i)%n < 2) cycle
               nuniq = nuniq + 1
               iuniq(:,nuniq) = (/j,i/)
            end do
         end do
         do k = c%ncel+1, f%ncpcel
            do j = 1, 2
               if (f%cpgp(j,f%cpcel(k)%idx)%n < 2) cycle
               ncell = ncell + 1
               icell(:,ncell) = (/j,k/)
            end do
         end do
       end associate
    end if
    if (nuniq > 0) then
       call iw_text("Gradient Paths",highlight=.true.)
       if (ncell > nuniq) then
          ! which paths: symmetry-unique or cell (the combo items are
          ! in the order of the tablecell values)
          itable = int(w%rep%gpaths%tablecell,c_int)
          call iw_combo_simple("Path list##gpathstablecombo","Symmetry-unique" // c_null_char //&
             "Cell" // c_null_char,itable,changed=ch,tooltips=&
             "The paths of the symmetry-unique critical points (an edit applies to all their&
             & copies in the cell)" // c_null_char //&
             "The paths of every critical point in the cell" // c_null_char,ttshown=ttshown)
          if (ch) w%rep%gpaths%tablecell = int(itable)
          call draw_path_table(w%rep%gpaths%tablecell == 1,merge(ncell,nuniq,w%rep%gpaths%tablecell == 1))
       else
          call draw_path_table(.false.,nuniq)
       end if

       ! show/hide all paths
       pshown1(1:size(w%rep%gpaths%pshown)) => w%rep%gpaths%pshown
       if (showhide_buttons(pshown1,"gpaths","gradient paths",ttshown)) changed = .true.

       ! show only the bond paths within or between molecules
       if (iw_button("Intramolecular##gpathsintramol")) call select_paths(1)
       call iw_tooltip("Show the bond paths between atoms of the same molecule and hide all others",ttshown)
       if (iw_button("Intermolecular##gpathsintermol",sameline=.true.)) call select_paths(2)
       call iw_tooltip("Show the bond paths between atoms of different molecules and hide all others",ttshown)

       ! show only the bond paths with a field value at the BCP above
       ! or below a threshold
       if (iw_button("Above##gpathsabove")) call select_paths(3)
       call iw_tooltip("Show the bond paths whose BCP has a field value above the threshold "//&
          "and hide all others",ttshown)
       if (iw_button("Below##gpathsbelow",sameline=.true.)) call select_paths(4)
       call iw_tooltip("Show the bond paths whose BCP has a field value below the threshold "//&
          "and hide all others",ttshown)
       ldum = iw_inputfloat("Field value at the BCP##gpathsfthr",w%rep%gpaths%fthr,decimal=6,width=12,&
          sameline=.true.)
       call iw_tooltip("Threshold for the Above and Below buttons (field value at the bond critical point, "//&
          "atomic units)",ttshown)
    end if

    ! the atoms at the ends of the bond paths, with the colors and radii
    ! of this object (the table below, shown only when they are drawn)
    ch = iw_checkbox("Show atoms at bond path ends##gpathsshowends",w%rep%gpaths%showends)
    changed = changed .or. ch
    call iw_tooltip("Also draw the atoms at the ends of the bond paths, with the colors and radii "//&
       "in the table below.",ttshown)
    if (w%rep%gpaths%showends) then
       ch = atom_table_widget(isys,w%rep%atoms%style%type,typechanged,ihighlight,highlight_type,&
          rgb=w%rep%atoms%style%rgb,rad=w%rep%atoms%style%rad)
       if (typechanged) call w%rep%atoms%style%reset(w%rep)
       changed = changed .or. ch .or. typechanged
    end if

    ! the path under the mouse in the table is highlighted in the view
    ! (a rebuild, not an edit of the object)
    if (any(ihover /= w%rep%gpaths%ihover)) then
       w%rep%gpaths%ihover = ihover
       if (iview > 0) win(iview)%sc%forcebuildlists = .true.
    end if

    ! process transient highlights
    if (ihighlight > 0) then
       call sysc(isys)%highlight_atoms(.true.,(/ihighlight/),highlight_type,&
          reshape(ColorHighlightScene,(/4,1/)))
    end if

  contains
    !> A label in the label column; the widget follows at xcol.
    subroutine label(str)
      character(len=*), intent(in) :: str

      call iw_text(str,alignframe=.true.)
      call igSameLine(xcol,-1._c_float)

    end subroutine label

    !> The table of the paths of the symmetry-unique CPs (cell = .false.)
    !> or of the cell CPs (cell = .true.), with n rows.
    subroutine draw_path_table(cell,n)
      logical, intent(in) :: cell
      integer, intent(in) :: n

      integer :: k
      integer(c_int) :: flags

      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      flags = ior(flags,ImGuiTableFlags_ScrollY)
      if (cell) then
         str1 = "##tablegpathscell" // c_null_char
      else
         str1 = "##tablegpathsuniq" // c_null_char
      end if
      sz%x = 0._c_float
      sz%y = iw_calcheight(min(10,n+1),0,.false.)
      if (igBeginTable(c_loc(str1),5,flags,sz,0._c_float)) then
         ncol = -1
         call iw_table_column("Show",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Color",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Radius (Å)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("End 1",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("End 2",icol=ncol,flags=ImGuiTableColumnFlags_WidthStretch)
         call iw_table_headers_row(freezetop=.true.,autofit=.true.)

         clipper = ImGuiListClipper_ImGuiListClipper()
         call ImGuiListClipper_Begin(clipper,n,-1._c_float)
         do while(ImGuiListClipper_Step(clipper))
            call c_f_pointer(clipper,clipper_f)
            do k = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
               call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
               if (cell) then
                  call draw_path_row(k,icell(1,k),0,icell(2,k))
               else
                  call draw_path_row(k,iuniq(1,k),iuniq(2,k),0)
               end if
            end do
         end do
         call ImGuiListClipper_End(clipper)
         call ImGuiListClipper_destroy(clipper)
         call igEndTable()
      end if

    end subroutine draw_path_table

    !> Row k of a table: path j of the symmetry-unique CP i (icp = 0;
    !> the values of its first cell copy, and an edit sets all copies),
    !> or of the cell CP icp (i = 0).
    subroutine draw_path_row(k,j,i,icp)
      integer, intent(in) :: k, j, i, icp

      integer :: iend, id, iu
      integer, allocatable :: kc(:)
      integer(c_int) :: idx(4)
      real(c_float) :: rgb(3), rad, xpos
      logical :: have, ldum, lsh
      character(len=:), allocatable :: lbl

      associate(c => sys(isys)%c, f => sys(isys)%f(w%rep%gpaths%ifield), g => w%rep%gpaths)
        ! the cell paths of the row (the first one gives the values
        ! shown), and its unique CP
        if (icp > 0) then
           kc = (/icp/)
           iu = f%cpcel(icp)%idx
           suffix = "_c" // string(k)
        else
           kc = g%pcopy(g%pfirst(i):g%pfirst(i+1)-1)
           iu = i
           suffix = "_u" // string(k)
        end if
        lsh = g%pshown(j,kc(1))
        rgb = g%prgb(:,j,kc(1))
        rad = g%prad(j,kc(1))

        ! show, and the row selectable that highlights the path in the
        ! view: first in the row, so the highlight holds while the mouse
        ! is on the checkbox too
        if (igTableSetColumnIndex(0_c_int)) then
           xpos = igGetCursorPosX()
           call iw_text("",alignframe=.true.)
           if (iw_highlight_selectable("##gpathsrowsel" // suffix)) ihover = (/j,i,icp/)
           call igSameLine(0._c_float,0._c_float)
           call igSetCursorPosX(xpos)
           if (iw_checkbox("##gpathsrowshown" // suffix,lsh)) then
              g%pshown(j,kc) = lsh
              changed = .true.
           end if
           call iw_tooltip("Show this gradient path",ttshown)
           call mixed_note(any(g%pshown(j,kc) .neqv. g%pshown(j,kc(1))))
        end if

        ! color and radius of this path
        if (igTableSetColumnIndex(1_c_int)) then
           if (iw_coloredit("##gpathsrowcolor" // suffix,rgb=rgb)) then
              g%prgb(1,j,kc) = rgb(1)
              g%prgb(2,j,kc) = rgb(2)
              g%prgb(3,j,kc) = rgb(3)
              changed = .true.
           end if
           call iw_tooltip("Color of this gradient path",ttshown)
           call mixed_note(any(g%prgb(:,j,kc) /= spread(g%prgb(:,j,kc(1)),2,size(kc))))
        end if
        if (igTableSetColumnIndex(2_c_int)) then
           if (iw_dragfloat_realc("##gpathsrowrad" // suffix,x1=rad,speed=0.002_c_float,&
              min=0.005_c_float,max=1._c_float,scale=real(bohrtoa,c_float),decimal=3,&
              flags=ImGuiSliderFlags_AlwaysClamp)) then
              g%prad(j,kc) = rad
              changed = .true.
           end if
           call iw_tooltip("Radius of this gradient path in Å",ttshown)
           call mixed_note(any(g%prad(j,kc) /= g%prad(j,kc(1))))
        end if

        ! the ends: the CP the path starts from (a BCP)...
        if (igTableSetColumnIndex(3_c_int)) then
           lbl = trim(f%cp(iu)%name)
           if (icp > 0) lbl = lbl // " " // string(icp)
           call cp_badge(iview,isys,w%rep%gpaths%ifield,iu,lbl // "##gpathsend1" // suffix)
        end if

        ! ... and the one it ends at (a nucleus, in the color this
        ! object draws it with; or a non-nuclear attractor; none if it
        ! leaves the molecule)
        if (igTableSetColumnIndex(4_c_int)) then
           if (icp > 0) then
              ! cell path: the cell CP at the end, and its lattice vector
              iend = f%cpcel(icp)%ipath(j)
              if (iend >= 1 .and. iend <= c%ncel) then
                 call g%path_end(isys,icp,j,idx(1),idx(2:4))
                 have = (idx(1) > 0)
                 rgb = 0._c_float
                 if (have) then
                    id = sysc(isys)%attype_celatom_to_id(w%rep%atoms%style%type,idx(1))
                    have = (id >= 1 .and. id <= w%rep%atoms%style%ntype)
                    if (have) rgb = w%rep%atoms%style%rgb(:,id)
                 end if
                 ldum = iw_atom_button(anchor_label(isys,idx,"?",species=.true.) // "##gpathsend2" // suffix,&
                    rgb,havergb=have,inert=.true.)
              elseif (iend > c%ncel .and. iend <= f%ncpcel) then
                 call cp_badge(iview,isys,w%rep%gpaths%ifield,f%cpcel(iend)%idx,&
                    trim(f%cp(f%cpcel(iend)%idx)%name) // " " // string(iend) //&
                    lvec_str(f%cpcel(icp)%ilvec(:,j)) // "##gpathsend2" // suffix)
              else
                 ! the cell list marks no end; the unique list tells
                 ! whether the path leaves the molecule
                 call no_end(f%cp(iu)%ipath(j))
              end if
           else
              ! symmetry-unique path: the symmetry-unique CP at the end
              iend = f%cp(i)%ipath(j)
              if (iend >= 1 .and. iend <= c%nneq) then
                 id = sysc(isys)%attype_type_id_to_id(atlisttype_nneq,iend,w%rep%atoms%style%type)
                 have = (id >= 1 .and. id <= w%rep%atoms%style%ntype)
                 rgb = 0._c_float
                 if (have) rgb = w%rep%atoms%style%rgb(:,id)
                 ldum = iw_atom_button(trim(c%at(iend)%name) // "##gpathsend2" // suffix,rgb,&
                    havergb=have,inert=.true.)
              elseif (iend > c%nneq .and. iend <= f%ncp) then
                 call cp_badge(iview,isys,w%rep%gpaths%ifield,iend,trim(f%cp(iend)%name) //&
                    "##gpathsend2" // suffix)
              else
                 call no_end(iend)
              end if
           end if
        end if
      end associate

    end subroutine draw_path_row

    !> Show only the bond paths of the BCPs that are intramolecular
    !> (mode 1: both ends are atoms of the same molecule),
    !> intermolecular (2: atoms of different molecules), or whose
    !> field value is above (3) or below (4) the threshold; hide all
    !> other paths.
    subroutine select_paths(mode)
      integer, intent(in) :: mode

      integer :: kk, a1, a2, l1(3), l2(3)
      logical :: ok

      associate(c => sys(isys)%c, f => sys(isys)%f(w%rep%gpaths%ifield), g => w%rep%gpaths)
        ! the molecules need the molecular data of the crystal
        if (mode <= 2 .and. .not.allocated(c%idatcelmol)) return
        g%pshown = .false.
        do kk = c%ncel+1, f%ncpcel
           if (mode <= 2) then
              ! both ends must be nuclei (atoms)
              call g%path_end(isys,kk,1,a1,l1)
              call g%path_end(isys,kk,2,a2,l2)
              if (a1 == 0 .or. a2 == 0) cycle
              ok = c%in_same_molecule(a1,l1,a2,l2)
              if (mode == 2) ok = .not.ok
           elseif (mode == 3) then
              ok = (f%cp(f%cpcel(kk)%idx)%s%f > g%fthr)
           else
              ok = (f%cp(f%cpcel(kk)%idx)%s%f < g%fthr)
           end if
           g%pshown(:,kk) = ok
        end do
      end associate
      changed = .true.

    end subroutine select_paths

    !> "mixed" after a widget of a row whose cell copies have
    !> different values for it (ismixed).
    subroutine mixed_note(ismixed)
      logical, intent(in) :: ismixed

      if (.not.ismixed) return
      call iw_text("mixed",disabled=.true.,sameline=.true.)
      call iw_tooltip("The cell copies of this path have different values (Cell tab); "//&
         "an edit here sets all of them",ttshown)

    end subroutine mixed_note

    !> The end cell of a path that does not end at a CP.
    subroutine no_end(iend)
      integer, intent(in) :: iend

      if (iend == -1) then
         call iw_text("(leaves the molecule)",disabled=.true.)
      else
         call iw_text("?",disabled=.true.)
      end if

    end subroutine no_end

  end function draw_editrep_gpaths

  !> Draw the overlay of the blocking job of the edit representation
  !> window: sampling the field of an isosurface (run_editrep).
  module subroutine block_editrep(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_wait_overlay
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    character(len=:), allocatable :: info
    integer :: isys

    info = ""
    isys = w%isys
    if (ok_system(isys,sys_init)) then
       info = "System: " // string(isys) // ": " // trim(sysc(isys)%seed%name)
       if (associated(w%rep)) then
          if (sys(isys)%goodfield(w%rep%iso%ifield)) &
             info = info // newline // "Field:  " // string(w%rep%iso%ifield) // ": " //&
             trim(sys(isys)%f(w%rep%iso%ifield)%name)
       end if
    end if
    info = info // newline // "Grid:   " // string(w%editrep_pending_n(1)) // " x " //&
       string(w%editrep_pending_n(2)) // " x " // string(w%editrep_pending_n(3))
    call iw_wait_overlay("Sampling the isosurface...",info,.true.)

  end subroutine block_editrep

  !> Run the blocking job of the edit representation window, called by
  !> the main loop after the frame with the overlay: sample the field of
  !> the isosurface on the staged grid and region, and apply them only
  !> if the sampling completes. A cancel (Esc) leaves the isosurface as
  !> it was.
  module subroutine run_editrep(w)
    use representations, only: reptype_isosurface
    use systems, only: sys, sys_init, ok_system
    use gui_main, only: begin_cancellable, end_cancellable
    class(window), intent(inout), target :: w

    integer :: isys, iview
    logical :: cancelled, outdomain, ok
    real*8, allocatable :: ff(:,:,:)

    ! the representation, its field, and its view must still be there
    isys = w%isys
    iview = w%anchor_view()
    if (iview == 0 .or. .not.associated(w%rep) .or. .not.ok_system(isys,sys_init)) return
    if (.not.win(iview)%isopen .or. .not.associated(win(iview)%sc)) return
    if (win(iview)%isys /= isys) return
    if (.not.w%rep%isinit .or. w%rep%type /= reptype_isosurface) return
    if (.not.sys(isys)%goodfield(w%rep%iso%ifield)) return

    associate(iso => w%rep%iso)
      call begin_cancellable()
      call iso%sample(isys,w%editrep_pending_n,iso%iregion,iso%rgn_x,iso%imosel,iso%imoidx,&
         ff,outdomain,ok)
      cancelled = end_cancellable()
      if (cancelled) then
         w%errmsg = "The sampling was cancelled: the grid is unchanged"
      else
         w%errmsg = ""
         call iso%apply_grid(w%editrep_pending_n,iso%iregion,iso%rgn_x)
         if (ok) call iso%set_samples(isys,ff,outdomain)
         win(iview)%sc%forcebuildlists = .true.
         call win(iview)%sc%undo_note("Calculate grid",.false.)
      end if
    end associate

  end subroutine run_editrep

  !xx! private procedures

  !> The prompt for tool itool (objtool_*, or objtool_kind0 + planarkind_*) of the 2D drawing editor
  !> in the view bar: short enough for it (vmbar_maxlen); the tooltip in
  !> the toolbar (planar_tool_hint) has the rest.
  function planar_tool_prompt(itool) result(str)
    use representations, only: planarkind_ellipse, planarkind_rect, planarkind_arrow,&
       planarkind_freehand, planarkind_curve, planarkind_NUM
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    integer :: ik

    ik = itool - objtool_kind0
    if (itool == objtool_select) then
       str = "Click a shape to select it; drag to edit it"
    elseif (itool == objtool_remove) then
       str = "Click a shape to remove it"
    elseif (ik == planarkind_ellipse) then
       str = "Drag to draw an ellipse"
    elseif (ik == planarkind_rect) then
       str = "Drag to draw a rectangle"
    elseif (ik == planarkind_arrow) then
       str = "Drag from the tail to the tip of the arrow"
    elseif (ik == planarkind_curve) then
       str = "Drag from the tail to the tip of the curve"
    elseif (ik == planarkind_freehand) then
       str = "Drag to draw a freehand line"
    elseif (ik >= 1 .and. ik <= planarkind_NUM) then
       str = "Click to add points; double-click to finish"
    else
       str = ""
    end if

  end function planar_tool_prompt

  !> The tooltip for tool itool (objtool_*, or objtool_kind0 + planarkind_*) of the 2D drawing
  !> editor in its toolbar. The keys are those bound in the 2D Drawing
  !> mode.
  function planar_tool_hint(itool) result(str)
    use representations, only: planarkind_ellipse, planarkind_rect, planarkind_arrow,&
       planarkind_freehand, planarkind_NUM, planarkind_name, planarkind_curve
    use keybindings, only: BIND_OBJEDIT_DRAW, BIND_OBJEDIT_EXIT, BIND_OBJEDIT_CONSTRAIN,&
       BIND_OBJEDIT_DELETE, BIND_OBJEDIT_FINISH, BIND_OBJEDIT_DELPOINT
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    integer :: ik

    ik = itool - objtool_kind0
    if (itool == objtool_select) then
       str = "Select: click a shape (" // kn(BIND_OBJEDIT_DRAW) // ") to select it, drag it to move it, " //&
          "and drag its handles to resize, rotate, or move its points (" //&
          kn(BIND_OBJEDIT_CONSTRAIN) // ": constrain). " // kn(BIND_OBJEDIT_DELETE) //&
          " removes the selected shape"
    elseif (itool == objtool_remove) then
       str = "Remove: click a shape (" // kn(BIND_OBJEDIT_DRAW) // ") to remove it"
    elseif (ik == planarkind_ellipse) then
       str = "Ellipse: drag to draw (" // kn(BIND_OBJEDIT_CONSTRAIN) // ": a circle)"
    elseif (ik == planarkind_rect) then
       str = "Rectangle: drag to draw (" // kn(BIND_OBJEDIT_CONSTRAIN) // ": a square)"
    elseif (ik == planarkind_arrow) then
       str = "Arrow or line: drag from the tail to the tip (" // kn(BIND_OBJEDIT_CONSTRAIN) //&
          ": snap to 45°)"
    elseif (ik == planarkind_curve) then
       str = "Curved arrow or line: drag from the tail to the tip (" // kn(BIND_OBJEDIT_CONSTRAIN) //&
          ": snap to 45°); with Select, drag the middle handle to bend it"
    elseif (ik == planarkind_freehand) then
       str = "Freehand: drag to draw"
    elseif (ik >= 1 .and. ik <= planarkind_NUM) then
       str = trim(planarkind_name(ik)) // ": click (" // kn(BIND_OBJEDIT_DRAW) // ") to add points; " //&
          "double-click, " // kn(BIND_OBJEDIT_FINISH) // ", or " // kn(BIND_OBJEDIT_EXIT) //&
          " to finish (" // kn(BIND_OBJEDIT_DELPOINT) // ": remove the last point)"
    else
       str = ""
    end if
    if (len(str) > 0) str = str // ". " // exit_hint()

  end function planar_tool_hint

  !> The prompt for tool itool (objtool_*, or objtool_kind0 + shapekind_*)
  !> of the 3D shapes editor in the view bar: short enough for it
  !> (vmbar_maxlen); the tooltip in the toolbar (shapes_tool_hint) has the
  !> rest.
  function shapes_tool_prompt(itool) result(str)
    use representations, only: shapekind_sphere, shapekind_box, shapekind_arrow,&
       shapekind_cone, shapekind_cylinder
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    integer :: ik

    ik = itool - objtool_kind0
    if (itool == objtool_select) then
       str = "Click a shape to select it; drag to edit it"
    elseif (itool == objtool_remove) then
       str = "Click a shape to remove it"
    elseif (ik == shapekind_sphere) then
       str = "Drag from the center out to the radius"
    elseif (ik == shapekind_box) then
       str = "Drag between two opposite corners"
    elseif (ik == shapekind_arrow) then
       str = "Drag from the tail to the tip of the arrow"
    elseif (ik == shapekind_cone) then
       str = "Drag from the base to the apex of the cone"
    elseif (ik == shapekind_cylinder) then
       str = "Drag from one end to the other"
    else
       str = ""
    end if

  end function shapes_tool_prompt

  !> The tooltip for tool itool (objtool_*, or objtool_kind0 + shapekind_*)
  !> of the 3D shapes editor in its toolbar. The keys are those bound in
  !> the object editing mode.
  function shapes_tool_hint(itool) result(str)
    use representations, only: shapekind_sphere, shapekind_box, shapekind_arrow,&
       shapekind_cone, shapekind_cylinder
    use keybindings, only: BIND_OBJEDIT_DRAW, BIND_OBJEDIT_CONSTRAIN, BIND_OBJEDIT_DELETE,&
       BIND_OBJEDIT_NOSNAP
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    integer :: ik
    character(len=:), allocatable :: snap

    ik = itool - objtool_kind0
    snap = ". The points snap to the atom under the mouse (" // kn(BIND_OBJEDIT_NOSNAP) //&
       ": do not snap)"
    if (itool == objtool_select) then
       str = "Select: click a shape (" // kn(BIND_OBJEDIT_DRAW) // ") to select it, drag it to " //&
          "move it, and drag its handles to resize it or move its ends (" //&
          kn(BIND_OBJEDIT_CONSTRAIN) // ": along an axis). " // kn(BIND_OBJEDIT_DELETE) //&
          " removes the selected shape" // snap
    elseif (itool == objtool_remove) then
       str = "Remove: click a shape (" // kn(BIND_OBJEDIT_DRAW) // ") to remove it"
    elseif (ik == shapekind_sphere) then
       str = "Sphere: drag from the center out to the radius" // snap
    elseif (ik == shapekind_box) then
       str = "Box: drag between two opposite corners, with edges along the cell (crystal) " //&
          "or cartesian (molecule) axes; the depth is the shorter edge (" //&
          kn(BIND_OBJEDIT_CONSTRAIN) // ": a cube)" // snap
    elseif (ik == shapekind_arrow) then
       str = "Arrow: drag from the tail to the tip (" // kn(BIND_OBJEDIT_CONSTRAIN) //&
          ": along an axis)" // snap
    elseif (ik == shapekind_cone) then
       str = "Cone: drag from the center of the base to the apex (" //&
          kn(BIND_OBJEDIT_CONSTRAIN) // ": along an axis)" // snap
    elseif (ik == shapekind_cylinder) then
       str = "Cylinder: drag from one end to the other (" // kn(BIND_OBJEDIT_CONSTRAIN) //&
          ": along an axis)" // snap
    else
       str = ""
    end if
    if (len(str) > 0) str = str // ". " // exit_hint()

  end function shapes_tool_hint

  !> The prompt for tool itool (objtool_*, or objtool_kind0 + 1 +
  !> textpos_*) of the text editor in the view bar: short enough for it
  !> (vmbar_maxlen); the tooltip in the toolbar (text_tool_hint) has the
  !> rest.
  function text_tool_prompt(itool) result(str)
    use representations, only: textpos_screen, textpos_point, textpos_atom, textpos_bond
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    integer :: ipl

    ipl = itool - objtool_kind0 - 1
    if (itool == objtool_select) then
       str = "Click a text to select it; drag to move it"
    elseif (itool == objtool_remove) then
       str = "Click a text to remove it"
    elseif (ipl == textpos_screen) then
       str = "Click to place an on-screen text"
    elseif (ipl == textpos_point) then
       str = "Click to place a text at a 3D point"
    elseif (ipl == textpos_atom) then
       str = "Click an atom to place a text on it"
    elseif (ipl == textpos_bond) then
       str = "Click a bond to place a text on it"
    else
       str = ""
    end if

  end function text_tool_prompt

  !> The tooltip for tool itool (objtool_*, or objtool_kind0 + 1 +
  !> textpos_*) of the text editor in its toolbar. The keys are those
  !> bound in the object editing mode.
  function text_tool_hint(itool) result(str)
    use representations, only: textpos_screen, textpos_point, textpos_atom, textpos_bond
    use keybindings, only: BIND_OBJEDIT_DRAW, BIND_OBJEDIT_CONSTRAIN, BIND_OBJEDIT_DELETE,&
       BIND_OBJEDIT_NOSNAP
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    integer :: ipl
    character(len=:), allocatable :: typing

    ipl = itool - objtool_kind0 - 1
    typing = ". Type the text right away: the text box of its object editor takes the keyboard"
    if (itool == objtool_select) then
       str = "Select: click a text (" // kn(BIND_OBJEDIT_DRAW) // ") to select it and drag it to " //&
          "move it: an on-screen text in the window, a text at a 3D point in the plane facing " //&
          "the camera (onto the atom under the mouse; " // kn(BIND_OBJEDIT_NOSNAP) // ": do not " //&
          "snap; " // kn(BIND_OBJEDIT_CONSTRAIN) // ": along an axis), and a text on an atom or " //&
          "bond away from it. Double-click a text to edit it. " // kn(BIND_OBJEDIT_DELETE) //&
          " removes the selected text"
    elseif (itool == objtool_remove) then
       str = "Remove: click a text (" // kn(BIND_OBJEDIT_DRAW) // ") to remove it"
    elseif (ipl == textpos_screen) then
       str = "On-screen text: click to place a text that stays at that position of the window" //&
          typing
    elseif (ipl == textpos_point) then
       str = "Text at a 3D point: click an atom to place the text at its position, or empty " //&
          "space for the point on the plane through the scene center (" //&
          kn(BIND_OBJEDIT_NOSNAP) // ": do not snap to atoms)" // typing
    elseif (ipl == textpos_atom) then
       str = "Text on an atom: click an atom to tie a text to it" // typing
    elseif (ipl == textpos_bond) then
       str = "Text on a bond: click a bond to tie a text to its midpoint" // typing
    else
       str = ""
    end if
    if (len(str) > 0) str = str // ". " // exit_hint()

  end function text_tool_hint

  !> The prompt for tool itool (objtool_*, or objtool_kind0 + number of
  !> atoms - 1) of the measurement editor in the view bar: short enough
  !> for it (vmbar_maxlen); the tooltip in the toolbar
  !> (measure_tool_hint) has the rest.
  function measure_tool_prompt(itool) result(str)
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    integer :: n

    n = itool - objtool_kind0 + 1
    if (itool == objtool_select) then
       str = "Click a measurement to select it; drag to edit it"
    elseif (itool == objtool_remove) then
       str = "Click a measurement to remove it"
    elseif (n == 2) then
       str = "Click two atoms, or a bond"
    elseif (n == 3) then
       str = "Click three atoms; the second is the vertex"
    elseif (n == 4) then
       str = "Click four atoms; the 2nd and 3rd are the axis"
    else
       str = ""
    end if

  end function measure_tool_prompt

  !> The tooltip for tool itool (objtool_*, or objtool_kind0 + number of
  !> atoms - 1) of the measurement editor in its toolbar. The keys are
  !> those bound in the object editing mode.
  function measure_tool_hint(itool) result(str)
    use keybindings, only: BIND_OBJEDIT_DRAW, BIND_OBJEDIT_EXIT, BIND_OBJEDIT_DELETE
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    integer :: n
    character(len=:), allocatable :: picks

    n = itool - objtool_kind0 + 1
    picks = ". The picked atoms are numbered in the view, and the values show at the " //&
       "mouse; click a picked atom again to drop it, " // kn(BIND_OBJEDIT_EXIT) //&
       " drops them all. Critical points can be picked too. Measuring the same atoms " //&
       "again selects the measurement"
    if (itool == objtool_select) then
       str = "Select: click a measurement (" // kn(BIND_OBJEDIT_DRAW) // "), its label or its " //&
          "lines, to select it; drag its label to move it, or one of its atoms onto another " //&
          "atom to measure from there. " // kn(BIND_OBJEDIT_DELETE) // " removes the selected " //&
          "measurement"
    elseif (itool == objtool_remove) then
       str = "Remove: click a measurement (" // kn(BIND_OBJEDIT_DRAW) // ") to remove it"
    elseif (n == 2) then
       str = "Distance: click two atoms (" // kn(BIND_OBJEDIT_DRAW) // "), or one bond" // picks
    elseif (n == 3) then
       str = "Angle: click three atoms (" // kn(BIND_OBJEDIT_DRAW) // "); the second one is " //&
          "the vertex" // picks
    elseif (n == 4) then
       str = "Dihedral: click four atoms (" // kn(BIND_OBJEDIT_DRAW) // "); the second and " //&
          "third ones are the axis" // picks
    else
       str = ""
    end if
    if (len(str) > 0) str = str // ". " // exit_hint()

  end function measure_tool_hint

  !> The prompt for tool itool (objtool_select, or objtool_kind0 + 1 +
  !> axplace_*) of the axes editor in the view bar: short enough for it
  !> (vmbar_maxlen); the tooltip in the toolbar (axes_tool_hint) has the
  !> rest.
  function axes_tool_prompt(itool) result(str)
    use representations, only: axplace_scene, axplace_window
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    if (itool == objtool_select) then
       str = "Drag the axes, an arrow tip, or a label"
    elseif (itool == objtool_kind0 + 1 + axplace_scene) then
       str = "Click an atom or a point to place the axes"
    elseif (itool == objtool_kind0 + 1 + axplace_window) then
       str = "Click to place the axes in the window"
    else
       str = ""
    end if

  end function axes_tool_prompt

  !> The tooltip for tool itool (objtool_select, or objtool_kind0 + 1 +
  !> axplace_*) of the axes editor in its toolbar. The keys are those
  !> bound in the object editing mode.
  function axes_tool_hint(itool) result(str)
    use representations, only: axplace_scene, axplace_window
    use keybindings, only: BIND_OBJEDIT_DRAW, BIND_OBJEDIT_CONSTRAIN, BIND_OBJEDIT_NOSNAP
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    if (itool == objtool_select) then
       str = "Select: drag (" // kn(BIND_OBJEDIT_DRAW) // ") the axes by their origin or " //&
          "shafts to move them (in the scene: onto the atom under the mouse; " //&
          kn(BIND_OBJEDIT_NOSNAP) // ": do not snap; " // kn(BIND_OBJEDIT_CONSTRAIN) //&
          ": along an axis), an arrow tip to resize them, or a label to move it"
    elseif (itool == objtool_kind0 + 1 + axplace_scene) then
       str = "Place in the scene: click an atom to put the origin of the axes on it, or empty " //&
          "space for the point on the plane through the scene center (" //&
          kn(BIND_OBJEDIT_NOSNAP) // ": do not snap to atoms)"
    elseif (itool == objtool_kind0 + 1 + axplace_window) then
       str = "Place in the window: click to anchor the axes at that position of the window, " //&
          "where they stay as the view moves"
    else
       str = ""
    end if
    if (len(str) > 0) str = str // ". " // exit_hint()

  end function axes_tool_hint

  !> The prompt for atom tool itool (objtool_kind0 + atomtool_*) in the
  !> view bar: short enough for it (vmbar_maxlen); the tooltip in the
  !> toolbar (atom_tool_hint) has the rest.
  function atom_tool_prompt(itool) result(str)
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    select case (itool - objtool_kind0)
    case (atomtool_paint)
       str = "Click or drag over atoms to paint them"
    case (atomtool_enlarge)
       str = "Click an atom to make it larger"
    case (atomtool_shrink)
       str = "Click an atom to make it smaller"
    case (atomtool_hide)
       str = "Click an atom to hide it"
    case (atomtool_poly)
       str = "Click an atom to toggle its polyhedron"
    case (atomtool_label)
       str = "Click an atom to show or hide its label"
    case (atomtool_shift)
       str = "Click an atom to shift it one cell"
    case default
       str = ""
    end select

  end function atom_tool_prompt

  !> The tooltip for atom tool itool (objtool_kind0 + atomtool_*) in the
  !> toolbar. The keys are those bound in the object editing mode.
  function atom_tool_hint(itool) result(str)
    use keybindings, only: BIND_OBJEDIT_DRAW
    integer, intent(in) :: itool
    character(len=:), allocatable :: str

    character(len=:), allocatable :: sel

    sel = ". Clicking a selected atom acts on the whole selection; if the atoms are not all of a " //&
       "kind, the object switches to one row per atom in its editor"
    select case (itool - objtool_kind0)
    case (atomtool_paint)
       str = "Paint: click (" // kn(BIND_OBJEDIT_DRAW) // ") atoms, or drag over them, to give " //&
          "them the paint color (chosen in the arrow menu of this button)" // sel
    case (atomtool_enlarge)
       str = "Enlarge: click (" // kn(BIND_OBJEDIT_DRAW) // ") an atom to make it 25% larger" // sel
    case (atomtool_shrink)
       str = "Shrink: click (" // kn(BIND_OBJEDIT_DRAW) // ") an atom to make it 25% smaller" // sel
    case (atomtool_hide)
       str = "Hide: click (" // kn(BIND_OBJEDIT_DRAW) // ") an atom to hide it in the view " //&
          "(Display Settings, or Reset in the arrow menu of this button, show it again)" // sel
    case (atomtool_poly)
       str = "Polyhedron: click (" // kn(BIND_OBJEDIT_DRAW) // ") an atom to show or hide the " //&
          "coordination polyhedron centered on it" // sel
    case (atomtool_label)
       str = "Label: click (" // kn(BIND_OBJEDIT_DRAW) // ") an atom to show or hide its label. " //&
          "The labels go by the label type of the labels object (atom names if it is made " //&
          "here): with atom names, all the symmetry-equivalent atoms of a crystal change " //&
          "together. Clicking a selected atom acts on the whole selection"
    case (atomtool_shift)
       str = "Lattice shift: click (" // kn(BIND_OBJEDIT_DRAW) // ") atoms to draw them one " //&
          "cell along the last lattice direction armed with the Lattice Shift buttons (+a, " //&
          "-a, +b, -b, +c, -c) in the arrow menu of this button (crystals only); the shifts " //&
          "add up, so the opposite direction undoes one. The Shift column of the Display " //&
          "Settings window shows and edits them" // sel
    case default
       str = ""
    end select
    if (len(str) > 0) str = str // ". " // exit_hint()

  end function atom_tool_hint

  !> The end of the toolbar tooltips of the object editors: the keys
  !> that turn the tool off.
  function exit_hint() result(str)
    use keybindings, only: BIND_OBJEDIT_EXIT, BIND_CANCEL
    character(len=:), allocatable :: str

    str = kn(BIND_OBJEDIT_EXIT) // " or " // kn(BIND_CANCEL) // ": back to navigation"

  end function exit_hint

  !> The name of the key bound to bind, trimmed for use in a message.
  function kn(bind)
    use keybindings, only: get_bind_keyname
    integer, intent(in) :: bind
    character(len=:), allocatable :: kn

    kn = trim(get_bind_keyname(bind))

  end function kn

end submodule editrep
