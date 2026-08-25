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

  !xx! private procedures
  ! function atom_selection_widget() result(changed)
  ! function periodicity_widgets(w,ttshown,tooltip_automatic,highlight) result(changed)

contains

  !> Update tasks for the edit representation window, before the
  !> window is created.
  module subroutine update_editrep(w)
    use systems, only: sys_init, ok_system
    class(window), intent(inout), target :: w

    integer :: isys
    logical :: doquit

    ! This window edits one object living in its anchor view's scene, so unlike
    ! the other view tools it does not follow the view to another system: the
    ! object belongs to the scene it was created in. Check the anchor, the
    ! system and the object are all still there, and close otherwise.
    isys = w%isys
    doquit = (w%anchor_view() == 0)
    if (.not.doquit) doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = .not.associated(w%rep)
    if (.not.doquit) doquit = .not.w%rep%isinit
    if (.not.doquit) doquit = (w%rep%type <= 0)

    ! if they aren't, quit the window
    if (doquit) call w%end()

  end subroutine update_editrep

  !> Draw the edit represenatation window.
  module subroutine draw_editrep(w)
    use representations, only: representation, reptype_atoms, reptype_unitcell, reptype_axes,&
       reptype_symelem, reptype_text, reptype_measure, reptype_isosurface
    use windows, only: win
    use keybindings, only: is_bind_event, BIND_OK_FOCUSED_DIALOG
    use systems, only: sysc, sys_init, ok_system
    use utils, only: iw_text, iw_tooltip, iw_button, iw_calcheight, iw_checkbox,&
       iw_inputtext, iw_close_event, iw_setpos_bottomright
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: isys, iview
    logical :: doquit, ok, ldum
    logical :: changed

    logical, save :: ttshown = .false. ! tooltip flag

    ! check the anchor, the system and the object are all still there (see
    ! update_editrep); the last test catches the anchor view moving to another
    ! system, which leaves this object behind in the old scene
    isys = w%isys
    iview = w%anchor_view()
    doquit = (iview == 0)
    if (.not.doquit) doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = .not.associated(w%rep)
    if (.not.doquit) doquit = .not.w%rep%isinit
    if (.not.doquit) doquit = (w%rep%type <= 0)
    if (.not.doquit) doquit = .not.win(iview)%isopen
    if (.not.doquit) doquit = .not.associated(win(iview)%sc)
    if (.not.doquit) doquit = (win(iview)%isys /= isys)

    if (.not.doquit) then
       ! whether the rep has changed
       changed = .false.

       ! system
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! name block (the representation type is fixed at creation and cannot be
       ! changed here)
       call iw_text("Name",highlight=.true.,alignframe=.true.)
       ldum = iw_inputtext("##nametextinput",bufsize=1023,texta=w%rep%name,width=30,sameline=.true.)
       call iw_tooltip("Name of this object",ttshown)

       ! shown checkbox
       changed = changed .or. iw_checkbox("Show",w%rep%shown,sameline=.true.)
       call iw_tooltip("Toggle show/hide this object",ttshown)

       ! type-dependent items
       if (w%rep%type == reptype_atoms) then
          changed = changed .or. w%draw_editrep_atoms(ttshown)
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
       end if

       ! rebuild draw lists if necessary
       if (changed) win(iview)%sc%forcebuildlists = .true.

       ! right-align and bottom-align for the rest of the contents
       call iw_setpos_bottomright(10,2)

       ! reset button
       if (iw_button("Reset",danger=.true.)) then
          call w%rep%set_defaults(0)
          win(iview)%sc%forcebuildlists = .true.
       end if

       ! close button
       ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
       ok = ok .or. iw_button("Close",sameline=.true.)
       if (ok) doquit = .true.
    end if

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_editrep

  !> Draw the editrep (Object) window, atoms class. Returns true if
  !> the scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_atoms(w,ttshown) result(changed)
    use representations, only: representation
    use systems, only: sys, sysc, atlisttype_species, atlisttype_ncel_ang, atlisttype_nmol,&
       atlisttype_nneq, atlisttype_ncel_frac
    use gui_main, only: ColorHighlightScene, ColorElement
    use tools_io, only: string
    use utils, only: iw_text, iw_tooltip, iw_helpermark, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_radiobutton, iw_calcheight, iw_clamp_color3, iw_checkbox, iw_coloredit,&
       iw_highlight_selectable, iw_dragfloat_real8, iw_inputtext, iw_inputint, iw_table_column,&
       iw_begintabitem
    use param, only: atmcov, atmvdw, atmcov0, newline, jmlcol, jmlcol2, bohrtoa
    use global, only: bondfactor_def, bonddelta_def
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: ispc, isys, iz
    character(kind=c_char,len=:), allocatable, target :: str1, str3, suffix
    real*8 :: x0(3)
    logical :: ch, ldum
    integer(c_int) :: lst, flags, nspcpair
    integer :: i, j, k, intable, nrow, is, ncol, ihighlight, highlight_type
    integer :: itype_combo, newtype
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    integer, allocatable :: indi(:), indj(:)
    type(ImVec2) :: sz

    integer(c_int), parameter :: lsttrans(0:7) = (/0,1,2,2,2,3,4,5/)
    integer(c_int), parameter :: lsttransi(0:5) = (/0,1,2,5,6,7/)

    integer(c_int), parameter :: ic_sp1 = 0
    integer(c_int), parameter :: ic_sp2 = 1
    integer(c_int), parameter :: ic_shown = 2

    ! initialize
    ihighlight = 0
    highlight_type = atlisttype_species
    changed = .false.
    isys = w%isys

    ! update representation to respond to changes in number of atoms and molecules
    call w%rep%update()

    ! row of display options
    changed = changed .or. iw_checkbox("Atoms##atomsglobaldisplay",w%rep%atoms%display,highlight=.true.)
    call iw_tooltip("Display atoms in the scene",ttshown)

    changed = changed .or. iw_checkbox("Bonds##bondsglobaldisplay",w%rep%bonds%display,sameline=.true.,highlight=.true.)
    call iw_tooltip("Display bonds in the scene",ttshown)

    changed = changed .or. iw_checkbox("Labels##labelsglobaldisplay",w%rep%labels%display,sameline=.true.,highlight=.true.)
    call iw_tooltip("Display atomic labels in the scene",ttshown)

    changed = changed .or. iw_checkbox("Polyhedra##polyglobaldisplay",w%rep%poly%display,sameline=.true.,highlight=.true.)
    call iw_tooltip("Display coordination polyhedra around the center atoms in the scene",ttshown)

    str1 = "##editrepatomstabbar" // string(w%isys) // c_null_char
    flags = ImGuiTabBarFlags_Reorderable
    flags = ior(flags,ImGuiTabBarFlags_AutoSelectNewTabs)
    if (igBeginTabBar(c_loc(str1),flags)) then

       !!!!! Selection tab !!!!!
       if (iw_begintabitem("Selection##editrepatoms_selectiontab")) then
          ! filter
          call iw_text("Filter",highlight=.true.,alignframe=.true.)
          call iw_helpermark("Show the atom if the filter expression evaluates to non-zero (true) at the atomic position. &
             &Structural variables are very useful for filters. Examples:"//newline//&
             "- '@x < 3' = all atoms with x lower than 3"//newline//&
             "- 'log($0) > 1' = log of the promolecular density higher than 1"//newline//&
             "- 'abs(@x) < 2 && abs(@y) < 2 && abs(@z) < 2' = atoms in the (-2,2) box"//newline//&
             "Click the Help button for more info.")
          if (iw_button("Help##helpfilter",sameline=.true.)) then
             str3 = "https://aoterodelaroza.github.io/critic2/manual/arithmetics" // c_null_char
             call openLink(c_loc(str3))
          end if
          call iw_tooltip("Open the manual page regarding arithmetic expressions.",ttshown)

          ! filter text input
          if (iw_inputtext("##filtertext",bufsize=1023,texta=w%rep%sel%filter,notlive=.true.)) then
             ! test the filter
             if (sys(isys)%c%ncel > 0) then
                x0 = sys(isys)%c%atcel(1)%r
             else
                x0 = 0d0
             end if
             changed = .true.
             w%rep%sel%errfilter = ""
          end if
          if (len_trim(w%rep%sel%filter) == 0) w%rep%sel%errfilter = ""
          call iw_tooltip("Apply this filter to the atoms in the system. Atoms are represented if non-zero.",&
             ttshown)
          if (iw_button("Clear",sameline=.true.)) then
             w%rep%sel%filter = ""
             w%rep%sel%errfilter = ""
             changed = .true.
          end if
          call iw_tooltip("Clear the filter",ttshown)
          if (len_trim(w%rep%sel%errfilter) > 0) &
             call iw_text("Error: " // trim(w%rep%sel%errfilter),danger=.true.)

          ! periodicity (evaluate first: the widgets must be drawn even if
          ! changed is already true, and .or. is allowed to short-circuit)
          if (.not.sys(isys)%c%ismolecule) then
             ch = periodicity_widgets(w,ttshown,&
                "Number of periodic cells controlled by the +/- options in the 'Scene' button of the view window")
             changed = changed .or. ch

             ! checkbox for molecular motif
             changed = changed .or. iw_checkbox("Show connected molecules",w%rep%sel%onemotif)
             call iw_tooltip("Translate atoms to display whole molecules",ttshown)

             ! checkbox for border
             changed = changed .or. iw_checkbox("Show atoms at cell edges",w%rep%sel%border,sameline=.true.)
             call iw_tooltip("Display atoms near the unit cell edges",ttshown)
          end if

          ! origin of the atoms
          if (.not.sys(isys)%c%ismolecule) then
             ! origin translation
             changed = changed .or. iw_dragfloat_real8("Translate Origin (fractional)##originatom",x3=w%rep%sel%origin,&
                speed=0.001d0,decimal=5)
             call iw_tooltip("Translation vector for the contents of the unit cell.",ttshown)

             ! origin shift
             changed = changed .or. iw_dragfloat_real8("Cell Origin Shift (fractional)##origincell",x3=w%rep%sel%tshift,&
                speed=0.001d0,decimal=5)
             call iw_tooltip("Displace the origin of the cell being represented.",ttshown)
          end if

          ! show the corner atoms of coordination polyhedra
          changed = changed .or. iw_checkbox("Show atoms at polyhedra corners##polyshowcorners",&
             w%rep%poly%showcorners)
          call iw_tooltip("When coordination polyhedra are displayed and atoms are shown, also "//&
             "draw the atoms at the polyhedra corners, even if they fall outside the current selection.",ttshown)

          ! draw the atom selection widget
          changed = changed .or. atom_selection_widget(isys,w%rep,.true.,.false.,ihighlight,highlight_type)

          call igEndTabItem()
       end if ! begin tab item (selection)

       !!!!! Atoms tab !!!!!
       if (w%rep%atoms%display) then
          if (iw_begintabitem("Atoms##editrepatoms_atomstab")) then
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

             ! draw the atom selection widget
             changed = changed .or. atom_selection_widget(isys,w%rep,&
                .false.,.true.,ihighlight,highlight_type)

             call igEndTabItem()
          end if ! begin tab item (atoms)
       end if

       !!!!! Bonds tab !!!!!
       if (w%rep%bonds%display) then
          if (iw_begintabitem("Bonds##editrepatoms_bondstab")) then
             !! bonds display !!

             !! global options !!
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
                "Single color"//c_null_char//"Two colors"//c_null_char,w%rep%bonds%color_style,sameline=.true.,changed=ch)
             call iw_tooltip("Use a single color for the bond, or two colors from the bonded atoms",ttshown)

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

             !! atom selection block !!
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

                call igTableSetupScrollFreeze(0, 1) ! top row always visible

                ! draw the header
                call igTableHeadersRow()
                call igTableSetColumnWidthAutoAll(igGetCurrentTable())

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
                   call igTableSetupScrollFreeze(0,1)
                   call igTableHeadersRow()

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

             call igEndTabItem()
          end if ! begin tab item (bonds)
       end if

       !!!!! Labels tab !!!!!
       if (w%rep%labels%display) then
          if (iw_begintabitem("Labels##editrepatoms_labelstab")) then
             !! labels display !!

             ! label styles
             !! global options !!
             call iw_text("Global Options",highlight=.true.,alignframe=.true.)
             if (iw_button("Reset##resetglobal",sameline=.true.,danger=.true.)) then
                w%rep%labels%type = 0
                call w%rep%set_defaults(3)
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
             sz%y = iw_calcheight(min(8,nrow+1),0,.false.)
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

                call igTableSetupScrollFreeze(0, 1) ! top row always visible

                ! draw the header
                call igTableHeadersRow()
                call igTableSetColumnWidthAutoAll(igGetCurrentTable())

                ! start the clipper
                clipper = ImGuiListClipper_ImGuiListClipper()
                call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)

                ! draw the rows
                do while(ImGuiListClipper_Step(clipper))
                   call c_f_pointer(clipper,clipper_f)
                   do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd

                      ! set up the next row
                      call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                      suffix = "_" // string(i)
                      ncol = -1

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

             ! style buttons: show/hide
             if (iw_button("Show All##showalllabels")) then
                w%rep%labels%style%shown = .true.
                changed = .true.
             end if
             call iw_tooltip("Show all labels",ttshown)
             if (iw_button("Hide All##hidealllabels",sameline=.true.)) then
                w%rep%labels%style%shown = .false.
                changed = .true.
             end if
             call iw_tooltip("Hide all labels",ttshown)
             if (iw_button("Toggle Show/Hide##togglealllabels",sameline=.true.)) then
                w%rep%labels%style%shown = .not.w%rep%labels%style%shown
                changed = .true.
             end if
             call iw_tooltip("Toggle the show/hide status for all bonds",ttshown)

             call igEndTabItem()
          end if ! begin tab item (labels)
       end if

       !!!!! Polyhedra tab !!!!!
       if (w%rep%poly%display) then
          if (iw_begintabitem("Polyhedra##editrepatoms_polyhedratab")) then
             ! make sure the style is initialized and matches the current system
             if (.not.w%rep%poly%style%isinit) then
                call w%rep%poly%style%reset(w%rep)
             elseif (w%rep%poly%style%ntype /= sysc(isys)%attype_number(w%rep%poly%style%type) .or.&
                size(w%rep%poly%style%corner,1) /= sys(isys)%c%nspc) then
                call w%rep%poly%style%reset(w%rep)
             end if

             ! center type selector
             call iw_text("Centers and Corners",highlight=.true.,alignframe=.true.)
             call iw_helpermark("In the table, each row corresponds to a polyhedron center. &
                &For each center, the columns show the atom name, whether the polyhedron is shown (Show),&
                & the distance range to the corners (Min/Max), and which atomic species are allowed as corners.",sameline=.true.)
             itype_combo = 0
             if (w%rep%poly%style%type == atlisttype_nneq) itype_combo = 1
             if (w%rep%poly%style%type == atlisttype_ncel_frac) itype_combo = 2
             call iw_combo_simple("Centers##polycentertype","Species" // c_null_char //&
                "Non-equivalent atoms" // c_null_char // "Cell atoms" // c_null_char,itype_combo)
             call iw_tooltip("How to group the atoms that act as polyhedra centers",ttshown)
             newtype = atlisttype_species
             if (itype_combo == 1) newtype = atlisttype_nneq
             if (itype_combo == 2) newtype = atlisttype_ncel_frac
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
                call igTableHeadersRow()

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
             call igEndTabItem()
          end if ! begin tab item (polyhedra)
       end if
       call igEndTabBar()
    end if ! begin tab bar

    ! process transient highlighs
    if (ihighlight > 0) then
       call sysc(isys)%highlight_atoms(.true.,(/ihighlight/),highlight_type,&
          reshape(ColorHighlightScene,(/4,1/)))
    end if

  end function draw_editrep_atoms

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

    ! periodicity (evaluate first: the widgets must be drawn even if changed is
    ! already true, and .or. is allowed to short-circuit)
    ch = periodicity_widgets(w,ttshown,&
       "Number of periodic cells controlled by the +/- options in the view menu")
    changed = changed .or. ch

    !! styles
    call iw_text("Style",highlight=.true.)
    changed = changed .or. iw_checkbox("Color crystallographic axes",w%rep%uc%coloraxes)
    call iw_tooltip("Represent crystallographic axes with colors (a=red,b=green,c=blue)",ttshown)
    changed = changed .or. iw_checkbox("Hide axes in vacuum directions",w%rep%uc%vaccutsticks)
    call iw_tooltip("In systems with vacuum direction(s), do not show the unit cell in the vacuum region",ttshown)

    changed = changed .or. iw_dragfloat_real8("Radius (Å)##outer",x1=w%rep%uc%radius,speed=0.005d0,&
       min=0d0,max=5d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Radii of the unit cell edges",ttshown)

    ch = iw_coloredit("Color",rgb=w%rep%uc%rgb,sameline=.true.)
    call iw_tooltip("Color of the unit cell edges",ttshown)
    if (ch) then
       w%rep%uc%rgb = min(w%rep%uc%rgb,1._c_float)
       w%rep%uc%rgb = max(w%rep%uc%rgb,0._c_float)
       changed = .true.
    end if

    !! inner divisions
    call iw_text("Inner Divisions",highlight=.true.)
    changed = changed .or. iw_checkbox("Display inner divisions",w%rep%uc%inner)
    call iw_tooltip("Represent the inner divisions inside a supercell",ttshown)
    if (w%rep%uc%inner) then
       changed = changed .or. iw_dragfloat_real8("Radius (Å)##inner",x1=w%rep%uc%radiusinner,speed=0.005d0,&
          min=0d0,max=5d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Radii of the inner unit cell edges",ttshown)

       changed = changed .or. iw_checkbox("Use dashed lines",w%rep%uc%innerstipple)
       call iw_tooltip("Use dashed lines for the inner cell divisions",ttshown)

       if (w%rep%uc%innerstipple) then
          changed = changed .or. iw_dragfloat_real8("Dash length (Å)",x1=w%rep%uc%innersteplen,speed=0.1d0,&
             min=0.001d0,max=100d0,scale=bohrtoa,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Length of the dashed lines for the inner cell divisions (in Å)",ttshown)
       end if
    end if

    ! origin of the unit cell
    call iw_text("Origin Shift",highlight=.true.)
    changed = changed .or. iw_dragfloat_real8("##originucx",x3=w%rep%sel%origin,speed=0.001d0,decimal=5)
    call iw_tooltip("Coordinates for the origin shift of the unit cell",ttshown)

  end function draw_editrep_unitcell

  !> Draw the editrep window, cartesian axes class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_axes(w,ttshown) result(changed)
    use windows, only: win
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
    ch = iw_radiobutton("Fixed in Window",int=w%rep%axes%placement,intval=1_c_int)
    if (ch) w%rep%axes%scale_auto = .true. ! re-size the gizmo for the scene
    changed = changed .or. ch
    call iw_tooltip("Anchor the axes at a fixed position in the window",ttshown)
    changed = changed .or. iw_radiobutton("At Position",int=w%rep%axes%placement,intval=0_c_int,sameline=.true.)
    call iw_tooltip("Place the axes at the cartesian origin",ttshown)
    if (w%rep%axes%placement == 0) then
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
       if (w%rep%axes%placement == 0) then
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

  !xx! private procedures

  !> Draw the atom selection table for crystal c on representation r
  !> and return whether any item has been changed.  idparent = ID of
  !> the parent window who owns the correpsonding scene.
  !> showselection = show the selection tab columns in the tables.
  !> showdrawopts = show the atoms tab columns (draw) in the tables.
  !> dohighlight = return true if highlight has been done
  function atom_selection_widget(isys,r,showselection,showdrawopts,ihighlight,highlight_type) &
     result(changed)
    use systems, only: sys, sysc, atlisttype_species, atlisttype_nneq, atlisttype_ncel_ang,&
       atlisttype_nmol, atlisttype_ncel_frac
    use representations, only: atom_geom_style, mol_geom_style
    use utils, only: iw_text, iw_combo_simple, iw_tooltip, iw_calcheight, iw_checkbox, iw_clamp_color3,&
       iw_calcwidth, iw_button, iw_coloredit, iw_highlight_selectable, iw_dragfloat_real8,&
       iw_table_column
    use crystalmod, only: crystal
    use global, only: iunit_ang, dunit0
    use tools_io, only: string, ioj_right
    use param, only: bohrtoa
    integer, intent(in) :: isys
    type(representation), intent(inout) :: r
    logical, intent(in) :: showselection
    logical, intent(in) :: showdrawopts
    integer, intent(out) :: ihighlight
    integer, intent(out) :: highlight_type
    logical :: changed

    logical :: domol, docoord
    logical :: ch
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2, suffix
    real*8 :: x0(3)
    type(ImVec2) :: sz0, szero
    integer :: ispc, i, iz, ncol, icol
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f

    logical, save :: ttshown = .false. ! tooltip flag

    integer, parameter :: atlisttype_allowed_crys(3) = (/atlisttype_species,atlisttype_nneq,&
       atlisttype_ncel_frac/)
    integer, parameter :: atlisttype_allowed_mol(2) = (/atlisttype_species,atlisttype_ncel_ang/)
    integer, allocatable :: atlisttype_allowed(:)

    ! initialize
    ihighlight = 0
    highlight_type = 0
    szero%x = 0
    szero%y = 0
    if (showselection) then
       call iw_text("Atom Selection",highlight=.true.)
    elseif (showdrawopts) then
       call iw_text("Atom Style",highlight=.true.)
    end if
    if (sys(isys)%c%ismolecule) then
       atlisttype_allowed = atlisttype_allowed_mol
    else
       atlisttype_allowed = atlisttype_allowed_crys
    end if

    ! selector and reset
    changed = .false.

    ! occupancy sectors
    if (showdrawopts .and. sys(isys)%c%haveocc) then
       changed = changed .or. iw_checkbox("Occupancy sectors",r%atoms%occ_sectors)
       call iw_tooltip("Draw partially occupied sites as spheres with a filled "//&
          "sector proportional to the occupancy",ttshown)
       if (r%atoms%occ_sectors) then
          changed = changed .or. iw_coloredit("Vacancy color",rgb=r%atoms%occ_empty_rgb,sameline=.true.)
          call iw_tooltip("Color of the sector corresponding to a vacancy drawn on partially occupied atoms",ttshown)
       end if
    end if

    ch = sysc(isys)%attype_combo_simple("Atom types##atomtypeselection",r%atoms%style%type,&
       atlisttype_allowed,units=.false.)
    call iw_tooltip("Group atoms by these categories",ttshown)
    if (ch) then
       call r%atoms%style%reset(r)
       changed = .true.
    end if

    ! whether to do the molecule column and the coordinates
    domol = (r%atoms%style%type == atlisttype_ncel_ang)
    docoord = (r%atoms%style%type == atlisttype_nneq .or. r%atoms%style%type == atlisttype_ncel_ang .or.&
       r%atoms%style%type == atlisttype_ncel_frac)
    ncol = 3
    if (showselection) ncol = ncol + 1 ! show
    if (showdrawopts) ncol = ncol + 2 ! col, radius
    if (domol) ncol = ncol + 1 ! mol
    if (docoord) ncol = ncol + 1 ! coordinates

    ! atom style table, for atoms
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tableatomstyles" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(5,r%atoms%style%ntype)+1,0,.false.)
    if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
       icol = -1

       ! header setup
       call iw_table_column("Id",icol=icol)

       call iw_table_column("Atom",icol=icol)

       call iw_table_column("Z ",icol=icol)

       if (showselection) then
          call iw_table_column("Show",icol=icol)
       end if

       if (showdrawopts) then
          call iw_table_column("Col",icol=icol)

          call iw_table_column("Radius",icol=icol)
       end if

       if (domol) then
          call iw_table_column("Mol",icol=icol)
       end if

       if (docoord) then
          if (r%atoms%style%type == atlisttype_ncel_ang) then
             str2 = "Coordinates (Å)"
          else
             str2 = "Coordinates (fractional)"
          end if
          call iw_table_column(str2,icol=icol,flags=ImGuiTableColumnFlags_WidthStretch)
       end if

       ! draw the header
       call igTableSetupScrollFreeze(0, 1) ! top row always visible
       call igTableHeadersRow()
       call igTableSetColumnWidthAutoAll(igGetCurrentTable())

       ! start the clipper
       clipper = ImGuiListClipper_ImGuiListClipper()
       call ImGuiListClipper_Begin(clipper,r%atoms%style%ntype,-1._c_float)

       ! draw the rows
       do while(ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
             suffix = "_" // string(i)
             icol = -1

             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             ispc = sysc(isys)%attype_species(r%atoms%style%type,i)
             iz = sys(isys)%c%spc(ispc)%z

             ! id
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) then
                call iw_text(string(i),alignframe=.true.)

                ! the highlight selectable
                if (iw_highlight_selectable("##selectablemoltable" // suffix)) then
                   ihighlight = i
                   highlight_type = r%atoms%style%type
                end if
             end if

             ! name
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) call iw_text(sysc(isys)%attype_name(r%atoms%style%type,i))

             ! Z
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) call iw_text(string(iz))

             ! shown
             if (showselection) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   changed = changed .or. iw_checkbox("##tableshown" // suffix ,r%atoms%style%shown(i))
                   call iw_tooltip("Toggle display of the atom/bond/label associated to this atom",ttshown)
                end if
             end if

             ! color
             if (showdrawopts) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   ch = iw_coloredit("##tablecolor" // suffix,rgb=r%atoms%style%rgb(:,i))
                   call iw_tooltip("Atom color",ttshown)
                   if (ch) then
                      r%atoms%style%rgb(:,i) = min(r%atoms%style%rgb(:,i),1._c_float)
                      r%atoms%style%rgb(:,i) = max(r%atoms%style%rgb(:,i),0._c_float)
                      changed = .true.
                   end if
                end if

                ! radius
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   ch = iw_dragfloat_real8("##tableradius" // string(i),x1=r%atoms%style%rad(i),speed=0.01d0,&
                      min=0d0,max=5d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
                   call iw_tooltip("Radius of the sphere representing the atom",ttshown)
                   if (ch) then
                      r%atoms%style%rad(i) = max(r%atoms%style%rad(i),0d0)
                      changed = .true.
                   end if
                end if
             end if

             ! molecule
             if (domol) then
                icol = icol + 1
                ! i is a complete list index in this case
                if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%idatcelmol(1,i)))
             end if

             ! rest of info
             if (docoord) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   s = ""
                   if (r%atoms%style%type > 0) then
                      x0 = sysc(isys)%attype_coordinates(r%atoms%style%type,i)
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

    if (showselection) then
       ! style buttons: show/hide
       if (iw_button("Show All##showallatoms")) then
          r%atoms%style%shown = .true.
          changed = .true.
       end if
       call iw_tooltip("Show all atoms/bonds/labels in the system",ttshown)
       if (iw_button("Hide All##hideallatoms",sameline=.true.)) then
          r%atoms%style%shown = .false.
          changed = .true.
       end if
       call iw_tooltip("Hide all atoms/bonds/labels in the system",ttshown)
       if (iw_button("Toggle Show/Hide##toggleallatoms",sameline=.true.)) then
          do i = 1, r%atoms%style%ntype
             r%atoms%style%shown(i) = .not.r%atoms%style%shown(i)
          end do
          changed = .true.
       end if
       call iw_tooltip("Toggle the show/hide status for all atoms/bonds/labels",ttshown)
    end if

    ! molecule selection
    ! initialized and more than one molecule
    if (r%mols%style%isinit .and. r%mols%style%ntype > 1) then
       if (showselection) then
          call iw_text("Molecule Selection",highlight=.true.)
       elseif (showdrawopts) then
          call iw_text("Molecule Style",highlight=.true.)
       end if

       ncol = 2
       if (showselection) ncol = ncol + 1 ! show
       if (showdrawopts) ncol = ncol + 2 ! tint, scale
       ncol = ncol + 1 ! center of mass

       ! molecule style table, for molecules
       flags = ImGuiTableFlags_None
       flags = ior(flags,ImGuiTableFlags_Resizable)
       flags = ior(flags,ImGuiTableFlags_Reorderable)
       flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       str1="##tablemolstyles" // c_null_char
       sz0%x = 0
       sz0%y = iw_calcheight(min(5,r%mols%style%ntype)+1,0,.false.)
       if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
          icol = -1

          ! header setup
          call iw_table_column("Id",icol=icol)

          call iw_table_column("Nat",icol=icol)

          if (showselection) then
             call iw_table_column("Show",icol=icol)
          end if

          if (showdrawopts) then
             call iw_table_column("Tint",icol=icol)

             call iw_table_column("Scale",icol=icol)
          end if

          if (sys(isys)%c%ismolecule) then
             str2 = "Center of mass (Å)"
          else
             str2 = "Center of mass (fractional)"
          end if
          call iw_table_column(str2,icol=icol,flags=ImGuiTableColumnFlags_WidthStretch)

          ! draw the header
          call igTableSetupScrollFreeze(0, 1) ! top row always visible
          call igTableHeadersRow()
          call igTableSetColumnWidthAutoAll(igGetCurrentTable())

          ! start the clipper
          clipper = ImGuiListClipper_ImGuiListClipper()
          call ImGuiListClipper_Begin(clipper,r%mols%style%ntype,-1._c_float)

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
                   if (iw_highlight_selectable("##selectableatomtable" // suffix)) then
                      ihighlight = i
                      highlight_type = atlisttype_nmol
                   end if
                end if

                ! nat
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%mol(i)%nat))

                ! shown
                if (showselection) then
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      changed = changed .or. iw_checkbox("##tablemolshown" // string(i) ,r%mols%style%shown(i))
                      call iw_tooltip("Toggle display of all atoms in this molecule",ttshown)
                   end if
                end if

                ! color
                if (showdrawopts) then
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      ch = iw_coloredit("##tablemolcolor" // string(i),rgb=r%mols%style%tint_rgb(:,i))
                      call iw_tooltip("Molecule color tint",ttshown)
                      if (ch) then
                         r%mols%style%tint_rgb(:,i) = min(r%mols%style%tint_rgb(:,i),1._c_float)
                         r%mols%style%tint_rgb(:,i) = max(r%mols%style%tint_rgb(:,i),0._c_float)
                         changed = .true.
                      end if
                   end if

                   ! radius
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      changed = changed .or. iw_dragfloat_real8("##tablemolradius" // string(i),&
                         x1=r%mols%style%scale_rad(i),speed=0.005d0,min=0d0,max=5d0,decimal=3,&
                         flags=ImGuiSliderFlags_AlwaysClamp)
                      call iw_tooltip("Scale factor for the atomic radii in this molecule",ttshown)
                   end if
                end if

                ! rest of info
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

       if (showselection) then
          ! style buttons: show/hide
          if (iw_button("Show All##showallmolecules")) then
             r%mols%style%shown = .true.
             changed = .true.
          end if
          call iw_tooltip("Show all molecules in the system",ttshown)
          if (iw_button("Hide All##hideallmolecules",sameline=.true.)) then
             r%mols%style%shown = .false.
             changed = .true.
          end if
          call iw_tooltip("Hide all molecules in the system",ttshown)
          if (iw_button("Toggle Show/Hide##toggleallmolecules",sameline=.true.)) then
             do i = 1, r%mols%style%ntype
                r%mols%style%shown(i) = .not.r%mols%style%shown(i)
             end do
             changed = .true.
          end if
          call iw_tooltip("Toggle the show/hide status for all molecules",ttshown)
       end if
    end if

  end function atom_selection_widget

  !> Draw the editrep (Object) window, symmetry-elements class. Returns true if
  !> the scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_symelem(w,ttshown) result(changed)
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_combo_simple,&
       iw_button, iw_calcheight, iw_table_column
    use systems, only: sys
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch
    integer :: i, icoord, nop
    integer(c_int) :: flags
    type(ImVec2) :: sz0
    character(kind=c_char,len=:), allocatable, target :: str1

    changed = .false.

    ! refresh the symmetry-element style (snapshot + visibility) if the geometry
    ! changed since the last reset
    call w%rep%update()
    nop = w%rep%symelem%style%nop

    !! origin
    call iw_text("Origin",highlight=.true.)
    if (sys(w%isys)%c%ismolecule) then
       ! molecules: cartesian options only (coordtype 1/2), referred to the
       ! molecular center; the combo is 0-based, hence the +/-1 offset
       icoord = max(w%rep%symelem%coordtype,1) - 1
       call iw_combo_simple("Coordinates##symelemcoord","Cartesian (Å)" // c_null_char // &
          "Cartesian (bohr)" // c_null_char,icoord,changed=ch)
       if (ch .or. (icoord+1) /= w%rep%symelem%coordtype) changed = .true.
       w%rep%symelem%coordtype = icoord + 1
    else
       icoord = w%rep%symelem%coordtype
       call iw_combo_simple("Coordinates##symelemcoord","Crystallographic" // c_null_char // &
          "Cartesian (Å)" // c_null_char // "Cartesian (bohr)" // c_null_char,icoord,changed=ch)
       if (ch) changed = .true.
       w%rep%symelem%coordtype = icoord
    end if
    call iw_tooltip("Coordinate system in which the origin is given",ttshown)
    changed = changed .or. iw_dragfloat_real8("##originsymelem",x3=w%rep%symelem%origin,speed=0.001d0,decimal=5)
    call iw_tooltip("Point the symmetry elements pass through (and, for crystals, every lattice point)",ttshown)

    !! color
    call iw_text("Color",highlight=.true.)
    changed = changed .or. iw_checkbox("Custom color##symelemcustomrgb",w%rep%symelem%usecustomrgb)
    call iw_tooltip("Color all elements with a single custom color. If off, mirror/glide planes use &
       &the default color and rotation axes are colored by rotation order.",ttshown)
    if (w%rep%symelem%usecustomrgb) &
       changed = changed .or. iw_coloredit("##symelemrgb",rgb=w%rep%symelem%rgb,sameline=.true.)

    !! operations
    call iw_text("Operations",highlight=.true.)
    if (w%rep%symelem%style%isinit) then
       ! all / none / toggle
       if (iw_button("All##symelemall")) then
          w%rep%symelem%style%shown = .true.
          changed = .true.
       end if
       call iw_tooltip("Show all operations",ttshown)
       if (iw_button("None##symelemnone",sameline=.true.)) then
          w%rep%symelem%style%shown = .false.
          changed = .true.
       end if
       call iw_tooltip("Hide all operations",ttshown)
       if (iw_button("Toggle##symelemtoggle",sameline=.true.)) then
          w%rep%symelem%style%shown = .not.w%rep%symelem%style%shown
          changed = .true.
       end if
       call iw_tooltip("Toggle the operation selection",ttshown)

       ! per-operation table
       flags = ImGuiTableFlags_None
       flags = ior(flags,ImGuiTableFlags_RowBg)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       str1 = "##symelemtable" // c_null_char
       sz0%x = 0
       sz0%y = iw_calcheight(min(nop,10)+1,0,.false.)
       if (igBeginTable(c_loc(str1),3,flags,sz0,0._c_float)) then
          call iw_table_column("Show",id=0)
          call iw_table_column("#",id=1)
          call iw_table_column("Symbol",id=2)
          call igTableSetupScrollFreeze(0,1)
          call igTableHeadersRow()

          do i = 1, nop
             call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
             if (igTableSetColumnIndex(0)) then
                if (w%rep%symelem%style%kind(i) == 0) then
                   ! identity/inversion: nothing to draw, show a disabled checkbox
                   call igBeginDisabled(.true._c_bool)
                   ch = iw_checkbox("##symelemshow" // string(i),w%rep%symelem%style%shown(i))
                   call igEndDisabled()
                else
                   if (iw_checkbox("##symelemshow" // string(i),w%rep%symelem%style%shown(i))) changed = .true.
                end if
             end if
             if (igTableSetColumnIndex(1)) call iw_text(string(i))
             if (igTableSetColumnIndex(2)) call iw_text(trim(w%rep%symelem%style%label(i)))
          end do
          call igEndTable()
       end if
    end if

  end function draw_editrep_symelem

  !> Draw the editrep window, text annotations. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_text(w,ttshown) result(changed)
    use interfaces_glfw, only: glfwGetTime
    use representations, only: text_item, textpos_screen, textpos_point, textpos_atom,&
       textpos_bond
    use gui_main, only: ColorLabel_def
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_combo_simple,&
       iw_button, iw_calcheight, iw_inputtext, iw_close_button, iw_highlight_selectable, iw_radiobutton,&
       iw_table_column
    use systems, only: sys, sysc
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch, ok, ldum
    integer :: i, k, iview, isel, idel, ipl
    integer(c_int) :: flags
    type(ImVec2) :: sz0
    character(kind=c_char,len=:), allocatable, target :: str1
    type(text_item), allocatable :: taux(:)

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
          w%editrep_pick_slot = 0
       elseif (win(iview)%viewmode >= 0) then
          ! the pick finished
          i = w%editrep_pick_item
          if (win(iview)%vmdata%idx(1) > 0) then
             if (w%editrep_pick_slot == 1) then
                if (w%rep%text%t(i)%placement == textpos_bond) then
                   ! bond anchor: stage this atom and chain the second pick
                   call w%editrep_pick%stage(win(iview)%vmdata%idx(1:4))
                   w%editrep_pick_slot = 2
                   win(iview)%vmdata%idx = 0
                   call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the second atom to bond",w%id)
                else
                   ! atom anchor: complete
                   w%rep%text%t(i)%idx1 = win(iview)%vmdata%idx(1:4)
                   changed = .true.
                   w%editrep_pick_item = 0
                   w%editrep_pick_slot = 0
                end if
             else
                ! second bond atom: complete, unless it repeats the first atom
                ! (a degenerate bond; treated as a cancel)
                if (.not.w%editrep_pick%same(win(iview)%vmdata%idx(1:4))) then
                   w%rep%text%t(i)%idx1 = w%editrep_pick%idx
                   w%rep%text%t(i)%idx2 = win(iview)%vmdata%idx(1:4)
                   changed = .true.
                end if
                w%editrep_pick_item = 0
                w%editrep_pick_slot = 0
             end if
          else
             ! the user cancelled the pick
             w%editrep_pick_item = 0
             w%editrep_pick_slot = 0
          end if
          if (w%editrep_pick_item == 0) win(iview)%vmdata%idx = 0
       end if
    end if

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
       call igTableSetupScrollFreeze(0,1)
       call igTableHeadersRow()

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
             if (w%rep%text%t(i)%placement == textpos_screen) then
                str1 = "on-screen"
             elseif (w%rep%text%t(i)%placement == textpos_point) then
                str1 = "3D position"
             elseif (w%rep%text%t(i)%placement == textpos_atom) then
                str1 = "atom"
             else
                str1 = "bond"
             end if
             call iw_text(str1,alignframe=.true.)
             ldum = iw_highlight_selectable("##textsel" // string(i),clicked=ch,&
                selected=(i == w%lastselected))
             if (ch) w%lastselected = i
          end if
       end do
       call igEndTable()
    end if

    ! add button
    if (iw_button("Add##textadd")) then
       allocate(taux(w%rep%text%ntext+1))
       if (w%rep%text%ntext > 0) taux(1:w%rep%text%ntext) = w%rep%text%t(1:w%rep%text%ntext)
       call move_alloc(taux,w%rep%text%t)
       w%rep%text%ntext = w%rep%text%ntext + 1
       w%rep%text%t(w%rep%text%ntext)%str = "Text"
       w%rep%text%t(w%rep%text%ntext)%rgb = ColorLabel_def
       w%lastselected = w%rep%text%ntext
       changed = .true.
    end if
    call iw_tooltip("Add a new text",ttshown)

    ! process a deletion
    if (idel > 0) then
       do k = idel, w%rep%text%ntext-1
          w%rep%text%t(k) = w%rep%text%t(k+1)
       end do
       w%rep%text%ntext = w%rep%text%ntext - 1
       if (w%editrep_pick_item == idel) then
          ! cancel a pick pending on the deleted item
          call win(iview)%viewmode_release_forced(w%id)
          w%editrep_pick_item = 0
          w%editrep_pick_slot = 0
       elseif (w%editrep_pick_item > idel) then
          w%editrep_pick_item = w%editrep_pick_item - 1
       end if
       changed = .true.
    end if

    !! section for changing the selected text
    if (w%rep%text%ntext > 0) then
       isel = min(max(w%lastselected,1),w%rep%text%ntext)
       w%lastselected = isel

       ! the text, edited in a multiline box
       call iw_text("Text",highlight=.true.)
       if (iw_inputtext("##textstredit",bufsize=4095,texta=w%rep%text%t(isel)%str,nlines=3)) &
          changed = .true.
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
          changed = changed .or. iw_dragfloat_real8("Position##textwinpos",x2=w%rep%text%t(isel)%winpos,&
             speed=0.005d0,min=0d0,max=1d0,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Position of the text in the viewport, as fractions of the window &
             &size from the left and bottom borders",ttshown)
          changed = changed .or. iw_radiobutton("In front##textfront",bool=w%rep%text%t(isel)%infront,&
             boolval=.true.)
          call iw_tooltip("Draw the text in front of the scene",ttshown)
          changed = changed .or. iw_radiobutton("Behind##textbehind",bool=w%rep%text%t(isel)%infront,&
             boolval=.false.,sameline=.true.)
          call iw_tooltip("Draw the text behind the scene",ttshown)
       elseif (ipl == textpos_point) then
          changed = changed .or. iw_dragfloat_real8("Position##textpos",x3=w%rep%text%t(isel)%pos,&
             speed=0.001d0,decimal=5)
          if (sys(w%isys)%c%ismolecule) then
             call iw_tooltip("Position of the text (Cartesian coordinates, Å)",ttshown)
          else
             call iw_tooltip("Position of the text (fractional coordinates)",ttshown)
          end if
       elseif (ipl == textpos_atom) then
          if (iw_button("Pick atom##textpickatom",disabled=(w%editrep_pick_item > 0))) then
             w%editrep_pick_item = isel
             w%editrep_pick_slot = 1
             call w%editrep_pick%arm()
             call win(iview)%viewmode_set_forced(vm_pick_atom,&
                "Pick the atom to anchor the text to",w%id)
          end if
          call iw_tooltip("Click, then pick the anchor atom in the view window",ttshown)
          call iw_text("Anchor: " // anchor_string(w%rep%text%t(isel)%idx1),sameline=.true.)
       else ! textpos_bond
          if (iw_button("Pick bond atoms##textpickbond",disabled=(w%editrep_pick_item > 0))) then
             w%editrep_pick_item = isel
             w%editrep_pick_slot = 1
             call w%editrep_pick%arm()
             call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the first atom to bond",w%id)
          end if
          call iw_tooltip("Click, then pick the two atoms of the bond in the view window",ttshown)
          call iw_text("Anchor: " // anchor_string(w%rep%text%t(isel)%idx1) // " - " //&
             anchor_string(w%rep%text%t(isel)%idx2),sameline=.true.)
       end if

       ! 3D placements: on-screen offset from the anchor and depth toggle
       if (ipl /= textpos_screen) then
          changed = changed .or. iw_dragfloat_real8("Offset (Å)##textoffset",x2=w%rep%text%t(isel)%offset,&
             speed=0.01d0,decimal=2)
          call iw_tooltip("Offset of the text from its anchor, in on-screen (in-plane) coordinates",ttshown)
          changed = changed .or. iw_checkbox("Depth##textdepth",w%rep%text%t(isel)%depth)
          call iw_tooltip("The text is hidden by objects in front of it. If unchecked, the &
             &text is drawn on top of all objects.",ttshown)
       end if

       ! style
       call iw_text("Style",highlight=.true.)
       changed = changed .or. iw_dragfloat_real8("Size##textscale",x1=w%rep%text%t(isel)%scale,&
          speed=0.01d0,min=0.05d0,max=10d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Size of the text",ttshown)
       changed = changed .or. iw_coloredit("Color##textcolor",rgb=w%rep%text%t(isel)%rgb,sameline=.true.)
       call iw_tooltip("Color of the text",ttshown)
       changed = changed .or. iw_checkbox("Scale with zoom##textzoom",w%rep%text%t(isel)%scalewithzoom)
       call iw_tooltip("Whether the text size scales when zooming in and out, or stays &
          &at a constant on-screen size",ttshown)
    end if

  contains
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

    !> Short description of an atom anchor: name + cell index (+ lattice vector).
    function anchor_string(idx) result(s)
      integer(c_int), intent(in) :: idx(4)
      character(len=:), allocatable :: s
      s = anchor_label(w%isys,idx,"(not set)")
    end function anchor_string

  end function draw_editrep_text

  !> Draw the editrep (Object) window, measurements class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_measure(w,ttshown) result(changed)
    use representations, only: measurement_item
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_button,&
       iw_atom_button, iw_calcheight, iw_close_button, iw_intstepper, iw_highlight_selectable,&
       iw_helpermark, iw_table_column, iw_begintabitem
    use keybindings, only: get_bind_keyname, BIND_NAV_MEASURE, BIND_NAV_MEASURE_TOGGLE
    use systems, only: sys, sysc, atlisttype_ncel_frac
    use tools_io, only: string
    use param, only: pi, bohrtoa
    use interfaces_glfw, only: glfwGetTime
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer(c_int) :: flags
    integer :: idel, k, iview, i
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
          if (win(iview)%vmdata%idx(1) > 0 .and. k >= 1 .and. k <= w%rep%measure%item(i)%n) then
             w%rep%measure%item(i)%idx(:,k) = win(iview)%vmdata%idx(1:4)
             changed = .true.
          end if
          w%editrep_pick_item = 0
          w%editrep_pick_slot = 0
          win(iview)%vmdata%idx = 0
       end if
    end if

    ! usage hint
    call iw_text("Measurements",highlight=.true.)
    call iw_helpermark("Measurements are made in the view window. Select the atoms with "//&
       trim(get_bind_keyname(BIND_NAV_MEASURE))//": two for a distance, three for an angle "//&
       "(the second atom is the vertex), four for a dihedral. Then "//&
       trim(get_bind_keyname(BIND_NAV_MEASURE_TOGGLE))//" on the last atom, or on empty space, "//&
       "adds the measurement to this table; the same click again removes it.")
    ! one tab per category: a selectable table, then the options for the
    ! selected row underneath
    str1 = "##editrepmeasuretabbar" // c_null_char
    flags = ImGuiTabBarFlags_None
    if (igBeginTabBar(c_loc(str1),flags)) then
       if (iw_begintabitem("Distances##editrepmeasure_disttab")) then
          call cat_table(2,"dist")
          call item_options(2)
          call igEndTabItem()
       end if
       if (iw_begintabitem("Angles##editrepmeasure_angtab")) then
          call cat_table(3,"ang")
          call item_options(3)
          call igEndTabItem()
       end if
       if (iw_begintabitem("Dihedrals##editrepmeasure_dihtab")) then
          call cat_table(4,"dih")
          call item_options(4)
          call igEndTabItem()
       end if
       call igEndTabBar()
    end if

    ! process a deletion (global item index) and keep the selection consistent
    if (idel > 0) then
       do k = idel, w%rep%measure%nitem-1
          w%rep%measure%item(k) = w%rep%measure%item(k+1)
       end do
       w%rep%measure%nitem = w%rep%measure%nitem - 1
       if (w%rep%measure%isel == idel) then
          w%rep%measure%isel = 0
       elseif (w%rep%measure%isel > idel) then
          w%rep%measure%isel = w%rep%measure%isel - 1
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

    !> Draw the table of measurement items of the given category (ncat
    !> = number of atoms: 2 distances, 3 angles, 4 dihedrals). Every
    !> atom is a colored button that commands the view into pick mode
    !> to replace it. Every row has the value and is selectable (sets
    !> that editor).
    subroutine cat_table(ncat,tabidn)
      integer, intent(in) :: ncat
      character(len=*), intent(in) :: tabidn

      integer(c_int) :: tflags
      integer :: i, k, nrow, ncol, icvalue
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
         call igTableSetupScrollFreeze(0,1)
         call igTableHeadersRow()

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
         call igEndTable()
      end if
    end subroutine cat_table

    !> Draw the per-item options for the selected measurement, underneath its
    !> table, but only when the selected item belongs to this category (ncat).
    subroutine item_options(ncat)
      integer, intent(in) :: ncat

      integer :: is, j
      integer(c_int) :: idec

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
         ! line style (distances): solid or dashed, with the dash length
         if (ncat == 2) then
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
         ! apply these style options to every measurement in this tab
         if (iw_button("Apply to All##measureitemapply",danger=.true.)) then
            do j = 1, w%rep%measure%nitem
               if (w%rep%measure%item(j)%n /= ncat .or. j == is) cycle
               call w%rep%measure%item(j)%copy_style(it)
            end do
            changed = .true.
         end if
         call iw_tooltip("Apply these style options to all measurements in this tab",ttshown)
      end associate
    end subroutine item_options

    !> Draw a button for atom islot of measurement item iitem (cell
    !> atom id + lattice vector), tinted with its color in the parent
    !> view. Label = species name + cell id ("?" if the atom is
    !> stale). Clicking it commands the parent view into pick mode.
    subroutine atom_button(iitem,islot,idn)
      integer, intent(in) :: iitem, islot
      character(len=*), intent(in) :: idn

      integer :: cid
      integer(c_int) :: idxfull(4)
      real(c_float) :: rgb(3)
      logical :: havergb, clicked
      character(len=:), allocatable :: lbl

      idxfull = w%rep%measure%item(iitem)%idx(:,islot)
      cid = idxfull(1)
      havergb = .false.
      rgb = 0._c_float
      if (cid >= 1 .and. cid <= sys(w%isys)%c%ncel) then
         lbl = anchor_label(w%isys,idxfull,"?",species=.true.)
         havergb = atom_view_rgb(iview,w%isys,atlisttype_ncel_frac,cid,rgb)
      else
         lbl = "?"
      end if

      clicked = iw_atom_button(lbl // idn,rgb,havergb=havergb,&
         disabled=(w%editrep_pick_item > 0))
      call iw_tooltip("Atom " // lbl // " (click to pick a replacement in the view)",ttshown)

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
      integer :: ia, idd
      logical :: okk
      okk = .true.
      do ia = 1, item%n
         idd = item%idx(1,ia)
         if (idd < 1 .or. idd > sys(w%isys)%c%ncel) then
            okk = .false.
            exit
         end if
         f(:,ia) = sys(w%isys)%c%atcel(idd)%x + item%idx(2:4,ia)
      end do
      if (.not.okk) then
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

  !> Draw the editrep window, isosurface class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_isosurface(w,ttshown) result(changed)
    use systems, only: sys
    use representations, only: iso_grid_size, iso_isgridfield, iso_level_custom,&
       iso_npts_custom_min, iso_npts_custom_max, iso_level_optstr_custom, iso_defaultlevel,&
       iso_region_cell, iso_region_frac, iso_region_ortho, iso_region_parallel,&
       iso_region_simplebox, iso_region_cube, iso_region_name, iso_region_modes_cry,&
       iso_region_modes_mol, iso_region_to_box, iso_region_seed,&
       iso_region_point_from_cart, iso_estimate_cost
    use utils, only: iw_text, iw_tooltip, iw_coloredit, iw_dragfloat_real8,&
       iw_calcwidth, iw_combo_simple, iw_button, iw_intstepper, duration_string
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: i, isys, iview, nstage(3), nshow(3), nrow, nmode, istat, ilevprev, ipad
    integer :: rmodes(size(iso_region_modes_mol))
    integer(c_int) :: ncus(3)
    real*8 :: box(3,0:3), prev0(3), prevv(3,3), flo(3), fhi(3), xpick(3)
    logical :: ch, ldum, goodf, isgrid, navail, okbox, capped, lapply
    logical(c_bool) :: is_selected
    real(c_float) :: rgba(4)
    real*8 :: speed
    character(kind=c_char,len=:), allocatable, target :: str1, str2
    type(ImVec2) :: szero

    ! initialize (the caller has already checked the anchor view and its
    ! scene are alive; the edited object lives in that scene, so the
    ! transient preview and the cell counts use it, not the main view's)
    changed = .false.
    isys = w%isys
    iview = w%anchor_view()
    szero%x = 0
    szero%y = 0

    ! handle a pending region-point pick commanded to the anchor view
    ! (armed by the Pick buttons in the Region section). The pick is
    ! self-validating: if the staged region mode no longer matches the
    ! one it was armed under -- any field change also resets the mode --
    ! or the field is gone, the picked point has no row to land in and
    ! the pick is cancelled, wherever the change came from.
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
    str1 = "##isofieldcombo" // c_null_char
    if (goodf) then
       str2 = string(w%rep%iso%ifield) // ": " // trim(sys(isys)%f(w%rep%iso%ifield)%name) // c_null_char
    else
       str2 = "<field not available>" // c_null_char
    end if
    call igSameLine(0._c_float,-1._c_float)
    call igSetNextItemWidth(iw_calcwidth(30,1))
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 0, sys(isys)%nf
          if (.not.sys(isys)%f(i)%isinit) cycle
          is_selected = (w%rep%iso%ifield == i)
          str2 = string(i) // ": " // trim(sys(isys)%f(i)%name) // c_null_char
          if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) then
             if (w%rep%iso%ifield /= i) then
                ! set the field with its default isovalue and grid level,
                ! applied immediately
                call w%rep%iso%set_field(isys,i)
                changed = .true.
             end if
          end if
          if (is_selected) &
             call igSetItemDefaultFocus()
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Field whose isosurface is displayed (selecting a field&
       & resets the isovalue to a default for that field)",ttshown)
    if (.not.goodf) then
       call iw_text("The selected field is not available in this system",danger=.true.)
       return
    end if

    ! field policy repair, before the sections below read it: the field
    ! may have stopped being a grid (e.g. a different field reloaded into
    ! the same slot), which invalidates the native level
    isgrid = iso_isgridfield(isys,w%rep%iso%ifield)
    if (.not.isgrid .and. w%rep%iso%ilevel == 0) then
       call w%rep%iso%set_field(isys,w%rep%iso%ifield)
       changed = .true.
    end if

    ! region section: the box over which the isosurface is calculated;
    ! staged like the level, committed by the same Calculate grid button.
    ! The modes offered depend on the system kind (iso_region_modes_*:
    ! no cell fractions for molecules, whose cell is an artifact; the
    ! molecule-centered simple box and cube are molecule-only). The
    ! staged box is shown in the view while this editor is open (posted
    ! at the end of the section).
    call iw_text("Region",highlight=.true.)
    if (sys(isys)%c%ismolecule) then
       nmode = size(iso_region_modes_mol)
       rmodes(1:nmode) = iso_region_modes_mol
    else
       nmode = size(iso_region_modes_cry)
       rmodes(1:nmode) = iso_region_modes_cry
    end if
    ! reset a staged mode this system kind does not offer
    if (all(rmodes(1:nmode) /= w%rep%iso%iregion)) then
       w%rep%iso%iregion = iso_region_cell
       call iso_region_seed(isys,w%rep%iso%iregion,w%rep%iso%rgn_x)
    end if
    str1 = "Mode##isoregion" // c_null_char
    str2 = trim(iso_region_name(w%rep%iso%iregion)) // c_null_char
    call igSetNextItemWidth(iw_calcwidth(len(iso_region_name)+4,0))
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nmode
          is_selected = (w%rep%iso%iregion == rmodes(i))
          str2 = trim(iso_region_name(rmodes(i))) // c_null_char
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
    call iw_tooltip("Region over which the isosurface is calculated: the whole&
       & unit cell, a cell-aligned box between two fractional points (crystals),&
       & an axis-aligned box between two Cartesian corners, an arbitrary&
       & parallelepiped given by its origin and the endpoints of its three&
       & edges, an axis-aligned box given by its center and half-lengths&
       & (molecules), or a cube given by its center and half-length (molecules).&
       & The region may extend beyond the unit cell.",ttshown)

    ! periodicity of the whole-cell isosurface: number of unit cells the
    ! mesh is replicated over, same as the other representations. Shown
    ! when the built mesh is periodic (whole cell of a crystal) -- the
    ! only case that replicates -- so the control is visible exactly
    ! while its state is live, regardless of the staged region. It
    ! applies immediately, unlike the staged region widgets around it.
    ! (Evaluate first: the widgets must be drawn even if changed is
    ! already true, and .or. is allowed to short-circuit.)
    if (w%rep%iso%per0_built) then
       ch = periodicity_widgets(w,ttshown,&
          "Number of periodic cells controlled by the +/- options in the 'Scene' button of the view window",&
          highlight=.false.)
       changed = changed .or. ch
    end if

    ! staged region coordinates (results discarded, Calculate grid commits);
    ! every point row gets a Pick button that fills it from the view (the
    ! half-length rows are lengths, not points)
    if (w%rep%iso%iregion == iso_region_simplebox .or.&
       w%rep%iso%iregion == iso_region_cube) then
       ldum = iw_dragfloat_real8("x0 (Å)##isorgn0",x3=w%rep%iso%rgn_x(:,0),&
          speed=0.05d0,decimal=3)
       call iw_tooltip(region_point_tooltip(w%rep%iso%iregion,0),ttshown)
       call region_pick_button(0)
       if (w%rep%iso%iregion == iso_region_cube) then
          ldum = iw_dragfloat_real8("Half-length (Å)##isorgn1",x1=w%rep%iso%rgn_x(1,1),&
             speed=0.05d0,decimal=3,min=0d0,flags=ImGuiSliderFlags_AlwaysClamp)
          ! the single half-length is stored replicated across rgn_x(:,1)
          w%rep%iso%rgn_x(2:3,1) = w%rep%iso%rgn_x(1,1)
       else
          ldum = iw_dragfloat_real8("Half-lengths (Å)##isorgn1",x3=w%rep%iso%rgn_x(:,1),&
             speed=0.05d0,decimal=3,min=0d0,flags=ImGuiSliderFlags_AlwaysClamp)
       end if
       call iw_tooltip(region_point_tooltip(w%rep%iso%iregion,1),ttshown)
    elseif (w%rep%iso%iregion /= iso_region_cell) then
       nrow = 1 ! last drag row: the far corner, or the third edge endpoint
       if (w%rep%iso%iregion == iso_region_parallel) nrow = 3
       do i = 0, nrow
          if (w%rep%iso%iregion == iso_region_frac) then
             str1 = "x" // string(i) // "##isorgn" // string(i)
             ldum = iw_dragfloat_real8(str1,x3=w%rep%iso%rgn_x(:,i),speed=0.01d0,decimal=4)
          else
             str1 = "x" // string(i) // " (Å)##isorgn" // string(i)
             ldum = iw_dragfloat_real8(str1,x3=w%rep%iso%rgn_x(:,i),speed=0.05d0,decimal=3)
          end if
          call iw_tooltip(region_point_tooltip(w%rep%iso%iregion,i),ttshown)
          call region_pick_button(i)
       end do
    end if

    ! the staged region box; a degenerate one blocks the grid below
    call iso_region_to_box(isys,w%rep%iso%iregion,w%rep%iso%rgn_x,box,okbox)
    if (.not.okbox) &
       call iw_text("The region is degenerate (zero volume)",danger=.true.,wrap=.true.)

    ! grid section: coarseness of the grid that supports the isosurface
    ! over the region above, the resulting dimensions, and the buttons
    ! that commit the staged grid and region
    call iw_text("Grid",highlight=.true.)

    ! coarseness level combo ("Native grid" only for grid fields over the
    ! whole cell)
    navail = isgrid .and. (w%rep%iso%iregion == iso_region_cell)
    if (.not.navail .and. w%rep%iso%ilevel == 0) &
       w%rep%iso%ilevel = iso_defaultlevel
    if (navail) then
       str1 = "Native grid" // c_null_char // iso_level_optstr_custom
    else
       str1 = iso_level_optstr_custom
    end if
    ilevprev = w%rep%iso%ilevel
    call iw_combo_simple("Level##isolevel",str1,w%rep%iso%ilevel,startsatone=.not.navail,&
       changed=ch)
    call iw_tooltip("Coarseness of the grid supporting the isosurface: the native grid&
       & of the field, a named density of points per angstrom, or a custom number of&
       & points",ttshown)

    ! custom level: the dimensions are given by the user, seeded with
    ! those of the level that was showing when Custom was selected
    ! (staged: committed by Calculate grid)
    if (w%rep%iso%ilevel == iso_level_custom) then
       if (ch .and. ilevprev /= iso_level_custom) then
          if (ilevprev == 0) then
             w%rep%iso%nptscustom = sys(isys)%f(w%rep%iso%ifield)%grid%n
          elseif (okbox) then
             w%rep%iso%nptscustom = level_grid_size(ilevprev)
          end if
       end if
       ncus = w%rep%iso%nptscustom
       ipad = ceiling(log10(max(maxval(ncus),1) + 0.1))
       ldum = iw_intstepper("n1##isonpts1",ncus(1),label="n1:",minval=int(iso_npts_custom_min,c_int),&
          maxval=int(iso_npts_custom_max,c_int),ndigit=ipad,notlive=.true.)
       ldum = iw_intstepper("n2##isonpts2",ncus(2),label="n2:",minval=int(iso_npts_custom_min,c_int),&
          maxval=int(iso_npts_custom_max,c_int),ndigit=ipad,notlive=.true.,sameline=.true.)
       ldum = iw_intstepper("n3##isonpts3",ncus(3),label="n3:",minval=int(iso_npts_custom_min,c_int),&
          maxval=int(iso_npts_custom_max,c_int),ndigit=ipad,notlive=.true.,sameline=.true.)
       call iw_tooltip("Number of grid points along each axis of the sampled region",ttshown)
       w%rep%iso%nptscustom = ncus
    end if

    ! staged grid dimensions
    lapply = .false.
    nstage = 0
    if (okbox) then
       nstage = level_grid_size(w%rep%iso%ilevel)
       nshow = nstage
       str2 = ""
       if (all(nstage == 0)) then
          nshow = sys(isys)%f(w%rep%iso%ifield)%grid%n
          str2 = " (native)"
       elseif (capped) then
          str2 = " (capped)"
       end if
       ! the custom level shows its dimensions in the steppers above, so
       ! the line is only needed there when the cap trimmed them
       if (w%rep%iso%ilevel /= iso_level_custom .or. capped) then
          call iw_text("Grid: " // string(nshow(1)) // " x " // string(nshow(2)) // " x " //&
             string(nshow(3)) // str2)
          call iw_tooltip("Dimensions of the grid that will support the isosurface",ttshown)
       end if
       lapply = .not.w%rep%iso%grid_isapplied(nstage,w%rep%iso%iregion,w%rep%iso%rgn_x)
    end if
    if (iw_button("Calculate grid",danger=.true.,disabled=.not.lapply)) then
       call w%rep%iso%apply_grid(nstage,w%rep%iso%iregion,w%rep%iso%rgn_x)
       changed = .true.
    end if
    call iw_tooltip("Use the selected grid and region for the isosurface (may require&
       & an expensive recalculation of the field samples)",ttshown)
    ! cost estimate for the staged options (nstage is all-zero for the
    ! native grid or a degenerate region: nothing to sample). The
    ! benchmark measures seconds per point, so the displayed total
    ! follows the staged options live
    if (iw_button("Estimate cost",sameline=.true.,disabled=all(nstage == 0))) &
       w%rep%iso%costest = iso_estimate_cost(isys,w%rep%iso%ifield,w%rep%iso%iregion,&
          w%rep%iso%rgn_x,nstage)
    call iw_tooltip("Estimate the time needed to sample the field with the selected&
       & options (evaluates the field at random points of the sampled box and&
       & extrapolates to the full grid)",ttshown)
    if (w%rep%iso%costest >= 0d0 .and. any(nstage /= 0)) &
       call iw_text("~" // duration_string(w%rep%iso%costest * product(real(nstage,8))),&
          sameline=.true.)
    ! red state warnings; a degenerate staged region already has its own
    ! message and a disabled button, so it is excluded
    if (okbox .and. .not.w%rep%iso%isgenerated(isys)) then
       call iw_text("The isosurface has not been generated yet (press Calculate grid)",&
          danger=.true.,wrap=.true.)
    elseif (lapply) then ! false for a degenerate region
       call iw_text("The grid and region options do not match the displayed isosurface&
          & (press Calculate grid)",danger=.true.,wrap=.true.)
    end if
    if (w%rep%iso%outdomain) &
       call iw_text("The field could not be evaluated in part of the region (zero there)",&
          danger=.true.,wrap=.true.)

    ! show the staged region in the view for as long as this editor is
    ! open (a transient box, rearmed every frame) -- except for a
    ! crystal in whole-cell mode, where the region is the cell itself
    ! and the box would only duplicate the unit-cell representation. In
    ! whole-cell mode with a molecular grid field the box is the valid
    ! window of the grid's own domain (matching the grid dimensions
    ! above), not the artificial cell; otherwise it is the staged
    ! region box. show_transient_box takes the user frame (with molx0
    ! for molecules); box and get_domain are cell-frame.
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
       call win(iview)%sc%show_transient_box(w%id,1,prev0,prevv,region_edgerad,&
          region_rgb,region_alpha)
    end if

    ! isovalue, with a drag speed proportional to the current value;
    ! notlive so the re-triangulation runs on commit, not every drag frame.
    ! The selectable values are capped by the range of the data backing the
    ! mesh (stamped by the last build; unclamped until the first build).
    call iw_text("Isosurface",highlight=.true.)
    speed = max(0.01d0 * abs(w%rep%iso%isoval(1)),1d-4)
    if (w%rep%iso%frange(1) <= w%rep%iso%frange(2)) then
       changed = changed .or. iw_dragfloat_real8("Isovalue",x1=w%rep%iso%isoval(1),speed=speed,&
          decimal=6,notlive=.true.,min=w%rep%iso%frange(1),max=w%rep%iso%frange(2),&
          flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Value of the field on the isosurface, in atomic units&
          & (grid range: " // string(w%rep%iso%frange(1),'e',decimal=4) // " to " //&
          string(w%rep%iso%frange(2),'e',decimal=4) // ")",ttshown)
    else
       changed = changed .or. iw_dragfloat_real8("Isovalue",x1=w%rep%iso%isoval(1),speed=speed,&
          decimal=6,notlive=.true.)
       call iw_tooltip("Value of the field on the isosurface (atomic units)",ttshown)
    end if

    ! color and opacity
    rgba(1:3) = w%rep%iso%rgb(:,1)
    rgba(4) = w%rep%iso%alpha(1)
    ch = iw_coloredit("Color",rgba=rgba,sameline=.true.)
    call iw_tooltip("Color and opacity of the isosurface",ttshown)
    if (ch) then
       w%rep%iso%rgb(:,1) = rgba(1:3)
       w%rep%iso%alpha(1) = rgba(4)
       changed = .true.
    end if

  contains
    !> Dimensions of the grid that coarseness level ilev gives over the
    !> staged region (host-associated box, valid only when okbox); also
    !> reports through capped whether the total-points cap trimmed them.
    function level_grid_size(ilev) result(nn)
      integer, intent(in) :: ilev
      integer :: nn(3)

      if (w%rep%iso%iregion == iso_region_cell) then
         nn = iso_grid_size(isys,ilev,w%rep%iso%nptscustom,capped,ifield=w%rep%iso%ifield)
      else
         nn = iso_grid_size(isys,ilev,w%rep%iso%nptscustom,capped,box=box)
      end if

    end function level_grid_size

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

  !> Draw the periodicity widgets shared by the atoms, unit-cell, and
  !> isosurface object editors: the radio buttons for the periodicity type
  !> and, when the type is manual, the a/b/c cell-number steppers plus a
  !> Reset button. tooltip_automatic = tooltip for the "Automatic" button,
  !> which points at a different control in each editor. ttshown = the
  !> tooltip flag. highlight = draw the "Periodicity" label as a section
  !> header (default; false where it is an item inside another section).
  !> Returns true if the periodicity changed.
  function periodicity_widgets(w,ttshown,tooltip_automatic,highlight) result(changed)
    use utils, only: iw_text, iw_tooltip, iw_radiobutton, iw_button, iw_intstepper
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    character(len=*,kind=c_char), intent(in) :: tooltip_automatic
    logical, intent(in), optional :: highlight
    logical :: changed

    logical :: ldum, highlight_
    integer :: ipad
    integer(c_int) :: nc(3)

    highlight_ = .true.
    if (present(highlight)) highlight_ = highlight

    changed = .false.
    call iw_text("Periodicity",highlight=highlight_,alignframe=.true.)

    ! radio buttons for the periodicity type
    changed = changed .or. iw_radiobutton("None",int=w%rep%sel%pertype,intval=0_c_int,sameline=.true.)
    call iw_tooltip("This object is represented only in the main cell and not repeated by translation",ttshown)
    changed = changed .or. iw_radiobutton("Automatic",int=w%rep%sel%pertype,intval=1_c_int,sameline=.true.)
    call iw_tooltip(tooltip_automatic,ttshown)
    changed = changed .or. iw_radiobutton("Manual",int=w%rep%sel%pertype,intval=2_c_int,sameline=.true.)
    call iw_tooltip("Manually set the number of periodic cells",ttshown)

    ! number of periodic cells, if manual (all three fields share a digit width)
    if (w%rep%sel%pertype == 2_c_int) then
       ipad = ceiling(log10(max(maxval(w%rep%sel%ncell),1) + 0.1))
       nc = w%rep%sel%ncell
       ldum = iw_intstepper("aaxis",nc(1),label="a:",minval=1_c_int,ndigit=ipad,notlive=.true.)
       ldum = iw_intstepper("baxis",nc(2),label="b:",minval=1_c_int,ndigit=ipad,notlive=.true.,sameline=.true.)
       ldum = iw_intstepper("caxis",nc(3),label="c:",minval=1_c_int,ndigit=ipad,notlive=.true.,sameline=.true.)
       if (any(nc /= w%rep%sel%ncell)) then
          w%rep%sel%ncell = nc
          changed = .true.
       end if

       if (iw_button("Reset",sameline=.true.)) then
          w%rep%sel%ncell = 1
          changed = .true.
       end if
    end if

  end function periodicity_widgets

end submodule editrep
