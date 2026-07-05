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

contains

  !> Update tasks for the edit representation window, before the
  !> window is created.
  module subroutine update_editrep(w)
    use windows, only: nwin, win, wintype_view
    use systems, only: sys_init, ok_system
    class(window), intent(inout), target :: w

    integer :: isys
    logical :: doquit

    ! check the system and representation are still active
    isys = w%isys
    doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = .not.associated(w%rep)
    if (.not.doquit) doquit = .not.w%rep%isinit
    if (.not.doquit) doquit = (w%rep%type <= 0)
    if (.not.doquit) doquit = .not.(w%idparent > 0 .and. w%idparent <= nwin)
    if (.not.doquit) doquit = .not.(win(w%idparent)%isinit)
    if (.not.doquit) doquit = win(w%idparent)%type /= wintype_view

    ! if they aren't, quit the window
    if (doquit) call w%end()

  end subroutine update_editrep

  !> Draw the edit represenatation window.
  module subroutine draw_editrep(w)
    use representations, only: representation, reptype_atoms, reptype_unitcell, reptype_axes,&
       reptype_symelem
    use windows, only: nwin, win, wintype_view
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG,&
       BIND_CLOSE_ALL_DIALOGS
    use systems, only: sysc, sys_init, ok_system
    use gui_main, only: g
    use utils, only: iw_text, iw_tooltip, iw_button, iw_calcwidth,&
       iw_calcheight, iw_checkbox, iw_inputtext
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: isys
    logical :: doquit, ok, ldum
    logical :: changed
    type(ImVec2) :: szavail

    logical, save :: ttshown = .false. ! tooltip flag

    ! check the system and representation are still active
    isys = w%isys
    doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = .not.associated(w%rep)
    if (.not.doquit) doquit = .not.w%rep%isinit
    if (.not.doquit) doquit = (w%rep%type <= 0)
    if (.not.doquit) doquit = .not.(w%idparent > 0 .and. w%idparent <= nwin)
    if (.not.doquit) doquit = .not.(win(w%idparent)%isinit)
    if (.not.doquit) doquit = .not.associated(win(w%idparent)%sc)
    if (.not.doquit) doquit = win(w%idparent)%type /= wintype_view
    if (.not.doquit) doquit = (win(w%idparent)%view_selected /= isys)

    if (.not.doquit) then
       ! whether the rep has changed
       changed = .false.

       ! system
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! name block (the representation type is fixed at creation and cannot be
       ! changed here)
       call igAlignTextToFramePadding()
       call iw_text("Name",highlight=.true.)
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
       end if

       ! rebuild draw lists if necessary
       if (changed) win(w%idparent)%sc%forcebuildlists = .true.

       ! right-align and bottom-align for the rest of the contents
       call igGetContentRegionAvail(szavail)
       call igSetCursorPosX(iw_calcwidth(10,2,from_end=.true.) - g%Style%ScrollbarSize)
       if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
          call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

       ! reset button
       if (iw_button("Reset",danger=.true.)) then
          call w%rep%set_defaults(0)
          win(w%idparent)%sc%forcebuildlists = .true.
       end if

       ! close button
       ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
       ok = ok .or. iw_button("Close",sameline=.true.)
       if (ok) doquit = .true.
    end if

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_editrep

  !> Draw the editrep (Object) window, atoms class. Returns true if
  !> the scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_atoms(w,ttshown) result(changed)
    use representations, only: representation
    use systems, only: sys, sysc, atlisttype_species, atlisttype_ncel_ang, atlisttype_nmol,&
       atlisttype_nneq, atlisttype_ncel_frac
    use gui_main, only: g, ColorHighlightScene, ColorElement
    use tools_io, only: string
    use utils, only: iw_text, iw_tooltip, iw_helpermark, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_radiobutton, iw_calcheight, iw_clamp_color3, iw_checkbox, iw_coloredit,&
       iw_highlight_selectable, iw_dragfloat_real8, iw_inputtext, iw_inputint
    use param, only: atmcov, atmvdw, atmcov0, newline, jmlcol, jmlcol2, bohrtoa
    use global, only: bondfactor_def, bonddelta_def
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    integer :: ispc, isys, iz, ipad
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3, suffix
    real*8 :: x0(3)
    logical :: ch, ldum
    integer(c_int) :: nc(3), lst, flags, nspcpair
    real(c_float) :: sqw
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
       str1 = "Selection##editrepatoms_selectiontab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
          ! filter
          call igAlignTextToFramePadding()
          call iw_text("Filter",highlight=.true.)
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

          ! periodicity
          if (.not.sys(isys)%c%ismolecule) then
             call igAlignTextToFramePadding()
             call iw_text("Periodicity",highlight=.true.)

             ! radio buttons for the periodicity type
             changed = changed .or. iw_radiobutton("None",int=w%rep%sel%pertype,intval=0_c_int,sameline=.true.)
             call iw_tooltip("This object is represented only in the main cell and not repeated by translation",ttshown)
             changed = changed .or. iw_radiobutton("Automatic",int=w%rep%sel%pertype,intval=1_c_int,sameline=.true.)
             call iw_tooltip("Number of periodic cells controlled by the +/- options &
                &in the 'Scene' button of the view window",ttshown)
             changed = changed .or. iw_radiobutton("Manual",int=w%rep%sel%pertype,intval=2_c_int,sameline=.true.)
             call iw_tooltip("Manually set the number of periodic cells",ttshown)

             ! number of periodic cells, if manual
             if (w%rep%sel%pertype == 2_c_int) then
                ! calculate widths
                ipad = ceiling(log10(max(maxval(w%rep%sel%ncell),1) + 0.1))
                sqw = max(iw_calcwidth(1,1),igGetTextLineHeightWithSpacing())
                call igPushItemWidth(sqw)

                nc = w%rep%sel%ncell
                call igAlignTextToFramePadding()
                call iw_text("a:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##aaxis")) w%rep%sel%ncell(1) = max(w%rep%sel%ncell(1)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                ldum = iw_inputint("##aaxis",w%rep%sel%ncell(1),width=ipad,notlive=.true.)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##aaxis")) w%rep%sel%ncell(1) = w%rep%sel%ncell(1)+1

                call igSameLine(0._c_float,-1._c_float)
                call iw_text("b:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##baxis")) w%rep%sel%ncell(2) = max(w%rep%sel%ncell(2)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                ldum = iw_inputint("##baxis",w%rep%sel%ncell(2),width=ipad,notlive=.true.)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##baxis")) w%rep%sel%ncell(2) = w%rep%sel%ncell(2)+1

                call igSameLine(0._c_float,-1._c_float)
                call iw_text("c:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##caxis")) w%rep%sel%ncell(3) = max(w%rep%sel%ncell(3)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                ldum = iw_inputint("##caxis",w%rep%sel%ncell(3),width=ipad,notlive=.true.)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##caxis")) w%rep%sel%ncell(3) = w%rep%sel%ncell(3)+1
                w%rep%sel%ncell = max(w%rep%sel%ncell,1)
                if (any(nc /= w%rep%sel%ncell)) changed = .true.
                call igPopItemWidth()

                if (iw_button("Reset",sameline=.true.)) then
                   w%rep%sel%ncell = 1
                   changed = .true.
                end if
             end if

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
          str1 = "Atoms##editrepatoms_atomstab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
             ! global options for atoms
             call igAlignTextToFramePadding()
             call iw_text("Global Options",highlight=.true.)
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
          str1 = "Bonds##editrepatoms_bondstab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
             !! bonds display !!

             !! global options !!
             call igAlignTextToFramePadding()
             call iw_text("Global Options",highlight=.true.)
             if (iw_button("Reset##resetglobal",sameline=.true.,danger=.true.)) then
                call w%rep%set_defaults(2)
                changed = .true.
             end if
             call iw_tooltip("Reset to the covalent bonding for this system and the default settings")

             ! rest of the options (record changes)
             ch = .false.
             call igAlignTextToFramePadding()
             call iw_text("Style")
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
             call igAlignTextToFramePadding()
             call iw_text("Color")
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
             call igAlignTextToFramePadding()
             call iw_text("Classify H-bonds")
             ch = ch .or. iw_checkbox("##hbondclassify",w%rep%bonds%hbond_classify,sameline=.true.)
             call iw_tooltip("Color each contact by its Jeffrey-Steiner hydrogen-bond strength &
                &(strong/moderate/weak), combining the H...A distance and the D-H...A angle information.",ttshown)

             if (w%rep%bonds%hbond_classify) then
                ch = ch .or. iw_coloredit("Strong##hbstrongcolor",rgb=w%rep%bonds%hbond_rgb(:,1))
                ch = ch .or. iw_coloredit("Moderate##hbmodcolor",rgb=w%rep%bonds%hbond_rgb(:,2),sameline=.true.)
                ch = ch .or. iw_coloredit("Weak##hbweakcolor",rgb=w%rep%bonds%hbond_rgb(:,3),sameline=.true.)

                call igAlignTextToFramePadding()
                call iw_text("H...A Distance (Å)")
                ch = ch .or. iw_dragfloat_real8("strong|moderate##hbdist1",x1=w%rep%bonds%hbond_dist(1),speed=0.01d0,&
                   min=0d0,max=5d0,scale=bohrtoa,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
                call iw_tooltip("H...A distance class boundaries",ttshown)
                ch = ch .or. iw_dragfloat_real8("moderate|weak##hbdist2",x1=w%rep%bonds%hbond_dist(2),speed=0.01d0,&
                   min=0d0,max=5d0,scale=bohrtoa,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
                call iw_tooltip("H...A distance class boundaries",ttshown)

                call igAlignTextToFramePadding()
                call iw_text("D-H...A Angle (°)")
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
                str2 = "Atom 1" // c_null_char
                flags = ImGuiTableColumnFlags_WidthFixed
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_sp1)

                str2 = "Atom 2" // c_null_char
                flags = ImGuiTableColumnFlags_WidthFixed
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_sp2)

                str2 = "Show" // c_null_char
                flags = ImGuiTableColumnFlags_WidthFixed
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_shown)

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
                         call igAlignTextToFramePadding()
                         call iw_text(trim(sys(isys)%c%spc(i)%name))
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
                   str2 = "Atom" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,0)
                   str2 = "Z" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,1)
                   str2 = "Radius (Å)" // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,2)
                   call igTableSetupScrollFreeze(0,1)
                   call igTableHeadersRow()

                   do i = 1, sys(isys)%c%nspc
                      iz = sys(isys)%c%spc(i)%z
                      if (iz <= 0) cycle
                      call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                      if (igTableSetColumnIndex(0)) then
                         call igAlignTextToFramePadding()
                         call iw_text(trim(sys(isys)%c%spc(i)%name))
                      end if
                      if (igTableSetColumnIndex(1)) then
                         call igAlignTextToFramePadding()
                         call iw_text(string(iz))
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
          str1 = "Labels##editrepatoms_labelstab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
             !! labels display !!

             ! label styles
             !! global options !!
             call igAlignTextToFramePadding()
             call iw_text("Global Options",highlight=.true.)
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
                str2 = "Id" // c_null_char
                ncol = ncol + 1
                flags = ImGuiTableColumnFlags_WidthFixed
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)

                if (intable /= atlisttype_nmol) then
                   str2 = "Atom" // c_null_char
                   ncol = ncol + 1
                   flags = ImGuiTableColumnFlags_WidthFixed
                   call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)

                   str2 = "Z " // c_null_char
                   ncol = ncol + 1
                   flags = ImGuiTableColumnFlags_WidthFixed
                   call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)
                end if

                str2 = "Show" // c_null_char
                ncol = ncol + 1
                flags = ImGuiTableColumnFlags_WidthFixed
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)

                str2 = "Text" // c_null_char
                ncol = ncol + 1
                flags = ImGuiTableColumnFlags_WidthStretch
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)

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
                         call igAlignTextToFramePadding()
                         call iw_text(string(i))

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
          str1 = "Polyhedra##editrepatoms_polyhedratab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
             ! make sure the style is initialized and matches the current system
             if (.not.w%rep%poly%style%isinit) then
                call w%rep%poly%style%reset(w%rep)
             elseif (w%rep%poly%style%ntype /= sysc(isys)%attype_number(w%rep%poly%style%type) .or.&
                size(w%rep%poly%style%corner,1) /= sys(isys)%c%nspc) then
                call w%rep%poly%style%reset(w%rep)
             end if

             ! center type selector
             call igAlignTextToFramePadding()
             call iw_text("Centers and Corners",highlight=.true.)
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
                str2 = "Id" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,0)
                str2 = "Atom" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,1)
                str2 = "Show" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,2)
                str2 = "Min (Å)" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,3)
                str2 = "Max (Å)" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,4)
                do j = 1, sys(isys)%c%nspc
                   str2 = trim(sys(isys)%c%spc(j)%name) // c_null_char
                   call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0.0_c_float,4+j)
                end do
                call igTableSetupScrollFreeze(1,1)
                call igTableHeadersRow()

                do i = 1, w%rep%poly%style%ntype
                   ispc = sysc(isys)%attype_species(w%rep%poly%style%type,i)
                   iz = sys(isys)%c%spc(ispc)%z
                   if (iz <= 0) cycle
                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   if (igTableSetColumnIndex(0)) then
                      call igAlignTextToFramePadding()
                      call iw_text(string(i))
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
    use gui_main, only: g
    use utils, only: iw_text, iw_tooltip, iw_calcwidth, iw_radiobutton, iw_button,&
       iw_clamp_color3, iw_checkbox, iw_coloredit, iw_dragfloat_real8, iw_inputint
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed, ldum
    integer :: ipad
    real(c_float) :: sqw
    integer(c_int) :: nc(3)

    logical :: ch

    ! initialize
    changed = .false.

    ! periodicity
    call igAlignTextToFramePadding()
    call iw_text("Periodicity",highlight=.true.)

    ! radio buttons for the periodicity type
    changed = changed .or. iw_radiobutton("None",int=w%rep%sel%pertype,intval=0_c_int,sameline=.true.)
    call iw_tooltip("This object is represented only in the main cell and not repeated by translation",ttshown)
    changed = changed .or. iw_radiobutton("Automatic",int=w%rep%sel%pertype,intval=1_c_int,sameline=.true.)
    call iw_tooltip("Number of periodic cells controlled by the +/- options in the view menu",ttshown)
    changed = changed .or. iw_radiobutton("Manual",int=w%rep%sel%pertype,intval=2_c_int,sameline=.true.)
    call iw_tooltip("Manually set the number of periodic cells",ttshown)

    ! number of periodic cells, if manual
    if (w%rep%sel%pertype == 2_c_int) then
       ! calculate widths
       ipad = ceiling(log10(max(maxval(w%rep%sel%ncell),1) + 0.1))
       sqw = max(iw_calcwidth(1,1),igGetTextLineHeightWithSpacing())
       call igPushItemWidth(sqw)

       nc = w%rep%sel%ncell
       call igAlignTextToFramePadding()
       call iw_text("a:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##aaxis")) w%rep%sel%ncell(1) = max(w%rep%sel%ncell(1)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       ldum = iw_inputint("##aaxis",w%rep%sel%ncell(1),width=ipad,notlive=.true.)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##aaxis")) w%rep%sel%ncell(1) = w%rep%sel%ncell(1)+1

       call igSameLine(0._c_float,-1._c_float)
       call iw_text("b:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##baxis")) w%rep%sel%ncell(2) = max(w%rep%sel%ncell(2)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       ldum = iw_inputint("##baxis",w%rep%sel%ncell(2),width=ipad,notlive=.true.)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##baxis")) w%rep%sel%ncell(2) = w%rep%sel%ncell(2)+1

       call igSameLine(0._c_float,-1._c_float)
       call iw_text("c:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##caxis")) w%rep%sel%ncell(3) = max(w%rep%sel%ncell(3)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       ldum = iw_inputint("##caxis",w%rep%sel%ncell(3),width=ipad,notlive=.true.)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##caxis")) w%rep%sel%ncell(3) = w%rep%sel%ncell(3)+1
       w%rep%sel%ncell = max(w%rep%sel%ncell,1)
       if (any(nc /= w%rep%sel%ncell)) changed = .true.
       call igPopItemWidth()

       if (iw_button("Reset",sameline=.true.)) then
          w%rep%sel%ncell = 1
          changed = .true.
       end if
    end if

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
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8,&
       iw_inputtext, iw_radiobutton, iw_combo_simple
    use systems, only: sys
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch
    integer :: icoord
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
          if (associated(win(w%idparent)%sc)) then
             zf = real(win(w%idparent)%sc%gizmo_zoom_factor(),8)
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
    use utils, only: iw_text, iw_combo_simple, iw_tooltip, iw_calcheight, iw_checkbox,&
       iw_clamp_color3, iw_calcwidth, iw_button, iw_coloredit, iw_highlight_selectable,&
       iw_dragfloat_real8
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

       if (showselection) then
          icol = icol + 1
          str2 = "Show" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
       end if

       if (showdrawopts) then
          icol = icol + 1
          str2 = "Col" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)

          icol = icol + 1
          str2 = "Radius" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
       end if

       if (domol) then
          icol = icol + 1
          str2 = "Mol" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
       end if

       if (docoord) then
          icol = icol + 1
          if (r%atoms%style%type == atlisttype_ncel_ang) then
             str2 = "Coordinates (Å)" // c_null_char
          else
             str2 = "Coordinates (fractional)" // c_null_char
          end if
          flags = ImGuiTableColumnFlags_WidthStretch
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
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
                call igAlignTextToFramePadding()
                call iw_text(string(i))

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
          icol = icol + 1
          str2 = "Id" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)

          icol = icol + 1
          str2 = "nat" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)

          if (showselection) then
             icol = icol + 1
             str2 = "Show" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
          end if

          if (showdrawopts) then
             icol = icol + 1
             str2 = "Tint" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)

             icol = icol + 1
             str2 = "Scale" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)
          end if

          icol = icol + 1
          if (sys(isys)%c%ismolecule) then
             str2 = "Center of mass (Å)" // c_null_char
          else
             str2 = "Center of mass (fractional)" // c_null_char
          end if
          flags = ImGuiTableColumnFlags_WidthStretch
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,icol)

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
                   call igAlignTextToFramePadding()
                   call iw_text(string(i))

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
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_real8,&
       iw_combo_simple, iw_button, iw_calcheight
    use systems, only: sys
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch
    integer :: i, icoord, nop
    integer(c_int) :: flags
    type(ImVec2) :: sz0
    character(kind=c_char,len=:), allocatable, target :: str1, str2

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
          str2 = "show" // c_null_char
          call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,0)
          str2 = "#" // c_null_char
          call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,1)
          str2 = "Symbol" // c_null_char
          call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_None,0._c_float,2)
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

end submodule editrep
