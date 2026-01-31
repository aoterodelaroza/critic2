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
    use representations, only: representation, reptype_atoms, reptype_unitcell
    use windows, only: nwin, win, wintype_view
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG,&
       BIND_CLOSE_ALL_DIALOGS
    use systems, only: sysc, sys_init, ok_system
    use gui_main, only: g
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_calcheight, iw_checkbox
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: isys, ll, itype
    logical :: doquit, ok
    logical(c_bool) :: changed
    character(kind=c_char,len=:), allocatable, target :: str1
    character(kind=c_char,len=1024), target :: txtinp
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

       ! name and type block
       call iw_text("Type and Name",highlight=.true.)

       ! the representation type
       itype = w%rep%type - 1
       call iw_combo_simple("##reptype","Atoms" // c_null_char // "Unit cell" // c_null_char,itype)
       if (w%rep%type /= itype + 1) changed = .true.
       w%rep%type = itype + 1
       call iw_tooltip("Type of object",ttshown)

       ! name text input
       str1 = "##nametextinput"
       txtinp = trim(adjustl(w%rep%name)) // c_null_char
       call igSameLine(0._c_float,-1._c_float)
       call igPushItemWidth(iw_calcwidth(30,1))
       if (igInputText(c_loc(str1),c_loc(txtinp),1023_c_size_t,ImGuiInputTextFlags_None,c_null_funptr,c_null_ptr)) then
          ll = index(txtinp,c_null_char)
          w%rep%name = txtinp(1:ll-1)
       end if
       call igPopItemWidth()
       call iw_tooltip("Name of this object",ttshown)

       ! shown checkbox
       changed = changed .or. iw_checkbox("Show",w%rep%shown,sameline=.true.)
       call iw_tooltip("Toggle show/hide this object",ttshown)

       ! type-dependent items
       if (w%rep%type == reptype_atoms) then
          changed = changed .or. w%draw_editrep_atoms(ttshown)
       elseif (w%rep%type == reptype_unitcell) then
          changed = changed .or. w%draw_editrep_unitcell(ttshown)
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
          call w%rep%set_defaults(win(w%idparent)%sc%style,0)
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
       atlisttype_nneq
    use gui_main, only: g, ColorHighlightScene, ColorElement
    use tools_io, only: string
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_radiobutton, iw_calcheight, iw_clamp_color3, iw_checkbox, iw_coloredit,&
       iw_highlight_selectable, iw_dragfloat_real8, iw_dragfloat_realc
    use param, only: atmcov, atmvdw, newline, jmlcol, jmlcol2, bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical(c_bool) :: changed

    integer :: ispc, isys, iz, ll, ipad
    character(kind=c_char,len=1024), target :: txtinp
    character(kind=c_char,len=33), target :: txtinp2
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3, suffix
    real*8 :: x0(3)
    logical(c_bool) :: ch, ldum
    integer(c_int) :: nc(3), lst, flags, nspcpair
    real(c_float) :: sqw, raux
    integer :: i, j, k, intable, nrow, is, ncol, ihighlight, highlight_type
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
    changed = changed .or. iw_checkbox("Atoms##atomsglobaldisplay",w%rep%atoms_display,highlight=.true.)
    call iw_tooltip("Display atoms in the scene",ttshown)

    changed = changed .or. iw_checkbox("Bonds##bondsglobaldisplay",w%rep%bonds_display,sameline=.true.,highlight=.true.)
    call iw_tooltip("Display bonds in the scene",ttshown)

    changed = changed .or. iw_checkbox("Labels##labelsglobaldisplay",w%rep%labels_display,sameline=.true.,highlight=.true.)
    call iw_tooltip("Display atomic labels in the scene",ttshown)

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
          call iw_text("(?)",sameline=.true.)
          call iw_tooltip("Show the atom if the filter expression evaluates to non-zero (true) at the atomic position. &
             &Structural variables are very useful for filters. Examples:"//newline//&
             "- '@x < 3' = all atoms with x lower than 3"//newline//&
             "- 'log($0) > 1' = log of the promolecular density higher than 1"//newline//&
             "- 'abs(@x) < 2 && abs(@y) < 2 && abs(@z) < 2' = atoms in the (-2,2) box"//newline//&
             "Click on the Help button for more info.")
          if (iw_button("Help##helpfilter",sameline=.true.)) then
             str3 = "https://aoterodelaroza.github.io/critic2/manual/arithmetics" // c_null_char
             call openLink(c_loc(str3))
          end if
          call iw_tooltip("Open the manual page about arithmetic expressions."&
             &"The 'basic usage' and 'structural variables' sections are relevant.",ttshown)

          ! filter text input
          str1 = "##filtertext" // c_null_char
          txtinp = trim(adjustl(w%rep%filter)) // c_null_char
          if (igInputText(c_loc(str1),c_loc(txtinp),1023_c_size_t,ImGuiInputTextFlags_EnterReturnsTrue,&
             c_null_funptr,c_null_ptr)) then
             ll = index(txtinp,c_null_char)
             w%rep%filter = txtinp(1:ll-1)

             ! test the filter
             if (sys(isys)%c%ncel > 0) then
                x0 = sys(isys)%c%atcel(1)%r
             else
                x0 = 0d0
             end if
             changed = .true.
             w%rep%errfilter = ""
          end if
          if (len_trim(w%rep%filter) == 0) w%rep%errfilter = ""
          call iw_tooltip("Apply this filter to the atoms in the system. Atoms are represented if non-zero.",&
             ttshown)
          if (iw_button("Clear",sameline=.true.)) then
             w%rep%filter = ""
             w%rep%errfilter = ""
             changed = .true.
          end if
          call iw_tooltip("Clear the filter",ttshown)
          if (len_trim(w%rep%errfilter) > 0) &
             call iw_text("Error: " // trim(w%rep%errfilter),danger=.true.)

          ! periodicity
          if (.not.sys(isys)%c%ismolecule) then
             call igAlignTextToFramePadding()
             call iw_text("Periodicity",highlight=.true.)

             ! radio buttons for the periodicity type
             changed = changed .or. iw_radiobutton("None",int=w%rep%pertype,intval=0_c_int,sameline=.true.)
             call iw_tooltip("This object is represented only in the main cell and not repeated by translation",ttshown)
             changed = changed .or. iw_radiobutton("Automatic",int=w%rep%pertype,intval=1_c_int,sameline=.true.)
             call iw_tooltip("Number of periodic cells controlled by the +/- options &
                &in the 'Scene' button of the view window",ttshown)
             changed = changed .or. iw_radiobutton("Manual",int=w%rep%pertype,intval=2_c_int,sameline=.true.)
             call iw_tooltip("Manually set the number of periodic cells",ttshown)

             ! number of periodic cells, if manual
             if (w%rep%pertype == 2_c_int) then
                ! calculate widths
                ipad = ceiling(log10(max(maxval(w%rep%ncell),1) + 0.1))
                sqw = max(iw_calcwidth(1,1),igGetTextLineHeightWithSpacing())
                call igPushItemWidth(sqw)

                nc = w%rep%ncell
                call igAlignTextToFramePadding()
                call iw_text("a:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##aaxis")) w%rep%ncell(1) = max(w%rep%ncell(1)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                str2 = "##aaxis" // c_null_char
                call igPushItemWidth(iw_calcwidth(ipad,1))
                ldum = igInputInt(c_loc(str2),w%rep%ncell(1),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
                call igPopItemWidth()
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##aaxis")) w%rep%ncell(1) = w%rep%ncell(1)+1

                call igSameLine(0._c_float,-1._c_float)
                call iw_text("b:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##baxis")) w%rep%ncell(2) = max(w%rep%ncell(2)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                str2 = "##baxis" // c_null_char
                call igPushItemWidth(iw_calcwidth(ipad,1))
                ldum = igInputInt(c_loc(str2),w%rep%ncell(2),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
                call igPopItemWidth()
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##baxis")) w%rep%ncell(2) = w%rep%ncell(2)+1

                call igSameLine(0._c_float,-1._c_float)
                call iw_text("c:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##caxis")) w%rep%ncell(3) = max(w%rep%ncell(3)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                str2 = "##caxis" // c_null_char
                call igPushItemWidth(iw_calcwidth(ipad,1))
                ldum = igInputInt(c_loc(str2),w%rep%ncell(3),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
                call igPopItemWidth()
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##caxis")) w%rep%ncell(3) = w%rep%ncell(3)+1
                w%rep%ncell = max(w%rep%ncell,1)
                if (any(nc /= w%rep%ncell)) changed = .true.
                call igPopItemWidth()

                if (iw_button("Reset",sameline=.true.)) then
                   w%rep%ncell = 1
                   changed = .true.
                end if
             end if

             ! checkbox for molecular motif
             changed = changed .or. iw_checkbox("Show connected molecules",w%rep%onemotif)
             call iw_tooltip("Translate atoms to display whole molecules",ttshown)

             ! checkbox for border
             changed = changed .or. iw_checkbox("Show atoms at cell edges",w%rep%border,sameline=.true.)
             call iw_tooltip("Display atoms near the unit cell edges",ttshown)
          end if

          ! origin of the atoms
          call igPushItemWidth(iw_calcwidth(21,3))
          if (.not.sys(isys)%c%ismolecule) then
             ! origin translation
             changed = changed .or. iw_dragfloat_real8("Translate Origin (fractional)##originatom",x3=w%rep%origin,&
                speed=0.001d0,sformat="%.5f")
             call iw_tooltip("Translation vector for the contents of the unit cell.",ttshown)

             ! origin shift
             changed = changed .or. iw_dragfloat_real8("Cell Origin Shift (fractional)##origincell",x3=w%rep%tshift,&
                speed=0.001d0,sformat="%.5f")
             call iw_tooltip("Displace the origin of the cell being represented.",ttshown)
          end if
          call igPopItemWidth()

          ! draw the atom selection widget
          changed = changed .or. atom_selection_widget(isys,w%rep,&
             .true.,.false.,ihighlight,highlight_type)

          call igEndTabItem()
       end if ! begin tab item (selection)

       !!!!! Atoms tab !!!!!
       if (w%rep%atoms_display) then
          str1 = "Atoms##editrepatoms_atomstab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
             ! global options for atoms
             call igAlignTextToFramePadding()
             call iw_text("Global Options",highlight=.true.)
             if (iw_button("Reset##resetglobalatoms",sameline=.true.,danger=.true.)) then
                call w%rep%set_defaults(win(w%idparent)%sc%style,1)
                changed = .true.
             end if
             call iw_tooltip("Reset to the default settings for the atom representation")
             call iw_combo_simple("Radii ##atomradiicombo","Covalent"//c_null_char//"Van der Waals"//c_null_char//&
                "Constant"//c_null_char,w%rep%atom_radii_type,changed=ch)
             call iw_tooltip("Set atomic radii to the tabulated values of this type",ttshown)

             call igSameLine(0._c_float,-1._c_float)
             call igPushItemWidth(iw_calcwidth(5,1))
             if (w%rep%atom_radii_type == 2) then
                ! constant size
                raux = w%rep%atom_radii_value * real(bohrtoa,c_float)
                ch = ch .or. iw_dragfloat_realc("Value##atomradii",x1=raux,speed=0.01_c_float,&
                   min=0._c_float,max=5._c_float,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
                call iw_tooltip("Atomic radii (Å)",ttshown)

                if (ch) then
                   w%rep%atom_radii_value = raux / real(bohrtoa,c_float)
                   w%rep%atom_style%rad(1:w%rep%atom_style%ntype) =w%rep%atom_radii_value
                   changed = .true.
                end if
             else
                ! variable size
                ch = ch .or. iw_dragfloat_realc("Scale##atomradiiscale",x1=w%rep%atom_radii_scale,speed=0.01_c_float,&
                   min=0._c_float,max=5._c_float,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
                call iw_tooltip("Scale factor for the tabulated atomic radii",ttshown)

                if (ch) then
                   do i = 1, w%rep%atom_style%ntype
                      ispc = sysc(isys)%attype_species(w%rep%atom_style%type,i)
                      iz = sys(isys)%c%spc(ispc)%z
                      if (w%rep%atom_radii_type == 0) then
                         w%rep%atom_style%rad(i) = real(atmcov(iz),c_float)
                      else
                         w%rep%atom_style%rad(i) = real(atmvdw(iz),c_float)
                      end if
                      w%rep%atom_style%rad(i) = w%rep%atom_style%rad(i) * w%rep%atom_radii_scale
                   end do
                   changed = .true.
                end if
             end if
             call igPopItemWidth()

             ! style buttons: set color
             call iw_combo_simple("Colors ##atomcolorselect","Current defaults" // c_null_char //&
                "jmol (light)" // c_null_char // "jmol2 (dark)" // c_null_char,&
                w%rep%atom_color_type,changed=ch)
             call iw_tooltip("Set the color of all atoms to the tabulated values",ttshown)
             if (ch) then
                do i = 1, w%rep%atom_style%ntype
                   ispc = sysc(isys)%attype_species(w%rep%atom_style%type,i)
                   iz = sys(isys)%c%spc(ispc)%z
                   if (w%rep%atom_color_type == 0) then
                      w%rep%atom_style%rgb(:,i) = ColorElement(:,iz)
                   elseif (w%rep%atom_color_type == 1) then
                      w%rep%atom_style%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
                   else
                      w%rep%atom_style%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
                   end if
                end do
                changed = .true.
             end if

             ! border size
             call igPushItemWidth(iw_calcwidth(5,1))
             changed = changed .or. iw_dragfloat_realc("Border Size (Å)",x1=w%rep%atom_border_size,&
                speed=0.002_c_float,min=0._c_float,max=1._c_float,scale=real(bohrtoa,c_float),&
                sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the thickness of the atom borders",ttshown)
             call igPopItemWidth()

             ! color
             changed = changed .or. iw_coloredit("Border Color",rgb=w%rep%atom_border_rgb,sameline=.true.)
             call iw_tooltip("Color of the border for the atoms",ttshown)

             ! draw the atom selection widget
             changed = changed .or. atom_selection_widget(isys,w%rep,&
                .false.,.true.,ihighlight,highlight_type)

             call igEndTabItem()
          end if ! begin tab item (atoms)
       end if

       !!!!! Bonds tab !!!!!
       if (w%rep%bonds_display) then
          str1 = "Bonds##editrepatoms_bondstab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
             !! bonds display !!

             !! global options !!
             call igAlignTextToFramePadding()
             call iw_text("Global Options",highlight=.true.)
             if (iw_button("Reset##resetglobal",sameline=.true.,danger=.true.)) then
                call w%rep%set_defaults(win(w%idparent)%sc%style,2)
                changed = .true.
             end if
             call iw_tooltip("Reset to the covalent bonding for this system and the default settings")

             ! rest of the options (record changes)
             ch = .false.
             call igAlignTextToFramePadding()
             call iw_text("Style")
             call iw_combo_simple("##tablebondstyleglobalselect",&
                "Single color"//c_null_char//"Two colors"//c_null_char,w%rep%bond_color_style,sameline=.true.,changed=ch)
             call iw_tooltip("Use a single color for the bond, or two colors from the bonded atoms",ttshown)

             call iw_text(" Radius (Å)",sameline=.true.)
             call igPushItemWidth(iw_calcwidth(5,1))
             call igSameLine(0._c_float,-1._c_float)
             ch = ch .or. iw_dragfloat_realc("##radiusbondtableglobal",x1=w%rep%bond_rad,speed=0.005_c_float,&
                min=0._c_float,max=2._c_float,scale=real(bohrtoa,c_float),sformat="%.3f",&
                flags=ImGuiSliderFlags_AlwaysClamp)
             call igPopItemWidth()
             call iw_tooltip("Radius of the bonds",ttshown)

             ! border size
             call igPushItemWidth(iw_calcwidth(5,1))
             ch = ch .or. iw_dragfloat_realc("Border Size (Å)",x1=w%rep%bond_border_size,speed=0.002_c_float,&
                min=0._c_float,max=1._c_float,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the thickness of the bond borders",ttshown)
             call igPopItemWidth()

             ! color
             ch = ch .or. iw_coloredit("Border Color",rgb=w%rep%bond_border_rgb,sameline=.true.)
             call iw_tooltip("Color of the border for the bonds",ttshown)

             ! color
             call igAlignTextToFramePadding()
             call iw_text("Color")
             ch = ch .or. iw_coloredit("##colorbondtableglobal",rgb=w%rep%bond_rgb,sameline=.true.)
             call iw_tooltip("Color of the bonds",ttshown)

             ! order
             call iw_text(" Order",sameline=.true.)
             call iw_combo_simple("##tablebondorderselectglobal",&
                "Dashed"//c_null_char//"Single"//c_null_char//"Double"//c_null_char//"Triple"//c_null_char,&
                w%rep%bond_order,sameline=.true.,changed=ldum)
             ch = ch .or. ldum
             call iw_tooltip("Bond order (dashed, single, double, etc.)",ttshown)

             ! both atoms
             call iw_text(" Both Atoms",sameline=.true.)
             ch = ch .or. iw_checkbox("##bothatomstableglobal",w%rep%bond_bothends,sameline=.true.)
             call iw_tooltip("Represent a bond if both end-atoms are in the scene (checked) or if only &
                &one end-atom is in the scene (unchecked)",ttshown)

             !! distance block !!
             call igAlignTextToFramePadding()
             call iw_text("Distances",highlight=.true.)
             call iw_text(" (",highlight=.true.,sameline_nospace=.true.)
             call iw_combo_simple("##tablebondglobaldistcombo","Factor"//c_null_char//"Range"//c_null_char,&
                w%rep%bond_distancetype,sameline_nospace=.true.)
             call iw_tooltip("Draw bonds whose lengths are a factor of the sum of atomic&
                & radii (Factor) or give bond distance range (Range)",ttshown)
             call iw_text(")",highlight=.true.,sameline_nospace=.true.)
             if (iw_button("Apply##applyglobal",sameline=.true.,danger=.true.)) then
                call w%rep%bond_style%generate_neighstars(w%rep)
                w%rep%bond_style%use_sys_nstar = .false.
                changed = .true.
             end if
             call iw_tooltip("Recalculate and draw bonds using the selected distance criteria",ttshown)

             if (w%rep%bond_distancetype == 0) then
                ! factor
                call igAlignTextToFramePadding()
                call iw_text("Between")
                call igSameLine(0._c_float,-1._c_float)
                call igPushItemWidth(iw_calcwidth(5,1))
                ldum = iw_dragfloat_real8("##bondtableglobalbfmin",x1=w%rep%bond_bfmin,speed=0.01d0,min=0.0d0,&
                   max=9.999d0,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
                call igPopItemWidth()
                call iw_tooltip("Bonds with length below this factor times the radii are not shown",ttshown)

                call iw_text("times",sameline=.true.)
                call iw_combo_simple("##bondtableglobalradtypemin","cov."//c_null_char//"vdw"//c_null_char,&
                   w%rep%bond_radtype(1),sameline=.true.)
                call iw_tooltip("Choose the atomic radii (covalent or van der Waals)",ttshown)
                call iw_text("radii",sameline=.true.)

                call igAlignTextToFramePadding()
                call iw_text("... and")
                call igSameLine(0._c_float,-1._c_float)
                call igPushItemWidth(iw_calcwidth(5,1))
                str2 = "##bondtableglobalbfmax" // c_null_char
                ldum = iw_dragfloat_real8("##bondtableglobalbfmax",x1=w%rep%bond_bfmax,speed=0.01d0,min=0.0d0,&
                   max=9.999d0,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
                call igPopItemWidth()
                call iw_tooltip("Bonds with length above this factor times the radii are not shown",ttshown)

                call iw_text("times",sameline=.true.)
                call iw_combo_simple("##bondtableglobalradtypemax","cov."//c_null_char//"vdw"//c_null_char,&
                   w%rep%bond_radtype(2),sameline=.true.)
                call iw_tooltip("Choose the atomic radii (covalent or van der Waals)",ttshown)
                call iw_text("radii",sameline=.true.)
             else
                ! range
                call igAlignTextToFramePadding()
                call iw_text("Between")
                call igSameLine(0._c_float,-1._c_float)
                call igPushItemWidth(iw_calcwidth(5,1))
                ldum = iw_dragfloat_real8("##bondtableglobaldmin",x1=w%rep%bond_dmin,speed=0.01d0,min=0.0d0,&
                   max=9.999d0,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
                call igPopItemWidth()
                call iw_text("Å",sameline=.true.)
                call iw_tooltip("Bonds with length below this factor times the radii are not shown",ttshown)

                call igAlignTextToFramePadding()
                call iw_text("... and")
                call igSameLine(0._c_float,-1._c_float)
                call igPushItemWidth(iw_calcwidth(5,1))
                ldum = iw_dragfloat_real8("##bondtableglobaldmax",x1=w%rep%bond_dmax,speed=0.01d0,min=0.0d0,&
                   max=9.999d0,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
                call igPopItemWidth()
                call iw_tooltip("Bonds with length above this factor times the radii are not shown",ttshown)
                call iw_text("Å",sameline=.true.)
             end if
             call igAlignTextToFramePadding()
             call iw_text("Intra/Inter-molecular")
             call iw_combo_simple("##tablebondimolselectglobal",&
                "any"//c_null_char//"intra"//c_null_char//"inter"//c_null_char,&
                w%rep%bond_imol,sameline=.true.)
             call iw_tooltip("Draw any bonds (any), only intramolecular (intra), or only intermolecular (inter)",&
                ttshown)

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
                         if (iw_checkbox("##bondtableshown" // suffix,w%rep%bond_style%shown(i,j))) then
                            ch = .true.
                            w%rep%bond_style%shown(j,i) = w%rep%bond_style%shown(i,j)
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
                w%rep%bond_style%shown = .true.
                ch = .true.
             end if
             call iw_tooltip("Show all bonds in the system",ttshown)
             if (iw_button("Hide All##hideallbonds",sameline=.true.)) then
                w%rep%bond_style%shown = .false.
                ch = .true.
             end if
             call iw_tooltip("Hide all bonds in the system",ttshown)
             if (iw_button("Toggle Show/Hide##toggleallatoms",sameline=.true.)) then
                do i = 1, sys(isys)%c%nspc
                   do j = i, sys(isys)%c%nspc
                      w%rep%bond_style%shown(j,i) = .not.w%rep%bond_style%shown(j,i)
                      w%rep%bond_style%shown(i,j) = w%rep%bond_style%shown(j,i)
                   end do
                end do
                ch = .true.
             end if
             call iw_tooltip("Toggle the show/hide status for all bonds",ttshown)

             ! immediately update if non-distances have changed
             if (ch) changed = .true.

             call igEndTabItem()
          end if ! begin tab item (bonds)
       end if

       !!!!! Labels tab !!!!!
       if (w%rep%labels_display) then
          str1 = "Labels##editrepatoms_labelstab" // c_null_char
          flags = ImGuiTabItemFlags_None
          if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
             !! labels display !!

             ! label styles
             !! global options !!
             call igAlignTextToFramePadding()
             call iw_text("Global Options",highlight=.true.)
             if (iw_button("Reset##resetglobal",sameline=.true.,danger=.true.)) then
                w%rep%label_type = 0
                call w%rep%set_defaults(win(w%idparent)%sc%style,3)
                changed = .true.
             end if
             call iw_tooltip("Reset to the labels to the default settings")

             if (sys(isys)%c%ismolecule) then
                lst = lsttrans(w%rep%label_type)
                call iw_combo_simple("Text##labelcontentselect","Atomic symbol"//c_null_char//&
                   "Atom name"// c_null_char//"Atom ID"// c_null_char//&
                   "Species ID"// c_null_char// "Atomic number"// c_null_char// "Molecule ID"// c_null_char,&
                   lst,changed=ch)
                w%rep%label_type = lsttransi(lst)
             else
                call iw_combo_simple("Text##labelcontentselect","Atomic symbol"//c_null_char//&
                   "Atom name"//c_null_char//"Cell atom ID"//c_null_char//&
                   "Cell atom ID + lattice vector"//c_null_char//"Symmetry-unique atom ID"//c_null_char//&
                   "Species ID"//c_null_char//"Atomic number"//c_null_char//"Molecule ID"//c_null_char//&
                   "Wyckoff position"//c_null_char,&
                   w%rep%label_type,changed=ch)
             end if
             if (ch) call w%rep%label_style%reset(w%rep)
             call iw_tooltip("Text to display in the atom labels",ttshown)
             changed = changed .or. ch

             ! scale, constant size, color
             call igPushItemWidth(iw_calcwidth(4,1))
             changed = changed .or. iw_dragfloat_realc("Scale##labelscale",x1=w%rep%label_scale,speed=0.01_c_float,&
                min=0._c_float,max=10._c_float,sformat="%.2f",flags=ImGuiSliderFlags_AlwaysClamp)
             call igPopItemWidth()
             call iw_tooltip("Scale factor for the atom labels",ttshown)

             changed = changed .or. iw_checkbox("Constant size##labelconstsize",&
                w%rep%label_const_size,sameline=.true.)
             call iw_tooltip("Labels have constant size (on) or labels scale with the&
                & size of the associated atom (off)",ttshown)

             changed = changed .or. iw_coloredit("Color##labelcolor",rgb=w%rep%label_rgb,sameline=.true.)
             call iw_tooltip("Color of the atom labels",ttshown)

             ! offset
             call igPushItemWidth(iw_calcwidth(21,3))
             changed = changed .or. iw_dragfloat_realc("Offset (Å)",x3=w%rep%label_offset,&
                speed=0.001_c_float,min=99.999_c_float,max=99.999_c_float,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Offset the position of the labels relative to the atom center",ttshown)
             call igPopItemWidth()

             ! table for label selection
             call iw_text("Label Selection",highlight=.true.)

             ! number of entries in the table
             select case(w%rep%label_type)
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
                         changed = changed .or. iw_checkbox("##labeltableshown" // suffix,w%rep%label_style%shown(i))
                         call iw_tooltip("Toggle display of labels for these atoms/molecules",ttshown)
                      end if

                      ! text
                      ncol = ncol + 1
                      if (igTableSetColumnIndex(ncol)) then
                         str1 = "##labeltabletext" // suffix // c_null_char
                         txtinp2 = trim(w%rep%label_style%str(i)) // c_null_char
                         call igPushItemWidth(iw_calcwidth(15,1))
                         if (igInputText(c_loc(str1),c_loc(txtinp2),32_c_size_t,ImGuiInputTextFlags_None,&
                            c_null_funptr,c_null_ptr)) then
                            ll = index(txtinp2,c_null_char)
                            w%rep%label_style%str(i) = txtinp2(1:ll-1)
                            changed = .true.
                         end if
                         call igPopItemWidth()
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
                w%rep%label_style%shown = .true.
                changed = .true.
             end if
             call iw_tooltip("Show all labels",ttshown)
             if (iw_button("Hide All##hidealllabels",sameline=.true.)) then
                w%rep%label_style%shown = .false.
                changed = .true.
             end if
             call iw_tooltip("Hide all labels",ttshown)
             if (iw_button("Toggle Show/Hide##togglealllabels",sameline=.true.)) then
                w%rep%label_style%shown = .not.w%rep%label_style%shown
                changed = .true.
             end if
             call iw_tooltip("Toggle the show/hide status for all bonds",ttshown)

             call igEndTabItem()
          end if ! begin tab item (labels)
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
       iw_clamp_color3, iw_checkbox, iw_coloredit, iw_dragfloat_realc, iw_dragfloat_real8
    use param, only: bohrtoa
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical(c_bool) :: changed, ldum
    integer :: ipad
    real(c_float) :: sqw
    integer(c_int) :: nc(3)

    character(kind=c_char,len=:), allocatable, target :: str2
    logical :: ch

    ! initialize
    changed = .false.

    ! periodicity
    call igAlignTextToFramePadding()
    call iw_text("Periodicity",highlight=.true.)

    ! radio buttons for the periodicity type
    changed = changed .or. iw_radiobutton("None",int=w%rep%pertype,intval=0_c_int,sameline=.true.)
    call iw_tooltip("This object is represented only in the main cell and not repeated by translation",ttshown)
    changed = changed .or. iw_radiobutton("Automatic",int=w%rep%pertype,intval=1_c_int,sameline=.true.)
    call iw_tooltip("Number of periodic cells controlled by the +/- options in the view menu",ttshown)
    changed = changed .or. iw_radiobutton("Manual",int=w%rep%pertype,intval=2_c_int,sameline=.true.)
    call iw_tooltip("Manually set the number of periodic cells",ttshown)

    ! number of periodic cells, if manual
    if (w%rep%pertype == 2_c_int) then
       ! calculate widths
       ipad = ceiling(log10(max(maxval(w%rep%ncell),1) + 0.1))
       sqw = max(iw_calcwidth(1,1),igGetTextLineHeightWithSpacing())
       call igPushItemWidth(sqw)

       nc = w%rep%ncell
       call igAlignTextToFramePadding()
       call iw_text("a:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##aaxis")) w%rep%ncell(1) = max(w%rep%ncell(1)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       str2 = "##aaxis" // c_null_char
       call igPushItemWidth(iw_calcwidth(ipad,1))
       ldum = igInputInt(c_loc(str2),w%rep%ncell(1),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
       call igPopItemWidth()
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##aaxis")) w%rep%ncell(1) = w%rep%ncell(1)+1

       call igSameLine(0._c_float,-1._c_float)
       call iw_text("b:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##baxis")) w%rep%ncell(2) = max(w%rep%ncell(2)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       str2 = "##baxis" // c_null_char
       call igPushItemWidth(iw_calcwidth(ipad,1))
       ldum = igInputInt(c_loc(str2),w%rep%ncell(2),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
       call igPopItemWidth()
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##baxis")) w%rep%ncell(2) = w%rep%ncell(2)+1

       call igSameLine(0._c_float,-1._c_float)
       call iw_text("c:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##caxis")) w%rep%ncell(3) = max(w%rep%ncell(3)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       str2 = "##caxis" // c_null_char
       call igPushItemWidth(iw_calcwidth(ipad,1))
       ldum = igInputInt(c_loc(str2),w%rep%ncell(3),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
       call igPopItemWidth()
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##caxis")) w%rep%ncell(3) = w%rep%ncell(3)+1
       w%rep%ncell = max(w%rep%ncell,1)
       if (any(nc /= w%rep%ncell)) changed = .true.
       call igPopItemWidth()

       if (iw_button("Reset",sameline=.true.)) then
          w%rep%ncell = 1
          changed = .true.
       end if
    end if

    !! styles
    call iw_text("Style",highlight=.true.)
    changed = changed .or. iw_checkbox("Color crystallographic axes",w%rep%uc_coloraxes)
    call iw_tooltip("Represent crystallographic axes with colors (a=red,b=green,c=blue)",ttshown)
    changed = changed .or. iw_checkbox("Hide axes in vacuum directions",w%rep%uc_vaccutsticks)
    call iw_tooltip("In systems with vacuum direction(s), do not show the unit cell in the vacuum region",ttshown)

    call igPushItemWidth(iw_calcwidth(5,1))
    changed = changed .or. iw_dragfloat_realc("Radius (Å)##outer",x1=w%rep%uc_radius,speed=0.005_c_float,&
       min=0._c_float,max=5._c_float,scale=real(bohrtoa,c_float),sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
    call igPopItemWidth()
    call iw_tooltip("Radii of the unit cell edges",ttshown)

    ch = iw_coloredit("Color",rgb=w%rep%uc_rgb,sameline=.true.)
    call iw_tooltip("Color of the unit cell edges",ttshown)
    if (ch) then
       w%rep%uc_rgb = min(w%rep%uc_rgb,1._c_float)
       w%rep%uc_rgb = max(w%rep%uc_rgb,0._c_float)
       changed = .true.
    end if

    !! inner divisions
    call iw_text("Inner Divisions",highlight=.true.)
    changed = changed .or. iw_checkbox("Display inner divisions",w%rep%uc_inner)
    call iw_tooltip("Represent the inner divisions inside a supercell",ttshown)
    if (w%rep%uc_inner) then
       call igPushItemWidth(iw_calcwidth(5,1))
       changed = changed .or. iw_dragfloat_realc("Radius (Å)##inner",x1=w%rep%uc_radiusinner,speed=0.005_c_float,&
          min=0._c_float,max=5._c_float,scale=real(bohrtoa,c_float),sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
       call igPopItemWidth()
       call iw_tooltip("Radii of the inner unit cell edges",ttshown)

       changed = changed .or. iw_checkbox("Use dashed lines",w%rep%uc_innerstipple)
       call iw_tooltip("Use dashed lines for the inner cell divisions",ttshown)

       if (w%rep%uc_innerstipple) then
          call igPushItemWidth(iw_calcwidth(5,1))
          changed = changed .or. iw_dragfloat_realc("Dash length (Å)",x1=w%rep%uc_innersteplen,speed=0.1_c_float,&
             min=0.001_c_float,max=100._c_float,scale=real(bohrtoa,c_float),sformat="%.1f",&
             flags=ImGuiSliderFlags_AlwaysClamp)
          call igPopItemWidth()
          call iw_tooltip("Length of the dashed lines for the inner cell divisions (in Å)",ttshown)
       end if
    end if

    ! origin of the unit cell
    call iw_text("Origin Shift",highlight=.true.)
    call igPushItemWidth(iw_calcwidth(21,3))
    changed = changed .or. iw_dragfloat_real8("##originucx",x3=w%rep%origin,speed=0.001d0,sformat="%.5f")
    call iw_tooltip("Coordinates for the origin shift of the unit cell",ttshown)
    call igPopItemWidth()

  end function draw_editrep_unitcell

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
       iw_dragfloat_realc
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
    logical(c_bool) :: ch
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2, str3, suffix
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
    ch = sysc(isys)%attype_combo_simple("Atom types##atomtypeselection",r%atom_style%type,&
       atlisttype_allowed,units=.false.)
    call iw_tooltip("Group atoms by these categories",ttshown)
    if (ch) then
       call r%atom_style%reset(r)
       changed = .true.
    end if

    ! whether to do the molecule column and the coordinates
    domol = (r%atom_style%type == atlisttype_ncel_ang)
    docoord = (r%atom_style%type == atlisttype_nneq .or. r%atom_style%type == atlisttype_ncel_ang .or.&
       r%atom_style%type == atlisttype_ncel_frac)
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
    sz0%y = iw_calcheight(min(5,r%atom_style%ntype)+1,0,.false.)
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
          if (r%atom_style%type == atlisttype_ncel_ang) then
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
       call ImGuiListClipper_Begin(clipper,r%atom_style%ntype,-1._c_float)

       ! draw the rows
       do while(ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
             suffix = "_" // string(i)
             icol = -1

             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             ispc = sysc(isys)%attype_species(r%atom_style%type,i)
             iz = sys(isys)%c%spc(ispc)%z

             ! id
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) then
                call igAlignTextToFramePadding()
                call iw_text(string(i))

                ! the highlight selectable
                if (iw_highlight_selectable("##selectablemoltable" // suffix)) then
                   ihighlight = i
                   highlight_type = r%atom_style%type
                end if
             end if

             ! name
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) call iw_text(string(sys(isys)%c%spc(ispc)%name))

             ! Z
             icol = icol + 1
             if (igTableSetColumnIndex(icol)) call iw_text(string(iz))

             ! shown
             if (showselection) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   changed = changed .or. iw_checkbox("##tableshown" // suffix ,r%atom_style%shown(i))
                   call iw_tooltip("Toggle display of the atom/bond/label associated to this atom",ttshown)
                end if
             end if

             ! color
             if (showdrawopts) then
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   ch = iw_coloredit("##tablecolor" // suffix,rgb=r%atom_style%rgb(:,i))
                   call iw_tooltip("Atom color",ttshown)
                   if (ch) then
                      r%atom_style%rgb(:,i) = min(r%atom_style%rgb(:,i),1._c_float)
                      r%atom_style%rgb(:,i) = max(r%atom_style%rgb(:,i),0._c_float)
                      changed = .true.
                   end if
                end if

                ! radius
                icol = icol + 1
                if (igTableSetColumnIndex(icol)) then
                   call igPushItemWidth(iw_calcwidth(5,1))
                   ch = iw_dragfloat_realc("##tableradius" // string(i),x1=r%atom_style%rad(i),speed=0.01_c_float,&
                      min=0._c_float,max=5._c_float,scale=real(bohrtoa,c_float),sformat=str3,flags=ImGuiSliderFlags_AlwaysClamp)
                   call iw_tooltip("Radius of the sphere representing the atom",ttshown)
                   if (ch) then
                      r%atom_style%rad(i) = max(r%atom_style%rad(i),0._c_float)
                      changed = .true.
                   end if
                   call igPopItemWidth()
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
                   if (r%atom_style%type > 0) then
                      x0 = sysc(isys)%attype_coordinates(r%atom_style%type,i)
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
          r%atom_style%shown = .true.
          changed = .true.
       end if
       call iw_tooltip("Show all atoms/bonds/labels in the system",ttshown)
       if (iw_button("Hide All##hideallatoms",sameline=.true.)) then
          r%atom_style%shown = .false.
          changed = .true.
       end if
       call iw_tooltip("Hide all atoms/bonds/labels in the system",ttshown)
       if (iw_button("Toggle Show/Hide##toggleallatoms",sameline=.true.)) then
          do i = 1, r%atom_style%ntype
             r%atom_style%shown(i) = .not.r%atom_style%shown(i)
          end do
          changed = .true.
       end if
       call iw_tooltip("Toggle the show/hide status for all atoms/bonds/labels",ttshown)
    end if

    ! molecule selection
    ! initialized and more than one molecule
    if (r%mol_style%isinit .and. r%mol_style%ntype > 1) then
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
       sz0%y = iw_calcheight(min(5,r%mol_style%ntype)+1,0,.false.)
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
          call ImGuiListClipper_Begin(clipper,r%mol_style%ntype,-1._c_float)

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
                      changed = changed .or. iw_checkbox("##tablemolshown" // string(i) ,r%mol_style%shown(i))
                      call iw_tooltip("Toggle display of all atoms in this molecule",ttshown)
                   end if
                end if

                ! color
                if (showdrawopts) then
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      ch = iw_coloredit("##tablemolcolor" // string(i),rgb=r%mol_style%tint_rgb(:,i))
                      call iw_tooltip("Molecule color tint",ttshown)
                      if (ch) then
                         r%mol_style%tint_rgb(:,i) = min(r%mol_style%tint_rgb(:,i),1._c_float)
                         r%mol_style%tint_rgb(:,i) = max(r%mol_style%tint_rgb(:,i),0._c_float)
                         changed = .true.
                      end if
                   end if

                   ! radius
                   icol = icol + 1
                   if (igTableSetColumnIndex(icol)) then
                      call igPushItemWidth(iw_calcwidth(5,1))
                      ch = iw_dragfloat_realc("##tablemolradius" // string(i),x1=r%mol_style%scale_rad(i),&
                         speed=0.005_c_float,min=0._c_float,max=5._c_float,&
                         sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)
                      call iw_tooltip("Scale factor for the atomic radii in this molecule",ttshown)
                      if (ch) then
                         r%mol_style%scale_rad(i) = max(r%mol_style%scale_rad(i),0._c_float)
                         changed = .true.
                      end if
                      call igPopItemWidth()
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
             r%mol_style%shown = .true.
             changed = .true.
          end if
          call iw_tooltip("Show all molecules in the system",ttshown)
          if (iw_button("Hide All##hideallmolecules",sameline=.true.)) then
             r%mol_style%shown = .false.
             changed = .true.
          end if
          call iw_tooltip("Hide all molecules in the system",ttshown)
          if (iw_button("Toggle Show/Hide##toggleallmolecules",sameline=.true.)) then
             do i = 1, r%mol_style%ntype
                r%mol_style%shown(i) = .not.r%mol_style%shown(i)
             end do
             changed = .true.
          end if
          call iw_tooltip("Toggle the show/hide status for all molecules",ttshown)
       end if
    end if

  end function atom_selection_widget

end submodule editrep
