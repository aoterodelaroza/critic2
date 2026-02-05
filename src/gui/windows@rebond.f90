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

! Routines for the rebond window.
submodule (windows) rebond
  use interfaces_cimgui
  implicit none
contains

  !> Update tasks for the rebond window, before the window is
  !> created.
  module subroutine update_rebond(w)
    use tools_io, only: string
    class(window), intent(inout), target :: w

    if (w%firstpass.or.w%tied_to_tree) then
       w%name = "Recalculate Bonds###rebond"  // string(w%id) // c_null_char
    else
       w%name = "Recalculate Bonds [detached]###rebond"  // string(w%id) // c_null_char
    end if

  end subroutine update_rebond

  !> Draw the contents of the rebond window.
  module subroutine draw_rebond(w)
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use systems, only: nsys, sysc, sys, sys_init, ok_system, lastchange_rebond
    use gui_main, only: g
    use utils, only: iw_text, iw_tooltip, iw_calcwidth, iw_button, iw_calcheight,&
       iw_dragfloat_real8
    use global, only: bondfactor_def
    use tools_io, only: string, nameguess
    use param, only: atmcov0, maxzat0, bohrtoa, newline
    class(window), intent(inout), target :: w

    logical :: doquit, ch
    integer :: i, iz, isys, natused
    type(ImVec2) :: szavail, szero, sz0
    integer(c_int) :: flags
    real(c_float) :: combowidth
    logical(c_bool) :: is_selected
    character(len=:,kind=c_char), allocatable, target :: str1, str2
    integer, allocatable :: iat(:)
    logical :: atused(maxzat0)

    integer, parameter :: ic_name = 0
    integer, parameter :: ic_z = 1
    integer, parameter :: ic_radius = 2

    ! window state
    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    ! first pass
    if (w%firstpass) then
       w%tied_to_tree = (w%isys == win(iwin_tree)%tree_selected)
    end if

    ! if tied to tree, update the isys
    if (w%tied_to_tree .and. (w%isys /= win(iwin_tree)%tree_selected)) &
       w%isys = win(iwin_tree)%tree_selected

    ! check if the system still exists
    if (.not.ok_system(w%isys,sys_init)) then
       ! this dialog does not make sense anymore, close it and exit
       call w%end()
       return
    end if
    isys = w%isys

    ! system combo
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
                w%isys = i
                isys = w%isys
                w%tied_to_tree = w%tied_to_tree .and. (w%isys == win(iwin_tree)%tree_selected)
             end if
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Recalculate the bonds in this system",ttshown)

    ! determine which atoms are to be used
    atused = .false.
    do i = 1, sys(isys)%c%nspc
       atused(sys(isys)%c%spc(i)%z) = .true.
    end do
    allocate(iat(count(atused)))
    natused = 0
    do i = 1, maxzat0
       if (atused(i)) then
          natused = natused + 1
          iat(natused) = i
       end if
    end do

    ! the radii table
    call iw_text("Atomic Radii",highlight=.true.)
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tableatomrcov" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(5,natused)+1,0,.false.)
    if (igBeginTable(c_loc(str1),3,flags,sz0,0._c_float)) then
       ! header setup
       str2 = "Atom" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_name)

       str2 = "Z" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_z)

       str2 = "Radius (Å)" // c_null_char
       flags = ImGuiTableColumnFlags_WidthStretch
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_radius)
       call igTableSetupScrollFreeze(0, 1) ! top row always visible

       ! draw the header
       call igTableHeadersRow()
       call igTableSetColumnWidthAutoAll(igGetCurrentTable())

       ! draw the rows
       do i = 1, natused
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
          iz = iat(i)

          ! name
          if (igTableSetColumnIndex(ic_name)) then
             call igAlignTextToFramePadding()
             call iw_text(string(nameguess(iz,.true.)))
          end if

          ! Z
          if (igTableSetColumnIndex(ic_z)) then
             call igAlignTextToFramePadding()
             call iw_text(string(iz))
          end if

          ! radius
          if (igTableSetColumnIndex(ic_radius)) then
             ch = iw_dragfloat_real8("##tableradius" // string(i),x1=sysc(isys)%atmcov(iz),speed=0.01d0,&
                min=0d0,max=2.65d0,scale=bohrtoa,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
          end if
       end do
       call igEndTable()
    end if

    ! bond factor
    call iw_text("Bond factor",highlight=.true.)
    call igSameLine(0._c_float,-1._c_float)
    ch = iw_dragfloat_real8("##bondfactor",x1=sysc(isys)%bondfactor,speed=0.001d0,min=1d0,max=2d0,&
       decimal=4,flags=ImGuiSliderFlags_AlwaysClamp)
    write (*,*) "xx bond = ", sysc(isys)%bondfactor
    call iw_tooltip("Bond factor parameter for connectivity calculation",ttshown)

    ! explanation message
    call iw_text("Atoms i and j are bonded if:"//newline//&
       "  d < (ri+rj)*(bond factor)"//newline//&
       "where d = atom distance; ri, rj = atom radii")

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(34,4,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! reset to the default bonds
    if (iw_button("Reset",danger=.true.)) then
       sysc(isys)%atmcov = atmcov0
       sysc(isys)%bondfactor = bondfactor_def
    end if
    call iw_tooltip("Reset the default bonding parameters and recalculate the system bonds",ttshown)

    ! apply the changes to all systems
    if (iw_button("Apply to All Systems",sameline=.true.)) then
       do i = 1, nsys
          if (sysc(i)%status >= sys_init) then
             sysc(i)%atmcov = sysc(isys)%atmcov
             sysc(i)%bondfactor = sysc(isys)%bondfactor
             call sys(i)%c%find_asterisms(sys(i)%c%nstar,sysc(i)%atmcov,sysc(i)%bondfactor)
             call sys(i)%c%fill_molecular_fragments()
             call sys(i)%c%calculate_molecular_equivalence()
             call sys(i)%c%calculate_periodicity()
             call sysc(i)%post_event(lastchange_rebond)
          end if
       end do
    end if

    ! apply the changes
    if (iw_button("Apply",sameline=.true.)) then
       ! find the atomic connectivity and the molecular fragments
       call sys(isys)%c%find_asterisms(sys(isys)%c%nstar,sysc(isys)%atmcov,sysc(isys)%bondfactor)
       call sys(isys)%c%fill_molecular_fragments()
       call sys(isys)%c%calculate_molecular_equivalence()
       call sys(isys)%c%calculate_periodicity()
       call sysc(isys)%post_event(lastchange_rebond)
    end if
    call iw_tooltip("Recalculate the system bonds with the selected parameters",ttshown)

    ! close button
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.
    doquit = doquit .or. iw_button("Close",sameline=.true.)

    ! quit the window
    if (doquit) &
       call w%end()

  end subroutine draw_rebond

end submodule rebond
