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

! Display window: edit the display object of the attached view
submodule (windows) display
  use interfaces_cimgui
  implicit none

contains

  !> Draw the contents of the display window.
  module subroutine draw_display(w)
    use systems, only: sys, sysc, ok_system, sys_init
    use gui_main, only: ColorHighlightScene
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_setpos_bottomright,&
       iw_checkbox, iw_dragfloat_real8, iw_inputtext, iw_arith_help, iw_arith_help_button,&
       iw_helpermark, iw_periodicity_widget
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged, changed, ch, typechanged
    integer :: isys, iview, ihighlight, highlight_type
    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the Display of its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)
    if (goodsys) goodsys = associated(win(iview)%sc)
    if (goodsys) goodsys = (win(iview)%sc%isinit > 0)
    if (w%firstpass) w%errmsg = ""

    changed = .false.
    ihighlight = 0
    highlight_type = 0
    if (goodsys) then
       ! the system
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       associate (disp => win(iview)%sc%disp)
         ! the Show masks must describe this system before they are edited
         call disp%update(isys)

         ! filter
         call iw_text("Filter",highlight=.true.,alignframe=.true.)
         call iw_helpermark("Show the atom if the filter expression evaluates to non-zero (true) at&
            & the atomic position; structural variables are very useful for filters."//newline//&
            iw_arith_help//newline//"Click the Help button for more info.")
         call iw_arith_help_button("##helpfilter",ttshown)
         if (iw_inputtext("##filtertext",bufsize=1023,texta=disp%filter,notlive=.true.)) then
            disp%errfilter = ""
            changed = .true.
         end if
         if (len_trim(disp%filter) == 0) disp%errfilter = ""
         call iw_tooltip("Apply this filter to the atoms in the system. Atoms are represented if non-zero.",&
            ttshown)
         if (iw_button("Clear",sameline=.true.)) then
            disp%filter = ""
            disp%errfilter = ""
            changed = .true.
         end if
         call iw_tooltip("Clear the filter",ttshown)
         if (len_trim(disp%errfilter) > 0) &
            call iw_text("Error: " // trim(disp%errfilter),danger=.true.,wrap=.true.)

         if (.not.sys(isys)%c%ismolecule) then
            ! periodicity (evaluate first: the widgets must be drawn even if
            ! changed is already true, and .or. is allowed to short-circuit)
            ch = iw_periodicity_widget(disp%ncell,ttshown)
            changed = changed .or. ch

            ! whole molecules and cell borders
            changed = changed .or. iw_checkbox("Show connected molecules",disp%onemotif)
            call iw_tooltip("Translate atoms to display whole molecules",ttshown)
            changed = changed .or. iw_checkbox("Show atoms at cell edges",disp%border,sameline=.true.)
            call iw_tooltip("Display atoms near the unit cell edges",ttshown)

            ! origin translation and cell origin shift
            changed = changed .or. iw_dragfloat_real8("Translate Origin (fractional)##originatom",&
               x3=disp%origin,speed=0.001d0,decimal=5)
            call iw_tooltip("Translation vector for the contents of the unit cell.",ttshown)
            changed = changed .or. iw_dragfloat_real8("Cell Origin Shift (fractional)##origincell",&
               x3=disp%tshift,speed=0.001d0,decimal=5)
            call iw_tooltip("Displace the origin of the cell being represented.",ttshown)
         end if

         ! the atom and molecule Show tables (a change of grouping remakes
         ! the masks, and the table is drawn next frame)
         ch = atom_table_widget(isys,disp%atype,typechanged,ihighlight,highlight_type,shown=disp%ashown)
         if (typechanged) call disp%reset_shown(isys)
         changed = changed .or. ch
         ch = mol_table_widget(isys,ihighlight,highlight_type,shown=disp%mshown)
         changed = changed .or. ch

         ! reset to the defaults for this system
         if (iw_button("Reset",danger=.true.)) then
            call disp%init(isys)
            changed = .true.
         end if
         call iw_tooltip("Reset the selection to the defaults for this system",ttshown)
       end associate

       ! the row under the mouse, highlighted in the scene
       if (ihighlight > 0) then
          call sysc(isys)%highlight_atoms(.true.,(/ihighlight/),highlight_type,&
             reshape(ColorHighlightScene,(/4,1/)))
       end if

       ! rebuild the draw lists if necessary
       if (changed) win(iview)%sc%forcebuildlists = .true.
    end if

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(5,1)

    ! final buttons: close
    if (iw_button("Close")) doquit = .true.
    call iw_tooltip("Close this window",ttshown)

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_display

end submodule display
