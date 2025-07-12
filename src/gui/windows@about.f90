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

! Routines for the about window.
submodule (windows) about
  use interfaces_cimgui
  implicit none
contains

  !> Draw the about window
  module subroutine draw_about(w)
    use config, only: istring_version, getstring
    use gui_main, only: g
    use utils, only: iw_text, iw_button, iw_calcwidth
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    class(window), intent(inout), target :: w

    real(c_float) :: wwidth, twidth
    type(ImVec2) :: szavail

    call iw_text("  ---- critic2 GUI ----",danger=.true.,centered=.true.)
    call iw_text("version " // getstring(istring_version),centered=.true.)
    call igNewLine()
    call iw_text("Critic2 is a program for visualization and analysis",centered=.true.)
    call iw_text("of structures and data in computational chemistry",centered=.true.)
    call igNewLine()
    call iw_text("Copyright (c) Alberto Otero de la Roza",centered=.true.)
    call iw_text("Distributed under GNU/GPL license, version 3",centered=.true.)
    call iw_text("Contact: aoterodelaroza@gmail.com",highlight=.true.,centered=.true.)
    call igNewLine()

    ! bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)
    wwidth = igGetWindowWidth()
    twidth = iw_calcwidth(5,1)
    call igSetCursorPosX((wwidth - twidth) * 0.5_c_float)
    if (iw_button("Close")) w%isopen = .false.

    ! exit if focused and received the close keybinding
    if ((w%focused() .and. (is_bind_event(BIND_CLOSE_FOCUSED_DIALOG) .or. is_bind_event(BIND_OK_FOCUSED_DIALOG))) .or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS))&
       w%isopen = .false.

  end subroutine draw_about

end submodule about
