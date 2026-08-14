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
    use utils, only: iw_text, iw_button, iw_setpos_bottomright, iw_close_event
    class(window), intent(inout), target :: w

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

    ! bottom-align and center the rest of the contents
    call iw_setpos_bottomright(5,1,centered=.true.)
    if (iw_button("Close")) w%isopen = .false.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) w%isopen = .false.

  end subroutine draw_about

end submodule about
