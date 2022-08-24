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

! Some utilities for building the GUI (e.g. wrappers around ImGui routines).
submodule (gui_utils) proc
  use iso_c_binding
  implicit none

contains

  ! Returns true if the last item has been hovered for at least thr
  ! seconds. If already_shown (the tooltip has already been displayed),
  ! do not use the delay.
  module function igIsItemHovered_delayed(flags,thr,already_shown)
    use gui_main, only: g
    use gui_interfaces_cimgui, only: igIsItemHovered
    integer(c_int), value :: flags
    real(c_float), intent(in) :: thr
       logical, intent(in) :: already_shown
    logical(c_bool) :: igIsItemHovered_delayed

    igIsItemHovered_delayed = igIsItemHovered(flags) .and. (already_shown .or. g%HoveredIdTimer >= thr)

  end function igIsItemHovered_delayed

end submodule proc
