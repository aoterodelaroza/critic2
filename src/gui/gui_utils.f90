! Copyright (c) 2019 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
module gui_utils
  use iso_c_binding
  implicit none

  private

  public :: igIsItemHovered_delayed

  ! module procedure interfaces
  interface
     module function igIsItemHovered_delayed(flags,thr,already_shown)
       integer(c_int), value :: flags
       real(c_float), intent(in) :: thr
       logical, intent(in) :: already_shown
       logical(c_bool) :: igIsItemHovered_delayed
     end function igIsItemHovered_delayed
  end interface

end module gui_utils
