! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.

! Public interface for the critic2 C/C++ library.
module libcritic2
  use iso_c_binding, only: c_ptr
  implicit none

  public

  ! module procedure interfaces
  interface
     module function create_structure_from_file(file) bind(c,name="create_structure_from_file")
       type(c_ptr), value, intent(in) :: file
       type(c_ptr) :: create_structure_from_file
     end function create_structure_from_file
     module subroutine destroy_structure(cr) bind(c,name="destroy_structure")
       type(c_ptr), value, intent(in) :: cr
     end subroutine destroy_structure
  end interface

end module libcritic2
