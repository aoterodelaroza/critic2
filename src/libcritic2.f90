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
  use iso_c_binding, only: c_ptr, c_int, c_double
  implicit none

  public

  ! module procedure interfaces
  interface
     module function c2_crystal_from_file(file) bind(c,name="c2_crystal_from_file")
       type(c_ptr), value, intent(in) :: file
       type(c_ptr) :: c2_crystal_from_file
     end function c2_crystal_from_file
     module function c2_crystal_from_lattice(natom,lattice,position,zat) bind(c,name="c2_crystal_from_lattice")
       integer(c_int), intent(in), value :: natom
       real(c_double), intent(in) :: lattice(3,3)
       real(c_double), intent(in) :: position(3,natom)
       integer(c_int), intent(in) :: zat(natom)
       type(c_ptr) :: c2_crystal_from_lattice
     end function c2_crystal_from_lattice
     module function c2_crystal_from_cellpar(natom,cel,ang,position,zat) bind(c,name="c2_crystal_from_cellpar")
       integer(c_int), intent(in), value :: natom
       real(c_double), intent(in) :: cel(3)
       real(c_double), intent(in) :: ang(3)
       real(c_double), intent(in) :: position(3,natom)
       integer(c_int), intent(in) :: zat(natom)
       type(c_ptr) :: c2_crystal_from_cellpar
     end function c2_crystal_from_cellpar
     module subroutine c2_describe_crystal(cr) bind(c,name="c2_describe_crystal")
       type(c_ptr), value, intent(in) :: cr
     end subroutine c2_describe_crystal
     module subroutine c2_destroy_crystal(cr) bind(c,name="c2_destroy_crystal")
       type(c_ptr), value, intent(in) :: cr
     end subroutine c2_destroy_crystal
  end interface

end module libcritic2
