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
  use iso_c_binding, only: c_ptr, c_int, c_double, c_bool
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
     module function c2_crystal_get_natom(cr) bind(c,name="c2_crystal_get_natom")
       type(c_ptr), value, intent(in) :: cr
       integer(c_int) :: c2_crystal_get_natom
     end function c2_crystal_get_natom
     module subroutine c2_crystal_get_structure(cr,natom,cel,ang,lattice,position,zat) &
        bind(c,name="c2_crystal_get_structure")
       type(c_ptr), value, intent(in) :: cr
       integer(c_int) :: natom
       type(c_ptr), value :: cel
       type(c_ptr), value :: ang
       type(c_ptr), value :: lattice
       type(c_ptr), value :: position
       type(c_ptr), value :: zat
     end subroutine c2_crystal_get_structure
     module subroutine c2_write_crystal(cr,file) bind(c,name="c2_write_crystal")
       type(c_ptr), value, intent(in) :: cr
       type(c_ptr), value, intent(in) :: file
     end subroutine c2_write_crystal
     module subroutine c2_destroy_crystal(cr) bind(c,name="c2_destroy_crystal")
       type(c_ptr), value, intent(in) :: cr
     end subroutine c2_destroy_crystal
     module function c2_peaks_from_crystal(cr,th2ini,th2end,lambda,fpol) bind(c,name="c2_peaks_from_crystal")
       type(c_ptr), value, intent(in) :: cr
       real(c_double), value :: th2ini, th2end, lambda, fpol
       type(c_ptr) :: c2_peaks_from_crystal
     end function c2_peaks_from_crystal
     module function c2_peaks_from_file(file) bind(c,name="c2_peaks_from_file")
       type(c_ptr), value, intent(in) :: file
       type(c_ptr) :: c2_peaks_from_file
     end function c2_peaks_from_file
     module subroutine c2_write_peaks(pk,file) bind(c,name="c2_write_peaks")
       type(c_ptr), value, intent(in) :: pk
       type(c_ptr), value, intent(in) :: file
     end subroutine c2_write_peaks
     module subroutine c2_destroy_peaks(pk) bind(c,name="c2_destroy_peaks")
       type(c_ptr), value, intent(in) :: pk
     end subroutine c2_destroy_peaks
     module function c2_compare_gpwdf(c1,p2,alpha,lambda,fpol) bind(c,name="c2_compare_gpwdf")
       type(c_ptr), value, intent(in) :: c1
       type(c_ptr), value, intent(in) :: p2
       real(c_double), value :: alpha, lambda, fpol
       real(c_double) :: c2_compare_gpwdf
     end function c2_compare_gpwdf
     module function c2_compare_vcgpwdf(c1,p2,crout,global,verbose,alpha,lambda,fpol,&
        maxfeval,besteps,max_elong,max_ang) bind(c,name="c2_compare_vcgpwdf")
       type(c_ptr), value, intent(in) :: c1
       type(c_ptr), value, intent(in) :: p2
       type(c_ptr) :: crout
       logical(c_bool), value :: global, verbose
       real(c_double), value :: alpha, lambda, fpol, besteps, max_elong, max_ang
       integer(c_int), value :: maxfeval
       real(c_double) :: c2_compare_vcgpwdf
     end function c2_compare_vcgpwdf
     module function c2_compare_vcgpwdf_global_safe(c1,p2,crout,verbose,lambda,fpol) &
        bind(c,name="c2_compare_vcgpwdf_global_safe")
       type(c_ptr), value, intent(in) :: c1
       type(c_ptr), value, intent(in) :: p2
       type(c_ptr) :: crout
       logical(c_bool), value :: verbose
       real(c_double), value :: lambda, fpol
       real(c_double) :: c2_compare_vcgpwdf_global_safe
     end function c2_compare_vcgpwdf_global_safe
     module function c2_compare_vcgpwdf_global_quick(c1,p2,crout,verbose,lambda,fpol) &
        bind(c,name="c2_compare_vcgpwdf_global_quick")
       type(c_ptr), value, intent(in) :: c1
       type(c_ptr), value, intent(in) :: p2
       type(c_ptr) :: crout
       logical(c_bool), value :: verbose
       real(c_double), value :: lambda, fpol
       real(c_double) :: c2_compare_vcgpwdf_global_quick
     end function c2_compare_vcgpwdf_global_quick
  end interface

end module libcritic2
