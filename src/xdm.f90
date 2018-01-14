! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! This module contains code for the calculation of the dispersion
! energy (and derivatives) using the exchange-hole dipole moment
! (XDM) model.
module xdm
  implicit none

  private

  public :: xdm_driver
  private :: xdm_grid
  private :: xdm_qe
  private :: xdm_wfn
  private :: write_cube
  private :: free_volume
  private :: frevol
  private :: calc_edisp
  private :: calc_coefs
  
  interface
     module subroutine xdm_driver(line)
       character*(*), intent(inout) :: line
     end subroutine xdm_driver
     module subroutine xdm_grid(line)
       character*(*), intent(inout) :: line
     end subroutine xdm_grid
     module subroutine xdm_qe(line0)
       character*(*), intent(inout) :: line0
     end subroutine xdm_qe
     module subroutine xdm_excitation(line0)
       character*(*), intent(inout) :: line0
     end subroutine xdm_excitation
     module subroutine xdm_excitation_readpostg(file,haveit,v,vfree,mm,lvec,inow)
       use systemmod, only: sy
       character*(*), intent(in) :: file
       logical, intent(inout) :: haveit(sy%c%ncel,0:1)
       real*8, intent(inout) :: v(sy%c%ncel,0:1)
       real*8, intent(inout) :: vfree(sy%c%ncel,0:1)
       real*8, intent(inout) :: mm(3,sy%c%ncel,0:1)
       integer, intent(inout) :: lvec(3,sy%c%ncel,0:1)
       integer, intent(in) :: inow
     end subroutine xdm_excitation_readpostg
     module subroutine xdm_wfn(a1o,a2o,chf)
       real*8, intent(in) :: a1o, a2o
       real*8, intent(in) :: chf
     end subroutine xdm_wfn
     module subroutine write_cube(file,line1,line2,n,c)
       character*(*), intent(in) :: file, line1, line2
       integer, intent(in) :: n(3)
       real*8, intent(in) :: c(0:n(1)-1,0:n(2)-1,0:n(3)-1)
     end subroutine write_cube
     module function free_volume(iz) result(afree)
       integer, intent(in) :: iz
       real*8 :: afree
     end function free_volume
     module function frevol(z,chf)
       integer, intent(in) :: z
       real*8, intent(in) :: chf
       real*8 :: frevol
     endfunction frevol
     module subroutine calc_edisp(c6,c8,c10,rvdw)
       use systemmod, only: sy
       real*8, intent(in) :: c6(sy%c%ncel,sy%c%ncel), c8(sy%c%ncel,sy%c%ncel), c10(sy%c%ncel,sy%c%ncel)
       real*8, intent(in) :: rvdw(sy%c%ncel,sy%c%ncel)
     end subroutine calc_edisp
     module function calc_edisp_from_mv(a1,a2,v,vfree,mm,lvec,i0,i1)
       use systemmod, only: sy
       real*8, intent(in) :: a1, a2
       real*8, intent(in) :: v(sy%c%ncel,0:1), vfree(sy%c%ncel,0:1)
       real*8, intent(in) :: mm(3,sy%c%ncel,0:1)
       integer, intent(in) :: lvec(3,sy%c%ncel,0:1)
       integer, intent(in) :: i0, i1
       real*8 :: calc_edisp_from_mv
     end function calc_edisp_from_mv
     module subroutine calc_coefs(a1,a2,chf,v,mm,c6,c8,c10,rvdw)
       use systemmod, only: sy
       real*8, intent(in) :: a1, a2, chf
       real*8, intent(in) :: v(sy%c%ncel), mm(3,sy%c%ncel)
       real*8, intent(out) :: c6(sy%c%ncel,sy%c%ncel), c8(sy%c%ncel,sy%c%ncel)
       real*8, intent(out) :: c10(sy%c%ncel,sy%c%ncel), rvdw(sy%c%ncel,sy%c%ncel)
     end subroutine calc_coefs
     module subroutine taufromelf(ielf,irho,itau)
       integer, intent(in) :: ielf, irho, itau
     end subroutine taufromelf
  end interface

end module xdm
