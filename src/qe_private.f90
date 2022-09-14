! Module with procedures for the Quantum ESPRESSO structure reader.
! Most routines in this module were adapted from QE.
!
! Copyright (C) 2006-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Attribution of specific routines at the top of the corresponding
! procedures.

! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Quantum ESPRESSO densities and structures.
module qe_private
  implicit none

  private

  public :: qe_latgen
  public :: sup_spacegroup

  integer, parameter :: dp = selected_real_kind(14,200)

  interface
     module subroutine qe_latgen(ibrav,celldm,a1,a2,a3,errmsg)
       integer, parameter :: dp = selected_real_kind(14,200)
       integer, intent(in) :: ibrav
       real(DP), intent(inout) :: celldm(6)
       real(DP), intent(out) :: a1(3), a2(3), a3(3)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine qe_latgen
     module subroutine sup_spacegroup(tau,ityp,extfor,if_pos,space_group_number,not_eq,&
        uniqueb,rhombohedral,choice,ibrav,nattot,tautot,ityptot)
       integer, parameter :: dp = selected_real_kind(14,200)
       INTEGER, INTENT(IN) :: space_group_number, choice
       LOGICAL, INTENT (IN) :: uniqueb, rhombohedral
       INTEGER, INTENT (INOUT) ::  not_eq
       INTEGER, INTENT(OUT) :: ibrav
       REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: tau, extfor
       INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: ityp
       INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: if_pos
       integer, intent(out) :: nattot
       real(dp), allocatable, intent(inout) :: tautot(:,:)
       integer, allocatable, intent(inout) :: ityptot(:)
     end subroutine sup_spacegroup
  end interface

end module qe_private
