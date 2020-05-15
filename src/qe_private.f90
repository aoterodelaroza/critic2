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

module qe_private
  implicit none

  private
  
  public :: qe_latgen
  public :: sup_spacegroup
  public :: nattot, tautot, extfortot, ityptot, if_postot
  
  integer, parameter :: dp = selected_real_kind(14,200)
  integer :: nattot
  real(dp), allocatable :: tautot(:,:), extfortot(:,:)
  integer, allocatable :: ityptot(:), if_postot(:,:)

  interface
     module subroutine qe_latgen(ibrav,celldm,a1,a2,a3,errmsg)
       integer, parameter :: dp = selected_real_kind(14,200)
       integer, intent(in) :: ibrav
       real(DP), intent(inout) :: celldm(6)
       real(DP), intent(out) :: a1(3), a2(3), a3(3)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine qe_latgen
     module subroutine sup_spacegroup(tau,ityp,extfor,if_pos,space_group_number,not_eq,&
        uniqueb,rhombohedral,choice,ibrav)
       INTEGER, INTENT(IN) :: space_group_number, choice
       LOGICAL, INTENT (IN) :: uniqueb, rhombohedral
       INTEGER, INTENT (INOUT) ::  not_eq
       INTEGER, INTENT(OUT) :: ibrav
       REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: tau, extfor
       INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: ityp
       INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: if_pos
     end subroutine sup_spacegroup
  end interface

end module qe_private
