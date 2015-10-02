! J.L. Casals Sainz contributed this mod to elk 2.1.25

! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoplot
! !INTERFACE:
subroutine rhoplot
! !USES:
use modmain
! !DESCRIPTION:
!   Outputs the charge density, read in from {\tt STATE.OUT}, for 1D, 2D or 3D
!   plotting.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! initialise universal variables
call init0
! read density from file
call readstate
! write the density plot to file
select case(task)
case(31)
  open(50,file='RHO1D.OUT',action='WRITE',form='FORMATTED')
  open(51,file='RHOLINES.OUT',action='WRITE',form='FORMATTED')
  call plot1d(50,51,1,lmaxvr,lmmaxvr,rhomt,rhoir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(rhoplot):")')
  write(*,'(" 1D density plot written to RHO1D.OUT")')
  write(*,'(" vertex location lines written to RHOLINES.OUT")')
case(32)
  open(50,file='RHO2D.OUT',action='WRITE',form='FORMATTED')
  call plot2d(50,1,lmaxvr,lmmaxvr,rhomt,rhoir)
  close(50)
  write(*,*)
  write(*,'("Info(rhoplot): 2D density plot written to RHO2D.OUT")')
case(33)
  open(50,file='RHO3D.OUT',action='WRITE',form='FORMATTED')
  call plot3d(50,1,lmaxvr,lmmaxvr,rhomt,rhoir)
  close(50)
  write(*,*)
  write(*,'("Info(rhoplot): 3D density plot written to RHO3D.OUT")')
case(34)
  open(50,file='RHO.OUT',action='WRITE',form='UNFORMATTED')
  write(50) lmmaxvr, nrmtmax, natmtot, ngtot
  write(50) rhomt, rhoir
  close(50)
  write(*,*)
  write(*,'("Info(elfplot): ELF mt+kpts description written to ELF.OUT")')
end select
return
end subroutine
!EOC

