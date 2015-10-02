
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potplot
! !INTERFACE:
subroutine potplot
! !USES:
use modmain
! !DESCRIPTION:
!   Outputs the exchange, correlation and Coulomb potentials, read in from
!   {\tt STATE.OUT}, for 1D, 2D or 3D plotting.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! initialise universal variables
call init0
! read the density and potentials from file
call readstate
! write the potential plots to file
select case(task)
case(41)
  open(50,file='VCL1D.OUT',action='WRITE',form='FORMATTED')
  open(51,file='VLINES.OUT',action='WRITE',form='FORMATTED')
  call plot1d(50,51,1,lmaxvr,lmmaxvr,vclmt,vclir)
  close(50)
  close(51)
  open(50,file='VXC1D.OUT',action='WRITE',form='FORMATTED')
  open(51,file='VLINES.OUT',action='WRITE',form='FORMATTED')
  call plot1d(50,51,1,lmaxvr,lmmaxvr,vxcmt,vxcir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(potplot):")')
  write(*,'(" 1D Coulomb potential plot written to VCL1D.OUT")')
  write(*,'(" 1D exchange-correlation potential plot written to VXC1D.OUT")')
  write(*,'(" vertex location lines written to VLINES.OUT")')
case(42)
  open(50,file='VCL2D.OUT',action='WRITE',form='FORMATTED')
  call plot2d(50,1,lmaxvr,lmmaxvr,vclmt,vclir)
  close(50)
  open(50,file='VXC2D.OUT',action='WRITE',form='FORMATTED')
  call plot2d(50,1,lmaxvr,lmmaxvr,vxcmt,vxcir)
  close(50)
  write(*,*)
  write(*,'("Info(potplot):")')
  write(*,'(" 2D Coulomb potential plot written to VCL2D.OUT")')
  write(*,'(" 2D exchange-correlation potential plot written to VXC2D.OUT")')
case(43)
  open(50,file='VCL3D.OUT',action='WRITE',form='FORMATTED')
  call plot3d(50,1,lmaxvr,lmmaxvr,vclmt,vclir)
  close(50)
  open(50,file='VXC3D.OUT',action='WRITE',form='FORMATTED')
  call plot3d(50,1,lmaxvr,lmmaxvr,vxcmt,vxcir)
  close(50)
  write(*,*)
  write(*,'("Info(potplot):")')
  write(*,'(" 3D Coulomb potential plot written to VCL3D.OUT")')
  write(*,'(" 3D exchange-correlation potential plot written to VXC3D.OUT")')
case(44)
  open(50,file='VCL.OUT',action='WRITE',form='UNFORMATTED')
  write(50) lmmaxvr, nrmtmax, natmtot, ngrtot
  write(50) vclmt, vclir
  close(50)
  open(50,file='VXC.OUT',action='WRITE',form='UNFORMATTED')
  write(50) lmmaxvr, nrmtmax, natmtot, ngrtot
  write(50) vxcmt, vxcir
  close(50)
  write(*,*)
  write(*,'("Info(potplot):")')
  write(*,'(" 3D Coulomb potential plot written to VCL.OUT")')
  write(*,'(" 3D exchange-correlation potential plot written to VXC.OUT")')
end select
write(*,*)
return
end subroutine
!EOC

