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

!> Integration and plotting of basins through bisection.
module bisect
  use surface, only: minisurf
  use systemmod, only: sy
  implicit none

  private 

  public :: basinplot
  public :: bundleplot
  public :: sphereintegrals
  public :: sphereintegrals_gauleg
  public :: sphereintegrals_lebedev
  public :: integrals
  
  interface
     module subroutine basinplot(line)
       character*(*), intent(in) :: line
     end subroutine basinplot
     module subroutine bundleplot(line)
       character*(*), intent(in) :: line
     end subroutine bundleplot
     module subroutine sphereintegrals(line)
       character*(*), intent(in) :: line
     end subroutine sphereintegrals
     module subroutine sphereintegrals_gauleg(x0,rad,ntheta,nphi,sprop,abserr,neval,meaneval)
       real*8, intent(in) :: x0(3), rad
       integer, intent(in) :: ntheta, nphi
       real*8, intent(out) :: sprop(sy%npropi) 
       real*8, intent(out) :: abserr
       integer, intent(out) :: neval, meaneval
     end subroutine sphereintegrals_gauleg
     module subroutine sphereintegrals_lebedev(x0,rad,nleb,sprop,abserr,neval,meaneval)
       real*8, intent(in) :: x0(3), rad
       integer, intent(in) :: nleb
       real*8, intent(out) :: sprop(sy%npropi) 
       real*8, intent(out) :: abserr
       integer, intent(out) :: neval, meaneval
     end subroutine sphereintegrals_lebedev
     module subroutine integrals(line)
       character*(*), intent(in) :: line
     end subroutine integrals
  end interface

end module bisect

