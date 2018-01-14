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
  implicit none

  private 

  private :: lim_surf
  private :: lim_bundle
  private :: bisect_msurface
  private :: bundle_msurface
  public :: basinplot
  public :: bundleplot
  public :: sphereintegrals_gauleg
  public :: sphereintegrals_lebedev
  public :: sphereintegrals
  private :: integrals_gauleg
  private :: integrals_lebedev
  public :: integrals
  private :: integrals_header
  private :: minisurf_write3dmodel
  private :: minisurf_writebasin
  private :: minisurf_writedbasin
  private :: minisurf_transform
  
  interface
     module subroutine lim_surf (cpid, xin, xfin, delta, xmed, tstep, nwarn)
       integer, intent(in) :: cpid
       real*8, dimension(3), intent(in) :: xin, xfin
       real*8, intent(in) :: delta
       real*8, dimension(3), intent(out) :: xmed
       integer, intent(out) :: tstep
       integer, intent(inout) :: nwarn
     end subroutine lim_surf
     module subroutine lim_bundle(xup, xdn, xin, xfin, delta, xmed, tstep, nwarn)
       real*8, intent(in) :: xup(3), xdn(3)
       real*8, intent(inout) :: xin(3), xfin(3)
       real*8, intent(in) :: delta
       real*8, intent(out) :: xmed(3)
       integer, intent(out) :: tstep
       integer, intent(inout) :: nwarn
     end subroutine lim_bundle
     module subroutine bisect_msurface(srf,cpid,prec,verbose)
       use surface, only: minisurf
       type(minisurf), intent(inout) :: srf
       integer, intent(in) :: cpid
       real*8, intent(in) :: prec
       logical, intent(in) :: verbose
     end subroutine bisect_msurface
     module subroutine bundle_msurface(srf,prec,verbose)
       use surface, only: minisurf
       type(minisurf), intent(inout) :: srf
       real*8, intent(in) :: prec
       logical, intent(in) :: verbose
     end subroutine bundle_msurface
     module subroutine basinplot(line)
       character*(*), intent(in) :: line
     end subroutine basinplot
     module subroutine bundleplot(line)
       character*(*), intent(in) :: line
     end subroutine bundleplot
     module subroutine sphereintegrals_gauleg(x0,rad,ntheta,nphi,sprop,abserr,neval,meaneval)
       use systemmod, only: sy
       real*8, intent(in) :: x0(3), rad
       integer, intent(in) :: ntheta, nphi
       real*8, intent(out) :: sprop(sy%npropi) 
       real*8, intent(out) :: abserr
       integer, intent(out) :: neval, meaneval
     end subroutine sphereintegrals_gauleg
     module subroutine sphereintegrals_lebedev(x0,rad,nleb,sprop,abserr,neval,meaneval)
       use systemmod, only: sy
       real*8, intent(in) :: x0(3), rad
       integer, intent(in) :: nleb
       real*8, intent(out) :: sprop(sy%npropi) 
       real*8, intent(out) :: abserr
       integer, intent(out) :: neval, meaneval
     end subroutine sphereintegrals_lebedev
     module subroutine sphereintegrals(line)
       character*(*), intent(in) :: line
     end subroutine sphereintegrals
     module subroutine integrals_gauleg(atprop,n1,n2,cpid,usefiles,verbose)
       use systemmod, only: sy
       real*8, intent(out) :: atprop(sy%npropi)
       integer, intent(in) :: n1, n2, cpid
       logical, intent(in) :: usefiles
       logical, intent(in) :: verbose
     end subroutine integrals_gauleg
     module subroutine integrals_lebedev(atprop,nleb,cpid,usefiles,verbose)
       use systemmod, only: sy
       real*8, intent(out) :: atprop(sy%npropi)
       integer, intent(in) :: nleb, cpid
       logical, intent(in) :: usefiles
       logical, intent(in) :: verbose
     end subroutine integrals_lebedev
     module subroutine integrals(line)
       character*(*), intent(in) :: line
     end subroutine integrals
     module subroutine integrals_header(meth,ntheta,nphi,np,cpid,usefiles,pname)
       integer, intent(in) :: meth, ntheta, nphi, np, cpid
       logical, intent(in) :: usefiles
       character(10), intent(in) :: pname
     end subroutine integrals_header
     module subroutine minisurf_write3dmodel(s,fmt,file,expr)
       use surface, only: minisurf
       type(minisurf), intent(inout) :: s
       character*3, intent(in) :: fmt
       character*(*), intent(in) :: file
       character*(*), intent(in), optional :: expr
     end subroutine minisurf_write3dmodel
     module subroutine minisurf_writebasin(s,offfile,doprops)
       use surface, only: minisurf
       type(minisurf), intent(inout) :: s
       character*(*), intent(in) :: offfile
       logical, intent(in) :: doprops
     end subroutine minisurf_writebasin
     module subroutine minisurf_writedbasin(s,npoint,offfile)
       use surface, only: minisurf
       type(minisurf), intent(inout) :: s
       integer, intent(in) :: npoint
       character*(*), intent(in) :: offfile
     end subroutine minisurf_writedbasin
     module subroutine minisurf_transform(s,op,tvec)
       use surface, only: minisurf
       type(minisurf) :: s
       integer, intent(in) :: op
       real*8, intent(in) :: tvec(3)
     end subroutine minisurf_transform
  end interface

end module bisect

