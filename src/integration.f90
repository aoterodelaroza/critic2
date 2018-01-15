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

!>  The following applies to the lebedev subroutines of this file:
!>  This subroutine is part of a set of subroutines that generate
!>  Lebedev grids [1-6] for integration on a sphere. The original
!>  C-code [1] was kindly provided by Dr. Dmitri N. Laikov and
!>  translated into fortran by Dr. Christoph van Wuellen.
!>  This subroutine was translated from C to fortran77 by hand.
!>  Users of this code are asked to include reference [1] in their
!>  publications, and in the user- and programmers-manuals
!>  describing their codes.
!>  This code was distributed through CCL (http://www.ccl.net/).
!>  [1] V.I. Lebedev, and D.N. Laikov
!>      "A quadrature formula for the sphere of the 131st
!>       algebraic order of accuracy"
!>      Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!>
!>  [2] V.I. Lebedev
!>      "A quadrature formula for the sphere of 59th algebraic
!>       order of accuracy"
!>      Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.
!>
!>  [3] V.I. Lebedev, and A.L. Skorokhodov
!>      "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!>      Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.
!>
!>  [4] V.I. Lebedev
!>      "Spherical quadrature formulas exact to orders 25-29"
!>      Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
!>
!>  [5] V.I. Lebedev
!>      "Quadratures on a sphere"
!>      Computational Mathematics and Mathematical Physics, Vol. 16,
!>      1976, pp. 10-24.
!>
!>  [6] V.I. Lebedev
!>      "Values of the nodes and weights of ninth to seventeenth
!>       order Gauss-Markov quadrature formulae invariant under the
!>       octahedron group with inversion"
!>      Computational Mathematics and Mathematical Physics, Vol. 15,
!>      1975, pp. 44-51.

!> Quadrature schemes and integration-related tools.
module integration
  use fieldmod, only: field
  use surface, only: minisurf
  use systemmod, only: sy
  implicit none

  private

  ! integration properties
  public :: intgrid_driver
  private :: intgrid_multipoles
  private :: intgrid_deloc_wfn
  private :: intgrid_deloc_wannier
  public :: int_radialquad
  public :: gauleg_mquad
  public :: lebedev_mquad
  private :: quadpack_f
  public :: int_output
  private :: int_output_multipoles
  private :: int_output_deloc_wfn
  private :: int_output_deloc_wannier
  private :: assign_strings
  public :: int_reorder_gridout
  private :: int_gridbasins
  private :: unpackidx
  private :: packidx
  
  interface
     module subroutine intgrid_driver(line)
       character*(*), intent(in) :: line
     end subroutine intgrid_driver
     module subroutine intgrid_multipoles(fint,idprop,natt,xgatt,idg,imtype,luw,mpole)
       real*8, intent(in) :: fint(:,:,:,:)
       integer, intent(in) :: idprop(:)
       integer, intent(in) :: natt
       real*8, intent(in) :: xgatt(3,natt)
       integer, intent(in) :: idg(:,:,:)
       integer, intent(in) :: imtype
       integer, intent(in) :: luw
       real*8, allocatable, intent(inout) :: mpole(:,:,:)
     end subroutine intgrid_multipoles
     module subroutine intgrid_deloc_wfn(natt,xgatt,idg,imtype,luw,sij)
       integer, intent(in) :: natt
       real*8, intent(in) :: xgatt(3,natt)
       integer, intent(in) :: idg(:,:,:)
       integer, intent(in) :: imtype
       integer, intent(in) :: luw
       real*8, intent(inout), allocatable :: sij(:,:,:,:,:)
     end subroutine intgrid_deloc_wfn
     module subroutine intgrid_deloc_wannier(natt,xgatt,idg,imtype,luw,sij)
       integer, intent(in) :: natt
       real*8, intent(in) :: xgatt(3,natt)
       integer, intent(in) :: idg(:,:,:)
       integer, intent(in) :: imtype
       integer, intent(in) :: luw
       complex*16, allocatable, intent(inout) :: sij(:,:,:,:,:)
     end subroutine intgrid_deloc_wannier
     module subroutine int_radialquad(x,theta,phi,r0,rend,lprop,abserr,neval,iaserr,ierr)
       real*8, intent(in) :: x(3)  !< The center of the basin (cartesian coords)
       real*8, intent(in) :: theta !< Polar angle of the ray
       real*8, intent(in) :: phi   !< Azimuthal angle or the ray
       real*8, intent(in) :: r0    !< Left limit of the radial interval
       real*8, intent(in) :: rend  !< Right limit of the radial interval
       real*8, intent(out) :: lprop(sy%npropi) !< The integrated properties
       real*8, intent(out) :: abserr !< Estimated absolute error
       integer, intent(out) :: neval !< Number of evaluations of grdall
       real*8, intent(out) :: iaserr(sy%npropi) !< Estimated IAS precision error
       integer, intent(out) :: ierr
     end subroutine int_radialquad
     module subroutine gauleg_mquad(srf,ntheta,nphi,rbeta,lprop,abserr,neval,iaserr)
       type(minisurf), intent(inout) :: srf !< Surface representing the basin
       integer, intent(in) :: ntheta !< Number of polar points
       integer, intent(in) :: nphi   !< Number of azimuthal points
       real*8, intent(in) :: rbeta   !< Beta-spehre radius
       real*8, intent(out) :: lprop(sy%npropi) !< Properties vector
       real*8, intent(out) :: abserr !< Integrated radial quad. error
       integer, intent(out) :: neval !< Number of evaluations
       real*8, intent(out) :: iaserr(sy%npropi) !< Integrated IAS precision error
     end subroutine gauleg_mquad
     module subroutine lebedev_mquad(srf,npts,rbeta,lprop,abserr,neval,iaserr)
       type(minisurf), intent(inout) :: srf !< Surface representing the basin
       integer, intent(in) :: npts   !< Number of points
       real*8, intent(in) :: rbeta   !< Beta-spehre radius
       real*8, intent(out) :: lprop(sy%npropi) !< Properties vector
       real*8, intent(out) :: abserr !< Integrated radial quad. error
       integer, intent(out) :: neval !< Number of evaluations
       real*8, intent(out) :: iaserr(sy%npropi) !< Integrated IAS precision error
     end subroutine lebedev_mquad
     module function quadpack_f(x,unit,xnuc) result(res)
       real*8, intent(in) :: x
       real*8, intent(in) :: unit(3)
       real*8, intent(in) :: xnuc(3)
       real*8 :: res(sy%npropi)
     end function quadpack_f
     module subroutine int_output(pmask,reason,nattr,icp,xattr,aprop,usesym,sij,mpole)
       logical, intent(in) :: pmask(sy%npropi)
       character*(*), intent(in) :: reason(sy%npropi)
       integer, intent(in) :: nattr
       integer, intent(in) :: icp(nattr)
       real*8, intent(in) :: xattr(3,nattr)
       real*8, intent(in) :: aprop(sy%npropi,nattr)
       logical, intent(in) :: usesym
       real*8, intent(in), allocatable, optional :: sij(:,:,:,:,:)
       real*8, intent(in), allocatable, optional :: mpole(:,:,:)
     end subroutine int_output
     module subroutine int_output_multipoles(nattr,icp,mpole)
       integer, intent(in) :: nattr
       integer, intent(in) :: icp(nattr)
       real*8, intent(in), allocatable, optional :: mpole(:,:,:)
     end subroutine int_output_multipoles
     module subroutine int_output_deloc_wfn(nattr,icp,sij)
       integer, intent(in) :: nattr
       integer, intent(in) :: icp(nattr)
       real*8, intent(in), allocatable, optional :: sij(:,:,:,:,:)
     end subroutine int_output_deloc_wfn
     module subroutine int_output_deloc_wannier(natt,icp,xgatt,sij)
       integer, intent(in) :: natt
       integer, intent(in) :: icp(natt)
       real*8, intent(in) :: xgatt(3,natt)
       complex*16, intent(in), allocatable, optional :: sij(:,:,:,:,:)
     end subroutine int_output_deloc_wannier
     module subroutine assign_strings(i,icp,usesym,scp,sncp,sname,smult,sz)
       integer, intent(in) :: i
       integer, intent(in) :: icp
       logical, intent(in) :: usesym
       character(len=:), allocatable, intent(out) :: sncp, scp, sname, sz, smult
     end subroutine assign_strings
     module subroutine int_reorder_gridout(ff,nattr,xgatt,idg,atexist,ratom,luw,icp)
       type(field), intent(inout) :: ff
       integer, intent(inout) :: nattr
       real*8, intent(inout), allocatable :: xgatt(:,:)
       integer, intent(inout), allocatable :: idg(:,:,:)
       logical, intent(in) :: atexist
       real*8, intent(in) :: ratom
       integer, intent(inout) :: luw
       integer, intent(inout), allocatable :: icp(:)
     end subroutine int_reorder_gridout
     module subroutine int_gridbasins(fmt,nattr,icp,xgatt,idg,imtype,luw)
       character*3, intent(in) :: fmt
       integer, intent(in) :: nattr
       integer, intent(in), allocatable :: icp(:)
       real*8, intent(in) :: xgatt(3,nattr)
       integer, intent(in) :: idg(:,:,:)
       integer, intent(in) :: imtype
       integer, intent(in) :: luw
     end subroutine int_gridbasins
     module subroutine unpackidx(idx,io,jo,ko,bo,nmo,nbnd,nwan)
       integer, intent(in) :: idx, nmo, nbnd, nwan(3)
       integer, intent(out) :: io, jo, ko, bo
     end subroutine unpackidx
     module subroutine packidx(io,jo,ko,bo,idx,nmo,nbnd,nwan)
       integer, intent(in) :: io, jo, ko, bo, nmo, nbnd, nwan(3)
       integer, intent(out) :: idx
     end subroutine packidx
  end interface

end module integration
