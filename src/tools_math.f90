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

! Math functions and subroutines
module tools_math
  implicit none

  private

  public :: crosscorr_triangle
  public :: crys2car_from_cellpar
  public :: car2crys_from_cellpar
  public :: factorial
  public :: genrlm_real
  public :: genylm
  public :: tosphere
  public :: ylmderiv
  public :: radial_derivs
  public :: ep
  public :: gcd
  public :: eig
  public :: eigns
  public :: rsindex
  public :: cross
  public :: mixed
  public :: mnorm2
  public :: detsym
  public :: det
  public :: matinv
  public :: plane_scale_extend
  public :: assign_ziso
  public :: comb
  public :: nchoosek
  public :: rmsd_walker
  public :: gauleg
  public :: good_lebedev
  public :: select_lebedev

  ! types of contour level choosing strategies
  integer, parameter, public :: niso_manual = 0
  integer, parameter, public :: niso_lin = 1
  integer, parameter, public :: niso_log = 2
  integer, parameter, public :: niso_atan = 3
  integer, parameter, public :: niso_bader = 4
  
  interface
     !xx! proc submodule
     module function crosscorr_triangle(h,f,g,l) result(dfg)
       use tools_io, only: ferror, faterr
       real*8, intent(in) :: h, l
       real*8, intent(in) :: f(:), g(:)
       real*8 :: dfg
     end function crosscorr_triangle
     module function crys2car_from_cellpar(aal,bbl) result(mat)
       real*8, intent(in) :: aal(3),bbl(3)
       real*8 :: mat(3,3)
     end function crys2car_from_cellpar
     module function car2crys_from_cellpar(aal,bbl) result(mat)
       real*8, intent(in) :: aal(3),bbl(3)
       real*8 :: mat(3,3)
     end function car2crys_from_cellpar
     module function factorial(n) result(f)
       real*8 :: f
       integer, intent(in) :: n
     end function factorial
     module subroutine genrlm_real(lmax,r,tp,rrlm)
       integer, intent(in) :: lmax
       real*8, intent(in) :: r
       real*8, intent(in) :: tp(2)
       real*8, intent(out) :: rrlm((lmax+1)*(lmax+1))
     end subroutine genrlm_real
     module subroutine genylm(lmax,tp,ylm)
       integer, intent(in) :: lmax
       real*8, intent(in) :: tp(2)
       complex*16, intent(out) :: ylm(*)
     end subroutine genylm
     module subroutine tosphere(v,r,tp)
       real*8, intent(in) :: v(3)
       real*8, intent(out) :: r
       real*8, intent(out) :: tp(2)
     end subroutine tosphere
     module subroutine ylmderiv(yl,r,l,m,c,cp,cpp,grad,hess)
       complex*16, dimension(:), intent(in) :: yl
       real*8, intent(in) :: r
       integer, intent(in) :: l, m
       real*8, intent(in) :: c, cp, cpp
       complex*16, intent(out) :: grad(3), hess(6)
     end subroutine ylmderiv
     module subroutine radial_derivs (rlm,rho,rho1,rho2,r0,a,b)
       real*8, dimension(:), intent(in) :: rlm
       real*8, intent(out) :: rho, rho1, rho2
       real*8, intent(in) :: r0
       real*8, intent(in) :: a, b
     end subroutine radial_derivs
     module function ep(x,i)
       integer, intent(in) :: i
       real*8, intent(in) ::  x
       real*8 :: ep
     end function ep
     module function gcd(n,num)
       integer :: n(num)
       integer :: num
       integer :: gcd 
     end function gcd
     module subroutine eig(mat,eval)
       real*8, intent(inout) :: mat(:,:)
       real*8, intent(out), optional :: eval(:)
     end subroutine eig
     module subroutine eigns(mat,eval,evali)
       real*8, intent(inout) :: mat(:,:)
       real*8, intent(out), optional :: eval(:)
       real*8, intent(out), optional :: evali(:)
     end subroutine eigns
     module subroutine rsindex(mat,ehess,r,s,eps)
       real*8, intent(inout)  :: mat(3,3)
       real*8, intent(out) :: ehess(3)
       integer, intent(out) :: r, s
       real*8, intent(in) :: eps
     end subroutine rsindex
     module function cross(v1,v2) result (vx)
       real*8, intent(in) :: v1(3)
       real*8, intent(in) :: v2(3)
       real*8 :: vx(3)
     end function cross
     module function mixed(v1,v2,v3)
       real*8, intent(in) :: v1(3)
       real*8, intent(in) :: v2(3)
       real*8, intent(in) :: v3(3)
       real*8 :: mixed
     end function mixed
     module function norm(v)
       real*8, intent(in) :: v(3)
       real*8 :: norm
     end function norm
     module function mnorm2(a)
       real*8, intent(in) :: a(3,3)
       real*8 :: mnorm2
     end function mnorm2
     module function detsym(m)
       real*8, intent(in) :: m(3,3)
       real*8 :: detsym
     end function detsym
     module function det(m)
       real*8, intent(in) :: m(3,3)
       real*8 :: det
     end function det
     module function matinv(m) result(mo)
       real*8, intent(in) :: m(3,3)
       real*8 :: mo(3,3)
     end function matinv
     module subroutine plane_scale_extend(x0,x1,x2,sxi,syi,zx0i,zx1i,zy0i,zy1i)
       real*8, intent(inout) :: x0(3), x1(3), x2(3)
       real*8, intent(in), optional :: sxi, syi, zx0i, zx1i, zy0i, zy1i
     end subroutine plane_scale_extend
     module subroutine assign_ziso(niso_type,niso,ziso,lin0,lin1,fmax,fmin)
       integer, intent(in) :: niso_type
       integer, intent(inout) :: niso
       real*8, allocatable, intent(inout) :: ziso(:)
       real*8, intent(in) :: lin0, lin1, fmax, fmin
     end subroutine assign_ziso
     module subroutine comb(n, p, l, c)
       integer :: n, p, l, c(p)
     end subroutine comb
     module function nchoosek(M, N)
       integer :: m, n
       integer :: nchoosek
     end function nchoosek
     module function rmsd_walker(x1o,x2o) result(rmsd)
       real*8, intent(in) :: x1o(:,:), x2o(:,:)
       real*8 :: rmsd
     end function rmsd_walker
     module subroutine gauleg (x1,x2,x,w,n)
       real*8, intent(in) :: x1
       real*8, intent(in) :: x2
       real*8, dimension(n), intent(out) :: x
       real*8, dimension(n), intent(out) :: w
       integer, intent(in) :: n
     end subroutine gauleg
     !xx! lebedev submodule
     module subroutine good_lebedev(npts)
       integer, intent(inout) :: npts
     end subroutine good_lebedev
     module subroutine select_lebedev(npts,xleb,yleb,zleb,wleb)
       integer, intent(in) :: npts
       real*8, intent(out) :: xleb(:), yleb(:), zleb(:), wleb(:)
     end subroutine select_lebedev
  end interface

end module tools_math
