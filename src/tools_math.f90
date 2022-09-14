! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

  !xx! proc submodule !xx!
  public :: crosscorr_triangle
  public :: m_x2c_from_cellpar
  public :: m_c2x_from_cellpar
  public :: factorial
  public :: genrlm_real
  public :: genylm
  public :: tosphere
  public :: ylmderiv
  public :: radial_derivs
  public :: ep
  public :: gcd
  public :: lcm
  public :: rational_approx
  public :: lattice_direction
  public :: eigsym
  public :: eig
  public :: rsindex
  public :: cross
  public :: mixed
  public :: mnorm2
  public :: det3sym
  public :: det3
  public :: matinv
  public :: matinvsym
  public :: plane_scale_extend
  public :: assign_ziso
  public :: comb
  public :: nchoosek
  public :: rmsd_walker
  public :: gauleg
  public :: bhole
  public :: xlnorm
  public :: fdamp_tt
  public :: fdamp_bj
  public :: munkres
  public :: umeyama_graph_matching
  public :: ullmann_graph_matching
  public :: invert_permutation
  !xx! lebedev submodule !xx!
  public :: good_lebedev
  public :: select_lebedev
  !xx! sobol submodule !xx!
  public :: sobol_sequence ! NOT THREAD-SAFE

  ! types of contour level choosing strategies
  integer, parameter, public :: niso_manual = 0
  integer, parameter, public :: niso_lin = 1
  integer, parameter, public :: niso_log = 2
  integer, parameter, public :: niso_atan = 3
  integer, parameter, public :: niso_bader = 4

  ! overloaded functions
  interface gcd
     module procedure gcd2
     module procedure gcdn
     module procedure gcd2_i8
     module procedure gcdn_i8
  end interface gcd
  interface lcm
     module procedure lcm2
     module procedure lcmn
     module procedure lcm2_i8
     module procedure lcmn_i8
  end interface lcm

  interface
     !xx! proc submodule
     module function crosscorr_triangle(h,f,g,l) result(dfg)
       use tools_io, only: ferror, faterr
       real*8, intent(in) :: h, l
       real*8, intent(in) :: f(:), g(:)
       real*8 :: dfg
     end function crosscorr_triangle
     module function m_x2c_from_cellpar(aal,bbl) result(mat)
       real*8, intent(in) :: aal(3),bbl(3)
       real*8 :: mat(3,3)
     end function m_x2c_from_cellpar
     module function m_c2x_from_cellpar(aal,bbl) result(mat)
       real*8, intent(in) :: aal(3),bbl(3)
       real*8 :: mat(3,3)
     end function m_c2x_from_cellpar
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
     module subroutine ylmderiv(yl,r,l,m,c,cp,cpp,nder,grad,hess)
       complex*16, dimension(:), intent(in) :: yl
       real*8, intent(in) :: r
       integer, intent(in) :: l, m
       real*8, intent(in) :: c, cp, cpp
       integer, intent(in) :: nder
       complex*16, intent(out), optional :: grad(3), hess(6)
     end subroutine ylmderiv
     module subroutine radial_derivs(rlm,a,b,r0,nder,rho,rho1,rho2)
       real*8, dimension(:), intent(in) :: rlm
       integer, intent(in) :: nder
       real*8, intent(in) :: r0
       real*8, intent(in) :: a, b
       real*8, intent(out) :: rho
       real*8, intent(out), optional :: rho1, rho2
     end subroutine radial_derivs
     module function ep(x,i)
       integer, intent(in) :: i
       real*8, intent(in) ::  x
       real*8 :: ep
     end function ep
     module function gcdn(n,num)
       integer, intent(in) :: n(num)
       integer, intent(in) :: num
       integer :: gcdn
     end function gcdn
     module function gcd2(m,n)
       integer, intent(in) :: m, n
       integer :: gcd2
     end function gcd2
     module function lcmn(n,num)
       integer, intent(in) :: n(num)
       integer, intent(in) :: num
       integer :: lcmn
     end function lcmn
     module function lcm2(m,n)
       integer, intent(in) :: m, n
       integer :: lcm2
     end function lcm2
     module function gcdn_i8(n,num)
       integer*8, intent(in) :: n(num)
       integer, intent(in) :: num
       integer*8 :: gcdn_i8
     end function gcdn_i8
     module function gcd2_i8(m,n)
       integer*8, intent(in) :: m, n
       integer*8 :: gcd2_i8
     end function gcd2_i8
     module function lcmn_i8(n,num)
       integer*8, intent(in) :: n(num)
       integer, intent(in) :: num
       integer*8 :: lcmn_i8
     end function lcmn_i8
     module function lcm2_i8(m,n)
       integer*8, intent(in) :: m, n
       integer*8 :: lcm2_i8
     end function lcm2_i8
     module subroutine rational_approx(x0,q,r,eps)
       real*8, intent(in) :: x0
       integer*8, intent(out):: q, r
       real*8, intent(in) :: eps
     end subroutine rational_approx
     module function lattice_direction(x0,allowr) result (yy)
       real*8, intent(in) :: x0(3)
       logical, intent(in) :: allowr
       integer :: yy(3)
     end function lattice_direction
     module subroutine eigsym(mat,n0,eval)
       integer, intent(in) :: n0
       real*8, intent(inout) :: mat(n0,n0)
       real*8, intent(out) :: eval(n0)
     end subroutine eigsym
     module subroutine eig(mat,n0,eval,evali)
       integer, intent(in) :: n0
       real*8, intent(inout) :: mat(n0,n0)
       real*8, intent(out) :: eval(n0)
       real*8, intent(out) :: evali(n0)
     end subroutine eig
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
     module function det3sym(m)
       real*8, intent(in) :: m(3,3)
       real*8 :: det3sym
     end function det3sym
     module function det3(m)
       real*8, intent(in) :: m(3,3)
       real*8 :: det3
     end function det3
     module subroutine matinv(m,n0,ier)
       integer, intent(in) :: n0
       real*8, intent(inout) :: m(n0,n0)
       integer, intent(out), optional :: ier
     end subroutine matinv
     module subroutine matinvsym(m,n0,ier)
       integer, intent(in) :: n0
       real*8, intent(inout) :: m(n0,n0)
       integer, intent(out), optional :: ier
     end subroutine matinvsym
     module subroutine plane_scale_extend(x0,x1,x2,sxi,syi,zx0i,zx1i,zy0i,zy1i)
       real*8, intent(inout) :: x0(3), x1(3), x2(3)
       real*8, intent(in), optional :: sxi, syi, zx0i, zx1i, zy0i, zy1i
     end subroutine plane_scale_extend
     module subroutine assign_ziso(niso_type,niso,ziso,fmin,fmax)
       integer, intent(in) :: niso_type
       integer, intent(inout) :: niso
       real*8, allocatable, intent(inout) :: ziso(:)
       real*8, intent(in) :: fmin, fmax
     end subroutine assign_ziso
     module subroutine comb(n, p, l, c)
       integer :: n, p, l, c(p)
     end subroutine comb
     module function nchoosek(M, N)
       integer :: m, n
       integer :: nchoosek
     end function nchoosek
     module function rmsd_walker(x1o,x2o,mrot) result(rmsd)
       real*8, intent(in) :: x1o(:,:), x2o(:,:)
       real*8, intent(out), optional :: mrot(3,3)
       real*8 :: rmsd
     end function rmsd_walker
     module subroutine gauleg (x1,x2,x,w,n)
       real*8, intent(in) :: x1
       real*8, intent(in) :: x2
       real*8, dimension(n), intent(out) :: x
       real*8, dimension(n), intent(out) :: w
       integer, intent(in) :: n
     end subroutine gauleg
     module subroutine bhole(rho,quad,hnorm,b,alf,prefac)
       real*8, intent(in) :: rho
       real*8, intent(in) :: quad
       real*8, intent(in) :: hnorm
       real*8, intent(out) :: b
       real*8, intent(out) :: alf
       real*8, intent(out) :: prefac
     end subroutine bhole
     module subroutine xlnorm(rho,quad,uxpos,xlnrm)
       real*8, intent(in) :: rho, quad, uxpos
       real*8, intent(out) :: xlnrm
     end subroutine xlnorm
     module function fdamp_tt(r,b,n)
       real*8, intent(in) :: r, b
       integer, intent(in) :: n
       real*8 :: fdamp_tt
     end function fdamp_tt
     module function fdamp_bj(r,rvdw,n)
       real*8, intent(in) :: r, rvdw
       integer, intent(in) :: n
       real*8 :: fdamp_bj
     end function fdamp_bj
     module subroutine munkres(n,c,as,cost)
       integer, intent(in) :: n
       real*8, intent(in) :: c(n,n)
       integer, intent(out) :: as(n)
       real*8, intent(out), optional :: cost
     end subroutine munkres
     module subroutine umeyama_graph_matching(n,ag,ah,perm)
       integer, intent(in) :: n
       real*8, intent(inout) :: ag(n,n)
       real*8, intent(inout) :: ah(n,n)
       integer, intent(out) :: perm(n)
     end subroutine umeyama_graph_matching
     module subroutine ullmann_graph_matching(iz1,ncon1,idcon1,iz2,ncon2,idcon2,nlist,list)
       integer, intent(in) :: iz1(:)
       integer, intent(in) :: ncon1(:)
       integer, intent(in) :: idcon1(:,:)
       integer, intent(in) :: iz2(:)
       integer, intent(in) :: ncon2(:)
       integer, intent(in) :: idcon2(:,:)
       integer, intent(out) :: nlist
       integer, allocatable, intent(inout) :: list(:,:)
     end subroutine ullmann_graph_matching
     module function invert_permutation(iperm)
       integer, intent(in) :: iperm(:)
       integer :: invert_permutation(size(iperm,1))
     end function invert_permutation
     !xx! lebedev submodule
     module subroutine good_lebedev(npts)
       integer, intent(inout) :: npts
     end subroutine good_lebedev
     module subroutine select_lebedev(npts,xleb,yleb,zleb,wleb)
       integer, intent(in) :: npts
       real*8, intent(out) :: xleb(:), yleb(:), zleb(:), wleb(:)
     end subroutine select_lebedev
     !xx! sobol submodule (not thread safe)
     module subroutine sobol_sequence(n,x,seed)
       integer, intent(in) :: n
       real*8, intent(out) :: x(n)
       integer*8, optional :: seed
     end subroutine sobol_sequence
  end interface

end module tools_math
