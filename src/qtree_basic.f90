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

!> Qtree, basic routines.
module qtree_basic
  use global, only: mneq
  implicit none

  private

  public :: qtree_initialize
  public :: qtree_checksymmetry
  public :: qtree_cleanup
  public :: find_beta_rodriguez
  public :: map_ode_pointers
  public :: cindex
  public :: crys2convex
  public :: locate_tetrah
  public :: neargp
  private :: find_tetrah_contacts
  private :: inverse_operation
  public :: get_tlengths
  private :: presplit_ws

  ! logical units for differences and stick tess files
  integer, public :: ludif
  integer, public :: lustick(0:20)

  ! nuclear and nonnuclear maxima
  integer, public :: nnuc

  ! eps 
  real*8, parameter, public :: eps_v_warn = 1d-10        ! volume assert, warn
  real*8, parameter, public :: eps_v_periodic = 1d-6     ! volume assert, deact. periodic
  real*8, parameter, public :: eps_v_strict = 1d-3       ! volume assert, strict
  real*8, parameter, public :: eps_tetrah_contact = 1d-5 ! tetrah contact
  real*8, parameter, public :: tetra_eps_strict = 1d-2   ! x in tetrah, strict
  real*8, parameter, public :: tetra_eps_warn = 1d-5     ! x in tetrah, warn
  real*8, public :: crys2convex_eps, crys2convex_eps1 

  ! symmetry operations of the (0,0,0) local group
  integer, public :: leqv
  real*8, public :: leqvf
  real*8, public :: lrotm(3,3,48)

  ! number of initial tetrahedra
  integer, public :: nt_orig 

  ! origin, edge vectors and volume of the initial tetrahedra (cryst. coords.)
  real*8, allocatable, public :: tvol(:), cmat(:,:,:), dmat(:,:,:)
  real*8, allocatable, public :: torig(:,:), tvec(:,:,:)
  real*8, allocatable, public :: borig(:,:), bvec(:,:,:)
  integer, allocatable, public :: tcontact(:,:)
  integer, parameter, public :: tcontact_void = 999999
  real*8, public :: maxlen
  real*8, public :: minlen
  logical, public :: periodic

  ! level of recursive subdivision
  integer, public :: maxl

  ! grd call count
  integer, public :: nterm
  integer, public :: ngrd_term, ngrd_int
  integer, public :: nlocate, nlocate_sloppy

  ! recursive algorithm term list
  ! atom ids for the terms.:
  ! 1..nneq -- the term. is the corresponding atom
  ! -1..-nneq -- the same, but the point is inside the beta-sphere
  ! 0 -- unassigned
  ! mneq+1 -- unknown (it is a non-nuclear CP)
  integer, parameter, public :: qtreei = selected_int_kind(2)
  integer, parameter, public :: qtreeidx = selected_int_kind(14)
  integer, parameter, public :: qtreer = selected_real_kind(14)
  logical, public :: savefgr, savelapgr
  integer, public :: nder

  ! integration
  real*8, allocatable, public :: atprop(:,:)
  integer*4, public :: korder(10)
  real*8, public :: kxyz(3,45,10), kw(45,10)
  logical, public :: intcorner_deferred

  ! integration spheres
  real*8, public :: r_betaint(mneq)
  real*8, public :: r_betagp(mneq)

  ! gradient modes, comparison
  integer, public :: ndiff
  integer, public :: ngrd1, ngrd2

  ! gradient integration pointers
  integer, pointer, public :: ode_o, ode_o2, ode_type
  real*8, pointer, public :: ode_a(:,:)
  real*8, pointer, public :: ode_b(:), ode_b2(:)
  logical, public :: ode_fsal
  real*8, parameter, public :: safety = 0.9d0
  real*8, public :: qinv, q1inv

  ! Butcher tableaus of ODE solvers
  ! c_1 |
  ! c_2 | a_21
  ! c_3 | a_31  a_32
  ! ... | ...   ...
  ! c_o | a_o1  a_o2  ...  a_oo-1
  !------------------------------
  !     | b_1   b_2        b_o
  !
  ! k_i = f(t_n+c_i*h, y_n+ a_i1 * h * k_1+...+ a_ii-1 * h * k_i-1)
  ! y_n+1 = y_n + h * sum_i=1..o (  b_i * k_i )
  !
  ! The step is calculated with the b (not b2) formula
  !
  !! Euler method 
  integer, target, public :: euler_o = 1
  real*8, target, public :: euler_c(1) = (/0d0/)
  real*8, target, public :: euler_a(1,1) = reshape((/1d0/),shape(euler_a))
  real*8, target, public :: euler_b(1) = (/1d0/)
  integer, target, public :: euler_type = 0
  !! Heun's method
  integer, target, public :: heun_o = 2
  real*8, target, public :: heun_c(2) = (/0d0, 1d0/)
  real*8, target, public :: heun_a(2,2) = reshape((/0d0, 1d0,&
     1d0, 0d0/),shape(heun_a))
  real*8, target, public :: heun_b(2) = (/0.5d0, 0.5d0/)
  integer, target, public :: heun_type = 0
  !! Kutta's method
  integer, target, public :: kutta_o = 3
  real*8, target, public :: kutta_c(3) = (/0d0, 0.5d0, 1d0/)
  real*8, target, public :: kutta_a(3,3) = reshape((/0.0d0, 0.5d0, -1.0d0,&
     0.5d0, 0.0d0,  2.0d0,&
     -1.0d0, 2.0d0,  0.0d0/),shape(kutta_a))
  real*8, target, public :: kutta_b(3) = (/1d0/6d0, 2d0/3d0, 1d0/6d0/)
  integer, target, public :: kutta_type = 0
  !! RK4 method
  integer, target, public :: rk4_o = 4
  real*8, target, public :: rk4_c(4) = (/0d0, 0.5d0, 0.5d0, 1d0/)
  real*8, target, public :: rk4_a(4,4) = reshape((/0.0d0, 0.5d0, 0.0d0, 0.0d0,&
     0.5d0, 0.0d0, 0.5d0, 0.0d0,&
     0.0d0, 0.5d0, 0.0d0, 1.0d0,&
     0.0d0, 0.0d0, 1.0d0, 0.0d0/),shape(rk4_a))
  real*8, target, public :: rk4_b(4) = (/1d0/6d0, 1d0/3d0, 1d0/3d0, 1d0/6d0/)
  integer, target, public :: rk4_type = 0
  !! Heun-Euler embedded method
  !! 1st order, 2nd order error est., 2 evals
  integer, target, public :: heul_o = 2
  integer, target, public :: heul_o2 = 1
  real*8, target, public :: heul_c(2) = (/0d0, 1d0/)
  real*8, target, public :: heul_a(2,2) = reshape((/0d0, 0d0,&
     1d0, 0d0/),shape(heul_a))
  real*8, target, public :: heul_b2(2) = (/0.5d0, 0.5d0/)
  real*8, target, public :: heul_b(2) = (/1.0d0, 0.0d0/) 
  integer, target, public :: heul_type = 1
  !! Bogacki-Shampine embedded method 
  !! 3th order, 5th order error est., 3(4) evals
  !! Local extrapolation, fsal
  integer, target, public :: bs_o = 4
  integer, target, public :: bs_o2 = 2
  real*8, target, public :: bs_c(4) = (/0d0, 0.5d0, 0.75d0, 1d0/)
  real*8, target, public :: bs_a(4,4) = reshape(&
     (/0.0d0,0d0,0d0,0d0,&
     1d0/2d0,        0.0d0,0d0,0d0,&
     0.0d0,          3d0/4d0,     0.0d0,0d0,&
     2d0/9d0,        1d0/3d0,     4d0/9d0,       0.0d0&
     /),shape(bs_a))
  real*8, target, public :: bs_b(4) = (/2d0/9d0, 1d0/3d0, 4d0/9d0, 0d0/)
  real*8, target, public :: bs_b2(4) = (/7d0/24d0, 1d0/4d0, 1d0/3d0, 1d0/8d0/)
  integer, target, public :: bs_type = 1
  !! Runge-Kutta-Cash-Karp embedded method 
  !! 4th order, 5th order error est., 6 evals
  integer, target, public :: rkck_o = 6
  integer, target, public :: rkck_o2 = 4
  real*8, target, public :: rkck_c(6) = (/0d0, 1d0/5d0, 3d0/10d0, 3d0/5d0, 1d0, 7d0/8d0/)
  real*8, target, public :: rkck_a(6,6) = reshape(&
     (/0.0d0,0d0,0d0,0d0,0d0,0d0,&
     1d0/5d0,        0.0d0,0d0,0d0,0d0,0d0,&
     3d0/40d0,       9d0/40d0,    0.0d0,0d0,0d0,0d0,&
     3d0/10d0,      -9d0/10d0,    6d0/5d0,       0.0d0,0d0,0d0,&
     -11d0/54d0,      5d0/2d0,    -70d0/27d0,     35d0/27d0,        0.0d0,0d0,&
     1631d0/55296d0, 175d0/512d0, 575d0/13824d0, 44275d0/110592d0, 253d0/4096d0, 0.0d0&
     /),shape(rkck_a))
  real*8, target :: rkck_b(6) = (/2825d0/27648d0, 0d0, 18575d0/48384d0, 13525d0/55296d0,&
     277d0/14336d0, 1d0/4d0 /)
  real*8, target :: rkck_b2(6) = (/ 37d0/378d0, 0d0, 250d0/621d0, 125d0/594d0, 0d0, 512d0/1771d0 /)
  integer, target :: rkck_type = 1
  !! Dormand-Prince embedded method 
  !! 4th order, 5th order error est., 6(7) evals
  !! Local extrapolation, fsal
  integer, target, public :: dp_o = 7
  integer, target, public :: dp_o2 = 4
  real*8, target, public :: dp_c(7) = (/0d0, 1d0/5d0, 3d0/10d0, 4d0/5d0, 8d0/9d0, 1d0, 1d0/)
  real*8, target, public :: dp_a(7,7) = reshape(&
     (/0.0d0,  0d0,0d0,0d0,0d0,0d0,0d0,&
     1d0/5d0,         0.0d0,0d0,0d0,0d0,0d0,0d0,&
     3d0/40d0,        9d0/40d0,       0.0d0,0d0,0d0,0d0,0d0,&
     44d0/45d0,      -56d0/15d0,      32d0/9d0,        0d0,0d0,0d0,0d0,&
     19372d0/6561d0, -25360d0/2187d0, 64448d0/6561d0, -212d0/729d0,  0d0,0d0,0d0,&
     9017d0/3168d0,  -355d0/33d0,     46732d0/5247d0,  49d0/176d0,  -5103d0/18656d0, 0d0,0d0,&
     35d0/384d0,      0d0,            500d0/1113d0,    125d0/192d0, -2187d0/6784d0,  11d0/84d0,      0d0&
     /),shape(dp_a))
  real*8, target, public :: dp_b2(7) = (/5179d0/57600d0, 0d0, 7571d0/16695d0, 393d0/640d0,&
     -92097d0/339200d0, 187d0/2100d0, 1d0/40d0/)
  real*8, target, public :: dp_b(7) = (/ 35d0/384d0, 0d0, 500d0/1113d0, 125d0/192d0, &
     -2187d0/6784d0, 11d0/84d0, 0d0 /)
  integer, target, public :: dp_type = 1

  ! 3! permutations
  integer, parameter, public :: perm3(3,3,6) = reshape((/&
     1, 0, 0,&  ! 1 2 3
     0, 1, 0,&
     0, 0, 1,&
     0, 1, 0,&  ! 2 3 1
     0, 0, 1,&
     1, 0, 0,&
     0, 0, 1,&  ! 3 1 2
     1, 0, 0,&
     0, 1, 0,&
     1, 0, 0,&  ! 1 3 2
     0, 0, 1,&
     0, 1, 0,&
     0, 0, 1,&  ! 3 2 1
     0, 1, 0,&
     1, 0, 0,&
     0, 1, 0,&  ! 2 1 3
     1, 0, 0,&
     0, 0, 1 /),shape(perm3))
  
  interface
     module subroutine qtree_initialize(lvl,plvl,acum_atprop,trm,fgr,lapgr,vgr,verbose)
       integer, intent(in) :: lvl, plvl
       logical, intent(in) :: verbose
       integer(qtreei), allocatable, intent(out) :: trm(:,:)
       real(qtreer), allocatable, intent(out) :: fgr(:,:), lapgr(:,:), vgr(:)
       real*8, allocatable, intent(out) :: acum_atprop(:,:)
     end subroutine qtree_initialize
     module subroutine qtree_checksymmetry()
     end subroutine qtree_checksymmetry
     module subroutine qtree_cleanup()
     end subroutine qtree_cleanup
     module subroutine find_beta_rodriguez(nuc,rbeta)
       integer, intent(in) :: nuc
       real*8, intent(inout) :: rbeta
     end subroutine find_beta_rodriguez
     module subroutine map_ode_pointers(odem)
       integer, intent(in) :: odem
     end subroutine map_ode_pointers
     module function cindex(i,lvl)
       integer, intent(in) :: i(3)
       integer, intent(in) :: lvl
       integer(qtreeidx) :: cindex
     end function cindex
     module subroutine crys2convex(x,base_t,rver)
       real*8, intent(in) :: x(3)
       integer, intent(out) :: base_t
       real*8, intent(out) :: rver(3)
     end subroutine crys2convex
     module subroutine locate_tetrah(x,base_t,rver,lrot)
       real*8, intent(inout) :: x(3)
       integer, intent(out) :: base_t
       real*8, intent(out) :: rver(3)
       integer, intent(out) :: lrot
     end subroutine locate_tetrah
     module subroutine neargp(xp,base_t,lrot,idx,dist)
       real*8, intent(inout) :: xp(3)
       integer, intent(inout) :: base_t
       integer, intent(inout) :: lrot
       integer(qtreeidx), intent(out) :: idx
       real*8, intent(out) :: dist
     end subroutine neargp
     module subroutine find_tetrah_contacts()
     end subroutine find_tetrah_contacts
     module subroutine inverse_operation(p,op,c,invp,invop,invc)
       integer, intent(in) :: p, op, c
       integer, intent(out) :: invp, invop, invc
     end subroutine inverse_operation
     module subroutine get_tlengths(minlen,maxlen)
       real*8, intent(out) :: minlen, maxlen
     end subroutine get_tlengths
     module subroutine presplit_ws(plvl,ntetrag,tetrag)
       integer, intent(in) :: plvl
       integer, intent(inout) :: ntetrag
       real*8, intent(inout), allocatable :: tetrag(:,:,:)
     end subroutine presplit_ws
  end interface

end module qtree_basic
