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

!> Qtree, basic routines.
module qtree_basic
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
  public :: get_tlengths

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
  integer, public :: minl

  ! qtree options
  real*8, allocatable, public :: sphfactor(:)
  real*8, allocatable, public :: sphintfactor(:)
  integer, public :: gradient_mode
  integer, public :: qtree_ode_mode
  real*8, public :: stepsize
  real*8, public :: ode_abserr
  integer, public :: integ_mode(20)
  integer, public :: integ_scheme
  integer, public :: keastnum
  integer, public :: plot_mode
  integer, public :: prop_mode
  integer, public :: mpstep
  real*8, public :: qtreefac
  real*8, public :: cub_abs
  real*8, public :: cub_rel
  integer, public :: cub_mpts
  logical, public :: docontacts
  real*8, public :: ws_origin(3)
  real*8, public :: ws_scale
  logical, public :: killext
  integer, public :: autosph
  logical, public :: checkbeta
  logical, public :: plotsticks
  integer, public :: color_allocate
  integer, public :: setsph_lvl
  real*8, public  :: vcutoff

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
  real*8, allocatable, public :: r_betaint(:)
  real*8, allocatable, public :: r_betagp(:)

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
     module subroutine get_tlengths(minlen,maxlen)
       real*8, intent(out) :: minlen, maxlen
     end subroutine get_tlengths
  end interface

end module qtree_basic
