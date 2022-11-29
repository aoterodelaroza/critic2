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

! Global variables, parameters and basic initialization subroutines.
! Contains global variables and parameters, and initialization
! procedures run at the beginning of the execution.
module global
  use param, only: bohrtoa, maxzat0
  implicit none

  public

  interface eval_next
     module procedure eval_next_real
     module procedure eval_next_int
  end interface eval_next

  ! Distance cutoffs
  real*8, parameter :: atomeps = 1d-2 !< Minimum distance to consider two atoms different
  real*8, parameter :: atomeps2 = atomeps*atomeps

  ! Beta-sphere default radius (bohr)
  real*8, parameter :: rbetadef = 0.1d0

  ! Environment variables
  character(len=:), allocatable :: critic_home !< CRITIC_HOME directory (with trailing "/")
  character(len=:), allocatable :: clib_file !< The path to the crystal library file
  character(len=:), allocatable :: mlib_file !< The path to the molecular library file

  ! Cutoff radius for 1d-12 atomic densities (max of r_LDA, r_PBE). Up
  ! to Pu (93). Values are bohr. This distance is usually enough to
  ! converge any density or orbital (dftb+, pi) contribution from an
  ! atom.
  real*8, parameter :: cutrad(maxzat0) = (/&
     2.149886192475d+01, 1.169139170668d+01, 3.430831385801d+01,& !001-003 H,  He, Li
     2.502075396007d+01, 2.801001395722d+01, 2.167675592180d+01,& !004-006 Be, B,  C
     1.749805708313d+01, 1.465173060207d+01, 1.263885024136d+01,& !007-009 N,  O,  F
     1.110521599057d+01, 3.523728402162d+01, 2.763528367271d+01,& !010-012 Ne, Na, Mg
     3.395549316507d+01, 2.847261278601d+01, 2.487715217494d+01,& !013-015 Al, Si, P
     2.222930087269d+01, 2.022231415676d+01, 1.857150175607d+01,& !016,018 S,  Cl, Ar
     3.884428729523d+01, 3.144587224767d+01, 2.970981796151d+01,& !019-021 K,  Ca, Sc
     2.864438811442d+01, 2.784088946336d+01, 2.925194799711d+01,& !022-024 Ti, V,  Cr
     2.660566532177d+01, 2.609916866690d+01, 2.558901439004d+01,& !025-027 Mn, Fe, Co
     2.517359152887d+01, 2.691554955610d+01, 2.435659411320d+01,& !028-030 Ni, Cu, Zn
     3.467478603212d+01, 2.914443825602d+01, 2.572575006996d+01,& !031-033 Ga, Ge, As
     2.323452863278d+01, 2.134146595122d+01, 1.981582897591d+01,& !034-036 Se, Br, Kr
     3.976877622180d+01, 3.266858263171d+01, 3.027851405458d+01,& !037-039 Rb, Sr, Y
     2.899491720657d+01, 2.967865003580d+01, 2.914637014504d+01,& !040-042 Zr, Nb, Mo
     2.697201600611d+01, 2.844039136970d+01, 2.814409350112d+01,& !043-045 Tc, Ru, Rh
     1.659926809140d+01, 2.771163603049d+01, 2.519886588880d+01,& !046-048 Pd, Ag, Cd
     3.538116802480d+01, 3.026454800251d+01, 2.695514982633d+01,& !049-051 In, Sn, Sb
     2.460024202780d+01, 2.277601677390d+01, 2.130658017554d+01,& !052-054 Te, I,  Xe
     4.137458546886d+01, 3.442036204804d+01, 3.242561614450d+01,& !055-057 Cs, Ba, La
     3.212250868201d+01, 3.306457792690d+01, 3.284026775197d+01,& !058-060 Ce, Pr, Nd
     3.262654222620d+01, 3.242292112974d+01, 3.222895504883d+01,& !061-063 Pm, Sm, Eu
     3.128431696800d+01, 3.180465714311d+01, 3.157435555252d+01,& !064-066 Gd, Tb, Dy
     3.135291924508d+01, 3.120231512615d+01, 3.099709515455d+01,& !067-069 Ho, Er, Tm
     3.079969409503d+01, 3.160515129459d+01, 2.709458469010d+01,& !070-072 Yb, Lu, Hf
     2.614193052742d+01, 2.548104664032d+01, 2.489113924347d+01,& !073-075 Ta, W,  Re
     2.441668377017d+01, 2.405143298004d+01, 2.466268008529d+01,& !076-078 Os, Ir, Pt
     2.439924398342d+01, 2.305709117567d+01, 3.643576493190d+01,& !079-081 Au, Hg, Tl
     3.110226831614d+01, 2.780342993946d+01, 2.541102668192d+01,& !082-084 Pb, Bi, Po
     2.360240806573d+01, 2.210165966674d+01, 4.053200388132d+01,& !085-087 At, Rn, Fr
     3.407838067822d+01, 3.585071927373d+01, 3.175945034367d+01,& !088-090 Ra, Ac, Th
     3.478340806986d+01, 3.489038964505d+01, 3.514212336660d+01,& !091-093 Pa, U,  Np
     3.120895952111d+01, 37d0, 37d0,&                             !094-096 Pu, Am, Cm
     37d0, 37d0, 37d0,&                                           !097-099 Bk, Cf, Es
     37d0, 37d0, 37d0,&                                           !100-102 Fm, Md, No
     37d0, 37d0, 37d0,&                                           !103-105 Lr, Rf, Db
     37d0, 37d0, 37d0,&                                           !106-108 Sg, Bh, Hs
     37d0, 37d0, 37d0,&                                           !109-111 Mt
     37d0, 37d0, 37d0,&                                           !112-114
     37d0, 37d0, 37d0,&                                           !115-117
     37d0,&                                                       !118
     0d0, 0d0, 0d0, 0d0, 0d0&                                     !119-123
     /) !< Cutoff radius for 1d-12 atomic densities (max of r_LDA, r_PBE)

  ! global control for critic
  character(len=:), allocatable :: fileroot !< file prefix
  logical :: testing = .false.
  logical :: quiet = .false.
  logical :: precisecube = .true.

  ! units
  integer :: iunit
  logical :: iunit_isdef

  integer, parameter :: iunit_bohr = 1
  integer, parameter :: iunit_ang = 2
  character*4, parameter :: iunitname0(2) = (/"bohr","ang_"/)
  real*8, parameter :: dunit0(2) = (/1d0,bohrtoa/)

  ! covalent bond factor
  real*8 :: bondfactor
  real*8, parameter :: bondfactor_def = 1.4d0

  ! default border for a molecular unit cell
  real*8, parameter :: rborder_def = 10d0 / bohrtoa

  ! guess and symmetry option (-1 = only for small systems, 0 = no, 1 = full)
  integer :: doguess = -1

  ! symmetry precision (spglib)
  real*8 :: symprec = 1d-2

  ! A crystal is considered small if it has less than this number of
  ! atoms in the unit cell.
  integer :: crsmall = 2000

  ! navigation options
  real*8 :: NAV_step !< gradient path step length (bohr)
  real*8 :: NAV_gradeps !< gradient path gradient mod termination threshold
  real*8 :: NAV_maxerr !< maximum error in gradient path tracing (bohr)
  integer :: NAV_stepper  !< the stepper in gp tracing
  integer, parameter :: NAV_stepper_euler = 1 !< Euler stepper (1 eval), poor-man's adaptive step
  integer, parameter :: NAV_stepper_rkck  = 2 !< Runge-Kutta-Cash-Karp embedded 4(5)-order, local extrapolation (6 eval), with error estimate
  integer, parameter :: NAV_stepper_dp    = 3 !< Dormand-prince embedded 4(5)-order, local extrapolation (7 eval), with error estimate
  integer, parameter :: NAV_stepper_bs    = 4 !< Bogacki-Shampine embedded 2(3)-order method, (5-1=4 eval, fsal), with error estimate
  integer, parameter :: NAV_stepper_heun  = 5 !< Heun stepper (2 eval), poor-man's adaptive step
  real*8 :: prunedist
  real*8, parameter :: gcpchange = 0.1d0 !< displacement from bcp/rcp when tracing bond/ring path

  ! critical points
  real*8 :: CP_hdegen = 1d-8 !< a CP is degenerate if any Hessian element is less than this value

  ! radial integration
  integer :: INT_radquad_type !< type of radial integration
  integer :: INT_radquad_nr !< number of nodes in radial integration
  real*8 :: INT_radquad_abserr !< req'd absolute error in radial integration
  real*8 :: INT_radquad_relerr !< req'd relative error in radial integration
  integer :: INT_radquad_errprop !< adaptive integration driver property
  logical :: INT_radquad_errprop_default !< has the property been set by the user?
  real*8 :: INT_iasprec !< precision of the IAS for radial integration
  integer, parameter :: INT_gauleg = 1 !< gauss legendre
  integer, parameter :: INT_qags = 2 !< quadpack qags
  integer, parameter :: INT_qng = 3 !< quadpack qng
  integer, parameter :: INT_qag = 4 !< quadpack qag
  integer, parameter :: INT_lebedev = 5 !< lebedev

  ! mesh type and quality for molecular integrations
  integer :: mesh_type !< type of mesh for molecular integrations (see meshmod)
  integer :: mesh_level !< level of mesh for molecular integrations (see meshmod)

  interface
     module subroutine critic_main()
     end subroutine critic_main
     module subroutine global_init(ghome,datadir)
       character*(*) :: ghome, datadir
     end subroutine global_init
     module subroutine global_set_defaults()
     end subroutine global_set_defaults
     module subroutine initial_banner()
     end subroutine initial_banner
     module subroutine help_me()
     end subroutine help_me
     module subroutine config_write()
     end subroutine config_write
     module subroutine critic_setvariables(line,lp)
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp
     end subroutine critic_setvariables
     module subroutine critic_clearvariable(line)
       character*(*), intent(in) :: line
     end subroutine critic_clearvariable
     module function eval_next_real(res,line,lp0)
       logical :: eval_next_real
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp0
       real*8, intent(out) :: res
     end function eval_next_real
     module function eval_next_int(res,line,lp0)
       logical :: eval_next_int
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp0
       integer, intent(out) :: res
     end function eval_next_int
     module subroutine list_radii()
     end subroutine list_radii
  end interface

end module global
