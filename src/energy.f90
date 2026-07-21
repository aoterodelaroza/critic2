! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> Backend-agnostic energy/force/stress calculator for a crystal or molecular
!> geometry.
module energy
  use iso_c_binding, only: c_ptr, c_null_ptr
  use param, only: bohrtoa, kcal2ha
  implicit none

  private

  public :: calculator
  public :: ff_backend_applicable

  ! energy backends
  integer, parameter, public :: ff_uff = 0 !< built-in Universal Force Field (UFF)
  integer, parameter, public :: ff_tblite = 1 !< tblite library (GFN2/GFN1-xTB)
  integer, parameter, public :: ff_tip4p = 2 !< built-in TIP4P water model (molecules only)

  ! tblite methods (used only by the tblite backend)
  integer, parameter, public :: tbm_gfn2 = 1 !< GFN2-xTB
  integer, parameter, public :: tbm_gfn1 = 2 !< GFN1-xTB

  ! TIP4P water model (Jorgensen et al., J. Chem. Phys. 79 (1983) 926).
  ! The harmonic intramolecular restraints are NOT part of TIP4P (a
  ! rigid model): they keep the monomers together during
  ! MD/relaxation, and are zero at the reference geometry. Force
  ! constants are the harmonic limit of q-TIP4P/F (Habershon et al.,
  ! J. Chem. Phys. 131 (2009) 024501).
  real*8, parameter, public :: tip4p_qh = 0.52d0 !< H charge (electrons); M charge = -2*qh
  real*8, parameter, public :: tip4p_doh = 0.9572d0 / bohrtoa !< reference OH distance, bohr
  real*8, parameter, public :: tip4p_hoh = 104.52d0 !< reference HOH angle, degrees
  real*8, parameter, public :: tip4p_dom = 0.15d0 / bohrtoa !< O-M site distance, bohr
  real*8, parameter, public :: tip4p_lja = 6d5 * kcal2ha / bohrtoa**12 !< LJ A (kcal A^12/mol -> au)
  real*8, parameter, public :: tip4p_ljc = 610d0 * kcal2ha / bohrtoa**6 !< LJ C (kcal A^6/mol -> au)
  real*8, parameter, public :: tip4p_kbond = 2d0 * 116.09d0 * 2.287d0**2 * kcal2ha * bohrtoa**2 !< OH restraint (au)
  real*8, parameter, public :: tip4p_kang = 87.85d0 * kcal2ha !< HOH restraint (hartree/rad^2)

  !> UFF atom-type parameters, in native UFF units (Rappe et al., JACS 114
  !> (1992) 10024). The full table is hardwired in uff_atom_params
  !> (src/energy@proc.F90); an empty label means "no UFF parameters".
  type uffparam
     character(len=8) :: label = ""
     real*8 :: r1 = 0d0     ! valence bond radius (angstrom)
     real*8 :: theta0 = 0d0 ! equilibrium valence angle (degrees)
     real*8 :: x1 = 0d0     ! van der Waals distance (angstrom)
     real*8 :: dwell = 0d0  ! van der Waals well depth (kcal/mol)
     real*8 :: zstar = 0d0  ! effective charge
     real*8 :: vtor = 0d0   ! sp3 torsional barrier Vi (kcal/mol)
     real*8 :: utor = 0d0   ! sp2 torsional barrier Uj (kcal/mol)
     real*8 :: chi = 0d0    ! GMP electronegativity (eV)
     real*8 :: jhard = 0d0  ! GMP hardness / idempotential (eV)
     real*8 :: qrad = 0d0   ! QEq orbital radius (angstrom, for the Slater kernel)
  end type uffparam

  !> Harmonic bond term (central pair i-j; j in image lvec).
  type bondterm
     integer :: i, j, lvec(3)
     real*8 :: r0 !< natural bond length (bohr)
     real*8 :: k !< force constant (hartree/bohr^2)
  end type bondterm

  !> Angle term (i-j-k, vertex j; i in image li, k in image lk). The bending
  !> energy is kk * P(u) with u = cos(theta) and P a polynomial of degree <= 4,
  !> which represents both the general cosine-Fourier form and GULP's periodic
  !> special cases (theta0 = 0, 90, 120, 180 degrees).
  type angleterm
     integer :: i, j, k, li(3), lk(3)
     real*8 :: p(0:4) !< coefficients of P(u), u = cos(theta)
     real*8 :: kk !< force constant k_IJK (hartree)
  end type angleterm

  !> Van der Waals (Lennard-Jones 12-6) nonbonded pair (i-j; j in image lvec).
  type vdwterm
     integer :: i, j, lvec(3)
     real*8 :: rmin !< vdW minimum distance x_IJ (bohr)
     real*8 :: deps !< vdW well depth D_IJ (hartree)
  end type vdwterm

  !> Torsion term (i-j-k-l, central bond j-k; i,k,l in images li,lk,ll rel. to j).
  type torsterm
     integer :: i, j, k, l, li(3), lk(3), ll(3)
     real*8 :: v !< barrier (hartree)
     integer :: n !< periodicity
     real*8 :: cosf !< phase term cos(n*phi0), +1 or -1
  end type torsterm

  !> Inversion (out-of-plane) term (central c + 3 neighbors n1/n2/n3 in images l1/l2/l3).
  type invterm
     integer :: c, n1, n2, n3, l1(3), l2(3), l3(3)
     real*8 :: k !< force constant per apex (hartree)
     real*8 :: c0, c1, c2 !< inversion coefficients
  end type invterm

  !> Energy/force/stress calculator
  type calculator
     integer :: backend = ff_uff !< selected backend (ff_*)
     integer :: method = tbm_gfn2 !< tblite method (tbm_*), only for ff_tblite
     logical :: ready = .false. !< true once %init has run successfully
     !! UFF
     integer :: nat = 0 !< number of atoms (= c%ncel)
     type(uffparam), allocatable :: attype(:) !< UFF parameters for each atom
     ! term lists
     integer :: nbond = 0, nang = 0, nnb = 0, ntor = 0, ninv = 0
     type(bondterm), allocatable :: bond(:)
     type(angleterm), allocatable :: ang(:)
     type(vdwterm), allocatable :: nb(:)
     type(torsterm), allocatable :: tor(:)
     type(invterm), allocatable :: inv(:)
     ! QEq partial charges
     logical :: do_elec = .false. !< whether the Coulomb/QEq term is active
     real*8, allocatable :: qeq(:) !< equilibrated partial charge per atom (electrons)
     real*8, allocatable :: qzeta(:) !< QEq Slater orbital exponents (1/bohr) at the solved charges
     integer, allocatable :: qpqn(:) !< QEq principal quantum numbers
     ! cached Ewald parameters
     real*8 :: qeta = 0d0, qrcut = 0d0, qhcut = 0d0 !< Ewald width and real/recip cutoffs
     integer :: qlrmax(3) = 0, qlhmax(3) = 0 !< real/recip cell ranges
     ! cached electrostatics neighbour list (CSR), rebuilt on the vdW-list schedule
     integer :: nqnb = 0 !< total cached electrostatic pairs
     integer, allocatable :: qnbptr(:) !< CSR row pointers (nat+1)
     integer, allocatable :: qnbj(:) !< neighbour cell-atom index
     integer, allocatable :: qnblv(:,:) !< neighbour lattice vector (3, nqnb)
     !! TIP4P
     integer :: nwat = 0 !< number of water molecules
     integer, allocatable :: iwat(:,:) !< water atoms (3,nwat): O, H1, H2 cell indices
     !! tblite
     type(c_ptr) :: tb_ctx = c_null_ptr
     type(c_ptr) :: tb_mol = c_null_ptr
     type(c_ptr) :: tb_calc = c_null_ptr
     type(c_ptr) :: tb_res = c_null_ptr
   contains
     procedure :: init => calc_init
     procedure :: evaluate => calc_evaluate
     procedure :: update_geometry => calc_update_geometry
     procedure :: free => calc_free
  end type calculator

  interface
     module subroutine calc_init(cl,c,backend,method,elec,errmsg)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
       integer, intent(in), optional :: backend
       integer, intent(in), optional :: method
       logical, intent(in), optional :: elec
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine calc_init
     module subroutine calc_evaluate(cl,c,ene,grad,stress,errmsg)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
       real*8, intent(out) :: ene
       real*8, intent(out) :: grad(:,:)
       real*8, intent(out), optional :: stress(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine calc_evaluate
     module subroutine calc_update_geometry(cl,c)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
     end subroutine calc_update_geometry
     module subroutine calc_free(cl)
       class(calculator), intent(inout) :: cl
     end subroutine calc_free
     module function ff_backend_applicable(backend,c) result(ok)
       use crystalmod, only: crystal
       integer, intent(in) :: backend
       class(crystal), intent(in) :: c
       logical :: ok
     end function ff_backend_applicable
     ! tblite backend (implemented in energy@proc.F90; stubs when built without tblite)
     module subroutine calc_init_tblite(cl,c,errmsg)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine calc_init_tblite
     module subroutine calc_eval_tblite(cl,c,ene,grad,stress,errmsg)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
       real*8, intent(out) :: ene
       real*8, intent(out) :: grad(:,:)
       real*8, intent(out), optional :: stress(3,3)
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine calc_eval_tblite
     module subroutine calc_update_tblite(cl,c)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
     end subroutine calc_update_tblite
     module subroutine calc_free_tblite(cl)
       class(calculator), intent(inout) :: cl
     end subroutine calc_free_tblite
  end interface

end module energy
