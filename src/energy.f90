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
!> geometry. This module is deliberately independent of any dynamics or GUI
!> code: it is meant to be reused by the interactive dynamics window, by
!> geometry relaxation, and by crystal-structure-prediction (CSP) searches.
!>
!> Two backends are provided behind a single interface:
!>  - ff_builtin: the built-in Universal Force Field (UFF, Rappe et al. 1992),
!>    with atom types derived from the connectivity graph and parameters read
!>    from dat/uff.dat. Covers bonds, angles, and van der Waals terms.
!>  - ff_tblite: an optional interface to the tblite library (GFN-FF / GFN2-xTB),
!>    only available when compiled with -DHAVE_TBLITE.
!>
!> All quantities are in atomic units (bohr, hartree).
module energy
  use iso_c_binding, only: c_ptr, c_null_ptr
  implicit none

  private

  public :: calculator

  ! energy backends
  integer, parameter, public :: ff_builtin = 0 !< built-in Universal Force Field (UFF)
  integer, parameter, public :: ff_tblite = 1 !< tblite library (GFN-FF/xTB)

  ! tblite methods (used only by the tblite backend)
  integer, parameter, public :: tbm_gfnff = 0 !< GFN-FF general force field
  integer, parameter, public :: tbm_gfn2 = 1 !< GFN2-xTB
  integer, parameter, public :: tbm_gfn1 = 2 !< GFN1-xTB

  !> UFF harmonic bond term (central pair i-j; j in image lvec).
  type bondterm
     integer :: i, j, lvec(3)
     real*8 :: r0 !< natural bond length (bohr)
     real*8 :: k !< force constant (hartree/bohr^2)
  end type bondterm

  !> UFF angle term (i-j-k, vertex j; i in image li, k in image lk).
  type angleterm
     integer :: i, j, k, li(3), lk(3)
     real*8 :: c0, c1, c2 !< precomputed Fourier coefficients
     real*8 :: kk !< force constant k_IJK (hartree)
  end type angleterm

  !> UFF van der Waals (Lennard-Jones 12-6) nonbonded pair (i-j; j in image lvec).
  type vdwterm
     integer :: i, j, lvec(3)
     real*8 :: rmin !< vdW minimum distance x_IJ (bohr)
     real*8 :: deps !< vdW well depth D_IJ (hartree)
     real*8 :: qscr2 = 0d0 !< squared Coulomb screening length a_IJ^2 (bohr^2), QEq electrostatics
     real*8 :: qij = 0d0 !< product of QEq charges q_i q_j (electrons^2)
  end type vdwterm

  !> UFF torsion term (i-j-k-l, central bond j-k; i,k,l in images li,lk,ll rel. to j).
  type torsterm
     integer :: i, j, k, l, li(3), lk(3), ll(3)
     real*8 :: v !< barrier (hartree)
     integer :: n !< periodicity
     real*8 :: cosf !< phase term cos(n*phi0), +1 or -1
  end type torsterm

  !> UFF inversion (out-of-plane) term (central c + 3 neighbors n1/n2/n3 in images l1/l2/l3).
  type invterm
     integer :: c, n1, n2, n3, l1(3), l2(3), l3(3)
     real*8 :: k !< force constant per apex (hartree)
     real*8 :: c0, c1, c2 !< inversion coefficients
  end type invterm

  !> Energy/force/stress calculator. Initialize with %init for a given crystal,
  !> then call %evaluate as many times as needed (positions are re-read from the
  !> crystal on each call). Release with %free.
  type calculator
     integer :: backend = ff_builtin !< selected backend (ff_*)
     integer :: method = tbm_gfnff !< tblite method (tbm_*), only for ff_tblite
     logical :: ready = .false. !< true once %init has run successfully
     integer :: nat = 0 !< number of atoms (= c%ncel)
     integer, allocatable :: attype(:) !< UFF atom-type index for each atom
     ! built-in UFF term lists
     integer :: nbond = 0, nang = 0, nnb = 0, ntor = 0, ninv = 0
     type(bondterm), allocatable :: bond(:)
     type(angleterm), allocatable :: ang(:)
     type(vdwterm), allocatable :: nb(:)
     type(torsterm), allocatable :: tor(:)
     type(invterm), allocatable :: inv(:)
     ! built-in UFF: QEq partial charges (electrostatics)
     logical :: do_elec = .false. !< whether the Coulomb/QEq term is active
     real*8, allocatable :: qeq(:) !< equilibrated partial charge per atom (electrons)
     ! tblite handles (only used by the tblite backend)
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
       character(len=:), allocatable, intent(out), optional :: errmsg
     end subroutine calc_init
     module subroutine calc_evaluate(cl,c,ene,grad,stress,errmsg)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
       real*8, intent(out) :: ene
       real*8, intent(out) :: grad(:,:)
       real*8, intent(out), optional :: stress(3,3)
       character(len=:), allocatable, intent(out), optional :: errmsg
     end subroutine calc_evaluate
     module subroutine calc_update_geometry(cl,c)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
     end subroutine calc_update_geometry
     module subroutine calc_free(cl)
       class(calculator), intent(inout) :: cl
     end subroutine calc_free
     ! tblite backend (implemented in energy@proc.F90; stubs when built without tblite)
     module subroutine calc_init_tblite(cl,c,errmsg)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
       character(len=:), allocatable, intent(out), optional :: errmsg
     end subroutine calc_init_tblite
     module subroutine calc_eval_tblite(cl,c,ene,grad,stress,errmsg)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
       real*8, intent(out) :: ene
       real*8, intent(out) :: grad(:,:)
       real*8, intent(out), optional :: stress(3,3)
       character(len=:), allocatable, intent(out), optional :: errmsg
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
