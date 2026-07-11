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
!>  - ff_builtin: a dependency-free classical force field (harmonic bonds and
!>    angles from the connectivity graph plus a Lennard-Jones nonbonded term).
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
  integer, parameter, public :: ff_builtin = 0 !< built-in classical force field
  integer, parameter, public :: ff_tblite = 1 !< tblite library (GFN-FF/xTB)

  ! tblite methods (used only by the tblite backend)
  integer, parameter, public :: tbm_gfnff = 0 !< GFN-FF general force field
  integer, parameter, public :: tbm_gfn2 = 1 !< GFN2-xTB
  integer, parameter, public :: tbm_gfn1 = 2 !< GFN1-xTB

  !> Energy/force/stress calculator. Initialize with %init for a given crystal,
  !> then call %evaluate as many times as needed (positions are re-read from the
  !> crystal on each call). Release with %free.
  type calculator
     integer :: backend = ff_builtin !< selected backend (ff_*)
     integer :: method = tbm_gfnff !< tblite method (tbm_*), only for ff_tblite
     logical :: ready = .false. !< true once %init has run successfully
     integer :: nat = 0 !< number of atoms (= c%ncel)
     ! built-in FF: harmonic bond terms
     integer :: nbond = 0
     integer, allocatable :: bond_i(:), bond_j(:) !< atcel indices of the pair
     integer, allocatable :: bond_lvec(:,:) !< lattice vector of atom j (3,nbond)
     real*8, allocatable :: bond_r0(:) !< equilibrium length (bohr)
     real*8, allocatable :: bond_k(:) !< force constant (hartree/bohr^2)
     ! built-in FF: harmonic angle terms (i-j-k, j is the vertex)
     integer :: nang = 0
     integer, allocatable :: ang_i(:), ang_j(:), ang_k(:) !< atcel indices
     integer, allocatable :: ang_li(:,:), ang_lk(:,:) !< lattice vectors of i and k
     real*8, allocatable :: ang_t0(:) !< equilibrium angle (radians)
     real*8, allocatable :: ang_kk(:) !< force constant (hartree/radian^2)
     ! built-in FF: Lennard-Jones nonbonded pairs
     integer :: nnb = 0
     integer, allocatable :: nb_i(:), nb_j(:) !< atcel indices
     integer, allocatable :: nb_lvec(:,:) !< lattice vector of atom j (3,nnb)
     real*8, allocatable :: nb_sig(:) !< LJ sigma (bohr)
     real*8, allocatable :: nb_eps(:) !< LJ epsilon (hartree)
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
     module subroutine calc_init(cl,c,backend,method,errmsg)
       use crystalmod, only: crystal
       class(calculator), intent(inout) :: cl
       class(crystal), intent(inout) :: c
       integer, intent(in), optional :: backend
       integer, intent(in), optional :: method
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
