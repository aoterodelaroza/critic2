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

!> Real-time molecular dynamics and geometry relaxations.
module dynamics
  use energy, only: calculator
  implicit none

  private

  public :: mdrun

  ! run modes
  integer, parameter, public :: md_dynamics = 0 !< constant-temperature dynamics
  integer, parameter, public :: md_relax = 1 !< geometry relaxation (FIRE)

  !> State of a running molecular-dynamics / relaxation simulation.
  type mdrun
     type(calculator) :: cl !< the energy/force backend
     logical :: ready = .false. !< true once %init has succeeded
     integer :: nat = 0 !< number of atoms
     integer :: mode = md_dynamics !< md_dynamics or md_relax
     real*8, allocatable :: r(:,:) !< positions (3,nat), bohr
     real*8, allocatable :: v(:,:) !< velocities (3,nat), bohr/a.u.
     real*8, allocatable :: f(:,:) !< forces (3,nat), hartree/bohr
     real*8, allocatable :: mass(:) !< atomic masses (nat), electron-mass units
     real*8 :: dt = 20d0 !< timestep, a.u. (~0.48 fs)
     real*8 :: temperature = 300d0 !< target temperature, K
     real*8 :: gamma = 2d-3 !< Langevin friction, 1/a.u.
     real*8 :: ekin = 0d0 !< last kinetic energy, hartree
     real*8 :: epot = 0d0 !< last potential energy, hartree
     ! interactive drag (harmonic spring atom -> cursor)
     integer :: drag_iat = 0 !< dragged atom (cell index), 0 = none
     real*8 :: drag_target(3) = 0d0 !< cursor target position, bohr
     real*8 :: drag_k = 0.02d0 !< drag spring constant, hartree/bohr^2
     ! FIRE relaxation state
     real*8 :: fire_alpha = 0.1d0
     real*8 :: fire_dt = 20d0
     integer :: fire_npos = 0
   contains
     procedure :: init => md_init
     procedure :: step => md_step
     procedure :: init_velocities => md_init_velocities
     procedure :: temperature_now => md_temperature
     procedure :: free => md_free
  end type mdrun

  interface
     module subroutine md_init(md,c,backend,method,temperature,dt,mode,errmsg)
       use crystalmod, only: crystal
       class(mdrun), intent(inout) :: md
       class(crystal), intent(inout) :: c
       integer, intent(in), optional :: backend
       integer, intent(in), optional :: method
       real*8, intent(in), optional :: temperature
       real*8, intent(in), optional :: dt
       integer, intent(in), optional :: mode
       character(len=:), allocatable, intent(out), optional :: errmsg
     end subroutine md_init
     module subroutine md_step(md,c)
       use crystalmod, only: crystal
       class(mdrun), intent(inout) :: md
       class(crystal), intent(inout) :: c
     end subroutine md_step
     module subroutine md_init_velocities(md)
       class(mdrun), intent(inout) :: md
     end subroutine md_init_velocities
     module function md_temperature(md) result(t)
       class(mdrun), intent(in) :: md
       real*8 :: t
     end function md_temperature
     module subroutine md_free(md)
       class(mdrun), intent(inout) :: md
     end subroutine md_free
  end interface

end module dynamics
