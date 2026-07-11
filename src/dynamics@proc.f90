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

!> Implementation of the real-time molecular-dynamics / relaxation engine.
submodule (dynamics) proc
  use param, only: kboltz, amu2au
  use tools_math, only: gauss_random
  implicit none

  ! FIRE relaxation parameters
  real*8, parameter :: fire_finc = 1.1d0
  real*8, parameter :: fire_fdec = 0.5d0
  real*8, parameter :: fire_falpha = 0.99d0
  real*8, parameter :: fire_alpha0 = 0.1d0
  integer, parameter :: fire_nmin = 5

contains

  !> Initialize an MD/relaxation run for crystal c. Optional arguments
  !> set the energy backend, method, target temperature (K), timestep
  !> (a.u.), and run mode (md_dynamics or md_relax). On error, errmsg
  !> is and non-zero.
  module subroutine md_init(md,c,backend,method,temperature,dt,mode,errmsg)
    use crystalmod, only: crystal
    use param, only: atmass
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c
    integer, intent(in), optional :: backend
    integer, intent(in), optional :: method
    real*8, intent(in), optional :: temperature
    real*8, intent(in), optional :: dt
    integer, intent(in), optional :: mode
    character(len=:), allocatable, intent(out), optional :: errmsg

    integer :: i

    call md%free()
    call random_seed()
    if (present(temperature)) md%temperature = temperature
    if (present(dt)) md%dt = dt
    if (present(mode)) md%mode = mode

    call md%cl%init(c,backend=backend,method=method,errmsg=errmsg)
    if (present(errmsg)) then
       if (allocated(errmsg)) return
    end if

    md%nat = c%ncel
    allocate(md%r(3,md%nat),md%v(3,md%nat),md%f(3,md%nat),md%mass(md%nat))
    do i = 1, md%nat
       md%r(:,i) = c%atcel(i)%r
       md%mass(i) = atmass(c%spc(c%atcel(i)%is)%z) * amu2au
    end do

    ! initial velocities and forces
    md%fire_dt = md%dt
    md%fire_alpha = fire_alpha0
    md%fire_npos = 0
    if (md%mode == md_dynamics) then
       call md%init_velocities()
    else
       md%v = 0d0
    end if
    call compute_forces(md,c)
    call kinetic(md)
    md%ready = .true.

  end subroutine md_init

  !> Advance the simulation by one step and write the new geometry back into c.
  module subroutine md_step(md,c)
    use crystalmod, only: crystal
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c

    if (.not.md%ready) return
    if (md%mode == md_relax) then
       call step_fire(md,c)
    else
       call step_langevin(md,c)
    end if
    call kinetic(md)

  end subroutine md_step

  !> (Re)draw velocities from the Maxwell-Boltzmann distribution at the target
  !> temperature and remove the center-of-mass drift.
  module subroutine md_init_velocities(md)
    class(mdrun), intent(inout) :: md

    integer :: i, k
    real*8 :: sigma, pcom(3), mtot

    do i = 1, md%nat
       sigma = sqrt(kboltz*md%temperature/md%mass(i))
       do k = 1, 3
          md%v(k,i) = sigma * gauss_random()
       end do
    end do

    ! remove center-of-mass momentum
    pcom = 0d0
    mtot = 0d0
    do i = 1, md%nat
       pcom = pcom + md%mass(i)*md%v(:,i)
       mtot = mtot + md%mass(i)
    end do
    if (mtot > 0d0) then
       do i = 1, md%nat
          md%v(:,i) = md%v(:,i) - pcom/mtot
       end do
    end if

  end subroutine md_init_velocities

  !> Instantaneous temperature (K) from the current kinetic energy.
  module function md_temperature(md) result(t)
    class(mdrun), intent(in) :: md
    real*8 :: t
    integer :: ndof
    t = 0d0
    if (md%nat <= 0) return
    ndof = max(3*md%nat - 3,1)
    t = 2d0*md%ekin/(real(ndof,8)*kboltz)
  end function md_temperature

  !> Release the run and its calculator.
  module subroutine md_free(md)
    class(mdrun), intent(inout) :: md
    call md%cl%free()
    md%ready = .false.
    md%nat = 0
    md%ekin = 0d0
    md%epot = 0d0
    if (allocated(md%r)) deallocate(md%r)
    if (allocated(md%v)) deallocate(md%v)
    if (allocated(md%f)) deallocate(md%f)
    if (allocated(md%mass)) deallocate(md%mass)
  end subroutine md_free

  !xx! private procedures

  !> Sync the crystal to the current positions, evaluate the force field, and
  !> store forces (including the interactive drag spring) and potential energy.
  subroutine compute_forces(md,c)
    use crystalmod, only: crystal
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c

    real*8, allocatable :: grad(:,:)

    call c%update_positions(md%r)
    allocate(grad(3,md%nat))
    call md%cl%evaluate(c,md%epot,grad)
    md%f = -grad
    if (md%drag_iat > 0 .and. md%drag_iat <= md%nat) &
       md%f(:,md%drag_iat) = md%f(:,md%drag_iat) + md%drag_k*(md%drag_target - md%r(:,md%drag_iat))

  end subroutine compute_forces

  !> One BAOAB Langevin velocity-Verlet step at the target temperature.
  subroutine step_langevin(md,c)
    use crystalmod, only: crystal
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c

    integer :: i, k
    real*8 :: c1, c2, sig

    c1 = exp(-md%gamma*md%dt)

    ! B: half kick
    do i = 1, md%nat
       md%v(:,i) = md%v(:,i) + 0.5d0*md%dt*md%f(:,i)/md%mass(i)
    end do
    ! A: half drift
    md%r = md%r + 0.5d0*md%dt*md%v
    ! O: Ornstein-Uhlenbeck thermostat
    c2 = sqrt(1d0-c1*c1)
    do i = 1, md%nat
       sig = c2 * sqrt(kboltz*md%temperature/md%mass(i))
       do k = 1, 3
          md%v(k,i) = c1*md%v(k,i) + sig*gauss_random()
       end do
    end do
    ! A: half drift
    md%r = md%r + 0.5d0*md%dt*md%v
    ! force evaluation at the new positions
    call compute_forces(md,c)
    ! B: half kick
    do i = 1, md%nat
       md%v(:,i) = md%v(:,i) + 0.5d0*md%dt*md%f(:,i)/md%mass(i)
    end do

    ! remove the center-of-mass drift injected by the thermostat (keeps the
    ! system from slowly translating out of view and makes 3N-3 the correct
    ! number of degrees of freedom)
    call remove_com(md)

  end subroutine step_langevin

  !> Subtract the center-of-mass velocity so the whole system does not drift.
  subroutine remove_com(md)
    class(mdrun), intent(inout) :: md
    integer :: i
    real*8 :: pcom(3), mtot
    pcom = 0d0
    mtot = 0d0
    do i = 1, md%nat
       pcom = pcom + md%mass(i)*md%v(:,i)
       mtot = mtot + md%mass(i)
    end do
    if (mtot <= 0d0) return
    do i = 1, md%nat
       md%v(:,i) = md%v(:,i) - pcom/mtot
    end do
  end subroutine remove_com

  !> One FIRE relaxation step (damped velocity-Verlet driving forces to zero).
  subroutine step_fire(md,c)
    use crystalmod, only: crystal
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c

    integer :: i
    real*8 :: power, fnorm, vnorm, dtmax

    dtmax = 5d0*md%dt

    ! velocity-Verlet with the adaptive FIRE timestep
    do i = 1, md%nat
       md%v(:,i) = md%v(:,i) + 0.5d0*md%fire_dt*md%f(:,i)/md%mass(i)
    end do
    md%r = md%r + md%fire_dt*md%v
    call compute_forces(md,c)
    do i = 1, md%nat
       md%v(:,i) = md%v(:,i) + 0.5d0*md%fire_dt*md%f(:,i)/md%mass(i)
    end do

    ! FIRE velocity mixing
    power = sum(md%f*md%v)
    fnorm = sqrt(sum(md%f*md%f))
    vnorm = sqrt(sum(md%v*md%v))
    if (fnorm > 1d-30) &
       md%v = (1d0-md%fire_alpha)*md%v + md%fire_alpha*(vnorm/fnorm)*md%f

    if (power > 0d0) then
       md%fire_npos = md%fire_npos + 1
       if (md%fire_npos > fire_nmin) then
          md%fire_dt = min(md%fire_dt*fire_finc,dtmax)
          md%fire_alpha = md%fire_alpha*fire_falpha
       end if
    else
       md%fire_npos = 0
       md%fire_dt = md%fire_dt*fire_fdec
       md%fire_alpha = fire_alpha0
       md%v = 0d0
    end if

  end subroutine step_fire

  !> Update the stored kinetic energy from the current velocities.
  subroutine kinetic(md)
    class(mdrun), intent(inout) :: md
    integer :: i
    md%ekin = 0d0
    do i = 1, md%nat
       md%ekin = md%ekin + 0.5d0*md%mass(i)*sum(md%v(:,i)**2)
    end do
  end subroutine kinetic

end submodule proc
