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
  real*8, parameter :: fire_drmax = 0.5d0 ! trust-region cap on the step displacements
  real*8, parameter :: fire_dt0 = 10d0 ! initial FIRE timestep (a.u.)
  real*8, parameter :: fire_dtmax = 40d0 ! ceiling for the adaptive FIRE timestep (a.u.)
  real*8, parameter :: fire_dtmin = 1d0 ! floor for the adaptive timestep (a.u.); stable for stiff O-H bonds

  ! enlarge the molecular cell when any atom comes within this fractional distance of a cell face.
  real*8, parameter :: md_enlarge_margin = 0.1d0

contains

  !> Initialize an MD/relaxation run for crystal c. Optional arguments
  !> set the energy backend, method, target temperature (K), timestep
  !> (a.u.), and run mode (md_dynamics or md_relax). errmsg is empty on success
  !> and holds the error message on failure.
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
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i

    errmsg = ""
    call md%free()
    call random_seed()
    if (present(temperature)) md%temperature = temperature
    if (present(dt)) md%dt = dt
    if (present(mode)) md%mode = mode

    call md%cl%init(c,backend=backend,method=method,errmsg=errmsg)
    if (len_trim(errmsg) > 0) return

    md%nat = c%ncel
    allocate(md%r(3,md%nat),md%r0(3,md%nat),md%v(3,md%nat),md%f(3,md%nat),md%mass(md%nat))
    do i = 1, md%nat
       md%r(:,i) = c%atcel(i)%r
       md%mass(i) = atmass(c%spc(c%atcel(i)%is)%z) * amu2au
    end do
    md%r0 = md%r

    ! initial velocities and forces
    md%fire_dt = fire_dt0
    md%fire_alpha = fire_alpha0
    md%fire_npos = 0
    md%istep = 0
    md%drag_iat = 0
    md%interacting = .false.
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
    ! keep the dragged atom clamped to the cursor with no residual velocity
    if (md%drag_iat > 0 .and. md%drag_iat <= md%nat) then
       md%r(:,md%drag_iat) = md%drag_target
       md%v(:,md%drag_iat) = 0d0
    end if
    ! recenter the molecule and grow its box if an atom nears a cell face
    call keep_in_cell(md,c)
    call kinetic(md)

  end subroutine md_step

  !> For molecules embedded in a display cell: translate the whole
  !> system so the center of mass sits at the cell center, and enlarge
  !> the molecular cell in place when any atom comes within
  !> md_enlarge_margin of a cell face.
  subroutine keep_in_cell(md,c)
    use crystalmod, only: crystal
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c

    integer :: i
    real*8 :: xcom(3), shift(3), xf(3)
    logical :: doenlarge

    if (.not.c%ismolecule .or. md%nat == 0) return

    ! while the user is dragging an atom or rigidly moving/rotating a molecule,
    ! do not recenter or resize (that would fight the manipulation)
    if (md%drag_iat > 0 .or. md%interacting) then
       call c%update_positions(md%r)
       return
    end if

    ! recenter the center of mass on the cell center (removes slow drift and
    ! keeps the molecule framed)
    xcom = 0d0
    do i = 1, md%nat
       xcom = xcom + md%mass(i) * md%r(:,i)
    end do
    xcom = xcom / sum(md%mass(1:md%nat))
    shift = c%x2c((/0.5d0,0.5d0,0.5d0/)) - xcom
    do i = 1, md%nat
       md%r(:,i) = md%r(:,i) + shift
    end do

    ! enlarge the cell if any atom sits within the margin of a face
    doenlarge = .false.
    do i = 1, md%nat
       xf = c%c2x(md%r(:,i))
       if (any(xf < md_enlarge_margin) .or. any(xf > 1d0 - md_enlarge_margin)) then
          doenlarge = .true.
          exit
       end if
    end do
    call c%update_positions(md%r)
    if (doenlarge) then
       call c%recompute_molecular_cell()
       do i = 1, md%nat
          md%r(:,i) = c%atcel(i)%r
       end do
    end if

  end subroutine keep_in_cell

  !> Restore the initial geometry and zero the velocities, writing the geometry
  !> back into c.
  module subroutine md_reset(md,c)
    use crystalmod, only: crystal
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c

    if (.not.md%ready) return
    md%r = md%r0
    md%v = 0d0
    md%drag_iat = 0
    md%interacting = .false.
    md%fire_dt = fire_dt0
    md%fire_alpha = fire_alpha0
    md%fire_npos = 0
    md%istep = 0
    call compute_forces(md,c)
    call kinetic(md)

  end subroutine md_reset

  !> (Re)draw velocities from the Maxwell-Boltzmann distribution at the target
  !> temperature and remove the center-of-mass drift.
  module subroutine md_init_velocities(md)
    class(mdrun), intent(inout) :: md

    integer :: i, k
    real*8 :: sigma

    do i = 1, md%nat
       sigma = sqrt(kboltz*md%temperature/md%mass(i))
       do k = 1, 3
          md%v(k,i) = sigma * gauss_random()
       end do
    end do

    ! remove the center-of-mass drift
    call remove_com(md)

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
    ! all buffers are allocated together in md_init
    if (allocated(md%r)) deallocate(md%r,md%r0,md%v,md%f,md%mass)
  end subroutine md_free

  !xx! private procedures

  !> Sync the crystal to the current positions, evaluate the force field, and
  !> store forces (including the interactive drag spring) and potential energy.
  subroutine compute_forces(md,c)
    use crystalmod, only: crystal
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c

    character(len=:), allocatable :: errmsg

    ! the dragged atom is clamped to the cursor position (a hard constraint, not
    ! a spring), so pin it before evaluating: neighbors then feel it exactly at
    ! the cursor and respond via the dynamics
    if (md%drag_iat > 0 .and. md%drag_iat <= md%nat) md%r(:,md%drag_iat) = md%drag_target

    ! evaluate the gradient straight into the persistent force buffer (no
    ! per-step allocation), then negate to get the forces
    call c%update_positions(md%r)
    ! periodically refresh the backend's neighbor list so contacts created by
    ! motion/dragging are captured (no-op for the tblite backend)
    md%istep = md%istep + 1
    if (md%nblist_every > 0 .and. mod(md%istep,md%nblist_every) == 0) &
       call md%cl%update_geometry(c)
    call md%cl%evaluate(c,md%epot,md%f,errmsg=errmsg)
    if (len_trim(errmsg) > 0) then
       md%ready = .false.
       return
    end if
    md%f = -md%f

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
    real*8 :: power, fnorm, vnorm, dmax
    logical :: capped

    ! velocity-Verlet with the adaptive FIRE timestep
    do i = 1, md%nat
       md%v(:,i) = md%v(:,i) + 0.5d0*md%fire_dt*md%f(:,i)/md%mass(i)
    end do
    ! trust region: scale the velocity so no atom moves more than fire_drmax this step
    dmax = 0d0
    do i = 1, md%nat
       dmax = max(dmax,norm2(md%v(:,i)))
    end do
    dmax = dmax * md%fire_dt
    capped = (dmax > fire_drmax)
    if (capped) md%v = md%v * (fire_drmax/dmax)
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

    if (power > 0d0 .and. .not.capped) then
       md%fire_npos = md%fire_npos + 1
       if (md%fire_npos > fire_nmin) then
          md%fire_dt = min(md%fire_dt*fire_finc,fire_dtmax)
          md%fire_alpha = md%fire_alpha*fire_falpha
       end if
    else
       md%fire_npos = 0
       md%fire_dt = max(md%fire_dt*fire_fdec,fire_dtmin)
       md%fire_alpha = fire_alpha0
       if (power <= 0d0) md%v = 0d0
    end if

    ! remove any net drift (relevant when a drag applies a net external force;
    ! a no-op for purely internal forces, which already sum to zero)
    call remove_com(md)

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
