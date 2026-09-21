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
  real*8, parameter :: md_drmax = 0.5d0 ! trust-region cap on per-step displacement shared by FIRE and Langevin
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
  module subroutine md_init(md,c,backend,method,temperature,dt,mode,eamfile,errmsg)
    use crystalmod, only: crystal
    use param, only: atmass
    class(mdrun), intent(inout) :: md
    class(crystal), intent(inout) :: c
    integer, intent(in), optional :: backend
    integer, intent(in), optional :: method
    real*8, intent(in), optional :: temperature
    real*8, intent(in), optional :: dt
    integer, intent(in), optional :: mode
    character(len=*), intent(in), optional :: eamfile
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i

    errmsg = ""
    call md%free()
    md%autostop = .true.
    if (present(temperature)) md%temperature = temperature
    if (present(dt)) md%dt = dt
    if (present(mode)) md%mode = mode

    call md%cl%init(c,backend=backend,method=method,eamfile=eamfile,errmsg=errmsg)
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
    md%simtime = 0d0
    md%stress = 0d0
    md%drag_iat = 0
    md%interacting = .false.
    if (md%mode == md_dynamics) then
       call md%init_velocities()
    else
       md%v = 0d0
    end if
    call compute_forces(md,c)
    if (allocated(md%errmsg)) then
       if (len_trim(md%errmsg) > 0) then
          errmsg = "force-field evaluation failed: " // md%errmsg
          return
       end if
    end if
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
       md%simtime = md%simtime + md%dt ! advance the clock by the step's timestep
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
    real*8 :: xcom(3), shift(3), xf(3), ri(3)
    logical :: doenlarge

    if (.not.c%ismolecule .or. md%nat == 0) return

    ! while the user is dragging an atom or rigidly moving/rotating a molecule,
    ! do not recenter or resize (that would fight the manipulation)
    if (md%drag_iat > 0 .or. md%interacting) then
       call c%update_positions(md%r)
       return
    end if

    ! recenter the center of mass on the cell center (removes slow drift and
    ! keeps the molecule framed). Held atoms are fixed in the cell frame, so
    ! a system that has any must not be translated under them
    if (.not.allocated(md%frozen)) then
       xcom = 0d0
       do i = 1, md%nat
          xcom = xcom + md%mass(i) * md%r(:,i)
       end do
       xcom = xcom / sum(md%mass(1:md%nat))
       shift = c%x2c((/0.5d0,0.5d0,0.5d0/)) - xcom
       do i = 1, md%nat
          md%r(:,i) = md%r(:,i) + shift
       end do
    end if

    ! enlarge the cell if any atom sits within the margin of a face
    doenlarge = .false.
    do i = 1, md%nat
       ! the ri temporary works around a gfortran miscompile of an array
       ! section of a polymorphic object passed to a type-bound function;
       ! see tools/check-gfortran-bug.md. Do NOT fold into c%c2x(md%r(:,i)).
       ri = md%r(:,i)
       xf = c%c2x(ri)
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

  !> Hold a subset of the atoms fixed for the rest of the run: frozen(i)
  !> true leaves atom i where it is and lets the others relax around it.
  !> Called with no argument, or with a mask that does not match the run
  !> or holds nothing, it releases every atom.
  module subroutine md_set_frozen(md,frozen)
    class(mdrun), intent(inout) :: md
    logical, intent(in), optional :: frozen(:)

    if (allocated(md%frozen)) deallocate(md%frozen)
    if (.not.present(frozen)) return
    if (size(frozen,1) /= md%nat) return
    if (.not.any(frozen)) return
    md%frozen = frozen
    call hold_frozen(md)

  end subroutine md_set_frozen

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
    md%simtime = 0d0
    md%stress = 0d0
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
    call hold_frozen(md)
    call remove_com(md)

  end subroutine md_init_velocities

  !> Rebase the run on the current geometry of c: the positions there become
  !> the restore point (md%r0) and the velocities are reseeded at the target
  !> temperature. Use before resuming a run that is still ready but whose
  !> restore point belongs to an older geometry.
  module subroutine md_rebase(md,c)
    use crystalmod, only: crystal
    class(mdrun), intent(inout) :: md
    class(crystal), intent(in) :: c

    integer :: i

    if (.not.md%ready .or. md%nat /= c%ncel) return
    do i = 1, md%nat
       md%r(:,i) = c%atcel(i)%r
    end do
    md%r0 = md%r
    md%istep = 0
    md%simtime = 0d0
    call md%init_velocities()

  end subroutine md_rebase

  !> Build a failure message for this run: base, followed by the dynamics'
  !> own error message if it has one.
  module function md_fail_message(md,base) result(msg)
    class(mdrun), intent(in) :: md
    character(len=*), intent(in) :: base
    character(len=:), allocatable :: msg

    msg = base
    if (allocated(md%errmsg)) then
       if (len_trim(md%errmsg) > 0) msg = msg // ": " // md%errmsg
    end if

  end function md_fail_message

  !> Instantaneous temperature (K) from the current kinetic energy.
  module function md_temperature(md) result(t)
    class(mdrun), intent(in) :: md
    real*8 :: t
    integer :: ndof
    t = 0d0
    if (md%nat <= 0) return
    ! free atoms only; the -3 removes the center of mass, which is fixed
    ! by the held atoms instead when there are any
    if (allocated(md%frozen)) then
       ndof = 3*count(.not.md%frozen)
    else
       ndof = 3*md%nat - 3
    end if
    ndof = max(ndof,1)
    t = 2d0*md%ekin/(real(ndof,8)*kboltz)
  end function md_temperature

  !> Local crystalline order of every atom, from the bond-orientational
  !> order (Steinhardt q6) of its neighbours within rcut (bohr): two
  !> neighbouring atoms share a "solid bond" when their normalized q6
  !> vectors overlap by more than 0.5 (ten Wolde and Frenkel), and
  !> order(i) is the fraction of solid bonds of atom i (0 = liquid-like, 1
  !> = crystalline). An atom with fewer than md_order_minnb neighbours
  !> has order 0. molten is the fraction of the atoms that are not held
  !> fixed with order below 0.5. Runs with the EAM backend only: the
  !> neighbours come from the pair list that backend caches (refreshed on
  !> the run's schedule), so the cost is just the spherical harmonics.
  !> errmsg is set if there is no such list or its cutoff is too short.
  module subroutine md_local_order(md,c,rcut,order,molten,errmsg)
    use crystalmod, only: crystal
    use tools_math, only: tosphere, genylm
    class(mdrun), intent(in) :: md
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rcut
    real*8, allocatable, intent(out) :: order(:)
    real*8, intent(out) :: molten
    character(len=:), allocatable, intent(out) :: errmsg

    integer, parameter :: md_order_minnb = 4 ! fewer neighbours than this = no local order
    real*8, parameter :: bond_thr = 0.5d0 ! q6 overlap above which a bond counts as solid
    real*8, parameter :: liquid_thr = 0.5d0 ! order below which an atom counts as liquid
    integer, parameter :: l6 = 6*7+1 ! genylm index of Y_6^0 (l*(l+1)+m+1)

    integer :: i, j, p, nsolid, nfree, nliq
    real*8 :: d(3), rcut2, r, tp(2), m(3,3)
    complex*16 :: y(49)
    complex*16, allocatable :: q(:,:)
    integer, allocatable :: nb(:)
    logical, allocatable :: keep(:)

    errmsg = ""
    molten = 0d0
    allocate(order(max(md%nat,0)))
    order = 0d0
    if (.not.md%ready .or. md%nat == 0) return
    if (.not.allocated(md%cl%enbptr) .or. .not.allocated(md%cl%enbj)) then
       errmsg = "local_order: no cached neighbour list (EAM backend only)"
       return
    end if
    if (size(md%cl%enbptr,1) < md%nat+1) then
       errmsg = "local_order: stale neighbour list"
       return
    end if
    if (md%cl%eam%rcut < rcut) then
       errmsg = "local_order: the potential cutoff is shorter than the neighbour cutoff"
       return
    end if
    rcut2 = rcut*rcut
    m = c%m_x2c

    ! q6 vector of every atom over its neighbours within rcut; keep marks
    ! the pairs of the cached list that are within rcut, for the bond pass
    allocate(q(-6:6,md%nat),nb(md%nat),keep(max(md%cl%nenb,1)))
    q = (0d0,0d0)
    nb = 0
    keep = .false.
    do i = 1, md%nat
       do p = md%cl%enbptr(i), md%cl%enbptr(i+1)-1
          j = md%cl%enbj(p)
          d = md%r(:,j) - md%r(:,i)
          if (.not.c%ismolecule) then
             if (any(md%cl%enblv(:,p) /= 0)) d = d + matmul(m,real(md%cl%enblv(:,p),8))
          end if
          if (dot_product(d,d) >= rcut2) cycle
          keep(p) = .true.
          nb(i) = nb(i) + 1
          call tosphere(d,r,tp)
          call genylm(6,tp,y)
          q(:,i) = q(:,i) + y(l6-6:l6+6)
       end do
       if (nb(i) >= md_order_minnb) then
          r = sqrt(sum(real(q(:,i)*conjg(q(:,i)),8)))
          if (r > 0d0) q(:,i) = q(:,i) / r
       end if
    end do

    ! solid bonds
    nfree = 0
    nliq = 0
    do i = 1, md%nat
       if (nb(i) >= md_order_minnb) then
          nsolid = 0
          do p = md%cl%enbptr(i), md%cl%enbptr(i+1)-1
             if (.not.keep(p)) cycle
             j = md%cl%enbj(p)
             if (nb(j) < md_order_minnb) cycle
             if (sum(real(q(:,i)*conjg(q(:,j)),8)) > bond_thr) nsolid = nsolid + 1
          end do
          order(i) = real(nsolid,8) / real(nb(i),8)
       end if
       ! molten fraction over the free atoms
       if (allocated(md%frozen)) then
          if (md%frozen(i)) cycle
       end if
       nfree = nfree + 1
       if (order(i) < liquid_thr) nliq = nliq + 1
    end do
    if (nfree > 0) molten = real(nliq,8) / real(nfree,8)

  end subroutine md_local_order

  !> Instantaneous pressure (GPa) of a periodic system: the kinetic term
  !> (2/3 E_kin/V) plus the virial term (-1/3 tr stress), from the stress cached
  !> by the last force evaluation. ok=.false. for a molecule (no cell).
  module function md_pressure(md,c,ok) result(p)
    use crystalmod, only: crystal
    use param, only: autogpa
    class(mdrun), intent(in) :: md
    class(crystal), intent(in) :: c
    logical, intent(out) :: ok
    real*8 :: p

    p = 0d0
    ok = (.not.c%ismolecule) .and. (c%omega > 0d0)
    if (.not.ok) return
    p = (2d0/3d0*md%ekin/c%omega - (md%stress(1,1)+md%stress(2,2)+md%stress(3,3))/3d0) * autogpa

  end function md_pressure

  !> Largest atomic force magnitude of the last evaluation, in
  !> eV/angstrom (0 if no forces are available).
  module function md_maxforce(md) result(f)
    use param, only: hartoev, bohrtoa
    class(mdrun), intent(in) :: md
    real*8 :: f

    f = 0d0
    if (md%nat > 0 .and. allocated(md%f)) &
       f = maxval(norm2(md%f,1)) * hartoev / bohrtoa

  end function md_maxforce

  !> Whether a relaxation run has converged: the maximum force is at
  !> or below the fconv threshold. Always false for dynamics runs.
  module function md_converged(md) result(conv)
    class(mdrun), intent(in) :: md
    logical :: conv

    conv = md%ready .and. md%mode == md_relax .and. md%nat > 0
    if (conv) conv = md%maxforce() <= md%fconv

  end function md_converged

  !> Release the run and its calculator.
  module subroutine md_free(md)
    class(mdrun), intent(inout) :: md
    if (allocated(md%frozen)) deallocate(md%frozen)
    call md%cl%free()
    md%ready = .false.
    md%errmsg = ""
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
    ! refresh the backend's neighbor list so contacts created by motion/dragging are captured
    md%istep = md%istep + 1
    if (md%nblist_every > 0 .and. mod(md%istep,md%nblist_every) == 0) then
       call c%update_env_after_move()
       call md%cl%update_geometry(c)
    end if
    ! cache the virial stress for crystals, on-screen pressure
    if (c%ismolecule) then
       call md%cl%evaluate(c,md%epot,md%f,errmsg=errmsg)
    else
       call md%cl%evaluate(c,md%epot,md%f,stress=md%stress,errmsg=errmsg)
    end if
    if (len_trim(errmsg) > 0) then
       md%errmsg = errmsg
       md%ready = .false.
       return
    end if
    md%errmsg = ""
    md%f = -md%f
    ! atoms held fixed feel no force (and carry no velocity)
    call hold_frozen(md)

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
    ! trust region: clamp per-atom speed so no atom drifts more than md_drmax this step
    ! step. Stops an atom dragged far outside the cell from blowing the integrator up
    call cap_velocity(md)
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
    ! the thermostat noise must not move the atoms held fixed
    call hold_frozen(md)
    ! A: half drift
    md%r = md%r + 0.5d0*md%dt*md%v
    ! force evaluation at the new positions
    call compute_forces(md,c)
    ! B: half kick
    do i = 1, md%nat
       md%v(:,i) = md%v(:,i) + 0.5d0*md%dt*md%f(:,i)/md%mass(i)
    end do
    ! clamp again so a pathological force cannot leave the reported kinetic
    ! energy / temperature (and the next step's drift) blown up
    call cap_velocity(md)

    ! remove the center-of-mass drift injected by the thermostat (keeps the
    ! system from slowly translating out of view and makes 3N-3 the correct
    ! number of degrees of freedom)
    call remove_com(md)

  end subroutine step_langevin

  !> Clamp each atom's speed so its drift over one timestep cannot exceed md_drmax.
  subroutine cap_velocity(md)
    class(mdrun), intent(inout) :: md
    integer :: i
    real*8 :: vmax, vn
    vmax = md_drmax / max(md%dt,1d-30)
    do i = 1, md%nat
       vn = norm2(md%v(:,i))
       if (vn > vmax) md%v(:,i) = md%v(:,i) * (vmax/vn)
    end do
  end subroutine cap_velocity

  !> Subtract the center-of-mass velocity so the whole system does not
  !> drift. Atoms held fixed pin the system already (and shifting every
  !> velocity would move them), so a run with a mask is left alone.
  subroutine remove_com(md)
    class(mdrun), intent(inout) :: md
    integer :: i
    real*8 :: pcom(3), mtot
    if (allocated(md%frozen)) return
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
    capped = (dmax > md_drmax)
    if (capped) md%v = md%v * (md_drmax/dmax)
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

    ! remove any net drift
    call remove_com(md)

  end subroutine step_fire

  !> Zero the force and the velocity of every atom held fixed.
  subroutine hold_frozen(md)
    class(mdrun), intent(inout) :: md

    integer :: i

    if (.not.allocated(md%frozen)) return
    do i = 1, md%nat
       if (.not.md%frozen(i)) cycle
       if (allocated(md%f)) md%f(:,i) = 0d0
       if (allocated(md%v)) md%v(:,i) = 0d0
    end do

  end subroutine hold_frozen

  !> Update the stored kinetic energy from the current velocities.
  subroutine kinetic(md)
    class(mdrun), intent(inout) :: md
    integer :: i
    md%ekin = 0d0
    do i = 1, md%nat
       md%ekin = md%ekin + 0.5d0*md%mass(i)*sum(md%v(:,i)**2)
    end do
  end subroutine kinetic

  !> Measure how expensive a molecular dynamics is on crystal c with
  !> the given force field: set up a run on a copy (so c is untouched)
  !> and take steps until nstepmax of them, or until tmax seconds have
  !> gone by, whichever comes first. Returns the number of steps timed
  !> (nstep), the setup time (tinit) and the average wall time per
  !> step (tstep), all in seconds.  Used to estimate the cost of a
  !> long run before committing to it.
  module subroutine md_benchmark(c,backend,method,temp,dt,nstepmax,tmax,nstep,tinit,tstep,errmsg,&
     eamfile)
    use crystalmod, only: crystal
    use energy, only: ff_backend_applicable, ff_backend_label
    type(crystal), intent(in) :: c
    integer, intent(in) :: backend, method, nstepmax
    real*8, intent(in) :: temp, dt, tmax
    integer, intent(out) :: nstep
    real*8, intent(out) :: tinit, tstep
    character(len=:), allocatable, intent(out) :: errmsg
    character(len=*), intent(in), optional :: eamfile

    type(mdrun) :: md
    type(crystal) :: cmd
    integer :: i, c0, c1, c2, rate

    nstep = 0
    tinit = 0d0
    tstep = 0d0
    errmsg = ""
    if (nstepmax < 1) then
       errmsg = "no steps to time"
       return
    end if
    if (.not.ff_backend_applicable(backend,c)) then
       errmsg = "the " // trim(ff_backend_label(backend,method)) //&
          " force field is not available for this system"
       return
    end if
    if (c%ncel == 0) then
       errmsg = "there are no atoms for the MD run"
       return
    end if

    ! time the setup and the steps on a copy, leaving c alone
    cmd = c
    call system_clock(count=c0,count_rate=rate)
    call md%init(cmd,backend=backend,method=method,temperature=temp,dt=dt,&
       mode=md_dynamics,eamfile=eamfile,errmsg=errmsg)
    call system_clock(count=c1)
    if (len_trim(errmsg) > 0) then
       call md%free()
       return
    end if

    ! stop early once the time budget is spent: a step is cheap with a force
    ! field and slow with a semiempirical hamiltonian, and the caller is
    ! waiting on this
    do i = 1, nstepmax
       call md%step(cmd)
       if (.not.md%ready) then
          errmsg = md%fail_message("force-field evaluation failed")
          call md%free()
          return
       end if
       nstep = i
       call system_clock(count=c2)
       if (rate > 0) then
          if (real(c2-c1,8) / real(rate,8) > tmax) exit
       end if
    end do
    call md%free()

    if (rate > 0 .and. nstep > 0) then
       tinit = real(c1-c0,8) / real(rate,8)
       tstep = real(c2-c1,8) / real(rate,8) / real(nstep,8)
    end if

  end subroutine md_benchmark

  !> Run an NVT (Langevin) molecular dynamics on a copy of crystal c
  !> with the given force field, temperature (K) and time step (a.u.),
  !> and append to seed(:) a seed for every snapshot: ngen of them,
  !> one every nstride steps, after nini equilibration steps. nseed is
  !> updated. On failure, errmsg.
  module subroutine bulk_md_seeds(c,backend,method,temp,dt,nini,ngen,nstride,seed,nseed,errmsg)
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed, realloc_crystalseed
    use energy, only: ff_backend_applicable, ff_backend_label
    use tools_io, only: uout, string
    use param, only: autofs
    type(crystal), intent(in) :: c
    integer, intent(in) :: backend, method, nini, ngen, nstride
    real*8, intent(in) :: temp, dt
    type(crystalseed), allocatable, intent(inout) :: seed(:)
    integer, intent(inout) :: nseed
    character(len=:), allocatable, intent(out) :: errmsg

    type(mdrun) :: md
    type(crystal) :: cmd
    character(len=:), allocatable :: ffname, snaplabel
    integer :: istep, ncoll, nsteptot, nprint
    logical :: collect

    errmsg = ""
    ffname = ff_backend_label(backend,method)
    if (.not.ff_backend_applicable(backend,c)) then
       errmsg = "the " // trim(ffname) // " force field is not available for this system"
       return
    end if
    if (c%ncel == 0) then
       errmsg = "there are no atoms for the MD run"
       return
    end if

    if (.not.allocated(seed)) allocate(seed(max(ngen,1)))

    ! run the MD on a copy so the current system geometry is preserved
    cmd = c
    call md%init(cmd,backend=backend,method=method,temperature=temp,dt=dt,&
       mode=md_dynamics,errmsg=errmsg)
    if (len_trim(errmsg) > 0) then
       call md%free() ! release any partial calculator state (e.g. tblite/xtb handles)
       return
    end if

    write (uout,'("+ Bulk configuration sampling by NVT molecular dynamics (WRITE BULK MD)")')
    write (uout,'("  force field: ",A)') trim(ffname)
    write (uout,'("  temperature (K): ",A)') string(temp,'f',decimal=2)
    write (uout,'("  time step (a.u. / fs): ",A," / ",A)') &
       string(dt,'f',decimal=2), string(dt*autofs,'f',decimal=4)
    write (uout,'("  equilibration steps: ",A)') string(nini)
    write (uout,'("  generation snapshots: ",A," (one every ",A," steps)")') &
       string(ngen), string(nstride)
    write (uout,'(/"   step         energy (Ha)        T (K)     snapshot")')

    nsteptot = nini + ngen*nstride
    nprint = max(1,nsteptot/50) ! keep the equilibration echo to ~50 lines
    ncoll = 0
    do istep = 1, nsteptot
       call md%step(cmd)
       if (.not.md%ready) then
          errmsg = md%fail_message("force-field evaluation failed during the MD run")
          call md%free()
          return
       end if

       ! collect a snapshot every nstride steps once equilibration is done
       collect = (istep > nini .and. mod(istep-nini,nstride) == 0 .and. ncoll < ngen)
       snaplabel = "-"
       if (collect) then
          ncoll = ncoll + 1
          if (nseed+1 > size(seed,1)) call realloc_crystalseed(seed,2*(nseed+1))
          call cmd%makeseed(seed(nseed+1),.false.)
          nseed = nseed + 1
          snaplabel = string(ncoll)
       end if

       if (collect .or. mod(istep,nprint) == 0 .or. istep == nsteptot) then
          write (uout,'(2X,I6,3X,A,3X,A,4X,A)') istep, &
             string(md%epot,'f',length=18,decimal=10), &
             string(md%temperature_now(),'f',length=8,decimal=1), snaplabel
       end if
    end do
    write (uout,'(/"  collected ",A," configurations")') string(ncoll)
    write (uout,*)

    call md%free()

  end subroutine bulk_md_seeds

end submodule proc
