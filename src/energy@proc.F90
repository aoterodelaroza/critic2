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

!> Implementation of the energy/force/stress calculator. Contains the built-in
!> Universal Force Field (UFF; Rappe et al. 1992) with connectivity-derived atom
!> types and parameters read from dat/uff.dat, and the tblite (GFN-FF / xTB) backend.
!> The tblite routines are compiled with real functionality only when built with
!> -DHAVE_TBLITE; otherwise they are stubs and the calculator falls back to the
!> built-in force field.
submodule (energy) proc
  use iso_c_binding
  use param, only: hartokjmol, hartoev
  implicit none

  ! UFF (Universal Force Field, Rappe et al., JACS 114 (1992) 10024) constants
  real*8, parameter :: rcut_lj = 12d0 !< van der Waals cutoff radius, bohr
  real*8, parameter :: uff_lambda = 0.1332d0 !< bond-order correction coefficient
  real*8, parameter :: uff_kpre = 664.12d0 !< force-constant prefactor, kcal.A/mol
  real*8, parameter :: kcal2ha = 4.184d0 / hartokjmol !< kcal/mol -> hartree

  !> One row of the UFF atom-type parameter table, in native UFF units.
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
  end type uffparam
  type(uffparam), allocatable, save :: uffpar(:) ! parameter table (from dat/uff.dat)
  integer, save :: nuffpar = 0
  logical, save :: uff_loaded = .false.

#ifdef HAVE_TBLITE
  ! Interfaces to the tblite C API (https://tblite.readthedocs.io/en/latest/api/c.html).
  ! Opaque handles are passed as C pointers. Array/optional pointers are passed
  ! by value as type(c_ptr) so that a NULL (c_null_ptr) or an actual address
  ! (c_loc(...)) can be given through the same argument.
  interface
     function c_tblite_new_error() bind(c,name="tblite_new_error") result(h)
       import c_ptr
       type(c_ptr) :: h
     end function c_tblite_new_error
     subroutine c_tblite_delete_error(h) bind(c,name="tblite_delete_error")
       import c_ptr
       type(c_ptr), intent(inout) :: h
     end subroutine c_tblite_delete_error
     function c_tblite_new_context() bind(c,name="tblite_new_context") result(h)
       import c_ptr
       type(c_ptr) :: h
     end function c_tblite_new_context
     subroutine c_tblite_delete_context(h) bind(c,name="tblite_delete_context")
       import c_ptr
       type(c_ptr), intent(inout) :: h
     end subroutine c_tblite_delete_context
     function c_tblite_new_structure(err,natoms,numbers,positions,charge,uhf,lattice,periodic) &
        bind(c,name="tblite_new_structure") result(h)
       import c_ptr, c_int
       type(c_ptr), value :: err
       integer(c_int), value :: natoms
       type(c_ptr), value :: numbers, positions, charge, uhf, lattice, periodic
       type(c_ptr) :: h
     end function c_tblite_new_structure
     subroutine c_tblite_update_structure_geometry(err,mol,positions,lattice) &
        bind(c,name="tblite_update_structure_geometry")
       import c_ptr
       type(c_ptr), value :: err, mol, positions, lattice
     end subroutine c_tblite_update_structure_geometry
     subroutine c_tblite_delete_structure(h) bind(c,name="tblite_delete_structure")
       import c_ptr
       type(c_ptr), intent(inout) :: h
     end subroutine c_tblite_delete_structure
     function c_tblite_new_gfn2_calculator(ctx,mol,config) &
        bind(c,name="tblite_new_gfn2_calculator") result(h)
       import c_ptr
       type(c_ptr), value :: ctx, mol, config
       type(c_ptr) :: h
     end function c_tblite_new_gfn2_calculator
     function c_tblite_new_gfn1_calculator(ctx,mol,config) &
        bind(c,name="tblite_new_gfn1_calculator") result(h)
       import c_ptr
       type(c_ptr), value :: ctx, mol, config
       type(c_ptr) :: h
     end function c_tblite_new_gfn1_calculator
     subroutine c_tblite_delete_calculator(h) bind(c,name="tblite_delete_calculator")
       import c_ptr
       type(c_ptr), intent(inout) :: h
     end subroutine c_tblite_delete_calculator
     function c_tblite_new_result() bind(c,name="tblite_new_result") result(h)
       import c_ptr
       type(c_ptr) :: h
     end function c_tblite_new_result
     subroutine c_tblite_delete_result(h) bind(c,name="tblite_delete_result")
       import c_ptr
       type(c_ptr), intent(inout) :: h
     end subroutine c_tblite_delete_result
     subroutine c_tblite_get_singlepoint(ctx,mol,calc,res) bind(c,name="tblite_get_singlepoint")
       import c_ptr
       type(c_ptr), value :: ctx, mol, calc, res
     end subroutine c_tblite_get_singlepoint
     subroutine c_tblite_get_result_energy(err,res,energy) bind(c,name="tblite_get_result_energy")
       import c_ptr, c_double
       type(c_ptr), value :: err, res
       real(c_double), intent(out) :: energy
     end subroutine c_tblite_get_result_energy
     subroutine c_tblite_get_result_gradient(err,res,gradient) bind(c,name="tblite_get_result_gradient")
       import c_ptr, c_double
       type(c_ptr), value :: err, res
       real(c_double), intent(out) :: gradient(*)
     end subroutine c_tblite_get_result_gradient
     subroutine c_tblite_get_result_virial(err,res,sigma) bind(c,name="tblite_get_result_virial")
       import c_ptr, c_double
       type(c_ptr), value :: err, res
       real(c_double), intent(out) :: sigma(*)
     end subroutine c_tblite_get_result_virial
  end interface
#endif

contains

  !> Initialize the calculator for crystal c. backend selects ff_builtin (default)
  !> or ff_tblite; method selects the tblite method. elec toggles the built-in
  !> UFF QEq electrostatics (default on). On error, errmsg is allocated.
  module subroutine calc_init(cl,c,backend,method,elec,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    integer, intent(in), optional :: backend
    integer, intent(in), optional :: method
    logical, intent(in), optional :: elec
    character(len=:), allocatable, intent(out), optional :: errmsg

    call cl%free()
    if (present(backend)) cl%backend = backend
    if (present(method)) cl%method = method
    cl%do_elec = .true.
    if (present(elec)) cl%do_elec = elec
    cl%nat = c%ncel

    if (cl%backend == ff_builtin) then
       call uff_setup(cl,c,errmsg)
       cl%ready = .not.(present(errmsg) .and. allocated(errmsg))
    else if (cl%backend == ff_tblite) then
       call calc_init_tblite(cl,c,errmsg)
       cl%ready = .not.(present(errmsg) .and. allocated(errmsg))
    else
       if (present(errmsg)) errmsg = "unknown energy backend"
    end if

  end subroutine calc_init

  !> Evaluate the energy (ene, hartree), energy gradient (grad(3,nat),
  !> hartree/bohr; force = -grad) and, optionally, the stress tensor
  !> (stress(3,3), hartree/bohr^3) at the current geometry of c.
  module subroutine calc_evaluate(cl,c,ene,grad,stress,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    real*8, intent(out) :: ene
    real*8, intent(out) :: grad(:,:)
    real*8, intent(out), optional :: stress(3,3)
    character(len=:), allocatable, intent(out), optional :: errmsg

    ene = 0d0
    grad = 0d0
    if (present(stress)) stress = 0d0

    if (.not.cl%ready) then
       if (present(errmsg)) errmsg = "calculator not initialized"
       return
    end if

    if (cl%backend == ff_builtin) then
       call uff_evaluate(cl,c,ene,grad,stress)
    else if (cl%backend == ff_tblite) then
       call calc_eval_tblite(cl,c,ene,grad,stress,errmsg)
    else
       if (present(errmsg)) errmsg = "unknown energy backend"
    end if

  end subroutine calc_evaluate

  !> Refresh backend caches after the geometry of c changed appreciably. For the
  !> built-in FF this rebuilds the nonbonded neighbor list; bonds and angles
  !> (which define the topology and equilibrium values) are kept.
  module subroutine calc_update_geometry(cl,c)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c

    if (.not.cl%ready) return
    if (cl%backend == ff_builtin) then
       call uff_nonbonded(cl,c)
    else if (cl%backend == ff_tblite) then
       call calc_update_tblite(cl,c)
    end if

  end subroutine calc_update_geometry

  !> Release all calculator resources.
  module subroutine calc_free(cl)
    class(calculator), intent(inout) :: cl

    if (cl%backend == ff_tblite) call calc_free_tblite(cl)
    cl%ready = .false.
    cl%nat = 0
    cl%nbond = 0
    cl%nang = 0
    cl%nnb = 0
    cl%ntor = 0
    cl%ninv = 0
    cl%do_elec = .false.
    if (allocated(cl%attype)) deallocate(cl%attype)
    if (allocated(cl%qeq)) deallocate(cl%qeq)
    if (allocated(cl%bond)) deallocate(cl%bond)
    if (allocated(cl%ang)) deallocate(cl%ang)
    if (allocated(cl%nb)) deallocate(cl%nb)
    if (allocated(cl%tor)) deallocate(cl%tor)
    if (allocated(cl%inv)) deallocate(cl%inv)

  end subroutine calc_free

  !xx! private procedures

  !> Build the UFF topology (bonds, angles, van der Waals pairs) from the
  !> connectivity graph. Atom types are assigned from the connectivity and all
  !> parameters converted to atomic units. On error (missing dat/uff.dat or an
  !> unparameterized element) errmsg is allocated.
  subroutine uff_setup(cl,c,errmsg)
    use crystalmod, only: crystal
    use types, only: neighstar
    use param, only: bohrtoa, rad
    use tools_io, only: string
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out), optional :: errmsg

    type(neighstar), allocatable :: ns(:)
    integer :: i, j, k, ka, kb, ni, nk, nb, na, pass
    integer :: lv(3), li(3), lk(3)
    real*8 :: rij, rjk, rik, t0, cost0, st0sq

    ! load the parameter table (dat/uff.dat) once
    call uff_load(errmsg)
    if (present(errmsg)) then
       if (allocated(errmsg)) return
    end if

    ! connectivity graph (use the crystal's own if available)
    if (allocated(c%nstar)) then
       ns = c%nstar
    else
       call c%find_asterisms(ns)
    end if

    ! assign a UFF atom type to every atom
    allocate(cl%attype(cl%nat))
    do i = 1, cl%nat
       cl%attype(i) = uff_type_index(c,i,ns)
       if (cl%attype(i) <= 0) then
          if (present(errmsg)) errmsg = "no UFF parameters for atomic number " // &
             string(c%spc(c%atcel(i)%is)%z)
          return
       end if
    end do

    ! ---- bonds ----  (each undirected bond stored once, canonically)
    do pass = 1, 2
       nb = 0
       do i = 1, cl%nat
          do k = 1, ns(i)%ncon
             j = ns(i)%idcon(k)
             lv = ns(i)%lcon(:,k)
             if (.not.bond_canonical(i,j,lv)) cycle
             nb = nb + 1
             if (pass == 2) then
                cl%bond(nb)%i = i
                cl%bond(nb)%j = j
                cl%bond(nb)%lvec = lv
                rij = uff_natural_length(cl%attype(i),cl%attype(j),bond_order(ns(i),k)) ! angstrom
                cl%bond(nb)%r0 = rij / bohrtoa                                          ! bohr
                ! k_IJ = 664.12 Z*_I Z*_J / r_IJ^3 [kcal/mol/A^2] -> hartree/bohr^2
                cl%bond(nb)%k = uff_kpre * uffpar(cl%attype(i))%zstar * uffpar(cl%attype(j))%zstar &
                   / rij**3 * kcal2ha * bohrtoa*bohrtoa
             end if
          end do
       end do
       if (pass == 1) then
          cl%nbond = nb
          allocate(cl%bond(nb))
       end if
    end do

    ! ---- angles (i-j-k, vertex j) ----  the Fourier coefficients depend only on
    ! theta0 and are precomputed here (linear centers fold into c0=1,c1=1,c2=0)
    na = 0
    do j = 1, cl%nat
       na = na + ns(j)%ncon*(ns(j)%ncon-1)/2
    end do
    cl%nang = na
    allocate(cl%ang(na))
    na = 0
    do j = 1, cl%nat
       t0 = uffpar(cl%attype(j))%theta0 * rad
       cost0 = cos(t0)
       st0sq = 1d0 - cost0*cost0
       do ka = 1, ns(j)%ncon
          ni = ns(j)%idcon(ka)
          li = ns(j)%lcon(:,ka)
          do kb = ka+1, ns(j)%ncon
             nk = ns(j)%idcon(kb)
             lk = ns(j)%lcon(:,kb)
             na = na + 1
             cl%ang(na)%i = ni
             cl%ang(na)%j = j
             cl%ang(na)%k = nk
             cl%ang(na)%li = li
             cl%ang(na)%lk = lk
             if (st0sq < 1d-8) then
                ! linear center: E = k (1 + cos theta) = k (c0 + c1 u + c2(2u^2-1))
                cl%ang(na)%c0 = 1d0; cl%ang(na)%c1 = 1d0; cl%ang(na)%c2 = 0d0
             else
                cl%ang(na)%c2 = 1d0/(4d0*st0sq)
                cl%ang(na)%c1 = -4d0*cl%ang(na)%c2*cost0
                cl%ang(na)%c0 = cl%ang(na)%c2*(2d0*cost0*cost0 + 1d0)
             end if
             ! natural bond lengths of the two arms (angstrom); r_IK by law of cosines
             rij = uff_natural_length(cl%attype(ni),cl%attype(j),bond_order(ns(j),ka))
             rjk = uff_natural_length(cl%attype(nk),cl%attype(j),bond_order(ns(j),kb))
             rik = sqrt(max(rij*rij + rjk*rjk - 2d0*rij*rjk*cost0,1d-10))
             ! k_IJK = 664.12 (Z*_I Z*_K/r_IK^5) r_IJ r_JK [3 r_IJ r_JK (1-cos^2 t0) - r_IK^2 cos t0]
             cl%ang(na)%kk = uff_kpre * uffpar(cl%attype(ni))%zstar * uffpar(cl%attype(nk))%zstar &
                / rik**5 * rij * rjk * (3d0*rij*rjk*(1d0-cost0*cost0) - rik*rik*cost0) * kcal2ha
          end do
       end do
    end do

    ! ---- torsions and inversions ----
    call uff_build_torsions(cl,c,ns)
    call uff_build_inversions(cl,c,ns)

    ! ---- QEq partial charges (before the pair list, which caches q_i q_j) ----
    if (cl%do_elec) call uff_qeq(cl,c)

    ! ---- van der Waals + electrostatics pairs ----
    call uff_nonbonded(cl,c)

  end subroutine uff_setup

  !> Build the UFF van der Waals (LJ 12-6) pair list within rcut_lj, excluding
  !> 1-2 (bonded) and 1-3 (share a common neighbor) pairs. Uses the atom types
  !> stored in cl%attype (set by uff_setup).
  subroutine uff_nonbonded(cl,c)
    use crystalmod, only: crystal
    use types, only: neighstar
    use param, only: icrd_crys, bohrtoa
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c

    type(neighstar), allocatable :: nslocal(:)
    integer :: i, j, m, nat, nnb, pass
    integer, allocatable :: eid(:), lvec(:,:)
    real*8, allocatable :: dist(:)
    integer :: lv(3)
    real*8 :: xij, dij
    logical :: havens, is13, haveq

    if (.not.allocated(cl%attype)) return

    ! read the connectivity in place (this routine is re-run periodically during
    ! MD, so avoid deep-copying the whole neighbor-star array each time)
    havens = allocated(c%nstar)
    if (.not.havens) call c%find_asterisms(nslocal)
    haveq = cl%do_elec .and. allocated(cl%qeq) ! charges available -> cache q_i q_j

    if (allocated(cl%nb)) deallocate(cl%nb)

    do pass = 1, 2
       nnb = 0
       do i = 1, cl%nat
          call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.false.,nat,eid=eid,dist=dist,&
             lvec=lvec,up2d=rcut_lj,nozero=.true.)
          do m = 1, nat
             j = eid(m)
             lv = lvec(:,m)
             if (.not.bond_canonical(i,j,lv)) cycle
             if (havens) then
                is13 = excluded_13(c%nstar,i,j,lv)
             else
                is13 = excluded_13(nslocal,i,j,lv)
             end if
             if (is13) cycle
             nnb = nnb + 1
             if (pass == 2) then
                cl%nb(nnb)%i = i
                cl%nb(nnb)%j = j
                cl%nb(nnb)%lvec = lv
                ! UFF combining rules: geometric mean of distance and well depth
                xij = sqrt(uffpar(cl%attype(i))%x1 * uffpar(cl%attype(j))%x1)     ! angstrom
                dij = sqrt(uffpar(cl%attype(i))%dwell * uffpar(cl%attype(j))%dwell) ! kcal/mol
                cl%nb(nnb)%rmin = xij / bohrtoa ! bohr
                cl%nb(nnb)%deps = dij * kcal2ha ! hartree
                ! electrostatics: cache the screening length^2 and the charge product
                cl%nb(nnb)%qscr2 = qeq_ascreen(cl%attype(i),cl%attype(j))**2
                if (haveq) then
                   cl%nb(nnb)%qij = cl%qeq(i)*cl%qeq(j)
                else
                   cl%nb(nnb)%qij = 0d0
                end if
             end if
          end do
       end do
       if (pass == 1) then
          cl%nnb = nnb
          allocate(cl%nb(nnb))
       end if
    end do

  end subroutine uff_nonbonded

  !> Equilibrate partial charges by electronegativity equalization (QEq/EEM):
  !> solve chi_i + sum_j J_ij q_j = mu with sum q = 0, where J_ii is the atomic
  !> hardness and J_ij an Ohno-screened Coulomb kernel (summed over periodic
  !> images within a cutoff). Charges are stored in cl%qeq; on a singular system
  !> the electrostatics is disabled.
  subroutine uff_qeq(cl,c)
    use crystalmod, only: crystal
    use param, only: icrd_crys
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c

    integer :: n, i, j, m, nat, info
    integer, allocatable :: eid(:), ipiv(:)
    real*8, allocatable :: dist(:), amat(:,:), bvec(:)
    real*8 :: ascr, rr
    real*8, parameter :: qcut = 20d0 ! charge-equilibration cutoff (bohr)

    n = cl%nat
    if (n <= 0) then
       cl%do_elec = .false.
       return
    end if

    allocate(amat(n+1,n+1),bvec(n+1),ipiv(n+1))
    amat = 0d0
    bvec = 0d0
    ! diagonal idempotential J_ii = 2*hardness and the electronegativity
    ! right-hand side (atomic units)
    do i = 1, n
       amat(i,i) = 2d0*uffpar(cl%attype(i))%jhard / hartoev
       bvec(i) = -uffpar(cl%attype(i))%chi / hartoev
    end do
    ! off-diagonal (and self-image) screened Coulomb, summed over images in qcut
    do i = 1, n
       call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.false.,nat,eid=eid,dist=dist,&
          up2d=qcut,nozero=.true.)
       do m = 1, nat
          j = eid(m)
          rr = dist(m)
          if (rr < 1d-10) cycle
          ascr = qeq_ascreen(cl%attype(i),cl%attype(j))
          amat(i,j) = amat(i,j) + 1d0/sqrt(rr*rr + ascr*ascr)
       end do
    end do
    ! charge-conservation constraint (sum q = 0) via a Lagrange multiplier
    do i = 1, n
       amat(i,n+1) = -1d0
       amat(n+1,i) = 1d0
    end do

    call dgesv(n+1,1,amat,n+1,ipiv,bvec,n+1,info)
    if (info /= 0) then
       ! singular (e.g. all-hard-zero); skip electrostatics
       cl%do_elec = .false.
       return
    end if

    if (allocated(cl%qeq)) deallocate(cl%qeq)
    allocate(cl%qeq(n))
    cl%qeq(1:n) = bvec(1:n)

  end subroutine uff_qeq

  !> Evaluate the UFF energy, gradient, and (optionally) stress, all in atomic
  !> units. Terms: harmonic bond stretch, cosine angle bend, Lennard-Jones 12-6
  !> van der Waals, screened Coulomb (QEq charges), torsions, and inversions.
  subroutine uff_evaluate(cl,c,ene,grad,stress)
    use crystalmod, only: crystal
    use tools_math, only: cross
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    real*8, intent(out) :: ene
    real*8, intent(out) :: grad(:,:)
    real*8, intent(out), optional :: stress(3,3)

    integer :: n, i, j, ii, jj, kk
    real*8 :: d(3), r, dr, gmag, g(3)
    real*8 :: ri(3), rj(3), rk(3), a(3), b(3), na, nb, u, edu
    real*8 :: dudi(3), dudk(3), dudj(3), gi(3), gj(3), gk(3)
    real*8 :: sr, sr2, sr6, sr12, fmag
    real*8 :: rr, rinv, gq(3)
    real*8 :: vir(3,3)
    logical :: dostress
    ! torsion locals
    integer :: t, nn, ti, tj, tk, tl
    real*8 :: rl(3), tv1(3), tv2(3), tv3(3), cp1(3), cp2(3)
    real*8 :: n1, n2, cphi, cc2, tn, dtn, dedc, vt, ct
    real*8 :: dva(3), dvb(3), dc1(3), dc2(3), dc3(3), gl(3)
    ! inversion locals
    integer :: inv, ap, ib1, ib2, iidx(3)
    real*8 :: rc(3), rnb(3,3), cosy, siny, edcy
    real*8 :: gc(3), gg1(3), gg2(3), gg3(3)

    dostress = present(stress)
    ene = 0d0
    grad = 0d0
    vir = 0d0

    ! bond stretch: E = 1/2 k (r-r_IJ)^2
    do n = 1, cl%nbond
       i = cl%bond(n)%i
       j = cl%bond(n)%j
       d = c%atcel(i)%r - atpos(c,j,cl%bond(n)%lvec)
       r = norm2(d)
       if (r < 1d-10) cycle
       dr = r - cl%bond(n)%r0
       ene = ene + 0.5d0 * cl%bond(n)%k * dr*dr
       gmag = cl%bond(n)%k * dr / r
       g = gmag * d
       grad(:,i) = grad(:,i) + g
       grad(:,j) = grad(:,j) - g
       if (dostress) call add_virial(vir,d,-g)
    end do

    ! angle bend, expressed in u = cos(theta) (no 1/sin(theta) singularity) with
    ! the Fourier coefficients precomputed at setup: E = k [C0 + C1 u + C2 (2u^2-1)]
    do n = 1, cl%nang
       ii = cl%ang(n)%i
       jj = cl%ang(n)%j
       kk = cl%ang(n)%k
       ri = atpos(c,ii,cl%ang(n)%li)
       rj = c%atcel(jj)%r
       rk = atpos(c,kk,cl%ang(n)%lk)
       a = ri - rj
       b = rk - rj
       na = norm2(a)
       nb = norm2(b)
       if (na < 1d-10 .or. nb < 1d-10) cycle
       u = dot_product(a,b)/(na*nb)
       u = max(min(u,1d0),-1d0)
       dudi = b/(na*nb) - (u/(na*na))*a
       dudk = a/(na*nb) - (u/(nb*nb))*b
       dudj = -(dudi+dudk)
       ene = ene + cl%ang(n)%kk * (cl%ang(n)%c0 + cl%ang(n)%c1*u + cl%ang(n)%c2*(2d0*u*u - 1d0))
       edu = cl%ang(n)%kk * (cl%ang(n)%c1 + 4d0*cl%ang(n)%c2*u)
       gi = edu*dudi
       gj = edu*dudj
       gk = edu*dudk
       grad(:,ii) = grad(:,ii) + gi
       grad(:,jj) = grad(:,jj) + gj
       grad(:,kk) = grad(:,kk) + gk
       if (dostress) then
          call add_virial_site(vir,ri,-gi)
          call add_virial_site(vir,rj,-gj)
          call add_virial_site(vir,rk,-gk)
       end if
    end do

    ! van der Waals (UFF 12-6): E = D [(x/r)^12 - 2 (x/r)^6]
    do n = 1, cl%nnb
       i = cl%nb(n)%i
       j = cl%nb(n)%j
       d = c%atcel(i)%r - atpos(c,j,cl%nb(n)%lvec)
       r = norm2(d)
       if (r < 1d-10 .or. r > rcut_lj) cycle
       sr = cl%nb(n)%rmin/r
       sr2 = sr*sr
       sr6 = sr2*sr2*sr2
       sr12 = sr6*sr6
       ene = ene + cl%nb(n)%deps*(sr12 - 2d0*sr6)
       ! dE/dr = 12 D (sr6 - sr12)/r ; gradient wrt d = (dE/dr) d/r
       fmag = 12d0*cl%nb(n)%deps*(sr6 - sr12)/r
       g = fmag * d / r
       grad(:,i) = grad(:,i) + g
       grad(:,j) = grad(:,j) - g
       if (dostress) call add_virial(vir,d,-g)
       ! screened Coulomb (QEq charges): E = q_i q_j / sqrt(r^2 + a_IJ^2), with the
       ! charge product and screening length^2 cached on the pair
       if (cl%do_elec) then
          rr = sqrt(r*r + cl%nb(n)%qscr2)
          rinv = 1d0/rr
          ene = ene + cl%nb(n)%qij*rinv
          gq = (-cl%nb(n)%qij*rinv*rinv*rinv) * d
          grad(:,i) = grad(:,i) + gq
          grad(:,j) = grad(:,j) - gq
          if (dostress) call add_virial(vir,d,-gq)
       end if
    end do

    ! torsions: E = V/2 (1 - cosTerm cos(n phi)), written via c = cos(phi) with
    ! the Chebyshev expansion so no sin/atan2 sign convention is needed
    do t = 1, cl%ntor
       ti = cl%tor(t)%i
       tj = cl%tor(t)%j
       tk = cl%tor(t)%k
       tl = cl%tor(t)%l
       ri = atpos(c,ti,cl%tor(t)%li)
       rj = c%atcel(tj)%r
       rk = atpos(c,tk,cl%tor(t)%lk)
       rl = atpos(c,tl,cl%tor(t)%ll)
       tv1 = rj - ri
       tv2 = rk - rj
       tv3 = rl - rk
       cp1 = cross(tv1,tv2)
       cp2 = cross(tv2,tv3)
       n1 = norm2(cp1)
       n2 = norm2(cp2)
       if (n1 < 1d-8 .or. n2 < 1d-8) cycle
       cphi = dot_product(cp1,cp2)/(n1*n2)
       cphi = max(min(cphi,1d0),-1d0)
       vt = cl%tor(t)%v
       nn = cl%tor(t)%n
       ct = cl%tor(t)%cosf
       cc2 = cphi*cphi
       select case(nn)
       case(2)
          tn = 2d0*cc2 - 1d0
          dtn = 4d0*cphi
       case(3)
          tn = cphi*(4d0*cc2 - 3d0)
          dtn = 12d0*cc2 - 3d0
       case(6)
          tn = 32d0*cc2**3 - 48d0*cc2*cc2 + 18d0*cc2 - 1d0
          dtn = cphi*(192d0*cc2*cc2 - 192d0*cc2 + 36d0)
       case default
          tn = cos(nn*acos(cphi))
          dtn = 0d0
       end select
       ene = ene + 0.5d0*vt*(1d0 - ct*tn)
       dedc = -0.5d0*vt*ct*dtn
       ! d(cos phi)/d(cp1,cp2), chained to the bond vectors and then to atoms
       dva = cp2/(n1*n2) - cphi*cp1/(n1*n1)
       dvb = cp1/(n1*n2) - cphi*cp2/(n2*n2)
       dc1 = cross(tv2,dva)
       dc2 = cross(dva,tv1) + cross(tv3,dvb)
       dc3 = cross(dvb,tv2)
       gi = dedc*(-dc1)
       gj = dedc*(dc1 - dc2)
       gk = dedc*(dc2 - dc3)
       gl = dedc*dc3
       grad(:,ti) = grad(:,ti) + gi
       grad(:,tj) = grad(:,tj) + gj
       grad(:,tk) = grad(:,tk) + gk
       grad(:,tl) = grad(:,tl) + gl
       if (dostress) then
          call add_virial_site(vir,ri,-gi)
          call add_virial_site(vir,rj,-gj)
          call add_virial_site(vir,rk,-gk)
          call add_virial_site(vir,rl,-gl)
       end if
    end do

    ! inversions (out-of-plane) at 3-coordinate centers, averaged over the 3
    ! apex choices: E = K [C0 + C1 sin(Y) + C2 cos(2 Y)]
    do inv = 1, cl%ninv
       ti = cl%inv(inv)%c
       iidx = (/cl%inv(inv)%n1,cl%inv(inv)%n2,cl%inv(inv)%n3/)
       rc = c%atcel(ti)%r
       rnb(:,1) = atpos(c,cl%inv(inv)%n1,cl%inv(inv)%l1)
       rnb(:,2) = atpos(c,cl%inv(inv)%n2,cl%inv(inv)%l2)
       rnb(:,3) = atpos(c,cl%inv(inv)%n3,cl%inv(inv)%l3)
       do ap = 1, 3
          if (ap == 1) then
             ib1 = 2; ib2 = 3
          else if (ap == 2) then
             ib1 = 1; ib2 = 3
          else
             ib1 = 1; ib2 = 2
          end if
          call uff_wilson_cosy(rc,rnb(:,ib1),rnb(:,ib2),rnb(:,ap),cosy,gc,gg1,gg2,gg3)
          siny = sqrt(max(1d0 - cosy*cosy,1d-12))
          ene = ene + cl%inv(inv)%k*(cl%inv(inv)%c0 + cl%inv(inv)%c1*siny &
             + cl%inv(inv)%c2*(2d0*siny*siny - 1d0))
          ! dE/d(cosY) = -K (C1 + 4 C2 sinY) cosY / sinY
          edcy = -cl%inv(inv)%k*(cl%inv(inv)%c1 + 4d0*cl%inv(inv)%c2*siny)*cosy/siny
          grad(:,ti) = grad(:,ti) + edcy*gc
          grad(:,iidx(ib1)) = grad(:,iidx(ib1)) + edcy*gg1
          grad(:,iidx(ib2)) = grad(:,iidx(ib2)) + edcy*gg2
          grad(:,iidx(ap))  = grad(:,iidx(ap))  + edcy*gg3
          if (dostress) then
             call add_virial_site(vir,rc,-edcy*gc)
             call add_virial_site(vir,rnb(:,ib1),-edcy*gg1)
             call add_virial_site(vir,rnb(:,ib2),-edcy*gg2)
             call add_virial_site(vir,rnb(:,ap),-edcy*gg3)
          end if
       end do
    end do

    if (dostress) then
       if (c%ismolecule .or. c%omega < 1d-10) then
          stress = 0d0
       else
          stress = -vir / c%omega
       end if
    end if

  end subroutine uff_evaluate

  !> Load the UFF parameter table from dat/uff.dat (cached after the first call).
  subroutine uff_load(errmsg)
    use global, only: critic_home
    use param, only: dirsep
    character(len=:), allocatable, intent(out), optional :: errmsg

    character(len=:), allocatable :: file
    character(len=256) :: line
    integer :: lu, ios, n, pass
    logical :: ex
    character(len=8) :: lbl
    real*8 :: r1, th, x1, dw, ze, z1, vi, uj, xi, hd, rd

    if (uff_loaded) return

    file = trim(critic_home) // dirsep // "force_fields" // dirsep // "uff.dat"
    inquire(file=file,exist=ex)
    if (.not.ex) then
       if (present(errmsg)) errmsg = "could not find the UFF parameter file: " // file
       return
    end if

    open(newunit=lu,file=file,status="old",action="read",iostat=ios)
    if (ios /= 0) then
       if (present(errmsg)) errmsg = "could not open the UFF parameter file: " // file
       return
    end if
    do pass = 1, 2
       if (pass == 2) rewind(lu)
       n = 0
       do
          read(lu,'(A)',iostat=ios) line
          if (ios /= 0) exit
          line = adjustl(line)
          if (len_trim(line) == 0) cycle
          if (line(1:1) == "#") cycle
          n = n + 1
          if (pass == 2) then
             read(line,*,iostat=ios) lbl, r1, th, x1, dw, ze, z1, vi, uj, xi, hd, rd
             if (ios /= 0) then
                close(lu)
                if (present(errmsg)) errmsg = "error parsing the UFF parameter file: " // file
                return
             end if
             uffpar(n)%label = lbl
             uffpar(n)%r1 = r1
             uffpar(n)%theta0 = th
             uffpar(n)%x1 = x1
             uffpar(n)%dwell = dw
             uffpar(n)%zstar = z1
             uffpar(n)%vtor = vi
             uffpar(n)%utor = uj
             uffpar(n)%chi = xi
             uffpar(n)%jhard = hd
          end if
       end do
       if (pass == 1) then
          nuffpar = n
          if (allocated(uffpar)) deallocate(uffpar)
          allocate(uffpar(n))
       end if
    end do
    close(lu)
    uff_loaded = .true.

  end subroutine uff_load

  !> Index of the UFF parameter row with the given atom-type label (0 if absent).
  integer function uff_find(lbl)
    character(len=*), intent(in) :: lbl
    integer :: i
    uff_find = 0
    do i = 1, nuffpar
       if (trim(uffpar(i)%label) == trim(lbl)) then
          uff_find = i
          return
       end if
    end do
  end function uff_find

  !> Assign a UFF atom type to atom iat from its element and connectivity
  !> (coordination number, bond orders, aromaticity). Returns the parameter-table
  !> index, or 0 if the element has no UFF parameters.
  function uff_type_index(c,iat,ns) result(ip)
    use crystalmod, only: crystal
    use types, only: neighstar
    class(crystal), intent(in) :: c
    integer, intent(in) :: iat
    type(neighstar), intent(in) :: ns(:)
    integer :: ip

    integer :: z, ncon
    logical :: arom
    character(len=8) :: lbl

    z = c%spc(c%atcel(iat)%is)%z
    ncon = ns(iat)%ncon
    arom = ns(iat)%isaromatic

    select case(z)
    case(1) ! H
       lbl = "H_"
    case(5) ! B
       if (ncon >= 4) then; lbl = "B_3"; else; lbl = "B_2"; end if
    case(6) ! C
       if (arom) then; lbl = "C_R"; elseif (ncon >= 4) then; lbl = "C_3"
       elseif (ncon == 3) then; lbl = "C_2"; else; lbl = "C_1"; end if
    case(7) ! N
       if (arom) then; lbl = "N_R"; elseif (ncon >= 3) then; lbl = "N_3"
       elseif (ncon == 2) then; lbl = "N_2"; else; lbl = "N_1"; end if
    case(8) ! O
       if (arom) then; lbl = "O_R"; elseif (ncon >= 2) then; lbl = "O_3"; else; lbl = "O_2"; end if
    case(14) ! Si
       lbl = "Si3"
    case(15) ! P
       if (ncon >= 4) then; lbl = "P_3+5"; else; lbl = "P_3+3"; end if
    case(16) ! S
       if (arom) then; lbl = "S_R"; elseif (ncon >= 4) then; lbl = "S_3+6"
       elseif (ncon == 3) then; lbl = "S_3+4"; else; lbl = "S_3+2"; end if
    case(9) ! F
       lbl = "F_"
    case(17) ! Cl
       lbl = "Cl"
    case(35) ! Br
       lbl = "Br"
    case(53) ! I
       lbl = "I_"
    case default
       lbl = uff_default_label(z)
    end select

    ip = uff_find(lbl)

  end function uff_type_index

  !> Canonical (default) UFF atom-type label for an element, used for atoms whose
  !> type does not depend on hybridization (noble gases, metals, ...). Returns an
  !> empty string for elements outside the UFF table.
  function uff_default_label(z) result(lbl)
    integer, intent(in) :: z
    character(len=8) :: lbl
    select case(z)
    case(1);   lbl = "H_"
    case(2);   lbl = "He4+4"
    case(3);   lbl = "Li"
    case(4);   lbl = "Be3+2"
    case(5);   lbl = "B_3"
    case(6);   lbl = "C_3"
    case(7);   lbl = "N_3"
    case(8);   lbl = "O_3"
    case(9);   lbl = "F_"
    case(10);  lbl = "Ne4+4"
    case(11);  lbl = "Na"
    case(12);  lbl = "Mg3+2"
    case(13);  lbl = "Al3"
    case(14);  lbl = "Si3"
    case(15);  lbl = "P_3+3"
    case(16);  lbl = "S_3+2"
    case(17);  lbl = "Cl"
    case(18);  lbl = "Ar4+4"
    case(19);  lbl = "K_"
    case(20);  lbl = "Ca6+2"
    case(21);  lbl = "Sc3+3"
    case(22);  lbl = "Ti6+4"
    case(23);  lbl = "V_3+5"
    case(24);  lbl = "Cr6+3"
    case(25);  lbl = "Mn6+2"
    case(26);  lbl = "Fe6+2"
    case(27);  lbl = "Co6+3"
    case(28);  lbl = "Ni4+2"
    case(29);  lbl = "Cu3+1"
    case(30);  lbl = "Zn3+2"
    case(31);  lbl = "Ga3+3"
    case(32);  lbl = "Ge3"
    case(33);  lbl = "As3+3"
    case(34);  lbl = "Se3+2"
    case(35);  lbl = "Br"
    case(36);  lbl = "Kr4+4"
    case(37);  lbl = "Rb"
    case(38);  lbl = "Sr6+2"
    case(39);  lbl = "Y_3+3"
    case(40);  lbl = "Zr3+4"
    case(41);  lbl = "Nb3+5"
    case(42);  lbl = "Mo6+6"
    case(43);  lbl = "Tc6+5"
    case(44);  lbl = "Ru6+2"
    case(45);  lbl = "Rh6+3"
    case(46);  lbl = "Pd4+2"
    case(47);  lbl = "Ag1+1"
    case(48);  lbl = "Cd3+2"
    case(49);  lbl = "In3+3"
    case(50);  lbl = "Sn3"
    case(51);  lbl = "Sb3+3"
    case(52);  lbl = "Te3+2"
    case(53);  lbl = "I_"
    case(54);  lbl = "Xe4+4"
    case(55);  lbl = "Cs"
    case(56);  lbl = "Ba6+2"
    case(57);  lbl = "La3+3"
    case(58);  lbl = "Ce6+3"
    case(59);  lbl = "Pr6+3"
    case(60);  lbl = "Nd6+3"
    case(61);  lbl = "Pm6+3"
    case(62);  lbl = "Sm6+3"
    case(63);  lbl = "Eu6+3"
    case(64);  lbl = "Gd6+3"
    case(65);  lbl = "Tb6+3"
    case(66);  lbl = "Dy6+3"
    case(67);  lbl = "Ho6+3"
    case(68);  lbl = "Er6+3"
    case(69);  lbl = "Tm6+3"
    case(70);  lbl = "Yb6+3"
    case(71);  lbl = "Lu6+3"
    case(72);  lbl = "Hf3+4"
    case(73);  lbl = "Ta3+5"
    case(74);  lbl = "W_6+6"
    case(75);  lbl = "Re6+5"
    case(76);  lbl = "Os6+6"
    case(77);  lbl = "Ir6+3"
    case(78);  lbl = "Pt4+2"
    case(79);  lbl = "Au4+3"
    case(80);  lbl = "Hg1+2"
    case(81);  lbl = "Tl3+3"
    case(82);  lbl = "Pb3"
    case(83);  lbl = "Bi3+3"
    case(84);  lbl = "Po3+2"
    case(85);  lbl = "At"
    case(86);  lbl = "Rn4+4"
    case(87);  lbl = "Fr"
    case(88);  lbl = "Ra6+2"
    case(89);  lbl = "Ac6+3"
    case(90);  lbl = "Th6+4"
    case(91);  lbl = "Pa6+4"
    case(92);  lbl = "U_6+4"
    case(93);  lbl = "Np6+4"
    case(94);  lbl = "Pu6+4"
    case(95);  lbl = "Am6+4"
    case(96);  lbl = "Cm6+3"
    case(97);  lbl = "Bk6+3"
    case(98);  lbl = "Cf6+3"
    case(99);  lbl = "Es6+3"
    case(100); lbl = "Fm6+3"
    case(101); lbl = "Md6+3"
    case(102); lbl = "No6+3"
    case(103); lbl = "Lw6+3"
    case default; lbl = ""
    end select
  end function uff_default_label

  !> Real bond order of the connection k of a neighbor star (aromatic -> 1.5,
  !> unset/dashed -> 1).
  pure function bond_order(nsa,k) result(bo)
    use types, only: neighstar
    type(neighstar), intent(in) :: nsa
    integer, intent(in) :: k
    real*8 :: bo
    integer :: o
    o = 1
    if (allocated(nsa%ordcon)) o = nsa%ordcon(k)
    if (o == -1) then
       bo = 1.5d0
    else if (o <= 0) then
       bo = 1d0
    else
       bo = real(o,8)
    end if
  end function bond_order

  !> UFF natural bond length r_IJ (angstrom) for atom types iti/itj and bond
  !> order bo, including bond-order and electronegativity corrections.
  function uff_natural_length(iti,itj,bo) result(rij)
    integer, intent(in) :: iti, itj
    real*8, intent(in) :: bo
    real*8 :: rij, ri, rj, ci, cj, rbo, ren

    ri = uffpar(iti)%r1
    rj = uffpar(itj)%r1
    ci = uffpar(iti)%chi
    cj = uffpar(itj)%chi
    rbo = -uff_lambda*(ri+rj)*log(bo)
    ren = 0d0
    if (ci*ri + cj*rj > 1d-10) &
       ren = ri*rj*(sqrt(ci)-sqrt(cj))**2 / (ci*ri + cj*rj)
    rij = ri + rj + rbo - ren

  end function uff_natural_length

  !> Ohno screening length a_IJ (bohr) for the QEq electrostatics between atom
  !> types iti and itj: a_IJ = 2/(J_ii+J_jj) with idempotential J = 2*hardness,
  !> i.e. hartoev/(hardness_i+hardness_j). Shared by the QEq matrix and the
  !> electrostatic energy so the two cannot drift.
  pure function qeq_ascreen(iti,itj) result(a)
    integer, intent(in) :: iti, itj
    real*8 :: a
    a = hartoev / max(uffpar(iti)%jhard + uffpar(itj)%jhard,1d-10)
  end function qeq_ascreen

  !> Hybridization class of a UFF atom type from its geometry code: 1 = linear
  !> (sp), 2 = trigonal/resonant (sp2, includes R), 3 = tetrahedral (sp3);
  !> 0 = other (linear terminal, ionic, square-planar/octahedral metal, ...).
  pure function uff_hyb(ip) result(h)
    integer, intent(in) :: ip
    integer :: h
    integer :: i
    character :: ch
    h = 0
    do i = 1, len_trim(uffpar(ip)%label)
       ch = uffpar(ip)%label(i:i)
       select case(ch)
       case("1"); h = 1; return
       case("2","R"); h = 2; return
       case("3"); h = 3; return
       case("4","5","6"); h = 0; return
       end select
    end do
  end function uff_hyb

  !> True if z is a group-16 (oxygen-column) element.
  pure function isgroup6(z) result(ok)
    integer, intent(in) :: z
    logical :: ok
    ok = (z==8 .or. z==16 .or. z==34 .or. z==52 .or. z==84)
  end function isgroup6

  !> UFF group-16 sp3 torsional barrier (kcal/mol): 2.0 for O, 6.8 otherwise.
  pure function group6_v(z) result(v)
    integer, intent(in) :: z
    real*8 :: v
    if (z == 8) then
       v = 2d0
    else
       v = 6.8d0
    end if
  end function group6_v

  !> Build the UFF torsion list (dihedrals about each central bond j-k whose
  !> atoms are both sp2 or sp3), assigning the barrier V, periodicity n, and
  !> phase from the two central hybridizations and the bond order.
  subroutine uff_build_torsions(cl,c,ns)
    use crystalmod, only: crystal
    use types, only: neighstar
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    type(neighstar), intent(in) :: ns(:)

    integer :: jj, kk, a, b, ni, nl, ka, hj, hk, nt, pass, norder
    integer :: ljk(3), li(3), lkl(3)
    real*8 :: v, costerm
    logical :: ok

    do pass = 1, 2
       nt = 0
       do jj = 1, cl%nat
          hj = uff_hyb(cl%attype(jj))
          if (hj /= 2 .and. hj /= 3) cycle
          do ka = 1, ns(jj)%ncon
             kk = ns(jj)%idcon(ka)
             ljk = ns(jj)%lcon(:,ka)
             if (.not.bond_canonical(jj,kk,ljk)) cycle ! each central bond once
             hk = uff_hyb(cl%attype(kk))
             if (hk /= 2 .and. hk /= 3) cycle
             call uff_torsion_params(cl,c,ns,jj,kk,bond_order(ns(jj),ka),v,norder,costerm,ok)
             if (.not.ok) cycle
             ! UFF distributes the barrier over all torsions about the central bond
             v = v / real(max((ns(jj)%ncon-1)*(ns(kk)%ncon-1),1),8)
             do a = 1, ns(jj)%ncon
                ni = ns(jj)%idcon(a)
                li = ns(jj)%lcon(:,a)
                if (ni == kk .and. all(li == ljk)) cycle ! skip the central bond
                do b = 1, ns(kk)%ncon
                   nl = ns(kk)%idcon(b)
                   lkl = ns(kk)%lcon(:,b)
                   ! L position relative to jj's cell = ljk + lkl; skip if it is jj itself
                   if (nl == jj .and. all(ljk+lkl == 0)) cycle
                   nt = nt + 1
                   if (pass == 2) then
                      cl%tor(nt)%i = ni
                      cl%tor(nt)%j = jj
                      cl%tor(nt)%k = kk
                      cl%tor(nt)%l = nl
                      cl%tor(nt)%li = li
                      cl%tor(nt)%lk = ljk
                      cl%tor(nt)%ll = ljk + lkl
                      cl%tor(nt)%v = v * kcal2ha
                      cl%tor(nt)%n = norder
                      cl%tor(nt)%cosf = costerm
                   end if
                end do
             end do
          end do
       end do
       if (pass == 1) then
          cl%ntor = nt
          allocate(cl%tor(nt))
       end if
    end do

  end subroutine uff_build_torsions

  !> UFF torsion parameters for the central bond jj-kk (both sp2/sp3): barrier v
  !> (kcal/mol), periodicity norder, and phase term costerm (+-1). ok=.false. if
  !> no torsion is defined.
  subroutine uff_torsion_params(cl,c,ns,jj,kk,bojk,v,norder,costerm,ok)
    use crystalmod, only: crystal
    use types, only: neighstar
    class(calculator), intent(in) :: cl
    class(crystal), intent(in) :: c
    type(neighstar), intent(in) :: ns(:)
    integer, intent(in) :: jj, kk
    real*8, intent(in) :: bojk
    real*8, intent(out) :: v, costerm
    integer, intent(out) :: norder
    logical, intent(out) :: ok

    integer :: hj, hk, zj, zk, is3, is2, z3, z2, m
    logical :: sp2sp2

    v = 0d0; norder = 1; costerm = 1d0; ok = .true.
    hj = uff_hyb(cl%attype(jj))
    hk = uff_hyb(cl%attype(kk))
    zj = c%spc(c%atcel(jj)%is)%z
    zk = c%spc(c%atcel(kk)%is)%z

    if (hj == 3 .and. hk == 3) then
       ! both sp3
       if (abs(bojk-1d0) < 1d-3 .and. isgroup6(zj) .and. isgroup6(zk)) then
          v = sqrt(group6_v(zj)*group6_v(zk)); norder = 2; costerm = -1d0
       else
          v = sqrt(uffpar(cl%attype(jj))%vtor * uffpar(cl%attype(kk))%vtor)
          norder = 3; costerm = -1d0
       end if
    else if (hj == 2 .and. hk == 2) then
       ! both sp2: bond-order dependent barrier (Rappe eq. 17)
       v = 5d0*sqrt(uffpar(cl%attype(jj))%utor * uffpar(cl%attype(kk))%utor) &
          *(1d0 + 4.18d0*log(bojk))
       norder = 2; costerm = 1d0
    else
       ! sp2-sp3 mixed
       if (hj == 3) then; is3 = jj; is2 = kk; else; is3 = kk; is2 = jj; end if
       z3 = c%spc(c%atcel(is3)%is)%z
       z2 = c%spc(c%atcel(is2)%is)%z
       if (isgroup6(z3) .and. .not.isgroup6(z2)) then
          v = 5d0*sqrt(uffpar(cl%attype(is3))%utor * uffpar(cl%attype(is2))%utor) &
             *(1d0 + 4.18d0*log(bojk))
          norder = 2; costerm = -1d0
       else
          ! sp3-sp2-sp2 conjugation: the sp2 center bonded to another sp2
          sp2sp2 = .false.
          do m = 1, ns(is2)%ncon
             if (ns(is2)%idcon(m) == is3) cycle
             if (uff_hyb(cl%attype(ns(is2)%idcon(m))) == 2) sp2sp2 = .true.
          end do
          if (sp2sp2) then
             v = 2d0; norder = 3; costerm = -1d0
          else
             v = 1d0; norder = 6; costerm = 1d0
          end if
       end if
    end if
    ! a zero barrier contributes nothing; drop it
    if (abs(v) < 1d-12) ok = .false.

  end subroutine uff_torsion_params

  !> Build the UFF inversion (out-of-plane) list for 3-coordinate C (sp2), N
  !> (sp2), and group-15 (P/As/Sb/Bi, pyramidal) centers.
  subroutine uff_build_inversions(cl,c,ns)
    use crystalmod, only: crystal
    use types, only: neighstar
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    type(neighstar), intent(in) :: ns(:)

    integer :: ic, ni, pass
    real*8 :: kk, c0, c1, c2
    logical :: ok

    do pass = 1, 2
       ni = 0
       do ic = 1, cl%nat
          if (ns(ic)%ncon /= 3) cycle
          call uff_inversion_params(cl,c,ns,ic,kk,c0,c1,c2,ok)
          if (.not.ok) cycle
          ni = ni + 1
          if (pass == 2) then
             cl%inv(ni)%c = ic
             cl%inv(ni)%n1 = ns(ic)%idcon(1)
             cl%inv(ni)%n2 = ns(ic)%idcon(2)
             cl%inv(ni)%n3 = ns(ic)%idcon(3)
             cl%inv(ni)%l1 = ns(ic)%lcon(:,1)
             cl%inv(ni)%l2 = ns(ic)%lcon(:,2)
             cl%inv(ni)%l3 = ns(ic)%lcon(:,3)
             cl%inv(ni)%k = kk * kcal2ha / 3d0 ! averaged over the 3 apex choices
             cl%inv(ni)%c0 = c0
             cl%inv(ni)%c1 = c1
             cl%inv(ni)%c2 = c2
          end if
       end do
       if (pass == 1) then
          cl%ninv = ni
          allocate(cl%inv(ni))
       end if
    end do

  end subroutine uff_build_inversions

  !> UFF inversion coefficients and force constant (kcal/mol) for a 3-coordinate
  !> center ic. ok=.false. if the element has no UFF inversion term.
  subroutine uff_inversion_params(cl,c,ns,ic,kk,c0,c1,c2,ok)
    use crystalmod, only: crystal
    use types, only: neighstar
    use param, only: rad
    class(calculator), intent(in) :: cl
    class(crystal), intent(in) :: c
    type(neighstar), intent(in) :: ns(:)
    integer, intent(in) :: ic
    real*8, intent(out) :: kk, c0, c1, c2
    logical, intent(out) :: ok

    integer :: z, h, m
    real*8 :: w0, cw
    logical :: boundo

    ok = .true.; kk = 0d0; c0 = 0d0; c1 = 0d0; c2 = 0d0
    z = c%spc(c%atcel(ic)%is)%z
    h = uff_hyb(cl%attype(ic))

    select case(z)
    case(6) ! C: planar sp2 (50 kcal/mol if bonded to O, else 6)
       if (h /= 2) then; ok = .false.; return; end if
       c0 = 1d0; c1 = -1d0; c2 = 0d0
       boundo = .false.
       do m = 1, ns(ic)%ncon
          if (c%spc(c%atcel(ns(ic)%idcon(m))%is)%z == 8) boundo = .true.
       end do
       if (boundo) then; kk = 50d0; else; kk = 6d0; end if
    case(7) ! N: planar sp2 (amide/aromatic)
       if (h /= 2) then; ok = .false.; return; end if
       c0 = 1d0; c1 = -1d0; c2 = 0d0; kk = 6d0
    case(15,33,51,83) ! P, As, Sb, Bi: pyramidal
       select case(z)
       case(15); w0 = 84.4339d0
       case(33); w0 = 86.9735d0
       case(51); w0 = 87.7047d0
       case default; w0 = 90d0
       end select
       w0 = w0 * rad
       cw = cos(w0)
       c2 = 1d0
       c1 = -4d0*cw
       c0 = -(c1*cw + c2*(2d0*cw*cw - 1d0))
       kk = 22d0
    case default
       ok = .false.
    end select

  end subroutine uff_inversion_params

  !> Wilson out-of-plane cosine for the inversion term. rc is the central atom,
  !> rj/rk the two base atoms, rl the apex. Returns cosy = uL . n where n is the
  !> unit normal of the (rj,rk) plane through rc, and the gradients of cosy with
  !> respect to the four atom positions.
  pure subroutine uff_wilson_cosy(rc,rj,rk,rl,cosy,gc,gj,gk,gl)
    use tools_math, only: cross
    real*8, intent(in) :: rc(3), rj(3), rk(3), rl(3)
    real*8, intent(out) :: cosy, gc(3), gj(3), gk(3), gl(3)

    real*8 :: ej(3), ek(3), el(3), cp(3), ncp, nl, trip
    real*8 :: ej2, ek2, ejk

    ej = rj - rc
    ek = rk - rc
    el = rl - rc
    cp = cross(ej,ek)
    ncp = norm2(cp)
    nl = norm2(el)
    if (ncp < 1d-10 .or. nl < 1d-10) then
       cosy = 0d0; gc = 0d0; gj = 0d0; gk = 0d0; gl = 0d0
       return
    end if
    trip = dot_product(el,cp)
    cosy = trip/(nl*ncp)
    ej2 = dot_product(ej,ej)
    ek2 = dot_product(ek,ek)
    ejk = dot_product(ej,ek)
    ! d(cosy)/d(el), d(ej), d(ek)  (see module notes)
    gl = cp/(nl*ncp) - trip*el/(nl**3*ncp)
    gj = cross(ek,el)/(nl*ncp) - trip*(ek2*ej - ejk*ek)/(nl*ncp**3)
    gk = cross(el,ej)/(nl*ncp) - trip*(ej2*ek - ejk*ej)/(nl*ncp**3)
    ! central atom: e* = r* - rc, so d/drc = -(sum of the three)
    gc = -(gj + gk + gl)

  end subroutine uff_wilson_cosy

  !> Cartesian position (bohr) of cell atom ia translated by integer lattice
  !> vector lv.
  pure function atpos(c,ia,lv) result(r)
    use crystalmod, only: crystal
    class(crystal), intent(in) :: c
    integer, intent(in) :: ia, lv(3)
    real*8 :: r(3)
    r = c%atcel(ia)%r + c%x2c(real(lv,8))
  end function atpos

  !> True if the directed bond (i -> j with lattice vector lv) is the canonical
  !> representative of its undirected pair (avoids double-counting).
  pure function bond_canonical(i,j,lv) result(ok)
    integer, intent(in) :: i, j, lv(3)
    logical :: ok
    if (i < j) then
       ok = .true.
    else if (i > j) then
       ok = .false.
    else
       ok = lv_positive(lv)
    end if
  end function bond_canonical

  !> Lexicographic positivity test for a lattice vector (used to pick one of
  !> lv / -lv when i == j).
  pure function lv_positive(lv) result(ok)
    integer, intent(in) :: lv(3)
    logical :: ok
    if (lv(1) /= 0) then
       ok = lv(1) > 0
    else if (lv(2) /= 0) then
       ok = lv(2) > 0
    else
       ok = lv(3) > 0
    end if
  end function lv_positive

  !> True if atoms (i, j+lv) are 1-2 (directly bonded) or 1-3 (both bonded to a
  !> common atom), and should therefore be excluded from the nonbonded term.
  function excluded_13(ns,i,j,lv) result(ex)
    use types, only: neighstar
    type(neighstar), intent(in) :: ns(:)
    integer, intent(in) :: i, j, lv(3)
    logical :: ex

    integer :: ka, kb, a, la(3), b, lab(3)

    ex = .false.
    ! 1-2: i directly bonded to (j,lv)
    do ka = 1, ns(i)%ncon
       if (ns(i)%idcon(ka) == j .and. all(ns(i)%lcon(:,ka) == lv)) then
          ex = .true.
          return
       end if
    end do
    ! 1-3: some neighbor a of i is bonded to (j,lv)
    do ka = 1, ns(i)%ncon
       a = ns(i)%idcon(ka)
       la = ns(i)%lcon(:,ka)
       do kb = 1, ns(a)%ncon
          b = ns(a)%idcon(kb)
          lab = la + ns(a)%lcon(:,kb)
          if (b == j .and. all(lab == lv)) then
             ex = .true.
             return
          end if
       end do
    end do

  end function excluded_13

  !> Accumulate the virial contribution of a pair: separation vector d = ri - rj
  !> and the force fi on atom i.
  pure subroutine add_virial(vir,d,fi)
    real*8, intent(inout) :: vir(3,3)
    real*8, intent(in) :: d(3), fi(3)
    integer :: p, q
    do q = 1, 3
       do p = 1, 3
          vir(p,q) = vir(p,q) + d(p)*fi(q)
       end do
    end do
  end subroutine add_virial

  !> Accumulate the virial contribution of a single (imaged) site at position r
  !> feeling force f. Summing over all sites of a cluster is translationally
  !> invariant because the cluster forces sum to zero.
  pure subroutine add_virial_site(vir,r,f)
    real*8, intent(inout) :: vir(3,3)
    real*8, intent(in) :: r(3), f(3)
    integer :: p, q
    do q = 1, 3
       do p = 1, 3
          vir(p,q) = vir(p,q) + r(p)*f(q)
       end do
    end do
  end subroutine add_virial_site

  ! ---- tblite (GFN-FF / xTB) backend ----

  module subroutine calc_init_tblite(cl,c,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out), optional :: errmsg

#ifdef HAVE_TBLITE
    type(c_ptr) :: err
    integer(c_int), allocatable, target :: numbers(:)
    real(c_double), allocatable, target :: coords(:,:)
    real(c_double), target :: lattice(3,3)
    logical(c_bool), target :: periodic(3)
    integer :: i

    ! atomic numbers and positions (bohr)
    allocate(numbers(c%ncel),coords(3,c%ncel))
    do i = 1, c%ncel
       numbers(i) = int(c%spc(c%atcel(i)%is)%z,c_int)
       coords(:,i) = c%atcel(i)%r
    end do

    err = c_tblite_new_error()
    if (c%ismolecule) then
       cl%tb_mol = c_tblite_new_structure(err,int(c%ncel,c_int),c_loc(numbers),c_loc(coords),&
          c_null_ptr,c_null_ptr,c_null_ptr,c_null_ptr)
    else
       lattice = c%m_x2c
       periodic = .true._c_bool
       cl%tb_mol = c_tblite_new_structure(err,int(c%ncel,c_int),c_loc(numbers),c_loc(coords),&
          c_null_ptr,c_null_ptr,c_loc(lattice),c_loc(periodic))
    end if
    cl%tb_ctx = c_tblite_new_context()

    ! GFN1/GFN2-xTB are topology-free (need only elements and coordinates)
    if (cl%method == tbm_gfn1) then
       cl%tb_calc = c_tblite_new_gfn1_calculator(cl%tb_ctx,cl%tb_mol,c_null_ptr)
    else
       cl%tb_calc = c_tblite_new_gfn2_calculator(cl%tb_ctx,cl%tb_mol,c_null_ptr)
    end if
    cl%tb_res = c_tblite_new_result()
    call c_tblite_delete_error(err)

    if (.not.(c_associated(cl%tb_mol).and.c_associated(cl%tb_calc).and.&
       c_associated(cl%tb_ctx).and.c_associated(cl%tb_res))) then
       if (present(errmsg)) errmsg = "tblite calculator initialization failed"
    end if
#else
    if (present(errmsg)) errmsg = "critic2 was compiled without tblite support"
#endif

  end subroutine calc_init_tblite

  module subroutine calc_eval_tblite(cl,c,ene,grad,stress,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    real*8, intent(out) :: ene
    real*8, intent(out) :: grad(:,:)
    real*8, intent(out), optional :: stress(3,3)
    character(len=:), allocatable, intent(out), optional :: errmsg

#ifdef HAVE_TBLITE
    type(c_ptr) :: err
    real(c_double), allocatable, target :: coords(:,:), grad_c(:)
    real(c_double), target :: lattice(3,3)
    real(c_double) :: e_c, sigma(9)
    integer :: i

    ene = 0d0
    grad = 0d0
    if (present(stress)) stress = 0d0

    allocate(coords(3,c%ncel),grad_c(3*c%ncel))
    do i = 1, c%ncel
       coords(:,i) = c%atcel(i)%r
    end do

    err = c_tblite_new_error()
    ! push the current geometry into the tblite structure
    if (c%ismolecule) then
       call c_tblite_update_structure_geometry(err,cl%tb_mol,c_loc(coords),c_null_ptr)
    else
       lattice = c%m_x2c
       call c_tblite_update_structure_geometry(err,cl%tb_mol,c_loc(coords),c_loc(lattice))
    end if

    ! single point
    call c_tblite_get_singlepoint(cl%tb_ctx,cl%tb_mol,cl%tb_calc,cl%tb_res)
    call c_tblite_get_result_energy(err,cl%tb_res,e_c)
    ene = e_c
    call c_tblite_get_result_gradient(err,cl%tb_res,grad_c)
    grad = reshape(grad_c,(/3,c%ncel/))
    if (present(stress) .and. .not.c%ismolecule .and. c%omega > 1d-10) then
       call c_tblite_get_result_virial(err,cl%tb_res,sigma)
       ! tblite returns the virial (dE/d(strain)); convert to stress
       stress = reshape(sigma,(/3,3/)) / c%omega
    end if
    call c_tblite_delete_error(err)
#else
    ene = 0d0
    grad = 0d0
    if (present(stress)) stress = 0d0
    if (present(errmsg)) errmsg = "critic2 was compiled without tblite support"
#endif

  end subroutine calc_eval_tblite

  module subroutine calc_update_tblite(cl,c)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    ! tblite reads the geometry on every evaluate (via update_structure_geometry),
    ! so nothing needs to be cached here
  end subroutine calc_update_tblite

  module subroutine calc_free_tblite(cl)
    class(calculator), intent(inout) :: cl
#ifdef HAVE_TBLITE
    if (c_associated(cl%tb_res)) call c_tblite_delete_result(cl%tb_res)
    if (c_associated(cl%tb_calc)) call c_tblite_delete_calculator(cl%tb_calc)
    if (c_associated(cl%tb_mol)) call c_tblite_delete_structure(cl%tb_mol)
    if (c_associated(cl%tb_ctx)) call c_tblite_delete_context(cl%tb_ctx)
    cl%tb_res = c_null_ptr
    cl%tb_calc = c_null_ptr
    cl%tb_mol = c_null_ptr
    cl%tb_ctx = c_null_ptr
#endif
  end subroutine calc_free_tblite

end submodule proc
