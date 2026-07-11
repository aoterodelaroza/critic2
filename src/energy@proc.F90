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
!> classical force field (harmonic bonds and angles from the connectivity graph
!> plus a Lennard-Jones nonbonded term) and the tblite (GFN-FF / xTB) backend.
!> The tblite routines are compiled with real functionality only when built with
!> -DHAVE_TBLITE; otherwise they are stubs and the calculator falls back to the
!> built-in force field.
submodule (energy) proc
  use iso_c_binding
  implicit none

  ! default built-in force-field parameters (atomic units)
  real*8, parameter :: kbond_def = 0.30d0 !< bond force constant, hartree/bohr^2 (per unit bond order)
  real*8, parameter :: kang_def = 0.10d0 !< angle force constant, hartree/rad^2
  real*8, parameter :: eps_lj = 2d-4 !< Lennard-Jones well depth, hartree
  real*8, parameter :: rcut_lj = 12d0 !< Lennard-Jones cutoff radius, bohr
  real*8, parameter :: sin_min = 1d-4 !< regularization for angles near 0/pi

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
  !> or ff_tblite; method selects the tblite method. On error, errmsg is allocated.
  module subroutine calc_init(cl,c,backend,method,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    integer, intent(in), optional :: backend
    integer, intent(in), optional :: method
    character(len=:), allocatable, intent(out), optional :: errmsg

    call cl%free()
    if (present(backend)) cl%backend = backend
    if (present(method)) cl%method = method
    cl%nat = c%ncel

    if (cl%backend == ff_builtin) then
       call builtin_setup(cl,c)
       cl%ready = .true.
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
       call builtin_evaluate(cl,c,ene,grad,stress)
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
       call builtin_nonbonded(cl,c)
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
    ! each group of arrays is allocated together, so a single guard suffices
    if (allocated(cl%bond_i)) deallocate(cl%bond_i,cl%bond_j,cl%bond_lvec,cl%bond_r0,cl%bond_k)
    if (allocated(cl%ang_i)) deallocate(cl%ang_i,cl%ang_j,cl%ang_k,cl%ang_li,cl%ang_lk,cl%ang_t0,cl%ang_kk)
    if (allocated(cl%nb_i)) deallocate(cl%nb_i,cl%nb_j,cl%nb_lvec,cl%nb_sig,cl%nb_eps)

  end subroutine calc_free

  !xx! private procedures

  !> Build the built-in force-field topology (bonds, angles, nonbonded pairs)
  !> from the connectivity graph and the current (equilibrium) geometry of c.
  subroutine builtin_setup(cl,c)
    use crystalmod, only: crystal
    use types, only: neighstar
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c

    type(neighstar), allocatable :: ns(:)
    integer :: i, j, k, ka, kb, ni, nk
    integer :: lv(3), li(3), lk(3)
    integer :: nb, na
    real*8 :: d(3), r0, a(3), b(3), t0
    integer :: ord

    ! connectivity graph (use the crystal's own if available)
    if (allocated(c%nstar)) then
       ns = c%nstar
    else
       call c%find_asterisms(ns)
    end if

    ! ---- bonds ----
    ! first pass: count (each undirected bond is stored once, canonically)
    nb = 0
    do i = 1, cl%nat
       do k = 1, ns(i)%ncon
          j = ns(i)%idcon(k)
          lv = ns(i)%lcon(:,k)
          if (.not.bond_canonical(i,j,lv)) cycle
          nb = nb + 1
       end do
    end do
    cl%nbond = nb
    allocate(cl%bond_i(nb),cl%bond_j(nb),cl%bond_lvec(3,nb),cl%bond_r0(nb),cl%bond_k(nb))

    ! second pass: fill
    nb = 0
    do i = 1, cl%nat
       do k = 1, ns(i)%ncon
          j = ns(i)%idcon(k)
          lv = ns(i)%lcon(:,k)
          if (.not.bond_canonical(i,j,lv)) cycle
          nb = nb + 1
          cl%bond_i(nb) = i
          cl%bond_j(nb) = j
          cl%bond_lvec(:,nb) = lv
          d = atpos(c,j,lv) - c%atcel(i)%r
          r0 = norm2(d)
          cl%bond_r0(nb) = r0
          ord = 1
          if (allocated(ns(i)%ordcon)) ord = max(ns(i)%ordcon(k),1)
          cl%bond_k(nb) = kbond_def * real(ord,8)
       end do
    end do

    ! ---- angles (i-j-k, vertex j) ----
    ! count: every unordered pair of neighbors of each vertex is an angle
    na = 0
    do j = 1, cl%nat
       na = na + ns(j)%ncon*(ns(j)%ncon-1)/2
    end do
    cl%nang = na
    allocate(cl%ang_i(na),cl%ang_j(na),cl%ang_k(na),cl%ang_li(3,na),cl%ang_lk(3,na),&
       cl%ang_t0(na),cl%ang_kk(na))

    ! second pass: fill
    na = 0
    do j = 1, cl%nat
       do ka = 1, ns(j)%ncon
          ni = ns(j)%idcon(ka)
          li = ns(j)%lcon(:,ka)
          do kb = ka+1, ns(j)%ncon
             nk = ns(j)%idcon(kb)
             lk = ns(j)%lcon(:,kb)
             na = na + 1
             cl%ang_i(na) = ni
             cl%ang_j(na) = j
             cl%ang_k(na) = nk
             cl%ang_li(:,na) = li
             cl%ang_lk(:,na) = lk
             a = atpos(c,ni,li) - c%atcel(j)%r
             b = atpos(c,nk,lk) - c%atcel(j)%r
             t0 = angle_between(a,b)
             cl%ang_t0(na) = t0
             cl%ang_kk(na) = kang_def
          end do
       end do
    end do

    ! ---- nonbonded pairs ----
    call builtin_nonbonded(cl,c)

  end subroutine builtin_setup

  !> Build the Lennard-Jones nonbonded pair list within rcut_lj, excluding 1-2
  !> (bonded) and 1-3 (share a common neighbor) pairs.
  subroutine builtin_nonbonded(cl,c)
    use crystalmod, only: crystal
    use types, only: neighstar
    use param, only: icrd_crys, atmvdw
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c

    type(neighstar), allocatable :: ns(:)
    integer :: i, j, m, nat, nnb, pass
    integer, allocatable :: eid(:), lvec(:,:)
    real*8, allocatable :: dist(:)
    integer :: lv(3)
    real*8 :: sig
    real*8, parameter :: onesixth_inv = 2d0**(1d0/6d0)

    if (allocated(c%nstar)) then
       ns = c%nstar
    else
       call c%find_asterisms(ns)
    end if

    if (allocated(cl%nb_i)) deallocate(cl%nb_i,cl%nb_j,cl%nb_lvec,cl%nb_sig,cl%nb_eps)

    do pass = 1, 2
       nnb = 0
       do i = 1, cl%nat
          call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.false.,nat,eid=eid,dist=dist,&
             lvec=lvec,up2d=rcut_lj,nozero=.true.)
          do m = 1, nat
             j = eid(m)
             lv = lvec(:,m)
             if (.not.bond_canonical(i,j,lv)) cycle
             if (excluded_13(ns,i,j,lv)) cycle
             nnb = nnb + 1
             if (pass == 2) then
                cl%nb_i(nnb) = i
                cl%nb_j(nnb) = j
                cl%nb_lvec(:,nnb) = lv
                sig = 0.5d0 * (atmvdw(c%spc(c%atcel(i)%is)%z) + atmvdw(c%spc(c%atcel(j)%is)%z))
                cl%nb_sig(nnb) = sig / onesixth_inv
                cl%nb_eps(nnb) = eps_lj
             end if
          end do
       end do
       if (pass == 1) then
          cl%nnb = nnb
          allocate(cl%nb_i(nnb),cl%nb_j(nnb),cl%nb_lvec(3,nnb),cl%nb_sig(nnb),cl%nb_eps(nnb))
       end if
    end do

  end subroutine builtin_nonbonded

  !> Evaluate the built-in force field: energy, gradient, and (optionally) stress.
  subroutine builtin_evaluate(cl,c,ene,grad,stress)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    real*8, intent(out) :: ene
    real*8, intent(out) :: grad(:,:)
    real*8, intent(out), optional :: stress(3,3)

    integer :: n, i, j, ii, jj, kk
    real*8 :: d(3), r, dr, gmag, g(3)
    real*8 :: ri(3), rj(3), rk(3), a(3), b(3), na, nb, u, theta, sint, coef
    real*8 :: dudi(3), dudk(3), dudj(3), gi(3), gj(3), gk(3)
    real*8 :: sr, sr2, sr6, sr12, fmag
    real*8 :: vir(3,3)
    logical :: dostress

    dostress = present(stress)
    ene = 0d0
    grad = 0d0
    vir = 0d0

    ! harmonic bonds: E = 1/2 k (r-r0)^2
    do n = 1, cl%nbond
       i = cl%bond_i(n)
       j = cl%bond_j(n)
       d = c%atcel(i)%r - atpos(c,j,cl%bond_lvec(:,n))
       r = norm2(d)
       if (r < 1d-10) cycle
       dr = r - cl%bond_r0(n)
       ene = ene + 0.5d0 * cl%bond_k(n) * dr*dr
       gmag = cl%bond_k(n) * dr / r
       g = gmag * d
       grad(:,i) = grad(:,i) + g
       grad(:,j) = grad(:,j) - g
       if (dostress) call add_virial(vir,d,-g)
    end do

    ! harmonic angles: E = 1/2 k (theta-theta0)^2
    do n = 1, cl%nang
       ii = cl%ang_i(n)
       jj = cl%ang_j(n)
       kk = cl%ang_k(n)
       ri = atpos(c,ii,cl%ang_li(:,n))
       rj = c%atcel(jj)%r
       rk = atpos(c,kk,cl%ang_lk(:,n))
       a = ri - rj
       b = rk - rj
       na = norm2(a)
       nb = norm2(b)
       if (na < 1d-10 .or. nb < 1d-10) cycle
       u = dot_product(a,b)/(na*nb)
       u = max(min(u,1d0),-1d0)
       theta = acos(u)
       sint = sqrt(max(1d0-u*u,sin_min*sin_min))
       ! du/dri and du/drk (see module notes)
       dudi = b/(na*nb) - (u/(na*na))*a
       dudk = a/(na*nb) - (u/(nb*nb))*b
       dudj = -(dudi+dudk)
       ene = ene + 0.5d0 * cl%ang_kk(n) * (theta-cl%ang_t0(n))**2
       coef = -cl%ang_kk(n) * (theta-cl%ang_t0(n)) / sint
       gi = coef*dudi
       gj = coef*dudj
       gk = coef*dudk
       grad(:,ii) = grad(:,ii) + gi
       grad(:,jj) = grad(:,jj) + gj
       grad(:,kk) = grad(:,kk) + gk
       if (dostress) then
          call add_virial_site(vir,ri,-gi)
          call add_virial_site(vir,rj,-gj)
          call add_virial_site(vir,rk,-gk)
       end if
    end do

    ! Lennard-Jones nonbonded: E = 4 eps [(s/r)^12 - (s/r)^6]
    do n = 1, cl%nnb
       i = cl%nb_i(n)
       j = cl%nb_j(n)
       d = c%atcel(i)%r - atpos(c,j,cl%nb_lvec(:,n))
       r = norm2(d)
       if (r < 1d-10 .or. r > rcut_lj) cycle
       sr = cl%nb_sig(n)/r
       sr2 = sr*sr
       sr6 = sr2*sr2*sr2
       sr12 = sr6*sr6
       ene = ene + 4d0*cl%nb_eps(n)*(sr12-sr6)
       ! dE/dr = 4 eps (-12 sr12 + 6 sr6)/r ; gradient wrt d = (dE/dr) d/r
       fmag = 4d0*cl%nb_eps(n)*(-12d0*sr12+6d0*sr6)/r
       g = fmag * d / r
       grad(:,i) = grad(:,i) + g
       grad(:,j) = grad(:,j) - g
       if (dostress) call add_virial(vir,d,-g)
    end do

    if (dostress) then
       if (c%ismolecule .or. c%omega < 1d-10) then
          stress = 0d0
       else
          stress = -vir / c%omega
       end if
    end if

  end subroutine builtin_evaluate

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

  !> Angle (radians) between vectors a and b.
  pure function angle_between(a,b) result(t)
    real*8, intent(in) :: a(3), b(3)
    real*8 :: t, u, na, nb
    na = norm2(a)
    nb = norm2(b)
    if (na < 1d-10 .or. nb < 1d-10) then
       t = 0d0
       return
    end if
    u = dot_product(a,b)/(na*nb)
    u = max(min(u,1d0),-1d0)
    t = acos(u)
  end function angle_between

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
