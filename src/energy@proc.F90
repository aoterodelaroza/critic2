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

!> Energy/force/stress calculator.
!>
!> The built-in UFF (valence + van der Waals terms) and the QEq
!> electrostatics were adapted from GULP by Julian Gale (the General
!> Utility Lattice Program, J. D. Gale, J. Chem. Soc. Faraday
!> Trans. 93 (1997) 629; J. D. Gale and A. L. Rohl, Mol. Simul. 29
!> (2003) 291). The functional forms, atom-type parameters, and unit
!> conventions follow GULP's so that the two codes agree, and the
!> two-center Slater Coulomb integrals used by QEq (qeq_css,
!> qeq_coeffs, qeq_caintgs, qeq_cbintgs) are ported from GULP's
!> gamma.F90.
submodule (energy) proc
  use iso_c_binding
  use param, only: hartoev, bohrtoa, fact, kcal2ha
  implicit none

  ! UFF (Universal Force Field, Rappe et al., JACS 114 (1992) 10024) constants
  real*8, parameter :: rcut_lj = 10d0 / bohrtoa !< van der Waals cutoff radius, bohr
  real*8, parameter :: uff_lambda = 0.1332d0 !< bond-order correction coefficient
  real*8, parameter :: uff_kpre = 664.12d0 !< force-constant prefactor, kcal.A/mol
  ! QEq (Rappe-Goddard) electrostatics, GULP convention
  real*8, parameter :: qeq_rqeq = 15d0 / bohrtoa ! Slater short-range correction cutoff
  real*8, parameter :: qeq_rqeqtaper = 12d0 / bohrtoa ! Slater taper
  real*8, parameter :: qeq_lambda = 0.5d0 !< orbital exponent scale factor
  integer, parameter :: qeq_itermax = 30 !< max QEq iterations (charge-dependent H)
  real*8, parameter :: qeq_scfcrit = 1d-10 !< QEq charge convergence (electrons)

  ! TIP4P water-model parameters (tip4p_qh, tip4p_doh, tip4p_hoh, ...) are
  ! declared public in the parent module (src/energy.f90) and reach this
  ! submodule by host association.

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
     subroutine c_tblite_set_context_verbosity(ctx,verbosity) &
        bind(c,name="tblite_set_context_verbosity")
       import c_ptr, c_int
       type(c_ptr), value :: ctx
       integer(c_int), value :: verbosity
     end subroutine c_tblite_set_context_verbosity
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

  !> Initialize the calculator for crystal c. backend selects ff_uff (default)
  !> or ff_tblite; method selects the tblite method. elec toggles the built-in
  !> UFF QEq electrostatics (default on). errmsg is empty on success and holds
  !> the error message on failure.
  module subroutine calc_init(cl,c,backend,method,elec,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    integer, intent(in), optional :: backend
    integer, intent(in), optional :: method
    logical, intent(in), optional :: elec
    character(len=:), allocatable, intent(out) :: errmsg

    errmsg = ""
    call cl%free()
    if (present(backend)) cl%backend = backend
    if (present(method)) cl%method = method
    cl%do_elec = .true.
    if (present(elec)) cl%do_elec = elec
    cl%nat = c%ncel

    if (cl%backend == ff_uff) then
       call uff_setup(cl,c,errmsg)
    else if (cl%backend == ff_tblite) then
       call calc_init_tblite(cl,c,errmsg)
    else if (cl%backend == ff_tip4p) then
       call tip4p_setup(cl,c,errmsg)
    else
       errmsg = "unknown energy backend"
    end if
    cl%ready = (len_trim(errmsg) == 0)

  end subroutine calc_init

  !> Whether energy the backend (ff_*) can be used for system c
  module function ff_backend_applicable(backend,c) result(ok)
    use crystalmod, only: crystal
    integer, intent(in) :: backend
    class(crystal), intent(in) :: c
    logical :: ok

    integer :: i, nh, no

    ok = .false.
    select case (backend)
    case (ff_uff)
       ok = .true.
    case (ff_tblite)
#ifdef HAVE_TBLITE
       ok = .true.
#else
       ok = .false.
#endif
    case (ff_tip4p)
       if (.not.c%ismolecule) return
       nh = 0
       no = 0
       do i = 1, c%ncel
          if (c%spc(c%atcel(i)%is)%z == 1) then
             nh = nh + 1
          else if (c%spc(c%atcel(i)%is)%z == 8) then
             no = no + 1
          else
             return ! a non-H/O atom: not all-water
          end if
       end do
       ok = (no > 0 .and. nh == 2*no)
    end select

  end function ff_backend_applicable

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
    character(len=:), allocatable, intent(out) :: errmsg

    errmsg = ""
    ene = 0d0
    grad = 0d0
    if (present(stress)) stress = 0d0

    if (.not.cl%ready) then
       errmsg = "calculator not initialized"
       return
    end if

    if (cl%backend == ff_uff) then
       call uff_evaluate(cl,c,ene,grad,stress)
    else if (cl%backend == ff_tblite) then
       call calc_eval_tblite(cl,c,ene,grad,stress,errmsg)
    else if (cl%backend == ff_tip4p) then
       call tip4p_evaluate(cl,c,ene,grad,stress)
    else
       errmsg = "unknown energy backend"
    end if

  end subroutine calc_evaluate

  !> Refresh backend caches after the geometry of c changes.
  module subroutine calc_update_geometry(cl,c)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c

    if (.not.cl%ready) return
    if (cl%backend == ff_uff) then
       call uff_nonbonded(cl,c)
       if (cl%do_elec) call uff_qeq_pairs(cl,c)
    else if (cl%backend == ff_tblite) then
       call calc_update_tblite(cl,c)
    end if
    ! ff_tip4p: nothing to update

  end subroutine calc_update_geometry

  !> Terminate a calculator.
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
    if (allocated(cl%qzeta)) deallocate(cl%qzeta)
    if (allocated(cl%qpqn)) deallocate(cl%qpqn)
    if (allocated(cl%qnbptr)) deallocate(cl%qnbptr)
    if (allocated(cl%qnbj)) deallocate(cl%qnbj)
    if (allocated(cl%qnblv)) deallocate(cl%qnblv)
    cl%nqnb = 0
    if (allocated(cl%bond)) deallocate(cl%bond)
    if (allocated(cl%ang)) deallocate(cl%ang)
    if (allocated(cl%nb)) deallocate(cl%nb)
    if (allocated(cl%tor)) deallocate(cl%tor)
    if (allocated(cl%inv)) deallocate(cl%inv)
    cl%nwat = 0
    if (allocated(cl%iwat)) deallocate(cl%iwat)

  end subroutine calc_free

  !xx! private procedures

  !xx! UFF (Universal Force Field, Rappe et al., JACS 114 (1992) 10024) constants

  !> Build the UFF topology (bonds, angles, van der Waals pairs) from the
  !> connectivity graph. Atom types are assigned from the connectivity and all
  !> parameters converted to atomic units. On error (an unparameterized element)
  !> errmsg is set (empty on success).
  subroutine uff_setup(cl,c,errmsg)
    use crystalmod, only: crystal
    use types, only: neighstar
    use param, only: rad
    use tools_io, only: string
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out) :: errmsg

    type(neighstar), allocatable :: ns(:)
    integer :: i, j, k, ka, kb, ni, nk, nb, na, ipass
    integer :: jj, kk, a, b, nl, hj, hk, nt, norder, ic
    integer :: lv(3), li(3), lk(3), ljk(3), lkl(3)
    real*8 :: rij, rjk, rik, t0, t0deg, cost0, st0sq, cc2
    real*8 :: v, costerm, c0, c1, c2, kinv
    logical :: ok

    errmsg = ""

    ! connectivity graph (use the crystal's own if available)
    if (allocated(c%nstar)) then
       ns = c%nstar
    else
       call c%find_asterisms(ns)
    end if

    ! assign the UFF parameters to every atom from its element and connectivity
    allocate(cl%attype(cl%nat))
    do i = 1, cl%nat
       call uff_atom_params(c,i,ns,cl%attype(i))
       if (len_trim(cl%attype(i)%label) == 0) then
          errmsg = "no UFF parameters for atomic number " // string(c%spc(c%atcel(i)%is)%z)
          return
       end if
    end do

    ! ---- bonds ----
    do ipass = 1, 2
       nb = 0
       do i = 1, cl%nat
          do k = 1, ns(i)%ncon
             j = ns(i)%idcon(k)
             lv = ns(i)%lcon(:,k)
             if (.not.bond_canonical(i,j,lv)) cycle
             nb = nb + 1
             if (ipass == 2) then
                cl%bond(nb)%i = i
                cl%bond(nb)%j = j
                cl%bond(nb)%lvec = lv
                rij = uff_natural_length(cl%attype(i),cl%attype(j),bond_order(ns(i),k))
                cl%bond(nb)%r0 = rij / bohrtoa
                ! k_IJ = 664.12 Z*_I Z*_J / r_IJ^3 [kcal/mol/A^2] -> hartree/bohr^2
                cl%bond(nb)%k = uff_kpre * cl%attype(i)%zstar * cl%attype(j)%zstar &
                   / rij**3 * kcal2ha * bohrtoa*bohrtoa
             end if
          end do
       end do
       if (ipass == 1) then
          cl%nbond = nb
          allocate(cl%bond(nb))
       end if
    end do

    ! ---- angles (i-j-k, vertex j) ----
    ! the bending energy is kk*P(u), u=cos(theta); P is precomputed from theta0
    na = 0
    do j = 1, cl%nat
       na = na + ns(j)%ncon*(ns(j)%ncon-1)/2
    end do
    cl%nang = na
    allocate(cl%ang(na))
    na = 0
    do j = 1, cl%nat
       t0deg = cl%attype(j)%theta0
       t0 = t0deg * rad
       cost0 = cos(t0)
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
             ! natural bond lengths of the two arms (angstrom); r_IK by law of cosines
             rij = uff_natural_length(cl%attype(ni),cl%attype(j),bond_order(ns(j),ka))
             rjk = uff_natural_length(cl%attype(nk),cl%attype(j),bond_order(ns(j),kb))
             rik = sqrt(max(rij*rij + rjk*rjk - 2d0*rij*rjk*cost0,1d-10))
             ! k_IJK = 664.12 (Z*_I Z*_K/r_IK^5) [3 r_IJ r_JK (1-cos^2 t0) - r_IK^2 cos t0]
             ! (GULP/Rappe convention: no extra r_IJ r_JK prefactor)
             cl%ang(na)%kk = uff_kpre * cl%attype(ni)%zstar * cl%attype(nk)%zstar &
                / rik**5 * (3d0*rij*rjk*(1d0-cost0*cost0) - rik*rik*cost0) * kcal2ha
             ! coefficients of P(u): general cosine-Fourier form, or GULP's periodic
             ! special cases selected by theta0 (same 1e-6 deg threshold GULP uses)
             cl%ang(na)%p = 0d0
             if (abs(t0deg-180d0) < 1d-6) then
                ! 1 + cos(theta)
                cl%ang(na)%p(0) = 1d0
                cl%ang(na)%p(1) = 1d0
             else if (abs(t0deg) < 1d-6) then
                ! 1 - cos(theta)
                cl%ang(na)%p(0) = 1d0
                cl%ang(na)%p(1) = -1d0
             else if (abs(t0deg-120d0) < 1d-6) then
                ! (1 - cos3theta)/9
                cl%ang(na)%p(0) = 1d0/9d0
                cl%ang(na)%p(1) = 1d0/3d0
                cl%ang(na)%p(3) = -4d0/9d0
             else if (abs(t0deg-90d0) < 1d-6) then
                ! (1 - cos4theta)/16
                cl%ang(na)%p(2) = 0.5d0
                cl%ang(na)%p(4) = -0.5d0
             else
                ! general: C0 + C1 u + C2 (2u^2-1) with C2=1/(4 sin^2 t0)
                st0sq = 1d0 - cost0*cost0
                cc2 = 1d0/(4d0*st0sq)
                cl%ang(na)%p(0) = cc2*(2d0*cost0*cost0 + 1d0) - cc2
                cl%ang(na)%p(1) = -4d0*cc2*cost0
                cl%ang(na)%p(2) = 2d0*cc2
             end if
          end do
       end do
    end do

    ! ---- torsions (dihedrals about each central bond j-k, both sp2/sp3) ----
    do ipass = 1, 2
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
                   if (ipass == 2) then
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
       if (ipass == 1) then
          cl%ntor = nt
          allocate(cl%tor(nt))
       end if
    end do

    ! ---- inversions (out-of-plane) for 3-coordinate C/N (sp2) and group-15 ----
    do ipass = 1, 2
       ni = 0
       do ic = 1, cl%nat
          if (ns(ic)%ncon /= 3) cycle
          call uff_inversion_params(cl,c,ns,ic,kinv,c0,c1,c2,ok)
          if (.not.ok) cycle
          ni = ni + 1
          if (ipass == 2) then
             cl%inv(ni)%c = ic
             cl%inv(ni)%n1 = ns(ic)%idcon(1)
             cl%inv(ni)%n2 = ns(ic)%idcon(2)
             cl%inv(ni)%n3 = ns(ic)%idcon(3)
             cl%inv(ni)%l1 = ns(ic)%lcon(:,1)
             cl%inv(ni)%l2 = ns(ic)%lcon(:,2)
             cl%inv(ni)%l3 = ns(ic)%lcon(:,3)
             cl%inv(ni)%k = kinv * kcal2ha / 3d0 ! averaged over the 3 apex choices
             cl%inv(ni)%c0 = c0
             cl%inv(ni)%c1 = c1
             cl%inv(ni)%c2 = c2
          end if
       end do
       if (ipass == 1) then
          cl%ninv = ni
          allocate(cl%inv(ni))
       end if
    end do

    ! ---- QEq partial charges (also builds the cached electrostatics list) ----
    if (cl%do_elec) call uff_qeq(cl,c)

    ! ---- van der Waals + electrostatics pairs ----
    call uff_nonbonded(cl,c)

  contains

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

      v = 0d0
      norder = 1
      costerm = 1d0
      ok = .true.
      hj = uff_hyb(cl%attype(jj))
      hk = uff_hyb(cl%attype(kk))
      zj = c%spc(c%atcel(jj)%is)%z
      zk = c%spc(c%atcel(kk)%is)%z

      if (hj == 3 .and. hk == 3) then
         ! both sp3
         if (abs(bojk-1d0) < 1d-3 .and. isgroup6(zj) .and. isgroup6(zk)) then
            v = sqrt(group6_v(zj)*group6_v(zk))
            norder = 2
            costerm = -1d0
         else
            v = sqrt(cl%attype(jj)%vtor * cl%attype(kk)%vtor)
            norder = 3
            costerm = -1d0
         end if
      else if (hj == 2 .and. hk == 2) then
         ! both sp2: bond-order dependent barrier (Rappe eq. 17)
         v = 5d0*sqrt(cl%attype(jj)%utor * cl%attype(kk)%utor) &
            *(1d0 + 4.18d0*log(bojk))
         norder = 2
         costerm = 1d0
      else
         ! sp2-sp3 mixed
         if (hj == 3) then
            is3 = jj
            is2 = kk
         else
            is3 = kk
            is2 = jj
         end if
         z3 = c%spc(c%atcel(is3)%is)%z
         z2 = c%spc(c%atcel(is2)%is)%z
         if (isgroup6(z3) .and. .not.isgroup6(z2)) then
            v = 5d0*sqrt(cl%attype(is3)%utor * cl%attype(is2)%utor) &
               *(1d0 + 4.18d0*log(bojk))
            norder = 2
            costerm = -1d0
         else
            ! sp3-sp2-sp2 conjugation: the sp2 center bonded to another sp2
            sp2sp2 = .false.
            do m = 1, ns(is2)%ncon
               if (ns(is2)%idcon(m) == is3) cycle
               if (uff_hyb(cl%attype(ns(is2)%idcon(m))) == 2) sp2sp2 = .true.
            end do
            if (sp2sp2) then
               v = 2d0
               norder = 3
               costerm = -1d0
            else
               v = 1d0
               norder = 6
               costerm = 1d0
            end if
         end if
      end if
      ! a zero barrier contributes nothing; drop it
      if (abs(v) < 1d-12) ok = .false.

    end subroutine uff_torsion_params

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

      ok = .true.
      kk = 0d0
      c0 = 0d0
      c1 = 0d0
      c2 = 0d0
      z = c%spc(c%atcel(ic)%is)%z
      h = uff_hyb(cl%attype(ic))

      select case(z)
      case(6) ! C: planar sp2 (50 kcal/mol if bonded to O, else 6)
         if (h /= 2) then
            ok = .false.
            return
         end if
         c0 = 1d0
         c1 = -1d0
         c2 = 0d0
         boundo = .false.
         do m = 1, ns(ic)%ncon
            if (c%spc(c%atcel(ns(ic)%idcon(m))%is)%z == 8) boundo = .true.
         end do
         if (boundo) then
            kk = 50d0
         else
            kk = 6d0
         end if
      case(7) ! N: planar sp2 (amide/aromatic)
         if (h /= 2) then
            ok = .false.
            return
         end if
         c0 = 1d0
         c1 = -1d0
         c2 = 0d0
         kk = 6d0
      case(15,33,51,83) ! P, As, Sb, Bi: pyramidal
         select case(z)
         case(15)
            w0 = 84.4339d0
         case(33)
            w0 = 86.9735d0
         case(51)
            w0 = 87.7047d0
         case default
            w0 = 90d0
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
      type(uffparam), intent(in) :: iti, itj
      real*8, intent(in) :: bo
      real*8 :: rij, ri, rj, ci, cj, rbo, ren

      ri = iti%r1
      rj = itj%r1
      ci = iti%chi
      cj = itj%chi
      rbo = -uff_lambda*(ri+rj)*log(bo)
      ren = 0d0
      if (ci*ri + cj*rj > 1d-10) &
         ren = ri*rj*(sqrt(ci)-sqrt(cj))**2 / (ci*ri + cj*rj)
      rij = ri + rj + rbo - ren

    end function uff_natural_length

    !> Hybridization class of a UFF atom type from its geometry code: 1 = linear
    !> (sp), 2 = trigonal/resonant (sp2, includes R), 3 = tetrahedral (sp3);
    !> 0 = other (linear terminal, ionic, square-planar/octahedral metal, ...).
    pure function uff_hyb(ip) result(h)
      type(uffparam), intent(in) :: ip
      integer :: h
      integer :: i
      character :: ch
      h = 0
      do i = 1, len_trim(ip%label)
         ch = ip%label(i:i)
         select case(ch)
         case("1")
            h = 1
            return
         case("2","R")
            h = 2
            return
         case("3")
            h = 3
            return
         case("4","5","6")
            h = 0
            return
         end select
      end do
    end function uff_hyb

    !> True if z is a group-16 element.
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

  end subroutine uff_setup

  !> Evaluate the UFF energy, gradient, and (optionally) stress, all
  !> in atomic units.
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
    real*8 :: ri(3), rj(3), rk(3), a(3), b(3), na, nb, u, edu, pu, dpu
    real*8 :: dudi(3), dudk(3), dudj(3), gi(3), gj(3), gk(3)
    real*8 :: sr, sr2, sr6, sr12, fmag
    real*8 :: vir(3,3)
    real*8 :: esofar, ecoul
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
    !$omp parallel do private(n,i,j,d,r,dr,gmag,g) reduction(+:ene,grad,vir)
    do n = 1, cl%nbond
       i = cl%bond(n)%i
       j = cl%bond(n)%j
       d = c%atcel(i)%r - c%atcel(j)%r - c%x2c(real(cl%bond(n)%lvec,8))
       r = norm2(d)
       if (r < 1d-10) cycle
       dr = r - cl%bond(n)%r0
       ene = ene + 0.5d0 * cl%bond(n)%k * dr*dr
       gmag = cl%bond(n)%k * dr / r
       g = gmag * d
       grad(:,i) = grad(:,i) + g
       grad(:,j) = grad(:,j) - g
       if (dostress) call virial(vir,d,-g)
    end do
    !$omp end parallel do
    esofar = ene

    ! angle bend, expressed in u = cos(theta) (no 1/sin(theta) singularity):
    ! E = kk * P(u) with the polynomial P precomputed at setup
    !$omp parallel do private(n,ii,jj,kk,ri,rj,rk,a,b,na,nb,u,dudi,dudk,dudj,pu,dpu,edu,gi,gj,gk) &
    !$omp    reduction(+:ene,grad,vir)
    do n = 1, cl%nang
       ii = cl%ang(n)%i
       jj = cl%ang(n)%j
       kk = cl%ang(n)%k
       ri = c%atcel(ii)%r + c%x2c(real(cl%ang(n)%li,8))
       rj = c%atcel(jj)%r
       rk = c%atcel(kk)%r + c%x2c(real(cl%ang(n)%lk,8))
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
       ! P(u) and dP/du (Horner), u = cos(theta)
       pu = cl%ang(n)%p(0) + u*(cl%ang(n)%p(1) + u*(cl%ang(n)%p(2) + u*(cl%ang(n)%p(3) + u*cl%ang(n)%p(4))))
       dpu = cl%ang(n)%p(1) + u*(2d0*cl%ang(n)%p(2) + u*(3d0*cl%ang(n)%p(3) + u*4d0*cl%ang(n)%p(4)))
       ene = ene + cl%ang(n)%kk * pu
       edu = cl%ang(n)%kk * dpu
       gi = edu*dudi
       gj = edu*dudj
       gk = edu*dudk
       grad(:,ii) = grad(:,ii) + gi
       grad(:,jj) = grad(:,jj) + gj
       grad(:,kk) = grad(:,kk) + gk
       if (dostress) then
          call virial(vir,ri,-gi)
          call virial(vir,rj,-gj)
          call virial(vir,rk,-gk)
       end if
    end do
    !$omp end parallel do
    esofar = ene

    ! van der Waals (UFF 12-6): E = D [(x/r)^12 - 2 (x/r)^6]
    !$omp parallel do private(n,i,j,d,r,sr,sr2,sr6,sr12,fmag,g) reduction(+:ene,grad,vir)
    do n = 1, cl%nnb
       i = cl%nb(n)%i
       j = cl%nb(n)%j
       d = c%atcel(i)%r - c%atcel(j)%r - c%x2c(real(cl%nb(n)%lvec,8))
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
       if (dostress) call virial(vir,d,-g)
    end do
    !$omp end parallel do
    esofar = ene

    ! QEq electrostatics (Ewald 1/r + self + Slater short-range correction)
    if (cl%do_elec) then
       call uff_electrostatics(cl,c,ene,grad,vir,dostress,ecoul)
    end if
    esofar = ene

    ! torsions: E = V/2 (1 - cosTerm cos(n phi)), written via c = cos(phi) with
    ! the Chebyshev expansion so no sin/atan2 sign convention is needed
    !$omp parallel do &
    !$omp    private(t,ti,tj,tk,tl,ri,rj,rk,rl,tv1,tv2,tv3,cp1,cp2,n1,n2,cphi,vt,nn,ct,cc2,tn,dtn,dedc,dva,dvb,dc1,dc2,dc3,gi,gj,gk,gl) &
    !$omp    reduction(+:ene,grad,vir)
    do t = 1, cl%ntor
       ti = cl%tor(t)%i
       tj = cl%tor(t)%j
       tk = cl%tor(t)%k
       tl = cl%tor(t)%l
       ri = c%atcel(ti)%r + c%x2c(real(cl%tor(t)%li,8))
       rj = c%atcel(tj)%r
       rk = c%atcel(tk)%r + c%x2c(real(cl%tor(t)%lk,8))
       rl = c%atcel(tl)%r + c%x2c(real(cl%tor(t)%ll,8))
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
          call virial(vir,ri,-gi)
          call virial(vir,rj,-gj)
          call virial(vir,rk,-gk)
          call virial(vir,rl,-gl)
       end if
    end do
    !$omp end parallel do
    esofar = ene

    ! inversions (out-of-plane) at 3-coordinate centers, averaged over the 3
    ! apex choices: E = K [C0 + C1 sin(Y) + C2 cos(2 Y)]
    !$omp parallel do private(inv,ti,iidx,rc,rnb,ap,ib1,ib2,cosy,siny,edcy,gc,gg1,gg2,gg3) &
    !$omp    reduction(+:ene,grad,vir)
    do inv = 1, cl%ninv
       ti = cl%inv(inv)%c
       iidx = (/cl%inv(inv)%n1,cl%inv(inv)%n2,cl%inv(inv)%n3/)
       rc = c%atcel(ti)%r
       rnb(:,1) = c%atcel(cl%inv(inv)%n1)%r + c%x2c(real(cl%inv(inv)%l1,8))
       rnb(:,2) = c%atcel(cl%inv(inv)%n2)%r + c%x2c(real(cl%inv(inv)%l2,8))
       rnb(:,3) = c%atcel(cl%inv(inv)%n3)%r + c%x2c(real(cl%inv(inv)%l3,8))
       do ap = 1, 3
          if (ap == 1) then
             ib1 = 2
             ib2 = 3
          else if (ap == 2) then
             ib1 = 1
             ib2 = 3
          else
             ib1 = 1
             ib2 = 2
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
             call virial(vir,rc,-edcy*gc)
             call virial(vir,rnb(:,ib1),-edcy*gg1)
             call virial(vir,rnb(:,ib2),-edcy*gg2)
             call virial(vir,rnb(:,ap),-edcy*gg3)
          end if
       end do
    end do
    !$omp end parallel do

    if (dostress) then
       if (c%ismolecule .or. c%omega < 1d-10) then
          stress = 0d0
       else
          stress = -vir / c%omega
       end if
    end if

  contains
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
         cosy = 0d0
         gc = 0d0
         gj = 0d0
         gk = 0d0
         gl = 0d0
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

    !> Add the QEq electrostatic energy (and, if grad present, forces and virial)
    !> at the current geometry with the fixed solved charges cl%qeq: the periodic
    !> Ewald 1/r energy, the QEq self-energy, and the short-range Slater
    !> correction. Atomic units. eelec returns the total electrostatic energy.
    subroutine uff_electrostatics(cl,c,ene,grad,vir,dostress,eelec)
      use crystalmod, only: crystal
      class(calculator), intent(inout) :: cl
      class(crystal), intent(inout) :: c
      real*8, intent(inout) :: ene
      real*8, intent(inout) :: grad(:,:)
      real*8, intent(inout) :: vir(3,3)
      logical, intent(in) :: dostress
      real*8, intent(out) :: eelec

      integer :: n, i, j, p, z
      real*8 :: ecoul, eself, eslat, gam, dgam, dzadum, qi, qj, chi, mu, zh0, gmag
      real*8 :: d(3), r

      n = cl%nat
      eelec = 0d0
      if (.not.cl%do_elec .or. .not.allocated(cl%qeq)) return

      ! keep fractional coordinates consistent with the cartesian ones: the Ewald
      ! sum and the neighbour search below key off %x, while the rest of the
      ! evaluator uses %r, and a caller (MD/relaxation) may move only %r
      do i = 1, n
         c%atcel(i)%x = c%c2x(c%atcel(i)%r)
      end do

      ! --- Coulomb 1/r: Ewald for periodic crystals; direct sum (below) for molecules ---
      ecoul = 0d0
      if (.not.c%ismolecule) &
         call c%ewald_energy_grad(cl%qeq,cl%qrcut,cl%qhcut,cl%qeta,cl%qlrmax,cl%qlhmax,ecoul,grad,vir,dostress,&
            nbptr=cl%qnbptr,nbj=cl%qnbj,nblv=cl%qnblv)

      ! --- QEq self energy (chi q + mu q^2; hydrogen higher-order) ---
      eself = 0d0
      !$omp parallel do private(i,z,qi,chi,mu,zh0) reduction(+:eself)
      do i = 1, n
         z = c%spc(c%atcel(i)%is)%z
         qi = cl%qeq(i)
         chi = cl%attype(i)%chi / hartoev
         mu = cl%attype(i)%jhard / hartoev
         if (z == 1) then
            zh0 = 0.75d0/(cl%attype(i)%qrad / bohrtoa)
            eself = eself + qi*(chi + qi*mu*(1d0 + 2d0*qi/(3d0*zh0)))
         else
            eself = eself + qi*(chi + qi*mu)
         end if
      end do
      !$omp end parallel do

      ! --- short-range: crystal = Slater correction (J-1/r) over images added to
      !     the Ewald 1/r; molecule = full direct Coulomb J over all finite pairs ---
      eslat = 0d0
      if (c%ismolecule) then
         !$omp parallel do private(i,j,qi,qj,d,r,gam,dgam,dzadum,gmag) reduction(+:ecoul,grad)
         do i = 1, n
            qi = cl%qeq(i)
            do j = 1, n
               if (i == j) cycle
               d = c%atcel(i)%r - c%atcel(j)%r
               r = norm2(d)
               if (r < 1d-10) cycle
               qj = cl%qeq(j)
               call qeq_slater(cl%qpqn(i),cl%qpqn(j),cl%qzeta(i),cl%qzeta(j),r,gam,dgam,dzadum)
               ! full two-centre integral J = (J-1/r) + 1/r and its r-derivative
               ecoul = ecoul + 0.5d0*qi*qj*(gam + 1d0/r)
               gmag = 0.5d0*qi*qj*(dgam - 1d0/(r*r))/r
               grad(:,i) = grad(:,i) + gmag*d
               grad(:,j) = grad(:,j) - gmag*d
            end do
         end do
         !$omp end parallel do
      else
         !$omp parallel do private(i,p,j,qi,qj,d,r,gam,dgam,dzadum,gmag) reduction(+:eslat,grad,vir)
         do i = 1, n
            qi = cl%qeq(i)
            do p = cl%qnbptr(i), cl%qnbptr(i+1)-1
               j = cl%qnbj(p)
               ! use the vector to the same image the force is taken along, so r and d agree
               d = c%atcel(i)%r - c%atcel(j)%r - c%x2c(real(cl%qnblv(:,p),8))
               r = norm2(d)
               if (r < 1d-10 .or. r > qeq_rqeq) cycle
               qj = cl%qeq(j)
               call qeq_slater(cl%qpqn(i),cl%qpqn(j),cl%qzeta(i),cl%qzeta(j),r,gam,dgam,dzadum)
               eslat = eslat + 0.5d0*qi*qj*gam
               ! force from the pairwise 0.5 q_i q_j J(r) term (double-counted over i,j)
               gmag = 0.5d0*qi*qj*dgam/r
               grad(:,i) = grad(:,i) + gmag*d
               grad(:,j) = grad(:,j) - gmag*d
               if (dostress) call virial(vir,d,-gmag*d)
            end do
         end do
         !$omp end parallel do
      end if

      eelec = ecoul + eself + eslat
      ene = ene + eelec

    end subroutine uff_electrostatics
  end subroutine uff_evaluate

  !> Equilibrate partial charges by the QEq method of Rappe and Goddard.
  subroutine uff_qeq(cl,c)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c

    integer :: n, i, info, iter, z
    integer, allocatable :: ipiv(:), npqn(:), zat(:)
    real*8, allocatable :: aew(:,:), amat0(:,:), amat(:,:), bvec(:), q(:), qold(:)
    real*8, allocatable :: zeta0(:), zeta(:), mu(:), chi(:), qrad(:), qfrac(:)
    real*8 :: zh0, qdiff, qsum
    real*8, parameter :: qmax = 8d0 ! reject unphysical charges (electrons)

    n = cl%nat
    if (n <= 0) then
       cl%do_elec = .false.
       return
    end if

    ! per-atom QEq parameters (atomic units)
    allocate(npqn(n),zat(n),zeta0(n),zeta(n),mu(n),chi(n),qrad(n))
    do i = 1, n
       z = c%spc(c%atcel(i)%is)%z
       zat(i) = z
       npqn(i) = qeq_pqn(z)
       chi(i) = cl%attype(i)%chi / hartoev            ! hartree
       mu(i) = cl%attype(i)%jhard / hartoev           ! hartree
       qrad(i) = cl%attype(i)%qrad / bohrtoa          ! bohr
       if (qrad(i) < 1d-6) then
          cl%do_elec = .false.
          return
       end if
       zeta0(i) = 0.5d0*qeq_lambda*dble(2*npqn(i)+1)/qrad(i)  ! 1/bohr
    end do

    ! periodic crystals: cache Ewald parameters, 1/r matrix, and neighbour list.
    ! finite molecules: no Ewald, the Coulomb is a direct sum.
    if (.not.c%ismolecule) then
       ! cache the Ewald parameters (width eta, real/reciprocal cutoffs and cell
       ! ranges); they depend only on the cell, so are reused for every evaluation
       ! during MD/relaxation, which also keeps eta fixed so the stress stays
       ! consistent. qfrac=1 makes the cutoffs charge-independent.
       allocate(qfrac(c%nspc))
       qfrac = 1d0
       call c%calculate_ewald_cutoffs(qfrac,cl%qrcut,cl%qhcut,cl%qeta,qsum,cl%qlrmax,cl%qlhmax)
       call uff_qeq_pairs(cl,c)
       allocate(aew(n,n))
       call c%ewald_matrix(aew,cl%qrcut,cl%qhcut,cl%qeta,cl%qlrmax,cl%qlhmax)
    end if

    ! charge-independent part of the QEq matrix, built once: the Ewald 1/r matrix
    ! plus the Slater terms of pairs whose exponents do not depend on charge
    ! (neither atom is hydrogen) plus the non-hydrogen hardness diagonal. Only the
    ! hydrogen-involving pairs and diagonal are rebuilt each iteration below.
    allocate(amat0(n,n))
    amat0 = 0d0
    if (.not.c%ismolecule) amat0 = aew
    call qeq_add_slater(cl,c,npqn,zat,zeta0,.false.,zeta0,amat0)  ! q unused when honly=.false.
    do i = 1, n
       if (zat(i) /= 1) amat0(i,i) = amat0(i,i) + 2d0*mu(i)
    end do

    ! iterative solve; only the hydrogen (charge-dependent) terms change
    allocate(amat(n+1,n+1),bvec(n+1),ipiv(n+1),q(n),qold(n))
    q = 0d0
    do iter = 1, qeq_itermax
       do i = 1, n
          zeta(i) = zeta0(i)
          if (zat(i) == 1) zeta(i) = zeta(i) + q(i)   ! r5 = 1 bohr^-1 in a.u.
       end do
       amat = 0d0
       amat(1:n,1:n) = amat0
       ! hydrogen-involving Slater pairs (charge-dependent exponents) + coupling
       call qeq_add_slater(cl,c,npqn,zat,zeta,.true.,q,amat(1:n,1:n))
       ! hydrogen hardness diagonal (charge-dependent); non-H diagonal is in amat0
       do i = 1, n
          if (zat(i) == 1) then
             zh0 = 0.75d0/qrad(i)   ! = 0.529177*0.75/R(ang) in atomic units
             amat(i,i) = amat(i,i) + 2d0*mu(i)*(1d0 + q(i)/zh0)
          end if
       end do
       ! charge-neutrality constraint via a Lagrange multiplier
       do i = 1, n
          amat(i,n+1) = -1d0
          amat(n+1,i) = 1d0
       end do
       bvec(1:n) = -chi(1:n)
       bvec(n+1) = 0d0
       call dgesv(n+1,1,amat,n+1,ipiv,bvec,n+1,info)
       if (info /= 0) then
          cl%do_elec = .false.
          return
       end if
       qold = q
       q = bvec(1:n)
       qdiff = maxval(abs(q - qold))
       if (qdiff < qeq_scfcrit) exit
    end do

    if (maxval(abs(q)) > qmax) then
       cl%do_elec = .false.
       return
    end if
    ! store charges and the converged orbital exponents / quantum numbers
    if (allocated(cl%qeq)) deallocate(cl%qeq)
    if (allocated(cl%qzeta)) deallocate(cl%qzeta)
    if (allocated(cl%qpqn)) deallocate(cl%qpqn)
    allocate(cl%qeq(n),cl%qzeta(n),cl%qpqn(n))
    cl%qeq = q
    cl%qpqn = npqn
    do i = 1, n
       cl%qzeta(i) = zeta0(i)
       if (zat(i) == 1) cl%qzeta(i) = cl%qzeta(i) + q(i)
    end do

  end subroutine uff_qeq

  !> Build the UFF van der Waals (LJ 12-6) pair list within rcut_lj, excluding
  !> 1-2 (bonded) and 1-3 (share a common neighbor) pairs. Uses the atom types
  !> stored in cl%attype (set by uff_setup).
  subroutine uff_nonbonded(cl,c)
    use crystalmod, only: crystal
    use types, only: neighstar
    use param, only: icrd_crys
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c

    type(neighstar), allocatable :: nslocal(:)
    integer :: i, j, m, nat, nnb, ipass
    integer, allocatable :: eid(:), lvec(:,:)
    real*8, allocatable :: dist(:)
    integer :: lv(3)
    real*8 :: xij, dij
    logical :: havens, is13

    if (.not.allocated(cl%attype)) return

    ! read the connectivity in place (this routine is re-run periodically during
    ! MD, so avoid deep-copying the whole neighbor-star array each time)
    havens = allocated(c%nstar)
    if (.not.havens) call c%find_asterisms(nslocal)

    if (allocated(cl%nb)) deallocate(cl%nb)

    do ipass = 1, 2
       nnb = 0
       do i = 1, cl%nat
          ! neighbour candidates within the vdW cutoff
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
             if (ipass == 2) then
                cl%nb(nnb)%i = i
                cl%nb(nnb)%j = j
                cl%nb(nnb)%lvec = lv
                ! UFF combining rules: geometric mean of distance and well depth
                xij = sqrt(cl%attype(i)%x1 * cl%attype(j)%x1)     ! angstrom
                dij = sqrt(cl%attype(i)%dwell * cl%attype(j)%dwell) ! kcal/mol
                cl%nb(nnb)%rmin = xij / bohrtoa ! bohr
                cl%nb(nnb)%deps = dij * kcal2ha ! hartree
             end if
          end do
       end do
       if (ipass == 1) then
          cl%nnb = nnb
          allocate(cl%nb(nnb))
       end if
    end do

  contains
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

  end subroutine uff_nonbonded

  !> Build the cached electrostatics neighbour list (CSR) out to the larger of
  !> the Slater cutoff and the Ewald real-space cutoff. Rebuilt periodically
  !> during MD (via calc_update_geometry) so the per-step evaluation does not
  !> repeat the spatial neighbour search. Molecules keep no list (direct sum).
  subroutine uff_qeq_pairs(cl,c)
    use crystalmod, only: crystal
    use param, only: icrd_crys
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    integer :: n, i, m, nat, tot, p
    integer, allocatable :: eid(:), lvln(:,:)
    real*8, allocatable :: dist(:)
    real*8 :: rmax

    cl%nqnb = 0
    if (c%ismolecule) return
    n = cl%nat
    rmax = max(qeq_rqeq,cl%qrcut)
    do i = 1, n
       c%atcel(i)%x = c%c2x(c%atcel(i)%r)
    end do
    if (allocated(cl%qnbptr)) deallocate(cl%qnbptr)
    allocate(cl%qnbptr(n+1))
    ! pass 1: count neighbours per atom
    tot = 0
    do i = 1, n
       call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.false.,nat,eid=eid,dist=dist,&
          lvec=lvln,up2d=rmax,nozero=.true.)
       cl%qnbptr(i) = tot + 1
       tot = tot + nat
    end do
    cl%qnbptr(n+1) = tot + 1
    cl%nqnb = tot
    if (allocated(cl%qnbj)) deallocate(cl%qnbj)
    if (allocated(cl%qnblv)) deallocate(cl%qnblv)
    allocate(cl%qnbj(max(tot,1)),cl%qnblv(3,max(tot,1)))
    ! pass 2: fill
    p = 0
    do i = 1, n
       call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.false.,nat,eid=eid,dist=dist,&
          lvec=lvln,up2d=rmax,nozero=.true.)
       do m = 1, nat
          p = p + 1
          cl%qnbj(p) = eid(m)
          cl%qnblv(:,p) = lvln(:,m)
       end do
    end do
  end subroutine uff_qeq_pairs

  !> UFF atom-type parameters for atom iat, resolved from its element Z and
  !> connectivity (number of neighbors ncon and aromaticity). The parameter
  !> table (Rappe et al., JACS 114 (1992) 10024) is hardwired here; only the
  !> atom types the connectivity-based typing can select are included. On return
  !> p%label == "" if the element has no UFF parameters.
  subroutine uff_atom_params(c,iat,ns,p)
    use crystalmod, only: crystal
    use types, only: neighstar
    class(crystal), intent(in) :: c
    integer, intent(in) :: iat
    type(neighstar), intent(in) :: ns(:)
    type(uffparam), intent(out) :: p

    integer :: z, ncon
    logical :: arom

    z = c%spc(c%atcel(iat)%is)%z
    ncon = ns(iat)%ncon
    arom = ns(iat)%isaromatic

    select case (z)
    case (5) ! B
       if (ncon >= 4) then
          p = uffparam("B_3", 0.838d0, 109.4712d0, 4.083d0, 0.180d0, 1.7550d0, 0.000d0, 2.000d0, 4.068d0, 4.750d0, 0.822d0)
       else
          p = uffparam("B_2", 0.828d0, 120.0d0, 4.083d0, 0.180d0, 1.7550d0, 0.000d0, 2.000d0, 4.068d0, 4.750d0, 0.822d0)
       end if
    case (6) ! C
       if (arom) then
          p = uffparam("C_R", 0.729d0, 120.0d0, 3.851d0, 0.105d0, 1.9120d0, 0.000d0, 2.000d0, 5.343d0, 5.063d0, 0.759d0)
       else if (ncon >= 4) then
          p = uffparam("C_3", 0.757d0, 109.4712d0, 3.851d0, 0.105d0, 1.9120d0, 2.119d0, 2.000d0, 5.343d0, 5.063d0, 0.759d0)
       else if (ncon == 3) then
          p = uffparam("C_2", 0.732d0, 120.0d0, 3.851d0, 0.105d0, 1.9120d0, 0.000d0, 2.000d0, 5.343d0, 5.063d0, 0.759d0)
       else
          p = uffparam("C_1", 0.706d0, 180.0d0, 3.851d0, 0.105d0, 1.9120d0, 0.000d0, 2.000d0, 5.343d0, 5.063d0, 0.759d0)
       end if
    case (7) ! N
       if (arom) then
          p = uffparam("N_R", 0.699d0, 120.0d0, 3.660d0, 0.069d0, 2.5438d0, 0.000d0, 2.000d0, 6.899d0, 5.880d0, 0.715d0)
       else if (ncon >= 3) then
          p = uffparam("N_3", 0.700d0, 106.7d0, 3.660d0, 0.069d0, 2.5438d0, 0.450d0, 2.000d0, 6.899d0, 5.880d0, 0.715d0)
       else if (ncon == 2) then
          p = uffparam("N_2", 0.685d0, 111.3d0, 3.660d0, 0.069d0, 2.5438d0, 0.000d0, 2.000d0, 6.899d0, 5.880d0, 0.715d0)
       else
          p = uffparam("N_1", 0.656d0, 180.0d0, 3.660d0, 0.069d0, 2.5438d0, 0.000d0, 2.000d0, 6.899d0, 5.880d0, 0.715d0)
       end if
    case (8) ! O
       if (arom) then
          p = uffparam("O_R", 0.680d0, 110.3d0, 3.500d0, 0.060d0, 2.2998d0, 0.000d0, 2.000d0, 8.741d0, 6.682d0, 0.669d0)
       else if (ncon >= 2) then
          p = uffparam("O_3", 0.658d0, 104.51d0, 3.500d0, 0.060d0, 2.2998d0, 0.018d0, 2.000d0, 8.741d0, 6.682d0, 0.669d0)
       else
          p = uffparam("O_2", 0.634d0, 120.0d0, 3.500d0, 0.060d0, 2.2998d0, 0.000d0, 2.000d0, 8.741d0, 6.682d0, 0.669d0)
       end if
    case (15) ! P
       if (ncon >= 4) then
          p = uffparam("P_3+5", 1.056d0, 109.4712d0, 4.147d0, 0.305d0, 2.8627d0, 2.400d0, 1.250d0, 5.463d0, 4.000d0, 1.102d0)
       else
          p = uffparam("P_3+3", 1.101d0, 93.8d0, 4.147d0, 0.305d0, 2.8627d0, 2.400d0, 1.250d0, 5.463d0, 4.000d0, 1.102d0)
       end if
    case (16) ! S
       if (arom) then
          p = uffparam("S_R", 1.077d0, 92.2d0, 4.035d0, 0.274d0, 2.7032d0, 0.000d0, 1.250d0, 6.928d0, 4.486d0, 1.047d0)
       else if (ncon >= 4) then
          p = uffparam("S_3+6", 1.027d0, 109.4712d0, 4.035d0, 0.274d0, 2.7032d0, 0.484d0, 1.250d0, 6.928d0, 4.486d0, 1.047d0)
       else if (ncon == 3) then
          p = uffparam("S_3+4", 1.049d0, 103.2d0, 4.035d0, 0.274d0, 2.7032d0, 0.484d0, 1.250d0, 6.928d0, 4.486d0, 1.047d0)
       else
          p = uffparam("S_3+2", 1.064d0, 92.1d0, 4.035d0, 0.274d0, 2.7032d0, 0.484d0, 1.250d0, 6.928d0, 4.486d0, 1.047d0)
       end if
    case (1) ! H
       p = uffparam("H_", 0.354d0, 180.0d0, 2.886d0, 0.044d0, 0.7120d0, 0.000d0, 0.000d0, 4.528d0, 6.9452d0, 0.371d0)
    case (2) ! He
       p = uffparam("He4+4", 0.849d0, 90.0d0, 2.362d0, 0.056d0, 0.0972d0, 0.000d0, 0.000d0, 9.660d0, 14.920d0, 1.300d0)
    case (3) ! Li
       p = uffparam("Li", 1.336d0, 180.0d0, 2.451d0, 0.025d0, 1.0255d0, 0.000d0, 2.000d0, 3.006d0, 2.386d0, 1.557d0)
    case (4) ! Be
       p = uffparam("Be3+2", 1.074d0, 109.4712d0, 2.745d0, 0.085d0, 1.5650d0, 0.000d0, 2.000d0, 4.877d0, 4.443d0, 1.240d0)
    case (9) ! F
       p = uffparam("F_", 0.668d0, 180.0d0, 3.364d0, 0.050d0, 1.7350d0, 0.000d0, 2.000d0, 10.874d0, 7.474d0, 0.706d0)
    case (10) ! Ne
       p = uffparam("Ne4+4", 0.920d0, 90.0d0, 3.243d0, 0.042d0, 0.1944d0, 0.000d0, 2.000d0, 11.040d0, 10.550d0, 1.768d0)
    case (11) ! Na
       p = uffparam("Na", 1.539d0, 180.0d0, 2.983d0, 0.030d0, 1.0809d0, 0.000d0, 1.250d0, 2.843d0, 2.296d0, 2.085d0)
    case (12) ! Mg
       p = uffparam("Mg3+2", 1.421d0, 109.4712d0, 3.021d0, 0.111d0, 1.7866d0, 0.000d0, 1.250d0, 3.951d0, 3.693d0, 1.500d0)
    case (13) ! Al
       p = uffparam("Al3", 1.244d0, 109.4712d0, 4.499d0, 0.505d0, 1.7924d0, 0.000d0, 1.250d0, 3.041d0, 3.590d0, 1.201d0)
    case (14) ! Si
       p = uffparam("Si3", 1.117d0, 109.4712d0, 4.295d0, 0.402d0, 2.3232d0, 1.225d0, 1.250d0, 4.168d0, 3.487d0, 1.176d0)
    case (17) ! Cl
       p = uffparam("Cl", 1.044d0, 180.0d0, 3.947d0, 0.227d0, 2.3484d0, 0.000d0, 1.250d0, 8.564d0, 4.946d0, 0.994d0)
    case (18) ! Ar
       p = uffparam("Ar4+4", 1.032d0, 90.0d0, 3.868d0, 0.185d0, 0.2994d0, 0.000d0, 1.250d0, 9.465d0, 6.355d0, 2.108d0)
    case (19) ! K
       p = uffparam("K_", 1.953d0, 180.0d0, 3.812d0, 0.035d0, 1.1645d0, 0.000d0, 0.700d0, 2.421d0, 1.920d0, 2.586d0)
    case (20) ! Ca
       p = uffparam("Ca6+2", 1.761d0, 90.0d0, 3.399d0, 0.238d0, 2.1414d0, 0.000d0, 0.700d0, 3.231d0, 2.880d0, 2.000d0)
    case (21) ! Sc
       p = uffparam("Sc3+3", 1.513d0, 109.4712d0, 3.295d0, 0.019d0, 2.5924d0, 0.000d0, 0.700d0, 3.395d0, 3.080d0, 1.750d0)
    case (22) ! Ti
       p = uffparam("Ti6+4", 1.412d0, 90.0d0, 3.175d0, 0.017d0, 2.6595d0, 0.000d0, 0.700d0, 3.470d0, 3.380d0, 1.607d0)
    case (23) ! V
       p = uffparam("V_3+5", 1.402d0, 109.4712d0, 3.144d0, 0.016d0, 2.6789d0, 0.000d0, 0.700d0, 3.650d0, 3.410d0, 1.470d0)
    case (24) ! Cr
       p = uffparam("Cr6+3", 1.345d0, 90.0d0, 3.023d0, 0.015d0, 2.4631d0, 0.000d0, 0.700d0, 3.415d0, 3.865d0, 1.402d0)
    case (25) ! Mn
       p = uffparam("Mn6+2", 1.382d0, 90.0d0, 2.961d0, 0.013d0, 2.4301d0, 0.000d0, 0.700d0, 3.325d0, 4.105d0, 1.533d0)
    case (26) ! Fe
       p = uffparam("Fe6+2", 1.335d0, 90.0d0, 2.912d0, 0.013d0, 2.4301d0, 0.000d0, 0.700d0, 3.760d0, 4.140d0, 1.393d0)
    case (27) ! Co
       p = uffparam("Co6+3", 1.241d0, 90.00d0, 2.872d0, 0.014d0, 2.430d0, 0.000d0, 0.700d0, 4.105d0, 4.175d0, 1.406d0)
    case (28) ! Ni
       p = uffparam("Ni4+2", 1.164d0, 90.0d0, 2.834d0, 0.015d0, 2.4301d0, 0.000d0, 0.700d0, 4.465d0, 4.205d0, 1.398d0)
    case (29) ! Cu
       p = uffparam("Cu3+1", 1.302d0, 109.4712d0, 3.495d0, 0.005d0, 1.7565d0, 0.000d0, 0.700d0, 3.729d0, 4.220d0, 1.434d0)
    case (30) ! Zn
       p = uffparam("Zn3+2", 1.193d0, 109.4712d0, 2.763d0, 0.124d0, 1.3084d0, 0.000d0, 0.700d0, 5.106d0, 4.285d0, 1.400d0)
    case (31) ! Ga
       p = uffparam("Ga3+3", 1.260d0, 109.4712d0, 4.383d0, 0.415d0, 1.8206d0, 0.000d0, 0.700d0, 2.999d0, 3.160d0, 1.211d0)
    case (32) ! Ge
       p = uffparam("Ge3", 1.197d0, 109.4712d0, 4.280d0, 0.379d0, 2.7888d0, 0.701d0, 0.700d0, 4.051d0, 3.438d0, 1.189d0)
    case (33) ! As
       p = uffparam("As3+3", 1.211d0, 92.1d0, 4.230d0, 0.309d0, 2.8640d0, 1.500d0, 0.700d0, 5.188d0, 3.809d0, 1.204d0)
    case (34) ! Se
       p = uffparam("Se3+2", 1.190d0, 90.6d0, 4.205d0, 0.291d0, 2.7645d0, 0.335d0, 0.700d0, 6.428d0, 4.131d0, 1.224d0)
    case (35) ! Br
       p = uffparam("Br", 1.192d0, 180.0d0, 4.189d0, 0.251d0, 2.5186d0, 0.000d0, 0.700d0, 7.790d0, 4.425d0, 1.141d0)
    case (36) ! Kr
       p = uffparam("Kr4+4", 1.147d0, 90.0d0, 4.141d0, 0.220d0, 0.4520d0, 0.000d0, 0.700d0, 8.505d0, 5.715d0, 2.270d0)
    case (37) ! Rb
       p = uffparam("Rb", 2.260d0, 180.0d0, 4.114d0, 0.040d0, 1.5922d0, 0.000d0, 0.200d0, 2.331d0, 1.846d0, 2.770d0)
    case (38) ! Sr
       p = uffparam("Sr6+2", 2.052d0, 90.0d0, 3.641d0, 0.235d0, 2.4486d0, 0.000d0, 0.200d0, 3.024d0, 2.440d0, 2.415d0)
    case (39) ! Y
       p = uffparam("Y_3+3", 1.698d0, 109.4712d0, 3.345d0, 0.072d0, 3.2573d0, 0.000d0, 0.200d0, 3.830d0, 2.810d0, 1.998d0)
    case (40) ! Zr
       p = uffparam("Zr3+4", 1.564d0, 109.4712d0, 3.124d0, 0.069d0, 3.6675d0, 0.000d0, 0.200d0, 3.400d0, 3.550d0, 1.758d0)
    case (41) ! Nb
       p = uffparam("Nb3+5", 1.473d0, 109.4712d0, 3.165d0, 0.059d0, 3.6179d0, 0.000d0, 0.200d0, 3.550d0, 3.380d0, 1.603d0)
    case (42) ! Mo
       p = uffparam("Mo6+6", 1.467d0, 90.0000d0, 3.052d0, 0.056d0, 3.4021d0, 0.000d0, 0.200d0, 3.465d0, 3.755d0, 1.530d0)
    case (43) ! Tc
       p = uffparam("Tc6+5", 1.322d0, 90.0000d0, 2.998d0, 0.048d0, 3.4021d0, 0.000d0, 0.200d0, 3.290d0, 3.990d0, 1.500d0)
    case (44) ! Ru
       p = uffparam("Ru6+2", 1.478d0, 90.0000d0, 2.963d0, 0.056d0, 3.4021d0, 0.000d0, 0.200d0, 3.575d0, 4.015d0, 1.500d0)
    case (45) ! Rh
       p = uffparam("Rh6+3", 1.332d0, 90.0000d0, 2.929d0, 0.053d0, 3.5081d0, 0.000d0, 0.200d0, 3.975d0, 4.005d0, 1.509d0)
    case (46) ! Pd
       p = uffparam("Pd4+2", 1.338d0, 90.0000d0, 2.899d0, 0.048d0, 3.2077d0, 0.000d0, 0.200d0, 4.320d0, 4.000d0, 1.544d0)
    case (47) ! Ag
       p = uffparam("Ag1+1", 1.386d0, 180.0000d0, 3.148d0, 0.036d0, 1.9557d0, 0.000d0, 0.200d0, 4.436d0, 3.134d0, 1.622d0)
    case (48) ! Cd
       p = uffparam("Cd3+2", 1.403d0, 109.4712d0, 2.848d0, 0.228d0, 1.6525d0, 0.000d0, 0.200d0, 5.034d0, 3.957d0, 1.600d0)
    case (49) ! In
       p = uffparam("In3+3", 1.459d0, 109.4712d0, 4.463d0, 0.599d0, 2.0704d0, 0.000d0, 0.200d0, 2.997d0, 2.896d0, 1.404d0)
    case (50) ! Sn
       p = uffparam("Sn3", 1.398d0, 109.4712d0, 4.392d0, 0.567d0, 2.9608d0, 0.199d0, 0.200d0, 3.987d0, 3.124d0, 1.354d0)
    case (51) ! Sb
       p = uffparam("Sb3+3", 1.407d0, 91.6000d0, 4.420d0, 0.449d0, 2.7042d0, 1.100d0, 0.200d0, 4.899d0, 3.342d0, 1.404d0)
    case (52) ! Te
       p = uffparam("Te3+2", 1.386d0, 90.2500d0, 4.470d0, 0.398d0, 2.8821d0, 0.300d0, 0.200d0, 5.816d0, 3.526d0, 1.380d0)
    case (53) ! I
       p = uffparam("I_", 1.382d0, 180.0000d0, 4.500d0, 0.339d0, 2.6537d0, 0.000d0, 0.200d0, 6.822d0, 3.762d0, 1.333d0)
    case (54) ! Xe
       p = uffparam("Xe4+4", 1.267d0, 90.0000d0, 4.404d0, 0.332d0, 0.5560d0, 0.000d0, 0.200d0, 7.595d0, 4.975d0, 2.459d0)
    case (55) ! Cs
       p = uffparam("Cs", 2.570d0, 180.0000d0, 4.517d0, 0.045d0, 1.5728d0, 0.000d0, 0.100d0, 2.183d0, 1.711d0, 2.984d0)
    case (56) ! Ba
       p = uffparam("Ba6+2", 2.277d0, 90.0000d0, 3.703d0, 0.364d0, 2.7266d0, 0.000d0, 0.100d0, 2.814d0, 2.396d0, 2.442d0)
    case (57) ! La
       p = uffparam("La3+3", 1.943d0, 109.4712d0, 3.522d0, 0.017d0, 3.3049d0, 0.000d0, 0.100d0, 2.8355d0, 2.7415d0, 2.071d0)
    case (58) ! Ce
       p = uffparam("Ce6+3", 1.841d0, 90.0000d0, 3.556d0, 0.013d0, 3.3049d0, 0.000d0, 0.100d0, 2.774d0, 2.692d0, 1.925d0)
    case (59) ! Pr
       p = uffparam("Pr6+3", 1.823d0, 90.0000d0, 3.606d0, 0.010d0, 3.3049d0, 0.000d0, 0.100d0, 2.858d0, 2.564d0, 2.007d0)
    case (60) ! Nd
       p = uffparam("Nd6+3", 1.816d0, 90.0000d0, 3.575d0, 0.010d0, 3.3049d0, 0.000d0, 0.100d0, 2.8685d0, 2.6205d0, 2.007d0)
    case (61) ! Pm
       p = uffparam("Pm6+3", 1.801d0, 90.0000d0, 3.547d0, 0.009d0, 3.3049d0, 0.000d0, 0.100d0, 2.881d0, 2.673d0, 2.000d0)
    case (62) ! Sm
       p = uffparam("Sm6+3", 1.780d0, 90.0000d0, 3.520d0, 0.008d0, 3.3049d0, 0.000d0, 0.100d0, 2.9115d0, 2.7195d0, 1.978d0)
    case (63) ! Eu
       p = uffparam("Eu6+3", 1.771d0, 90.0000d0, 3.493d0, 0.008d0, 3.3049d0, 0.000d0, 0.100d0, 2.8785d0, 2.7875d0, 2.227d0)
    case (64) ! Gd
       p = uffparam("Gd6+3", 1.735d0, 90.0000d0, 3.368d0, 0.009d0, 3.3049d0, 0.000d0, 0.100d0, 3.1665d0, 2.9745d0, 1.968d0)
    case (65) ! Tb
       p = uffparam("Tb6+3", 1.732d0, 90.0000d0, 3.451d0, 0.007d0, 3.3049d0, 0.000d0, 0.100d0, 3.018d0, 2.834d0, 1.954d0)
    case (66) ! Dy
       p = uffparam("Dy6+3", 1.710d0, 90.0000d0, 3.428d0, 0.007d0, 3.3049d0, 0.000d0, 0.100d0, 3.0555d0, 2.8715d0, 1.934d0)
    case (67) ! Ho
       p = uffparam("Ho6+3", 1.696d0, 90.0000d0, 3.409d0, 0.007d0, 3.4157d0, 0.000d0, 0.100d0, 3.127d0, 2.891d0, 1.925d0)
    case (68) ! Er
       p = uffparam("Er6+3", 1.673d0, 90.0000d0, 3.391d0, 0.007d0, 3.3049d0, 0.000d0, 0.100d0, 3.1865d0, 2.9145d0, 1.915d0)
    case (69) ! Tm
       p = uffparam("Tm6+3", 1.660d0, 90.0000d0, 3.374d0, 0.006d0, 3.3049d0, 0.000d0, 0.100d0, 3.2514d0, 2.9329d0, 2.000d0)
    case (70) ! Yb
       p = uffparam("Yb6+3", 1.637d0, 90.0000d0, 3.355d0, 0.228d0, 2.6177d0, 0.000d0, 0.100d0, 3.2889d0, 2.9650d0, 2.158d0)
    case (71) ! Lu
       p = uffparam("Lu6+3", 1.671d0, 90.0000d0, 3.640d0, 0.041d0, 3.2709d0, 0.000d0, 0.100d0, 2.9629d0, 2.4629d0, 1.896d0)
    case (72) ! Hf
       p = uffparam("Hf3+4", 1.611d0, 109.4712d0, 3.141d0, 0.072d0, 3.9212d0, 0.000d0, 0.100d0, 3.7000d0, 3.400d0, 1.759d0)
    case (73) ! Ta
       p = uffparam("Ta3+5", 1.511d0, 109.4712d0, 3.170d0, 0.081d0, 4.0748d0, 0.000d0, 0.100d0, 5.1000d0, 2.850d0, 1.605d0)
    case (74) ! W
       p = uffparam("W_6+6", 1.392d0, 90.0000d0, 3.069d0, 0.067d0, 3.6937d0, 0.000d0, 0.100d0, 4.6300d0, 3.310d0, 1.538d0)
    case (75) ! Re
       p = uffparam("Re6+5", 1.372d0, 90.0000d0, 2.954d0, 0.066d0, 3.6937d0, 0.000d0, 0.100d0, 3.9600d0, 3.920d0, 1.600d0)
    case (76) ! Os
       p = uffparam("Os6+6", 1.372d0, 90.0000d0, 3.120d0, 0.037d0, 3.6937d0, 0.000d0, 0.100d0, 5.1400d0, 3.630d0, 1.700d0)
    case (77) ! Ir
       p = uffparam("Ir6+3", 1.371d0, 90.0000d0, 2.840d0, 0.073d0, 3.7307d0, 0.000d0, 0.100d0, 5.0000d0, 4.000d0, 1.866d0)
    case (78) ! Pt
       p = uffparam("Pt4+2", 1.364d0, 90.0000d0, 2.754d0, 0.080d0, 3.3817d0, 0.000d0, 0.100d0, 4.7900d0, 4.430d0, 1.557d0)
    case (79) ! Au
       p = uffparam("Au4+3", 1.262d0, 90.0000d0, 3.293d0, 0.039d0, 2.6255d0, 0.000d0, 0.100d0, 4.8940d0, 2.586d0, 1.618d0)
    case (80) ! Hg
       p = uffparam("Hg1+2", 1.340d0, 180.0000d0, 2.705d0, 0.385d0, 1.7497d0, 0.000d0, 0.100d0, 6.2700d0, 4.160d0, 1.600d0)
    case (81) ! Tl
       p = uffparam("Tl3+3", 1.518d0, 120.0000d0, 4.347d0, 0.680d0, 2.0685d0, 0.000d0, 0.100d0, 3.2000d0, 2.900d0, 1.530d0)
    case (82) ! Pb
       p = uffparam("Pb3", 1.459d0, 109.4712d0, 4.297d0, 0.663d0, 2.8461d0, 0.100d0, 0.100d0, 3.9000d0, 3.530d0, 1.444d0)
    case (83) ! Bi
       p = uffparam("Bi3+3", 1.512d0, 90.0000d0, 4.370d0, 0.518d0, 2.4700d0, 1.000d0, 0.100d0, 4.6900d0, 3.740d0, 1.514d0)
    case (84) ! Po
       p = uffparam("Po3+2", 1.500d0, 90.0000d0, 4.709d0, 0.325d0, 2.3329d0, 0.300d0, 0.100d0, 4.2100d0, 4.210d0, 1.480d0)
    case (85) ! At
       p = uffparam("At", 1.545d0, 180.0000d0, 4.750d0, 0.284d0, 2.2357d0, 0.000d0, 0.100d0, 4.7500d0, 4.750d0, 1.470d0)
    case (86) ! Rn
       p = uffparam("Rn4+4", 1.420d0, 90.0000d0, 4.765d0, 0.248d0, 0.5832d0, 0.000d0, 0.100d0, 5.3700d0, 5.370d0, 2.200d0)
    case (87) ! Fr
       p = uffparam("Fr", 2.880d0, 180.0000d0, 4.900d0, 0.050d0, 1.8469d0, 0.000d0, 0.000d0, 2.0000d0, 2.000d0, 2.300d0)
    case (88) ! Ra
       p = uffparam("Ra6+2", 2.512d0, 90.0000d0, 3.677d0, 0.404d0, 2.9161d0, 0.000d0, 0.000d0, 2.8430d0, 2.434d0, 2.200d0)
    case (89) ! Ac
       p = uffparam("Ac6+3", 1.983d0, 90.0000d0, 3.478d0, 0.033d0, 3.8882d0, 0.000d0, 0.000d0, 2.8350d0, 2.835d0, 2.108d0)
    case (90) ! Th
       p = uffparam("Th6+4", 1.721d0, 90.0000d0, 3.396d0, 0.026d0, 4.2021d0, 0.000d0, 0.000d0, 3.1750d0, 2.905d0, 2.018d0)
    case (91) ! Pa
       p = uffparam("Pa6+4", 1.711d0, 90.0000d0, 3.424d0, 0.022d0, 3.8882d0, 0.000d0, 0.000d0, 2.9850d0, 2.905d0, 1.800d0)
    case (92) ! U
       p = uffparam("U_6+4", 1.684d0, 90.0000d0, 3.395d0, 0.022d0, 3.8882d0, 0.000d0, 0.000d0, 3.3410d0, 2.853d0, 1.713d0)
    case (93) ! Np
       p = uffparam("Np6+4", 1.666d0, 90.0000d0, 3.424d0, 0.019d0, 3.8882d0, 0.000d0, 0.000d0, 3.5490d0, 2.717d0, 1.800d0)
    case (94) ! Pu
       p = uffparam("Pu6+4", 1.657d0, 90.0000d0, 3.424d0, 0.016d0, 3.8882d0, 0.000d0, 0.000d0, 3.2430d0, 2.819d0, 1.840d0)
    case (95) ! Am
       p = uffparam("Am6+4", 1.660d0, 90.0000d0, 3.381d0, 0.014d0, 3.8882d0, 0.000d0, 0.000d0, 2.9895d0, 3.0035d0, 1.942d0)
    case (96) ! Cm
       p = uffparam("Cm6+3", 1.801d0, 90.0000d0, 3.326d0, 0.013d0, 3.8882d0, 0.000d0, 0.000d0, 2.8315d0, 3.1895d0, 1.900d0)
    case (97) ! Bk
       p = uffparam("Bk6+3", 1.761d0, 90.0000d0, 3.339d0, 0.013d0, 3.8882d0, 0.000d0, 0.000d0, 3.1935d0, 3.0355d0, 1.900d0)
    case (98) ! Cf
       p = uffparam("Cf6+3", 1.750d0, 90.0000d0, 3.313d0, 0.013d0, 3.8882d0, 0.000d0, 0.000d0, 3.1970d0, 3.101d0, 1.900d0)
    case (99) ! Es
       p = uffparam("Es6+3", 1.724d0, 90.0000d0, 3.299d0, 0.012d0, 3.8882d0, 0.000d0, 0.000d0, 3.3330d0, 3.089d0, 1.900d0)
    case (100) ! Fm
       p = uffparam("Fm6+3", 1.712d0, 90.0000d0, 3.286d0, 0.012d0, 3.8882d0, 0.000d0, 0.000d0, 3.4000d0, 3.100d0, 1.900d0)
    case (101) ! Md
       p = uffparam("Md6+3", 1.689d0, 90.0000d0, 3.274d0, 0.011d0, 3.8882d0, 0.000d0, 0.000d0, 3.4700d0, 3.110d0, 1.900d0)
    case (102) ! No
       p = uffparam("No6+3", 1.679d0, 90.0000d0, 3.248d0, 0.011d0, 3.8882d0, 0.000d0, 0.000d0, 3.4750d0, 3.175d0, 1.900d0)
    case (103) ! Lw
       p = uffparam("Lw6+3", 1.698d0, 90.0000d0, 3.236d0, 0.011d0, 3.8882d0, 0.000d0, 0.000d0, 3.5000d0, 3.200d0, 1.900d0)
    end select

  end subroutine uff_atom_params

  !xx! QEq (Rappe-Goddard) electrostatics -- ported from GULP's gamma.F90

  !> Binomial expansion coefficient used by the two-centre integral qeq_css
  !> (GULP coeffs()). When the product of the two orbitals' radial polynomials is
  !> rewritten in prolate-spheroidal coordinates, the term of "order" k carries
  !>   C_k(na,nb) = sum_i (na choose i) (nb choose l-i) (-1)^(l-i),  l = na+nb-k,
  !> the sum running over the i that keep both binomials in range.
  pure function qeq_coeffs(na,nb,k) result(co)
    use tools_math, only: nchoosek
    integer, intent(in) :: na, nb, k
    real*8 :: co
    integer :: l, ie, je, ia, il, i, j
    co = 0d0
    l = na + nb - k
    ie = min(l,na) + 1
    je = min(l,nb)
    ia = l - je + 1
    do il = ia, ie
       i = il - 1
       j = l - i
       co = co + nchoosek(na,i)*nchoosek(nb,j)*(-1d0)**j
    end do
  end function qeq_coeffs

  !> A-type auxiliary integrals of Roothaan/Mulliken (GULP caintgs),
  !>   A_n(x) = int_1^infinity t^n exp(-x t) dt,  n = 0..k,  x > 0,
  !> by the upward recursion A_n = (n A_{n-1} + exp(-x))/x with A_0 = exp(-x)/x.
  !> Returns a(1:k+1) with a(i) = A_{i-1}(x).
  pure subroutine qeq_caintgs(x,k,a)
    real*8, intent(in) :: x
    integer, intent(in) :: k
    real*8, intent(out) :: a(:)
    real*8 :: const, rx
    integer :: i
    const = exp(-x)
    rx = 1d0/x
    a(1) = const*rx
    do i = 1, k
       a(i+1) = (a(i)*dble(i) + const)*rx
    end do
  end subroutine qeq_caintgs

  !> B-type auxiliary integrals of Roothaan/Mulliken (GULP cbintgs),
  !>   B_n(x) = int_{-1}^{1} t^n exp(-x t) dt,  n = 0..k.
  !> For large |x| the upward recursion B_n = (n B_{n-1} + (-1)^n exp(x) -
  !> exp(-x))/x with B_0 = (exp(x)-exp(-x))/x is used; for small |x| it is
  !> numerically unstable, so a truncated Taylor series is used instead, and for
  !> x = 0 the exact B_n = (1-(-1)^(n+1))/(n+1). Returns b(1:k+1), b(i)=B_{i-1}(x).
  pure subroutine qeq_cbintgs(x,k,b)
    real*8, intent(in) :: x
    integer, intent(in) :: k
    real*8, intent(out) :: b(:)
    integer :: i, last, m
    real*8 :: absx, expx, expmx, rx, y, ytrm
    absx = abs(x)
    last = 0
    if (absx > 3d0) then
       last = -1
    else if (absx > 2d0) then
       if (k <= 10) then
          last = -1
       else
          last = 15
       end if
    else if (absx > 1d0) then
       if (k <= 7) then
          last = -1
       else
          last = 12
       end if
    else if (absx > 0.5d0) then
       if (k <= 5) then
          last = -1
       else
          last = 7
       end if
    else if (absx >= 1d-8) then
       last = 6
    else
       last = -2
    end if
    if (last == -1) then
       expx = exp(x)
       expmx = 1d0/expx
       rx = 1d0/x
       b(1) = (expx - expmx)*rx
       do i = 1, k
          b(i+1) = (dble(i)*b(i) + (-1d0)**i*expx - expmx)*rx
       end do
    else if (last == -2) then
       do i = 0, k
          b(i+1) = (1d0 - (-1d0)**(i+1))/dble(i+1)
       end do
    else
       do i = 0, k
          y = 0d0
          do m = 0, last
             ytrm = (-x)**(m-1)*(1d0 - (-1d0)**(m+i+1))/(fact(m)*dble(m+i+1))
             y = y + ytrm*(-x)
          end do
          b(i+1) = y
       end do
    end if
  end subroutine qeq_cbintgs

  !> Reduced two-centre integral over s-type Slater orbitals of principal
  !> quantum numbers n1 (centre A) and n2 (centre B), evaluated with the A/B
  !> auxiliary functions in prolate-spheroidal coordinates (Roothaan/Mulliken;
  !> GULP css):
  !>   s = 1/2 sum_{m=0}^{n1+n2} C_m(n1,n2) A_m(p) B_{n1+n2-m}(pt),
  !> with p = (alpha+beta)/2, pt = (alpha-beta)/2, alpha = 2 za r, beta = 2 zb r,
  !> A/B from qeq_caintgs/qeq_cbintgs and C from qeq_coeffs. Returns the value s
  !> and its derivatives with respect to r (deriv) and to the exponent za (dz1).
  subroutine qeq_css(nn1,nn2,alpha,beta,r,s,deriv,dz1)
    integer, intent(in) :: nn1, nn2
    real*8, intent(in) :: alpha, beta, r
    real*8, intent(out) :: s, deriv, dz1
    integer :: n1, n2, k, ulim, i1, nni1
    ! a/b hold the auxiliary integrals A_n/B_n filled up to index n1+n2+4; with the
    ! largest principal quantum number qeq_pqn returns (7, for Z>86) n1=n2=13, so
    ! the top index is exactly 30. Size 30 is the tight upper bound for real elements.
    real*8 :: p, pt, dpdr, dptdr, dpdz1, dptdz1, x, a(30), b(30), coff, dari, dbrn, da1i, db1n
    n1 = nn1
    n2 = nn2
    p = (alpha + beta)*0.5d0
    pt = (alpha - beta)*0.5d0
    dpdr = 0.5d0*(alpha + beta)/r
    dptdr = 0.5d0*(alpha - beta)/r
    dpdz1 = r
    dptdz1 = r
    if (n2 < n1) then
       k = n1
       n1 = n2
       n2 = k
       pt = -pt
       dptdr = -dptdr
       dptdz1 = -dptdz1
    end if
    s = 0d0
    deriv = 0d0
    dz1 = 0d0
    if (p > 86d0 .or. pt > 86d0) return
    call qeq_caintgs(p,n1+n2+3,a)
    call qeq_cbintgs(pt,n1+n2+3,b)
    ulim = n1 + n2 + 1
    x = 0d0
    do i1 = 1, ulim
       nni1 = n1 + n2 - i1 + 2
       coff = qeq_coeffs(n1,n2,i1-1)
       dari = -a(i1+1)*dpdr
       dbrn = -b(nni1+1)*dptdr
       da1i = -a(i1+1)*dpdz1
       db1n = -b(nni1+1)*dptdz1
       x = x + coff*a(i1)*b(nni1)
       deriv = deriv + coff*(dari*b(nni1) + a(i1)*dbrn)
       dz1 = dz1 + coff*(da1i*b(nni1) + a(i1)*db1n)
    end do
    s = x*0.5d0
    deriv = deriv*0.5d0
    dz1 = dz1*0.5d0
  end subroutine qeq_css

  !> QEq short-range kernel: J_ij(r) - 1/r and its r-derivative (atomic units,
  !> r in bohr, za/zb orbital exponents in 1/bohr), tapered to 0 at qeq_rqeq.
  !> Ported from GULP gammas() with nqeqjab=1 (energy + first derivative only).
  subroutine qeq_slater(na,nb,za,zb,r,gam,dgam,dza)
    integer, intent(in) :: na, nb
    real*8, intent(in) :: za, zb, r
    real*8, intent(out) :: gam, dgam, dza
    integer :: na2, nb2, i
    real*8 :: z2ra, z2rb, halfr, d2rtrm, drtrm, rtrm, ss, deriv, dz1
    real*8 :: ztrm, ctrm, trm1, trm2, trm3, tpfn, dtpfn, xt, rrd

    gam = 0d0
    dgam = 0d0
    dza = 0d0
    if (r < 1d-10 .or. r > qeq_rqeq) return
    z2ra = 2d0*za*r
    z2rb = 2d0*zb*r
    if (z2ra > 80d0 .or. z2rb > 80d0) return
    na2 = 2*na
    nb2 = 2*nb
    halfr = 0.5d0*r
    d2rtrm = halfr**(na2-2)
    drtrm = d2rtrm*halfr
    rtrm = drtrm*halfr
    ! first term (zb = 0)
    call qeq_css(na2-1,0,z2ra,0d0,r,ss,deriv,dz1)
    gam = rtrm*ss
    dgam = rtrm*deriv + dble(na)*drtrm*ss
    dza = rtrm*dz1
    ! sum over 2*nb
    rtrm = drtrm
    drtrm = d2rtrm
    ztrm = 0.5d0/(zb*dble(nb2))
    do i = nb2, 1, -1
       rtrm = rtrm*halfr
       drtrm = drtrm*halfr
       ztrm = ztrm*2d0*zb
       ctrm = ztrm/fact(nb2-i)
       call qeq_css(na2-1,nb2-i,z2ra,z2rb,r,ss,deriv,dz1)
       trm1 = dble(i)*ctrm
       trm2 = trm1*rtrm
       trm3 = trm1*dble(na2+nb2-i)*drtrm
       gam = gam - trm2*ss
       dgam = dgam - trm2*deriv - 0.5d0*trm3*ss
       dza = dza - trm2*dz1
    end do
    trm3 = ((2d0*za)**(na2+1))/fact(na2)
    gam = gam*trm3
    dgam = dgam*trm3
    ! dJ/dza; the extra term is d/dza of the (2 za)^(na2+1) prefactor times gam
    dza = dza*trm3 + dble(na2+1)/za*gam
    ! difference from the bare Coulomb 1/r (no za dependence)
    gam = gam - 1d0/r
    dgam = dgam + 1d0/(r*r)
    ! taper (1-x)^3 (1 + 3x + 6x^2) between qeq_rqeqtaper and qeq_rqeq
    if (r > qeq_rqeqtaper) then
       rrd = 1d0/(qeq_rqeq - qeq_rqeqtaper)
       xt = (r - qeq_rqeqtaper)*rrd
       tpfn = (1d0-xt)**3*(1d0 + 3d0*xt + 6d0*xt*xt)
       dtpfn = (-3d0*(1d0-xt)**2*(1d0 + 3d0*xt + 6d0*xt*xt) + (1d0-xt)**3*(3d0 + 12d0*xt))*rrd
       dgam = dgam*tpfn + gam*dtpfn
       gam = gam*tpfn
       dza = dza*tpfn
    end if
  end subroutine qeq_slater

  !> Principal quantum number used by QEq (GULP convention).
  pure function qeq_pqn(z) result(npqn)
    integer, intent(in) :: z
    integer :: npqn
    if (z <= 2) then
       npqn = 1
    else if (z <= 10) then
       npqn = 2
    else if (z <= 18) then
       npqn = 3
    else if (z <= 36) then
       npqn = 4
    else if (z <= 54) then
       npqn = 5
    else if (z <= 86) then
       npqn = 6
    else
       npqn = 7
    end if
  end function qeq_pqn

  !> Add the Slater short-range terms to the QEq matrix amat. honly selects the
  !> subset of pairs: .false. adds only pairs whose exponents are charge-
  !> independent (neither atom hydrogen), .true. adds only hydrogen-involving
  !> pairs (using the current exponents zeta and, for H, the q(dJ/dzeta) coupling
  !> in the coefficient of q_i). Crystal uses the cached neighbour list + Ewald
  !> (1/r already in amat); molecule uses the direct full J over all pairs.
  subroutine qeq_add_slater(cl,c,npqn,zat,zeta,honly,q,amat)
    use crystalmod, only: crystal
    class(calculator), intent(in) :: cl
    class(crystal), intent(in) :: c
    integer, intent(in) :: npqn(:), zat(:)
    real*8, intent(in) :: zeta(:), q(:)
    logical, intent(in) :: honly
    real*8, intent(inout) :: amat(:,:)

    integer :: n, i, j, p
    real*8 :: d(3), r, gam, dgam, dza
    logical :: ph

    n = cl%nat
    if (c%ismolecule) then
       do i = 1, n
          do j = 1, n
             if (i == j) cycle
             ph = (zat(i) == 1 .or. zat(j) == 1)
             if (ph .neqv. honly) cycle
             d = c%atcel(i)%r - c%atcel(j)%r
             r = norm2(d)
             if (r < 1d-10) cycle
             call qeq_slater(npqn(i),npqn(j),zeta(i),zeta(j),r,gam,dgam,dza)
             amat(i,j) = amat(i,j) + gam + 1d0/r
             if (honly .and. zat(i) == 1) amat(j,i) = amat(j,i) + q(i)*dza
          end do
       end do
    else
       do i = 1, n
          do p = cl%qnbptr(i), cl%qnbptr(i+1)-1
             j = cl%qnbj(p)
             ph = (zat(i) == 1 .or. zat(j) == 1)
             if (ph .neqv. honly) cycle
             d = c%atcel(i)%r - c%atcel(j)%r - c%x2c(real(cl%qnblv(:,p),8))
             r = norm2(d)
             if (r < 1d-10 .or. r > qeq_rqeq) cycle
             call qeq_slater(npqn(i),npqn(j),zeta(i),zeta(j),r,gam,dgam,dza)
             amat(i,j) = amat(i,j) + gam
             if (honly .and. zat(i) == 1) amat(j,i) = amat(j,i) + q(i)*dza
          end do
       end do
    end if
  end subroutine qeq_add_slater

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
       if (lv(1) /= 0) then
          ok = lv(1) > 0
       else if (lv(2) /= 0) then
          ok = lv(2) > 0
       else
          ok = lv(3) > 0
       end if
    end if
  end function bond_canonical

  !> Accumulate a virial contribution: the outer product of a vector v and a
  !> force f, vir(p,q) += v(p)*f(q). Called two equivalent ways: for a pair with
  !> the separation vector d = ri - rj and the pair force on i, or per site with
  !> an absolute (imaged) site position and the force on that site (summing a
  !> cluster's sites is translationally invariant because its forces sum to zero).
  pure subroutine virial(vir,v,f)
    real*8, intent(inout) :: vir(3,3)
    real*8, intent(in) :: v(3), f(3)
    integer :: p, q
    do q = 1, 3
       do p = 1, 3
          vir(p,q) = vir(p,q) + v(p)*f(q)
       end do
    end do
  end subroutine virial

  !xx! TIP4P water model (Jorgensen et al., J. Chem. Phys. 79 (1983) 926)

  !> Identify the water molecules for the TIP4P model from the connectivity
  !> graph and store them on the calculator. Molecules only; every atom must
  !> belong to a water (an O bonded to exactly two H). errmsg is empty on
  !> success and holds the error message on failure.
  subroutine tip4p_setup(cl,c,errmsg)
    use crystalmod, only: crystal
    use types, only: neighstar
    use tools_io, only: string
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out) :: errmsg

    type(neighstar), allocatable :: ns(:)
    logical, allocatable :: used(:)
    integer :: i, ih1, ih2

    errmsg = ""
    cl%do_elec = .false. ! the TIP4P charges are fixed; no QEq
    if (.not.c%ismolecule) then
       errmsg = "the TIP4P model is only available for molecules"
       return
    end if

    ! connectivity graph (use the crystal's own if available)
    if (allocated(c%nstar)) then
       ns = c%nstar
    else
       call c%find_asterisms(ns)
    end if

    ! find the waters: an O bonded to exactly two H, each used only once
    cl%nwat = 0
    if (allocated(cl%iwat)) deallocate(cl%iwat)
    allocate(cl%iwat(3,max(cl%nat/3,1)),used(cl%nat))
    used = .false.
    do i = 1, cl%nat
       if (c%spc(c%atcel(i)%is)%z /= 8) cycle
       if (ns(i)%ncon /= 2) then
          errmsg = "TIP4P requires an all-water system (O atom " // string(i) // &
             " does not have exactly two bonded atoms)"
          return
       end if
       ih1 = ns(i)%idcon(1)
       ih2 = ns(i)%idcon(2)
       if (c%spc(c%atcel(ih1)%is)%z /= 1 .or. c%spc(c%atcel(ih2)%is)%z /= 1 .or.&
          ih1 == ih2) then
          errmsg = "TIP4P requires an all-water system (O atom " // string(i) // &
             " is not bonded to two distinct H atoms)"
          return
       end if
       if (used(ih1) .or. used(ih2)) then
          errmsg = "TIP4P requires an all-water system (O atom " // string(i) // &
             " shares an H atom with another water)"
          return
       end if
       used(i) = .true.
       used(ih1) = .true.
       used(ih2) = .true.
       cl%nwat = cl%nwat + 1
       cl%iwat(1,cl%nwat) = i
       cl%iwat(2,cl%nwat) = ih1
       cl%iwat(3,cl%nwat) = ih2
    end do
    if (.not.all(used)) then
       errmsg = "TIP4P requires an all-water system (found " // string(cl%nwat) // &
          " waters for " // string(cl%nat) // " atoms)"
       return
    end if

  end subroutine tip4p_setup

  !> Evaluate the TIP4P energy, gradient, and (optionally) stress
  !> (always zero: molecules only), in atomic units. The monomers are
  !> held together by harmonic restraints on the OH bonds and the HOH
  !> angle, which vanish at the reference geometry (so rigid-geometry
  !> cluster energies match the standard TIP4P values).
  subroutine tip4p_evaluate(cl,c,ene,grad,stress)
    use crystalmod, only: crystal
    use param, only: rad
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    real*8, intent(out) :: ene
    real*8, intent(out) :: grad(:,:)
    real*8, intent(out), optional :: stress(3,3)

    integer :: w1, w2, k1, k2, io, jo, ih, ih1, ih2, k
    real*8 :: am, wo, theta0
    real*8 :: d(3), r, dr, gmag, g(3), r6, r12
    real*8 :: a(3), b(3), na, nb, u, st, theta, edu
    real*8 :: dudi(3), dudj(3), dudk(3)
    real*8 :: q(3)
    real*8, allocatable :: rs(:,:,:), gs(:,:,:)

    ene = 0d0
    grad = 0d0
    if (present(stress)) stress = 0d0

    ! M-site linear weight from the reference geometry
    am = tip4p_dom / (2d0 * tip4p_doh * cos(0.5d0 * tip4p_hoh * rad))
    wo = 1d0 - 2d0*am
    theta0 = tip4p_hoh * rad

    ! charge site positions (H1, H2, M) and charges
    q = (/ tip4p_qh, tip4p_qh, -2d0*tip4p_qh /)
    allocate(rs(3,3,cl%nwat),gs(3,3,cl%nwat))
    do w1 = 1, cl%nwat
       io = cl%iwat(1,w1)
       ih1 = cl%iwat(2,w1)
       ih2 = cl%iwat(3,w1)
       rs(:,1,w1) = c%atcel(ih1)%r
       rs(:,2,w1) = c%atcel(ih2)%r
       rs(:,3,w1) = wo * c%atcel(io)%r + am * (c%atcel(ih1)%r + c%atcel(ih2)%r)
    end do
    gs = 0d0

    ! intramolecular restraints: harmonic OH bonds and HOH angle (zero at the
    ! reference geometry; these keep the monomers together during MD/relaxation)
    do w1 = 1, cl%nwat
       io = cl%iwat(1,w1)
       do k = 2, 3
          ih = cl%iwat(k,w1)
          d = c%atcel(io)%r - c%atcel(ih)%r
          r = norm2(d)
          if (r < 1d-10) cycle
          dr = r - tip4p_doh
          ene = ene + 0.5d0 * tip4p_kbond * dr*dr
          g = tip4p_kbond * dr / r * d
          grad(:,io) = grad(:,io) + g
          grad(:,ih) = grad(:,ih) - g
       end do
       ih1 = cl%iwat(2,w1)
       ih2 = cl%iwat(3,w1)
       a = c%atcel(ih1)%r - c%atcel(io)%r
       b = c%atcel(ih2)%r - c%atcel(io)%r
       na = norm2(a)
       nb = norm2(b)
       if (na < 1d-10 .or. nb < 1d-10) cycle
       u = dot_product(a,b)/(na*nb)
       u = max(min(u,1d0),-1d0)
       theta = acos(u)
       st = sqrt(max(1d0 - u*u,1d-12))
       ene = ene + 0.5d0 * tip4p_kang * (theta - theta0)**2
       ! dE/du = -kang*(theta-theta0)/sin(theta)
       edu = -tip4p_kang * (theta - theta0) / st
       dudi = b/(na*nb) - (u/(na*na))*a
       dudk = a/(na*nb) - (u/(nb*nb))*b
       dudj = -(dudi+dudk)
       grad(:,ih1) = grad(:,ih1) + edu*dudi
       grad(:,ih2) = grad(:,ih2) + edu*dudk
       grad(:,io) = grad(:,io) + edu*dudj
    end do

    ! intermolecular: LJ between O atoms, Coulomb between charge sites
    do w1 = 1, cl%nwat
       io = cl%iwat(1,w1)
       do w2 = w1+1, cl%nwat
          ! LJ between the O atoms
          jo = cl%iwat(1,w2)
          d = c%atcel(io)%r - c%atcel(jo)%r
          r = norm2(d)
          if (r >= 1d-10) then
             r6 = r**6
             r12 = r6*r6
             ene = ene + tip4p_lja/r12 - tip4p_ljc/r6
             ! dE/dr = (-12 A/r^12 + 6 C/r^6)/r
             gmag = (-12d0*tip4p_lja/r12 + 6d0*tip4p_ljc/r6)/(r*r)
             g = gmag * d
             grad(:,io) = grad(:,io) + g
             grad(:,jo) = grad(:,jo) - g
          end if
          ! Coulomb between the charge sites
          do k1 = 1, 3
             do k2 = 1, 3
                d = rs(:,k1,w1) - rs(:,k2,w2)
                r = norm2(d)
                if (r < 1d-10) cycle
                ene = ene + q(k1)*q(k2)/r
                g = -q(k1)*q(k2)/(r*r*r) * d
                gs(:,k1,w1) = gs(:,k1,w1) + g
                gs(:,k2,w2) = gs(:,k2,w2) - g
             end do
          end do
       end do
    end do

    ! redistribute the site gradients onto the atoms: H sites are atoms; the
    ! M site spreads onto O/H1/H2 with the virtual-site weights
    do w1 = 1, cl%nwat
       io = cl%iwat(1,w1)
       ih1 = cl%iwat(2,w1)
       ih2 = cl%iwat(3,w1)
       grad(:,ih1) = grad(:,ih1) + gs(:,1,w1) + am*gs(:,3,w1)
       grad(:,ih2) = grad(:,ih2) + gs(:,2,w1) + am*gs(:,3,w1)
       grad(:,io) = grad(:,io) + wo*gs(:,3,w1)
    end do

  end subroutine tip4p_evaluate

  ! ---- tblite (GFN2/GFN1-xTB) backend ----

  module subroutine calc_init_tblite(cl,c,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out) :: errmsg

#ifdef HAVE_TBLITE
    type(c_ptr) :: err
    integer(c_int), allocatable, target :: numbers(:)
    real(c_double), allocatable, target :: coords(:,:)
    real(c_double), target :: lattice(3,3)
    logical(c_bool), target :: periodic(3)
    integer :: i

    errmsg = ""
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
    ! silence tblite's stdout printout (SCC iterations, summary, timings)
    if (c_associated(cl%tb_ctx)) call c_tblite_set_context_verbosity(cl%tb_ctx,0_c_int)

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
       errmsg = "tblite calculator initialization failed"
    end if
#else
    errmsg = "critic2 was compiled without tblite support"
#endif

  end subroutine calc_init_tblite

  module subroutine calc_eval_tblite(cl,c,ene,grad,stress,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    real*8, intent(out) :: ene
    real*8, intent(out) :: grad(:,:)
    real*8, intent(out), optional :: stress(3,3)
    character(len=:), allocatable, intent(out) :: errmsg

#ifdef HAVE_TBLITE
    type(c_ptr) :: err
    real(c_double), allocatable, target :: coords(:,:), grad_c(:)
    real(c_double), target :: lattice(3,3)
    real(c_double) :: e_c, sigma(9)
    integer :: i

    errmsg = ""
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
    errmsg = "critic2 was compiled without tblite support"
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
