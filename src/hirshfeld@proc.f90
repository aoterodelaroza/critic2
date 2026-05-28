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

! Hirshfeld integration
submodule (hirshfeld) proc
  use grid1mod, only: grid1
  implicit none

  integer, parameter :: mprops = 1

  ! ----- Hirshfeld-I module state (Bultinck 2007) -----
  ! Per-integration SCF state (qlo, qhi, frac, qfinal, isactive) lives on
  ! bas%hi_* in basindat (mirroring how bas%luw works for YT). Only the
  ! shared (Z, q) -> grid1 memoisation cache stays module-level here: it
  ! is conceptually like the agrid table in grid1mod and is safe to keep
  ! warm across multiple HIRSHFELD_I invocations.
  integer, parameter :: hi_qoff = 5      ! how many anion charges we admit
  integer, parameter :: hi_qcap = 100    ! cation cap (oversized to be safe)
  type(grid1), allocatable, target :: hi_cache(:,:)
  logical, allocatable :: hi_cache_tried(:,:)

contains

  !> Set the attractors for Hirshfeld integration. The actual
  !> integration is done using hirsh_weights.
  module subroutine hirsh_grid(s,bas)
    use systemmod, only: system
    use types, only: basindat
    use tools_io, only: ferror, faterr
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas

    integer :: i

    if (.not.s%isinit) &
       call ferror("hirsh_grid","system not initialized",faterr)
    if (.not.allocated(s%c)) &
       call ferror("hirsh_grid","system does not have crystal",faterr)

    ! Atoms are the attractors in this case
    allocate(bas%xattr(3,s%f(s%iref)%ncpcel))
    bas%xattr = 0d0
    if (bas%atexist) then
       bas%nattr = s%f(s%iref)%ncpcel
       do i = 1, s%f(s%iref)%ncpcel
          bas%xattr(:,i) = s%f(s%iref)%cpcel(i)%x
       end do
    else
       bas%nattr = 0
    end if

  end subroutine hirsh_grid

  !> For system s, calculate the hirshfeld weights for atom idb
  !> (complete list) on a grid and return it in w. The size of w bas%n
  !> determines the size of the grid and bas%f must contain the
  !> promolecular density. The size of w must be consistent with
  !> bas%n.
  module subroutine hirsh_weights(s,bas,idb,w)
    use systemmod, only: system
    use fragmentmod, only: fragment
    use types, only: basindat
    use param, only: VSMALL
    type(system), intent(inout) :: s
    type(basindat), intent(in) :: bas
    integer, intent(in) :: idb
    real*8, intent(out) :: w(:,:,:)

    type(fragment) :: fr
    real*8, allocatable :: faux(:,:,:)

    ! initialize
    w = 0d0
    fr%nat = 1
    fr%nspc = 1
    allocate(fr%at(1),fr%spc(1))

    ! Prepare a fragment with this atom
    fr%at(1)%x = s%c%atcel(idb)%x
    fr%at(1)%r = s%c%atcel(idb)%r
    fr%at(1)%is = 1
    fr%at(1)%cidx = idb
    fr%at(1)%idx = s%c%atcel(idb)%idx
    fr%at(1)%lvec = 0
    fr%spc(1) = s%c%spc(s%c%atcel(idb)%is)

    ! Calculate grid with the atomic density. The sum over cells
    ! turns into a sum over periodic copies of the atom.  call
    call s%c%promolecular_array3(faux,bas%n,fr=fr)

    ! hirshfeld weights
    w = faux / max(bas%f,VSMALL)

  end subroutine hirsh_weights

  !> Set the attractors for Voronoi integration and calculate
  !> the assignments of nodes to nuclei (bas%idg). The size
  !> of the grid is given by bas%n.
  module subroutine voronoi_grid(s,bas)
    use systemmod, only: system
    use tools_io, only: ferror, faterr
    use types, only: basindat
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas

    integer :: i

    if (.not.s%isinit) &
       call ferror("voronoi_grid","system not initialized",faterr)
    if (.not.allocated(s%c)) &
       call ferror("voronoi_grid","system does not have crystal",faterr)

    ! Atoms are the attractors in this case
    allocate(bas%xattr(3,s%f(s%iref)%ncpcel))
    bas%xattr = 0d0
    if (bas%atexist) then
       bas%nattr = s%f(s%iref)%ncpcel
       do i = 1, s%f(s%iref)%ncpcel
          bas%xattr(:,i) = s%f(s%iref)%cpcel(i)%x
       end do
    else
       bas%nattr = 0
    end if

    ! assign grid nodes to atoms
    call s%c%nearest_atom_grid(bas%n,bas%idg)

  end subroutine voronoi_grid

  ! Calculate hirshfeld populations and volumes using a mesh.
  module subroutine hirsh_nogrid()
    use fieldmod, only: field
    use meshmod, only: mesh
    use global, only: mesh_type, mesh_level
    use fragmentmod, only: fragment
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_center
    use param, only: VSMALL, im_rho

    type(field) :: fat
    type(mesh) :: m0, mat, mrho
    integer :: i, iat, cidx
    real*8 :: vsum, asum, ntotal, vtotal
    real*8, allocatable :: nat(:), vat(:)
    type(fragment) :: fr
    integer :: prop(1)

    ! header
    write (uout,'("* Hirshfeld atomic electron populations (using mesh integration)")')
    write (uout,'("  Reference field: ",A)') string(sy%iref)

    ! generate the meshes
    call m0%gen(sy%c,mesh_type,mesh_level)
    call mrho%gen(sy%c,mesh_type,mesh_level)
    call mat%gen(sy%c,mesh_type,mesh_level)

    ! report
    write (uout,'("+ Mesh details")')
    call m0%report()

    ! Initialize accumulation and fragment
    allocate(nat(sy%c%nneq),vat(sy%c%nneq))
    ntotal = 0d0
    vtotal = 0d0
    nat = 0d0
    vat = 0d0
    fr%nat = 1
    fr%nspc = 1
    allocate(fr%at(1),fr%spc(1))

    ! calculate the density and promolecular density on the mesh
    prop(1) = im_rho
    call mrho%fill(sy%f(sy%iref),prop,.not.sy%c%ismolecule)
    call m0%fill(sy%f(0),prop,.not.sy%c%ismolecule)

    do iat = 1, sy%c%nneq
       ! Fetch a representative with that iat
       cidx = 0
       do i = 1, sy%c%ncel
          if (sy%c%atcel(i)%idx == iat) then
             cidx = i
             exit
          end if
       end do

       ! Prepare a fragment with this atom
       fr%at(1)%x = sy%c%atcel(cidx)%x
       fr%at(1)%r = sy%c%atcel(cidx)%r
       fr%at(1)%is = 1
       fr%at(1)%idx = iat
       fr%at(1)%cidx = cidx
       fr%at(1)%lvec = 0
       fr%spc(1) = sy%c%spc(sy%c%at(iat)%is)

       ! Load the promolecular field with that fragment
       call fat%load_promolecular(sy%c,-1,"",fr)

       ! Calculate mesh with the atomic density. The sum over cells
       ! turns into a sum over periodic copies of the atom.
       prop(1) = im_rho
       call mat%fill(fat,prop,.not.sy%c%ismolecule)

       ! hirshfeld weights (with mesh weights)
       mat%f(:,1) = mat%f(:,1) / max(m0%f(:,1),VSMALL) * mrho%w

       ! accumulate
       asum = sum(mrho%f(:,1) * mat%f(:,1))
       vsum = sum(mat%f(:,1))
       nat(iat) = asum
       ntotal = ntotal + asum * sy%c%at(iat)%mult
       vat(iat) = vsum
       vtotal = vtotal + vsum * sy%c%at(iat)%mult
    end do

    ! write the results
    write (uout,'("+ Hirshfeld integration results")')
    write (uout,'("# N_hirsh = atomic electron populations")')
    write (uout,'("# V_hirsh = atomic volumes ")')
    write (uout,'("#nneq mult name   N_hirsh          V_hirsh")')
    do iat = 1, sy%c%nneq
       write (uout,'(5(A," "))') string(iat,length=4,justify=ioj_center), &
          string(sy%c%at(iat)%mult,length=4,justify=ioj_center), &
          string(sy%c%at(iat)%name,length=5,justify=ioj_center), &
          string(nat(iat),'f',length=16,decimal=10,justify=3), &
          string(vat(iat),'f',length=16,decimal=10,justify=3)
    end do
    write (uout,'("# total number of electrons: ",A)') string(ntotal,'e',decimal=10)
    write (uout,'("# total volume: ",A)') string(vtotal,'e',decimal=10)
    write (uout,*)

  end subroutine hirsh_nogrid

  !> Iterative Hirshfeld (Hirshfeld-I) SCF driver. On entry the
  !> attractors are already set up by hirsh_grid and bas%f holds the
  !> initial neutral promolecular density. On exit hirsh_isactive=.true.
  !> and the per-cellatom (qlo, qhi, frac) state is ready for use by
  !> hirsh_i_eval inside the field integrator.
  module subroutine hirsh_i_driver(s,bas)
    use systemmod, only: system
    use types, only: basindat
    use grid1mod, only: agrid
    use global, only: cutrad
    use tools_io, only: uout, string, ferror, faterr, warning
    use param, only: VSMALL, maxzat, icrd_crys
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas

    integer :: nca, iat, iter, i1, i2, i3, ii, iz, isp
    integer :: nat
    integer, allocatable :: nid(:), lvec(:,:)
    real*8, allocatable :: dist(:), rcutmax(:,:)
    real*8, allocatable :: nelec(:), vat(:), qprev(:)
    real*8 :: x0(3), xdelta(3,3), rho_i, rho_promol, dQ, dQmax, fac
    real*8 :: rho_lo, rho_hi
    logical :: converged

    if (.not.s%isinit) call ferror("hirsh_i_driver","system not initialized",faterr)
    if (.not.allocated(s%c)) call ferror("hirsh_i_driver","no crystal",faterr)

    nca = s%c%ncel

    ! allocate / reset SCF state on bas (mirrors how YT stores bas%luw)
    call hirsh_i_cleanup(bas)
    allocate(bas%hi_qlo(nca), bas%hi_qhi(nca), bas%hi_frac(nca), bas%hi_qfinal(nca))
    allocate(nelec(nca), vat(nca), qprev(nca))
    bas%hi_qlo = 0
    bas%hi_qhi = 0
    bas%hi_frac = 0d0
    bas%hi_qfinal = 0d0
    qprev = 0d0

    ! grid step vectors
    do ii = 1, 3
       xdelta(:,ii) = 0d0
       xdelta(ii,ii) = 1d0 / real(bas%n(ii),8)
    end do

    ! Cutoff per species using neutral grid extent (a reasonable bound:
    ! anion grids are slightly more diffuse but typically still within
    ! cutrad).
    if (allocated(rcutmax)) deallocate(rcutmax)
    allocate(rcutmax(s%c%nspc,2))
    rcutmax = 0d0
    do isp = 1, s%c%nspc
       iz = s%c%spc(isp)%z
       if (iz <= 0 .or. iz > maxzat) cycle
       if (agrid(iz)%isinit) then
          rcutmax(isp,2) = min(cutrad(iz),agrid(iz)%rmax)
       else
          call ferror('hirsh_i_driver','iterative hirshfeld requires atomic grids',faterr)
       end if
    end do

    ! header
    write (uout,'("* Hirshfeld-I (iterative Hirshfeld) SCF")')
    write (uout,'("  Bultinck, Van Alsenoy, Ayers, Carbo-Dorca,")')
    write (uout,'("    J. Chem. Phys. 126, 144111 (2007). (10.1063/1.2715563)")')
    write (uout,'("+ Tolerance (max |dQ|): ",A)') string(bas%hi_tol,'e',decimal=2)
    write (uout,'("+ Max iterations: ",A)') string(bas%hi_maxit)
    if (allocated(bas%hi_wfcdir)) then
       if (len_trim(bas%hi_wfcdir) > 0) &
          write (uout,'("+ User WFCDIR: ",A)') trim(bas%hi_wfcdir)
    end if
    write (uout,'("# iter  max|dQ|       sum(Q)        atom-charges")')

    converged = .false.
    iter = 0
    do
       iter = iter + 1

       ! preload cache for all (z, qlo) and (z, qhi) needed this iter
       ! (must happen serially before the OMP sweep below)
       call hirsh_i_preload(s,bas)

       ! activate evaluator so promol can be rebuilt using charged refs
       bas%hi_isactive = .true.

       ! single grid sweep: rebuild bas%f (promol) and accumulate per-atom N
       nelec = 0d0
       vat = 0d0
       !$omp parallel do private(x0,nat,rho_promol,rho_i,fac) firstprivate(nid,dist,lvec) &
       !$omp& reduction(+:nelec,vat)
       do i3 = 1, bas%n(3)
          do i2 = 1, bas%n(2)
             do i1 = 1, bas%n(1)
                x0 = (i1-1)*xdelta(:,1) + (i2-1)*xdelta(:,2) + (i3-1)*xdelta(:,3)
                call s%c%list_near_atoms(x0,icrd_crys,.false.,nat,nid,dist,lvec,up2dsp=rcutmax)
                ! promol
                rho_promol = 0d0
                do ii = 1, nat
                   call hirsh_i_eval(bas,nid(ii),dist(ii),rho_i)
                   rho_promol = rho_promol + rho_i
                end do
                bas%f(i1,i2,i3) = rho_promol
                fac = 1d0 / max(rho_promol,VSMALL)
                ! per-atom Hirshfeld weight w_A = rho_A/promol; integrate
                ! over reference density (only the population for now).
                do ii = 1, nat
                   call hirsh_i_eval(bas,nid(ii),dist(ii),rho_i)
                   nelec(nid(ii)) = nelec(nid(ii)) + fac * rho_i * s%f(s%iref)%grid%f(i1,i2,i3)
                   vat(nid(ii)) = vat(nid(ii)) + fac * rho_i
                end do
             end do
          end do
       end do
       !$omp end parallel do
       nelec = nelec * s%c%omega / product(bas%n)
       vat = vat * s%c%omega / product(bas%n)

       ! update per-cellatom Q and (qlo,qhi,frac)
       dQmax = 0d0
       do iat = 1, nca
          iz = s%c%spc(s%c%atcel(iat)%is)%z
          bas%hi_qfinal(iat) = real(iz,8) - nelec(iat)
          ! clamp Q into [-hi_qoff, iz) so qlo/qhi stay in cache range.
          ! qhi may equal iz (fully-stripped ion); interp returns 0 there,
          ! which is the physically correct empty-density limit.
          if (bas%hi_qfinal(iat) < -dble(hi_qoff)) bas%hi_qfinal(iat) = -dble(hi_qoff) + 1d-3
          if (bas%hi_qfinal(iat) > dble(iz) - 1d-3) bas%hi_qfinal(iat) = dble(iz) - 1d-3
          bas%hi_qlo(iat) = floor(bas%hi_qfinal(iat))
          bas%hi_qhi(iat) = bas%hi_qlo(iat) + 1
          bas%hi_frac(iat) = bas%hi_qfinal(iat) - dble(bas%hi_qlo(iat))
          dQ = abs(bas%hi_qfinal(iat) - qprev(iat))
          if (dQ > dQmax) dQmax = dQ
       end do

       ! status line
       write (uout,'(2X,A,2X,A,2X,A,2X,*(A," "))') &
          string(iter,length=4,justify=2), &
          string(dQmax,'e',decimal=3), &
          string(sum(bas%hi_qfinal),'f',decimal=5), &
          (trim(string(bas%hi_qfinal(iat),'f',decimal=4)), iat=1,nca)

       if (dQmax < bas%hi_tol) then
          converged = .true.
          exit
       end if
       if (iter >= bas%hi_maxit) exit
       qprev = bas%hi_qfinal
    end do

    if (.not.converged) then
       call ferror('hirsh_i_driver', &
          'Hirshfeld-I SCF did not converge within HIMAXIT iterations; using current state',warning)
    else
       write (uout,'("+ Hirshfeld-I SCF converged in ",A," iterations")') string(iter)
    end if
    write (uout,*)

    deallocate(nelec,vat,qprev)
    if (allocated(rcutmax)) deallocate(rcutmax)

  end subroutine hirsh_i_driver

  !> Evaluate the (interpolated) charged atomic density at distance r
  !> from cell-atom idcel using the SCF state stored on bas.
  module subroutine hirsh_i_eval(bas,idcel,dist,rho)
    use systemmod, only: sy
    use grid1mod, only: agrid
    use types, only: basindat
    type(basindat), intent(in) :: bas
    integer, intent(in) :: idcel
    real*8, intent(in) :: dist
    real*8, intent(out) :: rho

    integer :: iz, qlo, qhi
    real*8 :: rlo, rhi, rdum1, rdum2, fr
    type(grid1), pointer :: glo, ghi

    rho = 0d0
    if (idcel <= 0) return
    iz = sy%c%spc(sy%c%atcel(idcel)%is)%z

    if (.not.bas%hi_isactive) then
       ! fall back to neutral (used pre-SCF and when iterative is off)
       if (iz > 0) call agrid(iz)%interp(dist,rho,rdum1,rdum2)
       return
    end if

    qlo = bas%hi_qlo(idcel)
    qhi = bas%hi_qhi(idcel)
    fr  = bas%hi_frac(idcel)

    call hirsh_i_get_grid(iz,qlo,glo)
    call hirsh_i_get_grid(iz,qhi,ghi)
    rlo = 0d0; rhi = 0d0
    if (associated(glo)) call glo%interp(dist,rlo,rdum1,rdum2)
    if (associated(ghi)) call ghi%interp(dist,rhi,rdum1,rdum2)
    rho = (1d0 - fr) * rlo + fr * rhi

  end subroutine hirsh_i_eval

  !> Tear down the Hirshfeld-I per-integration state on bas and the
  !> shared (Z, q) memoisation cache. Allocatable components inside
  !> each cached grid1 are released by Fortran 2003+ deallocation
  !> semantics when the array is deallocated.
  module subroutine hirsh_i_cleanup(bas)
    use types, only: basindat
    type(basindat), intent(inout) :: bas
    if (allocated(hi_cache)) deallocate(hi_cache)
    if (allocated(hi_cache_tried)) deallocate(hi_cache_tried)
    if (allocated(bas%hi_qlo)) deallocate(bas%hi_qlo)
    if (allocated(bas%hi_qhi)) deallocate(bas%hi_qhi)
    if (allocated(bas%hi_frac)) deallocate(bas%hi_frac)
    if (allocated(bas%hi_qfinal)) deallocate(bas%hi_qfinal)
    bas%hi_isactive = .false.
  end subroutine hirsh_i_cleanup

  !> Thread-safe read-only lookup of a (z,q) radial density. Returns
  !> g => null() if the cache hasn't been populated yet (the caller can
  !> handle that by calling hirsh_i_preload before parallel regions).
  subroutine hirsh_i_get_grid(iz,q,g)
    use grid1mod, only: agrid
    use param, only: maxzat
    integer, intent(in) :: iz, q
    type(grid1), pointer :: g

    g => null()
    if (iz <= 0 .or. iz > maxzat) return

    ! neutral: use the agrid table directly
    if (q == 0) then
       if (agrid(iz)%isinit) g => agrid(iz)
       return
    end if

    if (.not.allocated(hi_cache)) return
    if (q < lbound(hi_cache,2) .or. q > ubound(hi_cache,2)) return

    if (hi_cache(iz,q)%isinit) then
       g => hi_cache(iz,q)
       return
    end if

    ! prior load failed or never attempted; fall back to neutral
    if (agrid(iz)%isinit) g => agrid(iz)
  end subroutine hirsh_i_get_grid

  !> Serial preload of all (z, q) atomic-density grids needed by the
  !> current SCF state. Must be called before any parallel evaluator
  !> sweep. Resolution order matches the documentation: (1) user WFCDIR
  !> file, (2) built-in cation table via read_db (q>=0), (3) anion
  !> extrapolation via read_db (q<0).
  subroutine hirsh_i_preload(s,bas)
    use systemmod, only: system
    use types, only: basindat
    use tools_io, only: lower, nameguess, ferror, warning, string
    use param, only: dirsep, maxzat
    type(system), intent(in) :: s
    type(basindat), intent(in) :: bas

    integer :: iat, iz
    character(len=:), allocatable :: fname
    logical :: exist, have_wfcdir

    if (.not.allocated(hi_cache)) then
       allocate(hi_cache(maxzat, -hi_qoff:hi_qcap))
       allocate(hi_cache_tried(maxzat, -hi_qoff:hi_qcap))
       hi_cache_tried = .false.
    end if

    have_wfcdir = .false.
    if (allocated(bas%hi_wfcdir)) then
       if (len_trim(bas%hi_wfcdir) > 0) have_wfcdir = .true.
    end if

    do iat = 1, s%c%ncel
       iz = s%c%spc(s%c%atcel(iat)%is)%z
       if (iz <= 0 .or. iz > maxzat) cycle
       call load_one(iz, bas%hi_qlo(iat))
       call load_one(iz, bas%hi_qhi(iat))
    end do

  contains
    subroutine load_one(iz,q)
      integer, intent(in) :: iz, q
      if (q == 0) return                                ! neutral lives in agrid
      if (q < lbound(hi_cache,2) .or. q > ubound(hi_cache,2)) return
      if (hi_cache(iz,q)%isinit) return
      if (hi_cache_tried(iz,q)) return
      hi_cache_tried(iz,q) = .true.

      ! (1) WFCDIR/<elem>_q<+|->N.wfc
      if (have_wfcdir) then
         if (q >= 0) then
            fname = trim(bas%hi_wfcdir) // dirsep // lower(nameguess(iz)) // &
                    "_q+" // trim(string(q)) // ".wfc"
         else
            fname = trim(bas%hi_wfcdir) // dirsep // lower(nameguess(iz)) // &
                    "_q" // trim(string(q)) // ".wfc"
         end if
         inquire(file=fname,exist=exist)
         if (exist) then
            call hirsh_i_load_userfile(hi_cache(iz,q),fname,iz,q)
            if (hi_cache(iz,q)%isinit) return
         end if
      end if

      ! (2),(3) Built-in tables / anion extrapolation via read_db.
      call hi_cache(iz,q)%read_db(iz,q)
      if (hi_cache(iz,q)%isinit) return

      call ferror('hirsh_i_preload', &
         'Could not load atomic density for z=' // trim(string(iz)) // &
         ' q=' // trim(string(q)) // '; using neutral',warning)
    end subroutine load_one
  end subroutine hirsh_i_preload

  !> Load a user-supplied .wfc file in critic2 format into a grid1
  !> using the public grid1_read_file wrapper.
  subroutine hirsh_i_load_userfile(g,fname,iz,q)
    use grid1mod, only: grid1_read_file
    type(grid1), intent(inout) :: g
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iz, q
    call grid1_read_file(g,fname,iz,q)
  end subroutine hirsh_i_load_userfile

end submodule proc
