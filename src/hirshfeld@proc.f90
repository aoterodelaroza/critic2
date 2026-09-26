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
  implicit none

  !xx! private procedures
  ! subroutine hirsh_register(c)
  ! function hirsh_irep(c) result(irep)
  ! subroutine hirsh_point(c,x0,icrd,xn,rcut,val,acc,nid,dist,lvec,rhoa)
  ! subroutine hirsh_accumulate(c,xn,acc,gridf,mx,mw,mf)
  ! subroutine hirsh_iterate(c,xn,tol,maxit,ncore,gridf,mx,mw,mf)
  ! subroutine hirsh_symmetrize(c,xn)
  ! subroutine hirsh_report(c,xn,nmin,nmax,pinned,conv,niter,dmax,ntotint)

contains

  !> Set the attractors for Hirshfeld integration, calculate the
  !> reference populations of the atoms (the neutral atoms, or the
  !> Hirshfeld-I populations if bas%hirsh_iter), and fill bas%f with
  !> the promolecular density for those populations. The integration
  !> itself is done in intgrid_hirshfeld_fields.
  module subroutine hirsh_grid(s,bas)
    use systemmod, only: system
    use types, only: basindat
    use tools_io, only: ferror, faterr, uout
    use param, only: icrd_crys
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas

    integer :: i, j, is, i1, i2, i3, nat
    real*8 :: x0(3), f
    real*8, allocatable :: rcut(:,:), ncore(:), dist(:)
    integer, allocatable :: nid(:), lvec(:,:)

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

    ! reference populations: the neutral atoms or Hirshfeld-I
    call hirsh_register(s%c)
    if (allocated(bas%hirsh_n)) deallocate(bas%hirsh_n)
    allocate(bas%hirsh_n(s%c%ncel))
    do i = 1, s%c%ncel
       bas%hirsh_n(i) = s%c%spc(s%c%atcel(i)%is)%z
    end do
    if (bas%hirsh_iter) then
       ! A pseudopotential grid holds only the valence density: the
       ! frozen core of each atom is added to its population, instead
       ! of integrating the cusped core density on the grid.
       allocate(ncore(s%c%ncel))
       ncore = 0d0
       if (s%f(s%iref)%usecore) then
          do i = 1, s%c%ncel
             is = s%c%atcel(i)%is
             if (s%f(s%iref)%zpsp(is) > 0) &
                ncore(i) = s%c%spc(is)%z - s%f(s%iref)%zpsp(is)
          end do
       end if
       call hirsh_iterate(s%c,bas%hirsh_n,bas%hirsh_tol,bas%hirsh_maxit,ncore=ncore,&
          gridf=s%f(s%iref)%grid%f)
       if (any(ncore > 0d0)) then
          write (uout,'("+ The Hirshfeld-I populations (N_HI) include the frozen cores (Z - ZPSP).")')
          write (uout,'("  The integrated properties below use the valence field only.")')
          write (uout,*)
       end if
    end if

    ! the promolecular density for the reference populations
    call hirsh_cutoffs(s%c,bas%hirsh_n,rcut)
    !$omp parallel do private(x0,nat,f) firstprivate(nid,dist,lvec) schedule(dynamic)
    do i3 = 1, bas%n(3)
       do i2 = 1, bas%n(2)
          do i1 = 1, bas%n(1)
             x0 = real((/i1,i2,i3/)-1,8) / real(bas%n,8)
             call s%c%list_near_atoms(x0,icrd_crys,.false.,nat,nid,dist,lvec,up2dsp=rcut)
             f = 0d0
             do j = 1, nat
                f = f + hirsh_rho(s%c,bas%hirsh_n,nid(j),dist(j))
             end do
             bas%f(i1,i2,i3) = f
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine hirsh_grid

  !> For system s, calculate the hirshfeld weights for atom idb
  !> (complete list) on a grid and return it in w. The size of w bas%n
  !> determines the size of the grid, bas%f must contain the
  !> promolecular density, and bas%hirsh_n the reference populations
  !> (see hirsh_grid). The size of w must be consistent with bas%n.
  module subroutine hirsh_weights(s,bas,idb,w)
    use systemmod, only: system
    use types, only: basindat
    use param, only: VSMALL, icrd_crys
    type(system), intent(inout) :: s
    type(basindat), intent(in) :: bas
    integer, intent(in) :: idb
    real*8, intent(out) :: w(:,:,:)

    integer :: j, i1, i2, i3, nat
    real*8 :: x0(3), f
    real*8, allocatable :: rcut(:,:), dist(:)
    integer, allocatable :: nid(:), lvec(:,:)

    ! the sum over the periodic copies of atom idb, divided by the promolecule
    call hirsh_cutoffs(s%c,bas%hirsh_n,rcut)
    !$omp parallel do private(x0,nat,f) firstprivate(nid,dist,lvec) schedule(dynamic)
    do i3 = 1, bas%n(3)
       do i2 = 1, bas%n(2)
          do i1 = 1, bas%n(1)
             x0 = real((/i1,i2,i3/)-1,8) / real(bas%n,8)
             call s%c%list_near_atoms(x0,icrd_crys,.false.,nat,nid,dist,lvec,up2dsp=rcut,id0=idb)
             f = 0d0
             do j = 1, nat
                f = f + hirsh_rho(s%c,bas%hirsh_n,idb,dist(j))
             end do
             w(i1,i2,i3) = f / max(bas%f(i1,i2,i3),VSMALL)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine hirsh_weights

  !> Reference density at distance r of the complete-cell atom i in
  !> crystal c, with the reference populations xn(1:ncel). Zero for
  !> atoms without atomic densities. The densities of the species
  !> must have been registered (hirsh_register).
  module function hirsh_rho(c,xn,i,r)
    use crystalmod, only: crystal
    use grid1mod, only: grid1_states_rho
    use param, only: maxzat
    type(crystal), intent(in) :: c
    real*8, intent(in) :: xn(:)
    integer, intent(in) :: i
    real*8, intent(in) :: r
    real*8 :: hirsh_rho

    integer :: iz

    iz = c%spc(c%atcel(i)%is)%z
    if (iz <= 0 .or. iz > maxzat) then
       hirsh_rho = 0d0
    else
       hirsh_rho = grid1_states_rho(iz,xn(i),r)
    end if

  end function hirsh_rho

  !> Cutoff radii for list_near_atoms (up2dsp) in rcut(1:nspc,2): for
  !> each species, the largest extent of the reference densities of
  !> its atoms, with populations xn(1:ncel) (see grid1_states_rcut).
  !> The densities must have been registered (hirsh_register).
  module subroutine hirsh_cutoffs(c,xn,rcut)
    use crystalmod, only: crystal
    use grid1mod, only: grid1_states_rcut
    use param, only: maxzat
    type(crystal), intent(in) :: c
    real*8, intent(in) :: xn(:)
    real*8, allocatable, intent(inout) :: rcut(:,:)

    integer :: i, is, iz

    if (allocated(rcut)) deallocate(rcut)
    allocate(rcut(c%nspc,2))
    rcut = 0d0
    do i = 1, c%ncel
       is = c%atcel(i)%is
       iz = c%spc(is)%z
       if (iz <= 0 .or. iz > maxzat) cycle
       rcut(is,2) = max(rcut(is,2),grid1_states_rcut(iz,xn(i)))
    end do

  end subroutine hirsh_cutoffs

  !> Parse one option of the HIRSHFELD keyword common to the grid and
  !> mesh paths: ITERATIVE, TOL tol.r, and MAXIT maxit.i. If word is
  !> one of them, read its value from line (at lp), set iter, tol, or
  !> maxit, and return found = .true.
  module subroutine hirsh_option(word,line,lp,iter,tol,maxit,found)
    use global, only: eval_next
    use tools_io, only: equal, isinteger, ferror, faterr
    character*(*), intent(in) :: word
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp
    logical, intent(inout) :: iter
    real*8, intent(inout) :: tol
    integer, intent(inout) :: maxit
    logical, intent(out) :: found

    logical :: ok

    found = .true.
    if (equal(word,"iterative")) then
       iter = .true.
    elseif (equal(word,"tol")) then
       ok = eval_next(tol,line,lp)
       if (.not.ok .or. tol <= 0d0) &
          call ferror("hirsh_option","Wrong TOL in HIRSHFELD",faterr,line,syntax=.true.)
    elseif (equal(word,"maxit")) then
       ok = isinteger(maxit,line,lp)
       if (.not.ok .or. maxit < 1) &
          call ferror("hirsh_option","Wrong MAXIT in HIRSHFELD",faterr,line,syntax=.true.)
    else
       found = .false.
    end if

  end subroutine hirsh_option

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

  !> Calculate hirshfeld populations and volumes using a mesh. The
  !> input line is HIRSHFELD [ITERATIVE] [TOL t.r] [MAXIT n.i].
  module subroutine hirsh_nogrid(line)
    use meshmod, only: mesh
    use global, only: mesh_type, mesh_level
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_center, lgetword, ferror, faterr
    use param, only: im_rho
    character*(*), intent(in) :: line

    type(mesh) :: m
    integer :: i, lp, iat, maxit
    real*8 :: tol
    logical :: iter, found
    character(len=:), allocatable :: word
    real*8, allocatable :: xn(:), acc(:,:)
    integer, allocatable :: irep(:)
    integer :: prop(1)

    ! parse the options
    iter = .false.
    tol = hirsh_tol_def
    maxit = hirsh_maxit_def
    lp = 1
    word = lgetword(line,lp)
    do while (.true.)
       word = lgetword(line,lp)
       call hirsh_option(word,line,lp,iter,tol,maxit,found)
       if (found) cycle
       if (len_trim(word) > 0) then
          call ferror("hirsh_nogrid","Unknown extra keyword (WCUBE, ONLY, and JSON need a grid): " &
             // word,faterr,line,syntax=.true.)
       else
          exit
       end if
    end do

    ! header
    write (uout,'("* Hirshfeld atomic electron populations (using mesh integration)")')
    write (uout,'("  Reference field: ",A)') string(sy%iref)

    ! generate the mesh
    call m%gen(sy%c,mesh_type,mesh_level,sy%f(sy%iref)%zpsp)
    write (uout,'("+ Mesh details")')
    call m%report()

    ! density on the mesh and reference populations
    prop(1) = im_rho
    call m%fill(sy%f(sy%iref),prop,.not.sy%c%ismolecule)
    call hirsh_register(sy%c)
    allocate(xn(sy%c%ncel))
    do i = 1, sy%c%ncel
       xn(i) = sy%c%spc(sy%c%atcel(i)%is)%z
    end do
    if (iter) &
       call hirsh_iterate(sy%c,xn,tol,maxit,mx=m%x,mw=m%w,mf=m%f(:,1))

    ! populations and volumes of the complete-cell atoms (including their periodic copies)
    allocate(acc(sy%c%ncel,2))
    call hirsh_accumulate(sy%c,xn,acc,mx=m%x,mw=m%w,mf=m%f(:,1))

    ! write the results
    irep = hirsh_irep(sy%c)
    write (uout,'("+ Hirshfeld integration results")')
    write (uout,'("# N_hirsh = atomic electron populations")')
    write (uout,'("# V_hirsh = atomic volumes ")')
    write (uout,'("#nneq mult name   N_hirsh          V_hirsh")')
    do iat = 1, sy%c%nneq
       write (uout,'(5(A," "))') string(iat,length=4,justify=ioj_center), &
          string(sy%c%at(iat)%mult,length=4,justify=ioj_center), &
          string(sy%c%at(iat)%name,length=5,justify=ioj_center), &
          string(acc(irep(iat),1),'f',length=16,decimal=10,justify=3), &
          string(acc(irep(iat),2),'f',length=16,decimal=10,justify=3)
    end do
    write (uout,'("# total number of electrons: ",A)') string(sum(acc(:,1)),'e',decimal=10)
    write (uout,'("# total volume: ",A)') string(sum(acc(:,2)),'e',decimal=10)
    write (uout,*)

  end subroutine hirsh_nogrid

  !xx! private procedures

  !> Register the charge-state densities of all species in crystal c.
  subroutine hirsh_register(c)
    use crystalmod, only: crystal
    use grid1mod, only: grid1_register_states
    use param, only: maxzat
    type(crystal), intent(in) :: c

    integer :: is

    do is = 1, c%nspc
       if (c%spc(is)%z <= 0 .or. c%spc(is)%z > maxzat) cycle
       call grid1_register_states(c%spc(is)%z)
    end do

  end subroutine hirsh_register

  !> A representative complete-cell atom for each non-equivalent atom
  !> in crystal c.
  function hirsh_irep(c) result(irep)
    use crystalmod, only: crystal
    type(crystal), intent(in) :: c
    integer :: irep(c%nneq)

    integer :: i

    irep = 0
    do i = c%ncel, 1, -1
       irep(c%atcel(i)%idx) = i
    end do

  end function hirsh_irep

  !> Hirshfeld partition at point x0 (icrd coordinates) of crystal
  !> c: val(:) times the weight of each atom in the environment,
  !> w_A = rho_A / sum_B rho_B, is added to acc(A,:), with A the
  !> complete-cell atom. The reference densities are those of the
  !> populations xn(1:ncel), cut off at rcut (hirsh_cutoffs). nid,
  !> dist, lvec, and rhoa are work space.
  subroutine hirsh_point(c,x0,icrd,xn,rcut,val,acc,nid,dist,lvec,rhoa)
    use crystalmod, only: crystal
    use types, only: realloc
    use param, only: VSMALL
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: icrd
    real*8, intent(in) :: xn(:)
    real*8, intent(in) :: rcut(:,:)
    real*8, intent(in) :: val(:)
    real*8, intent(inout) :: acc(:,:)
    integer, allocatable, intent(inout) :: nid(:), lvec(:,:)
    real*8, allocatable, intent(inout) :: dist(:), rhoa(:)

    integer :: i, nat
    real*8 :: pro

    call c%list_near_atoms(x0,icrd,.false.,nat,nid,dist,lvec,up2dsp=rcut)
    if (nat == 0) return
    if (.not.allocated(rhoa)) then
       allocate(rhoa(nat))
    elseif (size(rhoa) < nat) then
       call realloc(rhoa,nat)
    end if

    pro = 0d0
    do i = 1, nat
       rhoa(i) = hirsh_rho(c,xn,nid(i),dist(i))
       pro = pro + rhoa(i)
    end do
    if (pro < VSMALL) return
    do i = 1, nat
       acc(nid(i),:) = acc(nid(i),:) + val * (rhoa(i) / pro)
    end do

  end subroutine hirsh_point

  !> Integrate with the Hirshfeld weights of the reference
  !> populations xn(1:ncel) of crystal c: acc(A,1) = int w_A rho,
  !> and, if size(acc,2) > 1, acc(A,2) = int w_A (the volume). The
  !> density is either a grid in the unit cell (gridf) or its values
  !> (mf) on a mesh with points mx (Cartesian) and weights mw.
  subroutine hirsh_accumulate(c,xn,acc,gridf,mx,mw,mf)
    use crystalmod, only: crystal
    use param, only: icrd_crys, icrd_cart
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: xn(:)
    real*8, intent(inout) :: acc(:,:)
    real*8, intent(in), optional :: gridf(:,:,:)
    real*8, intent(in), optional :: mx(:,:), mw(:), mf(:)

    integer :: j, i1, i2, i3, n(3)
    real*8 :: x0(3), val(size(acc,2))
    real*8, allocatable :: rcut(:,:), rhoa(:), dist(:)
    integer, allocatable :: nid(:), lvec(:,:)

    call hirsh_cutoffs(c,xn,rcut)
    acc = 0d0
    if (present(gridf)) then
       n = shape(gridf)
       !$omp parallel do private(x0,val) firstprivate(nid,dist,lvec,rhoa) reduction(+:acc) &
       !$omp schedule(dynamic)
       do i3 = 1, n(3)
          do i2 = 1, n(2)
             do i1 = 1, n(1)
                x0 = real((/i1,i2,i3/)-1,8) / real(n,8)
                val = 1d0
                val(1) = gridf(i1,i2,i3)
                call hirsh_point(c,x0,icrd_crys,xn,rcut,val,acc,nid,dist,lvec,rhoa)
             end do
          end do
       end do
       !$omp end parallel do
       acc = acc * c%omega / real(product(n),8)
    else
       !$omp parallel do private(val) firstprivate(nid,dist,lvec,rhoa) reduction(+:acc) &
       !$omp schedule(dynamic)
       do j = 1, size(mw)
          val = mw(j)
          val(1) = mf(j) * mw(j)
          call hirsh_point(c,mx(:,j),icrd_cart,xn,rcut,val,acc,nid,dist,lvec,rhoa)
       end do
       !$omp end parallel do
    end if

  end subroutine hirsh_accumulate

  !> Iterative Hirshfeld (Hirshfeld-I) populations of the
  !> complete-cell atoms of crystal c. On input, xn(1:ncel) is the
  !> initial guess (usually the neutral atoms); on output, the
  !> converged populations. The iteration stops when max|dN| < tol or
  !> after maxit iterations. If present, ncore(1:ncel) is added to the
  !> integrated populations (frozen-core electrons of pseudopotential
  !> densities). The density is either a grid in the unit cell (gridf)
  !> or its values (mf) on a mesh with points mx (Cartesian) and
  !> weights mw.
  subroutine hirsh_iterate(c,xn,tol,maxit,ncore,gridf,mx,mw,mf)
    use crystalmod, only: crystal
    use grid1mod, only: sgrid
    use tools_io, only: uout, string, ioj_right
    use param, only: maxzat
    type(crystal), intent(inout) :: c
    real*8, intent(inout) :: xn(:)
    real*8, intent(in) :: tol
    integer, intent(in) :: maxit
    real*8, intent(in), optional :: ncore(:)
    real*8, intent(in), optional :: gridf(:,:,:)
    real*8, intent(in), optional :: mx(:,:), mw(:), mf(:)

    integer :: it, i, iz, niter
    real*8 :: dmax, ntotint
    real*8, allocatable :: xold(:), acc(:,:), nmin(:), nmax(:)
    logical, allocatable :: pinned(:)
    logical :: conv

    ! range of charge states of each atom
    allocate(nmin(c%ncel),nmax(c%ncel),pinned(c%ncel),xold(c%ncel),acc(c%ncel,1))
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       nmin(i) = xn(i)
       nmax(i) = xn(i)
       if (iz <= 0 .or. iz > maxzat) cycle
       nmin(i) = sgrid(iz)%nmin
       nmax(i) = sgrid(iz)%nmax
    end do

    write (uout,'("+ Iterative Hirshfeld (Hirshfeld-I) reference populations")')
    write (uout,'("  Please cite: ")')
    write (uout,'("    P. Bultinck, C. Van Alsenoy, P. W. Ayers, and R. Carbo-Dorca,")')
    write (uout,'("       J. Chem. Phys. 126, 144111 (2007). (10.1063/1.2715563)")')
    write (uout,'("  Reference densities interpolated between integer charge states")')
    write (uout,'("  Convergence: max|dN| < ",A," (maximum ",A," iterations)")') &
       string(tol,'e',decimal=2), string(maxit)
    if (present(ncore)) then
       if (any(ncore > 0d0)) &
          write (uout,'("  Frozen-core electrons (Z - ZPSP) added to the integrated valence populations")')
    end if
    write (uout,'("#   it      max|dN|")')

    xn = min(max(xn,nmin),nmax)
    pinned = .false.
    ntotint = sum(xn)
    conv = .false.
    niter = maxit
    dmax = 0d0
    do it = 1, maxit
       ! populations for the current reference densities
       xold = xn
       call hirsh_accumulate(c,xold,acc,gridf,mx,mw,mf)
       xn = acc(:,1)
       if (present(ncore)) xn = xn + ncore
       ntotint = sum(xn)

       ! symmetry-equivalent atoms share the population; clamp to the
       ! available charge states (the residual is measured after the
       ! clamp, so atoms pinned at the end of the range converge too)
       call hirsh_symmetrize(c,xn)
       pinned = (xn < nmin .or. xn > nmax)
       xn = min(max(xn,nmin),nmax)

       dmax = maxval(abs(xn - xold))
       write (uout,'(2X,A,4X,A)') string(it,length=4,justify=ioj_right), string(dmax,'e',decimal=4)
       if (dmax < tol) then
          conv = .true.
          niter = it
          exit
       end if
    end do

    call hirsh_report(c,xn,nmin,nmax,pinned,conv,niter,dmax,ntotint)

  end subroutine hirsh_iterate

  !> Average the populations xn(1:ncel) over symmetry-equivalent atoms.
  subroutine hirsh_symmetrize(c,xn)
    use crystalmod, only: crystal
    type(crystal), intent(in) :: c
    real*8, intent(inout) :: xn(:)

    integer :: i, k
    real*8 :: xs(c%nneq)
    integer :: ns(c%nneq)

    xs = 0d0
    ns = 0
    do i = 1, c%ncel
       k = c%atcel(i)%idx
       xs(k) = xs(k) + xn(i)
       ns(k) = ns(k) + 1
    end do
    do i = 1, c%ncel
       k = c%atcel(i)%idx
       xn(i) = xs(k) / ns(k)
    end do

  end subroutine hirsh_symmetrize

  !> Report the Hirshfeld-I populations xn(1:ncel), with ranges of
  !> charge states nmin and nmax, pinned atoms, convergence status,
  !> number of iterations, last residual, and the integrated total
  !> population in the last iteration (before clamping), ntotint.
  subroutine hirsh_report(c,xn,nmin,nmax,pinned,conv,niter,dmax,ntotint)
    use crystalmod, only: crystal
    use tools_io, only: uout, string, ioj_center, ioj_right, ferror, warning
    type(crystal), intent(in) :: c
    real*8, intent(in) :: xn(:), nmin(:), nmax(:)
    logical, intent(in) :: pinned(:)
    logical, intent(in) :: conv
    integer, intent(in) :: niter
    real*8, intent(in) :: dmax
    real*8, intent(in) :: ntotint

    ! integrated total and nuclear charge differing by more than this
    ! (a charged system differs by its charge) are reported
    real*8, parameter :: qwarn = 1.5d0

    integer :: i, iat, iz
    integer :: irep(c%nneq)
    real*8 :: ztot
    character(len=:), allocatable :: aux

    if (conv) then
       write (uout,'("+ Hirshfeld-I converged in ",A," iterations")') string(niter)
    else
       call ferror("hirsh_iterate","Hirshfeld-I not converged after " // string(niter) //&
          " iterations (max|dN| = " // string(dmax,'e',decimal=4) // ")",warning)
    end if

    irep = hirsh_irep(c)
    write (uout,'("# N_HI = Hirshfeld-I population; q_HI = Z - N_HI")')
    write (uout,'("# nmin/nmax = range of available charge states (number of electrons)")')
    write (uout,'("#nneq mult name    Z   nmin nmax      N_HI            q_HI")')
    do iat = 1, c%nneq
       i = irep(iat)
       iz = c%spc(c%at(iat)%is)%z
       aux = ""
       if (pinned(i)) aux = "  (pinned)"
       write (uout,'(8(A," "),A)') string(iat,length=4,justify=ioj_center), &
          string(c%at(iat)%mult,length=4,justify=ioj_center), &
          string(c%at(iat)%name,length=5,justify=ioj_center), &
          string(iz,length=4,justify=ioj_right), &
          string(nint(nmin(i)),length=4,justify=ioj_right), &
          string(nint(nmax(i)),length=4,justify=ioj_right), &
          string(xn(i),'f',length=16,decimal=10,justify=3), &
          string(iz-xn(i),'f',length=16,decimal=10,justify=3), aux
    end do

    ! the integrated total, not the clamped populations: clamping to
    ! the ends of the range hides a missing or double-counted core
    ztot = 0d0
    do i = 1, c%ncel
       ztot = ztot + c%spc(c%atcel(i)%is)%z
    end do
    write (uout,'("# total population: ",A,"  (nuclear charge: ",A,")")') &
       string(ntotint,'f',decimal=6), string(ztot,'f',decimal=1)
    if (any(pinned)) &
       call ferror("hirsh_iterate","Some atoms are pinned at the last available charge state",warning)
    if (ztot - ntotint > qwarn) then
       call ferror("hirsh_iterate","Total population smaller than the nuclear charge. &
          &If this is a pseudopotential valence density, use ZPSP.",warning)
    elseif (ntotint - ztot > qwarn) then
       call ferror("hirsh_iterate","Total population larger than the nuclear charge. &
          &Is ZPSP used with an all-electron density, or is the grid too coarse for the nuclear cusps?",warning)
    end if
    write (uout,*)

  end subroutine hirsh_report

end submodule proc
