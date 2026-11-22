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

! Routines for editing molecular and crystal structures
submodule (crystalmod) edit
  implicit none

contains

  !> Make a crystal seed (seed) from a crystal structure. If useabr,
  !> force the use of a particular value of useabr (1=aa and bb,
  !> 2=x2c, default is 2).
  module subroutine makeseed(c,seed,copysym,useabr,copybonding)
    use global, only: rborder_def
    use crystalseedmod, only: crystalseed
    class(crystal), intent(in) :: c
    type(crystalseed), intent(out) :: seed
    logical, intent(in) :: copysym
    integer, intent(in), optional :: useabr
    logical, intent(in), optional :: copybonding

    integer :: i
    integer :: useabr_
    logical :: copybonding_, wantsym, wantbonds, useneq

    ! initialize
    useabr_ = 2
    if (present(useabr)) useabr_ = useabr
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! What to copy: the symmetry (wantsym) and/or the bonding (wantbonds). The
    ! non-equivalent (asymmetric unit) atom list is used only when copying the
    ! symmetry without the bonding; with bonding (nstar is indexed over the
    ! complete cell list) the complete list is used instead, optionally
    ! together with the symmetry.
    wantsym = (.not.c%ismolecule .and. copysym .and. c%spgavail)
    wantbonds = (copybonding_ .and. allocated(c%nstar))
    useneq = (wantsym .and. .not.wantbonds)

    ! general
    seed%isused = .true.
    seed%name = ""
    seed%file = c%file
    seed%libname = c%libname
    seed%isformat = c%isformat

    ! atoms
    if (useneq) then
       ! the non-equivalent (asymmetric unit) atom list
       seed%nat = c%nneq
       allocate(seed%x(3,c%nneq),seed%is(c%nneq),seed%atname(c%nneq))
       if (c%haveocc) allocate(seed%occ(c%nneq))
       if (c%haveocc) allocate(seed%mix(c%nneq))
       do i = 1, c%nneq
          seed%x(:,i) = c%at(i)%x
          seed%is(i) = c%at(i)%is
          seed%atname(i) = c%at(i)%name
          if (c%haveocc) seed%occ(i) = c%at(i)%occ
          if (c%haveocc) seed%mix(i) = c%at(i)%mix
       end do
    else
       ! the complete cell list
       seed%nat = c%ncel
       allocate(seed%x(3,c%ncel),seed%is(c%ncel),seed%atname(c%ncel))
       if (c%haveocc) allocate(seed%occ(c%ncel))
       if (c%haveocc) allocate(seed%mix(c%ncel))
       do i = 1, c%ncel
          if (useabr_ == 0 .and. c%ismolecule) then
             ! molecule + useabr==0: struct_new expects absolute Cartesian (bohr)
             ! coordinates here, not crystallographic; +molx0 keeps the absolute
             ! positions invariant through the rebuild (struct_new recomputes molx0)
             seed%x(:,i) = c%atcel(i)%r + c%molx0
          else
             seed%x(:,i) = c%atcel(i)%x
          end if
          seed%is(i) = c%atcel(i)%is
          seed%atname(i) = c%at(c%atcel(i)%idx)%name
          if (c%haveocc) seed%occ(i) = c%at(c%atcel(i)%idx)%occ
          if (c%haveocc) &
             seed%mix(i) = c%at(c%atcel(i)%idx)%mix
       end do
    end if

    ! bonding (nstar is indexed by ncel)
    if (wantbonds) then
       seed%havebonds = .true.
       seed%nstar = c%nstar
    else
       seed%havebonds = .false.
       if (allocated(seed%nstar)) deallocate(seed%nstar)
    end if

    ! species
    seed%nspc = c%nspc
    allocate(seed%spc(c%nspc))
    do i = 1, c%nspc
       seed%spc(i) = c%spc(i)
    end do

    ! cell
    seed%useabr = useabr_
    if (useabr_ == 0) then
       seed%aa = c%aa
       seed%bb = c%bb
       seed%m_x2c = c%m_x2c
    elseif (useabr_ == 1) then
       seed%aa = c%aa
       seed%bb = c%bb
       seed%m_x2c = 0d0
    else
       seed%aa = 0d0
       seed%bb = 0d0
       seed%m_x2c = c%m_x2c
    end if

    ! symmetry. neqlist is .true. only for the asymmetric-unit case (useneq);
    ! a complete list carried together with the symmetry has neqlist=.false.
    seed%findsym = -1
    seed%checkrepeats = .false.
    if (wantsym) then
       seed%havesym = 1
       seed%neqlist = useneq
       seed%neqv = c%neqv
       seed%ncv = c%ncv
       allocate(seed%rotm(3,4,c%neqv),seed%cen(3,c%ncv))
       seed%rotm = c%rotm(:,:,1:c%neqv)
       seed%cen = c%cen(:,1:c%ncv)
    else
       seed%havesym = 0
       seed%neqlist = .false.
       seed%neqv = 0
       seed%ncv = 0
    end if

    ! molecular fields
    seed%ismolecule = c%ismolecule
    seed%cubic = (abs(c%aa(1)-c%aa(2)) < 1d-5).and.(abs(c%aa(1)-c%aa(3)) < 1d-5).and.&
       all(abs(c%bb-90d0) < 1d-3)
    seed%border = rborder_def
    seed%havex0 = .true.
    seed%molx0 = c%molx0

  end subroutine makeseed

  !> Make a crystal seed (seed) from a crystal structure, modulating
  !> the structure with a phonon displacement. The displacement
  !> corresponds to q-point qpt (fractional coords.), with eigenvector
  !> evec, and the given amplitude and phase. A supercell is built
  !> such that the structure is periodic according with the supplied
  !> qpt.
  module subroutine makeseed_nudged(c,seed,qpt,evec,amplitude,phase)
    use crystalseedmod, only: crystalseed
    use tools_math, only: rational_approx
    use param, only: atmass, pi, img, tpi, isformat_r_derived
    class(crystal), intent(in) :: c
    type(crystalseed), intent(out) :: seed
    real*8, intent(in) :: qpt(3)
    complex*16, intent(in) :: evec(:,:)
    real*8, intent(in) :: amplitude, phase

    real*8, parameter :: rational_approx_eps = 1d-3

    complex*16 :: displ, xdelta(3)
    integer :: i, nc(3), ix, iy, iz, k
    integer*8 :: q, r
    real*8 :: smass, xx(3), xd(3)

    ! general. This is an in-memory system derived from another one, not
    ! backed by a file: mark it so that restoring uses the initial seed.
    seed%isused = .true.
    seed%name = trim(c%origin_label()) // " (nudge)"
    seed%file = ""
    seed%isformat = isformat_r_derived

    ! calculate number of cells
    nc = 1
    do i = 1, 3
       if (abs(qpt(i)) > rational_approx_eps) then
          call rational_approx(qpt(i),q,r,rational_approx_eps)
          nc(i) = int(r)
       end if
    end do

    ! drawlist_sph(nsph)%xdelta = cmplx(xdelta1,kind=c_float_complex)
    ! x = x + real(displ * s%drawlist_sph(i)%xdelta,c_float)

    ! atoms
    displ = amplitude * exp(0.5d0 * phase * pi * img)
    k = 0
    seed%nat = c%ncel * nc(1) * nc(2) * nc(3)
    allocate(seed%x(3,seed%nat),seed%is(seed%nat),seed%atname(seed%nat))
    if (c%haveocc) allocate(seed%occ(seed%nat))
    if (c%haveocc) allocate(seed%mix(seed%nat))
    do i = 1, c%ncel
       smass = sqrt(atmass(c%spc(c%atcel(i)%is)%z))
       do ix = 0, nc(1)-1
          do iy = 0, nc(2)-1
             do iz = 0, nc(3)-1
                k = k + 1

                xx = c%atcel(i)%x + real((/ix,iy,iz/),8)
                xdelta = evec(:,i) * exp(img * tpi * dot_product(xx,qpt)) / smass
                xd = real(displ * xdelta,8)

                seed%x(:,k) = (xx + c%c2x(xd)) / real(nc,8)
                seed%is(k) = c%atcel(i)%is
                seed%atname(k) = c%at(c%atcel(i)%idx)%name
                if (c%haveocc) seed%occ(k) = c%at(c%atcel(i)%idx)%occ
                if (c%haveocc) &
                   seed%mix(k) = c%at(c%atcel(i)%idx)%mix
             end do
          end do
       end do
    end do

    ! species
    seed%nspc = c%nspc
    allocate(seed%spc(c%nspc))
    do i = 1, c%nspc
       seed%spc(i) = c%spc(i)
    end do

    ! cell
    seed%useabr = 2
    seed%aa = 0d0
    seed%bb = 0d0
    do i = 1, 3
       seed%m_x2c(:,i) = c%m_x2c(:,i) * nc(i)
    end do

    ! symmetry
    seed%findsym = -1
    seed%checkrepeats = .false.
    seed%havesym = 0
    seed%neqv = 0
    seed%ncv = 0

    ! molecular fields
    seed%ismolecule = c%ismolecule
    seed%cubic = (abs(c%aa(1)-c%aa(2)) < 1d-5).and.(abs(c%aa(1)-c%aa(3)) < 1d-5).and.&
       all(abs(c%bb-90d0) < 1d-3)
    seed%border = 0d0
    seed%havex0 = .true.
    seed%molx0 = c%molx0

  end subroutine makeseed_nudged

  !> Given a crystal structure (c) and three lattice vectors in cryst.
  !> coords (x0(:,1), x0(:,2), x0(:,3)), build the same crystal
  !> structure using the unit cell given by those vectors. If nnew,
  !> xnew, and isnew are given, replace the nnew atoms with the
  !> new positions (xnew) and species (isnew). If noenv is present
  !> and true, do not load the atomic grids or the environments in
  !> the new cell.
  module subroutine newcell(c,x00,t0,nnew,xnew,isnew,noenv,errmsg,copybonding,ti)
    use crystalseedmod, only: crystalseed
    use tools_math, only: det3, matinv, mnorm2
    use tools_io, only: string, uout
    use types, only: realloc
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x00(3,3)
    real*8, intent(in), optional :: t0(3)
    integer, intent(in), optional :: nnew
    real*8, intent(in), optional :: xnew(:,:)
    integer, intent(in), optional :: isnew(:)
    logical, intent(in), optional :: noenv
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: copybonding
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: ncseed
    logical :: ok, atgiven, issimple, copybonding_
    real*8 :: x0(3,3), x0inv(3,3), fvol, dmax0
    real*8 :: dd, t(3), xshift(3)
    integer :: i, j, k, m
    integer :: nn
    integer :: nlat, nlat2, nvec(3), ntot
    integer, allocatable :: lvec(:,:), isrc(:)
    real*8, allocatable :: xlat(:,:)
    character(len=:), allocatable :: errmsg_

    real*8, parameter :: epszero = 1d-8
    real*8, parameter :: eps = 1d-2
    real*8, parameter :: eps2 = eps * eps

    errmsg = ""
    if (c%ismolecule) then
       errmsg_ = 'NEWCELL incompatible with molecules'
       goto 999
    end if

    ! initialize
    atgiven = (present(nnew).and.present(xnew).and.present(isnew))
    x0 = x00
    dd = det3(x0)

    ! Whether to carry the bonding through the
    ! transformation. Unavailable when the new positions are given
    ! directly (atgiven).
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    copybonding_ = copybonding_ .and. allocated(c%nstar) .and. .not.atgiven

    ! check the new lattice vectors are sane
    if (abs(dd) * max(c%ncv,1) < eps) then
       errmsg_ = 'Invalid input vectors'
       goto 999
    elseif (dd < 0d0) then
       ! flip the cell
       x0 = -x0
       dd = det3(x0)
    endif

    ! set the origin translation
    if (present(t0)) then
       t = t0
    else
       t = 0d0
    end if

    ! species list for the new seed
    allocate(ncseed%spc(c%nspc))
    ncseed%nspc = c%nspc
    ncseed%spc = c%spc

    ! metrics of the new cell
    ncseed%m_x2c = matmul(c%m_x2c,x0)
    ncseed%useabr = 2

    ! check if we have new atomic positions
    if (atgiven) then
       ! sanity check
       if (size(xnew,1) /= 3 .or. size(xnew,2) /= nnew .or. size(isnew) /= nnew) then
          errmsg_ = 'NEWCELL: error in size of new atomic positions array'
          goto 999
       end if

       ! build the atomic info in the new seed from the info in input
       ncseed%nat = nnew
       allocate(ncseed%x(3,nnew),ncseed%is(nnew),ncseed%atname(nnew))
       do i = 1, nnew
          ncseed%x(:,i) = xnew(:,i)
          ncseed%is(i) = isnew(i)
          ncseed%atname(i) = c%spc(isnew(i))%name
       end do
    else
       ! check if this is a "simple" transformation (only integers in the diagonal)
       issimple = abs(x0(1,2)) < epszero.and.abs(x0(1,3)) < epszero.and.abs(x0(2,3)) < epszero.and.&
          abs(x0(2,1)) < epszero.and.abs(x0(3,1)) < epszero.and.abs(x0(3,2)) < epszero
       issimple = issimple .and. x0(1,1) > 0d0 .and. (abs(x0(1,1) - nint(x0(1,1))) < epszero)
       issimple = issimple .and. x0(2,2) > 0d0 .and. (abs(x0(2,2) - nint(x0(2,2))) < epszero)
       issimple = issimple .and. x0(3,3) > 0d0 .and. (abs(x0(3,3) - nint(x0(3,3))) < epszero)

       if (issimple) then
          ! this is a simple transformation: translate the atoms a number of times

          ! allocate
          do i = 1, 3
             nvec(i) = nint(x0(i,i))
          end do
          ntot = product(nvec)
          allocate(ncseed%x(3,c%ncel * ntot),ncseed%is(c%ncel * ntot),ncseed%atname(c%ncel * ntot))
          if (c%haveocc) allocate(ncseed%occ(c%ncel * ntot))
          if (c%haveocc) allocate(ncseed%mix(c%ncel * ntot))
          if (copybonding_) allocate(isrc(c%ncel * ntot))

          nn = 0
          do i = 1, nvec(1)
             do j = 1, nvec(2)
                do k = 1, nvec(3)
                   xshift = real((/i,j,k/) - 1,8) / real(nvec,8)
                   do m = 1, c%ncel
                      nn = nn + 1
                      ncseed%is(nn) = c%atcel(m)%is
                      ncseed%x(:,nn) = (c%atcel(m)%x-t) / real(nvec,8) + xshift
                      ncseed%atname(nn) = c%at(c%atcel(m)%idx)%name
                      if (c%haveocc) ncseed%occ(nn) = c%at(c%atcel(m)%idx)%occ
                      if (c%haveocc) &
                         ncseed%mix(nn) = c%at(c%atcel(m)%idx)%mix
                      if (copybonding_) isrc(nn) = m
                   end do
                end do
             end do
          end do
          ncseed%nat = c%ncel * ntot
       else
          ! this is not a simple transformation, so we need to search for the lattice
          ! vectors and compare the resulting atoms to avoid repetitions

          ! check that the vectors are pure translations
          do i = 1, 3
             ok = .false.
             do j = 1, c%ncv
                ok = (c%are_lclose(x0(:,i),c%cen(:,j),1d-4))
                if (ok) exit
             end do
             if (.not.ok) then
                write (uout,'("! Error: cell vectors are not lattice translations. This can happen for")')
                write (uout,'("! several reasons but a common one is that the symmetry in your crystal")')
                write (uout,'("! was not correctly calculated, or read from an external file (a cif")')
                write (uout,'("! file) that does not have all the symmetry operations. If this is the")')
                write (uout,'("! case, please try to run SYM RECALC before before attempting the cell")')
                write (uout,'("! transformation.")')
                errmsg_ = "Cell vector number " // string(i) // " is not a pure translation"
                goto 999
             end if
          end do

          ! inverse matrix
          x0inv = transpose(x0)
          call matinv(x0inv,3)

          ! calculate maximum distance
          dmax0 = 0d0
          dmax0 = max(dmax0,norm2(c%x2c(matmul((/1,0,0/),transpose(x0)))))
          dmax0 = max(dmax0,norm2(c%x2c(matmul((/0,1,0/),transpose(x0)))))
          dmax0 = max(dmax0,norm2(c%x2c(matmul((/0,0,1/),transpose(x0)))))
          dmax0 = max(dmax0,norm2(c%x2c(matmul((/1,1,0/),transpose(x0)))))
          dmax0 = max(dmax0,norm2(c%x2c(matmul((/1,0,1/),transpose(x0)))))
          dmax0 = max(dmax0,norm2(c%x2c(matmul((/0,1,1/),transpose(x0)))))
          dmax0 = max(dmax0,norm2(c%x2c(matmul((/1,1,1/),transpose(x0)))))

          ! number of atoms in the new cell
          fvol = abs(det3(x0))
          nn = nint(c%ncel * fvol)
          if (abs(nn - (c%ncel*fvol)) > eps) then
             errmsg_ = 'Inconsistent number of atoms in newcell'
             goto 999
          end if
          allocate(ncseed%x(3,nn),ncseed%is(nn),ncseed%atname(nn),isrc(nn))
          if (c%haveocc) allocate(ncseed%occ(nn))
          if (c%haveocc) allocate(ncseed%mix(nn))

          ! build the new atom list from the atoms of the current cell
          ! moved by its lattice translations, until the new cell is
          ! full. First the translations closest to the origin, which
          ! are all that is needed in a supercell.
          ncseed%nat = 0
          nlat = 0
          allocate(xlat(3,max(nint(dd),1)))
          call add_translation((/0d0,0d0,0d0/))
          if (ncseed%nat < nn) then
             call c%list_near_lattice_points((/0d0,0d0,0d0/),icrd_crys,.true.,nlat2,&
                lvec=lvec,up2d=dmax0*2.0d0,nozero=.true.)
             do i = 1, nlat2
                if (ncseed%nat >= nn) exit
                call add_translation(matmul(real(lvec(:,i),8),x0inv))
             end do
          end if

          ! Then the rest, by adding the current cell vectors to the
          ! translations found until there are no new ones. This is
          ! for a cell that is neither a supercell of the current cell
          ! nor a primitive cell of the lattice (the hexagonal cell of
          ! a cubic F lattice), which can need translations far away.
          i = 0
          do while (ncseed%nat < nn .and. i < nlat)
             i = i + 1
             do j = 1, 3
                call add_translation(xlat(:,i) + x0inv(j,:))
             end do
          end do
          call realloc(ncseed%x,3,ncseed%nat)
          call realloc(ncseed%is,ncseed%nat)
          call realloc(ncseed%atname,ncseed%nat)
          if (c%haveocc) call realloc(ncseed%occ,ncseed%nat)
          if (c%haveocc) call realloc(ncseed%mix,ncseed%nat)
          deallocate(xlat)

          ! the new cell must have all its atoms
          if (ncseed%nat /= nn) then
             errmsg_ = 'Could not find all the atoms of the new cell'
             goto 999
          end if

          ! a folding transformation (repeated atoms were skipped) does
          ! not preserve the bonds: fall back to recomputing them
          if (ncseed%nat /= c%ncel * nlat) copybonding_ = .false.
       end if
    end if

    ! map the bonds of the old cell onto the new atom list
    if (copybonding_) call map_bonds()

    ! rest of the seed information
    ncseed%isused = .true.
    ncseed%file = c%file
    ncseed%libname = c%libname
    ncseed%isformat = c%isformat
    ncseed%havesym = 0
    ncseed%findsym = -1
    ncseed%ismolecule = .false.

    ! initialize the structure
    call c%struct_new(ncseed,errmsg,noenv,ti=ti)
    if (len_trim(errmsg) > 0) return

    return

    ! error handling: report via errmsg
999 continue
    errmsg = errmsg_

  contains
    ! Add the atoms of the current cell moved by the lattice
    ! translation xl0 (crystallographic coordinates of the new cell)
    ! to the new atom list, unless the translation is equivalent to
    ! one that was already used.
    subroutine add_translation(xl0)
      real*8, intent(in) :: xl0(3)

      real*8 :: xl(3), xa(3), dxa(3)
      integer :: ia, ja
      logical :: isnew

      ! check if we have this translation already
      xl = xl0 - floor(xl0)
      do ia = 1, nlat
         dxa = xlat(:,ia) - xl
         dxa = dxa - nint(dxa)
         if (all(abs(dxa) < eps)) return
      end do
      nlat = nlat + 1
      if (nlat > size(xlat,2)) call realloc(xlat,3,2*nlat)
      xlat(:,nlat) = xl

      do ia = 1, c%ncel
         ! candidate atom
         xa = matmul(c%atcel(ia)%x-t,x0inv) + xl
         xa = xa - floor(xa)

         ! check if we have it already: a copy of the same atom by a
         ! centering vector that is a translation of the new cell
         isnew = .true.
         do ja = 1, ncseed%nat
            if (c%atcel(isrc(ja))%idx /= c%atcel(ia)%idx) cycle
            dxa = xa - ncseed%x(:,ja)
            dxa = matmul(ncseed%m_x2c,dxa - nint(dxa))
            if (dot_product(dxa,dxa) < eps2) then
               isnew = .false.
               exit
            end if
         end do
         if (.not.isnew) cycle

         ! add it to the list
         ncseed%nat = ncseed%nat + 1
         if (ncseed%nat > size(ncseed%x,2)) then
            call realloc(ncseed%x,3,2*ncseed%nat)
            call realloc(ncseed%is,2*ncseed%nat)
            call realloc(ncseed%atname,2*ncseed%nat)
            if (c%haveocc) call realloc(ncseed%occ,2*ncseed%nat)
            if (c%haveocc) call realloc(ncseed%mix,2*ncseed%nat)
            call realloc(isrc,2*ncseed%nat)
         end if
         ncseed%x(:,ncseed%nat) = xa
         ncseed%is(ncseed%nat) = c%atcel(ia)%is
         ncseed%atname(ncseed%nat) = c%at(c%atcel(ia)%idx)%name
         if (c%haveocc) ncseed%occ(ncseed%nat) = c%at(c%atcel(ia)%idx)%occ
         if (c%haveocc) ncseed%mix(ncseed%nat) = c%at(c%atcel(ia)%idx)%mix
         isrc(ncseed%nat) = ia
      end do

    end subroutine add_translation

    ! Build ncseed%nstar from c%nstar
    subroutine map_bonds()
      use types, only: neighstar
      integer :: i, k, kk, a, b, j, n
      real*8 :: xinv(3,3), w(3), dn(3)
      integer, allocatable :: idst(:), inxt(:)
      logical :: found

      real*8, parameter :: epsint = 1d-4 ! integer-vector threshold (new frac)

      ! inverse transformation (old frac to new frac)
      xinv = x0
      call matinv(xinv,3)

      ! per-source linked lists of new atoms (idst = first, inxt = next)
      allocate(idst(c%ncel),inxt(ncseed%nat))
      idst = 0
      do i = ncseed%nat, 1, -1
         inxt(i) = idst(isrc(i))
         idst(isrc(i)) = i
      end do

      allocate(ncseed%nstar(ncseed%nat))
      do i = 1, ncseed%nat
         a = isrc(i)
         n = c%nstar(a)%ncon
         ncseed%nstar(i)%isaromatic = c%nstar(a)%isaromatic
         allocate(ncseed%nstar(i)%idcon(n),ncseed%nstar(i)%lcon(3,n),&
            ncseed%nstar(i)%ordcon(n),ncseed%nstar(i)%aromdir(3,n))
         kk = 0
         do k = 1, n
            b = c%nstar(a)%idcon(k)
            w = ncseed%x(:,i) + matmul(xinv,c%atcel(b)%x +&
               real(c%nstar(a)%lcon(:,k),8) - c%atcel(a)%x)

            ! find the image of b an integer vector away from w
            found = .false.
            j = idst(b)
            do while (j > 0)
               dn = w - ncseed%x(:,j)
               if (all(abs(dn - nint(dn)) < epsint)) then
                  found = .true.
                  exit
               end if
               j = inxt(j)
            end do
            if (.not.found) cycle

            kk = kk + 1
            ncseed%nstar(i)%idcon(kk) = j
            ncseed%nstar(i)%lcon(:,kk) = nint(dn)
            ncseed%nstar(i)%ordcon(kk) = c%nstar(a)%ordcon(k)
            ncseed%nstar(i)%aromdir(:,kk) = c%nstar(a)%aromdir(:,k)
         end do
         ncseed%nstar(i)%ncon = kk
      end do
      ncseed%havebonds = .true.

    end subroutine map_bonds
  end subroutine newcell

  !> Transform to the standard cell. If toprim, convert to the
  !> primitive. If doforce, force the transformation to the primitive
  !> even if it does not lead to a smaller cell. Refine = refine the
  !> symmetry positions. noenv = do not initialize the
  !> environment. keepcell = idealize the cell and atomic positions
  !> but keep the original cell choice instead of switching to the
  !> conventional cell (only meaningful with refine). Return the
  !> transformation matrix, or a matrix of zeros if no change was
  !> done.
  module function cell_standard(c,toprim,doforce,refine,noenv,errmsg,ti,keepcell) result(x0)
    use iso_c_binding, only: c_double
    use spglib, only: spg_standardize_cell, SpglibDataset
    use global, only: symprec
    use tools_math, only: det3, matinv
    use types, only: realloc
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in) :: toprim
    logical, intent(in) :: doforce
    logical, intent(in) :: refine
    logical, intent(in), optional :: noenv
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti
    logical, intent(in), optional :: keepcell
    real*8 :: x0(3,3)

    integer :: ntyp, nat
    integer :: i, nnew, iprim, inorefine
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: types_(:)
    real*8 :: rmat(3,3), pmat(3,3)
    logical :: keepcell_
    type(SpglibDataset) :: dset

    x0 = 0d0
    errmsg = ""
    keepcell_ = .false.
    if (present(keepcell)) keepcell_ = keepcell .and. refine

    ! ignore molecules
    if (c%ismolecule) return

    ! save the transformation matrix to the conventional cell
    if (keepcell_) then
       call c%spglib_wrap(dset,.false.,errmsg,ti=ti)
       if (len_trim(errmsg) > 0) return
       pmat = transpose(dset%transformation_matrix)
    end if

    ! use spglib transformation to the standard cell
    rmat = transpose(c%m_x2c)
    nat = c%ncel
    ntyp = c%nspc
    allocate(x(3,4*c%ncel),types_(4*c%ncel)) ! to make room for spglib
    do i = 1, c%ncel
       x(:,i) = c%atcel(i)%x
       types_(i) = c%atcel(i)%is
    end do

    ! parse options
    iprim = 0
    if (toprim) iprim = 1
    inorefine = 1
    if (refine) inorefine = 0

    nnew = spg_standardize_cell(rmat,x,types_,nat,iprim,inorefine,symprec)
    if (nnew == 0) then
       errmsg = "Could not find primitive cell"
       return
    end if
    call realloc(x,3,nnew)
    call realloc(types_,nnew)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det3(rmat) < 0d0) rmat = -rmat

    ! rmat = transpose(matinv(c%spg%transformation_matrix))
    if (refine) then
       call c%newcell(rmat,nnew=nnew,xnew=x,isnew=types_,noenv=noenv,ti=ti,errmsg=errmsg)
       if (len_trim(errmsg) > 0) return
       ! transform the idealized conventional cell back to the original
       if (keepcell_) &
          call c%newcell(pmat,noenv=noenv,ti=ti,errmsg=errmsg)
    else
       ! if a primitive is wanted but det is not less than 1, do not make the change
       if (all(abs(rmat - eye) < symprec)) return
       if (toprim .and. .not.(det3(rmat) < 1d0-symprec) .and..not.doforce) return
       call c%newcell(rmat,noenv=noenv,ti=ti,errmsg=errmsg)
    end if
    if (len_trim(errmsg) > 0) return
    x0 = rmat

  end function cell_standard

  !> Transform to the Niggli cell. Return the transformation matrix,
  !> or a matrix of zeros if no change was done.
  module function cell_niggli(c,noenv,errmsg,ti) result(x0)
    use spglib, only: spg_niggli_reduce
    use global, only: symprec
    use tools_math, only: det3
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: noenv
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti
    real*8 :: x0(3,3)

    real*8 :: rmat(3,3)
    integer :: id, i

    x0 = 0d0
    errmsg = ""

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib niggli reduction
    rmat = transpose(c%m_x2c)
    id = spg_niggli_reduce(rmat,symprec)
    if (id == 0) then
       errmsg = "Could not find Niggli reduction"
       return
    end if
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det3(rmat) < 0d0) rmat = -rmat

    ! transform
    if (any(abs(rmat - eye) > symprec)) then
       call c%newcell(rmat,noenv=noenv,ti=ti,errmsg=errmsg)
       if (len_trim(errmsg) > 0) return
       x0 = rmat
    end if

  end function cell_niggli

  !> Transform to the Delaunay cell. Return the transformation matrix,
  !> or a matrix of zeros if no change was done.
  module function cell_delaunay(c,noenv,errmsg,ti) result(x0)
    use spglib, only: spg_delaunay_reduce
    use global, only: symprec
    use tools_math, only: det3
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: noenv
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti
    real*8 :: x0(3,3)

    real*8 :: rmat(3,3)
    integer :: id, i

    x0 = 0d0
    errmsg = ""

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib delaunay reduction
    rmat = transpose(c%m_x2c)
    id = spg_delaunay_reduce(rmat,symprec)
    if (id == 0) then
       errmsg = "Could not find Delaunay reduction"
       return
    end if
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det3(rmat) < 0d0) rmat = -rmat

    ! transform
    if (any(abs(rmat - eye) > symprec)) then
       call c%newcell(rmat,noenv=noenv,ti=ti,errmsg=errmsg)
       if (len_trim(errmsg) > 0) return
       x0 = rmat
    end if

  end function cell_delaunay

  !> Search the supercells of the current cell of every size from nmin
  !> (default 1) to inice. The candidates of size n are the distinct
  !> sublattices of index n, enumerated through their Hermite normal
  !> forms so that each appears exactly once (Hart and Forcade, PRB
  !> 77, 224115 (2008)). Each is reported in the basis that maximizes
  !> the radius of the inscribed sphere (its niceness): the basis
  !> whose reciprocal basis is Niggli-reduced. All of them are
  !> returned, in one flat list: the nc(n) cells of size n are
  !> cand(ic0(n)+1:ic0(n)+nc(n)), sorted by decreasing niceness (the
  !> first is the nicest cell of that size), each with its reach (half
  !> the shortest lattice vector) and the number of operations of the
  !> crystal compatible with the supercell lattice. If doreach, the
  !> cells are instead sorted by decreasing reach, ties by decreasing
  !> niceness, and each is reported in its Delaunay-reduced basis.
  module subroutine cell_nice_list(c,inice,nc,ic0,cand,errmsg,nmin,doreach)
    use spglib, only: spg_niggli_reduce, spg_delaunay_reduce
    use global, only: symprec
    use tools, only: mergesort
    use tools_math, only: idet3
    class(crystal), intent(inout) :: c
    integer, intent(in) :: inice
    integer, allocatable, intent(out) :: nc(:), ic0(:)
    type(nice_cell), allocatable, intent(out) :: cand(:)
    character(len=:), allocatable, intent(out) :: errmsg
    integer, intent(in), optional :: nmin
    logical, intent(in), optional :: doreach

    integer :: n, n0, ia, ib, ic, id, ie, if_, i, j, k, leqv, ntot, nstab
    integer :: mm(3,3), irot(3,3)
    integer, allocatable :: irotm(:,:,:), iord(:), iord2(:)
    type(nice_cell), allocatable :: tmp(:)
    logical :: reach
    real*8 :: rhead

    allocate(nc(inice),ic0(inice),cand(256))
    nc = 0
    ic0 = 0
    ntot = 0
    errmsg = ""
    if (c%ismolecule) then
       errmsg = "supercells can only be searched for crystals"
       return
    end if
    n0 = 1
    if (present(nmin)) n0 = max(nmin,1)
    reach = .false.
    if (present(doreach)) reach = doreach

    ! the distinct rotations of the crystal's current symmetry
    ! operations, in crystallographic coordinates (integer matrices);
    ! the same operations disp_count and CREATE_DISPLACEMENTS use
    allocate(irotm(3,3,c%neqv))
    leqv = 0
    do i = 1, c%neqv
       irot = nint(c%rotm(1:3,1:3,i))
       do j = 1, leqv
          if (all(irotm(:,:,j) == irot)) exit
       end do
       if (j <= leqv) cycle
       leqv = leqv + 1
       irotm(:,:,leqv) = irot
    end do

    ! enumerate the sublattices of each size
    do n = n0, inice
       ic0(n) = ntot
       ! Hermite normal forms of determinant n, with the lattice vectors
       ! as columns: (a,d,e), (0,b,f), (0,0,c) with abc = n, 0 <= d < b,
       ! 0 <= e,f < c; each sublattice of index n has exactly one, and
       ! there are sum_{abc=n} b c^2 of them
       do ia = 1, n
          if (mod(n,ia) /= 0) cycle
          do ib = 1, n / ia
             if (mod(n/ia,ib) /= 0) cycle
             ic = n / (ia * ib)
             do id = 0, ib-1
                do ie = 0, ic-1
                   do if_ = 0, ic-1
                      mm(:,1) = (/ia,id,ie/)
                      mm(:,2) = (/0,ib,if_/)
                      mm(:,3) = (/0,0,ic/)
                      if (is_orbit_rep(mm,nstab)) call add_candidate(n,mm,nstab*c%ncv*n)
                   end do
                end do
             end do
          end do
       end do

       ! sort the cells of this size by decreasing niceness or reach
       ! (stable, so ties keep the order of discovery)
       if (nc(n) > 1) then
          allocate(iord(nc(n)),tmp(nc(n)))
          do i = 1, nc(n)
             iord(i) = i
          end do
          if (reach) then
             ! by decreasing reach; equal reaches come out of the
             ! reductions with rounding noise, so each run of reaches
             ! equal within rounding is then sorted by decreasing niceness
             call mergesort(-cand(ic0(n)+1:ic0(n)+nc(n))%reach,iord,1,nc(n))
             i = 1
             do while (i < nc(n))
                rhead = cand(ic0(n)+iord(i))%reach
                j = i
                do while (j < nc(n))
                   if (abs(cand(ic0(n)+iord(j+1))%reach - rhead) > 1d-8 * rhead) exit
                   j = j + 1
                end do
                if (j > i) then
                   allocate(iord2(j-i+1))
                   do k = 1, j-i+1
                      iord2(k) = k
                   end do
                   call mergesort(-cand(ic0(n)+iord(i:j))%r,iord2,1,j-i+1)
                   iord(i:j) = iord(i-1+iord2)
                   deallocate(iord2)
                end if
                i = j + 1
             end do
          else
             call mergesort(-cand(ic0(n)+1:ic0(n)+nc(n))%r,iord,1,nc(n))
          end if
          tmp = cand(ic0(n)+1:ic0(n)+nc(n))
          do i = 1, nc(n)
             cand(ic0(n)+i) = tmp(iord(i))
          end do
          deallocate(iord,tmp)
       end if
    end do
    cand = cand(1:max(ntot,1))

  contains
    ! Is the sublattice spanned by the columns of mm the first of its
    ! orbit under the point group in the enumeration order? The rotated
    ! sublattice R L has Hermite normal form hnf(R M); if that form
    ! comes earlier in the enumeration, the orbit has already been
    ! processed. Also returns the number of rotations that give back M,
    ! the stabilizer of the lattice (the subgroup compatible with the
    ! supercell).
    function is_orbit_rep(mm,nstab)
      integer, intent(in) :: mm(3,3)
      integer, intent(out) :: nstab
      logical :: is_orbit_rep

      integer :: k, i, j, rm(3,3), hh(3,3)

      is_orbit_rep = .true.
      nstab = 0
      do k = 1, leqv
         do j = 1, 3
            do i = 1, 3
               rm(i,j) = irotm(i,1,k) * mm(1,j) + irotm(i,2,k) * mm(2,j) + irotm(i,3,k) * mm(3,j)
            end do
         end do
         call hnf(rm,hh)
         if (all(hh == mm)) then
            nstab = nstab + 1
         elseif (hnf_before(hh,mm)) then
            is_orbit_rep = .false.
            return
         end if
      end do

    end function is_orbit_rep

    ! Bring the sublattice of size n spanned by the columns of mm, with
    ! nops compatible operations, to the basis with the largest
    ! inscribed sphere, compute its niceness, and append it to the list.
    subroutine add_candidate(n,mm,nops)
      integer, intent(in) :: n, mm(3,3), nops

      integer :: ier, mred(3,3), i1, i2, i3
      real*8 :: x2c(3,3), rmat(3,3), r, scal, dmin

      ! The inscribed sphere has radius 1/(2 max|g_i|) with g_i the
      ! reciprocal vectors of the basis, so the nicest basis of the
      ! lattice is the one whose reciprocal basis is Niggli-reduced
      ! (shortest longest vector). Reciprocal vectors are scaled to the
      ! size of the direct ones for the reduction tolerance.
      x2c = matmul(c%m_x2c,real(mm,8))
      scal = (n * c%omega)**(2d0/3d0)
      rmat = transpose(recip(x2c)) * scal
      ier = spg_niggli_reduce(rmat,symprec)
      if (ier == 0) return
      r = 0.5d0 * scal / max(norm2(rmat(1,:)),norm2(rmat(2,:)),norm2(rmat(3,:)))
      x2c = recip(transpose(rmat) / scal)

      ! back to integers in the cell setting, right-handed
      rmat = matmul(c%m_c2x,x2c)
      mred = nint(rmat)
      if (maxval(abs(rmat - real(mred,8))) > 1d-6) return
      if (idet3(mred) < 0) mred = -mred

      ! The reach is half the shortest lattice vector. spglib returns
      ! the three shortest independent vectors of the Delaunay set,
      ! which contains the shortest lattice vector, so it is among
      ! the small combinations. The cell is reported in the
      ! Delaunay-reduced basis only when the cells are ranked by reach;
      ! otherwise the reduction starts from the already reduced x2c
      ! (same lattice, fewer iterations than the skewed HNF basis).
      dmin = 0d0
      if (reach) then
         rmat = transpose(matmul(c%m_x2c,real(mm,8)))
      else
         rmat = transpose(x2c)
      end if
      ier = spg_delaunay_reduce(rmat,symprec)
      if (ier == 0 .and. reach) return
      if (ier /= 0) then
         dmin = huge(1d0)
         do i1 = -1, 1
            do i2 = -1, 1
               do i3 = -1, 1
                  if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
                  dmin = min(dmin,norm2(i1 * rmat(1,:) + i2 * rmat(2,:) + i3 * rmat(3,:)))
               end do
            end do
         end do

         if (reach) then
            rmat = matmul(c%m_c2x,transpose(rmat))
            mred = nint(rmat)
            if (maxval(abs(rmat - real(mred,8))) > 1d-6) return
            if (idet3(mred) < 0) mred = -mred
         end if
      end if

      ! append
      if (ntot == size(cand,1)) then
         allocate(tmp(2*ntot))
         tmp(1:ntot) = cand(1:ntot)
         call move_alloc(tmp,cand)
      end if
      ntot = ntot + 1
      nc(n) = nc(n) + 1
      cand(ntot) = nice_cell(r=r,reach=0.5d0*dmin,m=mred,nops=nops)

    end subroutine add_candidate

    ! Hermite normal form h of the lattice spanned by the columns of the
    ! integer matrix m0 (nonzero determinant): lower triangular, positive
    ! diagonal, 0 <= h(2,1) < h(2,2), 0 <= h(3,1),h(3,2) < h(3,3). Column
    ! operations only, so h spans the same lattice.
    pure subroutine hnf(m0,h)
      integer, intent(in) :: m0(3,3)
      integer, intent(out) :: h(3,3)

      integer :: q

      h = m0
      ! zero the first row except (1,1), then the second except (2,2)
      call zero_entry(h,1,1,2)
      call zero_entry(h,1,1,3)
      call zero_entry(h,2,2,3)
      ! positive diagonal
      if (h(1,1) < 0) h(:,1) = -h(:,1)
      if (h(2,2) < 0) h(:,2) = -h(:,2)
      if (h(3,3) < 0) h(:,3) = -h(:,3)
      ! reduce the entries below the diagonal
      q = (h(2,1) - modulo(h(2,1),h(2,2))) / h(2,2)
      h(:,1) = h(:,1) - q * h(:,2)
      q = (h(3,1) - modulo(h(3,1),h(3,3))) / h(3,3)
      h(:,1) = h(:,1) - q * h(:,3)
      q = (h(3,2) - modulo(h(3,2),h(3,3))) / h(3,3)
      h(:,2) = h(:,2) - q * h(:,3)

    end subroutine hnf

    ! Euclid's algorithm on row r of h between columns i and j: on exit
    ! h(r,j) = 0 and h(r,i) = +-gcd, the lattice unchanged.
    pure subroutine zero_entry(h,r,i,j)
      integer, intent(inout) :: h(3,3)
      integer, intent(in) :: r, i, j

      integer :: q, t(3)

      do while (h(r,j) /= 0)
         q = h(r,i) / h(r,j)
         h(:,i) = h(:,i) - q * h(:,j)
         t = h(:,i)
         h(:,i) = h(:,j)
         h(:,j) = t
      end do

    end subroutine zero_entry

    ! Does the Hermite normal form h come before m in the enumeration
    ! order (a, then b, d, e, f)?
    pure function hnf_before(h,m)
      integer, intent(in) :: h(3,3), m(3,3)
      logical :: hnf_before

      integer :: kh(5), km(5), i

      kh = (/h(1,1),h(2,2),h(2,1),h(3,1),h(3,2)/)
      km = (/m(1,1),m(2,2),m(2,1),m(3,1),m(3,2)/)
      hnf_before = .false.
      do i = 1, 5
         if (kh(i) /= km(i)) then
            hnf_before = (kh(i) < km(i))
            return
         end if
      end do

    end function hnf_before

    ! Reciprocal basis (columns) of the basis a (columns), without the
    ! 2 pi: recip(recip(a)) = a.
    pure function recip(a) result(g)
      use tools_math, only: cross, det3
      real*8, intent(in) :: a(3,3)
      real*8 :: g(3,3)

      real*8 :: v

      v = det3(a)
      g(:,1) = cross(a(:,2),a(:,3)) / v
      g(:,2) = cross(a(:,3),a(:,1)) / v
      g(:,3) = cross(a(:,1),a(:,2)) / v

    end function recip
  end subroutine cell_nice_list

  !> Choose among the candidate supercells of one size given by
  !> cell_nice_list (cand, sorted by decreasing niceness; the
  !> displacement counts are filled in here as needed with
  !> disp_count). icrit = 0: the nicest cell, ties (radii equal to
  !> rounding) going to the most symmetric and then to the fewest
  !> displacements. icrit = 1 (MINDISP): the nicest cell among those
  !> with the fewest displacements of this size, ties going to the most
  !> symmetric. Returns the index of the choice (0 if there are no
  !> candidates). If doreach, the reach replaces the niceness in all of
  !> the above (cand must then come from cell_nice_list with doreach).
  !> If axes, the displacements are counted as CREATE_DISPLACEMENTS
  !> AXES generates them (six per independent atom).
  module subroutine cell_nice_select(c,cand,icrit,ibest,errmsg,doreach,axes)
    class(crystal), intent(inout) :: c
    type(nice_cell), intent(inout) :: cand(:)
    integer, intent(in) :: icrit
    integer, intent(out) :: ibest
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: doreach
    logical, intent(in), optional :: axes

    integer :: k, klast, nopsmax, smat(3,3)
    logical :: rtie
    real*8, allocatable :: key(:)

    errmsg = ""
    ibest = 0
    if (size(cand,1) < 1) return

    ! the sorting key: niceness (inscribed radius) or reach
    key = cand%r
    if (present(doreach)) then
       if (doreach) key = cand%reach
    end if

    ! the eligible candidates: all of them (MINDISP) or the nicest ones,
    ! within rounding of the first, with the most operations
    klast = size(cand,1)
    nopsmax = 0
    if (icrit /= 1) then
       klast = 1
       do while (klast < size(cand,1))
          if (abs(key(klast+1) - key(1)) > 1d-8 * key(1)) exit
          klast = klast + 1
       end do
       nopsmax = maxval(cand(1:klast)%nops)
    end if

    do k = 1, klast
       if (cand(k)%nops < nopsmax) cycle
       if (cand(k)%ndisp < 0) then
          smat = transpose(cand(k)%m)
          call c%disp_count(smat,cand(k)%ndisp,cand(k)%nindep,errmsg,axes=axes)
          if (len_trim(errmsg) > 0) return
       end if
       ! the candidates come sorted by decreasing niceness, so the
       ! earlier of two is at least as nice: a later one wins only with
       ! fewer displacements, or the same number and more operations at
       ! the same niceness (within rounding)
       if (ibest == 0) then
          ibest = k
       else
          rtie = abs(key(k) - key(ibest)) < 1d-8 * key(ibest)
          if (cand(k)%ndisp < cand(ibest)%ndisp .or. (cand(k)%ndisp == cand(ibest)%ndisp .and. &
             rtie .and. cand(k)%nops > cand(ibest)%nops)) ibest = k
       end if
    end do

  end subroutine cell_nice_select


  !> Create a new structure by reordering the atoms in the current
  !> structure. iperm is the permutation vector (atom i in the new
  !> structure is iperm(i) in the old structure).
  module subroutine reorder_atoms(c,iperm,isnneq,errmsg,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iperm(:)
    logical, intent(in) :: isnneq
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    real*8, allocatable :: x(:,:)
    integer, allocatable :: is(:)
    integer :: i
    logical :: copysym

    errmsg = ""

    ! make the new seed
    copysym = isnneq .and. .not.c%ismolecule .and. c%spgavail
    call c%makeseed(seed,copysym)
    seed0 = seed ! this seed rebuilds the current structure

    ! return if the permutation is not correct
    if (size(iperm,1) /= seed%nat) return

    ! reorder the atoms
    allocate(x(3,seed%nat),is(seed%nat))
    do i = 1, seed%nat
       x(:,i) = seed%x(:,iperm(i))
    end do
    seed%x = x
    deallocate(x,is)
    seed%is = seed%is(iperm(:))
    seed%atname = seed%atname(iperm(:))
    if (allocated(seed%occ)) seed%occ = seed%occ(iperm(:))
    if (allocated(seed%mix)) seed%mix = seed%mix(iperm(:))

    ! reload the crystal
    call c%struct_new(seed,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine reorder_atoms

  !> Create a new structure by reordering the molecular fragments in the
  !> current structure. iperm is the permutation vector (molecule i in
  !> the new structure is iperm(i) in the old structure). The atoms are
  !> relabeled so that the cell atoms are grouped by molecule and the
  !> molecules appear in the new order; molecular numbering follows the
  !> order of the atoms, so the reloaded structure has the requested
  !> molecular order.
  module subroutine reorder_molecules(c,iperm,errmsg,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iperm(:)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    real*8, allocatable :: x(:,:)
    integer, allocatable :: atperm(:)
    integer :: i, k, im, nat

    errmsg = ""

    ! the molecular decomposition must be available and the permutation valid
    if (c%nmol == 0 .or. .not.allocated(c%idatcelmol)) return
    if (size(iperm,1) /= c%nmol) return

    ! make the new seed (no symmetry, so seed atoms follow the cell-atom order)
    call c%makeseed(seed,.false.)
    seed0 = seed ! this seed rebuilds the current structure
    if (seed%nat /= c%ncel) return

    ! build the atom permutation: visit molecules in the new order and, for
    ! each, append its cell atoms keeping their original relative order
    allocate(atperm(c%ncel))
    nat = 0
    do i = 1, c%nmol
       im = iperm(i)
       if (im < 1 .or. im > c%nmol) return
       do k = 1, c%ncel
          if (c%idatcelmol(1,k) == im) then
             nat = nat + 1
             atperm(nat) = k
          end if
       end do
    end do
    if (nat /= c%ncel) return

    ! reorder the atoms
    allocate(x(3,seed%nat))
    do i = 1, seed%nat
       x(:,i) = seed%x(:,atperm(i))
    end do
    seed%x = x
    deallocate(x)
    seed%is = seed%is(atperm(:))
    seed%atname = seed%atname(atperm(:))
    if (allocated(seed%occ)) seed%occ = seed%occ(atperm(:))
    if (allocated(seed%mix)) seed%mix = seed%mix(atperm(:))

    ! reload the crystal
    call c%struct_new(seed,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine reorder_molecules

  !> Create a new structure by reordering the species in the current
  !> structure. iperm is the permutation vector (atom i in the new
  !> structure is iperm(i) in the old structure).
  module subroutine reorder_species(c,iperm,errmsg,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iperm(:)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    integer :: i, j
    integer, allocatable :: iiperm(:)

    errmsg = ""

    ! make the new seed
    call c%makeseed(seed,.false.)
    seed0 = seed ! this seed rebuilds the current structure

    ! return if the permutation is not correct
    if (size(iperm,1) /= seed%nspc) return
    allocate(iiperm(size(iperm,1)))
    do i = 1, size(iperm,1)
       iiperm(iperm(i)) = i
    end do

    ! reorder the species (remap the species indices, including the co-occupant
    ! lists, which are stored as species indices)
    seed%spc = seed%spc(iperm(:))
    do i = 1, seed%nat
       seed%is(i) = iiperm(seed%is(i))
       if (allocated(seed%mix)) then
          do j = 1, seed%mix(i)%nocc
             seed%mix(i)%is(j) = iiperm(seed%mix(i)%is(j))
          end do
       end if
    end do

    ! reload the crystal
    call c%struct_new(seed,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine reorder_species

  !> Re-assign atomic types to have an asymmetric unit with whole molecules
  module subroutine wholemols(c,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use tools, only: qcksort
    use types, only: realloc
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical, allocatable :: ismap(:,:), ldone(:)
    integer :: i, j, k, l, id
    real*8 :: x0(3)
    type(crystalseed) :: ncseed
    logical, allocatable :: sgroup(:,:), agroup(:), isuse(:)
    integer, allocatable :: pr(:,:), iorder(:), idxi(:)
    real*8 :: rot(3,4), xd(3), xxd(3)
    logical :: found, ok, again
    integer :: ig, ngroup, ngroupold

    real*8, parameter :: eps = 1d-4

    errmsg = ""

    ! this is only appropriate for molecular crystals
    if (c%ismolecule) return
    if (.not.c%ismol3d) return

    ! compute the product table for the rotational parts of this group
    allocate(pr(c%neqv,c%neqv))
    do i = 1, c%neqv
       do j = 1, c%neqv
          if (i == 1) then
             pr(i,j) = j
          else if (j == 1) then
             pr(i,j) = i
          else
             ! compute the product of i and j
             rot(:,1:3) = matmul(c%rotm(1:3,1:3,i),c%rotm(1:3,1:3,j))
             rot(:,4) = matmul(c%rotm(1:3,1:3,i),c%rotm(1:3,4,j)) + c%rotm(1:3,4,i)

             ! identify the product operation
             found = .false.
             loopk: do k = 1, c%neqv
                ok = all(abs(rot(1:3,1:3) - c%rotm(1:3,1:3,k)) < eps)
                if (ok) then
                   xd = rot(1:3,4) - c%rotm(1:3,4,k)
                   do l = 1, c%ncv
                      xxd = xd - c%cen(:,l)
                      if (all(abs(xxd - nint(xxd)) < eps)) then
                         found = .true.
                         pr(i,j) = k
                         exit loopk
                      end if
                   end do
                end if
             end do loopk
             if (.not.found) then
                errmsg = "symmetry operations do not form a group"
                return
             end if
          end if
       end do
    end do

    ! initialize the subgroup list with the trivial subgroups
    allocate(sgroup(c%neqv,2),agroup(c%neqv))
    ngroup = 2
    sgroup = .true.
    sgroup(2:,2) = .false.

    ! construct all the other subgroups
    ngroupold = 1
    do while (ngroupold < ngroup)
       ! consider the next subgroup
       ngroupold = ngroupold + 1
       ig = ngroupold

       ! try adding operations one by one
       do i = 2, c%neqv
          agroup = sgroup(:,ig)
          if (agroup(i)) cycle
          agroup(i) = .true.

          ! check if it is a subset of a known non-trivial subgroup
          found = .false.
          do j = 3, ngroup
             if (all((agroup .and. sgroup(:,j)) .eqv. agroup)) then
                found = .true.
                exit
             end if
          end do
          if (found) cycle

          ! complete the subgroup
          again = .true.
          do while (again)
             again = .false.
             do k = 1, c%neqv
                if (.not.agroup(k)) cycle
                do j = 1, c%neqv
                   if (.not.agroup(j)) cycle
                   if (.not.agroup(pr(k,j))) then
                      agroup(pr(k,j)) = .true.
                      again = .true.
                   end if
                end do
             end do
          end do

          ! check it is not the full group
          if (all(agroup)) cycle

          ! add it to the list
          ngroup = ngroup + 1
          if (ngroup > size(sgroup,2)) call realloc(sgroup,c%neqv,2*ngroup)
          sgroup(:,ngroup) = agroup
       end do
    end do
    call realloc(sgroup,c%neqv,ngroup)
    deallocate(agroup,pr)

    ! sort the subgroups by order
    allocate(iorder(ngroup),idxi(ngroup))
    do i = 1, ngroup
       iorder(i) = count(sgroup(:,i))
       idxi(i) = i
    end do
    call qcksort(iorder,idxi,1,ngroup)
    sgroup = sgroup(:,idxi)
    deallocate(iorder,idxi)

    ! identify the largest known subgroup with whole molecules in the asymmetric unit
    allocate(ismap(c%ncv,c%neqv),ldone(c%nneq))
    ig = 0
    do l = ngroup,1,-1
       ismap = .false.
       do i = 2, c%neqv
          if (.not.sgroup(i,l)) cycle
          do j = 1, c%ncv

             ldone = .false.
             do k = 1, c%ncel
                if (ldone(c%atcel(k)%idx)) cycle

                x0 = matmul(c%rotm(1:3,1:3,i),c%atcel(k)%x) + c%rotm(:,4,i) + c%cen(:,j)
                id = c%identify_atom(x0,icrd_crys)

                if (id == 0) then
                   errmsg = 'error identifying rotated atom'
                   return
                end if

                if (c%idatcelmol(1,k) == c%idatcelmol(1,id)) &
                   ismap(j,i) = .true.
                ldone(c%atcel(k)%idx) = .true.
             end do
          end do
       end do

       if (all(.not.ismap)) then
          ig = l
          exit
       end if
    end do
    deallocate(ldone,ismap)

    ! calculate the asymmetric unit for the new group, write it down in the seed
    allocate(ncseed%x(3,c%ncel),ncseed%is(c%ncel),ncseed%atname(c%ncel))
    if (c%haveocc) allocate(ncseed%occ(c%ncel))
    if (c%haveocc) allocate(ncseed%mix(c%ncel))
    ncseed%nat = 0
    allocate(isuse(c%ncel))
    isuse = .false.
    do i = 1, c%ncel
       if (isuse(i)) cycle
       isuse(i) = .true.
       do j = 1, c%neqv
          if (.not.sgroup(j,ig)) cycle
          do k = 1, c%ncv
             x0 = matmul(c%rotm(1:3,1:3,j),c%atcel(i)%x) + c%rotm(:,4,j) + c%cen(:,k)
             id = c%identify_atom(x0,icrd_crys)
             if (id == 0) then
                errmsg = 'error identifying rotated atom'
                return
             end if
             isuse(id) = .true.
          end do
       end do

       ncseed%nat = ncseed%nat + 1
       ncseed%x(:,ncseed%nat) = c%atcel(i)%x
       ncseed%is(ncseed%nat) = c%atcel(i)%is
       ncseed%atname(ncseed%nat) = c%at(c%atcel(i)%idx)%name
       if (c%haveocc) ncseed%occ(ncseed%nat) = c%at(c%atcel(i)%idx)%occ
       if (c%haveocc) &
          ncseed%mix(ncseed%nat) = c%at(c%atcel(i)%idx)%mix
    end do
    deallocate(isuse)

    ! write down the new symmetry operations in the seed
    ncseed%ncv = c%ncv
    allocate(ncseed%cen(3,c%ncv))
    ncseed%cen = c%cen(1:3,1:c%ncv)
    ncseed%neqv = count(sgroup(:,ig))
    allocate(ncseed%rotm(3,4,ncseed%neqv))
    k = 0
    do i = 1, c%neqv
       if (sgroup(i,ig)) then
          k = k + 1
          ncseed%rotm(:,:,k) = c%rotm(:,:,i)
       end if
    end do
    deallocate(sgroup)

    ! write the rest of the seed
    allocate(ncseed%spc(c%nspc))
    ncseed%nspc = c%nspc
    ncseed%spc = c%spc
    ncseed%useabr = 2
    ncseed%m_x2c = c%m_x2c
    ncseed%ismolecule = c%ismolecule
    ncseed%isused = .true.
    ncseed%file = c%file
    ncseed%libname = c%libname
    ncseed%isformat = c%isformat
    ncseed%havesym = 1
    ncseed%findsym = 0
    ncseed%checkrepeats = .false.
    ncseed%neqlist = .true.

    ! build the new crystal
    call c%struct_new(ncseed,errmsg=errmsg,ti=ti)

  end subroutine wholemols

  !> Remove or merge the atoms with IDs in the array iat(1:nat) from
  !> the structure.
  module subroutine edit_atom_list(c,nat,iat,remove,merge,duplicate,errmsg,copybonding,ti)
    use crystalseedmod, only: crystalseed
    use types, only: realloc
    class(crystal), intent(inout) :: c
    integer, intent(in) :: nat
    integer, intent(in) :: iat(nat)
    logical, intent(in), optional :: remove, merge, duplicate
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: copybonding
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    logical, allocatable :: useatoms(:), dupatoms(:)
    integer, allocatable :: imap(:)
    integer :: i, ipres, izmax, mergespc, natnew
    real*8 :: mergex(3), x0(3), xd(3), rdum
    logical :: remove_, merge_, duplicate_, copybonding_

    errmsg = ""

    ! return if nothing to do
    if (nat == 0) return

    ! check input options
    remove_ = .false.
    merge_ = .false.
    duplicate_ = .false.
    ipres = 0
    if (present(remove)) then
       remove_ = remove
       if (remove_) ipres = ipres + 1
    end if
    if (present(merge)) then
       merge_ = merge
       if (merge_) ipres = ipres + 1
    end if
    if (present(duplicate)) then
       duplicate_ = duplicate
       if (duplicate_) ipres = ipres + 1
    end if
    if (ipres == 0) return
    if (ipres > 1) then
       errmsg = 'More than one of merge/remove/duplicate'
       return
    end if

    ! carrying the bonding through the edit is only defined for removal
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    if (copybonding_ .and. .not.remove_) then
       errmsg = 'copybonding in edit_atom_list is only available with remove'
       return
    end if
    copybonding_ = copybonding_ .and. allocated(c%nstar)

    ! make seed from this crystal
    call c%makeseed(seed,copysym=.false.)
    seed0 = seed ! this seed rebuilds the current structure

    ! calculate the position of the merged atom: the heaviest atom in
    ! the list at the average position
    if (merge_) then
       ! first pass, calculate approximate average position based on atom 1
       x0 = c%atcel(iat(1))%x
       mergespc = c%atcel(iat(1))%is
       izmax = c%spc(mergespc)%z
       mergex = x0
       do i = 2, nat
          xd = c%atcel(iat(i))%x - x0
          call c%shortest(xd,rdum)
          mergex = mergex + (x0 + c%c2x(xd))
          if (c%spc(c%atcel(iat(i))%is)%z > izmax) then
             mergespc = c%atcel(iat(i))%is
             izmax = c%spc(c%atcel(iat(i))%is)%z
          end if
       end do
       mergex = mergex / real(nat,8)
    end if

    ! flag the atoms and species to use
    allocate(useatoms(c%ncel),dupatoms(c%ncel))
    useatoms = .true.
    dupatoms = .false.
    if (merge_ .or. remove_) then
       do i = 1, nat
          useatoms(iat(i)) = .false.
       end do
    else
       do i = 1, nat
          dupatoms(iat(i)) = .true.
       end do
    end if

    ! reallocate the seed arrays
    if (remove_) then
       natnew = count(useatoms)
    elseif (merge_) then
       natnew = count(useatoms) + 1
    else
       natnew = c%ncel + nat
    end if
    call realloc(seed%x,3,natnew)
    call realloc(seed%is,natnew)
    call realloc(seed%atname,natnew)
    if (allocated(seed%occ)) call realloc(seed%occ,natnew)
    ! rebuild the co-occupant lists fresh (the atom set changes here)
    if (allocated(seed%mix)) deallocate(seed%mix)
    if (c%haveocc) allocate(seed%mix(natnew))

    ! re-do the atom and species info
    seed%nat = 0
    do i = 1, c%ncel
       if (useatoms(i)) then
          seed%nat = seed%nat + 1
          seed%x(:,seed%nat) = c%atcel(i)%x
          seed%atname(seed%nat) = c%at(c%atcel(i)%idx)%name
          seed%is(seed%nat) = c%atcel(i)%is
          if (allocated(seed%occ)) seed%occ(seed%nat) = c%at(c%atcel(i)%idx)%occ
          if (c%haveocc) &
             seed%mix(seed%nat) = c%at(c%atcel(i)%idx)%mix
          if (dupatoms(i)) then
             seed%nat = seed%nat + 1
             seed%x(:,seed%nat) = c%atcel(i)%x
             seed%atname(seed%nat) = c%at(c%atcel(i)%idx)%name
             seed%is(seed%nat) = c%atcel(i)%is
             if (allocated(seed%occ)) seed%occ(seed%nat) = c%at(c%atcel(i)%idx)%occ
             if (c%haveocc) &
                seed%mix(seed%nat) = c%at(c%atcel(i)%idx)%mix
          end if
       end if
    end do

    ! if merging, add the merged atom at the end
    if (merge_) then
       seed%nat = seed%nat + 1
       seed%x(:,seed%nat) = mergex
       seed%is(seed%nat) = mergespc
       seed%atname(seed%nat) = c%spc(mergespc)%name
       if (allocated(seed%occ)) seed%occ(seed%nat) = 1d0
    end if

    ! carry the bonding of the surviving atoms through the rebuild
    ! (removal only, so the atom map is a plain compaction)
    if (copybonding_) then
       allocate(imap(c%ncel))
       natnew = 0
       do i = 1, c%ncel
          imap(i) = 0
          if (useatoms(i)) then
             natnew = natnew + 1
             imap(i) = natnew
          end if
       end do
       call c%bonds_subset(imap,seed%nstar)
       seed%havebonds = .true.
    end if

    ! build the new crystal
    call c%struct_new(seed,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine edit_atom_list

  !> Change the atoms with IDs in the array iat(1:nat) from their
  !> current species to is, and reset the atom name to the
  !> corresponding species name.
  module subroutine change_atom_species(c,nat,iat,is,copybonding,errmsg,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: nat
    integer, intent(in) :: iat(nat)
    integer, intent(in) :: is
    logical, intent(in), optional :: copybonding
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    integer :: i
    logical :: copybonding_

    errmsg = ""

    ! return if nothing to do
    if (nat == 0) return
    if (is < 1 .or. is > c%nspc) return

    ! make seed from this crystal, preserving bonding if requested
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    call c%makeseed(seed,copysym=.false.,copybonding=copybonding_)
    seed0 = seed ! this seed rebuilds the current structure

    ! apply the change; a re-speciated atom is no longer a mixed site
    do i = 1, nat
       seed%is(iat(i)) = is
       seed%atname(iat(i)) = c%spc(is)%name
       if (allocated(seed%mix)) seed%mix(iat(i))%nocc = 0
    end do

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine change_atom_species

  !> Move atom with complete list ID idx to position x in units of
  !> iunit_l (see global). If isnneq, move all atoms that are
  !> equivalent by symmetry. If dorelative, the movement is
  !> relative to its current position.
  module subroutine move_atom(c,idx,x,iunit_l,isnneq,dorelative,copybonding,errmsg,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: idx
    real*8, intent(in) :: x(3)
    integer, intent(in) :: iunit_l
    logical, intent(in) :: isnneq
    logical, intent(in) :: dorelative
    logical, intent(in), optional :: copybonding
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    real*8 :: xx(3)
    logical :: copysym, copybonding_

    errmsg = ""

    ! whether to use symmetry
    copysym = isnneq .and. .not.c%ismolecule .and. c%spgavail
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! make seed from this crystal, preserving bonding if requested
    call makeseed_for_edit(c,seed,copysym,copybonding_)
    seed0 = seed ! this seed rebuilds the current structure

    ! interpret units: internal Cartesian for a molecule, crystallographic otherwise
    xx = edit_seed_coords(c,x,iunit_l)
    if (dorelative) then
       if (c%ismolecule) then
          xx = c%atcel(idx)%r + xx
       else
          xx = seed%x(:,idx) + xx
       end if
    end if

    ! the final new position (absolute Cartesian for molecules)
    if (c%ismolecule) then
       seed%x(:,idx) = xx + c%molx0
    else
       seed%x(:,idx) = xx
    end if

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine move_atom

  !> Displace cell atom icel by dx (Cartesian, bohr) in place. This is
  !> the cheap counterpart of move_atom, meant for real-time dragging:
  !> it moves the atom and refreshes the environment, but does not
  !> rebuild the symmetry, the connectivity, or the molecular
  !> fragments, and does not wrap the atom into the main cell. The
  !> non-equivalent atom list is not touched either. The caller must
  !> restore full consistency with rebuild_after_move once the
  !> interactive edit is over.
  module subroutine move_atom_inplace(c,icel,dx)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: icel
    real*8, intent(in) :: dx(3)

    if (icel < 1 .or. icel > c%ncel) return
    call set_cellatom_r(c,icel,c%atcel(icel)%r + dx)

    ! moving a single atom is not a rigid motion of its fragment
    call atomlist_drop_molcache(c,(/icel/))
    call refresh_after_inplace_move(c)

  end subroutine move_atom_inplace

  !> Rigidly translate molecular fragment imol by dx (Cartesian, bohr)
  !> in place. Same notes as move_atom_inplace. The caller must call
  !> rebuild_after_move when the interactive edit is over.
  module subroutine move_molecule_inplace(c,imol,dx)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: imol
    real*8, intent(in) :: dx(3)

    integer :: k, icel

    if (imol < 1 .or. imol > c%nmol .or. .not.allocated(c%mol)) return
    do k = 1, c%mol(imol)%nat
       icel = c%mol(imol)%at(k)%cidx
       call set_cellatom_r(c,icel,c%atcel(icel)%r + dx)
    end do

    ! the cached standard frame rides along: the orientation does not change,
    ! only the center of mass
    if (c%mol(imol)%axes_computed) c%mol(imol)%xcm = c%mol(imol)%xcm + dx
    call refresh_after_inplace_move(c)

  end subroutine move_molecule_inplace

  !> Rigidly rotate molecular fragment imol about its center of mass in
  !> place. drot is the incremental rotation to apply (Cartesian).
  !> Same notes as move_atom_inplace. The caller must call
  !> rebuild_after_move when the interactive edit is over.
  module subroutine rotate_molecule_inplace(c,imol,drot)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: imol
    real*8, intent(in) :: drot(3,3)

    integer :: k, icel
    real*8 :: xcm(3), rnew(3), rlvec(3)

    if (imol < 1 .or. imol > c%nmol .or. .not.allocated(c%mol)) return
    if (.not.c%mol(imol)%discrete) return

    ! the rotation is applied in the fragment frame (the whole molecule, with
    ! its atoms translated by the per-atom lattice vectors); the cell atoms
    ! are the fragment positions minus those lattice vectors
    xcm = c%mol(imol)%cmass()
    do k = 1, c%mol(imol)%nat
       icel = c%mol(imol)%at(k)%cidx
       rlvec = c%x2c(real(c%mol(imol)%at(k)%lvec,8))
       rnew = xcm + matmul(drot,c%mol(imol)%at(k)%r - xcm)
       call set_cellatom_r(c,icel,rnew - rlvec)
    end do
    call c%mol(imol)%rotate_std(drot)
    call refresh_after_inplace_move(c)

  end subroutine rotate_molecule_inplace

  !> Rigidly translate the distinct cell atoms in the list idx by dx
  !> (Cartesian, bohr) in place. Same notes as move_atom_inplace. The
  !> caller must call rebuild_after_move when the interactive edit is
  !> over.
  module subroutine move_atoms_inplace(c,idx,dx)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: idx(:)
    real*8, intent(in) :: dx(3)

    integer :: k

    if (.not.atomlist_ok(c,idx)) return
    do k = 1, size(idx)
       call set_cellatom_r(c,idx(k),c%atcel(idx(k))%r + dx)
    end do
    call atomlist_drop_molcache(c,idx)
    call refresh_after_inplace_move(c)

  end subroutine move_atoms_inplace

  !> Rigidly rotate the distinct cell atoms in the list idx about
  !> their center of mass in place, each at the image given by the
  !> lattice vectors lvec (see atoms_images), not as they are in the
  !> main cell. drot is the incremental rotation to apply (Cartesian).
  !> The same lvec serves for all the steps of an interactive edit: the
  !> in-place editors do not wrap the atoms into the main cell. Same
  !> notes as move_atom_inplace. The caller must call
  !> rebuild_after_move when the interactive edit is over.
  module subroutine rotate_atoms_inplace(c,idx,lvec,drot)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: idx(:)
    integer, intent(in) :: lvec(:,:)
    real*8, intent(in) :: drot(3,3)

    integer :: k
    real*8 :: xcm(3), rlvec(3)

    if (.not.atomlist_ok(c,idx) .or. size(lvec,2) < size(idx)) return
    xcm = c%atoms_center(idx,lvec)
    do k = 1, size(idx)
       rlvec = c%x2c(real(lvec(:,k),8))
       call set_cellatom_r(c,idx(k),xcm + matmul(drot,c%atcel(idx(k))%r + rlvec - xcm) - rlvec)
    end do
    call atomlist_drop_molcache(c,idx)
    call refresh_after_inplace_move(c)

  end subroutine rotate_atoms_inplace

  !> Lattice vectors lvec(:,k) that make a compact image of the
  !> distinct cell atoms in the list idx, at c%atcel(idx(k))%x +
  !> lvec(:,k). The atoms of the list bonded to each other (c%nstar)
  !> are contiguous. A group of bonded atoms is where its first atom
  !> is in its molecular fragment if that is discrete (so a whole
  !> molecule is where its fragment is), or in the main cell
  !> otherwise; in a crystal, each group after the first is then moved
  !> to the lattice image whose center is nearest the center of the
  !> atoms placed before it, if there is a nearer one, so groups that
  !> are next to each other across a cell face stay together. Zero for
  !> an empty or invalid list.
  module subroutine atoms_images(c,idx,lvec)
    class(crystal), intent(in) :: c
    integer, intent(in) :: idx(:)
    integer, allocatable, intent(out) :: lvec(:,:)

    integer :: n, k, j, im, nat, nplaced, lv0(3), ish(3)
    integer, allocatable :: imap(:), imem(:), id(:), lv(:,:)
    logical, allocatable :: skip(:)
    logical :: disc
    real*8 :: xrun(3), xg(3)

    n = size(idx)
    allocate(lvec(3,n))
    lvec = 0
    if (.not.atomlist_ok(c,idx) .or. .not.allocated(c%nstar)) return

    ! the walks stay inside the list
    allocate(imap(c%ncel),imem(c%ncel),skip(c%ncel))
    imap = 0
    skip = .true.
    do k = 1, n
       imap(idx(k)) = k
       skip(idx(k)) = .false.
    end do

    imem = 0
    nplaced = 0
    xrun = 0d0
    do k = 1, n
       if (imem(idx(k)) /= 0) cycle

       ! the group of this atom, from its image in the fragment
       lv0 = 0
       if (allocated(c%idatcelmol) .and. allocated(c%mol)) then
          im = c%idatcelmol(1,idx(k))
          if (im >= 1 .and. im <= c%nmol) then
             if (c%mol(im)%discrete) lv0 = c%mol(im)%at(c%idatcelmol(2,idx(k)))%lvec
          end if
       end if
       call c%walk_component(idx(k),lv0,imem,nat,id,lv,disc,skipatom=skip)
       if (nat == 0) cycle

       ! next to the atoms placed before it (it stays where it is if
       ! another image is just as near)
       xg = 0d0
       do j = 1, nat
          xg = xg + c%atcel(id(j))%x + real(lv(:,j),8)
       end do
       ish = 0
       if (nplaced > 0 .and. .not.c%ismolecule) &
          ish = nint((xg / real(nat,8) - xrun / real(nplaced,8)) * (1d0 - 1d-6))
       do j = 1, nat
          lvec(:,imap(id(j))) = lv(:,j) - ish
       end do
       xrun = xrun + xg - real(nat * ish,8)
       nplaced = nplaced + nat
    end do

  end subroutine atoms_images

  !> Center of mass (Cartesian, bohr) of the distinct cell atoms in
  !> the list idx, each at the image given by the lattice vectors lvec
  !> (see atoms_images); the plain average if they have no mass.
  !> Returns zero for an empty or invalid list.
  module function atoms_center(c,idx,lvec) result(xcm)
    use param, only: atmass
    class(crystal), intent(in) :: c
    integer, intent(in) :: idx(:)
    integer, intent(in) :: lvec(:,:)
    real*8 :: xcm(3)

    integer :: k, iz
    real*8 :: wk, mtot, x(3), xsum(3)

    xcm = 0d0
    if (.not.atomlist_ok(c,idx) .or. size(lvec,2) < size(idx)) return
    xsum = 0d0
    mtot = 0d0
    do k = 1, size(idx)
       x = c%atcel(idx(k))%x + real(lvec(:,k),8)
       iz = c%spc(c%atcel(idx(k))%is)%z
       wk = 0d0
       if (iz > 0) wk = atmass(iz)
       xsum = xsum + x
       xcm = xcm + wk * x
       mtot = mtot + wk
    end do
    if (mtot > 0d0) then
       xcm = c%x2c(xcm / mtot)
    else
       xcm = c%x2c(xsum / real(size(idx),8))
    end if

  end function atoms_center

  !> Overwrite the Cartesian positions of all atoms in the main cell
  !> with rnew (3,ncel, bohr), updating the fractional coordinates
  !> accordingly. This is a fast in-place update intended for
  !> real-time animation (molecular dynamics, dragging): it does NOT
  !> rebuild the environment, symmetry, molecular fragments, or
  !> connectivity. In particular the cached fragment coordinates in
  !> c%mol go stale (unlike in the *_inplace editors, which keep them
  !> in step). Callers that need those refreshed (e.g. after stopping
  !> an animation) must rebuild the crystal explicitly with
  !> rebuild_after_move.
  module subroutine update_positions(c,rnew)
    use tools_io, only: ferror, faterr
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: rnew(:,:)

    integer :: i

    if (size(rnew,2) < c%ncel) &
       call ferror('update_positions','rnew size inconsistent with the number of cell atoms',faterr)
    do i = 1, c%ncel
       c%atcel(i)%r = rnew(:,i)
       c%atcel(i)%x = c%c2x(rnew(:,i))
    end do

  end subroutine update_positions

  !> Refresh the environment after atoms have moved via
  !> update_positions. Much cheaper than rebuild_after_move (no
  !> symmetry, bonds, or fragments); intended for the periodic
  !> neighbor-list refreshes of an interactive dynamics run.
  module subroutine update_env_after_move(c)
    class(crystal), intent(inout) :: c

    integer :: i

    do i = 1, c%ncel
       c%atcel(i)%rxc = c%x2xr(c%atcel(i)%x)
       c%atcel(i)%rxc = c%atcel(i)%rxc - floor(c%atcel(i)%rxc)
       c%atcel(i)%rxc = c%xr2c(c%atcel(i)%rxc)
    end do
    call c%build_env()

  end subroutine update_env_after_move

  !> Rebuild the crystal from the current cell positions (as left by
  !> update_positions), re-deriving the environment, connectivity,
  !> molecular fragments, and symmetry. Use this once after a run of
  !> in-place position updates (e.g. when interactive dynamics stops)
  !> to restore full consistency of the non-equivalent atom list and
  !> symmetry with the moved atoms. If copybonding, keep the current
  !> bond connectivity through the rebuild instead of re-deriving it
  !> from the atomic radii.
  module subroutine rebuild_after_move(c,copybonding,errmsg,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: copybonding
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    logical :: copybonding_

    errmsg = ""

    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! copysym=.false. builds the seed from the full cell (the moved atcel), so
    ! symmetry is re-derived from scratch (it is generally broken after motion)
    call makeseed_for_edit(c,seed,.false.,copybonding_)
    seed0 = seed ! this seed rebuilds the current structure
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine rebuild_after_move

  !> Delete the symmetry operations flagged in del and rebuild the crystal with
  !> the largest subgroup of operations that contains none of them (the identity
  !> is never removable).
  module subroutine reduce_symmetry(c,del,errmsg,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    logical, intent(in) :: del(:)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    logical, allocatable :: lkeep(:), deluse(:)
    integer :: i, k, nkeep

    errmsg = ""
    if (.not.c%isinit .or. c%ismolecule .or. c%neqv <= 1) return

    ! the identity (operation 1) is never removable
    allocate(deluse(c%neqv))
    deluse = del(1:c%neqv)
    deluse(1) = .false.
    if (.not.any(deluse)) then
       errmsg = "no removable operations were selected"
       return
    end if

    ! largest subgroup of the operations that avoids all the deleted ones
    allocate(lkeep(c%neqv))
    call c%symop_subgroup(deluse,lkeep)
    nkeep = count(lkeep)
    if (nkeep >= c%neqv) return

    ! create seed, then attach the reduced symmetry
    call c%makeseed(seed,copysym=.false.)
    seed0 = seed ! this seed rebuilds the current structure
    seed%havesym = 1
    seed%findsym = 0
    seed%checkrepeats = .false.
    seed%neqlist = .false.
    seed%neqv = nkeep
    seed%ncv = c%ncv
    if (allocated(seed%rotm)) deallocate(seed%rotm)
    if (allocated(seed%cen)) deallocate(seed%cen)
    allocate(seed%rotm(3,4,nkeep),seed%cen(3,c%ncv))
    seed%cen = c%cen(:,1:c%ncv)
    k = 0
    do i = 1, c%neqv
       if (lkeep(i)) then
          k = k + 1
          seed%rotm(:,:,k) = c%rotm(:,:,i)
       end if
    end do

    ! rebuild the crystal
    call c%struct_new(seed,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

    ! Try to recover the reduced space group
    call c%spglib_wrap(c%spg,.false.,errmsg,useidx=.true.,ti=ti)
    if (len_trim(errmsg) == 0 .and. c%spg%n_atoms > 0 .and. &
       c%spg%n_operations == c%neqv*c%ncv) then
       c%spgavail = .true.
       call c%spgtowyc()      ! clear the stale full-group letters, then fill
       call c%spgtowyc(c%spg)
    else
       errmsg = ""
       c%spgavail = .false.
       c%spg%n_atoms = 0
       call c%spgtowyc()
    end if

  end subroutine reduce_symmetry

  !> Recalculate the unit cell of a molecule to fit the current atomic
  !> positions, growing it as needed so the molecule (e.g. after a
  !> rigid rotation or translation in place) stays inside the
  !> cell.
  module subroutine recompute_molecular_cell(c)
    use tools_math, only: m_x2c_from_cellpar, matinv
    use global, only: rborder_def
    class(crystal), intent(inout) :: c

    integer :: i, j, k
    real*8 :: smin(3), smax(3), border, sh(3), saux(3)

    if (.not.c%ismolecule .or. c%ncel == 0) return

    ! calculate the bounding box
    border = max(rborder_def,1d-6)
    smin = c%atcel(1)%r + c%molx0
    smax = smin
    do k = 2, c%ncel
       saux = c%atcel(k)%r + c%molx0
       do j = 1, 3
          smin(j) = min(smin(j),saux(j))
          smax(j) = max(smax(j),saux(j))
       end do
    end do
    smin = smin - border
    smax = smax + border

    ! set up the new orthogonal cell
    sh = c%molx0 - smin
    c%aa = smax - smin
    c%bb = 90d0
    c%m_x2c = m_x2c_from_cellpar(c%aa,c%bb)
    c%m_c2x = c%m_x2c
    call matinv(c%m_c2x,3)
    c%molx0 = smin

    ! refresh the derived cell metrics and WS reduced cell for the new box so
    ! omega/gtensor/reduced matrices stay consistent
    call set_cell_metrics(c)

    ! molecular cell used for display (same formula as struct_new)
    c%molborder = max(border - max(2d0,0.8d0 * border),0d0) / c%aa

    ! shift the cell atoms
    do k = 1, c%ncel
       c%atcel(k)%r = c%atcel(k)%r + sh
       c%atcel(k)%x = c%c2x(c%atcel(k)%r)
       i = c%atcel(k)%idx
       if (i >= 1 .and. i <= c%nneq) then
          c%at(i)%r = c%atcel(k)%r
          c%at(i)%x = c%atcel(k)%x
       end if
    end do

    ! shift the cached molecular fragments to keep them consistent
    if (allocated(c%mol)) then
       do i = 1, c%nmol
          c%mol(i)%m_x2c = c%m_x2c
          do j = 1, c%mol(i)%nat
             c%mol(i)%at(j)%r = c%mol(i)%at(j)%r + sh
             c%mol(i)%at(j)%x = c%c2x(c%mol(i)%at(j)%r)
          end do
          if (c%mol(i)%axes_computed) c%mol(i)%xcm = c%mol(i)%xcm + sh
       end do
    end if

    ! rebuild the block-grid neighbor environment for the new cell/positions
    ! (through update_env_after_move, so the reduced-cell coordinates the
    ! neighbor queries use are refreshed along with the block hash)
    call c%update_env_after_move()

    ! recompute the vacuum descriptors for the new cell
    call c%calc_vacuum_lengths()

  end subroutine recompute_molecular_cell

  !> Compute the derived cell metrics (reciprocal matrices, metric tensors,
  !> volume, norms) and the Wigner-Seitz reduced cell from the already-set
  !> c%m_x2c / c%m_c2x. If errmsg is present, a failure in the
  !> Wigner-Seitz cell is returned in it instead of stopping.
  module subroutine set_cell_metrics(c,errmsg)
    use tools, only: wscell
    use tools_math, only: matinv, det3, mnorm2
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out), optional :: errmsg

    integer :: i
    real*8 :: g(3,3)
    real*8, allocatable :: area(:)

    g = matmul(transpose(c%m_x2c),c%m_x2c)
    c%m_rc2rx = transpose(c%m_x2c)
    c%m_rx2rc = transpose(c%m_c2x)
    c%gtensor = g
    c%omega = sqrt(max(det3(g),0d0))
    c%grtensor = g
    call matinv(c%grtensor,3)
    do i = 1, 3
       c%ar(i) = sqrt(c%grtensor(i,i))
    end do
    c%n2_x2c = mnorm2(c%m_x2c)
    c%n2_c2x = mnorm2(c%m_c2x)

    call wscell(c%m_x2c,.not.c%ismolecule,c%ws_nf,c%ws_nv,c%ws_mnfv,c%ws_iside,c%ws_nside,&
       c%ws_x,c%ws_ineighc,c%ws_ineighx,area,c%isortho,c%m_xr2x,c%m_x2xr,c%m_xr2c,c%m_c2xr,&
       c%n2_xr2x,c%n2_x2xr,c%n2_xr2c,c%n2_c2xr,c%ws_ineighxr,c%isortho_del,errmsg)

  end subroutine set_cell_metrics

  !> Rigidly translate the molecular fragment imol so that its center
  !> of mass is at position x in units of iunit_l (see global). If
  !> dorelative, x is interpreted as a displacement of the center of
  !> mass relative to its current position. All atoms in the fragment
  !> are translated by the same vector; symmetry is not preserved (P1).
  module subroutine move_molecule(c,imol,x,iunit_l,dorelative,copybonding,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use global, only: iunit_ang, iunit_bohr
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    integer, intent(in) :: imol
    real*8, intent(in) :: x(3)
    integer, intent(in) :: iunit_l
    logical, intent(in) :: dorelative
    logical, intent(in), optional :: copybonding
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    real*8 :: xx(3), dx(3), dxc(3)
    integer :: k
    logical :: copybonding_

    errmsg = ""

    ! consistency checks
    if (imol < 1 .or. imol > c%nmol .or. .not.allocated(c%idatcelmol)) return
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! interpret units: target/displacement in fractional coordinates
    if (iunit_l == iunit_ang) then
       xx = c%c2x(x / bohrtoa)
    elseif (iunit_l == iunit_bohr) then
       xx = c%c2x(x)
    else
       xx = x
    end if

    ! displacement (fractional) to apply to every atom in the fragment
    if (dorelative) then
       dx = xx
    else
       dx = xx - c%c2x(c%mol(imol)%cmass())
    end if

    ! make seed from this crystal (no symmetry, so seed atoms follow cell-atom
    ! order), preserving bonding if requested; re-fits the molecular cell
    call makeseed_for_edit(c,seed,.false.,copybonding_)
    seed0 = seed ! this seed rebuilds the current structure
    if (seed%nat /= c%ncel) return

    ! translate the atoms belonging to this fragment (Cartesian for molecules,
    ! crystallographic otherwise, matching the seed coordinate convention)
    if (c%ismolecule) then
       dxc = c%x2c(dx)
       do k = 1, c%ncel
          if (c%idatcelmol(1,k) == imol) &
             seed%x(:,k) = seed%x(:,k) + dxc
       end do
    else
       do k = 1, c%ncel
          if (c%idatcelmol(1,k) == imol) &
             seed%x(:,k) = seed%x(:,k) + dx
       end do
    end if

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine move_molecule

  !> Rigidly rotate the molecular fragment imol about its center of mass
  !> so that its standard orientation becomes the target orientation. The
  !> target is given by exactly one of: euler (ZYZ Euler angles, radians),
  !> quat (unit quaternion w,x,y,z), or rmat (3x3 rotation matrix). The
  !> rotation uses the fragment atoms (whole molecule, with lattice
  !> translations applied), not the cell atoms, so it stays rigid even
  !> when the molecule is split across cell boundaries. Only applies to
  !> discrete fragments. Symmetry is not preserved (P1).
  module subroutine rotate_molecule(c,imol,euler,quat,rmat,copybonding,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use tools_math, only: euler2mat, mat2euler, quat2mat
    use tools_io, only: ferror, faterr
    class(crystal), intent(inout) :: c
    integer, intent(in) :: imol
    real*8, intent(in), optional :: euler(3)
    real*8, intent(in), optional :: quat(4)
    real*8, intent(in), optional :: rmat(3,3)
    logical, intent(in), optional :: copybonding
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    real*8 :: rnew(3), amat(3,3), rrot(3,3), xcm(3)
    integer :: j, k, lvec(3), npres
    logical :: copybonding_

    errmsg = ""

    ! consistency checks
    if (imol < 1 .or. imol > c%nmol .or. .not.allocated(c%idatcelmol)) return
    if (.not.c%mol(imol)%discrete) return
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! exactly one orientation representation must be given
    npres = 0
    if (present(euler)) npres = npres + 1
    if (present(quat)) npres = npres + 1
    if (present(rmat)) npres = npres + 1
    if (npres /= 1) &
       call ferror('rotate_molecule','exactly one of euler/quat/rmat is required',faterr)

    ! make sure the standard frame is available
    if (.not.c%mol(imol)%axes_computed) call c%mol(imol)%compute_std()

    ! target orientation as a rotation matrix, then the rigid rotation that
    ! brings the current standard frame to that target
    if (present(euler)) then
       amat = euler2mat(euler)
    elseif (present(quat)) then
       amat = quat2mat(quat)
    else
       amat = rmat
    end if
    rrot = matmul(amat,transpose(c%mol(imol)%m_std))
    xcm = c%mol(imol)%xcm

    ! make seed from this crystal (no symmetry, so seed atoms follow cell-atom
    ! order), preserving bonding if requested; re-fits the molecular cell
    call makeseed_for_edit(c,seed,.false.,copybonding_)
    seed0 = seed ! this seed rebuilds the current structure
    if (seed%nat /= c%ncel) return

    ! rotate the fragment atoms (whole molecule) and write them to the seed
    do j = 1, c%mol(imol)%nat
       k = c%mol(imol)%at(j)%cidx
       ! atom position rigidly rotated about the molecule's center of mass
       rnew = xcm + matmul(rrot,c%mol(imol)%at(j)%r - xcm)
       if (c%ismolecule) then
          seed%x(:,k) = rnew + c%molx0 ! absolute cartesian (+molx0)
       else
          lvec = nint(seed%x(:,k) - c%c2x(c%mol(imol)%at(j)%r))
          seed%x(:,k) = c%c2x(rnew) + lvec ! crystallographic
       end if
    end do

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

    ! set the Euler angles explicitly (useful for continuous drag in the GUI)
    ! only if the bonding is preserved - otherwise the molecules may change
    if (copybonding_) then
       call c%mol(imol)%compute_std()
       c%mol(imol)%euler_std = mat2euler(amat)
    end if

  end subroutine rotate_molecule

  !> Modify the unit cell by changing the parameter given by iaxis:
  !> 1=a, 2=b, 3=c, -1=alpha, -2=beta, -3=gamma, 0=volume. The
  !> parameter is changed to x in units of iunit_l. If dorelative, the
  !> change is relative to the current value. If dofraction, x is
  !> interpreted as the fractional change in the current value.
  module subroutine move_cell(c,iaxis,x,iunit_l,dorelative,dofraction,copybonding,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use global, only: iunit_ang
    use param, only: bohrtoa, third
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iaxis
    real*8, intent(in) :: x
    integer, intent(in) :: iunit_l
    logical, intent(in) :: dorelative, dofraction
    logical, intent(in), optional :: copybonding
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    real*8 :: xx, ref
    logical :: copybonding_

    errmsg = ""

    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! make seed from this crystal, preserving the connectivity if requested
    if (iaxis /= 0) then
       call c%makeseed(seed,copysym=.false.,useabr=1,copybonding=copybonding_)
       seed0 = seed ! this seed rebuilds the current structure
    else
       call c%makeseed(seed,copysym=.false.,useabr=2,copybonding=copybonding_)
    end if

    ! interpret units
    if (iunit_l == iunit_ang.and..not.dofraction) then
       if (iaxis /= 0) then
          xx = x / bohrtoa
       else
          xx = x / bohrtoa**3
       end if
    else
       xx = x
    end if

    ! get the reference value
    if (iaxis == 1) then
       ref = c%aa(1)
    elseif (iaxis == 2) then
       ref = c%aa(2)
    elseif (iaxis == 3) then
       ref = c%aa(3)
    elseif (iaxis == -1) then
       ref = c%bb(1)
    elseif (iaxis == -2) then
       ref = c%bb(2)
    elseif (iaxis == -3) then
       ref = c%bb(3)
    elseif (iaxis == 0) then
       ref = c%omega
    else
       return
    end if

    ! transform
    if (dofraction) then
       ref = ref * xx
    elseif (dorelative) then
       ref = ref + xx
    else
       ref = xx
    end if

    ! put back the value
    if (iaxis == 1) then
       seed%aa(1) = ref
    elseif (iaxis == 2) then
       seed%aa(2) = ref
    elseif (iaxis == 3) then
       seed%aa(3) = ref
    elseif (iaxis == -1) then
       seed%bb(1) = ref
    elseif (iaxis == -2) then
       seed%bb(2) = ref
    elseif (iaxis == -3) then
       seed%bb(3) = ref
    elseif (iaxis == 0) then
       seed%m_x2c = seed%m_x2c * (ref / c%omega)**third
    else
       return
    end if

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine move_cell

  !> Modify the unit cell by changing the cell lengths and axes to the
  !> given values.
  module subroutine move_cell_all(c,aa,bb,copybonding,errmsg,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: aa(3), bb(3)
    logical, intent(in), optional :: copybonding
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    logical :: copybonding_

    errmsg = ""

    ! make seed from this crystal, preserving bonding if requested
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    call c%makeseed(seed,copysym=.false.,useabr=1,copybonding=copybonding_)
    seed0 = seed ! this seed rebuilds the current structure

    ! set the new axes
    seed%aa = aa
    seed%bb = bb

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

  end subroutine move_cell_all

  !> Find a unit cell with one of its axes oriented relative to the
  !> lattice. If isdir, the axis iaxis (1, 2 or 3) of the new cell is
  !> the shortest lattice vector along the direction ivec = [uvw], and
  !> the other two are the shortest lattice vectors that complete the
  !> cell (the shortest that can be an axis next to the first, then the
  !> shortest that goes with both; among vectors of the same length, the
  !> one that allows the shortest third axis and then the closest to
  !> perpendicular to the direction). Otherwise, ivec = (hkl) are the
  !> Miller indices of a family of lattice planes: the other two axes of
  !> the new cell are the shortest lattice vectors in the plane and the
  !> axis iaxis is the shortest lattice vector perpendicular to the
  !> plane, on the side the (hkl) normal points to. If the lattice has
  !> no vector perpendicular to the plane, it crosses too many planes
  !> (nplmax) or a primitive cell is requested (prim), it is the lattice
  !> vector that goes from one plane to the next, as close to the normal
  !> as possible. ivec is referred to the current cell. The search uses
  !> the whole lattice (centering vectors included), so the new cell is
  !> primitive, except if the axis perpendicular to a plane crosses
  !> several planes (never, if prim). On output, x0 has the lattice
  !> vectors of the new cell in the columns, in crystallographic
  !> coordinates of the current cell, ready to be passed to newcell.
  module subroutine cell_orient(c,ivec,isdir,iaxis,prim,x0,errmsg)
    use tools, only: primitive_lattice
    use tools_math, only: matinv, gcd, det3, iadj3
    use param, only: eye
    class(crystal), intent(in) :: c
    integer, intent(in) :: ivec(3)
    logical, intent(in) :: isdir
    integer, intent(in) :: iaxis
    logical, intent(in) :: prim
    real*8, intent(out) :: x0(3,3)
    character(len=:), allocatable, intent(out) :: errmsg

    real*8, parameter :: epsint = 1d-4 ! integer threshold
    real*8, parameter :: epslen = 1d-10 ! relative threshold for vector lengths
    ! largest index in the primitive basis: the integer products of
    ! three of them must not overflow
    integer, parameter :: ivmax = 500
    ! most planes crossed by the axis perpendicular to a plane
    integer, parameter :: nplmax = 100
    ! most lattice vectors of the same length considered for an axis
    integer, parameter :: mtie = 12

    real*8 :: mp(3,3), p(3,3), pinv(3,3), g(3,3), gp(3,3), ga(3), aga
    real*8 :: w(3), ru(3,3), det2, alpha, beta, q, qbest, q1, q2
    integer :: iv(3), u(3,3), ub(3,2), iw(3), ig, i, j, i0, j0, ibest, jbest
    integer :: nv, nw, xv(3,mtie), xw(3,mtie)
    logical :: ok

    errmsg = ""
    x0 = eye
    if (c%ismolecule) then
       errmsg = "The cell of a molecule cannot be transformed"
       return
    end if
    if (iaxis < 1 .or. iaxis > 3) then
       errmsg = "Invalid axis"
       return
    end if
    if (all(ivec == 0)) then
       errmsg = "The indices cannot all be zero"
       return
    end if

    ! a primitive basis of the lattice: Cartesian (mp) and in
    ! crystallographic coordinates of the current cell (p), and its
    ! metric tensor (g)
    call primitive_lattice(c%m_x2c,c%ncv,c%cen,mp,ok)
    if (.not.ok) mp = c%m_x2c
    p = matmul(c%m_c2x,mp)
    pinv = p
    call matinv(pinv,3)
    g = matmul(transpose(mp),mp)

    ! the direction or the plane, as a primitive integer vector in
    ! the primitive basis (direct or reciprocal). The current lattice
    ! vectors are integer combinations of the primitive ones; the
    ! products of hkl with the primitive vectors are multiples of 1/ncv
    if (isdir) then
       w = matmul(pinv,real(ivec,8))
    else
       w = matmul(transpose(p),real(ivec,8)) * real(c%ncv,8)
    end if
    if (any(abs(w) > 0.5d0 * real(huge(1),8))) then
       errmsg = "The indices are too large"
       return
    end if
    if (any(abs(w - nint(w)) > epsint)) then
       errmsg = "Inconsistent centering vectors: recalculate the symmetry"
       return
    end if
    iv = nint(w)
    ig = gcd(abs(iv),3)
    iv = iv / ig
    if (any(abs(iv) > ivmax)) then
       errmsg = "The indices are too large"
       return
    end if

    ! integer matrix with determinant +1 and the last column equal to
    ! iv: the last axis is the direction
    call complete_basis(iv,u)

    if (isdir) then
       ! metric of the lattice projected on the plane perpendicular to the direction
       ga = matmul(g,real(u(:,3),8))
       aga = dot_product(real(u(:,3),8),ga)
       do i = 1, 3
          do j = 1, 3
             gp(i,j) = g(i,j) - ga(i) * ga(j) / aga
          end do
       end do

       ! a first pair: the shortest projections, each as close to
       ! perpendicular to the direction as possible
       call gauss_reduce(u(:,1),u(:,2),gp)
       do i = 1, 2
          u(:,i) = u(:,i) - nint(dot_product(real(u(:,i),8),ga) / aga) * u(:,3)
       end do

       ! The shortest lattice vectors that complete the cell are not
       ! longer than those. Of the shortest vectors that can be the
       ! first, take the one with the shortest partner; if no partner
       ! is within reach, look farther. The first pair stays if the
       ! search is too long (very large indices).
       ub = u(:,1:2)
       q1 = qform(u(:,1),u(:,1),g)
       q2 = qform(u(:,2),u(:,2),g)
       call shortest_basis_vectors(u,1,g,min(q1,q2),nv,xv)
       ok = .false.
       q = max(q1,q2)
       do j = 1, 4
          qbest = q
          do i = 1, nv
             u(:,1) = xv(:,i)
             call shortest_basis_vectors(u,2,g,qbest,nw,xw)
             if (nw == 0) cycle
             q2 = qform(xw(:,1),xw(:,1),g)
             if (ok .and. q2 > qbest * (1d0 - epslen)) cycle
             ok = .true.
             qbest = q2
             ub(:,1) = xv(:,i)
             ub(:,2) = xw(:,1)
          end do
          if (ok .or. nv == 0) exit
          q = 4d0 * q
       end do
       u(:,1:2) = ub

       ! the two at an obtuse angle, and at acute angles with the
       ! direction if possible
       if (qform(u(:,1),u(:,2),g) > epslen * aga) u(:,2) = -u(:,2)
       if (dot_product(real(u(:,1) + u(:,2),8),ga) < -epslen * aga) u(:,1:2) = -u(:,1:2)
    else
       ! the first two axes are in the plane and the last goes to the
       ! next plane, for now: the transpose of the inverse
       u = transpose(iadj3(u))

       ! shortest vectors in the plane, at an obtuse angle
       call gauss_reduce(u(:,1),u(:,2),g)

       ! shortest vector to the next plane, from the metric of the
       ! three vectors (the third plus alpha and beta times the first
       ! two is perpendicular to the plane)
       ru = real(u,8)
       gp = matmul(transpose(ru),matmul(g,ru))
       det2 = gp(1,1) * gp(2,2) - gp(1,2) * gp(1,2)
       alpha = -(gp(1,3) * gp(2,2) - gp(2,3) * gp(1,2)) / det2
       beta = -(gp(2,3) * gp(1,1) - gp(1,3) * gp(1,2)) / det2
       i0 = nint(alpha)
       j0 = nint(beta)
       qbest = huge(1d0)
       ibest = i0
       jbest = j0
       do i = i0-1, i0+1
          do j = j0-1, j0+1
             iw = u(:,3) + i * u(:,1) + j * u(:,2)
             q = qform(iw,iw,g)
             if (q < qbest * (1d0 - epslen)) then
                qbest = q
                ibest = i
                jbest = j
             end if
          end do
       end do
       u(:,3) = u(:,3) + ibest * u(:,1) + jbest * u(:,2)
       alpha = alpha - ibest
       beta = beta - jbest

       ! the shortest lattice vector perpendicular to the plane, if a
       ! multiple of the perpendicular vector is in the lattice and
       ! the cell does not have to be primitive
       if (.not.prim) then
          do i = 1, nplmax
             if (abs(i * alpha - nint(i * alpha)) < epsint .and. abs(i * beta - nint(i * beta)) < epsint) then
                u(:,3) = i * u(:,3) + nint(i * alpha) * u(:,1) + nint(i * beta) * u(:,2)
                exit
             end if
          end do
       end if
    end if

    ! the three vectors in crystallographic coordinates of the current
    ! cell; exchange the first two if they are left-handed (the
    ! primitive basis can be)
    ru = matmul(p,real(u,8))
    if (det3(ru) < 0d0) ru(:,1:2) = ru(:,2:1:-1)

    ! the axes of the new cell, with the last vector in iaxis
    do i = 1, 3
       j = mod(iaxis + i - 1,3) + 1
       x0(:,j) = ru(:,i)
    end do

  end subroutine cell_orient

  !> Find the translation of the origin of the cell given by itype:
  !> origin_atom (the cell atom id), origin_cm (the center of mass),
  !> origin_center (the point that puts the atoms in the middle of the
  !> cell, along the axes that have vacuum) or origin_inversion (the
  !> inversion center closest to the current origin). Along the axes
  !> that have vacuum, the center of mass uses the atoms gathered on
  !> one side of it, so it does not depend on whether the atoms
  !> straddle the cell boundary; along the other axes, it uses the
  !> atoms as they are in the cell. Returns the new origin (t0) in
  !> crystallographic coordinates, to be passed to newcell.
  module subroutine cell_origin(c,itype,id,t0,errmsg)
    use param, only: atmass, eye
    class(crystal), intent(in) :: c
    integer, intent(in) :: itype
    integer, intent(in) :: id
    real*8, intent(out) :: t0(3)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, j, ia, l1, l2, l3
    real*8 :: mass, mtot, x(3), xc(3), xs(3), d, dbest

    real*8, parameter :: eps = 1d-5 ! comparison of rotation matrices
    real*8, parameter :: epsvac = 1d-6 ! an atom on the cell face, as in calc_vacuum_lengths (bohr)

    errmsg = ""
    t0 = 0d0
    if (c%ismolecule) then
       errmsg = "The origin of a molecule cannot be moved"
       return
    end if

    if (itype == origin_atom) then
       if (id < 1 .or. id > c%ncel) then
          errmsg = "Atom not found"
          return
       end if
       t0 = c%atcel(id)%x
    elseif (itype == origin_center) then
       if (.not.any(c%vaclength > iperiod_vacthr)) then
          errmsg = "There is no vacuum along any of the axes"
          return
       end if
       ! with the origin in the middle of the vacuum, the atoms are in
       ! the middle of the cell
       where (c%vaclength > iperiod_vacthr) t0 = 0.5d0 * (c%vacbot + c%vactop)
    elseif (itype == origin_cm) then
       mtot = 0d0
       do i = 1, c%ncel
          ! along the axes with vacuum, the atom in the period that
          ! starts at the top of the vacuum
          x = c%atcel(i)%x
          do ia = 1, 3
             if (c%vaclength(ia) > iperiod_vacthr) &
                x(ia) = x(ia) - floor(x(ia) - c%vactop(ia) + epsvac / c%aa(ia))
          end do
          mass = atmass(c%spc(c%atcel(i)%is)%z)
          mtot = mtot + mass
          t0 = t0 + mass * x
       end do
       t0 = t0 / max(mtot,1d-40)
    elseif (itype == origin_inversion) then
       ! an operation -x + t has inversion centers at (t + lattice vector)/2
       if (c%havesym < 1) then
          errmsg = "The symmetry of the crystal is not available"
          return
       end if
       dbest = huge(1d0)
       do i = 1, c%neqv
          if (any(abs(c%rotm(1:3,1:3,i) + eye) > eps)) cycle
          do j = 1, c%ncv
             x = c%rotm(1:3,4,i) + c%cen(:,j)
             do l1 = 0, 1
                do l2 = 0, 1
                   do l3 = 0, 1
                      ! the image of the center closest to the origin
                      xc = 0.5d0 * (x + real((/l1,l2,l3/),8))
                      xs = xc
                      call c%shortest(xs,d)
                      if (d < dbest - eps) then
                         dbest = d
                         t0 = xc + nint(c%c2x(xs) - xc)
                      end if
                   end do
                end do
             end do
          end do
          exit
       end do
       if (dbest == huge(1d0)) errmsg = "The space group has no inversion center"
    else
       errmsg = "Unknown origin type"
    end if

  end subroutine cell_origin

  !> Calculate the vacuum and the thickness of the contents of the
  !> cell along the crystallographic axis iaxis (1, 2 or 3), seen as a
  !> slab parallel to the other two axes. The cell is cut at the
  !> fractional coordinate cut along the axis or, if cut is not
  !> present, at the largest empty gap along it; with keepmol, the
  !> molecules are kept whole and the cut divides their centers. If
  !> cutat is present and a cell atom, the cell is cut at that atom
  !> instead (with keepmol, at the center of its molecule), which is
  !> then at the bottom of the slab.
  !> Returns the distance between the two faces of the slab (thick)
  !> and between the slab and its next periodic image (vac), both
  !> measured along the perpendicular to the other two axes, in bohr.
  !> With keepmol, vac is negative if the molecules on both sides of
  !> the cut interlock. If present, cutpos is the fractional
  !> coordinate at which the cell was cut, and ubot is the fractional
  !> coordinate along the axis of the lowest atom of the slab as it is
  !> laid out (see vacuum_layout).
  module subroutine vacuum_info(c,iaxis,keepmol,vac,thick,cut,cutpos,ubot,cutat)
    class(crystal), intent(in) :: c
    integer, intent(in) :: iaxis
    logical, intent(in) :: keepmol
    real*8, intent(out) :: vac, thick
    real*8, intent(in), optional :: cut
    real*8, intent(out), optional :: cutpos
    real*8, intent(out), optional :: ubot
    integer, intent(in), optional :: cutat

    integer, allocatable :: lsh(:)
    real*8 :: umin, umax, h

    vac = 0d0
    thick = 0d0
    if (present(cutpos)) cutpos = 0d0
    if (present(ubot)) ubot = 0d0
    if (c%ismolecule .or. c%ncel == 0 .or. iaxis < 1 .or. iaxis > 3) return

    ! height of the cell along the perpendicular to the other two axes
    h = 1d0 / c%ar(iaxis)
    call vacuum_layout(c,iaxis,keepmol,lsh,umin,umax,cut,cutpos,cutat)
    thick = (umax - umin) * h
    vac = h - thick
    if (present(ubot)) ubot = umin

  end subroutine vacuum_info

  !> Set the vacuum along the crystallographic axis iaxis (1, 2 or 3)
  !> to vac (bohr). The contents of the cell are cut at the fractional
  !> coordinate cut along the axis or, if cut is not present, at the
  !> largest empty gap along it (see vacuum_info), and the cell is
  !> lengthened or shortened along that axis until the distance
  !> between the slab and its next periodic image, measured along the
  !> perpendicular to the other two axes, is vac. The positions of the
  !> atoms in the slab relative to each other do not change. ipos is
  !> the position of the slab in the new cell: vacpos_center (in the
  !> middle of the cell) or vacpos_bottom (lowest atom at fractional
  !> coordinate zero, all the vacuum above). The new axis has the
  !> direction of the old one or, if orthogonal, it is perpendicular to
  !> the other two. If keepmol, the molecules are kept whole. If
  !> copybonding, keep the bond connectivity except for the bonds that
  !> cross the cut, which are kept only if the two sides do not move
  !> relative to each other (the same crystal with another origin). If
  !> the cell and the position of the slab in it are already the ones
  !> asked for, the crystal is left as it is; changed, if present, says
  !> whether it was rebuilt. (The outline of the new cell in the
  !> geometry window of the GUI, vacuum_newcell, follows the same rules
  !> for the axis and the position of the slab.)
  module subroutine set_vacuum(c,iaxis,vac,ipos,orthogonal,keepmol,cut,copybonding,errmsg,changed,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iaxis
    real*8, intent(in) :: vac
    integer, intent(in) :: ipos
    logical, intent(in) :: orthogonal, keepmol
    real*8, intent(in), optional :: cut
    logical, intent(in), optional :: copybonding
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(out), optional :: changed
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    logical :: copybonding_, samecell
    integer, allocatable :: lsh(:), imap(:), lshift(:,:)
    integer :: i, idrop
    real*8 :: umin, umax, h, hnew, z0, u, q, xn(3), tilt(3), x(3)

    ! shortest cell height allowed (bohr): the rebuild of a cell much
    ! flatter than an atom does not finish in a reasonable time
    real*8, parameter :: hmin = 1d0
    real*8, parameter :: epstilt = 1d-12 ! the axis is already perpendicular (fractional)
    real*8, parameter :: epssame = 1d-8 ! the height of the cell does not change (bohr)

    errmsg = ""
    if (present(changed)) changed = .false.
    if (c%ismolecule) then
       errmsg = "The vacuum cannot be set in a molecule"
       return
    end if
    if (iaxis < 1 .or. iaxis > 3) then
       errmsg = "Invalid axis"
       return
    end if
    if (ipos /= vacpos_center .and. ipos /= vacpos_bottom) then
       errmsg = "Invalid slab position"
       return
    end if
    if (vac < 0d0 .and. .not.keepmol) then
       errmsg = "The vacuum cannot be negative"
       return
    end if
    if (c%ncel == 0) return

    ! the unit vector perpendicular to the other two axes (along the
    ! reciprocal lattice vector) and the height of the cell along it
    xn = c%m_c2x(iaxis,:) / c%ar(iaxis)
    h = 1d0 / c%ar(iaxis)

    ! lay out the atoms as a slab
    call vacuum_layout(c,iaxis,keepmol,lsh,umin,umax,cut)

    ! the height of the new cell
    hnew = (umax - umin) * h + vac
    if (hnew < hmin) then
       errmsg = "The vacuum is too short: the cell would be less than 1 bohr high along this axis"
       return
    end if

    ! height of the bottom of the slab in the new cell
    if (ipos == vacpos_bottom) then
       z0 = 0d0
    else
       z0 = 0.5d0 * vac
    end if

    ! the tilt of the axis: the old lattice vector is h * xn minus
    ! tilt(j) times the other two
    tilt = 0d0
    if (orthogonal) then
       tilt = matmul(c%m_c2x,h * xn)
       tilt(iaxis) = 0d0
       where (abs(tilt) < epstilt) tilt = 0d0
    end if

    ! nothing to do if the cell does not change and the slab is where
    ! it is asked to be (the atoms would move by lattice vectors)
    samecell = (abs(hnew - h) < epssame .and. all(tilt == 0d0))
    q = umin - z0 / h
    if (samecell .and. abs(q - nint(q)) * h < epssame) return

    ! make seed from this crystal, preserving the connectivity if requested
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    call c%makeseed(seed,copysym=.false.,useabr=2,copybonding=copybonding_)
    seed0 = seed ! this seed rebuilds the current structure

    ! the new lattice vector. If it keeps its direction, so do the
    ! fractional coordinates along the other two axes. If it is made
    ! perpendicular to them, the atoms keep their Cartesian positions:
    ! an atom at height u goes with the tilt of the old axis.
    if (orthogonal) then
       seed%m_x2c(:,iaxis) = hnew * xn
    else
       seed%m_x2c(:,iaxis) = c%m_x2c(:,iaxis) * (hnew / h)
    end if

    ! the atomic positions: the height of an atom over the plane of
    ! the other two axes is its fractional coordinate times h
    do i = 1, c%ncel
       x = c%atcel(i)%x
       u = x(iaxis) + real(lsh(i),8) - umin
       x = x - tilt * u
       x(iaxis) = (u * h + z0) / hnew
       seed%x(:,i) = x
    end do

    ! the bonds: refer them to the shifted atoms and drop those across the cut
    if (seed%havebonds) then
       allocate(imap(c%ncel),lshift(3,c%ncel))
       lshift = 0
       do i = 1, c%ncel
          imap(i) = i
          lshift(iaxis,i) = lsh(i)
       end do
       idrop = iaxis
       if (samecell) idrop = 0
       call c%bonds_subset(imap,seed%nstar,lshift,dropaxis=idrop)
    end if

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if
    if (present(changed)) changed = .true.

  end subroutine set_vacuum

  !> Add atom with species is and position x in units of iunit_l (see
  !> global). If is <= 0, add a new species with Z = abs(is) to the
  !> system. If isnneq, replicate the atom by symmetry.
  !> Add one atom of species is at x (in units iunit_l; isnneq marks it
  !> as a new non-equivalent atom). copybonding keeps the bonding of the
  !> structure. The added atom is bonded to whatever a distance search
  !> finds, unless the caller names the atom it bonds to (bondto, a cell
  !> atom, at the position bondx if it is one of its periodic images):
  !> then that is the only bond made, and nothing else is re-bonded.
  module subroutine add_atom(c,is,x,iunit_l,isnneq,copybonding,bondto,bondx,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use types, only: realloc, siteocc
    use tools_io, only: nameguess
    class(crystal), intent(inout) :: c
    integer, intent(in) :: is
    real*8, intent(in) :: x(3)
    integer, intent(in) :: iunit_l
    logical, intent(in) :: isnneq
    logical, intent(in), optional :: copybonding
    integer, intent(in), optional :: bondto
    real*8, intent(in), optional :: bondx(3)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    type(siteocc), allocatable :: mixaux(:)
    real*8 :: xx(3)
    logical :: copysym, copybonding_
    integer :: is_, lvec(3)

    errmsg = ""

    ! whether to use symmetry; carrying the bonding requires the
    ! one-to-one atom list a symmetric add does not keep
    copysym = isnneq .and. .not.c%ismolecule .and. c%spgavail
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    copybonding_ = copybonding_ .and. .not.copysym .and. allocated(c%nstar)

    ! make seed from this crystal; re-fits the molecular cell to the added atom
    call makeseed_for_edit(c,seed,copysym,copybonding_)
    seed0 = seed ! this seed rebuilds the current structure

    ! interpret units: internal Cartesian for a molecule, crystallographic otherwise
    xx = edit_seed_coords(c,x,iunit_l)

    ! add the species
    if (is <= 0) then
       seed%nspc = seed%nspc + 1
       call realloc(seed%spc,seed%nspc)
       seed%spc(seed%nspc)%name = nameguess(abs(is),.true.)
       seed%spc(seed%nspc)%z = abs(is)
       is_ = seed%nspc
    else
       is_ = is
    end if

    ! add the atom
    seed%nat = seed%nat + 1
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    call realloc(seed%atname,seed%nat)
    if (copybonding_) call realloc(seed%nstar,seed%nat)
    if (allocated(seed%occ)) then
       call realloc(seed%occ,seed%nat)
       seed%occ(seed%nat) = 1d0
    end if
    ! the mixed-site data is indexed per atom by struct_new: grow it
    if (allocated(seed%mix)) then
       allocate(mixaux(seed%nat))
       mixaux(1:size(seed%mix,1)) = seed%mix
       call move_alloc(mixaux,seed%mix)
    end if
    ! the new atom position (absolute Cartesian for molecules)
    if (c%ismolecule) then
       seed%x(:,seed%nat) = xx + c%molx0
    else
       seed%x(:,seed%nat) = xx
    end if
    seed%is(seed%nat) = is_
    seed%atname(seed%nat) = seed%spc(is_)%name

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

    ! the bonds of the added atom: the one the caller asked for, or the
    ! ones a distance search finds
    if (copybonding_ .and. present(bondto)) then
       if (bondto >= 1 .and. bondto < c%ncel) then
          ! the caller may bond to an image, and the new atom may have
          ! been wrapped on the rebuild: recover the lattice vector from
          ! the offset the caller placed it against
          lvec = 0
          if (.not.c%ismolecule .and. present(bondx)) &
             lvec = nint((c%c2x(bondx) - xx) -&
             (c%atcel(bondto)%x - c%atcel(c%ncel)%x))
          call c%add_bond(c%ncel,bondto,lvec,1)
       end if
    elseif (copybonding_) then
       call bonds_attach_new(c,c%ncel-1,newnew=.true.,newold=.true.)
    end if

  end subroutine add_atom

  !> Add nat atoms with atomic numbers zat and Cartesian coordinates x
  !> to the crystal in a single rebuild.
  module subroutine add_fragment(c,nat,zat,x,copybonding,nstar0,newbonds,newbondx,errmsg,ti)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: nat
    integer, intent(in) :: zat(nat)
    real*8, intent(in) :: x(3,nat)
    logical, intent(in), optional :: copybonding
    type(neighstar), intent(in), optional :: nstar0(nat)
    integer, intent(in), optional :: newbonds(:,:)
    real*8, intent(in), optional :: newbondx(:,:)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    call c%replace_fragment(0,(/integer::/),nat,zat,x,copybonding,nstar0,&
       newbonds=newbonds,newbondx=newbondx,errmsg=errmsg,ti=ti)

  end subroutine add_fragment

  !> Delete the ndel cell atoms in idel and add nadd atoms with atomic
  !> numbers zat and Cartesian coordinates x to the crystal, in a
  !> single rebuild. With copybonding, the surviving atoms keep their
  !> bonds and the added atoms bring their internal bonding in nstar0.
  !> The bonds attaching the added atoms to the rest come from a distance
  !> search, unless newbonds is given: newbonds(:,i) = (added atom,
  !> pre-edit index of the atom it bonds to) is then the complete list of
  !> them, and an empty newbonds attaches the fragment to nothing.
  module subroutine replace_fragment(c,ndel,idel,nadd,zat,x,copybonding,nstar0,&
     newbonds,newbondx,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use types, only: realloc, siteocc, star_find
    use tools_io, only: nameguess
    class(crystal), intent(inout) :: c
    integer, intent(in) :: ndel
    integer, intent(in) :: idel(ndel)
    integer, intent(in) :: nadd
    integer, intent(in) :: zat(nadd)
    real*8, intent(in) :: x(3,nadd)
    logical, intent(in), optional :: copybonding
    type(neighstar), intent(in), optional :: nstar0(nadd)
    integer, intent(in), optional :: newbonds(:,:)
    real*8, intent(in), optional :: newbondx(:,:)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    type(crystalseed) :: seed0 ! pre-edit snapshot, for rollback
    character(len=:), allocatable :: errmsg0
    type(siteocc), allocatable :: mixaux(:)
    logical, allocatable :: keep(:)
    integer, allocatable :: imap(:)
    integer :: i, j, k, is, nat0, ia, ib, nb
    integer :: lvec(3)
    logical :: copybonding_, newnew, newold

    errmsg = ""

    if (ndel <= 0 .and. nadd <= 0) return
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    copybonding_ = copybonding_ .and. allocated(c%nstar)

    ! make seed from this crystal (atoms in atcel order); re-fits the
    ! molecular cell to the new fragment
    call makeseed_for_edit(c,seed,.false.,.false.)
    seed0 = seed ! this seed rebuilds the current structure

    ! the surviving atoms; they occupy 1..count(keep), in order, of the
    ! new atom list
    allocate(keep(seed%nat))
    keep = .true.
    do i = 1, ndel
       if (idel(i) >= 1 .and. idel(i) <= seed%nat) keep(idel(i)) = .false.
    end do

    ! Bonding of the surviving atoms. The added atoms start with their
    ! given internal bonding (nstar0, renumbered) or with no bonds;
    ! either way their bonds to the rest are found after the rebuild,
    ! unless the caller gave them (newbonds). imap (the map from the old
    ! numbering to the new one) is kept for the latter.
    if (copybonding_) then
       allocate(imap(seed%nat))
       nat0 = 0
       do i = 1, seed%nat
          imap(i) = 0
          if (keep(i)) then
             nat0 = nat0 + 1
             imap(i) = nat0
          end if
       end do
       call c%bonds_subset(imap,seed%nstar)
       call realloc(seed%nstar,nat0+nadd)
       if (present(nstar0)) then
          do i = 1, nadd
             seed%nstar(nat0+i) = nstar0(i)
             do k = 1, seed%nstar(nat0+i)%ncon
                seed%nstar(nat0+i)%idcon(k) = seed%nstar(nat0+i)%idcon(k) + nat0
             end do
          end do
       end if
       seed%havebonds = .true.
    end if

    ! compact out the deleted atoms
    if (ndel > 0) then
       nat0 = 0
       do i = 1, seed%nat
          if (.not.keep(i)) cycle
          nat0 = nat0 + 1
          if (nat0 == i) cycle
          seed%x(:,nat0) = seed%x(:,i)
          seed%is(nat0) = seed%is(i)
          seed%atname(nat0) = seed%atname(i)
          if (allocated(seed%occ)) seed%occ(nat0) = seed%occ(i)
          if (allocated(seed%mix)) seed%mix(nat0) = seed%mix(i)
       end do
       seed%nat = nat0
    end if

    ! refuse to build an empty structure
    if (seed%nat + nadd <= 0) return

    ! grow the atom arrays; mixed-site data must be grown too
    nat0 = seed%nat
    seed%nat = nat0 + nadd
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    call realloc(seed%atname,seed%nat)
    if (allocated(seed%occ)) then
       call realloc(seed%occ,seed%nat)
       seed%occ(nat0+1:seed%nat) = 1d0
    end if
    if (allocated(seed%mix)) then
       allocate(mixaux(seed%nat))
       mixaux(1:nat0) = seed%mix(1:nat0)
       call move_alloc(mixaux,seed%mix)
    end if

    do i = 1, nadd
       ! find or create the species for this atomic number
       is = 0
       do j = 1, seed%nspc
          if (seed%spc(j)%z == zat(i)) then
             is = j
             exit
          end if
       end do
       if (is == 0) then
          seed%nspc = seed%nspc + 1
          call realloc(seed%spc,seed%nspc)
          seed%spc(seed%nspc)%name = nameguess(zat(i),.true.)
          seed%spc(seed%nspc)%z = zat(i)
          is = seed%nspc
       end if

       ! the atom (absolute Cartesian for molecules)
       if (c%ismolecule) then
          seed%x(:,nat0+i) = x(:,i) + c%molx0
       else
          seed%x(:,nat0+i) = c%c2x(x(:,i))
       end if
       seed%is(nat0+i) = is
       seed%atname(nat0+i) = seed%spc(is)%name
    end do

    ! build the new crystal
    call c%struct_new(seed,errmsg=errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       ! the rebuild failed and struct_new leaves the crystal wiped:
       ! restore the pre-edit structure so the edit is a no-op
       call c%struct_new(seed0,errmsg0,ti=ti)
       return
    end if

    ! Bonds of the added atoms. Among themselves, a distance search finds
    ! them unless they brought their own (nstar0); to the rest, likewise
    ! unless the caller gave them (newbonds).
    if (copybonding_) then
       nat0 = seed%nat - nadd
       newnew = .not.present(nstar0)
       newold = .not.present(newbonds)
       if (newnew .or. newold) call bonds_attach_new(c,nat0,newnew,newold)
       if (present(newbonds)) then
          nb = 0
          do i = 1, size(newbonds,2)
             ia = newbonds(1,i)
             ib = newbonds(2,i)
             if (ia < 1 .or. ia > nadd .or. ib < 1 .or. ib > size(imap,1)) cycle
             if (imap(ib) == 0) cycle
             ! the atom the caller bonds to may be an image of the one in
             ! the cell, and both atoms may have been wrapped on the
             ! rebuild: recover the lattice vector by comparing the offset
             ! the caller placed the fragment against with the stored one
             lvec = 0
             if (.not.c%ismolecule .and. present(newbondx)) &
                lvec = nint(c%c2x(newbondx(:,i) - x(:,ia)) -&
                (c%atcel(imap(ib))%x - c%atcel(nat0+ia)%x))
             ia = nat0 + ia
             ib = imap(ib)
             if (ia == ib .and. all(lvec == 0)) cycle
             if (star_find(c%nstar(ia),ib,lvec) > 0) cycle
             call star_append(c%nstar(ia),ib,lvec,1)
             call star_append(c%nstar(ib),ia,-lvec,1)
             nb = nb + 1
          end do
          ! one pass over the whole structure, not one per bond
          if (nb > 0) call c%refresh_molecular_data()
       end if
    end if

  end subroutine replace_fragment

  ! Remove the bond between cell atoms iat1 and iat2 (iat2 at lattice vector
  ! lvec relative to iat1) by editing the neighbor stars in place. Removes both
  ! reciprocal entries (iat2 at lvec from iat1, iat1 at -lvec from iat2) and
  ! recomputes the molecular-fragment data; does NOT rebuild the structure.
  !> Index of the entry of iat1's neighbor star pointing at iat2 with
  !> lattice vector lvec, or 0 if there is no such bond (or the
  !> connectivity is not available).
  pure module function find_bond(c,iat1,iat2,lvec)
    use types, only: star_find
    class(crystal), intent(in) :: c
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)
    integer :: find_bond

    find_bond = 0
    if (.not.allocated(c%nstar)) return
    if (iat1 < 1 .or. iat1 > c%ncel .or. iat2 < 1 .or. iat2 > c%ncel) return
    find_bond = star_find(c%nstar(iat1),iat2,lvec)

  end function find_bond

  module subroutine remove_bond(c,iat1,iat2,lvec)
    use types, only: neighstar, star_find
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)

    if (.not.allocated(c%nstar)) return
    if (iat1 < 1 .or. iat1 > c%ncel .or. iat2 < 1 .or. iat2 > c%ncel) return

    ! remove iat2 (at lvec) from iat1's star and the reciprocal iat1 (at -lvec)
    ! from iat2's star
    call remove_one(c%nstar(iat1),iat2,lvec)
    call remove_one(c%nstar(iat2),iat1,-lvec)

    ! keep molecular fragments consistent with the new connectivity (in-place
    ! edit, no asterism recompute)
    call c%refresh_molecular_data()

  contains
    subroutine remove_one(ns,id,lv)
      type(neighstar), intent(inout) :: ns
      integer, intent(in) :: id, lv(3)
      integer :: k
      k = star_find(ns,id,lv)
      do while (k > 0)
         ns%idcon(k:ns%ncon-1) = ns%idcon(k+1:ns%ncon)
         ns%lcon(:,k:ns%ncon-1) = ns%lcon(:,k+1:ns%ncon)
         ns%ordcon(k:ns%ncon-1) = ns%ordcon(k+1:ns%ncon)
         ns%aromdir(:,k:ns%ncon-1) = ns%aromdir(:,k+1:ns%ncon)
         ns%ncon = ns%ncon - 1
         k = star_find(ns,id,lv)
      end do
    end subroutine remove_one
  end subroutine remove_bond

  ! Set the bond order (ordcon) of the bond between cell atoms iat1 and iat2
  ! (iat2 at lattice vector lvec relative to iat1). Updates both reciprocal
  ! entries in the neighbor stars in place; does not change the connectivity.
  ! Order convention: 0=dashed, 1=single, 2=double, 3=triple.
  module subroutine set_bond_order(c,iat1,iat2,lvec,order)
    use types, only: star_find
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)
    integer, intent(in) :: order

    integer :: k

    if (.not.allocated(c%nstar)) return
    if (iat1 < 1 .or. iat1 > c%ncel .or. iat2 < 1 .or. iat2 > c%ncel) return

    k = star_find(c%nstar(iat1),iat2,lvec)
    if (k > 0) c%nstar(iat1)%ordcon(k) = order
    k = star_find(c%nstar(iat2),iat1,-lvec)
    if (k > 0) c%nstar(iat2)%ordcon(k) = order

  end subroutine set_bond_order

  ! Add a bond between cell atoms iat1 and iat2 + lvec with the given
  ! bond order by editing the neighbor stars in place. Adds both
  ! reciprocal entries (iat2 at lvec from iat1, iat1 at -lvec from
  ! iat2) and recomputes the molecular-fragment data; does NOT rebuild
  ! the structure. Self-bonds and already-existing bonds are ignored.
  module subroutine add_bond(c,iat1,iat2,lvec,order)
    use types, only: star_find
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)
    integer, intent(in) :: order

    if (iat1 < 1 .or. iat1 > c%ncel .or. iat2 < 1 .or. iat2 > c%ncel) return
    if (iat1 == iat2 .and. all(lvec == 0)) return
    if (.not.allocated(c%nstar)) allocate(c%nstar(c%ncel))

    ! ignore the bond if it is already there
    if (star_find(c%nstar(iat1),iat2,lvec) > 0) return

    ! add iat2 (at lvec) to iat1's star and the reciprocal iat1 (at -lvec)
    ! to iat2's star
    call star_append(c%nstar(iat1),iat2,lvec,order)
    call star_append(c%nstar(iat2),iat1,-lvec,order)

    ! keep molecular fragments consistent with the new connectivity (in-place
    ! edit, no asterism recompute)
    call c%refresh_molecular_data()

  end subroutine add_bond

  !> Subset and renumber the neighbor stars (c%nstar) under the atom
  !> index map imap(1:ncel): imap(i) is the new index of cell atom i, or
  !> 0 if it is dropped. The result, in nstar(1:count(imap>0)), keeps
  !> only the bonds between two kept atoms, with idcon renumbered. If
  !> lshift(3,ncel) is given, atom i is understood to move by the
  !> lattice vector lshift(:,i) (e.g. a molecule unwrap) and the lcon
  !> vectors are corrected accordingly. With dropaxis (1-3; any other
  !> value is ignored), the bonds whose corrected lattice vector is
  !> nonzero along that axis are dropped. c%nstar must be allocated.
  module subroutine bonds_subset(c,imap,nstar,lshift,dropaxis)
    use types, only: neighstar, nstar_subset
    class(crystal), intent(in) :: c
    integer, intent(in) :: imap(:)
    type(neighstar), allocatable, intent(out) :: nstar(:)
    integer, intent(in), optional :: lshift(:,:)
    integer, intent(in), optional :: dropaxis

    call nstar_subset(c%nstar(1:c%ncel),imap,nstar,lshift,dropaxis=dropaxis)

  end subroutine bonds_subset

  !> Classified list of the bonded neighbors of cell atom icel, one
  !> entry per bond in the neighbor star (a neighbor may appear once
  !> per periodic image). Entries at degenerate distances (< 1d-10
  !> bohr) are skipped. Returns nsub = 0 if the connectivity is not
  !> available or icel is out of range; sub is (re)allocated to at
  !> least nsub.
  module subroutine substituents(c,icel,nsub,sub)
    class(crystal), intent(in) :: c
    integer, intent(in) :: icel
    integer, intent(out) :: nsub
    type(substituent), allocatable, intent(inout) :: sub(:)

    integer :: k, nb
    real*8 :: tv(3), v(3), d

    real*8, parameter :: eps_degen = 1d-10 ! coincident-atom threshold (bohr)

    nsub = 0
    if (.not.allocated(c%nstar)) return
    if (icel < 1 .or. icel > c%ncel) return

    if (allocated(sub)) then
       if (size(sub,1) < c%nstar(icel)%ncon) deallocate(sub)
    end if
    if (.not.allocated(sub)) allocate(sub(max(c%nstar(icel)%ncon,1)))

    do k = 1, c%nstar(icel)%ncon
       nb = c%nstar(icel)%idcon(k)
       tv = c%x2c(dble(c%nstar(icel)%lcon(:,k)))
       v = c%atcel(nb)%r + tv - c%atcel(icel)%r
       d = norm2(v)
       if (d < eps_degen) cycle
       nsub = nsub + 1
       sub(nsub)%id = nb
       sub(nsub)%kstar = k
       sub(nsub)%z = c%spc(c%atcel(nb)%is)%z
       sub(nsub)%terminal = (c%nstar(nb)%ncon == 1)
       sub(nsub)%u = v / d
       sub(nsub)%dist = d
       sub(nsub)%tv = tv
    end do

  end subroutine substituents

  !> Recompute the atomic connectivity (asterisms) and the derived
  !> molecular data (fragments, molecular equivalence, periodicity).
  module subroutine rebond(c,atmrad,bondfac,bonddelta,allowed)
    use param, only: maxzat0
    class(crystal), intent(inout) :: c
    real*8, intent(in), optional :: atmrad(0:maxzat0)
    real*8, intent(in), optional :: bondfac
    real*8, intent(in), optional :: bonddelta
    logical, intent(in), optional :: allowed(:,:)

    call c%find_asterisms(c%nstar,atmrad,bondfac,bonddelta=bonddelta,allowed=allowed)
    call c%refresh_molecular_data()

  end subroutine rebond

  !> Refresh the molecular data derived from the current connectivity: molecular
  !> fragments, molecular equivalence, and periodicity. This is the post-asterism
  !> tail shared by rebond, struct_new (when bonds are copied, not recomputed),
  !> and the in-place bond editors (add_bond/remove_bond). It assumes c%nstar is
  !> already populated; it does not recompute the asterisms.
  module subroutine refresh_molecular_data(c)
    class(crystal), intent(inout) :: c

    call c%fill_molecular_fragments()
    call c%calculate_molecular_equivalence()
    call c%calculate_periodicity()

  end subroutine refresh_molecular_data

  !xx! private procedures

  ! Lay out the atoms of the cell as a slab parallel to the other two
  ! axes than iaxis. The units of the layout are the atoms or, with
  ! keepmol, the molecules that are not periodic along the axis, taken
  ! whole in their connected placement. The cell is cut at the
  ! fractional coordinate cut along the
  ! axis: the units with their center below it go on top of the
  ! others, and a unit at the cut stays at the bottom. If cut is not
  ! present, the cell is cut in the middle of the largest gap between
  ! the centers of the units (the first one, if several have the same
  ! length, as in calc_vacuum_lengths). If cutat is present and a cell
  ! atom, the cell is cut at the center of the unit of that atom
  ! instead, whatever cut is. On output, lsh(i) is the
  ! number of lattice translations along the axis that put the cell
  ! atom i in the slab (up to a translation of the whole slab), umin
  ! and umax are the lowest and highest fractional coordinates of the
  ! shifted atoms and, if present, cutpos is the fractional coordinate
  ! of the cut, between 0 and 1.
  subroutine vacuum_layout(c,iaxis,keepmol,lsh,umin,umax,cut,cutpos,cutat)
    use tools, only: qcksort
    class(crystal), intent(in) :: c
    integer, intent(in) :: iaxis
    logical, intent(in) :: keepmol
    integer, allocatable, intent(out) :: lsh(:)
    real*8, intent(out) :: umin, umax
    real*8, intent(in), optional :: cut
    real*8, intent(out), optional :: cutpos
    integer, intent(in), optional :: cutat

    ! same criterion as calc_vacuum_lengths (bohr)
    real*8, parameter :: eps = 1d-6

    integer :: i, j, k, im, nunit, jbest
    integer, allocatable :: iunit(:), iord(:), nat(:), ish(:)
    real*8, allocatable :: ucen(:)
    real*8 :: u, epsf, gap, gapbest, ucut
    logical :: whole, atcut

    allocate(lsh(c%ncel))
    lsh = 0
    umin = 0d0
    umax = 0d0
    if (present(cutpos)) cutpos = 0d0
    if (c%ncel == 0) return
    epsf = eps / c%aa(iaxis)

    ! assign the atoms to units; the molecules are taken whole, in
    ! their connected placement (the lattice vectors of the fragment
    ! atoms do not gather the periodic core of a layer or a chain)
    allocate(iunit(c%ncel))
    iunit = 0
    nunit = 0
    if (keepmol .and. allocated(c%mol) .and. allocated(c%lvecmolc)) then
       do im = 1, c%nmol
          whole = c%mol(im)%discrete
          if (.not.whole) then
             whole = .true.
             if (allocated(c%mol(im)%lvec)) &
                whole = all(c%mol(im)%lvec(iaxis,1:c%mol(im)%nlvec) == 0)
          end if
          if (.not.whole) cycle
          nunit = nunit + 1
          do k = 1, c%mol(im)%nat
             i = c%mol(im)%at(k)%cidx
             iunit(i) = nunit
             lsh(i) = c%lvecmolc(iaxis,i)
          end do
       end do
    end if
    do i = 1, c%ncel
       if (iunit(i) == 0) then
          nunit = nunit + 1
          iunit(i) = nunit
       end if
    end do

    ! centers of the units, and the translations (ish) that bring them
    ! to the cell
    allocate(ucen(nunit),nat(nunit),ish(nunit))
    ucen = 0d0
    nat = 0
    do i = 1, c%ncel
       ucen(iunit(i)) = ucen(iunit(i)) + c%atcel(i)%x(iaxis) + real(lsh(i),8)
       nat(iunit(i)) = nat(iunit(i)) + 1
    end do
    ucen = ucen / real(nat,8)
    do j = 1, nunit
       ish(j) = -floor(ucen(j))
       if (ucen(j) + real(ish(j),8) > 1d0 - epsf) ish(j) = ish(j) - 1
       ucen(j) = ucen(j) + real(ish(j),8)
    end do

    ! the position of the cut
    atcut = .false.
    if (present(cutat)) atcut = (cutat >= 1 .and. cutat <= c%ncel)
    if (atcut) then
       ucut = max(ucen(iunit(cutat)),0d0)
    elseif (present(cut)) then
       ucut = cut - floor(cut)
       if (ucut > 1d0 - epsf) ucut = 0d0
    else
       ! the middle of the largest gap between consecutive centers; the
       ! last gap wraps around
       allocate(iord(nunit))
       do j = 1, nunit
          iord(j) = j
       end do
       call qcksort(ucen,iord,1,nunit)
       gapbest = -1d0
       jbest = nunit
       do j = 1, nunit
          if (j < nunit) then
             gap = ucen(iord(j+1)) - ucen(iord(j))
          else
             gap = ucen(iord(1)) + 1d0 - ucen(iord(j))
          end if
          if (gap > gapbest + epsf) then
             gapbest = gap
             jbest = j
          end if
       end do
       ucut = ucen(iord(jbest)) + 0.5d0 * gapbest
       ucut = ucut - floor(ucut)
    end if
    if (present(cutpos)) cutpos = ucut

    ! the units below the cut go on top of the others
    do j = 1, nunit
       if (ucen(j) < ucut - epsf) ish(j) = ish(j) + 1
    end do

    ! the limits of the slab
    umin = huge(1d0)
    umax = -huge(1d0)
    do i = 1, c%ncel
       lsh(i) = lsh(i) + ish(iunit(i))
       u = c%atcel(i)%x(iaxis) + real(lsh(i),8)
       umin = min(umin,u)
       umax = max(umax,u)
    end do

  end subroutine vacuum_layout

  ! Whether idx is a list of cell atoms: not empty, all in range.
  function atomlist_ok(c,idx)
    class(crystal), intent(in) :: c
    integer, intent(in) :: idx(:)
    logical :: atomlist_ok

    atomlist_ok = (size(idx) > 0)
    if (atomlist_ok) atomlist_ok = all(idx >= 1 .and. idx <= c%ncel)

  end function atomlist_ok

  ! The cell atoms in the list idx moved, in general not as their
  ! molecular fragments do: drop the lazily cached molecular data of
  ! those fragments (center of mass, inertia, standard frame, point
  ! group), to be recomputed on demand.
  subroutine atomlist_drop_molcache(c,idx)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: idx(:)

    integer :: k, im

    if (.not.allocated(c%idatcelmol) .or. .not.allocated(c%mol)) return
    do k = 1, size(idx)
       im = c%idatcelmol(1,idx(k))
       if (im >= 1 .and. im <= c%nmol) then
          c%mol(im)%axes_computed = .false.
          c%mol(im)%pg_computed = .false.
       end if
    end do

  end subroutine atomlist_drop_molcache

  ! Refresh the derived data after an in-place edit. For a molecule,
  ! grow the encompassing cell first if the edit pushed an atom out of
  ! it or up against its border: the block environment wraps
  ! coordinates into the cell, so an atom outside it would be binned
  ! against the opposite face and go missing from the neighbor queries
  ! around its true position.
  subroutine refresh_after_inplace_move(c)
    class(crystal), intent(inout) :: c

    integer :: k

    real*8, parameter :: margin = 1d-3 ! fraction of the cell edge

    if (c%ismolecule) then
       do k = 1, c%ncel
          if (any(c%atcel(k)%x < margin) .or. any(c%atcel(k)%x > 1d0-margin)) then
             ! this refreshes the environment for the new cell itself
             call c%recompute_molecular_cell()
             return
          end if
       end do
    end if
    call c%update_env_after_move()

  end subroutine refresh_after_inplace_move

  ! Set the position of cell atom icel to rnew (Cartesian, bohr),
  ! keeping the cached molecular fragment consistent with it (see the
  ! idatcelmol invariant in crystalmod.f90). Used by the in-place
  ! interactive editors; rnew is not wrapped into the main cell.
  subroutine set_cellatom_r(c,icel,rnew)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: icel
    real*8, intent(in) :: rnew(3)

    integer :: im, ia

    c%atcel(icel)%r = rnew
    c%atcel(icel)%x = c%c2x(rnew)
    if (.not.allocated(c%idatcelmol) .or. .not.allocated(c%mol)) return
    im = c%idatcelmol(1,icel)
    ia = c%idatcelmol(2,icel)
    c%mol(im)%at(ia)%x = c%atcel(icel)%x + c%mol(im)%at(ia)%lvec
    c%mol(im)%at(ia)%r = c%x2c(c%mol(im)%at(ia)%x)

  end subroutine set_cellatom_r

  ! Make a seed for an edit-and-rebuild operation. For a molecule, use
  ! an absolute-Cartesian seed (useabr=0) so struct_new re-fits the
  ! encompassing cell to the edited atoms instead of wrapping them into
  ! the old cell. A symmetric edit (copysym) uses the non-equivalent
  ! (nneq-indexed) seed so editing the representative expands to all its
  ! equivalents; bonds (nstar is ncel-indexed) cannot ride on that seed,
  ! so they are recomputed.
  subroutine makeseed_for_edit(c,seed,copysym,copybonding)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(in) :: c
    type(crystalseed), intent(out) :: seed
    logical, intent(in) :: copysym, copybonding

    if (c%ismolecule) then
       call c%makeseed(seed,copysym=.false.,useabr=0,copybonding=copybonding)
    else
       call c%makeseed(seed,copysym=copysym,copybonding=(copybonding .and. .not.copysym))
    end if

  end subroutine makeseed_for_edit

  ! Convert atomic position x in units iunit_l (iunit_ang, iunit_bohr,
  ! or fractional otherwise) to the coordinates used in an edit seed:
  ! internal Cartesian (bohr) for a molecule, crystallographic
  ! otherwise.
  function edit_seed_coords(c,x,iunit_l) result(xx)
    use global, only: iunit_ang, iunit_bohr
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x(3)
    integer, intent(in) :: iunit_l
    real*8 :: xx(3)

    if (iunit_l == iunit_ang) then
       xx = x / bohrtoa
       if (.not.c%ismolecule) xx = c%c2x(xx)
    elseif (iunit_l == iunit_bohr) then
       if (c%ismolecule) then
          xx = x
       else
          xx = c%c2x(x)
       end if
    else
       if (c%ismolecule) then
          xx = c%x2c(x)
       else
          xx = x
       end if
    end if

  end function edit_seed_coords

  !> Given an integer vector v whose components have no common
  !> divisor, find an integer matrix m with determinant +1 whose last
  !> column is v.
  subroutine complete_basis(v,m)
    use tools_math, only: idet3
    integer, intent(in) :: v(3)
    integer, intent(out) :: m(3,3)

    integer :: g, p, q, r, s, g2, w(3), wp(3)

    ! w = (v1,v2,0)/g and wp complete a basis of the xy plane
    call egcd(v(1),v(2),g,p,q)
    if (g == 0) then
       w = (/1,0,0/)
       wp = (/0,1,0/)
    else
       w = (/v(1)/g,v(2)/g,0/)
       wp = (/-q,p,0/)
    end if

    ! v = g*w + v3*z; t = -s*w + r*z completes the basis of the plane (w,z)
    call egcd(g,v(3),g2,r,s)
    m(:,1) = wp
    m(:,2) = -s * w
    m(3,2) = m(3,2) + r
    m(:,3) = v

    if (idet3(m) < 0) m(:,1) = -m(:,1)

  end subroutine complete_basis

  !> Extended Euclidean algorithm: the greatest common divisor g >= 0
  !> of a and b, and two integers p and q such that p*a + q*b = g.
  subroutine egcd(a,b,g,p,q)
    integer, intent(in) :: a, b
    integer, intent(out) :: g, p, q

    integer :: r0, r1, s0, s1, t0, t1, iq, iaux

    r0 = a
    r1 = b
    s0 = 1
    s1 = 0
    t0 = 0
    t1 = 1
    do while (r1 /= 0)
       iq = r0 / r1
       iaux = r0 - iq * r1
       r0 = r1
       r1 = iaux
       iaux = s0 - iq * s1
       s0 = s1
       s1 = iaux
       iaux = t0 - iq * t1
       t0 = t1
       t1 = iaux
    end do
    g = r0
    p = s0
    q = t0
    if (g < 0) then
       g = -g
       p = -p
       q = -q
    end if

  end subroutine egcd

  !> Gauss (Lagrange) reduction of the two-dimensional lattice spanned
  !> by the integer vectors a and b, with the scalar product given by
  !> the metric g: on output, a and b are the two shortest vectors that
  !> span the same lattice, a not longer than b and at an obtuse angle.
  subroutine gauss_reduce(a,b,g)
    integer, intent(inout) :: a(3), b(3)
    real*8, intent(in) :: g(3,3)

    integer :: iaux(3), it
    real*8 :: qa, qb, ab, ratio

    real*8, parameter :: eps = 1d-10
    integer, parameter :: maxit = 10000

    do it = 1, maxit
       qa = qform(a,a,g)
       qb = qform(b,b,g)
       if (qb < qa * (1d0 - eps)) then
          iaux = a
          a = b
          b = iaux
          qa = qb
       end if
       ab = qform(a,b,g)
       ratio = ab / qa
       if (abs(ratio) < 0.5d0 + eps) exit
       b = b - nint(ratio) * a
    end do
    if (qform(a,b,g) > eps * qform(a,a,g)) b = -b

  end subroutine gauss_reduce

  !> Find the shortest lattice vectors that can replace the column k
  !> (1 or 2) of the integer matrix u, whose columns are a basis of a
  !> lattice with metric g, and still give a basis. With k = 1, they
  !> are the shortest vectors that are part of a basis together with
  !> u(:,3); with k = 2, the shortest vectors that complete the basis
  !> started by u(:,1) and u(:,3). On output, xv(:,1:nv) are the
  !> vectors found, all of the same length, from the closest to
  !> perpendicular to u(:,3) to the farthest (at most size(xv,2), and
  !> one of each pair of opposite vectors). Only the vectors with
  !> squared length up to r2 are examined. nv is zero if there is none
  !> or there are too many to examine.
  subroutine shortest_basis_vectors(u,k,g,r2,nv,xv)
    use tools_math, only: matinv, gcd
    integer, intent(in) :: u(3,3)
    integer, intent(in) :: k
    real*8, intent(in) :: g(3,3), r2
    integer, intent(out) :: nv
    integer, intent(inout) :: xv(:,:)

    integer :: nx(3), x(3), ic(3), i1, i2, i3, i2lo, i3lo, j, l
    real*8 :: ginv(3,3), rn(3), ga(3), aga, q, qbest, cs, csv(size(xv,2))

    real*8, parameter :: eps = 1d-10
    real*8, parameter :: nmax = 2d5 ! most vectors to examine

    nv = 0

    ! the components of a vector of length r are at most r times the
    ! lengths of the reciprocal lattice vectors
    ginv = g
    call matinv(ginv,3)
    do i1 = 1, 3
       rn(i1) = sqrt(r2 * ginv(i1,i1)) * (1d0 + eps)
    end do
    if (product(2d0 * rn + 1d0) > nmax) return
    nx = floor(rn)

    ! one of each pair of opposite vectors: the first component that
    ! is not zero is positive
    ga = matmul(g,real(u(:,3),8))
    aga = dot_product(real(u(:,3),8),ga)
    qbest = r2
    do i1 = nx(1), 0, -1
       i2lo = -nx(2)
       if (i1 == 0) i2lo = 0
       do i2 = nx(2), i2lo, -1
          i3lo = -nx(3)
          if (i1 == 0 .and. i2 == 0) i3lo = 1
          do i3 = nx(3), i3lo, -1
             x = (/i1,i2,i3/)
             q = qform(x,x,g)
             if (q > qbest * (1d0 + eps)) cycle

             ! x and u(:,3) are part of a basis if the components of
             ! their vector product have no common divisor; with
             ! u(:,1), the three are a basis if the determinant is 1
             ic(1) = x(2) * u(3,3) - x(3) * u(2,3)
             ic(2) = x(3) * u(1,3) - x(1) * u(3,3)
             ic(3) = x(1) * u(2,3) - x(2) * u(1,3)
             if (k == 1) then
                if (gcd(abs(ic),3) /= 1) cycle
             else
                if (abs(dot_product(int(u(:,1),8),int(ic,8))) /= 1) cycle
             end if

             ! shorter: forget the ones found so far
             if (q < qbest * (1d0 - eps)) then
                nv = 0
                qbest = q
             end if

             ! its place in the list, by the angle with u(:,3)
             cs = abs(dot_product(real(x,8),ga))
             j = nv + 1
             do l = 1, nv
                if (cs < csv(l) - eps * aga) then
                   j = l
                   exit
                end if
             end do
             if (j > size(xv,2)) cycle
             nv = min(nv + 1,size(xv,2))
             do l = nv, j+1, -1
                xv(:,l) = xv(:,l-1)
                csv(l) = csv(l-1)
             end do
             xv(:,j) = x
             csv(j) = cs
          end do
       end do
    end do

  end subroutine shortest_basis_vectors

  !> Scalar product of the integer vectors a and b with the metric g.
  pure function qform(a,b,g)
    integer, intent(in) :: a(3), b(3)
    real*8, intent(in) :: g(3,3)
    real*8 :: qform

    qform = dot_product(real(a,8),matmul(g,real(b,8)))

  end function qform

  !> Append the entry (id, lattice vector lv, bond order ord) to the
  !> neighbor star ns, growing its arrays as needed.
  subroutine star_append(ns,id,lv,ord)
    use types, only: neighstar, realloc
    type(neighstar), intent(inout) :: ns
    integer, intent(in) :: id, lv(3), ord

    integer :: n

    n = ns%ncon + 1
    if (.not.allocated(ns%idcon)) then
       allocate(ns%idcon(n),ns%lcon(3,n),ns%ordcon(n),ns%aromdir(3,n))
    elseif (n > size(ns%idcon,1)) then
       call realloc(ns%idcon,n)
       call realloc(ns%lcon,3,n)
       call realloc(ns%ordcon,n)
       call realloc(ns%aromdir,3,n)
    end if
    ns%idcon(n) = id
    ns%lcon(:,n) = lv
    ns%ordcon(n) = ord
    ns%aromdir(:,n) = 0d0
    ns%ncon = n

  end subroutine star_append

  !> The atoms nkeep+1..ncel have just been added to the crystal
  !> carrying no bonds to the first nkeep: find their distance-based
  !> bonds and add them to c%nstar, then refresh the derived molecular
  !> data. Bonds between two added atoms are included only with newnew
  !> (they are already present when the fragment brought its own
  !> internal bonding), and bonds from an added atom to one of the first
  !> nkeep only with newold (a caller that knows how the fragment
  !> attaches makes those itself). Bonds among the first nkeep atoms are
  !> left untouched.
  subroutine bonds_attach_new(c,nkeep,newnew,newold)
    use global, only: bondfactor
    use types, only: neighstar
    use param, only: atmcov
    class(crystal), intent(inout) :: c
    integer, intent(in) :: nkeep
    logical, intent(in) :: newnew
    logical, intent(in) :: newold

    type(neighstar), allocatable :: ns(:)
    integer :: i, k, nb

    ! distance-based connectivity of the current structure, with the
    ! same radii and factor a full rebuild would use
    call c%find_asterisms(ns,atmcov,bondfactor)

    ! append the computed bonds involving the added atoms. A bond
    ! between two added atoms appears in both stars, so each direction
    ! is appended when the loop visits its origin atom; for a bond to
    ! one of the first nkeep atoms, the reciprocal entry is appended
    ! here as well.
    do i = nkeep+1, c%ncel
       do k = 1, ns(i)%ncon
          nb = ns(i)%idcon(k)
          if (nb <= nkeep) then
             if (.not.newold) cycle
             call star_append(c%nstar(i),nb,ns(i)%lcon(:,k),ns(i)%ordcon(k))
             call star_append(c%nstar(nb),i,-ns(i)%lcon(:,k),ns(i)%ordcon(k))
          elseif (newnew) then
             call star_append(c%nstar(i),nb,ns(i)%lcon(:,k),ns(i)%ordcon(k))
          end if
       end do
    end do

    call c%refresh_molecular_data()

  end subroutine bonds_attach_new

end submodule edit
