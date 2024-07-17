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

  !> Make a crystal seed (seed) from a crystal structure
  module subroutine makeseed(c,seed,copysym)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(in) :: c
    type(crystalseed), intent(out) :: seed
    logical, intent(in) :: copysym

    integer :: i

    ! general
    seed%isused = .true.
    seed%name = ""
    seed%file = c%file

    ! atoms
    if (copysym .and. c%spgavail) then
       seed%nat = c%nneq
       allocate(seed%x(3,c%nneq),seed%is(c%nneq))
       do i = 1, c%nneq
          seed%x(:,i) = c%at(i)%x
          seed%is(i) = c%at(i)%is
       end do
    else
       seed%nat = c%ncel
       allocate(seed%x(3,c%ncel),seed%is(c%ncel))
       do i = 1, c%ncel
          seed%x(:,i) = c%atcel(i)%x
          seed%is(i) = c%atcel(i)%is
       end do
    end if

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
    seed%m_x2c = c%m_x2c

    ! symmetry
    seed%findsym = -1
    seed%checkrepeats = .false.
    if (copysym .and. c%spgavail) then
       seed%havesym = 1
       seed%neqv = c%neqv
       seed%ncv = c%ncv
       allocate(seed%rotm(3,4,c%neqv),seed%cen(3,c%ncv))
       seed%rotm = c%rotm(:,:,1:c%neqv)
       seed%cen = c%cen(:,1:c%ncv)
    else
       seed%havesym = 0
       seed%neqv = 0
       seed%ncv = 0
    end if

    ! molecular fields
    seed%ismolecule = c%ismolecule
    seed%cubic = (abs(c%aa(1)-c%aa(2)) < 1d-5).and.(abs(c%aa(1)-c%aa(3)) < 1d-5).and.&
       all(abs(c%bb-90d0) < 1d-3)
    seed%border = 0d0
    seed%havex0 = .true.
    seed%molx0 = c%molx0

  end subroutine makeseed

  !> Given a crystal structure (c) and three lattice vectors in cryst.
  !> coords (x0(:,1), x0(:,2), x0(:,3)), build the same crystal
  !> structure using the unit cell given by those vectors. If nnew,
  !> xnew, and isnew are given, replace the nnew atoms with the
  !> new positions (xnew) and species (isnew). If noenv is present
  !> and true, do not load the atomic grids or the environments in
  !> the new cell.
  module subroutine newcell(c,x00,t0,nnew,xnew,isnew,noenv,ti)
    use crystalseedmod, only: crystalseed
    use tools_math, only: det3, matinv, mnorm2
    use tools_io, only: ferror, faterr, warning, string, uout
    use types, only: realloc
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x00(3,3)
    real*8, intent(in), optional :: t0(3)
    integer, intent(in), optional :: nnew
    real*8, intent(in), optional :: xnew(:,:)
    integer, intent(in), optional :: isnew(:)
    logical, intent(in), optional :: noenv
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: ncseed
    logical :: ok, found, atgiven, issimple
    real*8 :: x0(3,3), x0inv(3,3), fvol, dmax0
    real*8 :: x(3), dx(3), dd, t(3), xshift(3)
    integer :: i, j, k, m
    integer :: nn
    integer :: nlat, nlat2, nlatnew, nvec(3), ntot
    integer, allocatable :: lvec(:,:)
    real*8, allocatable :: xlat(:,:)

    real*8, parameter :: epszero = 1d-8
    real*8, parameter :: eps = 1d-2
    real*8, parameter :: eps2 = eps * eps

    if (c%ismolecule) &
       call ferror('newcell','NEWCELL incompatible with molecules',faterr)

    ! initialize
    atgiven = (present(nnew).and.present(xnew).and.present(isnew))
    x0 = x00
    dd = det3(x0)

    ! check the new lattice vectors are sane
    if (abs(dd) < eps) then
       call ferror('newcell','invalid input vectors',faterr)
    elseif (dd < 0d0) then
       ! flip the cell
       x0 = -x0
       dd = det3(x0)
       call ferror('newcell','det < 0; vectors flipped',warning)
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
       if (size(xnew,1) /= 3 .or. size(xnew,2) /= nnew .or. size(isnew) /= nnew) &
          call ferror('newcell','NEWCELL: error in size of new atomic positions array',faterr)

       ! build the atomic info in the new seed from the info in input
       ncseed%nat = nnew
       allocate(ncseed%x(3,nnew),ncseed%is(nnew))
       do i = 1, nnew
          ncseed%x(:,i) = xnew(:,i)
          ncseed%is(i) = isnew(i)
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
          allocate(ncseed%x(3,c%ncel * ntot),ncseed%is(c%ncel * ntot))

          nn = 0
          do i = 1, nvec(1)
             do j = 1, nvec(2)
                do k = 1, nvec(3)
                   xshift = real((/i,j,k/) - 1,8) / real(nvec,8)
                   do m = 1, c%ncel
                      nn = nn + 1
                      ncseed%is(nn) = c%atcel(m)%is
                      ncseed%x(:,nn) = (c%atcel(m)%x-t) / real(nvec,8) + xshift
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
                call ferror("newcell","Cell vector number " // string(i) // &
                   " is not a pure translation",faterr)
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
          call c%list_near_lattice_points((/0d0,0d0,0d0/),icrd_crys,.true.,nlat2,&
             lvec=lvec,up2d=dmax0*2.0d0,nozero=.true.)

          ! make lattice vector list
          nlatnew = max(nint(dd),1)
          allocate(xlat(3,nlatnew))
          nlat = 1
          xlat(:,1) = 0d0
          do i = 1, nlat2
             if (nlat == nlatnew) exit
             x = real(lvec(:,i),8)
             x = matmul(x,x0inv)
             x = x - floor(x)

             found = .false.
             do j = 1, nlat
                dx = xlat(:,j) - x
                dx = dx - nint(dx)
                if (all(abs(dx) < eps)) then
                   found = .true.
                   exit
                end if
             end do
             if (.not.found) then
                nlat = nlat + 1
                xlat(:,nlat) = x
             end if
          end do

          ! build the new atom list
          ncseed%nat = 0
          fvol = abs(det3(x0))
          nn = nint(c%ncel * fvol)
          if (abs(nn - (c%ncel*fvol)) > eps) &
             call ferror('newcell','inconsistent number of atoms in newcell',faterr)
          allocate(ncseed%x(3,nn),ncseed%is(nn))
          do i = 1, nlat
             do j = 1, c%ncel
                ! candidate atom
                x = matmul(c%atcel(j)%x-t,x0inv) + xlat(:,i)
                x = x - floor(x)

                ! check if we have it already
                ok = .true.
                do m = 1, ncseed%nat
                   dx = x - ncseed%x(:,m)
                   dx = matmul(ncseed%m_x2c,dx - nint(dx))
                   if (dot_product(dx,dx) < eps2) then
                      ok = .false.
                      exit
                   end if
                end do
                if (ok) then
                   ! add it to the list
                   ncseed%nat = ncseed%nat + 1
                   if (ncseed%nat > size(ncseed%x,2)) then
                      call realloc(ncseed%x,3,2*ncseed%nat)
                      call realloc(ncseed%is,2*ncseed%nat)
                   end if
                   ncseed%x(:,ncseed%nat) = x
                   ncseed%is(ncseed%nat) = c%atcel(j)%is
                end if
             end do
          end do
          call realloc(ncseed%x,3,ncseed%nat)
          call realloc(ncseed%is,ncseed%nat)
          deallocate(xlat)
       end if
    end if

    ! rest of the seed information
    ncseed%isused = .true.
    ncseed%file = c%file
    ncseed%havesym = 0
    ncseed%findsym = -1
    ncseed%ismolecule = .false.

    ! initialize the structure
    call c%struct_new(ncseed,.true.,noenv,ti=ti)

  end subroutine newcell

  !> Transform to the standard cell. If toprim, convert to the
  !> primitive standard cell. If doforce = .true., force the
  !> transformation to the primitive even if it does not lead to a
  !> smaller cell. Refine = refine the symmetry positions. noenv = do
  !> not initialize the environment. Return the transformation matrix,
  !> or a matrix of zeros if no change was done.
  module function cell_standard(c,toprim,doforce,refine,noenv,ti) result(x0)
    use iso_c_binding, only: c_double
    use spglib, only: spg_standardize_cell
    use global, only: symprec
    use tools_math, only: det3, matinv
    use tools_io, only: ferror, faterr
    use types, only: realloc
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in) :: toprim
    logical, intent(in) :: doforce
    logical, intent(in) :: refine
    logical, intent(in), optional :: noenv
    real*8 :: x0(3,3)
    type(thread_info), intent(in), optional :: ti

    integer :: ntyp, nat
    integer :: i, nnew, iprim, inorefine
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: types_(:)
    real*8 :: rmat(3,3)

    x0 = 0d0

    ! ignore molecules
    if (c%ismolecule) return

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
    if (nnew == 0) &
       call ferror("cell_standard","could not find primitive cell",faterr)
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
       call c%newcell(rmat,nnew=nnew,xnew=x,isnew=types_,noenv=noenv,ti=ti)
    else
       ! if a primitive is wanted but det is not less than 1, do not make the change
       if (all(abs(rmat - eye) < symprec)) return
       if (toprim .and. .not.(det3(rmat) < 1d0-symprec) .and..not.doforce) return
       call c%newcell(rmat,noenv=noenv,ti=ti)
    end if
    x0 = rmat

  end function cell_standard

  !> Transform to the Niggli cell. Return the transformation matrix,
  !> or a matrix of zeros if no change was done.
  module function cell_niggli(c,noenv,ti) result(x0)
    use spglib, only: spg_niggli_reduce
    use global, only: symprec
    use tools_io, only: ferror, faterr
    use tools_math, only: det3
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: noenv
    type(thread_info), intent(in), optional :: ti
    real*8 :: x0(3,3)

    real*8 :: rmat(3,3)
    integer :: id, i

    x0 = 0d0

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib niggli reduction
    rmat = transpose(c%m_x2c)
    id = spg_niggli_reduce(rmat,symprec)
    if (id == 0) &
       call ferror("cell_niggli","could not find Niggli reduction",faterr)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det3(rmat) < 0d0) rmat = -rmat

    ! transform
    if (any(abs(rmat - eye) > symprec)) then
       call c%newcell(rmat,noenv=noenv,ti=ti)
       x0 = rmat
    end if

  end function cell_niggli

  !> Transform to the Delaunay cell. Return the transformation matrix,
  !> or a matrix of zeros if no change was done.
  module function cell_delaunay(c,noenv,ti) result(x0)
    use spglib, only: spg_delaunay_reduce
    use global, only: symprec
    use tools_io, only: ferror, faterr
    use tools_math, only: det3
    use param, only: eye
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: noenv
    type(thread_info), intent(in), optional :: ti
    real*8 :: x0(3,3)

    real*8 :: rmat(3,3)
    integer :: id, i

    x0 = 0d0

    ! ignore molecules
    if (c%ismolecule) return

    ! use spglib delaunay reduction
    rmat = transpose(c%m_x2c)
    id = spg_delaunay_reduce(rmat,symprec)
    if (id == 0) &
       call ferror("cell_delaunay","could not find Delaunay reduction",faterr)
    rmat = transpose(rmat)
    do i = 1, 3
       rmat(:,i) = c%c2x(rmat(:,i))
    end do

    ! flip the cell?
    if (det3(rmat) < 0d0) rmat = -rmat

    ! transform
    if (any(abs(rmat - eye) > symprec)) then
       call c%newcell(rmat,noenv=noenv,ti=ti)
       x0 = rmat
    end if

  end function cell_delaunay

  !> Create a new structure by reordering the atoms in the current
  !> structure. iperm is the permutation vector (atom i in the new
  !> structure is iperm(i) in the old structure).
  module subroutine reorder_atoms(c,iperm,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iperm(:)
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8, allocatable :: x(:,:)
    integer, allocatable :: is(:)
    integer :: i

    ! make the new seed
    call c%makeseed(seed,.false.)

    ! reorder the atoms
    allocate(x(3,seed%nat),is(seed%nat))
    do i = 1, seed%nat
       x(:,i) = seed%x(:,iperm(i))
       is(i) = seed%is(iperm(i))
    end do
    seed%x = x
    seed%is = is
    deallocate(x,is)

    ! reload the crystal
    call c%struct_new(seed,.true.,ti=ti)

  end subroutine reorder_atoms

  !> Re-assign atomic types to have an asymmetric unit with whole molecules
  module subroutine wholemols(c,ti)
    use crystalseedmod, only: crystalseed
    use tools, only: qcksort
    use tools_io, only: ferror, faterr
    use types, only: realloc
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
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
             if (.not.found) &
                call ferror("wholemols","symmetry operations do not form a group",faterr)
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

                if (id == 0) &
                   call ferror('wholemols','error identifying rotated atom',faterr)

                if (c%idatcelmol(k) == c%idatcelmol(id)) &
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
    allocate(ncseed%x(3,c%ncel),ncseed%is(c%ncel))
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
             if (id == 0) &
                call ferror('wholemols','error identifying rotated atom',faterr)
             isuse(id) = .true.
          end do
       end do

       ncseed%nat = ncseed%nat + 1
       ncseed%x(:,ncseed%nat) = c%atcel(i)%x
       ncseed%is(ncseed%nat) = c%atcel(i)%is
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
    ncseed%havesym = 1
    ncseed%findsym = 0
    ncseed%checkrepeats = .false.

    ! build the new crystal
    call c%struct_new(ncseed,.true.,ti=ti)

  end subroutine wholemols

  !> Delete the atoms with IDs in the array iat(1:nat) from the
  !> structure.
  module subroutine delete_atoms(c,nat,iat,ti)
    use crystalseedmod, only: crystalseed
    use types, only: realloc
    class(crystal), intent(inout) :: c
    integer, intent(in) :: nat
    integer, intent(in) :: iat(nat)
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    logical, allocatable :: useatoms(:)
    integer, allocatable :: usespcs(:)
    integer :: i, nnspc

    ! make seed from this crystal
    call c%makeseed(seed,copysym=.false.)

    ! flag the atoms and species to use
    allocate(useatoms(c%ncel),usespcs(c%nspc))
    useatoms = .true.
    do i = 1, nat
       useatoms(iat(i)) = .false.
    end do
    usespcs = 0
    nnspc = 0
    do i = 1, c%ncel
       if (useatoms(i)) then
          if (usespcs(c%atcel(i)%is) == 0) then
             nnspc = nnspc + 1
             usespcs(c%atcel(i)%is) = nnspc
          end if
       end if
    end do

    ! re-do the atom and species info
    seed%nat = 0
    do i = 1, c%ncel
       if (useatoms(i)) then
          seed%nat = seed%nat + 1
          seed%x(:,seed%nat) = c%atcel(i)%x
          seed%is(seed%nat) = usespcs(c%atcel(i)%is)
       end if
    end do
    call realloc(seed%x,3,seed%nat)
    call realloc(seed%is,seed%nat)
    seed%nspc = nnspc
    do i = 1, c%nspc
       if (usespcs(i) > 0) seed%spc(usespcs(i)) = c%spc(i)
    end do
    call realloc(seed%spc,seed%nspc)

    ! build the new crystal
    call c%struct_new(seed,crashfail=.true.,ti=ti)

  end subroutine delete_atoms

  !> Move atom with complete list ID idx to position x in units of
  !> iunit_l (see global). If dorelative, the movement is relative to
  !> its current position.
  module subroutine move_atom(c,idx,x,iunit_l,dorelative,ti)
    use crystalseedmod, only: crystalseed
    use global, only: iunit_ang, iunit_bohr
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    integer, intent(in) :: idx
    real*8, intent(in) :: x(3)
    integer, intent(in) :: iunit_l
    logical, intent(in) :: dorelative
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8 :: xx(3)

    ! make seed from this crystal
    call c%makeseed(seed,copysym=.false.)

    ! interpret units
    if (iunit_l == iunit_ang) then
       xx = x / bohrtoa
       xx = c%c2x(xx)
    elseif (iunit_l == iunit_bohr) then
       xx = c%c2x(x)
    else
       xx = x
    end if

    if (dorelative) then
       seed%x(:,idx) = seed%x(:,idx) + xx
    else
       seed%x(:,idx) = xx
    end if

    ! build the new crystal
    call c%struct_new(seed,crashfail=.true.,ti=ti)

  end subroutine move_atom

end submodule edit
