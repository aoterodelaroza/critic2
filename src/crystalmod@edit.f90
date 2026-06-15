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
    use tools_io, only: ferror, faterr
    class(crystal), intent(in) :: c
    type(crystalseed), intent(out) :: seed
    logical, intent(in) :: copysym
    integer, intent(in), optional :: useabr
    logical, intent(in), optional :: copybonding

    integer :: i
    integer :: useabr_
    logical :: copybonding_

    ! initialize
    useabr_ = 2
    if (present(useabr)) useabr_ = useabr
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! copybonding and copysym are not compatible
    if (copybonding_ .and. copysym) &
       call ferror('makeseed','copybonding and copysym are incompatible',faterr)

    ! general
    seed%isused = .true.
    seed%name = ""
    seed%file = c%file
    seed%isformat = c%isformat

    ! atoms
    if (.not.c%ismolecule .and. copysym .and. c%spgavail) then
       ! crystals with symmetry
       seed%nat = c%nneq
       allocate(seed%x(3,c%nneq),seed%is(c%nneq),seed%atname(c%nneq))
       do i = 1, c%nneq
          seed%x(:,i) = c%at(i)%x
          seed%is(i) = c%at(i)%is
          seed%atname(i) = c%at(i)%name
       end do
    else
       ! crystals without symmetry
       seed%nat = c%ncel
       allocate(seed%x(3,c%ncel),seed%is(c%ncel),seed%atname(c%ncel))
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
       end do
    end if

    ! bonding (nstar is indexed by ncel)
    if (copybonding_ .and. allocated(c%nstar)) then
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
    use param, only: atmass, pi, img, tpi
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

    ! general
    seed%isused = .true.
    seed%name = trim(c%file) // " (nudge)"
    seed%file = c%file
    seed%isformat = c%isformat

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
  module subroutine newcell(c,x00,t0,nnew,xnew,isnew,noenv,errmsg,ti)
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

    ! check the new lattice vectors are sane
    if (abs(dd) < eps) then
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
          if (abs(nn - (c%ncel*fvol)) > eps) then
             errmsg_ = 'Inconsistent number of atoms in newcell'
             goto 999
          end if
          allocate(ncseed%x(3,nn),ncseed%is(nn),ncseed%atname(nn))
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
                      call realloc(ncseed%atname,2*ncseed%nat)
                   end if
                   ncseed%x(:,ncseed%nat) = x
                   ncseed%is(ncseed%nat) = c%atcel(j)%is
                   ncseed%atname(ncseed%nat) = c%at(c%atcel(j)%idx)%name
                end if
             end do
          end do
          call realloc(ncseed%x,3,ncseed%nat)
          call realloc(ncseed%is,ncseed%nat)
          call realloc(ncseed%atname,ncseed%nat)
          deallocate(xlat)
       end if
    end if

    ! rest of the seed information
    ncseed%isused = .true.
    ncseed%file = c%file
    ncseed%isformat = c%isformat
    ncseed%havesym = 0
    ncseed%findsym = -1
    ncseed%ismolecule = .false.

    ! initialize the structure
    call c%struct_new(ncseed,.false.,noenv,ti=ti)
    if (.not.c%isinit) errmsg = "Could not create the new crystal structure"

    return

    ! error handling: report via errmsg
999 continue
    errmsg = errmsg_

  end subroutine newcell

  !> Transform to the standard cell. If toprim, convert to the
  !> primitive standard cell. If doforce = .true., force the
  !> transformation to the primitive even if it does not lead to a
  !> smaller cell. Refine = refine the symmetry positions. noenv = do
  !> not initialize the environment. Return the transformation matrix,
  !> or a matrix of zeros if no change was done.
  module function cell_standard(c,toprim,doforce,refine,noenv,errmsg,ti) result(x0)
    use iso_c_binding, only: c_double
    use spglib, only: spg_standardize_cell
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
    real*8 :: x0(3,3)

    integer :: ntyp, nat
    integer :: i, nnew, iprim, inorefine
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: types_(:)
    real*8 :: rmat(3,3)

    x0 = 0d0
    errmsg = ""

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

  !> Search for "nice" supercells of the current cell, of increasing
  !> size from 1 up to inice times the input cell. For each size n,
  !> return in rmax(n) the radius of the largest sphere that fits in
  !> the nicest supercell of that size, and in mmax(:,:,n) the
  !> corresponding NEWCELL transformation matrix (lattice vectors of
  !> the supercell in the old setting, crystallographic coordinates).
  !> rmax(n) is zero if no supercell of that size was found.
  module subroutine cell_nice_list(c,inice,rmax,mmax)
    use tools_math, only: cross, det3
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    integer, intent(in) :: inice
    real*8, allocatable, intent(out) :: rmax(:)
    real*8, allocatable, intent(out) :: mmax(:,:,:)

    integer :: i, j, k, n, nat
    real*8 :: dmax0, mm(3,3), dd, x2c(3,3), r
    integer, allocatable :: lvec(:,:)

    allocate(rmax(inice),mmax(3,3,inice))
    rmax = 0d0
    mmax = 0d0
    if (c%ismolecule) return

    dmax0 = 1.2d0 * real(inice,8)**(1d0/3d0) * max(norm2(c%m_xr2c(:,1)),norm2(c%m_xr2c(:,2)),&
       norm2(c%m_xr2c(:,3))) + 1d-1
    call c%list_near_lattice_points((/0d0,0d0,0d0/),icrd_crys,.true.,nat,&
       lvec=lvec,up2d=dmax0,nozero=.true.)

    do i = 1, nat
       mm(:,1) = lvec(:,i)
       do j = i+1, nat
          mm(:,2) = lvec(:,j)
          do k = j+1, nat
             mm(:,3) = lvec(:,k)

             dd = det3(mm)
             if (dd < 0d0) mm = -mm
             dd = abs(dd)
             if (dd < 1d-2 .or. dd > (inice+0.5d0)) cycle
             n = nint(dd)
             if (n > inice) cycle

             x2c = matmul(c%m_x2c,mm)

             r = 0.5d0 * n * c%omega / max(norm2(cross(x2c(:,1),x2c(:,2))),&
                norm2(cross(x2c(:,1),x2c(:,3))),norm2(cross(x2c(:,2),x2c(:,3))))
             if (r > rmax(n)) then
                rmax(n) = r
                mmax(:,:,n) = mm
             end if
          end do
       end do
    end do

  end subroutine cell_nice_list

  !> Create a new structure by reordering the atoms in the current
  !> structure. iperm is the permutation vector (atom i in the new
  !> structure is iperm(i) in the old structure).
  module subroutine reorder_atoms(c,iperm,isnneq,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iperm(:)
    logical, intent(in) :: isnneq
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8, allocatable :: x(:,:)
    integer, allocatable :: is(:)
    integer :: i
    logical :: copysym

    ! make the new seed
    copysym = isnneq .and. .not.c%ismolecule .and. c%spgavail
    call c%makeseed(seed,copysym)

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

    ! reload the crystal
    call c%struct_new(seed,.true.,ti=ti)

  end subroutine reorder_atoms

  !> Create a new structure by reordering the molecular fragments in the
  !> current structure. iperm is the permutation vector (molecule i in
  !> the new structure is iperm(i) in the old structure). The atoms are
  !> relabeled so that the cell atoms are grouped by molecule and the
  !> molecules appear in the new order; molecular numbering follows the
  !> order of the atoms, so the reloaded structure has the requested
  !> molecular order.
  module subroutine reorder_molecules(c,iperm,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iperm(:)
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8, allocatable :: x(:,:)
    integer, allocatable :: atperm(:)
    integer :: i, k, im, nat

    ! the molecular decomposition must be available and the permutation valid
    if (c%nmol == 0 .or. .not.allocated(c%idatcelmol)) return
    if (size(iperm,1) /= c%nmol) return

    ! make the new seed (no symmetry, so seed atoms follow the cell-atom order)
    call c%makeseed(seed,.false.)
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

    ! reload the crystal
    call c%struct_new(seed,.true.,ti=ti)

  end subroutine reorder_molecules

  !> Create a new structure by reordering the species in the current
  !> structure. iperm is the permutation vector (atom i in the new
  !> structure is iperm(i) in the old structure).
  module subroutine reorder_species(c,iperm,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iperm(:)
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    integer :: i
    integer, allocatable :: iiperm(:)

    ! make the new seed
    call c%makeseed(seed,.false.)

    ! return if the permutation is not correct
    if (size(iperm,1) /= seed%nspc) return
    allocate(iiperm(size(iperm,1)))
    do i = 1, size(iperm,1)
       iiperm(iperm(i)) = i
    end do

    ! reorder the species
    seed%spc = seed%spc(iperm(:))
    do i = 1, seed%nat
       seed%is(i) = iiperm(seed%is(i))
    end do

    ! reload the crystal
    call c%struct_new(seed,.true.,ti=ti)

  end subroutine reorder_species

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
       ncseed%atname(ncseed%nat) = c%at(c%atcel(i)%idx)%name
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
    ncseed%isformat = c%isformat
    ncseed%havesym = 1
    ncseed%findsym = 0
    ncseed%checkrepeats = .false.

    ! build the new crystal
    call c%struct_new(ncseed,.true.,ti=ti)

  end subroutine wholemols

  !> Remove or merge the atoms with IDs in the array iat(1:nat) from
  !> the structure.
  module subroutine edit_atom_list(c,nat,iat,remove,merge,duplicate,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use types, only: realloc
    class(crystal), intent(inout) :: c
    integer, intent(in) :: nat
    integer, intent(in) :: iat(nat)
    logical, intent(in), optional :: remove, merge, duplicate
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    logical, allocatable :: useatoms(:), dupatoms(:)
    integer :: i, ipres, izmax, mergespc, natnew
    real*8 :: mergex(3), x0(3), xd(3), rdum
    logical :: remove_, merge_, duplicate_

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

    ! make seed from this crystal
    call c%makeseed(seed,copysym=.false.)

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

    ! re-do the atom and species info
    seed%nat = 0
    do i = 1, c%ncel
       if (useatoms(i)) then
          seed%nat = seed%nat + 1
          seed%x(:,seed%nat) = c%atcel(i)%x
          seed%atname(seed%nat) = c%at(c%atcel(i)%idx)%name
          seed%is(seed%nat) = c%atcel(i)%is
          if (dupatoms(i)) then
             seed%nat = seed%nat + 1
             seed%x(:,seed%nat) = c%atcel(i)%x
             seed%atname(seed%nat) = c%at(c%atcel(i)%idx)%name
             seed%is(seed%nat) = c%atcel(i)%is
          end if
       end if
    end do

    ! if merging, add the merged atom at the end
    if (merge_) then
       seed%nat = seed%nat + 1
       seed%x(:,seed%nat) = mergex
       seed%is(seed%nat) = mergespc
       seed%atname(seed%nat) = c%spc(mergespc)%name
    end if

    ! build the new crystal
    call c%struct_new(seed,crashfail=.false.,ti=ti)
    if (.not.c%isinit) errmsg = "Could not rebuild the structure after editing the atoms"

  end subroutine edit_atom_list

  !> Change the atoms with IDs in the array iat(1:nat) from their
  !> current species to is, and reset the atom name to the
  !> corresponding species name.
  module subroutine change_atom_species(c,nat,iat,is,copybonding,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    integer, intent(in) :: nat
    integer, intent(in) :: iat(nat)
    integer, intent(in) :: is
    logical, intent(in), optional :: copybonding
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    integer :: i
    logical :: copybonding_

    ! return if nothing to do
    if (nat == 0) return
    if (is < 1 .or. is > c%nspc) return

    ! make seed from this crystal, preserving bonding if requested
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    call c%makeseed(seed,copysym=.false.,copybonding=copybonding_)

    ! apply the change
    do i = 1, nat
       seed%is(iat(i)) = is
       seed%atname(iat(i)) = c%spc(is)%name
    end do

    ! build the new crystal
    call c%struct_new(seed,crashfail=.true.,ti=ti)

  end subroutine change_atom_species

  !> Move atom with complete list ID idx to position x in units of
  !> iunit_l (see global). If isnneq, move all atoms that are
  !> equivalent by symmetry. If dorelative, the movement is
  !> relative to its current position.
  module subroutine move_atom(c,idx,x,iunit_l,isnneq,dorelative,copybonding,ti)
    use crystalseedmod, only: crystalseed
    use global, only: iunit_ang, iunit_bohr
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    integer, intent(in) :: idx
    real*8, intent(in) :: x(3)
    integer, intent(in) :: iunit_l
    logical, intent(in) :: isnneq
    logical, intent(in) :: dorelative
    logical, intent(in), optional :: copybonding
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8 :: xx(3)
    logical :: copysym, copybonding_

    ! whether to use symmetry
    copysym = isnneq .and. .not.c%ismolecule .and. c%spgavail
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! make seed from this crystal, preserving bonding if requested. For a
    ! molecule, use an absolute-Cartesian seed (useabr=0) so struct_new re-fits
    ! the encompassing cell to the moved atom instead of wrapping it into the
    ! old cell. copybonding is incompatible with copysym (nneq-indexed seed).
    if (c%ismolecule) then
       call c%makeseed(seed,copysym=.false.,useabr=0,copybonding=copybonding_)
    else
       call c%makeseed(seed,copysym=copysym,copybonding=(copybonding_ .and. .not.copysym))
    end if

    ! interpret units. For a molecule, xx is kept as an internal Cartesian
    ! position (bohr); otherwise it is converted to crystallographic.
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
    call c%struct_new(seed,crashfail=.true.,ti=ti)

  end subroutine move_atom

  !> Recalculate the unit cell of a molecule to fit the current atomic
  !> positions, growing it as needed so the molecule (e.g. after a
  !> rigid rotation or translation in place) stays inside the
  !> cell. The absolute atomic positions (Cartesian + molx0) are
  !> preserved, so the rendered scene in the GUI does not move. The
  !> molecular cell (molborder) and the vacuum descriptors are
  !> recalculated. The neighbor star, molecular decomposition, and
  !> Wigner-Seitz data are NOT recomputed (no rebonding); they are
  !> reconciled by the next full rebuild (struct_new). This routine is
  !> only used for molecules, and only by the GUI (move_molecule and
  !> rotate_molecule with no rebonding).
  module subroutine recompute_molecular_cell(c)
    use tools_math, only: m_x2c_from_cellpar, matinv
    use global, only: rborder_def
    class(crystal), intent(inout) :: c

    integer :: i, j, k
    real*8 :: smin(3), smax(3), border, sh(3), saux(3)

    if (.not.c%ismolecule .or. c%ncel == 0) return

    ! bounding box of the absolute atomic positions (r + molx0), padded with
    ! the default border (same construction as struct_new's molecule branch)
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

    ! new orthogonal cell. The cell origin (molx0) sits at the lower padded
    ! corner, so the absolute positions (r + molx0) are left unchanged; sh is
    ! the uniform Cartesian shift applied to every internal Cartesian coordinate.
    sh = c%molx0 - smin
    c%aa = smax - smin
    c%bb = 90d0
    c%m_x2c = m_x2c_from_cellpar(c%aa,c%bb)
    c%m_c2x = c%m_x2c
    call matinv(c%m_c2x,3)
    c%molx0 = smin

    ! molecular cell used for display (same formula as struct_new)
    c%molborder = max(border - max(2d0,0.8d0 * border),0d0) / c%aa

    ! shift the cell atoms (and the matching non-equivalent atoms): the
    ! Cartesian positions move rigidly by sh, the fractional ones are
    ! recomputed in the new cell
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

    ! recompute the vacuum descriptors for the new cell
    call c%calc_vacuum_lengths()

  end subroutine recompute_molecular_cell

  !> Rigidly translate the molecular fragment imol so that its center
  !> of mass is at position x in units of iunit_l (see global). If
  !> dorelative, x is interpreted as a displacement of the center of
  !> mass relative to its current position. All atoms in the fragment
  !> are translated by the same vector; symmetry is not preserved (P1).
  module subroutine move_molecule(c,imol,x,iunit_l,dorelative,copybonding,ti)
    use crystalseedmod, only: crystalseed
    use global, only: iunit_ang, iunit_bohr
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    integer, intent(in) :: imol
    real*8, intent(in) :: x(3)
    integer, intent(in) :: iunit_l
    logical, intent(in) :: dorelative
    logical, intent(in), optional :: copybonding
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8 :: xx(3), dx(3), dxc(3)
    integer :: k
    logical :: copybonding_

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
    ! order), preserving bonding if requested. For a molecule, use an
    ! absolute-Cartesian seed (useabr=0) so struct_new re-fits the encompassing
    ! cell to the moved molecule.
    if (c%ismolecule) then
       call c%makeseed(seed,copysym=.false.,useabr=0,copybonding=copybonding_)
    else
       call c%makeseed(seed,copysym=.false.,copybonding=copybonding_)
    end if
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
    call c%struct_new(seed,crashfail=.true.,ti=ti)

  end subroutine move_molecule

  !> Rigidly rotate the molecular fragment imol about its center of mass
  !> so that the Euler angles (ZYZ, radians) of its standard orientation
  !> become euler. The rotation uses the fragment atoms (whole molecule,
  !> with lattice translations applied), not the cell atoms, so it stays
  !> rigid even when the molecule is split across cell boundaries. Only
  !> applies to discrete fragments. Symmetry is not preserved (P1).
  module subroutine rotate_molecule(c,imol,euler,copybonding,ti)
    use crystalseedmod, only: crystalseed
    use tools_math, only: euler2mat, mat2euler, mat2quat
    class(crystal), intent(inout) :: c
    integer, intent(in) :: imol
    real*8, intent(in) :: euler(3)
    logical, intent(in), optional :: copybonding
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8 :: rnew(3), rmat(3,3), rrot(3,3), xcm(3)
    integer :: j, k, lvec(3)
    logical :: copybonding_

    ! consistency checks
    if (imol < 1 .or. imol > c%nmol .or. .not.allocated(c%idatcelmol)) return
    if (.not.c%mol(imol)%discrete) return
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! make sure the standard frame is available
    if (.not.c%mol(imol)%axes_computed) call c%mol(imol)%compute_std()

    ! rigid rotation that brings the current standard frame to the target
    rmat = euler2mat(euler)
    rrot = matmul(rmat,transpose(c%mol(imol)%m_std))
    xcm = c%mol(imol)%xcm

    ! make seed from this crystal (no symmetry, so seed atoms follow cell-atom
    ! order), preserving bonding if requested. For a molecule, use an
    ! absolute-Cartesian seed (useabr=0) so struct_new re-fits the encompassing
    ! cell to the rotated molecule instead of wrapping atoms into the old cell.
    if (c%ismolecule) then
       call c%makeseed(seed,copysym=.false.,useabr=0,copybonding=copybonding_)
    else
       call c%makeseed(seed,copysym=.false.,copybonding=copybonding_)
    end if
    if (seed%nat /= c%ncel) return

    ! rotate the fragment atoms (whole molecule) and write them to the seed
    do j = 1, c%mol(imol)%nat
       k = c%mol(imol)%at(j)%cidx
       if (c%ismolecule) then
          seed%x(:,k) = c%molx0 + matmul(rrot,c%mol(imol)%at(j)%r - xcm) ! cartesian
       else
          lvec = nint(seed%x(:,k) - c%c2x(c%mol(imol)%at(j)%r))
          rnew = xcm + matmul(rrot,c%mol(imol)%at(j)%r - xcm)
          seed%x(:,k) = c%c2x(rnew) + lvec ! crystallographic
       end if
    end do

    ! build the new crystal
    call c%struct_new(seed,crashfail=.true.,ti=ti)

    ! set the Euler angles explicitly (useful for continuous drag in the GUI)
    ! only if the bonding is preserved - otherwise the molecules may change
    if (copybonding_) then
       call c%mol(imol)%compute_std()
       c%mol(imol)%euler_std = euler
    end if

  end subroutine rotate_molecule

  !> Modify the unit cell by changing the parameter given by iaxis:
  !> 1=a, 2=b, 3=c, -1=alpha, -2=beta, -3=gamma, 0=volume. The
  !> parameter is changed to x in units of iunit_l. If dorelative, the
  !> change is relative to the current value. If dofraction, x is
  !> interpreted as the fractional change in the current value.
  module subroutine move_cell(c,iaxis,x,iunit_l,dorelative,dofraction,ti)
    use crystalseedmod, only: crystalseed
    use global, only: iunit_ang
    use param, only: bohrtoa, third
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iaxis
    real*8, intent(in) :: x
    integer, intent(in) :: iunit_l
    logical, intent(in) :: dorelative, dofraction
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8 :: xx, ref

    ! make seed from this crystal
    if (iaxis /= 0) then
       call c%makeseed(seed,copysym=.false.,useabr=1)
    else
       call c%makeseed(seed,copysym=.false.,useabr=2)
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
    call c%struct_new(seed,crashfail=.true.,ti=ti)

  end subroutine move_cell

  !> Modify the unit cell by changing the cell lengths and axes to the
  !> given values.
  module subroutine move_cell_all(c,aa,bb,copybonding,ti)
    use crystalseedmod, only: crystalseed
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: aa(3), bb(3)
    logical, intent(in), optional :: copybonding
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    logical :: copybonding_

    ! make seed from this crystal, preserving bonding if requested
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    call c%makeseed(seed,copysym=.false.,useabr=1,copybonding=copybonding_)

    ! set the new axes
    seed%aa = aa
    seed%bb = bb

    ! build the new crystal
    call c%struct_new(seed,crashfail=.true.,ti=ti)

  end subroutine move_cell_all

  !> Add atom with species is and position x in units of iunit_l (see
  !> global). If is <= 0, add a new species with Z = abs(is) to the
  !> system. If isnneq, replicate the atom by symmetry.
  module subroutine add_atom(c,is,x,iunit_l,isnneq,ti)
    use crystalseedmod, only: crystalseed
    use types, only: realloc
    use global, only: iunit_ang, iunit_bohr
    use tools_io, only: nameguess
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    integer, intent(in) :: is
    real*8, intent(in) :: x(3)
    integer, intent(in) :: iunit_l
    logical, intent(in) :: isnneq
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: seed
    real*8 :: xx(3)
    logical :: copysym
    integer :: is_

    ! whether to use symmetry
    copysym = isnneq .and. .not.c%ismolecule .and. c%spgavail

    ! make seed from this crystal
    call c%makeseed(seed,copysym=copysym)

    ! interpret units
    if (iunit_l == iunit_ang) then
       xx = x / bohrtoa
       xx = c%c2x(xx)
    elseif (iunit_l == iunit_bohr) then
       xx = c%c2x(x)
    else
       xx = x
    end if

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
    seed%x(:,seed%nat) = xx
    seed%is(seed%nat) = is_
    seed%atname(seed%nat) = seed%spc(is_)%name

    ! build the new crystal
    call c%struct_new(seed,crashfail=.true.,ti=ti)

  end subroutine add_atom

  ! Remove the bond between cell atoms iat1 and iat2 (iat2 at lattice vector
  ! lvec relative to iat1) by editing the neighbor stars in place. Removes both
  ! reciprocal entries (iat2 at lvec from iat1, iat1 at -lvec from iat2) and
  ! recomputes the molecular-fragment data; does NOT rebuild the structure.
  module subroutine remove_bond(c,iat1,iat2,lvec)
    use types, only: neighstar, realloc
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)

    if (.not.allocated(c%nstar)) return
    if (iat1 < 1 .or. iat1 > c%ncel .or. iat2 < 1 .or. iat2 > c%ncel) return

    ! remove iat2 (at lvec) from iat1's star and the reciprocal iat1 (at -lvec)
    ! from iat2's star
    call remove_one(c%nstar(iat1),iat2,lvec)
    call remove_one(c%nstar(iat2),iat1,-lvec)

    ! keep molecular fragments consistent with the new connectivity (lightweight;
    ! same post-asterism steps the rebond window runs, no struct_new)
    call c%fill_molecular_fragments()
    call c%calculate_molecular_equivalence()
    call c%calculate_periodicity()

  contains
    subroutine remove_one(ns,id,lv)
      type(neighstar), intent(inout) :: ns
      integer, intent(in) :: id, lv(3)
      integer :: k, knew
      knew = 0
      do k = 1, ns%ncon
         if (ns%idcon(k) == id .and. all(ns%lcon(:,k) == lv)) cycle
         knew = knew + 1
         if (knew /= k) then
            ns%idcon(knew) = ns%idcon(k)
            ns%lcon(:,knew) = ns%lcon(:,k)
            ns%ordcon(knew) = ns%ordcon(k)
            ns%aromdir(:,knew) = ns%aromdir(:,k)
         end if
      end do
      if (knew /= ns%ncon) then
         ns%ncon = knew
         call realloc(ns%idcon,knew)
         call realloc(ns%lcon,3,knew)
         call realloc(ns%ordcon,knew)
         call realloc(ns%aromdir,3,knew)
      end if
    end subroutine remove_one
  end subroutine remove_bond

  ! Set the bond order (ordcon) of the bond between cell atoms iat1 and iat2
  ! (iat2 at lattice vector lvec relative to iat1). Updates both reciprocal
  ! entries in the neighbor stars in place; does not change the connectivity.
  ! Order convention: 0=dashed, 1=single, 2=double, 3=triple.
  module subroutine set_bond_order(c,iat1,iat2,lvec,order)
    use types, only: neighstar
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)
    integer, intent(in) :: order

    if (.not.allocated(c%nstar)) return
    if (iat1 < 1 .or. iat1 > c%ncel .or. iat2 < 1 .or. iat2 > c%ncel) return

    call set_one(c%nstar(iat1),iat2,lvec)
    call set_one(c%nstar(iat2),iat1,-lvec)

  contains
    subroutine set_one(ns,id,lv)
      type(neighstar), intent(inout) :: ns
      integer, intent(in) :: id, lv(3)
      integer :: k
      do k = 1, ns%ncon
         if (ns%idcon(k) == id .and. all(ns%lcon(:,k) == lv)) ns%ordcon(k) = order
      end do
    end subroutine set_one
  end subroutine set_bond_order

  ! Add a bond between cell atoms iat1 and iat2 + lvec with the given
  ! bond order by editing the neighbor stars in place. Adds both
  ! reciprocal entries (iat2 at lvec from iat1, iat1 at -lvec from
  ! iat2) and recomputes the molecular-fragment data; does NOT rebuild
  ! the structure. Self-bonds and already-existing bonds are ignored.
  module subroutine add_bond(c,iat1,iat2,lvec,order)
    use types, only: neighstar, realloc
    class(crystal), intent(inout) :: c
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)
    integer, intent(in) :: order

    integer :: k

    if (iat1 < 1 .or. iat1 > c%ncel .or. iat2 < 1 .or. iat2 > c%ncel) return
    if (iat1 == iat2 .and. all(lvec == 0)) return
    if (.not.allocated(c%nstar)) allocate(c%nstar(c%ncel))

    ! ignore the bond if it is already there
    do k = 1, c%nstar(iat1)%ncon
       if (c%nstar(iat1)%idcon(k) == iat2 .and. all(c%nstar(iat1)%lcon(:,k) == lvec)) return
    end do

    ! add iat2 (at lvec) to iat1's star and the reciprocal iat1 (at -lvec)
    ! to iat2's star
    call add_one(c%nstar(iat1),iat2,lvec)
    call add_one(c%nstar(iat2),iat1,-lvec)

    ! keep molecular fragments consistent with the new connectivity (lightweight;
    ! same post-asterism steps the rebond window runs, no struct_new)
    call c%fill_molecular_fragments()
    call c%calculate_molecular_equivalence()
    call c%calculate_periodicity()

  contains
    subroutine add_one(ns,id,lv)
      type(neighstar), intent(inout) :: ns
      integer, intent(in) :: id, lv(3)
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
      ns%ordcon(n) = order
      ns%aromdir(:,n) = 0d0
      ns%ncon = n
    end subroutine add_one
  end subroutine add_bond

end submodule edit
