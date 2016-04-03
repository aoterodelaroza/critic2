! Copyright (c) 2015 Alberto Otero de la Roza
! <alberto@fluor.quimica.uniovi.es>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.

! Find the symmetry operations of a finite molecule (from tessel)
module sympg
  implicit none

  private

  ! the only public routine is sym3d
  public :: sym3d

  ! public interface 
  public :: point_group
  public :: nopsym
  public :: oporder
  public :: opproper
  public :: opsym
  public :: nclas
  public :: ninclas
  public :: iinclas
  public :: opsymbol

  ! Symmetry operators and point group symmetry of the
  ! molecular motif (sym.inc).
  integer, parameter  :: MOPSYM = 480 !< Number of symmetry operations
  integer :: nopsym !< Number of symmetry operations.
  integer :: optype(MOPSYM) !< Symmetry operation types.
  integer :: oporder(MOPSYM) !< Rotation order (the n in C_n^m or S_n^m).
  integer :: opm(MOPSYM) !< The m in C_n^m or S_n^m.
  integer :: opinv(MOPSYM) !< Inverse of an operation.
  integer :: optable(MOPSYM,MOPSYM) !< Multiplication (Cayley) table.
  integer :: opmainord !< Order of the main axis.
  real*8  :: opsym(MOPSYM,3,3) !< The 3x3 matrix of the operation in the xyz rep.
  real*8  :: opaxis(MOPSYM,3) !< Rotation axis.
  real*8  :: opeuler(MOPSYM,3) !< Euler angles for the rotation axis.
  real*8  :: opangle(MOPSYM) !< Rotation angle (radians).
  character*(10) :: opsymbol(MOPSYM) !< Symbol for the operation.
  character*(8)  :: point_group !< Point group symbol.
  character*(8)  :: point_g2 !< Point group symbol (short version).
  logical :: opproper(MOPSYM) !< Proper/improper rotation (redundant).
  logical :: opisgener(MOPSYM) !< Is this a generator for the group?
  logical :: linear_mol !< True if the molecule is linear.

  !.Classes of symmetry operations:
  integer, parameter :: mclas = 12 !< classes of symmetry operations
  integer, parameter :: minclas = 20 !< classes of symmetry operations
  integer :: nclas !< number of classes
  integer :: opclas(MOPSYM)
  integer :: ninclas(MCLAS)
  integer :: iinclas(MCLAS,MINCLAS)

  ! Sym operator types:
  integer, parameter :: opt_identity     = 10 !< identity op.
  integer, parameter :: opt_rotation     = 11 !< rotation op.
  integer, parameter :: opt_inversion    = 12 !< inversion op.
  integer, parameter :: opt_sigma        = 13 !< plane op.
  integer, parameter :: opt_imp_rotation = 14 !< improper rotation op.
  integer, parameter :: opt_unknown      = 15 !< unknown op.

  !.Classification of atoms in orbits of equivalent-by-symmetry
  ! inidividuals:
  integer, parameter :: MMOL1 = 3000 !< max. orbits
  integer, allocatable :: orbmol(:) !< orbit to which an atom belongs.
  real*8  :: molradius !< maximum distance from the center to an atom.

  !.Tolerances that are required for the symmetry
  ! determination algorithm. Change the defaults with care.
  ! (symtol.inc)
  real*8, parameter :: TOLdist  = 3d-5 !< Check of a distance.
  real*8, parameter :: TOLeqvm  = 2d0 * TOLdist !< Comparison of matrices.
  real*8, parameter :: TOLsng   = 1d-12 !< Singularity of a 3x3 matrix.
  real*8, parameter :: TOLisint = 3d-5 !< Closeness to an integer.
  real*8, parameter :: TOLnull  = 3d-6 !< Null angles and projections.
  real*8, parameter :: TOLeigen = 3d-5 !< Eigenvalues of the orthogonal matrices.
  logical, parameter :: TOLdirty = .false. !< Activate the algorithm changes to deal with noisy data.

  !.The largest errors on the different tests.
  real*8 :: ERRsng !< Best nonsingular triplet found.
  real*8 :: ERReigen !< Error on the eigenvalues.

contains
  
  !> Check symmetry of nonlinear molecules.
  !> Adapted from TESSEL.
  subroutine sym3d (rmat, ax, ay, az, atZmol, nmol, verbose)
    use tools_io

    real*8, intent(in) :: rmat(3,3)
    integer   :: nmol, atZmol(nmol)
    real*8    :: ax(nmol), ay(nmol), az(nmol)
    logical, intent(in) :: verbose

    integer   :: k, ii, jj, i1, i2, i3, j1, j2, j3, ierror
    real*8    :: xmat(3,3), xm(3,3), xv(3,3), xop(3,3), xdet
    integer   :: ii1, ii2, ii3
    real*8    :: dbest
    integer*4 :: ntest1, ntest2

    !.Classify all atoms into orbits. All atoms in an orbit have the
    !     same atomic number and the same distance to the center of mass:
    call symorb (ax, ay, az, atZmol, nmol, verbose)

    !.The identity is always a sym operator:
    !
    linear_mol = .false.
    nopsym = 0
    call symmatfill (xmat, 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
    call symopadd (xmat)

    ! Is the inversion a sym op?
    !
    call symmatfill (xmat, -1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0)
    if (symopchk(rmat,xmat,ax,ay,az,atZmol,nmol)) call symopadd (xmat)

    !.Find a linearly independent triplet of atoms.
    ! Try first to find a very good triplet, or use the best available
    ! one anyway.
    ii1 = 0
    ii2 = 0
    ii3 = 0
    i1 = 0
    i2 = 0
    i3 = 0
    dbest = 0d0
    do i1 = 1, nmol
       xmat(1,1) = ax(i1)
       xmat(2,1) = ay(i1)
       xmat(3,1) = az(i1)
       do i2 = i1+1, nmol
          xmat(1,2) = ax(i2)
          xmat(2,2) = ay(i2)
          xmat(3,2) = az(i2)
          do i3 = i2+1, nmol
             xmat(1,3) = ax(i3)
             xmat(2,3) = ay(i3)
             xmat(3,3) = az(i3)
             call syminv (xmat, xv, xdet, ierror)
             if (ierror .eq. 0) then
                if (abs(xdet) .gt. 1d-1) goto 1001
                if (abs(xdet) .gt. dbest) then
                   dbest = abs(xdet)
                   ii1 = i1
                   ii2 = i2
                   ii3 = i3
                endif
             end if
          enddo ! i3
       enddo ! i2
    enddo ! i1
    if (ii1 == 0 .or. ii2 == 0 .or. ii3 == 0) then
       call ferror('sym3d','could not find atom triplet',faterr)
    end if
    xmat(1,1) = ax(ii1)
    xmat(2,1) = ay(ii1)
    xmat(3,1) = az(ii1)
    xmat(1,2) = ax(ii2)
    xmat(2,2) = ay(ii2)
    xmat(3,2) = az(ii2)
    xmat(1,3) = ax(ii3)
    xmat(2,3) = ay(ii3)
    xmat(3,3) = az(ii3)
    call syminv (xmat, xv, xdet, ierror)
    if (ierror .ne. 0) then
       call ferror('sym3d','singular triplet matrix',faterr)
    end if

    !.Run over all triplets that are compatible with the linear
    ! independent triplet already found. Each compatible triplet migth
    ! produce a new symmetry operator. To be compatible, both triplets
    ! must be formed by the same atoms in the same order.

1001 continue
    ERRsng = abs(xdet)
    ntest1 = 0
    ntest2 = 0
    do j1 = 1, nmol
       if (orbmol(i1).ne.orbmol(j1)) goto 1002
       xm(1,1) = ax(j1)
       xm(2,1) = ay(j1)
       xm(3,1) = az(j1)
       do j2 = 1, nmol
          if (orbmol(i2).ne.orbmol(j2)) goto 1003
          if (j1.eq.j2) goto 1003
          xm(1,2) = ax(j2)
          xm(2,2) = ay(j2)
          xm(3,2) = az(j2)
          do j3 = 1, nmol
             if (orbmol(i3).ne.orbmol(j3)) goto 1004
             if (j1.eq.j3 .or. j2.eq.j3) goto 1004
             ntest1 = ntest1 + 1
             xm(1,3) = ax(j3)
             xm(2,3) = ay(j3)
             xm(3,3) = az(j3)
             do ii = 1, 3
                do jj = 1, 3
                   xop(ii,jj) = 0d0
                   do k = 1, 3
                      xop(ii,jj) = xop(ii,jj) + xm(ii,k) * xv(k,jj)
                   enddo
                enddo
             enddo

             !.Check if this is a new sym operator:
             !
             if (symopchk(rmat,xop,ax,ay,az,atZmol,nmol)) then
                ntest2 = ntest2 + 1
                call symopadd (xop)
                if (TOLdirty) call symclosure(verbose)
             endif
1004         continue
          enddo  !j3
1003      continue
       enddo  !j2
1002   continue
    enddo  !j1

    if (verbose) &
       write (uout,600) ntest1, ntest2
600 format (&
       1x, 'ALG(sym3d) Triplet pairs tested:     ', i12/&
       1x, 'ALG(sym3d) Possible operators found: ', i12)

    !.Check the closure of the symmetry operators set:
    call symclosure(verbose)
    call symgetgroup(verbose)

    ! clean up
    if (allocated(orbmol)) deallocate(orbmol)

  end subroutine sym3d

  !> Check the closure of the symmetry operators set
  !> and get the properties of the symmetry operators.
  !> Adapted from TESSEL.
  subroutine symclosure (verbose)
    use tools_io

    logical, intent(in) :: verbose

    real*8         ::  xmat(3,3)
    integer        ::  i, j, ii, jj, kk, inum, old_nopsym, ini_nopsym
    logical        ::  newop, found, generated(MOPSYM)
    integer        ::  iop, jinv, iclas

    !.Initialize the Cayley table:
    !
    do i = 1, MOPSYM
       do j = 1, MOPSYM
          optable(i,j) = 0
       enddo
    enddo

    ! Build up the Cayley table and look for new operators at the same
    ! time:
    ini_nopsym = nopsym
    newop = .true.
    do while (newop)
       newop = .false.
       old_nopsym = nopsym
       do i = 1, old_nopsym
          do j = 1, old_nopsym
             if (optable(i,j) .le. 0) then

                !.Get the product of operators i x j and determine if
                ! this is a new or an old symmetry operator:
                !
                do ii = 1, 3
                   do jj = 1, 3
                      xmat(ii,jj) = 0d0
                      do kk = 1, 3
                         xmat(ii,jj) = xmat(ii,jj)&
                            + opsym(i,ii,kk) * opsym(j,kk,jj)
                      enddo
                   enddo
                enddo
                inum = symnumber(xmat)
                if (inum .le. 0) then
                   newop = .true.
                   optable(i,j) = nopsym + 1
                   call symopadd (xmat)
                else
                   optable(i,j) = inum
                endif
             endif
          enddo
       enddo
    enddo

    !.Determine the operation inverses:
    do i = 1, nopsym
       j = 1
       found = .false.
       do while (.not.found .and. j.le.nopsym)
          kk = optable(i,j)
          if (optype(kk) .eq. opt_identity) then
             opinv(i) = j
             found = .true.
          else
             j = j + 1
          endif
       enddo
       if (.not.found .and. verbose) write (uout,605) i
    enddo

    !.Check for a set of generators:
    !
    do i = 1, nopsym
       generated(i) = .false.
       opisgener(i) = .false.
    enddo
    do i = 1, nopsym
       if (.not.generated(i) .and. optype(i).ne.opt_identity) then
          opisgener(i) = .true.
       endif
       do j = 1, i
          kk = optable(i,j)
          generated(kk) = .true.
          kk = optable(j,i)
          generated(kk) = .true.
       enddo
    enddo

    !.Determine the classes of the symmetry operations:
    ! (Non-essential)
    ! A and B belong in the same class if there exist any operation X in
    ! the group such that: X^{-1} A X = B.
    !
    do i = 1, nopsym
       opclas(i) = -1
    enddo
    nclas = 0
    do i = 1, nopsym
       if (opclas(i) .le. 0) then

          !.This op belongs to a new class:
          !
          nclas = nclas + 1
          if (nclas.gt.MCLAS) call ferror ('symclosure',&
             'Too many classes. Increase MCLAS!', faterr)
          opclas(i) = nclas
          iclas = nclas
          ninclas(nclas) = 1
          iinclas(nclas,1) = i
       else

          !.The op belongs to a class already found:
          !
          iclas = opclas(i)
       endif

       !.Apply similarity transforms to get all ops in the same class:
       !
       do j = 1, nopsym
          jinv = opinv(j)
          iop = optable(optable(jinv, i), j)
          if (opclas(iop) .le. 0) then
             opclas(iop) = iclas
             ninclas(nclas) = ninclas(nclas) + 1
             if (ninclas(nclas).gt.MINCLAS) call ferror ('symclosure',&
                'Too many ops in a class. Increase MINCLAS!', faterr)
             iinclas(nclas,ninclas(nclas)) = iop
          else if (opclas(iop) .ne. iclas) then

             !.This is an error: iop and i should be in the same class,
             ! but they are not.
             if (verbose) write (uout,620) i, iop
          endif
       enddo
    enddo

    if (nopsym.gt.ini_nopsym .and. uout.ge.0 .and. verbose) then
       write (uout,600)
       write (uout,602) nopsym-ini_nopsym
    endif

600 format (/&
       1x, '++SYMCLOSURE: Complete the group by multiplying the',&
       1x, 'matrices already known.')
602 format (&
       1x, 'ALG(symclosure) ', i3,&
       ' new operations found because of the group closure')
605 format (&
       1x, 'DBG(symclosure): Inverse not found for opsym ', i3)
620 format (&
       1x, 'WARNING(symclosure): The ops', 2i5, ' should be in the',&
       1x, 'same class but they are not. We continue anyway!')

  end subroutine symclosure

  !> Classify all atoms into orbits. All atoms in an orbit
  !> have the same atomic number and the same distance to the center.
  !> Adapted from TESSEL.
  subroutine symorb (ax, ay, az, atZmol, nmol, verbose)
    use tools_io
    use types

    integer :: nmol, atZmol(nmol)
    real*8  :: ax(nmol), ay(nmol), az(nmol)
    logical, intent(in) :: verbose

    real*8  :: dis2
    integer :: i, j, iorb

    integer, parameter :: MORBIT = 20 !< max. orbits
    integer, parameter :: MATORB = 300 !< max. atomic orbits
    integer :: norbit !< number of orbits
    integer, allocatable :: natorb(:) !< number of atoms in an orbit.
    integer, allocatable :: iatorb(:,:) !< list of atoms in an orbit.
    integer, allocatable :: orbZ(:) !< atomic number of the atoms in this orbit.
    real*8, allocatable :: orbdis(:) !< distance to the center for all atoms in this orbit.

    ! allocate
    if (allocated(natorb)) deallocate(natorb)
    if (allocated(iatorb)) deallocate(iatorb)
    if (allocated(orbZ)) deallocate(orbZ)
    if (allocated(orbdis)) deallocate(orbdis)
    if (allocated(orbmol)) deallocate(orbmol)
    allocate(natorb(MORBIT),iatorb(MORBIT,MATORB),orbZ(MORBIT),orbdis(MORBIT))
    allocate(orbmol(MMOL1))

    !.Classify all atoms into orbits. All atoms in an orbit have the
    ! same atomic number and the same distance to the center of mass:
    norbit = 0
    molradius = 0d0
    do i = 1, nmol
       dis2 = sqrt(ax(i)*ax(i) + ay(i)*ay(i) + az(i)*az(i))
       if (dis2 .gt. molradius) molradius = dis2
       iorb = -1
       j = 1
       do while (j.le.norbit .and. iorb.lt.0)
          if (abs(dis2-orbdis(j)) .le. TOLdist .and. atZmol(i) .eq. orbZ(j)) then
             iorb = j
          endif
          j = j + 1
       enddo
       if (iorb .lt. 0) then
          norbit = norbit + 1
          if (norbit > MORBIT) then
             call realloc(natorb,2*norbit)
             call realloc(iatorb,2*norbit,size(iatorb,2))
             call realloc(orbZ,2*norbit)
             call realloc(orbdis,2*norbit)
          end if
          orbdis(norbit) = dis2
          orbZ(norbit) = atZmol(i)
          natorb(norbit) = 0
          iorb = norbit
       endif
       if (i > size(orbmol)) call realloc(orbmol,2*i)
       orbmol(i) = iorb
       natorb(iorb) = natorb(iorb) + 1
       if (natorb(iorb) > size(iatorb,2)) &
          call realloc(iatorb,size(iatorb,1),2*natorb(iorb))
       iatorb(iorb,natorb(iorb)) = i
    enddo
    if (verbose) write (uout,550) norbit
    if (verbose) then
       do i = 1, norbit
          write (uout,555) i, natorb(i), orbZ(i), orbdis(i),&
             (iatorb(i,j), j = 1, natorb(i))
       enddo
    end if

    ! cleanup
    deallocate(natorb,iatorb,orbZ,orbdis)

550 format (/&
       '+ SYMORB: Determining the atomic orbits.'/&
       '  ALG(symorb) Number of orbits: ', i5)
555 format (&
       '  ALG(symorb) Orbit ', i5/&
       '  ALG(symorb) Number of atoms: ', i5/&
       '  ALG(symorb) Atomic number and distance: ', i5, f18.9/&
       '  ALG(symorb) List of atoms in this orbit:'/&
       (1x, 10i6))

  end subroutine symorb

  !> If xmat is yet unknown add it to the list of sym ops.
  !> Adapted from TESSEL.
  subroutine symopadd (xmat)
    use tools_io, only: ferror, faterr, warning
    use param

    real*8            xmat(3,3)

    integer, parameter :: MOP = 20

    integer           i, j, k
    real*8            diff

    integer           ierror, ii, ibest
    logical           found
    real*8            xtrace, xdet, xnm, ang1, angle
    real*8            err, besterr
    real*8            xinv(3,3), xval(3), xvect(3,3), xeul(3)
    character*(40)    fmt

    do i = 1, nopsym
       diff = 0d0
       do j = 1, 3
          do k = 1, 3
             diff = diff + abs(xmat(j,k) - opsym(i,j,k))
          enddo
       enddo
       if (diff .le. TOLeqvm) return  ! This is not a new operator
    enddo
    nopsym = nopsym + 1
    if (nopsym .gt. MOPSYM) call ferror ('symopadd',&
       'Too many sym operators. Increase MOPSYM', faterr)
    do j = 1, 3
       do k = 1, 3
          opsym(nopsym,j,k) = xmat(j,k)
       enddo
    enddo

    ! Get the properties of the symmetry matrix:
    !
    xtrace = 0d0
    do i = 1, 3
       xtrace = xtrace + xmat(i,i)
    enddo
    call syminv (xmat, xinv, xdet, ierror)
    if (ierror .lt. 0) call ferror ('symopadd',&
       'Operator has a singular matrix!!!', faterr)

    !.Get the Euler angles for the rotation:
    !
    call symeuler (xmat, xeul)
    opeuler(nopsym,1) = xeul(1)
    opeuler(nopsym,2) = xeul(2)
    opeuler(nopsym,3) = xeul(3)

    ! Is this a proper or improper rotation?
    !
    if (abs(abs(xdet)-1d0) .gt. TOLisint) then
       call ferror ('symopadd', 'determinant is not +-1', warning)
    else if (xdet .gt. 0d0) then
       opproper(nopsym) = .true.
    else
       opproper(nopsym) = .false.
    endif

    !.Get the rotation axis except for E and i, for which it is
    ! fully degenerated:
    if (abs(abs(xtrace)-3d0) .le. TOLisint) then
       if (xtrace .gt. 0d0) then
          opsymbol(nopsym) = 'E'
          opaxis(nopsym,1) = 0d0
          opaxis(nopsym,2) = 0d0
          opaxis(nopsym,3) = 1d0
          oporder(nopsym) = 1
          opm(nopsym) = 1
          opangle(nopsym) = pi + pi
          optype(nopsym) = opt_identity
       else
          opsymbol(nopsym) = 'i'
          opaxis(nopsym,1) = 0d0
          opaxis(nopsym,2) = 0d0
          opaxis(nopsym,3) = 1d0
          oporder(nopsym) = 2
          opm(nopsym) = 1
          opangle(nopsym) = pi
          optype(nopsym) = opt_inversion
       endif
    else

       !.Get the rotation angle and the n and m of the C_n^m or S_n^m
       ! symmetry operator:
       !
       ang1 = 0.5d0*(xtrace-xdet)
       if (abs(ang1) .gt. 1d0) ang1 = sign(1d0, ang1)
       angle = acos(ang1)
       if (abs(angle) .le. TOLnull) then
          xnm = 1d0
          angle = pi + pi
       else
          xnm = 2d0 * pi / angle
       endif
       opangle(nopsym) = angle
       ii = 0
       found = .false.
       do while (.not. found .and. ii.le.MOP)
          ii = ii + 1
          found = abs(xnm*ii - nint(xnm*ii)) .le. TOLisint
       enddo
       if (found) then
          oporder(nopsym) = nint(xnm*ii)
          opm(nopsym) = ii
       else
          ibest = 1
          besterr = abs(xnm - nint(xnm))
          do ii = 2, MOP
             err = abs(xnm*ii - nint(xnm*ii))
             if (err .le. besterr) then
                besterr = err
                ibest = ii
             endif
          enddo
          oporder(nopsym) = nint(xnm*ibest)
          opm(nopsym) = ibest
       endif
       write (fmt,200) 1+int(log10(1d0*oporder(nopsym))), 1+int(log10(1d0*opm(nopsym)))
       !
       call symeigen (xmat, xval, xvect, xdet, ierror)
       if (ierror .ne. 0) then
          call ferror ('symopadd','Trouble finding the rotation axis', warning)
       endif
       opaxis(nopsym,1) = xvect(1,3)
       opaxis(nopsym,2) = xvect(2,3)
       opaxis(nopsym,3) = xvect(3,3)
       if (xdet .gt. 0d0) then
          write (opsymbol(nopsym),fmt) 'C', oporder(nopsym), opm(nopsym)
          optype(nopsym) = opt_rotation
       else if (abs(xtrace - 1d0) .le. TOLisint) then
          write (opsymbol(nopsym),210) 'sigma'
          optype(nopsym) = opt_sigma
       else if (oporder(nopsym).eq.1 .and. opm(nopsym).eq.1) then
          write (opsymbol(nopsym),210) 'sigma'
          optype(nopsym) = opt_sigma
       else
          write (opsymbol(nopsym),fmt) 'S', oporder(nopsym), opm(nopsym)
          optype(nopsym) = opt_imp_rotation
       endif
    endif
    !
200 format ('(a, "_", i', i2, ', "^", i', i2, ')')
210 format (a)
    !
  end subroutine symopadd

  !> Check if xmat is a sym operator.
  !> Adapted from TESSEL.
  function symopchk(rmat,xmat,ax,ay,az,atgroup,nmol)
    use tools_math

    logical :: symopchk
    real*8, intent(in) :: rmat(3,3)
    integer :: nmol, atgroup(nmol)
    real*8  :: xmat(3,3), ax(nmol), ay(nmol), az(nmol)

    integer :: i, j, k
    logical :: found
    real*8  :: xii, xij, xji, rmati(3,3)
    real*8  :: ximg(3), diff, xdif(3)

    !.First check: The matrix must be orthogonal:
    symopchk = .false.
    do i = 1, 3
       xii = xmat(1,i)*xmat(1,i) + xmat(2,i)*xmat(2,i) + xmat(3,i)*xmat(3,i)
       if (abs(xii-1d0) .gt. TOLeqvm) return
       do j = i+1, 3
          xij = xmat(1,i)*xmat(1,j) + xmat(2,i)*xmat(2,j) + xmat(3,i)*xmat(3,j)
          if (abs(xij) .gt. TOLeqvm) return
          xji = xmat(1,j)*xmat(1,i) + xmat(2,j)*xmat(2,i) + xmat(3,j)*xmat(3,i)
          if (abs(xji) .gt. TOLeqvm) return
       enddo
    enddo

    !.If the number of atoms is large, perhaps is better to check
    ! for repeated operations of symmetry before testing if it
    ! really transforms the molecule into itself:
    if (nmol .gt. nopsym+nopsym) then
       do i = 1, nopsym
          diff = 0d0
          do j = 1, 3
             do k = 1, 3
                diff = diff + abs(xmat(j,k) - opsym(i,j,k))
             enddo
          enddo
          if (diff .le. TOLeqvm) return  ! This is not a new operator
       enddo
    endif

    !.Transform every atom and check for the symmetric image:
    rmati = matinv(rmat)
    do i = 1, nmol
       ximg = (/ax(i), ay(i), az(i)/)
       ximg = matmul(xmat, ximg)
       found = .false.
       j = 1
       do while (.not.found .and. j.le.nmol)
          if (atgroup(i).eq.atgroup(j)) then
             xdif = ximg - (/ax(j),ay(j),az(j)/)
             xdif = matmul(xdif, rmati)
             xdif = xdif - nint(xdif)
             xdif = matmul(xdif, rmat)
             found = sqrt(dot_product(xdif,xdif)) .le. TOLdist
          endif
          j = j + 1
       enddo
       if (.not.found) return
    enddo
    symopchk = .true.

  end function symopchk

  !> Get the inverse of a 3x3 nonsingular matrix.
  !> The input matrix is not modified.
  !> Adapted from TESSEL.
  subroutine syminv (xmat, xinv, xdet, ierror)
    use tools_math

    real*8  :: xmat(3,3), xinv(3,3), xdet
    integer :: ierror

    integer :: i, j
    real*8  :: adj(3,3), det0, ddet, norm0

    norm0 = (norm(xmat(:,1)) + norm(xmat(:,2)) + norm(xmat(:,3)))/3d0
    ierror = 0
    adj(1,1) = xmat(2,2)*xmat(3,3) - xmat(2,3)*xmat(3,2)
    adj(1,2) = xmat(2,3)*xmat(3,1) - xmat(2,1)*xmat(3,3)
    adj(1,3) = xmat(3,2)*xmat(2,1) - xmat(2,2)*xmat(3,1)
    det0 = xmat(1,1)*adj(1,1) + xmat(1,2)*adj(1,2) + xmat(1,3)*adj(1,3)
    if (abs(det0/norm0) .le. TOLsng) then
       ierror = -1
       return
    endif
    adj(2,1) = xmat(1,3)*xmat(3,2) - xmat(1,2)*xmat(3,3)
    adj(2,2) = xmat(1,1)*xmat(3,3) - xmat(1,3)*xmat(3,1)
    adj(2,3) = xmat(1,2)*xmat(3,1) - xmat(1,1)*xmat(3,2)
    adj(3,1) = xmat(1,2)*xmat(2,3) - xmat(1,3)*xmat(2,2)
    adj(3,2) = xmat(1,3)*xmat(2,1) - xmat(1,1)*xmat(2,3)
    adj(3,3) = xmat(1,1)*xmat(2,2) - xmat(1,2)*xmat(2,1)
    ddet = 1d0/det0
    do i = 1, 3
       do j = 1, 3
          xinv(i,j) = adj(j,i) * ddet
       enddo
    enddo
    xdet = det0
    return

  end subroutine syminv

  !> Fill in the xmat matrix.
  !> Adapted from TESSEL.
  subroutine symmatfill (xmat,a11, a12, a13, a21, a22, a23, a31, a32, a33)

    real*8, intent(out) :: xmat(3,3)
    real*8, intent(in) :: a11, a12, a13, a21, a22, a23, a31, a32, a33

    xmat(1,1) = a11
    xmat(1,2) = a12
    xmat(1,3) = a13
    xmat(2,1) = a21
    xmat(2,2) = a22
    xmat(2,3) = a23
    xmat(3,1) = a31
    xmat(3,2) = a32
    xmat(3,3) = a33

  end subroutine symmatfill

  !> Determine the point group from the set of symmetry
  !> operators.
  !> Adapted from TESSEL.
  subroutine symgetgroup (verbose)
    use tools_io, only: uout, ferror, faterr, warning

    logical, intent(in) :: verbose

    integer, parameter :: MOP = 20

    integer       :: mainorder, mainmult, imain(MOP), mainimproper
    integer       :: binmult, ibin(MOP), nsigma, isigma(MOP), isigmah
    integer       :: nbinperp, i, j, j1
    logical       :: inversion, sigmah
    real*8        :: xprod
    character*(4) :: chmain, chimp

    !.Get the main axis, if any, and their multiplicity:
    ! Only the C_n^1 operations will be considered.
    ! Unfortunately, C_n^1 and C_n^(n-1) operations cannot be
    ! distinguished.
    ! Get also the list of binary axes, the sigma planes, and the
    ! order of the main improper axis (n > 2).
    mainorder = 0
    mainmult = 0
    mainimproper = 0
    binmult = 0
    inversion = .false.
    nsigma = 0
    do i = 1, nopsym
       if (optype(i) .eq. opt_rotation .and. opm(i) .eq. 1) then
          if (oporder(i) .gt. mainorder) then
             mainorder = oporder(i)
             mainmult = 1
             imain(mainmult) = i
          else if (oporder(i) .eq. mainorder) then
             mainmult = mainmult + 1
             if (mainmult .gt. MOP) call ferror ('symgetgroup',&
                'Main axis mult too high. Increase MOP!', faterr)
             imain(mainmult) = i
          endif
          if (oporder(i) .eq. 2) then
             binmult = binmult + 1
             if (binmult .gt. MOP) call ferror ('symgetgroup',&
                'Binary axis mult too high. Increase MOP!', faterr)
             ibin(binmult) = i
          endif
       else if (optype(i) .eq. opt_inversion) then
          inversion = .true.
       else if (optype(i) .eq. opt_sigma) then
          nsigma = nsigma + 1
          if (nsigma .gt. MOP) call ferror ('symgetgroup',&
             'Too many sigma planes. Increase MOP!', faterr)
          isigma(nsigma) = i
       else if (optype(i) .eq. opt_imp_rotation) then
          if (oporder(i) .gt. mainimproper) mainimproper = oporder(i)
       endif
    enddo
    if (mainorder .gt. 2) mainmult = mainmult/2

    ! If there is a single main order group, look for sigma_h planes:
    !
    sigmah = .false.
    isigmah = -1
    if (mainmult .eq. 1 .or. mainorder .eq. 2) then
       i = imain(1)
       j = 1
       do while (.not.sigmah .and. j .le. nsigma)
          j1 = isigma(j)
          xprod = opaxis(i,1)*opaxis(j1,1)&
             + opaxis(i,2)*opaxis(j1,2)&
             + opaxis(i,3)*opaxis(j1,3)
          if (abs(abs(xprod)-1d0) .le. TOLisint) then
             sigmah = .true.
             isigmah = j1
          endif
          j = j + 1
       enddo
    endif
    if (sigmah .and. isigmah > 0) then
       opsymbol(isigmah) = 'sigma_h'
    endif

    !.If there is a single main order group, look for n C2 perpendicular
    ! axis:
    nbinperp = 0
    if (mainmult .eq. 1 .or. mainorder .eq. 2) then
       i = imain(1)
       do j = 1, binmult
          j1 = ibin(j)
          xprod = opaxis(i,1)*opaxis(j1,1)&
             + opaxis(i,2)*opaxis(j1,2)&
             + opaxis(i,3)*opaxis(j1,3)
          if (abs(xprod) .le. TOLnull) then
             nbinperp = nbinperp + 1
          endif
       enddo
    endif

    !.Store as a character the main proper and improper orders:
    !
    if (mainorder .gt. 999) then
       write (chmain, '(i4)') mainorder
    else if (mainorder .gt. 99) then
       write (chmain, '(i3)') mainorder
    else if (mainorder .gt. 9) then
       write (chmain, '(i2)') mainorder
    else
       write (chmain, '(i1)') mainorder
    endif
    if (mainimproper .gt. 999) then
       write (chimp, '(i4)') mainimproper
    else if (mainimproper .gt. 99) then
       write (chimp, '(i3)') mainimproper
    else if (mainimproper .gt. 9) then
       write (chimp, '(i2)') mainimproper
    else
       write (chimp, '(i1)') mainimproper
    endif
    opmainord = mainorder

    !.Decision tree to get the point group. Start with the cubic-like
    ! groups:
    if (mainmult .gt. 1 .and. mainorder .gt. 2) then
       if (mainorder .eq. 5 .and. mainmult .ge. 6) then
          if (inversion) then
             point_group = 'Ih'
          else
             point_group = 'I'
          endif
          point_g2 = point_group
          opmainord = 5
       else if (mainorder .eq. 4 .and. mainmult .ge. 3) then
          if (inversion) then
             point_group = 'Oh'
          else
             point_group = 'O'
          endif
          point_g2 = point_group
          opmainord = 4
       else if (mainorder .eq. 3 .and. mainmult .ge. 4) then
          if (inversion) then
             point_group = 'Th'
          else if (nsigma .eq. 6) then
             point_group = 'Td'
          else
             point_group = 'T'
          endif
          point_g2 = point_group
          opmainord = 3
       else
          call ferror ('symgetgroup', 'Unknown cubic-like group',warning)
          if (verbose) &
             write (uout,700) mainorder, mainmult, inversion, binmult&
             , nbinperp, nsigma, sigmah, mainimproper
          point_group = '??'
       endif
    else if (mainorder .ge. 2) then
       if (mainorder .eq. nbinperp) then
          if (sigmah) then
             point_group = "D" // trim(chmain) // "h"
             point_g2 = 'Dnh'
          else if (mainorder .eq. nsigma) then
             point_group = "D" // trim(chmain) // "d"
             point_g2 = 'Dnd'
          else
             point_group = "D" // trim(chmain)
             point_g2 = 'Dn'
          endif
       else
          if (sigmah) then
             point_group = "C" // trim(chmain) // "h"
             point_g2 = 'Cnh'
          else if (mainorder .eq. nsigma) then
             point_group = "C" // trim(chmain) // "v"
             point_g2 = 'Cnv'
          else if (2*mainorder .eq. mainimproper) then
             point_group = "S" // trim(chimp)
             point_g2 = 'S2n'
             opmainord = mainorder
             opmainord = mainimproper
          else
             point_group = "C" // trim(chmain)
             point_g2 = 'Cn'
          endif
       endif
    else if (nopsym .eq. 2) then
       if (nsigma .eq. 1) then
          point_group = 'Cs'
          point_g2 = 'Cs'
       else if (inversion) then
          point_group = 'Ci'
          point_g2 = 'Ci'
       else
          call ferror ('symgetgroup', 'Unknown group', warning)
          if (verbose) &
             write (uout,700) mainorder, mainmult, inversion, binmult&
             , nbinperp, nsigma, sigmah, mainimproper
          point_group = '??'
       endif
    else if (nopsym .eq. 1) then
       point_group = 'C1'
       point_g2 = 'C1'
    else
       call ferror ('symgetgroup', 'Unknown group', warning)
       if (verbose) &
          write (uout,700) mainorder, mainmult, inversion, binmult&
          , nbinperp, nsigma, sigmah, mainimproper
       point_group = '??'
    endif

700 format (/&
       1x, '++SYMGETGROUP:'/&
       1x, 'DBG(symgetgroup) main order & mult.......:', 2i5/&
       1x, 'DBG(symgetgroup) inversion...............:', l5/&
       1x, 'DBG(symgetgroup) binary axes.............:', i5/&
       1x, 'DBG(symgetgroup) binary axes perp to main:', i5/&
       1x, 'DBG(symgetgroup) symmetry planes.........:', i5/&
       1x, 'DBG(symgetgroup) sigma_h plane...........:', l5/&
       1x, 'DBG(symgetgroup) main improper axis......:', i5)
    !
  end subroutine symgetgroup

  !> Get the (z,x,z) counterclockwise Euler angles that
  !> correspond to the rotation input "matrix".
  !> Adapted from TESSEL.
  subroutine symeuler (matrix, euler)
    use param

    real*8 :: matrix(3,3), euler(3)

    if (abs(matrix(3,3)) .gt. 1d0) then
       euler(2) = acos(sign(1d0,matrix(3,3)))
    else
       euler(2) = acos(matrix(3,3))
    endif
    if (abs(abs(matrix(3,3))-1d0) .le. 1d-6) then
       euler(3) = 0d0
       if (matrix(3,3) .gt. 0d0) then
          euler(1) = atan2(matrix(2,1), matrix(1,1))
       else
          euler(1) = atan2(-matrix(2,1), matrix(1,1))
       endif
    else
       euler(1) = atan2(matrix(3,1), matrix(3,2))
       euler(3) = atan2(matrix(1,3), -matrix(2,3))
    endif
    if (euler(1) .lt. 0d0) euler(1) = euler(1) + pi + pi
    if (euler(3) .lt. 0d0) euler(3) = euler(3) + pi + pi

  end subroutine symeuler

  !> Obtain the eigenvalues and eigenvectors of the
  !> orthogonal matrix xmat(,) by calling to the dgeev() routine in
  !> the lapack package. If xdet is +/-1, the routine uses this
  !> value as the determinant of xmat(,). Otherwise the determinant
  !> is obtained.
  !>
  !> The eigenvalues are sorted with xval(3) being the equivalent
  !> to the matrix determinant and xvect(:,3) being the normal vector
  !> corresponding to the rotation axis.
  !> The two other eigenvalues are meaningless (the imaginary part is
  !> not returned) but their eigenvectors define the plane normal to
  !> the rotation axis.
  !> The three eigenvectors are real and orthonormal.
  !> Adapted from TESSEL.
  subroutine symeigen (xmat, xval, xvect, xdet, ierror)

    integer           ierror
    real*8            xmat(3,3), xval(3), xvect(3,3), xdet
    !
    integer           i, j, ii, ij, ik
    real*8            ar(3,3), wr(3), wi(3), vr(3,3), vl(3,3)
    real*8            dif, difmax, xnorm
    integer, parameter :: lwork = 102
    real*8            work(lwork)

    !.Get the determinant if it is needed:
    !
    if (abs(abs(xdet)-1d0) .ge. TOLeigen) then
       xdet = xmat(1,1)*(xmat(2,2)*xmat(3,3) - xmat(2,3)*xmat(3,2))&
          + xmat(1,2)*(xmat(2,3)*xmat(3,1) - xmat(2,1)*xmat(3,3))&
          + xmat(1,3)*(xmat(2,1)*xmat(3,2) - xmat(2,2)*xmat(3,1))
    endif

    do i = 1, 3
       do j = 1, 3
          ar(i,j) = xmat(i,j)
       enddo
    enddo
    call dgeev ("N", "V", 3, ar, 3, wr, wi, vl, 3, vr, 3,&
       work, lwork, ierror)

    !.Get the eigenvalue identical to the determinant and return it as
    ! the third one:
    ik = 0
    difmax = 1d40
    do i = 1, 3
       dif = abs(wr(i)-xdet) + abs(wi(i))
       if (dif .lt. difmax) then
          ik = i
          difmax = dif
       endif
    enddo
    ii = mod(ik,3) + 1
    ij = mod(ii,3) + 1
    ERReigen = max(ERReigen, difmax)
    xval(3) = wr(ik)
    xvect(1,3) = vr(1,ik)
    xvect(2,3) = vr(2,ik)
    xvect(3,3) = vr(3,ik)
    xval(1) = wr(ii)
    xvect(1,1) = vr(1,ii)
    xvect(2,1) = vr(2,ii)
    xvect(3,1) = vr(3,ii)
    xval(2) = wr(ij)

    !.Get explicitely the v2 = v3 x v1 eigenvector:
    !
    xvect(1,2) = xvect(2,3) * xvect(3,1) - xvect(3,3) * xvect(2,1)
    xvect(2,2) = xvect(3,3) * xvect(1,1) - xvect(1,3) * xvect(3,1)
    xvect(3,2) = xvect(1,3) * xvect(2,1) - xvect(2,3) * xvect(1,1)

    !.Enforce the normalization:
    !
    do i = 1, 3
       xnorm = xvect(1,i)**2 + xvect(2,i)**2 + xvect(3,i)**2
       xnorm = 1d0 / sqrt(xnorm)
       xvect(1,i) = xvect(1,i) * xnorm
       xvect(2,i) = xvect(2,i) * xnorm
       xvect(3,i) = xvect(3,i) * xnorm
    enddo

  end subroutine symeigen

  !> Get the number of xmat in the symmetry operator list.
  !> A -1 value will be returned if xmat is not in the list.
  !> Adapted from TESSEL.
  function symnumber (xmat)

    integer :: symnumber

    real*8  :: xmat(3,3)

    integer :: i, j, k
    real*8  :: diff

    symnumber = -1
    do i = 1, nopsym
       diff = 0d0
       do j = 1, 3
          do k = 1, 3
             diff = diff + abs(xmat(j,k) - opsym(i,j,k))
          enddo
       enddo
       if (diff .le. TOLeqvm) then
          symnumber = i
          return
       endif
    enddo

  end function symnumber

end module sympg
