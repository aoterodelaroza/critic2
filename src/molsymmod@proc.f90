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

! Routines for handling molecular symmetry and related calculations
submodule (molsymmod) proc
  use types, only: molsymop, molsymop_identity, molsymop_inversion, molsymop_rotation,&
     molsymop_plane, molsymop_imp_rotation, molsymop_unknown
  use param, only: pi
  implicit none

  type orbit
     real*8 :: dist ! distance to the center of mass
     integer :: z ! atomic number
     integer :: nat ! number of atoms in the orbit
     integer, allocatable :: id(:) ! ids of the atoms in the orbit
  end type orbit

  real*8, parameter :: eps_dist = 3d-5 ! distance/coordinate comparisons: linearity, planarity, orbits, atom images (tessel: TOLdist)
  real*8, parameter :: eps_eqvm = 2d0 * eps_dist ! check whether two symmetry operations are equal (tessel: TOLeqvm = 2*TOLdist)
  real*8, parameter :: eps_sng = 1d-5 ! singularity checks: normal axis, triplet determinant (tessel: TOLsng)
  real*8, parameter :: eps_isint = 3d-5 ! eps for whether a number is an integer (tessel: TOLisint)
  real*8, parameter :: eps_null = 3d-6 ! eps for zero angle, in radians (tessel: TOLnull)

  integer, parameter :: norder_inf = 4 ! order of the infinite group (Cinfv, Dinfh)

  !xx! private procedures

contains

  !> Clear the information in a point group class.
  module subroutine point_group_clear(p)
    class(point_group), intent(inout) :: p

    p%avail = .false.
    p%isatom = .false.
    p%islinear = .false.
    p%isplanar = .false.
    p%xcm = 0d0
    p%nop = 0
    if (allocated(p%op)) deallocate(p%op)
    p%symbol = "??"

  end subroutine point_group_clear

  !> Print the point group operations to standard output.
  module subroutine point_group_report(p)
    use tools_io, only: uout, string, ioj_left
    class(point_group), intent(inout) :: p

    integer :: i, j

    write (uout,'("+ Molecular point group: ",A)') trim(p%symbol)
    write (uout,'("+ List of symmetry operations (",A,")")') string(p%nop)

    write (uout,'("#id  symbol                 axis/normal vector")')
    do i = 1, p%nop
       write (uout,'(99(A,X))') string(i,4,ioj_left),&
          string(p%op(i)%sym,10,ioj_left),&
          (string(p%op(i)%axis(j),'f',13,8,5),j=1,3)
    end do
    write (uout,*)

  end subroutine point_group_report

  !> Initialize the point group as C1
  module subroutine point_group_init_as_c1(p)
    use param, only: eye
    class(point_group), intent(inout) :: p

    call p%clear()
    p%nop = 0
    allocate(p%op(1))
    call add_symop(eye,p%nop,p%op)
    p%symbol = 'C1'
    p%avail = .true.

  end subroutine point_group_init_as_c1

  ! Calculate the symmetry operations and point group of a molecule given
  ! by nat atoms with Cartesian coordinates xin(3,nat) and atomic numbers
  ! z(nat). The result (operations and point-group symbol) is returned in
  ! pg; errmsg is non-empty on error.
  module subroutine calc_point_group(nat,xin,z,pg,errmsg)
    use tools_math, only: cross
    use param, only: atmass
    integer, intent(in) :: nat
    real*8, intent(in) :: xin(3,nat)
    integer, intent(in) :: z(nat)
    type(point_group), intent(inout) :: pg
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, iref, norbit
    real*8, allocatable :: x(:,:), w(:)
    integer, allocatable :: iorb(:)
    real*8 :: xref(3), xthis(3)
    type(orbit), allocatable :: orb(:)

    ! initialize; return if there are no atoms
    errmsg = ""
    call pg%clear()
    if (nat == 0) return
    pg%nop = 0
    allocate(pg%op(10))

    ! local (mutable) copy of the positions and atomic weights
    allocate(x(3,nat),w(nat))
    do i = 1, nat
       x(:,i) = xin(:,i)
       w(i) = atmass(z(i))
    end do

    ! calculate the center of mass
    pg%xcm = 0d0
    do i = 1, nat
       pg%xcm = pg%xcm + w(i) * x(:,i)
    end do
    pg%xcm = pg%xcm / sum(w)

    ! center the molecule
    do i = 1, nat
       x(:,i) = x(:,i) - pg%xcm
    end do

    ! whether the molecule is an atom
    pg%isatom = (nat == 1)

    ! whether the molecule is linear
    pg%islinear = .true.
    iref = -1
    if (nat > 2) then
       xref = x(:,2) - x(:,1)
       xref = xref / norm2(xref)
       do i = 3, nat
          xthis = x(:,i) - x(:,1)
          xthis = xthis / norm2(xthis)
          if (abs(abs(dot_product(xthis,xref)) - 1d0) > eps_dist) then
             iref = i
             pg%islinear = .false.
             exit
          end if
       end do
    end if

    ! whether the molecule is planar
    pg%isplanar = .true.
    if (.not.pg%islinear .and. nat > 3 .and. iref > 0) then
       ! vector normal to the plane
       xref = cross(x(:,2) - x(:,1), x(:,iref) - x(:,1))
       xref = xref / norm2(xref)

       ! check if all other vectors are in the same plane
       do i = 3, nat
          if (i == iref) cycle
          xthis = x(:,i) - x(:,1)
          xthis = xthis / norm2(xthis)
          if (abs(dot_product(xthis,xref)) > eps_dist) then
             pg%isplanar = .false.
             exit
          end if
       end do
    end if

    ! calculate the orbits for this system
    call calculate_orbits(nat,x,z,iorb,norbit,orb)

    ! calculate the symmetry operations
    if (pg%isatom) then
       call calcrotm_atom(pg%nop,pg%op)
    elseif (pg%islinear) then
       call calcrotm_linear(nat,x,norbit,orb,pg%nop,pg%op)
    elseif (pg%isplanar) then
       call calcrotm_planar(nat,x,iorb,norbit,orb,pg%nop,pg%op)
    else
       call calcrotm_general(nat,x,iorb,norbit,orb,pg%nop,pg%op)
    end if

    ! get the point group name
    if (pg%isatom) then
       pg%symbol = 'O3'
    elseif (pg%islinear) then
       pg%symbol = 'Cinfv'
       do i = 1, pg%nop
          if (pg%op(i)%type == molsymop_inversion) then
             pg%symbol = 'Dinfh'
             exit
          end if
       end do
    else
       call get_group_name(pg%nop,pg%op,pg%symbol)
    end if

    ! wrap up
    pg%avail = .true.

  end subroutine calc_point_group

  !xx! private procedures

  ! Calculate the orbits for the atoms in the molecule (nat atoms,
  ! x = coordinates, z = atomic numbers). Returns the orbit ID for each
  ! atom (iorb), the number of orbits (norbit), and the orbits themselves
  ! (orb).
  subroutine calculate_orbits(nat,x,z,iorb,norbit,orb)
    use types, only: realloc
    integer, intent(in) :: nat
    real*8, intent(in) :: x(3,nat)
    integer, intent(in) :: z(nat)
    integer, intent(inout), allocatable :: iorb(:)
    integer, intent(out) :: norbit
    type(orbit), allocatable, intent(inout) :: orb(:)

    integer :: i, j
    real*8 :: dist

    ! initialize
    norbit = 0
    if (allocated(iorb)) deallocate(iorb)
    if (allocated(orb)) deallocate(orb)
    allocate(iorb(nat),orb(10))

    main: do i = 1, nat
       dist = norm2(x(:,i))

       ! try to find if a current orbit fits
       do j = 1, norbit
          if (abs(orb(j)%dist - dist) < eps_dist .and. orb(j)%z == z(i)) then
             orb(j)%nat = orb(j)%nat + 1
             if (orb(j)%nat > size(orb(j)%id,1)) &
                call realloc(orb(j)%id,2*orb(j)%nat)
             orb(j)%id(orb(j)%nat) = i
             iorb(i) = j
             cycle main
          end if
       end do

       ! this must be a new orbit
       norbit = norbit + 1
       if (norbit > size(orb,1)) call realloc_orbit(2*norbit)
       orb(norbit)%nat = 1
       orb(norbit)%dist = dist
       orb(norbit)%z = z(i)
       allocate(orb(norbit)%id(10))
       orb(norbit)%id(1) = i
       iorb(i) = norbit
    end do main

    ! reallocate
    if (norbit /= size(orb,1)) call realloc_orbit(norbit)
    do i = 1, norbit
       call realloc(orb(i)%id,orb(i)%nat)
    end do

  contains
    subroutine realloc_orbit(n)
      integer, intent(in) :: n

      type(orbit), allocatable :: orbaux(:)

      allocate(orbaux(n))
      orbaux(1:min(n,size(orb,1))) = orb(1:min(n,size(orb,1)))
      call move_alloc(orbaux,orb)

    end subroutine realloc_orbit

  end subroutine calculate_orbits

  ! Calculate rotation matrices for an atom - just the identity and
  ! exit.
  subroutine calcrotm_atom(nrotm,mrotm)
    use tools_math, only: cross
    use param, only: eye
    integer, intent(inout) :: nrotm
    type(molsymop), intent(inout), allocatable :: mrotm(:)

    ! initialize the symmetry matrices with the identity matrix
    nrotm = 0
    call add_symop(eye,nrotm,mrotm)

  end subroutine calcrotm_atom

  ! Calculate rotation matrices of a linear molecule. The molecule is
  ! assumed to be centered at the origin. Input is the number of
  ! atoms, coordinates, atomic numbers, and orbits. Outputs
  ! rotation matrices.
  subroutine calcrotm_linear(nat,x,norbit,orb,nrotm,mrotm)
    use tools_math, only: cross, axisangle2mat
    use param, only: eye, tpi
    integer, intent(in) :: nat
    real*8, intent(in) :: x(3,nat)
    integer, intent(in) :: norbit
    type(orbit), intent(in) :: orb(norbit)
    integer, intent(inout) :: nrotm
    type(molsymop), intent(inout), allocatable :: mrotm(:)

    real*8 :: mat(3,3), xref(3), axis(3)

    ! initialize the symmetry matrices with the identity matrix
    nrotm = 0
    call add_symop(eye,nrotm,mrotm)

    ! see if the inversion is an operation
    mat = -eye
    if (is_symop(mat,nat,x,norbit,orb)) &
       call add_symop(mat,nrotm,mrotm)

    ! add the rotation along the axis
    axis = x(:,2) - x(:,1)
    axis = axis / norm2(axis)
    mat = axisangle2mat(axis,tpi / norder_inf)
    call add_symop(mat,nrotm,mrotm)

    ! add a reflection plane with normal perpendicular to the axis
    xref = (/0d0,0d0,1d0/)
    if (abs(abs(dot_product(axis,xref)) - 1d0) < eps_dist) &
       xref = (/0d0,1d0,0d0/)
    mat = reflection_matrix(cross(axis,xref))
    call add_symop(mat,nrotm,mrotm)

    ! complete the group
    call group_closure(nrotm,mrotm)

  end subroutine calcrotm_linear

  ! Calculate rotation matrices of a planar molecule. The molecule is
  ! assumed to be centered at the origin. Input is the number of
  ! atoms, coordinates, atomic numbers, and orbits. Outputs
  ! rotation matrices.
  subroutine calcrotm_planar(nat,x,iorb,norbit,orb,nrotm,mrotm)
    use tools_math, only: cross, matinv, det3
    use param, only: eye
    integer, intent(in) :: nat
    real*8, intent(in) :: x(3,nat)
    integer, intent(in) :: iorb(nat)
    integer, intent(in) :: norbit
    type(orbit), intent(in) :: orb(norbit)
    integer, intent(inout) :: nrotm
    type(molsymop), intent(inout), allocatable :: mrotm(:)

    integer :: i, i1, i2, j1, j2
    real*8 :: mat(3,3), mi(3,3), mop(3,3), axis(3), nn, dd

    ! initialize the symmetry matrices with the identity matrix
    nrotm = 0
    call add_symop(eye,nrotm,mrotm)

    ! see if the inversion is an operation
    mat = -eye
    if (is_symop(mat,nat,x,norbit,orb)) &
       call add_symop(mat,nrotm,mrotm)

    ! calculate the axis perpendicular to the molecular plane and add
    ! it as symop
    do i = 3, nat
       axis = cross(x(:,2) - x(:,1), x(:,i) - x(:,1))
       nn = norm2(axis)
       if (nn > eps_sng) then
          axis = axis / nn
          exit
       end if
    end do
    mat = reflection_matrix(axis)
    call add_symop(mat,nrotm,mrotm)

    ! run over unique pairs
    do i1 = 1, nat
       do i2 = i1+1, nat
          mi(:,1) = x(:,i1)
          mi(:,2) = x(:,i2)
          mi(:,3) = axis
          dd = det3(mi)
          if (abs(dd) < eps_sng) cycle
          call matinv(mi,3)

          ! run over pairs in the same orbits
          do j1 = 1, orb(iorb(i1))%nat
             do j2 = 1, orb(iorb(i2))%nat
                mat(:,1) = x(:,j1)
                mat(:,2) = x(:,j2)
                mat(:,3) = axis

                mop = matmul(mat,mi)
                if (is_symop(mop,nat,x,norbit,orb)) &
                   call add_symop(mop,nrotm,mrotm)
             end do
          end do
       end do
    end do

    ! complete the group
    call group_closure(nrotm,mrotm)

  end subroutine calcrotm_planar

  ! Calculate rotation matrices of a general molecule. The molecule is
  ! assumed to be centered at the origin. Input is the number of
  ! atoms, coordinates, atomic numbers, and orbits. Outputs rotation
  ! matrices.
  subroutine calcrotm_general(nat,x,iorb,norbit,orb,nrotm,mrotm)
    use tools_math, only: matinv, det3, cross
    use param, only: eye
    integer, intent(in) :: nat
    real*8, intent(in) :: x(3,nat)
    integer, intent(in) :: iorb(nat)
    integer, intent(in) :: norbit
    type(orbit), intent(in) :: orb(norbit)
    integer, intent(inout) :: nrotm
    type(molsymop), intent(inout), allocatable :: mrotm(:)

    integer :: i1, i2, i3, oi1, oi2, oi3, j1, j2, j3
    real*8 :: mat(3,3), mi(3,3), mop(3,3), dd, dbest

    ! initialize the symmetry matrices with the identity matrix
    nrotm = 0
    call add_symop(eye,nrotm,mrotm)

    ! see if the inversion is an operation
    mat = -eye
    if (is_symop(mat,nat,x,norbit,orb)) &
       call add_symop(mat,nrotm,mrotm)

    ! Choose a well-conditioned reference triplet (one not coplanar with
    ! the origin and as far from coplanar as possible). Using the first
    ! barely-acceptable triplet gives an ill-conditioned inverse and makes
    ! is_symop reject valid operations for highly symmetric/degenerate
    ! point sets (e.g. lattice-point clouds). Build it greedily in O(nat):
    ! the farthest atom, then the one most perpendicular to it, then the
    ! one maximizing the triple product.
    oi1 = 1
    do i1 = 2, nat
       if (norm2(x(:,i1)) > norm2(x(:,oi1))) oi1 = i1
    end do
    oi2 = 0
    dbest = -1d0
    do i2 = 1, nat
       if (i2 == oi1) cycle
       dd = norm2(cross(x(:,oi1),x(:,i2)))
       if (dd > dbest) then
          dbest = dd
          oi2 = i2
       end if
    end do
    oi3 = 0
    dbest = -1d0
    mi(:,1) = x(:,oi1)
    mi(:,2) = x(:,oi2)
    do i3 = 1, nat
       if (i3 == oi1 .or. i3 == oi2) cycle
       mi(:,3) = x(:,i3)
       dd = abs(det3(mi))
       if (dd > dbest) then
          dbest = dd
          oi3 = i3
       end if
    end do
    mi(:,1) = x(:,oi1)
    mi(:,2) = x(:,oi2)
    mi(:,3) = x(:,oi3)
    call matinv(mi,3)

    ! run over all other triplets with the same orbits
    do j1 = 1, nat
       if (iorb(j1) /= iorb(oi1)) cycle
       mat(:,1) = x(:,j1)
       do j2 = 1, nat
          if (iorb(j2) /= iorb(oi2)) cycle
          if (j2 == j1) cycle
          mat(:,2) = x(:,j2)
          do j3 = 1, nat
             if (iorb(j3) /= iorb(oi3)) cycle
             if (j3 == j1 .or. j3 == j2) cycle
             mat(:,3) = x(:,j3)

             mop = matmul(mat,mi)
             if (is_symop(mop,nat,x,norbit,orb)) &
                call add_symop(mop,nrotm,mrotm)
          end do
       end do
    end do

    ! complete the group
    call group_closure(nrotm,mrotm)

  end subroutine calcrotm_general

  ! Add symmety operation corresponding to mat to the list of rotation
  ! matrices (nrotm, mrotm), if it does not exist already.
  subroutine add_symop(mat,nrotm,mrotm)
    use tools_io, only: string
    use types, only: realloc
    use tools_math, only: det3, eig
    use param, only: tpi
    real*8, intent(in) :: mat(3,3)
    integer, intent(inout) :: nrotm
    type(molsymop), intent(inout), allocatable :: mrotm(:)

    integer :: i
    real*8 :: det, trace, angle, xx, xmin, xdif
    real*8 :: vp(3,3), eval(3), evali(3)

    ! trace and determinant of the operation
    trace = mat(1,1) + mat(2,2) + mat(3,3)
    det = det3(mat)

    ! skip if this is not a rotation
    if (abs(abs(det)-1d0) > eps_isint) return

    ! skip if the matrix is already known
    do i = 1, nrotm
       if (sum(abs(mat - mrotm(i)%m)) < eps_eqvm) then
          return
       end if
    end do

    ! reallocate if necessary and add the matrix
    nrotm = nrotm + 1
    if (nrotm > size(mrotm,1)) call realloc(mrotm,2*nrotm)
    mrotm(nrotm)%m = mat

    ! note whether the operation is proper/improper/wrong
    if (det > 0d0) then
       mrotm(nrotm)%proper = .true.
    else
       mrotm(nrotm)%proper = .false.
    end if

    ! identity and inversion treated separately
    if (abs(abs(trace)-3d0) < eps_isint) then
       if (trace > 0d0) then
          ! identity
          mrotm(nrotm)%type = molsymop_identity
          mrotm(nrotm)%opn = 1
          mrotm(nrotm)%opm = 1
          mrotm(nrotm)%axis = (/0d0,0d0,1d0/)
          mrotm(nrotm)%sym = "E"
       else
          ! inversion
          mrotm(nrotm)%type = molsymop_inversion
          mrotm(nrotm)%opn = 2
          mrotm(nrotm)%opm = 1
          mrotm(nrotm)%axis = (/0d0,0d0,1d0/)
          mrotm(nrotm)%sym = "i"
       end if
    else
       ! this is a C_n^m or S_n^m
       !! calculate the rotation angle
       angle = 0.5d0 * (trace - det)
       if (abs(angle) > 1d0) angle = sign(1d0, angle)
       angle = acos(angle)
       if (abs(angle) < eps_null) angle = tpi

       ! calculate the order of the rotation
       xmin = 1d40
       do i = 1, 20
          xx = i * tpi / angle
          xdif = abs(xx - nint(xx))
          if (xdif < eps_isint) then
             mrotm(nrotm)%opn = nint(xx)
             mrotm(nrotm)%opm = i
             exit
          elseif (xdif < xmin) then
             mrotm(nrotm)%opn = nint(xx)
             mrotm(nrotm)%opm = i
             xmin = xdif
          end if
       end do

       ! find the rotation axis
       vp = mat
       call eig(vp,3,eval,evali)
       xmin = 1d40
       do i = 1, 3
          xdif = abs(eval(i)-det) + abs(evali(i))
          if (xdif < xmin) then
             mrotm(nrotm)%axis = vp(:,i) / norm2(vp(:,i))
             xmin = xdif
          end if
       end do

       ! operation type
       if (det > 0d0) then
          mrotm(nrotm)%type = molsymop_rotation
          mrotm(nrotm)%sym = "C_" // string(mrotm(nrotm)%opn) // "^" // string(mrotm(nrotm)%opm)
       elseif (abs(trace - 1d0) < eps_isint) then
          mrotm(nrotm)%type = molsymop_plane
          mrotm(nrotm)%sym = "sigma"
       elseif (mrotm(nrotm)%opn == 2 .and. mrotm(nrotm)%opm == 1) then
          mrotm(nrotm)%type = molsymop_plane
          mrotm(nrotm)%sym = "sigma"
       else
          mrotm(nrotm)%type = molsymop_imp_rotation
          mrotm(nrotm)%sym = "S_" // string(mrotm(nrotm)%opn) // "^" // string(mrotm(nrotm)%opm)
       end if
    end if

  end subroutine add_symop

  ! Check whether the given matrix is a symmetry operation
  function is_symop(mat,nat,x,norbit,orb)
    real*8, intent(in) :: mat(3,3)
    integer, intent(in) :: nat
    real*8, intent(in) :: x(3,nat)
    integer, intent(in) :: norbit
    type(orbit), intent(in) :: orb(norbit)
    logical :: is_symop

    integer :: io, i, j
    real*8 :: x0(3), xt(3)
    logical :: found

    is_symop = .false.

    ! run over orbits
    do io = 1, norbit
       ! run over atoms
       do i = 1, orb(io)%nat
          ! calculate the transformed atom
          x0 = x(:,orb(io)%id(i))
          xt = matmul(mat,x0)

          ! check if the transformed atom is in the orbit
          found = .false.
          do j = 1, orb(io)%nat
             if (norm2(x(:,orb(io)%id(j)) - xt) < eps_dist) then
                found = .true.
                exit
             end if
          end do

          if (.not.found) return
       end do
    end do
    is_symop = .true.

  end function is_symop

  ! Find the missing operations to make the set of matrices a group
  ! under multiplication (closure). Updates the nrotm and mrotm
  ! arguments.
  subroutine group_closure(nrotm,mrotm)
    integer, intent(inout) :: nrotm
    type(molsymop), intent(inout), allocatable :: mrotm(:)

    logical :: again
    integer :: i, j, nhere, nignore
    real*8 :: mat(3,3)

    nignore = 0
    again = .true.
    do while (again)
       again = .false.
       nhere = nrotm
       do i = 2, nhere
          do j = 2, nhere
             if (i <= nignore .and. j <= nignore) cycle
             mat = matmul(mrotm(i)%m,mrotm(j)%m)
             call add_symop(mat,nrotm,mrotm)
          end do
       end do
       if (nrotm > nhere) then
          nignore = nhere
          again = .true.
       end if
    end do

  end subroutine group_closure

  ! Returns the matrix corresponding to a reflection given the vector
  ! normal to the plane. The plane passes through the origin.
  ! R = I - 2 * (n * nT) / (n * n)
  function reflection_matrix(normal) result(mat)
    implicit none
    real*8, intent(in)  :: normal(3)
    real*8 :: mat(3,3)

    real*8 :: nn(3)

    ! compute squared norm
    nn = normal / norm2(normal)

    ! reflection matrix
    mat(1,1) = 1.0d0 - 2.0d0*nn(1)*nn(1)
    mat(1,2) = -2.0d0*nn(1)*nn(2)
    mat(1,3) = -2.0d0*nn(1)*nn(3)

    mat(2,1) = -2.0d0*nn(2)*nn(1)
    mat(2,2) = 1.0d0 - 2.0d0*nn(2)*nn(2)
    mat(2,3) = -2.0d0*nn(2)*nn(3)

    mat(3,1) = -2.0d0*nn(3)*nn(1)
    mat(3,2) = -2.0d0*nn(3)*nn(2)
    mat(3,3) = 1.0d0 - 2.0d0*nn(3)*nn(3)

  end function reflection_matrix

  !> Calculate the point group name from the symmetry matrices. From
  !> tessel.  Input: the number and info about the symmetry
  !> operations. Returns the point group symbol and maybe changes some
  !> of the symbols (e.g. sigma_h).
  subroutine get_group_name(nrotm,mrotm,pg)
    use tools_io, only: string
    integer, intent(in) :: nrotm
    type(molsymop), intent(inout), allocatable :: mrotm(:)
    character(len=:), allocatable, intent(inout) :: pg

    integer :: mainorder, mainmult, binmult, nsigma, mainimproper, nbinperp
    integer, allocatable :: imain(:), ibin(:), isigma(:)
    logical :: inversion, sigmah
    integer :: i
    real*8 :: aref(3), axis(3), xprod

    ! initialize
    pg = ""
    mainorder = 0
    mainmult = 0
    mainimproper = 0
    binmult = 0
    nsigma = 0
    nbinperp = 0
    inversion = .false.
    allocate(imain(nrotm),ibin(nrotm),isigma(nrotm))

    ! make a note of the important parameters
    do i = 1, nrotm
       if (mrotm(i)%type == molsymop_rotation .and. mrotm(i)%opm == 1) then
          ! main order
          if (mrotm(i)%opn > mainorder) then
             mainorder = mrotm(i)%opn
             mainmult = 1
             imain(mainmult) = i
          elseif (mrotm(i)%opn == mainorder) then
             mainmult = mainmult + 1
             imain(mainmult) = i
          end if

          ! 2-fold rotations
          if (mrotm(i)%opn == 2) then
             binmult = binmult + 1
             ibin(binmult) = i
          end if
       elseif (mrotm(i)%type == molsymop_inversion) then
          inversion = .true.
       elseif (mrotm(i)%type == molsymop_plane) then
          nsigma = nsigma + 1
          isigma(nsigma) = i
       elseif (mrotm(i)%type == molsymop_imp_rotation) then
          mainimproper = max(mainimproper,mrotm(i)%opn)
       end if
    end do

    ! the opm = 1 cases are doubled if the order is greater than 2
    if (mainorder > 2) mainmult = mainmult/2

    ! if there is a single main order group, look for sigma_h planes:
    sigmah = .false.
    if (mainmult == 1 .or. mainorder == 2) then
       aref = mrotm(imain(1))%axis
       do i = 1, nsigma
          axis = mrotm(isigma(i))%axis
          xprod = dot_product(axis,aref)
          if (abs(abs(xprod)-1d0) < eps_isint) then
             sigmah = .true.
             mrotm(isigma(i))%sym = 'sigma_h'
          elseif (abs(xprod) < eps_isint) then
             mrotm(isigma(i))%sym = 'sigma_v'
          end if
       end do
    endif

    ! if there is a single main order group, look for n C2 perpendicular axis
    if (mainmult == 1 .or. mainorder == 2) then
       aref = mrotm(imain(1))%axis
       do i = 1, binmult
          axis = mrotm(ibin(i))%axis
          xprod = dot_product(axis,aref)
          if (abs(xprod) < eps_isint) &
             nbinperp = nbinperp + 1
       end do
    endif

    ! decision tree to get the point group
    if (mainmult > 1 .and. mainorder > 2) then
       ! cubic groups
       if (mainorder == 5 .and. mainmult >= 6) then
          if (inversion) then
             pg = 'Ih'
          else
             pg = 'I'
          endif
       else if (mainorder == 4 .and. mainmult >= 3) then
          if (inversion) then
             pg = 'Oh'
          else
             pg = 'O'
          endif
       else if (mainorder == 3 .and. mainmult >= 4) then
          if (inversion) then
             pg = 'Th'
          else if (nsigma == 6) then
             pg = 'Td'
          else
             pg = 'T'
          endif
       else
          pg = '??'
       endif

    else if (mainorder >= 2) then
       if (mainorder == nbinperp) then
          if (sigmah) then
             pg = "D" // string(mainorder) // "h"
          else if (mainorder .eq. nsigma) then
             pg = "D" // string(mainorder) // "d"
          else
             pg = "D" // string(mainorder)
          endif
       else
          if (sigmah) then
             pg = "C" // string(mainorder) // "h"
          else if (mainorder .eq. nsigma) then
             pg = "C" // string(mainorder) // "v"
          else if (2*mainorder .eq. mainimproper) then
             pg = "S" // string(mainimproper)
          else
             pg = "C" // string(mainorder)
          endif
       end if
    else if (nrotm == 2) then
       if (nsigma .eq. 1) then
          pg = 'Cs'
       else if (inversion) then
          pg = 'Ci'
       else
          pg = '??'
       endif
    else if (nrotm == 1) then
       pg = 'C1'
    else
       pg = '??'
    end if

  end subroutine get_group_name

end submodule proc
