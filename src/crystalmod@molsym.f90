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
submodule (crystalmod) molsym
  use param, only: pi
  implicit none

  type orbit
     real*8 :: dist ! distance to the center of mass
     integer :: z ! atomic number
     integer :: nat ! number of atoms in the orbit
     integer, allocatable :: id(:) ! ids of the atoms in the orbit
  end type orbit

  real*8, parameter :: eps_linear = 1d-5 ! check whether a molecule is linear
  real*8, parameter :: eps_planar = 1d-5 ! check whether a molecule is planar
  real*8, parameter :: eps_orbit = 1d-5 ! check to determine whether an atom is in an orbit
  real*8, parameter :: eps_equalop = 1d-5 ! check whether two symmetry operations are equal
  real*8, parameter :: eps_pos = 1d-5 ! whether to atomic positions are equal
  real*8, parameter :: eps_planar_axis = 1d-5 ! eps to find the axis normal to a plane
  real*8, parameter :: eps_planar_det3 = 1d-5 ! eps to find whether triplet is planar
  real*8, parameter :: eps_isinteger = 1d-5 ! eps for whether a number is an integer
  real*8, parameter :: eps_zeroangle = 1d0 * pi / 180d0  ! eps for zero angle (radians)

  integer, parameter :: norder_inf = 4 ! order of the infinite group (Cinfv, Dinfh)

  ! type for molecular symmetry operations
  type molsymop
     real*8 :: m(3,3) ! the rotation matrix
     integer :: type ! type of operation (rotation, plane,...)
     integer :: opn ! the n in C_n^m or S_n^m
     integer :: opm ! the m in C_n^m or S_n^m
     logical :: proper ! whether this is a proper symmetry operation
     real*8 :: axis(3) ! axis for rotations/normal for planes
     character(len=:), allocatable :: sym ! symbol
  end type molsymop

  ! identifiers for molecular symmetry operations
  integer, parameter :: molsymop_identity = 1
  integer, parameter :: molsymop_inversion = 2
  integer, parameter :: molsymop_rotation = 3
  integer, parameter :: molsymop_plane = 4
  integer, parameter :: molsymop_imp_rotation = 5
  integer, parameter :: molsymop_unknown = 6

  !xx! private procedures

contains

  ! Calculate the symmetry operations and point group of a molecule
  module subroutine calcmolsym(c,errmsg)
    use tools_math, only: cross
    use param, only: atmass
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, nat, iref, norbit
    real*8, allocatable :: x(:,:), w(:)
    integer, allocatable :: z(:), iorb(:)
    real*8 :: xcm(3), xref(3), xthis(3)
    logical :: islinear, isplanar
    type(orbit), allocatable :: orb(:)
    integer :: nrotm
    type(molsymop), allocatable :: mrotm(:)
    character(len=:), allocatable :: point_group

    ! initialize
    nrotm = 0
    allocate(mrotm(10))
    errmsg = ""

    ! allocate atomic positions and parameters
    nat = c%ncel
    allocate(x(3,nat),w(nat),z(nat))
    do i = 1, nat
       x(:,i) = c%atcel(i)%r
       z(i) = c%spc(c%atcel(i)%is)%z
       w(i) = atmass(z(i))
    end do

    ! calculate the center of mass
    xcm = 0d0
    do i = 1, nat
       xcm = xcm + w(i) * x(:,i)
    end do
    xcm = xcm / sum(w)

    ! center the molecule
    do i = 1, nat
       x(:,i) = x(:,i) - xcm
    end do

    ! whether the molecule is linear
    islinear = .true.
    iref = -1
    if (nat > 2) then
       xref = x(:,2) - x(:,1)
       xref = xref / norm2(xref)
       do i = 3, nat
          xthis = x(:,i) - x(:,1)
          xthis = xthis / norm2(xthis)
          if (abs(abs(dot_product(xthis,xref)) - 1d0) > eps_linear) then
             iref = i
             islinear = .false.
             exit
          end if
       end do
    end if

    ! whether the molecule is planar
    isplanar = .true.
    if (.not.islinear .and. nat > 3 .and. iref > 0) then
       ! vector normal to the plane
       xref = cross(x(:,2) - x(:,1), x(:,iref) - x(:,1))
       xref = xref / norm2(xref)

       ! check if all other vectors are in the same plane
       do i = 3, nat
          if (i == iref) cycle
          xthis = x(:,i) - x(:,1)
          xthis = xthis / norm2(xthis)
          if (abs(dot_product(xthis,xref)) > eps_planar) then
             isplanar = .false.
             exit
          end if
       end do
    end if

    ! calculate the orbits for this system
    call calculate_orbits(nat,x,z,iorb,norbit,orb)

    ! calculate the symmetry operations
    if (islinear) then
       call calcrotm_linear(nat,x,z,norbit,orb,nrotm,mrotm)
    elseif (isplanar) then
       call calcrotm_planar(nat,x,z,iorb,norbit,orb,nrotm,mrotm)
    else
       call calcrotm_general(nat,x,z,iorb,norbit,orb,nrotm,mrotm)
    end if

    ! get the point group name
    call get_group_name(nrotm,mrotm,point_group)

    ! ! print the point group information
    ! call print_group(point_group,nrotm,mrotm)

  end subroutine calcmolsym

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
    type(orbit), allocatable :: orbaux(:)

    ! initialize
    norbit = 0
    if (allocated(iorb)) deallocate(iorb)
    if (allocated(orb)) deallocate(orb)
    allocate(iorb(nat),orb(10))

    main: do i = 1, nat
       dist = norm2(x(:,i))

       ! try to find if a current orbit fits
       do j = 1, norbit
          if (abs(orb(j)%dist - dist) < eps_orbit .and. orb(j)%z == z(i)) then
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
       if (norbit > size(orb,1)) then
          allocate(orbaux(2*norbit))
          orbaux(1:size(orb,1)) = orb
          call move_alloc(orbaux,orb)
       end if
       orb(norbit)%nat = 1
       orb(norbit)%dist = dist
       orb(norbit)%z = z(i)
       allocate(orb(norbit)%id(10))
       orb(norbit)%id(1) = i
       iorb(i) = norbit
    end do main

    ! reallocate
    if (norbit /= size(orb,1)) then
       allocate(orbaux(norbit))
       orbaux = orb(norbit)
       call move_alloc(orbaux,orb)
    end if
    do i = 1, norbit
       call realloc(orb(i)%id,orb(i)%nat)
    end do

  end subroutine calculate_orbits

  ! Calculate rotation matrices of a linear molecule. The molecule is
  ! assumed to be centered at the origin. Input is the number of
  ! atoms, coordinates, atomic numbers, and orbits. Outputs
  ! rotation matrices.
  subroutine calcrotm_linear(nat,x,z,norbit,orb,nrotm,mrotm)
    use tools_math, only: cross
    use param, only: eye, tpi
    integer, intent(in) :: nat
    real*8, intent(in) :: x(3,nat)
    integer, intent(in) :: z(nat)
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
    mat = rotation_matrix(axis,tpi / norder_inf)
    call add_symop(mat,nrotm,mrotm)

    ! add a reflection plane with normal perpendicular to the axis
    xref = (/0d0,0d0,1d0/)
    if (abs(abs(dot_product(axis,xref)) - 1d0) < eps_linear) &
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
  subroutine calcrotm_planar(nat,x,z,iorb,norbit,orb,nrotm,mrotm)
    use tools_math, only: cross, matinv, det3
    use param, only: eye
    integer, intent(in) :: nat
    real*8, intent(in) :: x(3,nat)
    integer, intent(in) :: z(nat)
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
       if (nn > eps_planar_axis) then
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
          if (abs(dd) < eps_planar_det3) cycle
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
  subroutine calcrotm_general(nat,x,z,iorb,norbit,orb,nrotm,mrotm)
    use tools_math, only: matinv, det3
    use param, only: eye
    integer, intent(in) :: nat
    real*8, intent(in) :: x(3,nat)
    integer, intent(in) :: z(nat)
    integer, intent(in) :: iorb(nat)
    integer, intent(in) :: norbit
    type(orbit), intent(in) :: orb(norbit)
    integer, intent(inout) :: nrotm
    type(molsymop), intent(inout), allocatable :: mrotm(:)

    integer :: i1, i2, i3, oi1, oi2, oi3, j1, j2, j3
    real*8 :: mat(3,3), mi(3,3), mop(3,3), dd

    ! initialize the symmetry matrices with the identity matrix
    nrotm = 0
    call add_symop(eye,nrotm,mrotm)

    ! see if the inversion is an operation
    mat = -eye
    if (is_symop(mat,nat,x,norbit,orb)) &
       call add_symop(mat,nrotm,mrotm)

    ! find a triplet of atoms that is not coplanar with the origin
    oi1 = 0
    oi2 = 0
    oi3 = 0
    main: do i1 = 1, nat
       mi(:,1) = x(:,i1)
       do i2 = i1+1, nat
          mi(:,2) = x(:,i2)
          do i3 = i2+1, nat
             mi(:,3) = x(:,i3)
             dd = det3(mi)
             if (abs(dd) > eps_planar_det3) then
                oi1 = i1
                oi2 = i2
                oi3 = i3
                exit main
             end if
          end do
       end do
    end do main
    call matinv(mi,3)

    ! run over all other triplets with the same orbits
    do j1 = 1, nat
       if (iorb(j1) /= iorb(oi1)) cycle
       mat(:,1) = x(:,j1)
       do j2 = 1, nat
          if (iorb(j2) /= iorb(oi2)) cycle
          mat(:,2) = x(:,j2)
          do j3 = 1, nat
             if (iorb(j3) /= iorb(oi3)) cycle
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
    type(molsymop), allocatable :: aux(:)
    real*8 :: det, trace, angle, xx, xmin, xdif
    real*8 :: vp(3,3), eval(3), evali(3), axis(3)

    ! skip if the matrix is already known
    do i = 1, nrotm
       if (all(abs(mat - mrotm(i)%m) < eps_equalop)) then
          return
       end if
    end do

    ! reallocate if necessary and add the matrix
    nrotm = nrotm + 1
    if (nrotm > size(mrotm,1)) then
       allocate(aux(2*nrotm))
       aux(1:size(mrotm,1)) = mrotm
       call move_alloc(aux,mrotm)
    end if
    mrotm(nrotm)%m = mat

    ! trace and determinant of the operation
    trace = mat(1,1) + mat(2,2) + mat(3,3)
    det = det3(mat)

    ! whether the operation is proper
    if (abs(abs(det)-1d0) > eps_isinteger) then
       write (*,*) "error! (abs(abs(det)-1d0) > eps_isinteger)"
       stop 1
    elseif (det > 0d0) then
       mrotm(nrotm)%proper = .true.
    else
       mrotm(nrotm)%proper = .false.
    end if

    ! identity and inversion treated separately
    if (abs(abs(trace)-3d0) < eps_isinteger) then
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
       if (abs(angle) < eps_zeroangle) angle = tpi

       ! calculate the order of the rotation
       xmin = 1d40
       do i = 1, 20
          xx = i * tpi / angle
          xdif = abs(xx - nint(xx))
          if (xdif < eps_isinteger) then
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
       elseif (abs(trace - 1d0) < eps_isinteger) then
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
             if (norm2(x(:,orb(io)%id(j)) - xt) < eps_pos) then
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

  ! Returns the rotation matrix for a rotation around the given axis
  ! with the given angle in radians. Applies Rodrigues' rotation
  ! formula.
  function rotation_matrix(axis,angle) result(mat)
    use param, only: eye
    real*8, intent(in) :: axis(3), angle
    real*8 :: mat(3,3)

    real*8 :: nn, c, s, ax(3), temp(3)

    mat = eye
    nn = norm2(axis)
    if (nn < 1d-20) return

    c = cos(angle)
    s = sin(angle)
    ax = axis / nn
    temp = (1-c) * ax

    mat(1,1) = c + temp(1) * ax(1)
    mat(2,1) = temp(1) * ax(2) + s * ax(3)
    mat(3,1) = temp(1) * ax(3) - s * ax(2)

    mat(1,2) = temp(2) * ax(1) - s * ax(3)
    mat(2,2) = c + temp(2) * ax(2)
    mat(3,2) = temp(2) * ax(3) + s * ax(1)

    mat(1,3) = temp(3) * ax(1) + s * ax(2)
    mat(2,3) = temp(3) * ax(2) - s * ax(1)
    mat(3,3) = c + temp(3) * ax(3)

  end function rotation_matrix

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
          if (abs(abs(xprod)-1d0) < eps_isinteger) then
             sigmah = .true.
             mrotm(isigma(i))%sym = 'sigma_h'
          elseif (abs(xprod) < eps_isinteger) then
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
          if (abs(xprod) < eps_isinteger) &
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

  ! Print the point group operations to standard output.
  subroutine print_group(pg,nrotm,mrotm)
    use tools_io, only: uout, string, ioj_left
    character*(*), intent(in) :: pg
    integer, intent(in) :: nrotm
    type(molsymop), intent(in), allocatable :: mrotm(:)

    integer :: i, j

    write (uout,'("* Point group information")')
    write (uout,'("+ Point group: ",A)') trim(pg)
    write (uout,'("+ Number of symmetry operations: ",A)') string(nrotm)


    write (uout,'("#id  symbol                 axis/normal vector")')
    do i = 1, nrotm
       write (uout,'(99(A,X))') string(i,4,ioj_left),&
          string(mrotm(i)%sym,10,ioj_left),&
          (string(mrotm(i)%axis(j),'f',13,8,5),j=1,3)
    end do
    write (uout,*)

  end subroutine print_group

end submodule molsym
