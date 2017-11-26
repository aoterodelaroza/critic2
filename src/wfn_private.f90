!> Part of the following code has been adapted from postg
!> Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
!> Felix Kannemann <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
!> Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider <hs7@post.queensu.ca>,
!> and Axel D. Becke <axel.becke@dal.ca>

! Copyright (c) 2009-2017 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! molecular wavefunction readers and tools
module wfn_private
  implicit none
  
  private

  ! Molecular wavefunction type
  ! Order of the orbitals:
  ! Restricted wavefunctions (wfntyp = wfn_rhf):
  !
  !     alpha and beta           alpha and beta
  ! | ..................... | ..................... |
  ! 1                     nmoocc                  nmoall
  !
  ! Unrestricted wavefunctions (wfntyp = wfn_uhf):
  !
  !     alpha       beta        alpha       beta
  ! | ......... | ......... | ......... | ......... |
  ! 1         nalpha      nmoocc     nmoocc+      nmoall
  !                                nalpha_virt
  ! 
  ! Fractional wavefunctions (wfntyp = wfn_frac): 
  ! 
  ! | ............................................. |
  ! 1                                         nmoall = nmoocc
  ! 
  ! In rhf and uhf wavefunctions, the virtual orbitals (nmoocc to
  ! nmoall) are only available if hasvirtual = .true.
  type molwfn
     logical :: useecp !< this wavefunction comes from a calc with ECPs
     integer :: nmoall !< number of molecular orbitals (occupied and virtual)
     integer :: nmoocc !< number of occupied molecular orbitals
     integer :: nalpha !< number of alpha electrons
     integer :: nalpha_virt !< number of virtual alpha orbitals
     integer :: npri !< number of primitives
     integer :: wfntyp !< type of wavefunction (rhf, uhf, fractional occ)
     logical :: issto !< are the primitives GTOs or STOs?
     logical :: hasvirtual !< are the virtual orbitals known?
     integer :: ixmaxsto(4) !< maximum exponent for x, y, z, and r in STOs
     integer, allocatable :: icenter(:) !< primitive center
     integer, allocatable :: itype(:) !< primitive type (see li(:,:) array)
     real*8, allocatable :: d2ran(:) !< maximum d^2 (GTO) or d (STO) to discard the primitive
     real*8, allocatable :: e(:) !< primitive exponents
     real*8, allocatable :: occ(:) !< MO occupation numbers
     real*8, allocatable :: cmo(:,:) !< MO coefficients
     integer :: nedf !< number of EDFs (electron density functions - core density for ECPs)
     integer, allocatable :: icenter_edf(:) !< EDF centers
     integer, allocatable :: itype_edf(:) !< EDF types
     real*8, allocatable :: e_edf(:) !< EDF exponents
     real*8, allocatable :: c_edf(:) !< EDF coefficients
     ! atomic positions, internal copy
     integer :: nat !< number of atoms
     real*8, allocatable :: xat(:,:) !< atomic coordinates (Cartesian, bohr)
   contains
     procedure :: end => wfn_end !< deallocate all arrays in wfn object
     procedure :: read_wfn !< read wavefunction from a Gaussian wfn file
     procedure :: read_wfx !< read wavefunction from a Gaussian wfx file
     procedure :: read_fchk !< read wavefunction from a Gaussian formatted checkpoint file
     procedure :: read_molden !< read wavefunction from a molden file
     procedure :: register_struct !< pass the atomic number and positions to the object
     procedure :: rho2 !< calculate the density, derivatives, and other properties
     procedure :: calculate_mo !< calculate the MO values at a point (driver)
     procedure :: calculate_mo_sto !< calculate the MO values at a point (STO version)
     procedure :: calculate_mo_gto !< calculate the MO values at a point (GTO version)
  end type molwfn
  public :: molwfn

  public :: wfn_read_xyz_geometry
  public :: wfn_read_wfn_geometry
  public :: wfn_read_wfx_geometry
  public :: wfn_read_fchk_geometry
  public :: wfn_read_molden_geometry
  private :: gnorm
  public :: wfx_read_integers
  public :: wfx_read_reals1

  ! wfn type identifier
  integer, parameter, public :: wfn_rhf = 0
  integer, parameter, public :: wfn_uhf = 1
  integer, parameter, public :: wfn_frac = 2

  ! double factorials minus one
  integer, parameter :: dfacm1(0:8) = (/1,1,1,2,3,8,15,48,105/) 

  ! double factorials
  integer, parameter :: dfac(0:8) = (/1,1,2,3,8,15,48,105,384/) 

  ! number of angular components for shell type
  !                                    gs fs ds ps  s  p  d   f   g
  integer, parameter :: nshlt(-4:4) = (/9, 7, 5, 3, 1, 3, 6, 10, 15/) 

  ! initial and final types for cartesian shells
  integer, parameter :: jshl0(0:4) = (/1, 2, 5,  11, 21/) ! s, p, d, f, g
  integer, parameter :: jshl1(0:4) = (/1, 4, 10, 20, 35/) ! s, p, d, f, g

  ! named constants for the sph to car matrices below
  real*8, parameter :: s3 = sqrt(3d0)
  real*8, parameter :: s3_4 = sqrt(3d0/4d0)
  real*8, parameter :: s3_8 = sqrt(3d0/8d0)
  real*8, parameter :: s5_8 = sqrt(5d0/8d0)
  real*8, parameter :: s5_16 = sqrt(5d0/16d0)
  real*8, parameter :: s6 = sqrt(6d0)
  real*8, parameter :: s10 = sqrt(10d0)
  real*8, parameter :: s10_8 = sqrt(10d0/8d0)
  real*8, parameter :: s15 = sqrt(15d0)
  real*8, parameter :: s15_4 = sqrt(15d0/4d0)
  real*8, parameter :: s35_4 = sqrt(35d0/4d0)
  real*8, parameter :: s35_8 = sqrt(35d0/8d0)
  real*8, parameter :: s35_64 = sqrt(35d0/64d0)
  real*8, parameter :: s45 = sqrt(45d0)
  real*8, parameter :: s45_4 = sqrt(45d0/4d0)
  real*8, parameter :: s45_8 = sqrt(45d0/8d0)
  real*8, parameter :: s315_8 = sqrt(315d0/8d0)
  real*8, parameter :: s315_16 = sqrt(315d0/16d0)
  real*8, parameter :: d32 = 3d0/2d0
  real*8, parameter :: d34 = 3d0/4d0
  real*8, parameter :: d38 = 3d0/8d0

  ! -- Real solid harmonics r^l * Slm as a function of Cartesian products. -- 

  ! dsphcar: l = 2 
  !   spherical molden order: m = 0, 1, -1, 2, -2
  !   Cartesian molden order: xx, yy, zz, xy, xz, yz
  real*8 :: dsphcar(5,6) = reshape((/&
     !  0      1     -1      2     -2
     -0.5d0, 0.0d0, 0.0d0,  s3_4, 0.0d0,& ! xx
     -0.5d0, 0.0d0, 0.0d0, -s3_4, 0.0d0,& ! yy
      1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,& ! zz
      0.0d0, 0.0d0, 0.0d0, 0.0d0,    s3,& ! xy
      0.0d0,    s3, 0.0d0, 0.0d0, 0.0d0,& ! xz
      0.0d0, 0.0d0,    s3, 0.0d0, 0.0d0 & ! yz
     /),shape(dsphcar))

  ! fsphcar: l = 3 
  !   spherical molden order: m = 0, 1, -1, 2, -2, 3, -3
  !   Cartesian molden order: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
  real*8 :: fsphcar(7,10) = reshape((/&
     ! 0        1       -1         2       -2         3       -3 
     0.0d0,   -s3_8,   0.0d0,    0.0d0,   0.0d0,     s5_8,   0.0d0,& ! xxx 
     0.0d0,   0.0d0,   -s3_8,    0.0d0,   0.0d0,    0.0d0,   -s5_8,& ! yyy 
     1.0d0,   0.0d0,   0.0d0,    0.0d0,   0.0d0,    0.0d0,   0.0d0,& ! zzz 
     0.0d0,   -s3_8,   0.0d0,    0.0d0,   0.0d0,   -s45_8,   0.0d0,& ! xyy 
     0.0d0,   0.0d0,   -s3_8,    0.0d0,   0.0d0,    0.0d0,   s45_8,& ! xxy 
      -d32,   0.0d0,   0.0d0,    s15_4,   0.0d0,    0.0d0,   0.0d0,& ! xxz 
     0.0d0,      s6,   0.0d0,    0.0d0,   0.0d0,    0.0d0,   0.0d0,& ! xzz 
     0.0d0,   0.0d0,      s6,    0.0d0,   0.0d0,    0.0d0,   0.0d0,& ! yzz 
      -d32,   0.0d0,   0.0d0,   -s15_4,   0.0d0,    0.0d0,   0.0d0,& ! yyz 
     0.0d0,   0.0d0,   0.0d0,    0.0d0,     s15,    0.0d0,   0.0d0 & ! xyz 
     /),shape(fsphcar))

  ! gsphcar: l = 4
  !   spherical molden order: m = 0, 1, -1, 2, -2, 3, -3, 4, -4
  !   Cartesian molden order: xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz
  real*8 :: gsphcar(9,15) = reshape((/&
  !    0       1      -1       2      -2        3      -3         4      -4
       d38,  0.0d0,  0.0d0, -s5_16,  0.0d0,   0.0d0,  0.0d0,   s35_64,  0.0d0,& ! xxxx
       d38,  0.0d0,  0.0d0,  s5_16,  0.0d0,   0.0d0,  0.0d0,   s35_64,  0.0d0,& ! yyyy
       1d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! zzzz
     0.0d0,  0.0d0,  0.0d0,  0.0d0, -s10_8,   0.0d0,  0.0d0,    0.0d0,  s35_4,& ! xxxy
     0.0d0, -s45_8,  0.0d0,  0.0d0,  0.0d0,   s35_8,  0.0d0,    0.0d0,  0.0d0,& ! xxxz
     0.0d0,  0.0d0,  0.0d0,  0.0d0, -s10_8,   0.0d0,  0.0d0,    0.0d0, -s35_4,& ! xyyy
     0.0d0,  0.0d0, -s45_8,  0.0d0,  0.0d0,   0.0d0, -s35_8,    0.0d0,  0.0d0,& ! yyyz
     0.0d0,    s10,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xzzz
     0.0d0,  0.0d0,    s10,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! yzzz
       d34,  0.0d0,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0, -s315_16,  0.0d0,& ! xxyy
      -3d0,  0.0d0,  0.0d0,  s45_4,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xxzz
      -3d0,  0.0d0,  0.0d0, -s45_4,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! yyzz
     0.0d0,  0.0d0, -s45_8,  0.0d0,  0.0d0,   0.0d0, s315_8,    0.0d0,  0.0d0,& ! xxyz
     0.0d0, -s45_8,  0.0d0,  0.0d0,  0.0d0, -s315_8,  0.0d0,    0.0d0,  0.0d0,& ! xyyz
     0.0d0,  0.0d0,  0.0d0,  0.0d0,    s45,   0.0d0,  0.0d0,    0.0d0,  0.0d0 & ! xyzz
     /),shape(gsphcar))

  ! gsphcar: l = 4
  !   spherical fchk order: m = 0, 1, -1, 2, -2, 3, -3, 4, -4
  !   Cartesian fchk order:   zzzz yzzz yyzz yyyz yyyy xzzz xyzz xyyz xyyy xxzz xxyz xxyy xxxz xxxy xxxx
  ! This is just a permutation of gsphcar where the rows are arranged in fchk's primitive order
  real*8 :: gsphcar_fchk(9,15) = reshape((/&
  !    0       1      -1       2      -2        3      -3         4      -4
       1d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! zzzz
     0.0d0,  0.0d0,    s10,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! yzzz
      -3d0,  0.0d0,  0.0d0, -s45_4,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! yyzz
     0.0d0,  0.0d0, -s45_8,  0.0d0,  0.0d0,   0.0d0, -s35_8,    0.0d0,  0.0d0,& ! yyyz
       d38,  0.0d0,  0.0d0,  s5_16,  0.0d0,   0.0d0,  0.0d0,   s35_64,  0.0d0,& ! yyyy
     0.0d0,    s10,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xzzz
     0.0d0,  0.0d0,  0.0d0,  0.0d0,    s45,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xyzz
     0.0d0, -s45_8,  0.0d0,  0.0d0,  0.0d0, -s315_8,  0.0d0,    0.0d0,  0.0d0,& ! xyyz
     0.0d0,  0.0d0,  0.0d0,  0.0d0, -s10_8,   0.0d0,  0.0d0,    0.0d0, -s35_4,& ! xyyy
      -3d0,  0.0d0,  0.0d0,  s45_4,  0.0d0,   0.0d0,  0.0d0,    0.0d0,  0.0d0,& ! xxzz
     0.0d0,  0.0d0, -s45_8,  0.0d0,  0.0d0,   0.0d0, s315_8,    0.0d0,  0.0d0,& ! xxyz
       d34,  0.0d0,  0.0d0,  0.0d0,  0.0d0,   0.0d0,  0.0d0, -s315_16,  0.0d0,& ! xxyy
     0.0d0, -s45_8,  0.0d0,  0.0d0,  0.0d0,   s35_8,  0.0d0,    0.0d0,  0.0d0,& ! xxxz
     0.0d0,  0.0d0,  0.0d0,  0.0d0, -s10_8,   0.0d0,  0.0d0,    0.0d0,  s35_4,& ! xxxy
       d38,  0.0d0,  0.0d0, -s5_16,  0.0d0,   0.0d0,  0.0d0,   s35_64,  0.0d0 & ! xxxx
     /),shape(gsphcar))

  ! threshold for primitive exponential range
  real*8, parameter :: rprim_thres = 1d-12

contains

  !> Terminate the wfn arrays
  subroutine wfn_end(f)
    class(molwfn), intent(inout) :: f
    
    if (allocated(f%icenter)) deallocate(f%icenter)
    if (allocated(f%itype)) deallocate(f%itype)
    if (allocated(f%d2ran)) deallocate(f%d2ran)
    if (allocated(f%e)) deallocate(f%e)
    if (allocated(f%occ)) deallocate(f%occ)
    if (allocated(f%cmo)) deallocate(f%cmo)
    if (allocated(f%icenter_edf)) deallocate(f%icenter_edf)
    if (allocated(f%itype_edf)) deallocate(f%itype_edf)
    if (allocated(f%e_edf)) deallocate(f%e_edf)
    if (allocated(f%c_edf)) deallocate(f%c_edf)
    if (allocated(f%xat)) deallocate(f%xat)

  end subroutine wfn_end

  !> Read the molecular geometry from an xyz file
  subroutine wfn_read_xyz_geometry(file,n,x,z,name)
    use types, only: atom
    use tools_io, only: fopen_read, zatguess, isinteger, ferror, faterr, fclose
    use param, only: bohrtoa

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    real*8, allocatable, intent(inout) :: x(:,:) !< Coordinates (bohr)
    integer, allocatable, intent(inout) :: z(:) !< Atomic numbers
    character*(10), allocatable, intent(inout) :: name(:) !< Atomic names

    character*4 :: atsym
    integer :: lu
    integer :: i, lp
    logical :: ok

    ! deallocate
    if (allocated(x)) deallocate(x)
    if (allocated(z)) deallocate(z)
    if (allocated(name)) deallocate(name)

    ! read the number of atoms
    lu = fopen_read(file)
    read (lu,*) n
    read (lu,*)
    allocate(x(3,n),z(n),name(n))

    do i = 1, n
       read (lu,*) atsym, x(:,i)
       z(i) = zatguess(atsym)
       if (z(i) <= 0) then
          ! maybe it's a number
          lp = 1
          ok = isinteger(z(i),atsym,lp)
          if (.not.ok) &
             call ferror('wfn_read_xyz_geometry','could not determine atomic number',faterr)
       end if
       name(i) = trim(adjustl(atsym))
       x(:,i) = x(:,i) / bohrtoa
    end do
    call fclose(lu)

  end subroutine wfn_read_xyz_geometry

  !> Read the molecular geometry from a wfn file
  subroutine wfn_read_wfn_geometry(file,n,x,z,name)
    use types, only: atom
    use tools_io, only: fopen_read, ferror, faterr, zatguess, fclose

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    real*8, allocatable, intent(inout) :: x(:,:) !< Coordinates (bohr)
    integer, allocatable, intent(inout) :: z(:) !< Atomic numbers
    character*(10), allocatable, intent(inout) :: name(:) !< Atomic names

    character*4 :: atsym, orbtyp
    integer :: lu
    integer :: i, i1, i2
    real*8 :: zreal

    ! deallocate
    if (allocated(x)) deallocate(x)
    if (allocated(z)) deallocate(z)
    if (allocated(name)) deallocate(name)

    lu = fopen_read(file)

    ! read the number of atoms
    read (lu,*)
    read (lu,101) orbtyp, i1, i2, n
    if (n <= 0) &
       call ferror('wfn_read_wfn_geometry','wrong number of atoms',faterr)
    allocate(x(3,n),z(n),name(n))

    ! read the geometry
    do i = 1, n
       read(lu,106) atsym, x(:,i), zreal
       z(i) = zatguess(atsym)
       name(i) = trim(adjustl(atsym))
    end do
101 format (4X,A4,10X,3(I5,15X))
106 format(2X,A2,20X,3F12.8,10X,F5.1)

    call fclose(lu)

  end subroutine wfn_read_wfn_geometry

  !> Read the molecular geometry from a wfx file
  subroutine wfn_read_wfx_geometry(file,n,x,z,name)
    use types, only: atom
    use tools_io, only: fopen_read, getline_raw, ferror, faterr, nameguess, fclose

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    real*8, allocatable, intent(inout) :: x(:,:) !< Coordinates (bohr)
    integer, allocatable, intent(inout) :: z(:) !< Atomic numbers
    character*(10), allocatable, intent(inout) :: name(:) !< Atomic names

    integer :: lu
    character(len=:), allocatable :: line, line2
    integer :: i

    ! deallocate
    if (allocated(x)) deallocate(x)
    if (allocated(z)) deallocate(z)
    if (allocated(name)) deallocate(name)

    lu = fopen_read(file)

    ! read the number of atoms
    do while (getline_raw(lu,line,.true.))
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Number of Nuclei>") then
             read (lu,*) n
             exit
          endif
       endif
    enddo
    if (n == 0) &
       call ferror("wfn_read_wfx_geometry","Number of Nuclei tag not found",faterr)
    allocate(name(n),z(n),x(3,n))

    ! read the geometry
    rewind(lu)
    do while (.true.)
       read(lu,'(A)',end=20) line
       line2 = adjustl(line)
       line = line2
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Atomic Numbers>") then
             z = wfx_read_integers(lu,n)
             do i = 1, n
                name(i) = nameguess(z(i))
             end do
          elseif (trim(line) == "<Nuclear Cartesian Coordinates>") then
             x = reshape(wfx_read_reals1(lu,3*n),shape(x))
          endif
       endif
    enddo
20  continue

    call fclose(lu)

  end subroutine wfn_read_wfx_geometry

  !> Read the molecular geometry from a fchk file
  subroutine wfn_read_fchk_geometry(file,n,x,z,name)
    use types, only: atom
    use tools_io, only: fopen_read, getline_raw, isinteger, ferror, faterr, nameguess, fclose

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    real*8, allocatable, intent(inout) :: x(:,:) !< Coordinates (bohr)
    integer, allocatable, intent(inout) :: z(:) !< Atomic numbers
    character*(10), allocatable, intent(inout) :: name(:) !< Atomic names

    integer :: lu, lp, i, j
    character(len=:), allocatable :: line
    logical :: ok
    real*8, allocatable :: xat(:)

    ! deallocate
    if (allocated(x)) deallocate(x)
    if (allocated(z)) deallocate(z)
    if (allocated(name)) deallocate(name)

    lu = fopen_read(file)

    ! read the number of atoms
    n = 0
    do while (getline_raw(lu,line,.true.))
       lp = 45
       if (line(1:15) == "Number of atoms") then
          ok = isinteger(n,line,lp)
          if (.not.ok) call ferror('wfn_read_fchk_geometry','could not read number of atoms',faterr)
          exit
       endif
    enddo
    if (n == 0) &
       call ferror("wfn_read_fchk_geometry","error reading number of atoms",faterr)

    ! read the geometry
    allocate(z(n),x(3,n),name(n),xat(3*n))
    rewind(lu)
    do while (getline_raw(lu,line))
       lp = 45
       if (line(1:29) == "Current cartesian coordinates") then
          do i = 0, (3*n-1)/5
             read(lu,'(5E16.8)') (xat(5*i+j),j=1,min(5,3*n-5*i))
          enddo
       elseif (line(1:14) == "Atomic numbers") then
          do i = 0, (n-1)/6
             read(lu,'(6I12)') (z(6*i+j),j=1,min(6,n-6*i))
          enddo
       endif
    enddo
    do i = 1, n
       x(:,i) = xat((i-1)*3+1:(i-1)*3+3)
       name(i) = nameguess(z(i))
    end do
    
    ! clean up
    deallocate(xat)
    call fclose(lu)

  end subroutine wfn_read_fchk_geometry

  !> Read the molecular geometry from a molden file (only tested with
  !> new psi4 molden files).
  subroutine wfn_read_molden_geometry(file,n,x,z,name)
    use types, only: atom
    use tools_io, only: fopen_read, lower, getline_raw, lgetword, ferror, faterr, nameguess, fclose
    use param, only: bohrtoa

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    real*8, allocatable, intent(inout) :: x(:,:) !< Coordinates (bohr)
    integer, allocatable, intent(inout) :: z(:) !< Atomic numbers
    character*(10), allocatable, intent(inout) :: name(:) !< Atomic names

    integer :: lu, lp, idum, i
    character(len=:), allocatable :: line, keyword, word1, word2
    character*(1024) :: fixword
    logical :: ok, isang

    ! deallocate
    if (allocated(x)) deallocate(x)
    if (allocated(z)) deallocate(z)
    if (allocated(name)) deallocate(name)

    lu = fopen_read(file)

    ! read the number of atoms
    n = 0
    do while(next_keyword())
       if (trim(lower(keyword)) == "atoms") then
          n = 0
          ok = getline_raw(lu,line,.true.)
          do while(index(lower(line),"[") == 0 .and. len(trim(line)) > 0)
             n = n + 1
             ok = getline_raw(lu,line,.true.)
          end do
          exit
       end if
    end do
    if (n == 0) &
       call ferror("wfn_read_molden_geometry","error reading number of atoms",faterr)
    allocate(x(3,n),z(n),name(n))

    ! read the geometry
    rewind(lu)
    do while(next_keyword())
       if (trim(lower(keyword)) == "atoms") exit
    end do
    if (len_trim(keyword) == 0) &
       call ferror("wfn_read_molden_geometry","error reading geometry",faterr)

    ! geometry header -> detect the units for the geometry
    lp = 1
    word1 = lgetword(line,lp)
    word2 = lgetword(line,lp)
    isang = (trim(lower(word2)) == "(angs)".or.trim(lower(word2)) == "(ang)") 

    ! read the atomic numbers and positions
    do i = 1, n
       ok = getline_raw(lu,line,.true.)
       read(line,*) fixword, idum, z(i), x(:,i)
       if (isang) x(:,i) = x(:,i) / bohrtoa
       name(i) = nameguess(z(i))
    end do

    call fclose(lu)

  contains

    function next_keyword()

      integer :: istart, iend
      logical :: next_keyword

      keyword = ""
      next_keyword = .false.

      do while(.true.)
         ok = getline_raw(lu,line,.false.)
         if (.not.ok) return
         if (index(lower(line),"[") > 0) exit
      end do
      next_keyword =.true.
      istart = index(lower(line),"[") + 1
      iend = index(lower(line),"]") - 1
      keyword = line(istart:iend)
      return

    end function next_keyword

  end subroutine wfn_read_molden_geometry

  !> Read the wavefunction from a wfn file
  subroutine read_wfn(f,file)
    use tools_io, only: fopen_read, zatguess, ferror, faterr, fclose
    class(molwfn), intent(inout) :: f !< Output field
    character*(*), intent(in) :: file !< Input file

    integer :: luwfn, nat, iz, istat, i, j, num1, num2
    integer :: nalpha, ioc
    real*8 :: x(3), zreal, ene0, ene
    character*4 :: orbtyp
    character*2 :: dums, elem
    logical :: isfrac
    character*8 :: dum1

    f%useecp = .false.
    f%issto = .false.
    f%hasvirtual = .false.

    ! read number of atoms, primitives, orbitals
    luwfn = fopen_read(file)
    read (luwfn,*)
    read (luwfn,101) orbtyp, f%nmoocc, f%npri, nat
    f%nmoall = f%nmoocc
    f%nalpha_virt = 0
 
    ! atomic positions and numbers
    do i = 1, nat
       read(luwfn,106) elem, x, zreal
       iz = zatguess(elem)
       if (iz /= nint(zreal)) f%useecp = .true.
    end do
 
    ! center assignments, types of primitives
    if (allocated(f%icenter)) deallocate(f%icenter)
    allocate(f%icenter(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_wfn','could not allocate memory for icenter',faterr)
    if (allocated(f%itype)) deallocate(f%itype)
    allocate(f%itype(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_wfn','could not allocate memory for itype',faterr)
    read(luwfn,102) (f%icenter(i),i=1,f%npri)
    read(luwfn,102) (f%itype(i),i=1,f%npri)
    if (any(f%itype(1:f%npri) > 56)) then
       call ferror("read_wfn","primitive type not supported",faterr)
    endif

    ! primitive exponents
    if (allocated(f%e)) deallocate(f%e)
    allocate(f%e(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_wfn','could not allocate memory for exponents',faterr)
    read(luwfn,103) (f%e(i),i=1,f%npri)
 
    ! deal with ecps
    dums=""
    do while (dums.ne."MO")
       read (luwfn,'(A2)') dums
    enddo
    backspace(luwfn)
 
    ! occupations and orbital coefficients
    if (allocated(f%occ)) deallocate(f%occ)
    allocate(f%occ(f%nmoocc),stat=istat)
    if (istat /= 0) call ferror('read_wfn','could not allocate memory for occupations',faterr)
    if (allocated(f%cmo)) deallocate(f%cmo)
    allocate(f%cmo(f%nmoocc,f%npri),stat=istat)
    if (istat /= 0) call ferror('read_wfn','could not allocate memory for orbital coefficients',faterr)
    isfrac = .false.
    num1 = 0
    num2 = 0
    ene0 = -1d30
    nalpha = -1
    do i = 1, f%nmoocc
       read(luwfn,104) f%occ(i), ene
       read(luwfn,105) (f%cmo(i,j),j=1,f%npri)
       ioc = nint(f%occ(i))
       if (abs(ioc-f%occ(i)) > 1d-10) then
          isfrac = .true.
       else if (ioc == 1) then
          num1 = num1 + 1
       else if (ioc == 2) then
          num2 = num2 + 1
       endif
       if (ene < ene0-1d-3) nalpha = i-1
       ene0 = ene
    end do
    read(luwfn,*) dum1
 
    ! determine the type of wavefunction
    ! 0 - restricted, 1 - unrestricted, 2 - fractional
    if (isfrac) then
       f%wfntyp = wfn_frac
       f%nalpha = 0
    else if (num1 == 0) then
       f%wfntyp = wfn_rhf
       f%nalpha = nalpha
       f%nalpha_virt = 0
    else if (num2 == 0) then
       f%wfntyp = wfn_uhf
       f%nalpha = nalpha
       f%nalpha_virt = 0
    else
       call ferror("read_wfn","restricted-open wfn files not supported",faterr)
    endif
    call fclose(luwfn)

    ! calculate the range of each primitive (in distance^2)
    call calculate_d2ran(f)

101  format (4X,A4,10X,3(I5,15X))
102  format(20X,20I3)
103  format(10X,5E14.7)
104  format(35X,F12.7,15X,F12.6)
105  format(5(E16.8))
106  format(2X,A2,20X,3F12.8,10X,F5.1)

  end subroutine read_wfn

  !> Read the wavefunction from a wfx file
  subroutine read_wfx(f,file)
    use tools_io, only: fopen_read, getline_raw, ferror, faterr, fclose
    class(molwfn), intent(inout) :: f !< Output field
    character*(*), intent(in) :: file !< Input file

    integer :: luwfn, ncore, istat, i, num1, num2, ioc
    character(len=:), allocatable :: line, line2
    logical :: isfrac

    f%useecp = .false.
    f%issto = .false.
    f%hasvirtual = .false.
    f%nalpha_virt = 0

    ! first pass
    luwfn = fopen_read(file)
    f%nmoocc = 0
    ncore = 0
    f%npri = 0
    f%nedf = 0
    do while (getline_raw(luwfn,line))
       if (len_trim(line) < 1) exit
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Number of Occupied Molecular Orbitals>") then
             read (luwfn,*) f%nmoocc
             f%nmoall = f%nmoocc
          elseif (trim(line) == "<Number of Core Electrons>") then
             read (luwfn,*) ncore
          elseif (trim(line) == "<Number of Primitives>") then
             read(luwfn,*) f%npri
          elseif (trim(line) == "<Number of EDF Primitives>") then
             read(luwfn,*) f%nedf
          endif
       endif
    enddo
    
    if (f%nmoocc == 0) call ferror("read_wfx","Number of Occupied Molecular Orbitals tag not found",faterr)
    if (f%npri == 0) call ferror("read_wfx","Number of Primitives tag not found",faterr)
    if (ncore > 0) f%useecp = .true.

    ! allocate memory
    allocate(f%icenter(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_wfx','could not allocate memory for icenter',faterr)
    allocate(f%itype(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_wfx','could not allocate memory for itype',faterr)
    allocate(f%e(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_wfx','could not allocate memory for exponents',faterr)
    allocate(f%occ(f%nmoocc),stat=istat)
    if (istat /= 0) call ferror('read_wfx','could not allocate memory for occupations',faterr)
    allocate(f%cmo(f%nmoocc,f%npri),stat=istat)
    if (istat /= 0) call ferror('read_wfx','could not allocate memory for orbital coefficients',faterr)
    if (f%nedf > 0) then
       allocate(f%icenter_edf(f%nedf),stat=istat)
       if (istat /= 0) call ferror('read_wfx','could not allocate memory for icenter',faterr)
       allocate(f%itype_edf(f%nedf),stat=istat)
       if (istat /= 0) call ferror('read_wfx','could not allocate memory for itype',faterr)
       allocate(f%e_edf(f%nedf),stat=istat)
       if (istat /= 0) call ferror('read_wfx','could not allocate memory for exponents',faterr)
       allocate(f%c_edf(f%nedf),stat=istat)
       if (istat /= 0) call ferror('read_wfx','could not allocate memory for orbital coefficients',faterr)
    endif

    ! second pass
    rewind(luwfn)
    do while (.true.)
       read(luwfn,'(A)',end=20) line
       line2 = adjustl(line)
       line = line2
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Primitive Centers>") then
             f%icenter = wfx_read_integers(luwfn,f%npri)
          elseif (trim(line) == "<Primitive Types>") then
             f%itype = wfx_read_integers(luwfn,f%npri)
             if (any(f%itype(1:f%npri) > 56)) &
                call ferror("read_wfx","primitive type not supported",faterr)
          elseif (trim(line) == "<Primitive Exponents>") then
             f%e = wfx_read_reals1(luwfn,f%npri)
          elseif (trim(line) == "<Molecular Orbital Occupation Numbers>") then
             f%occ = wfx_read_reals1(luwfn,f%nmoocc)
          elseif (trim(line) == "<Molecular Orbital Primitive Coefficients>") then
             read(luwfn,*)
             do i = 1, f%nmoocc
                read(luwfn,*)
                read(luwfn,*)
                f%cmo(i,:) = wfx_read_reals1(luwfn,f%npri)
             enddo
          elseif (trim(line) == "<EDF Primitive Centers>") then
             f%icenter_edf = wfx_read_integers(luwfn,f%nedf)
          elseif (trim(line) == "<EDF Primitive Types>") then
             f%itype_edf = wfx_read_integers(luwfn,f%nedf)
             if (any(f%itype_edf(1:f%nedf) > 56)) &
                call ferror("read_wfx","primitive type not supported",faterr)
          elseif (trim(line) == "<EDF Primitive Exponents>") then
             f%e_edf = wfx_read_reals1(luwfn,f%nedf)
          elseif (trim(line) == "<EDF Primitive Coefficients>") then
             f%c_edf = wfx_read_reals1(luwfn,f%nedf)
          endif
       endif
    enddo
20  continue

    ! wavefunction type
    isfrac = .false.
    num1 = 0
    num2 = 0
    do i = 1, f%nmoocc
       ioc = nint(f%occ(i))
       if (abs(ioc-f%occ(i)) > 1d-10) then
          isfrac = .true.
       else if (ioc == 1) then
          num1 = num1 + 1
       else if (ioc == 2) then
          num2 = num2 + 1
       endif
    end do
    if (isfrac) then
       f%wfntyp = wfn_frac
    else if (num1 == 0) then
       f%wfntyp = wfn_rhf
    else if (num2 == 0) then
       f%wfntyp = wfn_uhf
    else
       call ferror("read_wfx","restricted-open wfx files not supported",faterr)
    endif
    
    ! calculate the range of each primitive (in distance^2)
    call calculate_d2ran(f)

    ! done
    call fclose(luwfn)

  end subroutine read_wfx

  !> Read the wavefunction from a Gaussian formatted checkpoint file (fchk)
  subroutine read_fchk(f,file,readvirtual)
    use tools_io, only: fopen_read, getline_raw, isinteger, faterr, ferror, fclose
    use param, only: pi
    class(molwfn), intent(inout) :: f !< Output field
    character*(*), intent(in) :: file !< Input file
    logical, intent(in) :: readvirtual !< Read the virtual orbitals

    character(len=:), allocatable :: line
    integer :: lp, i, j, k, k1, k2, nn, nm, nl, nc, ns, ncar, nsph
    integer :: luwfn, nelec, nalpha, nbassph, nbascar, nbeta, ncshel, nshel, lmax
    integer :: istat, ityp, nmoread, namoread, nmoalla, nmoallb
    logical :: ok, isecp
    real*8 :: norm, cons
    integer, allocatable :: ishlt(:), ishlpri(:), ishlat(:), itemp(:,:)
    real*8, allocatable :: exppri(:), ccontr(:), pccontr(:), motemp(:), mocoef(:,:)
    real*8, allocatable :: cnorm(:), cpri(:), rtemp(:,:)
    logical, allocatable :: icdup(:)

    ! translation between primitive ordering fchk -> critic2
    !           1   2 3 4    5  6  7  8  9 10    11  12  13  14  15  16  17  18  19  20
    ! fchk:     s   x y z   xx yy zz xy xz yz   xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz
    ! critic2:  s   x y z   xx yy zz xy xz yz   xxx yyy zzz xxy xxz yyz xyy xzz yzz xyz
    !
    !             21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   
    ! fchk:     zzzz yzzz yyzz yyyz yyyy xzzz xyzz xyyz xyyy xxzz xxyz xxyy xxxz xxxy xxxx
    ! critic2:  xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz 

    ! reorder the primitives within the same shell
    ! typtrans(fchk) = critic2
    integer, parameter :: typtrans(35) = (/&
       !1   2  3  4    5  6  7  8  9  10    11  12  13  14  15  16  17  18  19  20
       1,   2, 3, 4,   5, 6, 7, 8, 9, 10,   11, 12, 13, 17, 14, 15, 18, 19, 16, 20,&
       !21 22  23  24  25  26  27  28  29  30  31  32  33  34  35
       23, 29, 32, 27, 22, 28, 35, 34, 26, 31, 33, 30, 25, 24, 21&
       /)

    ! no ecps for now
    f%useecp = .false.
    f%issto = .false.

    ! first pass: dimensions
    luwfn = fopen_read(file)
    f%wfntyp = wfn_rhf

    ! first pass: dimensions
    isecp = .false.
    nelec = 0
    nalpha = 0
    nmoalla = 0
    nmoallb = 0
    do while (getline_raw(luwfn,line,.false.))
       lp = 50
       if (line(1:19) == "Number of electrons") then
          ok = isinteger(nelec,line,lp)
       elseif (line(1:25) == "Number of alpha electrons") then
          ok = isinteger(nalpha,line,lp)
       elseif (line(1:25) == "Number of basis functions") then
          ok = isinteger(nbassph,line,lp)
       elseif (line(1:24) == "Number of beta electrons") then
          ok = isinteger(nbeta,line,lp)
       elseif (line(1:27) == "Number of contracted shells") then
          ok = isinteger(ncshel,line,lp)
       elseif (line(1:26) == "Number of primitive shells") then
          ok = isinteger(nshel,line,lp)
       elseif (line(1:24) == "Highest angular momentum") then
          ok = isinteger(lmax,line,lp)
       elseif (line(1:22) == "Alpha Orbital Energies") then
          ok = isinteger(nmoalla,line,lp)
       elseif (line(1:21) == "Beta Orbital Energies") then
          ok = isinteger(nmoallb,line,lp)
          f%wfntyp = wfn_uhf
       elseif (line(1:8) == "ECP-LMax") then
          isecp = .true.
       endif
    enddo

    ! ECPs not implemented yet
    if (isecp) call ferror("read_fchk","ECPs not supported.",faterr)
    if (nelec == 0) call ferror("read_fchk","nelec = 0",faterr)
    if (nmoalla == 0) call ferror("read_fchk","nmoall = 0",faterr)

    ! Count the number of MOs
    if (f%wfntyp == wfn_rhf) then
       if (mod(nelec,2) == 1) call ferror("read_fchk","odd nelec but closed-shell wavefunction",faterr)
       f%nmoocc = nelec / 2
       f%nalpha = nelec / 2
       nmoread = nmoalla
       namoread = nmoalla
    else if (f%wfntyp == wfn_uhf) then
       f%nmoocc = nelec
       f%nalpha = nalpha
       nmoread = (nmoalla+nmoallb)
       namoread = nmoalla
    endif
    if (.not.readvirtual) then
       nmoread = f%nmoocc
       namoread = f%nalpha
       f%hasvirtual = .false.
       f%nmoall = nmoread
       f%nalpha_virt = 0
    else
       f%hasvirtual = (f%nmoall /= f%nmoocc)
       f%nmoall = nmoread
       f%nalpha_virt = nmoalla - nalpha
    end if

    ! allocate sutff
    allocate(f%occ(f%nmoocc),stat=istat)
    if (istat /= 0) call ferror('read_fchk','could not allocate memory for occ',faterr)
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),stat=istat)
    if (istat /= 0) call ferror('read_fchk','could not allocate memory for shell data',faterr)
    allocate(exppri(nshel),ccontr(nshel),pccontr(nshel),stat=istat)
    if (istat /= 0) call ferror('read_fchk','could not allocate memory for primitive data',faterr)

    ! type of wavefunction -> occupations
    if (f%wfntyp == wfn_uhf) then
       f%occ = 1
    else
       f%occ = 2
    end if

    ! rewind
    rewind(luwfn)

    ! second pass
    allocate(motemp(nbassph*nmoread),stat=istat)
    if (istat /= 0) call ferror('read_fchk','could not allocate memory for MO coefs',faterr)
    do while (getline_raw(luwfn,line,.false.))
       lp = 45
       if (line(1:11) == "Shell types") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)') (ishlt(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:30) == "Number of primitives per shell") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)') (ishlpri(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:30) == "Shell to atom map") then
          do i = 0, (ncshel-1)/6
             read(luwfn,'(6I12)') (ishlat(6*i+j),j=1,min(6,ncshel-6*i))
          enddo
       elseif (line(1:19) == "Primitive exponents") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)') (exppri(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:24) == "Contraction coefficients") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)') (ccontr(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:31) == "P(S=P) Contraction coefficients") then
          do i = 0, (nshel-1)/5
             read(luwfn,'(5E16.8)') (pccontr(5*i+j),j=1,min(5,nshel-5*i))
          enddo
       elseif (line(1:21) == "Alpha MO coefficients") then
          read(luwfn,'(5E16.8)') (motemp(i),i=1,namoread*nbassph)
       elseif (line(1:21) == "Beta MO coefficients") then
          read(luwfn,'(5E16.8)') (motemp(i),i=namoread*nbassph+1,nmoread*nbassph)
       endif
    enddo

    ! we are done with the file
    call fclose(luwfn)

    ! unfold sp shells next
    allocate(icdup(ncshel))
    nshel = 0
    ncshel = 0
    do i = 1, size(icdup)
       if (ishlt(i) == -1) then
          icdup(i) = .true.
          ncshel = ncshel + 2
          nshel = nshel + 2 * ishlpri(i)
       else
          icdup(i) = .false.
          ncshel = ncshel + 1
          nshel = nshel + ishlpri(i)
       endif
    end do

    ! transfer the information to the temporary sp-unfolded arrays
    allocate(itemp(ncshel,3),rtemp(nshel,2))
    nn = 0
    nm = 0
    nl = 0
    do i = 1, size(icdup)
       if (icdup(i)) then
          nn = nn + 1
          itemp(nn,1) = 0
          itemp(nn,2) = ishlpri(i)
          itemp(nn,3) = ishlat(i)
          rtemp(nl+1:nl+ishlpri(i),1) = exppri(nm+1:nm+ishlpri(i))
          rtemp(nl+1:nl+ishlpri(i),2) = ccontr(nm+1:nm+ishlpri(i))
          nl = nl + ishlpri(i)

          nn = nn + 1
          itemp(nn,1) = 1
          itemp(nn,2) = ishlpri(i)
          itemp(nn,3) = ishlat(i)
          rtemp(nl+1:nl+ishlpri(i),1) = exppri(nm+1:nm+ishlpri(i))
          rtemp(nl+1:nl+ishlpri(i),2) = pccontr(nm+1:nm+ishlpri(i))
          nl = nl + ishlpri(i)
       else
          nn = nn + 1
          itemp(nn,1) = ishlt(i)
          itemp(nn,2) = ishlpri(i)
          itemp(nn,3) = ishlat(i)
          rtemp(nl+1:nl+ishlpri(i),1) = exppri(nm+1:nm+ishlpri(i))
          rtemp(nl+1:nl+ishlpri(i),2) = ccontr(nm+1:nm+ishlpri(i))
          nl = nl + ishlpri(i)
       endif
       nm = nm + ishlpri(i)
    end do

    ! move the sp-unfolded information back and reallocate
    deallocate(ishlt,ishlpri,ishlat,exppri,ccontr,pccontr,icdup)
    allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),exppri(nshel),ccontr(nshel))
    ishlt = itemp(:,1)
    ishlpri = itemp(:,2)
    ishlat = itemp(:,3)
    exppri = rtemp(:,1)
    ccontr = rtemp(:,2)
    deallocate(itemp,rtemp)

    ! count the number of primitives and basis functions
    f%npri = 0
    nbascar = 0
    nbassph = 0
    do i = 1, ncshel
       ityp = ishlt(i)
       nbascar = nbascar + nshlt(abs(ityp))
       nbassph = nbassph + nshlt(ityp)
       f%npri = f%npri + nshlt(abs(ishlt(i))) * ishlpri(i)
    enddo
    
    ! convert spherical basis functions to Cartesian and build the mocoef
    ! deallocate the temporary motemp
    allocate(mocoef(nmoread,nbascar))
    nc = 0
    ns = 0
    do j = 1, ncshel
       nsph = nshlt(ishlt(j))
       ncar = nshlt(abs(ishlt(j)))
       if (nsph == ncar) then
          do i = 1, nmoread
             mocoef(i,nc+1:nc+ncar) = motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph)
          end do
       elseif (ishlt(j) == -2) then
          do i = 1, nmoread
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),dsphcar)
          end do
       elseif (ishlt(j) == -3) then
          do i = 1, nmoread
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),fsphcar)
          end do
       elseif (ishlt(j) == -4) then
          do i = 1, nmoread
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),gsphcar_fchk)
          end do
       else
          call ferror('read_fchk','h and higher primitives not supported yet',faterr)
       endif
       ns = ns + nsph
       nc = nc + ncar
    end do
    deallocate(motemp)

    ! Assign primitive center and type, exponents, etc.
    allocate(f%icenter(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_fchk','could not allocate memory for icenter',faterr)
    allocate(f%itype(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_fchk','could not allocate memory for itype',faterr)
    allocate(f%e(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_fchk','could not allocate memory for exponents',faterr)
    allocate(f%cmo(nmoread,f%npri),cpri(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_fchk','could not allocate memory for coeffs',faterr)

    ! normalize the primitive coefficients without the angular part
    allocate(cnorm(maxval(ishlpri)))
    nm = 0
    nn = 0
    do i = 1, ncshel
       do j = jshl0(abs(ishlt(i))), jshl1(abs(ishlt(i)))
          ityp = typtrans(j)
          ! primitive coefficients normalized
          do k = 1, ishlpri(i)
             cnorm(k) = ccontr(nm+k) * gnorm(ityp,exppri(nm+k)) 
          end do

          ! normalization constant for the basis function
          norm = 0d0
          do k1 = 1, ishlpri(i)
             do k2 = 1, ishlpri(i)
                norm = norm + cnorm(k1) * cnorm(k2) / (exppri(nm+k1)+exppri(nm+k2))**(abs(ishlt(i))+3d0/2d0)
             end do
          end do
          cons = pi**(3d0/2d0) * dfacm1(2*abs(ishlt(i))) / 2**(abs(ishlt(i)))
          norm = 1d0 / sqrt(norm * cons)

          ! gaussian fchk: multiply by sqrt((2lx-1)!! * (2ly-1)!! * (2lz-1)!! / (2l-1)!!)
          ! only for Cartesian primitives
          if (ishlt(i) >= 0) then
             if (ityp >= 8 .and. ityp <= 10) then
                norm = norm * sqrt(3d0)
             elseif (ityp >= 14 .and. ityp <= 19) then
                norm = norm * sqrt(5d0)
             elseif (ityp == 20) then
                norm = norm * sqrt(15d0)
             else if (ityp >= 24 .and. ityp <= 29) then
                norm = norm * sqrt(7d0)
             else if (ityp >= 30 .and. ityp <= 32) then
                norm = norm * sqrt(35d0/3d0)
             else if (ityp >= 33 .and. ityp <= 35) then
                norm = norm * sqrt(35d0)
             end if
          end if

          ! calculate and assign the normalized primitive coefficients
          do k = 1, ishlpri(i)
             nn = nn + 1
             cpri(nn) = cnorm(k) * norm
          end do
       end do
       nm = nm + ishlpri(i)
    end do
    deallocate(cnorm)

    ! build the wavefunction coefficients for the primitives
    nn = 0
    nm = 0
    nl = 0
    do i = 1, ncshel
       do j = jshl0(abs(ishlt(i))), jshl1(abs(ishlt(i)))
          ityp = typtrans(j)
          nl = nl + 1
          do k = 1, ishlpri(i)
             nn = nn + 1
             f%icenter(nn) = ishlat(i)
             f%itype(nn) = ityp
             f%e(nn) = exppri(nm+k)
             f%cmo(:,nn) = cpri(nn) * mocoef(:,nl)
          end do
       end do
       nm = nm + ishlpri(i)
    end do

    ! calculate the range of each primitive (in distance^2)
    call calculate_d2ran(f)

    ! clean up
    deallocate(ishlt,ishlpri,ishlat)
    deallocate(exppri,ccontr)
    deallocate(mocoef)

  end subroutine read_fchk

  !> Read the wavefunction from a molden file. Only molden files
  !> generated by psi4 (in the 2016 or later implementation) have been
  !> tested.
  subroutine read_molden(f,file,readvirtual)
    use tools_io, only: fopen_read, getline_raw, lower, ferror, faterr, lgetword, &
       isinteger, isreal, fclose, uout
    use param, only: pi
    class(molwfn), intent(inout) :: f !< Output field
    character*(*), intent(in) :: file !< Input file
    logical, intent(in) :: readvirtual !< Read the virtual orbitals

    character(len=:), allocatable :: line, keyword, word, word1, word2
    logical :: is5d, is7f, is9g, isalpha, ok, issto, isgto, isocc
    integer :: luwfn, istat, ityp
    integer :: i, j, k, k1, k2, ni, nj, nc, ns, nm, nn, nl, ncar, nsph
    integer :: nat, nelec, nalpha, nalphamo, nbetamo, ncshel, nshel, nbascar, nbassph
    integer :: idum, lp, lnmoa, lnmob, lnmo, lnmoav, lnmobv, ix, iy, iz, ir
    real*8 :: rdum, norm, cons
    integer, allocatable :: ishlt(:), ishlpri(:), ishlat(:)
    real*8, allocatable :: exppri(:), ccontr(:), motemp(:), cpri(:), mocoef(:,:), cnorm(:)

    ! Guide to the variables in this routine:
    !   issto and isgto = whether the wavefunction is given in GTOs or STOs.
    !                     Both can not be true at the same time.
    !
    !   ncshel = number of contraction shells
    !     -> ishlt(icshel) = angular momentum for this shell (s=0,p=1,d=2,dsph=-2,etc.)
    !     -> ishlpri(icshel) = number of primitives in this shell
    !     -> ishlat(icshel) = atom index where this shell is based
    !
    !   nshel = number of primitives (not counting the angular parts)
    !     -> exppri(ishel) = primitive exponent
    !     -> ccontr(ishel) = contraction coefficient for this primitive
    !
    !   m%npri = number of primitives (including angular parts)
    !
    !   nbassph = number of basis functions (spherical)
    !
    !   nbascar = number of basis functions (Cartesian)
    !
    !   lnmo = number of molecular orbitals, regardless of occupation
    !   m%nmo = number of occupied molecular orbitals
    !     -> m%occ(imo) = occupation of the molecular orbitals
    !
    !   m%nelec = number of electrons
    !   nalpha = number of alpha electrons
    !   m%wfntyp = type of wavefunction (0 = closed, 1 = open)
    !   m%charge = charge
    !   m%mult = multiplicity
    ! 
    ! For STOs, ncshel = nshel = m%npri = number of basis functions.
    ! The ishlt has a special meaning - the exponents for x, y, z, and r
    ! are packed with a 100 stride factor.

    ! translation between primitive ordering molden -> critic2
    !           1   2 3 4    5  6  7  8  9 10    11  12  13  14  15  16  17  18  19  20
    ! molden:   s   x y z   xx yy zz xy xz yz   xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz
    ! critic2:  s   x y z   xx yy zz xy xz yz   xxx yyy zzz xxy xxz yyz xyy xzz yzz xyz
    !
    !             21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   
    ! molden:   xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz
    ! critic2:  xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz 
    !
    ! h primitives and higher not supported in molden format.
    !
    !     typtrans(molden) = critic2
    integer, parameter :: typtrans(35) = (/&
       !1   2  3  4    5  6  7  8  9  10    11  12  13  14  15  16  17  18  19  20
       1,   2, 3, 4,   5, 6, 7, 8, 9, 10,   11, 12, 13, 17, 14, 15, 18, 19, 16, 20,&
       !21 22  23  24  25  26  27  28  29  30  31  32  33  34  35
       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35&
       /)

    ! initialize
    f%useecp = .false.
    f%issto = .false.
    is5d = .false.
    is7f = .false.
    is9g = .false.
    luwfn = fopen_read(file)
    nat = 0
    f%wfntyp = -1
    line = ""
    issto = .false.
    isgto = .false.

    ! parse the molden file, first pass -> read dimensions prior to allocation
    do while(next_keyword())
       if (trim(keyword) == "atoms") then
          ! read the number of atoms 
          ok = getline_raw(luwfn,line,.true.)
          do while(index(lower(line),"[") == 0 .and. len(trim(line)) > 0)
             nat = nat + 1
             ok = getline_raw(luwfn,line,.false.)
             if (.not.ok) exit
          end do
       elseif (trim(keyword) == "gto") then
          ! read the number of shells and primitives
          isgto = .true.
          ncshel = 0
          nshel = 0
          if (nat == 0) call ferror("read_molden","gto found but no atoms",faterr)
          do i = 1, nat
             ok = getline_raw(luwfn,line,.true.)
             ok = getline_raw(luwfn,line,.true.)
             do while (index(line,".") /= 0)
                word = ""
                read (line,*) word, idum, rdum
                ncshel = ncshel + 1
                nshel = nshel + idum
                do j = 1, idum
                   ok = getline_raw(luwfn,line,.true.)
                end do
                ok = getline_raw(luwfn,line,.false.)
                if (.not.ok) exit
             end do
             if (.not.ok) exit
          end do
       elseif (trim(keyword) == "sto") then
          ! read the number of shells and primitives
          issto = .true.
          f%npri = 0
          if (nat == 0) call ferror("read_molden","sto found but no atoms",faterr)
          do while (.true.)
             ok = getline_raw(luwfn,line,.true.)
             if (.not.ok) exit
             if (index(lower(line),"[") > 0) exit
             f%npri = f%npri + 1
          end do
       elseif (trim(keyword) == "mo") then
          ! read the number of MOs, number of electrons, and number of alpha electrons
          nelec = 0
          nalpha = 0
          nalphamo = 0
          nbetamo = 0
          do while(.true.)
             ok = getline_raw(luwfn,line,.false.)
             if (.not.ok) exit
             if (index(lower(line),"[") /= 0 .or. len_trim(line) == 0) exit
             if (index(lower(line),"ene=") /= 0) then
                ! spin
                ok = getline_raw(luwfn,line,.true.)
                isalpha = (index(lower(line),"alpha") > 0)
                if (isalpha) then
                   nalphamo = nalphamo + 1
                else
                   nbetamo = nbetamo + 1
                end if

                ! occupation
                ok = getline_raw(luwfn,line,.true.)
                word = ""
                read (line,*) word, rdum
                idum = nint(rdum)
                if (abs(rdum-idum) > 1d-6) then
                   write (uout,'("Fractional occupations are not supported yet for molden files.")')
                   write (uout,'("If you need this, please e-mail the critic2 developer.")')
                   call ferror('read_molden','Can not do fractional occupations with molden',faterr)
                end if
                if (idum == 1 .or. idum == 2) then
                   nelec = nelec + idum
                   if (isalpha) then
                      nalpha = nalpha + 1
                   endif
                   if (f%wfntyp < 0 .and. idum == 2) f%wfntyp = wfn_rhf
                   if (idum == 1) f%wfntyp = wfn_uhf
                elseif (idum == 0) then
                   continue
                else
                   call ferror('read_molden','wrong integer occupation',faterr)
                endif
             end if
          end do
       elseif (trim(keyword) == "5d" .or. trim(keyword) == "5d7f") then
          is5d = .true.
          is7f = .true.
          ok = getline_raw(luwfn,line,.false.)
       elseif (trim(keyword) == "5d10f") then
          is5d = .true.
          is7f = .false.
          ok = getline_raw(luwfn,line,.false.)
       elseif (trim(keyword) == "7f") then
          is5d = .false.
          is7f = .true.
          ok = getline_raw(luwfn,line,.false.)
       elseif (trim(keyword) == "9g") then
          is9g = .true.
          ok = getline_raw(luwfn,line,.false.)
       else 
          ! must be a keyword I don't know about -> skip
          ok = getline_raw(luwfn,line,.false.)
       end if
       if (.not.ok) exit
    end do
    
    if (isgto.and.issto) then
       call ferror('read_molden','Both [GTO] and [STO] blocks are present',faterr)
    else if (.not.isgto.and..not.issto) then
       call ferror('read_molden','No [GTO] or [STO] blocks present',faterr)
    else if (isgto) then
       f%issto = .false.
    else
       f%issto = .true.
    end if

    ! type of wavefunction -> number of MOs
    if (f%wfntyp == wfn_uhf) then
       f%nmoocc = nelec
    else
       if (mod(nelec,2) == 1) call ferror("read_molden","odd nelec but closed-shell wavefunction",faterr)
       f%nmoocc = nelec / 2
    end if
    f%nalpha = nalpha
    if (readvirtual) then
       f%nmoall = nalphamo + nbetamo
       f%hasvirtual = (f%nmoall /= f%nmoocc)
       f%nalpha_virt = nalphamo - nalpha
    else
       f%nalpha_virt = 0
       f%nmoall = f%nmoocc
       f%hasvirtual = .false.
    end if

    ! allocate stuff
    allocate(f%occ(f%nmoall),stat=istat)
    if (istat /= 0) call ferror('read_molden','alloc. memory for occ',faterr)
    f%occ = 0d0

    ! type of wavefunction -> occupations
    if (f%wfntyp == wfn_uhf) then
       f%occ(1:f%nmoocc) = 1d0
    else
       f%occ(1:f%nmoocc) = 2d0
    end if

    ! second pass
    rewind(luwfn)

    if (isgto) then
       !! GTO wavefunction !!
       allocate(ishlt(ncshel),ishlpri(ncshel),ishlat(ncshel),stat=istat)
       if (istat /= 0) call ferror('read_molden','alloc. memory for shell types',faterr)
       allocate(exppri(nshel),ccontr(nshel),stat=istat)
       if (istat /= 0) call ferror('read_molden','alloc. memory for prim. shells',faterr)

       ! Basis set specification 
       do while(getline_raw(luwfn,line,.true.))
          if (trim(lower(line)) == "[gto]") exit
       end do
       ni = 0
       nj = 0
       do i = 1, nat
          ok = getline_raw(luwfn,line,.true.)
          ok = getline_raw(luwfn,line,.true.)
          ! atomic blocks are separated by blank lines
          do while (len_trim(line) > 0)
             ni = ni + 1

             ! header -> l, npri, ???
             lp = 1
             word = lgetword(line,lp)
             ok = ok .and. isinteger(idum,line,lp)
             if (.not.ok) &
                call ferror("read_molden","error reading gto block",faterr)

             ! read exponents and coefficients
             do j = 1, idum
                nj = nj + 1
                ok = getline_raw(luwfn,line,.true.)
                read(line,*) exppri(nj), ccontr(nj)
             end do

             ! assign atomic centers, number of primitives, and shell types.
             ishlat(ni) = i
             ishlpri(ni) = idum
             if (lower(trim(word)) == "s") then
                ishlt(ni) = 0
             else if (lower(trim(word)) == "p") then
                ishlt(ni) = 1
             else if (lower(trim(word)) == "sp") then
                write (uout,'("SP shells are not supported yet for molden files.")')
                write (uout,'("If you need this, e-mail the critic2 developer.")')
                call ferror("read_molden","can't handle SP in molden format",faterr)
             else if (lower(trim(word)) == "d") then
                if (is5d) then
                   ishlt(ni) = -2
                else
                   ishlt(ni) = 2
                end if
             else if (lower(trim(word)) == "f") then
                if (is7f) then
                   ishlt(ni) = -3
                else
                   ishlt(ni) = 3
                end if
             else if (lower(trim(word)) == "g") then
                if (is9g) then
                   ishlt(ni) = -4
                else
                   ishlt(ni) = 4
                end if
             else
                write (uout,'("Shells with angular momentum higher than g are not supported yet for molden files.")')
                write (uout,'("If you need this, e-mail the critic2 developer.")')
                call ferror("read_molden","basis set type >g not supported in molden files",faterr)
             endif
             ok = getline_raw(luwfn,line,.false.)
             if (.not.ok) exit
          end do
       end do

       ! Count the number of primitives and basis functions
       f%npri = 0
       nbascar = 0
       nbassph = 0
       do i = 1, ncshel
          ityp = ishlt(i)
          nbascar = nbascar + nshlt(abs(ityp))
          nbassph = nbassph + nshlt(ityp)
          f%npri = f%npri + nshlt(abs(ityp)) * ishlpri(i)
       enddo
    else
       !! STO wavefunction !!
       allocate(f%itype(f%npri),f%icenter(f%npri),stat=istat)
       if (istat /= 0) call ferror('read_molden','alloc. memory for shell types',faterr)
       allocate(f%e(f%npri),ccontr(f%npri),stat=istat)
       if (istat /= 0) call ferror('read_molden','alloc. memory for prim. shells',faterr)

       do while(getline_raw(luwfn,line,.true.))
          if (trim(lower(line)) == "[sto]") exit
       end do

       ! Basis set specification 
       f%ixmaxsto = 0
       do i = 1, f%npri
          ok = getline_raw(luwfn,line,.true.)
          read(line,*) f%icenter(i), ix, iy, iz, ir, f%e(i), ccontr(i)
          f%itype(i) = ix + 100 * (iy + 100 * (iz + 100 * ir))
          f%ixmaxsto(1) = max(f%ixmaxsto(1),ix)
          f%ixmaxsto(2) = max(f%ixmaxsto(2),iy)
          f%ixmaxsto(3) = max(f%ixmaxsto(3),iz)
          f%ixmaxsto(4) = max(f%ixmaxsto(4),ir)
       end do
       nbassph = f%npri ! for the allocation
    end if

    ! Allocate space to read the the molecular orbital information
    allocate(motemp(nbassph*f%nmoall),stat=istat)
    if (istat /= 0) call ferror('read_molden','alloc. memory for MO coefs',faterr)
    motemp = 0d0

    ! advance to the MO coefficients
    do while(getline_raw(luwfn,line,.true.))
       if (trim(lower(line)) == "[mo]") exit
    end do

    ! read the MO coefficients
    lnmoa = 0
    lnmob = nalpha
    lnmoav = f%nmoocc
    lnmobv = f%nmoocc + f%nalpha_virt
    main: do while(.true.)
       ok = getline_raw(luwfn,line,.false.)
       if (.not.ok) exit main

       ! is this an alpha electron?
       if (index(lower(line),"spin=") /= 0 .and. f%wfntyp /= wfn_rhf) then
          lp = 1
          isalpha = (index(lower(line),"alpha") > 0)
       end if

       ! read the orbital if it is occupied or readvirtual is true
       ! assumes that the MO coefficients come right after the "Occup" tag
       if (index(lower(line),"occup=") /= 0) then
          lp = 1
          word1 = lgetword(line,lp)
          ok = isreal(rdum,line,lp)
          if (.not.ok) & 
             call ferror("read_molden","error reading mo block",faterr)
          idum = nint(rdum)
          isocc = (idum > 0)
          if (isocc .or. readvirtual) then
             if (isalpha.and.isocc) then
                lnmoa = lnmoa + 1
                lnmo = lnmoa
             elseif (isalpha.and..not.isocc) then
                lnmoav = lnmoav + 1
                lnmo = lnmoav
             else if (.not.isalpha.and.isocc) then
                lnmob = lnmob + 1
                lnmo = lnmob
             elseif (.not.isalpha.and..not.isocc) then
                lnmobv = lnmobv + 1
                lnmo = lnmobv
             end if

             do while (.true.)
                ok = getline_raw(luwfn,line,.false.)
                if (.not.ok) exit main
                if (index(lower(line),"sym=") > 0) exit
                read(line,*,end=99) idum, rdum
                motemp((lnmo-1)*nbassph+idum) = rdum
             end do
          end if
       end if
    end do main

    ! we are done with the file
    call fclose(luwfn)

    ! check we read the correct number of MOs
    if (lnmoa /= nalpha) &
       call ferror('read_molden','inconsistent number of MOs (lnmoa) in the second pass',faterr)
    if (lnmob /= f%nmoocc) &
       call ferror('read_molden','inconsistent number of MOs (lnmob) in the second pass',faterr)
    if (lnmoav /= f%nmoocc+f%nalpha_virt) &
       call ferror('read_molden','inconsistent number of MOs (lnmoav) in the second pass',faterr)
    if (lnmobv /= f%nmoall) &
       call ferror('read_molden','inconsistent number of MOs (lnmobv) in the second pass',faterr)

    ! STOs: build the MO coefficients and exit
    if (issto) then
       allocate(f%cmo(f%nmoall,f%npri),stat=istat)
       if (istat /= 0) call ferror('read_molden','could not allocate memory for coeffs',faterr)

       ! STOs: just copy the molecular orbital coefficients
       do j = 1, f%npri
          do i = 1, f%nmoall
             f%cmo(i,j) = ccontr(j) * motemp((i-1)*f%npri+j)
          end do
       end do
       deallocate(motemp,ccontr)
       
       ! calculate the range of each primitive (in distance^2)
       call calculate_d2ran(f)

       return
    end if

    !! From this point onward, only GTO wavefunctions !!
    ! convert spherical basis functions to Cartesian and build the mocoef
    ! deallocate the temporary motemp
    allocate(mocoef(f%nmoall,nbascar),stat=istat)
    if (istat /= 0) call ferror('read_molden','alloc. memory for mocoef',faterr)
    nc = 0
    ns = 0
    do j = 1, ncshel
       nsph = nshlt(ishlt(j))
       ncar = nshlt(abs(ishlt(j)))
       if (nsph == ncar) then
          do i = 1, f%nmoall
             mocoef(i,nc+1:nc+ncar) = motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph)
          end do
       elseif (ishlt(j) == -2) then
          do i = 1, f%nmoall
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),dsphcar)
          end do
       elseif (ishlt(j) == -3) then
          do i = 1, f%nmoall
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),fsphcar)
          end do
       elseif (ishlt(j) == -4) then
          do i = 1, f%nmoall
             mocoef(i,nc+1:nc+ncar) = matmul(motemp((i-1)*nbassph+ns+1:(i-1)*nbassph+ns+nsph),gsphcar)
          end do
       else
          call ferror('read_molden','h and higher primitives not supported in molden format',faterr)
       endif
       ns = ns + nsph
       nc = nc + ncar
    end do
    deallocate(motemp)

    ! from this point on, only Cartesians
    ishlt = abs(ishlt)

    ! Allocate molecule arrays for primitive center, type, exponents, etc.
    allocate(f%icenter(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_molden','could not allocate memory for icenter',faterr)
    allocate(f%itype(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_molden','could not allocate memory for itype',faterr)
    allocate(f%e(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_molden','could not allocate memory for exponents',faterr)
    allocate(f%cmo(f%nmoall,f%npri),cpri(f%npri),stat=istat)
    if (istat /= 0) call ferror('read_molden','could not allocate memory for coeffs',faterr)

    ! normalize the primitive coefficients without the angular part
    allocate(cnorm(maxval(ishlpri)))
    nm = 0
    nn = 0
    do i = 1, ncshel
       do j = jshl0(ishlt(i)), jshl1(ishlt(i))
          ityp = typtrans(j)

          ! primitive coefficients normalized
          do k = 1, ishlpri(i)
             cnorm(k) = ccontr(nm+k) * gnorm(ityp,exppri(nm+k)) 
          end do

          ! normalization constant for the basis function
          norm = 0d0
          do k1 = 1, ishlpri(i)
             do k2 = 1, ishlpri(i)
                norm = norm + cnorm(k1) * cnorm(k2) / (exppri(nm+k1)+exppri(nm+k2))**(ishlt(i)+3d0/2d0)
             end do
          end do
          cons = pi**(3d0/2d0) * dfacm1(2*ishlt(i)) / 2**(ishlt(i))
          norm = 1d0 / sqrt(norm * cons)

          ! calculate and assign the normalized primitive coefficients
          do k = 1, ishlpri(i)
             nn = nn + 1
             cpri(nn) = cnorm(k) * norm
          end do
       end do
       nm = nm + ishlpri(i)
    end do
    deallocate(cnorm)

    ! build the wavefunction coefficients for the primitives
    nn = 0
    nm = 0
    nl = 0
    do i = 1, ncshel
       do j = jshl0(ishlt(i)), jshl1(ishlt(i))
          ityp = typtrans(j)
          nl = nl + 1
          do k = 1, ishlpri(i)
             nn = nn + 1
             f%icenter(nn) = ishlat(i)
             f%itype(nn) = ityp
             f%e(nn) = exppri(nm+k)
             f%cmo(:,nn) = cpri(nn) * mocoef(:,nl)
          end do
       end do
       nm = nm + ishlpri(i)
    end do

    ! calculate the range of each primitive (in distance^2)
    call calculate_d2ran(f)

    ! clean up and exit
    deallocate(ishlt,ishlpri,ishlat)
    deallocate(exppri,ccontr)
    deallocate(mocoef,cpri)

    return
99  continue
    call ferror("read_molden","error reading molden file",faterr)
  contains

    function next_keyword()

      integer :: istart, iend
      logical :: next_keyword

      keyword = ""
      next_keyword = .false.

      do while(.true.)
         if (index(lower(line),"[") > 0) exit
         ok = getline_raw(luwfn,line,.false.)
         if (.not.ok) return
      end do
      next_keyword =.true.
      istart = index(lower(line),"[") + 1
      iend = index(lower(line),"]") - 1
      keyword = lower(line(istart:iend))
      return

    end function next_keyword

  end subroutine read_molden

  !> Register structural information.
  subroutine register_struct(f,ncel,atcel)
    use types, only: celatom
    class(molwfn), intent(inout) :: f
    integer, intent(in) :: ncel
    type(celatom), intent(in) :: atcel(:)

    integer :: i

    f%nat = ncel
    if (allocated(f%xat)) deallocate(f%xat)
    allocate(f%xat(3,f%nat))
    do i = 1, f%nat
       f%xat(:,i) = atcel(i)%r
    end do
    
  end subroutine register_struct

  !> Calculate the normalization factor of a primitive of a given type
  !> with exponent a.
  function gnorm(type,a) result(N)
    use tools_io, only: ferror, faterr
    use param, only: pi
    integer, intent(in) :: type
    real*8, intent(in) :: a
    real*8 :: N

    if (type == 1) then
       ! 1
       ! x
       N = 2**(3d0/4d0) * a**(3d0/4d0) / pi**(3d0/4d0)
    else if (type >= 2 .and. type <= 4) then
       ! 2 3 4
       ! x y z
       N = 2**(7d0/4d0) * a**(5d0/4d0) / pi**(3d0/4d0)
    else if (type >= 5 .and. type <= 7) then
       ! 5  6  7
       ! xx yy zz
       N = 2**(11d0/4d0) * a**(7d0/4d0) / pi**(3d0/4d0) / sqrt(3d0)
    else if (type >= 8 .and. type <= 10) then
       ! 7  8  9
       ! xy xz yz
       N = 2**(11d0/4d0) * a**(7d0/4d0) / pi**(3d0/4d0) 
    else if (type >= 11 .and. type <= 13) then
       ! 11  12  13
       ! xxx yyy zzz
       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0) / sqrt(15d0)
    else if (type >= 14 .and. type <= 19) then
       ! 14  15  16  17  18  19  
       ! xyy xxy xxz xzz yzz yyz 
       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0) / sqrt(3d0)
    else if (type == 20) then
       ! 20
       ! xyz
       N = 2**(15d0/4d0) * a**(9d0/4d0) / pi**(3d0/4d0)
    else if (type >= 21 .and. type <= 23) then
       ! 21   22   23   
       ! xxxx yyyy zzzz 
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0) / sqrt(35d0)
    else if (type >= 24 .and. type <= 29) then
       ! 24   25   26   27   28   29  
       ! xxxy xxxz xyyy yyyz xzzz yzzz
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0) / sqrt(5d0)
    else if (type >= 30 .and. type <= 32) then
       ! 30   31   32  
       ! xxyy xxzz yyzz
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0) * sqrt(3d0)
    else if (type >= 33 .and. type <= 35) then
       ! 33   34   35   
       ! xxyz xyyz xyzz 
       N = 2**(19d0/4d0) * a**(11d0/4d0) / pi**(3d0/4d0)
    else
       call ferror("gnorm","fixme: primitive type not supported",faterr)
    endif

  endfunction gnorm

  !> Determine the density and derivatives at a given target point.
  !> xpos is in cartesian coordiantes and assume that the molecule has
  !> been displaced to the center of a big cube. Same transformation
  !> as in load xyz/wfn/wfx. This routine is thread-safe.
  subroutine rho2(f,xpos,nder,rho,grad,h,gkin,vir,stress,xmo)
    use tools_io, only: ferror, faterr
    class(molwfn), intent(in) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    integer, intent(in) :: nder  !< Number of derivatives
    real*8, intent(out) :: rho !< Density
    real*8, intent(out) :: grad(3) !< Gradient
    real*8, intent(out) :: h(3,3) !< Hessian 
    real*8, intent(out) :: gkin !< G(r), kinetic energy density
    real*8, intent(out) :: vir !< Virial field
    real*8, intent(out) :: stress(3,3) !< Schrodinger stress tensor
    real*8, allocatable, intent(out), optional :: xmo(:) !< Values of the MO

    integer, parameter :: imax(0:2) = (/1,4,10/)

    real*8 :: al, ex, xl(3,0:2), xl2
    integer :: iat, ityp, i, imo1, l(3)
    integer :: imo, ix
    real*8 :: phi(1:f%nmoocc,imax(nder))
    real*8 :: dd(3,f%nat), d2(f%nat)
    real*8 :: hh(6), aocc
    
    real*8, parameter :: stoeps = 1d-40
    integer, parameter :: li(3,56) = reshape((/&
       0,0,0, & ! s
       1,0,0, 0,1,0, 0,0,1, & ! p
       2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1, & !d
       3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 0,2,1, &
       1,2,0, 1,0,2, 0,1,2, 1,1,1,& ! f
       4,0,0, 0,4,0, 0,0,4, 3,1,0, 3,0,1, 1,3,0, 0,3,1, 1,0,3,&
       0,1,3, 2,2,0, 2,0,2, 0,2,2, 2,1,1, 1,2,1, 1,1,2,& ! g
       0,0,5, 0,1,4, 0,2,3, 0,3,2, 0,4,1, 0,5,0, 1,0,4, 1,1,3,&
       1,2,2, 1,3,1, 1,4,0, 2,0,3, 2,1,2, 2,2,1, 2,3,0, 3,0,2,&
       3,1,1, 3,2,0, 4,0,1, 4,1,0, 5,0,0/),shape(li)) ! h

    ! calculate the MO values and derivatives at the point
    if (f%issto) then
       call f%calculate_mo_sto(xpos,phi,1,1,f%nmoocc,nder)
    else
       call f%calculate_mo_gto(xpos,phi,1,1,f%nmoocc,nder)
    end if

    ! save the MO values
    if (present(xmo)) then
       if (allocated(xmo)) then
          if (size(xmo) /= imo1) deallocate(xmo)
       end if
       if (.not.allocated(xmo)) allocate(xmo(imo1))
       xmo(1:imo1) = phi(1:imo1,1)
    end if

    ! calculate the density, etc.
    rho = 0d0
    grad = 0d0
    hh = 0d0
    gkin = 0d0
    vir = 0d0
    stress = 0d0
    do imo = 1, f%nmoocc
       aocc = f%occ(imo) 
       rho = rho + aocc * phi(imo,1) * phi(imo,1)
       if (nder>0) then
          grad = grad + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
          gkin = gkin + aocc * (phi(imo,2)*phi(imo,2) + phi(imo,3)*phi(imo,3) + phi(imo,4)*phi(imo,4))
       endif
       if (nder>1) then
          hh(1:3) = hh(1:3) + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
          hh(4) = hh(4) + 2 * aocc * (phi(imo,1)*phi(imo,8)+phi(imo,2)*phi(imo,3))
          hh(5) = hh(5) + 2 * aocc * (phi(imo,1)*phi(imo,9)+phi(imo,2)*phi(imo,4))
          hh(6) = hh(6) + 2 * aocc * (phi(imo,1)*phi(imo,10)+phi(imo,3)*phi(imo,4))
          stress(1,1) = stress(1,1) + aocc * (phi(imo,1) * phi(imo,5)  - phi(imo,2)*phi(imo,2))
          stress(1,2) = stress(1,2) + aocc * (phi(imo,1) * phi(imo,8)  - phi(imo,2)*phi(imo,3))
          stress(1,3) = stress(1,3) + aocc * (phi(imo,1) * phi(imo,9)  - phi(imo,2)*phi(imo,4))
          stress(2,2) = stress(2,2) + aocc * (phi(imo,1) * phi(imo,6)  - phi(imo,3)*phi(imo,3))
          stress(2,3) = stress(2,3) + aocc * (phi(imo,1) * phi(imo,10) - phi(imo,3)*phi(imo,4))
          stress(3,3) = stress(3,3) + aocc * (phi(imo,1) * phi(imo,7)  - phi(imo,4)*phi(imo,4))
       endif
    enddo
    stress(2,1) = stress(1,2)
    stress(3,1) = stress(1,3)
    stress(3,2) = stress(2,3)
    stress = stress * 0.5d0
    gkin = 0.5d0 * gkin
    vir = stress(1,1)+stress(2,2)+stress(3,3)

    ! core contribution to the density and its derivatives
    if (f%nedf > 0) then
       do i = 1, f%nedf
          ityp = f%itype_edf(i)
          iat = f%icenter_edf(i)
          al = f%e_edf(i)
          ex = exp(-al * d2(iat))

          l = li(1:3,ityp)
          do ix = 1, 3
             if (l(ix) == 0) then
                xl(ix,0) = 1d0
                xl(ix,1) = 0d0
                xl(ix,2) = 0d0
             else if (l(ix) == 1) then
                xl(ix,0) = dd(ix,iat)
                xl(ix,1) = 1d0
                xl(ix,2) = 0d0
             else if (l(ix) == 2) then
                xl(ix,0) = dd(ix,iat) * dd(ix,iat)
                xl(ix,1) = 2d0 * dd(ix,iat)
                xl(ix,2) = 2d0
             else if (l(ix) == 3) then
                xl(ix,0) = dd(ix,iat) * dd(ix,iat) * dd(ix,iat)
                xl(ix,1) = 3d0 * dd(ix,iat) * dd(ix,iat)
                xl(ix,2) = 6d0 * dd(ix,iat)
             else if (l(ix) == 4) then
                xl2 = dd(ix,iat) * dd(ix,iat)
                xl(ix,0) = xl2 * xl2
                xl(ix,1) = 4d0 * xl2 * dd(ix,iat)
                xl(ix,2) = 12d0 * xl2
             else if (l(ix) == 5) then
                xl2 = dd(ix,iat) * dd(ix,iat)
                xl(ix,0) = xl2 * xl2 * dd(ix,iat)
                xl(ix,1) = 5d0 * xl2 * xl2
                xl(ix,2) = 20d0 * xl2 * dd(ix,iat)
             else
                call ferror('wfn_rho2','power of L not supported',faterr)
             end if
          end do
          
          rho = rho + f%c_edf(i) * xl(1,0)*xl(2,0)*xl(3,0)*ex
          if (nder > 0) then
             grad(1) = grad(1) + f%c_edf(i) * (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * xl(2,0)*xl(3,0)*ex
             grad(2) = grad(2) + f%c_edf(i) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(1,0)*xl(3,0)*ex
             grad(3) = grad(3) + f%c_edf(i) * (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * xl(1,0)*xl(2,0)*ex
             if (nder > 1) then
                hh(1) = hh(1) + f%c_edf(i) * (xl(1,2)-2*al*(2*l(1)+1)*xl(1,0) + 4*al*al*dd(1,iat)**(l(1)+2)) * xl(2,0)*xl(3,0)*ex
                hh(2) = hh(2) + f%c_edf(i) * (xl(2,2)-2*al*(2*l(2)+1)*xl(2,0) + 4*al*al*dd(2,iat)**(l(2)+2)) * xl(3,0)*xl(1,0)*ex
                hh(3) = hh(3) + f%c_edf(i) * (xl(3,2)-2*al*(2*l(3)+1)*xl(3,0) + 4*al*al*dd(3,iat)**(l(3)+2)) * xl(1,0)*xl(2,0)*ex
                hh(4) = hh(4) + f%c_edf(i) * (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(3,0)*ex
                hh(5) = hh(5) + f%c_edf(i) * (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * xl(2,0)*ex
                hh(6) = hh(6) + f%c_edf(i) * (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(1,0)*ex
             endif
          endif
       end do
    end if

    ! re-order the hessian
    do i = 1, 3
       h(i,i) = hh(i)
    end do
    h(1,2) = hh(4)
    h(2,1) = h(1,2)
    h(1,3) = hh(5)
    h(3,1) = h(1,3)
    h(2,3) = hh(6)
    h(3,2) = h(2,3)

  end subroutine rho2

  !> Calculate the MO values at position xpos (Cartesian).
  subroutine calculate_mo(f,xpos,phi,philb,imo0,imo1,nder)
    class(molwfn), intent(in) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    real*8, intent(inout) :: phi(:,:) !< array for the final values
    integer, intent(in) :: philb !< lower bound for phi
    integer, intent(in) :: imo0 !< first MO
    integer, intent(in) :: imo1 !< last MO
    integer, intent(in) :: nder !< number of derivatives
    
    if (f%issto) then
       ! STO wavefunction
       call f%calculate_mo_sto(xpos,phi,philb,imo0,imo1,nder)
    else
       ! GTO wavefunction
       call f%calculate_mo_gto(xpos,phi,philb,imo0,imo1,nder)
    end if

  end subroutine calculate_mo

  !> Calculate the MO values at position xpos (Cartesian). STO version.
  subroutine calculate_mo_sto(f,xpos,phi,philb,imo0,imo1,nder)
    class(molwfn), intent(in) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    real*8, intent(inout) :: phi(:,:) !< array for the final values
    integer, intent(in) :: philb !< lower bound for phi
    integer, intent(in) :: imo0 !< first MO
    integer, intent(in) :: imo1 !< last MO
    integer, intent(in) :: nder !< number of derivatives

    integer :: i, j
    real*8 :: xx(4), al, ex, xratio(3)
    real*8 :: dx(4,-2:maxval(f%ixmaxsto),f%nat)
    integer :: ipri, imo, phimo, ityp, ixx(4), iat
    real*8 :: fprod(-2:0,-2:0,-2:0), f0r, f1r, f2r

    real*8, parameter :: stoeps = 1d-40

    ! calculate distances and their powers
    dx = 1d0
    dx(:,-2:-1,:) = 0d0
    do iat = 1, f%nat
       xx(1:3) = xpos - f%xat(:,iat)
       xx(4) = sqrt(xx(1)*xx(1)+xx(2)*xx(2)+xx(3)*xx(3))
       dx(:,1,iat) = xx
       ! positive powers
       do i = 1, 4
          do j = 2, f%ixmaxsto(i)
             dx(i,j,iat) = dx(i,j-1,iat) * xx(i)
          end do
       end do
    enddo

    ! build the MO values at the point
    phi = 0d0
    do ipri = 1, f%npri
       if (dx(4,1,f%icenter(ipri)) > f%d2ran(ipri)) cycle
       do imo = imo0, imo1
          phimo = imo - philb + 1

          ! unpack the exponents
          ityp = f%itype(ipri)
          ixx(1) = modulo(ityp,100)
          ityp = ityp / 100
          ixx(2) = modulo(ityp,100)
          ityp = ityp / 100
          ixx(3) = modulo(ityp,100)
          ixx(4) = ityp / 100

          ! calculate the exponentials
          iat = f%icenter(ipri)
          al = f%e(ipri)
          ex = exp(-al * dx(4,1,iat))

          ! MO value
          fprod = 0d0
          f0r = f%cmo(imo,ipri) * dx(4,ixx(4),iat) * ex
          fprod(0,0,0) = dx(1,ixx(1),iat) * dx(2,ixx(2),iat) * dx(3,ixx(3),iat)
          phi(phimo,1) = phi(phimo,1) + fprod(0,0,0) * f0r 
          if (nder > 0) then
             fprod(-1,0,0) = dx(1,ixx(1)-1,iat) * dx(2,ixx(2),iat) * dx(3,ixx(3),iat)
             fprod(0,-1,0) = dx(1,ixx(1),iat) * dx(2,ixx(2)-1,iat) * dx(3,ixx(3),iat)
             fprod(0,0,-1) = dx(1,ixx(1),iat) * dx(2,ixx(2),iat) * dx(3,ixx(3)-1,iat)
             f1r = (-al * dx(4,ixx(4),iat) + ixx(4) * dx(4,ixx(4)-1,iat)) * f%cmo(imo,ipri) * ex

             xratio(:) = dx(1:3,1,iat) / max(dx(4,1,iat),stoeps)
             phi(phimo,2) = phi(phimo,2) + ixx(1) * fprod(-1,0,0) * f0r + xratio(1) * fprod(0,0,0) * f1r
             phi(phimo,3) = phi(phimo,3) + ixx(2) * fprod(0,-1,0) * f0r + xratio(2) * fprod(0,0,0) * f1r
             phi(phimo,4) = phi(phimo,4) + ixx(3) * fprod(0,0,-1) * f0r + xratio(3) * fprod(0,0,0) * f1r

             if (nder > 1) then
                f2r = (al * al * dx(4,ixx(4),iat) - 2d0 * al * ixx(4) * dx(4,ixx(4)-1,iat) &
                   + ixx(4) * (ixx(4)-1) * dx(4,ixx(4)-2,iat)) * f%cmo(imo,ipri) * ex
                fprod(-2,0,0) = dx(1,ixx(1)-2,iat) * dx(2,ixx(2),iat) * dx(3,ixx(3),iat)
                fprod(0,-2,0) = dx(1,ixx(1),iat) * dx(2,ixx(2)-2,iat) * dx(3,ixx(3),iat)
                fprod(0,0,-2) = dx(1,ixx(1),iat) * dx(2,ixx(2),iat) * dx(3,ixx(3)-2,iat)
                fprod(-1,-1,0) = dx(1,ixx(1)-1,iat) * dx(2,ixx(2)-1,iat) * dx(3,ixx(3),iat)
                fprod(-1,0,-1) = dx(1,ixx(1)-1,iat) * dx(2,ixx(2),iat) * dx(3,ixx(3)-1,iat)
                fprod(0,-1,-1) = dx(1,ixx(1),iat) * dx(2,ixx(2)-1,iat) * dx(3,ixx(3)-1,iat)

                ! xx=5, yy=6, zz=7
                phi(phimo,5) = phi(phimo,5) + ixx(1) * (ixx(1)-1) * fprod(-2,0,0) * f0r &
                   + 2d0 * ixx(1) * fprod(-1,0,0) * xratio(1) * f1r & 
                   + (1 - xratio(1)*xratio(1)) * fprod(0,0,0) * f1r / max(dx(4,1,iat),stoeps) &
                   + fprod(0,0,0) * xratio(1)*xratio(1) * f2r
                phi(phimo,6) = phi(phimo,6) + ixx(2) * (ixx(2)-1) * fprod(0,-2,0) * f0r &
                   + 2d0 * ixx(2) * fprod(0,-1,0) * xratio(2) * f1r & 
                   + (1 - xratio(2)*xratio(2)) * fprod(0,0,0) * f1r / max(dx(4,1,iat),stoeps) &
                   + fprod(0,0,0) * xratio(2)*xratio(2) * f2r
                phi(phimo,7) = phi(phimo,7) + ixx(3) * (ixx(3)-1) * fprod(0,0,-2) * f0r &
                   + 2d0 * ixx(3) * fprod(0,0,-1) * xratio(3) * f1r & 
                   + (1 - xratio(3)*xratio(3)) * fprod(0,0,0) * f1r / max(dx(4,1,iat),stoeps) &
                   + fprod(0,0,0) * xratio(3)*xratio(3) * f2r

                ! xy=8, xz=9, yz=10
                phi(phimo,8) = phi(phimo,8) + ixx(1) * ixx(2) * fprod(-1,-1,0) * f0r &
                   + ixx(1) * fprod(-1,0,0) * xratio(2) * f1r & 
                   + ixx(2) * fprod(0,-1,0) * xratio(1) * f1r & 
                   + fprod(0,0,0) * xratio(1)*xratio(2) * (f2r - f1r / dx(4,1,iat))
                phi(phimo,9) = phi(phimo,9) + ixx(1) * ixx(3) * fprod(-1,0,-1) * f0r &
                   + ixx(1) * fprod(-1,0,0) * xratio(3) * f1r & 
                   + ixx(3) * fprod(0,0,-1) * xratio(1) * f1r & 
                   + fprod(0,0,0) * xratio(1)*xratio(3) * (f2r - f1r / dx(4,1,iat))
                phi(phimo,10) = phi(phimo,10) + ixx(2) * ixx(3) * fprod(0,-1,-1) * f0r &
                   + ixx(2) * fprod(0,-1,0) * xratio(3) * f1r & 
                   + ixx(3) * fprod(0,0,-1) * xratio(2) * f1r & 
                   + fprod(0,0,0) * xratio(2)*xratio(3) * (f2r - f1r / dx(4,1,iat))
             endif
          endif
       enddo
    enddo
    
  end subroutine calculate_mo_sto

  !> Calculate the MO values at position xpos (Cartesian). GTO version.
  subroutine calculate_mo_gto(f,xpos,phi,philb,imo0,imo1,nder)
    use tools_io, only: ferror, faterr
    class(molwfn), intent(in) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    real*8, intent(inout) :: phi(:,:) !< array for the final values
    integer, intent(in) :: philb !< lower bound for phi
    integer, intent(in) :: imo0 !< first MO
    integer, intent(in) :: imo1 !< last MO
    integer, intent(in) :: nder !< number of derivatives

    integer, parameter :: imax(0:2) = (/1,4,10/)

    integer :: l(3)
    integer :: iat, ityp, ipri, imo, ix, phimo
    real*8 :: dd(3,f%nat), d2(f%nat), al, ex, xl(3,0:2), xl2
    logical :: ldopri(f%npri)
    real*8 :: chi(f%npri,imax(nder))
    
    ! real*8, parameter :: stoeps = 1d-40
    integer, parameter :: li(3,56) = reshape((/&
       0,0,0, & ! s
       1,0,0, 0,1,0, 0,0,1, & ! p
       2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1, & !d
       3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 0,2,1, &
       1,2,0, 1,0,2, 0,1,2, 1,1,1,& ! f
       4,0,0, 0,4,0, 0,0,4, 3,1,0, 3,0,1, 1,3,0, 0,3,1, 1,0,3,&
       0,1,3, 2,2,0, 2,0,2, 0,2,2, 2,1,1, 1,2,1, 1,1,2,& ! g
       0,0,5, 0,1,4, 0,2,3, 0,3,2, 0,4,1, 0,5,0, 1,0,4, 1,1,3,&
       1,2,2, 1,3,1, 1,4,0, 2,0,3, 2,1,2, 2,2,1, 2,3,0, 3,0,2,&
       3,1,1, 3,2,0, 4,0,1, 4,1,0, 5,0,0/),shape(li)) ! h

    ! calculate distances
    do iat = 1, f%nat
       dd(:,iat) = xpos - f%xat(:,iat)
       d2(iat) = dd(1,iat)*dd(1,iat)+dd(2,iat)*dd(2,iat)+dd(3,iat)*dd(3,iat)
    enddo

    ! primitive values and their derivatives
    ldopri = .true.
    do ipri = 1, f%npri
       if (d2(f%icenter(ipri)) > f%d2ran(ipri)) then
          ldopri(ipri) = .false.
          cycle
       end if
       ityp = f%itype(ipri)
       iat = f%icenter(ipri)
       al = f%e(ipri)
       ex = exp(-al * d2(iat))
       
       l = li(1:3,ityp)
       do ix = 1, 3
          if (l(ix) == 0) then
             xl(ix,0) = 1d0
             xl(ix,1) = 0d0
             xl(ix,2) = 0d0
          else if (l(ix) == 1) then
             xl(ix,0) = dd(ix,iat)
             xl(ix,1) = 1d0
             xl(ix,2) = 0d0
          else if (l(ix) == 2) then
             xl(ix,0) = dd(ix,iat) * dd(ix,iat)
             xl(ix,1) = 2d0 * dd(ix,iat)
             xl(ix,2) = 2d0
          else if (l(ix) == 3) then
             xl(ix,0) = dd(ix,iat) * dd(ix,iat) * dd(ix,iat)
             xl(ix,1) = 3d0 * dd(ix,iat) * dd(ix,iat)
             xl(ix,2) = 6d0 * dd(ix,iat)
          else if (l(ix) == 4) then
             xl2 = dd(ix,iat) * dd(ix,iat)
             xl(ix,0) = xl2 * xl2
             xl(ix,1) = 4d0 * xl2 * dd(ix,iat)
             xl(ix,2) = 12d0 * xl2
          else if (l(ix) == 5) then
             xl2 = dd(ix,iat) * dd(ix,iat)
             xl(ix,0) = xl2 * xl2 * dd(ix,iat)
             xl(ix,1) = 5d0 * xl2 * xl2
             xl(ix,2) = 20d0 * xl2 * dd(ix,iat)
          else
             call ferror('wfn_rho2','power of L not supported',faterr)
          end if
       end do

       chi(ipri,1) = xl(1,0)*xl(2,0)*xl(3,0)*ex
       if (nder > 0) then
          chi(ipri,2) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * xl(2,0)*xl(3,0)*ex
          chi(ipri,3) = (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(1,0)*xl(3,0)*ex
          chi(ipri,4) = (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * xl(1,0)*xl(2,0)*ex
          if (nder > 1) then
             chi(ipri,5) = (xl(1,2)-2*al*(2*l(1)+1)*xl(1,0) + 4*al*al*dd(1,iat)**(l(1)+2)) * xl(2,0)*xl(3,0)*ex
             chi(ipri,6) = (xl(2,2)-2*al*(2*l(2)+1)*xl(2,0) + 4*al*al*dd(2,iat)**(l(2)+2)) * xl(3,0)*xl(1,0)*ex
             chi(ipri,7) = (xl(3,2)-2*al*(2*l(3)+1)*xl(3,0) + 4*al*al*dd(3,iat)**(l(3)+2)) * xl(1,0)*xl(2,0)*ex
             chi(ipri,8) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(3,0)*ex
             chi(ipri,9) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * xl(2,0)*ex
             chi(ipri,10)= (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(1,0)*ex
          endif
       endif
    enddo ! ipri = 1, npri

    ! build the MO values at the point
    phi = 0d0
    do ix = 1, imax(nder)
       do ipri = 1, f%npri
          if (.not.ldopri(ipri)) cycle
          do imo = imo0, imo1
             phimo = imo - philb + 1
             phi(phimo,ix) = phi(phimo,ix) + f%cmo(imo,ipri)*chi(ipri,ix)
          enddo
       enddo
    enddo
    
  end subroutine calculate_mo_gto

  !> Read a list of n integers from a logical unit
  function wfx_read_integers(lu,n) result(x)
    use tools_io, only: getline_raw, isinteger, ferror, faterr
    integer, intent(in) :: lu, n
    integer :: x(n)

    integer :: kk, lp, idum
    character(len=:), allocatable :: line
    logical :: ok

    kk = 0
    lp = 1
    ok = getline_raw(lu,line,.true.)
    do while(.true.)
       if (.not.isinteger(idum,line,lp)) then
          lp = 1
          ok = getline_raw(lu,line)
          if (.not.ok .or. line(1:2) == "</") exit
       else
          kk = kk + 1
          if (kk > n) call ferror("wfx_read_integers","exceeded size of the array",faterr)
          x(kk) = idum
       endif
    enddo

  endfunction wfx_read_integers

  !> Read a list of n reals from a logical unit
  function wfx_read_reals1(lu,n) result(x)
    use tools_io, only: getline_raw, isreal, ferror, faterr
    integer, intent(in) :: lu, n
    real*8 :: x(n)

    integer :: kk, lp
    real*8 :: rdum
    character(len=:), allocatable :: line
    logical :: ok

    kk = 0
    lp = 1
    ok = getline_raw(lu,line,.true.)
    do while(.true.)
       if (.not.isreal(rdum,line,lp)) then
          lp = 1
          ok = getline_raw(lu,line)
          if (.not.ok .or. line(1:1) == "<") exit
       else
          kk = kk + 1
          if (kk > n) call ferror("wfx_read_reals1","exceeded size of the array",faterr)
          x(kk) = rdum
       endif
    enddo

  endfunction wfx_read_reals1

  !> Calculate the distance range of each primitive
  subroutine calculate_d2ran(f)
    use tools_io, only: ferror, faterr
    class(molwfn), intent(inout) :: f

    integer :: istat, i

    if (allocated(f%d2ran)) deallocate(f%d2ran)
    allocate(f%d2ran(f%npri),stat=istat)
    if (istat /= 0) call ferror('calculate_d2ran','could not allocate memory for d2ran',faterr)
    do i = 1, f%npri
       f%d2ran(i) = -log(rprim_thres) / f%e(i)
    end do

  end subroutine calculate_d2ran
end module wfn_private
