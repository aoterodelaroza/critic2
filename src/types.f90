! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> User-defined types and overloaded reallocation procedures.
module types
  implicit none

  private
  public :: atom
  public :: celatom
  public :: anyatom
  public :: cp_type
  public :: scalar_value
  public :: scalar_value_noalloc
  public :: integrable
  public :: pointpropable
  public :: neighstar
  public :: realloc
  public :: gpathp

  ! overloaded functions
  interface realloc
     module procedure realloc_pointpropable
     module procedure realloc_integrable
     module procedure realloc_atom
     module procedure realloc_celatom
     module procedure realloc_anyatom
     module procedure realloc_cp
     module procedure realloc_gpathp
     module procedure realloc1l
     module procedure realloc1r
     module procedure realloc2r
     module procedure realloc3r
     module procedure realloc4r
     module procedure realloc5r
     module procedure realloc1i
     module procedure realloc2i
     module procedure realloc1c
     module procedure realloc1cmplx4
     module procedure realloc2cmplx4
     module procedure realloc4cmplx4
     module procedure realloc1cmplx8
     module procedure realloc5cmplx8
  end interface

  !> Non-equivalent atom list type (nneq)
  type atom
     real*8 :: x(3)   !< coordinates (crystallographic)
     real*8 :: r(3)   !< coordinates (cartesian)
     character*(10) :: name = "" !< name
     integer :: z = 0 !< atomic number
     real*8 :: qat = 0d0 !< ionic charge for promolecular densities (integer) and Ewald (fractional)
     integer :: mult  !< multiplicity
     real*8 :: rnn2   !< half the nearest neighbor distance
  end type atom
  
  !> Equivalent atom list (type)
  type celatom
     real*8 :: x(3) !< coordinates (crystallographic)
     real*8 :: r(3) !< coordinates (cartesian)
     integer :: idx !< corresponding atom from the non-equivalent atom list
     integer :: cidx !< corresponding equivalent atom from the complete atom list
     integer :: ir  !< rotation matrix to the representative equivalent atom
     integer :: ic  !< translation vector to the representative equivalent atom
     integer :: lvec(3) !< lattice vector to the representative equivalent atom
     integer :: lenv(3) !< lattice vector to the main cell atom (for environments)
  end type celatom

  !> Any atom in the crystal (not necessarily in the main cell), type
  type anyatom
     real*8 :: x(3) !< coordinates (crystallographic)
     real*8 :: r(3) !< coordinates (cartesian)
     integer :: idx !< corresponding atom from the non-equivalent atom list
     integer :: cidx !< corresponding atom from the complete atom list
     integer :: lvec(3) !< lattice vector to the atom in the complete atom list
     integer :: z !< atomic number
  end type anyatom

  !> Result of the evaluation of a scalar field
  type scalar_value
     ! basic
     real*8 :: f, fval, gf(3), hf(3,3), gfmod, gfmodval, del2f, del2fval
     ! kinetic energy density
     real*8 :: gkin
     ! schrodinger stress tensor
     real*8 :: stress(3,3)
     ! electronic potential energy density, virial field
     real*8 :: vir
     ! additional local properties
     real*8 :: hfevec(3,3), hfeval(3)
     integer :: r, s
     ! is it a nuclear position?
     logical :: isnuc
     ! molecular orbital values
     real*8, allocatable :: mo(:)
  end type scalar_value

  !> Result of the evaluation of a scalar field
  type scalar_value_noalloc
     ! basic
     real*8 :: f, fval, gf(3), hf(3,3), gfmod, gfmodval, del2f, del2fval
     ! kinetic energy density
     real*8 :: gkin
     ! schrodinger stress tensor
     real*8 :: stress(3,3)
     ! electronic potential energy density, virial field
     real*8 :: vir
     ! additional local properties
     real*8 :: hfevec(3,3), hfeval(3)
     integer :: r, s
     ! is it a nuclear position?
     logical :: isnuc
  end type scalar_value_noalloc

  !> Critical point type
  type cp_type
     ! Position and type
     real*8 :: x(3) !< Position (cryst. coords.)
     real*8 :: r(3) !< Position (cart. coords.)
     integer :: typ  !< Type of CP (-3,-1,1,3)
     integer :: typind !< Same as type, but mapped to (0,1,2,3)
     integer :: mult !< Multiplicity
     character*3 :: pg !< Point group symbol
     character*(10) :: name !< Name
     logical :: isdeg !< Is it a degenerate CP?
     logical :: isnuc !< Is it a nucleus?
     logical :: isnnm !< Is it a non-nuclear attractor/repulsor?

     ! beta-sphere 
     real*8 :: rbeta !< beta-sphere radius

     ! Properties at the CP
     type(scalar_value_noalloc) :: s  !< scalar value - evaluation of the reference field at the CP

     ! BCP and RCP ias properties
     integer :: ipath(2) !< Associated attractor (bcp) or repulsor (rcp), complete list
     integer :: ilvec(3,2) !< Lattice vector to shift the cp_(ipath) position of the actual attractor
     real*8 :: brdist(2) !< If b or r, distance to attractor/repulsor
     real*8 :: brpathlen(2) !< If b or r, distance to attractor/repulsor
     real*8 :: brang !< If b or r, angle wrt attractors/repulsors
     real*8 :: brvec(3) !< If b or r, the eigenvector along the bond (ring) path

     ! Complete list -> reduced CP list index and conversion
     integer :: idx !< Complete to non-equivalent list index
     integer :: ir !< Rotation matrix to the neq cp list
     integer :: ic !< Translation vector to the neq cp list 
     integer :: lvec(3) !< Lattice vector to the neq cp list 
  end type cp_type

  !> Information about an integrable field
  type integrable
     logical :: used = .false.
     integer :: itype
     integer :: fid
     character*(10) :: prop_name
     character*(2048) :: expr
     integer :: lmax
     real*8 :: x0(3)
  end type integrable

  !> Information about an point-property field
  type pointpropable
     character*(10) :: name
     character*(2048) :: expr
     integer :: ispecial ! 0 = expression, 1 = stress
     integer :: nf 
     integer, allocatable :: fused(:)
  end type pointpropable

  !> Neighbor star -> to analyze the connectivity in a molecular crystal
  type neighstar
     integer :: ncon = 0 !< number of neighbor for this atom
     integer, allocatable :: idcon(:) !< id (atcel) of the connected atom
     integer, allocatable :: lcon(:,:) !< lattice vector of the connected atom
  end type neighstar

  !> Point along a gradient path 
  type gpathp
     integer :: i
     real*8 :: x(3)
     real*8 :: r(3)
     real*8 :: f
     real*8 :: gf(3)
     real*8 :: hf(3,3)
  end type gpathp

contains
  
  !> Adapt the size of an allocatable 1D type(pointpropable) array
  subroutine realloc_pointpropable(a,nnew)
    type(pointpropable), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(pointpropable), allocatable :: temp(:)
    integer :: l1, u1

    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    if (u1 == nnew) return
    allocate(temp(l1:nnew))

    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc_pointpropable

  !> Adapt the size of an allocatable 1D type(integrable) array
  subroutine realloc_integrable(a,nnew)
    type(integrable), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(integrable), allocatable :: temp(:)
    integer :: l1, u1

    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    if (u1 == nnew) return
    allocate(temp(l1:nnew))

    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc_integrable

  !> Adapt the size of an allocatable 1D type(atom) array
  subroutine realloc_atom(a,nnew)
    type(atom), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(atom), allocatable :: temp(:)
    integer :: nold, i

    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)
    do i = nold+1,nnew
       a(i)%name = ""
       a(i)%z = 0
       a(i)%qat = 0d0
       a(i)%rnn2 = 0d0
    end do

  end subroutine realloc_atom

  !> Adapt the size of an allocatable 1D type(celatom) array
  subroutine realloc_celatom(a,nnew)
    type(celatom), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(celatom), allocatable :: temp(:)
    integer :: nold

    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_celatom

  !> Adapt the size of an allocatable 1D type(celatom) array
  subroutine realloc_anyatom(a,nnew)
    type(anyatom), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(anyatom), allocatable :: temp(:)
    integer :: nold

    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_anyatom

  !> Adapt the size of an allocatable 1D type(atom) array
  subroutine realloc_cp(a,nnew)
    type(cp_type), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew
    
    type(cp_type), allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_cp

  !> Adapt the size of an allocatable 1D type(atom) array
  subroutine realloc_gpathp(a,nnew)
    type(gpathp), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew
    
    type(gpathp), allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_gpathp

  !> Adapt the size of an allocatable 1D logical array
  subroutine realloc1l(a,nnew)
    logical, intent(inout), allocatable :: a(:) !< Input array, logical, 1D
    integer, intent(in) :: nnew !< new dimension
    
    logical, allocatable :: temp(:)
    integer :: l1, u1
    
    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    if (u1 == nnew) return
    allocate(temp(l1:nnew))
    
    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc1l

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc1r(a,nnew)
    real*8, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    real*8, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1r

  !> Adapt the size of an allocatable 2D real*8 array
  subroutine realloc2r(a,n1,n2)
    real*8, intent(inout), allocatable :: a(:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    real*8, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) then
       allocate(a(1:n1,1:n2))
       return
    end if
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    if (nold(1) == n1 .and. nold(2) == n2) return
    allocate(temp(n1,n2))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc3r(a,n1,n2,n3)
    real*8, intent(inout), allocatable :: a(:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3 !< new dimension
    
    real*8, allocatable :: temp(:,:,:)
    integer :: nold(3)
    
    if (.not.allocated(a)) then
       allocate(a(1:n1,1:n2,1:n3))
       return
    end if
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3) return
    allocate(temp(n1,n2,n3))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)))
    call move_alloc(temp,a)

  end subroutine realloc3r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc4r(a,n1,n2,n3,n4)
    real*8, intent(inout), allocatable :: a(:,:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3, n4 !< new dimension
    
    real*8, allocatable :: temp(:,:,:,:)
    integer :: nold(4)
    
    if (.not.allocated(a)) then
       allocate(a(1:n1,1:n2,1:n3,1:n4))
       return
    end if
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    nold(4) = size(a,4)
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3 .and.&
        nold(4) == n4) return
    allocate(temp(n1,n2,n3,n4))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)))
    call move_alloc(temp,a)

  end subroutine realloc4r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc5r(a,n1,n2,n3,n4,n5)
    real*8, intent(inout), allocatable :: a(:,:,:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3, n4, n5 !< new dimension
    
    real*8, allocatable :: temp(:,:,:,:,:)
    integer :: nold(5)
    
    if (.not.allocated(a)) then
       allocate(a(1:n1,1:n2,1:n3,1:n4,1:n5))
       return
    end if
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    nold(4) = size(a,4)
    nold(5) = size(a,5)
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3 .and.&
        nold(4) == n4 .and. nold(5) == n5) return
    allocate(temp(n1,n2,n3,n4,n5))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)),1:min(n5,nold(5))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)),1:min(n5,nold(5)))
    call move_alloc(temp,a)

  end subroutine realloc5r

  !> Adapt the size of an allocatable 1D integer array
  subroutine realloc1i(a,nnew)
    integer, intent(inout), allocatable :: a(:) !< Input array, integer, 1D
    integer, intent(in) :: nnew !< New dimension
    
    integer, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1i

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc2i(a,n1,n2)
    integer, intent(inout), allocatable :: a(:,:) !< Input array, integer, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    integer, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) then
       allocate(a(1:n1,1:n2))
       return
    end if
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    if (nold(1) == n1 .and. nold(2) == n2) return
    allocate(temp(n1,n2))
    
    temp = 0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2i

  !> Adapt the size of an allocatable 1D character array
  subroutine realloc1c(a,nnew)
    character*(*), intent(inout), allocatable :: a(:) !< Input array, character, 1D
    integer, intent(in) :: nnew !< New dimension
    
    character*(len(a)), allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1c

  !> Adapt the size of an allocatable 1D complex*8 array
  subroutine realloc1cmplx4(a,nnew)
    complex*8, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    complex*8, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1cmplx4

  !> Adapt the size of an allocatable 2D complex*8 array
  subroutine realloc2cmplx4(a,n1,n2)
    complex*8, intent(inout), allocatable :: a(:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    complex*8, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) then
       allocate(a(1:n1,1:n2))
       return
    end if
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    if (nold(1) == n1 .and. nold(2) == n2) return
    allocate(temp(n1,n2))
    
    temp = cmplx(0,8)
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2cmplx4

  !> Adapt the size of an allocatable 2D complex*8 array
  subroutine realloc4cmplx4(a,n1,n2,n3,n4)
    complex*8, intent(inout), allocatable :: a(:,:,:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2, n3, n4 !< new dimension
    
    complex*8, allocatable :: temp(:,:,:,:)
    integer :: nold(4), i
    
    if (.not.allocated(a)) then
       allocate(a(1:n1,1:n2,1:n3,1:n4))
       return
    end if
    do i = 1, 4
       nold(i) = size(a,i)
    end do
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3 .and.&
        nold(4) == n4) return
    allocate(temp(n1,n2,n3,n4))
    
    temp = cmplx(0,8)
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)))
    call move_alloc(temp,a)

  end subroutine realloc4cmplx4

  !> Adapt the size of an allocatable 1D complex*16 array
  subroutine realloc1cmplx8(a,nnew)
    complex*16, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    complex*16, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) then
       allocate(a(1:nnew))
       return
    end if
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1cmplx8

  !> Adapt the size of an allocatable 5D complex*8 array
  subroutine realloc5cmplx8(a,n1,n2,n3,n4,n5)
    complex*8, intent(inout), allocatable :: a(:,:,:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3, n4, n5 !< new dimension
    
    complex*8, allocatable :: temp(:,:,:,:,:)
    integer :: nold(5)
    
    if (.not.allocated(a)) then
       allocate(a(1:n1,1:n2,1:n3,1:n4,1:n5))
       return
    end if
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    nold(4) = size(a,4)
    nold(5) = size(a,5)
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3 .and.&
        nold(4) == n4 .and. nold(5) == n5) return
    allocate(temp(n1,n2,n3,n4,n5))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)),1:min(n5,nold(5))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)),1:min(n5,nold(5)))
    call move_alloc(temp,a)

  end subroutine realloc5cmplx8

end module types
