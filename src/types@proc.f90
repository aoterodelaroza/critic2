! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (types) proc
  implicit none

contains

  !> Clear of values and flags of a scalar value type.
  module subroutine scalar_value_clear(s)
    class(scalar_value), intent(inout) :: s

    s%f = 0d0
    s%fval = 0d0
    s%gf = 0d0
    s%hf = 0d0
    s%gfmod = 0d0
    s%gfmodval = 0d0
    s%del2f = 0d0
    s%del2fval = 0d0
    s%gkin = 0d0
    s%stress = 0d0
    s%vir = 0d0
    s%hfevec = 0d0
    s%hfeval = 0d0
    s%r = 0
    s%s = 0
    s%fup = 0d0
    s%fdn = 0d0
    s%fspin = 0d0
    s%fspc = 0d0
    s%isnuc = .false.
    s%avail_der1 = .false.
    s%avail_der2 = .false.
    s%avail_gkin = .false.
    s%avail_stress = .false.
    s%avail_vir = .false.
    s%avail_spin = .false.

  end subroutine scalar_value_clear

  !> Adapt the size of an allocatable 1D type(pointpropable) array
  module subroutine realloc_pointpropable(a,nnew)
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
  module subroutine realloc_integrable(a,nnew)
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

  !> Adapt the size of an allocatable 1D type(species) array
  module subroutine realloc_species(a,nnew)
    type(species), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(species), allocatable :: temp(:)
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
    end do

  end subroutine realloc_species

  !> Adapt the size of an allocatable 1D type(atom) array
  module subroutine realloc_basicatom(a,nnew)
    type(basicatom), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(basicatom), allocatable :: temp(:)
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

  end subroutine realloc_basicatom

  !> Adapt the size of an allocatable 1D type(atom) array
  module subroutine realloc_neqatom(a,nnew)
    type(neqatom), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(neqatom), allocatable :: temp(:)
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

  end subroutine realloc_neqatom

  !> Adapt the size of an allocatable 1D type(celatom) array
  module subroutine realloc_celatom(a,nnew)
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
  module subroutine realloc_anyatom(a,nnew)
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
  module subroutine realloc_cp(a,nnew)
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
  module subroutine realloc_gpathp(a,nnew)
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
  module subroutine realloc1l(a,nnew)
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
  module subroutine realloc1r(a,nnew)
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
  module subroutine realloc2r(a,n1,n2)
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
  module subroutine realloc3r(a,n1,n2,n3)
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
  module subroutine realloc4r(a,n1,n2,n3,n4)
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
  module subroutine realloc5r(a,n1,n2,n3,n4,n5)
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
  module subroutine realloc1i(a,nnew)
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
  module subroutine realloc2i(a,n1,n2)
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
  module subroutine realloc1c(a,nnew)
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
  module subroutine realloc1cmplx4(a,nnew)
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
  module subroutine realloc2cmplx4(a,n1,n2)
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
  module subroutine realloc4cmplx4(a,n1,n2,n3,n4)
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
  module subroutine realloc1cmplx8(a,nnew)
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
  module subroutine realloc5cmplx8(a,n1,n2,n3,n4,n5)
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

end submodule proc
