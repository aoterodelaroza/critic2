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

! Methods for the fragment class.
module fragmentmod
  use types
  implicit none

  private
  public :: fragment_init
  public :: fragment_merge_array
  public :: fragment_cmass

contains

  !> Initialize a fragment
  subroutine fragment_init(fr)
    type(fragment), intent(inout) :: fr
    
    if (allocated(fr%at)) deallocate(fr%at)
    allocate(fr%at(1))
    fr%nat = 0
    
  end subroutine fragment_init

  !> Merge two or more fragments, delete repeated atoms
  function fragment_merge_array(fra) result(fr)

    type(fragment), intent(in) :: fra(:)
    type(fragment) :: fr
    
    real*8, parameter :: eps = 1d-10

    integer :: i, j, k, nat0, nat1
    real*8 :: d2, x(3)
    logical :: found

    call fragment_init(fr)
    do i = 1, size(fra)
       nat0 = fr%nat
       nat1 = fr%nat + fra(i)%nat
       if (nat1 > size(fr%at)) call realloc(fr%at,2*nat1)
       do j = 1, fra(i)%nat
          found = .false.
          do k = 1, nat0
             x = abs(fra(i)%at(j)%r - fr%at(k)%r)
             found = all(x < eps)
             if (found) exit
          end do
          if (.not.found) then
             nat0 = nat0 + 1
             fr%at(nat0) = fra(i)%at(j)
          end if
       end do
       fr%nat = nat0
    end do
    call realloc(fr%at,fr%nat)
    
  end function fragment_merge_array

  !> Returns the center of mass (in cartesian coordinates).  If
  !> weight0 is false, then all atoms have the same weight.
  function fragment_cmass(fr,weight0) result (x)
    use param

    type(fragment), intent(in) :: fr
    logical, intent(in), optional :: weight0
    real*8 :: x(3)

    integer :: i
    logical :: weight
    real*8 :: sum

    weight = .true.
    if (present(weight0)) weight = weight0

    x = 0d0
    sum = 0d0
    if (weight) then
       do i = 1, fr%nat
          x = x + atmass(fr%at(i)%z) * fr%at(i)%r
          sum = sum + atmass(fr%at(i)%z)
       end do
    else
       do i = 1, fr%nat
          x = x + fr%at(i)%r
          sum = sum + 1d0
       end do
    end if
    x = x / sum

  end function fragment_cmass

end module fragmentmod
