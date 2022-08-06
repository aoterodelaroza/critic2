! Copyright (c) 2015 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>,
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

! Molecular fragment class.
module fragmentmod
  use types, only: anyatom, species
  implicit none

  private
  public :: fragment
  public :: realloc_fragment

  !> Type for a fragment of the crystal
  type fragment
     integer :: nat !< Number of atoms in the fragment
     type(anyatom), allocatable :: at(:) !< Atoms in the fragment
     integer :: nspc !< Number of species in the fragment
     type(species), allocatable :: spc(:) !< Species in the fragment
     logical :: discrete = .true. !< is this fragment discrete?
     integer :: nlvec = 0 !< number of connecting lattice vectors for non-discrete fragments
     integer, allocatable :: lvec(:,:) !< connecting lattice vectors
   contains
     procedure :: init => fragment_init ! initialize a fragment
     procedure :: append ! append a fragment to this fragment
     procedure :: merge_array ! merge several fragments
     procedure :: cmass ! calculate center of mass
     procedure :: writexyz ! write an xyz file
     procedure :: writecml ! write a cml file
     procedure :: writegjf ! write a gjf file
  end type fragment

  interface
     module subroutine fragment_init(fr)
       class(fragment), intent(inout) :: fr
     end subroutine fragment_init
     module subroutine merge_array(fr,fra,add)
       class(fragment), intent(inout) :: fr
       type(fragment), intent(in) :: fra(:)
       logical, intent(in), optional :: add
     end subroutine merge_array
     module subroutine append(fr,fra)
       class(fragment), intent(inout) :: fr
       class(fragment), intent(in) :: fra
     end subroutine append
     module function cmass(fr,weight0) result (x)
       class(fragment), intent(in) :: fr
       logical, intent(in), optional :: weight0
       real*8 :: x(3)
     end function cmass
     module subroutine writexyz(fr,file,usenames)
       class(fragment), intent(in) :: fr
       character*(*), intent(in) :: file
       logical, intent(in) :: usenames
     end subroutine writexyz
     module subroutine writecml(fr,file,r,luout)
       class(fragment), intent(in) :: fr
       character*(*), intent(in) :: file
       real*8, intent(in), optional :: r(3,3)
       integer, intent(out), optional :: luout
     end subroutine writecml
     module subroutine writegjf(fr,file)
       class(fragment), intent(in) :: fr
       character*(*), intent(in) :: file
     end subroutine writegjf
     module subroutine realloc_fragment(a,nnew)
       type(fragment), intent(inout), allocatable :: a(:)
       integer, intent(in) :: nnew
     end subroutine realloc_fragment
  end interface

end module fragmentmod
