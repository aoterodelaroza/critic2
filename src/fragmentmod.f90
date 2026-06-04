! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
  use types, only: anyatom, species, thread_info
  use molsymmod, only: point_group
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
     real*8 :: m_x2c(3,3) = 0d0 !< crystallographic -> Cartesian matrix (saved at build time)
     ! standard orientation (canonical body frame), computed lazily; see fragment_compute_std
     logical :: axes_computed = .false. !< whether the standard frame has been computed/cached
     logical :: isatom = .false.   !< whether this is... a single atom
     logical :: islinear = .false. !< ... a linear molecule
     logical :: isplanar = .false. !< ... a planar molecule
     real*8 :: xcm(3) = 0d0 !< center of mass (Cartesian)
     real*8 :: inertia(3) = 0d0 !< principal moments of inertia (ascending)
     real*8 :: m_std(3,3) = 0d0 !< standard orientation (det=+1); columns are the principal axes in Cartesian coordinates; canonical coords: matmul(transpose(m_std),at(k)%r-xcm)
     real*8 :: quat_std(4) = (/1d0,0d0,0d0,0d0/) !< unit quaternion (w,x,y,z) of the standard orientation = mat2quat(m_std)
     real*8 :: euler_std(3) = 0d0 !< ZYZ Euler angles (alpha,beta,gamma, radians) of the standard orientation = mat2euler(m_std)
     ! molecular point group, computed lazily; see fragment_pgsymbol
     logical :: pg_computed = .false. !< whether the point group has been computed/cached
     type(point_group) :: pg !< molecular point group of the fragment
   contains
     procedure :: init => fragment_init ! initialize a fragment
     procedure :: build => fragment_build ! build a fragment from atomic data
     procedure :: build_atom => fragment_build_atom ! build a fragment from a single atom
     procedure :: translate => fragment_translate ! translate the fragment by a lattice vector
     procedure :: translate_to_main_cell => fragment_translate_to_main_cell ! move the fragment to the main cell
     procedure :: standard_axes => fragment_standard_axes ! accessor for the standard (canonical) frame
     procedure :: compute_std => fragment_compute_std ! compute and cache the standard frame
     procedure :: pgsymbol => fragment_pgsymbol ! accessor for the point-group symbol (lazy)
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
     module subroutine fragment_build(fr,nspc,spc,nat,xyz,icrd,is,idx,cidx,m_x2c,&
        lvec,discrete,nlvec,mollvec)
       class(fragment), intent(inout) :: fr
       integer, intent(in) :: nspc
       type(species), intent(in) :: spc(nspc)
       integer, intent(in) :: nat
       real*8, intent(in) :: xyz(3,nat)
       integer, intent(in) :: icrd
       integer, intent(in) :: is(nat)
       integer, intent(in) :: idx(nat)
       integer, intent(in) :: cidx(nat)
       real*8, intent(in) :: m_x2c(3,3)
       integer, intent(in), optional :: lvec(3,nat)
       logical, intent(in), optional :: discrete
       integer, intent(in), optional :: nlvec
       integer, intent(in), optional :: mollvec(:,:)
     end subroutine fragment_build
     module subroutine fragment_build_atom(fr,spc,xyz,icrd,idx,cidx,m_x2c)
       class(fragment), intent(inout) :: fr
       type(species), intent(in) :: spc
       real*8, intent(in) :: xyz(3)
       integer, intent(in) :: icrd
       integer, intent(in) :: idx
       integer, intent(in) :: cidx
       real*8, intent(in) :: m_x2c(3,3)
     end subroutine fragment_build_atom
     module subroutine fragment_translate(fr,lvec)
       class(fragment), intent(inout) :: fr
       integer, intent(in) :: lvec(3)
     end subroutine fragment_translate
     module subroutine fragment_translate_to_main_cell(fr)
       class(fragment), intent(inout) :: fr
     end subroutine fragment_translate_to_main_cell
     module subroutine fragment_standard_axes(fr,m_std,inertia,xcm,isatom,islinear,isplanar,quat,euler)
       class(fragment), intent(inout) :: fr
       real*8, intent(out), optional :: m_std(3,3)
       real*8, intent(out), optional :: inertia(3)
       real*8, intent(out), optional :: xcm(3)
       logical, intent(out), optional :: isatom
       logical, intent(out), optional :: islinear
       logical, intent(out), optional :: isplanar
       real*8, intent(out), optional :: quat(4)
       real*8, intent(out), optional :: euler(3)
     end subroutine fragment_standard_axes
     module subroutine fragment_compute_std(fr)
       class(fragment), intent(inout) :: fr
     end subroutine fragment_compute_std
     module function fragment_pgsymbol(fr) result(symbol)
       class(fragment), intent(inout) :: fr
       character(len=:), allocatable :: symbol
     end function fragment_pgsymbol
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
     module subroutine writexyz(fr,file,usenames,ti)
       class(fragment), intent(in) :: fr
       character*(*), intent(in) :: file
       logical, intent(in) :: usenames
       type(thread_info), intent(in), optional :: ti
     end subroutine writexyz
     module subroutine writecml(fr,file,r,luout,ti)
       class(fragment), intent(in) :: fr
       character*(*), intent(in) :: file
       real*8, intent(in), optional :: r(3,3)
       integer, intent(out), optional :: luout
       type(thread_info), intent(in), optional :: ti
     end subroutine writecml
     module subroutine writegjf(fr,file,ti)
       class(fragment), intent(in) :: fr
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine writegjf
     module subroutine realloc_fragment(a,nnew)
       type(fragment), intent(inout), allocatable :: a(:)
       integer, intent(in) :: nnew
     end subroutine realloc_fragment
  end interface

end module fragmentmod
