! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Scene object and GL rendering utilities
module scenes
  use iso_c_binding
  implicit none

  private

  integer, parameter :: reptype_none = 0
  integer, parameter :: reptype_atoms = 1
  integer, parameter :: reptype_bonds = 2
  integer, parameter :: reptype_unitcell = 3
  integer, parameter :: reptype_labels = 4

  !> Representation: objects to draw on the scene
  type representation
     logical :: isinit = .false. ! true if the representation has been initialized
     logical(c_bool) :: shown = .false. ! true if the representation is currently shown
     integer :: type = reptype_none ! type of representation (atoms, bonds,...)
     integer :: id ! system ID
     integer :: idwin = 0 ! edit representation window ID
     character(kind=c_char,len=:), allocatable :: name ! name of the representation
     ! global parameters
     integer :: ncell(3) ! number of unit cells drawn (or zero if controlled by the global +/-)
     logical :: border = .true. ! draw atoms at the border of the unit cell
     logical :: onemotif = .false. ! draw connected molecules
     ! atoms
     real*8 :: atom_scale = 1d0 ! atomic radius scaling factor
     ! bonds
     real*8 :: bond_scale = 1d0 ! bond scaling factor
   contains
     procedure :: init => representation_init
     procedure :: end => representation_end
     procedure :: draw => representation_draw
     procedure :: draw_atoms
     procedure :: draw_bonds
     procedure :: draw_unitcell
     procedure :: draw_labels
  end type representation

  !> Scene: objects from the system to be drawn and plot settings
  type scene
     logical :: isinit = .false. ! whether the scene has been initialized
     integer :: id ! system ID
     real(c_float) :: scenerad = 1d0 ! scene radius
     real(c_float) :: scenecenter(3) ! scene center
     ! phong shader settings
     real(c_float) :: lightpos(3) ! light position
     real(c_float) :: lightcolor(3) ! light color
     real(c_float) :: ambient ! ambient light coefficent
     real(c_float) :: diffuse ! diffuse light coefficent
     real(c_float) :: specular ! specular light coefficent
     integer(c_int) :: shininess ! shininess parameter
     ! scene transformation matrices
     real(c_float) :: ortho_fov ! orthographic field of view
     real(c_float) :: persp_fov ! perspective field of view
     real(c_float) :: znear ! position of the near plane
     real(c_float) :: zfar  ! position of the far plane
     real(c_float) :: campos(3) ! position of the camera
     real(c_float) :: world(4,4) ! world transform matrix
     real(c_float) :: view(4,4) ! view transform matrix
     real(c_float) :: projection(4,4) ! projection transform matrix
     ! list of representations
     integer :: nrep = 0
     type(representation), allocatable :: rep(:)
   contains
     procedure :: init => scene_init
     procedure :: end => scene_end
     procedure :: reset => scene_reset
     procedure :: render => scene_render
     procedure :: representation_menu
     procedure :: get_new_representation_id
     procedure :: update_projection_matrix
     procedure :: update_view_matrix
  end type scene
  public :: scene

  ! module procedure interfaces
  interface
     ! scene
     module subroutine scene_init(s,isys)
       class(scene), intent(inout), target :: s
       integer, intent(in) :: isys
     end subroutine scene_init
     module subroutine scene_end(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_end
     module subroutine scene_reset(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_reset
     module subroutine scene_render(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_render
     module function representation_menu(s) result(changed)
       class(scene), intent(inout), target :: s
       logical :: changed
     end function representation_menu
     module function get_new_representation_id(s) result(id)
       class(scene), intent(inout), target :: s
       integer :: id
     end function get_new_representation_id
     module subroutine update_projection_matrix(s)
       class(scene), intent(inout), target :: s
     end subroutine update_projection_matrix
     module subroutine update_view_matrix(s)
       class(scene), intent(inout), target :: s
     end subroutine update_view_matrix
     ! representation
     module subroutine representation_init(r,isys)
       class(representation), intent(inout), target :: r
       integer, intent(in) :: isys
     end subroutine representation_init
     module subroutine representation_end(r)
       class(representation), intent(inout), target :: r
     end subroutine representation_end
     module subroutine representation_draw(r,xmin,xmax)
       class(representation), intent(inout), target :: r
       real*8, optional, intent(inout) :: xmin(3)
       real*8, optional, intent(inout) :: xmax(3)
     end subroutine representation_draw
     module subroutine draw_atoms(r,xmin,xmax)
       class(representation), intent(inout), target :: r
       real*8, optional, intent(inout) :: xmin(3)
       real*8, optional, intent(inout) :: xmax(3)
     end subroutine draw_atoms
     module subroutine draw_bonds(r,xmin,xmax)
       class(representation), intent(inout), target :: r
       real*8, optional, intent(inout) :: xmin(3)
       real*8, optional, intent(inout) :: xmax(3)
     end subroutine draw_bonds
     module subroutine draw_unitcell(r,xmin,xmax)
       class(representation), intent(inout), target :: r
       real*8, optional, intent(inout) :: xmin(3)
       real*8, optional, intent(inout) :: xmax(3)
     end subroutine draw_unitcell
     module subroutine draw_labels(r,xmin,xmax)
       class(representation), intent(inout), target :: r
       real*8, optional, intent(inout) :: xmin(3)
       real*8, optional, intent(inout) :: xmax(3)
     end subroutine draw_labels
  end interface

end module scenes

