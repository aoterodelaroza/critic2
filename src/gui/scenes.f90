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

  !> draw style for atoms
  type draw_style_atom
     logical(c_bool) :: shown
     real(c_float) :: rgba(4) ! color
     real(c_float) :: rad ! radius
  end type draw_style_atom

  integer, parameter, public :: reptype_none = 0
  integer, parameter, public :: reptype_atoms = 1
  integer, parameter, public :: reptype_unitcell = 2
  integer, parameter, public :: reptype_NUM = 2

  !> Representation: objects to draw on the scene
  type representation
     logical :: isinit = .false. ! true if the representation has been initialized
     logical(c_bool) :: shown = .false. ! true if the representation is currently shown
     integer :: type = reptype_none ! type of representation (atoms, cell,...)
     integer :: id ! system ID
     integer :: idrep ! representation ID
     integer :: idwin = 0 ! edit representation window ID
     integer :: iord = 0 ! representation order integer in menu
     character(kind=c_char,len=:), allocatable :: name ! name of the representation
     ! global parameters
     character(kind=c_char,len=:), allocatable :: filter ! filter for the representation
     logical :: goodfilter ! true if the filter is not in error
     integer(c_int) :: pertype = 1 ! periodicity control: 0=none, 1=auto, 2=manual
     integer(c_int) :: ncell(3) ! number of unit cells drawn
     ! atoms & bonds
     integer(c_int) :: global_display = 0 ! 0=atoms&bonds, 1=atoms, 2=bonds, 3=none
     logical(c_bool) :: border = .true. ! draw atoms at the border of the unit cell
     logical(c_bool) :: onemotif = .false. ! draw connected molecules
     integer(c_int) :: atom_style_type = 0 ! atom style type: 0=species, 1=nneq, 2=cell
     integer :: natom_style = 0 ! number of atom styles
     integer(c_int) :: atom_radii_reset_type = 0 ! option to reset radii: 0=covalent, 1=vdw
     real(c_float) :: atom_radii_reset_scale = 0.7_c_float ! reset radii, scale factor
     integer(c_int) :: atom_color_reset_type = 0 ! option to reset colors: 0=jmlcol, 1=jmlcol2
     integer(c_int) :: atom_res = 3 ! ball resolution for atoms
     type(draw_style_atom), allocatable :: atom_style(:) ! atom styles
   contains
     procedure :: init => representation_init
     procedure :: end => representation_end
     procedure :: reset_atom_style
     procedure :: draw => representation_draw
     procedure :: draw_atoms
     procedure :: draw_unitcell
  end type representation
  public :: representation

  !> Scene: objects from the system to be drawn and plot settings
  type scene
     logical :: isinit = .false. ! whether the scene has been initialized
     integer :: id ! system ID
     integer, allocatable :: iord(:) ! the representation order
     logical :: forcesort = .false. ! force sort the representations
     real(c_float) :: scenerad = 1d0 ! scene radius
     real(c_float) :: scenecenter(3) ! scene center
     integer(c_int) :: nc(3) ! number of unit cells drawn (global +/-)
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
     integer :: nrep = 0 ! number of representation
     type(representation), allocatable :: rep(:) ! representations
     integer, allocatable :: icount(:) ! last rep counter, for unique names
   contains
     procedure :: init => scene_init
     procedure :: end => scene_end
     procedure :: reset => scene_reset
     procedure :: render => scene_render
     procedure :: representation_menu
     procedure :: get_new_representation_id
     procedure :: update_projection_matrix
     procedure :: update_view_matrix
     procedure :: align_view_axis
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
     module function representation_menu(s,idcaller) result(changed)
       class(scene), intent(inout), target :: s
       integer(c_int), intent(in) :: idcaller
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
     module subroutine align_view_axis(s,iaxis)
       class(scene), intent(inout), target :: s
       integer, intent(in) :: iaxis
     end subroutine align_view_axis
     ! representation
     module subroutine representation_init(r,isys,irep,itype)
       class(representation), intent(inout), target :: r
       integer, intent(in) :: isys
       integer, intent(in) :: irep
       integer, intent(in), optional :: itype
     end subroutine representation_init
     module subroutine representation_end(r)
       class(representation), intent(inout), target :: r
     end subroutine representation_end
     module subroutine reset_atom_style(r)
       class(representation), intent(inout), target :: r
     end subroutine reset_atom_style
     module subroutine representation_draw(r,nc,xmin,xmax)
       class(representation), intent(inout), target :: r
       integer, intent(in) :: nc(3)
       real*8, optional, intent(inout) :: xmin(3)
       real*8, optional, intent(inout) :: xmax(3)
     end subroutine representation_draw
     module subroutine draw_atoms(r,nc,xmin,xmax)
       class(representation), intent(inout), target :: r
       integer, intent(in) :: nc(3)
       real*8, optional, intent(inout) :: xmin(3)
       real*8, optional, intent(inout) :: xmax(3)
     end subroutine draw_atoms
     module subroutine draw_unitcell(r,xmin,xmax)
       class(representation), intent(inout), target :: r
       real*8, optional, intent(inout) :: xmin(3)
       real*8, optional, intent(inout) :: xmax(3)
     end subroutine draw_unitcell
  end interface

end module scenes

