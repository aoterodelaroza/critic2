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
  use representations, only: representation, dl_sphere, dl_cylinder, dl_string
  use types, only: neighstar
  implicit none

  private

  !> Zoom minimum and maximum
  real(c_float), parameter, public :: min_zoom = 1._c_float
  real(c_float), parameter, public :: max_zoom = 100._c_float

  !> Animation default parameters
  real(c_float), parameter, public :: anim_speed_default = 4._c_float
  real(c_float), parameter, public :: anim_amplitude_default = 5._c_float
  real(c_float), parameter, public :: anim_amplitude_max = 15._c_float
  real(c_float), parameter, public :: anim_speed_max = 50._c_float

  ! scene style
  integer(c_int), parameter, public :: style_simple = 0
  integer(c_int), parameter, public :: style_phong = 1

  !> Scene: objects from the system to be drawn and plot settings
  type scene
     integer :: isinit = 0 ! 0=uninitialized, 1=initialized but not built, 2=init and built
     integer :: id ! system ID
     integer, allocatable :: iord(:) ! the representation order
     logical :: forcesort = .false. ! force sort the representations
     logical :: forceresetcam = .false. ! force reset of the camera
     logical :: forcebuildlists ! force rebuild of lists
     real*8 :: timelastrender = 0d0 ! time when the view was last rendered
     real*8 :: timelastbuild = 0d0 ! time of the last build
     real*8 :: timelastcamchange = 0d0 ! time the camera was last changed
     real*8 :: timerefanimation = 0d0 ! reference time for the animation
     real(c_float) :: scenerad = 1d0 ! scene radius
     real(c_float) :: scenecenter(3) ! scene center (world coords)
     real(c_float) :: scenexmin(3) ! scene xmin (world coords)
     real(c_float) :: scenexmax(3) ! scene xmax (world coords)
     integer(c_int) :: nc(3) ! number of unit cells drawn (global +/-)
     ! object resolutions
     integer(c_int) :: atom_res ! atom resolution
     integer(c_int) :: bond_res ! bond resolution
     integer(c_int) :: uc_res ! unit cell resolution
     ! scene/shader settings
     integer(c_int) :: style ! scene style (0=simple,1=phong)
     real(c_float) :: bgcolor(3) ! background color
     real(c_float) :: lightpos(3) ! light position
     real(c_float) :: lightcolor(3) ! light color
     real(c_float) :: ambient ! ambient light coefficent
     real(c_float) :: diffuse ! diffuse light coefficent
     real(c_float) :: specular ! specular light coefficent
     integer(c_int) :: shininess ! shininess parameter
     real(c_float) :: bordercolor(3) ! border color (simple shader)
     ! scene transformation matrices and camera options
     logical :: iscaminit = .false. ! true if the camera has been initialized
     real(c_float) :: camresetdist ! camera reset distance
     real(c_float) :: camratio ! window ratio for the camera
     real(c_float) :: ortho_fov ! orthographic field of view
     real(c_float) :: persp_fov ! perspective field of view
     real(c_float) :: campos(3) ! position of the camera (tworld coords)
     real(c_float) :: camfront(3) ! camera front vec (tworld coords)
     real(c_float) :: camup(3) ! camera up vec (tworld coords)
     real(c_float) :: world(4,4) ! world transform matrix
     real(c_float) :: view(4,4) ! view transform matrix
     real(c_float) :: projection(4,4) ! projection transform matrix
     integer :: lockedcam = 0 ! 0=no, -1=global locking group, n=locking group n (same as first system in the lock group)
     ! list of representations
     integer :: nrep = 0 ! number of representation
     type(representation), allocatable :: rep(:) ! representations
     integer, allocatable :: icount(:) ! last rep counter, for unique names
     ! measure atom sets
     integer :: nmsel
     integer :: msel(5,4) ! 1 is atom cell ID, 2:4 is lattice vector, 5 is sphere ID
     ! draw lists
     integer :: nsph ! number of spheres
     type(dl_sphere), allocatable :: drawlist_sph(:) ! sphere draw list
     integer :: ncyl ! number of cylinders
     type(dl_cylinder), allocatable :: drawlist_cyl(:) ! cylinder draw list
     integer :: ncylflat ! number of flat cylinders
     type(dl_cylinder), allocatable :: drawlist_cylflat(:) ! flat cylinder draw list
     integer :: nstring ! number of strings
     type(dl_string), allocatable :: drawlist_string(:) ! flat cylinder draw list
     ! vibration animation
     integer :: animation = 0 ! animate the scene? 0=off, 1=manual, 2=automatic
     integer :: iqpt_selected = 0 ! selected q-point in the vibrations window
     integer :: ifreq_selected = 0 ! selected frequency in the vibrations window
     real(c_float) :: anim_speed = anim_speed_default ! animation speed
     real(c_float) :: anim_amplitude = anim_amplitude_default ! animation amplitude
     real(c_float) :: anim_phase = 0._c_float ! animation phase (manual)
   contains
     procedure :: init => scene_init
     procedure :: end => scene_end
     procedure :: reset => scene_reset
     procedure :: reset_animation => scene_reset_animation
     procedure :: reset_atom_colors => scene_reset_atom_colors
     procedure :: build_lists => scene_build_lists
     procedure :: render => scene_render
     procedure :: renderpick => scene_render_pick
     procedure :: set_style_defaults => scene_set_style_defaults
     procedure :: cam_copy => scene_cam_copy
     procedure :: cam_zoom => scene_cam_zoom
     procedure :: cam_move => scene_cam_move
     procedure :: cam_rotate => scene_cam_rotate
     procedure :: representation_menu
     procedure :: get_new_representation_id
     procedure :: update_projection_matrix
     procedure :: update_view_matrix
     procedure :: align_view_axis
     procedure :: select_atom
     procedure :: add_representation
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
     module subroutine scene_reset_animation(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_reset_animation
     module subroutine scene_reset_atom_colors(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_reset_atom_colors
     module subroutine scene_build_lists(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_build_lists
     module subroutine scene_render(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_render
     module subroutine scene_render_pick(s)
       class(scene), intent(inout), target :: s
     end subroutine scene_render_pick
     module subroutine scene_set_style_defaults(s,style)
       class(scene), intent(inout), target :: s
       integer(c_int), intent(in), optional :: style
     end subroutine scene_set_style_defaults
     module subroutine scene_cam_copy(s,si)
       class(scene), intent(inout), target :: s
       type(scene), intent(in), target :: si
     end subroutine scene_cam_copy
     module subroutine scene_cam_zoom(s,ratio)
       class(scene), intent(inout), target :: s
       real(c_float), intent(in) :: ratio
     end subroutine scene_cam_zoom
     module subroutine scene_cam_move(s,xc)
       class(scene), intent(inout), target :: s
       real(c_float), intent(in) :: xc(3)
     end subroutine scene_cam_move
     module subroutine scene_cam_rotate(s,axis,ang)
       class(scene), intent(inout), target :: s
       real(c_float), intent(in) :: axis(3)
       real(c_float), intent(in) :: ang
     end subroutine scene_cam_rotate
     module function representation_menu(s,idparent) result(changed)
       class(scene), intent(inout), target :: s
       integer(c_int), intent(in) :: idparent
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
     module subroutine select_atom(s,idx)
       class(scene), intent(inout), target :: s
       integer, intent(in) :: idx(5)
     end subroutine select_atom
     module subroutine add_representation(s,itype,flavor)
       class(scene), intent(inout), target :: s
       integer, intent(in) :: itype
       integer, intent(in) :: flavor
     end subroutine add_representation
  end interface

end module scenes
