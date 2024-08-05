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

  !> Zoom minimum and maximum
  real(c_float), parameter, public :: min_zoom = 1._c_float
  real(c_float), parameter, public :: max_zoom = 100._c_float

  !> Animation default parameters
  real(c_float), parameter, public :: anim_speed_default = 4._c_float
  real(c_float), parameter, public :: anim_amplitude_default = 1._c_float

  !> spheres for the draw list
  type dl_sphere
     real(c_float) :: x(3) ! position
     real(c_float) :: r ! radius
     real(c_float) :: rgb(3) ! color
     integer(c_int) :: idx(4) ! atom ID (complete atom list) + lattice vector
     complex(c_float_complex) :: xdelta(3) ! delta-vector for vibration animations
  end type dl_sphere

  !> cylinders for the draw list
  type dl_cylinder
     real(c_float) :: x1(3) ! one end of the cylinder
     real(c_float) :: x2(3) ! other end of the cylinder
     real(c_float) :: r ! radius
     real(c_float) :: rgb(3) ! color
     complex(c_float_complex) :: x1delta(3) ! delta-vector for vibration animations (end 1)
     complex(c_float_complex) :: x2delta(3) ! delta-vector for vibration animations (end 2)
  end type dl_cylinder

  !> strings for the draw list
  type dl_string
     real(c_float) :: x(3) ! position
     real(c_float) :: r ! radius
     real(c_float) :: rgb(3) ! color
     real(c_float) :: scale ! scale (1.0 = radius)
     character(len=:), allocatable :: str ! string
     complex(c_float_complex) :: xdelta(3) ! delta-vector for vibration animations
  end type dl_string

  !> draw style for atoms
  type draw_style_atom
     logical(c_bool) :: shown
     real(c_float) :: rgb(3) ! color
     real(c_float) :: rad ! radius
  end type draw_style_atom

  integer, parameter, public :: reptype_none = 0
  integer, parameter, public :: reptype_atoms = 1
  integer, parameter, public :: reptype_unitcell = 2
  integer, parameter, public :: reptype_NUM = 2

  !> Representation: objects to draw on the scene
  type representation
     ! main variables
     logical :: isinit = .false. ! whether the representation has been initialized
     logical(c_bool) :: shown = .false. ! true if the representation is currently shown
     integer :: type = reptype_none ! type of representation (atoms, cell,...)
     integer :: id ! system ID
     integer :: idrep ! representation ID
     integer :: idwin = 0 ! edit representation window ID
     integer :: iord = 0 ! representation order integer in menu
     character(kind=c_char,len=:), allocatable :: name ! name of the representation
     ! global parameters
     integer(c_int) :: pertype = 1 ! periodicity control: 0=none, 1=auto, 2=manual
     integer(c_int) :: ncell(3) ! number of unit cells drawn
     real(c_float) :: origin(3) ! unit cell, origin shift
     real(c_float) :: tshift(3) ! origin of the unit cell display region
     ! atoms, bonds, labels
     character(kind=c_char,len=:), allocatable :: filter ! filter for the representation
     character(kind=c_char,len=:), allocatable :: errfilter ! filter error
     logical(c_bool) :: atoms_display = .true. ! whether to draw the atoms
     logical(c_bool) :: bonds_display = .true. ! whether to draw the bonds
     logical(c_bool) :: labels_display = .true. ! whether to draw the labels
     logical(c_bool) :: border = .true. ! draw atoms at the border of the unit cell
     logical(c_bool) :: onemotif = .false. ! draw connected molecules
     integer(c_int) :: atom_style_type = 0 ! atom style type: 0=species,1=nneq,2=cell
     integer :: natom_style = 0 ! number of atom styles
     integer(c_int) :: atom_radii_reset_type = 0 ! option to reset radii: 0=covalent, 1=vdw
     real(c_float) :: atom_radii_reset_scale = 0.7_c_float ! reset radii, scale factor
     integer(c_int) :: atom_color_reset_type = 0 ! option to reset colors: 0=jmlcol, 1=jmlcol2
     integer(c_int) :: bond_color_style ! bond color style: 0=single color, 1=half-and-half
     real(c_float) :: bond_rgb(3) ! bond color (single color style)
     real(c_float) :: bond_rad ! bond radius
     type(draw_style_atom), allocatable :: atom_style(:) ! atom styles
     integer(c_int) :: label_style ! 0=atom-symbol, 1=atom-name, 2=cel-atom, 3=cel-atom+lvec, 4=neq-atom, 5=spc, 6=Z, 7=mol, 8=wyckoff
     real(c_float) :: label_scale ! scale for the labels
     real(c_float) :: label_rgb(3) ! color of the labels
     logical(c_bool) :: label_const_size ! whether labels scale with objects or are constant size
     logical(c_bool) :: label_exclude_h ! whether to exclude hydrogen labels
     ! unit cell
     logical(c_bool) :: uc_inner ! unit cell, display inner cylinders
     logical(c_bool) :: uc_coloraxes ! unit cell, color the axes (x=red,y=green,z=blue)
     real(c_float) :: uc_radius ! unit cell cylinder radius
     real(c_float) :: uc_radiusinner ! unit cell cylinder radius (inner)
     real(c_float) :: uc_rgb(3) ! unit cell cylinder colors
     real(c_float) :: uc_innersteplen ! number of subdivisions for the inner sticks
     logical(c_bool) :: uc_innerstipple ! stippled lines for the inner lines
   contains
     procedure :: init => representation_init
     procedure :: end => representation_end
     procedure :: reset_atom_style
     procedure :: add_draw_elements
  end type representation
  public :: representation

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
     real(c_float) :: atomborder ! atom border (simple shader)
     real(c_float) :: bordercolor(3) ! border color (simple shader)
     ! scene transformation matrices
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
     integer :: lockedcam = 0 ! 0=no, n=locking group n (same as first system in the lock group)
     ! list of representations
     integer :: nrep = 0 ! number of representation
     type(representation), allocatable :: rep(:) ! representations
     integer, allocatable :: icount(:) ! last rep counter, for unique names
     ! measure atom sets
     integer :: nmsel
     integer :: msel(4,4)
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
   contains
     procedure :: init => scene_init
     procedure :: end => scene_end
     procedure :: reset => scene_reset
     procedure :: reset_animation => scene_reset_animation
     procedure :: build_lists => scene_build_lists
     procedure :: render => scene_render
     procedure :: renderpick => scene_render_pick
     procedure :: set_style_defaults => scene_set_style_defaults
     procedure :: copy_cam => scene_copy_cam
     procedure :: representation_menu
     procedure :: get_new_representation_id
     procedure :: update_projection_matrix
     procedure :: update_view_matrix
     procedure :: align_view_axis
     procedure :: select_atom
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
     recursive module subroutine scene_copy_cam(s,si,idx)
       class(scene), intent(inout), target :: s
       type(scene), intent(in), target, optional :: si
       integer, intent(in), optional :: idx
     end subroutine scene_copy_cam
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
     module subroutine select_atom(s,idx)
       class(scene), intent(inout), target :: s
       integer, intent(in) :: idx(4)
     end subroutine select_atom
     ! representation
     module subroutine representation_init(r,sc,isys,irep,itype,style)
       class(representation), intent(inout), target :: r
       type(scene), intent(inout), target :: sc
       integer, intent(in) :: isys
       integer, intent(in) :: irep
       integer, intent(in) :: itype
       integer, intent(in) :: style
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
     module subroutine add_draw_elements(r,nc,nsph,drawlist_sph,ncyl,drawlist_cyl,&
        ncylflat,drawlist_cylflat,nstring,drawlist_string,doanim,iqpt,ifreq)
       class(representation), intent(inout), target :: r
       integer, intent(in) :: nc(3)
       integer, intent(inout) :: nsph
       type(dl_sphere), intent(inout), allocatable :: drawlist_sph(:)
       integer, intent(inout) :: ncyl
       type(dl_cylinder), intent(inout), allocatable :: drawlist_cyl(:)
       integer, intent(inout) :: ncylflat
       type(dl_cylinder), intent(inout), allocatable :: drawlist_cylflat(:)
       integer, intent(inout) :: nstring
       type(dl_string), intent(inout), allocatable :: drawlist_string(:)
       logical, intent(in) :: doanim
       integer, intent(in) :: iqpt, ifreq
     end subroutine add_draw_elements
     module subroutine draw_text_direct(str,x0,siz,color,centered)
       character(len=*), intent(in) :: str
       real(c_float), intent(in) :: x0(2)
       real(c_float), intent(in) :: siz
       real(c_float), intent(in) :: color(3)
       logical, intent(in), optional :: centered
     end subroutine draw_text_direct
  end interface

end module scenes
