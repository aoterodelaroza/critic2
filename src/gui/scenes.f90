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
     integer(c_int) :: order ! order of the bond (0=dashed,1=single,2=double,3=triple)
     real(c_float) :: border ! border size
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

  !> Draw style for atoms
  type draw_style_atom
     logical :: isinit = .false. ! whether the style is intialized
     integer :: type ! atom style type: 0=species,1=nneq,2=cell
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     real(c_float), allocatable :: rgb(:,:) ! color (3,ntype)
     real(c_float), allocatable :: rad(:) ! radius (ntype)
   contains
     procedure :: reset => reset_atom_style
  end type draw_style_atom
  public :: draw_style_atom

  !> Draw style for molecules
  type draw_style_molecule
     logical :: isinit = .false. ! whether the style is intialized
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     real(c_float), allocatable :: tint_rgb(:,:) ! tint color (3,ntype)
     real(c_float), allocatable :: scale_rad(:) ! scale radius (ntype)
   contains
     procedure :: reset => reset_molecule_style
  end type draw_style_molecule
  public :: draw_style_molecule

  !> Draw style for bonds
  type draw_style_bond
     logical :: isinit = .false. ! whether the style is intialized
     logical :: isdef = .true. ! whether this is using the system's neighbor star
     ! temporary storage for the edit object/bonds tab, global options
     integer(c_int) :: distancetype_g ! selector for distance type (0=factor,1=range)
     real(c_float) :: dmin_g, dmax_g ! distance limits (angstrom)
     real(c_float) :: bfmin_g, bfmax_g ! bondfactor limits
     integer(c_int) :: radtype_g(2) ! radii type for min and max (0=covalent,1=vdw)
     integer(c_int) :: style_g ! bond style (0=single color, 1=two colors)
     real(c_float) :: rad_g ! radius
     real(c_float) :: border_g ! bond border
     real(c_float) :: rgb_g(3) ! color
     integer(c_int) :: order_g ! order (0=dashed,1=single,2=double,etc.)
     integer(c_int) :: imol_g ! molecular connections (0=any,1=intramol,2=intermol)
     logical :: bothends_g ! if true, both atoms need to be drawn to draw the bond
     logical, allocatable :: shown_g(:,:) ! by-species bond shown flags (nspc,nspc)
     ! the bond information
     type(neighstar), allocatable :: nstar(:) ! the neighbor star
   contains
     procedure :: reset => reset_bond_style
     procedure :: generate_neighstars_from_globals
  end type draw_style_bond
  public :: draw_style_bond

  ! types of representations
  integer, parameter, public :: reptype_none = 0
  integer, parameter, public :: reptype_atoms = 1
  integer, parameter, public :: reptype_unitcell = 2
  integer, parameter, public :: reptype_NUM = 2

  ! representation flavors
  integer, parameter, public :: repflavor_unknown = 0
  integer, parameter, public :: repflavor_atoms_basic = 1
  integer, parameter, public :: repflavor_atoms_vdwcontacts = 2
  integer, parameter, public :: repflavor_atoms_hbonds = 3
  integer, parameter, public :: repflavor_unitcell_basic = 4
  integer, parameter, public :: repflavor_NUM = 4

  !> Representation: objects to draw on the scene
  type representation
     ! main variables
     logical :: isinit = .false. ! whether the representation has been initialized
     logical :: shown = .false. ! true if the representation is currently shown
     integer :: type = reptype_none ! type of representation (atoms, cell,...)
     integer :: flavor = repflavor_unknown ! flavor of the representation
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
     logical :: atoms_display = .true. ! whether to draw the atoms
     logical :: bonds_display = .true. ! whether to draw the bonds
     logical :: labels_display = .true. ! whether to draw the labels
     logical :: border = .true. ! draw atoms at the border of the unit cell
     logical :: onemotif = .false. ! draw connected molecules
     integer(c_int) :: atom_radii_reset_type = 0 ! option to reset radii: 0=covalent, 1=vdw
     real(c_float) :: atom_radii_reset_scale = 0.7_c_float ! reset radii, scale factor
     integer(c_int) :: atom_color_reset_type = 0 ! option to reset colors: 0=jmlcol, 1=jmlcol2
     type(draw_style_atom) :: atom_style ! atom styles
     type(draw_style_molecule) :: mol_style ! molecule styles
     type(draw_style_bond) :: bond_style ! bond styles
     integer(c_int) :: label_style ! 0=atom-symbol, 1=atom-name, 2=cel-atom, 3=cel-atom+lvec, 4=neq-atom, 5=spc, 6=Z, 7=mol, 8=wyckoff
     real(c_float) :: label_scale ! scale for the labels
     real(c_float) :: label_rgb(3) ! color of the labels
     logical :: label_const_size ! whether labels scale with objects or are constant size
     logical :: label_exclude_h ! whether to exclude hydrogen labels
     ! unit cell
     logical :: uc_inner ! unit cell, display inner cylinders
     logical :: uc_coloraxes ! unit cell, color the axes (x=red,y=green,z=blue)
     real(c_float) :: uc_radius ! unit cell cylinder radius
     real(c_float) :: uc_radiusinner ! unit cell cylinder radius (inner)
     real(c_float) :: uc_rgb(3) ! unit cell cylinder colors
     real(c_float) :: uc_innersteplen ! number of subdivisions for the inner sticks
     logical :: uc_innerstipple ! stippled lines for the inner lines
   contains
     procedure :: init => representation_init
     procedure :: end => representation_end
     procedure :: update => update_structure
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
       integer, intent(in) :: idx(5)
     end subroutine select_atom
     ! draw_style_atom
     module subroutine reset_atom_style(d,isys,itype)
       class(draw_style_atom), intent(inout), target :: d
       integer, intent(in), value :: isys, itype
     end subroutine reset_atom_style
     ! draw_style_molecule
     module subroutine reset_molecule_style(d,isys)
       class(draw_style_molecule), intent(inout), target :: d
       integer, intent(in), value :: isys
     end subroutine reset_molecule_style
     ! draw_style_bond
     module subroutine reset_bond_style(d,isys,flavor)
       class(draw_style_bond), intent(inout), target :: d
       integer, intent(in), value :: isys
       integer, intent(in) :: flavor
     end subroutine reset_bond_style
     module subroutine generate_neighstars_from_globals(d,isys)
       class(draw_style_bond), intent(inout), target :: d
       integer, intent(in) :: isys
     end subroutine generate_neighstars_from_globals
     ! representation
     module subroutine representation_init(r,sc,isys,irep,itype,style,flavor)
       class(representation), intent(inout), target :: r
       type(scene), intent(inout), target :: sc
       integer, intent(in) :: isys
       integer, intent(in) :: irep
       integer, intent(in) :: itype
       integer, intent(in) :: style
       integer, intent(in) :: flavor
     end subroutine representation_init
     module subroutine representation_end(r)
       class(representation), intent(inout), target :: r
     end subroutine representation_end
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
     module subroutine update_structure(r)
       class(representation), intent(inout), target :: r
     end subroutine update_structure
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
