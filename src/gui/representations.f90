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
module representations
  use iso_c_binding
  use types, only: neighstar
  use shapes, only: dl_sphere, dl_cylinder, dl_cylinder_giz, dl_string, dl_string_giz,&
     dl_plane, dl_triangle, scene_objects, dl_append
  use param, only: bohrtoa, eye, maxzat0, atmcov0, mlen
  use global, only: bondfactor_def, bonddelta_def
  implicit none

  private

  ! default parameters for the representations (all distances in bohr)
  !--> atoms
  real*8, parameter, public :: atomborder_def = 0.05_c_float / bohrtoa ! atom border
  real*8, parameter, public :: atomborder_criticalpoints_def = 0.03_c_float / bohrtoa ! atom border (critical points)
  real*8, parameter, public :: atomborder_gradientpaths_def = 0.008_c_float / bohrtoa ! atom border (gradient paths)
  real*8, parameter, public :: atomcovradscale_def = 0.7_c_float ! atomic radius scale factor (covalent)
  real*8, parameter, public :: atomvdwradscale_def = 1.0_c_float ! atomic radius scale factor (vdw)
  real*8, parameter, public :: atomconstantrad_def = 0.4_c_float / bohrtoa ! atomic radius scale factor (vdw)
  real*8, parameter, public :: atomrad_licorice_def = 0.11_c_float / bohrtoa ! atomic radius value (licorice)
  real*8, parameter, public :: atomrad_criticalpoints_def = 0.13_c_float / bohrtoa ! atomic radius value (critical points)
  real*8, parameter, public :: atomrad_gradientpaths_def = 0.05_c_float / bohrtoa ! atomic radius value (gradient paths)
  !--> bonds
  real*8, parameter, public :: bondrad_def = 0.125d0 / bohrtoa ! bond radius
  real*8, parameter, public :: bondrad_licorice_def = 0.25d0 / bohrtoa ! bond radius (licorice)
  real*8, parameter, public :: bondrad_vdwcontacts_def = 0.15d0 / bohrtoa ! bond radius (vdw contacts)
  real*8, parameter, public :: bondfactor_vdwcontacts_def = 1.0d0 ! bond factor (vdw contacts: cutoff = sum of vdw radii)
  real*8, parameter, public :: bondfactor_hbonds_def = 1.2d0 ! bond factor for H-bonds (cutoff ~ weak H...A max, ~3.2 A for H...O)
  real*8, parameter, public :: hbond_dist_def(2) = (/1.5d0,2.2d0/) / bohrtoa ! H...A distance class boundaries (strong|mod, mod|weak)
  real*8, parameter, public :: hbond_ang_def(2) = (/130d0,170d0/) ! D-H...A angle class boundaries (weak|mod, mod|strong), degrees
  real*8, parameter, public :: hbond_angmin_def = 90d0 ! D-H...A angle below which a contact is not an H-bond (discarded)
  real*8, parameter, public :: bondborder_def = 0.03d0 / bohrtoa ! bond border
  real*8, parameter, public :: bondborder_stickflav_def = 0.025d0 / bohrtoa ! bond border (stick flavor)
  !--> unit cell
  real*8, parameter, public :: uc_radius_def = 0.08d0 / bohrtoa ! radius of sticks
  real*8, parameter, public :: uc_radiusinner_def = 0.08d0 / bohrtoa ! radius of inner sticks
  real*8, parameter, public :: uc_innersteplen_def = 1.0d0 / bohrtoa ! length of stipple
  !--> axes
  real*8, parameter, public :: axes_length_def = 1.0d0 / bohrtoa ! length of each axis
  real*8, parameter, public :: axes_radius_def = 0.06d0 / bohrtoa ! radius of the axis shafts
  real*8, parameter, public :: axes_conelength_def = 0.30d0 / bohrtoa ! length of the arrow head
  real*8, parameter, public :: axes_coneradius_def = 0.15d0 / bohrtoa ! base radius of the arrow head
  real*8, parameter, public :: axes_winpos_def(2) = (/0.10d0,0.15d0/) ! default axes position relative to window borders
  real*8, parameter, public :: axes_labeldistance_def = 0.2d0 / bohrtoa ! distance between label and arrow head
  real*8, parameter, public :: axes_labelscale_def = 0.3d0 ! label size
  real*8, parameter, public :: axes_winfrac_def = 0.15d0 ! window-anchored axes length as a fraction of the scene radius (auto-scale tuning knob)
  !--> rotation axis
  real*8, parameter, public :: rotaxis_radius_def = 0.05d0 / bohrtoa ! radius of the rotation-axis cylinder
  !--> symmetry elements
  real(c_float), parameter, public :: symelem_rgb_def(3) = (/0.85_c_float,0.10_c_float,0.85_c_float/) ! mirror-plane / default color
  real(c_float), parameter, public :: symelem_alpha = 0.3_c_float ! mirror-plane fill opacity (axes/frames are opaque)
  real*8, parameter, public :: symelem_margin = 1.1d0 ! symmetry element size factor
  real*8, parameter, public :: symelem_frame_radius = rotaxis_radius_def ! radius of the plane-border cylinders
  real*8, parameter, public :: symelem_axis_radius = 0.15d0 / bohrtoa ! radius of the axis cylinders
  ! per-order axis colors (used when nonzero; otherwise symelem_rgb_def)
  real(c_float), parameter, public :: symelem_rgb_order(3,2:6) = reshape((/&
     0.85_c_float,0.10_c_float,0.10_c_float,&   ! 2-fold: red
     0.10_c_float,0.65_c_float,0.10_c_float,&   ! 3-fold: green
     0.10_c_float,0.30_c_float,0.90_c_float,&   ! 4-fold: blue
     0.00_c_float,0.00_c_float,0.00_c_float,&   ! 5-fold: (unused, default)
     0.90_c_float,0.45_c_float,0.05_c_float/),& ! 6-fold: orange
     shape(symelem_rgb_order))

  !> Draw style for atoms (geometry-dependent parameters)
  type atom_geom_style
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer :: type ! atom style type (attlisttype_* in systems module)
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     real(c_float), allocatable :: rgb(:,:) ! color (3,ntype)
     real*8, allocatable :: rad(:) ! radius (ntype)
   contains
     procedure :: reset => atom_style_reset
     procedure :: reset_colors => atom_style_reset_colors
     procedure :: end => atom_style_end
  end type atom_geom_style
  public :: atom_geom_style

  !> Draw style for molecules (geometry-dependent parameters)
  type mol_geom_style
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     real(c_float), allocatable :: tint_rgb(:,:) ! tint color (3,ntype)
     real*8, allocatable :: scale_rad(:) ! scale radius (ntype)
   contains
     procedure :: reset => mol_style_reset
     procedure :: end => mol_style_end
  end type mol_geom_style
  public :: mol_geom_style

  !> Draw style for bonds (geometry-dependent parameters)
  type bond_geom_style
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     logical :: use_sys_nstar = .true. ! whether this is using the system's neighbor star
     logical, allocatable :: shown(:,:) ! by-species bond shown flags (nspc,nspc)
     type(neighstar), allocatable :: nstar(:) ! the neighbor star
   contains
     procedure :: generate_neighstars
     procedure :: copy_neighstars_from_system
     procedure :: reset => bond_style_reset
     procedure :: end => bond_style_end
  end type bond_geom_style
  public :: bond_geom_style

  !> Draw style for labels (geometry-dependent parameters)
  type label_geom_style
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     character*32, allocatable :: str(:) ! text
   contains
     procedure :: reset => label_style_reset
     procedure :: end => label_style_end
  end type label_geom_style
  public :: label_geom_style

  !> Draw style for coordination polyhedra (geometry-dependent parameters).
  !> Centers are enumerated by atom-list type (species/nneq/ncel/...). For each
  !> center type, the style records whether it is shown, which species act as
  !> polyhedron corners, and the min/max center-corner distance window.
  type coordpoly_geom_style
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer :: type ! center atom-list type (atlisttype_* in systems module)
     integer :: ntype = 0 ! number of center types
     logical, allocatable :: shown(:) ! draw polyhedra for this center type (ntype)
     logical, allocatable :: corner(:,:) ! species j is a corner of center type i (nspc,ntype)
     real*8, allocatable :: dmin(:) ! min center-corner distance per center type (ntype, bohr)
     real*8, allocatable :: dmax(:) ! max center-corner distance per center type (ntype, bohr)
   contains
     procedure :: reset => coordpoly_style_reset
     procedure :: end => coordpoly_style_end
  end type coordpoly_geom_style
  public :: coordpoly_geom_style

  !> Draw style for symmetry elements (geometry-dependent parameters). Holds a
  !> snapshot of the system's symmetry operations (kind/direction/order/label),
  !> refreshed from list_symops when the geometry changes, plus the per-operation
  !> visibility selected by the user.
  type symelem_style
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer :: nop = 0 ! number of symmetry operations
     logical, allocatable :: shown(:) ! per-op on/off (nop)
     integer, allocatable :: kind(:) ! per-op element kind (symop_kind_*; 0=none) (nop)
     real*8, allocatable :: dir(:,:) ! per-op cartesian unit dir/normal (3,nop)
     integer, allocatable :: order(:) ! per-op rotation order (nop)
     character(len=mlen), allocatable :: label(:) ! per-op label: HM symbol or molecular sym string (nop)
   contains
     procedure :: reset => symelem_style_reset
     procedure :: end => symelem_style_end
  end type symelem_style
  public :: symelem_style

  ! types of representations
  integer, parameter, public :: reptype_none = 0
  integer, parameter, public :: reptype_atoms = 1 ! atoms/bonds/labels
  integer, parameter, public :: reptype_unitcell = 2 ! unit cell
  integer, parameter, public :: reptype_axes = 3 ! cartesian/crystallographic axes gizmo
  integer, parameter, public :: reptype_rotaxis = 4 ! rotation axis for a molecule
  integer, parameter, public :: reptype_symelem = 5 ! symmetry element
  integer, parameter, public :: reptype_NUM = 5

  ! representation flavors
  integer, parameter, public :: repflavor_unknown = 0
  integer, parameter, public :: repflavor_atoms_ballandstick = 1
  integer, parameter, public :: repflavor_atoms_sticks = 2
  integer, parameter, public :: repflavor_atoms_licorice = 3
  integer, parameter, public :: repflavor_atoms_vdwcontacts = 4
  integer, parameter, public :: repflavor_atoms_hbonds = 5
  integer, parameter, public :: repflavor_atoms_criticalpoints = 6
  integer, parameter, public :: repflavor_atoms_gradientpaths = 7
  integer, parameter, public :: repflavor_atoms_polyhedra = 8
  integer, parameter, public :: repflavor_unitcell_basic = 9
  integer, parameter, public :: repflavor_axes = 10
  integer, parameter, public :: repflavor_rotaxis = 11
  integer, parameter, public :: repflavor_symelem = 12
  integer, parameter, public :: repflavor_NUM = 12

  !> Selection of the part of the system that is drawn: periodicity, origin
  !> shift, display region, and the atom filter (reptype_atoms; pertype, ncell
  !> and origin also control reptype_unitcell). Accessed as r%sel%...
  type rep_selection
     integer(c_int) :: pertype ! periodicity control: 0=none, 1=auto, 2=manual
     integer(c_int) :: ncell(3) ! number of unit cells drawn
     real*8 :: origin(3) ! origin shift of the representation
     real*8 :: tshift(3) ! origin of the unit cell display region
     logical :: border ! draw atoms at the border of the unit cell
     logical :: onemotif ! draw connected molecules
     character(kind=c_char,len=:), allocatable :: filter ! filter for the representation
     character(kind=c_char,len=:), allocatable :: errfilter ! filter error
  end type rep_selection
  public :: rep_selection

  !> Atom display options (reptype_atoms; accessed as r%atoms%...)
  type rep_atoms
     logical :: display ! whether to draw the atoms
     type(atom_geom_style) :: style ! atom styles (geometry-dependent)
     integer(c_int) :: radii_type ! option to reset radii: 0=covalent,1=vdw,2=constant
     real*8 :: radii_scale ! reset radii, scale factor
     real*8 :: radii_value ! reset radii, constant value (default: same as bond radius for licorice)
     integer(c_int) :: color_type ! option to reset colors: 0=current,1=jmlcol,2=jmlcol2
     real*8 :: border_size ! atom border size
     real(c_float) :: border_rgb(3) ! atom border color
  end type rep_atoms
  public :: rep_atoms

  !> Bond display options (reptype_atoms; accessed as r%bonds%...)
  type rep_bonds
     logical :: display ! whether to draw the bonds
     type(bond_geom_style) :: style ! bond styles (geometry-dependent)
     real*8 :: atmrad(0:maxzat0) = atmcov0 ! per-species covalent radii for bonding (bohr)
     real*8 :: bfactor = bondfactor_def ! bond factor for non-metal bonding
     real*8 :: bdelta = bonddelta_def ! bond delta for metal bonding (bohr)
     integer(c_int) :: color_style ! bond style (0=single color, 1=two colors)
     real*8 :: rad ! radius
     real*8 :: border_size ! bond border size
     real(c_float) :: border_rgb(3) ! bond border color
     real(c_float) :: rgb(3) ! color
     integer(c_int) :: order ! order (0=dashed,1=single,2=double,3=triple,4=ordcon value)
     integer(c_int) :: imol ! molecular connections (0=any,1=intramol,2=intermol)
     logical :: bothends ! if true, both atoms need to be drawn to draw the bond
     logical :: hbond_classify = .false. ! Jeffrey-Steiner hydrogen-bond strength classification
     real(c_float) :: hbond_rgb(3,3) = reshape((/& ! per-class colors
        0.00_c_float,0.00_c_float,1.00_c_float,&  ! strong: blue
        0.00_c_float,0.70_c_float,0.00_c_float,&  ! moderate: green
        1.00_c_float,0.00_c_float,0.00_c_float/),(/3,3/)) ! weak: red
     real*8 :: hbond_dist(2) = hbond_dist_def ! H...A distance class boundaries (strong|mod, mod|weak), bohr
     real*8 :: hbond_ang(2) = hbond_ang_def ! D-H...A angle class boundaries (weak|mod, mod|strong), degrees
  end type rep_bonds
  public :: rep_bonds

  !> Label display options (reptype_atoms; accessed as r%labels%...)
  type rep_labels
     logical :: display ! whether to draw the labels
     type(label_geom_style) :: style ! label styles (geometry-dependent)
     integer(c_int) :: type ! 0=atom-symbol, 1=atom-name, 2=cel-atom, 3=cel-atom+lvec, 4=neq-atom, 5=spc, 6=Z, 7=mol, 8=wyckoff
     real*8 :: scale ! scale for the labels
     real(c_float) :: rgb(3) ! color of the labels
     logical :: const_size ! whether labels scale with objects or are constant size
     real*8 :: offset(3) ! offset of the label
  end type rep_labels
  public :: rep_labels

  !> Unit cell display options (reptype_unitcell; accessed as r%uc%...)
  type rep_unitcell
     logical :: inner ! display inner cylinders
     logical :: coloraxes ! color the axes (x=red,y=green,z=blue)
     logical :: vaccutsticks ! cut sticks in systems with vacuum
     real*8 :: radius ! cylinder radius
     real*8 :: radiusinner ! cylinder radius (inner)
     real(c_float) :: rgb(3) ! cylinder colors
     real*8 :: innersteplen ! number of subdivisions for the inner sticks
     logical :: innerstipple ! stippled lines for the inner lines
  end type rep_unitcell
  public :: rep_unitcell

  !> Cartesian/crystallographic axes options (reptype_axes; accessed as r%axes%...)
  type rep_axes
     integer(c_int) :: kind ! 0 = cartesian, 1 = crystallographic
     real*8 :: rot(3,3) = eye ! orientation applied to the axis directions (columns are the axes); identity by default
     integer(c_int) :: placement ! 0 = at the origin, 1 = anchored at a fixed window position
     integer(c_int) :: coordtype ! origin coordinates: 0 = crystallographic, 1 = cartesian (angstrom), 2 = cartesian (bohr)
     real*8 :: origin(3) = 0d0 ! origin of the axes (coordinates per coordtype)
     real*8 :: winpos(2) ! window position (fractions from left and bottom) when window-anchored
     real*8 :: length ! length of each axis
     real*8 :: radius ! radius of the axis shafts
     real*8 :: conelength ! length of the arrowhead cones
     real*8 :: coneradius ! base radius of the arrowhead cones
     real(c_float) :: rgb(3,3) ! color of the x, y, z axes
     logical :: showlabels ! draw x/y/z labels at the axis tips
     real*8 :: labelscale ! scale for the axis labels
     logical :: labelconstsize ! whether the labels have constant size or scale with the arrowhead
     real(c_float) :: labelrgb(3) ! color of the axis labels
     character(kind=c_char,len=32) :: labelstr(3) ! text of the x, y, z axis labels
     real*8 :: labeldistance ! distance along the axis from the arrowhead to the label (all axes)
     real*8 :: labeloffset(3,3) ! per-axis (cartesian) offset of the label from that position
     real*8 :: scale ! global scale factor applied to the whole gizmo (arrows and labels)
     logical :: scale_auto ! auto-size the window-anchored gizmo from the scene radius
     logical :: scalewithzoom ! whether the window-anchored gizmo scales when the scene is zoomed
  end type rep_axes
  public :: rep_axes

  !> Rotation axis options (reptype_rotaxis; accessed as r%rotaxis%...)
  type rep_rotaxis
     real*8 :: origin(3) = 0d0 ! origin the axis line passes through (cartesian, bohr)
     real*8 :: dir(3) = (/0d0,0d0,1d0/) ! unit direction in cartesian (bohr); the axis line passes through origin
     real*8 :: length = 0d0 ! half-length: the cylinder spans origin +/- length*dir
     real*8 :: radius = rotaxis_radius_def ! radius of the rotation-axis cylinder
     real(c_float) :: rgb(3) = 0._c_float ! color of the rotation-axis cylinder
  end type rep_rotaxis
  public :: rep_rotaxis

  !> Symmetry element options (reptype_symelem; accessed as r%symelem%...).
  !> The transient fields describe a single hover/preview element (stamped by
  !> scene_show_symelems); size and cen are per-build inputs stamped by
  !> scene_build_lists for persistent sets too. The permanent fields hold the
  !> user-editable state of a persistent symmetry-element representation.
  type rep_symelem
     !! transient
     integer :: kind = 0 ! 0=none, 1=plane (mirror), 2=axis (rotation)
     real*8 :: origin_transient(3) = 0d0 ! transient element origin (cartesian, bohr)
     real*8 :: dir(3) = (/0d0,0d0,1d0/) ! axis direction or plane normal, unit, cartesian (bohr)
     real*8 :: size = 0d0 ! system bounding-sphere radius (bohr)
     real*8 :: cen(3) = 0d0 ! system center (bohr)
     integer :: order = 0 ! axis rotation order n (selects the axis color)
     real(c_float) :: rgb(3) = symelem_rgb_def ! color of the symmetry element
     !! permanent
     type(symelem_style) :: style ! operation snapshot + per-op visibility (geometry-dependent)
     real*8 :: origin(3) = 0d0 ! editable origin the elements pass through (coords per coordtype)
     integer(c_int) :: coordtype = 0 ! origin coords: 0=crystallographic, 1=cartesian (angstrom), 2=cartesian (bohr)
     logical :: usecustomrgb = .false. ! true: use rgb for all; false: per-order/default colors
  end type rep_symelem
  public :: rep_symelem

  !> Per-molecule display options (reptype_atoms; accessed as r%mols%...)
  type rep_mols
     type(mol_geom_style) :: style ! molecule styles (geometry-dependent)
  end type rep_mols
  public :: rep_mols

  !> Coordination polyhedra options (reptype_atoms; accessed as r%poly%...)
  type rep_poly
     logical :: display ! whether to draw the coordination polyhedra
     type(coordpoly_geom_style) :: style ! center/corner/distance geometry (geometry-dependent)
     real*8 :: alpha = 0.5d0 ! face opacity (1 = opaque)
     logical :: usecentercolor = .true. ! faces take the central atom color
     real(c_float) :: rgb(3) = 0._c_float ! face color when not using the central atom color
     real*8 :: edge_rad = 0.05d0 ! edge cylinder radius (bohr)
     real(c_float) :: edge_rgb(3) = 0._c_float ! edge cylinder color
     logical :: usecentercolor_edge = .true. ! edges take the central atom color
     real*8 :: coplanar_eps = 0.1d0 ! coplanarity tolerance for the planar-polygon path (bohr)
     logical :: showcorners = .true. ! also draw the corner atoms, even if outside the selection
  end type rep_poly
  public :: rep_poly

  !> Representation: objects to draw on the scene
  type representation
     ! main variables
     logical :: isinit = .false. ! whether the representation has been initialized
     logical :: shown = .false. ! true if the representation is currently shown
     integer :: type = reptype_none ! type of representation (atoms, cell,...)
     integer :: flavor = repflavor_unknown ! flavor of the representation
     integer :: id ! system ID
     integer :: idrep ! representation ID
     integer :: iord = 0 ! representation order integer in menu
     character(kind=c_char,len=:), allocatable :: name ! name of the representation
     ! per-object option groups
     type(rep_selection) :: sel ! which part of the system is drawn
     type(rep_atoms) :: atoms ! atom display options
     type(rep_bonds) :: bonds ! bond display options
     type(rep_labels) :: labels ! label display options
     type(rep_mols) :: mols ! per-molecule display options
     type(rep_unitcell) :: uc ! unit cell display options
     type(rep_axes) :: axes ! cartesian/crystallographic axes options
     type(rep_rotaxis) :: rotaxis ! rotation axis options
     type(rep_symelem) :: symelem ! symmetry element options
     type(rep_poly) :: poly ! coordination polyhedra options
   contains
     procedure :: init => representation_init
     procedure :: set_defaults => representation_set_defaults
     procedure :: end => representation_end
     procedure :: update => update_styles
     procedure :: add_draw_elements
     procedure :: reset_all_styles
  end type representation
  public :: representation

  ! module procedure interfaces
  interface
     module subroutine representation_init(r,isys,irep,itype,flavor,icount)
       class(representation), intent(inout) :: r
       integer, intent(in) :: isys
       integer, intent(in) :: irep
       integer, intent(in) :: itype
       integer, intent(in) :: flavor
       integer, intent(inout) :: icount(0:repflavor_NUM)
     end subroutine representation_init
     module subroutine representation_set_defaults(r,itype)
       class(representation), intent(inout) :: r
       integer, intent(in) :: itype
     end subroutine representation_set_defaults
     module subroutine representation_end(r)
       class(representation), intent(inout) :: r
     end subroutine representation_end
     module subroutine update_styles(r)
       class(representation), intent(inout) :: r
     end subroutine update_styles
     module subroutine add_draw_elements(r,nc,obj,doanim,iqpt,ifreq)
       class(representation), intent(inout) :: r
       integer, intent(in) :: nc(3)
       type(scene_objects), intent(inout) :: obj
       logical, intent(in) :: doanim
       integer, intent(in) :: iqpt, ifreq
     end subroutine add_draw_elements
     module subroutine reset_all_styles(r,itype)
       class(representation), intent(inout) :: r
       integer, intent(in) :: itype
     end subroutine reset_all_styles
     ! atom_geom_style
     module subroutine atom_style_reset(d,r)
       class(atom_geom_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine atom_style_reset
     module subroutine atom_style_reset_colors(d,r)
       class(atom_geom_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine atom_style_reset_colors
     module subroutine atom_style_end(d)
       class(atom_geom_style), intent(inout) :: d
     end subroutine atom_style_end
     ! mol_geom_style
     module subroutine mol_style_reset(d,r)
       class(mol_geom_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine mol_style_reset
     module subroutine mol_style_end(d)
       class(mol_geom_style), intent(inout) :: d
     end subroutine mol_style_end
     ! bond_geom_style
     module subroutine generate_neighstars(d,r)
       class(bond_geom_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine generate_neighstars
     module subroutine copy_neighstars_from_system(d,isys)
       class(bond_geom_style), intent(inout) :: d
       integer, intent(in) :: isys
     end subroutine copy_neighstars_from_system
     module subroutine bond_style_reset(d,r)
       class(bond_geom_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine bond_style_reset
     module subroutine bond_style_end(d)
       class(bond_geom_style), intent(inout) :: d
     end subroutine bond_style_end
     ! label_geom_style
     module subroutine label_style_reset(d,r)
       class(label_geom_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine label_style_reset
     module subroutine label_style_end(d)
       class(label_geom_style), intent(inout) :: d
     end subroutine label_style_end
     module subroutine coordpoly_style_reset(d,r)
       class(coordpoly_geom_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine coordpoly_style_reset
     module subroutine coordpoly_style_end(d)
       class(coordpoly_geom_style), intent(inout) :: d
     end subroutine coordpoly_style_end
     module subroutine symelem_style_reset(d,r)
       class(symelem_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine symelem_style_reset
     module subroutine symelem_style_end(d)
       class(symelem_style), intent(inout) :: d
     end subroutine symelem_style_end
  end interface

end module representations
