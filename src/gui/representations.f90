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
  use shapes, only: dl_sphere, dl_cylinder, dl_cylinder_over, dl_string, dl_string_over,&
     dl_plane, dl_triangle, dl_mesh, scene_objects, dl_append
  use param, only: bohrtoa, eye, maxzat0, atmcov0, mlen
  use grid3mod, only: hscale_num, hscale_linear, hscale_log, hscale_asinh
  use utils, only: iw_cmap_viridis, iw_cmap_rdbu, iw_colormap_lut
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
  logical, parameter, public :: occ_sectors_def = .true. ! render partial occupancies as sectors
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
  logical, parameter, public :: hbond_classify_def = .false. ! Jeffrey-Steiner H-bond strength classification
  !--> labels
  real*8, parameter, public :: label_scale_def = 0.5d0 ! label size
  real*8, parameter, public :: label_scale_criticalpoints_def = 0.3d0 ! label size (critical points)
  real*8, parameter, public :: label_offset_criticalpoints_def(3) = (/0d0,0.25d0,0d0/) ! label offset (critical points)
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
  !--> text annotations
  integer, parameter, public :: textpos_screen = 0 ! anchored to a viewport position
  integer, parameter, public :: textpos_point = 1 ! at a 3D position in the system
  integer, parameter, public :: textpos_atom = 2 ! tied to an atom (cell atom + lattice vector)
  integer, parameter, public :: textpos_bond = 3 ! tied to a bond (two atom images)
  real*8, parameter, public :: text_scale_def = label_scale_def ! text size (same semantics as labels)
  real*8, parameter, public :: text_winpos_def(2) = (/0.5d0,0.85d0/) ! default viewport position (near top-center)
  !--> coordination polyhedra
  real*8, parameter, public :: polyalpha_def = 0.75d0 ! face opacity (1 = opaque)
  real*8, parameter, public :: polyedgerad_def = 0.06d0 / bohrtoa ! edge cylinder radius
  real*8, parameter, public :: polycoplanar_def = 0.1d0 ! coplanarity tolerance for the planar-polygon path
  real(c_float), parameter, public :: polyface_rgb_def(3) = 0._c_float ! face color (when not using the central atom color)
  real(c_float), parameter, public :: polyedge_rgb_def(3) = 0._c_float ! edge color (when not using the central atom color)
  logical, parameter, public :: poly_usecentercolor_def = .true. ! faces take the central atom color
  logical, parameter, public :: poly_usecentercolor_edge_def = .false. ! edges take the central atom color
  logical, parameter, public :: poly_showcorners_def = .true. ! also draw the corner atoms outside the selection
  !--> measurements
  real*8, parameter, public :: measure_rad_def = 0.04d0 / bohrtoa ! radius of the measurement segments/edges
  real*8, parameter, public :: measure_sectorrad_def = 1.2d0 / bohrtoa ! radius of the angle/dihedral sectors
  real*8, parameter, public :: measure_planealpha_def = 0.4d0 ! opacity of the dihedral planes
  real*8, parameter, public :: measure_planeext_def = 0.1d0 ! fractional extension of the dihedral planes beyond the atoms
  real*8, parameter, public :: measure_sectoralpha_def = 0.5d0 ! opacity of the angle/dihedral sectors
  real*8, parameter, public :: measure_textscale_def = 0.75d0 ! size of the numeric labels
  real*8, parameter, public :: measure_dashlen_def = 0.15d0 / bohrtoa ! dash period for dashed segments
  integer, parameter, public :: measure_ndeclen_def = 3 ! decimal places for distances (angstrom)
  integer, parameter, public :: measure_ndecang_def = 1 ! decimal places for angles (degrees)
  real(c_float), parameter, public :: measure_rgb_dist_def(3) = (/0.90_c_float,0.65_c_float,0.10_c_float/) ! distance color (gold)
  real(c_float), parameter, public :: measure_rgb_ang_def(3) = (/0.10_c_float,0.65_c_float,0.85_c_float/) ! angle color (cyan)
  real(c_float), parameter, public :: measure_rgb_dih_def(3) = (/0.90_c_float,0.30_c_float,0.60_c_float/) ! dihedral color (pink)
  !--> isosurfaces
  integer, parameter, public :: iso_nlevel = 4 ! number of named coarseness levels (0 = native; 1..nlevel, set level)
  real*8, parameter, public :: iso_level_ptsang(iso_nlevel) = (/2d0,5d0,10d0,20d0/) ! points/ang of the named levels
  integer, parameter, public :: iso_level_custom = iso_nlevel + 1 ! level index for a custom points/ang
  integer, parameter, public :: iso_level_def = 2 ! default named level (medium)
  integer, public :: iso_defaultlevel = iso_level_def ! preference: level for newly created isosurfaces
  integer, parameter, public :: iso_npts_custom_min = 2 ! minimum points per axis (custom level)
  integer, parameter, public :: iso_npts_custom_max = 512 ! maximum points per axis (custom level)
  integer, parameter, public :: iso_maxpts_total = 4000000 ! total-points cap (~32 MB of cached samples); coarsen uniformly above it
  integer, parameter, public :: iso_npts_axmin = 16 ! per-axis minimum number of points
  character(len=*), parameter, public :: iso_level_optstr = &
     "Coarse (2 pts/Å)"//c_null_char//"Medium (5 pts/Å)"//c_null_char//&
     "Fine (10 pts/Å)"//c_null_char//"Very fine (20 pts/Å)"//c_null_char ! named levels
  character(len=*), parameter, public :: iso_level_optstr_custom = &
     iso_level_optstr//"Custom"//c_null_char ! named levels plus custom
  integer, parameter, public :: iso_region_cell = 0 ! region modes: whole unit cell (native level allowed; periodic in crystals)
  integer, parameter, public :: iso_region_frac = 1 ! cell-aligned box between fractional points x0 and x1 (crystals)
  integer, parameter, public :: iso_region_ortho = 2 ! axis-aligned Cartesian box between corners x0 and x1 (ang)
  integer, parameter, public :: iso_region_parallel = 3 ! origin x0 plus edge endpoints x1, x2, x3 (ang)
  integer, parameter, public :: iso_region_simplebox = 4 ! axis-aligned box: center x0, half-lengths x1 (ang; molecules)
  integer, parameter, public :: iso_region_cube = 5 ! cube: center x0, half-length x(1,1) (ang; molecules)
  integer, parameter, public :: iso_region_bbox = 6 ! bounding box of the atoms plus a buffer x(1,1) (ang; molecules)
  integer, parameter, public :: iso_region_NUM = 6 ! highest region mode
  integer, parameter, public :: iso_knd_cry = 1 ! display names, by mode and system kind (1 = crystal, 2 = molecule)
  integer, parameter, public :: iso_knd_mol = 2
  character(len=33), parameter, public :: iso_region_name(0:iso_region_NUM,2) = &
     reshape((/ character(len=33) :: &
     "Unit Cell", "Cell fractions", "Box (2 corners)",&
     "Parallelepiped (origin + 3 sides)", "", "", "",&
     "Cell", "", "Box (2 corners)", "Parallelepiped (origin + 3 sides)",&
     "Box (center + 3 half-lengths)", "Cube (center + half-length)",&
     "Expanded bounding box" /),(/iso_region_NUM+1,2/))
  ! one-line explanation of each mode, for the editor's tooltip. Kind-
  ! independent, unlike the names: every mode whose meaning differs between
  ! a crystal and a molecule is offered to only one of the two anyway
  character(len=64), parameter, public :: iso_region_desc(0:iso_region_NUM) = &
     (/ character(len=64) :: &
     "The whole cell (in a molecule, the box built around it)",&
     "A cell-aligned box between two points in fractional coordinates",&
     "An axis-aligned box between two Cartesian corners",&
     "An origin and the endpoints of the three edges",&
     "A center and a half-length along each axis",&
     "A center and a half-length",&
     "The box bounding the atoms, grown by a buffer on every side" /)
  integer, parameter, public :: iso_region_modes_cry(4) = (/iso_region_cell,&
     iso_region_frac,iso_region_ortho,iso_region_parallel/) ! region modes offered to each system kind
  integer, parameter, public :: iso_region_modes_mol(6) = (/iso_region_bbox,&
     iso_region_cell,iso_region_ortho,iso_region_parallel,iso_region_simplebox,&
     iso_region_cube/)
  real*8, parameter, public :: iso_bbox_buffer_def = 5d0 ! default buffer around the bounding box (ang)
  real*8, parameter, public :: iso_isoval_def = 0.1d0 ! default isovalue when no field statistics are available (a.u.)
  real*8, parameter, public :: iso_isoval_mo = 0.02d0 ! default isovalue for molecular orbitals (a.u.; the convention across molecular viewers)
  real*8, parameter, public :: iso_isoval_dens = 1d-3 ! default-isovalue policy: conventional molecular density contour (a.u.; same as SIGMAHOLE)
  real*8, parameter, public :: iso_qcharge_def = 0.9d0 ! fraction of the |f| integral enclosed by the default level
  real*8, parameter, public :: iso_spikeratio = 10d0 ! max|f|/mean|f| above which the field is spike-dominated
  integer, parameter, public :: iso_nhist = 96 ! bins of the value histogram shown in the editor
  character(len=8), parameter, public :: iso_hscale_name(hscale_num) = &
     (/ character(len=8) :: "Linear", "Log", "arcsinh" /) ! names of the histogram axis scales
  integer, parameter, public :: iso_hscale_y(2) = (/hscale_linear, hscale_log/) ! count axis scales
  real(c_float), parameter, public :: iso_alpha_def = 0.75_c_float ! default opacity
  integer, parameter, public :: iso_npalette = 8 ! number of palette colors
  real(c_float), parameter, public :: iso_rgb_palette(3,iso_npalette) = reshape((/&
     0.30_c_float,0.55_c_float,0.90_c_float,& ! blue
     0.90_c_float,0.45_c_float,0.15_c_float,& ! orange
     0.25_c_float,0.70_c_float,0.35_c_float,& ! green
     0.85_c_float,0.25_c_float,0.25_c_float,& ! red
     0.60_c_float,0.35_c_float,0.80_c_float,& ! purple
     0.20_c_float,0.75_c_float,0.80_c_float,& ! cyan
     0.90_c_float,0.75_c_float,0.20_c_float,& ! yellow
     0.90_c_float,0.45_c_float,0.70_c_float/),& ! pink
     (/3,iso_npalette/))
  real(c_float), parameter, public :: iso_rgb_def(3) = iso_rgb_palette(:,1) ! color of the first isosurface
  real(c_float), parameter, public :: iso_rgb_tol = 0.05_c_float ! two colors closer than this count as the same
  integer, parameter, public :: iso_map_color = 0 ! isosurface color: a single color for the whole surface
  integer, parameter, public :: iso_map_field = 1 ! the values of another field, through a colormap
  integer, parameter, public :: iso_map_expr = 2 ! the values of an arithmetic expression
  character(len=*), parameter, public :: iso_map_optstr = &
     "Color"//c_null_char//"Field map"//c_null_char//"Expression"//c_null_char ! in iso_map_* order
  integer, parameter, public :: iso_explen = 255 ! length of the expression that colors a surface
  integer, parameter, public :: iso_cmap_seq = iw_cmap_viridis ! colormap from utils (iw_cmap_*)
  integer, parameter, public :: iso_cmap_div = iw_cmap_rdbu
  integer, parameter, public :: iso_nlut = 256 ! entries of the colormap look-up table
  real(c_float), parameter, public :: iso_rgb_invalid(3) = 0.5_c_float ! color of a vertex outside the map field

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
     procedure :: alloc => symelem_style_alloc
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
  integer, parameter, public :: reptype_text = 6 ! user text annotations
  integer, parameter, public :: reptype_measure = 7 ! measurements (distances/angles/dihedrals)
  integer, parameter, public :: reptype_shapes = 8 ! list of geometric shapes
  integer, parameter, public :: reptype_isosurface = 9 ! isosurface of a scalar field
  integer, parameter, public :: reptype_NUM = 9

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
  integer, parameter, public :: repflavor_text = 13
  integer, parameter, public :: repflavor_measure = 14
  integer, parameter, public :: repflavor_shapes = 15
  integer, parameter, public :: repflavor_isosurface = 16
  integer, parameter, public :: repflavor_NUM = 16

  !> Selection of the part of the system that is drawn: periodicity, origin
  !> shift, display region, and the atom filter (reptype_atoms; pertype, ncell
  !> and origin also control reptype_unitcell; reptype_isosurface consumes
  !> pertype and ncell only -- the origin shift does not move the isosurface).
  !> Accessed as r%sel%...
  type rep_selection
     integer(c_int) :: pertype ! periodicity control: 0=none, 1=auto, 2=manual
     integer(c_int) :: ncell(3) ! number of unit cells drawn
     real*8 :: origin(3) ! origin shift of the representation
     real*8 :: tshift(3) ! origin of the unit cell display region
     logical :: border ! draw atoms at the border of the unit cell
     logical :: onemotif ! draw connected molecules
     character(kind=c_char,len=:), allocatable :: filter ! filter for the representation
     character(kind=c_char,len=:), allocatable :: errfilter ! filter error
   contains
     procedure :: ncells => rep_selection_ncells ! number of drawn cells from the periodicity control
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
     logical :: occ_sectors ! render partial occupancies as sectors
     real(c_float) :: occ_empty_rgb(3) ! color of the empty (unoccupied) sector
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
     logical :: hbond_classify = hbond_classify_def ! Jeffrey-Steiner hydrogen-bond strength classification
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

  ! shape kinds for the shapes representation
  integer, parameter, public :: shapekind_sphere = 1
  integer, parameter, public :: shapekind_box = 2

  !> A geometric shape in a shapes representation
  type rep_shape
     integer :: kind = shapekind_sphere ! shape kind (shapekind_*)
     real*8 :: x1(3) = 0d0 ! center/anchor (cartesian, bohr; molecules: absolute frame)
     real*8 :: v(3,3) = 0d0 ! box: the three edge vectors from x1 (cartesian, bohr)
     real*8 :: rad = 1d0 ! sphere: radius; box: thickness of the edges (bohr)
     real(c_float) :: rgb(3) = 0._c_float ! color
     real(c_float) :: alpha = 1._c_float ! sphere: opacity (1 = opaque); box: opacity of the faces (0 = wireframe only)
  end type rep_shape
  public :: rep_shape

  !> Shapes representation options (reptype_shapes; accessed as r%shapes%...)
  type rep_shapes
     integer :: nshape = 0 ! number of shapes in the list
     type(rep_shape), allocatable :: shape(:) ! the shapes
  end type rep_shapes
  public :: rep_shapes

  !> Symmetry element options (reptype_symelem; accessed as r%symelem%...).
  type rep_symelem
     type(symelem_style) :: style ! operation snapshot + per-op visibility (geometry-dependent)
     real*8 :: origin(3) = 0d0 ! origin the elements pass through (coords per coordtype)
     integer(c_int) :: coordtype = 0 ! origin coords: 0=crystallographic, 1=cartesian (angstrom), 2=cartesian (bohr)
     logical :: usecustomrgb = .false. ! true: use rgb for all; false: per-order/default colors
     real(c_float) :: rgb(3) = symelem_rgb_def ! color of the symmetry elements (usecustomrgb)
     real*8 :: size = 0d0 ! system bounding-sphere radius for sizing (bohr; stamped in build_lists)
     real*8 :: cen(3) = 0d0 ! system center for sizing (bohr; stamped in build_lists)
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
     real*8 :: alpha = polyalpha_def ! face opacity (1 = opaque)
     logical :: usecentercolor = poly_usecentercolor_def ! faces take the central atom color
     real(c_float) :: rgb(3) = polyface_rgb_def ! face color when not using the central atom color
     real*8 :: edge_rad = polyedgerad_def ! edge cylinder radius (bohr)
     real(c_float) :: edge_rgb(3) = polyedge_rgb_def ! edge cylinder color
     logical :: usecentercolor_edge = poly_usecentercolor_edge_def ! edges take the central atom color
     real*8 :: coplanar_eps = polycoplanar_def ! coplanarity tolerance for the planar-polygon path (bohr)
     logical :: showcorners = poly_showcorners_def ! also draw the corner atoms, even if outside the selection
  end type rep_poly
  public :: rep_poly

  !> User text annotation (reptype_text)
  type text_item
     logical :: shown = .true. ! whether this text is drawn
     character(kind=c_char,len=:), allocatable :: str ! the text (newlines allowed)
     integer(c_int) :: placement = textpos_screen ! placement mode (textpos_*)
     real*8 :: winpos(2) = text_winpos_def ! viewport position, fractions from left/bottom (textpos_screen)
     real*8 :: pos(3) = 0d0 ! 3D anchor: fractional (crystal) or angstrom (molecule) (textpos_point)
     integer(c_int) :: idx1(4) = 0 ! atom anchor/bond atom 1: cell atom + lattice vector; 0 = unset
     integer(c_int) :: idx2(4) = 0 ! bond atom 2: cell atom + lattice vector; 0 = unset
     real*8 :: scale = text_scale_def ! text size
     real(c_float) :: rgb(3) = 0._c_float ! text color
     logical :: scalewithzoom = .false. ! text scales with zoom (vs constant on-screen size)
     real*8 :: offset(2) = 0d0 ! in-plane on-screen offset (angstrom; 3D placements only)
     logical :: depth = .true. ! 3D placements: text is occluded by objects in front of it
     logical :: infront = .true. ! on-screen placement: text in front of (vs behind) the scene
  end type text_item
  public :: text_item

  !> Text annotation options (reptype_text; accessed as r%text%...)
  type rep_text
     integer :: ntext = 0 ! number of text items
     type(text_item), allocatable :: t(:) ! the text items
  end type rep_text
  public :: rep_text

  !> A measurement (reptype_measure). The number of atoms sets the
  !> kind: 2 = distance (Å), 3 = angle (deg, vertex = atom 2), 4 =
  !> dihedral (deg, axis = atoms 2-3).
  !> A measurement (reptype_measure). Style is per-item: each measurement
  !> carries its own color, radii, opacities, label size, and decimals.
  type measurement_item
     logical :: shown = .true. ! whether this measurement is drawn
     integer :: n = 0 ! number of atoms (2=distance, 3=angle, 4=dihedral)
     integer(c_int) :: idx(4,4) = 0 ! per-atom anchor (1:4,iatom): cell atom id + lattice vector
     real(c_float) :: rgb(3) = measure_rgb_dist_def ! color
     real*8 :: rad = measure_rad_def ! radius of the segment/arms/edges (bohr)
     real*8 :: sectorrad = measure_sectorrad_def ! radius of the angle/dihedral sector (bohr)
     real*8 :: sectoralpha = measure_sectoralpha_def ! opacity of the angle/dihedral sector
     real*8 :: planealpha = measure_planealpha_def ! opacity of the dihedral planes
     real*8 :: planeext = measure_planeext_def ! fractional extension of the dihedral planes beyond the atoms
     real*8 :: textscale = measure_textscale_def ! size of the numeric label
     integer :: ndec = measure_ndeclen_def ! decimal places for this item's value
     logical :: dashed = .false. ! draw the segment as a dashed (vs solid) line
     real*8 :: dashlen = measure_dashlen_def ! dash period for a dashed segment (bohr)
     logical :: orient = .true. ! orient the label along the segment (vs always horizontal)
     logical :: scalesystem = .true. ! label scales with the system size (vs constant on-screen size)
     real*8 :: offset(2) = 0d0 ! in-plane (screen-relative) offset of the label, angstrom
   contains
     procedure :: set_defaults => measurement_item_set_defaults
     procedure :: copy_style => measurement_item_copy_style
  end type measurement_item
  public :: measurement_item

  !> Measurement options (reptype_measure)
  type rep_measure
     integer :: nitem = 0 ! number of measurement items
     type(measurement_item), allocatable :: item(:) ! the measurement items
     integer :: isel = 0 ! selected item (for the per-item editor under the table)
  end type rep_measure
  public :: rep_measure

  !> One isosurface of the isosurface object: the level, how it is drawn,
  !> and the triangulation cached for it. The build stamp is on the
  !> isosurface itself, so isosurfaces can be added, removed, and
  !> re-leveled without disturbing the meshes of the others.
  type iso_slot
     real*8 :: isoval = 0d0 ! isovalue
     logical :: shown = .true. ! whether this isosurface is drawn in the scene
     real(c_float) :: rgb(3) = iso_rgb_def ! color (imap_mode = iso_map_color)
     real(c_float) :: alpha = iso_alpha_def ! opacity (1 = opaque)
     ! coloring by the values of another field (imap_mode = iso_map_field)
     integer :: imap_mode = iso_map_color ! where the color comes from (iso_map_*)
     integer :: imap = -1 ! field whose values color the surface (-1 = none chosen yet)
     character(len=iso_explen) :: mapexpr = "" ! expression whose values color the surface
     character(len=iso_explen) :: maperr = "" ! why that expression could not be used (empty = it could)
     integer :: icmap = iso_cmap_seq ! colormap (index into iw_cmap_name)
     logical :: icmap_auto = .true. ! colormap follows the values (sequential/diverging) until the user picks one
     real*8 :: maprange(2) = (/0d0,1d0/) ! values at the two ends of the colormap
     logical :: maprange_auto = .true. ! range follows the values on the surface until the user sets it
     real*8, allocatable :: mapval(:) ! values of the map field at the vertices (nv)
     logical(c_bool), allocatable :: mapok(:) ! whether each vertex is inside the map field's domain (nv)
     logical :: mapoutdomain = .false. ! some vertices fell outside the map field's domain
     type(dl_mesh) :: mesh ! the cached triangulation
     logical :: built = .false. ! whether mesh holds a triangulation
     real*8 :: isoval_built = 0d0 ! isovalue the mesh was built at
     integer :: imap_built = -1 ! map field the values were evaluated for (-1 = none)
     character(len=iso_explen) :: mapexpr_built = "" ! expression the values were evaluated for
     ! state the vertex colors were computed for; a change here recolors the mesh
     integer :: icmap_built = -1
     real*8 :: maprange_built(2) = 0d0
  end type iso_slot
  public :: iso_slot

  !> Isosurface display options (reptype_isosurface; accessed as r%iso%...)
  type rep_isosurface
     ! user options
     integer :: ifield = 0 ! field for the isosurface (index in sys(id)%f)
     integer :: imosel = 0 ! sample a molecular orbital of the field instead of the field itself:
                           ! an id_mo_* selector from the types module (all negative), or 0 = the field
     integer :: imoidx = 0 ! MO index accompanying imosel (id_mo_a, id_mo_b, id_mo_id)
     integer :: ilevel = iso_level_def ! staged grid coarseness level (0=native, 1..4=named, 5=custom)
     integer :: nptscustom(3) = 0 ! staged custom grid dimensions (ilevel = custom; seeded when it is selected)
     integer :: iregion = iso_region_cell ! staged region mode (iso_region_*)
     real*8 :: rgn_x(3,0:3) = 0d0 ! staged region coordinates: origin/corner/center (column 0) and far
                                  ! corner, edge endpoints, or half-lengths (columns 1-3); fractional or
                                  ! user-frame ang, per mode. The cube's single half-length is stored
                                  ! replicated across column 1, so all consumers read it generically.
     real*8 :: costest = -1d0 ! measured field-evaluation cost (seconds per sample point; <0 = none)
     integer :: nptsxyz(3) = 0 ! applied sampling grid; all-zero = native grid (grid fields
                               ! over the whole cell) or not yet generated (sampled fields)
     integer :: iregion_ap = iso_region_cell ! applied region mode (read only as whole-cell vs not)
     real*8 :: rgn_x_ap(3,0:3) = 0d0 ! applied region coordinates (user units per mode; the
                                     ! cell-frame box is derived at build time so it tracks cell edits)
     integer :: niso = 0 ! number of isosurfaces
     type(iso_slot), allocatable :: slot(:) ! the isosurfaces (niso of them)
     integer :: isel = 1 ! isosurface whose options are shown under the table in the editor
     real*8 :: timelastapply_grid = -1d0 ! time the grid + region were last applied (vs time_built)
     integer :: ifield_built = -1 ! field id when the meshes were built (-1 means never)
     integer :: imosel_built = 0 ! MO selector when the samples were taken
     integer :: imoidx_built = 0 ! MO index when the samples were taken
     integer :: fieldgen_built = -1 ! system field-set generation when the meshes were built
     logical :: per0_built = .false. ! whether the built meshes are periodic (whole cell of a crystal,
                                     ! non-partial data); gates the periodic replication and its editor UI
     logical :: outdomain = .false. ! some samples fell outside the field's domain (zeroed) in the last build
     real*8 :: time_built = -1d0 ! time the meshes were built (vs sysc timelastchange)
     real*8 :: frange(2) = (/1d0,-1d0/) ! min/max of the data backing the mesh; invalid if frange(1) > frange(2)
     ! histogram of the data backing the mesh, ready to plot as a staircase (2 points per bin):
     ! hist_x in field units, hist_y = number of points in the bin. Stamped with the build, like frange
     ! one staircase per axis scale (hscale_*): the bins are uniform in the transformed
     ! variable, so the bars match whichever axis the user selects
     real(c_double) :: hist_x(2*iso_nhist,hscale_num) = 0._c_double
     real(c_double) :: hist_y(2*iso_nhist,hscale_num) = 0._c_double
     logical :: hist_have(hscale_num) = .false. ! which scales the data allow
     real*8 :: hist_range(2,hscale_num) = 0d0 ! bin limits of each, in transformed units
     integer :: nhist = 0 ! number of staircase points (0 = no histogram yet)
     integer :: hist_xscale = 0 ! user: value-axis scale (0 = follow the suggested default)
     integer :: hist_yscale = hscale_log ! user: count-axis scale
     ! fraction of the integral of |f| (hist_q) and of the whole box volume (hist_v) at or above
     ! the lower edge of each bin, for the readout of what the current isovalue encloses. Binned
     ! logarithmically in |f| over hist_cumrange (log10 limits), independently of the staircase
     ! above: the display bins follow the axis on screen, and near zero those are too coarse
     real*8 :: hist_q(iso_nhist) = 0d0
     real*8 :: hist_v(iso_nhist) = 0d0
     real*8 :: hist_cumrange(2) = 0d0
     real*8, allocatable :: ff(:,:,:) ! cached field samples (non-grid fields; keyed by ifield_built and the applied-grid stamp)
   contains
     procedure :: set_field => iso_set_field ! select a field: default isovalue + grid level + applied dims
     procedure :: add_iso => iso_add_iso ! add an isosurface (isovalue and color chosen if not given)
     procedure :: del_iso => iso_del_iso ! remove an isosurface
     procedure :: apply_grid => iso_apply_grid ! commit a staged grid + region as the applied state
     procedure :: grid_isapplied => iso_grid_isapplied ! whether a staged grid + region is already applied
     procedure :: sampled_box => iso_sampled_box ! the box sampled by the applied state (domain policy)
     procedure :: isgenerated => iso_isgenerated ! whether the applied state describes a generated isosurface
     procedure :: stamp_histogram => iso_stamp_histogram ! recompute the value range and histogram from data
     procedure :: stamp_built => iso_stamp_built ! stamp the keys that say the samples are current
     procedure :: set_samples => iso_set_samples ! install externally supplied field samples
  end type rep_isosurface
  public :: rep_isosurface

  !> Representation: objects to draw on the scene
  type representation
     ! main variables
     logical :: isinit = .false. ! whether the representation has been initialized
     logical :: shown = .false. ! true if the representation is currently shown
     integer :: type = reptype_none ! type of representation (atoms, cell,...)
     integer :: flavor = repflavor_unknown ! flavor of the representation
     integer :: id ! system ID
     integer :: iord = 0 ! representation order integer in menu
     character(kind=c_char,len=:), allocatable :: name ! name of the representation
     ! transient-representation identity (owner > 0); unused in regular representations
     integer :: owner = 0 ! ID of the producer (owner) window; 0 = not transient
     integer :: itag = 0 ! producer-local content tag, for dedup/update-in-place
     logical :: armed = .false. ! re-armed this frame by the producer (otherwise reaped)
     ! per-object option groups
     type(rep_selection) :: sel ! which part of the system is drawn
     type(rep_atoms) :: atoms ! atom display options
     type(rep_bonds) :: bonds ! bond display options
     type(rep_labels) :: labels ! label display options
     type(rep_mols) :: mols ! per-molecule display options
     type(rep_unitcell) :: uc ! unit cell display options
     type(rep_axes) :: axes ! cartesian/crystallographic axes options
     type(rep_rotaxis) :: rotaxis ! rotation axis options
     type(rep_shapes) :: shapes ! geometric shapes options
     type(rep_symelem) :: symelem ! symmetry element options
     type(rep_poly) :: poly ! coordination polyhedra options
     type(rep_text) :: text ! text annotation options
     type(rep_measure) :: measure ! measurement options
     type(rep_isosurface) :: iso ! isosurface options
   contains
     procedure :: init => representation_init
     procedure :: set_defaults => representation_set_defaults
     procedure :: end => representation_end
     procedure :: update => update_styles
     procedure :: add_draw_elements
     procedure :: reset_all_styles
  end type representation
  public :: representation

  public :: iso_default_isovalue
  public :: iso_grid_size
  public :: iso_isgridfield
  public :: iso_region_to_box
  public :: iso_region_seed
  public :: iso_region_point_from_cart
  public :: iso_estimate_cost

  ! module procedure interfaces
  interface
     module function rep_selection_ncells(sel,nc) result(n)
       class(rep_selection), intent(in) :: sel
       integer, intent(in) :: nc(3)
       integer :: n(3)
     end function rep_selection_ncells
     module function iso_default_isovalue(isys,ifield) result(isoval)
       integer, intent(in) :: isys
       integer, intent(in) :: ifield
       real*8 :: isoval
     end function iso_default_isovalue
     module function iso_grid_size(isys,ilevel,ncustom,capped,ifield,box) result(n)
       integer, intent(in) :: isys
       integer, intent(in) :: ilevel
       integer, intent(in) :: ncustom(3)
       logical, intent(out), optional :: capped
       integer, intent(in), optional :: ifield
       real*8, intent(in), optional :: box(3,0:3)
       integer :: n(3)
     end function iso_grid_size
     module subroutine iso_region_to_box(isys,iregion,x,box,ok)
       integer, intent(in) :: isys
       integer, intent(in) :: iregion
       real*8, intent(in) :: x(3,0:3)
       real*8, intent(out) :: box(3,0:3)
       logical, intent(out), optional :: ok
     end subroutine iso_region_to_box
     module subroutine iso_region_seed(isys,iregion,x)
       integer, intent(in) :: isys
       integer, intent(in) :: iregion
       real*8, intent(out) :: x(3,0:3)
     end subroutine iso_region_seed
     module subroutine iso_region_point_from_cart(isys,iregion,xc,x)
       integer, intent(in) :: isys
       integer, intent(in) :: iregion
       real*8, intent(in) :: xc(3)
       real*8, intent(out) :: x(3)
     end subroutine iso_region_point_from_cart
     module function iso_estimate_cost(isys,ifield,iregion,x,n) result(secs)
       integer, intent(in) :: isys
       integer, intent(in) :: ifield
       integer, intent(in) :: iregion
       real*8, intent(in) :: x(3,0:3)
       integer, intent(in) :: n(3)
       real*8 :: secs
     end function iso_estimate_cost
     module subroutine iso_apply_grid(iso,n,iregion,x)
       class(rep_isosurface), intent(inout) :: iso
       integer, intent(in) :: n(3)
       integer, intent(in) :: iregion
       real*8, intent(in) :: x(3,0:3)
     end subroutine iso_apply_grid
     module function iso_grid_isapplied(iso,n,iregion,x) result(isap)
       class(rep_isosurface), intent(in) :: iso
       integer, intent(in) :: n(3)
       integer, intent(in) :: iregion
       real*8, intent(in) :: x(3,0:3)
       logical :: isap
     end function iso_grid_isapplied
     module subroutine iso_sample_domain(isys,ifield,iregion,x,n,xmat,x0c,cmat,per0,pereval,ok)
       integer, intent(in) :: isys
       integer, intent(in) :: ifield
       integer, intent(in) :: iregion
       real*8, intent(in) :: x(3,0:3)
       integer, intent(in) :: n(3)
       real*8, intent(out) :: xmat(3,3)
       real*8, intent(out) :: x0c(3)
       real*8, intent(out) :: cmat(3,3)
       logical, intent(out) :: per0
       logical, intent(out) :: pereval
       logical, intent(out) :: ok
     end subroutine iso_sample_domain
     module subroutine iso_sampled_box(iso,isys,xmat,x0c,cmat,per0,pereval,ok)
       class(rep_isosurface), intent(in) :: iso
       integer, intent(in) :: isys
       real*8, intent(out) :: xmat(3,3)
       real*8, intent(out) :: x0c(3)
       real*8, intent(out) :: cmat(3,3)
       logical, intent(out) :: per0
       logical, intent(out) :: pereval
       logical, intent(out) :: ok
     end subroutine iso_sampled_box
     module subroutine iso_stamp_histogram(iso,ff)
       class(rep_isosurface), intent(inout) :: iso
       real*8, intent(in) :: ff(:,:,:)
     end subroutine iso_stamp_histogram
     module subroutine iso_stamp_built(iso,isys)
       class(rep_isosurface), intent(inout) :: iso
       integer, intent(in) :: isys
     end subroutine iso_stamp_built
     module subroutine iso_set_samples(iso,isys,ff,outdomain)
       class(rep_isosurface), intent(inout) :: iso
       integer, intent(in) :: isys
       real*8, intent(in) :: ff(:,:,:)
       logical, intent(in) :: outdomain
     end subroutine iso_set_samples
     module function iso_isgenerated(iso,isys) result(gen)
       class(rep_isosurface), intent(in) :: iso
       integer, intent(in) :: isys
       logical :: gen
     end function iso_isgenerated
     module function iso_isgridfield(isys,ifield) result(isg)
       integer, intent(in) :: isys
       integer, intent(in) :: ifield
       logical :: isg
     end function iso_isgridfield
     module subroutine iso_set_field(iso,isys,ifield)
       class(rep_isosurface), intent(inout) :: iso
       integer, intent(in) :: isys
       integer, intent(in) :: ifield
     end subroutine iso_set_field
     module subroutine iso_add_iso(iso,isoval,rgb)
       class(rep_isosurface), intent(inout) :: iso
       real*8, intent(in), optional :: isoval
       real(c_float), intent(in), optional :: rgb(3)
     end subroutine iso_add_iso
     module subroutine iso_del_iso(iso,i)
       class(rep_isosurface), intent(inout) :: iso
       integer, intent(in) :: i
     end subroutine iso_del_iso
     module subroutine representation_init(r,isys,itype,flavor,icount)
       class(representation), intent(inout) :: r
       integer, intent(in) :: isys
       integer, intent(in) :: itype
       integer, intent(in) :: flavor
       integer, intent(inout) :: icount(0:repflavor_NUM)
     end subroutine representation_init
     module subroutine representation_set_defaults(r,itype)
       class(representation), intent(inout) :: r
       integer, intent(in) :: itype
     end subroutine representation_set_defaults
     module subroutine measurement_item_set_defaults(it,n)
       class(measurement_item), intent(inout) :: it
       integer, intent(in) :: n
     end subroutine measurement_item_set_defaults
     module subroutine measurement_item_copy_style(dst,src)
       class(measurement_item), intent(inout) :: dst
       type(measurement_item), intent(in) :: src
     end subroutine measurement_item_copy_style
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
     module subroutine symelem_style_alloc(d,nop)
       class(symelem_style), intent(inout) :: d
       integer, intent(in) :: nop
     end subroutine symelem_style_alloc
     module subroutine symelem_style_reset(d,r)
       class(symelem_style), intent(inout) :: d
       type(representation), intent(in) :: r
     end subroutine symelem_style_reset
     module subroutine symelem_style_end(d)
       class(symelem_style), intent(inout) :: d
     end subroutine symelem_style_end
  end interface

end module representations
