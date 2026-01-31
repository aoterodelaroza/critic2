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
  use shapes, only: dl_sphere, dl_cylinder, dl_string, scene_objects
  implicit none

  private

  ! default parameters for the representations (all distances in bohr)
  !--> atoms
  real(c_float), parameter, public :: atomborder_def = 0.1_c_float ! atom border
  real(c_float), parameter, public :: atomcovradscale_def = 0.7_c_float ! atomic radius scale factor (covalent)
  real(c_float), parameter, public :: atomvdwradscale_def = 1.0_c_float ! atomic radius scale factor (vdw)
  real(c_float), parameter, public :: atomrad_licorice_def = 0.23_c_float ! atomic radius value (licorice)
  !--> bonds
  real(c_float), parameter, public :: bondrad_def = 0.35_c_float ! bond radius
  real(c_float), parameter, public :: bondrad_licorice_def = 0.50_c_float ! bond radius (licorice)
  real(c_float), parameter, public :: bondborder_def = 0.1_c_float ! bond border
  real(c_float), parameter, public :: bondborder_stickflav_def = 0.05_c_float ! bond border (stick flavor)
  !--> unit cell
  real(c_float), parameter, public :: uc_radius_def = 0.15_c_float ! radius of sticks
  real(c_float), parameter, public :: uc_radiusinner_def = 0.15_c_float ! radius of inner sticks
  real(c_float), parameter, public :: uc_innersteplen_def = 2.0_c_float ! length of stipple

  !> Draw style for atoms (geometry-dependent parameters)
  type atom_geom_style
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer :: type ! atom style type (attlisttype_* in systems module)
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     real(c_float), allocatable :: rgb(:,:) ! color (3,ntype)
     real(c_float), allocatable :: rad(:) ! radius (ntype)
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
     real(c_float), allocatable :: scale_rad(:) ! scale radius (ntype)
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

  ! types of representations
  integer, parameter, public :: reptype_none = 0
  integer, parameter, public :: reptype_atoms = 1
  integer, parameter, public :: reptype_unitcell = 2
  integer, parameter, public :: reptype_NUM = 2

  ! representation flavors
  integer, parameter, public :: repflavor_unknown = 0
  integer, parameter, public :: repflavor_atoms_ballandstick = 1
  integer, parameter, public :: repflavor_atoms_sticks = 2
  integer, parameter, public :: repflavor_atoms_licorice = 3
  integer, parameter, public :: repflavor_atoms_vdwcontacts = 4
  integer, parameter, public :: repflavor_atoms_hbonds = 5
  integer, parameter, public :: repflavor_atoms_criticalpoints = 6
  integer, parameter, public :: repflavor_unitcell_basic = 7
  integer, parameter, public :: repflavor_NUM = 7

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
     ! global parameters
     integer(c_int) :: pertype ! periodicity control: 0=none, 1=auto, 2=manual
     integer(c_int) :: ncell(3) ! number of unit cells drawn
     real*8 :: origin(3) ! unit cell, origin shift
     real*8 :: tshift(3) ! origin of the unit cell display region
     ! atoms, bonds, labels
     character(kind=c_char,len=:), allocatable :: filter ! filter for the representation
     character(kind=c_char,len=:), allocatable :: errfilter ! filter error
     logical :: border ! draw atoms at the border of the unit cell
     logical :: onemotif ! draw connected molecules
     !--> atoms
     logical :: atoms_display ! whether to draw the atoms
     type(atom_geom_style) :: atom_style ! atom styles
     integer(c_int) :: atom_radii_type ! option to reset radii: 0=covalent,1=vdw,2=constant
     real(c_float) :: atom_radii_scale ! reset radii, scale factor
     real(c_float) :: atom_radii_value ! reset radii, constant value (default: same as bond radius for licorice)
     integer(c_int) :: atom_color_type ! option to reset colors: 0=current,1=jmlcol,2=jmlcol2
     real(c_float) :: atom_border_size ! atom border size
     real(c_float) :: atom_border_rgb(3) ! atom border color
     !--> bonds
     logical :: bonds_display ! whether to draw the bonds
     type(bond_geom_style) :: bond_style ! bond styles
     integer(c_int) :: bond_distancetype ! selector for distance type (0=factor,1=range)
     real*8 :: bond_dmin ! distance limits (angstrom)
     real*8 :: bond_dmax ! distance limits (angstrom)
     real*8 :: bond_bfmin ! bondfactor limits
     real*8 :: bond_bfmax ! bondfactor limits
     integer(c_int) :: bond_radtype(2) ! radii type for min and max (0=covalent,1=vdw)
     integer(c_int) :: bond_color_style ! bond style (0=single color, 1=two colors)
     real(c_float) :: bond_rad ! radius
     real(c_float) :: bond_border_size ! bond border size
     real(c_float) :: bond_border_rgb(3) ! bond color
     real(c_float) :: bond_rgb(3) ! color
     integer(c_int) :: bond_order ! order (0=dashed,1=single,2=double,etc.)
     integer(c_int) :: bond_imol ! molecular connections (0=any,1=intramol,2=intermol)
     logical :: bond_bothends ! if true, both atoms need to be drawn to draw the bond
     !--> labels
     logical :: labels_display ! whether to draw the labels
     type(label_geom_style) :: label_style ! bond styles
     integer(c_int) :: label_type ! 0=atom-symbol, 1=atom-name, 2=cel-atom, 3=cel-atom+lvec, 4=neq-atom, 5=spc, 6=Z, 7=mol, 8=wyckoff
     real(c_float) :: label_scale ! scale for the labels
     real(c_float) :: label_rgb(3) ! color of the labels
     logical :: label_const_size ! whether labels scale with objects or are constant size
     real(c_float) :: label_offset(3) ! offset of the label
     !--> molecules
     type(mol_geom_style) :: mol_style ! molecule styles
     ! unit cell
     logical :: uc_inner ! unit cell, display inner cylinders
     logical :: uc_coloraxes ! unit cell, color the axes (x=red,y=green,z=blue)
     logical :: uc_vaccutsticks ! unit cell, cut sticks in systems with vacuum
     real(c_float) :: uc_radius ! unit cell cylinder radius
     real(c_float) :: uc_radiusinner ! unit cell cylinder radius (inner)
     real(c_float) :: uc_rgb(3) ! unit cell cylinder colors
     real(c_float) :: uc_innersteplen ! number of subdivisions for the inner sticks
     logical :: uc_innerstipple ! stippled lines for the inner lines
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
     module subroutine representation_init(r,isys,irep,itype,style,flavor,icount)
       class(representation), intent(inout) :: r
       integer, intent(in) :: isys
       integer, intent(in) :: irep
       integer, intent(in) :: itype
       integer, intent(in) :: style
       integer, intent(in) :: flavor
       integer, intent(inout) :: icount(0:repflavor_NUM)
     end subroutine representation_init
     module subroutine representation_set_defaults(r,style,itype)
       class(representation), intent(inout) :: r
       integer, intent(in) :: style
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
  end interface

end module representations
