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

  ! default bond radius and atom border
  real(c_float), parameter, public :: bond_rad_def = 0.35_c_float
  real(c_float), parameter, public :: atomborder_def = 0.1_c_float

  !> Draw style for atoms
  type draw_style_atom
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer :: type ! atom style type: 0=species,1=nneq,2=cell
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     real(c_float), allocatable :: rgb(:,:) ! color (3,ntype)
     real(c_float), allocatable :: rad(:) ! radius (ntype)
     real(c_float) :: border_size = 0._c_float ! border size
     real(c_float) :: rgbborder(3) ! border color
   contains
     procedure :: reset => reset_atom_style
     procedure :: reset_colors => reset_colors_atom_style
  end type draw_style_atom
  public :: draw_style_atom

  !> Draw style for molecules
  type draw_style_molecule
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     real(c_float), allocatable :: tint_rgb(:,:) ! tint color (3,ntype)
     real(c_float), allocatable :: scale_rad(:) ! scale radius (ntype)
   contains
     procedure :: reset => reset_mol_style
  end type draw_style_molecule
  public :: draw_style_molecule

  !> Draw style for bonds
  type draw_style_bond
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     logical :: isdef = .true. ! whether this is using the system's neighbor star
     ! temporary storage for the edit object/bonds tab, global options
     integer(c_int) :: distancetype_g ! selector for distance type (0=factor,1=range)
     real(c_float) :: dmin_g, dmax_g ! distance limits (angstrom)
     real(c_float) :: bfmin_g, bfmax_g ! bondfactor limits
     integer(c_int) :: radtype_g(2) ! radii type for min and max (0=covalent,1=vdw)
     integer(c_int) :: style_g ! bond style (0=single color, 1=two colors)
     real(c_float) :: rad_g ! radius
     real(c_float) :: border_g ! bond border
     real(c_float) :: rgbborder_g(3) ! bond color
     real(c_float) :: rgb_g(3) ! color
     integer(c_int) :: order_g ! order (0=dashed,1=single,2=double,etc.)
     integer(c_int) :: imol_g ! molecular connections (0=any,1=intramol,2=intermol)
     logical :: bothends_g ! if true, both atoms need to be drawn to draw the bond
     logical, allocatable :: shown_g(:,:) ! by-species bond shown flags (nspc,nspc)
     ! the bond information
     type(neighstar), allocatable :: nstar(:) ! the neighbor star
   contains
     procedure :: generate_neighstars
     procedure :: copy_neighstars_from_system
     procedure :: reset => reset_bond_style
  end type draw_style_bond
  public :: draw_style_bond

  !> Draw style for labels
  type draw_style_label
     logical :: isinit = .false. ! whether the style is intialized
     real*8 :: timelastreset = 0d0 ! time the style was last reset
     integer(c_int) :: style = 0 ! 0=atom-symbol, 1=atom-name, 2=cel-atom, 3=cel-atom+lvec, 4=neq-atom, 5=spc, 6=Z, 7=mol, 8=wyckoff
     real(c_float) :: scale ! scale for the labels
     real(c_float) :: rgb(3) ! color of the labels
     logical :: const_size ! whether labels scale with objects or are constant size
     real(c_float) :: offset(3) ! offset of the label
     integer :: ntype = 0 ! number of entries in the style type (atoms or molecules)
     logical, allocatable :: shown(:) ! whether it is shown (ntype)
     character*32, allocatable :: str(:) ! text
   contains
     procedure :: reset => reset_label_style
  end type draw_style_label
  public :: draw_style_label

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
  integer, parameter, public :: repflavor_unitcell_basic = 6
  integer, parameter, public :: repflavor_NUM = 6

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
     integer(c_int) :: atom_color_reset_type = 0 ! option to reset colors: 0=current,1=jmlcol,2=jmlcol2
     type(draw_style_atom) :: atom_style ! atom styles
     type(draw_style_molecule) :: mol_style ! molecule styles
     type(draw_style_bond) :: bond_style ! bond styles
     type(draw_style_label) :: label_style ! bond styles
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
     procedure :: reset => representation_reset
     procedure :: end => representation_end
     procedure :: update => update_structure
     procedure :: add_draw_elements
     procedure :: reset_all_styles
  end type representation
  public :: representation

  ! module procedure interfaces
  interface
     module subroutine representation_init(r,isys,irep,itype,style,flavor,icount)
       class(representation), intent(inout), target :: r
       integer, intent(in) :: isys
       integer, intent(in) :: irep
       integer, intent(in) :: itype
       integer, intent(in) :: style
       integer, intent(in) :: flavor
       integer, intent(inout) :: icount(0:repflavor_NUM)
     end subroutine representation_init
     module subroutine representation_reset(r)
       class(representation), intent(inout), target :: r
     end subroutine representation_reset
     module subroutine representation_end(r)
       class(representation), intent(inout), target :: r
     end subroutine representation_end
     module subroutine update_structure(r)
       class(representation), intent(inout), target :: r
     end subroutine update_structure
     module subroutine add_draw_elements(r,nc,obj,doanim,iqpt,ifreq)
       class(representation), intent(inout), target :: r
       integer, intent(in) :: nc(3)
       type(scene_objects), intent(inout) :: obj
       logical, intent(in) :: doanim
       integer, intent(in) :: iqpt, ifreq
     end subroutine add_draw_elements
     module subroutine reset_all_styles(r)
       class(representation), intent(inout), target :: r
     end subroutine reset_all_styles
     ! draw_style_atom
     module subroutine reset_atom_style(d,isys)
       class(draw_style_atom), intent(inout), target :: d
       integer, intent(in) :: isys
     end subroutine reset_atom_style
     module subroutine reset_colors_atom_style(d,isys)
       class(draw_style_atom), intent(inout), target :: d
       integer, intent(in) :: isys
     end subroutine reset_colors_atom_style
     ! draw_style_molecule
     module subroutine reset_mol_style(d,isys)
       class(draw_style_molecule), intent(inout), target :: d
       integer, intent(in) :: isys
     end subroutine reset_mol_style
     ! draw_style_bond
     module subroutine generate_neighstars(d,isys)
       class(draw_style_bond), intent(inout), target :: d
       integer, intent(in) :: isys
     end subroutine generate_neighstars
     module subroutine copy_neighstars_from_system(d,isys)
       class(draw_style_bond), intent(inout), target :: d
       integer, intent(in) :: isys
     end subroutine copy_neighstars_from_system
     module subroutine reset_bond_style(d,isys,flavor)
       class(draw_style_bond), intent(inout), target :: d
       integer, intent(in) :: isys
       integer, intent(in), optional :: flavor
     end subroutine reset_bond_style
     ! draw_style_label
     module subroutine reset_label_style(d,isys)
       class(draw_style_label), intent(inout), target :: d
       integer, intent(in) :: isys
     end subroutine reset_label_style
  end interface

end module representations
