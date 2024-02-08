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

! Structure class and routines for basic crystallography computations
module crystalmod
  use environmod, only: environ
  use spglib, only: SpglibDataset
  use types, only: neqatom, celatom, neighstar, species
  use fragmentmod, only: fragment
  use param, only: maxzat0, mlen
  use types, only: thread_info
  implicit none

  private

  ! Holohedry identifier
  integer, parameter, public :: holo_unk = 0 ! unknown
  integer, parameter, public :: holo_tric = 1 ! triclinic
  integer, parameter, public :: holo_mono = 2 ! monoclinic
  integer, parameter, public :: holo_ortho = 3 ! orthorhombic
  integer, parameter, public :: holo_tetra = 4 ! tetragonal
  integer, parameter, public :: holo_trig = 5 ! trigonal
  integer, parameter, public :: holo_hex = 6 ! hexagonal
  integer, parameter, public :: holo_cub = 7 ! cubic
  character(len=12), parameter, public :: holo_string(0:7) = (/ &
     "unknown     ","triclinic   ","monoclinic  ","orthorhombic",&
     "tetragonal  ","trigonal    ","hexagonal   ","cubic       "/)

  ! Laue class identifier
  integer, parameter, public :: laue_unk = 0 ! unknown
  integer, parameter, public :: laue_1 = 1 ! -1
  integer, parameter, public :: laue_2m = 2 ! 2/m
  integer, parameter, public :: laue_mmm = 3 ! mmm
  integer, parameter, public :: laue_4m = 4 ! 4/m
  integer, parameter, public :: laue_4mmm = 5 ! 4/mmm
  integer, parameter, public :: laue_3 = 6 ! -3
  integer, parameter, public :: laue_3m = 7 ! -3m
  integer, parameter, public :: laue_6m = 8 ! 6/m
  integer, parameter, public :: laue_6mmm = 9 ! 6/mmm
  integer, parameter, public :: laue_m3 = 10 ! m-3
  integer, parameter, public :: laue_m3m = 11 ! m-3m
  character(len=12), parameter, public :: laue_string(0:11) = (/ &
     "unknown","-1     ","2/m    ","mmm    ","4/m    ","4/mmm  ",&
     "-3     ","-3m    ","6/m    ","6/mmm  ","m-3    ","m-3m   "/)

  !> The crystal class. A crystal contains the structural information for the
  !> system, and it can be an actual crystal (%ismolecule=.false.) or a molecule
  !> (%ismolecule=.true.) embedded in a large cell. When use in combination
  !> with the system class, the crystal needs to be passed as pointer to the
  !> system, and therefore needs to be declared as TARGET.
  type crystal
     ! Initialization flags
     logical :: isinit = .false. !< has the crystal structure been initialized?
     integer :: havesym = 0 !< was the symmetry determined? (0 - nosym, 1 - full)
     logical :: isewald = .false. !< do we have the data for ewald's sum?

     ! file name for the occasional critic2 trick
     character(len=mlen) :: file

     !! Initialization level: isinit !!
     ! species list
     integer :: nspc = 0 !< Number of species
     type(species), allocatable :: spc(:) !< Species
     ! non-equivalent atoms list
     integer :: nneq = 0 !< Number of non-equivalent atoms
     type(neqatom), allocatable :: at(:) !< Non-equivalent atom array
     ! complete atoms list
     integer :: ncel = 0 !< Number of atoms in the main cell
     type(celatom), allocatable :: atcel(:) !< List of atoms in the main cell
     ! cell and lattice metrics
     real*8 :: aa(3) !< cell lengths (bohr)
     real*8 :: bb(3) !< cell angles (degrees)
     real*8 :: omega !< unit cell volume
     real*8 :: gtensor(3,3) !< metric tensor (3,3)
     real*8 :: ar(3) !< reciprocal cell lengths
     real*8 :: grtensor(3,3) !< reciprocal metric tensor (3,3)
     ! crystallographic/cartesian conversion matrices and norm-2s
     real*8 :: m_x2c(3,3) !< input cell, crystallographic -> cartesian (m_x2c(:,i) are the lattice vectors)
     real*8 :: m_c2x(3,3) !< input cell, cartesian -> crystallographic
     real*8 :: m_xr2c(3,3) !< reduced cryst -> input cartesian
     real*8 :: m_c2xr(3,3) !< cartesian -> reduced cryst
     real*8 :: m_xr2x(3,3) !< reduced cryst -> input cryst
     real*8 :: m_x2xr(3,3) !< input cryst -> reduced cryst
     real*8 :: n2_x2c !< norm2 of input cell, crystallographic -> cartesian
     real*8 :: n2_c2x !< norm2 of input cell, cartesian -> crystallographic
     real*8 :: n2_xr2c !< norm2 of reduced cryst -> input cartesian
     real*8 :: n2_c2xr !< norm2 of cartesian -> reduced cryst
     real*8 :: n2_xr2x !< norm2 of reduced cryst -> input cryst
     real*8 :: n2_x2xr !< norm2 of input cryst -> reduced cryst
     ! space-group symmetry
     logical :: spgavail = .false. !< have spglib's symmetry?
     type(SpglibDataset) :: spg !< spglib's symmetry dataset
     integer :: neqv !< number of symmetry operations
     integer :: ncv  !< number of centering vectors
     real*8, allocatable :: cen(:,:) !< centering vectors
     real*8 :: rotm(3,4,48) !< symmetry operations
     ! variables for molecular systems
     logical :: ismolecule = .false. !< is it a molecule?
     real*8 :: molx0(3) !< centering vector for the molecule
     real*8 :: molborder(3) !< molecular cell border (cryst coords)
     ! wigner-seitz cell
     integer :: ws_nv !< number of vertices
     integer :: ws_nf !< number of facets
     integer :: ws_mnfv !< maximum number of vertices per facet
     integer, allocatable :: ws_ineighx(:,:) !< WS neighbor lattice points (cryst. coords.)
     real*8, allocatable :: ws_ineighc(:,:) !< WS neighbor lattice points (Cart. coords.)
     integer, allocatable :: ws_ineighxr(:,:) !< WS neighbor lattice points (del cell, cryst.)
     integer, allocatable :: ws_nside(:) !< number of sides of WS faces
     integer, allocatable :: ws_iside(:,:) !< sides of the WS faces
     real*8, allocatable :: ws_x(:,:) !< vertices of the WS cell (cryst. coords.)
     logical :: isortho !< is the cell orthogonal?
     logical :: isortho_del !< is the reduced cell orthogonal?
     ! core charges
     integer :: zpsp(maxzat0)

     !! Initialization level: isenv !!
     ! atomic environment of the cell
     type(environ) :: env
     integer :: nblock(3) ! number of environemt blocks
     integer,allocatable :: iblock0(:,:,:) ! starting atomic index for each block
     real*8 :: blockrmax ! radius of the largest sphere contained in a block
     real*8 :: blockomega ! volume of a block
     real*8 :: blockcv(3) ! cross products of the block lattice vectors

     !! Initialization level: isast !!
     ! asterisms
     type(neighstar), allocatable :: nstar(:) !< Neighbor stars
     integer :: nmol = 0 !< Number of molecules in the unit cell
     type(fragment), allocatable :: mol(:) !< Molecular fragments
     integer :: nlvac = 0 !< Number of vacuum lattice vectors
     integer :: lvac(3,2) !< Vacuum lattice vectors
     integer :: lcon(3,2) !< Connected lattice vectors
     real*8 :: vaclength(3) !< Vacuum length in bohr along the crystallographic axes (a,b,c)
     ! variables for 3d molecular crystals
     logical :: ismol3d !< Is this a 3d molecular crystal?
     integer, allocatable :: idxmol(:) !< -1: mol is fractional, 0: sym. unique, >0 index for nneq mol.
     integer, allocatable :: idatcelmol(:) !< cell atom i belongs to idatcelmol(i) molecule

     !! Initialization level: isewald !!
     ! ewald data
     real*8 :: rcut, hcut, eta, qsum, q2sum
     integer :: lrmax(3), lhmax(3)

   contains
     ! construction, destruction, initialization (proc)
     procedure :: init => struct_init !< Allocate arrays and nullify variables
     procedure :: end => struct_end !< Deallocate arrays and nullify variables
     procedure :: struct_new !< Initialize the structure from a crystal seed

     ! basic crystallographic operations (proc)
     procedure :: x2c !< Convert input cryst. -> cartesian
     procedure :: c2x !< Convert input cartesian -> cryst.
     procedure :: xr2c !< Convert reduced cryst. -> cartesian
     procedure :: c2xr !< Convert cartesian -> reduced cryst.
     procedure :: xr2x !< Convert reduced cryst. -> input cryst.
     procedure :: x2xr !< Convert input cryst. -> reduced cryst.
     procedure :: distance !< Distance between points in crystallographic coordinates
     procedure :: distmatrix !< calculate the distance matrix (molecules only)
     procedure :: eql_distance !< Shortest distance between lattice-translated vectors
     procedure :: shortest !< Gives the lattice-translated vector with shortest length
     procedure :: are_close !< True if a vector is at a distance less than eps of another
     procedure :: are_lclose !< True if a vector is at a distance less than eps of all latice translations of another
     procedure :: nearest_atom !< Calculate the atom nearest to a given point
     procedure :: nearest_atom_grid !< return the nearest atom IDs for a uniform grid of points
     procedure :: identify_atom !< Identify an atom in the unit cell
     procedure :: identify_spc !< Identify species in a structure from a string

     ! atomic environments and distance calculations
     procedure :: build_env
     procedure :: list_near_atoms
     procedure :: nearest_atom_env

     ! molecular environments and neighbors (mols)
     procedure :: identify_fragment !< Build an atomic fragment of the crystal
     procedure :: identify_fragment_from_xyz !< Build a crystal fragment from an xyz file
     procedure :: fill_molecular_fragments !< Find the molecular fragments in the crystal
     procedure :: calculate_molecular_equivalence !< Calculate symmetry relations between molecules
     procedure :: listatoms_cells !< List all atoms in n cells (maybe w border)
     procedure :: listatoms_sphcub !< List all atoms in a sphere or cube
     procedure :: listmolecules !< List all molecules in the crystal

     ! complex operations: ewald, promolecular, etc. (complex)
     procedure :: calculate_ewald_cutoffs !< Calculate the cutoffs for Ewald's sum
     procedure :: ewald_energy !< electrostatic energy (Ewald)
     procedure :: ewald_pot !< electrostatic potential (Ewald)
     procedure :: promolecular !< calculate core and promolecular densities
     procedure :: promolecular_grid !< calculate core and promolecular densities on a grid
     procedure :: get_pack_ratio !< Calculate the packing ratio
     procedure :: vdw_volume !< Calculate the van der waals volume
     procedure :: get_kpoints !< k-point grid for a given rklength

     ! changes to the crystal/molecular structure (edit)
     procedure :: makeseed !< make a crystal seed from a crystal
     procedure :: newcell !< Change the unit cell and rebuild the crystal
     procedure :: cell_standard !< Transform the the standard cell (possibly primitive)
     procedure :: cell_niggli !< Transform to the Niggli primitive cell
     procedure :: cell_delaunay !< Transform to the Delaunay primitive cell
     procedure :: reorder_atoms !< reorder the atoms in the crystal/molecule
     procedure :: wholemols !< Re-assign atomic types to have an asymmetric unit with whole molecules
     procedure :: delete_atoms !< Delete a list of atoms

     ! symmetry (symmetry)
     procedure :: sitesymm !< Determine the local-symmetry group symbol for a point
     procedure :: symeqv  !< Calculate the symmetry-equivalent positions of a point
     procedure :: get_mult !< Multiplicity of a point
     procedure :: spglib_wrap !< Get the spg from the crystal geometry
     procedure :: spgtowyc !< Copy the Wyckoff positions to a crystal from an spg
     procedure :: calcsym !< Calculate the symmetry operations from the crystal geometry
     procedure :: clearsym !< Clear symmetry info and transform to a P1
     procedure :: checkgroup !< Check that the space group operations are consistent
     procedure :: getiws !< Calculate the IWS and its tetrahedra partition around a point

     ! powder diffraction and related calcs (powderproc)
     procedure :: powder !< Calculate the powder diffraction pattern
     procedure :: rdf !< Calculate the radial distribution function
     procedure :: amd !< Calculate the average minimum distances (Widdowson et al.)

     ! output routines (write)
     procedure :: report => struct_report !< Write lots of information about the crystal structure to uout
     procedure :: struct_report_symmetry !< Write symmmetry information
     procedure :: struct_report_symxyz !< Write sym. ops. in crystallographic notation to uout
     procedure :: struct_write_json !< Write a json object containing the crystal structure info

     ! structure writers (write)
     procedure :: write_simple_driver
     procedure :: write_mol
     procedure :: write_3dmodel
     procedure :: write_espresso
     procedure :: write_vasp
     procedure :: write_abinit
     procedure :: write_elk
     procedure :: write_gaussian
     procedure :: write_tessel
     procedure :: write_critic
     procedure :: write_cif
     procedure :: write_d12
     procedure :: write_res
     procedure :: write_escher
     procedure :: write_db
     procedure :: write_gulp
     procedure :: write_lammps
     procedure :: write_siesta_fdf
     procedure :: write_siesta_in
     procedure :: write_dftbp_hsd
     procedure :: write_dftbp_gen
     procedure :: write_pyscf
     procedure :: write_fhi
     procedure :: write_tinkerfrac
     procedure :: write_castep_cell

     ! grid writers (write)
     procedure :: writegrid_cube
     procedure :: writegrid_vasp
     procedure :: writegrid_xsf
  end type crystal
  public :: crystal

  ! other crystallography tools that are crystal-independent (symmetry)
  public :: search_lattice
  public :: pointgroup_info

  ! module procedure interfaces
  interface
     !xx! proc submodule
     module subroutine struct_init(c)
       class(crystal), intent(inout) :: c
     end subroutine struct_init
     module subroutine struct_end(c)
       class(crystal), intent(inout) :: c
     end subroutine struct_end
     module subroutine struct_new(c,seed,crashfail,noenv,ti)
       use crystalseedmod, only: crystalseed
       class(crystal), intent(inout) :: c
       type(crystalseed), intent(in) :: seed
       logical, intent(in) :: crashfail
       logical, intent(in), optional :: noenv
       type(thread_info), intent(in), optional :: ti
     end subroutine struct_new
     pure module function x2c(c,xx) result(res)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: xx(3)
       real*8 :: res(3)
     end function x2c
     pure module function c2x(c,xx) result(res)
       class(crystal), intent(in) :: c
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function c2x
     pure module function xr2c(c,xx) result(res)
       class(crystal), intent(in) :: c
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function xr2c
     pure module function c2xr(c,xx) result(res)
       class(crystal), intent(in) :: c
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function c2xr
     pure module function xr2x(c,xx) result(res)
       class(crystal), intent(in) :: c
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function xr2x
     pure module function x2xr(c,xx) result(res)
       class(crystal), intent(in) :: c
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function x2xr
     pure module function distance(c,x1,x2)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: x1(3)
       real*8, intent(in) :: x2(3)
       real*8 :: distance
     end function distance
     module subroutine distmatrix(c,d,inverse,conn)
       class(crystal), intent(in) :: c
       real*8, allocatable, intent(inout) :: d(:,:)
       logical, intent(in), optional :: inverse
       logical, intent(in), optional :: conn
     end subroutine distmatrix
     pure module function eql_distance(c,x1,x2)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: x1(3)
       real*8, intent(in) :: x2(3)
       real*8 :: eql_distance
     end function eql_distance
     pure module subroutine shortest(c,x,dist)
       class(crystal), intent(in) :: c
       real*8, intent(inout) :: x(3)
       real*8, intent(out) :: dist
     end subroutine shortest
     module function are_close(c,x0,x1,eps,dd)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: x0(3), x1(3)
       real*8, intent(in) :: eps
       real*8, intent(out), optional :: dd
       logical :: are_close
     end function are_close
     module function are_lclose(c,x0,x1,eps,dd)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: x0(3), x1(3)
       real*8, intent(in) :: eps
       real*8, intent(out), optional :: dd
       logical :: are_lclose
     end function are_lclose
     module subroutine nearest_atom(c,xp,icrd,nid,dist,distmax,lvec,cidx0,idx0,is0,nozero)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: xp(3)
       integer, intent(in) :: icrd
       integer, intent(out) :: nid
       real*8, intent(out) :: dist
       real*8, intent(in), optional :: distmax
       integer, intent(out), optional :: lvec(3)
       integer, intent(in), optional :: cidx0
       integer, intent(in), optional :: idx0
       integer, intent(in), optional :: is0
       logical, intent(in), optional :: nozero
     end subroutine nearest_atom
     module subroutine nearest_atom_grid(c,n,idg)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: n(3)
       integer, allocatable, intent(inout) :: idg(:,:,:)
     end subroutine nearest_atom_grid
     module function identify_atom(c,x0,icrd,lvec,dist,distmax)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: x0(3)
       integer, intent(in) :: icrd
       integer, intent(out), optional :: lvec(3)
       real*8, intent(out), optional :: dist
       real*8, intent(in), optional :: distmax
       integer :: identify_atom
     end function identify_atom
     module function identify_spc(c,str) result(res)
       use crystalseedmod, only: crystalseed
       class(crystal), intent(inout) :: c
       character*(*), intent(in) :: str
       integer :: res
     end function identify_spc
     module subroutine build_env(c)
       class(crystal), intent(inout) :: c
     end subroutine build_env
     module subroutine list_near_atoms(c,xp,icrd,sorted,nat,eid,dist,lvec,ishell0,up2d,&
        up2dsp,up2dcidx,up2sh,up2n,nid0,id0,iz0,ispc0,nozero)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: xp(3)
       integer, intent(in) :: icrd
       logical, intent(in) :: sorted
       integer, intent(out) :: nat
       integer, allocatable, intent(inout), optional :: eid(:)
       real*8, allocatable, intent(inout), optional :: dist(:)
       integer, allocatable, intent(inout), optional :: lvec(:,:)
       integer, allocatable, intent(inout), optional :: ishell0(:)
       real*8, intent(in), optional :: up2d
       real*8, intent(in), optional :: up2dsp(1:c%nspc,2)
       real*8, intent(in), optional :: up2dcidx(1:c%ncel)
       integer, intent(in), optional :: up2sh
       integer, intent(in), optional :: up2n
       integer, intent(in), optional :: nid0
       integer, intent(in), optional :: id0
       integer, intent(in), optional :: iz0
       integer, intent(in), optional :: ispc0
       logical, intent(in), optional :: nozero
     end subroutine list_near_atoms
     module subroutine nearest_atom_env(c,xp,icrd,nid,dist,distmax,lvec,nid0,id0,iz0,ispc0,nozero)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: xp(3)
       integer, intent(in) :: icrd
       integer, intent(out) :: nid
       real*8, intent(out) :: dist
       real*8, intent(in), optional :: distmax
       integer, intent(out), optional :: lvec(3)
       integer, intent(in), optional :: nid0
       integer, intent(in), optional :: id0
       integer, intent(in), optional :: iz0
       integer, intent(in), optional :: ispc0
       logical, intent(in), optional :: nozero
     end subroutine nearest_atom_env
     module function identify_fragment(c,nat,x0) result(fr)
       class(crystal), intent(in) :: c
       integer, intent(in) :: nat
       real*8, intent(in) :: x0(3,nat)
       type(fragment) :: fr
     end function identify_fragment
     module function identify_fragment_from_xyz(c,file,errmsg,ti) result(fr)
       class(crystal), intent(in) :: c
       character*(*) :: file
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
       type(fragment) :: fr
     end function identify_fragment_from_xyz
     module subroutine fill_molecular_fragments(c)
       class(crystal), intent(inout) :: c
     end subroutine fill_molecular_fragments
     module subroutine calculate_molecular_equivalence(c)
       class(crystal), intent(inout) :: c
     end subroutine calculate_molecular_equivalence
     module function listatoms_cells(c,nx,doborder) result(fr)
       class(crystal), intent(in) :: c
       integer, intent(in) :: nx(3)
       logical, intent(in) :: doborder
       type(fragment) :: fr
     end function listatoms_cells
     module function listatoms_sphcub(c,rsph,xsph,rcub,xcub) result(fr)
       class(crystal), intent(in) :: c
       real*8, intent(in), optional :: rsph, xsph(3)
       real*8, intent(in), optional :: rcub, xcub(3)
       type(fragment) :: fr
     end function listatoms_sphcub
     module subroutine listmolecules(c,fri,nfrag,fr,isdiscrete)
       class(crystal), intent(inout) :: c
       type(fragment), intent(in) :: fri
       integer, intent(out) :: nfrag
       type(fragment), intent(out), allocatable :: fr(:)
       logical, intent(out), allocatable :: isdiscrete(:)
     end subroutine listmolecules
     module subroutine calculate_ewald_cutoffs(c)
       class(crystal), intent(inout) :: c
     end subroutine calculate_ewald_cutoffs
     module function ewald_energy(c) result(ewe)
       class(crystal), intent(inout) :: c
       real*8 :: ewe
     end function ewald_energy
     module function ewald_pot(c,x)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: x(3)
       real*8 :: ewald_pot
     end function ewald_pot
     module subroutine promolecular(c,x0,icrd,f,fp,fpp,nder,zpsp,fr)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: x0(3)
       integer, intent(in) :: icrd
       real*8, intent(out) :: f
       real*8, intent(out) :: fp(3)
       real*8, intent(out) :: fpp(3,3)
       integer, intent(in) :: nder
       integer, intent(in), optional :: zpsp(:)
       type(fragment), intent(in), optional :: fr
     end subroutine promolecular
     module subroutine promolecular_grid(c,f,n,zpsp,fr)
       use grid3mod, only: grid3
       class(crystal), intent(in) :: c
       type(grid3), intent(out) :: f
       integer, intent(in) :: n(3)
       integer, intent(in), optional :: zpsp(:)
       type(fragment), intent(in), optional :: fr
     end subroutine promolecular_grid
     module function get_pack_ratio(c) result (px)
       class(crystal), intent(inout) :: c
       real*8 :: px
     end function get_pack_ratio
     module function vdw_volume(c,relerr,rtable) result(vvdw)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: relerr
       real*8, intent(in), optional :: rtable(:)
       real*8 :: vvdw
     end function vdw_volume
     module subroutine get_kpoints(c,rk,nk)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: rk
       integer, intent(out) :: nk(3)
     end subroutine get_kpoints
     module subroutine makeseed(c,seed,copysym)
       use crystalseedmod, only: crystalseed
       class(crystal), intent(in) :: c
       type(crystalseed), intent(out) :: seed
       logical, intent(in) :: copysym
     end subroutine makeseed
     module subroutine newcell(c,x00,t0,nnew,xnew,isnew,noenv,ti)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: x00(3,3)
       real*8, intent(in), optional :: t0(3)
       integer, intent(in), optional :: nnew
       real*8, intent(in), optional :: xnew(:,:)
       integer, intent(in), optional :: isnew(:)
       logical, intent(in), optional :: noenv
       type(thread_info), intent(in), optional :: ti
     end subroutine newcell
     module function cell_standard(c,toprim,doforce,refine,noenv,ti) result(x0)
       class(crystal), intent(inout) :: c
       logical, intent(in) :: toprim
       logical, intent(in) :: doforce
       logical, intent(in) :: refine
       logical, intent(in), optional :: noenv
       type(thread_info), intent(in), optional :: ti
       real*8 :: x0(3,3)
     end function cell_standard
     module function cell_niggli(c,noenv,ti) result(x0)
       class(crystal), intent(inout) :: c
       logical, intent(in), optional :: noenv
       type(thread_info), intent(in), optional :: ti
       real*8 :: x0(3,3)
     end function cell_niggli
     module function cell_delaunay(c,noenv,ti) result(x0)
       class(crystal), intent(inout) :: c
       logical, intent(in), optional :: noenv
       type(thread_info), intent(in), optional :: ti
       real*8 :: x0(3,3)
     end function cell_delaunay
     module subroutine reorder_atoms(c,iperm,ti)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: iperm(:)
       type(thread_info), intent(in), optional :: ti
     end subroutine reorder_atoms
     module subroutine wholemols(c,ti)
       class(crystal), intent(inout) :: c
       type(thread_info), intent(in), optional :: ti
     end subroutine wholemols
     module subroutine delete_atoms(c,nat,iat,ti)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: nat
       integer, intent(in) :: iat(nat)
       type(thread_info), intent(in), optional :: ti
     end subroutine delete_atoms
     module function sitesymm(c,x0,eps0,leqv,lrotm)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: x0(3)
       real*8, intent(in), optional :: eps0
       character*3 :: sitesymm
       integer, optional :: leqv
       real*8, optional :: lrotm(3,3,48)
     end function sitesymm
     module subroutine symeqv(c,xp0,mmult,vec,irotm,icenv,eps0)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: xp0(3)
       integer, intent(out) :: mmult
       real*8, allocatable, intent(inout), optional :: vec(:,:)
       integer, allocatable, intent(inout), optional :: irotm(:)
       integer, allocatable, intent(inout), optional :: icenv(:)
       real*8, intent(in), optional :: eps0
     end subroutine symeqv
     module function get_mult(c,x0) result (mult)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: x0(3)
       integer :: mult
     end function get_mult
     module subroutine spglib_wrap(c,spg,usenneq,errmsg,ti)
       class(crystal), intent(in) :: c
       type(SpglibDataset), intent(inout) :: spg
       logical, intent(in) :: usenneq
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine spglib_wrap
     module subroutine spgtowyc(c,spg)
       class(crystal), intent(inout) :: c
       type(SpglibDataset), intent(inout), optional :: spg
     end subroutine spgtowyc
     module subroutine calcsym(c,usenneq,errmsg,ti)
       class(crystal), intent(inout) :: c
       logical, intent(in) :: usenneq
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine calcsym
     module subroutine clearsym(c,cel2neq,neq2cel)
       class(crystal), intent(inout) :: c
       logical, intent(in), optional :: cel2neq
       logical, intent(in), optional :: neq2cel
     end subroutine clearsym
     module subroutine checkgroup(c)
       class(crystal), intent(inout) :: c
     end subroutine checkgroup
     module subroutine getiws(c,xorigin,ntetrag,tetrag)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: xorigin(3)
       integer, intent(out), optional :: ntetrag
       real*8, allocatable, intent(inout), optional :: tetrag(:,:,:)
     end subroutine getiws
     module subroutine powder(c,mode,th2ini0,th2end0,lambda0,fpol,npts,sigma,ishard,&
        th2p,ip,hvecp,discardp,t,ih)
       class(crystal), intent(in) :: c
       integer, intent(in) :: mode
       real*8, intent(in) :: th2ini0, th2end0
       real*8, intent(in) :: lambda0
       real*8, intent(in) :: fpol
       integer, intent(in), optional :: npts
       real*8, intent(in), optional :: sigma
       logical, intent(in), optional :: ishard
       real*8, allocatable, intent(inout), optional :: th2p(:)
       real*8, allocatable, intent(inout), optional :: ip(:)
       integer, allocatable, intent(inout), optional :: hvecp(:,:)
       real*8, intent(in), optional :: discardp
       real*8, allocatable, intent(inout), optional :: t(:)
       real*8, allocatable, intent(inout), optional :: ih(:)
     end subroutine powder
     module subroutine rdf(c,rini,rend,sigma,ishard,npts,t,ih,npairs0,ipairs0,ihat,intpeak)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: rini
       real*8, intent(in) :: rend
       real*8, intent(in) :: sigma
       logical, intent(in) :: ishard
       integer, intent(in) :: npts
       real*8, allocatable, intent(inout) :: t(:)
       real*8, allocatable, intent(inout) :: ih(:)
       integer, intent(in), optional :: npairs0
       integer, intent(in), optional :: ipairs0(:,:)
       real*8, allocatable, intent(inout), optional :: ihat(:,:)
       real*8, intent(in), optional :: intpeak(:)
     end subroutine rdf
     module subroutine amd(c,imax,res)
       class(crystal), intent(in) :: c
       integer, intent(in) :: imax
       real*8, intent(out) :: res(imax)
     end subroutine amd
     module subroutine struct_report(c,lcrys,lq)
       class(crystal), intent(in) :: c
       logical, intent(in) :: lcrys
       logical, intent(in) :: lq
     end subroutine struct_report
     module subroutine struct_report_symmetry(c)
       class(crystal), intent(in) :: c
     end subroutine struct_report_symmetry
     module subroutine struct_report_symxyz(c,strfin,doaxes)
       class(crystal), intent(in) :: c
       character(len=mlen), intent(out), optional :: strfin(c%neqv*c%ncv)
       logical, intent(in), optional :: doaxes
     end subroutine struct_report_symxyz
     module subroutine struct_write_json(c,json,p)
       use json_module, only: json_value, json_core
       class(crystal), intent(in) :: c
       type(json_core), intent(inout) :: json
       type(json_value), pointer, intent(inout) :: p
     end subroutine struct_write_json
     module subroutine write_simple_driver(c,file,ti)
       class(crystal), intent(inout) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_simple_driver
     module subroutine write_mol(c,file,fmt,ix0,doborder0,onemotif0,molmotif0,&
        environ0,renv0,lnmer0,nmer0,rsph0,xsph0,rcub0,xcub0,usenames0,luout,ti)
       class(crystal), intent(inout) :: c
       character*(*), intent(in) :: file
       character*3, intent(in) :: fmt
       integer, intent(in), optional :: ix0(3)
       logical, intent(in), optional :: doborder0, onemotif0, molmotif0, environ0
       real*8, intent(in), optional :: renv0
       logical, intent(in), optional :: lnmer0
       integer, intent(in), optional :: nmer0
       real*8, intent(in), optional :: rsph0, xsph0(3)
       real*8, intent(in), optional :: rcub0, xcub0(3)
       logical, intent(in), optional :: usenames0
       integer, intent(out), optional :: luout
       type(thread_info), intent(in), optional :: ti
     end subroutine write_mol
     module subroutine write_3dmodel(c,file,fmt,ix0,doborder0,onemotif0,molmotif0,&
        docell0,domolcell0,rsph0,xsph0,rcub0,xcub0,gr0,ti)
       use graphics, only: grhandle
       class(crystal), intent(inout) :: c
       character*(*), intent(in) :: file
       character*3, intent(in) :: fmt
       integer, intent(in), optional :: ix0(3)
       logical, intent(in), optional :: doborder0, onemotif0, molmotif0
       logical, intent(in), optional :: docell0, domolcell0
       real*8, intent(in), optional :: rsph0, xsph0(3)
       real*8, intent(in), optional :: rcub0, xcub0(3)
       type(grhandle), intent(out), optional :: gr0
       type(thread_info), intent(in), optional :: ti
     end subroutine write_3dmodel
     module subroutine write_espresso(c,file,rklength,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       real*8, intent(in), optional :: rklength
       type(thread_info), intent(in), optional :: ti
     end subroutine write_espresso
     module subroutine write_vasp(c,file,verbose,append,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       logical, intent(in) :: verbose
       logical, intent(in), optional :: append
       type(thread_info), intent(in), optional :: ti
     end subroutine write_vasp
     module subroutine write_abinit(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_abinit
     module subroutine write_elk(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_elk
     module subroutine write_gaussian(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_gaussian
     module subroutine write_tessel(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_tessel
     module subroutine write_critic(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_critic
     module subroutine write_cif(c,file,usesym0,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       logical, intent(in) :: usesym0
       type(thread_info), intent(in), optional :: ti
     end subroutine write_cif
     module subroutine write_d12(c,file,dosym,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       logical, intent(in) :: dosym
       type(thread_info), intent(in), optional :: ti
     end subroutine write_d12
     module subroutine write_res(c,file,dosym,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       integer, intent(in) :: dosym
       type(thread_info), intent(in), optional :: ti
     end subroutine write_res
     module subroutine write_escher(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_escher
     module subroutine write_db(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_db
     module subroutine write_gulp(c,file,ti)
       class(crystal), intent(inout) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_gulp
     module subroutine write_lammps(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_lammps
     module subroutine write_siesta_fdf(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_siesta_fdf
     module subroutine write_siesta_in(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_siesta_in
     module subroutine write_dftbp_hsd(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_dftbp_hsd
     module subroutine write_dftbp_gen(c,file,lu0,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       integer, intent(in), optional :: lu0
       type(thread_info), intent(in), optional :: ti
     end subroutine write_dftbp_gen
     module subroutine write_pyscf(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_pyscf
     module subroutine write_fhi(c,file,frac,rklength,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       logical, intent(in) :: frac
       real*8, intent(in), optional :: rklength
       type(thread_info), intent(in), optional :: ti
     end subroutine write_fhi
     module subroutine write_tinkerfrac(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_tinkerfrac
     module subroutine write_castep_cell(c,file,rklength,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       real*8, intent(in), optional :: rklength
       type(thread_info), intent(in), optional :: ti
     end subroutine write_castep_cell
     module subroutine writegrid_cube(c,g,file,onlyheader,binary,xd0,x00,ishift0,ti)
       class(crystal), intent(in) :: c
       real*8, intent(in), allocatable :: g(:,:,:)
       character*(*), intent(in) :: file
       logical, intent(in) :: onlyheader
       logical, intent(in) :: binary
       real*8, intent(in), optional :: xd0(3,3)
       real*8, intent(in), optional :: x00(3)
       integer, intent(in), optional :: ishift0(3)
       type(thread_info), intent(in), optional :: ti
     end subroutine writegrid_cube
     module subroutine writegrid_vasp(c,g,file,onlyheader,ishift0,ti)
       class(crystal), intent(in) :: c
       real*8, intent(in), allocatable :: g(:,:,:)
       character*(*), intent(in) :: file
       logical, intent(in) :: onlyheader
       integer, intent(in), optional :: ishift0(3)
       type(thread_info), intent(in), optional :: ti
     end subroutine writegrid_vasp
     module subroutine writegrid_xsf(c,g,file,onlyheader,ishift0,ti)
       class(crystal), intent(in) :: c
       real*8, intent(in), allocatable :: g(:,:,:)
       character*(*), intent(in) :: file
       logical, intent(in) :: onlyheader
       integer, intent(in), optional :: ishift0(3)
       type(thread_info), intent(in), optional :: ti
     end subroutine writegrid_xsf
     module subroutine search_lattice(x2r,rmax,imax,jmax,kmax)
       real*8, intent(in) :: x2r(3,3), rmax
       integer, intent(out) :: imax, jmax, kmax
     end subroutine search_lattice
     module subroutine pointgroup_info(hmpg,schpg,holo,laue)
       character*(*), intent(in) :: hmpg
       character(len=3), intent(out) :: schpg
       integer, intent(out) :: holo
       integer, intent(out) :: laue
     end subroutine pointgroup_info
     !xx! write submodule
  end interface

end module crystalmod
