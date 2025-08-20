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
  use spglib, only: SpglibDataset
  use types, only: neqatom, celatom, neighstar, species, cp_type
  use fragmentmod, only: fragment
  use param, only: maxzat0, mlen
  use types, only: thread_info
  implicit none

  private

  ! system periodicity
  integer, parameter, public :: iperiod_3d_crystal = 0   ! 3d periodic
  integer, parameter, public :: iperiod_3d_layered = 1   ! 3d periodic, layered
  integer, parameter, public :: iperiod_3d_chain = 2     ! 3d periodic, chain
  integer, parameter, public :: iperiod_3d_molecular = 3 ! 3d periodic, molecular
  integer, parameter, public :: iperiod_2d = 4           ! slab
  integer, parameter, public :: iperiod_1d = 5           ! chain
  integer, parameter, public :: iperiod_0d = 6           ! molecule-in-a-box
  integer, parameter, public :: iperiod_mol_single = 7   ! a single molecule
  integer, parameter, public :: iperiod_mol_cluster = 8  ! a molecular cluster
  real*8,  parameter, public :: iperiod_vacthr = 15d0    ! threshold for vacuum detection

  ! holohedry identifier
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

  ! Defaults for powder X-ray diffraction routines
  real*8, parameter, public :: xrpd_lambda_def = 1.5406d0
  real*8, parameter, public :: xrpd_fpol_def = 0d0
  real*8, parameter, public :: xrpd_alpha_def = 0.5d0
  real*8, parameter, public :: xrpd_sigma_def = 0.05d0
  real*8, parameter, public :: xrpd_th2ini_def = 5d0
  real*8, parameter, public :: xrpd_th2end_def = 50d0
  real*8, parameter, public :: xrpd_besteps_def_safe = 1d-4
  real*8, parameter, public :: xrpd_besteps_def_quick = 1d-3
  integer, parameter, public :: xrpd_maxfeval_def_safe = 10000
  integer, parameter, public :: xrpd_maxfeval_def_quick = 5000
  real*8, parameter, public :: xrpd_max_elong_def_safe = 0.15d0
  real*8, parameter, public :: xrpd_max_elong_def_quick = 0.1d0
  real*8, parameter, public :: xrpd_max_ang_def_safe = 10d0
  real*8, parameter, public :: xrpd_max_ang_def_quick = 5d0

  !> Class for molecular or crystal vibrations
  type vibrations
     !! source file and format
     character(len=mlen) :: file ! source of vibration data
     integer :: ivformat ! format for the vibration data
     !! 2nd-order force constants
     logical :: hasfc2 = .false. ! true if FC2 is available
     real*8, allocatable :: fc2(:,:,:,:) ! 2nd-order FC matrix (3,3,nat,nat)
     integer :: fc2_acoustic(3) = -1 ! acoustic branches for vs calculation
     real*8 :: fc2_vs_delta = -1d0 ! delta for calculation of sound velocities (bohr-1)
     real*8 :: fc2_vs_center = -1d0 ! center for calculation of sound velocities (bohr-1)
     !! frequencies and eigenvectors
     logical :: hasvibs = .false. ! true if frequencies/eigenvectors are available
     integer :: nqpt ! number of q-points
     real*8, allocatable :: qpt(:,:) ! q-point coordinates (3,nqpt) (fractional)
     integer :: nfreq ! number of frequencies
     real*8, allocatable :: freq(:,:) ! frequencies (nfreq,nqpt) (cm-1)
     complex*16, allocatable :: vec(:,:,:,:) ! phonon eigenvector (3,nat,nfreq,nqpt)
     ! note: these are the vectors that come out of the dynamical
     ! matrix diagonalization and they must be orthonormal. For the
     ! displacements, divide by the sqrt(m_j).
   contains
     procedure :: end => vibrations_end !< terminate the vibrations object
     procedure :: print_summary => vibrations_print_summary !< print summary of vibs
     procedure :: print_fc2 => vibrations_print_fc2 !< print info about the FC2
     procedure :: print_freq => vibrations_print_freq !< print frequency info
     procedure :: print_eigenvector => vibrations_print_eigenvector !< print eigvec info
     procedure :: read_file => vibrations_read_file !< read a vib file, detect the format
     procedure :: apply_acoustic => vibrations_apply_acoustic !< apply acoustic sum rules to FC2
     procedure :: write_fc2 => vibrations_write_fc2 !< write FC2
     procedure :: calculate_q => vibrations_calculate_q !< calculate freqs and vec for a single q
     procedure :: calculate_vs => vibrations_calculate_vs !< calculate freqs and vec for a single q
     procedure :: calculate_vs_prepare => vibrations_calculate_vs_prepare !< prepare vs calculation
     procedure :: calculate_thermo => vibrations_calculate_thermo !< calculate thermodynamic properties
     procedure :: trim_fc2 => vibrations_trim_fc2
     procedure :: zero_fc2 => vibrations_zero_fc2
  end type vibrations
  public :: vibrations

  !> The crystal class. A crystal contains the structural information for the
  !> system, and it can be an actual crystal (%ismolecule=.false.) or a molecule
  !> (%ismolecule=.true.) embedded in a large cell. When use in combination
  !> with the system class, the crystal needs to be passed as pointer to the
  !> system, and therefore needs to be declared as TARGET.
  type crystal
     ! Initialization flags
     logical :: isinit = .false. !< has the crystal structure been initialized?
     integer :: havesym = 0 !< was the symmetry determined? (0 - nosym, 1 - full)

     ! file name and format for the occasional critic2 trick
     character(len=mlen) :: file
     integer :: isformat

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
     real*8 :: m_rx2rc(3,3) !< input cell, reciprocal cryst. -> reciprocal cartesian
     real*8 :: m_rc2rx(3,3) !< input cell, reciprocal cartesian -> reciprocal cryst.
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

     ! atomic environment of the cell
     integer :: nblock(3) ! number of environemt blocks
     integer,allocatable :: iblock0(:,:,:) ! starting atomic index for each block
     real*8 :: blockrmax ! radius of the largest sphere contained in a block
     real*8 :: blockomega ! volume of a block
     real*8 :: blockcv(3) ! cross products of the block lattice vectors

     ! asterisms
     type(neighstar), allocatable :: nstar(:) !< Neighbor stars
     integer :: nmol = 0 !< Number of molecules in the unit cell
     type(fragment), allocatable :: mol(:) !< Molecular fragments
     integer :: nlvac = 0 !< Number of vacuum lattice vectors
     integer :: lvac(3,2) !< Vacuum lattice vectors
     integer :: lcon(3,2) !< Connected lattice vectors
     real*8 :: vaclength(3) !< Vacuum length in bohr along the crystallographic axes (a,b,c)
     real*8 :: vacbot(3) !< Bottom of the vacuum region in each axis
     real*8 :: vactop(3) !< Top of the vacuum region in each axis
     ! variables for 3d molecular crystals
     integer :: iperiod !< periodicity (see iperiod_* constants)
     logical :: ismol3d !< Is this a 3d molecular crystal?
     integer, allocatable :: idxmol(:) !< -1: mol is fractional, 0: sym. unique, >0 index for nneq mol.
     integer, allocatable :: idatcelmol(:,:) !< cell atom i belongs to idatcelmol(1,i) molecule; idatcelmol(2,i) is the atom ID in the molecule
     ! --> These two are equal:
     ! c%mol(idatcelmol(1,i))%at(idatcelmol(2,i))%x
     ! c%atcel(i)%x + c%mol(idatcelmol(1))%at(idatcelmol(2))%lvec

     ! vibrations
     type(vibrations) :: vib !< molecular/crystal vibrations
   contains
     ! construction, destruction, initialization (proc)
     procedure :: init => struct_init !< Allocate arrays and nullify variables
     procedure :: end => struct_end !< Deallocate arrays and nullify variables
     procedure :: struct_new !< Initialize the structure from a crystal seed

     ! basic crystallographic operations (proc)
     procedure :: x2c !< Convert input cryst. -> cartesian
     procedure :: c2x !< Convert input cartesian -> cryst.
     procedure :: rx2rc !< Convert input reciprocal cryst. -> reciprocal cartesian
     procedure :: rc2rx !< Convert input reciprocal cartesian -> reciprocal cryst.
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
     procedure :: nearest_atom_grid !< return the nearest atom IDs for a uniform grid of points
     procedure :: identify_spc !< Identify species in a structure from a string

     ! atomic environments and distance calculations (env)
     procedure :: build_env
     procedure :: list_near_atoms
     procedure :: nearest_atom
     procedure :: get_rnn2
     procedure :: identify_atom
     procedure :: promolecular_atom
     procedure :: find_asterisms
     procedure :: list_near_lattice_points
     procedure :: nearest_lattice_point
     procedure :: in_same_molecule

     ! molecular environments and neighbors (mols)
     procedure :: identify_fragment !< Build an atomic fragment of the crystal
     procedure :: identify_fragment_from_xyz !< Build a crystal fragment from an xyz file
     procedure :: fill_molecular_fragments !< Find the molecular fragments in the crystal
     procedure :: calculate_molecular_equivalence !< Calculate symmetry relations between molecules
     procedure :: calculate_periodicity !< Calculate symmetry relations between molecules
     procedure :: listatoms_cells !< List all atoms in n cells (maybe w border)
     procedure :: listatoms_sphcub !< List all atoms in a sphere or cube
     procedure :: listmolecules !< List all molecules in the crystal

     ! complex operations: ewald, promolecular, etc. (complex)
     procedure :: calculate_ewald_cutoffs !< Calculate the cutoffs for Ewald's sum
     procedure :: ewald_energy !< electrostatic energy (Ewald)
     procedure :: ewald_pot !< electrostatic potential (Ewald)
     procedure :: promolecular_array3 !< calculate core and promolecular densities on a grid
     procedure :: get_pack_ratio !< Calculate the packing ratio
     procedure :: vdw_volume !< Calculate the van der waals volume
     procedure :: get_kpoints !< k-point grid for a given rklength

     ! changes to the crystal/molecular structure (edit)
     procedure :: makeseed !< make a crystal seed from a crystal
     procedure :: makeseed_nudged !< make a crystal seed, displaced by a phonon
     procedure :: newcell !< Change the unit cell and rebuild the crystal
     procedure :: cell_standard !< Transform the the standard cell (possibly primitive)
     procedure :: cell_niggli !< Transform to the Niggli primitive cell
     procedure :: cell_delaunay !< Transform to the Delaunay primitive cell
     procedure :: reorder_atoms !< reorder the atoms in the crystal/molecule
     procedure :: wholemols !< Re-assign atomic types to have an asymmetric unit with whole molecules
     procedure :: delete_atoms !< Delete a list of atoms
     procedure :: move_atom !< Move an atom
     procedure :: move_cell !< Move the unit cell

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
     procedure :: write_any_file
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
     procedure :: write_pdb
     procedure :: write_castep_cell
     procedure :: write_alamode

     ! grid writers (write)
     procedure :: writegrid_cube
     procedure :: writegrid_vasp
     procedure :: writegrid_xsf
  end type crystal
  public :: crystal

  !> Class for XRPD peak information
  type xrpd_peaklist
     integer :: npeak = 0 ! number of peaks
     logical :: haveth2limits = .false. ! whether th2ini and th2end are available
     logical :: haveradiation = .false. ! whether lambda and fpol are available
     logical :: havehvec = .false. ! whether hvec data (hkl indices) are available
     logical :: havegradients = .false. ! whether the gradients of th2 and ip are available
     logical :: havepeakshape = .false. ! whether the fwhm and gau/lor coefs are available
     real*8 :: lambda ! wavelength of the radiation (angstrom)
     real*8 :: fpol ! polarization correction factor (0 = unpolarized, 0.95 = syncrhotron)
     real*8 :: th2ini, th2end ! initial and final reflection angles (degrees)
     real*8, allocatable :: th2(:) ! reflection angles (2*theta, degrees)
     real*8, allocatable :: ip(:) ! peak intensities
     integer, allocatable :: hvec(:,:) ! reflection indices
     real*8, allocatable :: th2g(:,:) ! gradient of th2 wrt metric tensor (degrees)
     real*8, allocatable :: ipg(:,:) ! gradient of ip wrt metric tensor
     real*8, allocatable :: fwhm(:) ! peak full width at half maximum (fwhm, degrees)
     real*8, allocatable :: cgau(:) ! Gaussian/Lorentzian peak shape coefficient
   contains
     procedure :: end => xrpd_peaklist_end
     procedure :: from_crystal => xrpd_peaks_from_crystal_powder
     procedure :: from_peaks_file => xrpd_peaks_from_peaks_file
     procedure :: from_profile_file => xrpd_peaks_from_profile_file
     procedure :: write => xrpd_write_to_file
     procedure :: calculate_profile => xrpd_calculate_profile
  end type xrpd_peaklist
  public :: xrpd_peaklist

  ! other crystallography tools that are crystal-independent (symmetry)
  public :: search_lattice
  public :: pointgroup_info
  public :: crosscorr_gaussian
  public :: vcpwdf_compare
  public :: gaussian_compare
  public :: david_sivia_calculate_background
  public :: struct_detect_write_format

  ! module procedure interfaces
  interface
     !xx! crystal type
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
     pure module function rx2rc(c,xx) result(res)
       class(crystal), intent(in) :: c
       real*8, intent(in) :: xx(3)
       real*8 :: res(3)
     end function rx2rc
     pure module function rc2rx(c,xx) result(res)
       class(crystal), intent(in) :: c
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function rc2rx
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
     module subroutine nearest_atom_grid(c,n,idg)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: n(3)
       integer, allocatable, intent(inout) :: idg(:,:,:)
     end subroutine nearest_atom_grid
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
     module subroutine nearest_atom(c,xp,icrd,nid,dist,distmax,lvec,nid0,id0,iz0,ispc0,nozero)
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
     end subroutine nearest_atom
     module function get_rnn2(c,ineq)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: ineq
       real*8 :: get_rnn2
     end function get_rnn2
     module function identify_atom(c,x0,icrd,lvec,dist,distmax)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: x0(3)
       integer, intent(in) :: icrd
       integer, intent(out), optional :: lvec(3)
       real*8, intent(out), optional :: dist
       real*8, intent(in), optional :: distmax
       integer :: identify_atom
     end function identify_atom
     module subroutine promolecular_atom(c,x0,icrd,f,fp,fpp,nder,zpsp,fr)
       use fragmentmod, only: fragment
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: x0(3)
       integer, intent(in) :: icrd
       real*8, intent(out) :: f
       real*8, intent(out) :: fp(3)
       real*8, intent(out) :: fpp(3,3)
       integer, intent(in) :: nder
       integer, intent(in), optional :: zpsp(:)
       type(fragment), intent(in), optional :: fr
     end subroutine promolecular_atom
     module subroutine find_asterisms(c,nstar,atmrad,bondfac,rij)
       use param, only: maxzat0
       class(crystal), intent(inout) :: c
       type(neighstar), allocatable, intent(inout) :: nstar(:)
       real*8, intent(in), optional :: atmrad(0:maxzat0)
       real*8, intent(in), optional :: bondfac
       real*8, intent(in), optional :: rij(:,:,:)
     end subroutine find_asterisms
     module subroutine list_near_lattice_points(c,xp,icrd,sorted,nat,dist,lvec,ndiv,&
        up2d,up2n,nozero)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: xp(3)
       integer, intent(in) :: icrd
       logical, intent(in) :: sorted
       integer, intent(out) :: nat
       real*8, allocatable, intent(inout), optional :: dist(:)
       integer, allocatable, intent(inout), optional :: lvec(:,:)
       integer, intent(in), optional :: ndiv(3)
       real*8, intent(in), optional :: up2d
       integer, intent(in), optional :: up2n
       logical, intent(in), optional :: nozero
     end subroutine list_near_lattice_points
     module subroutine nearest_lattice_point(c,xp,icrd,dist,lvec,ndiv,nozero)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: xp(3)
       integer, intent(in) :: icrd
       real*8, intent(out) :: dist
       integer, intent(out), optional :: lvec(3)
       integer, intent(in), optional :: ndiv(3)
       logical, intent(in), optional :: nozero
     end subroutine nearest_lattice_point
     module function in_same_molecule(c,i,il,j,jl)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: i, j
       integer, intent(in) :: il(3), jl(3)
       logical :: in_same_molecule
     end function in_same_molecule
     module function identify_fragment(c,nat,x0) result(fr)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: nat
       real*8, intent(in) :: x0(3,nat)
       type(fragment) :: fr
     end function identify_fragment
     module function identify_fragment_from_xyz(c,file,errmsg,ti) result(fr)
       class(crystal), intent(inout) :: c
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
     module subroutine calculate_periodicity(c)
       class(crystal), intent(inout) :: c
     end subroutine calculate_periodicity
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
     module subroutine calculate_ewald_cutoffs(c,rcut,hcut,eta,qsum,lrmax,lhmax)
       class(crystal), intent(inout) :: c
       real*8, intent(out) :: rcut, hcut, eta, qsum
       integer, intent(out) :: lrmax(3), lhmax(3)
     end subroutine calculate_ewald_cutoffs
     module function ewald_energy(c) result(ewe)
       class(crystal), intent(inout) :: c
       real*8 :: ewe
     end function ewald_energy
     module function ewald_pot(c,x,rcut,hcut,eta,qsum,lrmax,lhmax)
       class(crystal), intent(inout) :: c
       real*8, intent(in) :: x(3)
       real*8, intent(out) :: rcut, hcut, eta, qsum
       integer, intent(out) :: lrmax(3), lhmax(3)
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
     module subroutine promolecular_array3(c,f,n,zpsp,fr)
       class(crystal), intent(inout) :: c
       real*8, intent(inout), allocatable :: f(:,:,:)
       integer, intent(in) :: n(3)
       integer, intent(in), optional :: zpsp(:)
       type(fragment), intent(in), optional :: fr
     end subroutine promolecular_array3
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
     module subroutine makeseed(c,seed,copysym,useabr)
       use crystalseedmod, only: crystalseed
       class(crystal), intent(in) :: c
       type(crystalseed), intent(out) :: seed
       logical, intent(in) :: copysym
       integer, intent(in), optional :: useabr
     end subroutine makeseed
     module subroutine makeseed_nudged(c,seed,qpt,evec,amplitude,phase)
       use crystalseedmod, only: crystalseed
       class(crystal), intent(in) :: c
       type(crystalseed), intent(out) :: seed
       real*8, intent(in) :: qpt(3)
       complex*16, intent(in) :: evec(:,:)
       real*8, intent(in) :: amplitude, phase
     end subroutine makeseed_nudged
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
     module subroutine move_atom(c,idx,x,iunit_l,dorelative,ti)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: idx
       real*8, intent(in) :: x(3)
       integer, intent(in) :: iunit_l
       logical, intent(in) :: dorelative
       type(thread_info), intent(in), optional :: ti
     end subroutine move_atom
     module subroutine move_cell(c,iaxis,x,iunit_l,dorelative,dofraction,ti)
       class(crystal), intent(inout) :: c
       integer, intent(in) :: iaxis
       real*8, intent(in) :: x
       integer, intent(in) :: iunit_l
       logical, intent(in) :: dorelative, dofraction
       type(thread_info), intent(in), optional :: ti
     end subroutine move_cell
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
       class(crystal), intent(inout) :: c
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
       class(crystal), intent(inout) :: c
       integer, intent(in) :: imax
       real*8, intent(out) :: res(imax)
     end subroutine amd
     module subroutine struct_report(c,lcrys,lq)
       class(crystal), intent(inout) :: c
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
       class(crystal), intent(inout) :: c
       type(json_core), intent(inout) :: json
       type(json_value), pointer, intent(inout) :: p
     end subroutine struct_write_json
     module subroutine write_any_file(c,file,errmsg,iwformat,ti)
       class(crystal), intent(inout) :: c
       character*(*), intent(in) :: file
       character(len=:), allocatable, intent(inout) :: errmsg
       integer, intent(in), optional :: iwformat
       type(thread_info), intent(in), optional :: ti
     end subroutine write_any_file
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
     module subroutine write_d12(c,file,dosym,doexternal,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       logical, intent(in) :: dosym
       logical, intent(in) :: doexternal
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
     module subroutine write_pdb(c,file,cp,cpcel,ixzassign,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(cp_type), intent(in), optional :: cp(:)
       type(cp_type), intent(in), optional :: cpcel(:)
       integer, intent(in), optional :: ixzassign(:)
       type(thread_info), intent(in), optional :: ti
     end subroutine write_pdb
     module subroutine write_castep_cell(c,file,rklength,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       real*8, intent(in), optional :: rklength
       type(thread_info), intent(in), optional :: ti
     end subroutine write_castep_cell
     module subroutine write_alamode(c,file,ti)
       class(crystal), intent(in) :: c
       character*(*), intent(in) :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine write_alamode
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
     !xx! vibrations type
     module subroutine vibrations_end(v,keepfc2)
       class(vibrations), intent(inout) :: v
       logical, intent(in), optional :: keepfc2
     end subroutine vibrations_end
     module subroutine vibrations_read_file(v,c,file,sline,ivformat,errmsg,ti)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       character*(*), intent(in) :: file, sline
       integer, intent(in) :: ivformat
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine vibrations_read_file
     module subroutine vibrations_print_summary(v)
       class(vibrations), intent(inout) :: v
     end subroutine vibrations_print_summary
     module subroutine vibrations_print_fc2(v,c,disteps,fc2eps,environ)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       real*8, intent(in), optional :: disteps, fc2eps
       logical, intent(in), optional :: environ
     end subroutine vibrations_print_fc2
     module subroutine vibrations_print_freq(v,id)
       class(vibrations), intent(inout) :: v
       integer, intent(in) :: id
     end subroutine vibrations_print_freq
     module subroutine vibrations_print_eigenvector(v,c,ifreq,idq,cartesian)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(in) :: c
       integer, intent(in) :: ifreq, idq
       logical, intent(in) :: cartesian
     end subroutine vibrations_print_eigenvector
     module subroutine vibrations_apply_acoustic(v,c,verbose)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       logical, intent(in), optional :: verbose
     end subroutine vibrations_apply_acoustic
     module subroutine vibrations_write_fc2(v,c,file,verbose)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       character(len=:), allocatable, intent(in), optional :: file
       logical, intent(in), optional :: verbose
     end subroutine vibrations_write_fc2
     module subroutine vibrations_calculate_q(v,c,q,freqo,veco)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       real*8, intent(in) :: q(3)
       real*8, intent(inout), allocatable, optional :: freqo(:)
       complex*16, intent(inout), allocatable, optional :: veco(:,:)
     end subroutine vibrations_calculate_q
     module subroutine vibrations_calculate_vs(v,c,q,vs)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       real*8, intent(in) :: q(3)
       real*8, intent(out) :: vs(3)
     end subroutine vibrations_calculate_vs
     module subroutine vibrations_calculate_vs_prepare(v,c,q,vs,verbose)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       real*8, intent(in) :: q(3)
       real*8, intent(out) :: vs(3)
       logical, intent(in) :: verbose
     end subroutine vibrations_calculate_vs_prepare
     module subroutine vibrations_calculate_thermo(v,t,zpe,fvib,svib,cv)
       class(vibrations), intent(inout) :: v
       real*8, intent(in) :: t
       real*8, intent(out) :: zpe, fvib, svib, cv
     end subroutine vibrations_calculate_thermo
     module subroutine vibrations_trim_fc2(v,c,dist,verbose)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       real*8, intent(in) :: dist
       logical, intent(in), optional :: verbose
     end subroutine vibrations_trim_fc2
     module subroutine vibrations_zero_fc2(v,c,eps0,verbose)
       class(vibrations), intent(inout) :: v
       type(crystal), intent(inout) :: c
       real*8, intent(in) :: eps0
       logical, intent(in), optional :: verbose
     end subroutine vibrations_zero_fc2
     !xx! xrpd_peaklist type
     module subroutine xrpd_peaklist_end(p)
       class(xrpd_peaklist), intent(inout) :: p
     end subroutine xrpd_peaklist_end
     module subroutine xrpd_peaks_from_crystal_powder(p,c,th2ini0,th2end0,lambda0,fpol,usehvecp,calcderivs,&
        errmsg,gg)
       class(xrpd_peaklist), intent(inout) :: p
       type(crystal), intent(in) :: c
       real*8, intent(in) :: th2ini0, th2end0
       real*8, intent(in) :: lambda0
       real*8, intent(in) :: fpol
       logical, intent(in) :: usehvecp
       logical, intent(in) :: calcderivs
       character(len=:), allocatable, intent(out) :: errmsg
       real*8, intent(in), optional :: gg(3,3)
     end subroutine xrpd_peaks_from_crystal_powder
     module subroutine xrpd_peaks_from_peaks_file(p,file,errmsg)
       class(xrpd_peaklist), intent(inout) :: p
       character*(*), intent(in) :: file
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine xrpd_peaks_from_peaks_file
     module subroutine xrpd_peaks_from_profile_file(p,xyfile,rms,errmsg,verbose0,&
        ymax_detect0,nadj0,pkinput,xorig,yorig,ycalc)
       class(xrpd_peaklist), intent(inout) :: p
       character*(*), intent(in) :: xyfile
       real*8, intent(out) :: rms
       character(len=:), allocatable, intent(out) :: errmsg
       logical, intent(in), optional :: verbose0
       real*8, intent(in), optional :: ymax_detect0
       integer, intent(in), optional :: nadj0
       type(xrpd_peaklist), intent(in), optional :: pkinput
       real*8, intent(inout), allocatable, optional :: xorig(:)
       real*8, intent(inout), allocatable, optional :: yorig(:)
       real*8, intent(inout), allocatable, optional :: ycalc(:)
     end subroutine xrpd_peaks_from_profile_file
     module subroutine xrpd_write_to_file(p,file)
       class(xrpd_peaklist), intent(in) :: p
       character*(*), intent(in) :: file
     end subroutine xrpd_write_to_file
     module subroutine xrpd_calculate_profile(p,n,x,y,errmsg,th2ini,th2end)
       class(xrpd_peaklist), intent(inout) :: p
       integer, intent(in) :: n
       real*8, allocatable, intent(inout) :: x(:), y(:)
       character(len=:), allocatable, intent(out) :: errmsg
       real*8, intent(in), optional :: th2ini, th2end
     end subroutine xrpd_calculate_profile
     !xx! independent procedures
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
     module subroutine crosscorr_gaussian(p1,p2,alpha,sigma,d12,errmsg,calcderivs,d12g)
       type(xrpd_peaklist), intent(in) :: p1, p2
       real*8, intent(in) :: alpha, sigma
       real*8, intent(out) :: d12
       character(len=:), allocatable, intent(out) :: errmsg
       logical, intent(in), optional :: calcderivs
       real*8, intent(out), optional :: d12g(6)
     end subroutine crosscorr_gaussian
     module subroutine vcpwdf_compare(c1,c2,diff,errmsg,max_elong,max_ang,max_vol,&
        powdiff_thr,c2out,verbose)
       use crystalseedmod, only: crystalseed
       type(crystal), intent(in) :: c1, c2
       real*8, intent(out) :: diff
       character(len=:), allocatable, intent(out) :: errmsg
       real*8, intent(in), optional :: max_elong
       real*8, intent(in), optional :: max_ang
       real*8, intent(in), optional :: max_vol
       real*8, intent(in), optional :: powdiff_thr
       type(crystal), intent(out), optional :: c2out
       logical, intent(in), optional :: verbose
     end subroutine vcpwdf_compare
     module subroutine gaussian_compare(c1,p2,imode,diff,errmsg,seedout,verbose0,alpha0,&
        lambda0,fpol0,maxfeval0,besteps0,max_elong_def0,max_ang_def0)
       use crystalseedmod, only: crystalseed
       type(crystal), intent(in) :: c1
       type(xrpd_peaklist), intent(in) :: p2
       integer, intent(in) :: imode
       real*8, intent(out) :: diff
       character(len=:), allocatable, intent(out) :: errmsg
       type(crystalseed), intent(out), optional :: seedout
       logical, intent(in), optional :: verbose0
       real*8, intent(in), optional :: alpha0
       real*8, intent(in), optional :: lambda0
       real*8, intent(in), optional :: fpol0
       integer, intent(in), optional :: maxfeval0
       real*8, intent(in), optional :: besteps0
       real*8, intent(in), optional :: max_elong_def0
       real*8, intent(in), optional :: max_ang_def0
     end subroutine gaussian_compare
     module function david_sivia_calculate_background(n,x,y,errmsg,nknot) result(yb)
       integer, intent(in) :: n
       real*8, intent(in) :: x(n)
       real*8, intent(in) :: y(n)
       character(len=:), allocatable, intent(out) :: errmsg
       integer, intent(in), optional :: nknot
       real*8 :: yb(n)
     end function david_sivia_calculate_background
     module subroutine struct_detect_write_format(file,isformat)
       character*(*), intent(in) :: file
       integer, intent(out) :: isformat
     end subroutine struct_detect_write_format
  end interface

end module crystalmod
