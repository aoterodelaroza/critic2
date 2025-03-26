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

! Crystal seed class, contains the structure file readers.
module crystalseedmod
  use fragmentmod, only: fragment
  use types, only: species, thread_info
  use param, only: mlen
  implicit none

  private

  !> Minimal amount of information to generate a crystal
  type crystalseed
     ! general
     logical :: isused = .false. !< Is the seed being used?
     character(len=mlen) :: file = "" !< Source file, if available
     character(len=mlen) :: name = "" !< Source file, if available
     integer :: isformat !< source file format
     ! atoms
     integer :: nat = 0 !< Number of atoms
     real*8, allocatable :: x(:,:) !< Atomic positions (crystal - fractional;molecule with useabr=0 - bohr)
     integer, allocatable :: is(:) !< Species of a given atom
     character*10, allocatable :: atname(:) !< atomic names
     ! species
     integer :: nspc = 0 !< Number of species
     type(species), allocatable :: spc(:) !< Species
     ! cell
     integer :: useabr = 0 !< 0 = uninit; 1 = use aa,bb; 2 = use m_x2c
     real*8 :: aa(3) !< Cell lengths (bohr)
     real*8 :: bb(3) !< Cell angles (degrees)
     real*8 :: m_x2c(3,3) !< Crystallographic to cartesian matrix
     ! symmetry
     integer :: havesym = 0 !< Have symmetry? 0 = no; 1 = yes
     integer :: findsym = -1 !< Find the symmetry? 0 = no; 1 = yes; -1 = only small
     logical :: checkrepeats = .false. !< Check if atoms are repeated on crystal build
     integer :: neqv = 0 !< Number of symmetry operations
     integer :: ncv = 0 !< Number ofcentering vectors
     real*8, allocatable :: cen(:,:) !< Centering vectors
     real*8, allocatable :: rotm(:,:,:) !< Space group operations
     ! molecular fields
     logical :: ismolecule = .false. !< Is this a molecule?
     logical :: cubic = .false. !< Use a cubic cell for the molecule
     real*8 :: border = 0d0 !< border of the molecular cell (bohr)
     logical :: havex0 = .false. !< an origin of the cell for molecules has bene given
     real*8 :: molx0(3) = 0d0 !< origin of the cell for molecules
     ! calculated properties
     real*8 :: energy = huge(1d0) !< energy (Hartree)
     real*8 :: pressure = huge(1d0) !< pressure (GPa)
   contains
     procedure :: end => seed_end !< Deallocate arrays and nullify variables
     procedure :: parse_crystal_env
     procedure :: parse_molecule_env
     procedure :: from_fragment
     procedure :: strip_hydrogens
     procedure :: read_library
     procedure :: read_any_file
     procedure :: read_cif
     procedure :: read_mol2
     procedure :: read_sdf
     procedure :: read_shelx
     procedure :: read_magres
     procedure :: read_f21
     procedure :: read_cube
     procedure :: read_bincube
     procedure :: read_wien
     procedure :: read_vasp
     procedure :: read_potcar
     procedure :: read_abinit
     procedure :: read_elk
     procedure :: read_mol
     procedure :: read_pdb
     procedure :: read_qeout
     procedure :: read_qein
     procedure :: read_crystalout
     procedure :: read_fploout
     procedure :: read_siesta
     procedure :: read_castep_cell
     procedure :: read_castep_geom
     procedure :: read_dmain
     procedure :: read_dftbp
     procedure :: read_xsf
     procedure :: read_pwc
     procedure :: read_axsf
     procedure :: read_aimsin
     procedure :: read_aimsout
     procedure :: read_tinkerfrac
     procedure :: read_xyz
     procedure :: report
  end type crystalseed
  public :: crystalseed

  ! assignment operator for the crystalseed class
  interface assignment(=)
     module subroutine assign_crystalseed(to,from)
       class(crystalseed), intent(out) :: to
       type(crystalseed), intent(in) :: from
     end subroutine assign_crystalseed
  end interface

  public :: realloc_crystalseed
  public :: struct_detect_format
  public :: struct_detect_ismol
  public :: read_seeds_from_file
  public :: read_alat_from_qeout

  ! module procedure interfaces
  interface
     module subroutine seed_end(seed)
       class(crystalseed), intent(inout) :: seed !< Crystal seed output
     end subroutine seed_end
     module subroutine parse_crystal_env(seed,lu,oksyn)
       class(crystalseed), intent(inout) :: seed
       integer, intent(in) :: lu
       logical, intent(out) :: oksyn
     end subroutine parse_crystal_env
     module subroutine parse_molecule_env(seed,lu,oksyn)
       class(crystalseed), intent(inout) :: seed
       integer, intent(in) :: lu
       logical, intent(out) :: oksyn
     end subroutine parse_molecule_env
     module subroutine read_library(seed,line,oksyn,mol,file,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: line
       logical, intent(out) :: oksyn
       logical, intent(in), optional :: mol
       character*(*), intent(in), optional :: file
       type(thread_info), intent(in), optional :: ti
     end subroutine read_library
     module subroutine from_fragment(seed,fr,order_by_cidx0)
       class(crystalseed), intent(inout) :: seed
       type(fragment), intent(in) :: fr
       logical, intent(in), optional :: order_by_cidx0
     end subroutine from_fragment
     module subroutine strip_hydrogens(seed)
       class(crystalseed), intent(inout) :: seed
     end subroutine strip_hydrogens
     module subroutine read_any_file(seed,file,mol0,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       integer, intent(in) :: mol0
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_any_file
     module subroutine read_cif(seed,file,dblock,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       character*(*), intent(in) :: dblock
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_cif
     module subroutine read_mol2(seed,file,rborder,docube,name,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character*(*), intent(in) :: name
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_mol2
     module subroutine read_sdf(seed,file,rborder,docube,errmsg,id,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character(len=:), allocatable, intent(out) :: errmsg
       integer, intent(in), optional :: id
       type(thread_info), intent(in), optional :: ti
     end subroutine read_sdf
     module subroutine read_shelx(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout)  :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_shelx
     module subroutine read_magres(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout)  :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_magres
     module subroutine read_f21(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_f21
     module subroutine read_cube(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_cube
     module subroutine read_bincube(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_bincube
     module subroutine read_wien(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_wien
     module subroutine read_vasp(seed,file,mol,hastypes,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       logical, intent(out) :: hastypes
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_vasp
     module subroutine read_abinit(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_abinit
     module subroutine read_elk(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_elk
     module subroutine read_mol(seed,file,fmt,rborder,docube,errmsg,alsovib,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       integer, intent(in) :: fmt
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character(len=:), allocatable, intent(out) :: errmsg
       logical, intent(out), optional :: alsovib
       type(thread_info), intent(in), optional :: ti
     end subroutine read_mol
     module subroutine read_pdb(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_pdb
     module subroutine read_qeout(seed,file,mol,istruct,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       integer, intent(in) :: istruct
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_qeout
     module subroutine read_qein(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_qein
     module subroutine read_crystalout(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_crystalout
     module subroutine read_fploout(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_fploout
     module subroutine read_siesta(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_siesta
     module subroutine read_castep_cell(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_castep_cell
     module subroutine read_castep_geom(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_castep_geom
     module subroutine read_dmain(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_dmain
     module subroutine read_dftbp(seed,file,rborder,docube,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_dftbp
     module subroutine read_xsf(seed,file,rborder,docube,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_xsf
     module subroutine read_pwc(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_pwc
     module subroutine read_axsf(seed,file,nread0,xnudge,rborder,docube,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       integer, intent(in) :: nread0
       real*8, intent(in) :: xnudge
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_axsf
     module subroutine read_aimsin(seed,file,mol,rborder,docube,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_aimsin
     module subroutine read_aimsout(seed,file,mol,rborder,docube,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_aimsout
     module subroutine read_tinkerfrac(seed,file,mol,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_tinkerfrac
     module subroutine read_xyz(seed,file,mol,rborder,docube,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_xyz
     module subroutine realloc_crystalseed(a,nnew)
       type(crystalseed), intent(inout), allocatable :: a(:)
       integer, intent(in) :: nnew
     end subroutine realloc_crystalseed
     module subroutine struct_detect_format(file,isformat,alsofield,ti)
       character*(*), intent(in) :: file
       integer, intent(out) :: isformat
       logical, intent(out), optional :: alsofield
       type(thread_info), intent(in), optional :: ti
     end subroutine struct_detect_format
     module subroutine struct_detect_ismol(file,isformat,ismol,ti)
       character*(*), intent(in) :: file
       integer, intent(in) :: isformat
       logical, intent(out) :: ismol
       type(thread_info), intent(in), optional :: ti
     end subroutine struct_detect_ismol
     module subroutine read_potcar(seed,file,errmsg,ti)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_potcar
     module subroutine read_seeds_from_file(file,mol0,isformat0,readlastonly,&
        nseed,seed,collapse,errmsg,iafield,iavib,ti)
       character*(*), intent(in) :: file
       integer, intent(in) :: mol0
       integer, intent(in) :: isformat0
       logical, intent(in) :: readlastonly
       integer, intent(out) :: nseed
       type(crystalseed), allocatable, intent(inout) :: seed(:)
       logical, intent(out) :: collapse
       character(len=:), allocatable, intent(out) :: errmsg
       integer, intent(out), optional :: iafield
       integer, intent(out), optional :: iavib
       type(thread_info), intent(in), optional :: ti
     end subroutine read_seeds_from_file
     module subroutine report(seed)
       class(crystalseed), intent(inout) :: seed
     end subroutine report
     module subroutine read_alat_from_qeout(file,alat,errmsg,ti)
       character*(*), intent(in) :: file
       real*8, intent(out) :: alat
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine read_alat_from_qeout
  end interface

end module crystalseedmod
