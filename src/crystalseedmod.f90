! Copyright (c) 2015 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>,
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

! Crystal seed class. External file readers.
module crystalseedmod
  use types, only: species
  implicit none

  private

  !> Minimal amount of information to generate a crystal
  type crystalseed
     ! general
     logical :: isused = .false. !< Is the seed being used?
     character(len=512) :: file = "" !< Source file, if available
     character(len=512) :: name = "" !< Source file, if available
     ! atoms
     integer :: nat = 0 !< Number of atoms
     real*8, allocatable :: x(:,:) !< Atomic positions (crystal - fractional;molecule with useabr=0 - bohr)
     integer, allocatable :: is(:) !< Species of a given atom
     ! species 
     integer :: nspc = 0 !< Number of species
     type(species), allocatable :: spc(:) !< Species
     ! cell
     integer :: useabr = 0 !< 0 = uninit; 1 = use aa,bb; 2 = use crys2car
     real*8 :: aa(3) !< Cell lengths (bohr)
     real*8 :: bb(3) !< Cell angles (degrees)
     real*8 :: crys2car(3,3) !< Crystallographic to cartesian matrix
     ! symmetry
     integer :: havesym = 0 !< Have symmetry? 0 = no; 1 = yes
     integer :: findsym = -1 !< Find the symmetry? 0 = no; 1 = yes; -1 = only small
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
   contains
     procedure :: parse_crystal_env 
     procedure :: parse_molecule_env 
     procedure :: read_library
     procedure :: read_cif 
     procedure :: read_shelx 
     procedure :: read_cube 
     procedure :: read_wien 
     procedure :: read_vasp 
     procedure :: read_potcar
     procedure :: read_abinit 
     procedure :: read_elk 
     procedure :: read_mol 
     procedure :: read_qeout 
     procedure :: read_qein 
     procedure :: read_crystalout 
     procedure :: read_siesta 
     procedure :: read_dftbp 
     procedure :: read_xsf 
     procedure :: spgs_wrap
  end type crystalseed
  public :: crystalseed

  public :: realloc_crystalseed
  public :: struct_detect_format
  public :: read_seeds_from_file
  public :: read_all_cif 
  public :: read_all_qeout
  private :: read_cif_items
  private :: is_espresso
  private :: qe_latgen

  ! module procedure interfaces
  interface
     module subroutine parse_crystal_env(seed,lu,oksyn)
       class(crystalseed), intent(inout) :: seed
       integer, intent(in) :: lu
       logical, intent(out) :: oksyn
     end subroutine
     module subroutine parse_molecule_env(seed,lu,oksyn)
       class(crystalseed), intent(inout) :: seed
       integer, intent(in) :: lu
       logical, intent(out) :: oksyn
     end subroutine parse_molecule_env
     module subroutine read_library(seed,line,mol,oksyn)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: line
       logical, intent(in) :: mol
       logical, intent(out) :: oksyn
     end subroutine read_library
     module subroutine read_cif(seed,file,dblock,mol,errmsg)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       character*(*), intent(in) :: dblock
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine read_cif
     module subroutine read_shelx(seed,file,mol,errmsg)
       class(crystalseed) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine read_shelx
     module subroutine read_cube(seed,file,mol)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_cube
     module subroutine read_wien(seed,file,mol)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_wien
     module subroutine read_vasp(seed,file,mol,hastypes)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       logical, intent(out), optional :: hastypes
     end subroutine read_vasp
     module subroutine read_abinit(seed,file,mol)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_abinit
     module subroutine read_elk(seed,file,mol)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_elk
     module subroutine read_mol(seed,file,fmt,rborder,docube)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       integer, intent(in) :: fmt
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
     end subroutine read_mol
     module subroutine read_qeout(seed,file,mol,istruct)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       integer, intent(in) :: istruct
     end subroutine read_qeout
     module subroutine read_qein(seed,file,mol)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_qein
     module subroutine read_crystalout(seed,file,mol)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_crystalout
     module subroutine read_siesta(seed,file,mol)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_siesta
     module subroutine read_dftbp(seed,file,molout,rborder,docube)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(out) :: molout
       real*8, intent(in) :: rborder
       logical, intent(in) :: docube
     end subroutine read_dftbp
     module subroutine read_xsf(seed,file,mol)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_xsf
     module subroutine realloc_crystalseed(a,nnew)
       type(crystalseed), intent(inout), allocatable :: a(:)
       integer, intent(in) :: nnew
     end subroutine realloc_crystalseed
     module subroutine struct_detect_format(file,isformat,ismol,alsofield)
       character*(*), intent(in) :: file
       integer, intent(out) :: isformat
       logical, intent(out) :: ismol
       logical, intent(out), optional :: alsofield
     end subroutine struct_detect_format
     module subroutine read_potcar(seed,file)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: file
     end subroutine read_potcar
     module subroutine read_seeds_from_file(file,mol0,nseed,seed,errmsg,iafield)
       character*(*), intent(in) :: file
       integer, intent(in) :: mol0
       integer, intent(out) :: nseed
       type(crystalseed), allocatable, intent(inout) :: seed(:)
       character(len=:), allocatable, intent(out) :: errmsg
       integer, intent(out), optional :: iafield
     end subroutine read_seeds_from_file
     module subroutine read_all_cif(nseed,seed,file,mol,errmsg)
       integer, intent(out) :: nseed
       type(crystalseed), intent(inout), allocatable :: seed(:)
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine read_all_cif
     module subroutine read_all_qeout(nseed,seed,file,mol)
       integer, intent(out) :: nseed
       type(crystalseed), intent(inout), allocatable :: seed(:)
       character*(*), intent(in) :: file
       logical, intent(in) :: mol
     end subroutine read_all_qeout
     module subroutine read_cif_items(seed,mol,errmsg)
       type(crystalseed), intent(inout) :: seed
       logical, intent(in) :: mol
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine read_cif_items
     module function is_espresso(file)
       logical :: is_espresso
       character*(*), intent(in) :: file
     end function is_espresso
     module subroutine qe_latgen(ibrav,celldm,a1,a2,a3)
       integer, parameter :: dp = selected_real_kind(14,200)
       integer, intent(in) :: ibrav
       real(DP), intent(inout) :: celldm(6)
       real(DP), intent(out) :: a1(3), a2(3), a3(3)
     end subroutine qe_latgen
     module subroutine spgs_wrap(seed,spg,usespgr)
       class(crystalseed), intent(inout) :: seed
       character*(*), intent(in) :: spg
       logical, intent(in) :: usespgr
     end subroutine spgs_wrap
  end interface

end module crystalseedmod
