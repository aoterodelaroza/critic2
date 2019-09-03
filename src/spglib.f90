! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

module spglib
  use iso_c_binding, only: c_char, c_int, c_double, c_ptr, c_null_char, c_f_pointer, c_associated

  implicit none

  private

  public :: SpglibSpaceGroupType
  public :: SpglibDataset
  public :: spg_get_major_version
  public :: spg_get_minor_version
  public :: spg_get_micro_version
  public :: spg_standardize_cell, spgat_standardize_cell
  public :: spg_get_symmetry, spgat_get_symmetry
  public :: spg_get_symmetry_with_collinear_spin, spgat_get_symmetry_with_collinear_spin
  public :: spg_get_multiplicity, spgat_get_multiplicity
  public :: spg_get_hall_number_from_symmetry
  public :: spg_delaunay_reduce
  public :: spg_niggli_reduce
  public :: spg_find_primitive, spgat_find_primitive
  public :: spg_get_international, spgat_get_international
  public :: spg_get_schoenflies, spgat_get_schoenflies
  public :: spg_get_pointgroup
  public :: spg_refine_cell, spgat_refine_cell
  public :: spg_get_error_code

  public :: spg_get_error_message
  public :: spg_get_dataset
  public :: spg_get_spacegroup_type
  public :: spg_get_hall_number_from_symbol
  public :: spg_get_symmetry_from_database
  public :: spg_list_spg

  enum, bind(C)
     enumerator ::  SPGLIB_SUCCESS = 0
     enumerator ::  SPGERR_SPACEGROUP_SEARCH_FAILED
     enumerator ::  SPGERR_CELL_STANDARDIZATION_FAILED
     enumerator ::  SPGERR_SYMMETRY_OPERATION_SEARCH_FAILED
     enumerator ::  SPGERR_ATOMS_TOO_CLOSE
     enumerator ::  SPGERR_POINTGROUP_NOT_FOUND
     enumerator ::  SPGERR_NIGGLI_FAILED
     enumerator ::  SPGERR_DELAUNAY_FAILED
     enumerator ::  SPGERR_ARRAY_SIZE_SHORTAGE
     enumerator ::  SPGERR_NONE
  end enum

  type :: SpglibSpaceGroupType
     integer(c_int) :: number
     character(len=11) :: international_short
     character(len=20) :: international_full
     character(len=32) :: international
     character(len=7) :: schoenflies
     character(len=17) :: hall_symbol
     character(len=6) :: choice
     character(len=6) :: pointgroup_international
     character(len=4) :: pointgroup_schoenflies
     integer(c_int) :: arithmetic_crystal_class_number
     character(len=7) :: arithmetic_crystal_class_symbol
  end type SpglibSpaceGroupType

  type :: SpglibDataset
     integer(c_int) :: spacegroup_number
     integer(c_int) :: hall_number
     character(len=11) :: international_symbol
     character(len=17) :: hall_symbol
     character(len=6) :: choice
     real(c_double)  :: transformation_matrix(3,3)
     real(c_double)  :: origin_shift(3)
     integer(c_int) :: n_operations
     integer(c_int), allocatable :: rotations(:,:,:)
     real(c_double), allocatable :: translations(:,:)
     integer(c_int) :: n_atoms
     integer(c_int), allocatable :: wyckoffs(:)
     character(len=7), allocatable :: site_symmetry_symbols(:)
     integer(c_int), allocatable :: equivalent_atoms(:) !Beware mapping refers to positions starting at 0
     integer(c_int), allocatable :: mapping_to_primitive(:)
     integer(c_int) :: n_std_atoms
     real(c_double) :: std_lattice(3,3)
     integer(c_int), allocatable :: std_types(:)
     real(c_double), allocatable :: std_positions(:,:)
     real(c_double)  :: std_rotation_matrix(3,3)
     integer(c_int), allocatable :: std_mapping_to_primitive(:)
     character(len=6) :: pointgroup_symbol
     integer(kind(SPGLIB_SUCCESS)) :: spglib_error
  end type SpglibDataset

  interface
     ! ## spglib.h
     ! int spg_get_major_version(void)
     !   Returns the major version number
     function spg_get_major_version() bind(c)
       import c_int
       integer(c_int) :: spg_get_major_version
     end function spg_get_major_version

     ! ## spglib.h
     ! int spg_get_minor_version(void)
     !   Returns the minor version number
     function spg_get_minor_version() bind(c)
       import c_int
       integer(c_int) :: spg_get_minor_version
     end function spg_get_minor_version

     ! ## spglib.h
     ! int spg_get_micro_version(void)
     !   Returns the micro version number
     function spg_get_micro_version() bind(c)
       import c_int
       integer(c_int) :: spg_get_micro_version
     end function spg_get_micro_version
     
     ! ## spglib.h
     ! int spg_standardize_cell(double lattice[3][3],double position[][3],int types[],
     !                          const int num_atom,const int to_primitive,const int no_idealize,
     !                          const double symprec);
     ! Returns the standardized unit cell, possibly transforming to a
     ! primitive and idealizing the structure. Returns 0 if failed;
     ! otherwise, returns the number of atoms in the standardized
     ! cell. to_primitve=1, return the standardized primitive.
     ! lattice, position, and types arguments are overwritten.
     ! no_idealize=1 idealizes the basis vectors and atomic positions.
     function spg_standardize_cell(lattice, position, types, num_atom, &
        to_primitive, no_idealize, symprec) bind(c)
       import c_int, c_double
       real(c_double), intent(inout) :: lattice(3,3)
       real(c_double), intent(in) :: position(3,*)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       integer(c_int), intent(in), value :: to_primitive
       integer(c_int), intent(in), value :: no_idealize
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_standardize_cell
     end function spg_standardize_cell

     function spgat_standardize_cell(lattice, position, types, num_atom, &
        to_primitive, no_idealize, symprec, angle_tolerance) bind(c)
       import c_int, c_double
       real(c_double), intent(inout) :: lattice(3,3)
       real(c_double), intent(in) :: position(3,*)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       integer(c_int), intent(in), value :: to_primitive
       integer(c_int), intent(in), value :: no_idealize
       real(c_double), intent(in), value :: symprec, angle_tolerance
       integer(c_int) :: spg_standardize_cell
     end function spgat_standardize_cell
     
     ! ## spglib.h
     ! int spg_get_symmetry(int rotation[][3][3], double translation[][3], const int max_size,
     !                      SPGCONST double lattice[3][3], SPGCONST double position[][3],
     !                      const int types[], const int num_atom, const double symprec);
     !
     ! Find the symmetry operations. The operations are stored in rotation and translation.
     ! The number of operations is returned as the return value. Rotations and translations
     ! are given in fractional coordinates.
     function spg_get_symmetry( rotation, translation, max_size, lattice, &
        & position, types, num_atom, symprec) bind(c)
       import c_int, c_double
       integer(c_int), intent(inout) :: rotation(3,3,*)
       real(c_double), intent(inout) :: translation(3,*)
       integer(c_int), intent(in), value :: max_size
       real(c_double), intent(in) :: lattice(3,3), position(3,*)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_get_symmetry
     end function spg_get_symmetry

     function spgat_get_symmetry( rotation, translation, max_size, lattice, &
        & position, types, num_atom, symprec, angle_tolerance) bind(c)
       import c_int, c_double
       integer(c_int), intent(inout) :: rotation(3,3,*)
       real(c_double), intent(inout) :: translation(3,*)
       integer(c_int), intent(in), value :: max_size
       real(c_double), intent(in) :: lattice(3,3), position(3,*)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec, angle_tolerance
       integer(c_int) :: spgat_get_symmetry
     end function spgat_get_symmetry

     
     !  int spg_get_symmetry_with_collinear_spin(int rotation[][3][3], double translation[][3], int equivalent_atoms[],
     !                                           const int max_size, SPGCONST double lattice[3][3], SPGCONST double position[][3],
     !                                           const int types[], const double spins[], const int num_atom, const double symprec);
     !
     ! Same as spg_get_symmetry but with collinear spins on atoms.
     function spg_get_symmetry_with_collinear_spin( rotation, translation, &
        & max_size, lattice, position, types, spins, num_atom, symprec) bind(c)
       import c_int, c_double
       integer(c_int), intent(inout) :: rotation(3,3,*)
       real(c_double), intent(inout) :: translation(3,*)
       integer(c_int), intent(in), value :: max_size
       real(c_double), intent(in) :: lattice(3,3), position(3,*)
       integer(c_int), intent(in) :: types(*)
       real(c_double), intent(in) :: spins(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_get_symmetry_with_collinear_spin
     end function spg_get_symmetry_with_collinear_spin

     function spgat_get_symmetry_with_collinear_spin( rotation, translation, &
        & max_size, lattice, position, types, spins, num_atom, symprec, angle_tolerance) bind(c)
       import c_int, c_double
       integer(c_int), intent(inout) :: rotation(3,3,*)
       real(c_double), intent(inout) :: translation(3,*)
       integer(c_int), intent(in), value :: max_size
       real(c_double), intent(in) :: lattice(3,3), position(3,*)
       integer(c_int), intent(in) :: types(*)
       real(c_double), intent(in) :: spins(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec, angle_tolerance
       integer(c_int) :: spgat_get_symmetry_with_collinear_spin
     end function spgat_get_symmetry_with_collinear_spin

     ! int spg_get_multiplicity(SPGCONST double lattice[3][3], SPGCONST double position[][3],
     !                          const int types[], const int num_atom, const double symprec);
     ! Return the number of symmetry operations. May be used in advance to allocate memory
     ! for the operations.
     function spg_get_multiplicity( lattice, position, types, num_atom, symprec) bind(c)
       import c_int, c_double
       real(c_double), intent(in) :: lattice(3,3), position(3,*)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_get_multiplicity
     end function spg_get_multiplicity

     function spgat_get_multiplicity( lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
       import c_int, c_double
       real(c_double), intent(in) :: lattice(3,3), position(3,*)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec, angle_tolerance
       integer(c_int) :: spgat_get_multiplicity
     end function spgat_get_multiplicity

     ! int spg_get_hall_number_from_symmetry(SPGCONST int rotation[][3][3],
     !                                       SPGCONST double translation[][3],
     !                                       const int num_operations,
     !                                       const double symprec);
     ! Calculate the hall number from the symmetry operations.
     function spg_get_hall_number_from_symmetry(rotation,translation,num_operations,symprec) bind(c)
       import c_int, c_double
       integer(c_int), intent(in) :: rotation(3,3,*)
       real(c_double), intent(in) :: translation(3,*)
       integer(c_int), intent(in), value :: num_operations
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_get_hall_number_from_symmetry
     end function spg_get_hall_number_from_symmetry
     
     ! int spg_delaunay_reduce(double lattice[3][3], const double symprec);
     ! Delaunay reduction for this cell. The result is overwritten in lattice.
     function spg_delaunay_reduce(lattice, symprec) bind(c)
       import c_int, c_double
       real(c_double), intent(inout) :: lattice(3,3)
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_delaunay_reduce
     end function spg_delaunay_reduce

     ! int spg_niggle_reduce(double lattice[3][3], const double symprec);
     ! Niggle reduction for this cell. The result is overwritten in lattice.
     function spg_niggli_reduce( lattice, symprec) bind(c)
       import c_int, c_double
       real(c_double), intent(inout) :: lattice(3,3)
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_niggli_reduce
     end function spg_niggli_reduce

     ! int spg_find_primitive(double lattice[3][3],double position[][3],int types[],
     !                        const int num_atom,const double symprec);
     ! Find a primitive cell and return the lattice vectors, positions, and types.
     ! The new number of atoms is returned, or 0 if not found.
     function spg_find_primitive(lattice, position, types, num_atom, symprec) bind(c)
       import c_int, c_double
       real(c_double), intent(inout) :: lattice(3,3), position(3,*)
       integer(c_int), intent(inout) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_find_primitive
     end function spg_find_primitive

     function spgat_find_primitive(lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
       import c_int, c_double
       real(c_double), intent(inout) :: lattice(3,3), position(3,*)
       integer(c_int), intent(inout) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec, angle_tolerance
       integer(c_int) :: spgat_find_primitive
     end function spgat_find_primitive
     
     ! int spg_get_international(char symbol[11],SPGCONST double lattice[3][3],SPGCONST double position[][3],
     !                           const int types[],const int num_atom,const double symprec);
     ! Determine the space group international table symbol from the crystal geometry.
     function spg_get_international(symbol, lattice, position, types, num_atom, symprec) bind(c)
       import c_char, c_int, c_double
       character(kind=c_char), intent(out) :: symbol(11)
       real(c_double), intent(in) :: lattice(3,3), position(3, *)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_get_international ! the number corresponding to 'symbol'. 0 on failure
     end function spg_get_international

     function spgat_get_international( symbol, lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
       import c_char, c_int, c_double
       character(kind=c_char), intent(out) :: symbol(11)
       real(c_double), intent(in) :: lattice(3,3), position(3, *)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec, angle_tolerance
       integer(c_int) :: spgat_get_international ! the number corresponding to 'symbol'. 0 on failure
     end function spgat_get_international

     ! int spg_get_schoenflies(char symbol[7],SPGCONST double lattice[3][3],SPGCONST double position[][3],
     !                         const int types[],const int num_atom,const double symprec);
     ! Determine the Schoenflies space group symbol from the crystal geometry.     
     function spg_get_schoenflies( symbol, lattice, position, types, num_atom, symprec) bind(c)
       import c_char, c_int, c_double
       character(kind=c_char), intent(out) :: symbol(7)
       real(c_double), intent(in) :: lattice(3,3), position(3, *)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_get_schoenflies ! the number corresponding to 'symbol'. 0 on failure
     end function spg_get_schoenflies

     function spgat_get_schoenflies( symbol, lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
       import c_char, c_int, c_double
       character(kind=c_char), intent(out) :: symbol(7)
       real(c_double), intent(in) :: lattice(3,3), position(3, *)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec, angle_tolerance
       integer(c_int) :: spgat_get_schoenflies ! the number corresponding to 'symbol'. 0 on failure
     end function spgat_get_schoenflies
     
     ! int spg_get_pointgroup(char symbol[6],int trans_mat[3][3],
     !                        SPGCONST int rotations[][3][3],const int num_rotations);
     ! Point group symbol from the rotational part of the symmetry operations.
     function spg_get_pointgroup(symbol, trans_mat, rotations, num_rotations) bind(c)
       import c_char, c_int, c_double
       character(kind=c_char) :: symbol(6)
       integer(c_int), intent(inout) :: trans_mat(3,3)
       integer(c_int), intent(in) :: rotations(3,3,*)
       integer(c_int), intent(in), value :: num_rotations
       integer(c_int) :: spg_get_pointgroup
     end function spg_get_pointgroup
     
     ! int spg_refine_cell(double lattice[3][3],double position[][3],int types[],
     !                     const int num_atom,const double symprec);
     ! Idealize atomic positions from input lattice vectors. The arrays need to be
     ! 4 times larger than those in the input cell.
     function spg_refine_cell( lattice, position, types, num_atom, symprec) bind(c)
       import c_int, c_double
       real(c_double), intent(inout) :: lattice(3,3), position(3,*)
       integer(c_int), intent(inout) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_refine_cell
     end function spg_refine_cell

     function spgat_refine_cell( lattice, position, types, num_atom, symprec, angle_tolerance) bind(c)
       import c_int, c_double
       real(c_double), intent(inout) :: lattice(3,3), position(3,*)
       integer(c_int), intent(inout) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec, angle_tolerance
       integer(c_int) :: spgat_refine_cell
     end function spgat_refine_cell

     ! SpglibError spg_get_error_code(void);
     ! Get the error code from the last operation
     function spg_get_error_code() bind(c, name='spg_get_error_code')
       integer(kind(SPGLIB_SUCCESS)) :: spg_get_error_code
     end function spg_get_error_code

     ! char *spg_get_error_message(SpglibError spglib_error);
     ! Returns the error message based on the id from the last operation.
     module function spg_get_error_message(spglib_error)
       integer(kind(SPGLIB_SUCCESS)) :: spglib_error
       character(len=32) :: spg_get_error_message
     end function spg_get_error_message
       
     ! SpglibDataset* spg_get_dataset(SPGCONST double lattice[3][3],SPGCONST double position[][3],
     !                                const int types[], const int num_atom, const double symprec);
     ! SpglibDataset* spg_get_dataset_with_hall_number(SPGCONST double lattice[3][3],SPGCONST double position[][3],
     !                                                 const int types[],const int num_atom,const int hall_number,
     !                                                 const double symprec);
     ! Get the dataset from the crystal geometry and precision. Optionally, indicate the
     ! hall number.
     module function spg_get_dataset(lattice, position, types, num_atom, symprec, hall_number) result(dset)
       use iso_c_binding, only: c_int, c_double
       real(c_double), intent(in) :: lattice(3,3)
       real(c_double), intent(in) :: position(3,*)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int), intent(in), value, optional :: hall_number
       type(SpglibDataset) :: dset
     end function spg_get_dataset

     ! SpglibSpacegroupType spg_get_spacegroup_type(const int hall_number);
     ! Get the space group tyep information from the hall symbol
     module function spg_get_spacegroup_type(hall_number) result(tp)
       integer, intent(in) :: hall_number
       type(SpglibSpaceGroupType) :: tp
     end function spg_get_spacegroup_type

     ! Return the hall number corresponding to a space group symbol.
     ! Returns -1 if not found.
     module function spg_get_hall_number_from_symbol(symbol0)
       character(len=*), intent(in) :: symbol0
       integer(c_int) :: spg_get_hall_number_from_symbol
     end function spg_get_hall_number_from_symbol
     
     ! Return the symmetry operations from the hall number.
     module subroutine spg_get_symmetry_from_database(hnum,nrot,ncv,rot,cv)
       integer, intent(in) :: hnum
       integer, intent(out) :: nrot
       integer, intent(out) :: ncv
       real*8, intent(inout), allocatable :: rot(:,:,:)
       real*8, intent(inout), allocatable :: cv(:,:)
     end subroutine spg_get_symmetry_from_database

     ! Write a list of all known space groups to the standard output
     module subroutine spg_list_spg()
     end subroutine spg_list_spg

  end interface

end module spglib
