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

  public :: SpglibDataset, SpglibSpaceGroupType
  public :: spg_get_major_version
  public :: spg_get_minor_version
  public :: spg_get_micro_version
  public :: spg_standardize_cell
  public :: spg_get_dataset
  public :: spg_get_symmetry, spgat_get_symmetry
  public :: spg_get_symmetry_with_collinear_spin, spgat_get_symmetry_with_collinear_spin
  public :: spg_get_multiplicity, spgat_get_multiplicity
  public :: spg_get_hall_number_from_symmetry
  public :: spg_delaunay_reduce, spg_niggli_reduce
  public :: spg_find_primitive, spgat_find_primitive
  public :: spg_get_international, spgat_get_international
  public :: spg_get_schoenflies, spgat_get_schoenflies
  public :: spg_get_pointgroup
  public :: spg_refine_cell, spgat_refine_cell
  public :: spg_get_ir_reciprocal_mesh
  public :: spg_get_stabilized_reciprocal_mesh
  public :: spg_get_error_code
  public :: spg_get_error_message

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
     
     ! int spg_get_symmetry_from_database(int rotations[192][3][3],double translations[192][3],
     !                                    const int hall_number);
     ! Return the symmetry operations from the hall number.
     function spg_get_symmetry_from_database(rotations,translations,hall_number) bind(c)
       import c_int, c_double
       integer(c_int), intent(inout) :: rotations(3,3,192)
       real(c_double), intent(inout) :: translations(3,192)
       integer(c_int) :: spg_get_symmetry_from_database
     end function spg_get_symmetry_from_database
     
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


     function spg_get_ir_reciprocal_mesh(grid_point, map, mesh, &
        & is_shift, is_time_reversal, lattice, position, types, num_atom, symprec) bind(c)
       import c_int, c_double
       !   Beware the map refers to positions starting at 0
       integer(c_int), intent(out) :: grid_point(3, *), map(*) ! size is product(mesh)
       integer(c_int), intent(in) :: mesh(3), is_shift(3)
       integer(c_int), intent(in), value :: is_time_reversal
       real(c_double), intent(in) :: lattice(3,3), position(3, *)
       integer(c_int), intent(in) :: types(*)
       integer(c_int), intent(in), value :: num_atom
       real(c_double), intent(in), value :: symprec
       integer(c_int) :: spg_get_ir_reciprocal_mesh ! the number of points in the reduced mesh
     end function spg_get_ir_reciprocal_mesh


     function spg_get_stabilized_reciprocal_mesh(grid_point, map, mesh, is_shift, &
        & is_time_reversal, lattice, num_rot, rotations, num_q, qpoints) bind(c)
       import c_int, c_double
       !   Beware the map refers to positions starting at 0
       integer(c_int), intent(inout) :: grid_point(3,*), map(*)
       integer(c_int), intent(in) :: mesh(3)
       integer(c_int), intent(in) :: is_shift(3)
       integer(c_int), intent(in), value :: is_time_reversal
       real(c_double), intent(in) :: lattice(3,3)
       integer(c_int), intent(in), value :: num_rot
       integer(c_int), intent(in) :: rotations(3,3,*)
       integer(c_int), intent(in), value :: num_q
       real(c_double), intent(in) :: qpoints(3,*)
       integer(c_int) :: spg_get_stabilized_reciprocal_mesh
     end function spg_get_stabilized_reciprocal_mesh

     ! SpglibError spg_get_error_code(void);
     ! Get the error code from the last operation
     function spg_get_error_code() bind(c, name='spg_get_error_code')
       integer(kind(SPGLIB_SUCCESS)) :: spg_get_error_code
     end function spg_get_error_code

  end interface


contains

  ! char *spg_get_error_message(SpglibError spglib_error);
  ! Returns the error message based on the id from the last operation.
  function spg_get_error_message(spglib_error)
    use c_interface_module, only: c_f_string
    integer(kind(SPGLIB_SUCCESS)) :: spglib_error
    character(len=32) :: spg_get_error_message
    integer :: i

    interface
       function spg_get_error_message_c(spglib_error) bind(c, name='spg_get_error_message')
         import c_ptr, SPGLIB_SUCCESS
         integer(kind(SPGLIB_SUCCESS)), value :: spglib_error
         type(c_ptr) :: spg_get_error_message_c
       end function spg_get_error_message_c
    end interface

    call c_f_string(spg_get_error_message_c(spglib_error),spg_get_error_message)

  end function spg_get_error_message

  ! SpglibDataset* spg_get_dataset(SPGCONST double lattice[3][3],SPGCONST double position[][3],
  !                                const int types[], const int num_atom, const double symprec);
  ! SpglibDataset* spg_get_dataset_with_hall_number(SPGCONST double lattice[3][3],SPGCONST double position[][3],
  !                                                 const int types[],const int num_atom,const int hall_number,
  !                                                 const double symprec);
  ! Get the dataset from the crystal geometry and precision. Optionally, indicate the
  ! hall number.
  function spg_get_dataset(lattice, position, types, num_atom, symprec, hall_number) result(dset)
    use c_interface_module, only: c_f_string
    real(c_double), intent(in) :: lattice(3,3)
    real(c_double), intent(in) :: position(3,*)
    integer(c_int), intent(in) :: types(*)
    integer(c_int), intent(in), value :: num_atom
    real(c_double), intent(in), value :: symprec
    integer(c_int), intent(in), value, optional :: hall_number
    type(SpglibDataset) :: dset

    type, bind(c) :: SpglibDataset_c
       integer(c_int) :: spacegroup_number
       integer(c_int) :: hall_number
       character(kind=c_char) :: international_symbol(11)
       character(kind=c_char) :: hall_symbol(17)
       character(kind=c_char) :: choice(6)
       real(c_double) :: transformation_matrix(3,3)
       real(c_double) :: origin_shift(3)
       integer(c_int) :: n_operations
       type(c_ptr) :: rotations
       type(c_ptr) :: translations
       integer(c_int) :: n_atoms
       type(c_ptr) :: wyckoffs
       type(c_ptr) :: site_symmetry_symbols
       type(c_ptr) :: equivalent_atoms
       type(c_ptr) :: mapping_to_primitive
       integer(c_int) :: n_std_atoms
       real(c_double) :: std_lattice(3,3)
       type(c_ptr) :: std_types
       type(c_ptr) :: std_positions
       real(c_double)  :: std_rotation_matrix(3,3)
       type(c_ptr) :: std_mapping_to_primitive
       character(kind=c_char) :: pointgroup_symbol(6)
    end type SpglibDataset_c

    interface
       function spg_get_dataset_c(lattice, position, types, num_atom, symprec) bind(c, name='spg_get_dataset')
         import c_int, c_double, c_ptr
         real(c_double), intent(in) :: lattice(3,3)
         real(c_double), intent(in) :: position(3,*)
         integer(c_int), intent(in) :: types(*)
         integer(c_int), intent(in), value :: num_atom
         real(c_double), intent(in), value :: symprec
         type(c_ptr) :: spg_get_dataset_c
       end function spg_get_dataset_c
       function spg_get_dataset_with_hall_number_c(lattice, position, types, num_atom, hall_number, symprec) bind(c, name='spg_get_dataset_with_hall_number')
         import c_int, c_double, c_ptr
         real(c_double), intent(in) :: lattice(3,3)
         real(c_double), intent(in) :: position(3,*)
         integer(c_int), intent(in) :: types(*)
         integer(c_int), intent(in), value :: num_atom
         integer(c_int), intent(in), value :: hall_number
         real(c_double), intent(in), value :: symprec
         type(c_ptr) :: spg_get_dataset_with_hall_number_c
       end function spg_get_dataset_with_hall_number_c

       subroutine spg_free_dataset_c(dataset) bind(c, name = 'spg_free_dataset')
         import SpglibDataset_c
         type(SpglibDataset_c), intent(inout) :: dataset
       end subroutine spg_free_dataset_c
    end interface

    type(SpglibDataset_c), pointer :: dset_c
    type(c_ptr) :: dataset_ptr_c
    integer(c_int) :: n_operations, n_atoms, n_std_atoms
    integer :: i
    integer(c_int) :: hall
    integer(kind(SPGLIB_SUCCESS)) :: SpglibErrcode
    real(c_double), pointer :: translations(:,:)
    integer(c_int), pointer :: rotations(:,:,:), wyckoffs(:), equivalent_atoms(:), std_types(:), std_positions(:,:)

    hall = 0
    if (present(hall_number)) hall = hall_number
    if (hall > 0) then
       dataset_ptr_c = spg_get_dataset_with_hall_number_c(lattice, position, types, num_atom, hall, symprec)
    else
       dataset_ptr_c = spg_get_dataset_c(lattice, position, types, num_atom, symprec)
    end if

    if (c_associated(dataset_ptr_c)) then
       dset%spglib_error = SPGLIB_SUCCESS

       call c_f_pointer(dataset_ptr_c, dset_c)

       dset%spacegroup_number     = dset_c%spacegroup_number
       dset%hall_number           = dset_c%hall_number
       dset%transformation_matrix = dset_c%transformation_matrix
       dset%origin_shift          = dset_c%origin_shift
       dset%n_operations          = dset_c%n_operations
       dset%n_atoms               = dset_c%n_atoms
       dset%n_std_atoms           = dset_c%n_std_atoms
       dset%std_lattice           = dset_c%std_lattice
       call c_f_string(dset_c%international_symbol,dset%international_symbol)
       call c_f_string(dset_c%hall_symbol,dset%hall_symbol)
       call c_f_string(dset_c%choice,dset%choice)
       call c_f_string(dset_c%pointgroup_symbol,dset%pointgroup_symbol)
       n_operations = dset_c%n_operations
       n_atoms      = dset_c%n_atoms
       n_std_atoms  = dset_c%n_std_atoms
       call c_f_pointer(dset_c%rotations, rotations, shape=[3, 3, n_operations])
       call c_f_pointer(dset_c%translations, translations, shape=[3, n_operations])
       call c_f_pointer(dset_c%wyckoffs, wyckoffs, shape=[n_atoms])
       call c_f_pointer(dset_c%equivalent_atoms, equivalent_atoms, shape=[n_atoms])
       call c_f_pointer(dset_c%std_types, std_types, shape=[n_std_atoms])
       call c_f_pointer(dset_c%std_positions, std_positions, shape=[3, n_std_atoms])
       allocate(dset%rotations(3, 3, n_operations))
       allocate(dset%translations(3, n_operations))
       allocate(dset%wyckoffs(n_atoms))
       allocate(dset%equivalent_atoms(n_atoms))
       allocate(dset%std_types(n_std_atoms))
       allocate(dset%std_positions(3, n_std_atoms))
       dset%rotations        = rotations
       dset%translations     = translations
       dset%wyckoffs         = wyckoffs
       dset%equivalent_atoms = equivalent_atoms
       dset%std_types        = std_types
       dset%std_positions    = std_positions
       call spg_free_dataset_c(dset_c)
    else
       dset%spglib_error = spg_get_error_code()
       dset%spacegroup_number = 0
       dset%hall_number = 0
       dset%international_symbol = ""
       dset%hall_symbol = ""
       dset%choice = ""
       dset%transformation_matrix = 0.0_c_double
       dset%origin_shift = 0.0_c_double
       dset%n_operations = 0
       dset%n_atoms = 0
       dset%n_std_atoms = 0
       dset%std_lattice = 0.0_c_double
       dset%pointgroup_symbol = ""
    end if

  end function spg_get_dataset

  ! SpglibSpacegroupType spg_get_spacegroup_type(const int hall_number);
  ! Get the space group tyep information from the hall symbol
  function spg_get_spacegroup_type(hall_number) result(tp)
    use c_interface_module, only: c_f_string
    integer, intent(in) :: hall_number
    type(SpglibSpaceGroupType) :: tp
    
    type, bind(c) :: SpglibSpaceGroupType_c
       integer(c_int) :: number
       character(kind=c_char) :: international_short(11)
       character(kind=c_char) :: international_full(20)
       character(kind=c_char) :: international(32)
       character(kind=c_char) :: schoenflies(7)
       character(kind=c_char) :: hall_symbol(17)
       character(kind=c_char) :: choice(6)
       character(kind=c_char) :: pointgroup_international(6)
       character(kind=c_char) :: pointgroup_schoenflies(4)
       integer(c_int) :: arithmetic_crystal_class_number
       character(kind=c_char) :: arithmetic_crystal_class_symbol(7)
    end type SpglibSpaceGroupType_c

    interface
       function spg_get_spacegroup_type_c(hall_number) bind(c, name='spg_get_spacegroup_type')
         import c_int, SpglibSpaceGroupType_c
         integer(c_int), intent(in), value :: hall_number
         type(SpglibSpaceGroupType_c) :: spg_get_spacegroup_type_c
       end function spg_get_spacegroup_type_c
    end interface

    type(SpglibSpaceGroupType_c) :: tpc
    
    tpc = spg_get_spacegroup_type_c(hall_number)
    tp%number = tpc%number
    tp%arithmetic_crystal_class_number = tpc%arithmetic_crystal_class_number
    call c_f_string(tpc%international_short,tp%international_short)
    call c_f_string(tpc%international_full,tp%international_full)
    call c_f_string(tpc%international,tp%international)
    call c_f_string(tpc%schoenflies,tp%schoenflies)
    call c_f_string(tpc%hall_symbol,tp%hall_symbol)
    call c_f_string(tpc%choice,tp%choice)
    call c_f_string(tpc%pointgroup_international,tp%pointgroup_international)
    call c_f_string(tpc%pointgroup_schoenflies,tp%pointgroup_schoenflies)
    call c_f_string(tpc%arithmetic_crystal_class_symbol,tp%arithmetic_crystal_class_symbol)
    
  end function spg_get_spacegroup_type

end module spglib
