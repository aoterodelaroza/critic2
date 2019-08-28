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

submodule (spglib) proc
  implicit none

contains

  ! char *spg_get_error_message(SpglibError spglib_error);
  ! Returns the error message based on the id from the last operation.
  module function spg_get_error_message(spglib_error)
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
  module function spg_get_dataset(lattice, position, types, num_atom, symprec, hall_number) result(dset)
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
  module function spg_get_spacegroup_type(hall_number) result(tp)
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

end submodule proc
