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

! interface to the spglib library
submodule (spglib) proc
  use hashmod, only: hash
  implicit none

  ! save the mapping between international symbol, etc.
  ! to hall number
  logical :: mapavail = .false.
  type(hash) :: ints, intf
  integer :: inthnum(230)

contains

  ! char *spg_get_error_message(SpglibError spglib_error);
  ! Returns the error message based on the id from the last operation.
  module function spg_get_error_message(spglib_error)
    use c_interface_module, only: c_f_string
    integer(kind(SPGLIB_SUCCESS)) :: spglib_error
    character(len=32) :: spg_get_error_message

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
    real(c_double), intent(inout) :: lattice(3,3)
    real(c_double), intent(inout) :: position(3,*)
    integer(c_int), intent(inout) :: types(*)
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
         real(c_double), intent(inout) :: lattice(3,3)
         real(c_double), intent(inout) :: position(3,*)
         integer(c_int), intent(inout) :: types(*)
         integer(c_int), intent(in), value :: num_atom
         real(c_double), intent(in), value :: symprec
         type(c_ptr) :: spg_get_dataset_c
       end function spg_get_dataset_c
       function spg_get_dataset_with_hall_number_c(lattice,position,types,num_atom,hall_number,symprec)&
          bind(c, name='spg_get_dataset_with_hall_number')
         import c_int, c_double, c_ptr
         real(c_double), intent(inout) :: lattice(3,3)
         real(c_double), intent(inout) :: position(3,*)
         integer(c_int), intent(inout) :: types(*)
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
    integer(c_int) :: hall
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

  ! Return the hall number corresponding to a space group symbol.
  ! Returns -1 if not found.
  module function spg_get_hall_number_from_symbol(symbol0) result(hnum)
    use tools_io, only: isinteger, lgetword, equal, lower, deblank
    character(len=*), intent(in) :: symbol0
    integer(c_int) :: hnum

    integer :: iaux, lp
    character(len=:), allocatable :: word, symbol

    call build_hall_mapping()

    hnum = -1
    lp = 1
    if (isinteger(iaux,symbol0,lp)) then
       word = lgetword(symbol0,lp)
       if (equal(word,'hm')) then
          if (iaux >= 1 .and. iaux <= 230) hnum = inthnum(iaux)
       else
          if (iaux >= 1 .and. iaux <= 530) hnum = iaux
       end if
    else
       symbol = trim(deblank(symbol0,todn=.true.))
       if (intf%iskey(symbol)) then
          hnum = intf%get(symbol,hnum)
       elseif (ints%iskey(symbol)) then
          hnum = ints%get(symbol,hnum)
       end if
    end if

  end function spg_get_hall_number_from_symbol

  !xx! private procedures

  ! Build the mapping that gives the Hall number from symbols, etc.
  subroutine build_hall_mapping()
    use tools_io, only: deblank, stripchar
    type(SpglibSpaceGroupType) :: sa
    integer :: i, iaux
    character(len=:), allocatable :: aux

    if (mapavail) return

    do i = 1, 530
       sa = spg_get_spacegroup_type(i)

       aux = deblank(sa%international_short,todn=.true.)
       aux = trim(stripchar(aux,"_"))
       if (.not.ints%iskey(aux)) call ints%put(aux,i)

       aux = deblank(sa%international_full,todn=.true.)
       aux = trim(stripchar(aux,"_"))
       if (.not.intf%iskey(aux)) call intf%put(aux,i)

       iaux = sa%number
       if (inthnum(iaux) == 0) inthnum(iaux) = i
    end do
    mapavail = .true.

  end subroutine build_hall_mapping

  ! Return the symmetry operations from the Hall number. If failed,
  ! return nrot = ncv = 0.
  module subroutine spg_get_symmetry_from_database(hnum,nrot,ncv,rot,cv)
    use types, only: realloc
    integer, intent(in) :: hnum
    integer, intent(out) :: nrot
    integer, intent(out) :: ncv
    real*8, intent(inout), allocatable :: rot(:,:,:)
    real*8, intent(inout), allocatable :: cv(:,:)

    integer(c_int) :: rot0(3,3,192)
    real(c_double) :: cv0(3,192)
    integer :: nop, i, j, k
    logical :: found
    real*8 :: xdif(3), xaux(3), xmaux(3,4)

    real*8, parameter :: eps = 1d-5
    integer, parameter :: eye(3,3) = reshape((/1,0,0,0,1,0,0,0,1/),shape(eye))

    interface
       ! int spg_get_symmetry_from_database(int rotations[192][3][3],double translations[192][3],
       !                                    const int hall_number);
       function spg_get_symmetry_from_database_c(rotations,translations,hall_number) bind(c,name='spg_get_symmetry_from_database')
         import c_int, c_double
         integer(c_int), intent(inout) :: rotations(3,3,192)
         real(c_double), intent(inout) :: translations(3,192)
         integer(c_int), value :: hall_number
         integer(c_int) :: spg_get_symmetry_from_database_c
       end function spg_get_symmetry_from_database_c
    end interface

    nrot = 0
    ncv = 0
    if (allocated(rot)) deallocate(rot)
    if (allocated(cv)) deallocate(cv)

    nop = spg_get_symmetry_from_database_c(rot0,cv0,hnum)
    if (nop <= 0) goto 999

    ! detect centering vectors
    allocate(cv(3,10))
    do i = 1, nop
       if (all(rot0(:,:,i) - eye == 0)) then
          ncv = ncv + 1
          if (ncv > size(cv,2)) call realloc(cv,3,2*ncv)
          cv(:,ncv) = cv0(:,i)
       end if
    end do
    call realloc(cv,3,ncv)

    ! allocate space for symmetry operations
    if (mod(nop,ncv) /= 0) goto 999

    ! reduce the number of operations using the centering vectors
    nrot = 0
    allocate(rot(3,4,nop/ncv))
    do i = 1, nop
       found = .false.
       jloop: do j = 1, nrot
          if (all(nint(rot(1:3,1:3,j)) == rot0(:,:,i))) then
             do k = 1, ncv
                xdif = (rot(1:3,4,j) + cv(:,k)) - cv0(:,i)
                xdif = xdif - nint(xdif)
                if (all(abs(xdif) < eps)) then
                   found = .true.
                   exit jloop
                end if
             end do
          end if
       end do jloop
       if (.not.found) then
          nrot = nrot + 1
          rot(1:3,1:3,nrot) = rot0(:,:,i)
          rot(1:3,4,nrot) = cv0(:,i)
       end if
    end do

    ! tranpose the rotation matrices
    do i = 1, nrot
       rot(1:3,1:3,i) = transpose(rot(1:3,1:3,i))
    end do

    ! check that the first centering is the zero vector
    if (.not.all(abs(cv(:,1)) < eps) .and. ncv > 1) then
       k = 0
       do j = 2, ncv
          if (all(abs(cv(:,j)) < eps)) then
             k = j
             exit
          end if
       end do
       if (k == 0) goto 999
       xaux = cv(:,k)
       cv(:,k) = cv(:,1)
       cv(:,1) = xaux
    end if

    ! check that the first operation is the identity
    if (.not.all(abs(rot(1:3,1:3,1) - eye) < eps) .and. all(abs(rot(1:3,4,1)) < eps) .and. nrot > 1) then
       k = 0
       do j = 2, nrot
          if (all(abs(rot(1:3,1:3,j) - eye) < eps) .and. all(abs(rot(1:3,4,j)) < eps)) then
             k = j
             exit
          end if
       end do
       if (k == 0) goto 999
       xmaux = rot(:,:,k)
       rot(:,:,k) = rot(:,:,1)
       rot(:,:,1) = xmaux
    end if

    return
999 continue
    ncv = 0
    nrot = 0
    if (allocated(cv)) deallocate(cv)
    if (allocated(rot)) deallocate(rot)

  end subroutine spg_get_symmetry_from_database

  ! Write a list of all known space groups to the standard output
  module subroutine spg_list_spg()
    use tools_io, only: deblank, stripchar, string, uout, ioj_left
    type(SpglibSpaceGroupType) :: sa
    integer :: i

    character(len=:), allocatable :: strs, strf, strcs

    write (uout,'("* LIST of space group types")')
    write (uout,'("# Hall = Hall space group number. ITA = International space group number.")')
    write (uout,'("# HM-short = short Hermann-Mauguin symbol. HM-long = long H-M symbol.")')
    write (uout,'("# choice = origin/setting choice. crys.-syst. = crystal system.")')
    write (uout,'("# Hall-symbol = Hall space group symbol.")')
    write (uout,'("#Hall ITA   HM-short HM-long       choice  crys.-syst.  Hall-symbol")')
    do i = 1, 530
       sa = spg_get_spacegroup_type(i)

       strs = deblank(sa%international_short)
       strs = trim(stripchar(strs,"_"))
       strf = deblank(sa%international_full)
       strf = trim(stripchar(strf,"_"))

       if (sa%number >= 1 .and. sa%number <= 2) then
          strcs = "triclinic"
       elseif (sa%number >= 3 .and. sa%number <= 15) then
          strcs = "monoclinic"
       elseif (sa%number >= 16 .and. sa%number <= 74) then
          strcs = "orthorhombic"
       elseif (sa%number >= 75 .and. sa%number <= 142) then
          strcs = "tetragonal"
       elseif (sa%number >= 143 .and. sa%number <= 167) then
          strcs = "trigonal"
       elseif (sa%number >= 168 .and. sa%number <= 194) then
          strcs = "hexagonal"
       elseif (sa%number >= 195 .and. sa%number <= 230) then
          strcs = "cubic"
       else
          strcs = "??"
       end if

       write (uout,'(6(A,X),"[",A,"]")') string(i,5,ioj_left), string(sa%number,5,ioj_left),&
          string(strs,8,ioj_left), string(strf,14,ioj_left), string(sa%choice,6,ioj_left), &
          string(strcs,12,ioj_left), string(sa%hall_symbol)
    end do
    write (uout,*)

  end subroutine spg_list_spg

end submodule proc
