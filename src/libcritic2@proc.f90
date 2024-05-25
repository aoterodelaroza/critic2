! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Implementation of the critic2 C/C++ library.
submodule (libcritic2) proc
  use iso_c_binding
  implicit none

  ! whether critic2 has been initialized
  logical :: critic2_init = .false.

  !xx! private procedures
  ! subroutine initialize_critic2()

contains

  !> Read a crystal structure from a file. Allocate space for the
  !> crystal structure and return a pointer to it or NULL if there was
  !> an error. Requires destruction of the object after its use.
  module function create_structure_from_file(file) bind(c,name="create_structure_from_file")
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use c_interface_module, only: c_f_string_alloc
    type(c_ptr), value, intent(in) :: file
    type(c_ptr) :: create_structure_from_file

    character(len=:), allocatable :: fname
    character(len=1,kind=c_char), dimension(:), pointer :: p_chars
    type(crystalseed) :: seed
    integer :: i
    character(len=:), allocatable :: errmsg
    type(crystal), pointer :: crf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()

    ! read seed from file
    create_structure_from_file = c_null_ptr
    call c_f_string_alloc(file,fname)
    call seed%read_any_file(fname,-1,errmsg)
    if (len_trim(errmsg) > 0) return

    ! output the crystal
    allocate(crf)
    call crf%struct_new(seed,.false.)
    if (.not.crf%isinit) return
    create_structure_from_file = c_loc(crf)

  end function create_structure_from_file

  !> Write the report about a crystal structure to standard output.
  module subroutine describe_structure(cr) bind(c,name="describe_structure")
    use crystalmod, only: crystal
    type(c_ptr), value, intent(in) :: cr

    type(crystal), pointer :: crf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(cr)) return
    call c_f_pointer(cr,crf)
    if (.not.associated(crf)) return

    ! write the report
    call crf%report(.true.,.false.)

  end subroutine describe_structure

  !> Destroy a crystal structure object.
  module subroutine destroy_structure(cr) bind(c,name="destroy_structure")
    use crystalmod, only: crystal
    type(c_ptr), value, intent(in) :: cr

    type(crystal), pointer :: crf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(cr)) return
    call c_f_pointer(cr,crf)
    if (.not.associated(crf)) return

    ! destroy the crystal object and free memory
    call crf%end()
    deallocate(crf)

  end subroutine destroy_structure

  !xx! private procedures

  !> Initialize critic2 when used as a library.
  subroutine initialize_critic2()
    use spglib, only: spg_build_hall_mapping
    use systemmod, only: systemmod_init
    use global, only: global_init, fileroot, initial_banner, config_write
    use tools_io, only: start_clock, stdargs, history_init
    use param, only: param_init
    use config, only: getstring, istring_datadir

    character(len=:), allocatable :: optv, ghome

    ! initialize
    call start_clock()
    call param_init()

    call stdargs(optv,ghome,fileroot)
    call history_init()

    call global_init(ghome,getstring(istring_datadir))
    call systemmod_init(1)
    call spg_build_hall_mapping() ! for spglib
    critic2_init = .true.

  end subroutine initialize_critic2

end submodule proc
