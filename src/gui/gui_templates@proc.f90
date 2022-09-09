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

! Some utilities for building the GUI (e.g. wrappers around ImGui routines).
submodule (gui_templates) proc
  use iso_c_binding
  use hashmod, only: hash
  use param, only: newline
  implicit none

  ! the keywords
  !- structural tools
  integer, parameter :: ikeyw_none = 0
  integer, parameter :: ikeyw_bz = 1          ! BZ
  integer, parameter :: ikeyw_econ = 2        ! ECON
  integer, parameter :: ikeyw_environ = 3     ! ENVIRON
  integer, parameter :: ikeyw_kpoints = 4     ! KPOINTS
  integer, parameter :: ikeyw_spg = 5         ! SPG
  integer, parameter :: ikeyw_sym = 6         ! SYM/SYMM/NOSYM/NOSYMM
  !- fields
  integer, parameter :: ikeyw_load = 7        ! LOAD
  integer, parameter :: ikeyw_reference = 8   ! REFERENCE
  integer, parameter :: ikeyw_setfield = 9    ! SETFIELD
  integer, parameter :: ikeyw_unload = 10     ! UNLOAD
  !- read and write files
  integer, parameter :: ikeyw_makemolsnc = 11 ! MAKEMOLSNC
  integer, parameter :: ikeyw_write = 12      ! WRITE
  !- misc
  integer, parameter :: ikeyw_libxc = 13      ! LIBXC
  integer, parameter :: ikeyw_molcell = 14    ! MOLCELL
  integer, parameter :: ikeyw_units = 15      ! UNITS
  integer, parameter :: ikeyw_zpsp = 16       ! ZPSP/Q/QAT/NOCORE
  !- variables
  integer, parameter :: ikeyw_list = 17       ! LIST
  integer, parameter :: ikeyw_clear = 18      ! CLEAR
  integer, parameter :: ikeyw_NUM = 18

  ! keyword sections (need to be sequential)
  integer, parameter :: isection_none = 0
  integer, parameter :: isection_structural_tools = 1 ! structural tools
  integer, parameter :: isection_fields = 2           ! fields
  integer, parameter :: isection_readwrite = 3        ! read & write files
  integer, parameter :: isection_variables = 4        ! variables
  integer, parameter :: isection_misc = 5             ! miscellaneous
  integer, parameter :: isection_NUM = 5
  integer, parameter :: ikeyw_section(ikeyw_NUM) = (/&
     isection_structural_tools,& ! BZ
     isection_structural_tools,& ! ECON
     isection_structural_tools,& ! ENVIRON
     isection_structural_tools,& ! KPOINTS
     isection_structural_tools,& ! SPG
     isection_structural_tools,& ! SYM
     isection_fields,&           ! LOAD
     isection_fields,&           ! SETFIELD
     isection_fields,&           ! UNLOAD
     isection_fields,&           ! REFERENCE
     isection_readwrite,&        ! MAKEMOLSNC
     isection_readwrite,&        ! WRITE
     isection_misc,&             ! LIBXC
     isection_misc,&             ! MOLCELL
     isection_misc,&             ! UNITS
     isection_misc,&             ! ZPSP/Q/QAT/NOCORE
     isection_variables,&        ! LIST
     isection_variables&         ! CLEAR
     /)

  ! keyword titles
  character(len=*,kind=c_char), parameter :: keyword_titles(ikeyw_NUM) = (/&
     "BZ (print Brillouin zone)              ",& ! BZ
     "ECON (effective coordination numbers)  ",& ! ECON
     "ENVIRON (calculate atomic environments)",& ! ENVIRON
     "KPOINTS (calculate k-point grid sizes) ",& ! KPOINTS
     "SPG (list space group types)           ",& ! SPG
     "SYM (symmetry analysis & refinement)   ",& ! SYM
     "LOAD (load a field)                    ",& ! LOAD
     "REFERENCE (set the reference field)    ",& ! REFERENCE
     "SETFIELD (set field options)           ",& ! SETFIELD
     "UNLOAD (unload a field)                ",& ! UNLOAD
     "MAKEMOLSNC (write DMACRYS mols file)   ",& ! MAKEMOLSNC
     "WRITE (write structure to a file)      ",& ! WRITE
     "LIBXC (list xc functionals)            ",& ! LIBXC
     "MOLCELL (change molecular cell)        ",& ! MOLCELL
     "UNITS (change distance units in output)",& ! UNITS
     "ZPSP/Q (atomic & pseudo charges)       ",& ! ZPSP/Q/QAT/NOCORE
     "LIST (list variables)                  ",& ! LIST
     "CLEAR (clear variables)                "& ! CLEAR
     /)

  ! section titles
  character(len=*,kind=c_char), parameter :: section_titles(isection_NUM) = (/&
     "Structural Tools  ",& ! structural tools
     "Load/Unload fields",& ! fields
     "Read/Write Files  ",& ! read/write files
     "Variables         ",& ! variables
     "Miscellaneous     "&  ! miscellaneous
     /)

  ! section ranges
  integer, parameter :: section_ranges(2,isection_NUM) = reshape((/&
     1,6,&   ! structural tools
     7,10,&  ! fields
     11,12,& ! read/write files
     13,16,& ! variables
     17,18&  ! miscellaneous
     /),shape(section_ranges))

  ! template (keyw) files
  character*(*), parameter :: template_file(ikeyw_NUM) = (/&
     "bz        ",& ! BZ
     "econ      ",& ! ECON
     "environ   ",& ! ENVIRON
     "kpoints   ",& ! KPOINTS
     "spg       ",& ! SPG
     "sym       ",& ! SYM
     "load      ",& ! LOAD
     "reference ",& ! REFERENCE
     "setfield  ",& ! SETFIELD
     "unload    ",& ! UNLOAD
     "makemolsnc",& ! WRITE
     "write     ",& ! WRITE
     "libxc     ",& ! LIBXC
     "molcell   ",& ! MOLCELL
     "units     ",& ! UNITS
     "zpsp      ",& ! ZPSP/Q/QAT/NOCORE
     "list      ",& ! LIST
     "clear     "& ! CLEAR
     /)

  ! documentation (md) files
  character*(*), parameter :: doclink(ikeyw_NUM) = (/&
     "structure/#c2-bz     ",& ! BZ
     "structure/#c2-econ   ",& ! ECON
     "structure/#c2-environ",& ! ENVIRON
     "structure/#c2-kpoints",& ! KPOINTS
     "crystal/#c2-spg      ",& ! SPG
     "crystal/#c2-symm     ",& ! SYM
     "fields/#c2-load      ",& ! LOAD
     "fields/#c2-reference ",& ! REFERENCE
     "fields/#c2-setfield  ",& ! SETFIELD
     "fields/#c2-unload    ",& ! UNLOAD
     "write/#c2-makemolsnc ",& ! MAKEMOLSNC
     "write/#c2-write      ",& ! WRITE
     "arithmetics/#libxc   ",& ! LIBXC
     "molecule/#c2-molcell ",& ! MOLCELL
     "inputoutput/#c2-units",& ! UNITS
     "crystal/#c2-charge   ",& ! ZPSP/Q/QAT/NOCORE
     "arithmetics/#c2-list ",& ! LIST
     "arithmetics/#c2-clear"&  ! CLEAR
     /)

  ! template hash
  type(hash) :: thash

  !xx! private procedures

contains

  !> Draw the keyword context menu in the templates and help buttons
  !> of the consolie input window. If textinsert is true, insert the
  !> tempalte, otherwise bring up the help.
  module subroutine draw_keyword_context_menu(textinsert)
    use gui_interfaces_cimgui
    logical, intent(in) :: textinsert

    character(kind=c_char,len=:), allocatable, target :: str1, str2

    integer :: i, j

    do i = 1, isection_NUM
       str1 = trim(section_titles(i)) // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          do j = section_ranges(1,i), section_ranges(2,i)
             str2 = trim(keyword_titles(j)) // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                call launch_keyword_action(textinsert,j)
             end if
          end do

          call igEndMenu()
       end if
    end do
    call igEndPopup()

  end subroutine draw_keyword_context_menu

  !xx! private procedures

  !> Launch the keyword action for the given keyword ikeyw. If
  !> textinsert, insert the template into the console
  !> input. Otherwise, open the documentation.
  subroutine launch_keyword_action(textinsert,ikeyw)
    use gui_window, only: win, iwin_console_input, stack_create_window
    use gui_interfaces_cimgui
    logical, intent(in) :: textinsert
    integer, intent(in) :: ikeyw

    character(kind=c_char,len=:), allocatable, target :: str

    if (textinsert) then
       str = get_keyword_template(ikeyw)
       call win(iwin_console_input)%fill_input_ci(str)
    else
       str = "https://aoterodelaroza.github.io/critic2/manual/" // trim(doclink(ikeyw)) //&
          c_null_char
       call openLink(c_loc(str))
    end if

  end subroutine launch_keyword_action

  !> Returns the template for the given keyword.
  function get_keyword_template(ikeyw)
    use global, only: critic_home
    use tools_io, only: fopen_read, fclose, getline_raw
    use param, only: dirsep, newline
    integer, intent(in) :: ikeyw
    character(kind=c_char,len=:), allocatable :: get_keyword_template

    character(kind=c_char,len=:), allocatable :: file, keyw, line
    logical :: exist
    integer :: lu

    keyw = trim(template_file(ikeyw))
    if (.not.thash%iskey(keyw)) then
       ! check the file exists
       file = trim(critic_home) // dirsep // "helpdoc" // dirsep // keyw // ".keyw"
       inquire(file=file,exist=exist)
       if (.not.exist) then
          get_keyword_template = &
             "## Template file not found. Please make sure the CRITIC_HOME is set or" // newline //&
             "## that critic2 is installed propertly."
          return
       end if

       ! read the file
       get_keyword_template = ""
       lu = fopen_read(file)
       do while (getline_raw(lu,line,.false.))
          get_keyword_template = get_keyword_template // line // newline
       end do
       call fclose(lu)

       ! save it to the hash
       call thash%put(keyw,get_keyword_template)
    else
       ! get it from the hash
       get_keyword_template = thash%get(keyw,'a')
    end if

  end function get_keyword_template

end submodule proc
