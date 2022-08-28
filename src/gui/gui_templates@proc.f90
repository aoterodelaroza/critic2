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
  integer, parameter :: ikeyw_none = 0
  integer, parameter :: ikeyw_kpoints = 1 ! KPOINTS
  integer, parameter :: ikeyw_NUM = 1

  ! keyword sections (need to be sequential)
  integer, parameter :: isection_none = 0
  integer, parameter :: isection_structural_tools = 1 ! structural tools
  integer, parameter :: isection_NUM = 1
  integer, parameter :: ikeyw_section(ikeyw_NUM) = (/&
     isection_structural_tools&
     /)

  ! keyword titles
  character(len=*,kind=c_char), parameter :: keyword_titles(ikeyw_NUM) = (/&
     "KPOINTS (Calculate k-Point Grid Sizes)"& ! KPOINTS
     /)

  ! section titles
  character(len=*,kind=c_char), parameter :: section_titles(isection_NUM) = (/&
     "Structural Tools"& ! structural tools
     /)

  ! section ranges
  integer, parameter :: section_ranges(2,1) = reshape((/&
     1,1& ! structural tools
     /),shape(section_ranges))

  ! template (keyw) files
  character*(*), parameter :: template_file(ikeyw_NUM) = (/&
     "kpoints"& ! KPOINTS
     /)

  ! documentation (md) files
  character*(*), parameter :: doclink(ikeyw_NUM) = (/&
     "structure/#c2-kpoints"& ! KPOINTS
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
       str1 = section_titles(i) // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          do j = section_ranges(1,i), section_ranges(2,i)
             str2 = keyword_titles(j) // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                call launch_keyword_action(textinsert,j)
             end if
          end do

          call igEndMenu()
       end if

       call igEndPopup()
    end do

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
       str = "https://aoterodelaroza.github.io/critic2/manual/" // doclink(ikeyw) //&
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

    keyw = template_file(ikeyw)
    if (.not.thash%iskey(keyw)) then
       ! check the file exists
       file = trim(critic_home) // dirsep // "helpdoc" // dirsep //&
          template_file(ikeyw) // ".keyw"
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
