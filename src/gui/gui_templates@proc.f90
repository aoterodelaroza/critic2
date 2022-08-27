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
  use param, only: newline
  implicit none

  character(len=*,kind=c_char), parameter :: template_kpoints = &
     "KPOINTS [rk.r] [RKMAX rkmax.r]" // newline

contains

  !> Draw the keyword context menu in the templates and help buttons
  !> of the consolie input window.
  module subroutine draw_keyword_context_menu()
    use gui_interfaces_cimgui
    use gui_window, only: win, iwin_console_input
    character(kind=c_char,len=:), allocatable, target :: str1, str2


    ! structural tools
    str1 = "Structural Tools" // c_null_char
    if (igBeginMenu(c_loc(str1),.true._c_bool)) then
       str2 = "KPOINTS (Calculate k-Point Grid Sizes)" // c_null_char
       if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
          call win(iwin_console_input)%fill_input_ci(template_kpoints)
       end if

       call igEndMenu()
    end if ! structural tools menu

    call igEndPopup()

  end subroutine draw_keyword_context_menu

end submodule proc
