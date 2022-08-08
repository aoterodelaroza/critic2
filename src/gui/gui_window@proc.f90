! Copyright (c) 2019 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
submodule (gui_window) proc
  use gui_interfaces_cimgui
  implicit none

  ! Count ID for keeping track of windows and widgets
  integer :: idcount = 0

contains

  !> Initialize a window of the given type. If isiopen, initialize it
  !> as open.
  module subroutine window_init(w,type,isopen)
    class(window), intent(inout) :: w
    integer, intent(in) :: type
    logical, intent(in) :: isopen

    w%isinit = .true.
    w%isopen = isopen
    w%type = type
    w%id = -1

  end subroutine window_init

  !> Draw an ImGui window.
  module subroutine window_draw(w)
    use tools_io, only: string
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str

    ! First pass: assign ID, name, and flags
    if (w%id < 0) then
       idcount = idcount + 1
       w%id = idcount
       if (w%type == wintype_tree) then
          ! w%name = "Tree###" // string(w%id) // c_null_char
          w%name = "Tree" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_view) then
          ! w%name = "View###" // string(w%id) // c_null_char
          w%name = "View" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_console) then
          ! w%name = "Console###" // string(w%id) // c_null_char
          w%name = "Console" // c_null_char
          w%flags = ImGuiWindowFlags_None
       end if
    end if

    if (w%isopen) then
       if (igBegin(c_loc(w%name),w%isopen,w%flags)) then
          ! assign the pointer ID
          w%ptr = igGetCurrentWindow()

          ! draw the window contents, depending on type
          if (w%type == wintype_tree) then
             str = "Hello Tree!"
             call igText(c_loc(str))
          elseif (w%type == wintype_view) then
             str = "Hello View!"
             call igText(c_loc(str))
          elseif (w%type == wintype_console) then
             str = "Hello Console!"
             call igText(c_loc(str))
          end if
       end if
       call igEnd()
    end if

  end subroutine window_draw

  !xx! private procedures

end submodule proc
