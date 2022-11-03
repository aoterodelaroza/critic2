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

! This module handles the key-bindings for the critic2 GUI.
module keybindings
  use iso_c_binding, only: c_int
  implicit none

  private

  public :: erase_bind
  public :: set_bind
  public :: set_default_keybindings
  public :: is_bind_event
  public :: get_bind_keyname
  public :: is_bind_mousescroll

  ! Public list of binds
  integer, parameter, public :: BIND_QUIT = 1 ! quit the program
  integer, parameter, public :: BIND_NEW = 2 ! create new systems
  integer, parameter, public :: BIND_OPEN = 3 ! open systems
  integer, parameter, public :: BIND_CLOSE_FOCUSED_DIALOG = 4 ! close focused dialog
  integer, parameter, public :: BIND_OK_FOCUSED_DIALOG = 5 ! OK focused dialog
  integer, parameter, public :: BIND_TREE_REMOVE_SYSTEM_FIELD = 6 ! tree: remove system or field
  integer, parameter, public :: BIND_TREE_MOVE_UP = 7 ! tree: move selection up
  integer, parameter, public :: BIND_TREE_MOVE_DOWN = 8 ! tree: move selection down
  integer, parameter, public :: BIND_INPCON_RUN = 9 ! tree: remove system
  integer, parameter, public :: BIND_VIEW_INC_NCELL = 10 ! view: increase number of cells
  integer, parameter, public :: BIND_VIEW_DEC_NCELL = 11 ! view: decrease number of cells
  integer, parameter, public :: BIND_NAV_ROTATE = 12 ! view: rotate the view
  integer, parameter, public :: BIND_NAV_TRANSLATE = 13 ! view: translate the view
  integer, parameter, public :: BIND_NAV_ZOOM = 14 ! view: zoom the view
  integer, parameter, public :: BIND_NAV_RESET = 15 ! view: reset the view
  integer, parameter, public :: BIND_NUM = 15 ! total number of binds
  ! #define BIND_CLOSE_ALL_DIALOGS 2  // Closes all windows
  ! #define BIND_VIEW_ALIGN_A_AXIS 3  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_B_AXIS 4  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_C_AXIS 5  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_X_AXIS 6  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_Y_AXIS 7  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_Z_AXIS 8  // Align view with a axis

  ! module procedure interfaces
  interface
     module subroutine erase_bind(key, mod, group)
       integer(c_int), intent(in) :: key, mod
       integer, intent(in) :: group
     end subroutine erase_bind
     module subroutine set_bind(bind, key, mod)
       use tools_io, only: ferror, faterr
       integer, intent(in) :: bind
       integer(c_int), intent(in) :: key, mod
     end subroutine set_bind
     module subroutine set_default_keybindings()
     end subroutine set_default_keybindings
     module function is_bind_event(bind,held)
       integer, intent(in) :: bind
       logical, intent(in), optional :: held
       logical :: is_bind_event
     end function is_bind_event
     module function get_bind_keyname(bind)
       integer, intent(in) :: bind
       character(len=:), allocatable :: get_bind_keyname
     end function get_bind_keyname
     module function is_bind_mousescroll(bind)
       integer, intent(in) :: bind
       logical :: is_bind_mousescroll
     end function is_bind_mousescroll
  end interface

end module keybindings
