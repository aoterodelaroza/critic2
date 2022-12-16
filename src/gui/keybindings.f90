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
  use interfaces_cimgui, only: ImGuiKey_COUNT
  use iso_c_binding, only: c_int
  implicit none

  private

  public :: erase_bind
  public :: set_bind
  public :: set_bind_from_user_input
  public :: set_default_keybindings
  public :: is_bind_event
  public :: get_bind_keyname
  public :: is_bind_mousescroll

  ! global flag for keybinding use
  logical, public :: use_keybindings = .true.

  ! mouse keybindings
  integer, parameter, public :: ImGuiKey_MouseLeft = ImGuiKey_COUNT + 1
  integer, parameter, public :: ImGuiKey_MouseLeftDouble = ImGuiKey_COUNT + 2
  integer, parameter, public :: ImGuiKey_MouseRight = ImGuiKey_COUNT + 3
  integer, parameter, public :: ImGuiKey_MouseRightDouble = ImGuiKey_COUNT + 4
  integer, parameter, public :: ImGuiKey_MouseMiddle = ImGuiKey_COUNT + 5
  integer, parameter, public :: ImGuiKey_MouseMiddleDouble = ImGuiKey_COUNT + 6
  integer, parameter, public :: ImGuiKey_MouseScroll = ImGuiKey_COUNT + 11

  ! Public list of binds
  integer, parameter, public :: BIND_QUIT = 1 ! quit the program
  integer, parameter, public :: BIND_NEW = 2 ! create new systems
  integer, parameter, public :: BIND_OPEN = 3 ! open systems
  integer, parameter, public :: BIND_CLOSE_ALL_DIALOGS = 4 ! close all open dialogs
  integer, parameter, public :: BIND_CLOSE_FOCUSED_DIALOG = 5 ! close focused dialog
  integer, parameter, public :: BIND_OK_FOCUSED_DIALOG = 6 ! OK focused dialog
  integer, parameter, public :: BIND_TREE_REMOVE_SYSTEM_FIELD = 7 ! tree: remove system or field
  integer, parameter, public :: BIND_TREE_MOVE_UP = 8 ! tree: move selection up
  integer, parameter, public :: BIND_TREE_MOVE_DOWN = 9 ! tree: move selection down
  integer, parameter, public :: BIND_INPCON_RUN = 10 ! tree: remove system
  integer, parameter, public :: BIND_VIEW_INC_NCELL = 11 ! view: increase number of cells
  integer, parameter, public :: BIND_VIEW_DEC_NCELL = 12 ! view: decrease number of cells
  integer, parameter, public :: BIND_VIEW_ALIGN_A_AXIS = 13 ! view: align view with a axis
  integer, parameter, public :: BIND_VIEW_ALIGN_B_AXIS = 14 ! view: align view with b axis
  integer, parameter, public :: BIND_VIEW_ALIGN_C_AXIS = 15 ! view: align view with c axis
  integer, parameter, public :: BIND_VIEW_ALIGN_X_AXIS = 16 ! view: align view with x axis
  integer, parameter, public :: BIND_VIEW_ALIGN_Y_AXIS = 17 ! view: align view with y axis
  integer, parameter, public :: BIND_VIEW_ALIGN_Z_AXIS = 18 ! view: align view with z axis
  integer, parameter, public :: BIND_NAV_ROTATE = 19 ! view: rotate the view
  integer, parameter, public :: BIND_NAV_TRANSLATE = 20 ! view: translate the view
  integer, parameter, public :: BIND_NAV_ZOOM = 21 ! view: zoom the view
  integer, parameter, public :: BIND_NAV_RESET = 22 ! view: reset the view
  integer, parameter, public :: BIND_NUM = 22 ! total number of binds

  ! Bind names
  character(len=31), parameter, public :: bindnames(BIND_NUM) = (/&
     "Quit                           ",& ! BIND_QUIT
     "New                            ",& ! BIND_NEW
     "Open file(s)                   ",& ! BIND_OPEN
     "Close all dialogs              ",& ! BIND_CLOSE_ALL_DIALOGS
     "Close focused dialog           ",& ! BIND_CLOSE_FOCUSED_DIALOG
     "OK in focused dialog           ",& ! BIND_OK_FOCUSED_DIALOG
     "Remove selected system or field",& ! BIND_TREE_REMOVE_SYSTEM_FIELD
     "Select previous system in tree ",& ! BIND_TREE_MOVE_UP
     "Select next system in tree     ",& ! BIND_TREE_MOVE_DOWN
     "Run the commands               ",& ! BIND_INPCON_RUN
     "Increase number of cells       ",& ! BIND_VIEW_INC_NCELL
     "Decrease number of cells       ",& ! BIND_VIEW_DEC_NCELL
     "Align with a axis              ",& ! BIND_VIEW_ALIGN_A_AXIS
     "Align with b axis              ",& ! BIND_VIEW_ALIGN_B_AXIS
     "Align with c axis              ",& ! BIND_VIEW_ALIGN_C_AXIS
     "Align with x axis              ",& ! BIND_VIEW_ALIGN_X_AXIS
     "Align with y axis              ",& ! BIND_VIEW_ALIGN_Y_AXIS
     "Align with z axis              ",& ! BIND_VIEW_ALIGN_Z_AXIS
     "Rotate the camera              ",& ! BIND_NAV_ROTATE
     "Translate the camera           ",& ! BIND_NAV_TRANSLATE
     "Camera zoom                    ",& ! BIND_NAV_ZOOM
     "Reset the camera               "/) ! BIND_NAV_RESET

  ! The key associated with each bind, bind -> key
  integer(c_int), public :: keybind(BIND_NUM)

  ! The modifiers associated with each bind, bind -> mod
  integer(c_int), parameter :: mod_none = 0
  integer(c_int), parameter :: mod_ctrl = 1
  integer(c_int), parameter :: mod_alt = 2
  integer(c_int), parameter :: mod_shift = 4
  integer(c_int), parameter :: mod_super = 8
  integer(c_int), public :: modbind(BIND_NUM)

  ! module procedure interfaces
  interface
     module subroutine erase_bind(key, mod, group)
       integer(c_int), intent(in) :: key, mod
       integer, intent(in) :: group
     end subroutine erase_bind
     module subroutine set_bind(bind, key, mod)
       integer, intent(in) :: bind
       integer(c_int), intent(in) :: key, mod
     end subroutine set_bind
     module function set_bind_from_user_input(bind)
       integer, intent(in) :: bind
       logical :: set_bind_from_user_input
     end function set_bind_from_user_input
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
