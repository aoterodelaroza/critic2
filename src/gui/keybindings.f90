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
  public :: is_mod_key

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
  integer, parameter, public :: BIND_CLOSE = 4 ! close system
  integer, parameter, public :: BIND_REOPEN = 5 ! reopen system
  integer, parameter, public :: BIND_GEOMETRY = 6 ! view/edit geometry
  integer, parameter, public :: BIND_SAVE = 7 ! save
  integer, parameter, public :: BIND_EXPORT_NOW = 8 ! export now
  integer, parameter, public :: BIND_CLOSE_ALL_DIALOGS = 9 ! close all open dialogs
  integer, parameter, public :: BIND_CLOSE_FOCUSED_DIALOG = 10 ! close focused dialog
  integer, parameter, public :: BIND_OK_FOCUSED_DIALOG = 11 ! OK focused dialog
  integer, parameter, public :: BIND_VIEWMODE_SELECT = 12 ! tree: remove system or field
  integer, parameter, public :: BIND_TREE_REMOVE_SYSTEM_FIELD = 13 ! tree: remove system or field
  integer, parameter, public :: BIND_TREE_MOVE_UP = 14 ! tree: move selection up
  integer, parameter, public :: BIND_TREE_MOVE_DOWN = 15 ! tree: move selection down
  integer, parameter, public :: BIND_INPCON_RUN = 16 ! tree: remove system
  integer, parameter, public :: BIND_VIEW_INC_NCELL = 17 ! view: increase number of cells
  integer, parameter, public :: BIND_VIEW_DEC_NCELL = 18 ! view: decrease number of cells
  integer, parameter, public :: BIND_VIEW_ALIGN_A_AXIS = 19 ! view: align view with a axis
  integer, parameter, public :: BIND_VIEW_ALIGN_B_AXIS = 20 ! view: align view with b axis
  integer, parameter, public :: BIND_VIEW_ALIGN_C_AXIS = 21 ! view: align view with c axis
  integer, parameter, public :: BIND_VIEW_ALIGN_X_AXIS = 22 ! view: align view with x axis
  integer, parameter, public :: BIND_VIEW_ALIGN_Y_AXIS = 23 ! view: align view with y axis
  integer, parameter, public :: BIND_VIEW_ALIGN_Z_AXIS = 24 ! view: align view with z axis
  integer, parameter, public :: BIND_VIEW_TOGGLE_ATOMS = 25 ! view: toggle atoms in the first rep
  integer, parameter, public :: BIND_VIEW_TOGGLE_BONDS = 26 ! view: toggle bonds in the first rep
  integer, parameter, public :: BIND_VIEW_CYCLE_LABELS = 27 ! view: cycle labels in the first rep
  integer, parameter, public :: BIND_VIEW_TOGGLE_CELL = 28 ! view: toggle cell in the first rep
  integer, parameter, public :: BIND_NAV_ROTATE = 29 ! navigation: rotate the view
  integer, parameter, public :: BIND_NAV_ROTATE_PERP = 30 ! navigation: rotate around axis perp. to screen
  integer, parameter, public :: BIND_NAV_TRANSLATE = 31 ! navigation: translate the view
  integer, parameter, public :: BIND_NAV_ZOOM = 32 ! navigation: zoom the view
  integer, parameter, public :: BIND_NAV_RESET = 33 ! navigation: reset the view
  integer, parameter, public :: BIND_NAV_MEASURE = 34 ! navigation: reset the view
  integer, parameter, public :: BIND_SELECT_ATOMS = 35 ! select: select single atoms
  integer, parameter, public :: BIND_SELECT_MOLECULES = 36 ! select: select molecules
  integer, parameter, public :: BIND_SELECT_MOLECULES_AND_DESELECT = 37 ! select: select molecules/deselect
  integer, parameter, public :: BIND_EDITGEOM_REMOVE = 38 ! edit geometry: remove atoms
  integer, parameter, public :: BIND_EDITGEOM_DESELECT = 39 ! edit geometry: deselect all atoms
  integer, parameter, public :: BIND_NUM = 39 ! total number of binds

  ! Bind names
  character(len=32), parameter, public :: bindnames(BIND_NUM) = (/&
     "Quit                            ",& ! BIND_QUIT
     "New                             ",& ! BIND_NEW
     "Open file(s)                    ",& ! BIND_OPEN
     "Close current system            ",& ! BIND_CLOSE
     "Reopen current system           ",& ! BIND_REOPEN
     "View/Edit Geometry              ",& ! BIND_GEOMETRY
     "Save                            ",& ! BIND_SAVE
     "Export to PNG                   ",& ! BIND_EXPORT_NOW
     "Close all dialogs               ",& ! BIND_CLOSE_ALL_DIALOGS
     "Close focused dialog            ",& ! BIND_CLOSE_FOCUSED_DIALOG
     "OK in focused dialog            ",& ! BIND_OK_FOCUSED_DIALOG
     "Select                          ",& ! BIND_VIEWMODE_SELECT
     "Remove selected system or field ",& ! BIND_TREE_REMOVE_SYSTEM_FIELD
     "Select previous system in tree  ",& ! BIND_TREE_MOVE_UP
     "Select next system in tree      ",& ! BIND_TREE_MOVE_DOWN
     "Run input commands              ",& ! BIND_INPCON_RUN
     "Increase number of cells shown  ",& ! BIND_VIEW_INC_NCELL
     "Decrease number of cells shown  ",& ! BIND_VIEW_DEC_NCELL
     "Align with a axis               ",& ! BIND_VIEW_ALIGN_A_AXIS
     "Align with b axis               ",& ! BIND_VIEW_ALIGN_B_AXIS
     "Align with c axis               ",& ! BIND_VIEW_ALIGN_C_AXIS
     "Align with x axis               ",& ! BIND_VIEW_ALIGN_X_AXIS
     "Align with y axis               ",& ! BIND_VIEW_ALIGN_Y_AXIS
     "Align with z axis               ",& ! BIND_VIEW_ALIGN_Z_AXIS
     "Toggle display of atoms         ",& ! BIND_VIEW_TOGGLE_ATOMS
     "Toggle display of bonds         ",& ! BIND_VIEW_TOGGLE_BONDS
     "Cycle through display of labels ",& ! BIND_VIEW_CYCLE_LABELS
     "Toggle display of unit cell     ",& ! BIND_VIEW_TOGGLE_CELL
     "Rotate the camera               ",& ! BIND_NAV_ROTATE
     "Rotate around perpendicular axis",& ! BIND_NAV_ROTATE_PERP
     "Translate the camera            ",& ! BIND_NAV_TRANSLATE
     "Camera zoom                     ",& ! BIND_NAV_ZOOM
     "Reset the camera                ",& ! BIND_NAV_RESET
     "Measure distances and angles    ",& ! BIND_NAV_MEASURE
     "Select atoms                    ",& ! BIND_SELECT_ATOMS
     "Select molecules                ",& ! BIND_SELECT_MOLECULES
     "Select molecules/Deselect       ",& ! BIND_SELECT_MOLECULES_AND_DESELECT
     "Remove selected atoms           ",& ! BIND_EDITGEOM_REMOVE
     "Deselect all atoms              "&  ! BIND_EDITGEOM_DESELECT
     /)

  ! The key associated with each bind, bind -> key
  integer(c_int), public :: keybind(BIND_NUM)

  ! The modifiers associated with each bind, bind -> mod
  integer(c_int), parameter :: mod_none = 0
  integer(c_int), parameter :: mod_shift = 1
  integer(c_int), parameter :: mod_ctrl = 2
  integer(c_int), parameter :: mod_alt = 4
  integer(c_int), parameter :: mod_super = 8
  integer(c_int), public :: modbind(BIND_NUM)

  ! The keybinding groups. The first group (1) must be the global.
  integer, parameter, public :: group_global = 1              ! keybindings that apply everywhere
  integer, parameter, public :: group_tree = 2                ! if the tree is active
  integer, parameter, public :: group_inpcon = 3              ! the input console is active
  integer, parameter, public :: group_dialog = 4              ! a dialog is active
  integer, parameter, public :: group_view = 5                ! if the view is active
  integer, parameter, public :: group_viewmode = 6            ! view mouse interaction modes
  integer, parameter, public :: group_viewmode_navigation = 7 ! view mouse interaction modes
  integer, parameter, public :: group_viewmode_select = 8     ! view mouse interaction modes
  integer, parameter, public :: group_editgeom = 9            ! if the edit geometry window is active
  integer, parameter, public :: group_NUM = 9                 ! total number of groups

  ! Names of the keybinding groups
  character(len=28), parameter, public :: groupnames(group_NUM) = (/&
     "Global                      ",&
     "Tree Window                 ",&
     "Input Window                ",&
     "Dialogs                     ",&
     "View Window                 ",&
     "View Mouse Interaction Modes",&
     "View Window: Navigation Mode",&
     "View Window: Select Mode    ",&
     "View/Edit Geometry          "/)

  ! Bind groups assignment
  integer, parameter, public :: groupbind(BIND_NUM) = (/&
     group_global,&              ! BIND_QUIT
     group_global,&              ! BIND_NEW
     group_global,&              ! BIND_OPEN
     group_global,&              ! BIND_CLOSE
     group_global,&              ! BIND_REOPEN
     group_global,&              ! BIND_GEOMETRY
     group_global,&              ! BIND_SAVE
     group_global,&              ! BIND_EXPORT_NOW
     group_global,&              ! BIND_CLOSE_ALL_DIALOGS
     group_dialog,&              ! BIND_CLOSE_FOCUSED_DIALOG
     group_dialog,&              ! BIND_OK_FOCUSED_DIALOG
     group_viewmode,&            ! BIND_VIEWMODE_SELECT
     group_tree,&                ! BIND_TREE_REMOVE_SYSTEM_FIELD
     group_tree,&                ! BIND_TREE_MOVE_UP
     group_tree,&                ! BIND_TREE_MOVE_DOWN
     group_inpcon,&              ! BIND_INPCON_RUN
     group_view,&                ! BIND_VIEW_INC_NCELL
     group_view,&                ! BIND_VIEW_DEC_NCELL
     group_view,&                ! BIND_VIEW_ALIGN_A_AXIS
     group_view,&                ! BIND_VIEW_ALIGN_B_AXIS
     group_view,&                ! BIND_VIEW_ALIGN_C_AXIS
     group_view,&                ! BIND_VIEW_ALIGN_X_AXIS
     group_view,&                ! BIND_VIEW_ALIGN_Y_AXIS
     group_view,&                ! BIND_VIEW_ALIGN_Z_AXIS
     group_view,&                ! BIND_VIEW_TOGGLE_ATOMS
     group_view,&                ! BIND_VIEW_TOGGLE_BONDS
     group_view,&                ! BIND_VIEW_CYCLE_LABELS
     group_view,&                ! BIND_VIEW_TOGGLE_CELL
     group_viewmode_navigation,& ! BIND_NAV_ROTATE
     group_viewmode_navigation,& ! BIND_NAV_ROTATE_PERP
     group_viewmode_navigation,& ! BIND_NAV_TRANSLATE
     group_viewmode_navigation,& ! BIND_NAV_ZOOM
     group_viewmode_navigation,& ! BIND_NAV_RESET
     group_viewmode_navigation,& ! BIND_NAV_MEASURE
     group_viewmode_select,&     ! BIND_SELECT_ATOMS
     group_viewmode_select,&     ! BIND_SELECT_MOLECULES
     group_viewmode_select,&     ! BIND_SELECT_MOLECULES_AND_DESELECT
     group_editgeom,&            ! BIND_EDITGEOM_REMOVE
     group_editgeom/)            ! BIND_EDITGEOM_DESELECT

  ! bindfull -> bindtype
  ! Binding type. If 0, requires pressing a key (not just a modifier)
  ! to trigger. If -1, requires only a modifier. If > 1, inherit
  ! the same modifier as the referenced key binding.
  integer, parameter, public :: bindtype(BIND_NUM) = (/&
     0,&  ! BIND_QUIT
     0,&  ! BIND_NEW
     0,&  ! BIND_OPEN
     0,&  ! BIND_CLOSE
     0,&  ! BIND_REOPEN
     0,&  ! BIND_GEOMETRY
     0,&  ! BIND_SAVE
     0,&  ! BIND_EXPORT_NOW
     0,&  ! BIND_CLOSE_ALL_DIALOGS
     0,&  ! BIND_CLOSE_FOCUSED_DIALOG
     0,&  ! BIND_OK_FOCUSED_DIALOG
     -1,& ! BOND_VIEWMODE_SELECT
     0,&  ! BIND_TREE_REMOVE_SYSTEM_FIELD
     0,&  ! BIND_TREE_MOVE_UP
     0,&  ! BIND_TREE_MOVE_DOWN
     0,&  ! BIND_INPCON_RUN
     0,&  ! BIND_VIEW_INC_NCELL
     0,&  ! BIND_VIEW_DEC_NCELL
     0,&  ! BIND_VIEW_ALIGN_A_AXIS
     0,&  ! BIND_VIEW_ALIGN_B_AXIS
     0,&  ! BIND_VIEW_ALIGN_C_AXIS
     0,&  ! BIND_VIEW_ALIGN_X_AXIS
     0,&  ! BIND_VIEW_ALIGN_Y_AXIS
     0,&  ! BIND_VIEW_ALIGN_Z_AXIS
     0,&  ! BIND_VIEW_TOGGLE_ATOMS
     0,&  ! BIND_VIEW_TOGGLE_BONDS
     0,&  ! BIND_VIEW_CYCLE_LABELS
     0,&  ! BIND_VIEW_TOGGLE_CELL
     0,&  ! BIND_NAV_ROTATE
     0,&  ! BIND_NAV_ROTATE_PERP
     0,&  ! BIND_NAV_TRANSLATE
     0,&  ! BIND_NAV_ZOOM
     0,&  ! BIND_NAV_RESET
     0,&  ! BIND_NAV_MEASURE
     BIND_VIEWMODE_SELECT,&  ! BIND_SELECT_ATOMS
     BIND_VIEWMODE_SELECT,&  ! BIND_SELECT_MOLECULES
     BIND_VIEWMODE_SELECT,&  ! BIND_SELECT_AND_DESELECT
     0,&  ! BIND_EDITGEOM_REMOVE
     0/)  ! BIND_EDITGEOM_DESELECT

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
       character(len=128) :: get_bind_keyname
     end function get_bind_keyname
     module function is_bind_mousescroll(bind)
       integer, intent(in) :: bind
       logical :: is_bind_mousescroll
     end function is_bind_mousescroll
     module function is_mod_key(key)
       integer, intent(in) :: key
       logical :: is_mod_key
     end function is_mod_key
  end interface

end module keybindings
