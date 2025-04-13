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
submodule (keybindings) proc
  use iso_c_binding
  use hashmod, only: hash
  implicit none

  ! Processing level for bind events. Right now 0 = all and 1 = none.
  ! Perhaps more will be added in the future.
  integer, parameter :: bindevent_level = 0

  ! Bind groups. The first group (1) must be the global.
  integer, parameter :: group_global = 1
  integer, parameter :: group_tree = 2   ! if the tree is active
  integer, parameter :: group_inpcon = 3 ! the input console is active
  integer, parameter :: group_dialog = 4 ! a dialog is active
  integer, parameter :: group_view = 5   ! if the view is active
  integer, parameter :: groupbind(BIND_NUM) = (/&
     group_global,& ! BIND_QUIT
     group_global,& ! BIND_NEW
     group_global,& ! BIND_GEOMETRY
     group_global,& ! BIND_OPEN
     group_global,& ! BIND_CLOSE
     group_global,& ! BIND_REOPEN
     group_global,& ! BIND_CLOSE_ALL_DIALOGS
     group_dialog,& ! BIND_CLOSE_FOCUSED_DIALOG
     group_dialog,& ! BIND_OK_FOCUSED_DIALOG
     group_tree,&   ! BIND_TREE_REMOVE_SYSTEM_FIELD
     group_tree,&   ! BIND_TREE_MOVE_UP
     group_tree,&   ! BIND_TREE_MOVE_DOWN
     group_inpcon,& ! BIND_INPCON_RUN
     group_view,&   ! BIND_VIEW_INC_NCELL
     group_view,&   ! BIND_VIEW_DEC_NCELL
     group_view,&   ! BIND_VIEW_ALIGN_A_AXIS
     group_view,&   ! BIND_VIEW_ALIGN_B_AXIS
     group_view,&   ! BIND_VIEW_ALIGN_C_AXIS
     group_view,&   ! BIND_VIEW_ALIGN_X_AXIS
     group_view,&   ! BIND_VIEW_ALIGN_Y_AXIS
     group_view,&   ! BIND_VIEW_ALIGN_Z_AXIS
     group_view,&   ! BIND_VIEW_TOGGLE_ATOMS
     group_view,&   ! BIND_VIEW_TOGGLE_BONDS
     group_view,&   ! BIND_VIEW_CYCLE_LABELS
     group_view,&   ! BIND_VIEW_TOGGLE_CELL
     group_view,&   ! BIND_NAV_ROTATE
     group_view,&   ! BIND_NAV_ROTATE_PERP
     group_view,&   ! BIND_NAV_TRANSLATE
     group_view,&   ! BIND_NAV_ZOOM
     group_view,&   ! BIND_NAV_RESET
     group_view/)   ! BIND_NAV_MEASURE

  integer, parameter :: ngroupbinds = 2

  ! Bind for a key, mod, and group combination
  type(hash) :: keymap

  character*18, parameter :: keynames(133) = (/&
     "Tab               ","Left Arrow        ","Right Arrow       ","Up Arrow          ","Down Arrow        ",&
     "Page Up           ","Page Down         ","Home              ","End               ","Insert            ",&
     "Delete            ","Backspace         ","Space             ","Enter             ","Escape            ",&
     "Left Ctrl         ","Left Shift        ","Left Alt          ","Left Super        ","Right Ctrl        ",&
     "Right Shift       ","Right Alt         ","Right Super       ","Menu              ","0                 ",&
     "1                 ","2                 ","3                 ","4                 ","5                 ",&
     "6                 ","7                 ","8                 ","9                 ","A                 ",&
     "B                 ","C                 ","D                 ","E                 ","F                 ",&
     "G                 ","H                 ","I                 ","J                 ","K                 ",&
     "L                 ","M                 ","N                 ","O                 ","P                 ",&
     "Q                 ","R                 ","S                 ","T                 ","U                 ",&
     "V                 ","W                 ","X                 ","Y                 ","Z                 ",&
     "F1                ","F2                ","F3                ","F4                ","F5                ",&
     "F6                ","F7                ","F8                ","F9                ","F10               ",&
     "F11               ","F12               ","'                 ",",                 ","-                 ",&
     ".                 ","/                 ",";                 ","=                 ","[                 ",&
     "\                 ","]                 ","`                 ","Caps Lock         ","Scroll Lock       ",&
     "Num Lock          ","Print Screen      ","Pause             ","Keypad 0          ","Keypad 1          ",&
     "Keypad 2          ","Keypad 3          ","Keypad 4          ","Keypad 5          ","Keypad 6          ",&
     "Keypad 7          ","Keypad 8          ","Keypad 9          ","Keypad .          ","Keypad /          ",&
     "Keypad *          ","Keypad -          ","Keypad +          ","Keypad Enter      ","Keypad =          ",&
     "GamepadStart      ","GamepadBack       ","GamepadFaceUp     ","GamepadFaceDown   ","GamepadFaceLeft   ",&
     "GamepadFaceRight  ","GamepadDpadUp     ","GamepadDpadDown   ","GamepadDpadLeft   ","GamepadDpadRight  ",&
     "GamepadL1         ","GamepadR1         ","GamepadL2         ","GamepadR2         ","GamepadL3         ",&
     "GamepadR3         ","GamepadLStickUp   ","GamepadLStickDown ","GamepadLStickLeft ","GamepadLStickRight",&
     "GamepadRStickUp   ","GamepadRStickDown ","GamepadRStickLeft ","GamepadRStickRight","Mod Ctrl          ",&
     "Mod Shift         ","Mod Alt           ","Mod Super         "/)

  !xx! private procedures
  ! function hkey(key,mod,group)

contains

  ! unbind the key/mod/group combination
  module subroutine erase_bind(key, mod, group)
    use interfaces_cimgui, only: ImGuiKey_None
    integer(c_int), intent(in) :: key, mod
    integer, intent(in) :: group

    character(len=:), allocatable :: hk
    integer :: oldbind

    hk = hkey(key,mod,group)
    if (keymap%iskey(hk)) then
       oldbind = keymap%get(hk,oldbind)
       modbind(oldbind) = mod_none
       keybind(oldbind) = ImGuiKey_None
       call keymap%delkey(hk)
    end if

  end subroutine erase_bind

  ! Associate a bind with a key/mod combination
  module subroutine set_bind(bind, key, mod)
    use interfaces_cimgui, only: ImGuiKey_None
    use tools_io, only: ferror, faterr
    integer, intent(in) :: bind
    integer(c_int), intent(in) :: key, mod

    integer :: group, i
    integer(c_int) :: oldkey, oldmod
    character(len=:), allocatable :: hk

    if (bind < 1 .or. bind > BIND_NUM) &
       call ferror('set_bind','BIND number out of range',faterr)

    group = groupbind(bind)

    ! erase the key+mod combination for this bind from the keymap
    oldkey = keybind(bind)
    oldmod = modbind(bind)
    if (oldkey /= ImGuiKey_None) then
       hk = hkey(oldkey,oldmod,group)
       if (keymap%iskey(hk)) call keymap%delkey(hk)
    end if
    ! unbind the previous owner of this key+mod combination in this group...
    call erase_bind(key,mod,group)
    if (group == group_global) then
       ! unbind in all other groups
       do i = 2, ngroupbinds
          call erase_bind(key,mod,i)
       end do
    else
       ! unbind from the global group
       call erase_bind(key,mod,group_global)
    end if

    ! make the new bind
    keybind(bind) = key
    modbind(bind) = mod
    call keymap%put(hkey(key,mod,group),bind)

  end subroutine set_bind

  ! Read user input and set key binding bind. Returns true if a key
  ! binding has been set.
  module function set_bind_from_user_input(bind)
    use interfaces_cimgui
    use gui_main, only: io
    integer, intent(in) :: bind
    logical :: set_bind_from_user_input

    integer :: i, key, mod

    ! get the current mod
    mod = get_current_mod()

    ! get the key
    key = -1
    do i = ImGuiKey_NamedKey_BEGIN, ImGuiKey_NamedKey_END-1
       if ((i==ImGuiKey_ModCtrl).or.(i==ImGuiKey_ModShift).or.(i==ImGuiKey_ModAlt).or.&
          (i==ImGuiKey_ModSuper).or.(i==ImGuiKey_LeftCtrl).or.(i==ImGuiKey_LeftShift).or.&
          (i==ImGuiKey_LeftAlt).or.(i==ImGuiKey_LeftSuper).or.(i==ImGuiKey_RightCtrl).or.&
          (i==ImGuiKey_RightShift).or.(i==ImGuiKey_RightAlt).or.(i==ImGuiKey_RightSuper)) cycle
       if (igIsKeyDown(i)) then
          key = i
          exit
       end if
    end do
    if (igIsMouseDown(ImGuiMouseButton_Left)) key = ImGuiKey_MouseLeft
    if (igIsMouseDown(ImGuiMouseButton_Right)) key = ImGuiKey_MouseRight
    if (igIsMouseDown(ImGuiMouseButton_Middle)) key = ImGuiKey_MouseMiddle
    if (abs(io%MouseWheel) > 1e-8_c_float) key = ImGuiKey_MouseScroll
    if (key /= -1) then
       set_bind_from_user_input = .true.
       call set_bind(bind,key,mod)
    else
       set_bind_from_user_input = .false.
    end if

  end function set_bind_from_user_input

  module subroutine set_default_keybindings()
    use interfaces_cimgui
    integer :: i

    ! initialize to no keys and modifiers
    do i = 1, BIND_NUM
       keybind(i) = ImGuiKey_None
       modbind(i) = mod_none
    end do
    call keymap%init()

    ! Default keybindings
    call set_bind(BIND_QUIT,ImGuiKey_Q,mod_ctrl)
    call set_bind(BIND_NEW,ImGuiKey_N,mod_ctrl)
    call set_bind(BIND_GEOMETRY,ImGuiKey_G,mod_none)
    call set_bind(BIND_OPEN,ImGuiKey_O,mod_ctrl)
    call set_bind(BIND_CLOSE,ImGuiKey_W,mod_ctrl)
    call set_bind(BIND_REOPEN,ImGuiKey_R,mod_ctrl)
    call set_bind(BIND_CLOSE_ALL_DIALOGS,ImGuiKey_Backspace,mod_none)
    call set_bind(BIND_CLOSE_FOCUSED_DIALOG,ImGuiKey_Escape,mod_none)
    call set_bind(BIND_OK_FOCUSED_DIALOG,ImGuiKey_Enter,mod_ctrl)
    call set_bind(BIND_TREE_REMOVE_SYSTEM_FIELD,ImGuiKey_Delete,mod_none)
    call set_bind(BIND_TREE_MOVE_UP,ImGuiKey_K,mod_none)
    call set_bind(BIND_TREE_MOVE_DOWN,ImGuiKey_J,mod_none)
    call set_bind(BIND_INPCON_RUN,ImGuiKey_Enter,mod_ctrl)
    call set_bind(BIND_VIEW_INC_NCELL,ImGuiKey_KeypadAdd,mod_none)
    call set_bind(BIND_VIEW_DEC_NCELL,ImGuiKey_KeypadSubtract,mod_none)
    call set_bind(BIND_VIEW_ALIGN_A_AXIS,ImGuiKey_A,mod_none)
    call set_bind(BIND_VIEW_ALIGN_B_AXIS,ImGuiKey_B,mod_none)
    call set_bind(BIND_VIEW_ALIGN_C_AXIS,ImGuiKey_C,mod_none)
    call set_bind(BIND_VIEW_ALIGN_X_AXIS,ImGuiKey_X,mod_none)
    call set_bind(BIND_VIEW_ALIGN_Y_AXIS,ImGuiKey_Y,mod_none)
    call set_bind(BIND_VIEW_ALIGN_Z_AXIS,ImGuiKey_Z,mod_none)
    call set_bind(BIND_VIEW_TOGGLE_ATOMS,ImGuiKey_Q,mod_none)
    call set_bind(BIND_VIEW_TOGGLE_BONDS,ImGuiKey_W,mod_none)
    call set_bind(BIND_VIEW_CYCLE_LABELS,ImGuiKey_E,mod_none)
    call set_bind(BIND_VIEW_TOGGLE_CELL,ImGuiKey_R,mod_none)
    call set_bind(BIND_NAV_ROTATE,ImGuiKey_MouseLeft,mod_none)
    call set_bind(BIND_NAV_ROTATE_PERP,ImGuiKey_MouseMiddle,mod_none)
    call set_bind(BIND_NAV_TRANSLATE,ImGuiKey_MouseRight,mod_none)
    call set_bind(BIND_NAV_ZOOM,ImGuiKey_MouseScroll,mod_none)
    call set_bind(BIND_NAV_RESET,ImGuiKey_MouseRightDouble,mod_none)
    call set_bind(BIND_NAV_MEASURE,ImGuiKey_MouseLeftDouble,mod_none)

  end subroutine set_default_keybindings

  ! Return whether the bind event is happening. If held (optional),
  ! the event happens only if the button is held down (for mouse).
  module function is_bind_event(bind,held)
    use gui_main, only: io
    use interfaces_cimgui
    integer, intent(in) :: bind
    logical, intent(in), optional :: held
    logical :: is_bind_event

    integer :: key, mod, modnow
    logical :: held_

    ! process options
    held_ = .false.
    if (present(held)) held_ = held

    ! some checks
    is_bind_event = .false.
    if (bindevent_level > 0) return
    if (bind < 1 .or. bind > BIND_NUM) return
    if (.not.use_keybindings) return

    ! get key and mod for this bind, and the current mod
    key = keybind(bind)
    mod = modbind(bind)
    modnow = get_current_mod()

    ! check if any bind is triggered
    if (key == ImGuiKey_None .or. mod /= modnow) then
       ! no key or the mod is not correct
       return
    elseif (key >= ImGuiKey_NamedKey_BEGIN .and. key < ImGuiKey_NamedKey_END .and..not.io%WantTextInput) then
       ! correct key ID and not keyboard captured or inputing text
       if (held_) then
          is_bind_event = igIsKeyDown(key)
       else
          is_bind_event = igIsKeyPressed(key,.true._c_bool)
       end if
    elseif (key == ImGuiKey_MouseLeft .and. held_) then
       is_bind_event = igIsMouseDown(ImGuiMouseButton_Left)
    elseif (key == ImGuiKey_MouseLeft .and. .not.held_) then
       is_bind_event = igIsMouseClicked(ImGuiMouseButton_Left,.false._c_bool)
    elseif (key == ImGuiKey_MouseRight .and. held_) then
       is_bind_event = igIsMouseDown(ImGuiMouseButton_Right)
    elseif (key == ImGuiKey_MouseRight .and. .not.held_) then
       is_bind_event = igIsMouseClicked(ImGuiMouseButton_Right,.false._c_bool)
    elseif (key == ImGuiKey_MouseMiddle .and. held_) then
       is_bind_event = igIsMouseDown(ImGuiMouseButton_Middle)
    elseif (key == ImGuiKey_MouseMiddle .and. .not.held_) then
       is_bind_event = igIsMouseClicked(ImGuiMouseButton_Middle,.false._c_bool)
    elseif (key == ImGuiKey_MouseLeftDouble .and. .not.held_) then
       is_bind_event = igIsMouseDoubleClicked(ImGuiMouseButton_Left)
    elseif (key == ImGuiKey_MouseRightDouble .and. .not.held_) then
       is_bind_event = igIsMouseDoubleClicked(ImGuiMouseButton_Right)
    elseif (key == ImGuiKey_MouseMiddleDouble .and. .not.held_) then
       is_bind_event = igIsMouseDoubleClicked(ImGuiMouseButton_Middle)
    elseif (key == ImGuiKey_MouseScroll) then
       is_bind_event = (abs(io%MouseWheel) > 1e-8_c_float)
    end if

  end function is_bind_event

  ! Return the key+mod combination for a given bind
  module function get_bind_keyname(bind)
    use interfaces_cimgui, only: ImGuiKey_None, ImGuiKey_NamedKey_BEGIN
    use c_interface_module, only: C_F_string_ptr_alloc
    use tools_io, only: lower
    integer, intent(in) :: bind
    character(len=128) :: get_bind_keyname

    integer :: group
    integer(c_int) :: key, mod

    get_bind_keyname = ""
    group = groupbind(bind)
    key = keybind(bind)
    mod = modbind(bind)

    if (key /= ImGuiKey_None) then
       get_bind_keyname = ""
       if (iand(mod,mod_ctrl)/=0)  get_bind_keyname = trim(get_bind_keyname) // "Ctrl+"
       if (iand(mod,mod_alt)/=0)   get_bind_keyname = trim(get_bind_keyname) // "Alt+"
       if (iand(mod,mod_shift)/=0) get_bind_keyname = trim(get_bind_keyname) // "Shift+"
       if (iand(mod,mod_super)/=0) get_bind_keyname = trim(get_bind_keyname) // "Super+"

       if (key == ImGuiKey_MouseLeft) then
          get_bind_keyname = trim(get_bind_keyname) // "Left Mouse"
       elseif (key == ImGuiKey_MouseLeftDouble) then
          get_bind_keyname = trim(get_bind_keyname) // "Double Left Mouse"
       elseif (key == ImGuiKey_MouseRight) then
          get_bind_keyname = trim(get_bind_keyname) // "Right Mouse"
       elseif (key == ImGuiKey_MouseRightDouble) then
          get_bind_keyname = trim(get_bind_keyname) // "Double Right Mouse"
       elseif (key == ImGuiKey_MouseMiddle) then
          get_bind_keyname = trim(get_bind_keyname) // "Middle Mouse"
       elseif (key == ImGuiKey_MouseMiddleDouble) then
          get_bind_keyname = trim(get_bind_keyname) // "Double Middle Mouse"
       elseif (key == ImGuiKey_MouseScroll) then
          get_bind_keyname = trim(get_bind_keyname) // "Mouse Wheel"
       else
          get_bind_keyname = trim(get_bind_keyname) // trim(keynames(key - ImGuiKey_NamedKey_BEGIN + 1))
       end if
    else
       get_bind_keyname = "<not bound>"
    end if

  end function get_bind_keyname

  ! Returns true if the bind corresponds to a mouse scroll
  module function is_bind_mousescroll(bind)
    integer, intent(in) :: bind
    logical :: is_bind_mousescroll

    is_bind_mousescroll = (keybind(bind) == ImGuiKey_MouseScroll)

  end function is_bind_mousescroll

  !xx! private procedures

  ! return a key for the keymap hash by combining key ID, mod, and group.
  function hkey(key,mod,group)
    use tools_io, only: string
    integer(c_int), intent(in) :: key, mod
    integer :: group
    character(len=:), allocatable :: hkey

    hkey = string(key) // "_" // string(mod) // "_" // string(group)

  end function hkey

  ! get current status of modifier keys
  function get_current_mod() result(mod)
    use interfaces_cimgui
    integer :: mod

    mod = 0
    if (igIsKeyDown(ImGuiKey_ModCtrl)) mod = mod + mod_ctrl
    if (igIsKeyDown(ImGuiKey_ModAlt)) mod = mod + mod_alt
    if (igIsKeyDown(ImGuiKey_ModShift)) mod = mod + mod_shift
    if (igIsKeyDown(ImGuiKey_ModSuper)) mod = mod + mod_super

  end function get_current_mod

end submodule proc
