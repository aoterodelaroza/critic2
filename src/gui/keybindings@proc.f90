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
  use interfaces_cimgui, only: ImGuiKey_COUNT
  use iso_c_binding
  use hashmod, only: hash
  implicit none

  ! mouse keybindings
  integer, parameter :: ImGuiKey_MouseLeft = ImGuiKey_COUNT + 1
  integer, parameter :: ImGuiKey_MouseLeftDouble = ImGuiKey_COUNT + 2
  integer, parameter :: ImGuiKey_MouseRight = ImGuiKey_COUNT + 3
  integer, parameter :: ImGuiKey_MouseRightDouble = ImGuiKey_COUNT + 4
  integer, parameter :: ImGuiKey_MouseMiddle = ImGuiKey_COUNT + 5
  integer, parameter :: ImGuiKey_MouseMiddleDouble = ImGuiKey_COUNT + 6
  integer, parameter :: ImGuiKey_MouseScroll = ImGuiKey_COUNT + 11

  ! Processing level for bind events. Right now 0 = all and 1 = none.
  ! Perhaps more will be added in the future.
  integer, parameter :: bindevent_level = 0

  ! Bind names
  character(len=31), parameter :: bindnames(BIND_NUM) = (/&
     "Quit                           ",& ! BIND_QUIT
     "New                            ",& ! BIND_NEW
     "Open file(s)                   ",& ! BIND_OPEN
     "Close focused dialog           ",& ! BIND_CLOSE_FOCUSED_DIALOG
     "OK in focused dialog           ",& ! BIND_OK_FOCUSED_DIALOG
     "Remove selected system or field",& ! BIND_TREE_REMOVE_SYSTEM_FIELD
     "Run the commands               ",& ! BIND_INPCON_RUN
     "Rotate the camera              ",& ! BIND_NAV_ROTATE
     "Translate the camera           ",& ! BIND_NAV_TRANSLATE
     "Camera zoom                    ",& ! BIND_NAV_ZOOM
     "Reset the camera               "/) ! BIND_NAV_RESET
  !   "Close all dialogs",
  !   "Align view with a axis",
  !   "Align view with b axis",
  !   "Align view with c axis",
  !   "Align view with x axis",
  !   "Align view with y axis",
  !   "Align view with z axis",

  ! Bind groups. The first group (1) must be the global.
  integer, parameter :: group_global = 1
  integer, parameter :: group_tree = 2   ! if the tree is active
  integer, parameter :: group_inpcon = 3 ! the input console is active
  integer, parameter :: group_dialog = 4 ! a dialog is active
  integer, parameter :: group_view = 5   ! if the view is active
  integer, parameter :: groupbind(BIND_NUM) = (/&
     group_global,& ! BIND_QUIT
     group_global,& ! BIND_NEW
     group_global,& ! BIND_OPEN
     group_dialog,& ! BIND_CLOSE_FOCUSED_DIALOG
     group_dialog,& ! BIND_OK_FOCUSED_DIALOG
     group_tree,&   ! BIND_TREE_REMOVE_SYSTEM_FIELD
     group_inpcon,& ! BIND_INPCON_RUN
     group_view,&   ! BIND_NAV_ROTATE
     group_view,&   ! BIND_NAV_TRANSLATE
     group_view,&   ! BIND_NAV_ZOOM
     group_view/)   ! BIND_NAV_RESET
  !   1, // close all dialogs
  !   2, // align view with a axis
  !   2, // align view with b axis
  !   2, // align view with c axis
  !   2, // align view with x axis
  !   2, // align view with y axis
  !   2, // align view with z axis

  integer, parameter :: ngroupbinds = 2

  ! The key associated with each bind, bind -> key
  integer(c_int) :: keybind(BIND_NUM)

  ! The modifiers associated with each bind, bind -> mod
  integer(c_int) :: modbind(BIND_NUM)

  ! Bind for a key, mod, and group combination
  type(hash) :: keymap

  !xx! private procedures
  ! function hkey(key,mod,group)

contains

  ! unbind the key/mod/group combination
  subroutine erase_bind(key, mod, group)
    use interfaces_cimgui, only: ImGuiKey_None
    integer(c_int), intent(in) :: key, mod
    integer, intent(in) :: group

    character(len=:), allocatable :: hk
    integer :: oldbind

    hk = hkey(key,mod,group)
    if (keymap%iskey(hk)) then
       oldbind = keymap%get(hk,oldbind)
       modbind(oldbind) = ImGuiKey_None
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

  module subroutine set_default_keybindings()
    use interfaces_cimgui
    integer :: i

    ! initialize to no keys and modifiers
    do i = 1, BIND_NUM
       keybind(i) = ImGuiKey_None
       modbind(i) = ImGuiKey_None
    end do
    call keymap%init()

    ! Default keybindings
    call set_bind(BIND_QUIT,ImGuiKey_Q,ImGuiKey_ModCtrl)
    call set_bind(BIND_NEW,ImGuiKey_N,ImGuiKey_ModCtrl)
    call set_bind(BIND_OPEN,ImGuiKey_O,ImGuiKey_ModCtrl)
    call set_bind(BIND_CLOSE_FOCUSED_DIALOG,ImGuiKey_Escape,ImGuiKey_None)
    call set_bind(BIND_OK_FOCUSED_DIALOG,ImGuiKey_Enter,ImGuiKey_ModCtrl)
    call set_bind(BIND_TREE_REMOVE_SYSTEM_FIELD,ImGuiKey_Delete,ImGuiKey_None)
    call set_bind(BIND_INPCON_RUN,ImGuiKey_Enter,ImGuiKey_ModCtrl)
    call set_bind(BIND_NAV_ROTATE,ImGuiKey_MouseLeft,ImGuiKey_None)
    call set_bind(BIND_NAV_TRANSLATE,ImGuiKey_MouseRight,ImGuiKey_None)
    call set_bind(BIND_NAV_ZOOM,ImGuiKey_MouseScroll,ImGuiKey_None)
    call set_bind(BIND_NAV_RESET,ImGuiKey_MouseLeftDouble,ImGuiKey_None)

    !   set_bind(BIND_CLOSE_ALL_DIALOGS,GLFW_KEY_DELETE,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_A_AXIS,GLFW_KEY_A,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_B_AXIS,GLFW_KEY_B,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_C_AXIS,GLFW_KEY_C,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_X_AXIS,GLFW_KEY_X,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_Y_AXIS,GLFW_KEY_Y,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_Z_AXIS,GLFW_KEY_Z,NOMOD);

  end subroutine set_default_keybindings

  ! Return whether the bind event is happening. If held (optional),
  ! the event happens only if the button is held down (for mouse).
  module function is_bind_event(bind,held)
    use gui_main, only: io
    use interfaces_cimgui
    integer, intent(in) :: bind
    logical, intent(in), optional :: held
    logical :: is_bind_event

    integer :: key, mod
    logical :: held_

    ! process options
    held_ = .false.
    if (present(held)) held_ = held

    ! some checks
    is_bind_event = .false.
    if (bindevent_level > 0) return
    if (bind < 1 .or. bind > BIND_NUM) return

    ! get current key and mod for this bind
    key = keybind(bind)
    mod = modbind(bind)

    if (key == ImGuiKey_None .or.(mod /= ImGuiKey_None.and..not.igIsKeyDown(mod))) then
       ! no key or the mod is not down
       return
    elseif (key >= ImGuiKey_NamedKey_BEGIN .and. key < ImGuiKey_NamedKey_END) then
       ! .and..not.io%WantCaptureKeyboard .and..not.io%WantTextInput ! this prevents using ESC with textinput
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
    use interfaces_cimgui, only: ImGuiKey_None, igGetKeyName, &
       ImGuiKey_ModCtrl, ImGuiKey_ModShift, ImGuiKey_ModAlt, ImGuiKey_ModSuper
    use c_interface_module, only: C_F_string_ptr_alloc
    use tools_io, only: lower
    integer, intent(in) :: bind
    character(len=:), allocatable :: get_bind_keyname

    integer :: group
    integer(c_int) :: key, mod
    type(c_ptr) :: name
    character(len=:), allocatable :: aux

    get_bind_keyname = ""
    group = groupbind(bind)
    key = keybind(bind)
    mod = modbind(bind)

    if (key /= ImGuiKey_None) then
       if (mod == ImGuiKey_ModCtrl) then
          get_bind_keyname = "Ctrl+"
       elseif (mod == ImGuiKey_ModShift) then
          get_bind_keyname = "Shift+"
       elseif (mod == ImGuiKey_ModAlt) then
          get_bind_keyname = "Alt+"
       elseif (mod == ImGuiKey_ModSuper) then
          get_bind_keyname = "Super+"
       end if
       name = igGetKeyName(key)
       call C_F_string_ptr_alloc(name,aux)
       get_bind_keyname = get_bind_keyname // lower(trim(aux))
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

end submodule proc
