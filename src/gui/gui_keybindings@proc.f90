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

! This module handles the key-bindings for the critic2 GUI.
submodule (gui_keybindings) proc
  use iso_c_binding
  use hashmod, only: hash
  implicit none

  ! Processing level for bind events. Right now 0 = all and 1 = none.
  ! Perhaps more will be added in the future.
  integer, parameter :: bindevent_level = 0

  ! Bind names
  character(len=15), parameter :: bindnames(BIND_NUM) = (/&
     "Quit           ",&
     "Tree: Move up  ",&
     "Tree: Move down"/)
  !   "Close last dialog",
  !   "Close all dialogs",
  !   "Align view with a axis",
  !   "Align view with b axis",
  !   "Align view with c axis",
  !   "Align view with x axis",
  !   "Align view with y axis",
  !   "Align view with z axis",
  !   "Camera rotate",
  !   "Camera pan",
  !   "Camera zoom",
  !   "Camera reset",

  ! Bind groups. The first group (1) must be the global.
  integer, parameter :: group_global = 1
  integer, parameter :: group_tree = 2 ! if the tree is active
  integer, parameter :: groupbind(BIND_NUM) = (/&
     group_global,& ! quit
     group_tree,& ! tree: move up
     group_tree/) ! tree: move down
  !   1, // close last dialog
  !   1, // close all dialogs
  !   2, // align view with a axis
  !   2, // align view with b axis
  !   2, // align view with c axis
  !   2, // align view with x axis
  !   2, // align view with y axis
  !   2, // align view with z axis
  !   2, // rotate camera (navigation)
  !   2, // pan camera (navigation)
  !   2, // zoom camera (navigation)
  !   2, // reset camera (navigation)
  integer, parameter :: ngroupbinds = 1

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
    use gui_interfaces_cimgui, only: ImGuiKey_None
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
    use gui_interfaces_cimgui, only: ImGuiKey_None
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
          call erase_bind(key,mod,i);
       end do
    else
       ! unbind from the global group
       call erase_bind(key,mod,group_global);
    end if

    ! make the new bind
    keybind(bind) = key
    modbind(bind) = mod
    call keymap%put(hkey(key,mod,group),bind)

  end subroutine set_bind

  module subroutine set_default_keybindings()
    use gui_interfaces_cimgui
    integer :: i

    ! initialize to no keys and modifiers
    do i = 1, BIND_NUM
       keybind(i) = ImGuiKey_None
       modbind(i) = ImGuiKey_None
    end do
    call keymap%init()

    ! Default keybindings
    call set_bind(BIND_QUIT,ImGuiKey_Q,ImGuiKey_ModCtrl)
    call set_bind(BIND_TREE_MOVE_UP,ImGuiKey_UpArrow,ImGuiKey_None)
    call set_bind(BIND_TREE_MOVE_DOWN,ImGuiKey_DownArrow,ImGuiKey_None)
    !   set_bind(BIND_CLOSE_LAST_DIALOG,GLFW_KEY_ESCAPE,NOMOD);
    !   set_bind(BIND_CLOSE_ALL_DIALOGS,GLFW_KEY_DELETE,NOMOD);
    !
    !   set_bind(BIND_VIEW_ALIGN_A_AXIS,GLFW_KEY_A,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_B_AXIS,GLFW_KEY_B,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_C_AXIS,GLFW_KEY_C,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_X_AXIS,GLFW_KEY_X,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_Y_AXIS,GLFW_KEY_Y,NOMOD);
    !   set_bind(BIND_VIEW_ALIGN_Z_AXIS,GLFW_KEY_Z,NOMOD);
    !
    !   // Default mouse bindings
    !   set_bind(BIND_NAV_ROTATE,GLFW_MOUSE_LEFT,NOMOD);
    !   set_bind(BIND_NAV_TRANSLATE,GLFW_MOUSE_RIGHT,NOMOD);
    !   set_bind(BIND_NAV_ZOOM,GLFW_MOUSE_SCROLL,NOMOD);
    !   set_bind(BIND_NAV_RESET,GLFW_MOUSE_LEFT_DOUBLE,NOMOD);
  end subroutine set_default_keybindings

  ! Return whether the bind event is happening. If held (optional),
  ! the event happens only if the button is held down (for mouse).
  module function is_bind_event(bind,held)
    use gui_interfaces_cimgui
    use gui_main, only: io
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
    elseif (key >= ImGuiKey_NamedKey_BEGIN .and. key < ImGuiKey_NamedKey_END .and.&
       .not.io%WantCaptureKeyboard .and..not.io%WantTextInput) then
       ! correct key ID and not keyboard captured or inputing text
       if (held_) then
          is_bind_event = igIsKeyDown(key)
       else
          is_bind_event = igIsKeyPressed(key,.true._c_bool)
       end if
    else
       ! mouse interaction
       !   else{
       !     if (key == GLFW_MOUSE_LEFT)
       !       if (!held)
       !       return IsMouseClicked(0);
       !       else
       !       return IsMouseDown(0);
       !     else if (key == GLFW_MOUSE_RIGHT)
       !       if (!held)
       !       return IsMouseClicked(1);
       !       else
       !       return IsMouseDown(1);
       !     else if (key == GLFW_MOUSE_MIDDLE)
       !       if (!held)
       !       return IsMouseClicked(2);
       !       else
       !       return IsMouseDown(2);
       !     else if (key == GLFW_MOUSE_BUTTON3)
       !       if (!held)
       !       return IsMouseClicked(3);
       !       else
       !       return IsMouseDown(3);
       !     else if (key == GLFW_MOUSE_BUTTON4)
       !       if (!held)
       !       return IsMouseClicked(4);
       !       else
       !       return IsMouseDown(4);
       !     else if (key == GLFW_MOUSE_LEFT_DOUBLE && !held)
       !       return IsMouseDoubleClicked(0);
       !     else if (key == GLFW_MOUSE_RIGHT_DOUBLE && !held)
       !       return IsMouseDoubleClicked(1);
       !     else if (key == GLFW_MOUSE_MIDDLE_DOUBLE && !held)
       !       return IsMouseDoubleClicked(2);
       !     else if (key == GLFW_MOUSE_BUTTON3_DOUBLE && !held)
       !       return IsMouseDoubleClicked(3);
       !     else if (key == GLFW_MOUSE_BUTTON4_DOUBLE && !held)
       !       return IsMouseDoubleClicked(4);
       !     else if (key == GLFW_MOUSE_SCROLL)
       !       return abs(GetCurrentContext()->IO.MouseWheel) > 1e-8;
       !     return false;
       !   }
    end if

  end function is_bind_event

  ! Return the key+mod combination for a given bind
  module function get_bind_keyname(bind)
    use gui_interfaces_cimgui, only: ImGuiKey_None, igGetKeyName, &
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
