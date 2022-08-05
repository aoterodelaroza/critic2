submodule (gui_keybindings) proc
  use iso_c_binding
  use hashmod, only: hash
  implicit none

  ! Processing level for bind events. Right now 0 = all and 1 = none.
  ! Perhaps more will be added in the future.
  integer, parameter :: bindevent_level = 0

  ! Bind names
  character(len=10), parameter :: bindnames(BIND_NUM) = (/&
     "Quit      "/)
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

  ! Bind groups. Group 1 is global.
  integer, parameter :: groupbind(BIND_NUM) = (/&
     1/) ! quit
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

contains

  module subroutine SetBind(bind, key, mod)
    use tools_io, only: ferror, faterr
    integer, intent(in) :: bind
    integer(c_int), intent(in) :: key, mod

    integer :: group
    integer(c_int) :: oldkey, oldmod

    if (bind < 1 .or. bind > BIND_NUM) &
       call ferror('SetBind','BIND number out of range',faterr)

    group = groupbind(bind)

    ! erase the key+mod combination for this bind from the keymap
    oldkey = keybind(bind)
    oldmod = modbind(bind)
    ! if (keymap.find(std::make_tuple(oldkey,oldmod,group)) != keymap.end())
    !   keymap.erase(std::make_tuple(oldkey,oldmod,group));

    ! unbind the previous owner of this key+mod combination in this group...
    !   EraseBind_(key,mod,group);

    !   if (group == 0){
    !     // ...and in all other groups
    !     for (int i = 1; i < ngroupbinds ; i++)
    !       EraseBind_(key,mod,i);
    !   } else {
    !     // ...and in the 0-group
    !     EraseBind_(key,mod,0);
    !   }

    ! make the new bind
    keybind(bind) = key
    modbind(bind) = mod
    ! keymap[std::make_tuple(key,mod,group)] = bind;

  end subroutine SetBind

  module subroutine SetDefaultKeyBindings()
    use gui_interfaces_cimgui

    integer :: i

    ! initialize to no keys and modifiers
    do i = 1, BIND_NUM
       keybind(i) = ImGuiKey_None
       modbind(i) = ImGuiKey_None
    end do

    ! Default keybindings
    call SetBind(BIND_QUIT,ImGuiKey_Q,ImGuiKey_ModCtrl)
!   SetBind(BIND_CLOSE_LAST_DIALOG,GLFW_KEY_ESCAPE,NOMOD);
!   SetBind(BIND_CLOSE_ALL_DIALOGS,GLFW_KEY_DELETE,NOMOD);
!
!   SetBind(BIND_VIEW_ALIGN_A_AXIS,GLFW_KEY_A,NOMOD);
!   SetBind(BIND_VIEW_ALIGN_B_AXIS,GLFW_KEY_B,NOMOD);
!   SetBind(BIND_VIEW_ALIGN_C_AXIS,GLFW_KEY_C,NOMOD);
!   SetBind(BIND_VIEW_ALIGN_X_AXIS,GLFW_KEY_X,NOMOD);
!   SetBind(BIND_VIEW_ALIGN_Y_AXIS,GLFW_KEY_Y,NOMOD);
!   SetBind(BIND_VIEW_ALIGN_Z_AXIS,GLFW_KEY_Z,NOMOD);
!
!   // Default mouse bindings
!   SetBind(BIND_NAV_ROTATE,GLFW_MOUSE_LEFT,NOMOD);
!   SetBind(BIND_NAV_TRANSLATE,GLFW_MOUSE_RIGHT,NOMOD);
!   SetBind(BIND_NAV_ZOOM,GLFW_MOUSE_SCROLL,NOMOD);
!   SetBind(BIND_NAV_RESET,GLFW_MOUSE_LEFT_DOUBLE,NOMOD);
  end subroutine SetDefaultKeyBindings

!   function IsModPressed(mod)
!     use gui_main, only: io
!     integer, intent(in) :: mod
!     logical :: IsModPressed

! !    iand(mod,GLFW_MOD_SHIFT)
! !   return (!(mod & GLFW_MOD_SHIFT) != io.KeyShift) &&
! !          (!(mod & GLFW_MOD_CONTROL) != io.KeyCtrl) &&
! !          (!(mod & GLFW_MOD_ALT) != io.KeyAlt) &&
! !        (!(mod & GLFW_MOD_SUPER) != io.KeySuper);
! ! }
!   end function IsModPressed

  ! Return whether the bind event is happening. If held (optional),
  ! the event happens only if the button is held down (for mouse).
  module function IsBindEvent(bind,held)
    use gui_interfaces_cimgui
    use gui_main, only: io
    integer, intent(in) :: bind
    logical, intent(in), optional :: held
    logical :: IsBindEvent

    integer :: key, mod
    logical :: held_

    ! process options
    held_ = .false.
    if (present(held)) held_ = held

    ! some checks
    IsBindEvent = .false.
    if (bindevent_level > 0) return
    if (bind < 1 .or. bind > BIND_NUM) return

    ! get current key and mod for this bind
    key = keybind(bind)
    mod = modbind(bind)

    if (key == ImGuiKey_None .or..not.igIsKeyDown(mod)) then
       ! no key or the mod is not down
       return
    elseif (key >= ImGuiKey_NamedKey_BEGIN .and. key < ImGuiKey_NamedKey_END .and.&
       .not.io%WantCaptureKeyboard .and..not.io%WantTextInput) then
       ! correct key ID and not keyboard captured or inputing text
       if (held_) then
          IsBindEvent = igIsKeyPressed(key,.true._c_bool)
       else
          IsBindEvent = igIsKeyDown(key)
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

  end function IsBindEvent

end submodule proc

! // ImGui settings, default settings and keybindings
! ImGuiIO& io = GetIO();
! io.IniFilename = nullptr;
! // io.MouseDrawCursor = true; // can't get anything other than the arrow otherwise
! DefaultSettings();
! SetDefaultKeyBindings();

! // Quit key binding
! if (IsBindEvent(BIND_QUIT,false))
!   glfwSetWindowShouldClose(rootwin, GLFW_TRUE);

! // key binds accessible to other files (to check for scroll)
! extern std::map<std::tuple<int,int,int>,int> keymap; // [key,mod,group] -> bind
!
! void SetDefaultKeyBindings();
!
! void RegisterCallback(int event,void *callback,void *data);
!
! void ProcessCallbacks();
!
! bool IsBindEvent(int bind,bool held);
!
! std::string BindKeyName(int event);
!
! void SetBind(int bind, int key, int mod);
!
! void SetBindEventLevel(int level=0);
!
! bool SetBindFromUserInput(int level);


! using namespace ImGui;
!
! std::map<std::tuple<int,int,int>,int> keymap = {}; // [key,mod] -> bind
!


! std::string BindKeyName(int bind){
!   std::string ckey = "";
!   if (keybind[bind] > 0){
!     const char *key = glfwGetKeyName(keybind[bind],0);
!     if (key){
!       ckey = std::string(key);
!       for (int i = 0; i < ckey.length(); i++)
!       ckey[i] = toupper(ckey[i]);
!     } else {
!       switch(keybind[bind]){
!       case GLFW_KEY_SPACE: ckey = "Space"; break;
!       case GLFW_KEY_WORLD_1: ckey = "World key 1"; break;
!       case GLFW_KEY_WORLD_2: ckey = "World key 2"; break;
!       case GLFW_KEY_ESCAPE: ckey = "Escape"; break;
!       case GLFW_KEY_ENTER: ckey = "Enter"; break;
!       case GLFW_KEY_TAB: ckey = "Tab"; break;
!       case GLFW_KEY_BACKSPACE: ckey = "Backspace"; break;
!       case GLFW_KEY_INSERT: ckey = "Insert"; break;
!       case GLFW_KEY_DELETE: ckey = "Delete"; break;
!       case GLFW_KEY_RIGHT: ckey = "Right"; break;
!       case GLFW_KEY_LEFT: ckey = "Left"; break;
!       case GLFW_KEY_DOWN: ckey = "Down"; break;
!       case GLFW_KEY_UP: ckey = "Up"; break;
!       case GLFW_KEY_PAGE_UP: ckey = "Page Up"; break;
!       case GLFW_KEY_PAGE_DOWN: ckey = "Page Down"; break;
!       case GLFW_KEY_HOME: ckey = "Home"; break;
!       case GLFW_KEY_END: ckey = "End"; break;
!       case GLFW_KEY_CAPS_LOCK: ckey = "Caps Lock"; break;
!       case GLFW_KEY_SCROLL_LOCK: ckey = "Scroll Lock"; break;
!       case GLFW_KEY_NUM_LOCK: ckey = "Num Lock"; break;
!       case GLFW_KEY_PRINT_SCREEN: ckey = "Print Screen"; break;
!       case GLFW_KEY_PAUSE: ckey = "Pause"; break;
!       case GLFW_KEY_F1: ckey = "F1"; break;
!       case GLFW_KEY_F2: ckey = "F2"; break;
!       case GLFW_KEY_F3: ckey = "F3"; break;
!       case GLFW_KEY_F4: ckey = "F4"; break;
!       case GLFW_KEY_F5: ckey = "F5"; break;
!       case GLFW_KEY_F6: ckey = "F6"; break;
!       case GLFW_KEY_F7: ckey = "F7"; break;
!       case GLFW_KEY_F8: ckey = "F8"; break;
!       case GLFW_KEY_F9: ckey = "F9"; break;
!       case GLFW_KEY_F10: ckey = "F10"; break;
!       case GLFW_KEY_F11: ckey = "F11"; break;
!       case GLFW_KEY_F12: ckey = "F12"; break;
!       case GLFW_KEY_F13: ckey = "F13"; break;
!       case GLFW_KEY_F14: ckey = "F14"; break;
!       case GLFW_KEY_F15: ckey = "F15"; break;
!       case GLFW_KEY_F16: ckey = "F16"; break;
!       case GLFW_KEY_F17: ckey = "F17"; break;
!       case GLFW_KEY_F18: ckey = "F18"; break;
!       case GLFW_KEY_F19: ckey = "F19"; break;
!       case GLFW_KEY_F20: ckey = "F20"; break;
!       case GLFW_KEY_F21: ckey = "F21"; break;
!       case GLFW_KEY_F22: ckey = "F22"; break;
!       case GLFW_KEY_F23: ckey = "F23"; break;
!       case GLFW_KEY_F24: ckey = "F24"; break;
!       case GLFW_KEY_F25: ckey = "F25"; break;
!       case GLFW_KEY_KP_0: ckey = "Numpad 0"; break;
!       case GLFW_KEY_KP_1: ckey = "Numpad 1"; break;
!       case GLFW_KEY_KP_2: ckey = "Numpad 2"; break;
!       case GLFW_KEY_KP_3: ckey = "Numpad 3"; break;
!       case GLFW_KEY_KP_4: ckey = "Numpad 4"; break;
!       case GLFW_KEY_KP_5: ckey = "Numpad 5"; break;
!       case GLFW_KEY_KP_6: ckey = "Numpad 6"; break;
!       case GLFW_KEY_KP_7: ckey = "Numpad 7"; break;
!       case GLFW_KEY_KP_8: ckey = "Numpad 8"; break;
!       case GLFW_KEY_KP_9: ckey = "Numpad 9"; break;
!       case GLFW_KEY_KP_DECIMAL: ckey = "Numpad ."; break;
!       case GLFW_KEY_KP_DIVIDE: ckey = "Numpad /"; break;
!       case GLFW_KEY_KP_MULTIPLY: ckey = "Numpad *"; break;
!       case GLFW_KEY_KP_SUBTRACT: ckey = "Numpad -"; break;
!       case GLFW_KEY_KP_ADD: ckey = "Numpad +"; break;
!       case GLFW_KEY_KP_ENTER: ckey = "Numpad Enter"; break;
!       case GLFW_KEY_KP_EQUAL: ckey = "Numpad ="; break;
!       case GLFW_KEY_MENU: ckey = "Menu"; break;
!
!       case GLFW_MOUSE_LEFT: ckey = "Left Mouse"; break;
!       case GLFW_MOUSE_LEFT_DOUBLE: ckey = "Double Left Mouse"; break;
!       case GLFW_MOUSE_RIGHT: ckey = "Right Mouse"; break;
!       case GLFW_MOUSE_RIGHT_DOUBLE: ckey = "Double Right Mouse"; break;
!       case GLFW_MOUSE_MIDDLE: ckey = "Middle Mouse"; break;
!       case GLFW_MOUSE_MIDDLE_DOUBLE: ckey = "Double Middle Mouse"; break;
!       case GLFW_MOUSE_BUTTON3: ckey = "Button3 Mouse"; break;
!       case GLFW_MOUSE_BUTTON3_DOUBLE: ckey = "Double Button3 Mouse"; break;
!       case GLFW_MOUSE_BUTTON4: ckey = "Button4 Mouse"; break;
!       case GLFW_MOUSE_BUTTON4_DOUBLE: ckey = "Double Button4 Mouse"; break;
!       case GLFW_MOUSE_SCROLL: ckey = "Mouse Wheel"; break;
!       }
!     }
!   }
!
!   return std::string((modbind[bind] & GLFW_MOD_SHIFT)?"Shift+":"") +
!     std::string((modbind[bind] & GLFW_MOD_CONTROL)?"Ctrl+":"") +
!     std::string((modbind[bind] & GLFW_MOD_ALT)?"Alt+":"") +
!     std::string((modbind[bind] & GLFW_MOD_SUPER)?"Super+":"") +
!     ckey;
! }
!
! static void EraseBind_(int key,int mod,int group){
!   if (keymap.find(std::make_tuple(key,mod,group)) != keymap.end()) {
!     int oldbind = keymap[std::make_tuple(key,mod,group)];
!     modbind[oldbind] = NOMOD;
!     keybind[oldbind] = NOKEY;
!     keymap.erase(std::make_tuple(key,mod,group));
!   }
! }

! void SetBindEventLevel(int level/*=0*/){
!   bindevent_level = level;
! }
!
! bool SetBindFromUserInput(int bind){
!   int mouse, key, newkey;
!   float scroll;
!   bool changed = true;
!
!   ImGui_ImplGlfwGL3_GetKeyMouseEvents(&mouse, &key, &scroll);
!
!   ImGuiIO& io = GetIO();
!   int mod = (io.KeyCtrl?GLFW_MOD_CONTROL:0x0000) | (io.KeyShift?GLFW_MOD_SHIFT:0x0000) |
!     (io.KeyAlt?GLFW_MOD_ALT:0x0000) | (io.KeySuper?GLFW_MOD_SUPER:0x0000);
!
!   if (scroll != 0.f){
!     newkey = GLFW_MOUSE_SCROLL;
!   } else if (mouse >= 0) {
!     switch (mouse){
!     case GLFW_MOUSE_BUTTON_LEFT:
!       newkey = GLFW_MOUSE_LEFT; break;
!     case GLFW_MOUSE_BUTTON_RIGHT:
!       newkey = GLFW_MOUSE_RIGHT; break;
!     case GLFW_MOUSE_BUTTON_MIDDLE:
!       newkey = GLFW_MOUSE_MIDDLE; break;
!     case GLFW_MOUSE_BUTTON_4:
!       newkey = GLFW_MOUSE_BUTTON3; break;
!     case GLFW_MOUSE_BUTTON_5:
!       newkey = GLFW_MOUSE_BUTTON4; break;
!     default:
!       changed = false;
!     }
!   } else if (key != NOKEY &&
!            key != GLFW_KEY_LEFT_SHIFT && key != GLFW_KEY_RIGHT_SHIFT &&
!            key != GLFW_KEY_LEFT_CONTROL && key != GLFW_KEY_RIGHT_CONTROL &&
!            key != GLFW_KEY_LEFT_ALT && key != GLFW_KEY_RIGHT_ALT &&
!            key != GLFW_KEY_LEFT_SUPER && key != GLFW_KEY_RIGHT_SUPER){
!     newkey = key;
!   } else {
!     changed = false;
!   }
!
!   if (changed)
!     SetBind(bind,newkey,mod);
!
!   return changed;
! }
!


