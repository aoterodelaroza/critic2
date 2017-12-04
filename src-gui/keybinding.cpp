/*
  Copyright (c) 2017 Alberto Otero de la Roza
  <aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
  Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
  Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
  <victor@fluor.quimica.uniovi.es>.

  critic2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at
  your option) any later version.

  critic2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "keybinding.h"
#include "imgui/imgui.h"
#include "imgui/imgui_internal.h"
#include <map>

using namespace std;
using namespace ImGui;

static int modbind[BIND_MAX]; // bind -> mod
int keybind[BIND_MAX]; // bind -> key
static map<pair<int,int>,int> keymap = {}; // [key,mod] -> bind

static bool IsModPressed(int mod){
  ImGuiIO& io = GetIO();
  return (!(mod & GLFW_MOD_SHIFT) != io.KeyShift) &&
         (!(mod & GLFW_MOD_CONTROL) != io.KeyCtrl) &&
         (!(mod & GLFW_MOD_ALT) != io.KeyAlt) &&
	 (!(mod & GLFW_MOD_SUPER) != io.KeySuper);
}

void RegisterDefaultBindings(){
  // Initialize to no keys and null callbacks
  for (int i = 0; i < BIND_MAX; i++){
    modbind[i] = NOMOD;
    keybind[i] = NOKEY;
  }

  // Default keybindings
  keybind[BIND_QUIT] = GLFW_KEY_Q;
  modbind[BIND_QUIT] = GLFW_MOD_CONTROL;
  keybind[BIND_CLOSE_LAST_DIALOG] = GLFW_KEY_ESCAPE;
  keybind[BIND_CLOSE_ALL_DIALOGS] = GLFW_KEY_DELETE;

  // Default mouse bindings
  keybind[BIND_NAV_ROTATE] = GLFW_MOUSE_LEFT;
  keybind[BIND_NAV_TRANSLATE] = GLFW_MOUSE_RIGHT;
  keybind[BIND_NAV_ZOOM] = GLFW_MOUSE_SCROLL;
  keybind[BIND_NAV_RESET] = GLFW_MOUSE_LEFT_DOUBLE;

  // Fill the keymap
  for (int i = 0; i < BIND_MAX; i++){
    if (keybind[i] != NOKEY && keybind[i] <= GLFW_KEY_LAST)
      keymap[make_pair(keybind[i],modbind[i])] = i;
  }
}

bool IsBindEvent(int bind, bool held){
  ImGuiIO& io = GetIO();
  int key = keybind[bind];
  int mod = modbind[bind];

  if (key == NOKEY || !IsModPressed(mod))
    return false;
  else if (key <= GLFW_KEY_LAST && !io.WantCaptureKeyboard && !io.WantTextInput)
    if (!held)
      return IsKeyPressed(key,false);
    else
      return IsKeyDown(key);
  else{
    if (key == GLFW_MOUSE_LEFT)
      if (!held)
	return IsMouseClicked(0);
      else
	return IsMouseDown(0);
    else if (key == GLFW_MOUSE_RIGHT)
      if (!held)
	return IsMouseClicked(1);
      else
	return IsMouseDown(1);
    else if (key == GLFW_MOUSE_MIDDLE)
      if (!held)
	return IsMouseClicked(2);
      else
	return IsMouseDown(2);
    else if (key == GLFW_MOUSE_BUTTON3)
      if (!held)
	return IsMouseClicked(3);
      else
	return IsMouseDown(3);
    else if (key == GLFW_MOUSE_BUTTON4)
      if (!held)
	return IsMouseClicked(4);
      else
	return IsMouseDown(4);
    else if (key == GLFW_MOUSE_LEFT_DOUBLE && !held)
      return IsMouseDoubleClicked(0);
    else if (key == GLFW_MOUSE_RIGHT_DOUBLE && !held)
      return IsMouseDoubleClicked(1);
    else if (key == GLFW_MOUSE_MIDDLE_DOUBLE && !held)
      return IsMouseDoubleClicked(2);
    else if (key == GLFW_MOUSE_BUTTON3_DOUBLE && !held)
      return IsMouseDoubleClicked(3);
    else if (key == GLFW_MOUSE_BUTTON4_DOUBLE && !held)
      return IsMouseDoubleClicked(4);
    else if (key == GLFW_MOUSE_SCROLL)
      return abs(GetCurrentContext()->IO.MouseWheel) > 1e-8;
    return false;
  }
}

string BindKeyName(int bind){
  string ckey = "";
  if (keybind[bind] > 0){
    const char *key = glfwGetKeyName(keybind[bind],0);
    if (key){
      ckey = string(key);
      for (int i = 0; i < ckey.length(); i++)
	ckey[i] = toupper(ckey[i]);
    } else {
      switch(keybind[bind]){
      case GLFW_MOUSE_LEFT: ckey = "Left Click"; break;
      case GLFW_MOUSE_LEFT_DOUBLE: ckey = "Double Click"; break;
      case GLFW_MOUSE_RIGHT: ckey = "Right Click"; break;
      case GLFW_MOUSE_RIGHT_DOUBLE: ckey = "Double Right Click"; break;
      case GLFW_MOUSE_MIDDLE: ckey = "Middle Click"; break;
      case GLFW_MOUSE_MIDDLE_DOUBLE: ckey = "Double Middle Click"; break;
      case GLFW_MOUSE_BUTTON3: ckey = "Button3 Click"; break;
      case GLFW_MOUSE_BUTTON3_DOUBLE: ckey = "Double Button3 Click"; break;
      case GLFW_MOUSE_BUTTON4: ckey = "Button4 Click"; break;
      case GLFW_MOUSE_BUTTON4_DOUBLE: ckey = "Double Button4 Click"; break;
      case GLFW_MOUSE_SCROLL: ckey = "Mouse Wheel"; break;
      }
    }
  }

  return string((modbind[bind] & GLFW_MOD_SHIFT)?"Shift+":"") +
    string((modbind[bind] & GLFW_MOD_CONTROL)?"Ctrl+":"") +
    string((modbind[bind] & GLFW_MOD_ALT)?"Alt+":"") +
    string((modbind[bind] & GLFW_MOD_SUPER)?"Super+":"") + 
    ckey;
}

