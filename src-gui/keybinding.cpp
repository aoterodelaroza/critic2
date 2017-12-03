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
static void *eventbind[BIND_MAX]; // bind -> event
static void *databind[BIND_MAX]; // bind -> data
static map<pair<int,int>,int> keymap = {}; // [key,mod] -> bind

void RegisterDefaultBindings(){
  // Initialize to no keys and null callbacks
  for (int i = 0; i < BIND_MAX; i++){
    modbind[i] = NOMOD;
    keybind[i] = NOKEY;
    eventbind[i] = nullptr;
    databind[i] = nullptr;
  }

  // Default keybindings
  keybind[BIND_QUIT] = GLFW_KEY_Q;
  modbind[BIND_QUIT] = GLFW_MOD_CONTROL;

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

void RegisterCallback(int bind,void *callback,void *data){
  eventbind[bind] = callback;
  databind[bind] = data;
}

void ProcessCallbacks(){
  ImGuiIO& io = GetIO();

  for (auto it=keymap.begin(); it != keymap.end(); it++){
    if (eventbind[it->second] &&
	IsKeyPressed(it->first.first,false) && 
	(!(it->first.second & GLFW_MOD_SHIFT) != io.KeyShift) &&
	(!(it->first.second & GLFW_MOD_CONTROL) != io.KeyCtrl) &&
	(!(it->first.second & GLFW_MOD_ALT) != io.KeyAlt) &&
	(!(it->first.second & GLFW_MOD_SUPER) != io.KeySuper)){
      ((void (*)(void *)) eventbind[it->second])(databind[it->second]);
    }
  }
}

bool IsBindEvent(int bind, bool held){
  int key = keybind[bind];

  if (key == NOKEY)
    return false;
  else if (key <= GLFW_KEY_LAST)
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
    ckey = string(glfwGetKeyName(keybind[bind],0));
    for (int i = 0; i < ckey.length(); i++)
      ckey[i] = toupper(ckey[i]);
  }

  return string((modbind[bind] & GLFW_MOD_SHIFT)?"Shift+":"") +
    string((modbind[bind] & GLFW_MOD_CONTROL)?"Ctrl+":"") +
    string((modbind[bind] & GLFW_MOD_ALT)?"Alt+":"") +
    string((modbind[bind] & GLFW_MOD_SUPER)?"Super+":"") + 
    ckey;
}

