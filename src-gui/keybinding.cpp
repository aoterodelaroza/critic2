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
#include <map>

using namespace std;
using namespace ImGui;

static int modbind[BIND_MAX]; // bind -> mod
static int keybind[BIND_MAX]; // bind -> key
static void *eventbind[BIND_MAX]; // bind -> event
static void *databind[BIND_MAX]; // bind -> data
static map<pair<int,int>,int> keymap = {}; // [key,mod] -> bind

void RegisterDefaultKeyBindings(){
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

  // Fill the keymap
  for (int i = 0; i < BIND_MAX; i++)
    keymap[make_pair(keybind[i],modbind[i])] = i;
}

void RegisterCallback(int event,void *callback,void *data){
  eventbind[event] = callback;
  databind[event] = data;
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

string EventKeyName(int event){
  string ckey = "";
  if (keybind[event] > 0){
    ckey = string(glfwGetKeyName(keybind[event],0));
    for (int i = 0; i < ckey.length(); i++)
      ckey[i] = toupper(ckey[i]);
  }

  return string((modbind[event] & GLFW_MOD_SHIFT)?"Shift+":"") +
    string((modbind[event] & GLFW_MOD_CONTROL)?"Ctrl+":"") +
    string((modbind[event] & GLFW_MOD_ALT)?"Alt+":"") +
    string((modbind[event] & GLFW_MOD_SUPER)?"Super+":"") + 
    ckey;
}

