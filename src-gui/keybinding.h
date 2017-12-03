// -*- c++ -*-
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

#ifndef KEYBINDING_H
#define KEYBINDING_H

#include <GLFW/glfw3.h>
#include <string>

using namespace std;

#define NOKEY 0
#define NOMOD 0x0000

#define BIND_QUIT 0
#define BIND_MAX 1 // number of BIND actions

void RegisterDefaultKeyBindings();

void RegisterCallback(int event,void *callback,void *data);

void ProcessCallbacks();

string EventKeyName(int event);

#endif
