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

#ifndef GLOBAL_H
#define GLOBAL_H

class Settings {
  public:

  // global signal for program termination
  bool want_quit;

  // show elements in the scene
  bool show_bonds;
  bool show_cps;
  bool show_atoms;
  bool show_cell;

  // characteristics of the elements in the scene
  char bondresolution;
  char atomresolution;
  float bondthickness;
  float atomsize;
  float cpsize;

  // camera position
  float cam_pos[3];
  float cam_target[3];
  float cam_up[3];

  // constructor
  Settings(bool ismolecule, float maxlen);

  // set the camera position
  void set_flags_and_cam(bool ismolecule, float maxlen, float maxclen);
  void set_flags(bool ismolecule);
  void set_cam_pos(float maxlen);
};

extern Settings settings;

#endif
