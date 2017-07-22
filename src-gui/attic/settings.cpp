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

#include "settings.h"

Settings settings = Settings(false,5.f);

Settings::Settings(bool ismolecule, float maxlen){
  want_quit = false;
  close_all_windows = false;
  preview_mode = false;

  bondresolution = 2;
  atomresolution = 1;

  bondthickness = 0.05;
  atomsize = 0.5;
  cpsize = 0.5;

  set_flags(ismolecule);
  set_cam_pos(maxlen);
}

void Settings::set_flags_and_cam(bool ismolecule, float maxlen, float maxclen){
  set_flags(ismolecule);
  if (ismolecule)
    set_cam_pos(maxlen);
  else
    set_cam_pos(maxclen>maxlen?maxclen:maxlen);
  rot.InitIdentity();
}

void Settings::set_flags(bool ismolecule){
  show_bonds = true;
  show_cps = true;
  show_atoms = true;
  show_cell = !ismolecule;
}

void Settings::set_cam_pos(float maxlen){
  cam_pos[0] = 0.f;
  cam_pos[1] = 0.f;
  cam_pos[2] = -2.f * maxlen;

  cam_target[0] = 0.f;
  cam_target[1] = 0.f;
  cam_target[2] = 1.f;

  cam_up[0] = 0.f;
  cam_up[1] = 1.f;
  cam_up[2] = 0.f;
}

