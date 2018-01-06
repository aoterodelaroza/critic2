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

#include <list>
#include <string>

#include "imgui/imgui.h"

#include "tree.h"
#include "critic2.h"

using namespace ImGui;

// 
// struct sscene{
//   int id;
// };
// struct sfile{
//   std::string filename;
//   std::list<sscene> sc;
// };
// static std::list<sfile> tree = {};

static int nfiles = 0;
static std::string *files = nullptr;
static std::string *paths = nullptr;

void UpdateTreeData(){
  nfiles = 0;
  if (files) delete files;
  if (paths) delete paths;

  nfiles = c2::nfiles;
  files = new std::string[nfiles];
  paths = new std::string[nfiles];

  for (int isc=0; isc<c2::nsc; isc++){
    c2::set_scene_pointers(isc+1);
    if (c2::isinit > 0){
      paths[c2::idfile-1] = files[c2::idfile-1] = std::string(c2::file);
      const size_t ll = files[c2::idfile-1].find_last_of("\\/");
      if (std::string::npos != ll)
       	files[c2::idfile-1].erase(0, ll + 1);
      if (std::string::npos != 0)
	paths[c2::idfile-1].erase(ll + 1,std::string::npos);
    }
  }
}

void ShowTree(){
  for (int i=0; i<c2::nfiles; i++){
    
    if (TreeNodeEx((void*)(intptr_t)i,0,files[i].c_str())){
      Text("Blah!");
      TreePop();
    }
  }
}
