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

#include "menu.h"
#include "keybinding.h"
#include "dialog.h"

#include "imgui/imgui.h"

using namespace ImGui;

void ShowMenu(GLFWwindow* rootwin){
  
  if (BeginMainMenuBar()){
    if (BeginMenu("File")){
      if (MenuItem("Quit",BindKeyName(BIND_QUIT).c_str()))
	glfwSetWindowShouldClose(rootwin, GLFW_TRUE);
      EndMenu();
    }
    if (BeginMenu("Edit")){
      if (MenuItem("Preferences..."))
	OpenDialog(DLG_Preferences);
      EndMenu();
    }
    if (BeginMenu("View")){
      if (MenuItem("Tree",NULL,dlgopen[DLG_Tree]))
	ToggleDialog(DLG_Tree);
      if (MenuItem("Preferences",NULL,dlgopen[DLG_Preferences]))
	ToggleDialog(DLG_Preferences);
      if (MenuItem("Structural Information",NULL,dlgopen[DLG_StructInfo]))
	ToggleDialog(DLG_StructInfo);
      EndMenu();
    }
    SameLine(0, GetContentRegionAvailWidth()-180.);
    Text("%.3f ms/frame (%.1f FPS)", 1000.0f / GetIO().Framerate, GetIO().Framerate);
  }
  EndMainMenuBar();
}
