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

#include "imgui.h"
#include "guiapps.h"
#include "critic2.h"
#include "global.h"

// Tooltip delay
static const float ttipdelay = 1.5;


static bool IsItemHoveredDelayed(float delay,float *time0,bool *reset);
static void AttachTooltip(const char* desc, float delay, float *time0, bool *reset);

void show_menu_bar(bool *want_quit){
  static float time0 = -1.;
  bool reset = true;

  // immediate actions
  if (ImGui::BeginMainMenuBar()){
    if (ImGui::BeginMenu("File")){
      if (ImGui::MenuItem("New","Ctrl+N")){structurenew_window_h = true;}
      AttachTooltip("Create a structure from scratch.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Open crystal","Ctrl+O")) {structureopen_window_h = 2;}
      AttachTooltip("Read the crystal structure from a file.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Open molecule","Ctrl+Alt+O")) {structureopen_window_h = 1;}
      AttachTooltip("Read the molecular structure from a file.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Open from library","Ctrl+L",false,false)) {}

      if (ImGui::MenuItem("Open recent",NULL,false,false)) {}

      ImGui::Separator();

      if (ImGui::MenuItem("Close","Ctrl+W")) {clear_scene(true);}
      AttachTooltip("Clear the current structure.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Quit","Ctrl+Q")){*want_quit = true;}
      AttachTooltip("Quit the program.\n",ttipdelay,&time0,&reset);

      ImGui::EndMenu();
      if (reset)
	time0 = -1.;
    }
    if (ImGui::BeginMenu("Calculate")) {
      if (ImGui::MenuItem("Generate Critical Points")) {call_auto();}
      AttachTooltip("Calculate the critical points.\nBleh and Blah!\n",ttipdelay,&time0,&reset);
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("View")) {
      if (ImGui::MenuItem("Toggle bonds","",show_bonds)) {show_bonds = !show_bonds;}
      AttachTooltip("Toggle show/hide bonds.\n",ttipdelay,&time0,&reset);
      
      if (ImGui::MenuItem("Toggle critical points","",show_cps)) {show_cps = !show_cps;}
      AttachTooltip("Toggle show/hide critical points.\n",ttipdelay,&time0,&reset);
      
      if (ImGui::MenuItem("Toggle atoms","",show_atoms)) {show_atoms = !show_atoms;}
      AttachTooltip("Toggle show/hide atoms.\n",ttipdelay,&time0,&reset);
      
      if (ImGui::MenuItem("Toggle cell","",show_cell)) {show_cell = !show_cell;}
      AttachTooltip("Toggle show/hide unit cell.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Show structure information","",structureinfo_window_h)) {structureinfo_window_h = !structureinfo_window_h;}
      AttachTooltip("Show information about the current structure.\n",ttipdelay,&time0,&reset);

      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  if (reset)
    time0 = -1.;
}

static bool IsItemHoveredDelayed(float delay,float *time0,bool *reset)
{
  float time = ImGui::GetTime();
  if (ImGui::IsItemHovered()){
    if (*time0 < 0.)
      *time0 = time;
    *reset = false;
  } 
  return (*time0 > 0.) && (time > *time0 + delay) && ImGui::IsItemHovered();
}

static void AttachTooltip(const char* desc, float delay, float *time0, bool *reset)
{
  if (IsItemHoveredDelayed(delay,time0,reset))
    {
      ImGui::BeginTooltip();
      ImGui::PushTextWrapPos(450.0f);
      ImGui::TextUnformatted(desc);
      ImGui::PopTextWrapPos();
      ImGui::EndTooltip();
    }
}

