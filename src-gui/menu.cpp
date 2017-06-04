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
#include "settings.h"

#include "imguifilesystem.h"

// Tooltip delay
static const float ttipdelay = 1.5;


static bool IsItemHoveredDelayed(float delay,float *time0,bool *reset);
static void AttachTooltip(const char* desc, float delay, float *time0, bool *reset);

void show_menu_bar(){
  static float time0 = -1.;
  bool reset = true;
  static ImGuiFs::Dialog fsopenfile;
  static bool firstpass = true;
  static int getfile = 0;

  // immediate actions
  if (ImGui::BeginMainMenuBar()){
    if (ImGui::BeginMenu("File")){
      if (ImGui::MenuItem("New","Ctrl+N",false,!settings.preview_mode))
	structurenew_window_h = true;
      AttachTooltip("Create a structure from scratch.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Open crystal","Ctrl+O",false,!settings.preview_mode)) 
	structureopen_window_h = 2;
      AttachTooltip("Read the crystal structure from a file.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Open molecule","Ctrl+Alt+O",false,!settings.preview_mode))
	structureopen_window_h = 1;
      AttachTooltip("Read the molecular structure from a file.\n",ttipdelay,&time0,&reset);

      if (ImGui::BeginMenu("Crystal library",!settings.preview_mode)){
	if (ImGui::MenuItem("Choose file"))
	  getfile = 1;
	ImGui::Separator();
	for (int i = 0; i < nlib_crys; i++){
	  if (ImGui::MenuItem(lib_crys[i])){
	    open_structure_from_library(i+1,0);
	    settings.set_flags_and_cam(false,box_xmaxlen,box_xmaxclen);
	  }
	}
	ImGui::EndMenu();
      }
      AttachTooltip("Read a crystal structure from the library file.\n",ttipdelay,&time0,&reset);

      if (ImGui::BeginMenu("Molecule library",!settings.preview_mode)){
	if (ImGui::MenuItem("Choose file"))
	  getfile = 2;
	ImGui::Separator();
	for (int i = 0; i < nlib_mol; i++){
	  if (ImGui::MenuItem(lib_mol[i])){
	    open_structure_from_library(i+1,1);
	    settings.set_flags_and_cam(true,box_xmaxlen,box_xmaxclen);
	  }
	}
	ImGui::EndMenu();
      }
      AttachTooltip("Read a molecular structure from the library file.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Open recent",NULL,false,false)){}

      ImGui::Separator();

      if (ImGui::MenuItem("Close","Ctrl+W",false,!settings.preview_mode)) 
	clear_scene(true);
      AttachTooltip("Clear the current structure.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Quit","Ctrl+Q")){settings.want_quit = true;}
      AttachTooltip("Quit the program.\n",ttipdelay,&time0,&reset);

      ImGui::EndMenu();
      if (reset)
	time0 = -1.;
    }
    if (ImGui::BeginMenu("Calculate")) {
      if (ImGui::MenuItem("Generate Critical Points",NULL,false,!settings.preview_mode)) 
	call_auto();
      AttachTooltip("Calculate the critical points.\nBleh and Blah!\n",ttipdelay,&time0,&reset);
      ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("View")) {
      if (ImGui::MenuItem("Toggle bonds","",settings.show_bonds)) {settings.show_bonds = !settings.show_bonds;}
      AttachTooltip("Toggle show/hide bonds.\n",ttipdelay,&time0,&reset);
      
      if (ImGui::MenuItem("Toggle critical points","",settings.show_cps)) {settings.show_cps = !settings.show_cps;}
      AttachTooltip("Toggle show/hide critical points.\n",ttipdelay,&time0,&reset);
      
      if (ImGui::MenuItem("Toggle atoms","",settings.show_atoms)) {settings.show_atoms = !settings.show_atoms;}
      AttachTooltip("Toggle show/hide atoms.\n",ttipdelay,&time0,&reset);
      
      if (ImGui::MenuItem("Toggle cell","",settings.show_cell)) {settings.show_cell = !settings.show_cell;}
      AttachTooltip("Toggle show/hide unit cell.\n",ttipdelay,&time0,&reset);

      if (ImGui::MenuItem("Show structure information","",structureinfo_window_h)) {structureinfo_window_h = !structureinfo_window_h;}
      AttachTooltip("Show information about the current structure.\n",ttipdelay,&time0,&reset);

      ImGui::Separator();
      if (ImGui::MenuItem("Close all windows","ESC")){settings.close_all_windows = true;}
      AttachTooltip("Close all open windows.\n",ttipdelay,&time0,&reset);

      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  // Library file dialog
  if (getfile){
    const char* filename = fsopenfile.chooseFileDialog(firstpass,"./",NULL);
    firstpass = false;
    if (fsopenfile.hasUserJustCancelledDialog() && strlen(filename) == 0){
      firstpass = true;
      getfile = 0;
    }
    if (strlen(filename) > 0){
      set_library_file(&filename, getfile); 
      firstpass = true;
      getfile = 0;
    }
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

