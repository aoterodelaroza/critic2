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
#include "draw.h"

#include "imguifilesystem.h"

// Variable definitions
bool structureinfo_window_h = false;
int structurenew_window_h = 0; // 0 - hidden, 1 - molecule 2 - crystal

// Process known handles for all windows
void guiapps_process_handles(){
  if (structureinfo_window_h) structureinfo_window(&structureinfo_window_h);
  if (structurenew_window_h > 0) structurenew_window(&structurenew_window_h);
}

// Open a layout showing the structural information for this 
// crystal/molecule
void structureinfo_window(bool *p_open){
  const char *title[] = {"Structure not available###structinfo", "Molecule information###structinfo", "Crystal information###structinfo"};

  const int elem_unselected = 0;
  const int elem_global_unavailable = 1;
  const int elem_global_molecule = 2;
  const int elem_global_crystal = 3;

  const int nelems = 4;
  const int elem_levels[nelems] = {0,0,1,2};

  const char *elem_title[] = {"", "", "General","General"};

  // Display different information depending on the type of structure (level)
  int level;
  if (!isinit)
    level = 0;
  else if (ismolecule)
    level = 1;
  else
    level = 2;

  static int selected = 0;

  ImGui::SetNextWindowSize(ImVec2(300, 300), ImGuiSetCond_Once);
  ImGui::SetNextWindowPos(ImVec2(0, 30), ImGuiSetCond_Once);
  if (ImGui::Begin(title[level], p_open)){
    ImGui::BeginGroup();

    // List of selectable elements on the left
    ImGui::BeginChild("left", ImVec2(150, 0), true);
    for (int i=2;i<nelems;i++){
      if (elem_levels[i] == level && ImGui::Selectable(elem_title[i], selected == i)){
	selected = i;
      }
    }
    ImGui::EndChild();
    ImGui::SameLine();
    if (elem_levels[selected] != level)
      selected = 0;

    ImGui::BeginChild("right", ImVec2(0, -ImGui::GetItemsLineHeightWithSpacing())); 
    ImGui::Text(elem_title[selected]);
    ImGui::Separator();
    if (selected == elem_global_molecule)
      ImGui::TextWrapped("This is a molecule!");
    else if (selected == elem_global_crystal)
      ImGui::TextWrapped("This is a crystal!");
    ImGui::EndChild();

    ImGui::EndGroup();
  }
  ImGui::End();
}

// Read the structure from an external file p_open = 0 - closed, 1 - molecule, 2 - crystal
void structurenew_window(int *p_open){
  static ImGuiFs::Dialog fsopenfile;
  static bool firstpass = true;

  const char* filename = fsopenfile.chooseFileDialog(firstpass,"./",NULL);
  firstpass = false;

  if (fsopenfile.hasUserJustCancelledDialog() && strlen(filename) == 0){
    // Dialog has been closed - set up for next time and prevent more calls for now
    firstpass = true;
    *p_open = 0;
  }

  if (strlen(filename) > 0){
    // Clean up previous and initialize the structure
    call_structure(&filename, *p_open == 1); 

    // Set default camera position, show cell if crystal, etc.
    draw_set_camera_pos(-1.);
    show_cell = (*p_open == 2);
    show_bonds = true;
    show_cps = true;
    show_atoms = true;

    // Close the dialog
    firstpass = true;
    *p_open = 0;
  }
}

