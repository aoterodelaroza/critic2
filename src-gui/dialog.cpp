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

#include "imgui/imgui_dock.h"

#include "dialog.h"
#include "view.h"
#include "settings.h"
#include "keybinding.h"

using namespace ImGui;

// Declarations for the dialog static functions (defined below)
static void DialogPreferences(bool *p_open);

// Dialog variables
static bool dlgopen[DLG_LAST] = {false};
static void (*dlgfun[DLG_LAST])(bool *) = {
  DialogPreferences,
};
static Dock *dlgdock[DLG_LAST] = {nullptr};
static Dialog_ dlglastopen = DLG_LAST;

void OpenDialog(Dialog_ dialog){
  if (!dlgopen[dialog])
    dlgopen[dialog] = true;
  else if (dlgdock[dialog] && dlgdock[dialog]->window){
    if (dlgdock[dialog]->status == Dock::Status_Docked){
      dlgdock[dialog]->parent->currenttab = dlgdock[dialog];
      dlgdock[dialog]->parent->focusContainer();      
    } else {
      FocusWindow(dlgdock[dialog]->window);
    }
  }
  dlglastopen = dialog;
}

void DialogDispatch(){
  // Detect the close-all and close-one
  if (IsBindEvent(BIND_CLOSE_LAST_DIALOG,false))
    CloseLastDialog();
  if (IsBindEvent(BIND_CLOSE_ALL_DIALOGS,false))
    CloseAllDialogs();

  // Process all the flags
  for (int i = 0; i < DLG_LAST; i++){
    if (dlgopen[i])
      dlgfun[i](&(dlgopen[i]));
  }
}

void CloseAllDialogs(){
  for (int i = 0; i < DLG_LAST; i++){
    dlgopen[i] = false;
    if (dlgdock[i])
      dlgdock[i]->closeDock();
  }
  dlglastopen = DLG_LAST;
}

void CloseLastDialog(){
  if (dlglastopen != DLG_LAST){
    dlgopen[dlglastopen] = false;
    if (dlgdock[dlglastopen])
      dlgdock[dlglastopen]->closeDock();
  }
  dlglastopen = DLG_LAST;
}

// Static functions for all dialogs known to the dispatcher //

// Preferences dialog //
static void DialogPreferences(bool *p_open){
  if (*p_open){
    SetNextWindowSize(ImVec2(500, 440), ImGuiSetCond_FirstUseEver);
    if (BeginDock("Preferences",p_open)){
      // Filter box
      ImGuiTextFilter Filter;
      Text("Filter");
      SameLine();
      Filter.Draw("",-1.f);

      // Selectable categories
      const int ncat = 5;
      char *catname[ncat] = {
	"General","Views","Key bindings","Interface","Fonts"
      };

      // Left panel
      static int catid = 0;
      BeginChild("leftpanel", ImVec2(150, 0), true);
      for (int i=0; i < ncat; i++){
	if (Selectable(catname[i], catid == i))
	  catid = i;
      }
      EndChild();
      SameLine();

      // Right panel
      BeginGroup();
      BeginChild("rightpanel", ImVec2(0,-GetItemsLineHeightWithSpacing()));
      Text(catname[catid]);
      Separator();
      if (catid == 0){
	// General
	if (TreeNode("Tooltips")){
	  Checkbox("Enable tooltips", &tooltip_enabled);
	  PushItemWidth(40); DragFloat("Tooltip delay (s)", &tooltip_delay, 0.1f, 0.0f, FLT_MAX, "%.1f", 1.0f);
	  PushItemWidth(40); DragFloat("Tooltip maximum width (pixel)", &tooltip_maxwidth, 5.0f, 0.0f, FLT_MAX, "%.1f", 1.0f);
	  TreePop();
	}
      } else if (catid == 1){
	// Views
	if (TreeNode("Lighting")){
	  // Views -> Lighting
	  bool changed = false, anychanged = false;
	  PushItemWidth(40 * 3); 
	  changed |= DragFloat3("Light position", &(view_lightpos[0]), 0.1f, -FLT_MAX, FLT_MAX, "%.1f", 1.0f); 
	  if (changed)
	    shader->setVec3("lightPos",value_ptr(view_lightpos));
	  anychanged |= changed;
	  
	  PushItemWidth(40 * 3); 
	  changed |= DragFloat3("Light color", &(view_lightcolor[0]), 0.01f, 0.0, 1.0f, "%.2f", 1.0f); 
	  if (changed)
	    shader->setVec3("lightColor",value_ptr(view_lightcolor));
	  anychanged |= changed;

	  changed = false;
	  PushItemWidth(40); 
	  changed |= DragFloat("Ambient light intensity", &view_ambient, 0.01f, 0.0f, 1.0f, "%.2f", 1.0f); 
	  if (changed)
	    shader->setFloat("ambient",view_ambient);
	  anychanged |= changed;

	  changed = false;
	  PushItemWidth(40); 
	  changed |= DragFloat("Diffuse light intensity", &view_diffuse, 0.01f, 0.0f, 1.0f, "%.2f", 1.0f); 
	  if (changed)
	    shader->setFloat("diffuse",view_diffuse);
	  anychanged |= changed;

	  changed = false;
	  PushItemWidth(40); 
	  changed |= DragFloat("Specular light intensity", &view_specular, 0.01f, 0.0f, 1.0f, "%.2f", 1.0f); 
	  if (changed)
	    shader->setFloat("specular",view_specular);
	  anychanged |= changed;

	  changed = false;
	  PushItemWidth(40); 
	  changed |= DragInt("Light shininess", &view_shininess, 1.0f, 0.0f, FLT_MAX, "%.2f"); 
	  if (changed)
	    shader->setInt("shininess",view_shininess);
	  anychanged |= changed;

	  if (anychanged)
	    ForceUpdateAllViews();
	  TreePop();
	}
	if (TreeNode("Mouse sensitivity")){
	  PushItemWidth(40); 
	  DragFloat("Rotation mouse sensitivity", &view_mousesens_rot, 0.1f, 0.0f, FLT_MAX, "%.1f", 1.0f); 
	  PushItemWidth(40); 
	  DragFloat("Zoom mouse sensitivity", &view_mousesens_zoom, 0.1f, 0.0f, FLT_MAX, "%.1f", 1.0f); 
	  TreePop();
	}
	if (TreeNode("Per-view settings (to be moved)")){
	  bool changed = false;
	  changed |= Checkbox("Wireframe rendering", &view_wireframe);
	  changed |= Checkbox("Orthgonal projection", &view_orthogonal);
	  PushItemWidth(40); 
	  changed |= DragFloat("Field of view (degrees)", &view_fov, 2.5f, 0.0f, 180.0f, "%.1f", 1.0f); 
	  PushItemWidth(40); 
	  changed |= DragFloat("Reset distance (scene radius)", &view_resetdistance, 0.05f, 0.0f, FLT_MAX, "%.2f", 1.0f); 
	  PushItemWidth(4*40); 
	  changed |= DragFloat4("Background color", view_bgrgb, 0.01f, 0.0, 1.0f, "%.2f", 1.0f); 
	  PushItemWidth(4*40); 
	  changed |= DragInt("Atom resolution", &view_isphres, 1.0f, 0.0f, 3.0f, "%.2f"); 
	  PushItemWidth(4*40); 
	  changed |= DragInt("Bond resolution", &view_icylres, 0.0f, 0.0f, 0.0f, "%.2f"); 
	  view_icylres = 0;

	  if (changed)
	    SetDefaultAllViews();
	  TreePop();
	}
      } else if (catid == 2){
	// Key bindings
      } else if (catid == 3){
	// Interface
      } else if (catid == 4){
	// Fonts
      }
      EndChild();

      // Line at the bottom
      BeginChild("buttons");
      if (Button("Reset to Defaults")){}
      SameLine();
      if (Button("Reset to File")){}
      SameLine();
      if (Button("OK")){
	*p_open = false;
      }
      EndChild();
      EndGroup();

      // if (Filter.IsActive()){
      // 	const char* buf_begin = Buf.begin();
      // 	const char* line = buf_begin;
      // 	for (int line_no = 0; line != NULL; line_no++)
      // 	  {
      // 	    const char* line_end = (line_no < LineOffsets.Size) ? buf_begin + LineOffsets[line_no] : NULL;
      // 	    if (Filter.PassFilter(line, line_end))
      // 	      TextUnformatted(line, line_end);
      // 	    line = line_end && line_end[1] ? line_end + 1 : NULL;
      // 	  }
      // }
      // else{
      // 	TextUnformatted(Buf.begin());
      // }

    }
    dlgdock[DLG_Preferences] = GetCurrentDock();
    EndDock();
  }
}
