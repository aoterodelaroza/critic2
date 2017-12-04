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

#include "dialog.h"
#include "imgui/imgui_dock.h"
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
      dlgdock[i]->CloseDock();
  }
  dlglastopen = DLG_LAST;
}

void CloseLastDialog(){
  if (dlglastopen != DLG_LAST){
    dlgopen[dlglastopen] = false;
    if (dlgdock[dlglastopen])
      dlgdock[dlglastopen]->CloseDock();
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
