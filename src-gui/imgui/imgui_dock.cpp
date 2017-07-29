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
// Rewritten from: git@github.com:vassvik/imgui_docking_minimal.git

// pass the container as an argument to begindock, to initialize the window as docked.
// save the settings to a config file
// strange size effect on docking

#include "imgui.h"
#define IMGUI_DEFINE_PLACEMENT_NEW
#include "imgui_internal.h"
#include "imgui_dock.h"
#include "imgui_widgets.h"
#include "imgui_impl_glfw.h"
#include <unordered_map>
#include <math.h>

using namespace ImGui;

// Dock context declarations
struct LabelHash{
  size_t operator()(const char *key) const{
    return (size_t) ImHash((const void *)key,0);
  };
};
static unordered_map<const char *,Dock*,LabelHash> dockht = {}; // global dock hash table
static Dock *getContainerAt(const ImVec2& pos); // get container at a given position

//xx// Dock context methods //xx//

static Dock *getContainerAt(const ImVec2& pos){
  Dock *dd;
  for (auto dpair : dockht){
    dd = dpair.second;
    if (dd->type != Dock::Type_Container) continue;
    if (dd->status != Dock::Status_Open) continue;
    if (IsMouseHoveringRect(dd->pos,dd->pos+dd->size,false))
      return dd;
  }
  return nullptr;
}

//xx// Dock methods //xx//

bool Dock::IsMouseHoveringTabBar(){
  const ImVec2 ytabcushion = ImVec2(0.f,10.);
  return !this->stack.empty() && IsMouseHoveringRect(this->tabbarrect.Min-ytabcushion,this->tabbarrect.Max+ytabcushion,false);
}

int Dock::getNearestTabBorder(){
  if (!this->IsMouseHoveringTabBar()) return -1;
  int ithis = this->tabsx.size()-1;
  float xpos = GetMousePos().x;
  for (int i = 0; i < this->tabsx.size(); i++){
    if (xpos < this->tabsx[i]){
      ithis = i;
      break;
    }
  }
  if (ithis > 0 && (xpos-this->tabsx[ithis-1]) < (this->tabsx[ithis]-xpos))
    ithis = ithis - 1;
  return ithis;
}

void Dock::showDropTargetFull(){
  SetNextWindowSize(ImVec2(0,0));
  Begin("##Overlay",nullptr,ImGuiWindowFlags_Tooltip|ImGuiWindowFlags_NoTitleBar|
        ImGuiWindowFlags_NoInputs|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_AlwaysAutoResize);
  ImDrawList* drawl = GetWindowDrawList();
  drawl->PushClipRectFullScreen();
  ImU32 docked_color = GetColorU32(ImGuiCol_FrameBg);
  docked_color = (docked_color & 0x00ffFFFF) | 0x80000000;
  drawl->AddRectFilled(this->pos, this->pos + this->size, docked_color);
  drawl->PopClipRect();
  End();
}

void Dock::showDropTargetOnTabBar(){
  const float triside = 10.f;

  int ithis = this->getNearestTabBorder();

  SetNextWindowSize(ImVec2(0,0));
  Begin("##Overlay",nullptr,ImGuiWindowFlags_Tooltip|ImGuiWindowFlags_NoTitleBar|
        ImGuiWindowFlags_NoInputs|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_AlwaysAutoResize);
  ImDrawList* drawl = GetWindowDrawList();
  drawl->PushClipRectFullScreen();
  ImU32 docked_color = GetColorU32(ImGuiCol_FrameBg);
  docked_color = (docked_color & 0x00ffFFFF) | 0x80000000;

  ImVec2 a, b, c;
  a.x = this->tabsx[ithis]; a.y = this->tabbarrect.Max.y;
  b.x = a.x + 0.5 * triside; b.y = a.y + triside * sqrt(3.)/2.;
  c.x = a.x - 0.5 * triside; c.y = b.y;
  drawl->AddTriangleFilled(a,b,c,docked_color);
  a.y = this->tabbarrect.Min.y;
  b.y = a.y - triside * sqrt(3.)/2.;
  c.y = b.y;
  drawl->AddTriangleFilled(a,b,c,docked_color);

  drawl->PopClipRect();
  End();
}

void Dock::newDock(Dock *dnew, int ithis /*=-1*/){
  if (this->sttype == Dock::Stack_None || this->sttype == Dock::Stack_Leaf){
    this->sttype = Dock::Stack_Leaf;
    dnew->parent = this;
    dnew->root = this->root;
    this->currenttab = dnew;
    if (ithis < 0 || ithis == this->stack.size())
      this->stack.push_back(dnew);
    else{
      int n = -1;
      for (auto it = this->stack.begin(); it != this->stack.end(); ++it){
	n++;
	if (n == ithis){
	  this->stack.insert(it,dnew);
	  break;
	}
      }
    }
  }
}

void Dock::drawContainer(bool allowresize){
    
  if (this->sttype == Dock::Stack_Leaf && this->stack.size() > 0){
    // Draw the tab
    this->drawTabBar();

    // Hide all tabs
    for (auto dd : this->stack) {
      dd->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoScrollbar|
        ImGuiWindowFlags_NoScrollWithMouse|ImGuiWindowFlags_NoCollapse|
        ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoInputs;
      dd->pos = ImVec2(0.,0.);
      dd->size = ImVec2(0.,0.);
      dd->hidden = true;
    }

    // Unhide current tab
    if (!this->hidden && !this->collapsed && this->currenttab){
      Dock *dd = this->currenttab;
      float topheight = this->tabbarrect.Max.y - this->pos.y;
      dd->pos = this->pos;
      dd->pos.y += topheight;
      dd->size = this->size;
      dd->size.y = this->size.y - topheight;
      dd->hidden = false;
      dd->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|
        ImGuiWindowFlags_NoCollapse|ImGuiWindowFlags_NoSavedSettings|
	ImGuiWindowFlags_NoBringToFrontOnFocus;
      if (!allowresize)
	dd->flags = dd->flags | ImGuiWindowFlags_NoResize;
    }
  } // if (this->sttype == Dock::Stack_Leaf && this->stack.size() > 0);
}

void Dock::clearContainer(){
  const float increment = 50.;

  ImVec2 pos = this->pos;
  for (auto dd : this->stack) {
    dd->status = Dock::Status_Open;
    dd->control_window_this_frame = true;
    dd->size = dd->size_saved;
    dd->collapsed = dd->collapsed_saved;
    dd->pos = pos;
    pos = pos + ImVec2(increment,increment);
  }
  this->currenttab = nullptr;
  this->stack.clear();
}

void Dock::drawTabBar(){
  ImGuiContext *g = GetCurrentContext();
  const float tabheight = GetTextLineHeightWithSpacing();
  const float barheight = 2 * GetTextLineHeightWithSpacing();
  const float maxtabwidth = 100.;
  ImVec4 text_color = g->Style.Colors[ImGuiCol_Text];
  text_color.w = 2.0 / g->Style.Alpha;

  // empty the list of tabs
  this->tabsx.Size = 0;

  // calculate the widths
  float tabwidth_long;
  if ((this->size.x - 2 * g->Style.ItemSpacing.x) >= this->stack.size() * maxtabwidth){
    tabwidth_long = maxtabwidth;
  } else{
    tabwidth_long = round(this->size.x - 2 * g->Style.ItemSpacing.x) / this->stack.size();
  }

  // the tabbar with alpha = 1.0, no spacing in the x
  PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.0,1.0));
  PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0.0,0.0));
  PushStyleColor(ImGuiCol_Text,text_color);

  // start placing the tabs
  char tmp[20];
  ImFormatString(tmp,IM_ARRAYSIZE(tmp),"tabs%d",(int)this->id);
  if (BeginChild(tmp,ImVec2(this->size.x,barheight),false)){
    bool active = false, hovered = false;
    Dock *dderase = nullptr, *ddlast = nullptr;
    ImVec2 center, pos0, pos1, pos1s, text_size;
    for (auto dd : this->stack) {
      SameLine();
      // make the x-button, update the container info
      bool dragged, dclicked, closeclicked;
      if (ButtonWithX(dd->label, ImVec2(tabwidth_long, tabheight), (dd == this->currenttab), false,
       		      dd->p_open, &dragged, &dclicked, &closeclicked, 2.f/g->Style.Alpha))
	this->currenttab = dd;
      this->tabsx.push_back(GetItemRectMax().x-tabwidth_long);
    
      // lift the tab using the main button
      if (dragged){
        dd->status = Dock::Status_Dragged;
        goto lift_this_tab;
      }
      // double click detaches the tab
      if (dclicked){
        dd->status = Dock::Status_Open;
        goto lift_this_tab;
      }
      // closed click kills the tab
      if (closeclicked){
	*(dd->p_open) = false;
	goto erase_this_tab;
      }

      ddlast = dd;
      continue;

    lift_this_tab:
      dd->control_window_this_frame = true;
      dd->size = dd->size_saved;
      dd->collapsed = dd->collapsed_saved;
      dd->pos = GetMousePos() - ImVec2(0.5*dd->size.x,barheight);

    erase_this_tab:
      dderase = dd;
      if (dd == this->currenttab){
        if (ddlast)
          this->currenttab = ddlast;
        else
          this->currenttab = nullptr;
      }
    } // dd in this->stack
    if (dderase){
      this->stack.remove(dderase);
      if (!this->currenttab && this->stack.size() > 0)
        this->currenttab = this->stack.front();
    }

    // If the tab (but not the buttons) is clicked, transfer to the
    // container. Raise the container & docked window to the top of
    // the window stack. 
    if (this->currenttab && !IsAnyItemActive()){
      if (g->HoveredWindow == GetCurrentWindow())
        this->SetContainerHoveredMovedActive(true);
      this->RaiseCurrentTab();
    }

    // last item in the tabsx
    this->tabsx.push_back(GetItemRectMax().x);
  } // BeginChild(tmp, ImVec2(this->size.x,barheight), true)
  this->tabbarrect = GetCurrentWindowRead()->DC.LastItemRect;
  this->tabbarrect.Min.x = this->pos.x;
  this->tabbarrect.Max.x = this->pos.x + this->size.x;
  EndChild();
  PopStyleVar();
  PopStyleVar();
  PopStyleColor();
}

void Dock::SetContainerHoveredMovedActive(bool setid){
  ImGuiContext *g = GetCurrentContext();

  if (g->IO.MouseClicked[0]){
    g->HoveredRootWindow = this->currenttab->window;
    g->HoveredWindow = this->currenttab->window;
    g->MovedWindow = this->window;
    g->MovedWindowMoveId = this->window->RootWindow->MoveId;
    if (setid)
      SetActiveID(g->MovedWindowMoveId, this->currenttab->window->RootWindow);
  }
}

void Dock::RaiseCurrentTab(){
  ImGuiContext *g = GetCurrentContext();

  if (this->currenttab){
    int ithis = -1, icont = -1;
    for (int i = 0; i < g->Windows.Size; i++){
      if (g->Windows[i] == this->currenttab->window)
        ithis = i;
      if (g->Windows[i] == this->window)
        icont = i;
    }
    if (icont >= 0 && ithis >= 0){
      g->Windows.erase(g->Windows.begin() + ithis);
      if (ithis < icont) icont--;
      if (icont >= g->Windows.size())
        g->Windows.push_back(this->currenttab->window);
      else
        g->Windows.insert(g->Windows.begin() + icont + 1,this->currenttab->window);
    }
  }
}

int Dock::getWindowStackPosition(){
  ImGuiContext *g = GetCurrentContext();

  int iexit = -1;
  for (int i = 0; i < g->Windows.Size; i++)
    if (g->Windows[i] == this->window){
      iexit = i;
      break;
    }
  return iexit;
}

//xx// Public interface //xx//

void ImGui::Container(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags extra_flags /*= 0*/){

  bool collapsed = true;
  ImGuiContext *g = GetCurrentContext();
  ImGuiWindowFlags flags = extra_flags;

  Dock *dd = dockht[label];
  if (!dd){
    dd = dockht[label] = new Dock;
    IM_ASSERT(dd);
  }

  // If the container has a window docked, set the minimum size
  if (dd->currenttab){
    ImVec2 size = dd->currenttab->window->SizeContents + dd->currenttab->window->WindowPadding;
    size.y += GetTextLineHeightWithSpacing() + dd->currenttab->window->WindowPadding.y +
      dd->window->TitleBarHeight();
    SetNextWindowSizeConstraints(size,ImVec2(FLT_MAX,FLT_MAX),nullptr,nullptr);
    flags = flags | ImGuiWindowFlags_NoResize;
  }

  // Render any container widgets in here
  if (Begin(label,p_open,flags)){
    collapsed = false;
  }

  // Fill the info for this dock
  dd->id = ImHash(label,0);
  dd->label = ImStrdup(label);
  dd->pos = GetWindowPos();
  dd->size = GetWindowSize();
  dd->size_saved = ImVec2(0.,0.);
  dd->type = Dock::Type_Container;
  dd->root = dd;
  dd->collapsed = collapsed;
  dd->window = GetCurrentWindow();
  dd->p_open = p_open;

  // Update the status
  if (g->ActiveId == GetCurrentWindow()->MoveId && g->IO.MouseDown[0]){
    // Dragging
    dd->status = Dock::Status_Dragged;
  } else {
    // Stationary -> open, closed, or collapsed
    if (!p_open || *p_open){
      if (collapsed){
        dd->status = Dock::Status_Collapsed;
      }
      else{
        dd->status = Dock::Status_Open;
      }
    } else {
      dd->status = Dock::Status_Closed;
    }
  }

  // If the container has just been closed, detach all docked windows
  if (dd->status == Dock::Status_Closed && !dd->stack.empty())
    dd->clearContainer();

  // Draw the container elements
  dd->drawContainer(!(extra_flags & ImGuiWindowFlags_NoResize));
  dd->size_saved.y = dd->tabbarrect.Max.y - dd->pos.y;

  // If the container is clicked, set the correct hovered/moved flags 
  // and raise container & docked window to the top of the stack.
  if (dd->currenttab){
    if (g->HoveredWindow == dd->window)
      dd->SetContainerHoveredMovedActive(!IsAnyItemActive());
    dd->RaiseCurrentTab();
  }

  End();
}

bool ImGui::BeginDock(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags flags /*= 0*/){

  bool collapsed;

  ImGuiContext *g = GetCurrentContext();

  Dock *dd = dockht[label];
  if (!dd) {
    dd = dockht[label] = new Dock;
    IM_ASSERT(dd);
  }

  if (dd->status == Dock::Status_Docked || dd->control_window_this_frame){
    // Docked or lifted: position and size are controlled
    SetNextWindowPos(dd->pos);
    SetNextWindowSize(dd->size);
    SetNextWindowCollapsed(dd->collapsed);
    if (dd->status == Dock::Status_Docked) {
      // Docked: flags and hidden controlled by the container, too
      flags = dd->flags;
      collapsed = dd->hidden;
      if (dd->hidden){
        Begin("",nullptr,dd->size,0.0,flags);
      } else {
        Begin(label,nullptr,flags);
      }
    } else if (dd->status == Dock::Status_Dragged) { 
      // the window has just been lifted from a container. Go back to
      // being a normal window with the new position and size; being
      // dragged.
      collapsed = !Begin(label,p_open,flags);
      dd->window = GetCurrentWindow();
      FocusWindow(dd->window);
      g->MovedWindow = dd->window;
      g->MovedWindowMoveId = dd->window->RootWindow->MoveId;
      SetActiveID(g->MovedWindowMoveId, dd->window->RootWindow);
    } else { 
      // the window has just been lifted, but not dragging
      collapsed = !Begin(label,p_open,flags);
      dd->window = GetCurrentWindow();
      FocusWindow(dd->window);
    }
  } else {
    // Floating window
    collapsed = !Begin(label,p_open,flags);
    dd->collapsed_saved = collapsed;
    if (!collapsed) dd->size_saved = dd->size;
  }

  // Fill the info for this dock
  dd->id = ImHash(label,0);
  dd->label = ImStrdup(label);
  dd->pos = GetWindowPos();
  dd->size = GetWindowSize();
  dd->type = Dock::Type_Dock;
  dd->flags = flags;
  dd->collapsed = collapsed;
  dd->window = GetCurrentWindow();
  dd->p_open = p_open;
  dd->control_window_this_frame = false;

  // Update the status
  Dock *ddest;
  if (g->ActiveId == GetCurrentWindow()->MoveId && g->IO.MouseDown[0]){
    // Dragging
    dd->status = Dock::Status_Dragged;
    ddest = getContainerAt(GetIO().MousePos);
  } else {
    int ithis = -1;
    if (dd->status == Dock::Status_Dragged && (ddest = getContainerAt(GetIO().MousePos)) &&
	(ddest->stack.empty() || ((ithis = ddest->getNearestTabBorder()) >= 0))){	 
      // Just stopped dragging and there is a container below
      dd->status = Dock::Status_Docked;
      ddest->newDock(dd,ithis);
    } else if (dd->status != Dock::Status_Docked){
      // Stationary -> open, closed, or collapsed
      if (!p_open || *p_open){
        if (collapsed){
          dd->status = Dock::Status_Collapsed;
        }
        else{
          dd->status = Dock::Status_Open;
        }
      } else {
        dd->status= Dock::Status_Closed;
      }
    }
  }

  // If dragged and hovering over a container, show the drop rectangles
  if (dd->status == Dock::Status_Dragged){
    if (ddest){
      if (ddest->stack.empty())
	ddest->showDropTargetFull();
      else if (ddest->IsMouseHoveringTabBar())
	ddest->showDropTargetOnTabBar();
    }
  }

  // If this window is docked, 
  if (dd->status == Dock::Status_Docked && !dd->hidden){
    // Transfer any unhandled clicks to the underlying container.
    if ((g->HoveredWindow == dd->window || g->HoveredWindow == dd->parent->window) &&
        IsMouseHoveringRect(dd->window->Pos,dd->window->Pos+dd->window->Size,false) &&
        g->IO.MouseClicked[0]){
      dd->parent->SetContainerHoveredMovedActive(!(IsMouseClicked(0) && IsAnyItemHovered()));
    }
    // Put it right on top of its container.
    dd->parent->RaiseCurrentTab();
    // If the resize grip is being used, resize the container too.
    if (g->ActiveId == dd->window->GetID("#RESIZE")){
      dd->parent->window->SizeFull = dd->window->Size + dd->parent->size_saved;
    }
  }

  return !collapsed;
}

void ImGui::EndDock() {
  End();
}

void ImGui::Print() {
  // for (auto dock : dockht){
  //   Text("key=%s id=%d label=%s\n", dock.first,dock.second->id,dock.second->label);
  //   Text("pos=(%f,%f) size=(%f,%f)\n",dock.second->pos.x,dock.second->pos.y,dock.second->size.x,dock.second->size.y);
  //   Text("type=%d status=%d\n", dock.second->type, dock.second->status);
  //   Text("sttype=%d list_size=%d\n", dock.second->sttype, dock.second->stack.size());
  //   if (dock.second->p_open)
  //     Text("p_open=%d\n", *(dock.second->p_open));
  // }

  for (auto dock : dockht){
    Text("key=%s id=%d label=%s\n", dock.first,dock.second->id,dock.second->label);
    Text("pos=(%f,%f) size=(%f,%f)\n",dock.second->pos.x,dock.second->pos.y,dock.second->size.x,dock.second->size.y);
    Text("type=%d status=%d\n", dock.second->type, dock.second->status);
    Text("sttype=%d list_size=%d\n", dock.second->sttype, dock.second->stack.size());
    if (dock.second->p_open)
      Text("p_open=%d\n", *(dock.second->p_open));
  }

  // ImGuiContext *g = GetCurrentContext();
  // for (int i = 0; i < g->Windows.Size; i++){
  //   Text("%d %s %p\n",i,g->Windows[i]->Name,g->Windows[i]);
  // }

  // where is g.hoveredrootwindow set?
  // void ImGui::SetHoveredID(ImGuiID id)
  //    g.HoveredId = id;
  //    g.HoveredIdAllowOverlap = false;
}

void ImGui::ShutdownDock(){
  for (auto dpair : dockht){
    if (dpair.second) delete dpair.second;
  }
  dockht.clear();
}

