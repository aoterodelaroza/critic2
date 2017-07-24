#include "imgui.h"
#define IMGUI_DEFINE_PLACEMENT_NEW
#include "imgui_internal.h"
#include "imgui_dock.h"

static ImVec2 operator+(ImVec2 lhs, ImVec2 rhs) {
    return ImVec2(lhs.x+rhs.x, lhs.y+rhs.y);
}
static ImVec2 operator-(ImVec2 lhs, ImVec2 rhs) {
    return ImVec2(lhs.x-rhs.x, lhs.y-rhs.y);
}
using namespace ImGui;

static DockContext g_dock;

//xx// DockContext methods //xx//

Dock *DockContext::getContainerAt(const ImVec2& pos){
  Dock *dd;
  for (auto dpair : this->dockht){
    dd = dpair.second;
    if (dd->type != Dock::Type_Container) continue;
    if (dd->status != Dock::Status_Open) continue;
    if (IsMouseHoveringRect(dd->pos,dd->pos+dd->size,false))
      return dd;
  }
  return nullptr;
}

//xx// Dock methods //xx//

void Dock::showDropTarget(){
  SetNextWindowSize(ImVec2(0,0));
  Begin("##Overlay",nullptr,ImGuiWindowFlags_Tooltip|ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|
        ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_AlwaysAutoResize);
  ImDrawList* canvas = GetWindowDrawList();
  canvas->PushClipRectFullScreen();
  ImU32 docked_color = GetColorU32(ImGuiCol_FrameBg);
  docked_color = (docked_color & 0x00ffFFFF) | 0x80000000;
  canvas->AddRectFilled(this->pos, this->pos + this->size, docked_color);
  canvas->PopClipRect();
  End();
}

void Dock::newDock(Dock *dnew){
  if (this->sttype == Dock::Stack_None){
    this->sttype = Dock::Stack_Leaf;
    this->stack.push_back(dnew);
    dnew->parent = this;
    dnew->root = this->root;
    this->currenttab = dnew;
  } else if (this->sttype == Dock::Stack_Leaf){
    this->stack.push_back(dnew);
    dnew->parent = this;
    dnew->root = this->root;
    this->currenttab = dnew;
  }
}

void Dock::drawContainer(){
    
  if (this->sttype == Dock::Stack_Leaf && this->stack.size() > 0){
    // Draw the tab
    this->drawTabBar();

    // Hide all tabs
    for (auto dd : this->stack) {
      dd->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar|
        ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoInputs;
      dd->pos = ImVec2(0.,0.);
      dd->size = ImVec2(0.,0.);
      dd->hidden = true;
    }

    // Unhide current tab
    if (!this->hidden && !this->collapsed){
      Dock *dd = this->currenttab;
      float h = 2 * GetTextLineHeightWithSpacing() + GetCurrentWindow()->TitleBarHeight();
      dd->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_ShowBorders|ImGuiWindowFlags_NoResize|
	ImGuiWindowFlags_NoCollapse|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoBringToFrontOnFocus;
      dd->pos = this->pos + ImVec2(0.,h);
      dd->size = this->size - ImVec2(0.,h);
      dd->hidden = false;
    }
  }
}

void Dock::drawTabBar(){
  float barheight = 2 * GetTextLineHeightWithSpacing();

  SetCursorScreenPos(this->pos + ImVec2(0.,GetCurrentWindow()->TitleBarHeight()));
  if (BeginChild("tabs", ImVec2(this->size.x,barheight), true)){

    // drawTabbarListButton(dock);
    for (auto dd : this->stack) {
      if (Button(dd->label))
	this->currenttab = dd;
      // currenttab and dragging
      // currenttab and close button
      SameLine();
    }
  }
  EndChild();
}

//xx// Public interface //xx//

void ImGui::Container(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags extra_flags /*= 0*/){

  bool collapsed = true;
  ImGuiWindowFlags flags = extra_flags | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoBringToFrontOnFocus;

  Dock *dd = g_dock.dockht[label];
  if (dd){
    // SetNextWindowPos(dd->pos);
    // SetNextWindowSize(dd->size);
  } else {
    dd = g_dock.dockht[label] = new Dock;
    IM_ASSERT(dd);
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
  dd->type = Dock::Type_Container;
  dd->root = dd;
  dd->collapsed = collapsed;
  dd->window = GetCurrentWindow();

  // Update the status
  ImGuiContext *g = GetCurrentContext();
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

  // Draw the container elements
  dd->drawContainer();

  End();
}

bool ImGui::BeginDock(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags extra_flags /*= 0*/){

  bool collapsed = true;
  ImGuiWindowFlags flags = extra_flags;

  Dock *dd = g_dock.dockht[label];
  if (!dd) {
    dd = g_dock.dockht[label] = new Dock;
    IM_ASSERT(dd);
  }

  if (dd->status == Dock::Status_Docked){
    // Docked: position, size, and flags controlled by the container
    SetNextWindowPos(dd->pos);
    SetNextWindowSize(dd->size);
    flags = dd->flags;
    collapsed = dd->hidden;
    if (dd->hidden){
      Begin("",nullptr,dd->size,0.0,flags);
    } else {
      Begin(label,nullptr,flags);
    }
  } else {
    // Floating window
    collapsed = !Begin(label,p_open,flags);
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

  // Update the status
  ImGuiContext *g = GetCurrentContext();
  Dock *ddest = g_dock.getContainerAt(GetIO().MousePos);
  if (g->ActiveId == GetCurrentWindow()->MoveId && g->IO.MouseDown[0]){
    // Dragging
    dd->status = Dock::Status_Dragged;
  } else {
    if (dd->status == Dock::Status_Dragged && ddest){
      // Just stopped dragging and there is a container below
      dd->status = Dock::Status_Docked;
      ddest->newDock(dd);
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
    Dock *ddest = g_dock.getContainerAt(GetIO().MousePos);
    if (ddest)
      ddest->showDropTarget();
  }

  return !collapsed;
}

void ImGui::EndDock() {
  End();
}

void ImGui::Print() {
  for (auto dock : g_dock.dockht){
    Text("key=%s id=%d label=%s\n", dock.first,dock.second->id,dock.second->label);
    Text("pos=(%f,%f) size=(%f,%f)\n",dock.second->pos.x,dock.second->pos.y,dock.second->size.x,dock.second->size.y);
    Text("type=%d status=%d\n", dock.second->type, dock.second->status);
    Text("sttype=%d list_size=%d\n", dock.second->sttype, dock.second->stack.size());
  }
}
