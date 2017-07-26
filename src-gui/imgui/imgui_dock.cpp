// relocate inside the bar tab
// improve look of the drawtabbar
// floating levels of containers and tabbed windows are screwed
// tabbar list too long and scrolling.

#include "imgui.h"
#define IMGUI_DEFINE_PLACEMENT_NEW
#include "imgui_internal.h"
#include "imgui_dock.h"
#include "imgui_impl_glfw.h"

static ImVec2 operator+(ImVec2 lhs, ImVec2 rhs) {
    return ImVec2(lhs.x+rhs.x, lhs.y+rhs.y);
}
static ImVec2 operator-(ImVec2 lhs, ImVec2 rhs) {
    return ImVec2(lhs.x-rhs.x, lhs.y-rhs.y);
}
ImVec2 operator*(ImVec2 lhs, float rhs) {
    return ImVec2(lhs.x*rhs, lhs.y*rhs);
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
    if (!this->hidden && !this->collapsed && this->currenttab){
      Dock *dd = this->currenttab;
      float h = 2 * GetTextLineHeightWithSpacing() + GetCurrentWindow()->TitleBarHeight();
      dd->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_ShowBorders|ImGuiWindowFlags_NoResize|
        ImGuiWindowFlags_NoCollapse|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoBringToFrontOnFocus;
      dd->pos = this->pos + ImVec2(0.,h);
      dd->size = this->size - ImVec2(0.,h);
      dd->hidden = false;
    }
  } // if (this->sttype == Dock::Stack_Leaf && this->stack.size() > 0);
}

void Dock::drawTabBar(){
  float barheight = 2 * GetTextLineHeightWithSpacing();

  SetCursorScreenPos(this->pos + ImVec2(0.,GetCurrentWindow()->TitleBarHeight()));
  char tmp[20];
  ImFormatString(tmp,IM_ARRAYSIZE(tmp),"tabs%d",(int)this->id);
  if (BeginChild(tmp,ImVec2(this->size.x,barheight), true)){
    ImDrawList* draw_list = GetWindowDrawList();
    ImU32 color = GetColorU32(ImGuiCol_FrameBg);
    ImU32 color_active = GetColorU32(ImGuiCol_FrameBgActive);
    ImU32 color_hovered = GetColorU32(ImGuiCol_FrameBgHovered);
    ImU32 text_color = GetColorU32(ImGuiCol_Text);
    ImU32 text_color_disabled = GetColorU32(ImGuiCol_TextDisabled);
    float line_height = GetTextLineHeightWithSpacing();
    float tab_base;

    // drawTabbarListButton(dock);

    bool active;
    Dock *dderase = nullptr, *ddlast = nullptr;
    char tmp2[20];
    ImVec2 pos, center;
    for (auto dd : this->stack) {
      SameLine(0, 15);
      const char* text_end = FindRenderedTextEnd(dd->label);
      ImVec2 size(CalcTextSize(dd->label, text_end).x, line_height);
      if (InvisibleButton(dd->label, size))
        this->currenttab = dd;
      active = (dd == this->currenttab);

      if (IsItemActive() && IsMouseDragging()){
        dd->status = Dock::Status_Dragged;
        goto lift_this_tab;
      }
      if (IsItemActive() && IsMouseDoubleClicked(0)){
        dd->status = Dock::Status_Open;
        goto lift_this_tab;
      }

      pos = GetItemRectMin();
      size.x += GetStyle().ItemSpacing.x;
      tab_base = pos.y;

      // Rectangle and text for the tab
      draw_list->AddRectFilled(pos+ImVec2(-8.0f, 0.0),pos+size,
                               IsItemHovered() ? color_hovered : (active ? color_active : color));
      draw_list->AddText(pos, text_color, dd->label, text_end);

      // The close button, if this window can be closed
      if (dd->p_open){
        // the close button itself
        SameLine();
        ImFormatString(tmp2,IM_ARRAYSIZE(tmp2),"##%d",(int)dd->id);
        if (Button(tmp2, ImVec2(16, 16))){
          *(dd->p_open) = false;
          goto erase_this_tab;
        }

        if (IsItemActive() && IsMouseDragging()){
          dd->status = Dock::Status_Dragged;
          goto lift_this_tab;
        }

        // the "x"
        SameLine();
        center = ((GetItemRectMin() + GetItemRectMax()) * 0.5f);
        draw_list->AddLine( center + ImVec2(-3.5f, -3.5f), center + ImVec2(3.5f, 3.5f), text_color_disabled);
        draw_list->AddLine( center + ImVec2(3.5f, -3.5f), center + ImVec2(-3.5f, 3.5f), text_color_disabled);
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
    ImVec2 cp(this->pos.x, tab_base + line_height);
    draw_list->AddLine(cp, cp + ImVec2(this->size.x, 0), color);
    if (dderase){
      this->stack.remove(dderase);
      if (!this->currenttab && this->stack.size() > 0)
        this->currenttab = this->stack.front();
    }
  } // BeginChild(tmp, ImVec2(this->size.x,barheight), true)
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
  dd->p_open = p_open;

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

  bool collapsed;
  ImGuiWindowFlags flags = extra_flags;

  ImGuiContext *g = GetCurrentContext();
  Dock *dd = g_dock.dockht[label];
  if (!dd) {
    dd = g_dock.dockht[label] = new Dock;
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

