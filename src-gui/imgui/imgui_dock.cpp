// do not Begin hidden windows
// tabbar list too long and scrolling.
// relocate inside the bar tab
// drag by grabbing inside the window
// compact styles ; handle styles 

#include "imgui.h"
#define IMGUI_DEFINE_PLACEMENT_NEW
#include "imgui_internal.h"
#include "imgui_dock.h"
#include "imgui_impl_glfw.h"
#include <unordered_map>

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

// Dock context declarations
struct LabelHash{
  size_t operator()(const char *key) const{
    return (size_t) ImHash((const void *)key,0);
  };
};
static unordered_map<const char *,Dock*,LabelHash> dockht = {}; // global dock hash table
static Dock *getContainerAt(const ImVec2& pos); // get container at a given position

//xx// DockContext methods //xx//

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
  ImGuiContext *g = GetCurrentContext();
  float barheight = g->FontSize + 6 * g->Style.ItemSpacing.y;
  float tabheight = g->FontSize + 2 * g->Style.ItemSpacing.y;
  float crossz = 0.3 * g->FontSize;

  SetCursorScreenPos(this->pos + ImVec2(0.,GetCurrentWindow()->TitleBarHeight()));
  char tmp[20];
  ImFormatString(tmp,IM_ARRAYSIZE(tmp),"tabs%d",(int)this->id);
  if (BeginChild(tmp,ImVec2(this->size.x,barheight), true)){
    ImDrawList* draw_list = GetWindowDrawList();
    ImU32 text_color = GetColorU32(ImGuiCol_Text);
    ImU32 text_color_disabled = GetColorU32(ImGuiCol_TextDisabled);
    ImU32 color = GetColorU32(ImGuiCol_FrameBg);
    ImU32 color_active = GetColorU32(ImGuiCol_FrameBgActive);
    ImU32 color_hovered = GetColorU32(ImGuiCol_FrameBgHovered);
    SetCursorScreenPos(this->pos + ImVec2(0.,GetCurrentWindow()->TitleBarHeight()));

    // drawTabbarListButton(dock);

    bool active;
    Dock *dderase = nullptr, *ddlast = nullptr;
    char tmp2[20];
    ImVec2 pos, center, pos0;
    float cursorx = 0.;
    bool hovered = false;
    for (auto dd : this->stack) {
      SameLine(cursorx, 0.);
      const char* text_end = FindRenderedTextEnd(dd->label);
      int tabwidth = CalcTextSize(dd->label, text_end).x + 2 * g->Style.ItemSpacing.x;
      ImVec2 size(tabwidth, tabheight);

      // main button
      if (InvisibleButton(dd->label, size))
        this->currenttab = dd;
      active = (dd == this->currenttab);
      pos0 = GetItemRectMin();

      hovered = IsItemHovered();
      // lift the tab using the main button
      if (IsItemActive() && IsMouseDragging()){
        dd->status = Dock::Status_Dragged;
        goto lift_this_tab;
      }
      // double click detaches the tab
      if (IsItemActive() && IsMouseDoubleClicked(0)){
        dd->status = Dock::Status_Open;
        goto lift_this_tab;
      }
      cursorx += size.x;
      if (cursorx == size.x) cursorx += g->Style.ItemSpacing.x;

      // The close button, if this window can be closed
      if (dd->p_open){
        // the close button itself
        SameLine(cursorx,0.);
        ImFormatString(tmp2,IM_ARRAYSIZE(tmp2),"##%d",(int)dd->id);
	ImVec2 size(tabheight, tabheight);
	tabwidth = tabwidth + tabheight;
	if (InvisibleButton(tmp2, size)){
          *(dd->p_open) = false;
          goto erase_this_tab;
        }
	hovered |= IsItemHovered();

	// tab being lifted using the x
        if (IsItemActive() && IsMouseDragging()){
          dd->status = Dock::Status_Dragged;
          goto lift_this_tab;
        }
	cursorx += size.x;
      }
      center = ((GetItemRectMin() + GetItemRectMax()) * 0.5f);
      cursorx += g->Style.ItemSpacing.x;

      // draw the rectangle and the text
      draw_list->AddRectFilled(pos0,pos0 + ImVec2(tabwidth,tabheight), hovered? color_hovered: (active? color_active: color));
      draw_list->AddText(pos0 + g->Style.ItemSpacing, text_color, dd->label, text_end);
      if (dd->p_open){
        // draw the "x"
	draw_list->AddLine( center + ImVec2(-crossz, -crossz), center + ImVec2( crossz, crossz), text_color_disabled);
	draw_list->AddLine( center + ImVec2( crossz, -crossz), center + ImVec2(-crossz, crossz), text_color_disabled);
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
  } // BeginChild(tmp, ImVec2(this->size.x,barheight), true)
  EndChild();
}

//xx// Public interface //xx//

void ImGui::Container(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags extra_flags /*= 0*/){

  bool collapsed = true;
  // ImGuiWindowFlags flags = extra_flags | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoBringToFrontOnFocus;
  ImGuiWindowFlags flags = extra_flags;

  Dock *dd = dockht[label];
  if (dd){
    // SetNextWindowPos(dd->pos);
    // SetNextWindowSize(dd->size);
  } else {
    dd = dockht[label] = new Dock;
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
  Dock *ddest = getContainerAt(GetIO().MousePos);
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
    Dock *ddest = getContainerAt(GetIO().MousePos);
    if (ddest)
      ddest->showDropTarget();
  }

  // If this window is docked, put it right on top of its container
  if (dd->status == Dock::Status_Docked && !dd->hidden){
    ImGuiContext *g = GetCurrentContext();
    int ithis = -1, icont = -1;
    for (int i = 0; i < g->Windows.Size; i++){
      if (g->Windows[i] == dd->window)
	ithis = i;
      if (g->Windows[i] == dd->parent->window)
	icont = i;
    }
    if (icont >= 0 && ithis >= 0){
      g->Windows.erase(g->Windows.begin() + ithis);
      if (icont >= g->Windows.size())
	g->Windows.push_back(dd->window);
      else
	g->Windows.insert(g->Windows.begin() + icont + 1,dd->window);
    }
  }

  return !collapsed;
}

void ImGui::EndDock() {
  End();
}

void ImGui::Print() {
  for (auto dock : dockht){
    Text("key=%s id=%d label=%s\n", dock.first,dock.second->id,dock.second->label);
    Text("pos=(%f,%f) size=(%f,%f)\n",dock.second->pos.x,dock.second->pos.y,dock.second->size.x,dock.second->size.y);
    Text("type=%d status=%d\n", dock.second->type, dock.second->status);
    Text("sttype=%d list_size=%d\n", dock.second->sttype, dock.second->stack.size());
  }
}

void ImGui::ShutdownDock(){
  for (auto dpair : dockht){
    if (dpair.second) delete dpair.second;
  }
  dockht.clear();
}

