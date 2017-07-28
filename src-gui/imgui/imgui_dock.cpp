// relocate inside the bar tab
// resize the container if there are tabbed windows; problem with the corner
// reorganize methods
// clearcontainer gives a cascade. some control over the window stack

#include "imgui.h"
#define IMGUI_DEFINE_PLACEMENT_NEW
#include "imgui_internal.h"
#include "imgui_dock.h"
#include "imgui_impl_glfw.h"
#include <unordered_map>
#include <math.h>

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
  Begin("##Overlay",nullptr,ImGuiWindowFlags_Tooltip|ImGuiWindowFlags_NoTitleBar|
	ImGuiWindowFlags_NoInputs|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_AlwaysAutoResize);
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
    float posymax = 0.;
    this->drawTabBar(&posymax);

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
      dd->pos = this->pos;
      dd->pos.y = posymax;
      dd->size = this->size;
      dd->size.y = this->size.y - (posymax - this->pos.y);
      dd->hidden = false;
      dd->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoResize|
        ImGuiWindowFlags_NoCollapse|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoBringToFrontOnFocus;
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

void Dock::drawTabBar(float *posymax){
  ImGuiContext *g = GetCurrentContext();
  const float tabheight = GetTextLineHeightWithSpacing();
  const float barheight = 2 * GetTextLineHeightWithSpacing();
  const float crossz = round(0.3 * g->FontSize);
  const float crosswidth = 2 * crossz + 6;
  const float maxtabwidth = 100.;
  const float mintabwidth = 2 * crosswidth + 1;
  const ImU32 text_color = GetColorU32(ImGuiCol_Text);
  const ImU32 text_color_disabled = GetColorU32(ImGuiCol_TextDisabled);
  const ImU32 color = GetColorU32(ImGuiCol_FrameBg);
  const ImU32 color_active = GetColorU32(ImGuiCol_FrameBgActive);
  const ImU32 color_hovered = GetColorU32(ImGuiCol_FrameBgHovered);

  // calculate the widths
  float tabwidth_long;
  if ((this->size.x - 2 * g->Style.ItemSpacing.x) >= this->stack.size() * maxtabwidth){
    tabwidth_long = maxtabwidth;
  } else{
    tabwidth_long = round(this->size.x - 2 * g->Style.ItemSpacing.x) / this->stack.size();
  }
  float tabwidth_short = tabwidth_long - crosswidth;

  // the tabbar with alpha = 1.0, no spacing in the x
  PushStyleVar(ImGuiStyleVar_Alpha, 1.0);
  PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.0,1.0));
  PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0.0,0.0));

  // start placing the tabs
  char tmp[20];
  ImFormatString(tmp,IM_ARRAYSIZE(tmp),"tabs%d",(int)this->id);
  if (BeginChild(tmp,ImVec2(this->size.x,barheight),false)){
    *posymax = GetItemRectMax().y;
    ImDrawList* draw_list = GetWindowDrawList();

    bool active = false, hovered = false;
    Dock *dderase = nullptr, *ddlast = nullptr;
    char tmp2[20];
    ImVec2 center, pos0, pos1, pos1s, text_size;
    ImRect clip_rect;
    for (auto dd : this->stack) {
      // size of the main button
      int tabwidth;
      if (dd->p_open && tabwidth_long >= mintabwidth)
	tabwidth = tabwidth_short;
      else
	tabwidth = tabwidth_long;
      ImVec2 size(tabwidth, tabheight);

      // main button
      SameLine();
      if (InvisibleButton(dd->label, size))
        this->currenttab = dd;
      active = (dd == this->currenttab);
      hovered = IsItemHovered();
      pos0 = GetItemRectMin();
      pos1s = GetItemRectMax();
      const char* text_end = FindRenderedTextEnd(dd->label);

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

      // draw the close button, if this window can be closed
      if (dd->p_open && tabwidth_long >= mintabwidth){
        // draw the close button itself
	SameLine();
        ImFormatString(tmp2,IM_ARRAYSIZE(tmp2),"##%d",(int)dd->id);
	ImVec2 size(crosswidth, tabheight);
	tabwidth = tabwidth + crosswidth;
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

        // draw the "x"
	center = ((GetItemRectMin() + GetItemRectMax()) * 0.5f);
      	draw_list->AddLine( center + ImVec2(-crossz, -crossz), center + ImVec2( crossz, crossz), text_color_disabled);
      	draw_list->AddLine( center + ImVec2( crossz, -crossz), center + ImVec2(-crossz, crossz), text_color_disabled);
      }
      pos1 = GetItemRectMax();

      // rectangle and text
      text_size = CalcTextSize(dd->label,text_end,true,false);
      clip_rect = ImRect(pos0,pos1s);
      draw_list->AddRectFilled(pos0,pos1,hovered?color_hovered:(active?color_active:color));
      draw_list->AddRect(pos0,pos1,color_active,0.0f,~0,1.0f);
      ImGui::RenderTextClipped(pos0,pos1s,dd->label,text_end,&text_size, ImVec2(0.5f,0.5f), &clip_rect);

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

  } // BeginChild(tmp, ImVec2(this->size.x,barheight), true)
  EndChild();
  PopStyleVar();
  PopStyleVar();
  PopStyleVar();
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

void ImGui::Container(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags flags /*= 0*/){

  bool collapsed = true;

  Dock *dd = dockht[label];
  if (!dd){
    dd = dockht[label] = new Dock;
    IM_ASSERT(dd);
  }

  // If the container has a window docked, set the minimum size
  // printf("%f %f\n",dd->window->SizeContents.x,dd->window->SizeContents.y);
  if (dd->currenttab){
    ImGuiContext *g = GetCurrentContext();
    ImVec2 size = dd->currenttab->window->SizeContents + dd->currenttab->window->WindowPadding;
    size.y += GetTextLineHeightWithSpacing() + dd->currenttab->window->WindowPadding.y +
      dd->window->TitleBarHeight();
    SetNextWindowSizeConstraints(size,ImVec2(FLT_MAX,FLT_MAX),nullptr,nullptr);
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

  // If the container has just been closed, detach all docked windows
  if (dd->status == Dock::Status_Closed && !dd->stack.empty()){
    dd->clearContainer();
  }

  // Draw the container elements
  dd->drawContainer();

  // If the container is clicked, set the correct hovered/moved flags 
  // and raise container & docked window to the top of the stack.
  if (dd->currenttab){
    ImGuiContext *g = GetCurrentContext();
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

  // If this window is docked, put it right on top of its container.
  // Transfer any unhandled clicks to the underlying container.
  if (dd->status == Dock::Status_Docked && !dd->hidden){
    if ((g->HoveredWindow == dd->window || g->HoveredWindow == dd->parent->window) &&
	IsMouseHoveringRect(dd->window->Pos,dd->window->Pos+dd->window->Size,false) &&
	g->IO.MouseClicked[0]){
      dd->parent->SetContainerHoveredMovedActive(!IsItemClicked(0));
    }
    dd->parent->RaiseCurrentTab();
  }

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
    // Text("type=%d status=%d\n", dock.second->type, dock.second->status);
    // Text("sttype=%d list_size=%d\n", dock.second->sttype, dock.second->stack.size());
    // if (dock.second->p_open)
    //   Text("p_open=%d\n", *(dock.second->p_open));
    // Text("stackpos=%d\n",dock.second->stackpos);
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

