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

//   xx rootcontainer xx
// fixing hovering over debug and over a non-dock window (visible)
// bring back the rootcontainer title and floating
// problem with several containers at the same time and rootcontainer
// problem docking containers to rootcontainers
// bar movement
// do not draw the container background when docked.
// autoresize of the container - manual and in the public interface
// deal with the small container size problem
// focus the whole rootcontainer on click
// eliminate a bar by rightclicking on it
// take out a container from a rootcontainer
// pass container or dock to rootcontainer. build rootcontainer tree from code
// resizable rootcontainer
// maybe: rootcontainer with its own window & control of the window stack 
// problem docking empty containers to a rootcontainer
// horizontal and vertical containers have active perpendicular edges?
// undock external containers
//   xx end xx
// triangles in the tabs; overlap
// clean up all methods
// see docking thread in imgui github

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
static unordered_map<string,Dock*> dockht = {}; // global dock hash table (string key)
static unordered_map<ImGuiWindow*,Dock*> dockwin = {}; // global dock hash table (window key)
static Dock *FindHoveredContainer(int type = -1); // find the container hovered by the mouse
static void placeWindow(ImGuiWindow* base,ImGuiWindow* moved,int idelta); // place a window relative to another in the window stack

//xx// Dock context methods //xx//

static Dock *FindHoveredDock(int type){
  ImGuiContext *g = GetCurrentContext();
  for (int i = g->Windows.Size-1; i >= 0; i--){
    ImGuiWindow *window = g->Windows[i];
    if (window->Flags & ImGuiWindowFlags_NoInputs)
      continue;
    if (window->Flags & ImGuiWindowFlags_ChildWindow)
      continue;
    ImRect bb(window->WindowRectClipped.Min - g->Style.TouchExtraPadding, window->WindowRectClipped.Max + g->Style.TouchExtraPadding);
    if (!bb.Contains(g->IO.MousePos))
      continue;
    Dock *dock = dockwin[window];
    if (!dock)
      continue;
    if (dock->collapsed)
      return nullptr;
    if (!dock->hoverable || dock->hidden)
      continue;
    if (type >= 0 && dock->type != type)
      return nullptr;
    return dock;
  }
  return nullptr;
}

static void placeWindow(ImGuiWindow* base,ImGuiWindow* moved,int idelta){
  ImGuiContext *g = GetCurrentContext();

  if (!base || !moved) return;

  int ibase = -1, imoved = -1;
  for (int i = 0; i < g->Windows.Size; i++){
    if (g->Windows[i] == base)
      ibase = i;
    if (g->Windows[i] == moved)
      imoved = i;
  }
  if (ibase < 0 || imoved < 0 || imoved == ibase + idelta) return;
  g->Windows.erase(g->Windows.begin() + imoved);
  if (imoved < ibase) ibase--;
  if (ibase + idelta > g->Windows.size())
    g->Windows.push_back(moved);
  else
    g->Windows.insert(g->Windows.begin() + max(ibase + idelta,0),moved);
}

//xx// Dock methods //xx//

bool Dock::IsMouseHoveringTabBar(){
  const ImVec2 ytabcushiondn = ImVec2(0.f,10.);
  const ImVec2 ytabcushionup = ImVec2(0.f,this->status==Dock::Status_Docked?0.:10.);
  return !this->stack.empty() && IsMouseHoveringRect(this->tabbarrect.Min-ytabcushionup,this->tabbarrect.Max+ytabcushiondn,false);
}

int Dock::IsMouseHoveringEdge(){
  // top, right, bottom, left
  const ImVec2 x0[8] = {{0.0,0.0}, {0.9,0.0}, {0.0,0.9}, {0.0,0.0}};
  const ImVec2 x1[8] = {{1.0,0.1}, {1.0,1.0}, {1.0,1.0}, {0.1,1.0}};
  const float minsize = 20.f;

  ImVec2 xmin, xmax;
  for (int i=0; i<4; i++){
    xmin.x = this->pos.x + x0[i].x * this->size.x;
    xmax.x = this->pos.x + x1[i].x * this->size.x;
    xmin.y = this->pos.y + x0[i].y * this->size.y;
    xmax.y = this->pos.y + x1[i].y * this->size.y;
    if (i == 0 && (xmax.y-xmin.y) < minsize)
      xmax.y = xmin.y + min(minsize,0.5f * this->size.y);
    else if (i == 2 && (xmax.y-xmin.y) < minsize)
      xmin.y = xmax.y - min(minsize,0.5f * this->size.y);
    else if (i == 3 && (xmax.x-xmin.x) < minsize)
      xmax.x = xmin.x + min(minsize,0.5f * this->size.x);
    else if (i == 1 && (xmax.x-xmin.x) < minsize)
      xmin.x = xmax.x - min(minsize,0.5f * this->size.x);

    if (IsMouseHoveringRect(xmin,xmax,false))
      return i+1;
  }
  return 0;
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
  Begin("##Drop",nullptr,ImGuiWindowFlags_Tooltip|ImGuiWindowFlags_NoTitleBar|
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
  Begin("##Drop",nullptr,ImGuiWindowFlags_Tooltip|ImGuiWindowFlags_NoTitleBar|
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

void Dock::showDropTargetEdge(int edge){
  // top, right, bottom, left
  const ImVec2 x0[4] = {{0.0,0.0}, {0.5,0.0}, {0.0,0.5}, {0.0,0.0}};
  const ImVec2 x1[4] = {{1.0,0.5}, {1.0,1.0}, {1.0,1.0}, {0.5,1.0}};

  if (edge > 0){
    ImVec2 xmin, xmax;
    xmin.x = this->pos.x + x0[edge-1].x * this->size.x;
    xmax.x = this->pos.x + x1[edge-1].x * this->size.x;
    xmin.y = this->pos.y + x0[edge-1].y * this->size.y;
    xmax.y = this->pos.y + x1[edge-1].y * this->size.y;

    SetNextWindowSize(ImVec2(0,0));
    Begin("##Drop",nullptr,ImGuiWindowFlags_Tooltip|ImGuiWindowFlags_NoTitleBar|
    	  ImGuiWindowFlags_NoInputs|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_AlwaysAutoResize);
    ImDrawList* drawl = GetWindowDrawList();
    drawl->PushClipRectFullScreen();
    ImU32 docked_color = GetColorU32(ImGuiCol_FrameBg);
    docked_color = (docked_color & 0x00ffFFFF) | 0x80000000;
    drawl->AddRectFilled(xmin,xmax,docked_color);
    drawl->PopClipRect();
    End();
  }
}

void Dock::newDock(Dock *dnew, int ithis /*=-1*/){
  if (!(this->type == Dock::Type_Container)) return;

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

void Dock::newDockRoot(Dock *dnew, int iedge){
  // 1:top, 2:right, 3:bottom, 4:left
  if (iedge == 0) return;

  Dock *dcont = nullptr;
  if (dnew->type == Dock::Type_Container)
    dcont = dnew;
  if (this->parent->type == Dock::Type_Root){
    Type_ type; 
    if (iedge == 1 || iedge == 3)
      type = Dock::Type_Horizontal;
    else
      type = Dock::Type_Vertical;
    dcont = this->OpRoot_ReplaceHV(type,iedge==1||iedge==4,dcont);
  } else if (this->parent->type == Dock::Type_Horizontal){
    if (iedge == 1 || iedge == 3){
      dcont = this->OpRoot_AddToHV(iedge==1,dcont);
    } else {
      dcont = this->OpRoot_ReplaceHV(Dock::Type_Vertical,iedge==4,dcont);
    }
  } else if (this->parent->type == Dock::Type_Vertical){
    if (iedge == 2 || iedge == 4){
      dcont = this->OpRoot_AddToHV(iedge==4,dcont);
    } else {
      dcont = this->OpRoot_ReplaceHV(Dock::Type_Horizontal,iedge==1,dcont);
    }
  }
  if (dnew->type != Dock::Type_Container)
    dcont->newDock(dnew);
}

Dock *Dock::OpRoot_ReplaceHV(Dock::Type_ type,bool before,Dock *dcont/*=nullptr*/){
  // 1:top, 2:right, 3:bottom, 4:left, 5:topleft, 6:topright, 7:bottomright, 8:bottomleft
  Dock *dpar = this->parent;
  Dock *root = dpar->root;
  if (!dcont){
    // new empty container
    char label1[strlen(root->label)+15];
    ImFormatString(label1,IM_ARRAYSIZE(label1),"%s__%d__",root->label,++(root->nchild));
    dcont = new Dock;
    IM_ASSERT(dcont);
    dcont->label = ImStrdup(label1);
    dockht[string(dcont->label)] = dcont;
    dcont->type = Dock::Type_Container;
    dcont->id = ImHash(label1,0);
    dcont->status == Dock::Status_Docked;
    dcont->automatic = true;
    dcont->hoverable = true;
  }

  // new horizontal or vertical container
  char label2[strlen(root->label)+15];
  ImFormatString(label2,IM_ARRAYSIZE(label2),"%s__%d__",root->label,++(root->nchild));
  Dock *dhv = new Dock;
  IM_ASSERT(dhv);
  dhv->label = ImStrdup(label2);
  dockht[string(dhv->label)] = dhv;
  dhv->type = type;
  dhv->id = ImHash(label2,0);
  dhv->status == Dock::Status_Docked;
  dhv->hoverable = false;
  dhv->automatic = true;

  // build the new horizontal/vertical
  if (before){
    dhv->stack.push_back(dcont);
    dhv->stack.push_back(this);
  } else {
    dhv->stack.push_back(this);
    dhv->stack.push_back(dcont);
  }

  // replace this with the new horizontal/vertical container in the parent's stack
  for(auto it = dpar->stack.begin(); it != dpar->stack.end(); it++){ 
    if (*it == this){
      dpar->stack.insert(dpar->stack.erase(it),dhv);
      break;
    }
  }

  // return the new container
  return dcont;
}

Dock *Dock::OpRoot_AddToHV(bool before,Dock *dcont/*=nullptr*/){
  // 1:top, 2:right, 3:bottom, 4:left, 5:topleft, 6:topright, 7:bottomright, 8:bottomleft
  Dock *dpar = this->parent;
  Dock *root = dpar->root;
  if (!dcont){
    // new empty container
    char label1[strlen(root->label)+15];
    ImFormatString(label1,IM_ARRAYSIZE(label1),"%s__%d__",root->label,++(root->nchild));
    dcont = new Dock;
    IM_ASSERT(dcont);
    dcont->label = ImStrdup(label1);
    dockht[string(dcont->label)] = dcont;
    dcont->type = Dock::Type_Container;
    dcont->id = ImHash(label1,0);
    dcont->status == Dock::Status_Docked;
    dcont->hoverable = true;
    dcont->automatic = true;
  }

  // add to the parent's stack
  for(auto it = dpar->stack.begin(); it != dpar->stack.end(); it++){ 
    if (*it == this){
      if (before)
	dpar->stack.insert(it,dcont);
      else
	dpar->stack.insert(++it,dcont);
      break;
    }
  }

  // return the new container
  return dcont;
}

void Dock::OpRoot_FillEmpty(){
  if (!this->stack.empty()) return;

  char tmp[strlen(this->label)+15];
  ImFormatString(tmp,IM_ARRAYSIZE(tmp),"%s__%d__",this->label,++(this->nchild));
  Dock *dcont = new Dock;
  IM_ASSERT(dcont);
  dcont->label = ImStrdup(tmp);
  dockht[string(dcont->label)] = dcont;
  dcont->type = Dock::Type_Container;
  dcont->id = ImHash(tmp,0);
  dcont->status == Dock::Status_Docked;
  dcont->hoverable = true;
  dcont->automatic = true;
  this->stack.push_back(dcont);
}

void Dock::killDock(Dock *parent/*=nullptr*/, Dock *replacement/*=nullptr*/){
  if (parent){
    for(auto it = parent->stack.begin(); it != parent->stack.end(); it++){ 
      if (*it == this){
	if (!replacement)
	  parent->stack.erase(it);
	else
	  parent->stack.insert(parent->stack.erase(it),replacement);
	break;
      }
    }
  }
  dockht.erase(string(this->label));
  dockwin.erase(this->window);
  if (this) delete this;
}

void Dock::drawContainer(bool noresize){
  if (!(this->type == Dock::Type_Container)) return;
    
  if (this->stack.size() > 0){
    // Draw the tab
    this->drawTabBar();

    // Hide all tabs
    for (auto dd : this->stack) 
      dd->hideTabWindow(this);

    // Unhide current tab
    if (!this->hidden && !this->collapsed && this->currenttab)
      this->currenttab->showTabWindow(this,noresize);

  } // if (this->stack.size() > 0);
}

void Dock::drawRootContainer(Dock *root){
  const float barwidth = 8.;

  this->root = root;
  if (this->type == Dock::Type_Root){
    for (auto dd : this->stack){
      dd->pos = this->pos;
      dd->size = this->size;
      dd->parent = this;
      dd->drawRootContainer(root);
    }
  } else if (this->type == Dock::Type_Horizontal || this->type == Dock::Type_Vertical) {
    int n = -1;
    int ntot = this->stack.size();
    float width1, width2;
    for (auto dd : this->stack){
      n++;
      if (n == 0){
	width1 = 0.f;
	width2 = 0.5f * barwidth;
      } else if (n == this->stack.size()-1){
	width1 = 0.5f * barwidth;
	width2 = 0.5f * barwidth;
      } else {
	width1 = 0.5f * barwidth;
	width2 = 1.0f * barwidth;
      }
      dd->pos = this->pos;
      if (this->type == Dock::Type_Horizontal)
	dd->pos.y += (n * this->size.y) / ntot + width1;
      else
	dd->pos.x += (n * this->size.x) / ntot + width1;
      dd->size = this->size;
      if (this->type == Dock::Type_Horizontal){
	dd->size.y /= ntot;
	dd->size.y -= width2;
      }
      else{
	dd->size.x /= ntot;
	dd->size.x -= width2;
      }
      dd->parent = this;
      dd->drawRootContainer(root);
    }
  } else if (this->type == Dock::Type_Container) {
    // Draw the docked container window
    bool transparentframe = this->currenttab;
    SetNextWindowPos(this->pos);
    SetNextWindowSize(this->size);
    this->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoResize|
	ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoCollapse|
	ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoBringToFrontOnFocus;
    if (transparentframe)
      PushStyleColor(ImGuiCol_WindowBg,TransparentColor(ImGuiCol_WindowBg));
    Begin(this->label,nullptr,this->flags);
    this->window = GetCurrentWindow();
    dockwin[this->window] = this;
    this->status = Dock::Status_Docked;
    this->hoverable = true;
    this->drawContainer(true);
    this->size_saved.y = this->tabbarrect.Max.y - this->pos.y;
    End();
    if (transparentframe)
      PopStyleColor();

    // rootcontainer -> container -> dock
    placeWindow(this->root->window,this->window,+1);
    if (this->currenttab)
      placeWindow(this->window,this->currenttab->window,+1);
  }
}

void Dock::drawRootContainerBars(Dock *root){
  const float barwidth = 6.;

  this->root = root;
  if (this->type == Dock::Type_Root){
    for (auto dd : this->stack)
      dd->drawRootContainerBars(root);
  } else if (this->type == Dock::Type_Horizontal || this->type == Dock::Type_Vertical) {
    int n = -1;
    int ntot = this->stack.size();

    ImVec2 pos, size, limits;
    int direction;

    for (auto dd : this->stack){
      n++;
      if (n != 0){
	pos = this->pos;
	size = this->size;
	if (this->type == Dock::Type_Horizontal){
	  pos.y += (n * this->size.y) / ntot - 0.5f * barwidth;
	  size.y = barwidth;
	  direction = 2;
	} else {
	  pos.x += (n * this->size.x) / ntot - 0.5f * barwidth;
	  size.x = barwidth;
	  direction = 1;
	}
	limits = ImVec2(pos.x,pos.x+size.x-barwidth);
	SlidingBar(root->window, &pos, size, limits, direction, 2.0f);
      }
      dd->drawRootContainerBars(root);
    }
  }
}

void Dock::drawTabBar(){
  ImGuiContext *g = GetCurrentContext();
  const float tabheight = GetTextLineHeightWithSpacing();
  const float barheight = 2 * GetTextLineHeightWithSpacing();
  const float maxtabwidth = 100.;
  ImVec4 text_color = g->Style.Colors[ImGuiCol_Text];
  text_color.w = 2.0 / g->Style.Alpha;
  bool raise = false;

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
  char tmp[strlen(this->label)+15];
  ImFormatString(tmp,IM_ARRAYSIZE(tmp),"%s__tab__",this->label);
  if (BeginChild(tmp,ImVec2(this->size.x,barheight),false)){
    bool active = false, hovered = false;
    Dock *dderase = nullptr, *ddlast = nullptr;
    ImVec2 center, pos0, pos1, pos1s, text_size;
    for (auto dd : this->stack) {
      SameLine();
      // make the x-button, update the container info
      bool dragged, dclicked, closeclicked;
      if (ButtonWithX(dd->label, ImVec2(tabwidth_long, tabheight), (dd == this->currenttab), false,
                      dd->p_open, &dragged, &dclicked, &closeclicked, 2.f/g->Style.Alpha)){
        this->currenttab = dd;
	dd->parent = this;
	dd->root = this->root;
      }
      this->tabsx.push_back(GetItemRectMax().x-tabwidth_long);
    
      // lift the tab using the main button
      if (dragged){
        dd->status = Dock::Status_Dragged;
	dd->hoverable = false;
        goto lift_this_tab;
      }
      // double click detaches the tab
      if (dclicked){
        dd->status = Dock::Status_Open;
	dd->hoverable = true;
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
      if (!this->currenttab && this->stack.size() > 0){
        this->currenttab = this->stack.front();
	this->currenttab->parent = this;
      }
    }
    // last item in the tabsx
    this->tabsx.push_back(GetItemRectMax().x);
  } // BeginChild(tmp, ImVec2(this->size.x,barheight), true)
  this->tabbarrect = GetCurrentWindowRead()->DC.LastItemRect;
  this->tabbarrect.Min.x = this->pos.x;
  EndChild();
  PopStyleVar();
  PopStyleVar();
  PopStyleColor();
}

void Dock::clearContainer(){
  const float increment = 50.;

  ImVec2 pos = this->pos;
  for (auto dd : this->stack) {
    dd->status = Dock::Status_Open;
    dd->hoverable = true;
    dd->control_window_this_frame = true;
    dd->size = dd->size_saved;
    dd->collapsed = dd->collapsed_saved;
    dd->pos = pos;
    pos = pos + ImVec2(increment,increment);
  }
  this->currenttab = nullptr;
  this->stack.clear();
}

void Dock::focusContainer(){
  ImGuiContext *g = GetCurrentContext();

  // Push the container and the docked window to the top of the stack
  for (int i = 0; i < g->Windows.Size; i++){
    if (g->Windows[i] == this->window)
      g->Windows.erase(g->Windows.begin() + i);
  }
  g->Windows.push_back(this->window);
  if (this->currenttab){
    for (int i = 0; i < g->Windows.Size; i++){
      if (g->Windows[i] == this->currenttab->window)
	g->Windows.erase(g->Windows.begin() + i);
    }
    g->Windows.push_back(this->currenttab->window);
  }

  // The docked window becomes focused, if possible. Otherwise, the container
  if (!this->currenttab){
    g->HoveredRootWindow = this->window;
    g->HoveredWindow = this->window;
  } else {
    g->HoveredRootWindow = this->currenttab->window;
    g->HoveredWindow = this->currenttab->window;
    if (g->ActiveIdWindow == this->currenttab->window && !g->ActiveIdIsAlive)
      ClearActiveID();
  }

  // The container is being moved
  if (!IsAnyItemActive() && !IsAnyItemHovered() && g->IO.MouseClicked[0]){
    g->MovedWindow = this->window;
    g->MovedWindowMoveId = this->window->RootWindow->MoveId;
    if (this->currenttab)
      SetActiveID(g->MovedWindowMoveId, this->currenttab->window->RootWindow);
  }
}

void Dock::killContainerMaybe(){
  // Only kill containers, horziontals, and verticals that were automatically generated
  if (!this || (this->type != Dock::Type_Container && this->type != Dock::Type_Horizontal && 
		this->type != Dock::Type_Vertical) || !(this->automatic))
    return;
  Dock *dpar = this->parent;
  if (!dpar) return;

  if (this->type == Dock::Type_Container && this->stack.empty()){
    // An empty container
    // Do not remove the last container from a root container, even if it's empty
    if (dpar->type == Dock::Type_Root && dpar->stack.size() <= 1) return;
    this->killDock(dpar);

    // Try to kill its parent
    dpar->killContainerMaybe();
  } else if ((this->type == Dock::Type_Vertical || this->type == Dock::Type_Horizontal) && this->stack.empty()){
    // If a horizontal or vertical is empty, turn it into a container
    this->type = Dock::Type_Container;
    // then to try kill it
    this->killContainerMaybe();
  } else if ((this->type == Dock::Type_Vertical || this->type == Dock::Type_Horizontal) && this->stack.size() == 1){
    // This vertical/horizontal container only has window -> eliminate
    // it and connect the single window to its parent
    this->killDock(dpar,this->stack.back());

    // Try to kill its parent
    dpar->killContainerMaybe();
  }
}

void Dock::hideTabWindow(Dock *dcont){
  if (!dcont) return;
  this->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoScrollbar|
    ImGuiWindowFlags_NoScrollWithMouse|ImGuiWindowFlags_NoCollapse|
    ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoInputs;
  this->pos = ImVec2(0.,0.);
  this->size = ImVec2(0.,0.);
  this->hidden = true;
}

void Dock::showTabWindow(Dock *dcont, bool noresize){
  if (!dcont) return;

  float topheight = dcont->tabbarrect.Max.y - dcont->pos.y;
  this->pos = dcont->pos;
  this->pos.y += topheight;
  this->size = dcont->size;
  this->size.y = dcont->size.y - topheight;
  this->hidden = false;
  this->flags = ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|
    ImGuiWindowFlags_NoCollapse|ImGuiWindowFlags_NoSavedSettings|
    ImGuiWindowFlags_NoBringToFrontOnFocus|ImGuiWindowFlags_HorizontalScrollbar;
  if (noresize)
    this->flags = this->flags | ImGuiWindowFlags_NoResize;
}

//xx// Public interface //xx//

Dock *ImGui::RootContainer(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags extra_flags /*= 0*/){
  bool collapsed;
  ImGuiContext *g = GetCurrentContext();
  ImGuiWindowFlags flags = extra_flags;

  Dock *dd = dockht[string(label)];
  if (!dd){
    dd = new Dock;
    IM_ASSERT(dd);
    dd->label = ImStrdup(label);
    dockht[string(dd->label)] = dd;
    dd->type = Dock::Type_Root;
    dd->id = ImHash(label,0);
    dd->nchild = 0;
  }

  // Initialize with a container if empty
  dd->OpRoot_FillEmpty();

  // xxxx // set minimum size? // xxxx //

  // Making an invisible window
  PushStyleColor(ImGuiCol_WindowBg,TransparentColor(ImGuiCol_WindowBg));
  flags = flags | ImGuiWindowFlags_NoResize;
  collapsed = !Begin(label,p_open,flags);

  // set the properties of the rootcontainer window
  dd->window = GetCurrentWindow();
  dd->pos = dd->window->Pos + ImVec2(0.f,dd->window->TitleBarHeight());
  dd->size = dd->window->Size - ImVec2(0.f,dd->window->TitleBarHeight());
  dd->size_saved = ImVec2(0.,0.);
  dd->type = Dock::Type_Root;
  dd->flags = extra_flags;
  dd->root = dd;
  dd->collapsed = collapsed;
  dockwin[dd->window] = dd;
  dd->p_open = p_open;
  dd->hoverable = false;

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

  // xxxx If the root container has just been closed, detach all docked windows

  // Traverse the tree and draw all the bars
  dd->drawRootContainerBars(dd);

  // End the root container window
  End();
  PopStyleColor();

  // Traverse the tree and draw all the containers
  dd->drawRootContainer(dd);

  return dd;
}

Dock *ImGui::Container(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags extra_flags /*= 0*/){

  bool collapsed = true;
  ImGuiContext *g = GetCurrentContext();
  ImGuiWindowFlags flags = extra_flags;

  Dock *dd = dockht[string(label)];
  if (!dd){
    dd = new Dock;
    IM_ASSERT(dd);
    dd->label = ImStrdup(label);
    dockht[string(dd->label)] = dd;
    dd->type = Dock::Type_Container;
    dd->id = ImHash(label,0);
  }

  // If docked, the root container takes care of everything
  if (dd->status == Dock::Status_Docked) 
    return dd;

  // If the container has a window docked, set the minimum size.
  if (dd->currenttab){
    SetNextWindowSizeConstraints(dd->size_saved,ImVec2(FLT_MAX,FLT_MAX),nullptr,nullptr);
    flags = flags | ImGuiWindowFlags_NoResize;
  } else {
    SetNextWindowSizeConstraints(ImVec2(0.,4*(GetTextLineHeightWithSpacing()+g->Style.ItemSpacing.y)),
					ImVec2(FLT_MAX,FLT_MAX),nullptr,nullptr);
  }

  // Render any container widgets in here
  bool transparentframe = dd->currenttab;
  if (transparentframe)
    PushStyleColor(ImGuiCol_WindowBg,TransparentColor(ImGuiCol_WindowBg));
  collapsed = !Begin(label,p_open,flags);

  // Fill the info for this dock
  dd->pos = GetWindowPos();
  dd->size = GetWindowSize();
  dd->size_saved = ImVec2(0.,0.);
  dd->type = Dock::Type_Container;
  dd->flags = extra_flags;
  dd->root = dd;
  dd->collapsed = collapsed;
  dd->window = GetCurrentWindow();
  dockwin[dd->window] = dd;
  dd->p_open = p_open;

  // Update the status
  Dock *ddest;
  if (g->ActiveId == GetCurrentWindow()->MoveId && g->IO.MouseDown[0]){
    // Dragging
    dd->status = Dock::Status_Dragged;
    dd->hoverable = false;
    ddest = FindHoveredDock(Dock::Type_Container);
  } else {
    int ithis = -1, iedge = 0;
    bool dropit = (dd->status == Dock::Status_Dragged && (ddest = FindHoveredDock(Dock::Type_Container)));
    if (dropit && ddest->stack.empty() && ddest->parent && ddest->parent->type == Dock::Type_Root){
      // drop it into the root container and replace it
      ddest->killDock(ddest->parent,dd);
      dd->status = Dock::Status_Docked;
      dd->hoverable = true;
    } else if (dropit && ddest->status == Dock::Status_Docked && ((iedge = ddest->IsMouseHoveringEdge()) > 0)){
      // drop into the edge
      dd->status = Dock::Status_Docked;
      dd->hoverable = true;
      ddest->newDockRoot(dd,iedge);
    } else {
      // Stationary -> open, closed, or collapsed
      if (!p_open || *p_open){
	if (collapsed){
	  dd->status = Dock::Status_Collapsed;
	  dd->hoverable = true;
	}
	else{
	  dd->status = Dock::Status_Open;
	  dd->hoverable = true;
	}
      } else {
	dd->status = Dock::Status_Closed;
	dd->hoverable = false;
      }
    }
  }

  // If dragged and hovering over a container, show the drop rectangles
  if (dd->status == Dock::Status_Dragged){
    if (ddest){
      if (ddest->stack.empty() && ddest->parent && ddest->parent->type == Dock::Type_Root)
        ddest->showDropTargetFull();
      else if (ddest->status == Dock::Status_Docked)
        ddest->showDropTargetEdge(ddest->IsMouseHoveringEdge());
    }
  }

  // If the container has just been closed, detach all docked windows
  if (dd->status == Dock::Status_Closed && !dd->stack.empty())
    dd->clearContainer();

  // Draw the container elements
  dd->drawContainer(extra_flags & ImGuiWindowFlags_NoResize);
  dd->size_saved.y = dd->tabbarrect.Max.y - dd->pos.y;

  // If the container is clicked, set the correct hovered/moved flags
  // and raise container & docked window to the top of the stack.
  if (dd->currenttab && g->IO.MouseClicked[0] && g->HoveredRootWindow == dd->window)
      dd->focusContainer();

  // Put the current tab on top of the current window
  if (dd->currenttab)
    placeWindow(dd->window,dd->currenttab->window,+1);

  End();
  if (transparentframe)
    PopStyleColor();

  return dd;
}

bool ImGui::BeginDock(const char* label, bool* p_open /*=nullptr*/, ImGuiWindowFlags flags /*= 0*/, Dock* oncedock /*=nullptr*/){
  bool collapsed;
  ImGuiContext *g = GetCurrentContext();

  // Create the entry in the dock context if it doesn't exist
  Dock *dd = dockht[string(label)];
  if (!dd) {
    dd = new Dock;
    IM_ASSERT(dd);
    dd->label = ImStrdup(label);
    dockht[string(dd->label)] = dd;
    dd->type = Dock::Type_Dock;
    dd->id = ImHash(label,0);

    // This is the first pass -> if oncedock exists, dock to that container
    if (oncedock){
      oncedock->newDock(dd,-1);
      dd->status = Dock::Status_Docked;
      dd->hoverable = false;
      dd->control_window_this_frame = true;
      dd->showTabWindow(oncedock,dd->flags & ImGuiWindowFlags_NoResize);
      dd->parent = oncedock;
    }
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
        Begin(label,nullptr,dd->size,0.0,flags);
      } else {
        Begin(label,nullptr,flags);
      }
    } else if (dd->status == Dock::Status_Dragged) { 
      // the window has just been lifted from a container. Go back to
      // being a normal window with the new position and size; being
      // dragged.
      collapsed = !Begin(label,p_open,flags);
      dd->window = GetCurrentWindow();
      dockwin[dd->window] = dd;
      FocusWindow(dd->window);
      g->MovedWindow = dd->window;
      g->MovedWindowMoveId = dd->window->RootWindow->MoveId;
      SetActiveID(g->MovedWindowMoveId, dd->window->RootWindow);
      dd->parent->killContainerMaybe();
      dd->parent = nullptr;
    } else { 
      // the window has just been lifted, but not dragging
      collapsed = !Begin(label,p_open,flags);
      dd->window = GetCurrentWindow();
      dockwin[dd->window] = dd;
      FocusWindow(dd->window);
      dd->parent->killContainerMaybe();
      dd->parent = nullptr;
    }
  } else {
    // Floating window
    collapsed = !Begin(label,p_open,flags);
    dd->collapsed_saved = collapsed;
    if (!collapsed) dd->size_saved = dd->size;
  }

  // Fill the info for this dock
  dd->pos = GetWindowPos();
  dd->size = GetWindowSize();
  dd->type = Dock::Type_Dock;
  dd->flags = flags;
  dd->collapsed = collapsed;
  dd->window = GetCurrentWindow();
  dockwin[dd->window] = dd;
  dd->p_open = p_open;
  dd->control_window_this_frame = false;

  // Update the status
  Dock *ddest;
  if (g->ActiveId == GetCurrentWindow()->MoveId && g->IO.MouseDown[0]){
    // Dragging
    dd->status = Dock::Status_Dragged;
    dd->hoverable = false;
    ddest = FindHoveredDock(Dock::Type_Container);
  } else {
    int ithis = -1, iedge = 0;
    bool dropit = (dd->status == Dock::Status_Dragged && (ddest = FindHoveredDock(Dock::Type_Container)));
    if (dropit && (ddest->stack.empty() || ((ithis = ddest->getNearestTabBorder()) >= 0))){
      // Just stopped dragging and there is a container below
      dd->status = Dock::Status_Docked;
      dd->hoverable = false;
      ddest->newDock(dd,ithis);
    } else if (dropit && ddest->status == Dock::Status_Docked && ((iedge = ddest->IsMouseHoveringEdge()) > 0)){
      // stopped dragging and there is a root container below
      dd->status = Dock::Status_Docked;
      dd->hoverable = false;
      ddest->newDockRoot(dd,iedge);
    } else if (dd->status != Dock::Status_Docked){
      // Stationary -> open, closed, or collapsed
      if (!p_open || *p_open){
        if (collapsed){
          dd->status = Dock::Status_Collapsed;
	  dd->hoverable = true;
        }
        else{
          dd->status = Dock::Status_Open;
	  dd->hoverable = true;
        }
      } else {
        dd->status= Dock::Status_Closed;
	dd->hoverable = false;
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
      else if (ddest->status == Dock::Status_Docked)
        ddest->showDropTargetEdge(ddest->IsMouseHoveringEdge());
    }
  }

  // If this window is being clicked or dragged, focus the container+dock
  if (dd->status == Dock::Status_Docked && !dd->hidden)
    if (g->HoveredWindow == dd->window && g->IO.MouseClicked[0])
      dd->parent->focusContainer();

  if (dd->status == Dock::Status_Docked && !dd->hidden){
    // Move the container right under this dock
    placeWindow(dd->window,dd->parent->window,-1);
    // If the resize grip is being used, resize the container too.
    if (g->ActiveId == dd->window->GetID("#RESIZE"))
      dd->parent->window->SizeFull = dd->window->Size + dd->parent->size_saved;
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

  // for (auto dock : dockht){
  //   Text("key=%s id=%d label=%s\n", dock.first.c_str(),dock.second->id,dock.second->label);
  //   Text("pos=(%f,%f) size=(%f,%f)\n",dock.second->pos.x,dock.second->pos.y,dock.second->size.x,dock.second->size.y);
  //   Text("type=%d status=%d list_size=%d\n", dock.second->type, dock.second->status, dock.second->stack.size());
  //   if (dock.second->p_open)
  //     Text("p_open=%d\n", *(dock.second->p_open));
  // }

  ImGuiContext *g = GetCurrentContext();
  for (int i = 0; i < g->Windows.Size; i++){
    Text("%d %s %p\n",i,g->Windows[i]->Name,g->Windows[i]);
  }

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
  dockwin.clear();
}
