// -*-c++-*-
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

#ifndef IMGUI_DOCK_H
#define IMGUI_DOCK_H

#include "imgui.h"
#include "imgui_internal.h"
#include <list>

using namespace std;

namespace ImGui{

  struct Dock{
    enum Type_{Type_None,Type_Container,Type_Dock};
    enum Status_{Status_None,Status_Open,Status_Collapsed,Status_Closed,
		 Status_Dragged,Status_Docked,Status_Lifted};
    enum Stack_{Stack_None,Stack_Vertical,Stack_Horizontal,Stack_Leaf};

    ImU32 id = 0; // id of the window
    char* label = nullptr; // window label
    ImVec2 pos; // position of the window
    ImVec2 size; // size of the window
    ImVec2 size_saved; // saved size (before docking for dockable window, of the titlebar+tab for container)
    ImRect tabbarrect = {}; // rectangle for the container tab bar
    ImVector<float> tabsx = {};
    ImGuiWindowFlags flags = 0; // flags for the window
    bool hidden = false; // whether a docked window is hidden
    bool collapsed = false; // whether a docked window is collapsed
    bool collapsed_saved = false; // saved collapsed (before docking)
    Type_ type = Type_None; // type of docking window
    Status_ status = Status_None; // status of the docking window
    Stack_ sttype = Stack_None; // stacking type when docked
    list<Dock *> stack = {}; // stack of docks at this level
    Dock *currenttab = nullptr; // currently selected tab
    Dock *parent = nullptr; // node immediately above current node
    Dock *root = nullptr; // root node of this tree
    ImGuiWindow* window = nullptr; // associated window
    bool *p_open = nullptr; // the calling routine open window bool
    bool control_window_this_frame = false; // the pos, size, etc. change window's attributes this frame

    Dock(): pos(0,0), size(-1,-1), size_saved(-1,-1) {};
    ~Dock(){ MemFree(label);}

    // Is the mouse hovering the tab bar of this dock?
    bool IsMouseHoveringTabBar();
    // Get the nearest tab border in the tab
    int getNearestTabBorder();
    // Show the drop targets for this window (full window)
    void showDropTargetFull();
    // Show the drop targets for this window (on the tab bar)
    void showDropTargetOnTabBar();
    // Add a new dock to a container (at position ithis)
    void newDock(Dock *dnew, int ithis = -1);
    // Draw container
    void drawContainer(bool allowresize);
    // Clear all docked windows from a container
    void clearContainer();
    // Draw the tab bar of a tabbed container
    void drawTabBar();
    // Set the hovered window, moved window, and active ID for a container/tab pair
    void SetContainerHoveredMovedActive(bool setid);
    // Raise the current displayed tab in a container in the window stack to 
    // just before the container
    void RaiseCurrentTab();
    // Get the position of this dock's window in the window stack
    int getWindowStackPosition();

  }; // struct Dock

  // Create a container with the given label. If p_open, with a close
  // button.  Extra window flags are passed to the container window.
  void Container(const char* label, bool* p_open=nullptr, ImGuiWindowFlags extra_flags=0);

  // Create/end a window that can be docked to a container. If p_open,
  // with a close button.  Extra window flags are passed to the
  // container window.
  bool BeginDock(const char* label, bool* p_open=nullptr, ImGuiWindowFlags extra_flags=0);
  void EndDock();

  // Free the memory occupied by the docks and the dock hash table
  void ShutdownDock();

  // Print information about the current known docks. Debug only.
  void Print();
  
} // namespace ImGui
#endif
