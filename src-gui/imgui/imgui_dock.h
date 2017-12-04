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
// Original code by vassvik (?) released as public domain.

// This code provides a collection of docking windows and containers
// for the immediate-mode graphical user interface (ImGui) library
// (https://github.com/ocornut/imgui) by Omar Cornut and
// collaborators. Three components are implemented:
//
// 1. Docks - Behave like normal windows, but can attach themselves to
// containers and root containers.
//
// 2. Containers - Windows that display the contents of docks attached
// to them. Containers have a tab bar at the top that shows all
// containers attached. At any time, one tab is active and the rest
// are hidden. Very similar to how a browser window (e.g. chrome,
// firefox) works.
//
// 3. Root containers - docks and containers can be attached to root
// container windows. The window space of a root container is split
// when a new window is added to it, depending on the position at
// which it is dropped. Containers inside root containers keep track
// of attached docks, and have an additional "lift grip" by which they
// can be pulled out of the root container.
//
// Besides the attaching/detaching, all three components behave like
// normal windows - they can be resized, auto-resized, moved,
// collapsed, and closed. The usual flags (NoResize,
// NoBringToFrontOnFocus,...) and placement commands
// (SetNextWindowPos,...) work as well. 
//
// The public interface contains the following:
//
// - The ImGui::Dock structure: docks, containers, and root containers
// are instances of this class.
//
// - ImGui::RootContainer: create a root container and return a
// pointer to the new Dock object.
//
// - ImGui::Container: create a container and return a pointer to the
// new Dock object. Neither root containers nor containers need to
// call an End() function (because no items are allowed inside them).
//
// - ImGui::BeginDock and ImGui::EndDock: the equivalent for docks to what
// ImGui::Begin() and ImGui::End() are for windows. They open and
// close a dock. If BeginDock returns true, interactive items can be
// added to the dock. In addition, the argument oncedock allows
// docking the dock to a container in the first pass. More convoluted
// ways to attach docks and containers on initialization are possible
// using the newDock and newDockRoot functions - an example will be
// provided elsewhere.
//
// - ImGui::GetCurrentDock: when used beteween BeginDock and EndDock,
// returns a pointer to the currently open dock. Otherwise, returns
// null.
// 
// - ImGui::ShutdownDock: deallocates memory for the dock hash
// table. Should be run once docks are no longer needed, or at the end
// of the program.
//
// Some notes:
//
// 1. Call RootContainer before any Container attached to it, and
// Container before any BeginDock/EndDock attached to it. Doing it the
// other way around works, but since root containers set the position
// of containers, and containers of docks, it will result in a lag of
// one or two frames when moving windows about.
//
// 2. Some widgets from imgui_widgets are used, so you will need that
// file and its header to use docks.
// 
// Have fun! -- Alberto

#ifndef IMGUI_DOCK_H
#define IMGUI_DOCK_H

#include "imgui.h"
#include "imgui_internal.h"
#include <list>

using namespace std;

namespace ImGui{

  typedef int DockFlags;

  struct Dock{
    // Enum for the type of docking windows
    enum Type_{Type_None,Type_Root,Type_Container,Type_Dock,Type_Horizontal,
               Type_Vertical};
    // Status of a single dock
    enum Status_{Status_None,Status_Open,Status_Collapsed,Status_Closed,
                 Status_Dragged,Status_Docked};
    // Dock flags
    enum DockFlags_{
      DockFlags_NoLiftContainer = 1 << 0, // A container is not allowed to be lifted
      DockFlags_Transparent = 1 << 1,     // This window is transparent (but still handles inputs)
    };

    char* label = nullptr; // dock and window label
    ImGuiWindow* window = nullptr; // associated window
    Type_ type = Type_None; // type of docking window
    Status_ status = Status_None; // status of the docking window
    ImVec2 pos = {}; // position of the window
    ImVec2 pos_saved = {}; // position of the window (before docking)
    ImVec2 size = {-1.f,-1.f}; // size of the window
    ImVec2 size_saved = {}; // saved size (before docking for dockable window)
    ImGuiWindowFlags flags = 0; // flags for the window
    ImGuiWindowFlags flags_saved = 0; // flags for the window (before docking)
    DockFlags dockflags = 0; // flags for this dock
    bool collapsed = false; // whether a docked window is collapsed
    bool collapsed_saved = false; // saved collapsed (before docking)
    ImRect tabbarrect = {}; // rectangle for the container tab bar
    float tabdz = 0.f; // z position for the end of the tab bar (container)
    ImGuiWindow* tabwin = nullptr; // pointer to the tab window (for cleaning up the window stack)
    ImVector<float> tabsx = {}; // tab positions for container; sliders for h/v-container
    ImVector<bool> tabsfixed = {}; // whether a slider is fixed.
    bool hidden = false; // whether a docked window is hidden
    bool hoverable = true; // whether a window responds to being hovered
    list<Dock *> stack = {}; // stack of docks at this level
    Dock *currenttab = nullptr; // currently selected tab (container)
    Dock *parent = nullptr; // immediate dock to which this is dock
    Dock *root = nullptr; // root container to which this is docked
    bool *p_open = nullptr; // the calling routine open window bool
    bool control_window_this_frame = false; // the pos, size, etc. change window's attributes this frame
    int nchild_ = 0; // number of children (to generate labels in rootcontainer)
    int nchild = 0; // number of children (to count for the last dock in rootcontainer)
    bool automatic = false; // whether this dock was automatically generated in a rootcontainer

    Dock(){};
    ~Dock(){ MemFree(label);}

    // Is the mouse hovering the tab bar of this dock? (no rectangle clipping)
    bool IsMouseHoveringTabBar();
    // Is the mouse hovering the drop edges of a container? (no rectangle clipping)
    // Returns the edge id.
    int IsMouseHoveringEdge();
    // Get the nearest tab border in the tab when hovering a
    // container. Returns the tab number or -1 if the tab bar is not
    // hovered or there are no tabs.
    int getNearestTabBorder();

    // Show a drop targets on this window that covers the whole window.
    void showDropTargetFull();
    // Show the drop targets on this window's tab bar
    void showDropTargetOnTabBar();
    // Show the drop targets on the edge of the container. edge is the
    // id for the edge (1:top, 2:right, 3:bottom, 4:left).
    void showDropTargetEdge(int edge);

    // Find the integer index of dock dthis in the stack of this
    // container or h/v-container.
    int OpStack_Find(Dock *dthis);
    // Insert the dock dnew in the stack of this at position ithis (or 
    // at the back if ithis == -1.
    void OpStack_Insert(Dock *dnew, int ithis=-1);
    // Replace the dock replaced with replacement in this dock's
    // stack. If erase, kill the replaced dock.
    void OpStack_Replace(Dock *replaced, Dock *replacement, bool erase);
    // Remove the dock dd from the stack of this container. If erase,
    // kill the dock as well.
    void OpStack_Remove(Dock *dd, bool erase);

    // Replace this dock (part of a root container tree) with a
    // horizontal (type==Type_Horizontal) or vertical
    // (type==Type_Vertical) container. The new container has the
    // current dock plus container dcont (if null, a new dcont is
    // allocated). The new container is placed before (before==true)
    // or after (false) the old one. Returns the new container.
    Dock *OpRoot_ReplaceHV(Type_ type,bool before,Dock *dcont=nullptr);
    // Add a new container (dcont) to the horizontal/vertical parent
    // of this dock.  If !dcont, a new container is allocated. The new
    // container is placed before (before==true) or after (false) the
    // old one, and is returned by this function.
    Dock *OpRoot_AddToHV(bool before,Dock *dcont=nullptr);
    // Fill an empty root container with at least one empty automatic
    // container. 
    void OpRoot_FillEmpty();

    // Raise this dock's window to the top of the window stack.
    void raiseDock();
    // Raise this dock's window to the top of the window
    // stack. Or sink it if it is NoBringToFrontOnFocus.
    void raiseOrSinkDock();
    // Focus a container and its child and parent. Sets the move ID
    // and the active ID.
    void focusContainer();
    // Lift this container from its rootcontainer
    void liftContainer();

    // Add a new dock (dnew) to a container (this) at position ithis
    // in the tab bar, and make it the current tab. If ithis == -1,
    // add it to the end of the tab bar.
    void newDock(Dock *dnew, int ithis = -1);
    // Add a new dock to a root container, or a container docked to a
    // root container. If this is a root or h-v container, the new
    // dock is added to the last container in the tree (bottom
    // right). If the added dock is a container, just add it - if it
    // is a normal dock, create an automatic container for it.  The
    // new dock is added at edge ithis of the container (1:top,
    // 2:right, 3:bottom, 4:left). edge = 5 is used to add to or
    // replace an automatic container in an empty root container.
    // Returns the new container.
    Dock *newDockRoot(Dock *dnew, int iedge);

    // Undock a container, restore its position, size, etc. flags, and
    // place it at the top or the bottom of the window stack. The
    // stack of the parent container is not modified.
    void unDock();
    // Clear all docked windows from a container
    void clearContainer();
    // Clear all docked windows from a root container
    void clearRootContainer();

    // Kill this automatic container if it is empty. If this
    // h/v-container is empty, convert it to a container. If it has
    // one window, kill it and connect its child to its parent.
    void killContainerMaybe();

    // Draw the tab bar of a container. On output, erased is true if a
    // tab in this container was closed.
    void drawTabBar(bool *erased=nullptr);
    // Hide this dock docked to a container on an inactive tab.
    void hideTabWindow();
    // Show this dock docked to a container (dcont) on an active
    // tab. If noresize, do not show the resize grip.
    void showTabWindow(Dock *dcont, bool noresize);
    // Draw the contents of a container. If noresize, do not show the
    // resize grip. On output, erased is true if a tab in this
    // container was closed.
    void drawContainer(bool noresize, bool *erased=nullptr);

    // Traverse the tree of this root container and return its minimum
    // size based on its contents. Recursive.
    void getMinSize(ImVec2 *minsize, ImVec2 *autosize);
    // Center all sliding bars in this root container.
    void resetRootContainerBars();
    // Sets the position of a sliding bar. The sliding bar is on edge
    // iedge (1:top, 2:right, 3:bottom, 4:left) of this container,
    // which is docked to a root container. xpos is the position of
    // the bar given as a fraction of the window size (between 0 and
    // 1).
    void setSlidingBarPosition(int iedge, float xpos);
    // Traverse the tree of this root container and draw all sliding
    // bars in it. Sets the tabsx vector containing the positions of
    // the bars. root is a pointer to the root container. Recursive.
    void drawRootContainerBars(Dock* root);
    // Traverse the tree of this root container and draw all
    // containers in it. Must be called after drawRootContainerBars to
    // have correct sliding bar positions. root is a pointer to the
    // root container. On output, lift contains a pointer to a
    // container to be lifted or null. erased points to the container
    // where one of the tabs was closed. count keeps count of the
    // number of docked windows in the root container. Recursive.
    void drawRootContainer(Dock* root, Dock **lift, Dock **erased, int *count = nullptr);

    // Sets the size of this dock/container in its detached
    // state. Useful when a dock/container is immediately attached in
    // the first pass and does not have the chance to save this
    // variable from the created window.
    void setDetachedDockSize(float x, float y);

    // Close a dock window. This function is used to kill a dock
    // externally. Unlike normal ImGui windows, making p_open false is
    // not enough if the window is docked, because the container has to
    // be modified, too. Therefore, CloseDock() should always be used
    // when setting p_open = false.
    void CloseDock();
  }; // struct Dock

  // Create a root container with the given label. If p_open, with a
  // close button (close status as *p_open). Extra window flags are
  // passed to the container window. Dock flags can also be passed
  // (see above). Returns a pointer to the root container dock
  // object. Root containers can hold containers and docks.
  Dock *RootContainer(const char* label, bool* p_open=nullptr, ImGuiWindowFlags extra_flags=0, DockFlags dock_flags=0);

  // Create a container with the given label. If p_open, with a close
  // button (close status as *p_open). Extra window flags are passed
  // to the container window. Dock flags can also be passed (see
  // above). Returns a pointer to the container dock
  // object. Containers contain docks and can be docked to root
  // containers.
  Dock *Container(const char* label, bool* p_open=nullptr, ImGuiWindowFlags extra_flags=0, DockFlags dock_flags=0);

  // Create/end a dock window. If p_open, with a close button. If
  // p_open, with a close button (close status as *p_open). Extra
  // window flags are passed to the window.  Dock flags can also be
  // passed (see above). If a pointer to a container is passed in
  // oncedock, dock to that container in the first pass. Returns true
  // if the window is open and accepts items (same as ImGui's
  // Begin). The EndDock() function needs to be used after all items
  // have been added. See CloseDock() note above for how to set
  // p_open externally to close a dock.
  bool BeginDock(const char* label, bool* p_open=nullptr, ImGuiWindowFlags extra_flags=0, 
                  DockFlags dock_flags=0, Dock *oncedock=nullptr);
  void EndDock();

  // GetCurrentDock() gives the pointer to the current open dock (same
  // as GetCurrentWindow(), but for docks). Must appear between
  // BeginDock and EndDock. Returns null if no dock is open.
  Dock *GetCurrentDock();

  // Free the memory occupied dock hash tables.
  void ShutdownDock();

  // Print information about the current known docks. For debug purposes.
  void PrintDock__();
  
} // namespace ImGui
#endif
