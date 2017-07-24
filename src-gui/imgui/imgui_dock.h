// -*-c++-*-
// Rewritten from: git@github.com:vassvik/imgui_docking_minimal.git

#ifndef IMGUI_DOCK_H
#define IMGUI_DOCK_H

#include "imgui.h"
#include "imgui_internal.h"
#include <unordered_map>
#include <list>

using namespace std;

namespace ImGui{

  struct Dock{
    enum Type_{Type_None,Type_Container,Type_Dock};
    enum Status_{Status_None,Status_Open,Status_Collapsed,Status_Closed,
		 Status_Dragged,Status_Docked};
    enum Stack_{Stack_None,Stack_Vertical,Stack_Horizontal,Stack_Leaf};

    ImU32 id = 0; // id of the window
    char* label = nullptr; // window label
    ImVec2 pos; // position of the window
    ImVec2 size; // size of the window
    ImGuiWindowFlags flags = 0; // flags for the window
    bool hidden = false; // whether a docked window is hidden
    bool collapsed = false; // whether a docked window is collapsed
    Type_ type = Type_None; // type of docking window
    Status_ status = Status_None; // status of the docking window
    Stack_ sttype = Stack_None; // stacking type when docked
    list<Dock *> stack = {}; // stack of docks at this level
    Dock *currenttab = nullptr; // currently selected tab
    Dock *parent = nullptr; // node immediately above current node
    Dock *root = nullptr; // root node of this tree
    ImGuiWindow* window = nullptr; // associated window

    Dock(): pos(0,0), size(-1,-1) {};
    ~Dock(){ MemFree(label);}

    // Show the drop targets for this window
    void showDropTarget();
    // Add a new dock to a container
    void newDock(Dock *dnew);
    // Draw container
    void drawContainer();
    // Draw the tab bar of a tabbed container
    void drawTabBar();

  }; // struct Dock

  struct DockContext{
    struct LabelHash{
      size_t operator()(const char *key) const{
	return (size_t) ImHash((const void *)key,0);
      };
    };
    unordered_map<const char *,Dock*,LabelHash> dockht = {};
    Dock* m_current = nullptr;
    DockContext(){};
    ~DockContext(){};
    Dock *getContainerAt(const ImVec2& pos);
  }; // struct DockContext

  // Create a container with the given label. If p_open, with a close
  // button.  Extra window flags are passed to the container window.
  void Container(const char* label, bool* p_open=nullptr, ImGuiWindowFlags extra_flags=0);

  // Create/end a window that can be docked to a container. If p_open,
  // with a close button.  Extra window flags are passed to the
  // container window.
  bool BeginDock(const char* label, bool* p_open=nullptr, ImGuiWindowFlags extra_flags=0);
  void EndDock();

  // Print information about the current known docks. Debug only.
  void Print();
  
} // namespace ImGui
#endif
