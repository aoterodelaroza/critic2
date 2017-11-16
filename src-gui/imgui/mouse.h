// -*-c++-*-

// Mouse state struct

#ifndef MOUSE_H
#define MOUSE_H

#include "imgui.h"
#include "imgui_internal.h"
#include <glm/glm.hpp>

// A mouse class to encapsulate all information about the state of the
// mouse when it is hovering a certain window region.
struct MouseState
{
  bool hover   = false;
  bool lclick  = false;
  bool ldown   = false;
  bool mclick  = false;
  bool mdown   = false;
  bool rclick  = false;
  bool rdown   = false;
  bool ldrag   = false;
  bool rdrag   = false;
  bool ldclick = false;
  float scroll = 0.f;
  glm::vec2 pos = {0.f,0.f};
  glm::vec2 ndpos = {0.f,0.f};

  void Fill(){
    ImGuiContext *g = ImGui::GetCurrentContext();

    hover = false;
    lclick = ImGui::IsMouseClicked(0);
    mclick = ImGui::IsMouseClicked(2);
    rclick = ImGui::IsMouseClicked(1);
    ldclick = ImGui::IsMouseDoubleClicked(0);

    ldrag = ImGui::IsMouseDragging(0);
    rdrag = ImGui::IsMouseDragging(1);

    ldown = g->IO.MouseDown[0];
    rdown = g->IO.MouseDown[1];
    mdown = g->IO.MouseDown[2];

    pos = {g->IO.MousePos.x,g->IO.MousePos.y};
    ndpos = {0.f,0.f};
  
    scroll = g->IO.MouseWheel;
  }
};
#endif
