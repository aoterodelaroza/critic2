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

// Several widgets independent from each other are provided by this
// file. See below for a description.

#ifndef IMGUI_WIDGETS_H
#define IMGUI_WIDGETS_H

#include "imgui_widgets.h"
#include "imgui.h"
#include "imgui_internal.h"
#include "mouse.h"

using namespace std;

// Named constants
const float small_alpha = 1e-15;

// Helper functions
static inline ImVec2 operator+(ImVec2 lhs, ImVec2 rhs) {
    return ImVec2(lhs.x+rhs.x, lhs.y+rhs.y);
}
static inline ImVec2 operator-(ImVec2 lhs, ImVec2 rhs) {
    return ImVec2(lhs.x-rhs.x, lhs.y-rhs.y);
}
static inline ImVec2 operator*(ImVec2 lhs, float rhs) {
    return ImVec2(lhs.x*rhs, lhs.y*rhs);
}
static inline ImVec4 OpaqueColor(ImGuiCol_ color, float newalpha){
  ImGuiContext *g = ImGui::GetCurrentContext();
  ImVec4 col = g->Style.Colors[color];
  col.w = newalpha;
}
static inline ImVec4 TransparentColor(ImGuiCol_ color){
  ImGuiContext *g = ImGui::GetCurrentContext();
  ImVec4 col = g->Style.Colors[color];
  col.w = small_alpha;
}

namespace ImGui{

  // Sliding bar for splits. label: used to calculate the ID. window:
  // window containing the bar. pos: position of the top left of the bar on
  // input and output. size: size of the bar. minx and maxx: minimum and maximum
  // positions in direction direction (1=x, 2=y).
  void SlidingBar(const char* label, ImGuiWindow* window, ImVec2 *pos, ImVec2 size, 
                  float minx, float maxx, int direction);

  // Button with a clickable "X" at the end. label: label for the
  // button (and generates the ID of the main button). size: size of
  // the button. activetab: whether the button uses the "active" or
  // "inactive" color. scrollbarcol: if true, use scroll bar grab
  // colors (otherwise, use framebg colors). p_open whether the "x" is
  // shown and close status. dragged: on output, true if the button is
  // being dragged. dclicked: on output, true if the button was double 
  // clicked. closeclicked: on output, true if the X has been clicked.
  // alphamul: alpha multiplier for all colors. Returns true if the
  // main part of the button (not the x) has been clicked.
  bool ButtonWithX(const char* label, const ImVec2& size, bool activetab, bool scrollbarcol,
                   bool *p_open, bool *dragged, bool *dclicked, float alpha = 1.f);

  // A resize grip drawn on window that controls the size of cwindow.
  // On output, dclicked is true if double-click (auto-resize)
  // happened.
  void ResizeGripOther(const char *label, ImGuiWindow* window, ImGuiWindow* cwindow, bool *dclicked=nullptr);

  // Lift grip. A grip the with button colors drawn on the bottom left
  // corner of the window. True if the grip is clicked.
  bool LiftGrip(const char *label, ImGuiWindow* window);

  // xxxx //
  bool ImageInteractive(ImTextureID texture, bool *hover, ImVec2 *ndpos);

  // xxxx //
  bool InvisibleButtonEx(const char* str_id, const ImVec2& size_arg, bool* hovered, bool *held);
} // namespace ImGui
#endif
