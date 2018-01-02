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

using namespace std;

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
  const float small_alpha = 1e-15;
  ImGuiContext *g = ImGui::GetCurrentContext();
  ImVec4 col = g->Style.Colors[color];
  col.w = small_alpha;
}

// Colors for the widgets
enum ImGuiColWidgets_ {
  ImGuiColWidgets_Slidingbar,
  ImGuiColWidgets_SlidingbarHovered,
  ImGuiColWidgets_SlidingbarActive,
  ImGuiColWidgets_Tab,
  ImGuiColWidgets_TabHovered,
  ImGuiColWidgets_TabPressed,
  ImGuiColWidgets_TabActive,
  ImGuiColWidgets_TabXFg,
  ImGuiColWidgets_TabXFgHovered,
  ImGuiColWidgets_TabXFgActive,
  ImGuiColWidgets_TabXBg,
  ImGuiColWidgets_TabXBgHovered,
  ImGuiColWidgets_TabXBgActive,
  ImGuiColWidgets_TabBorder,
  ImGuiColWidgets_LiftGrip,
  ImGuiColWidgets_LiftGripHovered,
  ImGuiColWidgets_LiftGripActive,
  ImGuiColWidgets_DropTarget,
  ImGuiColWidgets_DropTargetActive,
  ImGuiColWidgets_COUNT,
};

// Style for the widgets
struct ImGuiStyleWidgets_ {
  ImVec4 Colors[ImGuiColWidgets_COUNT];
  float TabRounding;
  float TabBorderSize;
  float DropTargetLooseness;
  float DropTargetMinsizeEdge;
  float DropTargetMaxsizeEdge;
  float DropTargetEdgeFraction;
  float DropTargetFullFraction;
  float TabHeight;
  float TabMaxWidth;
  float CascadeIncrement;
  float SlidingBarWidth;

  void DefaultStyle(){
    TabRounding = 7.0f;
    TabBorderSize = 0.0f;
    DropTargetLooseness = 4.0f;
    DropTargetMinsizeEdge = 20.f;
    DropTargetMaxsizeEdge = 20.f;
    DropTargetEdgeFraction = 0.1f;
    DropTargetFullFraction = 0.4f;
    TabHeight = 19.0f;
    TabMaxWidth = 100.f;
    CascadeIncrement = 25.f;
    SlidingBarWidth = 4.f;
  }

  void DefaultColors(){
    Colors[ImGuiColWidgets_Slidingbar]        = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
    Colors[ImGuiColWidgets_SlidingbarHovered] = ImVec4(0.60f, 0.60f, 0.70f, 1.00f);
    Colors[ImGuiColWidgets_SlidingbarActive]  = ImVec4(0.70f, 0.70f, 0.90f, 1.00f);
    Colors[ImGuiColWidgets_Tab]               = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
    Colors[ImGuiColWidgets_TabHovered]        = ImVec4(0.45f, 0.45f, 0.90f, 1.00f);
    Colors[ImGuiColWidgets_TabPressed]        = ImVec4(0.46f, 0.54f, 0.80f, 1.00f);
    Colors[ImGuiColWidgets_TabActive]         = ImVec4(0.53f, 0.53f, 0.87f, 1.00f);
    Colors[ImGuiColWidgets_TabXFg]            = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
    Colors[ImGuiColWidgets_TabXFgHovered]     = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
    Colors[ImGuiColWidgets_TabXFgActive]      = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
    Colors[ImGuiColWidgets_TabXBg]            = ImVec4(0.80f, 0.20f, 0.00f, 0.00f);
    Colors[ImGuiColWidgets_TabXBgHovered]     = ImVec4(0.80f, 0.20f, 0.00f, 1.00f);
    Colors[ImGuiColWidgets_TabXBgActive]      = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
    Colors[ImGuiColWidgets_TabBorder]         = ImVec4(0.50f, 0.50f, 0.50f, 0.50f);
    Colors[ImGuiColWidgets_LiftGrip]          = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
    Colors[ImGuiColWidgets_LiftGripHovered]   = ImVec4(0.80f, 0.40f, 0.20f, 1.00f);
    Colors[ImGuiColWidgets_LiftGripActive]    = ImVec4(1.00f, 0.40f, 0.20f, 1.00f);
    Colors[ImGuiColWidgets_DropTarget]        = ImVec4(0.43f, 0.43f, 0.43f, 0.43f);
    Colors[ImGuiColWidgets_DropTargetActive]  = ImVec4(0.80f, 0.80f, 0.80f, 0.80f);
  }

  // gettabheight gettabwidth getslidingbarwidth ...
  ImGuiStyleWidgets_(){
    DefaultStyle();
    DefaultColors();
  };
};
extern ImGuiStyleWidgets_ ImGuiStyleWidgets;

// Widgets added to ImGui
namespace ImGui{
  // Returns true if mouse is hovering the inside of a convex
  // polygon.
  bool IsMouseHoveringConvexPoly(const ImVec2* points, const int num_points);

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
  bool ButtonWithX(const char* label, const ImVec2& size, bool activetab,
                   bool *p_open, bool *dragged, bool *dclicked, float alpha = 1.f);

  // A resize grip drawn on window that controls the size of cwindow.
  // On output, dclicked is true if double-click (auto-resize)
  // happened.
  void ResizeGripOther(const char *label, ImGuiWindow* window, ImGuiWindow* cwindow, bool *dclicked=nullptr);

  // Lift grip. A grip the with button colors drawn on the bottom left
  // corner of the window. True if the grip is clicked.
  bool LiftGrip(const char *label, ImGuiWindow* window);

  // xxxx //
  bool ImageInteractive(ImTextureID texture, float a, bool *hover, ImRect *vrect);

  // xxxx //
  bool InvisibleButtonEx(const char* str_id, const ImVec2& size_arg, bool* hovered, bool *held);

  // xxxx //
  void AttachTooltip(const char* desc, float delay, float maxwidth, ImFont* font);

} // namespace ImGui

#endif

