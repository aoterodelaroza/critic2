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

#ifndef IMGUI_WIDGETS_H
#define IMGUI_WIDGETS_H

#include "imgui_widgets.h"
#include "imgui.h"
#include "imgui_internal.h"

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

namespace ImGui{

  // Sliding bar for root container splits
  void SlidingBar(ImGuiWindow* window, ImVec2 *pos, ImVec2 size, ImVec2 limits, int direction,
		  float alpha=-1.f);

  // Button with a clickable "X" at the end
  bool ButtonWithX(const char* label, const ImVec2& size, bool activetab, bool buttoncol,
		   bool *p_open, bool *dragged, bool *dclicked, bool *closeclicked, 
		   float alpha = 1.f);

} // namespace ImGui
#endif
