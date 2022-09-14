// Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
// Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
// <victor@fluor.quimica.uniovi.es>.
//
// critic2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at
// your option) any later version.
//
// critic2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <http://www.gnu.org/licenses/>.

#include "my_imgui.h"

using namespace ImGui;

bool my_CloseButton(const char* str_id, ImVec4 buttonbg) {
  float sz = GetFrameHeight();
  ImVec2 size = {sz,sz};
  ImGuiButtonFlags flags = ImGuiButtonFlags_None;

  ImGuiContext& g = *GImGui;
  ImGuiWindow* window = GetCurrentWindow();
  if (window->SkipItems)
    return false;

  const ImGuiID id = window->GetID(str_id);
  const ImRect bb(window->DC.CursorPos, window->DC.CursorPos + size);
  const float default_size = GetFrameHeight();
  ItemSize(size, (size.y >= default_size) ? g.Style.FramePadding.y : -1.0f);
  if (!ItemAdd(bb, id))
    return false;

  if (g.LastItemData.InFlags & ImGuiItemFlags_ButtonRepeat)
    flags |= ImGuiButtonFlags_Repeat;

  bool hovered, held;
  bool pressed = ButtonBehavior(bb, id, &hovered, &held, flags);

  // Render
  const ImU32 col = GetColorU32((held && hovered) ? ImGuiCol_ButtonActive : hovered ? ImGuiCol_ButtonHovered : GetColorU32(buttonbg));
  ImVec2 center = bb.GetCenter();
  window->DrawList->AddCircleFilled(center, ImMax(2.0f, g.FontSize * 0.5f + 1.0f), col, 12);

  float cross_extent = g.FontSize * 0.5f * 0.7071f - 1.0f;
  ImU32 cross_col = GetColorU32(ImGuiCol_Text);
  center -= ImVec2(0.5f, 0.5f);
  window->DrawList->AddLine(center + ImVec2(+cross_extent, +cross_extent), center + ImVec2(-cross_extent, -cross_extent), cross_col, 1.0f);
  window->DrawList->AddLine(center + ImVec2(+cross_extent, -cross_extent), center + ImVec2(-cross_extent, +cross_extent), cross_col, 1.0f);

  return pressed;
}
