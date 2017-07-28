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

#include "imgui_widgets.h"
#include "imgui.h"
#include "imgui_internal.h"

using namespace ImGui;

bool ImGui::ButtonWithX(const char* label, const ImVec2& size, bool activetab, bool buttoncol,
			bool *p_open, bool *dragged, bool *dclicked, bool *closeclicked, 
			float alphamul /*=1.f*/){
  // lengths and colors
  ImGuiContext *g = GetCurrentContext();
  const float crossz = round(0.3 * g->FontSize);
  const float crosswidth = 2 * crossz + 6;
  const float mintabwidth = 2 * crosswidth + 1;
  const ImU32 cross_color = GetColorU32(ImGuiCol_Text,alphamul);
  const ImU32 cross_color_hovered = GetColorU32(ImGuiCol_TextSelectedBg,alphamul);
  ImU32 color = GetColorU32(ImGuiCol_FrameBg,alphamul);
  ImU32 color_active = GetColorU32(ImGuiCol_FrameBgActive,alphamul);
  ImU32 color_hovered = GetColorU32(ImGuiCol_FrameBgHovered,alphamul);
  if (buttoncol){
    color = GetColorU32(ImGuiCol_Button,alphamul);
    color_active = GetColorU32(ImGuiCol_ButtonActive,alphamul);
    color_hovered = GetColorU32(ImGuiCol_ButtonHovered,alphamul);
  }

  // size of the main button
  ImVec2 mainsize = size;
  if (p_open && size.x >= mintabwidth)
    mainsize.x -= crosswidth;

  // main button
  bool clicked = InvisibleButton(label, mainsize);

  // some positions and other variables
  bool hovered = IsItemHovered();
  ImVec2 pos0 = GetItemRectMin();
  ImVec2 pos1s = GetItemRectMax();

  // set the output flags for the main button
  *dragged = IsItemActive() && IsMouseDragging();
  *dclicked = IsItemActive() && IsMouseDoubleClicked(0);
  *closeclicked = false;
  
  // draw the close button, if this window can be closed
  ImVec2 center;
  if (p_open && size.x >= mintabwidth){
    // draw the close button itself
    SameLine();
    char tmp2[strlen(label)+6];
    ImFormatString(tmp2,IM_ARRAYSIZE(tmp2),"%s__x__",label);
    ImVec2 smallsize(crosswidth, size.y);

    // cross button
    *closeclicked = (InvisibleButton(tmp2, smallsize));

    // update output flags and variables for drawing
    *dragged = *dragged | (!*closeclicked && IsItemActive() && IsMouseDragging());
    hovered |= IsItemHovered();
    center = ((GetItemRectMin() + GetItemRectMax()) * 0.5f);
  }
  ImVec2 pos1 = GetItemRectMax();

  // rectangle and text
  ImDrawList* draw_list = GetWindowDrawList();
  const char* text_end = FindRenderedTextEnd(label);
  ImVec2 text_size = CalcTextSize(label,text_end,true,false);
  ImRect clip_rect = ImRect(pos0,pos1s);
  draw_list->AddRectFilled(pos0,pos1,hovered?color_hovered:(activetab?color_active:color));
  draw_list->AddRect(pos0,pos1,color_active,0.0f,~0,1.0f);
  RenderTextClipped(pos0,pos1s,label,text_end,&text_size, ImVec2(0.5f,0.5f), &clip_rect);
  
  // draw the "x"
  if (p_open && size.x >= mintabwidth){
    draw_list->AddLine(center+ImVec2(-crossz,-crossz), center+ImVec2(crossz,crossz),IsItemHovered()?cross_color_hovered:cross_color);
    draw_list->AddLine(center+ImVec2( crossz,-crossz), center+ImVec2(-crossz,crossz),IsItemHovered()?cross_color_hovered:cross_color);
  }

  return clicked;
}
