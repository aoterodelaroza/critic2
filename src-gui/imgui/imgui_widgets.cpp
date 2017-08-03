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
#include <math.h>

using namespace ImGui;

void ImGui::SlidingBar(const char *label, ImGuiWindow* window, ImVec2 *pos, 
		       ImVec2 size, float minx, float maxx, 
		       int direction, float alpha/*=-1*/){
  ImDrawList* dl = window->DrawList;
  ImGuiContext *g = GetCurrentContext();
  bool hovered, held;
  const ImU32 color = GetColorU32(ImGuiCol_FrameBg);
  const ImU32 colorhovered = GetColorU32(ImGuiCol_TextDisabled);
  const ImU32 coloractive = GetColorU32(ImGuiCol_Text);
  
  const ImRect slidingrect(*pos,*pos+size);
  const ImGuiID slidingid = window->GetID(label);
  ButtonBehavior(slidingrect, slidingid, &hovered, &held);

  if (held){
    if (direction == 1)
      pos->x = fmax(fmin(g->IO.MousePos.x,maxx),minx);
    else
      pos->y = fmax(fmin(g->IO.MousePos.y,maxx),minx);
  }

  // draw the rectangle
  if (alpha >= 0.f)
    PushStyleVar(ImGuiStyleVar_Alpha,alpha);
  dl->PushClipRectFullScreen();
  dl->AddRectFilled(slidingrect.GetTL(),slidingrect.GetBR(),
		    held?coloractive:(hovered?colorhovered:color),
		    g->Style.ScrollbarRounding);
  dl->PopClipRect();
  if (alpha >= 0.f)
    PopStyleVar();

}

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
  ImDrawList* drawl = GetWindowDrawList();
  const char* text_end = FindRenderedTextEnd(label);
  ImVec2 text_size = CalcTextSize(label,text_end,true,false);
  ImRect clip_rect = ImRect(pos0,pos1s);
  drawl->AddRectFilled(pos0,pos1,hovered?color_hovered:(activetab?color_active:color));
  drawl->AddRect(pos0,pos1,color_active,0.0f,~0,1.0f);
  RenderTextClipped(pos0,pos1s,label,text_end,&text_size, ImVec2(0.5f,0.5f), &clip_rect);
  
  // draw the "x"
  if (p_open && size.x >= mintabwidth){
    drawl->AddLine(center+ImVec2(-crossz,-crossz), center+ImVec2(crossz,crossz),IsItemHovered()?cross_color_hovered:cross_color);
    drawl->AddLine(center+ImVec2( crossz,-crossz), center+ImVec2(-crossz,crossz),IsItemHovered()?cross_color_hovered:cross_color);
  }

  return clicked;
}

void ImGui::ResizeGripOther(const char *label, ImGuiWindow* window, ImGuiWindow* cwindow){
  ImGuiContext *g = GetCurrentContext();
//  const ImVec2 br = window->Rect().GetBR();
//  ImDrawList* dl = window->DrawList;
  const ImVec2 br = GetCurrentWindow()->Pos + GetCurrentWindow()->Size;
  ImDrawList* dl = GetWindowDrawList();
  const float resize_corner_size = ImMax(g->FontSize * 1.35f, g->Style.WindowRounding + 1.0f + g->FontSize * 0.2f);
  const ImRect resize_rect(br - ImVec2(resize_corner_size * 0.75f, resize_corner_size * 0.75f), br);
  const ImGuiID resize_id = window->GetID("#blehblah");

  // button behavior
  bool hovered, held;
  ButtonBehavior(resize_rect, resize_id, &hovered, &held, true);

  // resize grip (from imgui.cpp)
  ImU32 resize_col = GetColorU32(held ? ImGuiCol_ResizeGripActive : hovered ? ImGuiCol_ResizeGripHovered : ImGuiCol_ResizeGrip);
  dl->PathLineTo(br + ImVec2(-resize_corner_size, -window->BorderSize));
  dl->PathLineTo(br + ImVec2(-window->BorderSize, -resize_corner_size));
  dl->PathArcToFast(ImVec2(br.x - g->Style.WindowRounding - window->BorderSize, br.y - g->Style.WindowRounding - window->BorderSize), g->Style.WindowRounding, 0, 3);
  dl->PathFillConvex(resize_col);
}

