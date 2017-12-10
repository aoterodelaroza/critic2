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

#include <algorithm>
#include <unordered_map>

using namespace ImGui;

// Tooltip variables and behavior //
// Variables:
// - tooltip_shownid: ID for the tooltip element shown.
// - tooltip_hoveredid: ID for the tooltip element hovered last frame.
// - tooltip_thisframe: saves the current frame to detect when
//   AttachTooltip is called for the first time.
// - tooltip_washovered: set to true if any AttachTooltip element
//   has been hovered in the last frame.
// - tooltip_time: the last time when no item was hovered.
// - tooltip_lastactive - the last time when a tooltip was shown.
// How it works:
// The tooltip is shown if:
// - The mouse hovers over a single tooltip element for longer than
//   delay, which is passed as a function argument.
// - The mouse hovers over a tooltip element and it was already
//   showing a tooltip from a different element. In this case,
//   there is no delay.
// - The mouse hovers over a tooltip element, then hovers none for a
//   period of time less than delay, then overs over a new tooltip
//   element. In this case, the delay is calculated as the difference
//   between the time the mouse hovered over the new tooltip element
//   (tooltip_time) and the time the mouse left the previous tooltip
//   element (tooltip_lastactive), with a maximum equal to
//   delay. Thus, the delay is equal to the time the mouse spent
//   outside any tooltip elment.
static ImGuiID tooltip_shownid = 0;
static ImGuiID tooltip_hoveredid = 0;
static int tooltip_thisframe = -1;
static bool tooltip_washovered = false;
static float tooltip_time = 0.f;
static float tooltip_lastactive = 0.f;

// Variable for the widget styles
ImGuiStyleWidgets_ ImGuiStyleWidgets;

// Constructor for the style struct //
ImGuiStyleWidgets_::ImGuiStyleWidgets_(){
  // colors
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

  // style vars
  TabRounding = 7.0f;
  TabBorderSize = 0.0f;
  DropTargetLooseness = 4.0f;
  DropTargetMinsizeEdge = 30.f;
  DropTargetMaxsizeEdge = 90.f;
  DropTargetEdgeFraction = 0.1f;
  DropTargetFullFraction = 0.4f;
  TabHeight = 19.0f;
  TabMaxWidth = 100.f;
  CascadeIncrement = 25.f;
  SlidingBarWidth = 4.f;
};

// Function definitions //

bool ImGui::IsMouseHoveringConvexPoly(const ImVec2* points, const int num_points){
    ImGuiContext& g = *GImGui;

    float p0[num_points][2];
    
    // differences
    for (int i = 0; i < num_points; i++){
      p0[i][0] = points[i].x - g.IO.MousePos.x;
      p0[i][1] = points[i].y - g.IO.MousePos.y;
    }

    // first sign
    float a0 = p0[0][0] * p0[num_points-1][1] - p0[num_points-1][0] * p0[0][1];

    // compare to all other signs
    for (int i = 0; i < num_points-1; i++){
      float a = p0[i+1][0] * p0[i][1] - p0[i][0] * p0[i+1][1];
      if (a * a0 < 0.f) 
	return false;
    }

    return true;
}

void ImGui::SlidingBar(const char *label, ImGuiWindow* window, ImVec2 *pos, 
                       ImVec2 size, float minx, float maxx, int direction){
  ImDrawList* dl = window->DrawList;
  ImGuiContext *g = GetCurrentContext();
  bool hovered, held;
  const ImU32 color = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_Slidingbar]);
  const ImU32 coloractive = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarHovered]);
  const ImU32 colorhovered = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarActive]);
  
  const ImRect slidingrect(*pos,*pos+size);
  const ImGuiID slidingid = window->GetID(label);
  ButtonBehavior(slidingrect, slidingid, &hovered, &held);

  if (held){
    if (direction == 1)
      pos->x = max(min(g->IO.MousePos.x - 0.5f * size.x,maxx),minx);
    else
      pos->y = max(min(g->IO.MousePos.y - 0.5f * size.y,maxx),minx);
  }

  // draw the rectangle
  dl->PushClipRectFullScreen();
  dl->AddRectFilled(slidingrect.GetTL(),slidingrect.GetBR(),
                    held?coloractive:(hovered?colorhovered:color),
                    g->Style.ScrollbarRounding);
  dl->PopClipRect();
}

bool ImGui::ButtonWithX(const char* label, const ImVec2& size, bool activetab,
                        bool *p_open, bool *dragged, bool *dclicked, float alphamul /*=1.f*/){
  ImGuiContext *g = GetCurrentContext();
  const float crossz = round(0.3 * g->FontSize);
  const float crosswidth = 3.5f * crossz + 6;
  const float mintabwidth = 2 * crosswidth + 1;

  const ImU32 colorxfg = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFg]);
  const ImU32 colorxfg_hovered = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgHovered]);
  const ImU32 colorxfg_pressed = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgActive]);
  const ImU32 colorxbg = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBg]);
  const ImU32 colorxbg_hovered = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgHovered]);
  const ImU32 colorxbg_pressed = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgActive]);
  ImU32 color = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_Tab]);
  ImU32 color_active  = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabActive]);
  ImU32 color_pressed = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabPressed]);
  ImU32 color_hovered = GetColorU32(ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabHovered]);

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
  bool pressed = IsItemActive() && IsMouseDown(0);
  *dragged = IsItemActive() && IsMouseDragging();
  *dclicked = IsItemActive() && IsMouseDoubleClicked(0);
  
  // draw the close button, if this window can be closed
  ImVec2 center;
  bool xhovered = false, xpressed = false;
  if (p_open && size.x >= mintabwidth){
    // draw the close button itself
    SameLine();
    char tmp2[strlen(label)+6];
    ImFormatString(tmp2,IM_ARRAYSIZE(tmp2),"%s__x__",label);
    ImVec2 smallsize(crosswidth, size.y);

    // cross button
    *p_open = !(InvisibleButton(tmp2, smallsize));

    // update output flags and variables for drawing
    *dragged = *dragged | (p_open && IsItemActive() && IsMouseDragging());
    xhovered = IsItemHovered();
    hovered |= xhovered;
    xpressed = IsItemActive() && IsMouseDown(0);
    center = ((GetItemRectMin() + GetItemRectMax()) * 0.5f);
  }
  ImVec2 pos1 = GetItemRectMax();

  // rectangle and text
  ImDrawList* drawl = GetWindowDrawList();
  const char* text_end = FindRenderedTextEnd(label);
  ImVec2 text_size = CalcTextSize(label,text_end,true,false);
  ImRect clip_rect = ImRect(pos0,pos1s);
  drawl->AddRectFilled(pos0,pos1,activetab? color_active:
                       pressed? color_pressed:
                       hovered? color_hovered:
                       color,ImGuiStyleWidgets.TabRounding,ImDrawCornerFlags_TopLeft|ImDrawCornerFlags_TopRight);

  if (ImGuiStyleWidgets.TabBorderSize > 0.0f)
    drawl->AddRect(pos0,pos1,GetColorU32(g->Style.Colors[ImGuiColWidgets_TabBorder]),
		   ImGuiStyleWidgets.TabRounding,ImDrawCornerFlags_TopLeft|ImDrawCornerFlags_TopRight,1.0f);
  RenderTextClipped(pos0,pos1s,label,text_end,&text_size, ImVec2(0.5f,0.5f), &clip_rect);
  
  // draw the "x"
  if (p_open && size.x >= mintabwidth){
    drawl->AddCircleFilled(center,crossz * sqrt(2.f) * 1.25f,
                           xpressed?colorxbg_pressed:
                           xhovered?colorxbg_hovered:
                           colorxbg,36);
    drawl->AddLine(center+ImVec2(-crossz,-crossz), center+ImVec2(crossz,crossz),
                   xpressed?colorxfg_pressed:
                   xhovered?colorxfg_hovered:
                   colorxfg);
    drawl->AddLine(center+ImVec2( crossz,-crossz), center+ImVec2(-crossz,crossz),
                   xpressed?colorxfg_pressed:
                   xhovered?colorxfg_hovered:
                   colorxfg);
  }

  return clicked;
}

void ImGui::ResizeGripOther(const char *label, ImGuiWindow* window, ImGuiWindow* cwindow, bool *dclicked/*=nullptr*/){
  static bool first = true;
  if (dclicked) *dclicked = false;

  static ImVec2 pos_orig = {};
  static ImVec2 size_orig = {};
  static ImVec2 csize_orig = {};
  ImGuiContext *g = GetCurrentContext();
  const ImVec2 br = window->Rect().GetBR();
  ImDrawList* dl = window->DrawList;
  const float resize_corner_size = ImMax(g->FontSize * 1.35f, g->Style.WindowRounding + 1.0f + g->FontSize * 0.2f);
  const ImRect resize_rect(br - ImVec2(resize_corner_size * 0.75f, resize_corner_size * 0.75f), br);
  char tmp[strlen(label)+15];
  ImFormatString(tmp,IM_ARRAYSIZE(tmp),"%s__resize__",label);
  const ImGuiID resize_id = window->GetID(tmp);

  // no clipping; save previous clipping
  ImRect saverect = window->ClipRect;
  window->ClipRect = ImVec4(-FLT_MAX,-FLT_MAX,+FLT_MAX,+FLT_MAX);
  dl->PushClipRectFullScreen();

  // button behavior
  bool hovered, held;
  ButtonBehavior(resize_rect, resize_id, &hovered, &held, ImGuiButtonFlags_FlattenChilds);

  // update the static flags
  if (held){
    if (first) {
      pos_orig = window->Pos;
      size_orig = window->SizeFull;
      csize_orig = cwindow->SizeFull;
    }
    first = false;
  } else {
    first = true;
  }

  // mouse cursor
  if (hovered || held)
    g->MouseCursor = ImGuiMouseCursor_ResizeNWSE;

  // apply the size change
  if (g->HoveredWindow == window && held && g->IO.MouseDoubleClicked[0]){
    if (dclicked) *dclicked = true;
    ImVec2 size_auto_fit = ImClamp(cwindow->SizeContents + cwindow->WindowPadding, g->Style.WindowMinSize, 
                            ImMax(g->Style.WindowMinSize, g->IO.DisplaySize - g->Style.DisplaySafeAreaPadding));

    cwindow->SizeFull = size_auto_fit;
    ClearActiveID();
  } else if (held)
    cwindow->SizeFull = csize_orig + (g->IO.MousePos - g->ActiveIdClickOffset + resize_rect.GetSize() - pos_orig) - size_orig;
  cwindow->Size = cwindow->SizeFull;

  // resize grip (from imgui.cpp)
  ImU32 resize_col = GetColorU32(held ? ImGuiCol_ResizeGripActive : hovered ? ImGuiCol_ResizeGripHovered : ImGuiCol_ResizeGrip);
  dl->PathLineTo(br + ImVec2(-resize_corner_size, -window->WindowBorderSize));
  dl->PathLineTo(br + ImVec2(-window->WindowBorderSize, -resize_corner_size));
  dl->PathArcToFast(ImVec2(br.x - g->Style.WindowRounding - window->WindowBorderSize, br.y - g->Style.WindowRounding - window->WindowBorderSize), g->Style.WindowRounding, 0, 3);
  dl->PathFillConvex(resize_col);
  dl->PopClipRect();
  window->ClipRect = saverect;
}

bool ImGui::LiftGrip(const char *label, ImGuiWindow* window){
  ImGuiContext *g = GetCurrentContext();
  const ImVec2 bl = window->Rect().GetBL();
  ImDrawList* dl = window->DrawList;
  const float lift_corner_size = ImMax(g->FontSize * 1.35f, g->Style.WindowRounding + 1.0f + g->FontSize * 0.2f);
  const ImRect lift_rect(bl - ImVec2(0.f, lift_corner_size * 0.75f), bl + ImVec2(lift_corner_size * 0.75f, 0.f));
  char tmp[strlen(label)+15];
  ImFormatString(tmp,IM_ARRAYSIZE(tmp),"%s__lift__",label);
  const ImGuiID lift_id = window->GetID(tmp);

  // no clipping; save previous clipping
  ImRect saverect = window->ClipRect;
  window->ClipRect = ImVec4(-FLT_MAX,-FLT_MAX,+FLT_MAX,+FLT_MAX);
  dl->PushClipRectFullScreen();

  // button behavior
  bool hovered, held;
  ButtonBehavior(lift_rect, lift_id, &hovered, &held, ImGuiButtonFlags_FlattenChilds);

  // lift grip (from imgui.cpp's resize grip)
  ImU32 lift_col = GetColorU32(held?ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripActive]: 
                               hovered?ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripHovered]: 
                               ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGrip]);
  dl->PathLineTo(bl + ImVec2(window->WindowBorderSize, -lift_corner_size));
  dl->PathLineTo(bl + ImVec2(lift_corner_size, -window->WindowBorderSize));
  dl->PathArcToFast(ImVec2(bl.x + g->Style.WindowRounding + window->WindowBorderSize, bl.y - g->Style.WindowRounding - window->WindowBorderSize), g->Style.WindowRounding, 3, 6);
  dl->PathFillConvex(lift_col);
  dl->PopClipRect();
  window->ClipRect = saverect;
  
  return held && IsMouseDragging();
}

bool ImGui::ImageInteractive(ImTextureID texture, float a, bool *hover, ImRect *vrect){
  ImGuiWindow *win = GetCurrentWindow(); 
  if (win->SkipItems)
    return false;
  ImGuiContext *g = GetCurrentContext();

  PushID((void *)texture);
  const ImGuiID id = win->GetID("#imageinteractive");
  PopID();

  vrect->Min = win->DC.CursorPos;
  vrect->Max = win->DC.CursorPos + win->ContentsRegionRect.Max - (win->DC.CursorPos - win->Pos + ImVec2(1.f,1.f));

  if (!ItemAdd(*vrect, id))
      return false;

  float x = vrect->Max.x - vrect->Min.x;
  float y = vrect->Max.y - vrect->Min.y;
  float xratio = x/fmax(x,y);
  float yratio = y/fmax(x,y);
  float rx = 0.5f * (1.f - xratio) * a;
  float ry = 0.5f * (1.f - yratio) * a;

  bool held;
  bool pressed = ButtonBehavior(*vrect, id, hover, &held);
  win->DrawList->AddImage(texture,vrect->Min,vrect->Max,ImVec2(rx, a - ry),ImVec2(a - rx, ry));
}

bool ImGui::InvisibleButtonEx(const char* str_id, const ImVec2& size_arg, bool* hovered, bool *held){
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems)
        return false;

    const ImGuiID id = window->GetID(str_id);
    ImVec2 size = CalcItemSize(size_arg, 0.0f, 0.0f);
    const ImRect bb(window->DC.CursorPos, window->DC.CursorPos + size);
    ItemSize(bb);
    if (!ItemAdd(bb, id))
        return false;

    bool pressed = ButtonBehavior(bb, id, hovered, held);

    return pressed;
}

void ImGui::AttachTooltip(const char* desc, float delay, float maxwidth, ImFont* font){
  ImGuiContext *g = GetCurrentContext();
  ImGuiID id = g->CurrentWindow->DC.LastItemId;
  float time = GetTime();
  int thisframe = GetFrameCount();

  if (thisframe != tooltip_thisframe){
    // run once every frame, in the first call
    if (!tooltip_washovered){
      tooltip_time = time;
      tooltip_shownid = 0;
      tooltip_hoveredid = 0;
    }
    tooltip_thisframe = thisframe;
    tooltip_washovered = false;
  }

  if (g->HoveredId == id){
    // If no tooltip is being shown and the mouse moves from one tooltip element
    // to another, this is the same as if it moved from a zone without any
    // tooltip elements.
    if (id != tooltip_hoveredid && tooltip_shownid == 0)
      tooltip_time = time;

    if (tooltip_lastactive != 0.f)
      delay = fmin(delay,fmax(tooltip_time - tooltip_lastactive,0.f));

    tooltip_washovered = true;
    tooltip_hoveredid = id;

    if (time - tooltip_time > delay){
      tooltip_shownid = id;
      tooltip_lastactive = time;
      BeginTooltip();
      PushTextWrapPos(maxwidth);
      PushFont(font);
      TextUnformatted(desc);
      PopFont();
      PopTextWrapPos();
      EndTooltip();
    }
  }
}
