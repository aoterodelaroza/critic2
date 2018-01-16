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

#include "settings.h"
#include "message.h"
#include "imgui/imgui.h"
#include "imgui/imgui_internal.h"
#include "imgui/icon_glyphs.h"

#include <list>

using namespace ImGui;

static std::list<Message*> mlist = {};

void NewMessage(MessageType_ type,const char *text){
  Message *m = new Message;

  m->time = GetTime();
  m->message = ImStrdup(text);
  m->type = type;
  mlist.push_back(m);
}

void MessageDispatch(){
  ImGuiContext *g = GetCurrentContext();
  float time = GetTime();

  float lasty = g->IO.DisplaySize.y;
  int n = 0;
  auto it = mlist.begin();
  while (it != mlist.end()){
    float dt = time - (*it)->time;
    if ((*it) && (dt > ImGuiStyleUI.MessageExpire)){
      delete (*it);
      mlist.erase(it++);
      continue;
    }

    ImVec4 color;
    if ((*it)->type == Message_Info)
      color = ImGuiStyleUI.Colors[ImGuiColUI_MessageInfo];
    else if ((*it)->type == Message_Warning)
      color = ImGuiStyleUI.Colors[ImGuiColUI_MessageWarning];
    else if ((*it)->type == Message_Error)
      color = ImGuiStyleUI.Colors[ImGuiColUI_MessageError];
    float alphaf = 1.0f;
    if (dt > 0.75f * ImGuiStyleUI.MessageExpire)
      alphaf = 1.f - (dt - 0.75f * ImGuiStyleUI.MessageExpire) / (0.25f * ImGuiStyleUI.MessageExpire);
    color.w *= alphaf;
    PushStyleColor(ImGuiCol_WindowBg,GetColorU32(color));
    
    ++n;
    PushFont(fonticon);
    float messagewidth = ImGuiStyleUI.MessageWidth;
    ImVec2 size1 = CalcTextSize(ICON_SM_INFO,NULL,false,messagewidth);
    PopFont();
    ImVec2 size2 = CalcTextSize((*it)->message,NULL,false,messagewidth - size1.x - 2.f * g->Style.ItemSpacing.x);
    ImVec2 size = {messagewidth,std::max(size2.y + 0.5f * (ImGuiStyleUI.FontSizeIcon - ImGuiStyleUI.FontSize),size1.y) + 2.f * g->Style.WindowPadding.y};
    ImVec2 pos = {g->IO.DisplaySize.x - messagewidth - 2.f*g->Style.ItemSpacing.x,
		  lasty - 2.f*g->Style.ItemSpacing.y - size.y};
    SetNextWindowPos(pos);
    SetNextWindowSize(size);
    PushStyleVar(ImGuiStyleVar_WindowBorderSize,0.0f);

    char tmp[20];
    ImFormatString(tmp,IM_ARRAYSIZE(tmp),"##message__%d__",n);
    Begin(tmp,nullptr,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar|
	  ImGuiWindowFlags_NoCollapse|ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoSavedSettings|
	  ImGuiWindowFlags_AlwaysAutoResize);

    PushFont(fonticon);
    if ((*it)->type == Message_Info)
      Text(ICON_SM_INFO);
    else if ((*it)->type == Message_Warning)
      Text(ICON_SM_WARN);
    else if ((*it)->type == Message_Error)
      Text(ICON_SM_ERROR);
    SameLine();
    PopFont();

    GetCurrentWindow()->DC.CurrentLineTextBaseOffset = 0.5f * (ImGuiStyleUI.FontSizeIcon - ImGuiStyleUI.FontSize);
    TextWrapped((*it)->message);
    bool toerase = IsWindowHovered() && (IsMouseClicked(0) || IsMouseClicked(1));
    End();

    PopStyleColor(1);
    PopStyleVar(1);
    lasty = pos.y;
    if (toerase){
      delete (*it);
      mlist.erase(it++);
    } else {
      ++it;
    }
  }
}

