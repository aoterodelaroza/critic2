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
#include <GLFW/glfw3.h>

#include "imgui/additional_fonts.h"
#include "imgui/imgui_widgets.h"
#include "json/json.hpp"

#include <unistd.h>
#include <sys/stat.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace ImGui;
using json = nlohmann::json;

// Global variables: font pointers and bake size (see settings.h)
ImFont* fontdefault = nullptr;
ImFont* fonticon = nullptr;
const float fontsizebake = 36.0f;

// UI style
ImGuiStyleUI_ ImGuiStyleUI;

// Shader
Shader *shader;

// View defaults and settings
bool view_wireframe; // render objects in wireframe
bool view_orthogonal; // orthogonal/perspective projection
float view_fov; // field of view (degrees)
float view_resetdistance; // distance to scene on reset (scenerad)
float view_bgrgb[4]; // background color
glm::vec3 view_lightpos; // light position relative to the camera
glm::vec3 view_lightcolor; // light color
float view_ambient; // ambient light intensity
float view_diffuse; // diffuse light intensity
float view_specular; // specular light intensity
int view_shininess; // light shininess
int view_isphres; // resolution of the atom spheres
int view_icylres; // resolution of the bond cylinders
float view_mousesens_rot; // mouse rotation sensitivity (scale factor)
float view_mousesens_zoom; // mouse zoom sensitivity (scale factor)

// configuration file name
string conffile = "";

void DefaultSettings(){
  static bool firstpass = true;

  // default colors
  UIStyleColorsClassic();

  // default ImGui settings (imgui.h)
  ImGuiStyle& style = GetStyle();
  style.Alpha                   = 1.0f;
  style.WindowPadding           = ImVec2(8,8);
  style.WindowRounding          = 7.0f;
  style.WindowBorderSize        = 0.0f;
  style.WindowMinSize           = ImVec2(32,32);
  style.WindowTitleAlign        = ImVec2(0.0f,0.5f);
  style.ChildRounding           = 0.0f;
  style.ChildBorderSize         = 1.0f;
  style.PopupRounding           = 0.0f;
  style.PopupBorderSize         = 1.0f;
  style.FramePadding            = ImVec2(4,3);
  style.FrameRounding           = 0.0f;
  style.FrameBorderSize         = 0.0f;
  style.ItemSpacing             = ImVec2(8,4);
  style.ItemInnerSpacing        = ImVec2(4,4);
  style.TouchExtraPadding       = ImVec2(0,0);
  style.IndentSpacing           = 21.0f;
  style.ColumnsMinSpacing       = 6.0f;
  style.ScrollbarSize           = 16.0f;
  style.ScrollbarRounding       = 9.0f;
  style.GrabMinSize             = 10.0f;
  style.GrabRounding            = 0.0f;
  style.ButtonTextAlign         = ImVec2(0.5f,0.5f);
  style.DisplayWindowPadding    = ImVec2(22,22);
  style.DisplaySafeAreaPadding  = ImVec2(4,4);
  style.AntiAliasedLines        = true;
  style.AntiAliasedShapes       = true;
  style.CurveTessellationTol    = 1.25f;

  // default widget settings (imgui_widget.h)
  ImGuiStyleWidgets.DefaultColors();
  ImGuiStyleWidgets.DefaultStyle();

  // default UI settings (settings.h)
  ImGuiStyleUI.DefaultColors();
  ImGuiStyleUI.DefaultStyle();

  // default view settings (settings.h)
  view_wireframe = false;
  view_orthogonal = false;
  view_fov = 45.0f;
  view_resetdistance = 1.3f;
  view_bgrgb[0] = view_bgrgb[1] = view_bgrgb[2] = 0.f;
  view_bgrgb[3] = 1.0f;
  view_lightpos[0] = view_lightpos[1] = 20.f;
  view_lightpos[2] = 0.f;
  view_lightcolor[0] = view_lightcolor[1] = view_lightcolor[2] = 1.f;
  view_ambient = 0.2f;
  view_diffuse = 0.4f;
  view_specular = 0.6f;
  view_shininess = 8;
  view_isphres = 2;
  view_icylres = 0;
  view_mousesens_rot = 1.0f;
  view_mousesens_zoom = 1.0f;

  // set up shader and fill uniforms
  if (firstpass){
    shader = new Shader();
    shader->use();
  }
  shader->setVec3("lightPos",value_ptr(view_lightpos));
  shader->setVec3("lightColor",value_ptr(view_lightcolor));
  shader->setFloat("ambient",view_ambient);
  shader->setFloat("diffuse",view_diffuse);
  shader->setFloat("specular",view_specular);
  shader->setInt("shininess",view_shininess);

  ImGuiIO& io = GetIO();
  if (firstpass){
    // load standard fonts
    // const ImWchar* ImFontAtlas::GetGlyphRangesDefault()
    ImFontConfig fntconfig;
    fntconfig.MergeMode = false;
    fntconfig.PixelSnapH = true;
    fntconfig.OversampleH = fntconfig.OversampleV = 4;
    fntconfig.MergeMode = false;
    fntconfig.SizePixels = fontsizebake;

    static const ImWchar ranges[] = { 0x0020, 0x00FF, 0 }; // basic latin + latin supplement
    static const ImWchar rangesicon[] = { ICON_MIN_SM, ICON_MAX_SM, 0 }; // icons

    // Bake the first font - always the icons.
    strcpy(fntconfig.Name, "Icons");
    fonticon = io.Fonts->AddFontFromMemoryCompressedBase85TTF(smallicons_data_base85, fontsizebake, &fntconfig, rangesicon);

    // Bake the rest of the fonts
    strcpy(fntconfig.Name, "Proggy Clean");
    io.Fonts->AddFontFromMemoryCompressedBase85TTF(proggyclean_data_base85, fontsizebake, &fntconfig, ranges);
    strcpy(fntconfig.Name, "Cousine Regular");
    io.Fonts->AddFontFromMemoryCompressedBase85TTF(cousine_data_base85, fontsizebake, &fntconfig, ranges);
    strcpy(fntconfig.Name, "Droid Sans");
    io.Fonts->AddFontFromMemoryCompressedBase85TTF(droidsans_data_base85, fontsizebake, &fntconfig, ranges);
    strcpy(fntconfig.Name, "Karla Regular");
    io.Fonts->AddFontFromMemoryCompressedBase85TTF(karla_data_base85, fontsizebake, &fntconfig, ranges);
    strcpy(fntconfig.Name, "Proggy Tiny");
    io.Fonts->AddFontFromMemoryCompressedBase85TTF(proggytiny_data_base85, fontsizebake, &fntconfig, ranges);
    strcpy(fntconfig.Name, "Roboto Medium");
    io.Fonts->AddFontFromMemoryCompressedBase85TTF(roboto_data_base85, fontsizebake, &fntconfig, ranges);
  }

  // set the defualt sizes
  io.Fonts->Fonts[0]->Scale = ImGuiStyleUI.FontSizeIcon / fontsizebake;
  for (int i = 1; i < io.Fonts->Fonts.Size; i++)
    io.Fonts->Fonts[i]->Scale = ImGuiStyleUI.FontSize / fontsizebake;

  // set the default font
  fontdefault = io.Fonts->Fonts[1];
  io.FontDefault = fontdefault;

  firstpass = false;
}

// classic colors style
void UIStyleColorsClassic(){
  ImGuiStyle& style = GetStyle();
  style.Colors[ImGuiCol_Text]                   = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  style.Colors[ImGuiCol_TextDisabled]           = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  style.Colors[ImGuiCol_WindowBg]               = ImVec4(0.00f, 0.00f, 0.00f, 0.70f);
  style.Colors[ImGuiCol_ChildBg]                = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
  style.Colors[ImGuiCol_PopupBg]                = ImVec4(0.11f, 0.11f, 0.14f, 0.92f);
  style.Colors[ImGuiCol_Border]                 = ImVec4(0.50f, 0.50f, 0.50f, 0.50f);
  style.Colors[ImGuiCol_BorderShadow]           = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
  style.Colors[ImGuiCol_FrameBg]                = ImVec4(0.43f, 0.43f, 0.43f, 0.39f);
  style.Colors[ImGuiCol_FrameBgHovered]         = ImVec4(0.47f, 0.47f, 0.69f, 0.40f);
  style.Colors[ImGuiCol_FrameBgActive]          = ImVec4(0.42f, 0.41f, 0.64f, 0.69f);
  style.Colors[ImGuiCol_TitleBg]                = ImVec4(0.27f, 0.27f, 0.54f, 0.83f);
  style.Colors[ImGuiCol_TitleBgActive]          = ImVec4(0.32f, 0.32f, 0.63f, 0.87f);
  style.Colors[ImGuiCol_TitleBgCollapsed]       = ImVec4(0.40f, 0.40f, 0.80f, 0.20f);
  style.Colors[ImGuiCol_MenuBarBg]              = ImVec4(0.40f, 0.40f, 0.55f, 0.80f);
  style.Colors[ImGuiCol_ScrollbarBg]            = ImVec4(0.20f, 0.25f, 0.30f, 0.60f);
  style.Colors[ImGuiCol_ScrollbarGrab]          = ImVec4(0.40f, 0.40f, 0.80f, 0.30f);
  style.Colors[ImGuiCol_ScrollbarGrabHovered]   = ImVec4(0.40f, 0.40f, 0.80f, 0.40f);
  style.Colors[ImGuiCol_ScrollbarGrabActive]    = ImVec4(0.41f, 0.39f, 0.80f, 0.60f);
  style.Colors[ImGuiCol_CheckMark]              = ImVec4(0.90f, 0.90f, 0.90f, 0.50f);
  style.Colors[ImGuiCol_SliderGrab]             = ImVec4(1.00f, 1.00f, 1.00f, 0.30f);
  style.Colors[ImGuiCol_SliderGrabActive]       = ImVec4(0.41f, 0.39f, 0.80f, 0.60f);
  style.Colors[ImGuiCol_Button]                 = ImVec4(0.35f, 0.40f, 0.61f, 0.62f);
  style.Colors[ImGuiCol_ButtonHovered]          = ImVec4(0.40f, 0.48f, 0.71f, 0.79f);
  style.Colors[ImGuiCol_ButtonActive]           = ImVec4(0.46f, 0.54f, 0.80f, 1.00f);
  style.Colors[ImGuiCol_Header]                 = ImVec4(0.40f, 0.40f, 0.90f, 0.45f);
  style.Colors[ImGuiCol_HeaderHovered]          = ImVec4(0.45f, 0.45f, 0.90f, 0.80f);
  style.Colors[ImGuiCol_HeaderActive]           = ImVec4(0.53f, 0.53f, 0.87f, 0.80f);
  style.Colors[ImGuiCol_Separator]              = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
  style.Colors[ImGuiCol_SeparatorHovered]       = ImVec4(0.60f, 0.60f, 0.70f, 1.00f);
  style.Colors[ImGuiCol_SeparatorActive]        = ImVec4(0.70f, 0.70f, 0.90f, 1.00f);
  style.Colors[ImGuiCol_ResizeGrip]             = ImVec4(1.00f, 1.00f, 1.00f, 0.16f);
  style.Colors[ImGuiCol_ResizeGripHovered]      = ImVec4(0.78f, 0.82f, 1.00f, 0.60f);
  style.Colors[ImGuiCol_ResizeGripActive]       = ImVec4(0.78f, 0.82f, 1.00f, 0.90f);
  style.Colors[ImGuiCol_CloseButton]            = ImVec4(0.50f, 0.50f, 0.90f, 0.50f);
  style.Colors[ImGuiCol_CloseButtonHovered]     = ImVec4(0.70f, 0.70f, 0.90f, 0.60f);
  style.Colors[ImGuiCol_CloseButtonActive]      = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
  style.Colors[ImGuiCol_PlotLines]              = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
  style.Colors[ImGuiCol_PlotLinesHovered]       = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogram]          = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogramHovered]   = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TextSelectedBg]         = ImVec4(0.00f, 0.00f, 1.00f, 0.35f);
  style.Colors[ImGuiCol_ModalWindowDarkening]   = ImVec4(0.20f, 0.20f, 0.20f, 0.35f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Slidingbar]        = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarHovered] = ImVec4(0.60f, 0.60f, 0.70f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarActive]  = ImVec4(0.70f, 0.70f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Tab]               = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabHovered]        = ImVec4(0.45f, 0.45f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabPressed]        = ImVec4(0.46f, 0.54f, 0.80f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabActive]         = ImVec4(0.53f, 0.53f, 0.87f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFg]            = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgHovered]     = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgActive]      = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBg]            = ImVec4(0.80f, 0.20f, 0.00f, 0.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgHovered]     = ImVec4(0.80f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgActive]      = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabBorder]         = ImVec4(0.50f, 0.50f, 0.50f, 0.50f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGrip]          = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripHovered]   = ImVec4(0.80f, 0.40f, 0.20f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripActive]    = ImVec4(1.00f, 0.40f, 0.20f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTarget]        = ImVec4(0.43f, 0.43f, 0.43f, 0.39f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTargetActive]  = ImVec4(0.80f, 0.80f, 0.80f, 0.80f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIcon]         = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconHovered]  = ImVec4(0.8549f,0.6471f,0.1255f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconActive]   = ImVec4(0.7216f,0.5254,0.04314f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconInactive] = ImVec4(0.5f,0.5f,0.5f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageInfo] = ImVec4(0.0f,1.0f,0.0f,0.4f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageWarning] = ImVec4(1.0f,1.0f,0.0f,0.4f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageError] = ImVec4(1.0f,0.0f,0.0f,0.4f);
}

// dark colors style
void UIStyleColorsDark(){
  ImGuiStyle& style = GetStyle();
  style.Colors[ImGuiCol_Text]                   = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
  style.Colors[ImGuiCol_TextDisabled]           = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
  style.Colors[ImGuiCol_WindowBg]               = ImVec4(0.06f, 0.06f, 0.06f, 0.94f);
  style.Colors[ImGuiCol_ChildBg]                = ImVec4(1.00f, 1.00f, 1.00f, 0.00f);
  style.Colors[ImGuiCol_PopupBg]                = ImVec4(0.08f, 0.08f, 0.08f, 0.94f);
  style.Colors[ImGuiCol_Border]                 = ImVec4(0.43f, 0.43f, 0.50f, 0.50f);
  style.Colors[ImGuiCol_BorderShadow]           = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
  style.Colors[ImGuiCol_FrameBg]                = ImVec4(0.16f, 0.29f, 0.48f, 0.54f);
  style.Colors[ImGuiCol_FrameBgHovered]         = ImVec4(0.26f, 0.59f, 0.98f, 0.40f);
  style.Colors[ImGuiCol_FrameBgActive]          = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
  style.Colors[ImGuiCol_TitleBg]                = ImVec4(0.04f, 0.04f, 0.04f, 1.00f);
  style.Colors[ImGuiCol_TitleBgActive]          = ImVec4(0.16f, 0.29f, 0.48f, 1.00f);
  style.Colors[ImGuiCol_TitleBgCollapsed]       = ImVec4(0.00f, 0.00f, 0.00f, 0.51f);
  style.Colors[ImGuiCol_MenuBarBg]              = ImVec4(0.14f, 0.14f, 0.14f, 1.00f);
  style.Colors[ImGuiCol_ScrollbarBg]            = ImVec4(0.02f, 0.02f, 0.02f, 0.53f);
  style.Colors[ImGuiCol_ScrollbarGrab]          = ImVec4(0.31f, 0.31f, 0.31f, 1.00f);
  style.Colors[ImGuiCol_ScrollbarGrabHovered]   = ImVec4(0.41f, 0.41f, 0.41f, 1.00f);
  style.Colors[ImGuiCol_ScrollbarGrabActive]    = ImVec4(0.51f, 0.51f, 0.51f, 1.00f);
  style.Colors[ImGuiCol_CheckMark]              = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_SliderGrab]             = ImVec4(0.24f, 0.52f, 0.88f, 1.00f);
  style.Colors[ImGuiCol_SliderGrabActive]       = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_Button]                 = ImVec4(0.26f, 0.59f, 0.98f, 0.40f);
  style.Colors[ImGuiCol_ButtonHovered]          = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_ButtonActive]           = ImVec4(0.06f, 0.53f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_Header]                 = ImVec4(0.26f, 0.59f, 0.98f, 0.31f);
  style.Colors[ImGuiCol_HeaderHovered]          = ImVec4(0.26f, 0.59f, 0.98f, 0.80f);
  style.Colors[ImGuiCol_HeaderActive]           = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_Separator]              = ImVec4(0.61f, 0.61f, 0.61f, 1.00f);
  style.Colors[ImGuiCol_SeparatorHovered]       = ImVec4(0.10f, 0.40f, 0.75f, 0.78f);
  style.Colors[ImGuiCol_SeparatorActive]        = ImVec4(0.10f, 0.40f, 0.75f, 1.00f);
  style.Colors[ImGuiCol_ResizeGrip]             = ImVec4(0.26f, 0.59f, 0.98f, 0.25f);
  style.Colors[ImGuiCol_ResizeGripHovered]      = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
  style.Colors[ImGuiCol_ResizeGripActive]       = ImVec4(0.26f, 0.59f, 0.98f, 0.95f);
  style.Colors[ImGuiCol_CloseButton]            = ImVec4(0.41f, 0.41f, 0.41f, 0.50f);
  style.Colors[ImGuiCol_CloseButtonHovered]     = ImVec4(0.98f, 0.39f, 0.36f, 1.00f);
  style.Colors[ImGuiCol_CloseButtonActive]      = ImVec4(0.98f, 0.39f, 0.36f, 1.00f);
  style.Colors[ImGuiCol_PlotLines]              = ImVec4(0.61f, 0.61f, 0.61f, 1.00f);
  style.Colors[ImGuiCol_PlotLinesHovered]       = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogram]          = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogramHovered]   = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TextSelectedBg]         = ImVec4(0.26f, 0.59f, 0.98f, 0.35f);
  style.Colors[ImGuiCol_ModalWindowDarkening]   = ImVec4(0.80f, 0.80f, 0.80f, 0.35f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Slidingbar]        = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarHovered] = ImVec4(0.60f, 0.60f, 0.70f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarActive]  = ImVec4(0.70f, 0.70f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Tab]               = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabHovered]        = ImVec4(0.45f, 0.45f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabPressed]        = ImVec4(0.46f, 0.54f, 0.80f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabActive]         = ImVec4(0.53f, 0.53f, 0.87f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFg]            = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgHovered]     = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgActive]      = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBg]            = ImVec4(0.80f, 0.20f, 0.00f, 0.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgHovered]     = ImVec4(0.80f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgActive]      = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabBorder]         = ImVec4(0.50f, 0.50f, 0.50f, 0.50f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGrip]          = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripHovered]   = ImVec4(0.80f, 0.40f, 0.20f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripActive]    = ImVec4(1.00f, 0.40f, 0.20f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTarget]        = ImVec4(0.43f, 0.43f, 0.43f, 0.39f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTargetActive]  = ImVec4(0.80f, 0.80f, 0.80f, 0.80f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIcon]         = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconHovered]  = ImVec4(0.8549f,0.6471f,0.1255f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconActive]   = ImVec4(0.7216f,0.5254,0.04314f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconInactive] = ImVec4(0.5f,0.5f,0.5f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageInfo] = ImVec4(0.0f,1.0f,0.0f,0.4f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageWarning] = ImVec4(1.0f,1.0f,0.0f,0.4f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageError] = ImVec4(1.0f,0.0f,0.0f,0.4f);
}

// light colors style
void UIStyleColorsLight(){
  ImGuiStyle& style = GetStyle();
  style.Colors[ImGuiCol_Text]                   = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TextDisabled]           = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  style.Colors[ImGuiCol_WindowBg]               = ImVec4(0.94f, 0.94f, 0.94f, 1.00f);
  style.Colors[ImGuiCol_ChildBg]                = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
  style.Colors[ImGuiCol_PopupBg]                = ImVec4(1.00f, 1.00f, 1.00f, 0.98f);
  style.Colors[ImGuiCol_Border]                 = ImVec4(0.00f, 0.00f, 0.00f, 0.30f);
  style.Colors[ImGuiCol_BorderShadow]           = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
  style.Colors[ImGuiCol_FrameBg]                = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
  style.Colors[ImGuiCol_FrameBgHovered]         = ImVec4(0.26f, 0.59f, 0.98f, 0.40f);
  style.Colors[ImGuiCol_FrameBgActive]          = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
  style.Colors[ImGuiCol_TitleBg]                = ImVec4(0.96f, 0.96f, 0.96f, 1.00f);
  style.Colors[ImGuiCol_TitleBgActive]          = ImVec4(0.82f, 0.82f, 0.82f, 1.00f);
  style.Colors[ImGuiCol_TitleBgCollapsed]       = ImVec4(1.00f, 1.00f, 1.00f, 0.51f);
  style.Colors[ImGuiCol_MenuBarBg]              = ImVec4(0.86f, 0.86f, 0.86f, 1.00f);
  style.Colors[ImGuiCol_ScrollbarBg]            = ImVec4(0.98f, 0.98f, 0.98f, 0.53f);
  style.Colors[ImGuiCol_ScrollbarGrab]          = ImVec4(0.69f, 0.69f, 0.69f, 0.80f);
  style.Colors[ImGuiCol_ScrollbarGrabHovered]   = ImVec4(0.49f, 0.49f, 0.49f, 0.80f);
  style.Colors[ImGuiCol_ScrollbarGrabActive]    = ImVec4(0.49f, 0.49f, 0.49f, 1.00f);
  style.Colors[ImGuiCol_CheckMark]              = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_SliderGrab]             = ImVec4(0.26f, 0.59f, 0.98f, 0.78f);
  style.Colors[ImGuiCol_SliderGrabActive]       = ImVec4(0.46f, 0.54f, 0.80f, 0.60f);
  style.Colors[ImGuiCol_Button]                 = ImVec4(0.26f, 0.59f, 0.98f, 0.40f);
  style.Colors[ImGuiCol_ButtonHovered]          = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_ButtonActive]           = ImVec4(0.06f, 0.53f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_Header]                 = ImVec4(0.26f, 0.59f, 0.98f, 0.31f);
  style.Colors[ImGuiCol_HeaderHovered]          = ImVec4(0.26f, 0.59f, 0.98f, 0.80f);
  style.Colors[ImGuiCol_HeaderActive]           = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
  style.Colors[ImGuiCol_Separator]              = ImVec4(0.39f, 0.39f, 0.39f, 1.00f);
  style.Colors[ImGuiCol_SeparatorHovered]       = ImVec4(0.14f, 0.44f, 0.80f, 0.78f);
  style.Colors[ImGuiCol_SeparatorActive]        = ImVec4(0.14f, 0.44f, 0.80f, 1.00f);
  style.Colors[ImGuiCol_ResizeGrip]             = ImVec4(0.80f, 0.80f, 0.80f, 0.56f);
  style.Colors[ImGuiCol_ResizeGripHovered]      = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
  style.Colors[ImGuiCol_ResizeGripActive]       = ImVec4(0.26f, 0.59f, 0.98f, 0.95f);
  style.Colors[ImGuiCol_CloseButton]            = ImVec4(0.59f, 0.59f, 0.59f, 0.50f);
  style.Colors[ImGuiCol_CloseButtonHovered]     = ImVec4(0.98f, 0.39f, 0.36f, 1.00f);
  style.Colors[ImGuiCol_CloseButtonActive]      = ImVec4(0.98f, 0.39f, 0.36f, 1.00f);
  style.Colors[ImGuiCol_PlotLines]              = ImVec4(0.39f, 0.39f, 0.39f, 1.00f);
  style.Colors[ImGuiCol_PlotLinesHovered]       = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogram]          = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogramHovered]   = ImVec4(1.00f, 0.45f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TextSelectedBg]         = ImVec4(0.26f, 0.59f, 0.98f, 0.35f);
  style.Colors[ImGuiCol_ModalWindowDarkening]   = ImVec4(0.20f, 0.20f, 0.20f, 0.35f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Slidingbar]        = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarHovered] = ImVec4(0.60f, 0.60f, 0.70f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarActive]  = ImVec4(0.70f, 0.70f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Tab]               = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabHovered]        = ImVec4(0.45f, 0.45f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabPressed]        = ImVec4(0.46f, 0.54f, 0.80f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabActive]         = ImVec4(0.53f, 0.53f, 0.87f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFg]            = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgHovered]     = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgActive]      = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBg]            = ImVec4(0.80f, 0.20f, 0.00f, 0.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgHovered]     = ImVec4(0.80f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgActive]      = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabBorder]         = ImVec4(0.50f, 0.50f, 0.50f, 0.50f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGrip]          = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripHovered]   = ImVec4(0.80f, 0.40f, 0.20f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripActive]    = ImVec4(1.00f, 0.40f, 0.20f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTarget]        = ImVec4(0.43f, 0.43f, 0.43f, 0.39f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTargetActive]  = ImVec4(0.80f, 0.80f, 0.80f, 0.80f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIcon]         = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconHovered]  = ImVec4(0.8549f,0.6471f,0.1255f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconActive]   = ImVec4(0.7216f,0.5254,0.04314f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconInactive] = ImVec4(0.5f,0.5f,0.5f,1.0f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageInfo] = ImVec4(0.0f,1.0f,0.0f,0.4f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageWarning] = ImVec4(1.0f,1.0f,0.0f,0.4f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageError] = ImVec4(1.0f,0.0f,0.0f,0.4f);
}

// Find the configuration file or create a new one if it does not
// exist.
bool FindConfigurationFile(){
  // __MORPHOS__, _WIN32, __HAIKU__, __APPLE__

  // Find the path
  const char* xdghome = getenv("XDG_CONFIG_HOME");
  const char* home = getenv("HOME");
  string thisfile = {};

  const int npath = 4;
  string paths[npath] = {
    xdghome?string(xdghome)  + "critic2":"",
    home?string(home) + "/.config/critic2":"",
    home?string(home) + "/.critic2":"",
    "."
  };

  for (int i=0; i<npath; i++){
    if (paths[i].length() == 0)
      continue;

    // check if the directory exists
    bool isdir = (access(paths[i].c_str(),F_OK) == 0);

    // try to create directory
    bool dircreated = false;
    if (!isdir){
      mode_t oldmask = umask(0);
      int ok = mkdir(paths[i].c_str(), S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
      umask(oldmask);
      if (ok != 0)
	continue;
      dircreated = true;
    }

    // check if directory is readable and writable
    if (access(paths[i].c_str(),R_OK) != 0 || access(paths[i].c_str(),W_OK) != 0)
      goto cleanup;

    thisfile = paths[i] + "/critic2.cfg";

    // check if file exists
    if (access(thisfile.c_str(),F_OK) != 0)
      if (!WriteConfigurationFile(thisfile))
	goto cleanup;
      
    // check if the file is readable and writable
    if (access(thisfile.c_str(),R_OK) != 0 || access(thisfile.c_str(),W_OK) != 0)
      goto cleanup;

    conffile = thisfile;
    return true;

    cleanup:
    if (dircreated)
      rmdir(paths[i].c_str());
  }
  return false;
}


// JSON input and output //
static void to_json(json& j, const ImVec4& p) {
  j = json{p.x, p.y, p.z, p.w};
}
static void to_json(json& j, const ImVec2& p) {
  j = json{p.x, p.y};
}
static void from_json(const json& j, ImVec4& p) {
  p.x = j[0];
  p.y = j[1];
  p.z = j[2];
  p.w = j[3];
}
static void from_json(const json& j, ImVec2& p) {
  p.x = j[0];
  p.y = j[1];
}

bool WriteConfigurationFile(string file){
  ofstream fo(file);
  if (!fo.good())
    return false;

  json j;

  // UI colors (settings.h)
  j["Color_ViewIcon"] = ImGuiStyleUI.Colors[ImGuiColUI_ViewIcon];
  j["Color_ViewIconHovered"] = ImGuiStyleUI.Colors[ImGuiColUI_ViewIconHovered];
  j["Color_ViewIconActive"] = ImGuiStyleUI.Colors[ImGuiColUI_ViewIconActive];
  j["Color_ViewIconInactive"] = ImGuiStyleUI.Colors[ImGuiColUI_ViewIconInactive];
  j["Color_MessageInfo"] = ImGuiStyleUI.Colors[ImGuiColUI_MessageInfo];
  j["Color_MessageWarning"] = ImGuiStyleUI.Colors[ImGuiColUI_MessageWarning];
  j["Color_MessageError"] = ImGuiStyleUI.Colors[ImGuiColUI_MessageError];

  // UI styles (settings.h)
  j["MessageWidth"] = ImGuiStyleUI.MessageWidth;
  j["MessageExpire"] = ImGuiStyleUI.MessageExpire;
  j["FontSizeIcon"] = ImGuiStyleUI.FontSizeIcon;
  j["FontSize"] = ImGuiStyleUI.FontSize;
  j["TooltipEnabled"] = ImGuiStyleUI.TooltipEnabled;
  j["TooltipDelay"] = ImGuiStyleUI.TooltipDelay;
  j["TooltipMaxwidth"] = ImGuiStyleUI.TooltipMaxwidth;

  // Widget colors (imgui_widgets.h)
  j["Color_Slidingbar"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_Slidingbar];
  j["Color_SlidingbarHovered"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarHovered];
  j["Color_SlidingbarActive"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarActive];
  j["Color_Tab"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_Tab];
  j["Color_TabHovered"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabHovered];
  j["Color_TabPressed"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabPressed];
  j["Color_TabActive"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabActive];
  j["Color_TabXFg"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFg];
  j["Color_TabXFgHovered"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgHovered];
  j["Color_TabXFgActive"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgActive];
  j["Color_TabXBg"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBg];
  j["Color_TabXBgHovered"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgHovered];
  j["Color_TabXBgActive"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgActive];
  j["Color_TabBorder"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabBorder];
  j["Color_DropTarget"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTarget];
  j["Color_DropTargetActive"] = ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTargetActive];

  // Widget styles (imgui_widgets.h)
  j["TabRounding"] = ImGuiStyleWidgets.TabRounding;
  j["TabBorderSize"] = ImGuiStyleWidgets.TabBorderSize;
  j["DropTargetLooseness"] = ImGuiStyleWidgets.DropTargetLooseness;
  j["DropTargetMinsizeEdge"] = ImGuiStyleWidgets.DropTargetMinsizeEdge;
  j["DropTargetMaxsizeEdge"] = ImGuiStyleWidgets.DropTargetMaxsizeEdge;
  j["DropTargetEdgeFraction"] = ImGuiStyleWidgets.DropTargetEdgeFraction;
  j["DropTargetFullFraction"] = ImGuiStyleWidgets.DropTargetFullFraction;
  j["TabHeight"] = ImGuiStyleWidgets.TabHeight;
  j["TabMaxWidth"] = ImGuiStyleWidgets.TabMaxWidth;
  j["SlidingBarWidth"] = ImGuiStyleWidgets.SlidingBarWidth;

  // imgui colors
  ImGuiStyle& style = GetStyle();
  j["Color_Text"] = style.Colors[ImGuiCol_Text];
  j["Color_TextDisabled"] = style.Colors[ImGuiCol_TextDisabled];
  j["Color_WindowBg"] = style.Colors[ImGuiCol_WindowBg];
  j["Color_ChildBg"] = style.Colors[ImGuiCol_ChildBg];
  j["Color_PopupBg"] = style.Colors[ImGuiCol_PopupBg];
  j["Color_Border"] = style.Colors[ImGuiCol_Border];
  j["Color_BorderShadow"] = style.Colors[ImGuiCol_BorderShadow];
  j["Color_FrameBg"] = style.Colors[ImGuiCol_FrameBg];
  j["Color_FrameBgHovered"] = style.Colors[ImGuiCol_FrameBgHovered];
  j["Color_FrameBgActive"] = style.Colors[ImGuiCol_FrameBgActive];
  j["Color_TitleBg"] = style.Colors[ImGuiCol_TitleBg];
  j["Color_TitleBgActive"] = style.Colors[ImGuiCol_TitleBgActive];
  j["Color_TitleBgCollapsed"] = style.Colors[ImGuiCol_TitleBgCollapsed];
  j["Color_MenuBarBg"] = style.Colors[ImGuiCol_MenuBarBg];
  j["Color_ScrollbarBg"] = style.Colors[ImGuiCol_ScrollbarBg];
  j["Color_ScrollbarGrab"] = style.Colors[ImGuiCol_ScrollbarGrab];
  j["Color_ScrollbarGrabHovered"] = style.Colors[ImGuiCol_ScrollbarGrabHovered];
  j["Color_ScrollbarGrabActive"] = style.Colors[ImGuiCol_ScrollbarGrabActive];
  j["Color_CheckMark"] = style.Colors[ImGuiCol_CheckMark];
  j["Color_SliderGrab"] = style.Colors[ImGuiCol_SliderGrab];
  j["Color_SliderGrabActive"] = style.Colors[ImGuiCol_SliderGrabActive];
  j["Color_Button"] = style.Colors[ImGuiCol_Button];
  j["Color_ButtonHovered"] = style.Colors[ImGuiCol_ButtonHovered];
  j["Color_ButtonActive"] = style.Colors[ImGuiCol_ButtonActive];
  j["Color_Header"] = style.Colors[ImGuiCol_Header];
  j["Color_HeaderHovered"] = style.Colors[ImGuiCol_HeaderHovered];
  j["Color_HeaderActive"] = style.Colors[ImGuiCol_HeaderActive];
  j["Color_Separator"] = style.Colors[ImGuiCol_Separator];
  j["Color_SeparatorHovered"] = style.Colors[ImGuiCol_SeparatorHovered];
  j["Color_SeparatorActive"] = style.Colors[ImGuiCol_SeparatorActive];
  j["Color_ResizeGrip"] = style.Colors[ImGuiCol_ResizeGrip];
  j["Color_ResizeGripHovered"] = style.Colors[ImGuiCol_ResizeGripHovered];
  j["Color_ResizeGripActive"] = style.Colors[ImGuiCol_ResizeGripActive];
  j["Color_CloseButton"] = style.Colors[ImGuiCol_CloseButton];
  j["Color_CloseButtonHovered"] = style.Colors[ImGuiCol_CloseButtonHovered];
  j["Color_CloseButtonActive"] = style.Colors[ImGuiCol_CloseButtonActive];
  j["Color_PlotLines"] = style.Colors[ImGuiCol_PlotLines];
  j["Color_PlotLinesHovered"] = style.Colors[ImGuiCol_PlotLinesHovered];
  j["Color_PlotHistogram"] = style.Colors[ImGuiCol_PlotHistogram];
  j["Color_PlotHistogramHovered"] = style.Colors[ImGuiCol_PlotHistogramHovered];
  j["Color_TextSelectedBg"] = style.Colors[ImGuiCol_TextSelectedBg];
  j["Color_ModalWindowDarkening"] = style.Colors[ImGuiCol_ModalWindowDarkening];

  // imgui styles
  j["Alpha"] = style.Alpha;
  j["WindowPadding"] = style.WindowPadding;
  j["WindowRounding"] = style.WindowRounding;
  j["WindowBorderSize"] = style.WindowBorderSize;
  j["WindowMinSize"] = style.WindowMinSize;
  j["WindowTitleAlign"] = style.WindowTitleAlign;
  j["ChildRounding"] = style.ChildRounding;
  j["ChildBorderSize"] = style.ChildBorderSize;
  j["PopupRounding"] = style.PopupRounding;
  j["PopupBorderSize"] = style.PopupBorderSize;
  j["FramePadding"] = style.FramePadding;
  j["FrameRounding"] = style.FrameRounding;
  j["FrameBorderSize"] = style.FrameBorderSize;
  j["ItemSpacing"] = style.ItemSpacing;
  j["ItemInnerSpacing"] = style.ItemInnerSpacing;
  j["TouchExtraPadding"] = style.TouchExtraPadding;
  j["IndentSpacing"] = style.IndentSpacing;
  j["ColumnsMinSpacing"] = style.ColumnsMinSpacing;
  j["ScrollbarSize"] = style.ScrollbarSize;
  j["ScrollbarRounding"] = style.ScrollbarRounding;
  j["GrabMinSize"] = style.GrabMinSize;
  j["GrabRounding"] = style.GrabRounding;
  j["ButtonTextAlign"] = style.ButtonTextAlign;
  j["DisplayWindowPadding"] = style.DisplayWindowPadding;
  j["DisplaySafeAreaPadding"] = style.DisplaySafeAreaPadding;
  j["AntiAliasedLines"] = style.AntiAliasedLines;
  j["AntiAliasedShapes"] = style.AntiAliasedShapes;
  j["CurveTessellationTol"] = style.CurveTessellationTol;

  if (!fo.is_open())
    return false;

  fo << setw(4) << j << endl;
  fo.close();
  return true;
}

bool ReadConfigurationFile(string file){
  ifstream fi(file);
  if (!fi.good() || !fi.is_open())
    return false;

  json j;
  fi >> j;

  // UI colors (settings.h)
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIcon] = j["Color_ViewIcon"].get<ImVec4>();
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconHovered] = j["Color_ViewIconHovered"].get<ImVec4>();
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconActive] = j["Color_ViewIconActive"].get<ImVec4>();
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconInactive] = j["Color_ViewIconInactive"].get<ImVec4>();
  ImGuiStyleUI.Colors[ImGuiColUI_MessageInfo] = j["Color_MessageInfo"].get<ImVec4>();
  ImGuiStyleUI.Colors[ImGuiColUI_MessageWarning] = j["Color_MessageWarning"].get<ImVec4>();
  ImGuiStyleUI.Colors[ImGuiColUI_MessageError] = j["Color_MessageError"].get<ImVec4>();

  // UI styles (settings.h)
  ImGuiStyleUI.MessageWidth = j["MessageWidth"].get<float>();
  ImGuiStyleUI.MessageExpire = j["MessageExpire"].get<float>();
  ImGuiStyleUI.FontSizeIcon = j["FontSizeIcon"].get<float>();
  ImGuiStyleUI.FontSize = j["FontSize"].get<float>();
  ImGuiStyleUI.TooltipEnabled = j["TooltipEnabled"].get<bool>();
  ImGuiStyleUI.TooltipDelay = j["TooltipDelay"].get<float>();
  ImGuiStyleUI.TooltipMaxwidth = j["TooltipMaxwidth"].get<float>();

  // Widget colors (imgui_widgets.h)
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Slidingbar] = j["Color_Slidingbar"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarHovered] = j["Color_SlidingbarHovered"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarActive] = j["Color_SlidingbarActive"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Tab] = j["Color_Tab"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabHovered] = j["Color_TabHovered"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabPressed] = j["Color_TabPressed"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabActive] = j["Color_TabActive"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFg] = j["Color_TabXFg"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgHovered] = j["Color_TabXFgHovered"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgActive] = j["Color_TabXFgActive"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBg] = j["Color_TabXBg"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgHovered] = j["Color_TabXBgHovered"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgActive] = j["Color_TabXBgActive"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabBorder] = j["Color_TabBorder"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTarget] = j["Color_DropTarget"].get<ImVec4>();
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTargetActive] = j["Color_DropTargetActive"].get<ImVec4>();

  // Widget styles (imgui_widgets.h)
  ImGuiStyleWidgets.TabRounding = j["TabRounding"].get<float>();
  ImGuiStyleWidgets.TabBorderSize = j["TabBorderSize"].get<float>();
  ImGuiStyleWidgets.DropTargetLooseness = j["DropTargetLooseness"].get<float>();
  ImGuiStyleWidgets.DropTargetMinsizeEdge = j["DropTargetMinsizeEdge"].get<float>();
  ImGuiStyleWidgets.DropTargetMaxsizeEdge = j["DropTargetMaxsizeEdge"].get<float>();
  ImGuiStyleWidgets.DropTargetEdgeFraction = j["DropTargetEdgeFraction"].get<float>();
  ImGuiStyleWidgets.DropTargetFullFraction = j["DropTargetFullFraction"].get<float>();
  ImGuiStyleWidgets.TabHeight = j["TabHeight"].get<float>();
  ImGuiStyleWidgets.TabMaxWidth = j["TabMaxWidth"].get<float>();
  ImGuiStyleWidgets.SlidingBarWidth = j["SlidingBarWidth"].get<float>();

  // imgui colors
  ImGuiStyle& style = GetStyle();
  style.Colors[ImGuiCol_Text] = j["Color_Text"].get<ImVec4>();
  style.Colors[ImGuiCol_TextDisabled] = j["Color_TextDisabled"].get<ImVec4>();
  style.Colors[ImGuiCol_WindowBg] = j["Color_WindowBg"].get<ImVec4>();
  style.Colors[ImGuiCol_ChildBg] = j["Color_ChildBg"].get<ImVec4>();
  style.Colors[ImGuiCol_PopupBg] = j["Color_PopupBg"].get<ImVec4>();
  style.Colors[ImGuiCol_Border] = j["Color_Border"].get<ImVec4>();
  style.Colors[ImGuiCol_BorderShadow] = j["Color_BorderShadow"].get<ImVec4>();
  style.Colors[ImGuiCol_FrameBg] = j["Color_FrameBg"].get<ImVec4>();
  style.Colors[ImGuiCol_FrameBgHovered] = j["Color_FrameBgHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_FrameBgActive] = j["Color_FrameBgActive"].get<ImVec4>();
  style.Colors[ImGuiCol_TitleBg] = j["Color_TitleBg"].get<ImVec4>();
  style.Colors[ImGuiCol_TitleBgActive] = j["Color_TitleBgActive"].get<ImVec4>();
  style.Colors[ImGuiCol_TitleBgCollapsed] = j["Color_TitleBgCollapsed"].get<ImVec4>();
  style.Colors[ImGuiCol_MenuBarBg] = j["Color_MenuBarBg"].get<ImVec4>();
  style.Colors[ImGuiCol_ScrollbarBg] = j["Color_ScrollbarBg"].get<ImVec4>();
  style.Colors[ImGuiCol_ScrollbarGrab] = j["Color_ScrollbarGrab"].get<ImVec4>();
  style.Colors[ImGuiCol_ScrollbarGrabHovered] = j["Color_ScrollbarGrabHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_ScrollbarGrabActive] = j["Color_ScrollbarGrabActive"].get<ImVec4>();
  style.Colors[ImGuiCol_CheckMark] = j["Color_CheckMark"].get<ImVec4>();
  style.Colors[ImGuiCol_SliderGrab] = j["Color_SliderGrab"].get<ImVec4>();
  style.Colors[ImGuiCol_SliderGrabActive] = j["Color_SliderGrabActive"].get<ImVec4>();
  style.Colors[ImGuiCol_Button] = j["Color_Button"].get<ImVec4>();
  style.Colors[ImGuiCol_ButtonHovered] = j["Color_ButtonHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_ButtonActive] = j["Color_ButtonActive"].get<ImVec4>();
  style.Colors[ImGuiCol_Header] = j["Color_Header"].get<ImVec4>();
  style.Colors[ImGuiCol_HeaderHovered] = j["Color_HeaderHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_HeaderActive] = j["Color_HeaderActive"].get<ImVec4>();
  style.Colors[ImGuiCol_Separator] = j["Color_Separator"].get<ImVec4>();
  style.Colors[ImGuiCol_SeparatorHovered] = j["Color_SeparatorHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_SeparatorActive] = j["Color_SeparatorActive"].get<ImVec4>();
  style.Colors[ImGuiCol_ResizeGrip] = j["Color_ResizeGrip"].get<ImVec4>();
  style.Colors[ImGuiCol_ResizeGripHovered] = j["Color_ResizeGripHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_ResizeGripActive] = j["Color_ResizeGripActive"].get<ImVec4>();
  style.Colors[ImGuiCol_CloseButton] = j["Color_CloseButton"].get<ImVec4>();
  style.Colors[ImGuiCol_CloseButtonHovered] = j["Color_CloseButtonHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_CloseButtonActive] = j["Color_CloseButtonActive"].get<ImVec4>();
  style.Colors[ImGuiCol_PlotLines] = j["Color_PlotLines"].get<ImVec4>();
  style.Colors[ImGuiCol_PlotLinesHovered] = j["Color_PlotLinesHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_PlotHistogram] = j["Color_PlotHistogram"].get<ImVec4>();
  style.Colors[ImGuiCol_PlotHistogramHovered] = j["Color_PlotHistogramHovered"].get<ImVec4>();
  style.Colors[ImGuiCol_TextSelectedBg] = j["Color_TextSelectedBg"].get<ImVec4>();
  style.Colors[ImGuiCol_ModalWindowDarkening] = j["Color_ModalWindowDarkening"].get<ImVec4>();

  // imgui styles
  style.Alpha = j["Alpha"].get<float>();
  style.WindowPadding = j["WindowPadding"].get<ImVec2>();
  style.WindowRounding = j["WindowRounding"].get<float>();
  style.WindowBorderSize = j["WindowBorderSize"].get<float>();
  style.WindowMinSize = j["WindowMinSize"].get<ImVec2>();
  style.WindowTitleAlign = j["WindowTitleAlign"].get<ImVec2>();
  style.ChildRounding = j["ChildRounding"].get<float>();
  style.ChildBorderSize = j["ChildBorderSize"].get<float>();
  style.PopupRounding = j["PopupRounding"].get<float>();
  style.PopupBorderSize = j["PopupBorderSize"].get<float>();
  style.FramePadding = j["FramePadding"].get<ImVec2>();
  style.FrameRounding = j["FrameRounding"].get<float>();
  style.FrameBorderSize = j["FrameBorderSize"].get<float>();
  style.ItemSpacing = j["ItemSpacing"].get<ImVec2>();
  style.ItemInnerSpacing = j["ItemInnerSpacing"].get<ImVec2>();
  style.TouchExtraPadding = j["TouchExtraPadding"].get<ImVec2>();
  style.IndentSpacing = j["IndentSpacing"].get<float>();
  style.ColumnsMinSpacing = j["ColumnsMinSpacing"].get<float>();
  style.ScrollbarSize = j["ScrollbarSize"].get<float>();
  style.ScrollbarRounding = j["ScrollbarRounding"].get<float>();
  style.GrabMinSize = j["GrabMinSize"].get<float>();
  style.GrabRounding = j["GrabRounding"].get<float>();
  style.ButtonTextAlign = j["ButtonTextAlign"].get<ImVec2>();
  style.DisplayWindowPadding = j["DisplayWindowPadding"].get<ImVec2>();
  style.DisplaySafeAreaPadding = j["DisplaySafeAreaPadding"].get<ImVec2>();
  style.AntiAliasedLines = j["AntiAliasedLines"].get<bool>();
  style.AntiAliasedShapes = j["AntiAliasedShapes"].get<bool>();
  style.CurveTessellationTol = j["CurveTessellationTol"].get<float>();

  fi.close();
  return true;
}

// Some global callbacks //

void error_callback(int error, const char* description){
  fprintf(stderr, "Error %d: %s\n", error, description);
}
