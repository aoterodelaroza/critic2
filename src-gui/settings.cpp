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
#include "keybinding.h"
#include <GLFW/glfw3.h>

#include "imgui/additional_fonts.h"
#include "imgui/imgui_widgets.h"
#include "imgui/imgui_dock.h"
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
const float fontsizebake = 48.0f;

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
bool view_show_atoms; // show atoms?
float view_scale_atoms; // global scale atoms
int view_isphres; // resolution of the atom spheres
bool view_show_bonds; // show bonds?
float view_scale_bonds; // global scale bonds
int view_icylres; // resolution of the bond cylinders
float view_mousesens_rot; // mouse rotation sensitivity (scale factor)
float view_mousesens_zoom; // mouse zoom sensitivity (scale factor)

// configuration file name
string conffile = "";

// List of configuration variables and types
enum vartype_ {Type_None,Type_Bool,Type_Int,Type_Float,Type_ImVec2,Type_ImVec4};

struct confvar_ {
  char *name;
  void *variable;
  vartype_ vartype;
};
static confvar_ confvar[] = {
  // UI colors (settings.h)
  {"Color_BackDrop",&ImGuiStyleUI.Colors[ImGuiColUI_BackDrop],Type_ImVec4},
  {"Color_ViewIcon",&ImGuiStyleUI.Colors[ImGuiColUI_ViewIcon],Type_ImVec4},
  {"Color_ViewIconHovered",&ImGuiStyleUI.Colors[ImGuiColUI_ViewIconHovered],Type_ImVec4},
  {"Color_ViewIconActive",&ImGuiStyleUI.Colors[ImGuiColUI_ViewIconActive],Type_ImVec4},
  {"Color_ViewIconInactive",&ImGuiStyleUI.Colors[ImGuiColUI_ViewIconInactive],Type_ImVec4},
  {"Color_MessageInfo",&ImGuiStyleUI.Colors[ImGuiColUI_MessageInfo],Type_ImVec4},
  {"Color_MessageWarning",&ImGuiStyleUI.Colors[ImGuiColUI_MessageWarning],Type_ImVec4},
  {"Color_MessageError",&ImGuiStyleUI.Colors[ImGuiColUI_MessageError],Type_ImVec4},

  // UI styles (settings.h)
  {"MessageWidth",&ImGuiStyleUI.MessageWidth,Type_Float},
  {"MessageExpire",&ImGuiStyleUI.MessageExpire,Type_Float},
  {"FontSizeIcon",&ImGuiStyleUI.FontSizeIcon,Type_Float},
  {"FontSize",&ImGuiStyleUI.FontSize,Type_Float},
  {"FontSelected",&ImGuiStyleUI.FontSelected,Type_Int},
  {"TooltipEnabled",&ImGuiStyleUI.TooltipEnabled,Type_Bool},
  {"TooltipDelay",&ImGuiStyleUI.TooltipDelay,Type_Float},
  {"TooltipMaxwidth",&ImGuiStyleUI.TooltipMaxwidth,Type_Float},
  {"UIScaleFactor",&ImGuiStyleUI.UIScaleFactor,Type_Float},

  // View default options (settings.h)
  {"ViewWireframe",&view_wireframe,Type_Bool},
  {"ViewOrthogonal",&view_orthogonal,Type_Bool},
  {"ViewFov",&view_fov,Type_Float},
  {"ViewResetdistance",&view_resetdistance,Type_Float},
  {"ViewBgrgb[0]",&view_bgrgb[0],Type_Float},
  {"ViewBgrgb[1]",&view_bgrgb[1],Type_Float},
  {"ViewBgrgb[2]",&view_bgrgb[2],Type_Float},
  {"ViewBgrgb[3]",&view_bgrgb[3],Type_Float},
  {"ViewLightpos.x",&view_lightpos.x,Type_Float},
  {"ViewLightpos.y",&view_lightpos.y,Type_Float},
  {"ViewLightpos.z",&view_lightpos.z,Type_Float},
  {"ViewLightcolor.x",&view_lightcolor.x,Type_Float},
  {"ViewLightcolor.y",&view_lightcolor.y,Type_Float},
  {"ViewLightcolor.z",&view_lightcolor.z,Type_Float},
  {"ViewAmbient",&view_ambient,Type_Float},
  {"ViewDiffuse",&view_diffuse,Type_Float},
  {"ViewSpecular",&view_specular,Type_Float},
  {"ViewShininess",&view_shininess,Type_Int},
  {"ViewShowAtoms",&view_show_atoms,Type_Bool},
  {"ViewScaleAtoms",&view_scale_atoms,Type_Float},
  {"ViewIsphres",&view_isphres,Type_Int},
  {"ViewShowBonds",&view_show_bonds,Type_Bool},
  {"ViewScaleBonds",&view_scale_bonds,Type_Float},
  {"ViewIcylres",&view_icylres,Type_Int},
  {"ViewMousesensRot",&view_mousesens_rot,Type_Float},
  {"ViewMousesensZoom",&view_mousesens_zoom,Type_Float},

  // Widget colors (imgui_widgets.h)
  {"Color_Slidingbar",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_Slidingbar],Type_ImVec4},
  {"Color_SlidingbarHovered",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarHovered],Type_ImVec4},
  {"Color_SlidingbarActive",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarActive],Type_ImVec4},
  {"Color_Tab",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_Tab],Type_ImVec4},
  {"Color_TabHovered",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabHovered],Type_ImVec4},
  {"Color_TabPressed",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabPressed],Type_ImVec4},
  {"Color_TabActive",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabActive],Type_ImVec4},
  {"Color_TabXFg",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFg],Type_ImVec4},
  {"Color_TabXFgHovered",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgHovered],Type_ImVec4},
  {"Color_TabXFgActive",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgActive],Type_ImVec4},
  {"Color_TabXBg",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBg],Type_ImVec4},
  {"Color_TabXBgHovered",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgHovered],Type_ImVec4},
  {"Color_TabXBgActive",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgActive],Type_ImVec4},
  {"Color_TabBorder",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabBorder],Type_ImVec4},
  {"Color_DropTarget",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTarget],Type_ImVec4},
  {"Color_DropTargetActive",&ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTargetActive],Type_ImVec4},

  // Widget styles (imgui_widgets.h)
  {"TabRounding",&ImGuiStyleWidgets.TabRounding,Type_Float},
  {"TabBorderSize",&ImGuiStyleWidgets.TabBorderSize,Type_Float},
  {"DropTargetLooseness",&ImGuiStyleWidgets.DropTargetLooseness,Type_Float},
  {"DropTargetMinsizeEdge",&ImGuiStyleWidgets.DropTargetMinsizeEdge,Type_Float},
  {"DropTargetMaxsizeEdge",&ImGuiStyleWidgets.DropTargetMaxsizeEdge,Type_Float},
  {"DropTargetEdgeFraction",&ImGuiStyleWidgets.DropTargetEdgeFraction,Type_Float},
  {"DropTargetFullFraction",&ImGuiStyleWidgets.DropTargetFullFraction,Type_Float},
  {"TabHeight",&ImGuiStyleWidgets.TabHeight,Type_Float},
  {"TabMaxWidth",&ImGuiStyleWidgets.TabMaxWidth,Type_Float},
  {"SlidingBarWidth",&ImGuiStyleWidgets.SlidingBarWidth,Type_Float},

  // imgui colors (imgui.h)
  {"Color_Text",&GImGui->Style.Colors[ImGuiCol_Text],Type_ImVec4},
  {"Color_TextDisabled",&GImGui->Style.Colors[ImGuiCol_TextDisabled],Type_ImVec4},
  {"Color_WindowBg",&GImGui->Style.Colors[ImGuiCol_WindowBg],Type_ImVec4},
  {"Color_ChildBg",&GImGui->Style.Colors[ImGuiCol_ChildBg],Type_ImVec4},
  {"Color_PopupBg",&GImGui->Style.Colors[ImGuiCol_PopupBg],Type_ImVec4},
  {"Color_Border",&GImGui->Style.Colors[ImGuiCol_Border],Type_ImVec4},
  {"Color_BorderShadow",&GImGui->Style.Colors[ImGuiCol_BorderShadow],Type_ImVec4},
  {"Color_FrameBg",&GImGui->Style.Colors[ImGuiCol_FrameBg],Type_ImVec4},
  {"Color_FrameBgHovered",&GImGui->Style.Colors[ImGuiCol_FrameBgHovered],Type_ImVec4},
  {"Color_FrameBgActive",&GImGui->Style.Colors[ImGuiCol_FrameBgActive],Type_ImVec4},
  {"Color_TitleBg",&GImGui->Style.Colors[ImGuiCol_TitleBg],Type_ImVec4},
  {"Color_TitleBgActive",&GImGui->Style.Colors[ImGuiCol_TitleBgActive],Type_ImVec4},
  {"Color_TitleBgCollapsed",&GImGui->Style.Colors[ImGuiCol_TitleBgCollapsed],Type_ImVec4},
  {"Color_MenuBarBg",&GImGui->Style.Colors[ImGuiCol_MenuBarBg],Type_ImVec4},
  {"Color_ScrollbarBg",&GImGui->Style.Colors[ImGuiCol_ScrollbarBg],Type_ImVec4},
  {"Color_ScrollbarGrab",&GImGui->Style.Colors[ImGuiCol_ScrollbarGrab],Type_ImVec4},
  {"Color_ScrollbarGrabHovered",&GImGui->Style.Colors[ImGuiCol_ScrollbarGrabHovered],Type_ImVec4},
  {"Color_ScrollbarGrabActive",&GImGui->Style.Colors[ImGuiCol_ScrollbarGrabActive],Type_ImVec4},
  {"Color_CheckMark",&GImGui->Style.Colors[ImGuiCol_CheckMark],Type_ImVec4},
  {"Color_SliderGrab",&GImGui->Style.Colors[ImGuiCol_SliderGrab],Type_ImVec4},
  {"Color_SliderGrabActive",&GImGui->Style.Colors[ImGuiCol_SliderGrabActive],Type_ImVec4},
  {"Color_Button",&GImGui->Style.Colors[ImGuiCol_Button],Type_ImVec4},
  {"Color_ButtonHovered",&GImGui->Style.Colors[ImGuiCol_ButtonHovered],Type_ImVec4},
  {"Color_ButtonActive",&GImGui->Style.Colors[ImGuiCol_ButtonActive],Type_ImVec4},
  {"Color_Header",&GImGui->Style.Colors[ImGuiCol_Header],Type_ImVec4},
  {"Color_HeaderHovered",&GImGui->Style.Colors[ImGuiCol_HeaderHovered],Type_ImVec4},
  {"Color_HeaderActive",&GImGui->Style.Colors[ImGuiCol_HeaderActive],Type_ImVec4},
  {"Color_Separator",&GImGui->Style.Colors[ImGuiCol_Separator],Type_ImVec4},
  {"Color_SeparatorHovered",&GImGui->Style.Colors[ImGuiCol_SeparatorHovered],Type_ImVec4},
  {"Color_SeparatorActive",&GImGui->Style.Colors[ImGuiCol_SeparatorActive],Type_ImVec4},
  {"Color_ResizeGrip",&GImGui->Style.Colors[ImGuiCol_ResizeGrip],Type_ImVec4},
  {"Color_ResizeGripHovered",&GImGui->Style.Colors[ImGuiCol_ResizeGripHovered],Type_ImVec4},
  {"Color_ResizeGripActive",&GImGui->Style.Colors[ImGuiCol_ResizeGripActive],Type_ImVec4},
  {"Color_CloseButton",&GImGui->Style.Colors[ImGuiCol_CloseButton],Type_ImVec4},
  {"Color_CloseButtonHovered",&GImGui->Style.Colors[ImGuiCol_CloseButtonHovered],Type_ImVec4},
  {"Color_CloseButtonActive",&GImGui->Style.Colors[ImGuiCol_CloseButtonActive],Type_ImVec4},
  {"Color_PlotLines",&GImGui->Style.Colors[ImGuiCol_PlotLines],Type_ImVec4},
  {"Color_PlotLinesHovered",&GImGui->Style.Colors[ImGuiCol_PlotLinesHovered],Type_ImVec4},
  {"Color_PlotHistogram",&GImGui->Style.Colors[ImGuiCol_PlotHistogram],Type_ImVec4},
  {"Color_PlotHistogramHovered",&GImGui->Style.Colors[ImGuiCol_PlotHistogramHovered],Type_ImVec4},
  {"Color_TextSelectedBg",&GImGui->Style.Colors[ImGuiCol_TextSelectedBg],Type_ImVec4},
  {"Color_ModalWindowDarkening",&GImGui->Style.Colors[ImGuiCol_ModalWindowDarkening],Type_ImVec4},

  // imgui styles
  {"Alpha",&GImGui->Style.Alpha,Type_Float},
  {"WindowPadding",&GImGui->Style.WindowPadding,Type_ImVec2},
  {"WindowRounding",&GImGui->Style.WindowRounding,Type_Float},
  {"WindowBorderSize",&GImGui->Style.WindowBorderSize,Type_Float},
  {"WindowMinSize",&GImGui->Style.WindowMinSize,Type_ImVec2},
  {"WindowTitleAlign",&GImGui->Style.WindowTitleAlign,Type_ImVec2},
  {"ChildRounding",&GImGui->Style.ChildRounding,Type_Float},
  {"ChildBorderSize",&GImGui->Style.ChildBorderSize,Type_Float},
  {"PopupRounding",&GImGui->Style.PopupRounding,Type_Float},
  {"PopupBorderSize",&GImGui->Style.PopupBorderSize,Type_Float},
  {"FramePadding",&GImGui->Style.FramePadding,Type_ImVec2},
  {"FrameRounding",&GImGui->Style.FrameRounding,Type_Float},
  {"FrameBorderSize",&GImGui->Style.FrameBorderSize,Type_Float},
  {"ItemSpacing",&GImGui->Style.ItemSpacing,Type_ImVec2},
  {"ItemInnerSpacing",&GImGui->Style.ItemInnerSpacing,Type_ImVec2},
  {"TouchExtraPadding",&GImGui->Style.TouchExtraPadding,Type_ImVec2},
  {"IndentSpacing",&GImGui->Style.IndentSpacing,Type_Float},
  {"ColumnsMinSpacing",&GImGui->Style.ColumnsMinSpacing,Type_Float},
  {"ScrollbarSize",&GImGui->Style.ScrollbarSize,Type_Float},
  {"ScrollbarRounding",&GImGui->Style.ScrollbarRounding,Type_Float},
  {"GrabMinSize",&GImGui->Style.GrabMinSize,Type_Float},
  {"GrabRounding",&GImGui->Style.GrabRounding,Type_Float},
  {"ButtonTextAlign",&GImGui->Style.ButtonTextAlign,Type_ImVec2},
  {"DisplayWindowPadding",&GImGui->Style.DisplayWindowPadding,Type_ImVec2},
  {"DisplaySafeAreaPadding",&GImGui->Style.DisplaySafeAreaPadding,Type_ImVec2},
  {"AntiAliasedLines",&GImGui->Style.AntiAliasedLines,Type_Bool},
  {"AntiAliasedShapes",&GImGui->Style.AntiAliasedShapes,Type_Bool},
  {"CurveTessellationTol",&GImGui->Style.CurveTessellationTol,Type_Float},

  // keybinding
  {"Bind_QUIT_mod",&modbind[BIND_QUIT],Type_Int},
  {"Bind_QUIT_key",&keybind[BIND_QUIT],Type_Int},
  {"Bind_CLOSE_LAST_DIALOG_mod",&modbind[BIND_CLOSE_LAST_DIALOG],Type_Int},
  {"Bind_CLOSE_LAST_DIALOG_key",&keybind[BIND_CLOSE_LAST_DIALOG],Type_Int},
  {"Bind_CLOSE_ALL_DIALOGS_mod",&modbind[BIND_CLOSE_ALL_DIALOGS],Type_Int},
  {"Bind_CLOSE_ALL_DIALOGS_key",&keybind[BIND_CLOSE_ALL_DIALOGS],Type_Int},
  {"Bind_VIEW_ALIGN_A_AXIS_mod",&modbind[BIND_VIEW_ALIGN_A_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_A_AXIS_key",&keybind[BIND_VIEW_ALIGN_A_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_B_AXIS_mod",&modbind[BIND_VIEW_ALIGN_B_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_B_AXIS_key",&keybind[BIND_VIEW_ALIGN_B_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_C_AXIS_mod",&modbind[BIND_VIEW_ALIGN_C_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_C_AXIS_key",&keybind[BIND_VIEW_ALIGN_C_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_X_AXIS_mod",&modbind[BIND_VIEW_ALIGN_X_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_X_AXIS_key",&keybind[BIND_VIEW_ALIGN_X_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_Y_AXIS_mod",&modbind[BIND_VIEW_ALIGN_Y_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_Y_AXIS_key",&keybind[BIND_VIEW_ALIGN_Y_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_Z_AXIS_mod",&modbind[BIND_VIEW_ALIGN_Z_AXIS],Type_Int},
  {"Bind_VIEW_ALIGN_Z_AXIS_key",&keybind[BIND_VIEW_ALIGN_Z_AXIS],Type_Int},
  {"Bind_NAV_ROTATE_mod",&modbind[BIND_NAV_ROTATE],Type_Int},
  {"Bind_NAV_ROTATE_key",&keybind[BIND_NAV_ROTATE],Type_Int},
  {"Bind_NAV_TRANSLATE_mod",&modbind[BIND_NAV_TRANSLATE],Type_Int},
  {"Bind_NAV_TRANSLATE_key",&keybind[BIND_NAV_TRANSLATE],Type_Int},
  {"Bind_NAV_ZOOM_mod",&modbind[BIND_NAV_ZOOM],Type_Int},
  {"Bind_NAV_ZOOM_key",&keybind[BIND_NAV_ZOOM],Type_Int},
  {"Bind_NAV_RESET_mod",&modbind[BIND_NAV_RESET],Type_Int},
  {"Bind_NAV_RESET_key",&keybind[BIND_NAV_RESET],Type_Int},

  // sentinel
  {nullptr,nullptr,Type_None},
};

void DefaultSettings(){
  static bool firstpass = true;

  // default colors
  UIStyleColorsMutantOrange();

  // default ImGui settings (imgui.h)
  ImGuiStyle& style = GetStyle();
  style.Alpha                   = 1.0f;
  style.WindowPadding           = ImVec2(8,8);
  style.WindowRounding          = 7.0f;
  style.WindowBorderSize        = 1.0f;
  style.WindowMinSize           = ImVec2(32,32);
  style.WindowTitleAlign        = ImVec2(0.0f,0.5f);
  style.ChildRounding           = 7.0f;
  style.ChildBorderSize         = 1.0f;
  style.PopupRounding           = 7.0f;
  style.PopupBorderSize         = 1.0f;
  style.FramePadding            = ImVec2(4,3);
  style.FrameRounding           = 7.0f;
  style.FrameBorderSize         = 0.0f;
  style.ItemSpacing             = ImVec2(8,4);
  style.ItemInnerSpacing        = ImVec2(4,4);
  style.TouchExtraPadding       = ImVec2(0,0);
  style.IndentSpacing           = 21.0f;
  style.ColumnsMinSpacing       = 6.0f;
  style.ScrollbarSize           = 12.0f;
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
  ImGuiStyleWidgets.DefaultStyle();

  // default UI settings (settings.h)
  ImGuiStyleUI.DefaultStyle();

  // default view settings (settings.h)
  view_wireframe = false;
  view_orthogonal = false;
  view_fov = 45.0f;
  view_resetdistance = 1.1f;
  view_bgrgb[0] = view_bgrgb[1] = view_bgrgb[2] = 0.f;
  view_bgrgb[3] = 1.0f;
  view_lightpos[0] = view_lightpos[1] = 20.f;
  view_lightpos[2] = 0.f;
  view_lightcolor[0] = view_lightcolor[1] = view_lightcolor[2] = 1.f;
  view_ambient = 0.2f;
  view_diffuse = 0.4f;
  view_specular = 0.6f;
  view_shininess = 8;
  view_show_atoms = true;
  view_scale_atoms = 1.0f;
  view_isphres = 2;
  view_show_bonds = true;
  view_scale_bonds = 1.0f;
  view_icylres = 1;
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

  // set the defualt fonts and sizes
  SetUIFont(ImGuiStyleUI.FontSelected,ImGuiStyleUI.FontSize,ImGuiStyleUI.FontSizeIcon);

  firstpass = false;
}

void ScaleUI(float scale){
  ImGuiIO& io = GetIO();

  ImGuiStyleUI.MessageWidth *= scale;
  ImGuiStyleUI.FontSizeIcon *= scale;
  ImGuiStyleUI.FontSize *= scale;
  SetUIFont(-1,ImGuiStyleUI.FontSize,ImGuiStyleUI.FontSizeIcon);
  ImGuiStyleUI.TooltipMaxwidth *= scale;
  ImGuiStyleWidgets.TabRounding *= scale;
  ImGuiStyleWidgets.DropTargetLooseness *= scale;
  ImGuiStyleWidgets.DropTargetMinsizeEdge *= scale;
  ImGuiStyleWidgets.DropTargetMaxsizeEdge *= scale;
  ImGuiStyleWidgets.TabHeight *= scale;
  ImGuiStyleWidgets.TabMaxWidth *= scale;
  ImGuiStyleWidgets.SlidingBarWidth *= scale;
  GImGui->Style.WindowPadding.x *= scale;
  GImGui->Style.WindowPadding.y *= scale;
  GImGui->Style.WindowRounding *= scale;
  GImGui->Style.WindowMinSize.x *= scale;
  GImGui->Style.WindowMinSize.y *= scale;
  GImGui->Style.ChildRounding *= scale;
  GImGui->Style.PopupRounding *= scale;
  GImGui->Style.FramePadding.x *= scale;
  GImGui->Style.FramePadding.y *= scale;
  GImGui->Style.FrameRounding *= scale;
  GImGui->Style.ItemSpacing.x *= scale;
  GImGui->Style.ItemSpacing.y *= scale;
  GImGui->Style.ItemInnerSpacing.x *= scale;
  GImGui->Style.ItemInnerSpacing.y *= scale;
  GImGui->Style.TouchExtraPadding.x *= scale;
  GImGui->Style.TouchExtraPadding.y *= scale;
  GImGui->Style.IndentSpacing *= scale;
  GImGui->Style.ColumnsMinSpacing *= scale;
  GImGui->Style.ScrollbarSize *= scale;
  GImGui->Style.ScrollbarRounding *= scale;
  GImGui->Style.GrabMinSize *= scale;
  GImGui->Style.GrabRounding *= scale;
  GImGui->Style.DisplayWindowPadding.x *= scale;
  GImGui->Style.DisplayWindowPadding.y *= scale;
  GImGui->Style.DisplaySafeAreaPadding.x *= scale;
  GImGui->Style.DisplaySafeAreaPadding.y *= scale;
}

void SetUIFont(int ifont/*=-1*/,float size/*=-1.f*/,float sizeicon/*=-1.f*/){
  ImGuiIO& io = GetIO();

  if (ifont != -1){
    ImGuiStyleUI.FontSelected = ifont;
    fontdefault = io.Fonts->Fonts[ImGuiStyleUI.FontSelected];
    io.FontDefault = fontdefault;
  }
  if (size > 0.f){
    ImGuiStyleUI.FontSize = size;
    for (int i = 1; i < io.Fonts->Fonts.Size; i++)
      io.Fonts->Fonts[i]->Scale = ImGuiStyleUI.FontSize / fontsizebake;
  }
  if (sizeicon > 0.f){
    ImGuiStyleUI.FontSizeIcon = sizeicon;
    io.Fonts->Fonts[0]->Scale = ImGuiStyleUI.FontSizeIcon / fontsizebake;
  }
}

// classic colors style
void UIStyleColorsMutantOrange(){

  ImGuiStyle& style = GetStyle();
  style.Colors[ImGuiCol_Text]                   = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  style.Colors[ImGuiCol_TextDisabled]           = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  style.Colors[ImGuiCol_WindowBg]               = ImVec4(0.00f, 0.00f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_ChildBg]                = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
  style.Colors[ImGuiCol_PopupBg]                = ImVec4(0.14f, 0.14f, 0.14f, 0.90f);
  style.Colors[ImGuiCol_Border]                 = ImVec4(0.43f, 0.43f, 0.43f, 0.60f);
  style.Colors[ImGuiCol_BorderShadow]           = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
  style.Colors[ImGuiCol_FrameBg]                = ImVec4(0.43f, 0.43f, 0.43f, 0.40f);
  style.Colors[ImGuiCol_FrameBgHovered]         = ImVec4(0.60f, 0.60f, 0.60f, 0.40f);
  style.Colors[ImGuiCol_FrameBgActive]          = ImVec4(0.80f, 0.80f, 0.80f, 0.40f);
  style.Colors[ImGuiCol_TitleBg]                = ImVec4(0.40f, 0.20f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_TitleBgActive]          = ImVec4(0.75f, 0.50f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_TitleBgCollapsed]       = ImVec4(0.95f, 0.70f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_MenuBarBg]              = ImVec4(0.75f, 0.50f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_ScrollbarBg]            = ImVec4(0.43f, 0.43f, 0.43f, 0.40f);
  style.Colors[ImGuiCol_ScrollbarGrab]          = ImVec4(0.40f, 0.20f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_ScrollbarGrabHovered]   = ImVec4(0.75f, 0.50f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_ScrollbarGrabActive]    = ImVec4(0.95f, 0.70f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_CheckMark]              = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  style.Colors[ImGuiCol_SliderGrab]             = ImVec4(0.60f, 0.60f, 0.60f, 0.70f);
  style.Colors[ImGuiCol_SliderGrabActive]       = ImVec4(0.40f, 0.70f, 0.20f, 0.70f);
  style.Colors[ImGuiCol_Button]                 = ImVec4(0.40f, 0.20f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_ButtonHovered]          = ImVec4(0.75f, 0.50f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_ButtonActive]           = ImVec4(0.40f, 0.70f, 0.20f, 0.90f);
  style.Colors[ImGuiCol_Header]                 = ImVec4(0.40f, 0.70f, 0.20f, 0.40f);
  style.Colors[ImGuiCol_HeaderHovered]          = ImVec4(0.60f, 0.80f, 0.30f, 0.90f);
  style.Colors[ImGuiCol_HeaderActive]           = ImVec4(0.70f, 0.95f, 0.40f, 0.90f);
  style.Colors[ImGuiCol_Separator]              = ImVec4(0.43f, 0.43f, 0.43f, 1.00f);
  style.Colors[ImGuiCol_SeparatorHovered]       = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  style.Colors[ImGuiCol_SeparatorActive]        = ImVec4(0.80f, 0.80f, 0.80f, 1.00f);
  style.Colors[ImGuiCol_ResizeGrip]             = ImVec4(0.43f, 0.43f, 0.43f, 0.40f);
  style.Colors[ImGuiCol_ResizeGripHovered]      = ImVec4(0.60f, 0.60f, 0.60f, 0.40f);
  style.Colors[ImGuiCol_ResizeGripActive]       = ImVec4(0.80f, 0.80f, 0.80f, 0.40f);
  style.Colors[ImGuiCol_CloseButton]            = ImVec4(0.43f, 0.43f, 0.43f, 0.50f);
  style.Colors[ImGuiCol_CloseButtonHovered]     = ImVec4(0.80f, 0.20f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_CloseButtonActive]      = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_PlotLines]              = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
  style.Colors[ImGuiCol_PlotLinesHovered]       = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogram]          = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_PlotHistogramHovered]   = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
  style.Colors[ImGuiCol_TextSelectedBg]         = ImVec4(0.40f, 0.70f, 0.20f, 0.35f);
  style.Colors[ImGuiCol_ModalWindowDarkening]   = ImVec4(0.00f, 0.00f, 0.00f, 0.35f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Slidingbar]        = ImVec4(0.43f, 0.43f, 0.43f, 0.40f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarHovered] = ImVec4(0.40f, 0.20f, 0.00f, 0.90f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_SlidingbarActive]  = ImVec4(0.75f, 0.50f, 0.00f, 0.90f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_Tab]               = ImVec4(0.43f, 0.43f, 0.43f, 0.40f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabHovered]        = ImVec4(0.75f, 0.50f, 0.00f, 0.90f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabPressed]        = ImVec4(0.95f, 0.70f, 0.00f, 0.90f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabActive]         = ImVec4(0.40f, 0.70f, 0.20f, 0.40f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFg]            = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgHovered]     = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXFgActive]      = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBg]            = ImVec4(0.80f, 0.20f, 0.00f, 0.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgHovered]     = ImVec4(0.80f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabXBgActive]      = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_TabBorder]         = ImVec4(0.43f, 0.43f, 0.43f, 0.60f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGrip]          = ImVec4(0.60f, 0.20f, 0.00f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripHovered]   = ImVec4(0.80f, 0.40f, 0.20f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_LiftGripActive]    = ImVec4(1.00f, 0.40f, 0.20f, 1.00f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTarget]        = ImVec4(0.43f, 0.43f, 0.43f, 0.40f);
  ImGuiStyleWidgets.Colors[ImGuiColWidgets_DropTargetActive]  = ImVec4(0.80f, 0.80f, 0.80f, 0.80f);
  ImGuiStyleUI.Colors[ImGuiColUI_BackDrop]         = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIcon]         = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconHovered]  = ImVec4(0.95f, 0.70f, 0.00f, 0.90f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconActive]   = ImVec4(0.75f, 0.50f, 0.00f, 0.90f);
  ImGuiStyleUI.Colors[ImGuiColUI_ViewIconInactive] = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageInfo]      = ImVec4(0.0f,1.0f,0.0f,0.4f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageWarning]   = ImVec4(1.0f,1.0f,0.0f,0.4f);
  ImGuiStyleUI.Colors[ImGuiColUI_MessageError]     = ImVec4(1.0f,0.0f,0.0f,0.4f);
}

// classic colors style
void UIStyleColorsClassic(){
  ImGuiStyle& style = GetStyle();
  style.Colors[ImGuiCol_Text]                   = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  style.Colors[ImGuiCol_TextDisabled]           = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
  style.Colors[ImGuiCol_WindowBg]               = ImVec4(0.00f, 0.00f, 0.00f, 0.90f);
  style.Colors[ImGuiCol_ChildBg]                = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
  style.Colors[ImGuiCol_PopupBg]                = ImVec4(0.11f, 0.11f, 0.14f, 0.90f);
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
  ImGuiStyleUI.Colors[ImGuiColUI_BackDrop]         = ImVec4(0.20f, 0.30f, 0.30f, 1.00f);
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
  ImGuiStyleUI.Colors[ImGuiColUI_BackDrop]         = ImVec4(0.20f, 0.30f, 0.30f, 1.00f);
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
  ImGuiStyleUI.Colors[ImGuiColUI_BackDrop]         = ImVec4(0.20f, 0.30f, 0.30f, 1.00f);
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
namespace ImGui{
  static void to_json(json& j, const Dock& p) {
    j["label"] = p.label;
    int n = 0;
    for (auto dd: p.stack){
      ++n;
      j[to_string(n)] = *dd;
    }
  }
}

void WriteLayout(Dock* root){
  json j;
  
  j["root"] = (*root);
  cout << setw(4) << j << endl;
}

bool WriteConfigurationFile(string file){
  ofstream fo(file);
  if (!fo.good())
    return false;

  json j;

  for (int i = 0;;i++){
    if (!confvar[i].name)
      break;
    switch(confvar[i].vartype){
    case Type_Bool:
      j[confvar[i].name] = *((bool *) confvar[i].variable); break;
    case Type_Int:
      j[confvar[i].name] = *((int *) confvar[i].variable); break;
    case Type_Float:
      j[confvar[i].name] = *((float *) confvar[i].variable); break;
    case Type_ImVec2:
      j[confvar[i].name] = *((ImVec2 *) confvar[i].variable); break;
    case Type_ImVec4:
      j[confvar[i].name] = *((ImVec4 *) confvar[i].variable); break;
    }
  }

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
  try {
    fi >> j;
  } catch(json::parse_error) {
    return false;
  }

  for (int i = 0;;i++){
    if (!confvar[i].name)
      break;
    switch(confvar[i].vartype){
    case Type_Bool:
      if (j.count(confvar[i].name) != 0) *((bool *) confvar[i].variable) = j[confvar[i].name].get<bool>(); break;
    case Type_Int:
      if (j.count(confvar[i].name) != 0) *((int *) confvar[i].variable) = j[confvar[i].name].get<int>(); break;
    case Type_Float:
      if (j.count(confvar[i].name) != 0) *((float *) confvar[i].variable) = j[confvar[i].name].get<float>(); break;
    case Type_ImVec2:
      if (j.count(confvar[i].name) != 0) *((ImVec2 *) confvar[i].variable) = j[confvar[i].name].get<ImVec2>(); break;
    case Type_ImVec4:
      if (j.count(confvar[i].name) != 0) *((ImVec4 *) confvar[i].variable) = j[confvar[i].name].get<ImVec4>(); break;
    }
  }

  fi.close();

  // Set the default font just read
  SetUIFont(ImGuiStyleUI.FontSelected,ImGuiStyleUI.FontSize,ImGuiStyleUI.FontSizeIcon);

  return true;
}

// Some global callbacks //

void error_callback(int error, const char* description){
  fprintf(stderr, "Error %d: %s\n", error, description);
}
