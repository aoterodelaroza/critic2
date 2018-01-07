// -*-c++-*-

// global settings for the GUI (maybe to be changeable later on).

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <glm/glm.hpp>

#include "imgui/imgui.h"
#include "imgui/imgui_dock.h"

#include "shader.h"
#include "view.h"

// constants
#define PI 3.14159265358979323846
#define AUTOANG 0.52917720859

// UI style
enum ImGuiColUI_ {
  ImGuiColUI_BackDrop,
  ImGuiColUI_ViewIcon,
  ImGuiColUI_ViewIconHovered,
  ImGuiColUI_ViewIconActive,
  ImGuiColUI_ViewIconInactive,
  ImGuiColUI_MessageInfo,
  ImGuiColUI_MessageWarning,
  ImGuiColUI_MessageError,
  ImGuiColUI_COUNT,
};
struct ImGuiStyleUI_ {
  ImVec4 Colors[ImGuiColUI_COUNT];

  float MessageWidth;
  float MessageExpire;
  float FontSizeIcon;
  float FontSize;
  int FontSelected;
  bool TooltipEnabled;
  float TooltipDelay;
  float TooltipMaxwidth;
  float UIScaleFactor;

  void DefaultStyle(){
    MessageWidth = 400.f;
    MessageExpire = 5.f;
    FontSizeIcon = 24.0f;
    FontSize = 14.0f;
    FontSelected = 1;
    TooltipEnabled = true;
    TooltipDelay = 1.5f;
    TooltipMaxwidth = 450.f;
    UIScaleFactor = 1.0f;
  }
  void DefaultColors(){
    Colors[ImGuiColUI_BackDrop]         = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
    Colors[ImGuiColUI_ViewIcon]         = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
    Colors[ImGuiColUI_ViewIconHovered]  = ImVec4(0.95f, 0.70f, 0.00f, 0.90f);
    Colors[ImGuiColUI_ViewIconActive]   = ImVec4(0.75f, 0.50f, 0.00f, 0.90f);
    Colors[ImGuiColUI_ViewIconInactive] = ImVec4(0.60f, 0.60f, 0.60f, 1.00f);
    Colors[ImGuiColUI_MessageInfo]      = ImVec4(0.0f,1.0f,0.0f,0.4f);
    Colors[ImGuiColUI_MessageWarning]   = ImVec4(1.0f,1.0f,0.0f,0.4f);
    Colors[ImGuiColUI_MessageError]     = ImVec4(1.0f,0.0f,0.0f,0.4f);
  }
  ImGuiStyleUI_(){
    DefaultStyle();
    DefaultColors();
  };
};
extern ImGuiStyleUI_ ImGuiStyleUI;

// Shader
extern Shader *shader;

// Fonts (defined in main.cpp)
extern ImFont* fontdefault;
extern ImFont* fonticon;
extern const float fontsizebake;

// Color picker options
const int coloreditflags = ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_AlphaBar | 
  ImGuiColorEditFlags_AlphaPreviewHalf | ImGuiColorEditFlags_PickerHueBar | 
  ImGuiColorEditFlags_Uint8 | ImGuiColorEditFlags_RGB;

// Framebuffer texture default side length
const float FBO_tex_a = 1024.f;

// Mouse constants
const float znear = 0.1f; // znear for the camera
const float zfar = 1000.f; // zfar for the camera
const float min_zoom = 1.f; // minimum distance to origin (zoom, bohr)
const float max_zoom = 100.f; // maximum distance to origin (zoom, scenerad)

// View defaults and settings
extern bool view_wireframe;
extern bool view_orthogonal;
extern float view_fov;
extern float view_resetdistance; 
extern float view_bgrgb[4];
extern glm::vec3 view_lightpos;
extern glm::vec3 view_lightcolor;
extern float view_ambient;
extern float view_diffuse;
extern float view_specular;
extern int view_shininess;
extern bool view_show_atoms;
extern float view_scale_atoms;
extern int view_isphres;
extern bool view_show_bonds;
extern float view_scale_bonds;
extern int view_icylres;
const float mousesens_rot0 = 2.0f;
extern float view_mousesens_rot;
const float mousesens_zoom0 = 0.15f;
extern float view_mousesens_zoom;

// Configuration file
extern std::string conffile;

// Default settings
void DefaultSettings();
void ScaleUI(float scale);
void SetUIFont(int ifont=-1,float size=-1.f,float sizeicon=-1.f);
void UIStyleColorsMutantOrange();
void UIStyleColorsClassic();
void UIStyleColorsDark();
void UIStyleColorsLight();

// Configuration file
bool FindConfigurationFile();
void WriteLayout(ImGui::Dock* root);
bool WriteConfigurationFile(std::string file);
bool ReadConfigurationFile(std::string file);

// Callbacks
void error_callback(int error, const char* description);

#endif SETTINGS_H

