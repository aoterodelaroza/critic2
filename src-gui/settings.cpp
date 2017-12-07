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

using namespace ImGui;

// Global variables: fonts (see settings.h)
ImFont* fontdefault = nullptr;
ImFont* fonticon = nullptr;
const float fontsizebake = 36.0f;
const float fontsizeicon = 24.0f;
float fontsize = 14.0f;

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

// Tooltips
bool tooltip_enabled; // enabled tooltips
float tooltip_delay; // delay in seconds
float tooltip_maxwidth; // maxwidth of the tooltip

ImGuiStyleUI_::ImGuiStyleUI_(){
  Colors[ImGuiColUI_ViewIcon]         = ImVec4(0.90f, 0.90f, 0.90f, 1.00f);
  Colors[ImGuiColUI_ViewIconHovered]  = ImVec4(0.8549f,0.6471f,0.1255f,1.0f);
  Colors[ImGuiColUI_ViewIconActive]   = ImVec4(0.7216f,0.5254,0.04314f,1.0f);
  Colors[ImGuiColUI_ViewIconInactive] = ImVec4(0.5f,0.5f,0.5f,1.0f);
}

void DefaultSettings(){
  // views
  view_wireframe = false;
  view_orthogonal = false;
  view_fov = 45.0f;
  view_resetdistance = 1.25f; 
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

  // tooltips
  tooltip_enabled = true; // enabled tooltips
  tooltip_delay = 1.5f; // delay in seconds
  tooltip_maxwidth = 450.f; // maxwidth of the tooltip

  // set up shader and fill uniforms
  shader = new Shader();
  shader->use();
  shader->setVec3("lightPos",value_ptr(view_lightpos));
  shader->setVec3("lightColor",value_ptr(view_lightcolor));
  shader->setFloat("ambient",view_ambient);
  shader->setFloat("diffuse",view_diffuse);
  shader->setFloat("specular",view_specular);
  shader->setInt("shininess",view_shininess);

  // load standard fonts
  // const ImWchar* ImFontAtlas::GetGlyphRangesDefault()
  ImGuiIO& io = GetIO();
  ImFontConfig fntconfig; 
  fntconfig.MergeMode = false; 
  fntconfig.PixelSnapH = true;
  fntconfig.OversampleH = fntconfig.OversampleV = 4;
  fntconfig.MergeMode = false;
  fntconfig.SizePixels = fontsizebake;

  static const ImWchar ranges[] = { 0x0020, 0x00FF, 0 }; // basic latin + latin supplement
  static const ImWchar rangesicon[] = { 0x0021, 0x002c, 0 }; // icons

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

  // set the defualt sizes
  io.Fonts->Fonts[0]->Scale = fontsizeicon / fontsizebake;
  for (int i = 1; i < io.Fonts->Fonts.Size; i++)
    io.Fonts->Fonts[i]->Scale = fontsize / fontsizebake;
    
  // set the default font
  fontdefault = io.Fonts->Fonts[1];
  io.FontDefault = fontdefault;
}

// Some global callbacks //

void error_callback(int error, const char* description){
  fprintf(stderr, "Error %d: %s\n", error, description);
}

