// -*-c++-*-

// global settings for the GUI (maybe to be changeable later on).

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdio.h>
#include <stdlib.h>

#include <glm/glm.hpp>
#include "imgui/imgui.h"
#include "shader.h"

// constants
#define PI 3.14159265358979323846
#define AUTOANG 0.52917720859

// Shader
extern Shader *shader;

// Fonts (defined in main.cpp)
extern ImFont* fontdefault;
extern ImFont* fonticon;
extern const float fontsizebake;
extern const float fontsizeicon;
extern float fontsize;

// Color picker options
const int coloreditflags = ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_AlphaBar | 
  ImGuiColorEditFlags_AlphaPreview | ImGuiColorEditFlags_PickerHueBar | 
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
extern int view_isphres;
extern int view_icylres;
const float mousesens_rot0 = 2.0f;
extern float view_mousesens_rot;
const float mousesens_zoom0 = 0.15f;
extern float view_mousesens_zoom;

// Tooltips
extern bool tooltip_enabled;
extern float tooltip_delay;
extern float tooltip_maxwidth;

// Default settings
void DefaultSettings();

// Callbacks
void error_callback(int error, const char* description);

#endif SETTINGS_H

