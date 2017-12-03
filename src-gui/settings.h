// -*-c++-*-

// global settings for the GUI (maybe to be changeable later on).

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdio.h>
#include <stdlib.h>

#include <glm/glm.hpp>
#include "imgui/imgui.h"
#include "imgui/mouse.h"

// constants
#define PI 3.14159265358979323846
#define AUTOANG 0.52917720859

// Fonts (defined in main.cpp)
extern ImFont* fontdefault;
extern ImFont* fonticon;
const float fonticon_size = 24.0f;

// Mouse state
extern MouseState mstate;

// Framebuffer texture default side length
 const float FBO_tex_a = 1024.f;

// Background color for views
const float bgrgb[4] = {0.0f,0.0f,0.0f,1.0f};

// Light parameters
const glm::vec3 lightPos = {20.f,20.f,0.f};
const glm::vec3 lightColor = {1.f,1.f,1.f};
const float ambient = 0.2f;
const float diffuse = 0.4f;
const float specular = 0.6f;
const int shininess = 8;

// Atom & bond resolution: 0 -> nmaxsph-1,nmaxcyl-1
const int isphres = 2;
const int icylres = 0;

// Grid and Cartesian axes
const float radgrid = 0.02;

// Mouse constants
const float mousesens_rot = 2.0f; // Mouse rotate sensitivity
const float mousesens_zoom = 0.15f; // Mouse zoom sensitivity
const float zfov = 45.f; // fov for the perspective
const float znear = 0.1f; // znear for the camera
const float zfar = 1000.f; // zfar for the camera
const float min_zoom = 1.f; // minimum distance to origin (zoom)

// Tooltip delay
const float tooltip_delay = 1.5f; // in seconds
const float tooltip_maxwidth = 450.f; // in pixels

#endif SETTINGS_H

