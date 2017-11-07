// -*-c++-*-

// global settings for the GUI (maybe to be changeable later on).

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdio.h>
#include <stdlib.h>

#include <glm/glm.hpp>

// Framebuffer texture square side length
const int nmaxtex = 9;
const float FBO_tex_a[nmaxtex] = {256.f,362.f,512.f,724.f,1024.f,1448.f,2048.f,2896.f,4096.f};

// Background color for views
const float bgrgb[4] = {0.0f,0.0f,0.0f,1.0f};

// Light parameters
const glm::vec3 lightPos = {20.f,20.f,0.f};
const glm::vec3 lightColor = {1.f,1.f,1.f};
const float ambient = 0.2f;
const float diffuse = 0.4f;
const float specular = 0.6f;
const int shininess = 8;

// atom & bond resolution: 0 -> nmaxsph-1,nmaxcyl-1
const int isphres = 3;
const int icylres = 0;

// Mouse constants
const float mousesens_rot = 2.0f; // Mouse rotate sensitivity
const float mousesens_zoom = 0.15f; // Mouse zoom sensitivity
const float zfov = 45.f; // fov for the perspective
const float znear = 0.1f; // znear for the camera
const float zfar = 1000.f; // zfar for the camera
const float min_zoom = 1.f; // minimum distance to origin (zoom)

#endif SETTINGS_H

