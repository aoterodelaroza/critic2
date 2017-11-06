// -*-c++-*-

// global settings for the GUI (maybe to be changeable later on).

#ifndef SETTINGS_H
#define SETTINGS_H

// Framebuffer texture size
const float FBO_tex_x = 1024.f;
const float FBO_tex_y = 1024.f;

// Background color for views
const float bgrgb[4] = {0.0f,0.0f,0.0f,1.0f};

// atom & bond resolution
const int isphres = 0;

// Mouse constants
const float mousesens_rot = 0.001f; // Mouse rotate sensitivity
const float mousesens_zoom = 0.15f; // Mouse zoom sensitivity
const float zfov = 45.f; // fov for the perspective
const float znear = 0.1f; // znear for the camera
const float zfar = 1000.f; // zfar for the camera
const float min_zoom = 1.f; // minimum distance to origin (zoom)

#endif SETTINGS_H

