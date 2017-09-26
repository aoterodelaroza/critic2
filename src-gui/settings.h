// -*-c++-*-

// global settings for the GUI (maybe to be changeable later on).

#ifndef SETTINGS_H
#define SETTINGS_H

// Framebuffer texture size
const float FBO_tex_x = 1024.f;
const float FBO_tex_y = 1024.f;

// Mouse constants
const float mousesens_pan = 0.005f; // Mouse pan sensitivity
const float mousesens_rot = 0.001f; // Mouse rotate sensitivity
const float mousesens_zoom = 0.1f; // Mouse zoom sensitivity
const float znear = 0.1f; // znear for the camera and min. dist. to lookat

#endif SETTINGS_H

