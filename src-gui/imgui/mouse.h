// -*-c++-*-

// Mouse state struct

#ifndef MOUSE_H
#define MOUSE_H

#include <glm/glm.hpp>

// A mouse class to encapsulate all information about the state of the
// mouse when it is hovering a certain window region.
struct MouseState
{
  bool hover   = false;
  bool lclick  = false;
  bool ldown   = false;
  bool mclick  = false;
  bool mdown   = false;
  bool rclick  = false;
  bool rdown   = false;
  bool ldrag   = false;
  bool rdrag   = false;
  bool ldclick = false;
  float scroll = 0.f;
  glm::vec2 pos = {0.f,0.f};
  glm::vec2 ndpos = {0.f,0.f};
};
#endif
