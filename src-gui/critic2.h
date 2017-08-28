// Public interface to libcritic2.a

#ifndef CRITIC2_H
#define CRITIC2_H

namespace c2 {
  extern "C" void gui_initialize(void *);
  extern "C" void draw_scene();
  extern "C" void gui_end();
}

#endif
