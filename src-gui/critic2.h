// Public interface to libcritic2.a

#ifndef CRITIC2_H
#define CRITIC2_H

namespace c2 {
  extern "C" void gui_initialize(void *);
  extern "C" void open_file(const char **filename, int ismolecule); 
  extern "C" void draw_scene(int isc);
  extern "C" void gui_end();
}

#endif
