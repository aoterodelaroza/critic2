// Public interface to libcritic2.a

#ifndef CRITIC2_H
#define CRITIC2_H

namespace c2 {
  // atom type
  extern "C" struct c_atom {
    char label[11];
    int z;
    float r[3];
    float rad;
    float rgb[4];
  };

  // routines
  extern "C" void gui_initialize(void *);
  extern "C" void open_file(const char **filename, int ismolecule); 
  extern "C" void set_scene_pointers(int isc);
  extern "C" void gui_end();

  // pointers to the current scene
  extern "C" int nat;
  extern "C" struct c_atom *at;
  extern "C" float scenerad;
}

#endif
