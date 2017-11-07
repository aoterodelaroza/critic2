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

  // bond type
  extern "C" struct c_bond {
    float r1[3];
    float r2[3];
    float rgb1[4];
    float rgb2[4];
    float rad;
  };

  // routines
  extern "C" void gui_initialize(void *);
  extern "C" void open_file(const char **filename, int ismolecule); 
  extern "C" void set_scene_pointers(int isc);
  extern "C" void gui_end();

  // pointers to the current scene
  extern "C" int nat;
  extern "C" struct c_atom *at;
  extern "C" int nbond;
  extern "C" struct c_bond *bond;
  extern "C" float scenerad;
}

#endif
