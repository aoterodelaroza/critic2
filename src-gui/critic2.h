// Public interface to libcritic2.a

#ifndef CRITIC2_H
#define CRITIC2_H

namespace c2 {
  // atom type
  extern "C" struct c_atom {
    float x[3];
    float r[3];
    int is;
    int z;
    char name[11];
    int idx;
    int cidx;
    int flvec[3];
    int ifrag;
    float rad;
    float rgb[4];
    int ncon;
  };

  // routines
  extern "C" void gui_initialize();
  extern "C" int open_file(const char **filename, int ismolecule); 
  extern "C" void scene_initialize(int isc); 
  extern "C" void set_scene_pointers(int isc);
  extern "C" void gui_end();

  // number of scenes
  extern "C" int nsc;

  // pointers to the current scene
  extern "C" int isinit;
  extern "C" float scenerad;

  extern "C" int nat;
  extern "C" struct c_atom *at;

  extern "C" int mncon;
  extern "C" int *idcon;
  extern "C" int *lcon;

  extern "C" int nmol;
  extern "C" int *moldiscrete;

  extern "C" float *avec[3];
  extern "C" int ismolecule;
  extern "C" float *molx0;
  extern "C" float *molborder;
}

#endif
