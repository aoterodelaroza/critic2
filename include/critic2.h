// Public interface for src/gui_interface.f90

//xx// Variables made available through host association in 
//xx// gui_interface.f90. See gui_interface.f90 for documentation.
// Sticks
extern "C" struct c_stick {
  float r1[3];
  float r2[3];
  float rmid[3];
  float length;
  float thick;
  float rgb[3];
  float rot[4][4];
};

// Balls
extern "C" struct c_ball {
  float r[3];
  float rad;
  float rgb[3];
};

// Atoms
extern "C" struct c_atom {
  char name[11];
  int z;
  struct c_ball b;
};
extern "C" int nat;
extern "C" struct c_atom *at;

// Bonds
extern "C" struct c_bond {
  int i1;
  int i2;
  struct c_stick s;
};
extern "C" int nbond;
extern "C" struct c_bond *bond;

// Critical points
extern "C" struct c_critp {
  int type;
  char name[11];
  struct c_ball b;
};
extern "C" int ncritp;
extern "C" struct c_critp *critp;

// Unit cell
extern "C" float cell_x0[3];
extern "C" float cell_lat[3][3];
extern "C" int cell_nstick;
extern "C" struct c_stick cell_s[12];

// Bounding box
extern "C" float box_xmin[3];
extern "C" float box_xmax[3];
extern "C" float box_xcm[3];
extern "C" float box_xmaxlen;

//xx// Procedures made available in gui_interface.f90

// Initialize critic2
extern "C" void critic2_initialize();

// End the critic2 run
extern "C" void critic2_end();

// Read a new molecule/crystal from an external file
extern "C" void call_structure(const char **filename, int isMolecule); 

// Calculate critical points for the current field
extern "C" void call_auto(); 

// Update the data in the scene (fortran module variables) to make
// it available to the GUI 
extern "C" void update_scene(); 

// Clear the scene
extern "C" void clear_scene(bool unload); 

