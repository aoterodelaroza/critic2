/* Test program */

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <critic2.h>
#include <GL/gl3w.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw_gl3.h>
#include <matrix_math.h>
#include <tinyfiledialogs.h>
#include "geometry.h"

// #ifdef WIN32 //platform spisific sleep functions
// #include <synchapi.h>
// #endif // WIN32
// #if defined(LINUX) || defined(__APPLE__)
#include <unistd.h>
// #endif // LINUX || __APPLE__

using namespace std;

// Bond and atom resolutions (0 to 3)
static const char bondresolution = 2;
static const char atomresolution = 1;

static void ShowAppMainMenuBar(bool * show_bonds, bool * show_cps, bool * show_atoms);
static void ShowMenuFile();

///platform independent convertion vars into a string representation
string charConverter(float t) {
  char buffer[64];
  int len = sprintf(buffer, "%f", t);
  string val = buffer;
  return val;
}

string charConverter(size_t t) {
  char buffer[64];
  int len = sprintf(buffer, "%d", t);
  string val = buffer;
  return val;
}

string charConverter(int t) {
  char buffer[64];
  int len = sprintf(buffer, "%d", t);
  string val = buffer;
  return val;
}

//
//  Global Variables and Structs
//
///inpute vareables describing the current state of all used key and mouse inputs
struct {
  bool LeftMouseButton = 0;
  bool RightMouseButton = 0;
  bool MiddleMouseButton = 0;
  double ScrollYOffset = 0.f;
  double MPosX = 0.f;
  double MPosY = 0.f;
  double lastX = 0.f;
  double lastY = 0.f;
  double diffX = 0.f;
  double diffY = 0.f;
  bool ShiftKey = 0;
  bool CtrlKey = 0;
  bool AltKey = 0;
} input;

struct {
  Vector3f center;
  Vector3f dimensions;
} boundingCube;

///information about the current state of the camera
static CameraInfo cam;

// descriptions of the state of sevral Shader vareables
struct {
  GLuint gWorldLocation;
  GLuint gWVPLocation;
  GLuint vColorLocation;
  GLuint lColorLocation;
  GLuint lDirectionLocation;
  GLuint fAmbientIntensityLocation;
} ShaderVarLocations;

// all information regarding an atom
struct atom{
  string name = "";
  bool selected = false;
  int atomicNumber;
  Vector3f atomPosition;
  string atomTreeName; //must be saved to preserve imgui tree Id's
  int numberOfBonds;
  int * bonds;
  bool * neighCrystalBonds;
  int atomTreePosition;
};

// information about a bond
struct bond{
  atom * a1;
  atom * a2;
  Vector3f center;
  Matrix4f rotation;
  float length;
  bool neighCrystalBond;
};

/// information about a critical point
struct criticalPoint {
  Vector3f cpPosition;
  int type;
  string typeName = "";
  bool selected = false;
};

bond * bonds;
criticalPoint * loadedCriticalPoints;

atom * loadedAtoms;
int loadedAtomsAmount = 0; //number of loaded atoms (equal to what would be loadedAtoms.length)
int selectedAtom = 0; //atom selected by tree or search bar
int bondsAmount = 0; //the number of bonds in the structure
int loadedCPAmount = 0; //the number of critical points in the structure
int selectedCP = 0; //critical point selected by the tree

// Standard error print to console
static void error_callback(int error, const char* description)
{
  fprintf(stderr, "Error %d: %s\n", error, description);
}

/// add a shader to the gl program
static void AddShader(GLuint ShaderProgram, const char* pShaderText, GLenum ShaderType)
{
  GLuint ShaderObj = glCreateShader(ShaderType);

  if (ShaderObj == 0) {
    fprintf(stderr, "Error creating shader type %d\n", ShaderType);
    exit(0);
  }

  const GLchar * p[1];
  p[0] = pShaderText;
  GLint Lengths[1];
  Lengths[0] = strlen(pShaderText);
  glShaderSource(ShaderObj, 1, p, Lengths);
  glCompileShader(ShaderObj);
  GLint success;
  glGetShaderiv(ShaderObj, GL_COMPILE_STATUS, &success);
  if (!success) {
    GLchar InfoLog[1024];
    glGetShaderInfoLog(ShaderObj, 1024, NULL, InfoLog);
    fprintf(stderr, "Error compiling shader type %d: '%s'\n", ShaderType, InfoLog);
    exit(1);
  }
  glAttachShader(ShaderProgram, ShaderObj);
}

// create a shader programe based on a shader script
static GLuint LightingShader()
{

  GLuint ShaderProgram = glCreateProgram();
  if (ShaderProgram == 0){
    exit(1);
  }

  const char * vs = "#version 330 \n \
    uniform mat4 gWorld; \n \
    uniform mat4 gWVP; \n \
    layout (location = 0) in vec3 inPosition; \n \
    layout (location = 1) in vec3 inNormal; \n \
    smooth out vec3 vNormal; \n \
    void main() { \n \
      gl_Position = gWVP * vec4(inPosition, 1.0); \n \
      vNormal = (gWorld * vec4(inNormal, 0.0)).xyz; \n \
      }";

  const char * fs = "#version 330 \n \
    smooth in vec3 vNormal; \n \
    uniform vec4 vColor; \n \
    out vec4 outputColor; \n \
    uniform vec3 lColor; \n \
    uniform vec3 lDirection; \n \
    uniform float fAmbientIntensity; \n \
    void main() { \n \
      float fDiffuseIntensity = max(0.0, dot(normalize(vNormal), lDirection)); \n \
      outputColor = vColor; \n \
      }";

  //outputColor = vColor * vec4(lColor * (fAmbientIntensity+fDiffuseIntensity), 1.0);
  //      outputColor = vColor;

  AddShader(ShaderProgram, vs, GL_VERTEX_SHADER);
  AddShader(ShaderProgram, fs, GL_FRAGMENT_SHADER);

  GLint success = 0;

  glLinkProgram(ShaderProgram);
  glGetProgramiv(ShaderProgram, GL_LINK_STATUS, &success);
  if (success == 0) exit(1);

  glValidateProgram(ShaderProgram);
  glGetProgramiv(ShaderProgram, GL_VALIDATE_STATUS, &success);
  if (success == 0) exit(1);

  return ShaderProgram;
}

/// create a mesh object that can be reused for evrey object that uses the same mesh
void CreateAndFillBuffers(GLuint * VertexBuffer, GLuint * IndexBuffer,
                          GLfloat * Vertices, unsigned int * Indices,
                          unsigned int NumVertices, unsigned int NumIndices)
{
  glGenBuffers(1, VertexBuffer);
  glBindBuffer(GL_ARRAY_BUFFER, *VertexBuffer);
  glBufferData(GL_ARRAY_BUFFER, NumVertices*sizeof(GLfloat), Vertices, GL_STATIC_DRAW);

  glGenBuffers(1, IndexBuffer);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *IndexBuffer);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, NumIndices*sizeof(unsigned int), Indices, GL_STATIC_DRAW);

}

// get mouse scroll
void ScrollCallback(GLFWwindow * window, double xoffset, double yoffset)
{
  float camZoomFactor = boundingCube.dimensions.Length() * 0.2f;
  cam.Pos[2] += yoffset * camZoomFactor;
}

// rotate the bonds into place on screen based on the atoms specifiyed
void GenerateBondInfo(bond * b, atom * a1, atom * a2)
{
  b->a1 = a1;
  b->a2 = a2;

  Vector3f mid = (a1->atomPosition + a2->atomPosition)/2;
  Vector3f q = a1->atomPosition - a2->atomPosition;
  float d = q.Length()/2;

  b->length = d;
  b->center = mid;

  Vector3f n_q = Vector3f(q);
  n_q.Normalize();

  Vector3f z_axis = Vector3f(0, 0, 1);
  z_axis.Normalize();
  Vector3f axis = z_axis.Cross(n_q);
  axis.Normalize();
  float angle = acosf(z_axis.Dot(n_q));

  Matrix4f Rot;
  Rot.InitRotateAxisTransform(axis, angle);

  b->rotation = Rot;
}

// draw a bond between 2 atoms defined in the bond struct
void DrawBond(Pipeline * p, GLuint CylVB, GLuint CylIB, bond * b)
{
  float grey[3] = {.5, .5, .5};
  float white[3] = {1, 1, 1};

  p->Scale(0.05f, 0.05f, b->length);

  Vector3f pos = b->center - boundingCube.center;
  p->Translate(pos.x, pos.y, pos.z);
  p->SetRotationMatrix(b->rotation);

  float dir[3] = {cam.Target[0], cam.Target[1], cam.Target[2]};
  glUniformMatrix4fv(ShaderVarLocations.gWVPLocation, 1, GL_TRUE,
                     (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(ShaderVarLocations.gWorldLocation, 1, GL_TRUE,
                     (const GLfloat *)p->GetWorldTrans());
  glUniform4fv(ShaderVarLocations.vColorLocation, 1, (const GLfloat *)&grey);
  glUniform4fv(ShaderVarLocations.lColorLocation, 1, (const GLfloat *)&white);
  glUniform4fv(ShaderVarLocations.lDirectionLocation, 1, (const GLfloat *)&dir);
  glUniform1f(ShaderVarLocations.fAmbientIntensityLocation, 0.8);

  glBindBuffer(GL_ARRAY_BUFFER, CylVB);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, CylIB);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glDrawElements(GL_TRIANGLES, 240, GL_UNSIGNED_INT, 0);
}

//global vars for an atom's mesh
GLuint atomVB; //atom vertacies
GLuint atomIB; //atom indecies ~(direction of verts)
unsigned int numbIndeces;

///remove refrences to the currently loaded molecule
void destructLoadedMolecule(){
  if (loadedAtomsAmount > 0){
    for (int i=loadedAtomsAmount-1; i>=0; i--){
      delete loadedAtoms[i].bonds;
    }
    loadedAtomsAmount = 0;
  }
  if (bondsAmount > 0){
    delete bonds;
    bondsAmount = 0;
  }
}

///remove refrences to the currently loaded critical points
void destructCriticalPoints() {
  if (loadedCPAmount > 0) {
    loadedCriticalPoints = NULL;
    loadedCPAmount = 0;
  }
}


/// load atoms using the critic2 external C functions defined in interface.f90
void loadAtoms() {
  //fill loadedAtoms array
  int n;
  get_num_atoms(&n);

  loadedAtomsAmount = n;
  loadedAtoms = new atom[loadedAtomsAmount];

  float minX = FLT_MAX;
  float maxX = -FLT_MAX;
  float minY = FLT_MAX;
  float maxY = -FLT_MAX;
  float minZ = FLT_MAX;
  float maxZ = -FLT_MAX;

  for (int i=0; i < n; i++) {
    int atomicN;
    double x;
    double y;
    double z;

    get_atom_position(i+1, &atomicN, &x, &y, &z);

    loadedAtoms[i].atomicNumber = atomicN;
    loadedAtoms[i].atomPosition = Vector3f(x, y, z);


    // calculate boundary and of molecule
    if (x < minX){
      minX = x;
    } else if (x > maxX){
      maxX = x;
    }
    if (y < minY){
      minY = y;
    } else if (y > maxY){
      maxY = y;
    }
    if (z < minZ){
      minZ = z;
    } else if (z > maxZ){
      maxZ = z;
    }
  }

  //calculate center point of and size of molecule
  boundingCube.center = Vector3f((minX + maxX)/2, (minY + maxY)/2, (minZ + maxZ)/2);
  boundingCube.dimensions = Vector3f((maxX - minX), (maxY - minY), (maxZ - minZ));

  //tree names must be constant
  for (size_t x = 0; x < loadedAtomsAmount; x++) {
    std::string nodeName = "";
    nodeName += "Elem Name: ";
    nodeName += loadedAtoms[x].name;
    nodeName += "Atomic #:";
    nodeName += charConverter(loadedAtoms[x].atomicNumber);
    nodeName += "  ID: ";
    nodeName += charConverter(x);
    nodeName += "##TreeID = "; //extra info for imgui to find selection
    nodeName += charConverter(x);
    loadedAtoms[x].atomTreeName = nodeName;
    loadedAtoms[x].atomTreePosition = x;
  }

}

/// load bonds using the critic2 external C functions defined in interface.f90
void loadBonds() {

  for (int i = 0; i < loadedAtomsAmount; i++) {
    int nstarN;

    num_of_bonds(i+1, &nstarN);

    loadedAtoms[i].numberOfBonds = nstarN;
    loadedAtoms[i].bonds = new int[nstarN];
    loadedAtoms[i].neighCrystalBonds = new bool[nstarN];
    bondsAmount += nstarN;

    for (int j = 0; j < nstarN; j++) {
      int connected_atom;
      bool neighCrystal = false;

      get_atom_bond(i+1, j+1, &connected_atom, &neighCrystal);
      loadedAtoms[i].bonds[j] = connected_atom-1;
      loadedAtoms[i].neighCrystalBonds[j] = neighCrystal;
    }
  }

  bonds = new bond[bondsAmount];
  int bondidx = 0;
  for (int i=0; i<loadedAtomsAmount; i++){
    int numBonds = loadedAtoms[i].numberOfBonds;
    for (int j=0; j<numBonds; j++){
      bonds[bondidx].neighCrystalBond = loadedAtoms[i].neighCrystalBonds[j];
      GenerateBondInfo(&bonds[bondidx], &loadedAtoms[i], &loadedAtoms[loadedAtoms[i].bonds[j]]);
      bondidx += 1;
    }
  }
}

/// load critical points using the critic2 external C functions defined in interface.f90
void loadCriticalPoints() {
  int numCP;
  num_of_crit_points(&numCP);

  loadedCPAmount += (numCP - loadedAtomsAmount);
  loadedCriticalPoints = new criticalPoint[loadedCPAmount];

  for (int i = loadedAtomsAmount + 1; i <= numCP; i++) {
    int cpType;
    double x;
    double y;
    double z;

    get_cp_pos_type(i, &cpType, &x, &y, &z);

    loadedCriticalPoints[(i-(loadedAtomsAmount+1))].cpPosition = Vector3f(x, y, z);

    loadedCriticalPoints[(i-(loadedAtomsAmount+1))].type = cpType;

    if (cpType == -1) {
      loadedCriticalPoints[(i-(loadedAtomsAmount+1))].typeName += "bond";
    } else if (cpType == 1) {
      loadedCriticalPoints[(i-(loadedAtomsAmount+1))].typeName += "ring";
    } else if (cpType == 3) {
      loadedCriticalPoints[(i-(loadedAtomsAmount+1))].typeName += "cage";
    }
  }
}

// iterates though all bonds[] and draws them
void drawAllBonds(Pipeline * p, GLuint CylVB, GLuint CylIB)
{
  for (int i=0; i< bondsAmount; i++){
    if (!bonds[i].neighCrystalBond) {
      DrawBond(p, CylVB, CylIB, &bonds[i]);
    }
  }
}

/// setting color of meshes based on some input
///returns the color of an atom based on the atomic number
///and desired color Intesity (brightness)
Vector3f getAtomColor(int atomicNumber) {
  Vector3f color;
  if (atomicNumber == 1) {     // Hydrogen = white
    color = Vector3f(1, 1, 1);
  } else if (atomicNumber == 6) {   // Carbon = black
    color = Vector3f(0, 0, 0);
  } else if (atomicNumber == 7) {   // Nitrogen = dark blue
    color = Vector3f(0, 0, 0.5);
  } else if (atomicNumber == 8) {   // Oxygen = red
    color = Vector3f(1, 0, 0);
  } else if (atomicNumber == 9 || atomicNumber == 17) {   // Fluorine & Chlorine = green
    color = Vector3f(0, 1, 0);
  } else if (atomicNumber == 35) {  // Bromine = dark red
    color = Vector3f(0.5, 0, 0);
  } else if (atomicNumber == 53) {  // Iodine = dark violet
    color = Vector3f(0.5, 0, 0.5);
  } else if (atomicNumber == 2 ||    // noble gases (He, Ne, Ar, Kr, Xe, Rn) = cyan
	     atomicNumber == 10 ||
	     atomicNumber == 18 ||
	     atomicNumber == 36 ||
	     atomicNumber == 54 ||
	     atomicNumber == 86) {
    color = Vector3f(0, 1, 1);
  } else if (atomicNumber == 15) {  // Potassium = orange
    color = Vector3f(1, 0.5, 0);
  } else if (atomicNumber == 16) {  // sulfur = yellow
    color = Vector3f(1, 1, 0);
  } else if (atomicNumber == 3 ||   // alkali metals = violet
             atomicNumber == 11 ||
             atomicNumber == 19 ||
             atomicNumber == 37 ||
             atomicNumber == 55 ||
             atomicNumber == 87){
    color = Vector3f(1, 0, 1);
  } else if (atomicNumber == 4 ||   // alkaline earth metals = dark green
             atomicNumber == 12 ||
             atomicNumber == 20 ||
             atomicNumber == 38 ||
             atomicNumber == 56 ||
             atomicNumber == 88){
    color = Vector3f(0, 0.5, 0);
  } else if (atomicNumber == 81) {  // titaniam = gray
    color = Vector3f(0.75, 0.75, 0.75);
  } else if (atomicNumber == 26) {  // iron = dark orange
    color = Vector3f(0.75, 0.25, 0);
  } else if ((atomicNumber >= 21 && atomicNumber <= 30) ||  // transition metals = orange/pink
             (atomicNumber >= 39 && atomicNumber <= 48) ||
             (atomicNumber >= 57 && atomicNumber <= 80) ||
             (atomicNumber >= 89 && atomicNumber <= 112)){
    color = Vector3f(1, 0.5, 0.3);
  } else  {   // all other atoms = pink
    color = Vector3f(1, 0.25, 0.5);
  }
  return color;
}

///Set critical point color based on type
Vector3f getCritPointColor(int cpType) {
  Vector3f color;
  if (cpType == -1) {     // bond cp = yellow
    color = Vector3f(1, 1, 0);
  } else if (cpType == 1) {   // ring cp = light grey
    color = Vector3f(0.6588, 0.6588, 0.6588);
  } else  {   // cage cp = light purple
    color = Vector3f(0.94, 0.81, 0.99);
  }
  return color;
}

//describes what each atom color is
//basic color Legend
void atomLegend_displayColorBoxAndName(int colorNumber, string * colorNames) {
  Vector3f AtomicNumberColor = getAtomColor(colorNumber);
  ImVec4 color = ImColor(AtomicNumberColor.x, AtomicNumberColor.y, AtomicNumberColor.z, 1.0);
  ImGui::ColorButton(color); ImGui::SameLine();
  ImGui::Text(colorNames[colorNumber].c_str());
}

/// create a menu describing how atoms are colored
void atomColorLegend() {
  bool p_open = false;
  if (ImGui::CollapsingHeader("atom color legend")) {
    string * colorNames = new string[101];
    colorNames[1] = "Hydrogen";
    colorNames[2] = "Noble Gas";
    colorNames[3] = "Alkali Metals";
    colorNames[4] = "Akaline Earth";

    colorNames[6] = "Carbon";
    colorNames[7] = "Nitrogen";
    colorNames[8] = "Oxygen";
    colorNames[9] = "Flourine & Chlorine";

    colorNames[15] = "Potassium";
    colorNames[16] = "Sulfur";

    colorNames[21] = "transition Metal";

    colorNames[26] = "iron";

    colorNames[35] = "Bromine";

    colorNames[53] = "Iodine";

    colorNames[81] = "titaniam";

    colorNames[100] = "other";

    for (size_t i = 1; i < 5; i++) {
      atomLegend_displayColorBoxAndName(i, colorNames);
    }
    for (size_t i = 6; i < 10; i++) {
      atomLegend_displayColorBoxAndName(i, colorNames);
    }

    atomLegend_displayColorBoxAndName(15, colorNames);
    atomLegend_displayColorBoxAndName(16, colorNames);
    atomLegend_displayColorBoxAndName(21, colorNames);
    atomLegend_displayColorBoxAndName(26, colorNames);
    atomLegend_displayColorBoxAndName(35, colorNames);
    atomLegend_displayColorBoxAndName(53, colorNames);
    atomLegend_displayColorBoxAndName(81, colorNames);
    atomLegend_displayColorBoxAndName(100, colorNames);
  }
}

// ///will be used to draw atom number over the atom using imgui window
// float* getScreenPositionOfVertex(float *vertexLocation) {
// 	//TODO transfrom from vertex location to screen location
// 	return NULL;
// }

/// draw an atom using gl functions
void drawAtomInstance(int id, Vector3f posVector, Vector3f color,
                      Pipeline * p, GLuint SphereVB, GLuint SphereIB) {

  //if atom is selected, brighten it
  if (loadedAtoms[id].selected) {
    color = color * 1.5;
  }

  float scaleAmount = (float)loadedAtoms[id].atomicNumber;
  // float scaleAmount = 3;
  if (scaleAmount < 4.0f) {
    scaleAmount = 0.2f;
  } else {
    scaleAmount = 0.4f;
  }
  p->Scale(scaleAmount, scaleAmount, scaleAmount);

  Vector3f pos = posVector - boundingCube.center;
  p->Translate(pos.x, pos.y, pos.z);
  p->Rotate(0.f, 0.f, 0.f);

  glUniformMatrix4fv(ShaderVarLocations.gWVPLocation, 1, GL_TRUE,
                     (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(ShaderVarLocations.gWorldLocation, 1, GL_TRUE,
                     (const GLfloat *)p->GetWorldTrans());
  glBindBuffer(GL_ARRAY_BUFFER, SphereVB);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, SphereIB);
  glUniform4fv(ShaderVarLocations.vColorLocation, 1, (const GLfloat *)&color);
  glDrawElements(GL_TRIANGLES, 6144, GL_UNSIGNED_INT, 0);
  /*
  //TODO draw atom ID number
  ImGui::SetNextWindowSize(ImVec2(5, 5), ImGuiSetCond_Always);
  ImGui::SetNextWindowCollapsed(true);
  //float * winPos;
  //matrix math to transform posVector to pixel location of an atoms center

  //ImGui::SetNextWindowPos(ImVec2(winPos[0], winPos[1])); //TODO set location of identifying number
  ImGui::Begin(std::to_string(identifyer).c_str(), false);
  ImGui::End();
  */
}

/// draw a critical point using gl functions
void drawCritPointInstance(int identifier, Vector3f posVector, const GLfloat color[4],
			   Pipeline * p, GLuint SphereVB, GLuint SphereIB) {
  //selection start
  float inc = 1.f;
  if (loadedCriticalPoints[identifier].selected) { //selection is color based
    inc = 1.5f;
  }

  GLfloat n_Color[] = {color[0] * inc,color[1] * inc,color[2] * inc, color[3]};
  //selection end

  float scaleAmount = 0.1f;
  p->Scale(scaleAmount, scaleAmount, scaleAmount);
  Vector3f pos = posVector - boundingCube.center;
  p->Translate(pos.x, pos.y, pos.z);

  p->Rotate(0.f, 0.f, 0.f); //no rotation required
  glUniformMatrix4fv(ShaderVarLocations.gWVPLocation, 1, GL_TRUE, (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(ShaderVarLocations.gWorldLocation, 1, GL_TRUE, (const GLfloat *)p->GetWorldTrans());
  glBindBuffer(GL_ARRAY_BUFFER, SphereVB);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, SphereIB);
  glUniform4fv(ShaderVarLocations.vColorLocation, 1, (const GLfloat *)&n_Color);
  glDrawElements(GL_TRIANGLES, 6144, GL_UNSIGNED_INT, 0);
  /*
  //TODO draw crit point ID number or type?
  ImGui::SetNextWindowSize(ImVec2(5, 5), ImGuiSetCond_Always);
  ImGui::SetNextWindowCollapsed(true);
  //float * winPos;
  //matrix math to transform posVector to pixel location of an atoms center

  //ImGui::SetNextWindowPos(ImVec2(winPos[0], winPos[1])); //TODO set location of identifying number
  ImGui::Begin(charConverter(identifyer).c_str(), false);
  ImGui::End();
  */
}

bool flashAtoms = false; // toggle with selection toggles (in gui)
//selctedAtom from tree selection
//the number of frames the deselected atoms stay invisable
int framesMax = 15; // ~0.5 seconds
int framesLeft = 0;
bool otherAtomsVisable = true;


bool flashCP = false;
int framesMaxCP = 15;
int framesLeftCP = 0;
bool otherCriticalPointsVisable = true;


///draws all atoms in the loadedAtoms struct
void drawAllAtoms(Pipeline * p, GLuint SphereVB, GLuint SphereIB) {
  if (flashAtoms && loadedAtomsAmount > 0) { //flash mode
    if (framesLeft <= 0) {
      otherAtomsVisable = !otherAtomsVisable;
      framesLeft = framesMax;
    }
    if (otherAtomsVisable) { //all visible
      for (size_t x = 0; x < loadedAtomsAmount; x++) {
	Vector3f color = getAtomColor(loadedAtoms[x].atomicNumber);
	drawAtomInstance(x, loadedAtoms[x].atomPosition, color, p, SphereVB, SphereIB);
      }
    } else { // only selected atom visable
      Vector3f color = getAtomColor(loadedAtoms[selectedAtom].atomicNumber);
      drawAtomInstance(selectedAtom, loadedAtoms[selectedAtom].atomPosition, color, p, SphereVB, SphereIB);
    }
    framesLeft--;
  } else { // regular drawing
    for (size_t x = 0; x < loadedAtomsAmount; x++){
      Vector3f color = getAtomColor(loadedAtoms[x].atomicNumber);
      drawAtomInstance(x, loadedAtoms[x].atomPosition, color, p, SphereVB, SphereIB);
    }
  }
}

///draws all loaded critical points
void drawAllCPs(Pipeline * p, GLuint SphereVB, GLuint SphereIB) {
  //cp flashing
  if (flashCP) {
    if (framesLeftCP <= 0) {
      otherCriticalPointsVisable = !otherCriticalPointsVisable;
      framesLeftCP = framesMaxCP;
    }
    if(otherCriticalPointsVisable)
      for (int x = 0; x < loadedCPAmount; x++) {
	Vector3f color = getCritPointColor(loadedCriticalPoints[x].type);
	drawCritPointInstance(x, loadedCriticalPoints[x].cpPosition, color, p, SphereVB, SphereIB);
      }
    else {
      Vector3f color = getCritPointColor(loadedCriticalPoints[selectedCP].type);
      drawCritPointInstance(selectedCP, loadedCriticalPoints[selectedCP].cpPosition, color, p, SphereVB, SphereIB);
    }
    framesLeftCP--;
  }
  else {
    for (int x = 0; x < loadedCPAmount; x++) {
      Vector3f color = getCritPointColor(loadedCriticalPoints[x].type);
      drawCritPointInstance(x, loadedCriticalPoints[x].cpPosition, color, p, SphereVB, SphereIB);
    }
  }
}

/// moves cam over atom (alligned to z axis)
void lookAtAtom(int atomNumber) {
  cam.Pos[0] = loadedAtoms[atomNumber].atomPosition.x;
  cam.Pos[1] = loadedAtoms[atomNumber].atomPosition.y;
}

/// moves cam over crit point (alligned to z axis)
void lookAtCritPoint(int critPointNum) {
  cam.Pos[0] = loadedCriticalPoints[critPointNum].cpPosition.x;
  cam.Pos[1] = loadedCriticalPoints[critPointNum].cpPosition.y;
}


///information to display in the stats list
///selects an atom focusing the view and displaying additonal info
///in the
void selectAtom(int atomIndex) {
  loadedAtoms[atomIndex].selected = true;
  lookAtAtom(atomIndex);
  selectedAtom = atomIndex;
}

void selectCriticalPoint(int cpIndex) {
  loadedCriticalPoints[cpIndex].selected = true;
  lookAtCritPoint(cpIndex);
  selectedCP = cpIndex;
}


//This methods are used to display additonal information about
//a perticular critical point or atom
//number of displayVars is currently assumed to be 3
void displayCol(string * displayStats, int numberOfCol) {
  for (size_t i = 0; i < numberOfCol; i++) {
    //text is wraped if too large
    ImGui::TextWrapped(displayStats[i].c_str());
    ImGui::NextColumn();
  }
}

void atomBondAmountInfo(string * displayVars, int atomNumber) {
  displayVars[0] = "number of bonds";
  displayVars[1] = charConverter(loadedAtoms[atomNumber].numberOfBonds);
  displayVars[2] = "";
}

void atomAtomicNumberInfo(string * displayVars, int atomNumber) {
  displayVars[0] = "atomic number";
  displayVars[1] = charConverter(loadedAtoms[atomNumber].atomicNumber);
  displayVars[2] = "";
}

void criticalPointTypeInfo(string * displayVars, int criticalPointIndex) {
  displayVars[0] = "critical point Type";
  displayVars[1] = loadedCriticalPoints[criticalPointIndex].typeName;
  displayVars[2] = "";
}



/**
   display any information about the currently selected atom
   this is A menu bar item and must be called in a window
*/
void drawSelectedAtomStats() {
  if (loadedAtomsAmount == 0) {
    return;
  }
  if (ImGui::CollapsingHeader("Selected atom information")) {
    const int numberOfColumns = 3;

    ImGui::Columns(numberOfColumns, "mycolumns");
    ImGui::Separator();

    //----column names
    ImGui::Text("Info Type"); ImGui::NextColumn();
    ImGui::Text("Value1"); ImGui::NextColumn();
    ImGui::Text("Value2"); ImGui::NextColumn();
    //-----

    ImGui::Separator();

    string displayStats[numberOfColumns];

    //loading stats into string array and displaying them
    atomBondAmountInfo(displayStats, selectedAtom);
    displayCol(displayStats, numberOfColumns);

    atomAtomicNumberInfo(displayStats, selectedAtom);
    displayCol(displayStats, numberOfColumns);

    ImGui::Columns(1);
  }

}

void drawSelectedCPStats() {
  if (loadedCPAmount == 0 || selectedCP < 0) {
    return;
  }
  if (selectedCP >= loadedCPAmount) {
    // cout << "error in cp stat display: selected greater than number" << endl;
    return;
  }

  bool p_open = false;

  if (ImGui::CollapsingHeader("critical point information")) {
    const int numberOfColums = 3;

    ImGui::Columns(numberOfColums, "mycolumns");
    ImGui::Separator();
    ImGui::Text("Info Type"); ImGui::NextColumn();
    ImGui::Text("Value1"); ImGui::NextColumn();
    ImGui::Text("Value2"); ImGui::NextColumn();
    ImGui::Separator();

    string displayStats[numberOfColums];

    criticalPointTypeInfo(displayStats, selectedCP);
    displayCol(displayStats, numberOfColums);

    //CP stat 2
    //displayCol(displayStats, numberOfColums);
    ImGui::Columns(1);
  }
}

// /// print information about the current state of the camra
// void printCamStats() {
// 	bool p_open = false;
// 	ImGui::SetNextWindowSize(ImVec2(300, 75), ImGuiSetCond_Appearing);
// 	ImGui::Begin("cam stats", &p_open);
// 	string camPos = "cam pos: " + (string)charConverter(cam.Pos[0]) + "," + charConverter(cam.Pos[1]) + "," + charConverter(cam.Pos[2]);
// 	string camTarget = "cam target: " + (string)charConverter(cam.Target[0]) + "," + charConverter(cam.Target[1]) + "," + charConverter(cam.Target[2]);
// 	string camUp = "cam up: " + (string)charConverter(cam.Up[0]) + "," + charConverter(cam.Up[1]) + "," + charConverter(cam.Up[2]);
// 
// 	ImGui::Text(camPos.c_str());
// 	ImGui::Text(camTarget.c_str());
// 	ImGui::Text(camUp.c_str());
// 
// 
// 	ImGui::End();
// }
// 
// /// draw menu items in a toolbar, currently replaced with dropdowns
// void drawToolBar(int screen_w, int screen_h,
//                  bool * show_bonds, bool * show_cps, bool * show_atoms) {
// 	ImGui::SetNextWindowSize(ImVec2(50, screen_h),ImGuiSetCond_Once);
// 	ImGui::SetNextWindowPos(ImVec2(0, 0),ImGuiSetCond_Always);
//   ImGuiWindowFlags flags = 0;
//     flags |= ImGuiWindowFlags_AlwaysAutoResize;
//     flags |= ImGuiWindowFlags_NoResize;
//     flags |= ImGuiWindowFlags_NoMove;
//     flags |= ImGuiWindowFlags_NoTitleBar;
// 
// 	bool p_open = false;
// 	ImGui::Begin("ToolBar",&p_open, flags);
//   if (ImGui::Button("Load Molecule")){
//       char const * lTheOpenFileName = tinyfd_openFileDialog(
//     		"Select Molecule file",
//     		"../../examples/data/benzene.wfx",
//     		0,
//         NULL,
//     		NULL,
//     		0);
// 
//       if (lTheOpenFileName == NULL) {
//         return;
//       }
// 
//       init_struct();
//       call_structure(lTheOpenFileName, (int) strlen(lTheOpenFileName), 1);
//       destructLoadedMolecule();
//       destructCriticalPoints();
//       loadAtoms();
//       loadBonds();
//   }
//   if (ImGui::Button("Load Crystal")){
//       char const * lTheOpenFileName = tinyfd_openFileDialog(
//     		"Select Molecule file",
//     		"../../examples/data/ammonia.big.vel.cube",
//     		0,
//     		NULL,
//     		NULL,
//     		0);
// 
//       if (lTheOpenFileName == NULL) {
//         return;
//       }
// 
//       init_struct();
//       call_structure(lTheOpenFileName, (int) strlen(lTheOpenFileName), 0);
//       destructLoadedMolecule();
//       destructCriticalPoints();
//       loadAtoms();
//       loadBonds();
//   }
//   if (ImGui::Button("Generate Critical Points")){
//     auto_cp();
//   }
//   if (ImGui::Button("Load Critical Points")) {
//     destructCriticalPoints();
//     loadCriticalPoints();
//   }
//   if (ImGui::Button("Clear")){
//     destructLoadedMolecule();
//     destructCriticalPoints();
//   }
//   ImGui::Checkbox("Bonds", show_bonds);
//   ImGui::Checkbox("Crit Pts", show_cps);
//   ImGui::Checkbox("Atoms", show_atoms);
// 
//   if (ImGui::Checkbox("`elc.find", &flashAtoms)) { //set flashing to defults
// 	  framesMax = 15; // ~0.5 seconds
// 	  framesLeft = 0;
// 	  otherAtomsVisable = true;
//   }
// 
//   ImGui::End();
// }

/// This is the main imgui window describing atoms,bonds,critical points, and additonal information
void drawMainMenuTree(int screen_w, int screen_h) {
  ImGui::SetNextWindowSize(ImVec2(300, screen_h),ImGuiSetCond_Always);
  ImGui::SetNextWindowPos(ImVec2(screen_w-300, 0),ImGuiSetCond_Always);
  ImGuiWindowFlags flags = 0;
  //    flags |= ImGuiWindowFlags_AlwaysAutoResize;
  flags |= ImGuiWindowFlags_NoResize;
  flags |= ImGuiWindowFlags_NoMove;

  bool p_open = false;
  ImGui::Begin("Tree View",&p_open, flags);
  if (ImGui::CollapsingHeader("atoms and bonds")) {
    int closeOthers = -1;
    for (size_t x = 0; x < loadedAtomsAmount; x++) {
      if (ImGui::TreeNode(loadedAtoms[x].atomTreeName.c_str())) {
	if (loadedAtoms[x].selected == false) {
	  closeOthers = x;
	  selectAtom(x);
	}

	//selection based on atoms bonds
	for (size_t i = 0; i < bondsAmount; i++) {
	  if (bonds[i].a1->atomTreePosition == loadedAtoms[x].atomTreePosition) {
	    string bondName = "bondedTo" + charConverter(bonds[i].a2->atomTreePosition);
	    if (ImGui::TreeNode(bondName.c_str())) {
	      //select a2
	      closeOthers = bonds[i].a2->atomTreePosition;
	      // cout << "going to atom" + charConverter(closeOthers) << endl;
	      ImGui::TreePop();
	    }
	  }
	  else if (bonds[i].a2->atomTreePosition == loadedAtoms[x].atomTreePosition) {
	    string bondName = "bondedTo" + charConverter(bonds[i].a1->atomTreePosition);
	    if (ImGui::TreeNode(bondName.c_str())) {
	      //select a1
	      closeOthers = bonds[i].a1->atomTreePosition;
	      // cout << "going to atom" + charConverter(closeOthers) << endl;
	      ImGui::TreePop();
	    }
	  }
	}

	ImGui::TreePop();
      }
      else {
	loadedAtoms[x].selected = false;
      }
    }

    if (closeOthers != -1) {
      for (int i = 0; i < loadedAtomsAmount; i++){
	ImGui::GetStateStorage()->SetInt(ImGui::GetID(loadedAtoms[i].atomTreeName.c_str()), 0); // close all tabs
      }
      ImGui::GetStateStorage()->SetInt(ImGui::GetID(loadedAtoms[selectedAtom].atomTreeName.c_str()), 1);
      closeOthers = -1;
    }

  }
  if (loadedCPAmount != 0) {
    if (ImGui::CollapsingHeader("Critical Points")) {
      int closeOthers = -1;

      for (size_t i = 0; i < loadedCPAmount; i++) {
	if (ImGui::TreeNode((loadedCriticalPoints[i].typeName + ":" + charConverter(i)).c_str())) { //critical point tree node
	  if (loadedCriticalPoints[i].selected == false) {
	    loadedCriticalPoints[i].selected = true;
	    lookAtCritPoint(i);
	    closeOthers = i;
	    selectedCP = i;
	  }
	  ImGui::TreePop();
	}
	else {
	  loadedCriticalPoints[i].selected = false;
	}
      }

      if (closeOthers != -1) { //only one cp tab should be open at a time
	for (int i = 0; i < loadedAtomsAmount; i++) {
	  ImGui::GetStateStorage()->SetInt(ImGui::GetID((loadedCriticalPoints[i].typeName + ":" + charConverter(i)).c_str()), 0); // close all tabs
	}
	ImGui::GetStateStorage()->SetInt(ImGui::GetID((loadedCriticalPoints[closeOthers].typeName + ":" + charConverter(closeOthers)).c_str()), 1); // leave selected tab open
	closeOthers = -1;
      }
    }

    // start of cp info window
    drawSelectedCPStats();

  }
  drawSelectedAtomStats();
  atomColorLegend();


  ImGui::End();
}

///search bar to find atoms by atomic number
void createAtomSearchBar() {
  ImGui::SetNextWindowSize(ImVec2(200, 70), ImGuiSetCond_Appearing);
  ImGui::Begin("Atom Search by atomic #");
  ImGuiWindowFlags Flags;
  int atomAtomicNumber = 0;
  if (ImGui::InputInt("atomic#", &atomAtomicNumber, 0, 0, ImGuiInputTextFlags_EnterReturnsTrue)) { //search on enter
    if (atomAtomicNumber > 0 && atomAtomicNumber < loadedAtomsAmount) {
      selectAtom(atomAtomicNumber);
    }
  }

  //checkbox to flash atoms
  if (ImGui::Checkbox("Selc.find", &flashAtoms)) { //set flashing to defults
    framesMax = 15; // ~0.5 seconds
    framesLeft = 0;
    otherAtomsVisable = true;
  }

  ImGui::End();
}

void createCriticalPointSearchBar() {
  ImGui::SetNextWindowSize(ImVec2(200, 70), ImGuiSetCond_Appearing);
  ImGui::Begin("Critical Search by ID #");
  ImGuiWindowFlags Flags;
  int criticalpointID = 0;
  if (ImGui::InputInt("CP_ID#", &criticalpointID, 0, 0, ImGuiInputTextFlags_EnterReturnsTrue)) { //search on enter
    if (criticalpointID > 0 && criticalpointID < loadedCPAmount) {
      selectCriticalPoint(criticalpointID);
    }
  }

  //checkbox to flash atoms
  if (ImGui::Checkbox("Selc.findCP", &flashCP)) { //set flashing to defults
    framesMaxCP = 15; // ~0.5 seconds
    framesLeftCP = 0;
    otherCriticalPointsVisable = true;
  }

  ImGui::End();
}

static void showMenuFunctions(){
  if (ImGui::MenuItem("Generate Critical Points")) {
    auto_cp();
  }
  if (ImGui::MenuItem("Load Critical Points")) {
    destructCriticalPoints();
    loadCriticalPoints();
  }
}

static void showMenuVisuals(bool * show_bonds, bool * show_cps, bool * show_atoms) {
  if (ImGui::MenuItem("show/hide Bonds")) {
    *show_bonds = !*show_bonds;
  }
  if (ImGui::MenuItem("show/hide Crit Pts")) {
    *show_cps = !*show_cps;
  }
  if (ImGui::MenuItem("show/hide Atoms")) {
    *show_atoms = !*show_atoms;
  }

  if (ImGui::MenuItem("flash current selection")) {
    flashAtoms = !flashAtoms;
    framesMax = 15; // ~0.5 seconds
    framesLeft = 0;
    otherAtomsVisable = true;
  }
}

static void ShowAppMainMenuBar(bool * show_bonds, bool * show_cps, bool * show_atoms)
{
  if (ImGui::BeginMainMenuBar())
    {
      if (ImGui::BeginMenu("File")){
	ShowMenuFile();
	ImGui::EndMenu();
      }
      if (ImGui::BeginMenu("Functions")) {
	showMenuFunctions();
	ImGui::EndMenu();
      }
      if (ImGui::BeginMenu("Visuals")) {
	showMenuVisuals(show_bonds, show_cps, show_atoms);
	ImGui::EndMenu();
      }

      ImGui::EndMainMenuBar();
    }
}

static void ShowMenuFile()
{
  if (ImGui::MenuItem("Load Molecule")) {
    char const * lTheOpenFileName = tinyfd_openFileDialog(
							  "Select Molecule file",
							  "../../examples/data/benzene.wfx",
							  0,
							  NULL,
							  NULL,
							  0);

    if (lTheOpenFileName == NULL) {
      return;
    }

    init_struct();
    call_structure(lTheOpenFileName, (int)strlen(lTheOpenFileName), 1);
    destructLoadedMolecule();
    destructCriticalPoints();
    loadAtoms();
    loadBonds();
  }
  if (ImGui::MenuItem("Load Crystal")) {
    char const * lTheOpenFileName = tinyfd_openFileDialog(
							  "Select Molecule file",
							  "../../examples/data/ammonia.big.vel.cube",
							  0,
							  NULL,
							  NULL,
							  0);

    if (lTheOpenFileName == NULL) {
      return;
    }

    init_struct();
    call_structure(lTheOpenFileName, (int)strlen(lTheOpenFileName), 0);
    destructLoadedMolecule();
    destructCriticalPoints();
    loadAtoms();
    loadBonds();
  }
  if (ImGui::MenuItem("Clear")) {
    destructLoadedMolecule();
    destructCriticalPoints();
  }
}

int main(int argc, char *argv[])
{
  // Initialize the critic2 library
  critic2_initialize();

  // Setup window
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    return 1;
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_DECORATED, GL_TRUE);
  glfwWindowHint(GLFW_VISIBLE, GL_TRUE);

  // Create the window
  GLFWwindow* window = glfwCreateWindow(640, 480, "gcritic2", NULL, NULL);
  glfwMakeContextCurrent(window);
  gl3wInit();

  // Setup ImGui binding
  ImGui_ImplGlfwGL3_Init(window, true);

  // Event callbacks
  glfwSetScrollCallback(window, ScrollCallback);

  //Setup up OpenGL stuff
  GLuint VertexArray;
  glGenVertexArrays(1, &VertexArray);
  glBindVertexArray(VertexArray);

  //GLuint trishader = CompileShaders();
  //gWorldLocation = glGetUniformLocation(trishader, "gWorld");
  //mColorLocation = glGetUniformLocation(trishader, "mColor");
  GLuint lightshader = LightingShader();
  ShaderVarLocations.gWorldLocation = glGetUniformLocation(lightshader, "gWorld");
  ShaderVarLocations.gWVPLocation = glGetUniformLocation(lightshader, "gWVP");
  ShaderVarLocations.vColorLocation = glGetUniformLocation(lightshader, "vColor");
  ShaderVarLocations.lColorLocation = glGetUniformLocation(lightshader, "lColor");
  ShaderVarLocations.lDirectionLocation = glGetUniformLocation(lightshader, "lDirection");
 
  //glEnables
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);

  // initialize pipeline
  Pipeline p;

  // Load sphere mesh
  GLuint SphereIB, SphereVB;
  CreateAndFillBuffers(&SphereVB, &SphereIB, isphv[atomresolution], isphi[atomresolution], 
		       3*nsphv[atomresolution], 3*nsphi[atomresolution]);
 
  // Load cylinder mesh
  GLuint CylIB, CylVB;
  CreateAndFillBuffers(&CylVB, &CylIB, icycv[bondresolution], icyci[bondresolution], 
		       3*ncycv[bondresolution], 3*ncyci[bondresolution]);
 
 
  // Imgui static variables
  static bool show_bonds = true;
  static bool show_cps = true;
  static bool show_atoms = true;
 
  // input variables;
  // c means for current loop, l means last loop, p means last pressed
  static int cLMB;
  static int cRMB;
  static int lLMB;
  static int lRMB;
  static double cMPosX;
  static double cMPosY;
  static double lMPosX;
  static double lMPosY;
  static double pMPosX;
  static double pMPosY;
  static double scrollY;
 
  cam.Pos[0] = 0.f; cam.Pos[1] = 0.f; cam.Pos[2] = -10.f;
  cam.Target[0] = 0.f; cam.Target[1] = 0.f; cam.Target[2] = 1.f;
  cam.Up[0] = 0.f; cam.Up[1] = 1.f; cam.Up[2] = 0.f;
 
  time_t lastTime = time(0);
  time_t curTime = lastTime;
  double frameTime = 35.0;

  Vector3f curRotAxis = Vector3f(0, 0, 0);
  Vector3f lastRotAxis = Vector3f(0, 0, 0);
  Vector3f rotAxis = Vector3f(0, 0, 0);
 
  static float lastRotAng = 0;
  static float curRotAng = 0;
  static float rotAng = 0;
 
  static float diffX;
  static float diffY;
 
  Matrix4f lastRot;
  Matrix4f curRot;
  Matrix4f rot;
  lastRot.InitIdentity();
  curRot.InitIdentity();
  rot.InitIdentity();
  bool show_test_window = true;
  //
  // Main loop ------------------------------------------------------------------
  //
  while (!glfwWindowShouldClose(window)){
    curTime = time(0);
    if ((difftime(lastTime, curTime) < frameTime)) {
#ifdef WIN32
      Sleep(frameTime - difftime(lastTime, curTime));
#endif // WIN32
      //#if defined(LINUX) || defined(__APPLE__)
      usleep(frameTime - difftime(lastTime, curTime));
      //#endif // LINUX || __APPLE__
 
 
    }
    lastTime = curTime;
 
 
    glfwPollEvents();
    ImGui_ImplGlfwGL3_NewFrame();
 
    // Process input
    ImGuiIO& io = ImGui::GetIO();
    lLMB = cLMB;
    lRMB = cRMB;
    cLMB = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    cRMB = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
    lMPosX = cMPosX;
    lMPosY = cMPosY;
    glfwGetCursorPos(window, &cMPosX, &cMPosY);
 
    float camPanFactor = fabs(0.00115f * cam.Pos[2]);
    float camRotateFactor = 0.015f;
    if (!io.WantCaptureMouse) {
      if (cLMB == GLFW_PRESS){
	cam.Pos[0] -= camPanFactor * (cMPosX - lMPosX);
	cam.Pos[1] += camPanFactor * (cMPosY - lMPosY);
      }
      if (cRMB == GLFW_PRESS){
	if (lRMB != GLFW_PRESS){
	  pMPosX = cMPosX;
	  pMPosY = cMPosY;
 
	  lastRot = rot;
	} else {
 
	  diffX = (float)(cMPosX - pMPosX);
	  diffY = (float)(cMPosY - pMPosY);
 
	  curRotAxis = Vector3f(diffX, -diffY, 0);
	  curRotAxis = curRotAxis.Cross(Vector3f(0, 0, 1));
	  curRotAng = curRotAxis.Length() * camRotateFactor;
	  curRotAxis.Normalize();
 
	  curRot.InitRotateAxisTransform(curRotAxis, curRotAng);
	  rot = curRot * lastRot;
	}
      }
    }
 
    // Rendering
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
 
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.1f, 0.1f, 0.4f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 
    p.SetPersProjInfo(45, display_w, display_h, 1.f, 1000.f);
    p.SetOrthoProjInfo(-10.f, 10.f, -10.f, 10.f, -1000.f, 1000.f);
    p.SetPostRotationMatrix(rot);
    p.SetCamera(cam);
 
    glEnableVertexAttribArray(0);
 
    // molecule drawing
    if (show_bonds){
      drawAllBonds(&p, CylVB, CylIB);
    }
    if (show_atoms){
      drawAllAtoms(&p, SphereVB, SphereIB);
    }
    if (show_cps){
      drawAllCPs(&p, SphereVB, SphereIB);
    }
    //color legend
    // imgui overlays
    //    		printCamStats();
    //        ShowAppMainMenuBar();
 
    drawMainMenuTree(display_w, display_h-5);
    //drawToolBar(display_w, display_h, &show_bonds, &show_cps, &show_atoms);
    createAtomSearchBar();
    createCriticalPointSearchBar();
    ShowAppMainMenuBar(&show_bonds, &show_cps, &show_atoms);
 
    glDisableVertexAttribArray(0);
    glUseProgram(lightshader);
 
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    ImGui::Render();
    glfwSwapBuffers(window);
  }
 
  // Cleanup on program end
  ImGui_ImplGlfwGL3_Shutdown();
  glfwTerminate();

  // Terminate the critic2 run
  critic2_end();

  return 0;
}

