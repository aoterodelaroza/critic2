/*
Copyright (c) 2015 Alberto Otero de la Roza
<aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
<victor@fluor.quimica.uniovi.es>.

critic2 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

critic2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

#include <critic2.h>

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw_gl3.h"

#include "matrix_math.h"
#include "geometry.h"
#include "imguifilesystem.h"

// #ifdef WIN32 //platform spisific sleep functions
// #include <synchapi.h>
// #endif // WIN32
// #if defined(LINUX) || defined(__APPLE__)
#include <unistd.h>
// #endif // LINUX || __APPLE__

using namespace std;

// Forward declarations //
static void new_structure_dialog(bool *p_open, int ismolecule);
static void ShowAppMainMenuBar(bool* show_bonds, bool* show_cps, bool* show_atoms);

// Static GUI variables //
// Bond and atom resolutions (0 = coarse -> 3 = smooth)
static const char bondresolution = 2;
static const char atomresolution = 1;

// Bond thickness and atom/CP size
static const float bondthickness = 0.05;
static const float atomsize = 0.5;
static const float cpsize = 0.2;

// Current state of the camera
static CameraInfo cam;

// Shader and shader variables
static GLuint lightshader;
static struct {
  GLuint gWorldLocation;
  GLuint gWVPLocation;
  GLuint vColorLocation;
  GLuint lColorLocation;
  GLuint lDirectionLocation;
  GLuint fAmbientIntensityLocation;
} ShaderVarLocations;

// xxxx //

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
  float camZoomFactor = xmaxlen * 0.2f;
  cam.Pos[2] += yoffset * camZoomFactor;
}

// draw a bond between 2 atoms defined in the bond struct
void DrawBond(Pipeline * p, GLuint CylVB, GLuint CylIB, c_bond *b)
{
  float grey[3] = {.5, .5, .5};
  float white[3] = {1, 1, 1};

  p->Scale(bondthickness, bondthickness, b->length);
  p->Translate(b->r2[0]-xcm[0], b->r2[1]-xcm[1], b->r2[2]-xcm[2]);
  p->SetRotationMatrix(b->rot);

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

/// draw an atom using gl functions
void drawAtomInstance(int id, float posVector[3], float color[3],
                      Pipeline * p, GLuint SphereVB, GLuint SphereIB) {

  Vector3f center = Vector3f(xcm[0],xcm[1],xcm[2]);
  Vector3f position = Vector3f(posVector[0],posVector[1],posVector[2]);

  float rscal = atomsize * at[id].rad;
  p->Scale(rscal,rscal,rscal);

  Vector3f pos = position - center;
  p->Translate(pos.x, pos.y, pos.z);
  p->Rotate(0.f, 0.f, 0.f);

  glUniformMatrix4fv(ShaderVarLocations.gWVPLocation, 1, GL_TRUE,
                     (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(ShaderVarLocations.gWorldLocation, 1, GL_TRUE,
                     (const GLfloat *)p->GetWorldTrans());
  glBindBuffer(GL_ARRAY_BUFFER, SphereVB);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, SphereIB);
  glUniform4fv(ShaderVarLocations.vColorLocation, 1, color);
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
void drawCritPointInstance(int identifier, float posVector[3], float color[3],
			   Pipeline * p, GLuint SphereVB, GLuint SphereIB) {
  p->Scale(cpsize, cpsize, cpsize);
  p->Translate(posVector[0]-xcm[0],posVector[1]-xcm[1],posVector[2]-xcm[2]);
  p->Rotate(0.f, 0.f, 0.f); 
  glUniformMatrix4fv(ShaderVarLocations.gWVPLocation, 1, GL_TRUE, (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(ShaderVarLocations.gWorldLocation, 1, GL_TRUE, (const GLfloat *)p->GetWorldTrans());
  glBindBuffer(GL_ARRAY_BUFFER, SphereVB);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, SphereIB);
  glUniform4fv(ShaderVarLocations.vColorLocation, 1, color);
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
bool otherAtomsVisible = true;


bool flashCP = false;
int framesMaxCP = 15;
int framesLeftCP = 0;
bool otherCriticalPointsVisible = true;


/// moves cam over atom (alligned to z axis)
void lookAtAtom(int atomNumber) {
  cam.Pos[0] = at[atomNumber].r[0];
  cam.Pos[1] = at[atomNumber].r[1];
}

/// moves cam over crit point (alligned to z axis)
void lookAtCritPoint(int critPointNum) {
  cam.Pos[0] = critp[critPointNum].r[0];
  cam.Pos[1] = critp[critPointNum].r[1];
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
  // displayVars[1] = to_string(at[atomNumber].numberOfBonds);
  displayVars[1] = "";
  displayVars[2] = "";
}

void atomAtomicNumberInfo(string * displayVars, int atomNumber) {
  displayVars[0] = "atomic number";
  displayVars[1] = to_string(at[atomNumber].z);
  displayVars[2] = "";
}

void criticalPointTypeInfo(string * displayVars, int criticalPointIndex) {
  displayVars[0] = "critical point Type";
  // displayVars[1] = loadedCriticalPoints[criticalPointIndex].typeName;
  displayVars[2] = "";
}

static void showMenuFunctions(){
  if (ImGui::MenuItem("Generate Critical Points")) {
    call_auto();
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
    otherAtomsVisible = true;
  }
}

static void ShowAppMainMenuBar(bool * show_bonds, bool * show_cps, bool * show_atoms, bool * want_quit)
{
  // logical variables for persistent dialogs
  static bool show_new_structure_dialog = false;
  static int ismolecule;

  // persistent dialog calls
  if (show_new_structure_dialog) new_structure_dialog(&show_new_structure_dialog,ismolecule);

  // immediate actions
  if (ImGui::BeginMainMenuBar())
    {
      if (ImGui::BeginMenu("File")){
	if (ImGui::MenuItem("New","Ctrl+N")){}
	if (ImGui::MenuItem("Open crystal","Ctrl+O")){
	  show_new_structure_dialog = true;
	  ismolecule = 0;
	}
	if (ImGui::MenuItem("Open molecule","Ctrl+Alt+O")) {
	  show_new_structure_dialog = true;
	  ismolecule = 1;
	}
	if (ImGui::MenuItem("Open from library","Ctrl+L")) {}
	if (ImGui::MenuItem("Open recent")) {}
	ImGui::Separator();
	if (ImGui::MenuItem("Close","Ctrl+W")) {
	}
	if (ImGui::MenuItem("Quit","Ctrl+Q")) 
	  *want_quit = true;
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

  // Shader
  lightshader = LightingShader();
  ShaderVarLocations.gWorldLocation = glGetUniformLocation(lightshader, "gWorld");
  ShaderVarLocations.gWVPLocation = glGetUniformLocation(lightshader, "gWVP");
  ShaderVarLocations.vColorLocation = glGetUniformLocation(lightshader, "vColor");
  ShaderVarLocations.lColorLocation = glGetUniformLocation(lightshader, "lColor");
  ShaderVarLocations.lDirectionLocation = glGetUniformLocation(lightshader, "lDirection");
  ShaderVarLocations.fAmbientIntensityLocation = glGetUniformLocation(lightshader, "fAmbientIntensity");
 
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
  static bool want_quit = false;
 
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
      for (int i=0; i<nbond; i++){
	DrawBond(&p, CylVB, CylIB, &bond[i]);
      }
    }
    if (show_atoms){
      for (size_t x = 0; x < nat; x++){
	drawAtomInstance(x, at[x].r, at[x].rgb, &p, SphereVB, SphereIB);
      }
    }
    if (show_cps){
      for (int x = 0; x < ncritp; x++) {
	drawCritPointInstance(x, critp[x].r, critp[x].rgb, &p, SphereVB, SphereIB);
      }
    }
 
    ShowAppMainMenuBar(&show_bonds, &show_cps, &show_atoms, &want_quit);
    if (want_quit)
      glfwSetWindowShouldClose(window, GLFW_TRUE);
 
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

// 
static void new_structure_dialog(bool *p_open, int ismolecule){
  static ImGuiFs::Dialog fsopenfile;
  static bool firstpass = true;

  const char* filename = fsopenfile.chooseFileDialog(firstpass,"./",NULL,"bleh");
  firstpass = false;

  if (fsopenfile.hasUserJustCancelledDialog() || strlen(filename) > 0){
    // Dialog has been closed - set up for next time and prevent more calls for now
    firstpass = true;
    *p_open = false;
  }

  if (strlen(filename) > 0){
    // Clean up previous and initialize the structure
    call_structure(filename, (int)strlen(filename), ismolecule); 
    firstpass = true;
    *p_open = false;
  }
}

