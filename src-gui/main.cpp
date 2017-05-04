/*
Copyright (c) 2017 Alberto Otero de la Roza
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
#include "callback.h"
#include "guiapps.h"

// #ifdef WIN32 //platform spisific sleep functions
// #include <synchapi.h>
// #endif // WIN32
// #if defined(LINUX) || defined(__APPLE__)
#include <unistd.h>
// #endif // LINUX || __APPLE__

using namespace std;

// Forward declarations //
static void ShowAppMainMenuBar();

// GUI global variables (main.h) //
// Bond and atom resolutions (0 = coarse -> 3 = smooth)
const char bondresolution = 2;
const char atomresolution = 1;

// Bond thickness and atom/CP size
const float bondthickness = 0.05;
const float atomsize = 0.5;
const float cpsize = 0.5;

// Tooltipdelay
const float ttipdelay = 1.5;

// Show/hide elements of the interface
bool show_bonds = true;
bool show_cps = true;
bool show_atoms = true;
bool show_cell = true;

// Quit flag
bool want_quit = false;

// Current state of the camera
CameraInfo cam;

// xxxx not processed yet //

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

bool IsItemHoveredDelayed(float delay,float *time0,bool *reset)
{
  float time = ImGui::GetTime();
  if (ImGui::IsItemHovered()){
    if (*time0 < 0.)
      *time0 = time;
    *reset = false;
  } 
  return (*time0 > 0.) && (time > *time0 + delay) && ImGui::IsItemHovered();
}

static void AttachTooltip(const char* desc, float delay, float *time0, bool *reset)
{
  if (IsItemHoveredDelayed(delay,time0,reset))
    {
      ImGui::BeginTooltip();
      ImGui::PushTextWrapPos(450.0f);
      ImGui::TextUnformatted(desc);
      ImGui::PopTextWrapPos();
      ImGui::EndTooltip();
    }
}

// add a shader to the gl program
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

// draw a bond between 2 atoms defined in the bond struct
void drawstick(Pipeline *p, const c_stick *s)
{
  p->Scale(s->thick, s->thick, s->length);
  p->Translate(s->r2[0], s->r2[1], s->r2[2]);
  p->SetRotationMatrix(s->rot);

  // float dir[3] = {cam.Target[0], cam.Target[1], cam.Target[2]};
  glUniformMatrix4fv(ShaderVarLocations.gWVPLocation, 1, GL_TRUE, (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(ShaderVarLocations.gWorldLocation, 1, GL_TRUE, (const GLfloat *)p->GetWorldTrans());
  glUniform4fv(ShaderVarLocations.vColorLocation, 1, (const GLfloat *)&(s->rgb));
  // glUniform4fv(ShaderVarLocations.lColorLocation, 1, (const GLfloat *)&white);
  // glUniform4fv(ShaderVarLocations.lDirectionLocation, 1, (const GLfloat *)&dir);
  // glUniform1f(ShaderVarLocations.fAmbientIntensityLocation, 0.8);

  glBindBuffer(GL_ARRAY_BUFFER, bufcylv[bondresolution]);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufcyli[bondresolution]);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glDrawElements(GL_TRIANGLES, 3*ncyli[bondresolution], GL_UNSIGNED_INT, 0);
}

// Draw a ball, with optional scaling
void drawball(Pipeline *p, const c_ball *b, float scal = 1.0)
{
  p->Scale(b->rad * scal,b->rad * scal,b->rad * scal);
  p->Translate(b->r[0], b->r[1], b->r[2]);

  glUniformMatrix4fv(ShaderVarLocations.gWVPLocation, 1, GL_TRUE, (const GLfloat *)p->GetWVPTrans());
  glUniformMatrix4fv(ShaderVarLocations.gWorldLocation, 1, GL_TRUE, (const GLfloat *)p->GetWorldTrans());
  glUniform4fv(ShaderVarLocations.vColorLocation, 1, b->rgb);

  glBindBuffer(GL_ARRAY_BUFFER, bufsphv[atomresolution]);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufsphi[atomresolution]);
  glDrawElements(GL_TRIANGLES, 3*nsphi[atomresolution], GL_UNSIGNED_INT, 0);
}

static void showMenuFunctions(){
  static float time0 = -1.;
  bool reset = true;
  if (ImGui::MenuItem("Generate Critical Points")) {
    call_auto();
  }
  AttachTooltip("Calculate the critical points.\nBleh and Blah!\n",ttipdelay,&time0,&reset);
  if (reset)
    time0 = -1.;
}

static void showMenuVisuals() {
  static float time0 = -1.;
  bool reset = true;

  if (ImGui::MenuItem("Toggle bonds","",show_bonds)) {
    show_bonds = !show_bonds;
  }
  AttachTooltip("Toggle show/hide bonds.\n",ttipdelay,&time0,&reset);

  if (ImGui::MenuItem("Toggle critical points","",show_cps)) {
    show_cps = !show_cps;
  }
  AttachTooltip("Toggle show/hide critical points.\n",ttipdelay,&time0,&reset);

  if (ImGui::MenuItem("Toggle atoms","",show_atoms)) {
    show_atoms = !show_atoms;
  }
  AttachTooltip("Toggle show/hide atoms.\n",ttipdelay,&time0,&reset);

  if (ImGui::MenuItem("Toggle cell","",show_cell)) {
    show_cell = !show_cell;
  }
  AttachTooltip("Toggle show/hide unit cell.\n",ttipdelay,&time0,&reset);

  if (ImGui::MenuItem("Show structure information","",structureinfo_window_h))
    structureinfo_window_h = !structureinfo_window_h;
  AttachTooltip("Show information about the current structure.\n",ttipdelay,&time0,&reset);

  if (reset)
    time0 = -1.;
}

static void ShowAppMainMenuBar()
{
  // immediate actions
  if (ImGui::BeginMainMenuBar())
    {
      if (ImGui::BeginMenu("File")){
	static float time0 = -1.;
	bool reset = true;

	if (ImGui::MenuItem("New","Ctrl+N",false,false)){}

	if (ImGui::MenuItem("Open crystal","Ctrl+O"))
	  structurenew_window_h = 2;
	AttachTooltip("Read the crystal structure from a file.\n",ttipdelay,&time0,&reset);

	if (ImGui::MenuItem("Open molecule","Ctrl+Alt+O"))
	  structurenew_window_h = 1;
	AttachTooltip("Read the molecular structure from a file.\n",ttipdelay,&time0,&reset);

	if (ImGui::MenuItem("Open from library","Ctrl+L",false,false)) {}

	if (ImGui::MenuItem("Open recent",NULL,false,false)) {}

	ImGui::Separator();

	if (ImGui::MenuItem("Close","Ctrl+W")) {
	  clear_scene(true);
	}
	AttachTooltip("Clear the current structure.\n",ttipdelay,&time0,&reset);

	if (ImGui::MenuItem("Quit","Ctrl+Q")) 
	  want_quit = true;
	AttachTooltip("Quit the program.\n",ttipdelay,&time0,&reset);

	ImGui::EndMenu();

	if (reset)
	  time0 = -1.;
      }
      if (ImGui::BeginMenu("Calculate")) {
	showMenuFunctions();
	ImGui::EndMenu();
      }
      if (ImGui::BeginMenu("View")) {
	showMenuVisuals();
	ImGui::EndMenu();
      }

      ImGui::EndMainMenuBar();
    }

  // Process key bindings associated to the menu
  ImGuiIO& io = ImGui::GetIO();
  if (io.KeyCtrl && io.KeysDown[GLFW_KEY_Q])
    want_quit = true;
  if (io.KeyCtrl && io.KeysDown[GLFW_KEY_W])
    clear_scene(true);
  if (io.KeyCtrl && io.KeysDown[GLFW_KEY_O])
    structurenew_window_h = 2;
  if (io.KeyCtrl && io.KeyAlt && io.KeysDown[GLFW_KEY_O])
    structurenew_window_h = 1;
}

// 
int main(int argc, char *argv[])
{
  // Initialize the critic2 library
  critic2_initialize();

  // Create the window and connect callbacks; initialize glfw/gl3w
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    return 1;
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_DECORATED, GL_TRUE);
  glfwWindowHint(GLFW_VISIBLE, GL_TRUE);
  GLFWwindow* window = glfwCreateWindow(1024, 768, "gcritic2", NULL, NULL);
  glfwMakeContextCurrent(window);
  gl3wInit();

  // Set up ImGui binding
  ImGui_ImplGlfwGL3_Init(window, true);

  // Connect scroll callback
  glfwSetScrollCallback(window, scroll_callback);

  // Some default start-up values for imgui
  ImGui::GetIO().IniFilename = NULL; // no ini file pollution
  cam.Pos[0] = 0.f; cam.Pos[1] = 0.f; cam.Pos[2] = -10.f;
  cam.Target[0] = 0.f; cam.Target[1] = 0.f; cam.Target[2] = 1.f;
  cam.Up[0] = 0.f; cam.Up[1] = 1.f; cam.Up[2] = 0.f;

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

  // Load meshes
  CreateAndFillBuffers();

  // Initialize pipeline
  Pipeline p;
 
  // Concatenate the input arguments and pass them to critic2
  if (argc > 1){
    string argall = "";
    for(int i=1;i<argc;i++)
      argall = argall + argv[i] + " ";
    call_structure((const char **) &argall, -1);
    cam.Pos[0] = 0.f; cam.Pos[1] = 0.f; cam.Pos[2] = -2.*box_xmaxlen;
    show_cell = !ismolecule;
    structureinfo_window_h = true;
  }

  // Imgui static variables
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
 
    // Process mouse input
    lLMB = cLMB;
    lRMB = cRMB;
    cLMB = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    cRMB = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
    lMPosX = cMPosX;
    lMPosY = cMPosY;
    glfwGetCursorPos(window, &cMPosX, &cMPosY);
 
    float camPanFactor = fabs(0.00115f * cam.Pos[2]);
    float camRotateFactor = 0.015f;
    if (!ImGui::GetIO().WantCaptureMouse) {
      if (cRMB == GLFW_PRESS){
	cam.Pos[0] -= camPanFactor * (cMPosX - lMPosX);
	cam.Pos[1] += camPanFactor * (cMPosY - lMPosY);
      }
      if (cLMB == GLFW_PRESS){
	if (lLMB != GLFW_PRESS){
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
 
    // Draw the scene elements
    if (show_bonds){
      for (int i=0; i<nbond; i++){
	drawstick(&p, &(bond[i].s));
      }
    }
    if (show_atoms){
      for (int i=0; i<nat; i++){
	drawball(&p, &(at[i].b), atomsize);
      }
    }
    if (show_cps){
      for (int i=0; i<ncritp; i++) {
	drawball(&p, &(critp[i].b), cpsize);
      }
    }
    if (show_cell){
      for (int i=0; i<cell_nstick; i++) {
	drawstick(&p, &(cell_s[i]));
      }
    }

    // Process GUI elment handles
    guiapps_process_handles();

    // Menus
    ShowAppMainMenuBar();
    if (want_quit)
      glfwSetWindowShouldClose(window, GLFW_TRUE);
 
    // Dummy window -> arbitrary objects rendered on the screen
    // ImGui::Begin("",NULL,ImVec2(200.,200.),0.0,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoInputs);
    // ImDrawList* drawlist = ImGui::GetWindowDrawList();
    // ImFont *font = ImGui::GetFont();
    // drawlist->AddText(ImVec2(100.,100.),ImGui::GetColorU32(ImGuiCol_Text),"bleh!");
    // ImGui::End();

    // Render
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

