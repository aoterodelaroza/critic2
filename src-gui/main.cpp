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
#include <critic2.h>

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw_gl3.h"

#include "geometry.h"
#include "guiapps.h"
#include "shader.h"
#include "menu.h"
#include "draw.h"
#include "global.h"

// #ifdef WIN32 //platform spisific sleep functions
// #include <synchapi.h>
// #endif // WIN32
// #if defined(LINUX) || defined(__APPLE__)
#include <unistd.h>
// #endif // LINUX || __APPLE__

using namespace std;

// Global definitions
bool show_bonds = true;
bool show_cps = true;
bool show_atoms = true;
bool show_cell = true;

// Quit flag
bool want_quit = false;

// 
int main(int argc, char *argv[])
{
  // Initialize the critic2 library
  critic2_initialize();

  // Create the window and connect callbacks; initialize glfw/gl3w
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
  draw_set_camera_pos(-1.);  

  // Shader
  GLuint lightshader = LightingShader();
 
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
    draw_set_camera_pos(box_xmaxlen);
    show_cell = !ismolecule;
    structureinfo_window_h = true;
  }

  // 
  time_t lastTime = time(0);
  time_t curTime = lastTime;
  double frameTime = 35.0;

  // main loop
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
 
    static Matrix4f rot;
    process_mouse_input(window, &rot);

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
 
    // draw the scene elements
    draw_all_elements(&p,lightshader,window);

    // process GUI elment handles
    guiapps_process_handles();

    // menus
    show_menu_bar(&want_quit);

    // process key bindings
    ImGuiIO& io = ImGui::GetIO();
    if (io.KeyCtrl && io.KeysDown[GLFW_KEY_Q])
      want_quit = true;
    if (io.KeyCtrl && io.KeysDown[GLFW_KEY_W])
      clear_scene(true);
    if (io.KeyCtrl && io.KeysDown[GLFW_KEY_O] && !structurenew_window_h)
      structurenew_window_h = 2;
    if (io.KeyCtrl && io.KeyAlt && io.KeysDown[GLFW_KEY_O] && !structurenew_window_h)
      structurenew_window_h = 1;

    // handle quit signal
    if (want_quit)
      glfwSetWindowShouldClose(window, GLFW_TRUE);
 
    // Dummy window -> arbitrary objects rendered on the screen
    // ImGui::Begin("",NULL,ImVec2(200.,200.),0.0,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoInputs);
    // ImDrawList* drawlist = ImGui::GetWindowDrawList();
    // ImFont *font = ImGui::GetFont();
    // drawlist->AddText(ImVec2(100.,100.),ImGui::GetColorU32(ImGuiCol_Text),"bleh!");
    // ImGui::End();

    // render
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

