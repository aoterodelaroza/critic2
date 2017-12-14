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

#include "imgui/gl3w.h"
#include <GLFW/glfw3.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw_gl3.h"
#include "imgui/imgui_dock.h"
#include "imgui/imgui_widgets.h"

#include <stdio.h>
#include <stdlib.h>
#include "critic2.h"
#include "shader.h"
#include "shapes.h"
#include "view.h"
#include "dialog.h"
#include "message.h"
#include "keybinding.h"

using namespace std;
using namespace ImGui;

int main(int argc, char *argv[]){
  // Initialize
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()) exit(EXIT_FAILURE);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

  // Set up window
  GLFWwindow* rootwin = glfwCreateWindow(850, 720, "gcritic2", nullptr, nullptr);
  assert(rootwin!=nullptr);
  glfwMakeContextCurrent(rootwin);
  gl3wInit();
  glfwSetInputMode(rootwin, GLFW_STICKY_KEYS, 1);

   // Initialize critic2
  c2::gui_initialize((void *) rootwin);

  // Setup ImGui binding
  ImGui_ImplGlfwGL3_Init(rootwin, true);

  // ImGui settings, default settings and keybindings
  ImGuiIO& io = GetIO();
  io.IniFilename = nullptr;
  DefaultSettings();
  SetDefaultKeyBindings();

  // Find and read configuration file
  if (!FindConfigurationFile())
    NewMessage(Message_Error,"Could not open configuration file. Using defaults.");
  else if (!ReadConfigurationFile(conffile))
    NewMessage(Message_Error,"Could not read configuration file. Using defaults.");

  // Shader and opengl settings
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  // Concatenate the input arguments and pass them to critic2
  if (argc > 1){
    string argall = "";
    for(int i=1;i<argc;i++)
      argall = argall + argv[i] + " ";
    c2::open_file((const char **) &argall, -1);
  }

  // Create and fill vertex, element, and frame buffers (shapes.h)
  CreateAndFillBuffers();

  // Create the main view 
  View *mainview = CreateView("Main view",1);
  OpenDialog(DLG_Tree);

  // Main loop
  while (!glfwWindowShouldClose(rootwin)){
    // New frame
    ImGui_ImplGlfwGL3_ResetKeyMouseEvents();
    glfwPollEvents();
    ImGui_ImplGlfwGL3_NewFrame();
    ImGuiContext *g = GetCurrentContext();

    // Quit key binding
    if (IsBindEvent(BIND_QUIT,false))
      glfwSetWindowShouldClose(rootwin, GLFW_TRUE);

    // Background
    ImVec4 colorbd = ImGuiStyleUI.Colors[ImGuiColUI_BackDrop];
    glClearColor(colorbd.x,colorbd.y,colorbd.z,colorbd.w);
    glClear(GL_COLOR_BUFFER_BIT);

    // Main menu bar
    if (BeginMainMenuBar()){
      if (BeginMenu("File")){
        if (MenuItem("Quit",BindKeyName(BIND_QUIT).c_str()))
	  glfwSetWindowShouldClose(rootwin, GLFW_TRUE);
        EndMenu();
      }
      if (BeginMenu("Edit")){
        if (MenuItem("Preferences..."))
	  OpenDialog(DLG_Preferences);
        EndMenu();
      }
      SameLine(0, GetWindowSize().x-250.);
      Text("%.3f ms/frame (%.1f FPS)", 1000.0f / GetIO().Framerate, GetIO().Framerate);
    }
    EndMainMenuBar();

    // Root container
    PushStyleVar(ImGuiStyleVar_WindowRounding, 0.f);
    float menubarh = g->FontBaseSize + g->Style.FramePadding.y * 2.0f;
    SetNextWindowPos(ImVec2(0.,menubarh));
    SetNextWindowSize(ImVec2(g->IO.DisplaySize.x,g->IO.DisplaySize.y-menubarh));
    Dock *droot = RootContainer("critic2root",nullptr,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize|
                                ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoSavedSettings|
                                ImGuiWindowFlags_NoBringToFrontOnFocus);
    PopStyleVar();

    // Draw all views
    DrawAllViews();

    // Menu dialog dispatch
    DialogDispatch();

    // Message dispatch
    MessageDispatch();

    // Dock everything in the first pass
    static bool first = true;
    if (first){
      first = false;
      Dock *dcont= droot->newDockRoot(mainview->dock,Dock::Drop_Tab);
      dcont->newDockRoot(dlgdock[DLG_Tree],Dock::Drop_Left);
      dcont->setSlidingBarPosition(Dock::Drop_Left,0.2f);
    }

    // xxxx for imgui testing xxxx //
    // ShowTestWindow();

    // // xxxx for debugging the dock system xxxx //
    // PrintDock__();

    // xxxx for testing //
    // RootContainer("test rootcontainer");
    // Container("test container");

    // // xxxx //
    // if (IsBindEvent(BIND_QUIT,false)){
    //   WriteLayout(droot);
    // }

    // Render and swap
    Render();
    glfwSwapBuffers(rootwin);
  }

  // Cleanup
  c2::gui_end();
  DeleteBuffers();
  ShutdownDock();
  delete shader;
  ImGui_ImplGlfwGL3_Shutdown();
  glfwTerminate();

  return 0;
}

