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

#include "critic2.h"
#include "shapes.h"
#include "view.h"
#include "dialog.h"
#include "message.h"
#include "keybinding.h"
#include "menu.h"
#include "tree.h"
#include "text.h"

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

using namespace ImGui;

int main(int argc, char *argv[]){
  // Initialize critic2
  c2::gui_initialize();

  // Parse the command line
  if (argc > 1){
    int mol = -1;
    for(int i=1; i<argc; i++){
      if (!strcmp(argv[i],"-c")){
	mol = 0;
      } else if (!strcmp(argv[i],"-m")) {
	mol = 1;
      } else if (!strcmp(argv[i],"-h")) {
	printf(" gcritic2 [-c|-m|] file1 [-c|-m|] file2...\n");
	printf("-c file : read file as a crystal\n");
	printf("-m file : read file as a molecule\n");
	printf("   file : let critic2 decide\n");
	printf("I need to write this xxxx\n");
	exit(0);
      } else {
	if (!c2::open_file((const char **) &(argv[i]), mol)){
	  std::string message = "Could not open file: " + string(argv[i]) +"\nReason: " + c2::errmsg;
	  NewMessage(Message_Error,message.c_str());
	}
	mol = -1;
      }
    }
    if (c2::nsc > 0)
      c2::scene_initialize(1);
    UpdateTreeData();
  }

  // Initialize glfw
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()) exit(EXIT_FAILURE);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_SAMPLES, 4);

  // Set up window
  GLFWwindow* rootwin = glfwCreateWindow(850, 720, "gcritic2", nullptr, nullptr);
  assert(rootwin!=nullptr);
  glfwMakeContextCurrent(rootwin);
  gl3wInit();
  glfwSetInputMode(rootwin, GLFW_STICKY_KEYS, 1);

  // Setup ImGui binding
  ImGui_ImplGlfwGL3_Init(rootwin, true);

  // ImGui settings, default settings and keybindings
  ImGuiIO& io = GetIO();
  io.IniFilename = nullptr;
  // io.MouseDrawCursor = true; // can't get anything other than the arrow otherwise
  DefaultSettings();
  SetDefaultKeyBindings();

  // Initialize the freetype library
  InitFreetype();

  // Find and read configuration file
  if (!FindConfigurationFile())
    NewMessage(Message_Error,"Could not open configuration file. Using defaults.");
  else if (!ReadConfigurationFile(conffile))
    NewMessage(Message_Error,"Could not read configuration file. Using defaults.");

  // Opengl settings
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  // Create and fill vertex, element, and frame buffers (shapes.h)
  CreateAndFillBuffers();

  // Create the main view 
  mainview = CreateView("Main view",1);
  OpenDialog(DLG_Tree);
  dlglastopen = DLG_LAST;

  // Main loop
  while (!glfwWindowShouldClose(rootwin)){
    static bool first = true;

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

    // Root container
    PushStyleVar(ImGuiStyleVar_WindowRounding, 0.f);
    PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.f);
    float menubarh = g->FontBaseSize + g->Style.FramePadding.y * 2.0f;
    SetNextWindowPos(ImVec2(0.,menubarh));
    SetNextWindowSize(ImVec2(io.DisplaySize.x,io.DisplaySize.y-menubarh));
    Dock *droot = RootContainer("critic2root",nullptr,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize|
                                ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoSavedSettings|
                                ImGuiWindowFlags_NoBringToFrontOnFocus);
    PopStyleVar(2);

    // Main menu bar
    ShowMenu(rootwin);

    // Menu dialog dispatch (before views so it picks up the popup
    // when esc is pressed - fix this later).
    DialogDispatch();

    // Draw all views
    DrawAllViews();

    // Message dispatch
    MessageDispatch();

    // Dock everything in the first pass and set default sizes
    if (first){
      first = false;

      // default size and position of the detached main view
      mainview->dock->setDetachedDockPosition(0.25f*io.DisplaySize.x,0.25f*io.DisplaySize.y);
      mainview->dock->setDetachedDockSize(0.5f*io.DisplaySize.x,0.5f*io.DisplaySize.y);
      mainview->dock->setSplitWeight(4.f,4.f);

      // attach the main view and the tree to the root
      Dock *dcont= droot->newDockRoot(mainview->dock,Dock::Drop_Tab);
      dcont->newDockRoot(dlgdock[DLG_Tree],Dock::Drop_Left);
    }

    // xxxx for imgui testing xxxx //
    // ShowTestWindow();

    // xxxx for debugging the dock system xxxx //
    // PrintDock__();

    // // xxxx for testing //
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
  delete mainview;

  ImGui_ImplGlfwGL3_Shutdown();
  glfwTerminate();

  return 0;
}

