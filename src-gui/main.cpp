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

using namespace std;
using namespace ImGui;

static void error_callback(int error, const char* description){
  fprintf(stderr, "Error %d: %s\n", error, description);
}

int main(int argc, char *argv[]){
  // Initialize
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()) exit(EXIT_FAILURE);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

  // Set up window
  GLFWwindow* rootwin = glfwCreateWindow(1280, 720, "gcritic2", nullptr, nullptr);
  assert(rootwin!=nullptr);
  glfwMakeContextCurrent(rootwin);
  gl3wInit();

  // Initialize the critic2 library
  c2::gui_initialize((void *) rootwin);

  // Setup ImGui binding
  ImGui_ImplGlfwGL3_Init(rootwin, true);

  // GUI settings
  ImGuiIO& io = GetIO();
  io.IniFilename = nullptr;

  // Main loop
  while (!glfwWindowShouldClose(rootwin)){
    // New frame
    glfwPollEvents();
    ImGui_ImplGlfwGL3_NewFrame();
    ImGuiContext *g = GetCurrentContext();

    // Draw a basic example of the GUI
    // Main menu bar
    if (BeginMainMenuBar()){
      if (BeginMenu("File")){
        if (MenuItem("Quit","Ctrl+Q"))
          glfwSetWindowShouldClose(rootwin, GLFW_TRUE);
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

    // Containers
    Dock *dviewcont = nullptr;
    static bool show_viewcont = true;
    if (show_viewcont)
      dviewcont = ImGui::Container("Container##viewcontainer",&show_viewcont,0,Dock::DockFlags_NoLiftContainer | Dock::DockFlags_Transparent);

    // Docks
    static bool show_treedock = true;
    static bool show_infodock = true;
    Dock *dtreedock = nullptr, *dinfodock = nullptr, *dviewdock = nullptr;
    if (show_treedock){
      if (BeginDock("Tree view",&show_treedock)){
        Text("Tree View!");     
        if (Button("Tree view")){printf("Tree view\n");}
      }
      dtreedock = GetCurrentDock();
      EndDock();
    }
    if (show_infodock){
      if (BeginDock("Info view",&show_infodock)){
        Text("Info View!");     
        if (Button("Info view")){printf("Info view\n");}
      }
      dinfodock = GetCurrentDock();
      EndDock();
    }
    if (BeginDock("Main view",nullptr,0,Dock::DockFlags_NoLiftContainer,dviewcont)){
      c2::draw_scene();
    }
    dviewdock = GetCurrentDock();
    EndDock();

    // Dock everything in the first pass
    static bool first = true;
    if (first){
      first = false;
      droot->newDockRoot(dviewcont,5);

      Dock *dtmp = dviewcont->newDockRoot(dtreedock,4);
      dviewcont->setSlidingBarPosition(4,0.2f);

      dtmp->newDockRoot(dinfodock,3);
      dtmp->setSlidingBarPosition(3,0.7f);


      dviewdock->setDetachedDockSize(300.f,300.f);
      dtreedock->setDetachedDockSize(300.f,300.f);
      dinfodock->setDetachedDockSize(300.f,300.f);
    }

    // Draw the current scene

    // Render and swap
    Render();
    glfwSwapBuffers(rootwin);
  }

  // Cleanup
  ShutdownDock();
  ImGui_ImplGlfwGL3_Shutdown();
  glfwTerminate();

  return 0;
}

