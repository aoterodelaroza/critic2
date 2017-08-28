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
#include "camera.h"
#include "shapes.h"

using namespace std;
using namespace ImGui;

// xxxx //
const ImVec2 texsize = ImVec2(1024.f,1024.f);
// xxxx //

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

  // Initialize critic2
  c2::gui_initialize((void *) rootwin);

  // Setup ImGui binding
  ImGui_ImplGlfwGL3_Init(rootwin, true);

  // GUI settings
  ImGuiIO& io = GetIO();
  io.IniFilename = nullptr;

  // Shader and opengl settings
  Shader shader = {};
  shader.use();
  glEnable(GL_DEPTH_TEST); 

  // Vertex and element buffers (shapes.h)
  CreateAndFillBuffers();

  // Concatenate the input arguments and pass them to critic2
  if (argc > 1){
    string argall = "";
    for(int i=1;i<argc;i++)
      argall = argall + argv[i] + " ";
    c2::open_file((const char **) &argall, -1);
  }

  // xxxx //
  float vertices[] = {
    -1.0f, -1.0f, 0.0f,
    1.0f, -1.0f, 0.0f, 
    0.0f,  1.0f, 0.0f  
  }; 

  unsigned int VBO, VAO;
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);
  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 
  glBindVertexArray(0); 

  unsigned int fbo;
  glGenFramebuffers(1, &fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  unsigned int texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texsize.x, texsize.y, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);  
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0); 
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    exit(EXIT_FAILURE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  // xxxx //

  // Main loop
  while (!glfwWindowShouldClose(rootwin)){
    // New frame
    glfwPollEvents();
    ImGui_ImplGlfwGL3_NewFrame();
    ImGuiContext *g = GetCurrentContext();

    // Background
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

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
      dviewcont = Container("Container##viewcontainer",&show_viewcont,0,Dock::DockFlags_NoLiftContainer | Dock::DockFlags_Transparent);

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

    PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0,g->Style.WindowPadding.y));
    if (BeginDock("Main view",nullptr,0,Dock::DockFlags_NoLiftContainer,dviewcont)){
      glBindFramebuffer(GL_FRAMEBUFFER, fbo);
      glViewport(0.,0.,texsize.x,texsize.y);
      c2::draw_scene(1);
      ImGuiWindow *win = GetCurrentWindow(); // win->Pos win->Pos+Size
      GetWindowDrawList()->AddImage((void *) texture,win->Pos + ImVec2(0.f,win->TitleBarHeight()),
				    win->Pos + win->Size - ImVec2(0.f,win->TitleBarHeight()), 
				    ImVec2(0, 1), ImVec2(1, 0));
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
    dviewdock = GetCurrentDock();
    EndDock();
    PopStyleVar();

    // Begin("Style Editor");
    // ShowStyleEditor(); 
    // End();

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

  DeleteBuffers();

  // xxxx //
  glDeleteFramebuffers(1, &fbo); 

  // Terminate critic2
  c2::gui_end();

  // Cleanup
  ShutdownDock();
  ImGui_ImplGlfwGL3_Shutdown();
  glfwTerminate();

  return 0;
}

