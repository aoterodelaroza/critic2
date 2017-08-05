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

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_dock.h"
#include "imgui/imgui_widgets.h"

#include <GL/glu.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>
#include "critic2.h"

using namespace std;

static void error_callback(int error, const char* description){
  fprintf(stderr, "Error %d: %s\n", error, description);
}

int main(int argc, char *argv[]){
  // Initialize
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()) exit(EXIT_FAILURE);

  // Set up window
  GLFWwindow* window = glfwCreateWindow(1280, 720, "gcritic2", nullptr, nullptr);
  assert(window!=nullptr);
  glfwMakeContextCurrent(window);

  // Initialize the critic2 library
  c2::gui_initialize((void *) window);

  // Setup ImGui binding
  ImGui_ImplGlfw_Init(window, true);

  // GUI settings
  ImGuiIO& io = ImGui::GetIO();
  io.IniFilename = nullptr;

  // Main loop
  while (!glfwWindowShouldClose(window)){
    // New frame
    glfwPollEvents();
    ImGui_ImplGlfw_NewFrame();

    // Draw the GUI
    if (ImGui::BeginMainMenuBar()){
      if (ImGui::BeginMenu("File")){
        if (ImGui::MenuItem("Quit","Ctrl+Q"))
          glfwSetWindowShouldClose(window, GLFW_TRUE);
        ImGui::EndMenu();
      }
      ImGui::SameLine(0, ImGui::GetWindowSize().x-250.);
      ImGui::Text("%.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
      ImGui::EndMainMenuBar();
    }

    static bool show_scene1 = true;
    static bool show_scene2 = true;
    static bool show_scene3 = true;
    static bool show_scene4 = true;
    static bool show_scene5 = true;
    static bool show_scene6 = true;
    static bool show_scene7 = true;
    static bool show_scene8 = true;
    static bool show_scene9 = true;


    // // for testing the docking routines for save/load layout // //
    static bool first = true;
    ImGui::Dock *d1, *d2, *d3, *d4, *d5, *d6, *dcont1, *dcont2, *droot;

    if (show_scene9){
      ImGui::SetNextWindowPos(ImVec2(450,550),ImGuiSetCond_Once);
      ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_Once);
      droot = ImGui::RootContainer("rootcontain",&show_scene9);
    }

    if (show_scene7){
      ImGui::SetNextWindowPos(ImVec2(200,200),ImGuiSetCond_Once);
      ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
      dcont1 = ImGui::Container("contain7",&show_scene7);
    }
    if (show_scene8)
      dcont2 = ImGui::Container("contain8",&show_scene8);

    ImGui::Dock *dcont3 = ImGui::Container("contain9");
    ImGui::Dock *dcont4 = ImGui::Container("contain10");


    if (show_scene1){
      ImGui::SetNextWindowPos(ImVec2(450,550),ImGuiSetCond_Once);
      ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_Once);
      if (ImGui::BeginDock("Style Editor", &show_scene1)){
        ImGui::ShowStyleEditor();
      }
      d1 = ImGui::GetCurrentDock();
      ImGui::EndDock();
    }
    if (show_scene2){
      ImGui::SetNextWindowPos(ImVec2(450,550),ImGuiSetCond_Once);
      ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_Once);
      if (ImGui::BeginDock("Info2", &show_scene2)){
        ImGui::Text("Hello2!");
        if (ImGui::Button("Button2")){printf("Button2\n");}
      }
      d2 = ImGui::GetCurrentDock();
      ImGui::EndDock();
    }
    if (show_scene3){
      ImGui::SetNextWindowPos(ImVec2(450,550),ImGuiSetCond_Once);
      ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_Once);
      if (ImGui::BeginDock("Info3", &show_scene3)){
        ImGui::Text("Hello3!");
        if (ImGui::Button("Button3")){printf("Button3\n");}
      }
      d3 = ImGui::GetCurrentDock();
      ImGui::EndDock();
    }
    if (show_scene4){
      ImGui::SetNextWindowPos(ImVec2(450,550),ImGuiSetCond_Once);
      ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_Once);
      if (ImGui::BeginDock("Info4", &show_scene4)){
        ImGui::Text("Hello4!");
        if (ImGui::Button("Button4")){printf("Button4\n");}
      }
      d4 = ImGui::GetCurrentDock();
      ImGui::EndDock();
    }
    if (show_scene5){
      ImGui::SetNextWindowPos(ImVec2(450,550),ImGuiSetCond_Once);
      ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_Once);
      if (ImGui::BeginDock("Info5", &show_scene5)){
        ImGui::Text("Hello5!");
        if (ImGui::Button("Button5")){printf("Button5\n");}
      }
      d5 = ImGui::GetCurrentDock();
      ImGui::EndDock();
    }
    if (show_scene6){
      ImGui::SetNextWindowPos(ImVec2(450,550),ImGuiSetCond_Once);
      ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_Once);
      if (ImGui::BeginDock("Info6", &show_scene6)){
        ImGui::Text("Hello6!");
        if (ImGui::Button("Button6")){printf("Button6\n");}
      }
      d6 = ImGui::GetCurrentDock();
      ImGui::EndDock();
    }
    // if (first){
    //   first = false;
    //   dcont1->newDock(d1);
    //   dcont1->newDock(d2);

    //   dcont2->newDock(d3);
    //   dcont2->newDock(d4,0);

    //   ImGui::Dock *tmp;
    //   tmp = droot->newDockRoot(dcont1,5);
    //   tmp = tmp->newDockRoot(dcont3,1);
    //   dcont1->newDockRoot(dcont2,2);
    //   tmp->newDockRoot(dcont4,2);
    // }

    if (first){
      first = false;
      ImGui::Dock *tmp;
      tmp = droot->newDockRoot(d1,5);
      tmp = tmp->newDockRoot(d2,1);
      tmp = tmp->newDockRoot(d3,1);
      tmp = tmp->newDockRoot(dcont3,3);
      dcont3->newDock(d4);
      dcont3->newDock(d5,0);
      dcont3->newDock(d6,0);
    }

    // for testing the resize grip // //

    // // for testing the resize grip // //
    // ImGuiWindow *cwindow;
    // ImGui::SetNextWindowPos(ImVec2(40.f,40.f),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200.f,200.f),ImGuiSetCond_Once);
    // if (ImGui::Begin("onewindow",nullptr)){
    //   ImGui::Text("Blah!!");
    // }
    // cwindow = ImGui::GetCurrentWindow();
    // ImGui::End();

    // ImGui::SetNextWindowPos(ImVec2(340.f,40.f),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200.f,200.f),ImGuiSetCond_Once);
    // if (ImGui::Begin("anotherwindow",nullptr,ImGuiWindowFlags_NoResize)){
    //   ImGui::PrintDock__();
    //   ImGui::ResizeGripOther("anotherwindow", ImGui::GetCurrentWindow(), cwindow);
    // }
    // ImGui::End();
    // // for testing the resize grip // //

    // // for hovering checks // //
    // ImGui::SetNextWindowPos(ImVec2(20,20),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
    // ImGui::Container("contain1");
    // ImGui::SetNextWindowPos(ImVec2(250,20),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
    // ImGui::Container("contain2");
    // ImGui::SetNextWindowPos(ImVec2(500,20),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
    // ImGui::Container("contain3");
    // ImGui::SetNextWindowPos(ImVec2(20,250),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
    // ImGui::Container("contain4");
    // ImGui::SetNextWindowPos(ImVec2(250,250),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
    // ImGui::Container("contain5");
    // ImGui::SetNextWindowPos(ImVec2(500,250),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
    // ImGui::Container("contain6");
    // ImGui::SetNextWindowPos(ImVec2(750,500),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
    // if (ImGui::BeginDock("dock1")){}
    // ImGui::EndDock();
    // ImGui::SetNextWindowPos(ImVec2(750,500),ImGuiSetCond_Once);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_Once);
    // if (ImGui::BeginDock("dock2")){}
    // ImGui::EndDock();
    // ImGui::SetNextWindowPos(ImVec2(200,500),ImGuiSetCond_FirstUseEver);
    // ImGui::SetNextWindowSize(ImVec2(200,200),ImGuiSetCond_FirstUseEver);
    // ImGui::Dock *cont3 = ImGui::RootContainer("rootcontain");
    // // for hovering checks // //

    // // // for general checks // //
    // if (show_scene7) {
    //   ImGui::SetNextWindowPos(ImVec2(100,300),ImGuiSetCond_FirstUseEver);
    //   ImGui::SetNextWindowSize(ImVec2(400,400),ImGuiSetCond_FirstUseEver);
    //   ImGui::Dock *cont3 = ImGui::RootContainer("rootcontain",&show_scene7,ImGuiWindowFlags_NoBringToFrontOnFocus);
    // }
    // ImGui::Container("contain4");
    // ImGui::Container("contain5");
    // if (show_scene8) ImGui::Container("contain6",&show_scene8);
    // if (show_scene9) ImGui::Container("contain7",&show_scene9);

    // if (show_scene1){
    //   ImGui::SetNextWindowPos(ImVec2(500,200),ImGuiSetCond_FirstUseEver);
    //   ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_FirstUseEver);
    //   if (ImGui::BeginDock("Info",&show_scene1)){
    //     ImGui::Text("Hello!");
    //     if (ImGui::Button("Button##1")){printf("Button\n");}
    //     ImGui::SameLine();
    //     if (ImGui::Button("Button##2")){printf("Button\n");}
    //     ImGui::SameLine();
    //     if (ImGui::Button("Button##3")){printf("Button\n");}
    //     if (ImGui::Button("Button##4")){printf("Button\n");}
    //     ImGui::SameLine();
    //     if (ImGui::Button("Button##5")){printf("Button\n");}
    //     ImGui::SameLine();
    //     if (ImGui::Button("Button##6")){printf("Button\n");}

    //     static char command[2048] = "";
    //     static char command2[2048] = "";
    //     ImGui::Text("Input:");
    //     ImGui::SameLine();
    //     if (ImGui::InputText("###inputconsole", command, IM_ARRAYSIZE(command), ImGuiInputTextFlags_EnterReturnsTrue|ImGuiInputTextFlags_AutoSelectAll|ImGuiInputTextFlags_AlwaysInsertMode)){
    //     }
    //     if (ImGui::InputText("###inputconsole2", command2, IM_ARRAYSIZE(command2), ImGuiInputTextFlags_EnterReturnsTrue|ImGuiInputTextFlags_AutoSelectAll|ImGuiInputTextFlags_AlwaysInsertMode)){
    //     }
    //   }
    //   ImGui::EndDock();
    // }
    // if (show_scene2){
    //   ImGui::SetNextWindowPos(ImVec2(500,200),ImGuiSetCond_FirstUseEver);
    //   ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_FirstUseEver);
    //   if (ImGui::BeginDock("Info2", &show_scene2,ImGuiWindowFlags_NoTitleBar)){
    //     ImGui::Text("Hello2!");
    //     if (ImGui::Button("Button2")){printf("Button2\n");}
    //   }
    //   ImGui::EndDock();
    // }
    // if (show_scene3){
    //   ImGui::SetNextWindowPos(ImVec2(500,200),ImGuiSetCond_FirstUseEver);
    //   ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_FirstUseEver);
    //   if (ImGui::BeginDock("NoClose1")){
    //     ImGui::Text("Hello3!");
    //     if (ImGui::Button("Button3")){printf("Button3\n");}
    //   }
    //   ImGui::EndDock();
    // }
    // if (show_scene4){
    //   ImGui::SetNextWindowPos(ImVec2(500,200),ImGuiSetCond_FirstUseEver);
    //   ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_FirstUseEver);
    //   if (ImGui::BeginDock("NoClose2")){
    // 	ImGui::Text("Hello4!");
    // 	if (ImGui::Button("Button4")){printf("Button4\n");}
    //   }
    //   ImGui::EndDock();
    // }
    // if (show_scene5){
    //   ImGui::SetNextWindowPos(ImVec2(500,200),ImGuiSetCond_FirstUseEver);
    //   ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_FirstUseEver);
    //   if (ImGui::BeginDock("Info5", &show_scene5,ImGuiWindowFlags_NoBringToFrontOnFocus)){
    // 	ImGui::Text("Hello5!");
    // 	if (ImGui::Button("Button5")){printf("Button5\n");}
    //   }
    //   ImGui::EndDock();
    // }

    // // for general checks // //

    // // docking info // //
    ImGui::SetNextWindowPos(ImVec2(900,100),ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(300,300),ImGuiSetCond_FirstUseEver);
    if (ImGui::Begin("justawindow")){ImGui::PrintDock__();}
    ImGui::EndDock();
    // // docking info // //

    // Draw the current scene
    c2::draw_scene();

    // Render and swap
    ImGui::Render();
    glfwSwapBuffers(window);
  }

  // Cleanup
  ImGui::ShutdownDock();
  ImGui_ImplGlfw_Shutdown();
  glfwTerminate();

  return 0;
}

