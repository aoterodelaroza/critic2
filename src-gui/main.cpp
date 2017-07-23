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

#include <GL/glu.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>
#include "critic2.h"

using namespace std;

static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error %d: %s\n", error, description);
}

int main(int argc, char *argv[])
{
    // Initialize
    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) exit(EXIT_FAILURE);

    // Set up window
    GLFWwindow* window = glfwCreateWindow(1280, 720, "gcritic2", NULL, NULL);
    assert(window!=NULL);
    glfwMakeContextCurrent(window);

    // Initialize the critic2 library
    c2::gui_initialize((void *) window);

    // Setup ImGui binding
    ImGui_ImplGlfw_Init(window, true);

    bool show_another_window = false;
    ImVec4 clear_color = ImColor(114, 144, 154);

    // Main loop
    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        ImGui_ImplGlfw_NewFrame();

        // 1. Show a simple window
        // Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets appears in a window automatically called "Debug"
        {
            static float f = 0.0f;
            ImGui::Text("Hello, world!");
            ImGui::SliderFloat("float", &f, 0.0f, 1.0f);
            ImGui::ColorEdit3("clear color", (float*)&clear_color);
            if (ImGui::Button("Another Window")) show_another_window ^= 1;
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        }

        // 2. Show another simple window, this time using an explicit Begin/End pair
        if (show_another_window)
        {
            ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiSetCond_FirstUseEver);
            ImGui::Begin("Another Window", &show_another_window);
            ImGui::Text("Hello");
            ImGui::End();
        }

        // Draw the current scene
	c2::draw_scene();

	// Render and swap
        ImGui::Render();
        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplGlfw_Shutdown();
    glfwTerminate();

  return 0;
}

