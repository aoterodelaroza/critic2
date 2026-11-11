// Copyright (c) 2015 Alberto Otero de la Roza
// <aoterodelaroza@gmail.com>,
// Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
// <victor@fluor.quimica.uniovi.es>.
//
// critic2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at
// your option) any later version.
//
// critic2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <http://www.gnu.org/licenses/>.

#include "glfwwrapper.h"
#include "imgui.h"

#include <GLFW/glfw3.h>
#ifdef CRITIC2_GLFW_X11
#define GLFW_EXPOSE_NATIVE_X11
#include <GLFW/glfw3native.h>
#include <X11/Xutil.h>
#endif

// constants //
const int const_GLFW_TRUE = GLFW_TRUE;
const int const_GLFW_FALSE = GLFW_FALSE;
const int const_GLFW_SAMPLES = GLFW_SAMPLES;
const int const_GLFW_CONTEXT_VERSION_MAJOR = GLFW_CONTEXT_VERSION_MAJOR;
const int const_GLFW_CONTEXT_VERSION_MINOR = GLFW_CONTEXT_VERSION_MINOR;
const int const_GLFW_OPENGL_PROFILE = GLFW_OPENGL_PROFILE;
const int const_GLFW_OPENGL_FORWARD_COMPAT = GLFW_OPENGL_FORWARD_COMPAT;
const int const_GLFW_OPENGL_CORE_PROFILE = GLFW_OPENGL_CORE_PROFILE;
const int const_GLFW_STICKY_KEYS = GLFW_STICKY_KEYS;
const int const_GLFW_VISIBLE = GLFW_VISIBLE;
const int const_GLFW_KEY_ESCAPE = GLFW_KEY_ESCAPE;
const int const_GLFW_PRESS = GLFW_PRESS;

// functions //
// Nonzero if GLFW runs on Wayland, where window positions are not
// available. glfwGetPlatform is new in GLFW 3.4; older versions
// are assumed to run on X11.
int glfw_platform_is_wayland(){
#if GLFW_VERSION_MAJOR > 3 || (GLFW_VERSION_MAJOR == 3 && GLFW_VERSION_MINOR >= 4)
  return glfwGetPlatform() == GLFW_PLATFORM_WAYLAND;
#else
  return 0;
#endif
}

// Nonzero if key is pressed in any OS window of the GUI: the main
// window and, with detachable windows, every window of an imgui
// viewport (the keyboard follows the focused one). If sticky,
// GLFW_STICKY_KEYS is set first on each of them, so a press between
// two queries is not lost (the main window has it already).
int glfw_any_window_key(int key, int sticky){
  ImGuiPlatformIO& pio = ImGui::GetPlatformIO();
  int pressed = 0;
  for (int i = 0; i < pio.Viewports.Size; i++){
    GLFWwindow* w = (GLFWwindow*) pio.Viewports[i]->PlatformHandle;
    if (!w) continue;
    if (sticky) glfwSetInputMode(w, GLFW_STICKY_KEYS, GLFW_TRUE);
    if (glfwGetKey(w, key) == GLFW_PRESS) pressed = 1;
  }
  return pressed;
}

// Nonzero if GLFW runs on X11. glfwGetPlatform is new in GLFW 3.4;
// older versions built with X11 support are assumed to run on it.
static int glfw_platform_is_x11(){
#if GLFW_VERSION_MAJOR > 3 || (GLFW_VERSION_MAJOR == 3 && GLFW_VERSION_MINOR >= 4)
  return glfwGetPlatform() == GLFW_PLATFORM_X11;
#else
  return 1;
#endif
}

// Mark the (not yet shown) window as not accepting the keyboard focus:
// on X11, the ICCCM "No Input" model (WM_HINTS input = False and no
// WM_TAKE_FOCUS, which GLFW never sets), so the window manager never
// focuses it and the focus stays on the window that had it. Used for
// the OS windows of tooltips, popups and menus, whose keyboard input
// comes through the main window. A no-op on the other platforms (the
// native X11 calls would raise a fatal GLFW error there).
void glfw_window_set_no_input_focus(GLFWwindow *window){
#ifdef CRITIC2_GLFW_X11
  if (!glfw_platform_is_x11()) return;
  Display *dpy = glfwGetX11Display();
  Window win = glfwGetX11Window(window);
  if (!dpy || !win) return;
  XWMHints *hints = XGetWMHints(dpy, win);
  if (!hints) hints = XAllocWMHints();
  if (!hints) return;
  hints->flags |= InputHint;
  hints->input = False;
  XSetWMHints(dpy, win, hints);
  XFree(hints);
#else
  (void) window;
#endif
}

