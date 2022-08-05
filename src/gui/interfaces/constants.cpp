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

// Preprocessor-constants for Fortran.

#include "../imgui/imgui_impl_opengl3_loader.h"
#include "../imgui/imgui_impl_glfw.h"
#include "../imgui/imgui.h"
#include "../imgui/imgui_internal.h"
#include "../imgui/cimgui.h"
#include <GLFW/glfw3.h>
#include "constants.h"

// OpenGL
const int const_GL_COLOR_BUFFER_BIT = GL_COLOR_BUFFER_BIT;

// cimgui
//// enum ImGuiConfigFlags_
const int const_ImGuiConfigFlags_None = ImGuiConfigFlags_None;
const int const_ImGuiConfigFlags_NavEnableKeyboard = ImGuiConfigFlags_NavEnableKeyboard;
const int const_ImGuiConfigFlags_NavEnableGamepad = ImGuiConfigFlags_NavEnableGamepad;
const int const_ImGuiConfigFlags_NavEnableSetMousePos = ImGuiConfigFlags_NavEnableSetMousePos;
const int const_ImGuiConfigFlags_NavNoCaptureKeyboard = ImGuiConfigFlags_NavNoCaptureKeyboard;
const int const_ImGuiConfigFlags_NoMouse = ImGuiConfigFlags_NoMouse;
const int const_ImGuiConfigFlags_NoMouseCursorChange = ImGuiConfigFlags_NoMouseCursorChange;
const int const_ImGuiConfigFlags_DockingEnable = ImGuiConfigFlags_DockingEnable;
const int const_ImGuiConfigFlags_ViewportsEnable = ImGuiConfigFlags_ViewportsEnable;
const int const_ImGuiConfigFlags_DpiEnableScaleViewports= ImGuiConfigFlags_DpiEnableScaleViewports;
const int const_ImGuiConfigFlags_DpiEnableScaleFonts = ImGuiConfigFlags_DpiEnableScaleFonts;
const int const_ImGuiConfigFlags_IsSRGB = ImGuiConfigFlags_IsSRGB;
const int const_ImGuiConfigFlags_IsTouchScreen = ImGuiConfigFlags_IsTouchScreen;
//// enum ImGuiKey_
const int const_ImGuiKey_None = ImGuiKey_None;
const int const_ImGuiKey_Tab = ImGuiKey_Tab;
const int const_ImGuiKey_LeftArrow = ImGuiKey_LeftArrow;
const int const_ImGuiKey_RightArrow = ImGuiKey_RightArrow;
const int const_ImGuiKey_UpArrow = ImGuiKey_UpArrow;
const int const_ImGuiKey_DownArrow = ImGuiKey_DownArrow;
const int const_ImGuiKey_PageUp = ImGuiKey_PageUp;
const int const_ImGuiKey_PageDown = ImGuiKey_PageDown;
const int const_ImGuiKey_Home = ImGuiKey_Home;
const int const_ImGuiKey_End = ImGuiKey_End;
const int const_ImGuiKey_Insert = ImGuiKey_Insert;
const int const_ImGuiKey_Delete = ImGuiKey_Delete;
const int const_ImGuiKey_Backspace = ImGuiKey_Backspace;
const int const_ImGuiKey_Space = ImGuiKey_Space;
const int const_ImGuiKey_Enter = ImGuiKey_Enter;
const int const_ImGuiKey_Escape = ImGuiKey_Escape;
const int const_ImGuiKey_LeftCtrl = ImGuiKey_LeftCtrl;
const int const_ImGuiKey_LeftShift = ImGuiKey_LeftShift;
const int const_ImGuiKey_LeftAlt = ImGuiKey_LeftAlt;
const int const_ImGuiKey_LeftSuper = ImGuiKey_LeftSuper;
const int const_ImGuiKey_RightCtrl = ImGuiKey_RightCtrl;
const int const_ImGuiKey_RightShift = ImGuiKey_RightShift;
const int const_ImGuiKey_RightAlt = ImGuiKey_RightAlt;
const int const_ImGuiKey_RightSuper = ImGuiKey_RightSuper;
const int const_ImGuiKey_Menu = ImGuiKey_Menu;
const int const_ImGuiKey_0 = ImGuiKey_0;
const int const_ImGuiKey_1 = ImGuiKey_1;
const int const_ImGuiKey_2 = ImGuiKey_2;
const int const_ImGuiKey_3 = ImGuiKey_3;
const int const_ImGuiKey_4 = ImGuiKey_4;
const int const_ImGuiKey_5 = ImGuiKey_5;
const int const_ImGuiKey_6 = ImGuiKey_6;
const int const_ImGuiKey_7 = ImGuiKey_7;
const int const_ImGuiKey_8 = ImGuiKey_8;
const int const_ImGuiKey_9 = ImGuiKey_9;
const int const_ImGuiKey_A = ImGuiKey_A;
const int const_ImGuiKey_B = ImGuiKey_B;
const int const_ImGuiKey_C = ImGuiKey_C;
const int const_ImGuiKey_D = ImGuiKey_D;
const int const_ImGuiKey_E = ImGuiKey_E;
const int const_ImGuiKey_F = ImGuiKey_F;
const int const_ImGuiKey_G = ImGuiKey_G;
const int const_ImGuiKey_H = ImGuiKey_H;
const int const_ImGuiKey_I = ImGuiKey_I;
const int const_ImGuiKey_J = ImGuiKey_J;
const int const_ImGuiKey_K = ImGuiKey_K;
const int const_ImGuiKey_L = ImGuiKey_L;
const int const_ImGuiKey_M = ImGuiKey_M;
const int const_ImGuiKey_N = ImGuiKey_N;
const int const_ImGuiKey_O = ImGuiKey_O;
const int const_ImGuiKey_P = ImGuiKey_O;
const int const_ImGuiKey_Q = ImGuiKey_Q;
const int const_ImGuiKey_R = ImGuiKey_R;
const int const_ImGuiKey_S = ImGuiKey_S;
const int const_ImGuiKey_T = ImGuiKey_T;
const int const_ImGuiKey_U = ImGuiKey_U;
const int const_ImGuiKey_V = ImGuiKey_V;
const int const_ImGuiKey_W = ImGuiKey_W;
const int const_ImGuiKey_X = ImGuiKey_X;
const int const_ImGuiKey_Y = ImGuiKey_Y;
const int const_ImGuiKey_Z = ImGuiKey_Z;
const int const_ImGuiKey_F1 = ImGuiKey_F1;
const int const_ImGuiKey_F2 = ImGuiKey_F2;
const int const_ImGuiKey_F3 = ImGuiKey_F3;
const int const_ImGuiKey_F4 = ImGuiKey_F4;
const int const_ImGuiKey_F5 = ImGuiKey_F5;
const int const_ImGuiKey_F6 = ImGuiKey_F6;
const int const_ImGuiKey_F7 = ImGuiKey_F7;
const int const_ImGuiKey_F8 = ImGuiKey_F8;
const int const_ImGuiKey_F9 = ImGuiKey_F9;
const int const_ImGuiKey_F10 = ImGuiKey_F10;
const int const_ImGuiKey_F11 = ImGuiKey_F11;
const int const_ImGuiKey_F12 = ImGuiKey_F12;
const int const_ImGuiKey_Apostrophe = ImGuiKey_Apostrophe;
const int const_ImGuiKey_Comma = ImGuiKey_Comma;
const int const_ImGuiKey_Minus = ImGuiKey_Minus;
const int const_ImGuiKey_Period = ImGuiKey_Period;
const int const_ImGuiKey_Slash = ImGuiKey_Slash;
const int const_ImGuiKey_Semicolon = ImGuiKey_Semicolon;
const int const_ImGuiKey_Equal = ImGuiKey_Equal;
const int const_ImGuiKey_LeftBracket = ImGuiKey_LeftBracket;
const int const_ImGuiKey_Backslash = ImGuiKey_Backslash;
const int const_ImGuiKey_RightBracket = ImGuiKey_RightBracket;
const int const_ImGuiKey_GraveAccent = ImGuiKey_GraveAccent;
const int const_ImGuiKey_CapsLock = ImGuiKey_CapsLock;
const int const_ImGuiKey_ScrollLock = ImGuiKey_ScrollLock;
const int const_ImGuiKey_NumLock = ImGuiKey_NumLock;
const int const_ImGuiKey_PrintScreen = ImGuiKey_PrintScreen;
const int const_ImGuiKey_Pause = ImGuiKey_Pause;
const int const_ImGuiKey_Keypad0 = ImGuiKey_Keypad0;
const int const_ImGuiKey_Keypad1 = ImGuiKey_Keypad1;
const int const_ImGuiKey_Keypad2 = ImGuiKey_Keypad2;
const int const_ImGuiKey_Keypad3 = ImGuiKey_Keypad3;
const int const_ImGuiKey_Keypad4 = ImGuiKey_Keypad4;
const int const_ImGuiKey_Keypad5 = ImGuiKey_Keypad5;
const int const_ImGuiKey_Keypad6 = ImGuiKey_Keypad6;
const int const_ImGuiKey_Keypad7 = ImGuiKey_Keypad7;
const int const_ImGuiKey_Keypad8 = ImGuiKey_Keypad8;
const int const_ImGuiKey_Keypad9 = ImGuiKey_Keypad9;
const int const_ImGuiKey_KeypadDecimal = ImGuiKey_KeypadDecimal;
const int const_ImGuiKey_KeypadDivide = ImGuiKey_KeypadDivide;
const int const_ImGuiKey_KeypadMultiply = ImGuiKey_KeypadMultiply;
const int const_ImGuiKey_KeypadSubtract = ImGuiKey_KeypadSubtract;
const int const_ImGuiKey_KeypadAdd = ImGuiKey_KeypadAdd;
const int const_ImGuiKey_KeypadEnter = ImGuiKey_KeypadEnter;
const int const_ImGuiKey_KeypadEqual = ImGuiKey_KeypadEqual;
const int const_ImGuiKey_GamepadStart = ImGuiKey_GamepadStart;
const int const_ImGuiKey_GamepadBack = ImGuiKey_GamepadBack;
const int const_ImGuiKey_GamepadFaceUp = ImGuiKey_GamepadFaceUp;
const int const_ImGuiKey_GamepadFaceDown = ImGuiKey_GamepadFaceDown;
const int const_ImGuiKey_GamepadFaceLeft = ImGuiKey_GamepadFaceLeft;
const int const_ImGuiKey_GamepadFaceRight = ImGuiKey_GamepadFaceRight;
const int const_ImGuiKey_GamepadDpadUp = ImGuiKey_GamepadDpadUp;
const int const_ImGuiKey_GamepadDpadDown = ImGuiKey_GamepadDpadDown;
const int const_ImGuiKey_GamepadDpadLeft = ImGuiKey_GamepadDpadLeft;
const int const_ImGuiKey_GamepadDpadRight = ImGuiKey_GamepadDpadRight;
const int const_ImGuiKey_GamepadL1 = ImGuiKey_GamepadL1;
const int const_ImGuiKey_GamepadR1 = ImGuiKey_GamepadR1;
const int const_ImGuiKey_GamepadL2 = ImGuiKey_GamepadL2;
const int const_ImGuiKey_GamepadR2 = ImGuiKey_GamepadR2;
const int const_ImGuiKey_GamepadL3 = ImGuiKey_GamepadL3;
const int const_ImGuiKey_GamepadR3 = ImGuiKey_GamepadR3;
const int const_ImGuiKey_GamepadLStickUp = ImGuiKey_GamepadLStickUp;
const int const_ImGuiKey_GamepadLStickDown = ImGuiKey_GamepadLStickDown;
const int const_ImGuiKey_GamepadLStickLeft = ImGuiKey_GamepadLStickLeft;
const int const_ImGuiKey_GamepadLStickRight = ImGuiKey_GamepadLStickRight;
const int const_ImGuiKey_GamepadRStickUp = ImGuiKey_GamepadRStickUp;
const int const_ImGuiKey_GamepadRStickDown = ImGuiKey_GamepadRStickDown;
const int const_ImGuiKey_GamepadRStickLeft = ImGuiKey_GamepadRStickLeft;
const int const_ImGuiKey_GamepadRStickRight = ImGuiKey_GamepadRStickRight;
const int const_ImGuiKey_ModCtrl = ImGuiKey_ModCtrl;
const int const_ImGuiKey_ModShift = ImGuiKey_ModShift;
const int const_ImGuiKey_ModAlt = ImGuiKey_ModAlt;
const int const_ImGuiKey_ModSuper = ImGuiKey_ModSuper;
const int const_ImGuiKey_COUNT = ImGuiKey_COUNT;
const int const_ImGuiKey_NamedKey_BEGIN = ImGuiKey_NamedKey_BEGIN;
const int const_ImGuiKey_NamedKey_END = ImGuiKey_NamedKey_END;
const int const_ImGuiKey_NamedKey_COUNT = ImGuiKey_NamedKey_COUNT;
const int const_ImGuiKey_KeysData_SIZE = ImGuiKey_KeysData_SIZE;
const int const_ImGuiKey_KeysData_OFFSET = ImGuiKey_KeysData_OFFSET;
//// enum ImGuiNavInput_
const int const_ImGuiNavInput_Activate = ImGuiNavInput_Activate;
const int const_ImGuiNavInput_Cancel = ImGuiNavInput_Cancel;
const int const_ImGuiNavInput_Input = ImGuiNavInput_Input;
const int const_ImGuiNavInput_Menu = ImGuiNavInput_Menu;
const int const_ImGuiNavInput_DpadLeft = ImGuiNavInput_DpadLeft;
const int const_ImGuiNavInput_DpadRight = ImGuiNavInput_DpadRight;
const int const_ImGuiNavInput_DpadUp = ImGuiNavInput_DpadUp;
const int const_ImGuiNavInput_DpadDown = ImGuiNavInput_DpadDown;
const int const_ImGuiNavInput_LStickLeft = ImGuiNavInput_LStickLeft;
const int const_ImGuiNavInput_LStickRight = ImGuiNavInput_LStickRight;
const int const_ImGuiNavInput_LStickUp = ImGuiNavInput_LStickUp;
const int const_ImGuiNavInput_LStickDown = ImGuiNavInput_LStickDown;
const int const_ImGuiNavInput_FocusPrev = ImGuiNavInput_FocusPrev;
const int const_ImGuiNavInput_FocusNext = ImGuiNavInput_FocusNext;
const int const_ImGuiNavInput_TweakSlow = ImGuiNavInput_TweakSlow;
const int const_ImGuiNavInput_TweakFast = ImGuiNavInput_TweakFast;
const int const_ImGuiNavInput_KeyLeft_ = ImGuiNavInput_KeyLeft_;
const int const_ImGuiNavInput_KeyRight_ = ImGuiNavInput_KeyRight_;
const int const_ImGuiNavInput_KeyUp_ = ImGuiNavInput_KeyUp_;
const int const_ImGuiNavInput_KeyDown_ = ImGuiNavInput_KeyDown_;
const int const_ImGuiNavInput_COUNT = ImGuiNavInput_COUNT;

// GLFW
const int const_GLFW_SAMPLES = GLFW_SAMPLES;
const int const_GLFW_CONTEXT_VERSION_MAJOR = GLFW_CONTEXT_VERSION_MAJOR;
const int const_GLFW_CONTEXT_VERSION_MINOR = GLFW_CONTEXT_VERSION_MINOR;
const int const_GLFW_OPENGL_PROFILE = GLFW_OPENGL_PROFILE;
const int const_GLFW_OPENGL_FORWARD_COMPAT = GLFW_OPENGL_FORWARD_COMPAT;
const int const_GLFW_OPENGL_CORE_PROFILE = GLFW_OPENGL_CORE_PROFILE;
const int const_GLFW_STICKY_KEYS = GLFW_STICKY_KEYS;

