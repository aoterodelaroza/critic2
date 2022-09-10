// Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

#ifndef CIMPLOT_INCLUDED
#define CIMPLOT_INCLUDED

#include "cimgui.h"

// types
typedef int ImPlotFlags;
typedef int ImPlotAxisFlags;
typedef int ImPlotLineFlags;
typedef int ImPlotMarker;

// prototypes
CIMGUI_API ImPlotContext* ipCreateContext();
CIMGUI_API void ipDestroyContext(ImPlotContext* ctx);
CIMGUI_API bool ipBeginPlot(const char* title_id, const ImVec2 size, ImPlotFlags flags);
CIMGUI_API void ipEndPlot();
CIMGUI_API void ipSetupAxes(const char* x_label, const char* y_label, ImPlotAxisFlags x_flags, ImPlotAxisFlags y_flags);
CIMGUI_API void ipPlotLine(const char* label_id, const double* xs, const double* ys, int count, ImPlotLineFlags flags, int offset);
CIMGUI_API void SetNextMarkerStyle(ImPlotMarker marker, float size, const ImVec4 fill, float weight, const ImVec4 outline);
CIMGUI_API void ipShowDemoWindow(bool* p_open);

#endif //CIMPLOT_INCLUDED
