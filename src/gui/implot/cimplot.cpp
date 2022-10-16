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

#include "../imgui/imgui.h"
#include "../imgui/imgui_internal.h"
#include "implot.h"
#include "cimplot.h"

CIMGUI_API ImPlotContext* ipCreateContext(){
  return ImPlot::CreateContext();
}
CIMGUI_API void ipDestroyContext(ImPlotContext* ctx){
  return ImPlot::DestroyContext(ctx);
}
CIMGUI_API bool ipBeginPlot(const char* title_id, const ImVec2 size, ImPlotFlags flags){
  return ImPlot::BeginPlot(title_id,size,flags);
}
CIMGUI_API void ipEndPlot(){
  ImPlot::EndPlot();
}
CIMGUI_API void ipSetupAxisFormat(ImAxis axis, const char* fmt){
  ImPlot::SetupAxisFormat(axis,fmt);
}
CIMGUI_API void ipSetupAxisTicks(ImAxis axis, double v_min, double v_max, int n_ticks){
  ImPlot::SetupAxisTicks(axis,v_min,v_max,n_ticks);
}
CIMGUI_API void ipSetupAxes(const char* x_label, const char* y_label, ImPlotAxisFlags x_flags, ImPlotAxisFlags y_flags){
  ImPlot::SetupAxes(x_label,y_label,x_flags,y_flags);
}
CIMGUI_API void ipPlotLine(const char* label_id, const double* xs, const double* ys, int count, ImPlotLineFlags flags, int offset){
  ImPlot::PlotLine(label_id,xs,ys,count,flags,offset,sizeof(double));
}
CIMGUI_API void ipGetPlotCurrentLimits(double *xmin,double *xmax,double *ymin,double *ymax){
  ImPlotRect im = ImPlot::GetPlotLimits();
  *xmin = im.X.Min;
  *xmax = im.X.Max;
  *ymin = im.Y.Min;
  *ymax = im.Y.Max;
}
CIMGUI_API void ipSetNextMarkerStyle(ImPlotMarker marker, float size, const ImVec4 fill, float weight, const ImVec4 outline){
  ImPlot::SetNextMarkerStyle(marker,size,fill,weight,outline);
}
CIMGUI_API void ipShowDemoWindow(bool* p_open){
  return ImPlot::ShowDemoWindow(p_open);
}
