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
// the OpenGL3 backend (built with C linkage, see cimgui_impl.h)
extern "C" {
  bool ImGui_ImplOpenGL3_Init(const char* glsl_version);
  void ImGui_ImplOpenGL3_Shutdown(void);
  void ImGui_ImplOpenGL3_NewFrame(void);
}
#ifdef IMGUI_ENABLE_FREETYPE
#include "../imgui/imgui_freetype.h"
#endif

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
CIMGUI_API void ipSetupAxisScale(ImAxis axis, ImPlotScale scale){
  ImPlot::SetupAxisScale(axis,scale);
}
CIMGUI_API void ipSetupAxisLimits(ImAxis axis, double v_min, double v_max, ImPlotCond cond){
  ImPlot::SetupAxisLimits(axis,v_min,v_max,cond);
}
CIMGUI_API bool ipDragLineX(int id, double* x, const ImVec4 col, float thickness, ImPlotDragToolFlags flags){
  return ImPlot::DragLineX(id,x,col,thickness,flags);
}
CIMGUI_API void ipPlotToPixels(double x, double y, ImAxis x_axis, ImAxis y_axis, float *px, float *py){
  ImVec2 p = ImPlot::PlotToPixels(x,y,x_axis,y_axis);
  *px = p.x;
  *py = p.y;
}
CIMGUI_API void ipSampleColormap(float t, ImPlotColormap cmap, ImVec4 *out){
  *out = ImPlot::SampleColormap(t,cmap);
}
CIMGUI_API void ipSetupAxisTicksValues(ImAxis axis, const double* values, int n_ticks, bool keep_default){
  ImPlot::SetupAxisTicks(axis,values,n_ticks,NULL,keep_default);
}
CIMGUI_API void ipSetupAxes(const char* x_label, const char* y_label, ImPlotAxisFlags x_flags, ImPlotAxisFlags y_flags){
  ImPlot::SetupAxes(x_label,y_label,x_flags,y_flags);
}
CIMGUI_API void ipPlotLine(const char* label_id, const double* xs, const double* ys, int count, ImPlotLineFlags flags, int offset){
  ImPlot::PlotLine(label_id,xs,ys,count,flags,offset,sizeof(double));
}
CIMGUI_API void ipPlotStairs(const char* label_id, const double* xs, const double* ys, int count, ImPlotStairsFlags flags, int offset){
  ImPlot::PlotStairs(label_id,xs,ys,count,flags,offset,sizeof(double));
}
CIMGUI_API void ipGetPlotCurrentLimits(double *xmin,double *xmax,double *ymin,double *ymax){
  ImPlotRect im = ImPlot::GetPlotLimits();
  *xmin = im.X.Min;
  *xmax = im.X.Max;
  *ymin = im.Y.Min;
  *ymax = im.Y.Max;
}
CIMGUI_API void ipSetNextLineStyle(const ImVec4 col, float weight){
  ImPlot::SetNextLineStyle(col,weight);
}
CIMGUI_API void ipSetNextMarkerStyle(ImPlotMarker marker, float size, const ImVec4 fill, float weight, const ImVec4 outline){
  ImPlot::SetNextMarkerStyle(marker,size,fill,weight,outline);
}
CIMGUI_API void ipColormapScale(const char* label, double scale_min, double scale_max, const ImVec2 size, const char* format, ImPlotColormap cmap){
  ImPlot::ColormapScale(label,scale_min,scale_max,size,format,0,cmap);
}
CIMGUI_API void ipGetPlotMousePos(double *x, double *y){
  ImPlotPoint p = ImPlot::GetPlotMousePos();
  *x = p.x;
  *y = p.y;
}
CIMGUI_API bool ipIsPlotHovered(){
  return ImPlot::IsPlotHovered();
}
CIMGUI_API void ipShowDemoWindow(bool* p_open){
  return ImPlot::ShowDemoWindow(p_open);
}

// Plots for a publication: a second ImGui and ImPlot context, black
// on white, that draws into a window covering its display and renders
// it scale times larger (DisplayFramebufferScale), with the font (the
// base85-compressed TTF ttf85, glyphs in ranges) baked at that size,
// up to bakemax pixels. The contexts current before the call are
// restored by ipPrintContextDestroy.
struct ipPrintContext {
  ImGuiContext *gui0, *gui;
  ImPlotContext *plot0, *plot;
  float scale;
};
CIMGUI_API void* ipPrintContextCreate(float fontsize, float scale, float bakemax, const char* ttf85,
                                      const ImWchar* ranges, const char* glsl){
  const ImVec4 white(1.f,1.f,1.f,1.f), black(0.f,0.f,0.f,1.f), none(0.f,0.f,0.f,0.f);

  ipPrintContext *pc = IM_NEW(ipPrintContext)();
  pc->scale = scale;
  pc->gui0 = ImGui::GetCurrentContext();
  pc->plot0 = ImPlot::GetCurrentContext();
  pc->gui = ImGui::CreateContext();
  ImGui::SetCurrentContext(pc->gui);
  pc->plot = ImPlot::CreateContext();
  ImPlot::SetCurrentContext(pc->plot);

  ImGuiIO& io = ImGui::GetIO();
  io.IniFilename = NULL;
  io.LogFilename = NULL;
#ifdef IMGUI_ENABLE_FREETYPE
  io.Fonts->FontBuilderFlags = ImGuiFreeTypeBuilderFlags_Bold;
#endif
  float bake = ImMin(fontsize * scale, ImMax(bakemax, fontsize));
  ImFont* font = io.Fonts->AddFontFromMemoryCompressedBase85TTF(ttf85, bake, NULL, ranges);
  font->Scale = fontsize / bake;

  ImGui::StyleColorsLight();
  ImGuiStyle& st = ImGui::GetStyle();
  st.WindowBorderSize = 0.f;
  st.Colors[ImGuiCol_WindowBg] = white;
  st.Colors[ImGuiCol_FrameBg] = white;
  st.Colors[ImGuiCol_Text] = black;

  ImPlot::StyleColorsLight();
  ImPlotStyle& ps = ImPlot::GetStyle();
  ps.Colors[ImPlotCol_FrameBg] = white;
  ps.Colors[ImPlotCol_PlotBg] = white;
  ps.Colors[ImPlotCol_PlotBorder] = black;
  ps.Colors[ImPlotCol_AxisText] = black;
  ps.Colors[ImPlotCol_AxisTick] = black;
  ps.Colors[ImPlotCol_AxisGrid] = none;
  ps.Colors[ImPlotCol_AxisBg] = none;
  ps.Colors[ImPlotCol_AxisBgHovered] = none;
  ps.Colors[ImPlotCol_AxisBgActive] = none;

  ImGui_ImplOpenGL3_Init(glsl);
  return pc;
}
CIMGUI_API void ipPrintFrameBegin(void* p, float width, float height){
  ipPrintContext *pc = (ipPrintContext*) p;
  ImGuiIO& io = ImGui::GetIO();
  io.DisplaySize = ImVec2(width, height);
  io.DisplayFramebufferScale = ImVec2(pc->scale, pc->scale);
  io.DeltaTime = 1.f / 60.f;
  ImGui_ImplOpenGL3_NewFrame();
  ImGui::NewFrame();
  ImGui::SetNextWindowPos(ImVec2(0.f, 0.f));
  ImGui::SetNextWindowSize(ImVec2(width, height));
  ImGui::Begin("##print", NULL, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoSavedSettings |
               ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoNav);
  // the antialiasing fringe is one pixel of the image, not of the window
  ImGui::GetWindowDrawList()->_FringeScale = 1.f / pc->scale;
}
CIMGUI_API ImDrawData* ipPrintFrameEnd(){
  ImGui::End();
  ImGui::Render();
  return ImGui::GetDrawData();
}
CIMGUI_API void ipPrintContextDestroy(void* p){
  ipPrintContext *pc = (ipPrintContext*) p;
  ImGui_ImplOpenGL3_Shutdown();
  ImPlot::DestroyContext(pc->plot);
  ImGui::DestroyContext(pc->gui);
  ImGui::SetCurrentContext(pc->gui0);
  ImPlot::SetCurrentContext(pc->plot0);
  IM_DELETE(pc);
}
