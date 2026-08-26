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

#include "implotwrapper.h"

#include "implot.h"

//// implot
// enum ImAxis_
const int const_ImAxis_X1 = ImAxis_X1;
const int const_ImAxis_X2 = ImAxis_X2;
const int const_ImAxis_X3 = ImAxis_X3;
const int const_ImAxis_Y1 = ImAxis_Y1;
const int const_ImAxis_Y2 = ImAxis_Y2;
const int const_ImAxis_Y3 = ImAxis_Y3;
// enum ImPlotStairsFlags_
const int const_ImPlotStairsFlags_None = ImPlotStairsFlags_None;
const int const_ImPlotStairsFlags_Shaded = ImPlotStairsFlags_Shaded;
// enum ImPlotCond_
const int const_ImPlotCond_Always = ImPlotCond_Always;
const int const_ImPlotCond_Once = ImPlotCond_Once;
// enum ImPlotScale_
const int const_ImPlotScale_Linear = ImPlotScale_Linear;
const int const_ImPlotScale_Log10 = ImPlotScale_Log10;
const int const_ImPlotScale_SymLog = ImPlotScale_SymLog;
// enum ImPlotDragToolFlags_
const int const_ImPlotDragToolFlags_None = ImPlotDragToolFlags_None;
const int const_ImPlotDragToolFlags_NoInputs = ImPlotDragToolFlags_NoInputs;
// enum ImPlotColormap_ (the continuous ones)
const int const_ImPlotColormap_Viridis = ImPlotColormap_Viridis;
const int const_ImPlotColormap_Plasma = ImPlotColormap_Plasma;
const int const_ImPlotColormap_Hot = ImPlotColormap_Hot;
const int const_ImPlotColormap_Cool = ImPlotColormap_Cool;
const int const_ImPlotColormap_Jet = ImPlotColormap_Jet;
const int const_ImPlotColormap_RdBu = ImPlotColormap_RdBu;
const int const_ImPlotColormap_Spectral = ImPlotColormap_Spectral;
const int const_ImPlotColormap_Greys = ImPlotColormap_Greys;
// enum ImPlotFlags_
const int const_ImPlotFlags_None = ImPlotFlags_None;
const int const_ImPlotFlags_NoTitle = ImPlotFlags_NoTitle;
const int const_ImPlotFlags_NoLegend = ImPlotFlags_NoLegend;
const int const_ImPlotFlags_NoMouseText = ImPlotFlags_NoMouseText;
const int const_ImPlotFlags_NoInputs = ImPlotFlags_NoInputs;
const int const_ImPlotFlags_NoMenus = ImPlotFlags_NoMenus;
const int const_ImPlotFlags_NoBoxSelect = ImPlotFlags_NoBoxSelect;
const int const_ImPlotFlags_NoChild = ImPlotFlags_NoChild;
const int const_ImPlotFlags_NoFrame = ImPlotFlags_NoFrame;
const int const_ImPlotFlags_Equal = ImPlotFlags_Equal;
const int const_ImPlotFlags_Crosshairs = ImPlotFlags_Crosshairs;
const int const_ImPlotFlags_CanvasOnly = ImPlotFlags_CanvasOnly;
// enum ImPlotAxisFlags_
const int const_ImPlotAxisFlags_None = ImPlotAxisFlags_None;
const int const_ImPlotAxisFlags_NoLabel = ImPlotAxisFlags_NoLabel;
const int const_ImPlotAxisFlags_NoGridLines = ImPlotAxisFlags_NoGridLines;
const int const_ImPlotAxisFlags_NoTickMarks = ImPlotAxisFlags_NoTickMarks;
const int const_ImPlotAxisFlags_NoTickLabels = ImPlotAxisFlags_NoTickLabels;
const int const_ImPlotAxisFlags_NoInitialFit = ImPlotAxisFlags_NoInitialFit;
const int const_ImPlotAxisFlags_NoMenus = ImPlotAxisFlags_NoMenus;
const int const_ImPlotAxisFlags_NoSideSwitch = ImPlotAxisFlags_NoSideSwitch;
const int const_ImPlotAxisFlags_NoHighlight = ImPlotAxisFlags_NoHighlight;
const int const_ImPlotAxisFlags_Opposite = ImPlotAxisFlags_Opposite;
const int const_ImPlotAxisFlags_Foreground = ImPlotAxisFlags_Foreground;
const int const_ImPlotAxisFlags_Invert = ImPlotAxisFlags_Invert;
const int const_ImPlotAxisFlags_AutoFit = ImPlotAxisFlags_AutoFit;
const int const_ImPlotAxisFlags_RangeFit = ImPlotAxisFlags_RangeFit;
const int const_ImPlotAxisFlags_PanStretch = ImPlotAxisFlags_PanStretch;
const int const_ImPlotAxisFlags_LockMin = ImPlotAxisFlags_LockMin;
const int const_ImPlotAxisFlags_LockMax = ImPlotAxisFlags_LockMax;
const int const_ImPlotAxisFlags_Lock = ImPlotAxisFlags_Lock;
const int const_ImPlotAxisFlags_NoDecorations = ImPlotAxisFlags_NoDecorations;
const int const_ImPlotAxisFlags_AuxDefault = ImPlotAxisFlags_AuxDefault;
// enum ImPlotLineFlags_
const int const_ImPlotLineFlags_None = ImPlotLineFlags_None;
const int const_ImPlotLineFlags_Segments = ImPlotLineFlags_Segments;
const int const_ImPlotLineFlags_Loop = ImPlotLineFlags_Loop;
const int const_ImPlotLineFlags_SkipNaN = ImPlotLineFlags_SkipNaN;
const int const_ImPlotLineFlags_NoClip = ImPlotLineFlags_NoClip;
const int const_ImPlotLineFlags_Shaded = ImPlotLineFlags_Shaded;
// enum ImPlotMarker_
const int const_ImPlotMarker_None = ImPlotMarker_None;
const int const_ImPlotMarker_Circle = ImPlotMarker_Circle;
const int const_ImPlotMarker_Square = ImPlotMarker_Square;
const int const_ImPlotMarker_Diamond = ImPlotMarker_Diamond;
const int const_ImPlotMarker_Up = ImPlotMarker_Up;
const int const_ImPlotMarker_Down = ImPlotMarker_Down;
const int const_ImPlotMarker_Left = ImPlotMarker_Left;
const int const_ImPlotMarker_Right = ImPlotMarker_Right;
const int const_ImPlotMarker_Cross = ImPlotMarker_Cross;
const int const_ImPlotMarker_Plus = ImPlotMarker_Plus;
const int const_ImPlotMarker_Asterisk = ImPlotMarker_Asterisk;
// const int const_ImPlotMarker_COUNT

