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

// Wrapper to the IMPLOT library

#ifndef IMPLOTWRAPPER_H
#define IMPLOTWRAPPER_H

// enum ImAxis_
extern "C" const int const_ImAxis_X1;
extern "C" const int const_ImAxis_X2;
extern "C" const int const_ImAxis_X3;
extern "C" const int const_ImAxis_Y1;
extern "C" const int const_ImAxis_Y2;
extern "C" const int const_ImAxis_Y3;
// enum ImPlotFlags_
extern "C" const int const_ImPlotFlags_None;
extern "C" const int const_ImPlotFlags_NoTitle;
extern "C" const int const_ImPlotFlags_NoLegend;
extern "C" const int const_ImPlotFlags_NoMouseText;
extern "C" const int const_ImPlotFlags_NoInputs;
extern "C" const int const_ImPlotFlags_NoMenus;
extern "C" const int const_ImPlotFlags_NoBoxSelect;
extern "C" const int const_ImPlotFlags_NoChild;
extern "C" const int const_ImPlotFlags_NoFrame;
extern "C" const int const_ImPlotFlags_Equal;
extern "C" const int const_ImPlotFlags_Crosshairs;
extern "C" const int const_ImPlotFlags_CanvasOnly;
// enum ImPlotAxisFlags_
extern "C" const int const_ImPlotAxisFlags_None;
extern "C" const int const_ImPlotAxisFlags_NoLabel;
extern "C" const int const_ImPlotAxisFlags_NoGridLines;
extern "C" const int const_ImPlotAxisFlags_NoTickMarks;
extern "C" const int const_ImPlotAxisFlags_NoTickLabels;
extern "C" const int const_ImPlotAxisFlags_NoInitialFit;
extern "C" const int const_ImPlotAxisFlags_NoMenus;
extern "C" const int const_ImPlotAxisFlags_NoSideSwitch;
extern "C" const int const_ImPlotAxisFlags_NoHighlight;
extern "C" const int const_ImPlotAxisFlags_Opposite;
extern "C" const int const_ImPlotAxisFlags_Foreground;
extern "C" const int const_ImPlotAxisFlags_Invert;
extern "C" const int const_ImPlotAxisFlags_AutoFit;
extern "C" const int const_ImPlotAxisFlags_RangeFit;
extern "C" const int const_ImPlotAxisFlags_PanStretch;
extern "C" const int const_ImPlotAxisFlags_LockMin;
extern "C" const int const_ImPlotAxisFlags_LockMax;
extern "C" const int const_ImPlotAxisFlags_Lock;
extern "C" const int const_ImPlotAxisFlags_NoDecorations;
extern "C" const int const_ImPlotAxisFlags_AuxDefault;
// enum ImPlotLineFlags_
extern "C" const int const_ImPlotLineFlags_None;
extern "C" const int const_ImPlotLineFlags_Segments;
extern "C" const int const_ImPlotLineFlags_Loop;
extern "C" const int const_ImPlotLineFlags_SkipNaN;
extern "C" const int const_ImPlotLineFlags_NoClip;
extern "C" const int const_ImPlotLineFlags_Shaded;
// enum ImPlotMarker_
extern "C" const int const_ImPlotMarker_None;
extern "C" const int const_ImPlotMarker_Circle;
extern "C" const int const_ImPlotMarker_Square;
extern "C" const int const_ImPlotMarker_Diamond;
extern "C" const int const_ImPlotMarker_Up;
extern "C" const int const_ImPlotMarker_Down;
extern "C" const int const_ImPlotMarker_Left;
extern "C" const int const_ImPlotMarker_Right;
extern "C" const int const_ImPlotMarker_Cross;
extern "C" const int const_ImPlotMarker_Plus;
extern "C" const int const_ImPlotMarker_Asterisk;
// const int const_ImPlotMarker_COUNT

#endif
