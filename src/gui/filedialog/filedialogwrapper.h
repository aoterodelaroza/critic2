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

// Wrapper to filedialog

#ifndef FILEDIALOGWRAPPER_H
#define FILEDIALOGWRAPPER_H

//// filedialog
// enum ImGuiFileDialogFlags;
extern "C" const int const_ImGuiFileDialogFlags_None;
extern "C" const int const_ImGuiFileDialogFlags_ConfirmOverwrite;
extern "C" const int const_ImGuiFileDialogFlags_DontShowHiddenFiles;
extern "C" const int const_ImGuiFileDialogFlags_DisableCreateDirectoryButton;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnType;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnSize;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnDate;
extern "C" const int const_ImGuiFileDialogFlags_NoDialog;
extern "C" const int const_ImGuiFileDialogFlags_ReadOnlyFileNameField;
extern "C" const int const_ImGuiFileDialogFlags_CaseInsensitiveExtention;
extern "C" const int const_ImGuiFileDialogFlags_Modal;
// enum IGFD_FileStyle
extern "C" const int const_IGFD_FileStyle_None;
extern "C" const int const_IGFD_FileStyleByTypeFile;
extern "C" const int const_IGFD_FileStyleByTypeDir;
extern "C" const int const_IGFD_FileStyleByTypeLink;
extern "C" const int const_IGFD_FileStyleByExtention;
extern "C" const int const_IGFD_FileStyleByFullName;
extern "C" const int const_IGFD_FileStyleByContainedInFullName;

#endif
