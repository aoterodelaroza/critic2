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

#include "filedialogwrapper.h"

#include "ImGuiFileDialog.h"

//// filedialog
// enum ImGuiFileDialogFlags;
const int const_ImGuiFileDialogFlags_None = ImGuiFileDialogFlags_None;
const int const_ImGuiFileDialogFlags_ConfirmOverwrite = ImGuiFileDialogFlags_ConfirmOverwrite;
const int const_ImGuiFileDialogFlags_DontShowHiddenFiles = ImGuiFileDialogFlags_DontShowHiddenFiles;
const int const_ImGuiFileDialogFlags_DisableCreateDirectoryButton = ImGuiFileDialogFlags_DisableCreateDirectoryButton;
const int const_ImGuiFileDialogFlags_HideColumnType = ImGuiFileDialogFlags_HideColumnType;
const int const_ImGuiFileDialogFlags_HideColumnSize = ImGuiFileDialogFlags_HideColumnSize;
const int const_ImGuiFileDialogFlags_HideColumnDate = ImGuiFileDialogFlags_HideColumnDate;
const int const_ImGuiFileDialogFlags_NoDialog = ImGuiFileDialogFlags_NoDialog;
const int const_ImGuiFileDialogFlags_ReadOnlyFileNameField = ImGuiFileDialogFlags_ReadOnlyFileNameField;
const int const_ImGuiFileDialogFlags_CaseInsensitiveExtention = ImGuiFileDialogFlags_CaseInsensitiveExtention;
const int const_ImGuiFileDialogFlags_Modal = ImGuiFileDialogFlags_Modal;
// enum IGFD_FileStyle
const int const_IGFD_FileStyle_None = IGFD_FileStyle_None;
const int const_IGFD_FileStyleByTypeFile = IGFD_FileStyleByTypeFile;
const int const_IGFD_FileStyleByTypeDir = IGFD_FileStyleByTypeDir;
const int const_IGFD_FileStyleByTypeLink = IGFD_FileStyleByTypeLink;
const int const_IGFD_FileStyleByExtention = IGFD_FileStyleByExtention;
const int const_IGFD_FileStyleByFullName = IGFD_FileStyleByFullName;
const int const_IGFD_FileStyleByContainedInFullName = IGFD_FileStyleByContainedInFullName;
