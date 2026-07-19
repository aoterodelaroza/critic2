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

#if defined(_WIN32) || defined(__WIN32__) || defined(__WINDOWS__)
  /* Windows */
  #include <windows.h>
  #include <shellapi.h>
  #include <direct.h>
  #include <wchar.h>
  #define GETCWD _getcwd
#else
  /* Unix */
  #include <unistd.h>
  #define GETCWD getcwd
#endif
#include <stdio.h>
#include <stdlib.h>

// get the current working directory (portable?)
int getCurrentWorkDir(char *str, size_t siz){
  if (GETCWD(str, siz) == str)
    return 0;
  else
    return 1;
}

// Call the external browser to open a link (portable?)
void openLink(const char* link){
#if defined(_WIN32) || defined(__WIN32__) || defined(__WINDOWS__)
    // ShellExecute opens the link without spawning cmd.exe
    ShellExecuteA(NULL, "open", link, NULL, NULL, SW_SHOWNORMAL);
#else
    char command[1024];
  #if defined(__APPLE__)
    snprintf(command, 1024, "open \"%s\"", link);
  #else
    snprintf(command, 1024, "xdg-open \"%s\"", link);
  #endif
    system(command);
#endif
}

#if defined(_WIN32) || defined(__WIN32__) || defined(__WINDOWS__)
// Pop a modal error message box. The Windows GUI build is a GUI-subsystem
// (-mwindows) app with no console, so a fatal error written to stdout/stderr
// would be invisible; this makes it visible. user32 is already linked by the
// GUI executable (via GLFW).
void guiMessageBox(const char* title, const char* msg){
    MessageBoxA(NULL, msg ? msg : "", title ? title : "critic2",
                MB_OK | MB_ICONERROR | MB_SETFOREGROUND | MB_TOPMOST);
}

// --- OpenGL software/Mesa fallback (driven from gui_start, gui_main@proc.F90) ---
// critic2.exe statically imports opengl32.dll, so Windows loads it at process
// startup from the executable's own directory (then the system directory). The
// package uses that: the native build in bin\ has no opengl32.dll beside it (so
// it gets the system driver), while the software build in bin-mesa\ ships Mesa's
// opengl32.dll beside it (so it gets Mesa). This routine detects the latter case
// and forces Mesa's software rasterizer (llvmpipe), so the software build renders
// even on a machine with no GPU / no OpenGL 3.3. Call before the GL context is
// created; it is a no-op for the native build.
void guiForceSoftwareGLIfBundled(void){
    wchar_t path[MAX_PATH];
    DWORD n = GetModuleFileNameW(NULL, path, MAX_PATH);
    if (n == 0 || n >= MAX_PATH) return;
    wchar_t *slash = wcsrchr(path, L'\\');
    if (!slash) return;
    if ((size_t)(slash - path) + 1 + 12 >= MAX_PATH) return; /* + "opengl32.dll" + NUL */
    wcscpy(slash + 1, L"opengl32.dll");
    if (GetFileAttributesW(path) != INVALID_FILE_ATTRIBUTES)
        SetEnvironmentVariableA("GALLIUM_DRIVER", "llvmpipe");
}
#endif

