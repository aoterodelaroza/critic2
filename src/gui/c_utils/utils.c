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
// startup from the executable's own directory (then the system directory) --
// before any of our code runs, so we cannot redirect it at runtime. Instead the
// GUI is tried first with the system OpenGL; if creating the 3.3 context fails,
// we relaunch a SECOND copy of the executable that is installed in bin\mesa
// alongside the bundled Mesa opengl32.dll. That copy's static import therefore
// resolves to Mesa (a software/D3D12 OpenGL), giving OpenGL 3.3+ everywhere.

// 1 if CRITIC2_SOFTWARE_GL is set in the environment (i.e. this is the
// relaunched, bundled-Mesa instance)
int guiGLSoftwareRequested(void){
    return GetEnvironmentVariableA("CRITIC2_SOFTWARE_GL", NULL, 0) > 0 ? 1 : 0;
}

// prepend dir + ';' to the child's PATH so the bin\mesa copy finds the runtime
// DLLs (libgfortran, glfw3, ...) that live in the parent bin directory
static void guiPrependPath(const wchar_t *dir){
    DWORD have = GetEnvironmentVariableW(L"PATH", NULL, 0); /* chars incl. NUL, 0 if unset */
    size_t need = wcslen(dir) + 2 + (have ? have : 1);
    wchar_t *np = (wchar_t*)malloc(need * sizeof(wchar_t));
    if (!np) return;
    wcscpy(np, dir);
    wcscat(np, L";");
    if (have)
        GetEnvironmentVariableW(L"PATH", np + wcslen(np), (DWORD)(need - wcslen(np)));
    SetEnvironmentVariableW(L"PATH", np);
    free(np);
}

// relaunch the bundled-Mesa copy of this executable (bin\mesa\<name>). The
// original command line (its arguments) is preserved, critic_home is forwarded
// via CRITIC_HOME so the new process finds its data, and GALLIUM_DRIVER selects
// the software renderer. Returns 0 if the fallback is unavailable (no bundled
// Mesa) or the relaunch failed; on success it ends the process and does not
// return.
int guiRelaunchSoftwareGL(const char *critic_home){
    wchar_t exepath[MAX_PATH];
    DWORD n = GetModuleFileNameW(NULL, exepath, MAX_PATH);
    if (n == 0 || n >= MAX_PATH) return 0;
    wchar_t *slash = wcsrchr(exepath, L'\\');
    if (!slash) return 0;

    wchar_t exedir[MAX_PATH], mesaexe[MAX_PATH];
    size_t dlen = (size_t)(slash - exepath);
    if (dlen >= MAX_PATH) return 0;
    wcsncpy(exedir, exepath, dlen);
    exedir[dlen] = L'\0';
    if (_snwprintf(mesaexe, MAX_PATH, L"%s\\mesa\\%s", exedir, slash + 1) < 0) return 0;
    if (GetFileAttributesW(mesaexe) == INVALID_FILE_ATTRIBUTES) return 0; /* no bundled Mesa */

    SetEnvironmentVariableA("CRITIC2_SOFTWARE_GL", "1"); /* mark + loop guard */
    SetEnvironmentVariableA("GALLIUM_DRIVER", "llvmpipe"); /* force software rendering */
    if (critic_home && *critic_home)
        SetEnvironmentVariableA("CRITIC_HOME", critic_home);
    guiPrependPath(exedir);

    STARTUPINFOW si;
    PROCESS_INFORMATION pi;
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));
    /* run the bin\mesa copy but keep the original command line (its arguments) */
    wchar_t *cl = GetCommandLineW();
    size_t len = wcslen(cl) + 1;
    wchar_t *cmd = (wchar_t*)malloc(len * sizeof(wchar_t));
    if (!cmd) return 0;
    memcpy(cmd, cl, len * sizeof(wchar_t));
    BOOL ok = CreateProcessW(mesaexe, cmd, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi);
    free(cmd);
    if (!ok) return 0;
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    ExitProcess(0); /* relaunched successfully; end this process */
    return 1;       /* unreachable */
}
#endif

