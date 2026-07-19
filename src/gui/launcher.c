// Portable launcher for the Windows critic2 package. One of these sits at the
// package root; it runs an executable in bin/ or bin-mesa/ addressed relative to
// the launcher's OWN location, so it works wherever the package is unpacked --
// with no console window, because the launcher is a GUI-subsystem program (a
// .lnk shortcut cannot do this portably: it bakes in an absolute target path).
// LAUNCH_TARGET (the exe to run) and LAUNCH_WORKDIR (its working directory, so a
// bin-mesa exe finds the shared runtime DLLs in bin/) are relative to the
// launcher's folder and are baked in at compile time. Forward slashes are used
// throughout; CreateProcess accepts them.
#include <windows.h>
#include <stdio.h>

#ifndef LAUNCH_TARGET
#define LAUNCH_TARGET  L"bin/critic2-gui.exe"
#endif
#ifndef LAUNCH_WORKDIR
#define LAUNCH_WORKDIR L"bin"
#endif

int WINAPI wWinMain(HINSTANCE hInst, HINSTANCE hPrev, PWSTR lpCmd, int nShow){
    (void)hInst; (void)hPrev; (void)nShow;

    wchar_t dir[MAX_PATH], target[MAX_PATH], workdir[MAX_PATH];
    DWORD n = GetModuleFileNameW(NULL, dir, MAX_PATH);
    if (n == 0 || n >= MAX_PATH) return 1;
    wchar_t *slash = wcsrchr(dir, L'\\');
    if (!slash) return 1;
    *slash = L'\0';                               /* dir = launcher's folder */
    if (_snwprintf(target,  MAX_PATH, L"%s/%s", dir, LAUNCH_TARGET)  < 0) return 1;
    if (_snwprintf(workdir, MAX_PATH, L"%s/%s", dir, LAUNCH_WORKDIR) < 0) return 1;

    /* child command line: "target" <forwarded arguments> */
    static wchar_t cl[32768];
    _snwprintf(cl, 32768, L"\"%s\" %s", target, (lpCmd && *lpCmd) ? lpCmd : L"");

    DWORD flags = 0;
#ifdef LAUNCH_CONSOLE
    flags = CREATE_NEW_CONSOLE;                   /* the CLI gets its own console */
#endif

    STARTUPINFOW si;
    PROCESS_INFORMATION pi;
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));
    if (!CreateProcessW(target, cl, NULL, NULL, FALSE, flags, NULL, workdir, &si, &pi))
        return 1;
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    return 0;
}
