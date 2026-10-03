/*
  Copyright (c) 2007-2026 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
  Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
  <victor@fluor.quimica.uniovi.es>.

  critic2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at
  your option) any later version.

  critic2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Write the absolute path of the running executable into buf (at most
   n bytes, null-terminated) and return its length, or 0 if it cannot
   be determined. Unlike argv[0], this works when critic2 was started
   by name through the PATH. Called from global_init to locate the data
   files of a relocatable (unpacked) install. */
#include <string.h>
#if defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#include <stdlib.h>
#include <limits.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

int critic2_exepath(char *buf, int n) {
  if (n <= 1) return 0;
  buf[0] = '\0';
#if defined(_WIN32)
  DWORD len = GetModuleFileNameA(NULL, buf, (DWORD) n);
  if (len == 0 || len >= (DWORD) n) return 0;
  return (int) len;
#elif defined(__APPLE__)
  char raw[PATH_MAX], real[PATH_MAX];
  uint32_t size = sizeof(raw);
  if (_NSGetExecutablePath(raw, &size) != 0) return 0;
  if (!realpath(raw, real)) return 0;
  size_t len = strlen(real);
  if (len >= (size_t) n) return 0;
  memcpy(buf, real, len + 1);
  return (int) len;
#elif defined(__linux__)
  ssize_t len = readlink("/proc/self/exe", buf, (size_t) n - 1);
  if (len <= 0 || len >= (ssize_t) n - 1) return 0;
  buf[len] = '\0';
  return (int) len;
#else
  return 0;
#endif
}
