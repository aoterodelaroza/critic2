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

/* Create one directory (the parent must exist). Returns 0 on success
   or if the directory already exists, nonzero otherwise. Called from
   tools_io%mkpath, which creates the parents first. */
#include <errno.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#endif

int critic2_mkdir(const char *path) {
  int ier;
#ifdef _WIN32
  ier = _mkdir(path);
#else
  ier = mkdir(path, 0777);
#endif
  if (ier != 0 && errno == EEXIST) ier = 0;
  return ier;
}
