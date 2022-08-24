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
  #include <direct.h>
  #define GETCWD _getcwd
#else
  /* Unix */
  #include <unistd.h>
  #define GETCWD getcwd
#endif

int getCurrentWorkDir(char *str, size_t siz){
  if (GETCWD(str, siz) == str)
    return 0;
  else
    return 1;
}
