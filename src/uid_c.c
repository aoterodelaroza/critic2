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

/* The counter behind global%new_uid: a positive integer not returned
   before in this run, from any thread (the GUI loads fields in its
   initialization threads, which are not OpenMP threads, so an OpenMP
   atomic is not enough, nor available without OpenMP). */
static long long critic2_uid_last = 0;

long long critic2_next_uid(void) {
  return __atomic_add_fetch(&critic2_uid_last, 1, __ATOMIC_SEQ_CST);
}
