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

#include "tinycwrapper.h"
#include "tinycthread.h"
#include <stdlib.h>

// Allocate and deallocate mutex objects
void *allocate_mtx(){
  return (void *) malloc(sizeof(mtx_t));
}
void deallocate_mtx(void *p){
  free(p);
}

// Allocate and deallocate mutex threads
void *allocate_thrd(){
  return (void *) malloc(sizeof(thrd_t));
}
void deallocate_thrd(void *p){
  free(p);
}

