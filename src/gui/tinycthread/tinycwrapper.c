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
#include <stdlib.h>

/* Function return values */
const int const_thrd_error = thrd_error;
const int const_thrd_success = thrd_success;
const int const_thrd_timedout = thrd_timedout;
const int const_thrd_busy = thrd_busy;
const int const_thrd_nomem = thrd_nomem;
/* Mutex types */
const int const_mtx_plain = mtx_plain;
const int const_mtx_timed = mtx_timed;
const int const_mtx_recursive = mtx_recursive;

// Allocate and deallocate mutex objects
void *allocate_mtx(){
  return (void *) malloc(sizeof(mtx_t));
}
void deallocate_mtx(void *p){
  if (p) free(p);
}

// Wrapper for thread functions
int wrap_thrd_join(thrd_t *thr, int *res){
  return thrd_join(*thr,res);
}

// Full memory barrier: the writes before it are visible to other
// threads before the writes after it, and the reads after it are not
// done before the reads before it. Being an external call, it is also
// a compiler barrier for the Fortran caller.
void memory_fence(){
  __atomic_thread_fence(__ATOMIC_SEQ_CST);
}
