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

// Wrapper to the tinythread library

#ifndef TINYCWRAPPER_H
#define TINYCWRAPPER_H

#include "tinycthread.h"

// constants //
/* Function return values */
extern const int const_thrd_error;
extern const int const_thrd_success;
extern const int const_thrd_timedout;
extern const int const_thrd_busy;
extern const int const_thrd_nomem;
/* Mutex types */
extern const int const_mtx_plain;
extern const int const_mtx_timed;
extern const int const_mtx_recursive;

// wrapper functions //
// mutexes
void *allocate_mtx();
void deallocate_mtx(void *p);
// threads
int wrap_thrd_join(thrd_t *thr, int *res);

#endif
