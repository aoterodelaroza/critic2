/* 
Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

#include <stdio.h>
#include <stdlib.h>
#include <libqhull.h>

void doqhull(char *file1, char *file2, char *file3, int *ithr){

  // open the file with the star points
  FILE *fin;
  fin = fopen(file1,"r");

  // open the file with the vertices
  FILE *fout;
  fout = fopen(file2,"w");

  // initialize the options
  int nopts = 5;
  char *opts[] = {"qhull","v","Qbb","QV0","p"};
  qh_init_A(fin, fout, stderr, nopts, opts);  
  qh_initflags(qh qhull_command);

  // initialize the point array using the scratch file
  int dim, numpoints;
  boolT ismalloc;
  coordT *points;
  points = qh_readpoints(&numpoints, &dim, &ismalloc);

  qh_init_B(points, numpoints, dim, ismalloc);
  qh_qhull();
  qh_check_output();
  qh_produce_output();
  qh_freeqhull(True);

  fclose(fin);
  fclose(fout);

  // open the file with the vertices
  fin = fopen(file2,"r");
  
  // open the file with the vertices
  fout = fopen(file3,"w");

  // initialize the options
  // the OFF format orients the facets so I don't have to do it myself
  nopts = 4;
  char stra[10];
  sprintf(stra,"C-1e-%d",*ithr);
  char *opts2[] = {"qhull","Qs","o",stra};
  qh_init_A(fin, fout, stderr, nopts, opts2);  
  qh_initflags(qh qhull_command);

  // initialize the point array using the scratch file
  points = qh_readpoints(&numpoints, &dim, &ismalloc);

  qh_init_B(points, numpoints, dim, ismalloc);
  qh_qhull();
  qh_check_output();
  qh_produce_output();
  qh_freeqhull(True);

  fclose(fin);
  fclose(fout);
}

