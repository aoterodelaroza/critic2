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

static FILE *fidsave = NULL;

void runqhull1(int n, double xstar[n][3], int *nf, int *nv, int *mnfv){
  // write input file
  FILE *fid1 = tmpfile();
  fprintf(fid1,"3\n");
  fprintf(fid1,"%d\n", n+1);
  fprintf(fid1,"%.15e %.15e %.15e \n",0.,0.,0.);
  for (int i = 0; i < n; i++)
    fprintf(fid1,"%.15e %.15e %.15e \n",xstar[i][0],xstar[i][1],xstar[i][2]);
  rewind(fid1);

  // output file
  FILE *fid2 = tmpfile();

  // Use qhull
  // v + Qbb: voronoi
  // QV0: only for the first point (origin)
  // Fi: give the voronoi-relevant vectors (referred to as facet hyperplanes)
  // p: print vertices
  // Fv: print facet vertex indices
  // Fa: print areas
  int nopts = 6;
  char *opts[] = {"qhull","v","Qbb","QV0","Fv","p"};
  qh_init_A(fid1, fid2, stderr, nopts, opts);  
  qh_initflags(qh qhull_command);

  int dim, numpoints;
  boolT ismalloc;
  coordT *points;
  points = qh_readpoints(&numpoints, &dim, &ismalloc);

  qh_init_B(points, numpoints, dim, ismalloc);
  qh_qhull();
  qh_check_output();
  qh_produce_output();
  qh_freeqhull(True);
  int curlong, totlong;
  qh_memfreeshort (&curlong, &totlong);
  fclose(fid1);

  // read the file and write down the dimensions for the arrays
  *mnfv = 0;
  rewind(fid2);
  char buf[1024];
  fgets(buf, sizeof(buf), fid2);
  sscanf(buf,"%d",nf);
  for (int i=0; i<*nf; i++){
    fgets(buf, sizeof(buf), fid2);
    int kk;
    sscanf(buf,"%d",&kk);
    if (kk-2 > *mnfv)
      *mnfv = kk-2;
  }
  fgets(buf, sizeof(buf), fid2);
  fgets(buf, sizeof(buf), fid2);
  sscanf(buf,"%d",nv);
  fidsave = fid2;
}

void runqhull2(int nf, int nv, int mnfv, int ivws[nf], double xvws[nv][3], 
	       int nfvws[nf], int fvws[nf][mnfv]){

  rewind(fidsave);
  char buf[1024];

  // read the faces
  fgets(buf, sizeof(buf), fidsave);
  for (int i=0; i<nf; i++){
    int idum;
    fscanf(fidsave,"%d %d %d", &(nfvws[i]),&idum,&(ivws[i]));
    nfvws[i] -= 2;
    for (int j=0; j<nfvws[i]; j++)
      fscanf(fidsave,"%d", &(fvws[i][j]));
    for (int j=nfvws[i]; j<mnfv; j++)
      fvws[i][j] = 0;
  }
  fgets(buf, sizeof(buf), fidsave);

  // read the vertices
  fgets(buf, sizeof(buf), fidsave);
  fgets(buf, sizeof(buf), fidsave);
  for (int i=0; i<nv; i++){
    fgets(buf, sizeof(buf), fidsave);
    sscanf(buf,"%lf %lf %lf",&(xvws[i][0]),&(xvws[i][1]),&(xvws[i][2]));
  }

  fclose(fidsave);
}

