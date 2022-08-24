/*
Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
#include <math.h>
#include "libqhull_r.h"

// From a list of n vertices (xstar), calculate the Voronoi polyhedron
// of the first point in xstar and return its number of faces (nf),
// number of vertices (nv), and the maximum number of vertex per face (mnfv).
// The temporary file containing the vertex/edge/face information
// remains open until the user calls step2. The handle is saved in fidsave_voronoi.
void runqhull_voronoi_step1(int n, double xstar[n][3], int *nf, int *nv, int *mnfv, FILE **fid){
  qhT qhT_;

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
  qh_init_A(&qhT_, fid1, fid2, stderr, 0, NULL);
  qhT_.NOerrexit= False;
  char options [2000];
  sprintf(options,"qhull v Qbb QV0 Fv p Qs");
  qh_initflags(&qhT_, options);

  int dim, numpoints;
  boolT ismalloc;
  coordT *points;
  points = qh_readpoints(&qhT_, &numpoints, &dim, &ismalloc);

  qh_init_B(&qhT_, points, numpoints, dim, ismalloc);
  qh_qhull(&qhT_);
  qh_check_output(&qhT_);
  qh_produce_output(&qhT_);
  qh_freeqhull(&qhT_, True);
  int curlong, totlong;
  qh_memfreeshort(&qhT_, &curlong, &totlong);
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
  *fid = fid2;
}

// Read the Voronoi polyhedron calculated in step 1. Input: number of faces
// (nf), number of vertices (nv) and maximum number of vertices per face (mnfv).
// Returns: the neighbor vertex identifier for a given face (ivws), the vertices
// of the polyhedron (xvws), the number of vertices for each face (nfvws), and
// the list of vertices for each face (fvws).
void runqhull_voronoi_step2(int nf, int nv, int mnfv, int ivws[nf], double xvws[nv][3],
	       int nfvws[nf], int fvws[nf][mnfv], FILE *fid){

  rewind(fid);
  char buf[1024];

  // read the faces
  fgets(buf, sizeof(buf), fid);
  for (int i=0; i<nf; i++){
    int idum;
    fscanf(fid,"%d %d %d", &(nfvws[i]),&idum,&(ivws[i]));
    nfvws[i] -= 2;
    for (int j=0; j<nfvws[i]; j++)
      fscanf(fid,"%d", &(fvws[i][j]));
    for (int j=nfvws[i]; j<mnfv; j++)
      fvws[i][j] = 0;
  }
  fgets(buf, sizeof(buf), fid);

  // read the vertices
  fgets(buf, sizeof(buf), fid);
  fgets(buf, sizeof(buf), fid);
  for (int i=0; i<nv; i++){
    fgets(buf, sizeof(buf), fid);
    sscanf(buf,"%lf %lf %lf",&(xvws[i][0]),&(xvws[i][1]),&(xvws[i][2]));
  }

  fclose(fid);
}

// From a list of n vertices (xvert), and an inside point (x0),
// calculate the convex hull of the vertices by first projecting the
// input vertices onto the unit sphere, then running qhull. The number
// of faces in the convex hull (nf) is returned.  The temporary file
// containing the face information remains open until the user calls
// step2. The handle is saved in fidsave_basintri.
void runqhull_basintriangulate_step1(int n, double x0[3], double xvert[n][3], int *nf, FILE **fid){
  qhT qhT_;

  // write input file
  FILE *fid1 = tmpfile();
  fprintf(fid1,"3\n");
  fprintf(fid1,"%d\n", n);
  for (int i = 0; i < n; i++){
    double x[3] = {xvert[i][0] - x0[0], xvert[i][1] - x0[1], xvert[i][2] - x0[2]};
    double norm = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    fprintf(fid1,"%.15e %.15e %.15e \n",x[0]/norm,x[1]/norm,x[2]/norm);
  }
  rewind(fid1);

  // output file
  FILE *fid2 = tmpfile();

  // Use qhull
  // qhull (no options = convex hull)
  // Fv: print facet vertex indices
  // Qt: faces are triangles
  qh_init_A(&qhT_, fid1, fid2, stderr, 0, NULL);
  qhT_.NOerrexit= False;
  char options [2000];
  sprintf(options,"qhull Fv Qt");
  qh_initflags(&qhT_, options);

  int dim, numpoints;
  boolT ismalloc;
  coordT *points;
  points = qh_readpoints(&qhT_, &numpoints, &dim, &ismalloc);

  qh_init_B(&qhT_, points, numpoints, dim, ismalloc);
  qh_qhull(&qhT_);
  qh_check_output(&qhT_);
  qh_produce_output(&qhT_);
  qh_freeqhull(&qhT_, True);
  int curlong, totlong;
  qh_memfreeshort(&qhT_, &curlong, &totlong);
  fclose(fid1);

  // read the file and write down the dimensions for the arrays
  rewind(fid2);
  fscanf(fid2,"%d",nf);
  *fid = fid2;
}

// Read the convex hull faces calculated in step 1. Input: number of
// faces (nf). Returns: the list of vertex identifiers for each face
// (iface).
void runqhull_basintriangulate_step2(int nf, int iface[nf][3], FILE *fid){
  rewind(fid);
  char buf[1024];

  // read the faces
  fgets(buf, sizeof(buf), fid);
  for (int i=0; i<nf; i++){
    int idum;
    fscanf(fid,"%d %d %d %d", &idum, &(iface[i][0]), &(iface[i][1]), &(iface[i][2]));
    for (int j=0; j<3; j++)
      iface[i][j] += 1;
  }

  fclose(fid);
}
