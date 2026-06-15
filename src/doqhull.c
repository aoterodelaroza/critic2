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
#include "qhull_ra.h"

// This file drives the (reentrant) qhull library directly through its
// in-memory API: input points are passed as a coordT array and the results
// are read straight from qhull's facet/vertex data structures. No temporary
// files or text parsing are involved. Each calculation is split into a step1
// routine (run qhull, collect the result into a malloc'd context, return the
// array dimensions) and a step2 routine (copy the context into the
// Fortran-allocated arrays and free the context). qhull errors are trapped
// with setjmp and reported back through *ier instead of aborting the program.

/*=====================================================================
  Voronoi polyhedron (Wigner-Seitz cell)
  =====================================================================*/

// Collected Voronoi result for the cell of input site 0.
typedef struct {
  int nf;        // number of faces (Voronoi ridges of site 0)
  int cap;       // capacity of the per-face arrays
  int mnfv;      // maximum number of vertices over all faces
  int nv;        // number of Voronoi vertices (cell-0 centers)
  int *ivws;     // [nf] neighbor site id (index into xstar) for each face
  int *nfvws;    // [nf] number of vertices for each face
  int **fv;      // [nf][nfvws] 0-based Voronoi-vertex index for each face
  double *xvws;  // [nv][3] Voronoi vertex coordinates, indexed by id
  int err;       // nonzero if an unbounded ridge was found (not a finite cell)
} voronoi_ctx;

static void voronoi_ctx_grow(voronoi_ctx *c){
  if (c->nf < c->cap)
    return;
  int ncap = c->cap ? 2*c->cap : 16;
  c->ivws  = (int*)  realloc(c->ivws,  (size_t)ncap*sizeof(int));
  c->nfvws = (int*)  realloc(c->nfvws, (size_t)ncap*sizeof(int));
  c->fv    = (int**) realloc(c->fv,    (size_t)ncap*sizeof(int*));
  c->cap = ncap;
}

// printvridge-style callback (see printvridgeT in io_r.h). Called by
// qh_printvdiagram2 once per Voronoi ridge of the cell of site 0. vertexA is
// the neighboring input site; centers is the ordered set of Voronoi vertices
// (Delaunay facet circumcenters) bounding the ridge. The context pointer is
// smuggled in through the fp slot (qhull only forwards fp to this callback).
static void voronoi_collect(qhT *qh, FILE *fp, vertexT *vertex, vertexT *vertexA,
                            setT *centers, boolT unbounded){
  voronoi_ctx *c = (voronoi_ctx*) fp;
  facetT *facet, **facetp;
  int m, j;

  (void) vertex;
  if (unbounded){
    c->err = 1;
    return;
  }
  voronoi_ctx_grow(c);
  m = qh_setsize(qh, centers);
  c->ivws[c->nf]  = qh_pointid(qh, vertexA->point);
  c->nfvws[c->nf] = m;
  if (m > c->mnfv)
    c->mnfv = m;
  c->fv[c->nf] = (int*) malloc((size_t)(m > 0 ? m : 1)*sizeof(int));
  j = 0;
  // 1-based vertex id (visitid); xvws is stored so that Fortran's xvws(:,id)
  // selects the center numbered "id".
  FOREACHsetelement_(facetT, centers, facet)
    c->fv[c->nf][j++] = (int)facet->visitid;
  c->nf++;
}

// From a list of n vertices (xstar), calculate the Voronoi polyhedron of the
// origin (added as input site 0) and return its number of faces (nf), number
// of vertices (nv), maximum number of vertices per face (mnfv), and a context
// handle (ctx) holding the full result for step2. *ier is nonzero on failure.
void runqhull_voronoi_step1(int n, double xstar[n][3], int *nf, int *nv, int *mnfv,
                            void **ctx, int *ier){
  qhT qh_qh;
  qhT *qh = &qh_qh;
  coordT *points;
  voronoi_ctx *c;
  facetT *facet;
  unsigned int numfacets;
  int exitcode, curlong, totlong, i, idx, maxid;
  int numcenters;
  boolT isLower;
  setT *vertices;

  *ier = 0; *nf = 0; *nv = 0; *mnfv = 0; *ctx = NULL;

  // input points: the origin (site 0) followed by the n star vectors
  points = (coordT*) malloc((size_t)3*(n+1)*sizeof(coordT));
  points[0] = 0.; points[1] = 0.; points[2] = 0.;
  for (i = 0; i < n; i++){
    points[3*(i+1)+0] = xstar[i][0];
    points[3*(i+1)+1] = xstar[i][1];
    points[3*(i+1)+2] = xstar[i][2];
  }

  c = (voronoi_ctx*) calloc(1, sizeof(voronoi_ctx));

  // v + Qbb: Voronoi; QV0: cell of the first point (origin); Qs: search for
  // the initial simplex. (Fv/p only set internal output state; nothing is
  // printed because no output file is given.)
  qh_zero(qh, stderr);
  exitcode = qh_new_qhull(qh, 3, n+1, points, False, "qhull v Qbb QV0 Fv p Qs", NULL, stderr);
  if (!exitcode){
    exitcode = setjmp(qh->errexit);
    if (!exitcode){
      qh->NOerrexit = False;

      // number the Voronoi vertices (facet->visitid) exactly as qhull's text
      // output would, then compute their coordinates (facet->center).
      vertices = qh_markvoronoi(qh, qh->facet_list, NULL, !qh_ALL, &isLower, &numcenters);

      // walk the Voronoi ridges of site 0, collecting faces into the context.
      // Centers (facet->center) are computed afterwards, matching qhull's text
      // path (Fv is emitted before the 'p' format triggers qh_setvoronoi_all),
      // so qh_detvridge3 sees the same center state and produces the same order.
      qh_printvdiagram2(qh, (FILE*)c, voronoi_collect, vertices, qh_RIDGEall, True /*inorder*/);
      qh_settempfree(qh, &vertices);
      qh_setvoronoi_all(qh);

      // gather the Voronoi vertex coordinates by id (visitid-1)
      numfacets = (unsigned int) qh->num_facets;
      maxid = 0;
      FORALLfacets {
        if (facet->visitid && facet->visitid < numfacets && (int)facet->visitid > maxid)
          maxid = (int)facet->visitid;
      }
      c->nv = maxid;
      // calloc: any center id in [1,maxid] whose facet->center is unexpectedly
      // NULL leaves a deterministic 0 rather than uninitialized garbage.
      c->xvws = (double*) calloc((size_t)3*(maxid > 0 ? maxid : 1), sizeof(double));
      FORALLfacets {
        if (facet->visitid && facet->visitid < numfacets && facet->center){
          idx = (int)facet->visitid - 1;
          c->xvws[3*idx+0] = facet->center[0];
          c->xvws[3*idx+1] = facet->center[1];
          c->xvws[3*idx+2] = facet->center[2];
        }
      }
      if (c->err)
        exitcode = 1;
    }
    qh->NOerrexit = True;
  }
  *ier = exitcode;

  qh_freeqhull(qh, !qh_ALL);
  qh_memfreeshort(qh, &curlong, &totlong);
  free(points);

  *nf = c->nf;
  *nv = c->nv;
  *mnfv = c->mnfv;
  *ctx = (void*) c;
}

// Copy the Voronoi polyhedron collected in step1 into the Fortran arrays and
// free the context. Returns: the neighbor site id for each face (ivws), the
// vertices of the polyhedron (xvws), the number of vertices per face (nfvws),
// and the list of vertices for each face (fvws). The memory layout matches the
// Fortran declarations: xvws(3,nv) and fvws(mnfv,nf), iside on the Fortran side.
void runqhull_voronoi_step2(int nf, int nv, int mnfv, int ivws[nf], double xvws[nv][3],
                            int nfvws[nf], int fvws[nf][mnfv], void *ctx){
  voronoi_ctx *c = (voronoi_ctx*) ctx;
  int i, j;

  for (i = 0; i < nf; i++){
    ivws[i] = c->ivws[i];
    nfvws[i] = c->nfvws[i];
    for (j = 0; j < c->nfvws[i]; j++)
      fvws[i][j] = c->fv[i][j];
    for (j = c->nfvws[i]; j < mnfv; j++)
      fvws[i][j] = 0;
  }
  for (i = 0; i < nv; i++){
    xvws[i][0] = c->xvws[3*i+0];
    xvws[i][1] = c->xvws[3*i+1];
    xvws[i][2] = c->xvws[3*i+2];
  }

  for (i = 0; i < c->nf; i++)
    free(c->fv[i]);
  free(c->fv);
  free(c->ivws);
  free(c->nfvws);
  free(c->xvws);
  free(c);
}

/*=====================================================================
  Convex-hull triangulation of a basin surface
  =====================================================================*/

typedef struct {
  int nf;       // number of triangular faces
  int cap;      // capacity of the face array
  int (*f)[3];  // [nf][3] 1-based vertex indices into the original xvert list
} tri_ctx;

static void tri_ctx_grow(tri_ctx *c){
  if (c->nf < c->cap)
    return;
  int ncap = c->cap ? 2*c->cap : 64;
  c->f = (int (*)[3]) realloc(c->f, (size_t)ncap*sizeof(*c->f));
  c->cap = ncap;
}

// From a list of n vertices (xvert) and an inside point (x0), calculate the
// convex hull of the vertices by first projecting them onto the unit sphere
// (centered on x0), then running qhull. Returns the number of triangular faces
// (nf) and a context handle (ctx) holding their vertex indices for step2.
void runqhull_basintriangulate_step1(int n, double x0[3], double xvert[n][3],
                                     int *nf, void **ctx, int *ier){
  qhT qh_qh;
  qhT *qh = &qh_qh;
  const double thrsnorm = 1e-40;
  coordT *points;
  int *origidx;
  tri_ctx *c;
  facetT *facet;
  vertexT *vertex, **vertexp;
  int exitcode, curlong, totlong, i, nactual, k, v[3];

  *ier = 0; *nf = 0; *ctx = NULL;

  // center on x0, drop points coincident with the center, project the rest
  // onto the unit sphere, and remember the original index of each kept point.
  points = (coordT*) malloc((size_t)3*(n > 0 ? n : 1)*sizeof(coordT));
  origidx = (int*) malloc((size_t)(n > 0 ? n : 1)*sizeof(int));
  nactual = 0;
  for (i = 0; i < n; i++){
    double x[3] = {xvert[i][0]-x0[0], xvert[i][1]-x0[1], xvert[i][2]-x0[2]};
    double norm = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    if (norm > thrsnorm){
      points[3*nactual+0] = x[0]/norm;
      points[3*nactual+1] = x[1]/norm;
      points[3*nactual+2] = x[2]/norm;
      origidx[nactual] = i;
      nactual++;
    }
  }

  c = (tri_ctx*) calloc(1, sizeof(tri_ctx));

  // (no options = convex hull) Qt: triangulated output; Pp: no precision
  // warnings; QJ: joggle the input to avoid precision problems.
  qh_zero(qh, stderr);
  exitcode = qh_new_qhull(qh, 3, nactual, points, False, "qhull Qt Pp QJ", NULL, stderr);
  if (!exitcode){
    exitcode = setjmp(qh->errexit);
    if (!exitcode){
      qh->NOerrexit = False;
      // each hull facet is a triangle (Qt); read its three vertices, mapping
      // qhull point ids back to the original xvert indices (1-based).
      FORALLfacets {
        if (facet->upperdelaunay || !facet->normal)
          continue;
        k = 0;
        FOREACHvertex_(facet->vertices){
          if (k < 3)
            v[k] = origidx[qh_pointid(qh, vertex->point)] + 1;
          k++;
        }
        if (k == 3){
          tri_ctx_grow(c);
          c->f[c->nf][0] = v[0];
          c->f[c->nf][1] = v[1];
          c->f[c->nf][2] = v[2];
          c->nf++;
        }
      }
    }
    qh->NOerrexit = True;
  }
  *ier = exitcode;

  qh_freeqhull(qh, !qh_ALL);
  qh_memfreeshort(qh, &curlong, &totlong);
  free(points);
  free(origidx);

  *nf = c->nf;
  *ctx = (void*) c;
}

// Copy the convex-hull faces collected in step1 into the Fortran array and
// free the context. Returns the list of (1-based) vertex indices per face.
void runqhull_basintriangulate_step2(int nf, int iface[nf][3], void *ctx){
  tri_ctx *c = (tri_ctx*) ctx;
  int i;

  for (i = 0; i < nf; i++){
    iface[i][0] = c->f[i][0];
    iface[i][1] = c->f[i][1];
    iface[i][2] = c->f[i][2];
  }

  free(c->f);
  free(c);
}
