// Copyright (c) 2017 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
// Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
// <victor@fluor.quimica.uniovi.es>.
//
// critic2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// critic2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "shapes.h"
#include "settings.h"
#include <cmath>
#include <cstring>
#include <map>

# define PI 3.14159265358979323846

using namespace std;

// Some constants
static float tau = (1.0 + sqrt(5))/2.0;
static float rad0 = sqrt(3 - tau);
static float xico = (tau-1)/rad0;
static float zico = 1/rad0;

static float cylr = cosf(30.f * PI / 180.f);
static float cylnorm = sqrtf(0.5f + cylr * cylr);
static float cylrn = cylr / cylnorm;
static float halfn = 0.5f / cylnorm;
static float halfn2 = 2 * halfn;

// global variables
GLuint sphVAO[nmaxsph];
GLuint sphVBO;
GLuint sphEBO[nmaxsph];
const GLuint sphnve[nmaxsph]    = {     12,  42, 162,  642};
const GLuint sphnel[nmaxsph]    = {     20,  80, 320, 1280};
const GLuint sphneladd[nmaxsph+1] = {0, 20, 100, 420, 1700};
GLuint cylVAO[nmaxcyl];
GLuint cylVBO;
GLuint cylEBO[nmaxcyl];
const GLuint cylnve[nmaxcyl]    = {     14};
const GLuint cylnel[nmaxcyl]    = {     24};
const GLuint cylneladd[nmaxcyl+1] = {0, 24};

// icosahedron  vertices
static GLfloat *sphv;
static const GLfloat sphv0[] = {
  -xico,   0.0,  zico,  -xico,   0.0,  zico,
   xico,   0.0,  zico,   xico,   0.0,  zico,
  -xico,   0.0, -zico,  -xico,   0.0, -zico,
   xico,   0.0, -zico,   xico,   0.0, -zico,
    0.0,  zico,  xico,    0.0,  zico,  xico,
    0.0,  zico, -xico,    0.0,  zico, -xico,
    0.0, -zico,  xico,    0.0, -zico,  xico,
    0.0, -zico, -xico,    0.0, -zico, -xico,
   zico,  xico,   0.0,   zico,  xico,   0.0,
  -zico,  xico,   0.0,  -zico,  xico,   0.0,
   zico, -xico,   0.0,   zico, -xico,   0.0,
  -zico, -xico,   0.0,  -zico, -xico,   0.0,
};

// icosehedron faces
static GLuint *sphi;
static GLuint sphi0[] = {
  0, 1,  4,
  0, 4,  9,
  9, 4,  5,
  4, 8,  5,
  4, 1,  8,
  8, 1, 10,
  8, 10, 3,
  5, 8,  3,
  5, 3,  2,
  2, 3,  7,
  7, 3, 10,
  7, 10, 6,
  7, 6, 11,
 11, 6,  0,
  0, 6,  1,
  6, 10, 1,
  9, 11, 0,
  9, 2, 11,
  9, 5,  2,
  7, 11, 2,
};

// cylinder vertices (hexagonal prism)
static GLfloat *cylv;
static GLfloat cylv0[] = {
    0.0,  0.0, -0.5,     0.0,     0.0,   -1.0,
   cylr,  0.5, -0.5,   cylrn,   halfn, -halfn,
   cylr, -0.5, -0.5,   cylrn,  -halfn, -halfn,
    0.0, -1.0, -0.5,     0.0, -halfn2, -halfn,
  -cylr, -0.5, -0.5,  -cylrn,  -halfn, -halfn,
  -cylr,  0.5, -0.5,  -cylrn,   halfn, -halfn,
    0.0,  1.0, -0.5,     0.0,  halfn2, -halfn,
    0.0,  0.0,  0.5,     0.0,     0.0,    1.0,
   cylr,  0.5,  0.5,   cylrn,   halfn,  halfn,
   cylr, -0.5,  0.5,   cylrn,  -halfn,  halfn,
    0.0, -1.0,  0.5,     0.0, -halfn2,  halfn,
  -cylr, -0.5,  0.5,  -cylrn,  -halfn,  halfn,
  -cylr,  0.5,  0.5,  -cylrn,   halfn,  halfn,
    0.0,  1.0,  0.5,     0.0,  halfn2,  halfn,
};

// hexagonal prism faces
static GLuint *cyli;
static GLuint cyli0[] = {
  0,  1,   2,
  0,  2,   3,
  0,  3,   4,
  0,  4,   5,
  0,  5,   6,
  0,  6,   1,
  7,  9,   8,
  7, 10,   9,
  7, 11,  10,
  7, 12,  11,
  7, 13,  12,
  7,  8,  13, 
  1,  8,   2, 
  2,  8,   9,
  2,  9,   3,
  3,  9,  10,
  3, 10,   4,
  4, 10,  11,
  4, 11,   5,
  5, 11,  12,
  5, 12,   6,
  6, 12,  13,
  6, 13,   1,
  1, 13,   8,
};

// create the buffers for the sphere and cylinder objects
void CreateAndFillBuffers(){
  // allocate space for the icospheres and copy the icosahedron as the first icosphere
  sphv = new GLfloat [6*sphnve[nmaxsph-1]];
  sphi = new GLuint [3*sphneladd[nmaxsph]];
  memcpy(sphv,sphv0,sizeof(sphv0));
  memcpy(sphi,sphi0,sizeof(sphi0));

  // recursive subdivision for the other icospheres
  for (int i = 1; i < nmaxsph; i++){
    map<pair<int,int>,int> edgev = {};

    int n = sphnve[i-1] - 1;
    int nface = sphneladd[i]-1;
    for (int j = sphneladd[i-1]; j < sphneladd[i]; j++){
      int k1 = sphi[3*j+0];
      int k2 = sphi[3*j+1];
      int k3 = sphi[3*j+2];
      int nk1, nk2, nk3;

      // create the new vertices
      if (edgev[make_pair(k1,k2)]){
        nk1 = edgev[make_pair(k1,k2)];
      } else {
        edgev[make_pair(k1,k2)] = edgev[make_pair(k2,k1)] = nk1 = ++n;
        sphv[6*n+0] = sphv[6*n+3] = 0.5f * (sphv[6*k1+0] + sphv[6*k2+0]);
        sphv[6*n+1] = sphv[6*n+4] = 0.5f * (sphv[6*k1+1] + sphv[6*k2+1]);
        sphv[6*n+2] = sphv[6*n+5] = 0.5f * (sphv[6*k1+2] + sphv[6*k2+2]);
      }
      if (edgev[make_pair(k1,k3)]){
        nk2 = edgev[make_pair(k1,k3)];
      } else {
        edgev[make_pair(k1,k3)] = edgev[make_pair(k3,k1)] = nk2 = ++n;
        sphv[6*n+0] = sphv[6*n+3] = 0.5f * (sphv[6*k1+0] + sphv[6*k3+0]);
        sphv[6*n+1] = sphv[6*n+4] = 0.5f * (sphv[6*k1+1] + sphv[6*k3+1]);
        sphv[6*n+2] = sphv[6*n+5] = 0.5f * (sphv[6*k1+2] + sphv[6*k3+2]);
      }
      if (edgev[make_pair(k2,k3)]){
        nk3 = edgev[make_pair(k2,k3)];
      } else {
        edgev[make_pair(k2,k3)] = edgev[make_pair(k3,k2)] = nk3 = ++n;
        sphv[6*n+0] = sphv[6*n+3] = 0.5f * (sphv[6*k2+0] + sphv[6*k3+0]);
        sphv[6*n+1] = sphv[6*n+4] = 0.5f * (sphv[6*k2+1] + sphv[6*k3+1]);
        sphv[6*n+2] = sphv[6*n+5] = 0.5f * (sphv[6*k2+2] + sphv[6*k3+2]);
      }

      // create the new faces
      nface++;
      sphi[3*nface+0] = k1;  sphi[3*nface+1] = nk1; sphi[3*nface+2] = nk2; 
      nface++;
      sphi[3*nface+0] = nk1; sphi[3*nface+1] = nk3; sphi[3*nface+2] = nk2; 
      nface++;
      sphi[3*nface+0] = nk1; sphi[3*nface+1] = k2;  sphi[3*nface+2] = nk3; 
      nface++;
      sphi[3*nface+0] = nk2; sphi[3*nface+1] = nk3; sphi[3*nface+2] = k3; 
    }

    // normalize the vertices
    for (int j = sphnve[i-1]; j < sphnve[i]; j++){
      float anorm = sqrtf(sphv[6*j+0]*sphv[6*j+0]+sphv[6*j+1]*sphv[6*j+1]+sphv[6*j+2]*sphv[6*j+2]);
      sphv[6*j+0] = sphv[6*j+3] = sphv[6*j+0] / anorm;
      sphv[6*j+1] = sphv[6*j+4] = sphv[6*j+1] / anorm;
      sphv[6*j+2] = sphv[6*j+5] = sphv[6*j+2] / anorm;
    }
  }

  // build the buffers for the spheres
  glGenVertexArrays(nmaxsph, sphVAO);
  glGenBuffers(1, &sphVBO);
  glGenBuffers(nmaxsph, sphEBO);

  glBindBuffer(GL_ARRAY_BUFFER, sphVBO);
  glBufferData(GL_ARRAY_BUFFER, 6 * sphnve[nmaxsph-1] * sizeof(GLfloat), sphv, GL_STATIC_DRAW);

  for (int i=0;i<nmaxsph;i++){
    glBindVertexArray(sphVAO[i]);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphEBO[i]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * sphnel[i] * sizeof(GLuint), &(sphi[3*sphneladd[i]]), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);
  }

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  // allocate space for the icospheres and copy the icosahedron as the first icosphere
  cylv = new GLfloat [6*cylnve[nmaxcyl-1]];
  cyli = new GLuint [3*cylneladd[nmaxcyl]];
  memcpy(cylv,cylv0,sizeof(cylv0));
  memcpy(cyli,cyli0,sizeof(cyli0));

  // normalize cylinder vertices
  for (int j = 0; j < cylnve[0]; j++){
    float anorm = sqrtf(cylv[6*j+0]*cylv[6*j+0]+cylv[6*j+1]*cylv[6*j+1]+cylv[6*j+2]*cylv[6*j+2]);
    cylv[6*j+3] = cylv[6*j+0] / anorm;
    cylv[6*j+4] = cylv[6*j+1] / anorm;
    cylv[6*j+5] = cylv[6*j+2] / anorm;
  }

  // build the buffers for the cylinders
  glGenVertexArrays(nmaxcyl, cylVAO);
  glGenBuffers(1, &cylVBO);
  glGenBuffers(nmaxcyl, cylEBO);

  glBindBuffer(GL_ARRAY_BUFFER, cylVBO);
  glBufferData(GL_ARRAY_BUFFER, 6 * cylnve[nmaxcyl-1] * sizeof(GLfloat), cylv, GL_STATIC_DRAW);

  for (int i=0;i<nmaxcyl;i++){
    glBindVertexArray(cylVAO[i]);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cylEBO[i]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * cylnel[i] * sizeof(GLuint), &(cyli[3*cylneladd[i]]), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);
  }

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  
}

// delete the buffers for the sphere and cylinder objects
void DeleteBuffers(){
  glDeleteVertexArrays(nmaxsph, sphVAO);
  glDeleteBuffers(1, &sphVBO);
  glDeleteBuffers(nmaxsph, sphEBO);
  glDeleteVertexArrays(nmaxcyl, cylVAO);
  glDeleteBuffers(1, &cylVBO);
  glDeleteBuffers(nmaxcyl, cylEBO);
  delete sphi;
  delete sphv;
  delete cyli;
  delete cylv;
}
