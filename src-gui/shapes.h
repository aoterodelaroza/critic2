// -*-c++-*-
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

#ifndef SHAPES_H
#define SHAPES_H

// public vertex buffer and attribute objects
GLuint FBO;
GLuint FBOtex;
GLuint sphereVAO[1], sphereVBO[1], sphereEBO[1];
const GLuint spherenv[1] = {36};
GLuint cylVAO[1], cylVBO[1], cylEBO[1];
const GLuint cylnv[1] = {42};

// Texture size
const float FBO_tex_x = 1024.f;
const float FBO_tex_y = 1024.f;

// Some constants
static float tau = (1.0 + sqrt(5))/2.0;
static float rad0 = sqrt(3 - tau);
static float xico = (tau-1)/rad0;
static float zico = 1/rad0;

// icosahedron  vertices
static GLfloat icov[] = {
  -xico,   0.0,  zico,
   xico,   0.0,  zico,
  -xico,   0.0, -zico,
   xico,   0.0, -zico,
    0.0,  zico,  xico,
    0.0,  zico, -xico,
    0.0, -zico,  xico,
    0.0, -zico, -xico,
   zico,  xico,   0.0,
  -zico,  xico,   0.0,
   zico, -xico,   0.0,
  -zico, -xico,   0.0
};

// icosehedron faces
static GLuint icoi[] = {
  1,  4,  0,
  4,  9,  0,
  4,  5,  9,
  8,  5,  4,
  1,  8,  4,
  1, 10,  8,
  10, 3,  8,
  8,  3,  5,
  3,  2,  5,
  3,  7,  2,
  3, 10,  7,
  10, 6,  7,
  6, 11,  7,
  6,  0, 11,
  6,  1,  0,
  10, 1,  6,
  11, 0,  9,
  2, 11,  9,
  5,  2,  9,
  11, 2,  7
};

// hexagonal prism vertices
static GLfloat hpv[] = {
  -0.86602540,  0.50000000,  0.00000000,
  -0.86602540, -0.50000000,  0.00000000,
  -0.00000000, -1.00000000,  0.00000000,
   0.86602540, -0.50000000,  0.00000000,
   0.86602540,  0.50000000,  0.00000000,
   0.00000000,  1.00000000,  0.00000000,
  -0.86602540,  0.50000000,  1.00000000,
  -0.86602540, -0.50000000,  1.00000000,
  -0.00000000, -1.00000000,  1.00000000,
   0.86602540, -0.50000000,  1.00000000,
   0.86602540,  0.50000000,  1.00000000,
   0.00000000,  1.00000000,  1.00000000,
   0.00000000,  0.00000000,  0.00000000,
   0.00000000,  0.00000000,  1.00000000,
};

// hexagonal prism faces
static GLuint hpi[] = {
  12,  0,  1,
  13,  6,  7,
   0,  1,  6,
   7,  6,  1,
  12,  1,  2,
  13,  7,  8,
   1,  2,  7,
   8,  7,  2,
  12,  2,  3,
  13,  8,  9,
   2,  3,  8,
   9,  8,  3,
  12,  3,  4,
  13,  9, 10,
   3,  4,  9,
  10,  9,  4,
  12,  4,  5,
  13, 10, 11,
   4,  5, 10,
  11, 10,  5,
  12,  5,  0,
  13, 11,  6,
   5,  0, 11,
   6, 11,  0
};

// create the buffers for the sphere and cylinder objects
void CreateAndFillBuffers(){
  // setup up OpenGL stuff
  glGenVertexArrays(1, sphereVAO);
  glGenBuffers(1, sphereVBO);
  glGenBuffers(1, sphereEBO);
  glGenVertexArrays(1, cylVAO);
  glGenBuffers(1, cylVBO);
  glGenBuffers(1, cylEBO);

  glBindVertexArray(sphereVAO[0]);

  glBindBuffer(GL_ARRAY_BUFFER, sphereVBO[0]);
  glBufferData(GL_ARRAY_BUFFER, sizeof(icov), icov, GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphereEBO[0]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(icoi), icoi, GL_STATIC_DRAW);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  glBindVertexArray(cylVAO[0]);

  glBindBuffer(GL_ARRAY_BUFFER, cylVBO[0]);
  glBufferData(GL_ARRAY_BUFFER, sizeof(icov), icov, GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cylEBO[0]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(icoi), icoi, GL_STATIC_DRAW);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  glGenFramebuffers(1, &FBO);
  glGenTextures(1, &FBOtex);

  glBindFramebuffer(GL_FRAMEBUFFER, FBO);
  glBindTexture(GL_TEXTURE_2D, FBOtex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, FBO_tex_x, FBO_tex_y, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);  
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, FBOtex, 0); 
  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    exit(EXIT_FAILURE);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

// delete the buffers for the sphere and cylinder objects
void DeleteBuffers(){
  glDeleteVertexArrays(1, sphereVAO);
  glDeleteBuffers(1, sphereVBO);
  glDeleteBuffers(1, sphereEBO);
  glDeleteTextures(1, &FBOtex);
  glDeleteFramebuffers(1, &FBO); 
}
#endif
