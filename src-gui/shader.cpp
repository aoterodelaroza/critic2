/*
Copyright (c) 2017 Alberto Otero de la Roza
<aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
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
#include <string.h>
#include "shader.h"

GLuint lightshader;
GLuint gWorldLocation;
GLuint gWVPLocation;
GLuint vColorLocation;
GLuint lColorLocation;
GLuint lDirectionLocation;
GLuint fAmbientIntensityLocation;

// Add a shader to the gl program
void AddShader(GLuint ShaderProgram, const char* pShaderText, GLenum ShaderType)
{
  GLuint ShaderObj = glCreateShader(ShaderType);

  if (ShaderObj == 0) {
    fprintf(stderr, "Error creating shader type %d\n", ShaderType);
    exit(1);
  }

  const GLchar * p[1];
  p[0] = pShaderText;
  GLint Lengths[1];
  Lengths[0] = strlen(pShaderText);
  glShaderSource(ShaderObj, 1, p, Lengths);
  glCompileShader(ShaderObj);
  GLint success;
  glGetShaderiv(ShaderObj, GL_COMPILE_STATUS, &success);
  if (!success) {
    GLchar InfoLog[1024];
    glGetShaderInfoLog(ShaderObj, 1024, NULL, InfoLog);
    fprintf(stderr, "Error compiling shader type %d: '%s'\n", ShaderType, InfoLog);
    exit(1);
  }
  glAttachShader(ShaderProgram, ShaderObj);
}

// Create a shader program based on a shader script
GLuint LightingShader()
{

  GLuint ShaderProgram = glCreateProgram();
  if (ShaderProgram == 0){
    exit(1);
  }

  const char * vs = "#version 330 \n \
    uniform mat4 gWorld; \n \
    uniform mat4 gWVP; \n \
    layout (location = 0) in vec3 inPosition; \n \
    layout (location = 1) in vec3 inNormal; \n \
    smooth out vec3 vNormal; \n \
    void main() { \n \
      gl_Position = gWVP * vec4(inPosition, 1.0); \n \
      vNormal = (gWorld * vec4(inNormal, 0.0)).xyz; \n \
      }";

  const char * fs = "#version 330 \n \
    smooth in vec3 vNormal; \n \
    uniform vec4 vColor; \n \
    out vec4 outputColor; \n \
    uniform vec3 lColor; \n \
    uniform vec3 lDirection; \n \
    uniform float fAmbientIntensity; \n \
    void main() { \n \
      float fDiffuseIntensity = max(0.0, dot(normalize(vNormal), lDirection)); \n \
      outputColor = vColor; \n \
      }";

  AddShader(ShaderProgram, vs, GL_VERTEX_SHADER);
  AddShader(ShaderProgram, fs, GL_FRAGMENT_SHADER);

  GLint success = 0;

  glLinkProgram(ShaderProgram);
  glGetProgramiv(ShaderProgram, GL_LINK_STATUS, &success);
  if (success == 0) exit(1);

  glValidateProgram(ShaderProgram);
  glGetProgramiv(ShaderProgram, GL_VALIDATE_STATUS, &success);
  if (success == 0) exit(1);

  return ShaderProgram;
}

