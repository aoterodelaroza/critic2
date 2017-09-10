// -*-c++-*-

// based on the camera class from learnopengl.com

#ifndef CAMERA_H
#define CAMERA_H

#include "imgui/gl3w.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtx/transform.hpp>
#include "imgui/mouse.h"
#include "settings.h"
#include <stdio.h>

using namespace glm;

// An abstract camera class that processes input and calculates the corresponding Eular Angles, Vectors and Matrices for use in OpenGL
class Camera
{
private:
  mat4 view; // view matrix
  mat4 projection; // projection matrix
  MouseState pmstate; // mouse state for the previous frame
  vec2 posdrag; // saved position for drag offset
  mat4 viewdrag; // saved view matrix for drag offset
public:
  bool isortho = true; // true = orthographic; false = perspective
  float fov = 45.f; // field of view for perspective
  float zoom = 1.f; // zoom for orthographic
  float srad = 20.f; // scene radius

  Camera(vec3 position,vec3 lookat,vec3 up){
    view = lookAt(position, lookat, up);
    setProjection();
  }

  Camera(float posX, float posY, float posZ, 
	 float laX, float laY, float laZ,
	 float upX, float upY, float upZ){
    view = lookAt(vec3(posX, posY, posZ), vec3(laX, laY, laZ), 
		  vec3(upX, upY, upZ));
    setProjection();
  }

  void setProjection(){
    if (isortho)
      projection = ortho(-0.5f*srad*zoom,0.5f*srad*zoom,
			 -0.5f*srad*zoom,0.5f*srad*zoom,0.1f,100.f);
    else
      projection = infinitePerspective(radians(fov),1.0f,0.1f);
  }
  void processMouseEvents(MouseState *m);
  void applyMatrices(GLuint shader_id);
};
#endif
