// -*-c++-*-

// based on the camera class from learnopengl.com

#ifndef CAMERA_H
#define CAMERA_H

#include "imgui/gl3w.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "imgui/mouse.h"

using namespace glm;

// An abstract camera class that processes input and calculates the corresponding Eular Angles, Vectors and Matrices for use in OpenGL
class Camera
{
private:
  MouseState pmstate; // mouse state for the previous frame
  vec2 posdrag; // saved position for drag offset
public:
  vec3 Position; // camera position
  vec3 Up; // camera up direction (normal)
  vec3 LookAt; // camera look at position
  bool isortho = false; // true = orthographic; false = perspective
  float fov = 45.f; // field of view for perspective
  float orthoa = 10.f; // half the side of the ortho square
  float depth = 100.f; // depth of vision
  float mousesens_pan = 0.02; // Mouse pan sensitivity
  float mousesens_rot = 0.02; // Mouse rotate sensitivity
  float mousesens_zoom = 0.2; // Mouse zoom sensitivity

  Camera(vec3 position,vec3 up,vec3 lookat){
    Position = position;
    Up = up;
    LookAt = lookat;
  }

  Camera(float posX, float posY, float posZ, 
	 float upX, float upY, float upZ, 
	 float laX, float laY, float laZ){
    Position = vec3(posX, posY, posZ);
    Up = normalize(vec3(upX, upY, upZ));
    LookAt = vec3(laX, laY, laZ);
  }

  mat4 GetViewMatrix(){
    return lookAt(Position, LookAt, Up);
  }

  mat4 GetProjectionMatrix(){
    if (isortho)
      return ortho(-orthoa,orthoa,-orthoa,orthoa,0.1f,depth);
    else
      return perspective(radians(fov),1.0f,0.1f,depth);
  }

  void processMouseEvents(MouseState *m);
};
#endif
