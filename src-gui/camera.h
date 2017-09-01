// -*-c++-*-

// based on the camera class from learnopengl.com

#ifndef CAMERA_H
#define CAMERA_H

#include "imgui/gl3w.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;

// An abstract camera class that processes input and calculates the corresponding Eular Angles, Vectors and Matrices for use in OpenGL
class Camera
{
public:
  vec3 Position = vec3(0.0f, 0.0f, 20.0f); // camera position
  vec3 Up = vec3(0.0f, 1.0f, 0.0f); // camera up direction
  vec3 LookAt = vec3(0.0f, 0.0f, 0.0f); // camera look at position
  bool isortho = true; // true = orthographic; false = perspective
  float fov = 45.f; // field of view for perspective
  float orthoa = 10.f; // half the side of the ortho square
  float depth = 100.f; // depth of vision

  Camera(){}

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
};
#endif
