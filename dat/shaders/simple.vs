#version 330 core
uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

layout (location = 0) in vec3 inPosition;
layout (location = 1) in vec3 inNormal;

void main(){
  gl_Position = projection * view * world * model * vec4(inPosition, 1.0);
}
