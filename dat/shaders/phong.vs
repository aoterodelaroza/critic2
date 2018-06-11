#version 330 core
uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;
uniform mat3 normrot;

layout (location = 0) in vec3 inPosition;
layout (location = 1) in vec3 inNormal;

out vec3 Normal;
out vec4 fragPos;

void main(){
  fragPos = view * world * model * vec4(inPosition, 1.0);
  gl_Position = projection * fragPos;
  Normal = normalize(normrot * inNormal);
}
