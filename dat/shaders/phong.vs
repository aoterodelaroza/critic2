#version 330 core
uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

layout (location = 0) in vec3 inPosition;
layout (location = 1) in vec3 inNormal;

out vec3 Normal;
out vec4 fragPos;

void main(){
  mat3 normrot = transpose(inverse(mat3(view) * mat3(world) * mat3(model)));
  fragPos = view * world * model * vec4(inPosition, 1.0);
  gl_Position = projection * fragPos;
  Normal = normalize(normrot * inNormal);
}
