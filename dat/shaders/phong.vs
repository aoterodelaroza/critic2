#version 330 core
uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

layout (location = 0) in vec3 inPosition;

out vec3 Normal;
out vec4 fragPos;

void main(){
  vec3 nn = normalize(inPosition);
  mat4 vwm = view * world * model;
  fragPos = vwm * vec4(inPosition, 1.0);
  gl_Position = projection * fragPos;
  Normal = normalize(transpose(inverse(mat3(vwm))) * nn);
}
