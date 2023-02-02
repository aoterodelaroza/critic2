#version 330 core
uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

layout (location = 0) in vec3 inPosition;

out vec3 center;
out vec3 vertex;

void main(){
  mat4 vwm = view * world * model;
  vec4 center4 = vwm * vec4(0.,0.,0.,1.);
  vec4 vertex4 = vwm * vec4(inPosition, 1.0);
  center = center4.xyz / center4.w;
  vertex = vertex4.xyz / vertex4.w;
  gl_Position = projection * vertex4;
}
