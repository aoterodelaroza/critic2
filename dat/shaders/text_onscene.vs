#version 330 core
layout (location = 0) in vec3 x0;
layout (location = 1) in float rshift;
layout (location = 2) in vec2 vertex;
layout (location = 3) in vec2 tex;
out vec2 TexCoords;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;
uniform float depth;

void main(){
  vec4 x = view * world * vec4(x0, 1.0);
  vec4 x0 = projection * x;
  x /= x.w;
  x.z += rshift;
  x = projection * x;
  x /= x.w;
  x0 /= x0.w;

  gl_Position = vec4(x0.xy + vertex, x.z, 1.0);
  TexCoords = tex;
}
