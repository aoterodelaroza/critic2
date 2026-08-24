#version 330 core

// Indexed triangle mesh with per-vertex normals (isosurfaces). Vertices are
// in scene coordinates; the normal is carried to the fragment shader in view
// space for the headlight shading.

layout (location = 0) in vec3 inPosition; // vertex position
layout (location = 1) in vec3 inNormal;   // vertex normal

uniform mat4 view;
uniform mat4 world;
uniform mat4 projection;
uniform vec3 xshift; // offset of this copy of the mesh (scene coordinates)

out vec3 vNormal; // view-space normal

void main(){
  // view*world is a rotation plus a translation in this codebase, so the
  // normal transforms with its upper-left 3x3 block
  vNormal = mat3(view * world) * inNormal;
  gl_Position = projection * view * world * vec4(inPosition + xshift, 1.0);
}
