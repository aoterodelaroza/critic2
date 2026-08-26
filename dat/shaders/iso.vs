#version 330 core

// Indexed triangle mesh with per-vertex normals and colors (isosurfaces).
// Vertices are in scene coordinates; the normal is carried to the fragment
// shader in view space for the headlight shading. The color is per vertex so
// an isosurface can be colored by the values of another field; a flat-colored
// mesh simply carries the same color on every vertex.

layout (location = 0) in vec3 inPosition; // vertex position
layout (location = 1) in vec3 inNormal;   // vertex normal
layout (location = 2) in vec3 inColor;    // vertex color (flat meshes repeat one color)

uniform mat4 view;
uniform mat4 world;
uniform mat4 projection;
uniform vec3 xshift; // offset of this copy of the mesh (scene coordinates)

out vec3 vNormal; // view-space normal
out vec3 vColor;  // vertex color, interpolated across the triangle

void main(){
  // view*world is a rotation plus a translation in this codebase, so the
  // normal transforms with its upper-left 3x3 block
  vNormal = mat3(view * world) * inNormal;
  vColor = inColor;
  gl_Position = projection * view * world * vec4(inPosition + xshift, 1.0);
}
