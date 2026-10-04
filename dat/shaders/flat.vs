#version 330 core

// Planar (screen-space) shapes: triangles pre-tessellated on the CPU,
// given directly in the NDC of the render buffer, with a per-shape depth.

layout (location = 0) in vec3 inPosition; // NDC position (x,y) and depth (z)
layout (location = 1) in vec4 inColor;    // color and opacity

out vec4 fColor;

void main(){
  fColor = inColor;
  gl_Position = vec4(inPosition, 1.0);
}
