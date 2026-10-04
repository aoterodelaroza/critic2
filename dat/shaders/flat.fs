#version 330 core

// Planar (screen-space) shapes: per-vertex flat color.

in vec4 fColor;
out vec4 outColor;

void main(){
  outColor = fColor;
}
