#version 330 core

// Flat-colored instanced mesh (planes, polyhedra faces, cone arrowheads).

flat in vec4 fColor;
out vec4 outColor;

void main(){
  outColor = fColor;
}
