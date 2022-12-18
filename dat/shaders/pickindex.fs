#version 330 core

uniform vec4 idx;
out vec4 outputColor;

void main(){
  outputColor = idx;
}
