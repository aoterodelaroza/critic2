#version 330 core

uniform vec4 vColor;
uniform vec3 bordercolor;
uniform float rborder;

in vec3 center;
in vec3 vertex;
out vec4 outputColor;

void main(){
  if (rborder > 0){
    // radius and projected radius calculations in view space
    vec3 vx = vertex - center;
    vec3 n = vec3(0.,0.,-1.);
    float radius = length(vx);
    float rproj = length(vx - dot(vx,n) * n);

    if (radius - rproj < rborder)
      outputColor = vec4(bordercolor,vColor.a);
    else
      outputColor = vColor;
  } else {
    outputColor = vColor;
  }
}
