#version 330 core

uniform vec4 vColor;
uniform vec3 bordercolor;
uniform float rborder_sph;
uniform float rborder_cyl;

in vec3 center;
in vec3 vertex;
in vec3 up;
in vec3 side;
out vec4 outputColor;

void main(){

  if (rborder_sph > 0){
    // sphere: radius and projected radius calculations in view space
    vec3 vx = vertex - center;
    vec3 n = vec3(0.,0.,-1.);
    float radius = length(vx);
    float rproj = length(vx - dot(vx,n) * n);

    if (radius - rproj < rborder_sph)
      outputColor = vec4(bordercolor,vColor.a);
    else
      outputColor = vColor;
  } else if (rborder_cyl > 0){
    // sphere: projection calculations in view space
    vec3 upn = normalize(up-center);
    vec3 n = vec3(0.,0.,-1.);
    float radius = length(side-center);
    vec3 side2 = normalize(cross(n,upn));
    vec3 vx = vertex-center;
    vx = vx - dot(vx,upn) * upn;
    if (radius * (1 - abs(dot(normalize(vx),side2))) < rborder_cyl)
      outputColor = vec4(bordercolor,vColor.a);
    else
      outputColor = vColor;
  } else {
    outputColor = vColor;
  }
}
