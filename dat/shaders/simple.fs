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
    // cylinder calculations
    // half-height and radius of the cylinder
    float hheight = length(up-center);
    float radius = length(side-center);

    // up vector and vector pointing towards camera
    vec3 upn = normalize(up-center);
    vec3 n = vec3(0.,0.,-1.);

    // side2 is perpendicular to up and n
    vec3 side2 = normalize(cross(n,upn));

    // calculate the distance along the cylinder (-hheight to hheight) and project out the up
    vec3 vx = vertex-center;
    float ralong = dot(vx,upn);
    vx = vx - ralong * upn;

    // // discard to make dashed
    // int ndashes = 10;
    // if (mod(trunc(0.5 * (hheight + ralong) / hheight * ndashes),2) == 1) discard;

    // border
    if (radius * (1 - abs(dot(normalize(vx),side2))) < rborder_cyl)
      outputColor = vec4(bordercolor,vColor.a);
    else
      outputColor = vColor;
  } else {
    outputColor = vColor;
  }
}
