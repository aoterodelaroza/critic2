#version 330 core

uniform int object_type; // 0=sphere, 1=cyl, 2=cyl_flat
uniform float rborder; // border of the object
uniform vec3 bordercolor; // color of the border
uniform vec4 vColor; // color of the object
uniform int isortho; // 1=orthographic, 0=perspective projection

in vec3 center;
in vec3 vertex;
in vec3 up;
in vec3 side;
out vec4 outputColor;

void main(){

  // vector to the center
  vec3 vx = vertex - center;

  // view direction of this fragment: -z for orthographic and the ray
  // from the camera to the fragment in perspective
  vec3 ncam = (isortho != 0) ? vec3(0.,0.,-1.) : normalize(vertex);

  if (object_type == 0){
    // sphere //

    // project onto the screen plane passing through the center of the sphere.
    // radius = radius of the sphere, rproj = distance to the center
    float radius = length(vx);
    float rproj = length(vx - dot(vx,ncam) * ncam);

    // use the border color if close to the edge of the sphere
    if (radius - rproj < rborder)
      outputColor = vec4(bordercolor,vColor.a);
    else
      outputColor = vColor;
  } else if (object_type == 1) {
    // cylinder //

    // radius of the cylinder in view coordinates
    float radius = length(side-center);

    // up vector = along the cylinder, ncam = view direction at this fragment
    // (constant -z for orthographic, ray to the fragment for perspective)
    vec3 upn = normalize(up-center);
    vec3 ncam = (isortho != 0) ? vec3(0.,0.,-1.) : normalize(vertex);

    // vperp is perpendicular to upn and ncam
    vec3 vperp = normalize(cross(ncam,upn));

    // project out the component along the cylinder axis
    vec3 vx = vertex-center;
    vx = vx - dot(vx,upn) * upn;

    // border
    if (radius - radius * abs(dot(normalize(vx),vperp)) < rborder)
      outputColor = vec4(bordercolor,vColor.a);
    else
      outputColor = vColor;
  } else if (object_type == 2) {
    // flat cylinder //
    outputColor = vColor;
  }
}
