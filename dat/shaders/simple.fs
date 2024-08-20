#version 330 core

uniform int object_type; // 0=sphere, 1=cyl, 2=cyl_flat
uniform float rborder; // border of the object
uniform vec3 bordercolor; // color of the border
uniform vec4 vColor; // color of the object
uniform int ndash_cyl; // number of dashes in the cylinder

in vec3 center;
in vec3 vertex;
in vec3 up;
in vec3 side;
out vec4 outputColor;

void main(){

  if (object_type == 0){
    // sphere //

    // project onto the screen plane passing through the center of the sphere
    vec3 vx = vertex - center;
    vec3 ncam = vec3(0.,0.,-1.);

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

    // calculate half-height and radius of the cylinder in view coordinates
    float hheight = length(up-center);
    float radius = length(side-center);

    // up vector = along the cylinder, ncam = towards the camera
    vec3 upn = normalize(up-center);
    vec3 ncam = vec3(0.,0.,-1.);

    // vperp is perpendicular to upn and ncam
    vec3 vperp = normalize(cross(ncam,upn));

    // calculate the distance along the cylinder (-hheight to hheight) and project out the up
    vec3 vx = vertex-center;
    float ralong = dot(vx,upn);
    vx = vx - ralong * upn;

    // discard to make dashed
    if (ndash_cyl > 0 && mod(trunc(0.5 * (hheight + ralong) / hheight * ndash_cyl),2) == 1){
      discard;
    } else {
      // border
      if (radius - radius * abs(dot(normalize(vx),vperp)) < rborder)
	outputColor = vec4(bordercolor,vColor.a);
      else
	outputColor = vColor;
    }
  } else if (object_type == 2) {
    // flat cylinder //
    outputColor = vColor;
  }
}
