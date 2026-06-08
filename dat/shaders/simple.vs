#version 330 core
// transformation matrices
uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

uniform int object_type; // 0=sphere, 1=cyl, 2=cyl_flat
uniform float rborder; // border of the object
uniform vec3 bordercolor; // color of the border
uniform float delta_cyl; // displacement for the cylinder along screen plane
uniform vec3 bond_outward; // aromatic bonds: outward (ring-exterior) direction; 0 otherwise

// coordinates of the vertex
layout (location = 0) in vec3 inPosition;

// various positions in view coordinates passed to the fragment shader
out vec3 center;
out vec3 vertex;
out vec3 up;
out vec3 side;

void main(){
  mat4 vwm = view * world * model;

  vec4 vertex4 = vwm * vec4(inPosition, 1.0);

  if (object_type == 0) {
    // sphere -> pass vertex and center positions in view coordinates
    vertex = vertex4.xyz;
    center = (vwm * vec4(0.,0.,0.,1.)).xyz;
  } else if (object_type == 1) {
    // cylinder -> pass vertex, center, up, and side vectors
    vertex = vertex4.xyz;
    center = (vwm * vec4(0.,0.,0.,1.)).xyz;
    up = (vwm * vec4(0.,0.,0.5,1.)).xyz;
    side = (vwm * vec4(0.,0.5,0.,1.)).xyz;

    if (delta_cyl != 0.) {
      vec3 upn = normalize(up-center);
      vec3 ncam = vec3(0.,0.,-1.);
      vec3 vperp = normalize(cross(ncam,upn));

      float dd = delta_cyl;
      // aromatic bonds: orient the offset so +delta is on the ring-exterior
      // side (i.e. the inner, -delta sub-bond is the dashed one)
      if (dot(bond_outward,bond_outward) > 1e-10) {
        vec3 outv = mat3(view*world) * bond_outward;
        if (dot(outv,vperp) < 0.) dd = -dd;
      }
      vertex4 = vec4(vertex + dd * vperp,1.0);
    }
  } else {
    // flat cylinder
  }
  gl_Position = projection * vertex4;
}
