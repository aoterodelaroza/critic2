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
uniform int isortho; // 1=orthographic, 0=perspective projection
uniform int isgizmo; // 1=screen-fixed overlay (window-anchored gizmo)
uniform vec3 gizmo_ndc; // target NDC position for the gizmo (xy used)
uniform float gizmo_scale; // per-item zoom-compensation factor (gizf)

// coordinates of the vertex
layout (location = 0) in vec3 inPosition;

// information passed to the fragment shader
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
      // apply a delta in a direction perpendicular to the camera
      // +delta corresponds to the ring exterior
      vec3 upn = normalize(up-center);
      vec3 ncam = (isortho != 0) ? vec3(0.,0.,-1.) : normalize(center);
      vec3 vperp = normalize(cross(ncam,upn));

      float dd = delta_cyl;
      if (dot(bond_outward,bond_outward) > 1e-10) {
        vec3 outv = mat3(view*world) * bond_outward;
        if (dot(outv,vperp) < 0.) dd = -dd;
      }
      vertex4 = vec4(vertex + dd * vperp,1.0);
    }
  } else {
    // flat cylinder
  }
  if (isgizmo != 0) {
    // screen-fixed overlay: anchor the gizmo geometry (built around the local
    // origin) at gizmo_ndc with an orthographic projection scaled by
    // gizmo_scale. The local point includes the model placement (model*p); only
    // the rotation of view*world is applied (no eye translation), so each item
    // keeps its offset from the gizmo origin. The view-space outputs
    // (center/vertex/up/side) above are unchanged; the fragment shader only
    // uses their differences and ncam, so lighting/border are unaffected.
    vec3 local = (model * vec4(inPosition,1.0)).xyz;
    vec4 rc = projection * vec4(mat3(view * world) * local * gizmo_scale, 0.0);
    gl_Position = vec4(gizmo_ndc.xy + rc.xy, rc.z, 1.0);
  } else {
    gl_Position = projection * vertex4;
  }
}
