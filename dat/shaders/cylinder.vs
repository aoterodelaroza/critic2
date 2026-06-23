#version 330 core

// Instanced capped-cylinder shader. Each instance is a billboard quad
// that covers the screen projection of the segment; the fragment
// shader ray-casts the actual cylinder. Endpoints arrive already
// animated (CPU side).

// per-vertex: quad corner in [-1,1]^2
layout (location = 0) in vec2 a_corner;

// per-instance attributes
layout (location = 1) in vec3 a_x1;          // endpoint 1 (pre-world coords)
layout (location = 2) in vec3 a_x2;          // endpoint 2 (pre-world coords)
layout (location = 3) in float a_radius;
layout (location = 4) in vec4 a_color;
layout (location = 5) in float a_border;
layout (location = 6) in vec3 a_bordercolor;
layout (location = 7) in float a_delta;      // lateral (screen-perp) offset for multi-bonds
layout (location = 8) in vec3 a_outward;     // aromatic ring-exterior direction (0 otherwise)

uniform mat4 view;
uniform mat4 world;
uniform mat4 projection;
uniform int isortho;       // 1=orthographic, 0=perspective
uniform int isanchored;       // 1=window-anchored overlay, 0=normal scene
uniform vec3 anchored_ndc;    // target NDC position for the anchor (xy used)
uniform float anchored_scale; // anchor zoom-compensation factor

out vec3 fA;               // eye-space endpoint 1
out vec3 fB;               // eye-space endpoint 2
out vec3 fPosEye;          // billboard point in eye space
out float fR;
flat out vec4 fColor;
flat out float fBorder;
flat out vec3 fBorderColor;

void main(){
  // eye-space endpoints. For a window-anchored object only the rotation of
  // view*world is applied (no eye translation) and the geometry is scaled to a
  // constant on-screen size.
  mat4 vm = view * world;
  vec3 ae, be;
  float rad = a_radius;
  if (isanchored != 0){
    ae = mat3(vm) * (a_x1 * anchored_scale);
    be = mat3(vm) * (a_x2 * anchored_scale);
    rad = a_radius * anchored_scale;
  } else {
    ae = (vm * vec4(a_x1, 1.0)).xyz;
    be = (vm * vec4(a_x2, 1.0)).xyz;
  }
  vec3 mid = 0.5 * (ae + be);

  // multi-bond lateral offset, perpendicular to the view direction and the
  // cylinder axis (matches the legacy delta_cyl behaviour). +delta is the
  // aromatic ring exterior when a_outward is set.
  if (a_delta != 0.0){
    vec3 axe = be - ae;
    float al = length(axe);
    if (al > 1e-7){
      vec3 dax = axe / al;
      vec3 ncam = (isortho != 0) ? vec3(0.0,0.0,-1.0) : normalize(mid);
      vec3 vperp = normalize(cross(ncam, dax));
      float dd = a_delta;
      if (dot(a_outward, a_outward) > 1e-10){
        vec3 outv = mat3(vm) * a_outward;
        if (dot(outv, vperp) < 0.0) dd = -dd;
      }
      ae += dd * vperp;
      be += dd * vperp;
      mid += dd * vperp;
    }
  }

  vec3 axis = be - ae;
  float len = max(length(axis), 1e-7);
  vec3 d = axis / len;

  // perspective enlargement so the quad covers the silhouette (misses are
  // discarded in the fragment shader)
  float dist = max(length(mid), rad + 1e-4);
  float fac = (isortho != 0) ? 1.0 : dist / sqrt(max(dist*dist - rad*rad, 1e-6));
  fac = min(fac, 4.0);
  float rr = rad * fac;

  // screen-perpendicular to the axis
  vec3 e = (isortho != 0) ? vec3(0.0,0.0,-1.0) : normalize(mid);
  vec3 cr = cross(d, e);
  vec3 pos;
  if (length(cr) < 1e-4){
    // (near) end-on: a screen-facing disc of radius rr covers the projection
    pos = mid + rr * (a_corner.x * vec3(1.0,0.0,0.0) + a_corner.y * vec3(0.0,1.0,0.0));
  } else {
    vec3 side = normalize(cr);
    // span the axis (with an rr overhang for the caps) and the screen-perp
    pos = mid + d * (a_corner.x * (0.5 * len + rr)) + side * (a_corner.y * rr);
  }

  fA = ae;
  fB = be;
  fPosEye = pos;
  fR = rad;
  fColor = a_color;
  fBorder = a_border;
  fBorderColor = a_bordercolor;

  vec4 c = projection * vec4(pos, 1.0);
  gl_Position = (isanchored != 0) ? vec4(anchored_ndc.xy + c.xy, c.z, c.w) : c;
}
