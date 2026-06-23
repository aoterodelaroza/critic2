#version 330 core

// Instanced sphere shader. Each instance is a camera-facing quad
// sized to cover the sphere silhouette; the fragment shader ray-casts
// the actual sphere.

// per-vertex: quad corner in [-1,1]^2
layout (location = 0) in vec2 a_corner;

// per-instance attributes
layout (location = 1) in vec3 a_center;      // sphere center (pre-world coords)
layout (location = 2) in float a_radius;     // sphere radius
layout (location = 3) in vec4 a_color;       // rgba color
layout (location = 4) in float a_border;     // silhouette border width
layout (location = 5) in vec3 a_bordercolor; // border color
layout (location = 6) in vec3 a_xdelta_re;   // vibration delta (real part)
layout (location = 7) in vec3 a_xdelta_im;   // vibration delta (imag part)
layout (location = 8) in vec4 a_idx;         // picking index (bit-packed floats)

uniform mat4 view;
uniform mat4 world;
uniform mat4 projection;
uniform vec3 displ;       // complex vibration amplitude (re,im,unused); 0 when not animating
uniform int isortho;      // 1=orthographic, 0=perspective
uniform int isanchored;      // 1=window-anchored overlay, 0=normal scene
uniform vec3 anchored_ndc;   // target NDC position for the anchor (xy used)
uniform float anchored_scale;// anchor zoom-compensation factor

out vec3 fCenterEye;      // sphere center in eye space
out float fRadius;
out vec3 fPosEye;         // billboard point in eye space
flat out vec4 fColor;
flat out float fBorder;
flat out vec3 fBorderColor;
flat out vec4 fIdx;

void main(){
  // animated center (pre-world coords), then to eye space. For a window-anchored
  // object only the rotation of view*world is applied (no eye translation) and
  // the geometry is scaled to a constant on-screen size.
  mat4 vm = view * world;
  vec3 cpre = a_center + displ.x * a_xdelta_re - displ.y * a_xdelta_im;
  vec3 ec;
  float rad = a_radius;
  if (isanchored != 0){
    ec = mat3(vm) * (cpre * anchored_scale);
    rad = a_radius * anchored_scale;
  } else {
    ec = (vm * vec4(cpre, 1.0)).xyz;
  }

  // camera-facing quad in eye space (camera looks down -z). For perspective,
  // enlarge so the quad always covers the projected silhouette of the sphere.
  float d = max(length(ec), rad + 1e-4);
  float fac = (isortho != 0) ? 1.0 : d / sqrt(max(d*d - rad*rad, 1e-6));
  fac = min(fac, 4.0);
  vec3 pos = ec + rad * fac * vec3(a_corner, 0.0);

  fCenterEye = ec;
  fRadius = rad;
  fPosEye = pos;
  fColor = a_color;
  fBorder = a_border;
  fBorderColor = a_bordercolor;
  fIdx = a_idx;

  vec4 c = projection * vec4(pos, 1.0);
  gl_Position = (isanchored != 0) ? vec4(anchored_ndc.xy + c.xy, c.z, c.w) : c;
}
