#version 330 core

// Instanced plain mesh: a per-instance model matrix maps the unit mesh onto
// each object; the fragment shader outputs the per-instance flat color.
// Vibration animation (polyhedra triangles) displaces the vertices here:
// the per-instance corner deltas are interpolated with the barycentric
// weights of the reference triangle ((0,0,0),(1,0,0),(0,1,0)), read off
// inPosition.xy. Planes and cones pack zero deltas.

layout (location = 0) in vec3 inPosition; // mesh vertex
layout (location = 1) in vec4 mcol0;      // model matrix columns
layout (location = 2) in vec4 mcol1;
layout (location = 3) in vec4 mcol2;
layout (location = 4) in vec4 mcol3;
layout (location = 5) in vec4 a_color;
layout (location = 6) in vec3 a_d1re;     // vertex 1 vibration delta (real part)
layout (location = 7) in vec3 a_d1im;     // vertex 1 vibration delta (imag part)
layout (location = 8) in vec3 a_d2re;     // vertex 2
layout (location = 9) in vec3 a_d2im;
layout (location = 10) in vec3 a_d3re;    // vertex 3
layout (location = 11) in vec3 a_d3im;

uniform mat4 view;
uniform mat4 world;
uniform mat4 projection;
uniform vec3 displ;           // complex vibration amplitude (re,im,unused); 0 when not animating
uniform int isanchored;       // 1=window-anchored overlay, 0=normal scene
uniform vec3 anchored_ndc;    // target NDC position for the anchor (xy used)
uniform float anchored_scale; // anchor zoom-compensation factor

flat out vec4 fColor;

void main(){
  mat4 model = mat4(mcol0, mcol1, mcol2, mcol3);
  fColor = a_color;

  // vibration displacement: barycentric interpolation of the corner deltas
  // (exact at the corners, linear inside; identical to displacing the corners
  // on the CPU). The weights multiply zeros for planes and cones.
  float w2 = inPosition.x;
  float w3 = inPosition.y;
  float w1 = 1.0 - w2 - w3;
  vec3 dre = w1 * a_d1re + w2 * a_d2re + w3 * a_d3re;
  vec3 dim = w1 * a_d1im + w2 * a_d2im + w3 * a_d3im;
  vec3 local = (model * vec4(inPosition, 1.0)).xyz + displ.x * dre - displ.y * dim;

  if (isanchored != 0){
    // window-anchored: rotation-only of view*world, scaled to a constant
    // on-screen size, placed at anchored_ndc with the (orthographic) projection
    vec4 c = projection * vec4(mat3(view * world) * (local * anchored_scale), 1.0);
    gl_Position = vec4(anchored_ndc.xy + c.xy, c.z, c.w);
  } else {
    gl_Position = projection * view * world * vec4(local, 1.0);
  }
}
