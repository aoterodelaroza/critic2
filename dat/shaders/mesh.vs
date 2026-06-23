#version 330 core

// Instanced plain mesh: a per-instance model matrix maps the unit mesh onto
// each object; the fragment shader outputs the per-instance flat color.

layout (location = 0) in vec3 inPosition; // mesh vertex
layout (location = 1) in vec4 mcol0;      // model matrix columns
layout (location = 2) in vec4 mcol1;
layout (location = 3) in vec4 mcol2;
layout (location = 4) in vec4 mcol3;
layout (location = 5) in vec4 a_color;

uniform mat4 view;
uniform mat4 world;
uniform mat4 projection;
uniform int isanchored;       // 1=window-anchored overlay, 0=normal scene
uniform vec3 anchored_ndc;    // target NDC position for the anchor (xy used)
uniform float anchored_scale; // anchor zoom-compensation factor

flat out vec4 fColor;

void main(){
  mat4 model = mat4(mcol0, mcol1, mcol2, mcol3);
  fColor = a_color;
  if (isanchored != 0){
    // window-anchored: rotation-only of view*world, scaled to a constant
    // on-screen size, placed at anchored_ndc with the (orthographic) projection
    vec3 local = (model * vec4(inPosition, 1.0)).xyz;
    vec4 c = projection * vec4(mat3(view * world) * (local * anchored_scale), 1.0);
    gl_Position = vec4(anchored_ndc.xy + c.xy, c.z, c.w);
  } else {
    gl_Position = projection * view * world * model * vec4(inPosition, 1.0);
  }
}
