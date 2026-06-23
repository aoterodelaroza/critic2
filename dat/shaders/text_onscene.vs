#version 330 core
layout (location = 0) in vec3 x0;
layout (location = 1) in vec3 xshift;
layout (location = 2) in vec2 vertex;
layout (location = 3) in vec2 tex;
out vec2 TexCoords;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;
uniform float depth;
uniform int isgizmo; // 1=screen-fixed overlay (window-anchored gizmo)
uniform vec3 gizmo_ndc; // target NDC position for the gizmo (xy used)
uniform float gizmo_scale; // per-item zoom-compensation factor (gizf)

void main(){
  if (isgizmo != 0) {
    // screen-fixed overlay: place the (gizmo-local) anchor x0 at gizmo_ndc with
    // an orthographic projection scaled by gizmo_scale; the per-glyph offset
    // (vertex) is already in NDC. xshift.xy is an eye-space in-plane shift (0
    // for the gizmo); xshift.z is the eye-space depth offset.
    vec3 e = (mat3(view * world) * x0) * gizmo_scale;
    e.xy += xshift.xy;
    vec4 cxy = projection * vec4(e, 0.0);
    e.z += xshift.z;
    vec4 cz = projection * vec4(e, 0.0);
    gl_Position = vec4(gizmo_ndc.xy + cxy.xy + vertex, cz.z, 1.0);
  } else {
    vec4 x = view * world * vec4(x0, 1.0);
    x.xy += xshift.xy;
    vec4 x0 = projection * x;
    x.z += xshift.z;
    x = projection * x;
    x /= x.w;
    x0 /= x0.w;
    gl_Position = vec4(x0.xy + vertex, x.z, 1.0);
  }
  TexCoords = tex;
}
