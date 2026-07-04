#version 330 core
layout (location = 0) in vec3 x0;         // anchor (world coords, equilibrium)
layout (location = 1) in vec3 xshift;
layout (location = 2) in vec2 vertex;
layout (location = 3) in vec2 tex;
layout (location = 4) in vec3 x0delta_re; // anchor vibration delta (real part)
layout (location = 5) in vec3 x0delta_im; // anchor vibration delta (imag part)
out vec2 TexCoords;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;
uniform vec3 displ; // complex vibration amplitude (re,im,unused); 0 when not animating
uniform int isanchored; // 1=screen-fixed overlay (window-anchored gizmo)
uniform vec3 anchored_ndc; // target NDC position for the gizmo (xy used)
uniform float anchored_scale; // per-item zoom-compensation factor (gizf)

void main(){
  // animated anchor: vibration displacement applied here so the streamed and
  // cached glyph vertices hold the equilibrium anchor
  vec3 xa = x0 + displ.x * x0delta_re - displ.y * x0delta_im;
  if (isanchored != 0) {
    // screen-fixed overlay: place the anchor xa at anchored_ndc with
    // an orthographic projection scaled by anchored_scale; the
    // per-glyph offset (vertex) is already in NDC. xshift.xy is an
    // eye-space in-plane shift (0 for the gizmo); xshift.z is the
    // eye-space depth offset.
    vec3 e = (mat3(view * world) * xa) * anchored_scale;
    e.xy += xshift.xy;
    vec4 cxy = projection * vec4(e, 0.0);
    e.z += xshift.z;
    vec4 cz = projection * vec4(e, 0.0);
    gl_Position = vec4(anchored_ndc.xy + cxy.xy + vertex, cz.z, 1.0);
  } else {
    vec4 x = view * world * vec4(xa, 1.0);
    x.xy += xshift.xy;
    vec4 c0 = projection * x;
    x.z += xshift.z;
    x = projection * x;
    x /= x.w;
    c0 /= c0.w;
    gl_Position = vec4(c0.xy + vertex, x.z, 1.0);
  }
  TexCoords = tex;
}
