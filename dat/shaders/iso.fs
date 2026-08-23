#version 330 core

// Two-sided Blinn-Phong headlight shading for isosurface meshes; the mesh
// color and opacity come from the rgba uniform.

in vec3 vNormal;
out vec4 outColor;

uniform vec4 rgba;

void main(){
  // headlight: light and viewer at the camera; in view space the camera looks
  // down -z, so the lit side has positive normal z. Two-sided: flip the
  // normal towards the camera.
  vec3 n = normalize(vNormal);
  if (n.z < 0.0) n = -n;

  float diff = n.z;               // dot(n, l) with l = (0,0,1)
  float spec = pow(n.z, 32.0);    // halfway vector = l = viewer direction
  vec3 color = rgba.rgb * (0.25 + 0.75 * diff) + 0.15 * spec;
  outColor = vec4(color, rgba.a);
}
