#version 330 core

// Ray-cast sphere fragment shader. Solves the ray-sphere intersection
// in eye space, writes the true depth, and applies flat color + a
// screen-space silhouette border. In pick mode it outputs the
// per-instance index instead of a color.

in vec3 fCenterEye;
in float fRadius;
in vec3 fPosEye;
flat in vec4 fColor;
flat in float fBorder;
flat in vec3 fBorderColor;
flat in vec4 fIdx;

uniform mat4 projection;
uniform int isortho;   // 1=orthographic, 0=perspective
uniform int uPick;     // 1=output index, 0=output color

out vec4 outColor;

void main(){
  // view ray in eye space. For orthographic rays start the origin in front of
  // the sphere (window-anchored geometry can straddle z=0, so a z=0 origin would
  // put a sphere facing the camera behind the ray).
  vec3 ro, rd;
  if (isortho != 0){
    ro = vec3(fPosEye.xy, fCenterEye.z + fRadius + 1.0);
    rd = vec3(0.0, 0.0, -1.0);
  } else {
    ro = vec3(0.0);
    rd = normalize(fPosEye);
  }

  // ray-sphere intersection (center fCenterEye, radius fRadius)
  vec3 oc = ro - fCenterEye;
  float b = dot(oc, rd);
  float c = dot(oc, oc) - fRadius * fRadius;
  float disc = b * b - c;
  if (disc < 0.0) discard;
  float t = -b - sqrt(disc);
  vec3 hit = ro + t * rd;

  // write the true fragment depth
  vec4 clip = projection * vec4(hit, 1.0);
  gl_FragDepth = 0.5 * (clip.z / clip.w) + 0.5;

  if (uPick != 0){
    outColor = fIdx;
    return;
  }

  // flat color + silhouette border: border when the surface point is close to
  // the rim (its in-screen distance from the center approaches the radius)
  vec3 vx = hit - fCenterEye;                 // |vx| == fRadius
  float rproj = length(vx - dot(vx, rd) * rd);
  if (fRadius - rproj < fBorder)
    outColor = vec4(fBorderColor, fColor.a);
  else
    outColor = fColor;
}
