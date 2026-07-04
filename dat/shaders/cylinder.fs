#version 330 core

// Ray-cast capped-cylinder impostor. Intersects the view ray with a finite
// cylinder (side + two flat caps) in eye space, writes the true depth, and
// applies flat color + a screen-space silhouette border.

in vec3 fA;
in vec3 fB;
in vec3 fPosEye;
in float fR;
flat in vec4 fColor;
flat in float fBorder;
flat in vec3 fBorderColor;

uniform mat4 projection;
uniform int isortho;

out vec4 outColor;

void main(){
  // view ray in eye space. For orthographic rays start the origin in front of
  // the whole cylinder (window-anchored geometry can straddle z=0, so a z=0
  // origin would put the half facing the camera behind the ray).
  vec3 ro, rd;
  if (isortho != 0){
    ro = vec3(fPosEye.xy, max(fA.z, fB.z) + fR + 1.0);
    rd = vec3(0.0, 0.0, -1.0);
  } else {
    ro = vec3(0.0);
    rd = normalize(fPosEye);
  }

  vec3 pa = fA;
  vec3 axis = fB - fA;
  float len = max(length(axis), 1e-7);
  vec3 d = axis / len;

  // ray vs infinite cylinder around the axis (components perpendicular to d)
  vec3 dp = ro - pa;
  float rdd = dot(rd, d);
  float dpd = dot(dp, d);
  vec3 rdp = rd - rdd * d;
  vec3 dpp = dp - dpd * d;
  float A = dot(rdp, rdp);
  float B = 2.0 * dot(dpp, rdp);
  float C = dot(dpp, dpp) - fR * fR;

  float tbest = 1e30;
  float hbest = 0.0;
  bool hit = false;

  // side of the cylinder
  if (A > 1e-12){
    float disc = B * B - 4.0 * A * C;
    if (disc >= 0.0){
      float t = (-B - sqrt(disc)) / (2.0 * A);
      float h = dpd + t * rdd;
      if (t > 0.0 && h >= 0.0 && h <= len){
        tbest = t; hbest = h; hit = true;
      }
    }
  }

  // flat caps at h=0 and h=len
  if (abs(rdd) > 1e-12){
    float t0 = (0.0 - dpd) / rdd;
    if (t0 > 0.0 && t0 < tbest){
      vec3 v = (ro + t0 * rd) - pa;
      vec3 vr = v - dot(v, d) * d;
      if (dot(vr, vr) <= fR * fR){ tbest = t0; hbest = 0.0; hit = true; }
    }
    float t1 = (len - dpd) / rdd;
    if (t1 > 0.0 && t1 < tbest){
      vec3 v = (ro + t1 * rd) - (pa + len * d);
      vec3 vr = v - dot(v, d) * d;
      if (dot(vr, vr) <= fR * fR){ tbest = t1; hbest = len; hit = true; }
    }
  }

  if (!hit) discard;

  vec3 hitp = ro + tbest * rd;

  // true fragment depth
  vec4 clip = projection * vec4(hitp, 1.0);
  gl_FragDepth = 0.5 * (clip.z / clip.w) + 0.5;

  // flat color + silhouette border: the border shows where the surface radial
  // direction lines up with the screen-perpendicular of the axis
  vec3 cv = cross(rd, d);
  vec3 vperp = (length(cv) > 1e-6) ? cv / length(cv) : vec3(0.0);
  vec3 rad = (hitp - pa) - dot(hitp - pa, d) * d;
  float rl = length(rad);
  float proj = (rl > 1e-6) ? abs(dot(rad / rl, vperp)) : 0.0;
  if (fR - fR * proj < fBorder)
    outColor = vec4(fBorderColor, fColor.a);
  else
    outColor = fColor;
}
