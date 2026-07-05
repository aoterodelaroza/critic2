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
flat in float fOcc;
flat in vec3 fOccEmpty;
flat in vec3 fPieCum;    // cumulative sector boundaries t2,t3,ttot (mixed sites)
flat in vec3 fPieCol2;   // color of pie sector 2
flat in vec3 fPieCol3;   // color of pie sector 3
flat in vec3 fPieCol4;   // color of pie sector 4

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

  // Occupancy pie: up to four species sectors proportional to their occupancy,
  // and the empty-sector color (fOccEmpty) for any residual vacancy. Sector
  // boundaries are the cumulative fractions t1=fOcc, (t2,t3,ttot)=fPieCum,
  // measured from the +y axis so the 0.5-occupancy division is vertical. For a
  // full/single-species atom the extra sectors are zero-width.
  const float TAU = 6.2831853;
  vec4 base = fColor;
  // enter the pie for a partially-occupied atom, or a mixed site with a real
  // second sector (fPieCum.x = t2 > t1 = fOcc) even if the representative is ~full
  if (fOcc < 0.999 || fPieCum.x > fOcc + 1e-4){
    vec2 p = vx.xy;
    float pf = (atan(p.y, p.x) - 1.5707963) / TAU;  // fraction from +y, [0,1)
    if (pf < 0.0) pf += 1.0;

    float t1 = fOcc, t2 = fPieCum.x, t3 = fPieCum.y, ttot = fPieCum.z;
    vec3 col;
    if (pf < t1)        col = fColor.rgb;
    else if (pf < t2)   col = fPieCol2;
    else if (pf < t3)   col = fPieCol3;
    else if (pf < ttot) col = fPieCol4;
    else                col = fOccEmpty;
    base = vec4(col, fColor.a);

    // border along each sector-division ray, same thickness as the rim border
    float halfw = 0.5 * fBorder;
    float bnd[5];
    bnd[0] = 0.0; bnd[1] = t1; bnd[2] = t2; bnd[3] = t3; bnd[4] = ttot;
    for (int k = 0; k < 5; k++){
      vec2 d = vec2(-sin(bnd[k]*TAU), cos(bnd[k]*TAU));
      if (dot(p, d) >= 0.0 && abs(p.x*d.y - p.y*d.x) < halfw)
        base = vec4(fBorderColor, fColor.a);
    }
  }

  if (fRadius - rproj < fBorder)
    outColor = vec4(fBorderColor, fColor.a);
  else
    outColor = base;
}
