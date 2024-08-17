#version 330 core
uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;
uniform float rborder_sph;
uniform float rborder_cyl;

layout (location = 0) in vec3 inPosition;

out vec3 center;
out vec3 vertex;
out vec3 up;
out vec3 side;

void main(){
  mat4 vwm = view * world * model;

  vec4 vertex4 = vwm * vec4(inPosition, 1.0);

  if (rborder_sph > 0 || rborder_cyl > 0){
    vertex = vertex4.xyz / vertex4.w;

    vec4 center4 = vwm * vec4(0.,0.,0.,1.);
    center = center4.xyz / center4.w;

    if (rborder_cyl > 0){
      vec4 up4 = vwm * vec4(0.,0.,0.5,1.);
      up = up4.xyz / up4.w;
      vec4 side4 = vwm * vec4(0.,0.5,0.,1.);
      side = side4.xyz / side4.w;

      // vec3 upn = normalize(up-center);
      // vec3 n = vec3(0.,0.,-1.);
      // vec3 side2 = normalize(cross(n,upn));
      //
      // vertex = vertex + side2 * 0.1;
      // vertex4 = vec4(vertex,1.0);
    } else {
      up = vec3(0.,0.,0.);
    }
  } else {
    center = vec3(0.,0.,0.);
    vertex = vec3(0.,0.,0.);
    up = vec3(0.,0.,0.);
  }
  gl_Position = projection * vertex4;
}
