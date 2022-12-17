#version 330 core

uniform int uselighting;
uniform vec3 vColor;
uniform vec3 lightPos;
uniform vec3 lightColor;
uniform float ambient;
uniform float diffuse;
uniform float specular;
uniform int shininess;

in vec3 Normal;
in vec4 fragPos;

out vec4 outputColor;

void main(){
  if (uselighting != 0){
    vec3 viewdir = normalize(-vec3(fragPos));
    vec3 lightdir = normalize(lightPos - vec3(fragPos));
    vec3 reflectdir = reflect(-lightdir,Normal);
    float diff = diffuse * max(dot(Normal,lightdir),0.f);
    float spec = specular * pow(max(dot(viewdir,reflectdir),0.0f),shininess);
    outputColor = vec4((ambient + diff + spec) * lightColor,1.0f) * vec4(vColor,1.0f);
  } else {
    outputColor = vec4(vColor,1.0f);
  }
}
