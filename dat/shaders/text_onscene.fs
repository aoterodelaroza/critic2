#version 330 core
in vec2 TexCoords;
out vec4 color;

uniform sampler2D text;
uniform vec3 textColor;

void main(){
    color = vec4(textColor, 1.0) * texture(text, TexCoords);
    if (color.a > 0.25)
       color.a = 1.0;
    else
       color.a = 0.0;
}
