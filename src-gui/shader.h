// -*-c++-*-

// based on the shader class from learnopengl.com

#ifndef SHADER_H
#define SHADER_H

#include <iostream>
#include <string>

#include "imgui/gl3w.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Shader{
public:
  GLuint id;
  GLuint idtext;

  Shader(){
    id = glCreateProgram();
    if (!id)
      exit(EXIT_FAILURE);
    idtext = glCreateProgram();
    if (!idtext)
      exit(EXIT_FAILURE);

    const char *vs[] = {
      "#version 330 core                                                          \n"
      "uniform mat4 model;                                                        \n"
      "uniform mat4 projection;                                                   \n"
      "uniform mat4 view;                                                         \n"
      "uniform mat4 world;                                                        \n"
      "uniform mat3 normrot;                                                      \n"
      "                                                                           \n"
      "layout (location = 0) in vec3 inPosition;                                  \n"
      "layout (location = 1) in vec3 inNormal;                                    \n"
      "                                                                           \n"
      "out vec3 Normal;                                                           \n"
      "out vec4 fragPos;                                                          \n"
      "                                                                           \n"
      "void main() {                                                              \n"
      "  fragPos = view * world * model * vec4(inPosition, 1.0);                  \n"
      "  gl_Position = projection * fragPos;                                      \n"
      "  Normal = normalize(normrot * inNormal);                                  \n"
      "}\n"
    };

    const char *fs[] =  {
      "#version 330 core                                                             \n"
      "                                                                              \n"
      "uniform int uselighting;                                                      \n"
      "uniform vec4 vColor;                                                          \n"
      "uniform vec3 lightPos;                                                        \n"
      "uniform vec3 lightColor;                                                      \n"
      "uniform float ambient;                                                        \n"
      "uniform float diffuse;                                                        \n"
      "uniform float specular;                                                       \n"
      "uniform int shininess;                                                        \n"
      "                                                                              \n"
      "in vec3 Normal;                                                               \n"
      "in vec4 fragPos;                                                              \n"
      "                                                                              \n"
      "out vec4 outputColor;                                                         \n"
      "                                                                              \n"
      "void main() {                                                                 \n"
      "  if (uselighting != 0){                                                      \n"
      "    vec3 viewdir = normalize(-vec3(fragPos));                                 \n"
      "    vec3 lightdir = normalize(lightPos - vec3(fragPos));                      \n"
      "    vec3 reflectdir = reflect(-lightdir,Normal);                              \n"
      "    float diff = diffuse * max(dot(Normal,lightdir),0.f);                     \n"
      "    float spec = specular * pow(max(dot(viewdir,reflectdir),0.0f),shininess); \n"
      "    outputColor = vec4((ambient + diff + spec) * lightColor,1.0f) * vColor;   \n"
      "  } else {                                                                    \n"
      "    outputColor = vColor;                                                     \n"
      "  }                                                                           \n"
      "}\n"
    };

    const char *vstext[] = {
      "#version 330 core\n"
      "layout (location = 0) in vec4 vertex;\n"
      "out vec2 TexCoords;\n"
      "\n"
      "uniform mat4 projection;\n"
      "\n"
      "void main()\n"
      "{\n"
      "    gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);\n"
      "    TexCoords = vertex.zw;\n"
      "}  \n"
    };

    const char *fstext[] = {
      "#version 330 core\n"
      "in vec2 TexCoords;\n"
      "out vec4 color;\n"
      "\n"
      "uniform sampler2D text;\n"
      "uniform vec3 textColor;\n"
      "\n"
      "void main()\n"
      "{    \n"
      "    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);\n"
      "    color = vec4(textColor, 1.0) * sampled;\n"
      "}  \n"
    };

    // vertex shader
    GLuint vertex = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex, 1, vs, NULL);
    glCompileShader(vertex);
    checkCompileErrors(vertex, "VERTEX");
    // fragment shader
    GLuint fragment = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment, 1, fs, NULL);
    glCompileShader(fragment);
    checkCompileErrors(fragment, "FRAGMENT");
    // vertex shader (text)
    GLuint vertextext = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertextext, 1, vstext, NULL);
    glCompileShader(vertextext);
    checkCompileErrors(vertextext, "VERTEXTEXT");
    // fragment shader (text)
    GLuint fragmenttext = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmenttext, 1, fstext, NULL);
    glCompileShader(fragmenttext);
    checkCompileErrors(fragmenttext, "FRAGMENTTEXT");

    // shader program
    glAttachShader(id, vertex);
    glAttachShader(id, fragment);
    glLinkProgram(id);
    checkCompileErrors(id, "PROGRAM");
    glDeleteShader(vertex);
    glDeleteShader(fragment);

    // shader program (text)
    glAttachShader(idtext, vertextext);
    glAttachShader(idtext, fragmenttext);
    glLinkProgram(idtext);
    checkCompileErrors(idtext, "PROGRAM");
    glDeleteShader(vertextext);
    glDeleteShader(fragmenttext);

    // set the global variables
    glUseProgram(id); 
  }

  void use(){ 
    glUseProgram(id); 
  }
  void usetext(){ 
    glUseProgram(idtext); 
  }

  void setMat3(const char *name, const GLfloat * value) const{
    glUniformMatrix3fv(glGetUniformLocation(id,name), 1, GL_FALSE, value);
  }
  void setMat4(const char *name, const GLfloat * value) const{
    glUniformMatrix4fv(glGetUniformLocation(id,name), 1, GL_FALSE, value);
  }
  void setVec3(const char *name, const GLfloat * value) const{
    glUniform3fv(glGetUniformLocation(id,name), 1, value);
  }
  void setVec4(const char *name, const GLfloat * value) const{
    glUniform4fv(glGetUniformLocation(id,name), 1, value);
  }
  void setInt(const char *name, int value) const{ 
    glUniform1i(glGetUniformLocation(id,name), value); 
  }
  void setFloat(const char *name, float value) const{ 
    glUniform1f(glGetUniformLocation(id,name), value); 
  }
  void setTextProjection(float atex) const{ 
    glm::mat4 projection = glm::ortho(0.0f, atex, 0.0f, atex);
    glUniformMatrix4fv(glGetUniformLocation(idtext,"projection"), 1, GL_FALSE, glm::value_ptr(projection));
  }
  void setTextColor(const GLfloat *value) const{ 
    glUniform3fv(glGetUniformLocation(idtext,"textColor"), 1, value);
  }

private:
  // utility function for checking shader compilation/linking errors.
  void checkCompileErrors(unsigned int shader, std::string type) const{
    int success;
    char infoLog[1024];
    if (type != "PROGRAM"){
      glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
      if (!success){
	glGetShaderInfoLog(shader, 1024, NULL, infoLog);
	std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << std::endl << infoLog << std::endl;
	exit(EXIT_FAILURE);
      }
    } else {
      glGetProgramiv(shader, GL_LINK_STATUS, &success);
      if (!success){
	glGetProgramInfoLog(shader, 1024, NULL, infoLog);
	std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << std::endl << infoLog << std::endl;
	exit(EXIT_FAILURE);
      }
    }
  }
};
#endif
