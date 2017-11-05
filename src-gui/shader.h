// -*-c++-*-

// based on the shader class from learnopengl.com

#ifndef SHADER_H
#define SHADER_H

#include <glm/glm.hpp>
#include "imgui/gl3w.h"
#include <iostream>

using namespace std;
using namespace glm;

class Shader{
public:
  GLuint id;

  Shader(){
    id = glCreateProgram();
    if (!id)
      exit(EXIT_FAILURE);

    const char *vs[] =
      {
        "#version 330 core                                                          \n"
        "uniform mat4 model;                                                        \n"
        "uniform mat4 projection;                                                   \n"
        "uniform mat4 view;                                                         \n"
        "uniform mat4 world;                                                        \n"
        "layout (location = 0) in vec3 inPosition;                                  \n"
        "void main() {                                                              \n"
        "  gl_Position = projection * view * world * model * vec4(inPosition, 1.0); \n"
        "}\n"
      };

    const char *fs[] =  
      {
        "#version 330 core     \n"
        "uniform vec4 vColor;  \n"
        "out vec4 outputColor; \n"
        "void main() {         \n"
        "outputColor = vColor; \n"
        "}\n"
      };

    // vertex shader
    GLuint vertex = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex, 1, vs, NULL);
    glCompileShader(vertex);
    checkCompileErrors(vertex, "VERTEX");
    // fragment Shader
    GLuint fragment = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment, 1, fs, NULL);
    glCompileShader(fragment);
    checkCompileErrors(fragment, "FRAGMENT");
    // shader Program
    glAttachShader(id, vertex);
    glAttachShader(id, fragment);
    glLinkProgram(id);
    checkCompileErrors(id, "PROGRAM");
    glDeleteShader(vertex);
    glDeleteShader(fragment);
  }

  void use(){ 
    glUseProgram(id); 
  }

  void setMat4(const char *name, const GLfloat * value) const{
    glUniformMatrix4fv(glGetUniformLocation(id,name), 1, GL_FALSE, value);
  }
  void setVec4(const char *name, const GLfloat * value) const{
    glUniform4fv(glGetUniformLocation(id,name), 1, value);
  }

private:
  // utility function for checking shader compilation/linking errors.
  void checkCompileErrors(unsigned int shader, string type){
    int success;
    char infoLog[1024];
    if (type != "PROGRAM"){
      glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
      if (!success){
	glGetShaderInfoLog(shader, 1024, NULL, infoLog);
	cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << endl << infoLog << endl;
      }
    } else {
      glGetProgramiv(shader, GL_LINK_STATUS, &success);
      if (!success){
	glGetProgramInfoLog(shader, 1024, NULL, infoLog);
	cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << endl << infoLog << endl;
      }
    }
  }
};
#endif
