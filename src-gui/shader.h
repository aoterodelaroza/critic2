// -*-c++-*-

// based on the shader class from learnopengl.com

#ifndef SHADER_H
#define SHADER_H

#include "imgui/gl3w.h"
#include <iostream>

class Shader{
public:
  GLuint id;

  Shader(){
    id = glCreateProgram();
    if (!id)
      exit(EXIT_FAILURE);

    const char *vs = 
      "#version 330 core\n"
      "uniform mat4 model; \n"
      "uniform mat4 view; \n"
      "uniform mat4 projection; \n"
      "layout (location = 0) in vec3 inPosition; \n"
      "void main() { \n"
      "  gl_Position = vec4(inPosition, 1.0) * model * view * projection; \n"
      "}";    

//      "gl_Position = gWVP * vec4(inPosition, 1.0); \n"
	
    const char *fs = 
      "#version 330 core\n"
      "uniform vec4 vColor; \n"
      "out vec4 outputColor; \n"
      "void main() { \n"
      "outputColor = vColor; \n"
      "}";

    // vertex shader
    GLuint vertex = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex, 1, &vs, NULL);
    glCompileShader(vertex);
    checkCompileErrors(vertex, "VERTEX");
    // fragment Shader
    GLuint fragment = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment, 1, &fs, NULL);
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

  void setBool(const std::string &name, bool value) const{         
    glUniform1i(glGetUniformLocation(id, name.c_str()), (int)value); 
  }

  void setInt(const std::string &name, int value) const{ 
    glUniform1i(glGetUniformLocation(id, name.c_str()), value); 
  }

  void setFloat(const std::string &name, float value) const{ 
    glUniform1f(glGetUniformLocation(id, name.c_str()), value); 
  }

private:
  // utility function for checking shader compilation/linking errors.
  void checkCompileErrors(unsigned int shader, std::string type){
    int success;
    char infoLog[1024];
    if (type != "PROGRAM"){
      glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
      if (!success){
	glGetShaderInfoLog(shader, 1024, NULL, infoLog);
	std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << std::endl << infoLog << std::endl;
      }
    } else {
      glGetProgramiv(shader, GL_LINK_STATUS, &success);
      if (!success){
	glGetProgramInfoLog(shader, 1024, NULL, infoLog);
	std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << std::endl << infoLog << std::endl;
      }
    }
  }
};
#endif
