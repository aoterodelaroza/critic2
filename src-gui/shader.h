// -*-c++-*-

// based on the shader class from learnopengl.com

#ifndef SHADER_H
#define SHADER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "imgui/gl3w.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Shader{
public:
  GLuint id;

  // Constructor adapted from learnopengl.com by Joey de Vries (https://joeydevries.com).
  Shader(const char* vpath, const char* fpath){
    std::string vcode, fcode;
    std::ifstream vfile, ffile;

    // ensure ifstream objects can throw exceptions:
    vfile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
    ffile.exceptions (std::ifstream::failbit | std::ifstream::badbit);

    try{
      vfile.open(vpath);
      ffile.open(fpath);

      std::stringstream vstream, fstream;
      vstream << vfile.rdbuf();
      fstream << ffile.rdbuf();

      vfile.close();
      ffile.close();

      vcode = vstream.str();
      fcode = fstream.str();
    } catch (std::ifstream::failure e) {
      printf("Could not find shader file: %s or %s\n",vpath,fpath);
      exit(1);
    }
    const char* vcodec = vcode.c_str();
    const char* fcodec = fcode.c_str();

    unsigned int vertex = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex, 1, &vcodec, NULL);
    glCompileShader(vertex);
    checkCompileErrors(vertex, "VERTEX");

    unsigned int fragment = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment, 1, &fcodec, NULL);
    glCompileShader(fragment);
    checkCompileErrors(fragment, "FRAGMENT");

    id = glCreateProgram();
    glAttachShader(id, vertex);
    glAttachShader(id, fragment);
    glLinkProgram(id);

    glDeleteShader(vertex);
    glDeleteShader(fragment);
  }

  void use(){ 
    glUseProgram(id); 
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

private:
  void checkCompileErrors(unsigned int shader, std::string type) const{
    int success;
    char infoLog[1024];
    if (type != "PROGRAM"){
      glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
      if (!success){
	glGetShaderInfoLog(shader, 1024, NULL, infoLog);
	std::cout << "Shader compilation error of type: " << type << std::endl << infoLog << std::endl;
	exit(EXIT_FAILURE);
      }
    } else {
      glGetProgramiv(shader, GL_LINK_STATUS, &success);
      if (!success){
	glGetProgramInfoLog(shader, 1024, NULL, infoLog);
	std::cout << "Shader linking error of type: " << type << std::endl << infoLog << std::endl;
	exit(EXIT_FAILURE);
      }
    }
  }
};
#endif
