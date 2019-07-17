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
  // id of the shader program. glDeleteProgram(id): A value of id = 0 will be silently ignored.
  GLuint id = 0;

  Shader(const char* vpath, const char* fpath){
    loadShader(vpath,fpath,true);
  }

  // Adapted from learnopengl.com by Joey de Vries (https://joeydevries.com).
  void loadShader(const char* vpath, const char* fpath, bool crashiffail){
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

    GLuint vertex, fragment, idnew;
    
    vertex = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex, 1, &vcodec, NULL);
    glCompileShader(vertex);
    if (!checkCompileErrors(vertex, "VERTEX", crashiffail))
      goto fail;

    fragment = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment, 1, &fcodec, NULL);
    glCompileShader(fragment);
    if (!checkCompileErrors(fragment, "FRAGMENT", crashiffail))
      goto faildeletevertex;

    id = glCreateProgram();
    glAttachShader(id, vertex);
    glAttachShader(id, fragment);
    glLinkProgram(id);
    if (!checkCompileErrors(id, "PROGRAM", crashiffail))
      goto faildeletefragment;

    idnew = glCreateProgram();
    glAttachShader(idnew, vertex);
    glAttachShader(idnew, fragment);
    glLinkProgram(idnew);
    if (!checkCompileErrors(idnew, "PROGRAM", crashiffail))
      goto faildeletefragment;
    
    glDeleteProgram(id);
    id = idnew;

    faildeletefragment:
    glDeleteShader(fragment);

    faildeletevertex:
    glDeleteShader(vertex);

    fail:
    return;
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
  bool checkCompileErrors(unsigned int shader, std::string type, bool crashiffail) const{
    int success;
    char infoLog[1024];
    if (type != "PROGRAM"){
      glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
      if (!success){
	glGetShaderInfoLog(shader, 1024, NULL, infoLog);
	std::cout << "Shader compilation error of type: " << type << std::endl << infoLog << std::endl;
        if (crashiffail) 
          exit(EXIT_FAILURE);
        else
          return false;
      }
    } else {
      glGetProgramiv(shader, GL_LINK_STATUS, &success);
      if (!success){
	glGetProgramInfoLog(shader, 1024, NULL, infoLog);
	std::cout << "Shader linking error of type: " << type << std::endl << infoLog << std::endl;
        if (crashiffail) 
          exit(EXIT_FAILURE);
        else
          return false;
      }
    }
    return true;
  }
};
#endif
