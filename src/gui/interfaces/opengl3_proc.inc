! -*-F90-*-

!xx! OpenGL3 procedures

! -- Type equivalence chart (see glcorearb.h) --
! GLbitfield = GLenum = GLuint = GLint = GLsizei = c_int
! GLhalf = GLushort = GLshort = c_short
! GLclampf = GLfloat = c_float
! GLubyte = GLbyte = GLchar = GLboolean = c_char
! GLdouble = GLclampd = c_double
! GLsizeiptr = ptrdiff_t
! GLvoid = void

! void glAttachShader(GLuint program, GLuint shader);
subroutine glAttachShader(prog,shader) bind(c,name="glAttachShader")
  import c_int
  integer(c_int), value :: prog
  integer(c_int), value :: shader
end subroutine glAttachShader
! void glBindBuffer(GLenum target, GLuint buffer);
subroutine glBindBuffer(target_,buffer) bind(c,name="glBindBuffer")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: buffer
end subroutine glBindBuffer
! void glBindFramebuffer(GLenum target, GLuint framebuffer);
subroutine glBindFramebuffer(target_,framebuffer) bind(c,name="glBindFramebuffer")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: framebuffer
end subroutine glBindFramebuffer
! void glBindRenderbuffer(GLenum target, GLuint renderbuffer);
subroutine glBindRenderbuffer(target_,renderbuffer) bind(c,name="glBindRenderbuffer")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: renderbuffer
end subroutine glBindRenderbuffer
! void glBindTexture(GLenum target, GLuint texture);
subroutine glBindTexture(target_,texture) bind(c,name="glBindTexture")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: texture
end subroutine glBindTexture
! void glBufferData(GLenum target, GLsizeiptr size, const GLvoid *data, GLenum usage);
subroutine glBufferData(target_,size,data,usage) bind(c,name="glBufferData")
  import c_int, c_int64_t, c_ptr, c_ptrdiff_t
  integer(c_int), value :: target_
  integer(c_ptrdiff_t), value :: size
  type(c_ptr), intent(in), value :: data
  integer(c_int), value :: usage
end subroutine glBufferData
! void glBindVertexArray(GLuint array);
subroutine glBindVertexArray(array) bind(c,name="glBindVertexArray")
  import c_int
  integer(c_int), value :: array
end subroutine glBindVertexArray
! GLenum glCheckFramebufferStatus(GLenum target);
function glCheckFramebufferStatus(target_) bind(c,name="glCheckFramebufferStatus")
  import c_int
  integer(c_int), value :: target_
  integer(c_int) :: glCheckFramebufferStatus
end function glCheckFramebufferStatus
! void glClear (GLbitfield mask);
subroutine glClear (mask) bind(c,name="glClear")
  import c_int
  integer(c_int), value :: mask
end subroutine glClear
! void glClearColor(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha)
subroutine glClearColor(red, green, blue, alpha) bind(c,name="glClearColor")
  import c_float
  real(c_float), value :: red, green, blue, alpha
end subroutine glClearColor
! void glCompileShader(GLuint shader)
subroutine glCompileShader(shader) bind(c,name="glCompileShader")
  import c_int
  integer(c_int), value :: shader
end subroutine glCompileShader
! GLuint glCreateProgram(void);
function glCreateProgram() bind(c,name="glCreateProgram")
  import c_int
  integer(c_int) :: glCreateProgram
end function glCreateProgram
! GLuint glCreateShader(GLenum shaderType);
function glCreateShader(shaderType) bind(c,name="glCreateShader")
  import c_int
  integer(c_int), value :: shaderType
  integer(c_int) :: glCreateShader
end function glCreateShader
! void glDeleteBuffers(GLsizei n, const GLuint *buffers)
subroutine glDeleteBuffers(n,buffers) bind(c,name="glDeleteBuffers")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: buffers
end subroutine glDeleteBuffers
! void glDeleteFramebuffers(GLsizei n, GLuint *framebuffers);
subroutine glDeleteFramebuffers(n,framebuffers) bind(c,name="glDeleteFramebuffers")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: framebuffers
end subroutine glDeleteFramebuffers
! void glDeleteProgram(GLuint program);
subroutine glDeleteProgram(prog) bind(c,name="glDeleteProgram")
  import c_int
  integer(c_int), value :: prog
end subroutine glDeleteProgram
! void glDeleteRenderbuffers(GLsizei n, GLuint *renderbuffers);
subroutine glDeleteRenderbuffers(n,renderbuffers) bind(c,name="glDeleteRenderbuffers")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: renderbuffers
end subroutine glDeleteRenderbuffers
! void glDeleteShader(GLuint program);
subroutine glDeleteShader(shader) bind(c,name="glDeleteShader")
  import c_int
  integer(c_int), value :: shader
end subroutine glDeleteShader
! void glDeleteTextures(GLsizei n, const GLuint *textures);
subroutine glDeleteTextures(n,textures) bind(c,name="glDeleteTextures")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: textures
end subroutine glDeleteTextures
! void glDeleteVertexArrays(GLsizei n, const GLuint *arrays)
subroutine glDeleteVertexArrays(n,arrays) bind(c,name="glDeleteVertexArrays")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: arrays
end subroutine glDeleteVertexArrays
! void glDetachShader(GLuint program, GLuint shader)
subroutine glDetachShader(prog,shader) bind(c,name="glDetachShader")
  import c_int
  integer(c_int), value :: prog
  integer(c_int), value :: shader
end subroutine glDetachShader
! void glDrawArrays(GLenum mode, GLint first, GLsizei count);
subroutine glDrawArrays(mode,first,count) bind(c,name="glDrawArrays")
  import c_int
  integer(c_int), value :: mode
  integer(c_int), value :: first
  integer(c_int), value :: count
end subroutine glDrawArrays
! void glEnable(GLenum cap)
subroutine glEnable(cap) bind(c,name="glEnable")
  import c_int
  integer(c_int), value :: cap
end subroutine glEnable
! void glEnableVertexAttribArray(GLuint index);
subroutine glEnableVertexAttribArray(index) bind(c,name="glEnableVertexAttribArray")
  import c_int
  integer(c_int), value :: index
end subroutine glEnableVertexAttribArray
! void glFramebufferRenderbuffer(GLenum target, GLenum attachment, GLenum renderbuffertarget,
!                                GLuint renderbuffer);
subroutine glFramebufferRenderbuffer(target_,attachment,renderbuffertarget,renderbuffer) bind(c,name="glFramebufferRenderbuffer")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: attachment
  integer(c_int), value :: renderbuffertarget
  integer(c_int), value :: renderbuffer
end subroutine glFramebufferRenderbuffer
! void glGenBuffers(GLsizei n, GLuint* buffers);
subroutine glGenBuffers(n,buffers) bind(c,name="glGenBuffers")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), value :: buffers
end subroutine glGenBuffers
! void glGenFramebuffers(GLsizei n,GLuint *ids);
subroutine glGenFramebuffers(n,ids) bind(c,name="glGenFramebuffers")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), value :: ids
end subroutine glGenFramebuffers
! GLenum glGetError( void);
function glGetError() bind(c,name="glGetError")
  import c_int
  integer(c_int) :: glGetError
end function glGetError
! void glGetProgramInfoLog(GLuint program, GLsizei maxLength, GLsizei *length, GLchar *infoLog);
subroutine glGetProgramInfoLog(prog,maxLength,length,infoLog) bind(c,name="glGetProgramInfoLog")
  import c_int, c_ptr
  integer(c_int), value :: prog
  integer(c_int), value :: maxLength
  integer(c_int) :: length
  type(c_ptr), value :: infoLog
end subroutine glGetProgramInfoLog
! void glGetProgramiv(GLuint program, GLenum pname, GLint *params);
subroutine glGetProgramiv(prog,pname,params) bind(c,name="glGetProgramiv")
  import c_int
  integer(c_int), value :: prog
  integer(c_int), value :: pname
  integer(c_int) :: params
end subroutine glGetProgramiv
! void glGetShaderInfoLog(GLuint shader, GLsizei maxLength, GLsizei *length, GLchar *infoLog);
subroutine glGetShaderInfoLog(shader,maxLength,length,infoLog) bind(c,name="glGetShaderInfoLog")
  import c_int, c_ptr
  integer(c_int), value :: shader
  integer(c_int), value :: maxLength
  integer(c_int) :: length
  type(c_ptr), value :: infoLog
end subroutine glGetShaderInfoLog
! void glGetShaderiv(GLuint shader, GLenum pname, GLint *params);
subroutine glGetShaderiv(shader,pname,params) bind(c,name="glGetShaderiv")
  import c_int
  integer(c_int), value :: shader
  integer(c_int), value :: pname
  integer(c_int) :: params
end subroutine glGetShaderiv
! GLint glGetUniformLocation(GLuint program, const GLchar *name);
function glGetUniformLocation(prog_,name) bind(c,name="glGetUniformLocation")
  import c_int, c_ptr
  integer(c_int), value :: prog_
  type(c_ptr), intent(in), value :: name
  integer(c_int) :: glGetUniformLocation
end function glGetUniformLocation
! void glFramebufferTexture2D(GLenum target, GLenum attachment, GLenum textarget, GLuint texture,
!      GLint level);
subroutine glFramebufferTexture2D(target_,attachment,textarget,texture,level) bind(c,name="glFramebufferTexture2D")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: attachment
  integer(c_int), value :: textarget
  integer(c_int), value :: texture
  integer(c_int), value :: level
end subroutine glFramebufferTexture2D
! void glGenRenderbuffers(GLsizei n,GLuint *renderbuffers)
subroutine glGenRenderbuffers(n,renderbuffers) bind(c,name="glGenRenderbuffers")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), value :: renderbuffers
end subroutine glGenRenderbuffers
! void glGenTextures(GLsizei n,GLuint *textures);
subroutine glGenTextures(n,textures) bind(c,name="glGenTextures")
  import c_ptr, c_int
  integer(c_int), value :: n
  type(c_ptr), value :: textures
end subroutine glGenTextures
! void glGenVertexArrays(GLsizei n, GLuint *arrays);
subroutine glGenVertexArrays(n,arrays) bind(c,name="glGenVertexArrays")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), value :: arrays
end subroutine glGenVertexArrays
! void glLinkProgram(GLuint program)
subroutine glLinkProgram(prog) bind(c,name="glLinkProgram")
  import c_int
  integer(c_int), value :: prog
end subroutine glLinkProgram
! void glRenderbufferStorage(GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
subroutine glRenderbufferStorage(target_,internalformat,width,height) bind(c,name="glRenderbufferStorage")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: internalformat
  integer(c_int), value :: width
  integer(c_int), value :: height
end subroutine glRenderbufferStorage
! void glShaderSource(GLuint shader, GLsizei count, const GLchar **string, const GLint *length);
subroutine glShaderSource(shader,count,string,length) bind(c,name="glShaderSource")
  import c_int, c_ptr
  integer(c_int), value :: shader
  integer(c_int), value :: count
  type(c_ptr), intent(in) :: string
  type(c_ptr), intent(in), value :: length
end subroutine glShaderSource
! void glTexImage2D(GLenum target, GLint level, GLint internalFormat, GLsizei width, GLsizei height,
!                   GLint border, GLenum format, GLenum type, const GLvoid * data);
subroutine glTexImage2D(target_,level,internalFormat,width,height,border,format_,type,data) bind(c,name="glTexImage2D")
  import c_int, c_ptr
  integer(c_int), value :: target_
  integer(c_int), value :: level
  integer(c_int), value :: internalFormat
  integer(c_int), value :: width
  integer(c_int), value :: height
  integer(c_int), value :: border
  integer(c_int), value :: format_
  integer(c_int), value :: type
  type(c_ptr), intent(in), value :: data
end subroutine glTexImage2D
! void glTexParameterf(GLenum target_, GLenum pname, GLfloat param);
subroutine glTexParameterf(target_,pname,param) bind(c,name="glTexParameterf")
  import c_int, c_float
  integer(c_int), value :: target_
  integer(c_int), value :: pname
  real(c_float), value :: param
end subroutine glTexParameterf
! void glTexParameteri(GLenum target_, GLenum pname, GLint param);
subroutine glTexParameteri(target_,pname,param) bind(c,name="glTexParameteri")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: pname
  integer(c_int), value :: param
end subroutine glTexParameteri
! void glUniform1i(GLint location, GLint v0);
subroutine glUniform1i(location,v0) bind(c,name="glUniform1i")
  import c_int, c_float
  integer(c_int), value :: location
  integer(c_int), value :: v0
end subroutine glUniform1i
! void glUniform1f(GLint location, GLfloat v0)
subroutine glUniform1f(location,v0) bind(c,name="glUniform1f")
  import c_int, c_float
  integer(c_int), value :: location
  real(c_float), value :: v0
end subroutine glUniform1f
! void glUniform3fv(GLint location, GLsizei count, const GLfloat *value);
subroutine glUniform3fv(location,count,val) bind(c,name="glUniform3fv")
  import c_int, c_float, c_ptr
  integer(c_int), value :: location
  integer(c_int), value :: count
  type(c_ptr), intent(in), value :: val
end subroutine glUniform3fv
! void glUniform4fv(GLint location, GLsizei count, const GLfloat *value);
subroutine glUniform4fv(location,count,val) bind(c,name="glUniform4fv")
  import c_int, c_float, c_ptr
  integer(c_int), value :: location
  integer(c_int), value :: count
  type(c_ptr), intent(in), value :: val
end subroutine glUniform4fv
! void glUniformMatrix3fv(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
subroutine glUniformMatrix3fv(location,count,transpose,val) bind(c,name="glUniformMatrix3fv")
  import c_int, c_signed_char, c_ptr
  integer(c_int), value :: location
  integer(c_int), value :: count
  integer(c_signed_char), value :: transpose
  type(c_ptr), intent(in), value :: val
end subroutine glUniformMatrix3fv
! void glUniformMatrix4fv(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
subroutine glUniformMatrix4fv(location,count,transpose,val) bind(c,name="glUniformMatrix4fv")
  import c_int, c_signed_char, c_ptr
  integer(c_int), value :: location
  integer(c_int), value :: count
  integer(c_signed_char), value :: transpose
  type(c_ptr), intent(in), value :: val
end subroutine glUniformMatrix4fv
! void glUseProgram(GLuint program);
subroutine glUseProgram(prog_) bind(c,name="glUseProgram")
  import c_int
  integer(c_int), value :: prog_
end subroutine glUseProgram
! void glVertexAttribPointer(GLuint index, GLint size, GLenum type, GLboolean normalized,
!                            GLsizei stride, const GLvoid * pointer);
subroutine glVertexAttribPointer(index,size,type,normalized,stride,point) bind(c,name="glVertexAttribPointer")
  import c_int, c_signed_char, c_ptr
  integer(c_int), value :: index
  integer(c_int), value :: size
  integer(c_int), value :: type
  integer(c_signed_char), value :: normalized
  integer(c_int), value :: stride
  type(c_ptr), intent(in), value :: point
end subroutine glVertexAttribPointer
! void glViewport (GLint x, GLint y, GLsizei width, GLsizei height);
subroutine glViewport(x, y, width, height) bind(c,name="glViewport")
  import c_int
  integer(c_int), value :: x, y, width, height
end subroutine glViewport
