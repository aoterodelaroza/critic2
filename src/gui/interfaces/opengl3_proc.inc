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
! GLintptr = c_intptr_t

! void glActiveTexture(GLenum texture);
subroutine glActiveTexture(texture) bind(c,name="glActiveTexture_")
  import c_int
  integer(c_int), value :: texture
end subroutine glActiveTexture
! void glAttachShader(GLuint program, GLuint shader);
subroutine glAttachShader(prog,shader) bind(c,name="glAttachShader_")
  import c_int
  integer(c_int), value :: prog
  integer(c_int), value :: shader
end subroutine glAttachShader
! void glBindBuffer(GLenum target, GLuint buffer);
subroutine glBindBuffer(target_,buffer) bind(c,name="glBindBuffer_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: buffer
end subroutine glBindBuffer
! void glBindFramebuffer(GLenum target, GLuint framebuffer);
subroutine glBindFramebuffer(target_,framebuffer) bind(c,name="glBindFramebuffer_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: framebuffer
end subroutine glBindFramebuffer
! void glBlendEquation(GLenum mode);
subroutine glBlendEquation(mode) bind(c,name="glBlendEquation_")
  import c_int
  integer(c_int), value :: mode
end subroutine glBlendEquation
! void glBlendFunc(GLenum sfactor, GLenum dfactor);
subroutine glBlendFunc(sfactor,dfactor) bind(c,name="glBlendFunc_")
  import c_int
  integer(c_int), value :: sfactor
  integer(c_int), value :: dfactor
end subroutine glBlendFunc
! void glBlendFuncSeparate(GLenum srcRGB, GLenum dstRGB, GLenum srcAlpha, GLenum dstAlpha);
subroutine glBlendFuncSeparate(srcRGB, dstRGB, srcAlpha, dstAlpha) bind(c,name="glBlendFuncSeparate_")
  import c_int
  integer(c_int), value :: srcRGB
  integer(c_int), value :: dstRGB
  integer(c_int), value :: srcAlpha
  integer(c_int), value :: dstAlpha
end subroutine glBlendFuncSeparate
! void glBlitFramebuffer(GLint srcX0,GLint srcY0,GLint srcX1,GLint srcY1,GLint dstX0,GLint dstY0,
!                        GLint dstX1,GLint dstY1,GLbitfield mask,GLenum filter);
subroutine glBlitFramebuffer(srcX0,srcY0,srcX1,srcY1,dstX0,dstY0,dstX1,dstY1,mask,filter) bind(c,name="glBlitFramebuffer_")
  import c_int
  integer(c_int), value :: srcX0
  integer(c_int), value :: srcY0
  integer(c_int), value :: srcX1
  integer(c_int), value :: srcY1
  integer(c_int), value :: dstX0
  integer(c_int), value :: dstY0
  integer(c_int), value :: dstX1
  integer(c_int), value :: dstY1
  integer(c_int), value :: mask
  integer(c_int), value :: filter
end subroutine glBlitFramebuffer
! void glBindRenderbuffer(GLenum target, GLuint renderbuffer);
subroutine glBindRenderbuffer(target_,renderbuffer) bind(c,name="glBindRenderbuffer_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: renderbuffer
end subroutine glBindRenderbuffer
! void glBindTexture(GLenum target, GLuint texture);
subroutine glBindTexture(target_,texture) bind(c,name="glBindTexture_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: texture
end subroutine glBindTexture
! void glBufferData(GLenum target, GLsizeiptr size, const GLvoid *data, GLenum usage);
subroutine glBufferData(target_,size,data,usage) bind(c,name="glBufferData_")
  import c_int, c_int64_t, c_ptr, c_ptrdiff_t
  integer(c_int), value :: target_
  integer(c_ptrdiff_t), value :: size
  type(c_ptr), intent(in), value :: data
  integer(c_int), value :: usage
end subroutine glBufferData
! void glBufferSubData(GLenum target, GLintptr offset, GLsizeiptr size, const GLvoid *data);
subroutine glBufferSubData(target_,offset,size,data) bind(c,name="glBufferSubData_")
  import c_int, c_intptr_t, c_ptrdiff_t, c_ptr
  integer(c_int), value :: target_
  integer(c_intptr_t), value :: offset
  integer(c_ptrdiff_t), value :: size
  type(c_ptr), intent(in), value :: data
end subroutine glBufferSubData
! void glBindVertexArray(GLuint array);
subroutine glBindVertexArray(array) bind(c,name="glBindVertexArray_")
  import c_int
  integer(c_int), value :: array
end subroutine glBindVertexArray
! GLenum glCheckFramebufferStatus(GLenum target);
function glCheckFramebufferStatus(target_) bind(c,name="glCheckFramebufferStatus_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int) :: glCheckFramebufferStatus
end function glCheckFramebufferStatus
! void glClear (GLbitfield mask);
subroutine glClear (mask) bind(c,name="glClear_")
  import c_int
  integer(c_int), value :: mask
end subroutine glClear
! void glClearColor(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha)
subroutine glClearColor(red, green, blue, alpha) bind(c,name="glClearColor_")
  import c_float
  real(c_float), value :: red, green, blue, alpha
end subroutine glClearColor
! void glClearStencil(GLint s)
subroutine glClearStencil(s) bind(c,name="glClearStencil_")
  import c_int
  integer(c_int), value :: s
end subroutine glClearStencil
! void glCompileShader(GLuint shader)
subroutine glCompileShader(shader) bind(c,name="glCompileShader_")
  import c_int
  integer(c_int), value :: shader
end subroutine glCompileShader
! GLuint glCreateProgram(void);
function glCreateProgram() bind(c,name="glCreateProgram_")
  import c_int
  integer(c_int) :: glCreateProgram
end function glCreateProgram
! GLuint glCreateShader(GLenum shaderType);
function glCreateShader(shaderType) bind(c,name="glCreateShader_")
  import c_int
  integer(c_int), value :: shaderType
  integer(c_int) :: glCreateShader
end function glCreateShader
! void glDeleteBuffers(GLsizei n, const GLuint *buffers)
subroutine glDeleteBuffers(n,buffers) bind(c,name="glDeleteBuffers_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: buffers
end subroutine glDeleteBuffers
! void glDeleteFramebuffers(GLsizei n, GLuint *framebuffers);
subroutine glDeleteFramebuffers(n,framebuffers) bind(c,name="glDeleteFramebuffers_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: framebuffers
end subroutine glDeleteFramebuffers
! void glDeleteProgram(GLuint program);
subroutine glDeleteProgram(prog) bind(c,name="glDeleteProgram_")
  import c_int
  integer(c_int), value :: prog
end subroutine glDeleteProgram
! void glDeleteRenderbuffers(GLsizei n, GLuint *renderbuffers);
subroutine glDeleteRenderbuffers(n,renderbuffers) bind(c,name="glDeleteRenderbuffers_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: renderbuffers
end subroutine glDeleteRenderbuffers
! void glDeleteShader(GLuint program);
subroutine glDeleteShader(shader) bind(c,name="glDeleteShader_")
  import c_int
  integer(c_int), value :: shader
end subroutine glDeleteShader
! void glDeleteTextures(GLsizei n, const GLuint *textures);
subroutine glDeleteTextures(n,textures) bind(c,name="glDeleteTextures_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: textures
end subroutine glDeleteTextures
! void glDeleteVertexArrays(GLsizei n, const GLuint *arrays)
subroutine glDeleteVertexArrays(n,arrays) bind(c,name="glDeleteVertexArrays_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: arrays
end subroutine glDeleteVertexArrays
! void glDetachShader(GLuint program, GLuint shader)
subroutine glDetachShader(prog,shader) bind(c,name="glDetachShader_")
  import c_int
  integer(c_int), value :: prog
  integer(c_int), value :: shader
end subroutine glDetachShader
! void glDisable(GLenum cap)
subroutine glDisable(cap) bind(c,name="glDisable_")
  import c_int
  integer(c_int), value :: cap
end subroutine glDisable
! void glDrawArrays(GLenum mode, GLint first, GLsizei count);
subroutine glDrawArrays(mode,first,count) bind(c,name="glDrawArrays_")
  import c_int
  integer(c_int), value :: mode
  integer(c_int), value :: first
  integer(c_int), value :: count
end subroutine glDrawArrays
! void glDrawElements(GLenum mode, GLsizei count, GLenum type, const GLvoid *indices);
subroutine glDrawElements(mode,count,type,indices) bind(c,name="glDrawElements_")
  import c_int, c_ptr
  integer(c_int), value :: mode
  integer(c_int), value :: count
  integer(c_int), value :: type
  type(c_ptr), intent(in), value :: indices
end subroutine glDrawElements
! void glEnable(GLenum cap)
subroutine glEnable(cap) bind(c,name="glEnable_")
  import c_int
  integer(c_int), value :: cap
end subroutine glEnable
! void glEnableVertexAttribArray(GLuint index);
subroutine glEnableVertexAttribArray(index) bind(c,name="glEnableVertexAttribArray_")
  import c_int
  integer(c_int), value :: index
end subroutine glEnableVertexAttribArray
! void glFramebufferRenderbuffer(GLenum target, GLenum attachment, GLenum renderbuffertarget,
!                                GLuint renderbuffer);
subroutine glFramebufferRenderbuffer(target_,attachment,renderbuffertarget,renderbuffer) bind(c,name="glFramebufferRenderbuffer_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: attachment
  integer(c_int), value :: renderbuffertarget
  integer(c_int), value :: renderbuffer
end subroutine glFramebufferRenderbuffer
! void glGenBuffers(GLsizei n, GLuint* buffers);
subroutine glGenBuffers(n,buffers) bind(c,name="glGenBuffers_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), value :: buffers
end subroutine glGenBuffers
! void glGenerateMipmap(GLenum target);
subroutine glGenerateMipmap(targt) bind(c,name="glGenerateMipmap")
  import c_int
  integer(c_int), value :: targt
end subroutine glGenerateMipmap
! void glGenFramebuffers(GLsizei n,GLuint *ids);
subroutine glGenFramebuffers(n,ids) bind(c,name="glGenFramebuffers_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), value :: ids
end subroutine glGenFramebuffers
! GLenum glGetError( void);
function glGetError() bind(c,name="glGetError_")
  import c_int
  integer(c_int) :: glGetError
end function glGetError
! void glGetProgramInfoLog(GLuint program, GLsizei maxLength, GLsizei *length, GLchar *infoLog);
subroutine glGetProgramInfoLog(prog,maxLength,length,infoLog) bind(c,name="glGetProgramInfoLog_")
  import c_int, c_ptr
  integer(c_int), value :: prog
  integer(c_int), value :: maxLength
  integer(c_int) :: length
  type(c_ptr), value :: infoLog
end subroutine glGetProgramInfoLog
! void glGetProgramiv(GLuint program, GLenum pname, GLint *params);
subroutine glGetProgramiv(prog,pname,params) bind(c,name="glGetProgramiv_")
  import c_int
  integer(c_int), value :: prog
  integer(c_int), value :: pname
  integer(c_int) :: params
end subroutine glGetProgramiv
! void glGetShaderInfoLog(GLuint shader, GLsizei maxLength, GLsizei *length, GLchar *infoLog);
subroutine glGetShaderInfoLog(shader,maxLength,length,infoLog) bind(c,name="glGetShaderInfoLog_")
  import c_int, c_ptr
  integer(c_int), value :: shader
  integer(c_int), value :: maxLength
  integer(c_int) :: length
  type(c_ptr), value :: infoLog
end subroutine glGetShaderInfoLog
! void glGetShaderiv(GLuint shader, GLenum pname, GLint *params);
subroutine glGetShaderiv(shader,pname,params) bind(c,name="glGetShaderiv_")
  import c_int
  integer(c_int), value :: shader
  integer(c_int), value :: pname
  integer(c_int) :: params
end subroutine glGetShaderiv
! GLint glGetUniformLocation(GLuint program, const GLchar *name);
function glGetUniformLocation(prog_,name) bind(c,name="glGetUniformLocation_")
  import c_int, c_ptr
  integer(c_int), value :: prog_
  type(c_ptr), intent(in), value :: name
  integer(c_int) :: glGetUniformLocation
end function glGetUniformLocation
! void glFramebufferTexture2D(GLenum target, GLenum attachment, GLenum textarget, GLuint texture,
!      GLint level);
subroutine glFramebufferTexture2D(target_,attachment,textarget,texture,level) bind(c,name="glFramebufferTexture2D_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: attachment
  integer(c_int), value :: textarget
  integer(c_int), value :: texture
  integer(c_int), value :: level
end subroutine glFramebufferTexture2D
! void glGenRenderbuffers(GLsizei n,GLuint *renderbuffers)
subroutine glGenRenderbuffers(n,renderbuffers) bind(c,name="glGenRenderbuffers_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), value :: renderbuffers
end subroutine glGenRenderbuffers
! void glGenTextures(GLsizei n,GLuint *textures);
subroutine glGenTextures(n,textures) bind(c,name="glGenTextures_")
  import c_ptr, c_int
  integer(c_int), value :: n
  type(c_ptr), value :: textures
end subroutine glGenTextures
! void glGenVertexArrays(GLsizei n, GLuint *arrays);
subroutine glGenVertexArrays(n,arrays) bind(c,name="glGenVertexArrays_")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), value :: arrays
end subroutine glGenVertexArrays
! void glLinkProgram(GLuint program)
subroutine glLinkProgram(prog) bind(c,name="glLinkProgram_")
  import c_int
  integer(c_int), value :: prog
end subroutine glLinkProgram
! void glPixelStorei(GLenum pname,GLint param)
subroutine glPixelStorei(pname,param) bind(c,name="glPixelStorei_")
  import c_int
  integer(c_int), value :: pname, param
end subroutine glPixelStorei
! void glReadPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid * data);
subroutine glReadPixels(x,y,width,height,form,type,data) bind(c,name="glReadPixels_")
  import c_int, c_ptr
  integer(c_int), value :: x, y
  integer(c_int), value :: width, height
  integer(c_int), value :: form
  integer(c_int), value :: type
  type(c_ptr), intent(in), value :: data
end subroutine glReadPixels
! void glRenderbufferStorage(GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
subroutine glRenderbufferStorage(target_,internalformat,width,height) bind(c,name="glRenderbufferStorage_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: internalformat
  integer(c_int), value :: width
  integer(c_int), value :: height
end subroutine glRenderbufferStorage
! void glRenderbufferStorageMultisample(GLenum target,GLsizei samples,GLenum internalformat,GLsizei width,GLsizei height);
subroutine glRenderbufferStorageMultisample(target_,samples,internalformat,width,height) &
   bind(c,name="glRenderbufferStorageMultisample_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: samples
  integer(c_int), value :: internalformat
  integer(c_int), value :: width
  integer(c_int), value :: height
end subroutine glRenderbufferStorageMultisample
! void glShaderSource(GLuint shader, GLsizei count, const GLchar **string, const GLint *length);
subroutine glShaderSource(shader,count,string,length) bind(c,name="glShaderSource_")
  import c_int, c_ptr
  integer(c_int), value :: shader
  integer(c_int), value :: count
  type(c_ptr), intent(in) :: string
  type(c_ptr), intent(in), value :: length
end subroutine glShaderSource
! void glStencilFunc(GLenum func, GLint ref, GLuint mask);
subroutine glStencilFunc(func, ref, mask) bind(c,name="glStencilFunc_")
  import c_int
  integer(c_int), value :: func
  integer(c_int), value :: ref
  integer(c_int), value :: mask
end subroutine glStencilFunc
! void glStencilMask(GLuint mask);
subroutine glStencilMask(mask) bind(c,name="glStencilMask_")
  import c_int
  integer(c_int), value :: mask
end subroutine glStencilMask
! void glStencilOp(GLenum sfail, GLenum dpfail, GLenum dppass);
subroutine glStencilOp(sfail, dpfail, dppass) bind(c,name="glStencilOp_")
  import c_int
  integer(c_int), value :: sfail
  integer(c_int), value :: dpfail
  integer(c_int), value :: dppass
end subroutine glStencilOp
! void glTexImage2D(GLenum target, GLint level, GLint internalFormat, GLsizei width, GLsizei height,
!                   GLint border, GLenum format, GLenum type, const GLvoid * data);
subroutine glTexImage2D(target_,level,internalFormat,width,height,border,format_,type,data) bind(c,name="glTexImage2D_")
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
! void glTexImage2DMultisample(GLenum target,GLsizei samples,GLint internalformat,GLsizei width,
!                              GLsizei height,GLboolean fixedsamplelocations);
subroutine glTexImage2DMultisample(target_,samples,internalformat,width,height,fixedsamplelocations) &
   bind(c,name="glTexImage2DMultisample_")
  import c_int, c_signed_char
  integer(c_int), value :: target_
  integer(c_int), value :: samples
  integer(c_int), value :: internalformat
  integer(c_int), value :: width
  integer(c_int), value :: height
  integer(c_signed_char), value :: fixedsamplelocations
end subroutine glTexImage2DMultisample
! void glTexParameterf(GLenum target_, GLenum pname, GLfloat param);
subroutine glTexParameterf(target_,pname,param) bind(c,name="glTexParameterf_")
  import c_int, c_float
  integer(c_int), value :: target_
  integer(c_int), value :: pname
  real(c_float), value :: param
end subroutine glTexParameterf
! void glTexParameteri(GLenum target_, GLenum pname, GLint param);
subroutine glTexParameteri(target_,pname,param) bind(c,name="glTexParameteri_")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: pname
  integer(c_int), value :: param
end subroutine glTexParameteri
! void glUniform1i(GLint location, GLint v0);
subroutine glUniform1i(location,v0) bind(c,name="glUniform1i_")
  import c_int, c_float
  integer(c_int), value :: location
  integer(c_int), value :: v0
end subroutine glUniform1i
! void glUniform1f(GLint location, GLfloat v0)
subroutine glUniform1f(location,v0) bind(c,name="glUniform1f_")
  import c_int, c_float
  integer(c_int), value :: location
  real(c_float), value :: v0
end subroutine glUniform1f
! void glUniform3fv(GLint location, GLsizei count, const GLfloat *value);
subroutine glUniform3fv(location,count,val) bind(c,name="glUniform3fv_")
  import c_int, c_float, c_ptr
  integer(c_int), value :: location
  integer(c_int), value :: count
  type(c_ptr), intent(in), value :: val
end subroutine glUniform3fv
! void glUniform4fv(GLint location, GLsizei count, const GLfloat *value);
subroutine glUniform4fv(location,count,val) bind(c,name="glUniform4fv_")
  import c_int, c_float, c_ptr
  integer(c_int), value :: location
  integer(c_int), value :: count
  type(c_ptr), intent(in), value :: val
end subroutine glUniform4fv
! void glUniformMatrix3fv(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
subroutine glUniformMatrix3fv(location,count,transpose,val) bind(c,name="glUniformMatrix3fv_")
  import c_int, c_signed_char, c_ptr
  integer(c_int), value :: location
  integer(c_int), value :: count
  integer(c_signed_char), value :: transpose
  type(c_ptr), intent(in), value :: val
end subroutine glUniformMatrix3fv
! void glUniformMatrix4fv(GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
subroutine glUniformMatrix4fv(location,count,transpose,val) bind(c,name="glUniformMatrix4fv_")
  import c_int, c_signed_char, c_ptr
  integer(c_int), value :: location
  integer(c_int), value :: count
  integer(c_signed_char), value :: transpose
  type(c_ptr), intent(in), value :: val
end subroutine glUniformMatrix4fv
! void glUseProgram(GLuint program);
subroutine glUseProgram(prog_) bind(c,name="glUseProgram_")
  import c_int
  integer(c_int), value :: prog_
end subroutine glUseProgram
! void glVertexAttribPointer(GLuint index, GLint size, GLenum type, GLboolean normalized,
!                            GLsizei stride, const GLvoid * pointer);
subroutine glVertexAttribPointer(index,size,type,normalized,stride,point) bind(c,name="glVertexAttribPointer_")
  import c_int, c_signed_char, c_ptr
  integer(c_int), value :: index
  integer(c_int), value :: size
  integer(c_int), value :: type
  integer(c_signed_char), value :: normalized
  integer(c_int), value :: stride
  type(c_ptr), intent(in), value :: point
end subroutine glVertexAttribPointer
! void glViewport (GLint x, GLint y, GLsizei width, GLsizei height);
subroutine glViewport(x, y, width, height) bind(c,name="glViewport_")
  import c_int
  integer(c_int), value :: x, y, width, height
end subroutine glViewport
