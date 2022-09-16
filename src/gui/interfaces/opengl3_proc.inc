! -*-F90-*-

!xx! OpenGL3 procedures

! -- Type equivalence chart --
! GLenum = GLuint = GLint = GLsizei = c_int
! GLbitfield = c_int (unsigned int)
! GLfloat = c_float

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
! void glDeleteFramebuffers(GLsizei n, GLuint *framebuffers);
subroutine glDeleteFramebuffers(n,framebuffers) bind(c,name="glDeleteFramebuffers")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: framebuffers
end subroutine glDeleteFramebuffers
! void glDeleteRenderbuffers(GLsizei n, GLuint *renderbuffers);
subroutine glDeleteRenderbuffers(n,renderbuffers) bind(c,name="glDeleteRenderbuffers")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: renderbuffers
end subroutine glDeleteRenderbuffers
! void glDeleteTextures(GLsizei n, const GLuint *textures);
subroutine glDeleteTextures(n,textures) bind(c,name="glDeleteTextures")
  import c_int, c_ptr
  integer(c_int), value :: n
  type(c_ptr), intent(in), value :: textures
end subroutine glDeleteTextures
! void glEnable(GLenum cap)
subroutine glEnable(cap) bind(c,name="glEnable")
  import c_int
  integer(c_int), value :: cap
end subroutine glEnable
! void glFramebufferRenderbuffer(GLenum target, GLenum attachment, GLenum renderbuffertarget,
!                                GLuint renderbuffer);
subroutine glFramebufferRenderbuffer(target_,attachment,renderbuffertarget,renderbuffer) bind(c,name="glFramebufferRenderbuffer")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: attachment
  integer(c_int), value :: renderbuffertarget
  integer(c_int), value :: renderbuffer
end subroutine glFramebufferRenderbuffer
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
! void glRenderbufferStorage(GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
subroutine glRenderbufferStorage(target_,internalformat,width,height) bind(c,name="glRenderbufferStorage")
  import c_int
  integer(c_int), value :: target_
  integer(c_int), value :: internalformat
  integer(c_int), value :: width
  integer(c_int), value :: height
end subroutine glRenderbufferStorage
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
! void glViewport (GLint x, GLint y, GLsizei width, GLsizei height);
subroutine glViewport(x, y, width, height) bind(c,name="glViewport")
  import c_int
  integer(c_int), value :: x, y, width, height
end subroutine glViewport
