! -*-F90-*-

!xx! OpenGL3 procedures

! -- Type equivalence chart --
! GLint = GLsizei = c_int
! GLbitfield = c_int (unsigned int)
! GLfloat = c_float

! void glViewport (GLint x, GLint y, GLsizei width, GLsizei height);
subroutine glViewport(x, y, width, height) bind(c,name="glViewport")
  import c_int
  integer(c_int), value :: x, y, width, height
end subroutine glViewport
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
