  ! -*-F90-*-

  !xx! GLFW constants
  integer(c_int), bind(C,name="const_GLFW_SAMPLES") :: GLFW_SAMPLES
  integer(c_int), bind(C,name="const_GLFW_CONTEXT_VERSION_MAJOR") :: GLFW_CONTEXT_VERSION_MAJOR
  integer(c_int), bind(C,name="const_GLFW_CONTEXT_VERSION_MINOR") :: GLFW_CONTEXT_VERSION_MINOR
  integer(c_int), bind(C,name="const_GLFW_OPENGL_PROFILE") :: GLFW_OPENGL_PROFILE
  integer(c_int), bind(C,name="const_GLFW_OPENGL_FORWARD_COMPAT") :: GLFW_OPENGL_FORWARD_COMPAT
  integer(c_int), bind(C,name="const_GLFW_OPENGL_CORE_PROFILE") :: GLFW_OPENGL_CORE_PROFILE
  integer(c_int), bind(C,name="const_GLFW_STICKY_KEYS") :: GLFW_STICKY_KEYS
