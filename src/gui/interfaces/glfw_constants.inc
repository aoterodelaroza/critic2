! -*-F90-*-

!xx! GLFW constants
integer(c_int), parameter :: GLFW_SAMPLES                = int(z'0002100D',c_int)
integer(c_int), parameter :: GLFW_CONTEXT_VERSION_MAJOR  = int(z'00022002',c_int)
integer(c_int), parameter :: GLFW_CONTEXT_VERSION_MINOR  = int(z'00022003',c_int)
integer(c_int), parameter :: GLFW_OPENGL_PROFILE         = int(z'00022008',c_int)
integer(c_int), parameter :: GLFW_OPENGL_FORWARD_COMPAT  = int(z'00022006',c_int)
integer(c_int), parameter :: GLFW_OPENGL_CORE_PROFILE    = int(z'00032001',c_int)
integer(c_int), parameter :: GLFW_STICKY_KEYS            = int(z'00033002',c_int)
