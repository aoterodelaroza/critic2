! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! OpenGL shader programs
submodule (shaders) proc
  use iso_c_binding
  implicit none

  character*(*), parameter :: shader_file(shader_NUM) = (/&
     "test        ",& ! shader_test
     "simple      ",& ! shader_simple
     "phong       ",& ! shader_phong
     "text_direct ",& ! shader_text_direct
     "text_onscene",& ! shader_text_onscene
     "pickindex   "&  ! shader_text_onscene
     /)

  ! shader programs
  integer(c_int) :: ishad_prog(shader_NUM)

  ! current shader
  integer :: icur_shader = 0

contains

  !> Compile and link the shader programs
  module subroutine shaders_init()
    use interfaces_opengl3
    use global, only: critic_home
    use tools_io, only: fopen_read, fclose, getline_raw, ferror, faterr, uout
    use param, only: dirsep, newline

    integer :: i, j, lu
    integer(c_int) :: shad(2), success, length
    character(len=:), allocatable :: file, line
    logical :: exist
    character(kind=c_char,len=:), allocatable, target :: sprog
    character(kind=c_char,len=1024), target :: infolog

    character*2, parameter :: prefix(2) = (/"vs","fs"/)
    integer(c_int) :: shtype(2)

    shtype = (/GL_VERTEX_SHADER,GL_FRAGMENT_SHADER/)
    do i = 1, shader_NUM
       ! compile
       do j = 1, 2
          ! read the vertex shader file
          file = trim(critic_home) // dirsep // "shaders" // dirsep // trim(shader_file(i)) // "." // prefix(j)
          inquire(file=file,exist=exist)
          if (.not.exist) &
             call ferror('shaders_init','could not fine shader file: ' // file,faterr)
          sprog = ""
          lu = fopen_read(file)
          do while (getline_raw(lu,line))
             sprog = sprog // line // newline
          end do
          call fclose(lu)
          sprog = sprog // c_null_char

          ! compile the shader
          shad(j) = glCreateShader(shtype(j))
          call glShaderSource(shad(j), 1_c_int, c_loc(sprog), c_null_ptr)
          call glCompileShader(shad(j))
          call glGetShaderiv(shad(j), GL_COMPILE_STATUS, success)
          if (success == GL_FALSE) then
             call glGetShaderInfoLog(shad(j), 1023, length, c_loc(infolog))
             write (uout,'("Error: ", A)') trim(infolog)
             call ferror('shaders_init','error compiling shader',faterr)
          end if
       end do

       ! create the program
       ishad_prog(i) = glCreateProgram()
       call glAttachShader(ishad_prog(i), shad(1))
       call glAttachShader(ishad_prog(i), shad(2))
       call glLinkProgram(ishad_prog(i))
       call glGetProgramiv(ishad_prog(i), GL_LINK_STATUS, success)
       if (success == GL_FALSE) then
          call glGetProgramInfoLog(ishad_prog(i), 1023, length, c_loc(infolog))
          write (uout,'("Error: ", A)') trim(infolog)
          call ferror('shaders_init','error linking shader',faterr)
       end if

       ! detach and delete the shaders
       do j = 1, 2
          call glDetachShader(ishad_prog(i),shad(j))
          call glDeleteShader(shad(j))
       end do
    end do

  end subroutine shaders_init

  !> Delete the shader programs
  module subroutine shaders_end()
    use interfaces_opengl3
    integer :: i

    do i = 1, shader_NUM
       call glDeleteProgram(ishad_prog(i))
    end do

  end subroutine shaders_end

  !> Use shader with the given ID, if valid and not currently in use
  module subroutine useshader(id)
    use interfaces_opengl3, only: glUseProgram
    integer, intent(in) :: id

    if (id < 1 .or. id > shader_NUM) return

    if (id /= icur_shader) then
       call glUseProgram(ishad_prog(id))
       icur_shader = id
    end if

  end subroutine useshader

  !> Set an integer uniform.
  module subroutine setuniform_int(name,x)
    use interfaces_opengl3, only: glUniform1i, glGetUniformLocation
    character*(*), intent(in), target :: name
    integer(c_int), intent(in) :: x

    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: idx

    str = trim(name) // c_null_char
    idx = glGetUniformLocation(ishad_prog(icur_shader),c_loc(str))
    call glUniform1i(idx,x)

  end subroutine setuniform_int

  !> Set a float uniform.
  module subroutine setuniform_float(name,x)
    use interfaces_opengl3, only: glUniform1f, glGetUniformLocation
    character*(*), intent(in), target :: name
    real(c_float), intent(in) :: x

    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: idx

    str = trim(name) // c_null_char
    idx = glGetUniformLocation(ishad_prog(icur_shader),c_loc(str))
    call glUniform1f(idx,x)

  end subroutine setuniform_float

  !> Set a vec3 float uniform.
  module subroutine setuniform_vec3(name,x)
    use interfaces_opengl3, only: glUniform3fv, glGetUniformLocation
    character*(*), intent(in), target :: name
    real(c_float), intent(in), target :: x(3)

    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: idx

    str = trim(name) // c_null_char
    idx = glGetUniformLocation(ishad_prog(icur_shader),c_loc(str))
    call glUniform3fv(idx,1,c_loc(x))

  end subroutine setuniform_vec3

  !> Set a vec4 float uniform.
  module subroutine setuniform_vec4(name,x)
    use interfaces_opengl3, only: glUniform4fv, glGetUniformLocation
    character*(*), intent(in), target :: name
    real(c_float), intent(in), target :: x(4)

    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: idx

    str = trim(name) // c_null_char
    idx = glGetUniformLocation(ishad_prog(icur_shader),c_loc(str))
    call glUniform4fv(idx,1,c_loc(x))

  end subroutine setuniform_vec4

  !> Set a mat3 float uniform
  module subroutine setuniform_mat3(name,x)
    use interfaces_opengl3, only: glUniformMatrix3fv, glGetUniformLocation, GL_FALSE
    character*(*), intent(in), target :: name
    real(c_float), intent(in), target :: x(3,3)

    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: idx

    str = trim(name) // c_null_char
    idx = glGetUniformLocation(ishad_prog(icur_shader),c_loc(str))
    call glUniformMatrix3fv(idx, 1, int(GL_FALSE,c_signed_char), c_loc(x))

  end subroutine setuniform_mat3

  !> Set a mat4 float uniform
  module subroutine setuniform_mat4(name,x)
    use interfaces_opengl3, only: glUniformMatrix4fv, glGetUniformLocation, GL_FALSE
    character*(*), intent(in), target :: name
    real(c_float), intent(in), target :: x(4,4)

    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: idx

    str = trim(name) // c_null_char
    idx = glGetUniformLocation(ishad_prog(icur_shader),c_loc(str))
    call glUniformMatrix4fv(idx, 1, int(GL_FALSE,c_signed_char), c_loc(x))

  end subroutine setuniform_mat4

end submodule proc
