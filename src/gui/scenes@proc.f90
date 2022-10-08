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

! Scene object and GL rendering utilities
submodule (scenes) proc
  implicit none

  real*8 :: min_scenerad = 10.d0

  ! object resolution
  integer, parameter :: isphres = 3 ! sphere
  integer, parameter :: icylres = 3 ! cylinder

  ! some math parameters
  real(c_float), parameter :: zero = 0._c_float
  real(c_float), parameter :: one = 1._c_float
  real(c_float), parameter :: eye4(4,4) = reshape((/&
     one,zero,zero,zero,&
     zero,one,zero,zero,&
     zero,zero,one,zero,&
     zero,zero,zero,one/),shape(eye4))

contains

  !> Initialize a scene object associated with system isys.
  module subroutine scene_init(s,isys)
    use gui_main, only: nsys, sysc, sys_init, sys
    class(scene), intent(inout), target :: s
    integer, intent(in) :: isys

    real*8 :: xmin(3), xmax(3), x(3)

    if (isys < 1 .or. isys > nsys) return
    if (sysc(isys)%status /= sys_init) return

    ! basic variables
    s%id = isys
    s%showcell = .not.(sys(isys)%c%ismolecule)
    s%ncell = 1
    s%isinit = .true.

    ! default imotif
    if (sys(isys)%c%ismolecule) then
       s%imotif = 0
    elseif (sys(isys)%c%ismol3d) then
       s%imotif = 3
    else
       s%imotif = 1
    end if

    ! phong default settings
    s%lightpos = (/20._c_float,20._c_float,0._c_float/)
    s%lightcolor = (/1._c_float,1._c_float,1._c_float/)
    s%ambient = 0.2_c_float
    s%diffuse = 0.4_c_float
    s%specular = 0.6_c_float
    s%shininess = 8_c_int

    ! reset the camera
    call s%reset()

  end subroutine scene_init

  !> Terminate a scene object
  module subroutine scene_end(s)
    class(scene), intent(inout), target :: s

    s%isinit = .false.
    s%id = 0

  end subroutine scene_end

  !> Reset the camera position and direction. Sets scenerad,
  !> scenecenter, ortho_fov, persp_fov, v_front, v_up, v_pos, view,
  !> world, projection, and znear.
  module subroutine scene_reset(s)
    use utils, only: lookat
    use gui_main, only: sys
    use param, only: pi
    class(scene), intent(inout), target :: s

    real(c_float) :: pic, hw2
    real*8 :: xmin(3), xmax(3), x(3), center(3)
    integer :: isys
    integer :: i, i1, i2, i3

    ! initialize
    isys = s%id

    ! calculate the initial scene radius
    xmin = 0d0
    xmax = 0d0
    if (sys(isys)%c%ncel > 0) then
       xmin = sys(isys)%c%atcel(1)%r
       xmax = sys(isys)%c%atcel(1)%r
       do i = 2, sys(isys)%c%ncel
          xmax = max(sys(isys)%c%atcel(i)%r,xmax)
          xmin = min(sys(isys)%c%atcel(i)%r,xmin)
       end do

       if (.not.sys(isys)%c%ismolecule) then
          do i1 = 0, 1
             do i2 = 0, 1
                do i3 = 0, 1
                   x = real((/i1,i2,i3/),8)
                   x = sys(isys)%c%x2c(x)
                   xmax = max(x,xmax)
                   xmin = min(x,xmin)
                end do
             end do
          end do
       end if
    end if
    s%scenerad = max(norm2(xmax-xmin),min_scenerad)

    ! calculate the scene center
    center = (/0.5d0,0.5d0,0.5d0/)
    center = sys(isys)%c%x2c(center)

    ! default transformation matrices
    pic = real(pi,c_float)
    s%ortho_fov = 45._c_float
    s%persp_fov = 45._c_float

    s%world = eye4

    s%v_front = (/0._c_float,0._c_float,-1._c_float/)
    s%v_up = (/0._c_float,1._c_float,0._c_float/)
    s%v_pos = real(center,c_float)
    s%v_pos(3) = s%v_pos(3) + 1.1_c_float * s%scenerad * &
       sqrt(real(maxval(s%ncell),c_float)) / tan(0.5_c_float * s%ortho_fov * pic / 180._c_float)
    s%view = lookat(s%v_pos,s%v_pos+s%v_front,s%v_up)

    s%znear = 0.1_c_float
    call s%update_projection_matrix()

  end subroutine scene_reset

  !> Draw the scene
  module subroutine scene_render(s)
    use interfaces_opengl3
    use gui_main, only: sysc, sys_init, sys
    use shaders, only: shader_test, shader_phong, useshader, setuniform_int,&
       setuniform_float, setuniform_vec3, setuniform_vec4, setuniform_mat3,&
       setuniform_mat4
    use shapes, only: sphVAO, cylVAO
    use param, only: jmlcol2, atmcov
    class(scene), intent(inout), target :: s

    real(c_float), target :: vcolor(4), model(4,4), rgba(4), x0(3), x1(3), rad
    integer :: i, iz

    real*8, parameter :: uc(3,2,12) = reshape((/&
       0._c_float,0._c_float,0._c_float,  1._c_float,0._c_float,0._c_float,&
       0._c_float,0._c_float,0._c_float,  0._c_float,1._c_float,0._c_float,&
       0._c_float,0._c_float,0._c_float,  0._c_float,0._c_float,1._c_float,&
       1._c_float,1._c_float,1._c_float,  1._c_float,1._c_float,0._c_float,&
       1._c_float,1._c_float,1._c_float,  1._c_float,0._c_float,1._c_float,&
       1._c_float,1._c_float,1._c_float,  0._c_float,1._c_float,1._c_float,&
       1._c_float,0._c_float,0._c_float,  1._c_float,1._c_float,0._c_float,&
       0._c_float,1._c_float,0._c_float,  1._c_float,1._c_float,0._c_float,&
       0._c_float,1._c_float,0._c_float,  0._c_float,1._c_float,1._c_float,&
       0._c_float,0._c_float,1._c_float,  0._c_float,1._c_float,1._c_float,&
       0._c_float,0._c_float,1._c_float,  1._c_float,0._c_float,1._c_float,&
       1._c_float,0._c_float,0._c_float,  1._c_float,0._c_float,1._c_float/),shape(uc))

    if (.not.s%isinit) return
    if (.not.sysc(s%id)%status == sys_init) return

    ! set up the shader and the uniforms
    call useshader(shader_phong)
    call setuniform_int("uselighting",1_c_int)
    call setuniform_vec3("lightPos",s%lightpos)
    call setuniform_vec3("lightColor",s%lightcolor)
    call setuniform_float("ambient",s%ambient)
    call setuniform_float("diffuse",s%diffuse)
    call setuniform_float("specular",s%specular)
    call setuniform_int("shininess",s%shininess)
    call setuniform_mat4("world",s%world)
    call setuniform_mat4("view",s%view)
    call setuniform_mat4("projection",s%projection)

    ! draw the atoms
    call glBindVertexArray(sphVAO(isphres))
    do i = 1, sys(s%id)%c%ncel
       iz = sys(s%id)%c%spc(sys(s%id)%c%atcel(i)%is)%z
       rgba(1:3) = real(jmlcol2(:,iz),c_float) / 255._c_float
       rgba(4) = 1._c_float
       x0 = real(sys(s%id)%c%atcel(i)%r,c_float)
       rad = 0.7_c_float * real(atmcov(iz),c_float)
       call s%draw_sphere(x0,rad,rgba)
    end do
    call glBindVertexArray(0)

    ! draw the unit cell
    if (s%showcell) then
       call glBindVertexArray(cylVAO(icylres))
       rad = 0.2_c_float
       rgba = (/1._c_float,0._c_float,0._c_float,1._c_float/)
       do i = 1, 12
          x0 = real(sys(s%id)%c%x2c(uc(:,1,i)),c_float)
          x1 = real(sys(s%id)%c%x2c(uc(:,2,i)),c_float)
          call s%draw_cylinder(x0,x1,rad,rgba)
       end do
       call glBindVertexArray(0)
    end if

  end subroutine scene_render

  !> Update the projection matrix from the v_pos
  module subroutine update_projection_matrix(s)
    use utils, only: ortho
    use param, only: pi
    class(scene), intent(inout), target :: s

    real(c_float) :: pic, hw2

    pic = real(pi,c_float)
    hw2 = tan(0.5_c_float * s%ortho_fov * pic / 180._c_float) * s%v_pos(3)
    s%projection = ortho(-hw2,hw2,-hw2,hw2,s%znear,1000._c_float)

  end subroutine update_projection_matrix

  !> Draw a sphere with center x0, radius rad and color rgba. Requires
  !> having the sphere VAO bound.
  module subroutine draw_sphere(s,x0,rad,rgba)
    use interfaces_opengl3
    use shaders, only: setuniform_vec4, setuniform_mat4
    use shapes, only: sphnel
    class(scene), intent(inout), target :: s
    real(c_float), intent(in) :: x0(3)
    real(c_float), intent(in) :: rad
    real(c_float), intent(in) :: rgba(4)

    real(c_float) :: model(4,4)

    ! the model matrix: scale and translate
    model = eye4
    model(1:3,4) = x0
    model(1,1) = rad
    model(2,2) = rad
    model(3,3) = rad

    ! draw the sphere
    call setuniform_vec4("vColor",rgba)
    call setuniform_mat4("model",model)
    call glDrawElements(GL_TRIANGLES, int(3*sphnel(isphres),c_int), GL_UNSIGNED_INT, c_null_ptr)

  end subroutine draw_sphere

  !> Draw a cylinder from x1 to x2 with radius rad and color
  !> rgba. Requires having the cylinder VAO bound.
  module subroutine draw_cylinder(s,x1,x2,rad,rgba)
    use interfaces_opengl3
    use utils, only: cross_c
    use shaders, only: setuniform_vec4, setuniform_mat4
    use shapes, only: cylnel
    class(scene), intent(inout), target :: s
    real(c_float), intent(in) :: x1(3)
    real(c_float), intent(in) :: x2(3)
    real(c_float), intent(in) :: rad
    real(c_float), intent(in) :: rgba(4)

    real(c_float) :: xmid(3), xdif(3), up(3), crs(3), model(4,4), blen
    real(c_float) :: a, ca, sa, axis(3), temp(3)

    xmid = 0.5_c_float * (x1 + x2)
    xdif = x2 - x1
    blen = norm2(xdif)
    xdif = xdif / blen
    up = (/0._c_float,0._c_float,1._c_float/)
    crs = cross_c(up,xdif)

    ! the model matrix
    model = eye4
    model(1:3,4) = xmid ! translate(m_model,xmid);
    if (dot_product(crs,crs) > 1e-14_c_float) then
       ! m_model = m_model * rotate(acos(dot(xdif,up)),crs);
       a = acos(dot_product(xdif,up))
       ca = cos(a)
       sa = sin(a)
       axis = crs / norm2(crs)
       temp = (1._c_float - ca) * axis

       model(1,1) = ca + temp(1) * axis(1)
       model(2,1) = temp(1) * axis(2) + sa * axis(3)
       model(3,1) = temp(1) * axis(3) - sa * axis(2)

       model(1,2) = temp(2) * axis(1) - sa * axis(3)
       model(2,2) = ca + temp(2) * axis(2)
       model(3,2) = temp(2) * axis(3) + sa * axis(1)

       model(1,3) = temp(3) * axis(1) + sa * axis(2)
       model(2,3) = temp(3) * axis(2) - sa * axis(1)
       model(3,3) = ca + temp(3) * axis(3)
    end if
    ! m_model = scale(m_model,glm::vec3(rad,rad,blen));
    model(:,1) = model(:,1) * rad
    model(:,2) = model(:,2) * rad
    model(:,3) = model(:,3) * blen

    ! draw the cylinder
    call setuniform_vec4("vColor",rgba)
    call setuniform_mat4("model",model)
    call glDrawElements(GL_TRIANGLES, int(3*cylnel(icylres),c_int), GL_UNSIGNED_INT, c_null_ptr)

  end subroutine draw_cylinder

end submodule proc
