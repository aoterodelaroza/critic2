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

! Some math utilities adapted from the glm library
submodule (utils) math
  use iso_c_binding
  implicit none

contains

  !> Return the projection matrix for an infinite perspective with
  !> field of view fovy, aspect ratio aspect, and near plane znear
  module subroutine infiniteperspective(m,fovy,aspect,znear)
    real(c_float), intent(out) :: m(4,4)
    real(c_float), intent(in) :: fovy, aspect, znear

    real(c_float) :: range, left, right, bottom, top

    range = tan(fovy / 2._c_float) * znear
    left = -range * aspect
    right = range * aspect
    bottom = -range
    top = range

    m = 0._c_float
    m(1,1) = (2._c_float * znear) / (right - left)
    m(2,2) = (2._c_float * znear) / (top - bottom)
    m(3,3) = -1._c_float
    m(4,3) = -1._c_float
    m(3,4) = -2._c_float * znear

  end subroutine infiniteperspective

  !> Return the projection matrix for an orthographic perspective with
  !> fustrum limits left, right, bottom, top and clipping planes znear
  !> and zfar.
  module subroutine ortho(m,left,right,bottom,top,znear,zfar)
    real(c_float), intent(out) :: m(4,4)
    real(c_float), intent(in) :: left, right, bottom, top, znear, zfar

    m = 0._c_float
    m(1,1) = 2._c_float / (right - left)
    m(2,2) = 2._c_float / (top - bottom)
    m(3,3) = -2._c_float / (zfar - znear)
    m(1,4) = - (right + left) / (right - left)
    m(2,4) = - (top + bottom) / (top - bottom)
    m(3,4) = - (zfar + znear) / (zfar - znear)
    m(4,4) = 1._c_float

  end subroutine ortho

  !> Camera (view) transformation matrix for a camera at position eye
  !> pointing in the directon towards center and with the given up vector.
  module subroutine lookat(m,eye,center,up)
    use tools_math, only: cross_cfloat
    real(c_float), intent(out) :: m(4,4)
    real(c_float), intent(in) :: eye(3)
    real(c_float), intent(in) :: center(3)
    real(c_float), intent(in) :: up(3)

    real(c_float) :: f(3), s(3), u(3)

    f = center - eye
    f = f / norm2(f)
    s = cross_cfloat(f,up)
    s = s / norm2(s)
    u = cross_cfloat(s,f)
    u = u / norm2(u)
    m(:,1) = (/ s(1), u(1), -f(1), 0._c_float/)
    m(:,2) = (/ s(2), u(2), -f(2), 0._c_float/)
    m(:,3) = (/ s(3), u(3), -f(3), 0._c_float/)
    m(:,4) = (/-dot_product(s,eye),-dot_product(u,eye),dot_product(f,eye),1._c_float/)

  end subroutine lookat

  !> From 3D coordinates in pos (world, view, etc.) to window
  !> coordinates. model is the [world]*view*[model] matrix, proj
  !> is the projection matrix. viewport_a is the side of the current
  !> (square) viewport (0,0,a,a).
  module subroutine project(m,pos,mview,proj,viewport_a)
    real(c_float), intent(out) :: m(3)
    real(c_float), intent(in) :: pos(3)
    real(c_float), intent(in) :: mview(4,4)
    real(c_float), intent(in) :: proj(4,4)
    integer(c_int), intent(in) :: viewport_a

    real(c_float) :: tmp(4)

    tmp(1:3) = pos
    tmp(4) = 1._c_float
    tmp = matmul(mview,tmp)
    tmp = matmul(proj,tmp)
    tmp = tmp / tmp(4)
    tmp = tmp * 0.5_c_float + 0.5_c_float
    tmp(1) = tmp(1) * real(viewport_a,c_float)
    tmp(2) = tmp(2) * real(viewport_a,c_float)
    m = tmp(1:3)

  end subroutine project

  !> From 3D coordinates in window coordinates to world/view/etc.
  !> coordinates. model is the [world]*view*[model] matrix, proj
  !> is the projection matrix. viewport_a is the side of the current
  !> (square) viewport (0,0,a,a).
  module subroutine unproject(m,pos,mview,proj,viewport_a)
    use tools_math, only: matinv_cfloat
    real(c_float), intent(out) :: m(3)
    real(c_float), intent(in) :: pos(3)
    real(c_float), intent(in) :: mview(4,4)
    real(c_float), intent(in) :: proj(4,4)
    integer(c_int), intent(in) :: viewport_a

    real(c_float) :: tmp(4), minv(4,4)
    integer :: ier

    ! invert the projection * mview matrix
    minv = matmul(proj,mview)
    call matinv_cfloat(minv,4,ier)

    ! unproject
    tmp(1:3) = pos
    tmp(4) = 1._c_float
    tmp(1) = tmp(1) / real(viewport_a,c_float)
    tmp(2) = tmp(2) / real(viewport_a,c_float)
    tmp = tmp * 2._c_float - 1._c_float
    tmp = matmul(minv,tmp)
    tmp = tmp / tmp(4)
    m = tmp(1:3)

  end subroutine unproject

  ! Apply a translation matrix to m by vector v. Return the translated matrix.
  module subroutine translate(m,mat,v)
    real(c_float), intent(out) :: m(4,4)
    real(c_float), intent(in) :: mat(4,4)
    real(c_float), intent(in) :: v(3)

    m = mat
    m(:,4) = mat(:,1) * v(1) + mat(:,2) * v(2) + mat(:,3) * v(3) + mat(:,4)

  end subroutine translate

  ! Apply a rotation matrix to m by axis and angle. Return the rotated matrix.
  module subroutine rotate(m,mat,angle,axis)
    real(c_float), intent(out) :: m(4,4)
    real(c_float), intent(in) :: mat(4,4)
    real(c_float), intent(in) :: angle
    real(c_float), intent(in) :: axis(3)

    real(c_float) :: a, c, s, ax(3), temp(3), mm(3,3), res(4,4), nn

    nn = norm2(axis)
    if (nn < 1d-20) then
       m = mat
       return
    end if

    a = angle
    c = cos(angle)
    s = sin(angle)
    ax = axis / nn
    temp = (1-c) * ax

    mm(1,1) = c + temp(1) * ax(1)
    mm(2,1) = temp(1) * ax(2) + s * ax(3)
    mm(3,1) = temp(1) * ax(3) - s * ax(2)

    mm(1,2) = temp(2) * ax(1) - s * ax(3)
    mm(2,2) = c + temp(2) * ax(2)
    mm(3,2) = temp(2) * ax(3) + s * ax(1)

    mm(1,3) = temp(3) * ax(1) + s * ax(2)
    mm(2,3) = temp(3) * ax(2) - s * ax(1)
    mm(3,3) = c + temp(3) * ax(3)

    res(:,1) = mat(:,1) * mm(1,1) + mat(:,2) * mm(2,1) + mat(:,3) * mm(3,1)
    res(:,2) = mat(:,1) * mm(1,2) + mat(:,2) * mm(2,2) + mat(:,3) * mm(3,2)
    res(:,3) = mat(:,1) * mm(1,3) + mat(:,2) * mm(2,3) + mat(:,3) * mm(3,3)
    res(:,4) = mat(:,4)
    m = res

  end subroutine rotate

  ! Calculate m * v + t, where t is the translation in the 4x4 matrix.
  ! Returns the resulting 3-vector.
  module subroutine mult(m,mat,v,notrans)
    real(c_float), intent(out) :: m(3)
    real(c_float), intent(in) :: mat(4,4)
    real(c_float), intent(in) :: v(3)
    logical, intent(in), optional :: notrans

    logical :: notrans_
    real(c_float) :: vx(4), m3(3,3)

    notrans_ = .false.
    if (present(notrans)) notrans_ = notrans

    if (notrans_) then
       m3 = mat(1:3,1:3)
       m = matmul(m3,v)
    else
       vx(1:3) = v
       vx(4) = 1._c_float
       vx = matmul(mat,vx)
       m = vx(1:3) / vx(4)
    end if

  end subroutine mult

  ! Calculate inv(m) * v + t, where t is the translation in the
  ! 4x4 matrix. Returns the resulting 3-vector.
  module subroutine invmult(m,mat,v,notrans)
    use tools_math, only: matinv_cfloat
    real(c_float), intent(out) :: m(3)
    real(c_float), intent(in) :: mat(4,4)
    real(c_float), intent(in) :: v(3)
    logical, intent(in), optional :: notrans

    integer :: ier
    logical :: notrans_
    real(c_float) :: vx(4), mx(4,4), m3(3,3)

    notrans_ = .false.
    if (present(notrans)) notrans_ = notrans

    mx = mat
    call matinv_cfloat(mx,4,ier)

    if (notrans_) then
       m3 = mx(1:3,1:3)
       vx(1:3) = matmul(m3,v)
       m = vx(1:3)
    else
       vx(1:3) = v
       vx(4) = 1._c_float
       vx = matmul(mx,vx)
       m = vx(1:3) / vx(4)
    end if

  end subroutine invmult

end submodule math
