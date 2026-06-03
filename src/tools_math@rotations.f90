! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Utilities for representing and converting 3D rotations: rotation
! matrices, quaternions, etc. Quaternion convention throughout: q =
! (w,x,y,z), scalar part first, unit norm, active right-handed rotation.
submodule (tools_math) rotations
  implicit none

contains

  !> Convert a proper rotation matrix r (det=+1) into the corresponding
  !> unit quaternion q = (w,x,y,z) (scalar part first), using Shepperd's
  !> method (the branch with the largest pivot is chosen for numerical
  !> stability, so all rotations including 180 degrees are handled). The
  !> sign is fixed deterministically so that the first significant
  !> component in (w,x,y,z) order is positive, removing the q/-q
  !> double-cover ambiguity. Inverse of quat2mat.
  module function mat2quat(r) result(q)
    real*8, intent(in) :: r(3,3)
    real*8 :: q(4)

    integer :: i
    real*8 :: tr, s, nn

    tr = r(1,1) + r(2,2) + r(3,3)
    if (tr > r(1,1) .and. tr > r(2,2) .and. tr > r(3,3)) then
       s = sqrt(tr + 1d0) * 2d0       ! s = 4*w
       q(1) = 0.25d0 * s
       q(2) = (r(3,2) - r(2,3)) / s
       q(3) = (r(1,3) - r(3,1)) / s
       q(4) = (r(2,1) - r(1,2)) / s
    else if (r(1,1) > r(2,2) .and. r(1,1) > r(3,3)) then
       s = sqrt(1d0 + r(1,1) - r(2,2) - r(3,3)) * 2d0   ! s = 4*x
       q(1) = (r(3,2) - r(2,3)) / s
       q(2) = 0.25d0 * s
       q(3) = (r(1,2) + r(2,1)) / s
       q(4) = (r(1,3) + r(3,1)) / s
    else if (r(2,2) > r(3,3)) then
       s = sqrt(1d0 + r(2,2) - r(1,1) - r(3,3)) * 2d0   ! s = 4*y
       q(1) = (r(1,3) - r(3,1)) / s
       q(2) = (r(1,2) + r(2,1)) / s
       q(3) = 0.25d0 * s
       q(4) = (r(2,3) + r(3,2)) / s
    else
       s = sqrt(1d0 + r(3,3) - r(1,1) - r(2,2)) * 2d0   ! s = 4*z
       q(1) = (r(2,1) - r(1,2)) / s
       q(2) = (r(1,3) + r(3,1)) / s
       q(3) = (r(2,3) + r(3,2)) / s
       q(4) = 0.25d0 * s
    end if

    ! normalize
    nn = sqrt(dot_product(q,q))
    if (nn > 1d-40) q = q / nn

    ! deterministic sign: first significant component (w,x,y,z) positive
    do i = 1, 4
       if (abs(q(i)) > 1d-12) then
          if (q(i) < 0d0) q = -q
          exit
       end if
    end do

  end function mat2quat

  !> Convert a quaternion q = (w,x,y,z) (scalar part first) into the
  !> corresponding rotation matrix. The quaternion is normalized
  !> internally; a (near) zero quaternion yields the identity. Inverse of
  !> mat2quat.
  module function quat2mat(q) result(r)
    real*8, intent(in) :: q(4)
    real*8 :: r(3,3)

    real*8 :: qq(4), w, x, y, z, nn

    nn = sqrt(dot_product(q,q))
    if (nn < 1d-40) then
       r = 0d0
       r(1,1) = 1d0
       r(2,2) = 1d0
       r(3,3) = 1d0
       return
    end if
    qq = q / nn
    w = qq(1); x = qq(2); y = qq(3); z = qq(4)

    r(1,1) = 1d0 - 2d0*(y*y + z*z)
    r(1,2) = 2d0*(x*y - w*z)
    r(1,3) = 2d0*(x*z + w*y)
    r(2,1) = 2d0*(x*y + w*z)
    r(2,2) = 1d0 - 2d0*(x*x + z*z)
    r(2,3) = 2d0*(y*z - w*x)
    r(3,1) = 2d0*(x*z - w*y)
    r(3,2) = 2d0*(y*z + w*x)
    r(3,3) = 1d0 - 2d0*(x*x + y*y)

  end function quat2mat

  !> Convert a proper rotation matrix r (det=+1) into the ZYZ intrinsic
  !> Euler angles e = (alpha,beta,gamma) in radians, such that
  !> r = Rz(alpha) Ry(beta) Rz(gamma) (active, right-handed rotations).
  !> The ranges are alpha,gamma in (-pi,pi] and beta in [0,pi]. At the
  !> gimbal-lock poles (beta = 0 or pi, where only alpha+/-gamma is
  !> determined) the convention is gamma = 0. Inverse of euler2mat.
  module function mat2euler(r) result(e)
    real*8, intent(in) :: r(3,3)
    real*8 :: e(3)

    real*8, parameter :: eps = 1d-12
    real*8 :: sb

    sb = sqrt(r(1,3)*r(1,3) + r(2,3)*r(2,3))   ! |sin(beta)|
    e(2) = atan2(sb, r(3,3))                    ! beta in [0,pi]
    if (sb > eps) then
       e(1) = atan2(r(2,3), r(1,3))             ! alpha
       e(3) = atan2(r(3,2), -r(3,1))            ! gamma
    else
       ! gimbal lock: beta = 0 or pi, only alpha (+/-) gamma is defined
       e(3) = 0d0
       if (r(3,3) > 0d0) then
          e(1) = atan2(r(2,1), r(1,1))          ! beta = 0:  R = Rz(alpha)
       else
          e(1) = atan2(-r(2,1), -r(1,1))        ! beta = pi: R = Rz(alpha) Ry(pi)
       end if
    end if

  end function mat2euler

  !> Convert ZYZ intrinsic Euler angles e = (alpha,beta,gamma) (radians)
  !> into the corresponding rotation matrix r = Rz(alpha) Ry(beta)
  !> Rz(gamma) (active, right-handed). Inverse of mat2euler.
  module function euler2mat(e) result(r)
    real*8, intent(in) :: e(3)
    real*8 :: r(3,3)

    real*8 :: ca, sa, cb, sb, cg, sg

    ca = cos(e(1)); sa = sin(e(1))
    cb = cos(e(2)); sb = sin(e(2))
    cg = cos(e(3)); sg = sin(e(3))

    r(1,1) = ca*cb*cg - sa*sg
    r(1,2) = -ca*cb*sg - sa*cg
    r(1,3) = ca*sb
    r(2,1) = sa*cb*cg + ca*sg
    r(2,2) = -sa*cb*sg + ca*cg
    r(2,3) = sa*sb
    r(3,1) = -sb*cg
    r(3,2) = sb*sg
    r(3,3) = cb

  end function euler2mat

end submodule rotations
