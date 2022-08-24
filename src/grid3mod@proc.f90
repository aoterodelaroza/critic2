! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Class for 3d grids and related tools.
submodule (grid3mod) proc
  implicit none

  !xx! private procedures
  ! subroutine grinterp_nearest(f,x0,y)
  ! subroutine grinterp_trilinear(f,x0,y,yp)
  ! subroutine grinterp_trispline(f,x0,y,yp,ypp)
  ! subroutine grinterp_tricubic(f,xi,y,yp,ypp)
  ! subroutine grinterp_smr(f,xi,y,yp,ypp,i0ref)
  ! function grid_near(f,x) result(res)
  ! function grid_floor(f,x,main,shift)
  ! function grid_ceiling(f,x,main,shift)
  ! function euclidean_near(f,x)
  ! subrotine init_geometry(f,x2c,n,env)
  ! subrotine copy_geometry(f,g)
  ! subroutine init_trispline(f)
  ! subroutine init_smr(f)
  ! subroutine smr_kernelfun(r,k,f,fp,fpp)

  ! Notes about the 3D-FFT order of the array elements.
  ! - cfftnd accepts a 1D array representing the 3D array in Fortran
  !   memory order (i.e. reshape can be used).
  !
  ! - An element of this array is:
  !     igfft = j3*n(2)*n(1)+j2*n(1)+j1+1
  !   where jx: 1 -> nx.
  !
  ! - The jx: 1 -> nx are recentered to ix: -nx+nx/2+1 -> nx/2 via
  !   the operation:
  !   jx = modulo(ix,n(x))
  !   i.e. the jx > nx/2 have nx subtracted from them to make ix.
  !
  ! - The reciprocal-space lattice vector for a given set of ix is:
  !     vgc = i1 * (a*) + i2 * (b*) + i3 * (c*)
  !   where
  !     a* = 2 * pi * cross(c,b) / vol
  !     b* = 2 * pi * cross(a,c) / vol
  !     c* = 2 * pi * cross(b,a) / vol
  !     with vol = abs(det3(x2c)) = a * cross(b,c)
  !
  ! - To run over all ix and jx (ntot = n1 * n2 * n3):
  !
  !   do ig = 1, ntot
  !     i3 = mod(ig-1,n(3)) + n(3)/2 - n(3) + 1
  !     iaux = (ig-1 - (i3-1 + n(3) - n(3)/2)) / n(3)
  !     i2 = mod(iaux,n(2)) + n(2)/2 - n(2) + 1
  !     iaux = (iaux - (i2-1 + n(2) - n(2)/2)) / n(2)
  !     i1 = iaux + n(1)/2 - n(1) + 1
  !     j1 = modulo(i1,n(1))
  !     j2 = modulo(i2,n(2))
  !     j3 = modulo(i3,n(3))
  !   end do
  !
  ! - To transform the ix back to the ig index:
  !   ig = i3-1+n(3)-n(3)/2 + (i2-1+n(2)-n(2)/2) * n(3) + (i1-1+n(1)-n(1)/2) * n(3) * n(2) + 1
  !

  ! The 64x64 matrix for tricubic interpolation
  real*8, parameter :: c(64,64) = reshape((/&                      ! values for c(i,j), with...  (i,  j)
     1d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&  ! 1-16, 1
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 1
     -3d0,  0d0,  9d0, -6d0,  0d0,  0d0,  0d0,  0d0,  9d0,  0d0, -27d0,  18d0, -6d0,  0d0,  18d0, -12d0,& ! 33-48, 1
     2d0,  0d0, -6d0,  4d0,  0d0,  0d0,  0d0,  0d0, -6d0,  0d0,  18d0, -12d0,  4d0,  0d0, -12d0,   8d0,&  ! 49-64, 1
     0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&  ! 1-16, 2
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 2
     0d0,  0d0, -9d0,  6d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  27d0, -18d0,  0d0,  0d0, -18d0,  12d0,&  ! 33-48, 2
     0d0,  0d0,  6d0, -4d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -18d0,  12d0,  0d0,  0d0,  12d0,  -8d0,&  ! 49-64, 2
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&  ! 1-16, 3
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 3
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -9d0,  0d0,  27d0, -18d0,  6d0,  0d0, -18d0,  12d0,&  ! 33-48, 3
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0,  0d0, -18d0,  12d0, -4d0,  0d0,  12d0,  -8d0,&  ! 49-64, 3
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&  ! 1-16, 4
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 4
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -27d0,  18d0,  0d0,  0d0,  18d0, -12d0,&  ! 33-48, 4
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  18d0, -12d0,  0d0,  0d0, -12d0,   8d0,&  ! 49-64, 4
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 5
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 5
     3d0,  0d0, -9d0,  6d0,  0d0,  0d0,  0d0,  0d0, -9d0,  0d0,  27d0, -18d0,  6d0,  0d0, -18d0,  12d0,&  ! 33-48, 5
     -2d0,  0d0,  6d0, -4d0,  0d0,  0d0,  0d0,  0d0,  6d0,  0d0, -18d0,  12d0, -4d0,  0d0,  12d0,  -8d0,& ! 49-64, 5
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 6
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 6
     0d0,  0d0,  9d0, -6d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -27d0,  18d0,  0d0,  0d0,  18d0, -12d0,&  ! 33-48, 6
     0d0,  0d0, -6d0,  4d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  18d0, -12d0,  0d0,  0d0, -12d0,   8d0,&  ! 49-64, 6
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 7
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 7
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  9d0,  0d0, -27d0,  18d0, -6d0,  0d0,  18d0, -12d0,&  ! 33-48, 7
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  0d0,  18d0, -12d0,  4d0,  0d0, -12d0,   8d0,&  ! 49-64, 7
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 8
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 8
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  27d0, -18d0,  0d0,  0d0, -18d0,  12d0,&  ! 33-48, 8
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -18d0,  12d0,  0d0,  0d0,  12d0,  -8d0,&  ! 49-64, 8
     0d0,  1d0, -2d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&  ! 1-16, 9
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 9
     0d0, -3d0,  6d0, -3d0,  0d0,  0d0,  0d0,  0d0,  0d0,  9d0, -18d0,   9d0,  0d0, -6d0,  12d0,  -6d0,&  ! 33-48, 9
     0d0,  2d0, -4d0,  2d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  12d0,  -6d0,  0d0,  4d0,  -8d0,   4d0,&  ! 49-64, 9
     0d0,  0d0, -1d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&  ! 1-16, 10
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 10
     0d0,  0d0,  3d0, -3d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   9d0,  0d0,  0d0,   6d0,  -6d0,&  ! 33-48, 10
     0d0,  0d0, -2d0,  2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -6d0,  0d0,  0d0,  -4d0,   4d0,&  ! 49-64, 10
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&  ! 1-16, 11
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 11
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -9d0,  18d0,  -9d0,  0d0,  6d0, -12d0,   6d0,&  ! 33-48, 11
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0, -12d0,   6d0,  0d0, -4d0,   8d0,  -4d0,&  ! 49-64, 11
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&  ! 1-16, 12
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 12
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -9d0,  0d0,  0d0,  -6d0,   6d0,&  ! 33-48, 12
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -6d0,   6d0,  0d0,  0d0,   4d0,  -4d0,&  ! 49-64, 12
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 13
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 13
     0d0,  3d0, -6d0,  3d0,  0d0,  0d0,  0d0,  0d0,  0d0, -9d0,  18d0,  -9d0,  0d0,  6d0, -12d0,   6d0,&  ! 33-48, 13
     0d0, -2d0,  4d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0, -12d0,   6d0,  0d0, -4d0,   8d0,  -4d0,&  ! 49-64, 13
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 14
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 14
     0d0,  0d0, -3d0,  3d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -9d0,  0d0,  0d0,  -6d0,   6d0,&  ! 33-48, 14
     0d0,  0d0,  2d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -6d0,   6d0,  0d0,  0d0,   4d0,  -4d0,&  ! 49-64, 14
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 15
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 15
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  9d0, -18d0,   9d0,  0d0, -6d0,  12d0,  -6d0,&  ! 33-48, 15
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  12d0,  -6d0,  0d0,  4d0,  -8d0,   4d0,&  ! 49-64, 15
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 16
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 16
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   9d0,  0d0,  0d0,   6d0,  -6d0,&  ! 33-48, 16
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -6d0,  0d0,  0d0,  -4d0,   4d0,&  ! 49-64, 16
     0d0,  0d0,  0d0,  0d0,  1d0,  0d0, -3d0,  2d0, -2d0,  0d0,   6d0,  -4d0,  1d0,  0d0,  -3d0,   2d0,&  ! 1-16, 17
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 17
     0d0,  0d0,  0d0,  0d0, -3d0,  0d0,  9d0, -6d0,  6d0,  0d0, -18d0,  12d0, -3d0,  0d0,   9d0,  -6d0,&  ! 33-48, 17
     0d0,  0d0,  0d0,  0d0,  2d0,  0d0, -6d0,  4d0, -4d0,  0d0,  12d0,  -8d0,  2d0,  0d0,  -6d0,   4d0,&  ! 49-64, 17
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   3d0,  -2d0,&  ! 1-16, 18
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 18
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -9d0,  6d0,  0d0,  0d0,  18d0, -12d0,  0d0,  0d0,  -9d0,   6d0,&  ! 33-48, 18
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0, -4d0,  0d0,  0d0, -12d0,   8d0,  0d0,  0d0,   6d0,  -4d0,&  ! 49-64, 18
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  0d0,   3d0,  -2d0,  1d0,  0d0,  -3d0,   2d0,&  ! 1-16, 19
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 19
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -3d0,  0d0,   9d0,  -6d0,&  ! 33-48, 19
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  0d0,   6d0,  -4d0,  2d0,  0d0,  -6d0,   4d0,&  ! 49-64, 19
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   2d0,  0d0,  0d0,   3d0,  -2d0,&  ! 1-16, 20
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 20
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -9d0,   6d0,&  ! 33-48, 20
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   6d0,  -4d0,&  ! 49-64, 20
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 21
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 21
     0d0,  0d0,  0d0,  0d0,  3d0,  0d0, -9d0,  6d0, -6d0,  0d0,  18d0, -12d0,  3d0,  0d0,  -9d0,   6d0,&  ! 33-48, 21
     0d0,  0d0,  0d0,  0d0, -2d0,  0d0,  6d0, -4d0,  4d0,  0d0, -12d0,   8d0, -2d0,  0d0,   6d0,  -4d0,&  ! 49-64, 21
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 22
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 22
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  9d0, -6d0,  0d0,  0d0, -18d0,  12d0,  0d0,  0d0,   9d0,  -6d0,&  ! 33-48, 22
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  4d0,  0d0,  0d0,  12d0,  -8d0,  0d0,  0d0,  -6d0,   4d0,&  ! 49-64, 22
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 23
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 23
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  3d0,  0d0,  -9d0,   6d0,&  ! 33-48, 23
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0,  0d0,  -6d0,   4d0, -2d0,  0d0,   6d0,  -4d0,&  ! 49-64, 23
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 24
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 24
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   9d0,  -6d0,&  ! 33-48, 24
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -4d0,  0d0,  0d0,  -6d0,   4d0,&  ! 49-64, 24
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 25
     1d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&  ! 17-32, 25
     -2d0,  0d0,  6d0, -4d0,  0d0,  0d0,  0d0,  0d0,  6d0,  0d0, -18d0,  12d0, -4d0,  0d0,  12d0,  -8d0,& ! 33-48, 25
     1d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&  ! 49-64, 25
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 26
     0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&  ! 17-32, 26
     0d0,  0d0, -6d0,  4d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  18d0, -12d0,  0d0,  0d0, -12d0,   8d0,&  ! 33-48, 26
     0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&  ! 49-64, 26
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 27
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&  ! 17-32, 27
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  0d0,  18d0, -12d0,  4d0,  0d0, -12d0,   8d0,&  ! 33-48, 27
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&  ! 49-64, 27
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 28
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&  ! 17-32, 28
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -18d0,  12d0,  0d0,  0d0,  12d0,  -8d0,&  ! 33-48, 28
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&  ! 49-64, 28
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 29
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 29
     -1d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,& ! 33-48, 29
     1d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&  ! 49-64, 29
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 30
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 30
     0d0,  0d0, -3d0,  2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&  ! 33-48, 30
     0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&  ! 49-64, 30
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 31
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 31
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  0d0,   9d0,  -6d0,  2d0,  0d0,  -6d0,   4d0,&  ! 33-48, 31
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  0d0,  -9d0,   6d0, -2d0,  0d0,   6d0,  -4d0,&  ! 49-64, 31
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 32
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 32
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -9d0,   6d0,  0d0,  0d0,   6d0,  -4d0,&  ! 33-48, 32
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   9d0,  -6d0,  0d0,  0d0,  -6d0,   4d0,&  ! 49-64, 32
     0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -2d0,  1d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  1d0,  -2d0,   1d0,&  ! 1-16, 33
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 33
     0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  6d0, -3d0,  0d0,  6d0, -12d0,   6d0,  0d0, -3d0,   6d0,  -3d0,&  ! 33-48, 33
     0d0,  0d0,  0d0,  0d0,  0d0,  2d0, -4d0,  2d0,  0d0, -4d0,   8d0,  -4d0,  0d0,  2d0,  -4d0,   2d0,&  ! 49-64, 33
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  1d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -1d0,   1d0,&  ! 1-16, 34
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 34
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -3d0,  0d0,  0d0,  -6d0,   6d0,  0d0,  0d0,   3d0,  -3d0,&  ! 33-48, 34
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  2d0,  0d0,  0d0,   4d0,  -4d0,  0d0,  0d0,  -2d0,   2d0,&  ! 49-64, 34
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,   2d0,  -1d0,  0d0,  1d0,  -2d0,   1d0,&  ! 1-16, 35
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 35
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -3d0,   6d0,  -3d0,&  ! 33-48, 35
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  2d0,  -4d0,   2d0,&  ! 49-64, 35
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   1d0,  -1d0,  0d0,  0d0,  -1d0,   1d0,&  ! 1-16, 36
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 36
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   3d0,  -3d0,&  ! 33-48, 36
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -2d0,   2d0,&  ! 49-64, 36
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 37
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 37
     0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -6d0,  3d0,  0d0, -6d0,  12d0,  -6d0,  0d0,  3d0,  -6d0,   3d0,&  ! 33-48, 37
     0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  4d0, -2d0,  0d0,  4d0,  -8d0,   4d0,  0d0, -2d0,   4d0,  -2d0,&  ! 49-64, 37
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 38
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 38
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  3d0,  0d0,  0d0,   6d0,  -6d0,  0d0,  0d0,  -3d0,   3d0,&  ! 33-48, 38
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0, -2d0,  0d0,  0d0,  -4d0,   4d0,  0d0,  0d0,   2d0,  -2d0,&  ! 49-64, 38
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 39
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 39
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  3d0,  -6d0,   3d0,&  ! 33-48, 39
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0,  -4d0,   2d0,  0d0, -2d0,   4d0,  -2d0,&  ! 49-64, 39
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 40
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 40
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -3d0,   3d0,&  ! 33-48, 40
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -2d0,   2d0,  0d0,  0d0,   2d0,  -2d0,&  ! 49-64, 40
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 41
     0d0,  1d0, -2d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&  ! 17-32, 41
     0d0, -2d0,  4d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  6d0, -12d0,   6d0,  0d0, -4d0,   8d0,  -4d0,&  ! 33-48, 41
     0d0,  1d0, -2d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&  ! 49-64, 41
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 42
     0d0,  0d0, -1d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&  ! 17-32, 42
     0d0,  0d0,  2d0, -2d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -6d0,   6d0,  0d0,  0d0,   4d0,  -4d0,&  ! 33-48, 42
     0d0,  0d0, -1d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&  ! 49-64, 42
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 43
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&  ! 17-32, 43
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  12d0,  -6d0,  0d0,  4d0,  -8d0,   4d0,&  ! 33-48, 43
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&  ! 49-64, 43
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 44
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&  ! 17-32, 44
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -6d0,  0d0,  0d0,  -4d0,   4d0,&  ! 33-48, 44
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&  ! 49-64, 44
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 45
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 45
     0d0, -1d0,  2d0, -1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&  ! 33-48, 45
     0d0,  1d0, -2d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&  ! 49-64, 45
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 46
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 46
     0d0,  0d0,  1d0, -1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&  ! 33-48, 46
     0d0,  0d0, -1d0,  1d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&  ! 49-64, 46
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 47
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 47
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,   6d0,  -3d0,  0d0,  2d0,  -4d0,   2d0,&  ! 33-48, 47
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0,  -6d0,   3d0,  0d0, -2d0,   4d0,  -2d0,&  ! 49-64, 47
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 48
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 48
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -3d0,  0d0,  0d0,  -2d0,   2d0,&  ! 33-48, 48
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   3d0,  0d0,  0d0,   2d0,  -2d0,&  ! 49-64, 48
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 49
     0d0,  0d0,  0d0,  0d0,  1d0,  0d0, -3d0,  2d0, -2d0,  0d0,   6d0,  -4d0,  1d0,  0d0,  -3d0,   2d0,&  ! 17-32, 49
     0d0,  0d0,  0d0,  0d0, -2d0,  0d0,  6d0, -4d0,  4d0,  0d0, -12d0,   8d0, -2d0,  0d0,   6d0,  -4d0,&  ! 33-48, 49
     0d0,  0d0,  0d0,  0d0,  1d0,  0d0, -3d0,  2d0, -2d0,  0d0,   6d0,  -4d0,  1d0,  0d0,  -3d0,   2d0,&  ! 49-64, 49
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 50
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   3d0,  -2d0,&  ! 17-32, 50
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -6d0,  4d0,  0d0,  0d0,  12d0,  -8d0,  0d0,  0d0,  -6d0,   4d0,&  ! 33-48, 50
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   3d0,  -2d0,&  ! 49-64, 50
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 51
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  0d0,   3d0,  -2d0,  1d0,  0d0,  -3d0,   2d0,&  ! 17-32, 51
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0,  0d0,  -6d0,   4d0, -2d0,  0d0,   6d0,  -4d0,&  ! 33-48, 51
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  0d0,   3d0,  -2d0,  1d0,  0d0,  -3d0,   2d0,&  ! 49-64, 51
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 52
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   2d0,  0d0,  0d0,   3d0,  -2d0,&  ! 17-32, 52
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   6d0,  -4d0,  0d0,  0d0,  -6d0,   4d0,&  ! 33-48, 52
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   2d0,  0d0,  0d0,   3d0,  -2d0,&  ! 49-64, 52
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 53
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 53
     0d0,  0d0,  0d0,  0d0, -1d0,  0d0,  3d0, -2d0,  2d0,  0d0,  -6d0,   4d0, -1d0,  0d0,   3d0,  -2d0,&  ! 33-48, 53
     0d0,  0d0,  0d0,  0d0,  1d0,  0d0, -3d0,  2d0, -2d0,  0d0,   6d0,  -4d0,  1d0,  0d0,  -3d0,   2d0,&  ! 49-64, 53
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 54
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 54
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -3d0,  2d0,  0d0,  0d0,   6d0,  -4d0,  0d0,  0d0,  -3d0,   2d0,&  ! 33-48, 54
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  3d0, -2d0,  0d0,  0d0,  -6d0,   4d0,  0d0,  0d0,   3d0,  -2d0,&  ! 49-64, 54
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 55
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 55
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  1d0,  0d0,  -3d0,   2d0, -1d0,  0d0,   3d0,  -2d0,&  ! 33-48, 55
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  0d0,   3d0,  -2d0,  1d0,  0d0,  -3d0,   2d0,&  ! 49-64, 55
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 56
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 56
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   3d0,  -2d0,  0d0,  0d0,  -3d0,   2d0,&  ! 33-48, 56
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -3d0,   2d0,  0d0,  0d0,   3d0,  -2d0,&  ! 49-64, 56
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 57
     0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -2d0,  1d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  1d0,  -2d0,   1d0,&  ! 17-32, 57
     0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  4d0, -2d0,  0d0,  4d0,  -8d0,   4d0,  0d0, -2d0,   4d0,  -2d0,&  ! 33-48, 57
     0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -2d0,  1d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  1d0,  -2d0,   1d0,&  ! 49-64, 57
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 58
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  1d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -1d0,   1d0,&  ! 17-32, 58
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0, -2d0,  0d0,  0d0,  -4d0,   4d0,  0d0,  0d0,   2d0,  -2d0,&  ! 33-48, 58
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  1d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -1d0,   1d0,&  ! 49-64, 58
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 59
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,   2d0,  -1d0,  0d0,  1d0,  -2d0,   1d0,&  ! 17-32, 59
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  2d0,  -4d0,   2d0,  0d0, -2d0,   4d0,  -2d0,&  ! 33-48, 59
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,   2d0,  -1d0,  0d0,  1d0,  -2d0,   1d0,&  ! 49-64, 59
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 60
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   1d0,  -1d0,  0d0,  0d0,  -1d0,   1d0,&  ! 17-32, 60
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -2d0,   2d0,  0d0,  0d0,   2d0,  -2d0,&  ! 33-48, 60
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   1d0,  -1d0,  0d0,  0d0,  -1d0,   1d0,&  ! 49-64, 60
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 61
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 61
     0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  2d0, -1d0,  0d0,  2d0,  -4d0,   2d0,  0d0, -1d0,   2d0,  -1d0,&  ! 33-48, 61
     0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -2d0,  1d0,  0d0, -2d0,   4d0,  -2d0,  0d0,  1d0,  -2d0,   1d0,&  ! 49-64, 61
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 62
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 62
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  1d0, -1d0,  0d0,  0d0,  -2d0,   2d0,  0d0,  0d0,   1d0,  -1d0,&  ! 33-48, 62
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  1d0,  0d0,  0d0,   2d0,  -2d0,  0d0,  0d0,  -1d0,   1d0,&  ! 49-64, 62
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 63
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 63
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  1d0,  -2d0,   1d0,  0d0, -1d0,   2d0,  -1d0,&  ! 33-48, 63
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,   2d0,  -1d0,  0d0,  1d0,  -2d0,   1d0,&  ! 49-64, 63
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 1-16, 64
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   0d0,   0d0,  0d0,  0d0,   0d0,   0d0,&  ! 17-32, 64
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  -1d0,   1d0,  0d0,  0d0,   1d0,  -1d0,&  ! 33-48, 64
     0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,  0d0,   1d0,  -1d0,  0d0,  0d0,  -1d0,   1d0&   ! 49-64, 64
     /),shape(c))

  integer, parameter :: smr_kkern = 3

contains

  !> Build a 3d grid of dimension n using an arithmetic expression
  !> (expr). sptr = C pointer to the associated system.
  module subroutine new_eval(f,sptr,n,expr,x2c,env)
    use tools_math, only: matinv
    use arithmetic, only: eval_grid
    use types, only: realloc
    use iso_c_binding, only: c_ptr
    class(grid3), intent(inout) :: f
    type(c_ptr), intent(in) :: sptr
    integer, intent(in) :: n(3)
    character(*), intent(in) :: expr
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env

    logical :: iok

    call f%end()
    call init_geometry(f,x2c,n,env)
    f%mode = mode_default
    f%isinit = .true.
    allocate(f%f(n(1),n(2),n(3)))
    call eval_grid(n,expr,sptr,f%f,iok)
    if (.not.iok) call f%end()

  end subroutine new_eval

  module subroutine grid_end(f)
    class(grid3), intent(inout) :: f

    f%mode = mode_default
    f%isinit = .false.
    f%isqe = .false.
    f%iswan = .false.
    f%smr_nlist = 0
    if (allocated(f%f)) deallocate(f%f)
    if (allocated(f%c2)) deallocate(f%c2)
    if (allocated(f%smr_ilist)) deallocate(f%smr_ilist)
    if (allocated(f%smr_xlist)) deallocate(f%smr_xlist)
    if (allocated(f%smr_phiinv)) deallocate(f%smr_phiinv)
    if (allocated(f%smr_rho0)) deallocate(f%smr_rho0)
    if (allocated(f%qe%kpt)) deallocate(f%qe%kpt)
    if (allocated(f%qe%wk)) deallocate(f%qe%wk)
    if (allocated(f%qe%ek)) deallocate(f%qe%ek)
    if (allocated(f%qe%occ)) deallocate(f%qe%occ)
    if (allocated(f%qe%ngk)) deallocate(f%qe%ngk)
    if (allocated(f%qe%igk_k)) deallocate(f%qe%igk_k)
    if (allocated(f%qe%nl)) deallocate(f%qe%nl)
    if (allocated(f%qe%nlm)) deallocate(f%qe%nlm)
    if (allocated(f%qe%center)) deallocate(f%qe%center)
    if (allocated(f%qe%spread)) deallocate(f%qe%spread)
    if (allocated(f%qe%u)) deallocate(f%qe%u)

  end subroutine grid_end

  !> Set the interpolation mode. The possible modes arenearest,
  !> trilinear, trispline, and tricubic (lowercase).
  module subroutine setmode(f,mode)
    use tools_io, only: equal, lower
    use param, only: icrd_crys, VSMALL
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: mode

    character(len=:), allocatable :: lmode
    integer :: i, j, k
    real*8 :: xdelta(3,3), x(3), rho, rhof(3), rhoff(3,3)

    ! parse the mode
    lmode = lower(mode)
    if (equal(lmode,'smoothrho')) then
       f%mode = mode_smr
       f%smr_nlist = 0
       f%smr_nenv = 8**3
       f%smr_fdmax = 1.75d0
       if (allocated(f%smr_ilist)) deallocate(f%smr_ilist)
       if (allocated(f%smr_xlist)) deallocate(f%smr_xlist)
       if (allocated(f%smr_phiinv)) deallocate(f%smr_phiinv)
       if (allocated(f%smr_rho0)) deallocate(f%smr_rho0)
    else if (equal(lmode,'tricubic')) then
       f%mode = mode_tricubic
    else if (equal(lmode,'trispline')) then
       f%mode = mode_trispline
       if (allocated(f%c2)) deallocate(f%c2)
    else if (equal(lmode,'trilinear')) then
       f%mode = mode_trilinear
    else if (equal(lmode,'nearest')) then
       f%mode = mode_nearest
    else if (equal(lmode,'default')) then
       f%mode = mode_default
    end if

    ! if test interpolation, calculate the rho0 and derivatives here
    if (f%mode == mode_smr) then
       if (allocated(f%smr_rho0)) deallocate(f%smr_rho0)
       allocate(f%smr_rho0(f%n(1),f%n(2),f%n(3)))

       do i = 1, 3
          xdelta(:,i) = 0d0
          xdelta(i,i) = 1d0 / real(f%n(i),8)
       end do
       !$omp parallel do private(x,rho,rhof,rhoff)
       do k = 1, f%n(3)
          do j = 1, f%n(2)
             do i = 1, f%n(1)
                x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)
                call f%atenv%promolecular(x,icrd_crys,rho,rhof,rhoff,0)

                !$omp critical(write)
                f%smr_rho0(i,j,k) = max(rho,VSMALL)
                !$omp end critical(write)
             end do
          end do
       end do
       !$omp end parallel do
    end if

  end subroutine setmode

  !> Normalize the grid to a given value. omega is the cell volume.
  module subroutine normalize(f,norm,omega)
    class(grid3), intent(inout) :: f
    real*8, intent(in) :: norm, omega

    f%f = f%f / (sum(f%f) * omega / real(product(f%n),8)) * norm
    if (allocated(f%c2)) deallocate(f%c2)

  end subroutine normalize

  !> Build a grid field from a three-dimensional array
  module subroutine from_array3(f,g,x2c,env)
    use tools_math, only: matinv
    use tools_io, only: ferror, faterr
    class(grid3), intent(inout) :: f
    real*8, intent(in) :: g(:,:,:)
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env

    integer :: istat, n(3)

    call f%end()
    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default
    n = ubound(g) - lbound(g) + 1
    call init_geometry(f,x2c,n,env)
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('from_array3','Error allocating grid',faterr)
    f%f = g

  end subroutine from_array3

  !> Read a grid in Gaussian CUBE format
  module subroutine read_cube(f,file,x2c,env,ti)
    use tools_io, only: fopen_read, ferror, faterr, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    integer :: luc
    integer :: nat
    integer :: istat, n(3), i, j, k
    logical :: ismo

    call f%end()
    luc = fopen_read(file,ti=ti)

    read (luc,*)
    read (luc,*)
    read (luc,*,iostat=istat) nat
    ismo = (nat < 0)
    nat = abs(nat)

    if (istat /= 0) &
       call ferror('read_cube','Error reading nat',faterr,file)
    do i = 1, 3
       read (luc,*,iostat=istat) n(i)
       if (istat /= 0) &
          call ferror('read_cube','Error reading nx, ny, nz',faterr,file)
    end do
    do i = 1, nat
       read (luc,*)
    end do
    if (ismo) &
       read (luc,*)

    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default
    call init_geometry(f,x2c,n,env)
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_cube','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),k=1,n(3)),j=1,n(2)),i=1,n(1))
    if (istat /= 0) &
       call ferror('read_cube','Error reading grid',faterr,file)

    call fclose(luc)

  end subroutine read_cube

  !> Read a grid in binary CUBE format
  module subroutine read_bincube(f,file,x2c,env,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, ferror, faterr, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    integer :: luc
    integer :: nat
    integer :: istat, n(3), i, iz
    logical :: ismo
    real*8 :: x0(3), xd(3,3), rdum

    call f%end()
    luc = fopen_read(file,form="unformatted",ti=ti)

    read (luc,iostat=istat) nat, x0
    ismo = (nat < 0)
    nat = abs(nat)

    if (istat /= 0) &
       call ferror('read_cube','Error reading nat',faterr,file)
    read (luc,iostat=istat) n, xd
    if (istat /= 0) &
       call ferror('read_cube','Error reading nx, ny, nz',faterr,file)
    do i = 1, nat
       read (luc) iz, rdum, x0
    end do

    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default
    call init_geometry(f,x2c,n,env)
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_cube','Error allocating grid',faterr,file)
    read(luc,iostat=istat) f%f
    if (istat /= 0) &
       call ferror('read_cube','Error reading grid',faterr,file)

    call fclose(luc)

  end subroutine read_bincube

  !> Read a grid in siesta RHO format
  module subroutine read_siesta(f,file,x2c,env,ti)
    use tools_io, only: fopen_read, faterr, ferror, fclose
    use tools_math, only: matinv
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    integer :: luc, nspin, istat
    integer :: i, iy, iz, n(3)
    real*8 :: r(3,3)
    real*4, allocatable :: g(:)

    ! initialize
    call f%end()
    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default

    ! open file
    luc = fopen_read(file,'unformatted',ti=ti)

    ! assume unformatted
    read (luc) r
    read (luc) n, nspin
    call init_geometry(f,x2c,n,env)

    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_siesta','Error allocating grid',faterr,file)
    allocate(g(n(1)),stat=istat)
    if (istat /= 0) &
       call ferror('read_siesta','Error allocating auxiliary grid',faterr,file)
    f%f = 0d0
    do i = 1, nspin
       do iz = 1, n(3)
          do iy = 1, n(2)
             read (luc) g
             f%f(:,iy,iz) = f%f(:,iy,iz) + g
          end do
       end do
    end do
    deallocate(g)

    call fclose(luc)

  end subroutine read_siesta

  !> Read a grid in abinit format
  module subroutine read_abinit(f,file,x2c,env,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, ferror, faterr, fclose
    use abinit_private, only: hdr_type, hdr_io
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: errmsg
    integer :: luc
    integer :: fform0, istat, n(3)
    type(hdr_type) :: hdr
    real*8, allocatable :: g(:,:,:)

    call f%end()
    luc = fopen_read(file,'unformatted',ti=ti)

    ! read the header
    call hdr_io(fform0,hdr,1,luc,errmsg)
    if (len_trim(errmsg) > 0) &
       call ferror('read_abinit',errmsg,faterr,file)

    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default
    n = hdr%ngfft(:)
    call init_geometry(f,x2c,n,env)
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_abinit','Error allocating grid',faterr,file)
    allocate(g(n(1),n(2),n(3)))
    read(luc,iostat=istat) g
    f%f = g
    deallocate(g)
    if (istat /= 0) &
       call ferror('read_abinit','Error reading grid',faterr,file)

    call fclose(luc)

  end subroutine read_abinit

  !> Read a grid in VASP (CHGCAR/CHG) format from a file. Omega is the
  !> cell volume used to divide the grid values (the CHGCAR contains
  !> density * omega). In CHGCAR containing more than one grid block,
  !> ibl can be used to choose which block to read (density, spin
  !> density, etc.). If vscal, scale by volume.
  module subroutine read_vasp(f,file,x2c,vscal,ibl,env,ti)
    use tools_math, only: det3, matinv
    use tools_io, only: fopen_read, getline_raw, faterr, ferror, fclose, string, &
       isinteger
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: x2c(3,3)
    logical, intent(in) :: vscal
    integer, intent(in), optional :: ibl
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    integer :: luc
    integer :: istat, n(3), i, j, k, ibcur, nnew(3)
    character(len=:), allocatable :: line
    logical :: ok

    call f%end()
    luc = fopen_read(file,ti=ti)

    do while(.true.)
       ok = getline_raw(luc,line,.true.)
       if (len(trim(line)) == 0) exit
    end do

    read (luc,*,iostat=istat) n
    if (istat /= 0) &
       call ferror('read_vasp','Error reading nx, ny, nz',faterr,file)

    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default
    call init_geometry(f,x2c,n,env)
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_vasp','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('read_vasp','Error reading grid',faterr,file)

    if (present(ibl)) then
       ibcur = 1
       do while (ibcur < ibl)
          nnew = 0
          do while (all(nnew /= n))
             ok = getline_raw(luc,line,.false.)
             if (.not.ok) &
                call ferror('read_vasp','Error reading grid block '//string(ibcur)//' < requested '//string(ibl),faterr,file)
             ok = isinteger(nnew(1),line)
             ok = ok .and. isinteger(nnew(2),line)
             ok = ok .and. isinteger(nnew(3),line)
          end do

          read(luc,*,iostat=istat) (((f%f(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
          if (istat /= 0) &
             call ferror('read_vasp','Error reading grid in block '//string(ibcur),faterr,file)
          ibcur = ibcur + 1
       end do
    end if
    if (vscal) &
       f%f(:,:,:) = f%f(:,:,:) / det3(x2c)
    call fclose(luc)

  end subroutine read_vasp

  !> Read a grid in aimpac qub format
  module subroutine read_qub(f,file,x2c,env,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, ferror, faterr, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    integer :: luc
    integer :: istat, n(3), i, j, k

    call f%end()
    luc = fopen_read(file,ti=ti)

    read (luc,*,iostat=istat) n
    if (istat /= 0) &
       call ferror('read_qub','Error reading n1, n2, n3',faterr,file)

    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default
    call init_geometry(f,x2c,n,env)
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_qub','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('read_qub','Error reading grid',faterr,file)

    call fclose(luc)

  end subroutine read_qub

  !> Read a grid in xcrysden xsf format -- only first 3d grid in first 3d block
  module subroutine read_xsf(f,file,x2c,env,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, getline_raw, lgetword, equal, ferror, faterr, &
       fclose
    use types, only: realloc
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    integer :: luc
    integer :: istat, n(3), lp, i, j, k
    character(len=:), allocatable :: line, word
    logical :: found, ok
    real*8, dimension(3) :: x0, x1, x2, x3
    real*8 :: pmat(3,3)
    real*8, allocatable :: ggloc(:,:,:)

    call f%end()

    ! open file for reading
    luc = fopen_read(file,ti=ti)

    ! position at the beginning of the first grid, ignore the rest
    found = .false.
    ok = .false.
    do while (getline_raw(luc,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,'primvec'))then
          ok = .true.
          read(luc,*,iostat=istat) pmat
          if (istat /= 0) &
             call ferror('read_xsf','Error PRIMVEC',faterr,file)
       else if (equal(word,'begin_block_datagrid_3d').or.equal(word,'begin_block_datagrid3d'))then
          found = .true.
          exit
       end if
    end do
    if (.not.found) call ferror('read_xsf','BEGIN_BLOCK_DATAGRID_3D not found',faterr,file)
    if (.not.ok) call ferror('read_xsf','PRIMVEC not found',faterr,file)

    found = .false.
    do while (getline_raw(luc,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word(1:min(17,len(word))),'begin_datagrid_3d') .or. equal(word(1:min(11,len(word))),'datagrid_3d').or.&
           equal(word(1:min(16,len(word))),'begin_datagrid3d') .or. equal(word(1:min(10,len(word))),'datagrid3d')) then
          found = .true.
          exit
       end if
    end do
    if (.not.found) call ferror('read_xsf','BEGIN_DATAGRID_3D... not found',faterr,file)

    ! grid dimension
    read (luc,*,iostat=istat) n
    if (istat /= 0) &
       call ferror('read_xsf','Error reading n1, n2, n3',faterr,file)
    call init_geometry(f,x2c,n-1,env)

    ! origin and edge vectors
    read (luc,*,iostat=istat) x0, x1, x2, x3

    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default
    allocate(ggloc(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('read_xsf','Error allocating grid',faterr,file)
    if (istat /= 0) &
       call ferror('read_xsf','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((ggloc(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('read_xsf','Error reading grid',faterr,file)

    allocate(f%f(f%n(1),f%n(2),f%n(3)),stat=istat)
    f%f = ggloc(1:n(1)-1,1:n(2)-1,1:n(3)-1)
    deallocate(ggloc)

    call fclose(luc)

  end subroutine read_xsf

  !> Read pwc file created by pw2critic.x in Quantum
  !> ESPRESSO. Contains the Bloch states, k-points, and structural
  !> info. Calculates the electron density from the Bloch states.
  !> ispin = 0 (all-electron density), 1 (spin-up), 2 (spin-down).
  !> ikpt = use only the indicated k-points. ibnd = use only the
  !> indicated bands. emin,emax: only the bands in the energy range.
  module subroutine read_pwc(f,fpwc,ispin,ikpt,ibnd,emin,emax,x2c,env,ti)
    use tools_math, only: det3, matinv
    use tools_io, only: fopen_read, fclose, ferror, faterr
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: fpwc
    integer, intent(in) :: ispin
    integer, intent(in), allocatable :: ikpt(:)
    integer, intent(in), allocatable :: ibnd(:)
    real*8, intent(in) :: emin, emax
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    integer :: i, is, ik, ib, iver, n(3)
    integer :: luc
    integer :: npwx, ngms, nkstot, nsp, nat
    real*8 :: at(3,3), fspin, alat
    complex*16, allocatable :: raux(:,:,:), rseq(:), evc(:)
    logical :: ok1, ok2

    real*8, parameter :: epsocc = 1d-6

    ! initialize
    call f%end()
    f%qe%fpwc = fpwc
    f%isqe = .true.
    f%iswan = .false.

    ! open file
    luc = fopen_read(fpwc,form="unformatted",ti=ti)

    ! header and lattice vectors
    read (luc) iver ! version
    if (iver < 2) &
       call ferror('read_pwc','This pwc file is too old. Please update your QE and regenerate it.',faterr)

    read (luc) nsp, nat, alat ! nsp, nat, alat
    read (luc) ! atm
    read (luc) ! ityp
    read (luc) ! tau
    read (luc) at
    at = at * alat

    ! read the dimensions for the arrays
    read (luc) f%qe%nks, f%qe%nbnd, f%qe%nspin, f%qe%gamma_only
    if (f%qe%nspin == 1) then
       fspin = 2d0
    else
       fspin = 1d0
    end if
    read (luc) f%qe%nk(1), f%qe%nk(2), f%qe%nk(3)
    read (luc) n
    read (luc) npwx, ngms
    nkstot = f%qe%nspin * f%qe%nks
    call init_geometry(f,x2c,n,env)

    ! read k-point info
    if (allocated(f%qe%kpt)) deallocate(f%qe%kpt)
    allocate(f%qe%kpt(3,f%qe%nks))
    if (allocated(f%qe%wk)) deallocate(f%qe%wk)
    allocate(f%qe%wk(f%qe%nks))
    if (allocated(f%qe%ek)) deallocate(f%qe%ek)
    allocate(f%qe%ek(f%qe%nbnd,f%qe%nks,f%qe%nspin))
    if (allocated(f%qe%occ)) deallocate(f%qe%occ)
    allocate(f%qe%occ(f%qe%nbnd,f%qe%nks,f%qe%nspin))
    read (luc) f%qe%kpt
    read (luc) f%qe%wk
    read (luc) f%qe%ek
    read (luc) f%qe%occ

    ! read k-point mapping
    if (allocated(f%qe%ngk)) deallocate(f%qe%ngk)
    allocate(f%qe%ngk(f%qe%nks))
    if (allocated(f%qe%igk_k)) deallocate(f%qe%igk_k)
    allocate(f%qe%igk_k(npwx,f%qe%nks))
    if (allocated(f%qe%nl)) deallocate(f%qe%nl)
    allocate(f%qe%nl(ngms))
    read (luc) f%qe%ngk
    read (luc) f%qe%igk_k
    read (luc) f%qe%nl
    if (f%qe%gamma_only) then
       if (allocated(f%qe%nlm)) deallocate(f%qe%nlm)
       allocate(f%qe%nlm(ngms))
       read (luc) f%qe%nlm
    end if

    ! convert k-point coordinates to reciprocal crystallographic and
    ! band energies to Hartree
    do i = 1, f%qe%nks
       f%qe%kpt(:,i) = matmul(f%qe%kpt(:,i),at) / alat
    end do
    f%qe%ek = 0.5d0 * f%qe%ek

    ! calculate the electron density
    if (allocated(f%f)) deallocate(f%f)
    allocate(f%f(f%n(1),f%n(2),f%n(3)))
    allocate(raux(f%n(1),f%n(2),f%n(3)),rseq(f%n(1)*f%n(2)*f%n(3)))
    allocate(evc(maxval(f%qe%ngk(1:f%qe%nks))))
    f%f = 0d0

    do is = 1, f%qe%nspin
       do ik = 1, f%qe%nks
          do ib = 1, f%qe%nbnd
             rseq = 0d0
             read (luc) evc(1:f%qe%ngk(ik))

             ! range checks
             if (is /= ispin .and. ispin /= 0) cycle
             if (f%qe%ek(ib,ik,is) < emin .or. f%qe%ek(ib,ik,is) > emax) cycle
             if (allocated(ikpt)) then
                if (all(ik /= ikpt)) cycle
             end if
             if (allocated(ibnd)) then
                if (all(ib /= ibnd)) cycle
             end if

             rseq(f%qe%nl(f%qe%igk_k(1:f%qe%ngk(ik),ik))) = evc(1:f%qe%ngk(ik))
             if (f%qe%gamma_only) &
                rseq(f%qe%nlm(f%qe%igk_k(1:f%qe%ngk(ik),ik))) = conjg(evc(1:f%qe%ngk(ik)))
             raux = reshape(rseq,shape(raux))
             call cfftnd(3,f%n,+1,raux)
             f%f = f%f + f%qe%occ(ib,ik,is) * real(conjg(raux) * raux,8)
          end do
       end do
    end do
    f%f = f%f * fspin / (det3(at) * sum(f%qe%wk))

    ! assign nbndw
    if (f%qe%nspin == 1) then
       f%qe%nbndw = f%qe%nbnd
    else
       f%qe%nbndw = f%qe%nbnd

       ok1 = .true.
       ok2 = .true.
       do i = f%qe%nbnd, 1, -1
          ok1 = ok1 .and. all(abs(f%qe%occ(i,1:f%qe%nks,1)) < epsocc)
          ok2 = ok2 .and. all(abs(f%qe%occ(i,1:f%qe%nks,2)) < epsocc)
          if (.not.ok1.and..not.ok2) then
             exit
          else
             if (ok1) f%qe%nbndw(1) = i-1
             if (ok2) f%qe%nbndw(2) = i-1
          end if
       end do
    end if

    ! close and clean up
    call fclose(luc)
    f%isinit = .true.
    f%mode = mode_default

    ! wannier transformation, if given
    if (allocated(f%qe%center)) deallocate(f%qe%center)
    if (allocated(f%qe%spread)) deallocate(f%qe%spread)
    if (allocated(f%qe%u)) deallocate(f%qe%u)

  end subroutine read_pwc

  !> Read a grid in elk format -- only first 3d grid in first 3d block
  module subroutine read_elk(f,file,x2c,env,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, ferror, faterr, fclose
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: file !< Input file
    real*8, intent(in) :: x2c(3,3)
    type(environ), intent(in), target :: env
    type(thread_info), intent(in), optional :: ti

    integer :: luc, ios
    integer :: n(3), i, j, k
    real*8 :: dum(3)

    call f%end()

    ! open file for reading
    luc = fopen_read(file,ti=ti)

    ! grid dimension
    read (luc,*,iostat=ios) n
    if (ios /= 0) &
       call ferror('read_elk','Error reading n1, n2, n3',faterr,file)

    call init_geometry(f,x2c,n,env)
    allocate(f%f(n(1),n(2),n(3)),stat=ios)
    if (ios /= 0) &
       call ferror('read_elk','Error allocating grid',faterr,file)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             read (luc,*) dum, f%f(i,j,k)
          end do
       end do
    end do

    call fclose(luc)

    n = f%n
    f%isinit = .true.
    f%isqe = .false.
    f%iswan = .false.
    f%mode = mode_default

  end subroutine read_elk

  !> Reads information from a wannier90 checkpoint file (.chk). Sets
  !> the wannier information in the qe field only (center, spread, u).
  !> Specs from wannier90, 2.0.1 (works, too: 2.1.0)
  module subroutine read_wannier_chk(f,fileup,filedn,ti)
    use tools_math, only: matinv
    use tools_io, only: faterr, ferror, uout, fopen_read, fclose
    use param, only: bohrtoa
    class(grid3), intent(inout) :: f
    character*(*), intent(in) :: fileup
    character*(*), intent(in), optional :: filedn
    type(thread_info), intent(in), optional :: ti

    integer :: lu(2)
    character(len=33) :: header
    integer :: i, j, k, nspin, is
    integer :: nbnd, jbnd, idum, nks, nk(3), ik1, ik2, ik3
    real*8 :: rlatt(3,3), rclatt(3,3), rlatti(3,3)
    character(len=20) :: chkpt1
    logical :: have_disentangled
    real*8, allocatable :: kpt(:,:)

    ! check qe is available
    if (.not.f%isinit) then
       call ferror("read_wannier_chk","cannot read wannier data with non-initialized grid",faterr)
    end if
    if (.not.f%isqe) then
       call ferror("read_wannier_chk","cannot read wannier data without qe data",faterr)
    end if

    ! check nspin consistency
    if (present(filedn)) then
       if (f%qe%nspin /= 2) &
          call ferror("read_wannier_chk","two chk files but nspin = 1",faterr)
    else
       if (f%qe%nspin /= 1) &
          call ferror("read_wannier_chk","one chk file but nspin = 2",faterr)
    end if
    nspin = f%qe%nspin

    ! open files and initialize
    lu(1) = fopen_read(fileup,form="unformatted",ti=ti)
    if (nspin == 2) &
       lu(2) = fopen_read(filedn,form="unformatted",ti=ti)

    ! header and number of bands
    do is = nspin, 1, -1
       read(lu(is)) header
       read(lu(is)) nbnd
       read(lu(is)) jbnd
       if (jbnd > 0) &
          call ferror("read_wannier_chk","number of excluded bands must be 0",faterr)
       if (nbnd /= f%qe%nbnd .and. nspin == 1) &
          call ferror("read_wannier_chk","number of bands different in wannier and qe",faterr)
       read(lu(is)) (idum,i=1,jbnd)

       ! real and reciprocal lattice
       read(lu(is)) ((rlatt(i,j),i=1,3),j=1,3)
       read(lu(is)) ((rclatt(i,j),i=1,3),j=1,3)

       ! number of k-points
       read(lu(is)) nks
       read(lu(is)) nk
       if (nks == 0 .or.any(nk == 0) .or. nks/=(nk(1)*nk(2)*nk(3))) &
          call ferror("read_wannier_chk","error in number of k-points (wannier)",faterr)
       if (nks /= f%qe%nks) &
          call ferror("read_wannier_chk","number of k-points from wannier different than qe",faterr)
    end do

    ! k-points
    allocate(kpt(3,nks))
    do is = nspin, 1, -1
       read(lu(is)) ((kpt(i,j),i=1,3),j=1,nks)
       do i = 1, nks
          ik1 = nint(kpt(1,i) * nk(1))
          ik2 = nint(kpt(2,i) * nk(2))
          ik3 = nint(kpt(3,i) * nk(3))
          if (abs(kpt(1,i) * nk(1) - ik1) > 1d-5 .or.abs(kpt(2,i) * nk(2) - ik2) > 1d-5 .or.&
             abs(kpt(3,i) * nk(3) - ik3) > 1d-5) then
             write (uout,*) kpt(:,i)
             write (uout,*) kpt(1,i)*nk(1),kpt(1,i)*nk(2),kpt(1,i)*nk(3)
             write (uout,*) ik1, ik2, ik3
             call ferror("read_wannier_chk","not a (uniform) monkhorst-pack grid or shifted grid",faterr)
          end if
          if (any(abs(kpt(:,i) - f%qe%kpt(:,i)) > 1d-5)) then
             write (uout,*) i
             write (uout,*) kpt(:,i)
             write (uout,*) f%qe%kpt(:,i)
             call ferror("read_wannier_chk","inconsistent wannier/qe k-point coordinates",faterr)
          end if
       end do
    end do

    ! overwrite the k-point in the qe data (if this was a nscf calculation, the nk(3)
    ! from qe is zero)
    f%qe%nk = nk
    f%qe%nbndw = 0

    do is = nspin, 1, -1
       read(lu(is)) idum ! number of nearest k-point neighbours
       read(lu(is)) jbnd ! number of wannier functions
       f%qe%nbndw(is) = jbnd

       ! checkpoint position and disentanglement
       read(lu(is)) chkpt1
       read(lu(is)) have_disentangled
       if (have_disentangled) &
          call ferror("read_wannier_chk","cannot handle disentangled wannier functions",faterr)
    end do

    ! rest of the file: u and m matrices, wannier centers and spreads
    jbnd = maxval(f%qe%nbndw(1:nspin))
    if (allocated(f%qe%u)) deallocate(f%qe%u)
    if (allocated(f%qe%center)) deallocate(f%qe%center)
    if (allocated(f%qe%spread)) deallocate(f%qe%spread)
    allocate(f%qe%u(jbnd,jbnd,nks,nspin))
    allocate(f%qe%center(3,jbnd,nspin),f%qe%spread(jbnd,nspin))
    f%qe%center = 0d0
    f%qe%spread = 0d0
    do is = nspin, 1, -1
       read(lu(is)) (((f%qe%u(i,j,k,is),i=1,f%qe%nbndw(is)),j=1,f%qe%nbndw(is)),k=1,nks)
       read(lu(is)) ! m matrix
       read(lu(is)) ((f%qe%center(i,j,is),i=1,3),j=1,f%qe%nbndw(is))
       read(lu(is)) (f%qe%spread(i,is),i=1,f%qe%nbndw(is))
       call fclose(lu(is))
    end do

    ! convert centers to crystallographic and spread to bohr
    rlatti = rlatt
    call matinv(rlatti,3)
    do is = 1, nspin
       do i = 1, f%qe%nbndw(is)
          f%qe%center(:,i,is) = matmul(f%qe%center(:,i,is),rlatti)
          do j = 1, 3
             if (f%qe%center(j,i,is) > nk(j)) &
                f%qe%center(j,i,is) = f%qe%center(j,i,is) - nk(j)
             if (f%qe%center(j,i,is) < 0d0) &
                f%qe%center(j,i,is) = f%qe%center(j,i,is) + nk(j)
          end do
          if (f%qe%spread(i,is) < 0d0) then
             f%qe%spread(i,is) = -sqrt(abs(f%qe%spread(i,is))) / bohrtoa
          else
             f%qe%spread(i,is) = sqrt(f%qe%spread(i,is)) / bohrtoa
          end if
       end do
    end do

    ! clean up
    f%iswan = .true.

  end subroutine read_wannier_chk

  !> Interpolate the function value, first and second derivative at
  !> point x0 (crystallographic coords.) using the grid g.  This
  !> routine is thread-safe.
  module subroutine interp(f,xi,y,yp,ypp)
    class(grid3), intent(inout) :: f !< Grid to interpolate
    real*8, intent(in) :: xi(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    real*8 :: x0(3)

    x0 = xi - floor(xi)
    y = 0d0
    yp = 0d0
    ypp = 0d0

    if (f%mode == mode_nearest) then
       call grinterp_nearest(f,x0,y)
       yp = 0d0
       ypp = 0d0
    else if (f%mode == mode_trilinear) then
       call grinterp_trilinear(f,x0,y,yp)
       ypp = 0d0
    else if (f%mode == mode_trispline) then
       call grinterp_trispline(f,x0,y,yp,ypp)
    else if (f%mode == mode_tricubic) then
       call grinterp_tricubic(f,x0,y,yp,ypp)
    else if (f%mode == mode_smr) then
       call grinterp_smr(f,x0,y,yp,ypp)
    end if

    ! convert the gradient and Hessian to Cartesian coordinates
    if (f%mode /= mode_smr) then
       yp = matmul(transpose(f%c2x),yp)
       ypp = matmul(matmul(transpose(f%c2x),ypp),f%c2x)
    end if

  end subroutine interp

  !> Given the grid field in frho calculate the laplacian or Hessian
  !> derivatives of frho with FFT and save them in flap. ix can be
  !> one of 0 (lap), 1 (hxx), 2 (hyy), 3 (hzz).
  module subroutine laplacian_hxx(flap,frho,ix)
    use tools_io, only: ferror, faterr
    use tools_math, only: cross, det3
    use param, only: pi
    class(grid3), intent(inout) :: flap
    type(grid3), intent(in) :: frho
    integer, intent(in) :: ix

    integer :: n(3), i1, i2, i3, ntot, igfft
    complex*16, allocatable :: zaux(:)
    real*8 :: bvec(3,3)
    real*8, allocatable :: vgc(:,:)

    call flap%end()
    if (.not.frho%isinit) &
       call ferror('grid_laplacian','no density grid',faterr)

    n = frho%n
    call copy_geometry(flap,frho)
    flap%isinit = .true.
    flap%mode = mode_default
    ntot = n(1) * n(2) * n(3)

    ! reciprocal lattice vectors
    bvec(:,1) = cross(frho%x2c(:,3),frho%x2c(:,2))
    bvec(:,2) = cross(frho%x2c(:,1),frho%x2c(:,3))
    bvec(:,3) = cross(frho%x2c(:,2),frho%x2c(:,1))
    bvec = 2d0 * pi / abs(det3(frho%x2c)) * bvec

    ! prepare the vectors
    allocate(vgc(3,ntot))
    do i1 = n(1)/2-n(1)+1, n(1)/2
       do i2 = n(2)/2-n(2)+1, n(2)/2
          do i3 = n(3)/2-n(3)+1, n(3)/2
             igfft = modulo(i3,n(3))*n(2)*n(1)+modulo(i2,n(2))*n(1)+modulo(i1,n(1))+1
             vgc(:,igfft)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
          end do
       end do
    end do

    ! allocate the work array and do the FT
    allocate(zaux(ntot))
    zaux = 0d0
    zaux = reshape(frho%f,shape(zaux))
    call cfftnd(3,n,-1,zaux)

    ! calculate the derivative in reciprocal space
    if (ix == 0) then
       do igfft = 1, ntot
          zaux(igfft) = -dot_product(vgc(:,igfft),vgc(:,igfft)) * zaux(igfft)
       end do
    else
       do igfft = 1, ntot
          zaux(igfft) = -vgc(ix,igfft) * vgc(ix,igfft) * zaux(igfft)
       end do
    end if
    deallocate(vgc)

    ! FT back to real space
    call cfftnd(3,n,+1,zaux)
    allocate(flap%f(n(1),n(2),n(3)))
    flap%f = real(reshape(real(zaux,8),shape(flap%f)),8)
    deallocate(zaux)

  end subroutine laplacian_hxx

  !> Calculate the gradient norm of a scalar field using FFT.
  module subroutine gradrho(fgrho,frho)
    use tools_io, only: faterr, ferror
    use tools_math, only: det3, cross
    use param, only: pi
    class(grid3), intent(inout) :: fgrho
    type(grid3), intent(in) :: frho

    integer :: n(3), i, i1, i2, i3, ntot, igfft
    complex*16, allocatable :: zaux(:)
    real*8 :: bvec(3,3)
    real*8, allocatable :: vgc(:,:)

    call fgrho%end()
    if (.not.frho%isinit) &
       call ferror('grid_gradgrho','no density grid',faterr)

    ! allocate slot
    n = frho%n
    call copy_geometry(fgrho,frho)
    fgrho%isinit = .true.
    fgrho%mode = mode_default
    ntot = n(1) * n(2) * n(3)

    ! reciprocal lattice vectors
    bvec(:,1) = cross(frho%x2c(:,3),frho%x2c(:,2))
    bvec(:,2) = cross(frho%x2c(:,1),frho%x2c(:,3))
    bvec(:,3) = cross(frho%x2c(:,2),frho%x2c(:,1))
    bvec = 2d0 * pi / abs(det3(frho%x2c)) * bvec

    ! prepare the vectors
    allocate(vgc(3,ntot))
    do i1 = n(1)/2-n(1)+1, n(1)/2
       do i2 = n(2)/2-n(2)+1, n(2)/2
          do i3 = n(3)/2-n(3)+1, n(3)/2
             igfft = modulo(i3,n(3))*n(2)*n(1)+modulo(i2,n(2))*n(1)+modulo(i1,n(1))+1
             vgc(:,igfft)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
          end do
       end do
    end do

    ! allocate work space and final grid
    allocate(zaux(ntot),fgrho%f(n(1),n(2),n(3)))
    zaux = 0d0
    fgrho%f = 0d0

    do i = 1, 3
       ! FT the source grid to reciprocal space
       zaux = reshape(frho%f,shape(zaux))
       call cfftnd(3,n,-1,zaux)

       ! calculate the derivative in reciprocal space
       do igfft = 1, ntot
          zaux(igfft) = vgc(i,igfft) * cmplx(-aimag(zaux(igfft)),dble(zaux(igfft)),8)
       end do

       ! FT back to real space and accumulate
       call cfftnd(3,n,+1,zaux)
       fgrho%f = fgrho%f + (real(reshape(real(zaux,8),shape(fgrho%f)),8))**2
    end do

    ! wrap up
    fgrho%f = sqrt(fgrho%f)
    deallocate(zaux)

  end subroutine gradrho

  !> Given the grid field in frho, calculate the electrostatic
  !> (Hartree) potential generated by it in fpot using FFT. The k=0
  !> coefficient of the potential Fourier expansion (and thus the
  !> average potential) is set to zero. If isry is .true., then the
  !> output potential is in Ry; otherwise in Hartree.
  module subroutine pot(fpot,frho,isry)
    use tools_io, only: ferror, faterr
    use tools_math, only: cross, det3
    use param, only: pi
    class(grid3), intent(inout) :: fpot
    type(grid3), intent(in) :: frho
    logical, intent(in) :: isry

    integer :: n(3), i1, i2, i3, ntot, igfft
    complex*16, allocatable :: zaux(:)
    real*8 :: bvec(3,3), vgc2
    real*8, allocatable :: vgc(:,:)

    call fpot%end()
    if (.not.frho%isinit) &
       call ferror('grid_pot','no density grid',faterr)

    ! allocate slot
    n = frho%n
    call copy_geometry(fpot,frho)
    fpot%isinit = .true.
    fpot%mode = mode_default
    ntot = n(1) * n(2) * n(3)

    ! reciprocal lattice vectors
    bvec(:,1) = cross(frho%x2c(:,3),frho%x2c(:,2))
    bvec(:,2) = cross(frho%x2c(:,1),frho%x2c(:,3))
    bvec(:,3) = cross(frho%x2c(:,2),frho%x2c(:,1))
    bvec = 2d0 * pi / abs(det3(frho%x2c)) * bvec

    ! prepare the vectors
    allocate(vgc(3,ntot))
    do i1 = n(1)/2-n(1)+1, n(1)/2
       do i2 = n(2)/2-n(2)+1, n(2)/2
          do i3 = n(3)/2-n(3)+1, n(3)/2
             igfft = modulo(i3,n(3))*n(2)*n(1)+modulo(i2,n(2))*n(1)+modulo(i1,n(1))+1
             vgc(:,igfft)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
          end do
       end do
    end do

    ! allocate the work array and do the FT
    allocate(zaux(ntot))
    zaux = 0d0
    zaux = reshape(frho%f,shape(zaux))
    call cfftnd(3,n,-1,zaux)

    ! calculate the derivative in reciprocal space
    do igfft = 1, ntot
       vgc2 = dot_product(vgc(:,igfft),vgc(:,igfft))
       if (vgc2 < 1d-12) then
          zaux(igfft) = 0d0
       else
          zaux(igfft) = -zaux(igfft) / vgc2
       end if
    end do

    ! FT back to real space
    call cfftnd(3,n,+1,zaux)
    allocate(fpot%f(n(1),n(2),n(3)))
    if (isry) then
       fpot%f = -8d0 * pi * real(reshape(real(zaux,8),shape(fpot%f)),8)
    else
       fpot%f = -4d0 * pi * real(reshape(real(zaux,8),shape(fpot%f)),8)
    end if

    deallocate(vgc,zaux)

  end subroutine pot

  !> Resample an input grid field (frho) to a new size given by n2
  !> using FFT.
  module subroutine resample(frs,frho,n2)
    use tools_io, only: ferror, faterr
    class(grid3), intent(inout) :: frs
    type(grid3), intent(in) :: frho
    integer, intent(in) :: n2(3)

    complex*16, allocatable :: zaux(:), zaux2(:)
    integer :: ig, i1, i2, i3, iaux
    integer :: n(3), ntot, ntot2, igfft, igfft2

    call frs%end()
    if (.not.frho%isinit) &
       call ferror('resample','no density grid',faterr)

    ! allocate slot
    call init_geometry(frs,frho%x2c,n2,frho%env)
    frs%isinit = .true.
    frs%mode = mode_default

    ! FT the source density to reciprocal space
    n = frho%n
    ntot = n(1) * n(2) * n(3)
    allocate(zaux(ntot))
    zaux = 0d0
    zaux = reshape(frho%f,shape(zaux))
    call cfftnd(3,frho%n,-1,zaux)

    ! allocate space for the resampled grid
    ntot2 = n2(1) * n2(2) * n2(3)
    allocate(zaux2(ntot2))
    zaux2 = 0d0

    ! place the elements of the original density into the new array
    do ig = 1, ntot
       i3 = mod(ig-1,n(3)) + n(3)/2 - n(3) + 1
       iaux = (ig-1 - (i3-1 + n(3) - n(3)/2)) / n(3)
       i2 = mod(iaux,n(2)) + n(2)/2 - n(2) + 1
       iaux = (iaux - (i2-1 + n(2) - n(2)/2)) / n(2)
       i1 = iaux + n(1)/2 - n(1) + 1

       igfft = modulo(i3,n(3))*n(2)*n(1)+modulo(i2,n(2))*n(1)+modulo(i1,n(1))+1
       igfft2 = modulo(i3,n2(3))*n2(2)*n2(1)+modulo(i2,n2(2))*n2(1)+modulo(i1,n2(1))+1
       zaux2(igfft2) = zaux(igfft)
    end do
    deallocate(zaux)

    ! FT back to real space
    call cfftnd(3,n2,+1,zaux2)
    allocate(frs%f(n2(1),n2(2),n2(3)))
    frs%f = real(reshape(real(zaux2,8),shape(frs%f)),8)
    deallocate(zaux2)

  end subroutine resample

  !> Read the pwc file from QE and write temporary, perhaps rotated
  !> versions of the same file. This operation pre-arranges the data
  !> in a manner that makes it easy for get_qe_wnr to work in
  !> parallel. Returns the LUs of the spin-up and spin-down rotated
  !> pwc files, which remain open. luevc_ibnd keeps track of the next
  !> band in the LU (returns 1 for both spins). If useu is false, do
  !> not rotate the coefficients (but still write the scratch files).
  module subroutine rotate_qe_evc(f,luevc,luevc_ibnd,useu,ti)
    use tools_io, only: fopen_scratch, string, fopen_read, fclose
    class(grid3), intent(inout) :: f
    integer, intent(out) :: luevc(2)
    integer, intent(out) :: luevc_ibnd(2)
    logical, intent(in) :: useu
    type(thread_info), intent(in), optional :: ti

    integer :: i, ik, ibnd, jbnd, luc, ireg
    complex*16, allocatable :: evc(:), evcnew(:)

    luevc = -1
    luevc_ibnd = 0

    luc = fopen_read(f%qe%fpwc,form="unformatted",ti=ti)

    ! skip the header of the pwc file
    if (f%qe%gamma_only) then
       ireg = 18
    else
       ireg = 17
    end if

    allocate(evc(maxval(f%qe%ngk(1:f%qe%nks))),evcnew(maxval(f%qe%ngk(1:f%qe%nks))))

    do i = 1, f%qe%nspin
       luevc(i) = fopen_scratch(form="unformatted",ti=ti)
    end do
    do ibnd = 1, f%qe%nbnd
       rewind(luc)
       do i = 1, ireg
          read (luc)
       end do
       do i = 1, f%qe%nspin
          do ik = 1, f%qe%nks
             evcnew = 0d0
             do jbnd = 1, f%qe%nbnd
                read (luc) evc(1:f%qe%ngk(ik))
                if (ibnd > f%qe%nbndw(i) .or. jbnd > f%qe%nbndw(i)) cycle
                if (useu) then
                   evcnew(1:f%qe%ngk(ik)) = evcnew(1:f%qe%ngk(ik)) + f%qe%u(jbnd,ibnd,ik,i) * evc(1:f%qe%ngk(ik))
                elseif (ibnd == jbnd) then
                   evcnew(1:f%qe%ngk(ik)) = evc(1:f%qe%ngk(ik))
                end if
             end do
             if (ibnd > f%qe%nbndw(i)) cycle
             write (luevc(i)) evcnew(1:f%qe%ngk(ik))
          end do
       end do
    end do
    do i = 1, f%qe%nspin
       rewind(luevc(i))
       luevc_ibnd(i) = 1
    end do

    deallocate(evc,evcnew)
    call fclose(luc)

  end subroutine rotate_qe_evc

  !> Build the Wannier functions for band ibnd and spin ispin from QEs
  !> Bloch coefficients, parallel version. Returns the Wannier
  !> functions in supercell grid fout. This version uses the logical
  !> units for the scratch files containing the rotated Bloch
  !> coefficients (luevc, one for each spin) and a pointer that tracks
  !> the current band index in those files (luevc_ibnd). The latter is
  !> updated by calls to this routine. omega is the cell volume (used
  !> for normalization).
  module subroutine get_qe_wnr(f,omega,ibnd,ispin,luevc,luevc_ibnd,fout)
    use tools_io, only: uout, ferror, faterr
    use param, only: tpi, img
    class(grid3), intent(in) :: f
    real*8, intent(in) :: omega
    integer, intent(in) :: ibnd
    integer, intent(in) :: ispin
    integer, intent(in) :: luevc(2)
    integer, intent(inout) :: luevc_ibnd(2)
    complex*16, intent(out) :: fout(:,:,:,:)

    integer :: i, j, k
    integer :: ik, ik1, ik2, ik3, ikk, ilat, ikg, ik0
    complex*16, allocatable :: evc(:), rseq(:), raux(:,:,:), raux2(:,:,:)
    real*8 :: xkpt(3)

    ! some checks
    if (f%n(1) /= size(fout,1).or.f%n(2) /= size(fout,2).or.f%n(3) /= size(fout,3)) &
       call ferror("get_qe_wnr","inconsistent grid size",faterr)
    if (f%qe%nk(1)*f%qe%nk(2)*f%qe%nk(3) /= size(fout,4)) &
       call ferror("get_qe_wnr","inconsistent number of k-points",faterr)
    if (f%qe%nk(1)*f%qe%nk(2)*f%qe%nk(3) /= f%qe%nks) &
       call ferror("get_qe_wnr","inconsistent number of k-points",faterr)

    allocate(evc(maxval(f%qe%ngk(1:f%qe%nks))))

    if (luevc_ibnd(ispin) > ibnd) then
       rewind(luevc(ispin))
       luevc_ibnd(ispin) = 1
    end if
    if (luevc_ibnd(ispin) /= ibnd) then
       do i = 1, ibnd-luevc_ibnd(ispin)
          do ik = 1, f%qe%nks
             read (luevc(ispin)) evc(1:f%qe%ngk(ik))
          end do
       end do
       luevc_ibnd(ispin) = ibnd
    end if

    ! allocate auxiliary arrays
    allocate(rseq(f%n(1)*f%n(2)*f%n(3)))
    allocate(raux(f%n(1),f%n(2),f%n(3)))
    allocate(raux2(f%n(1),f%n(2),f%n(3)))
    rseq = 0d0
    raux2 = 0d0
    raux = 0d0

    ! run over k-points
    fout = 0d0
    ikg = 0
    !$omp parallel do private(ik1,ik2,ik3,ilat,ik0,xkpt) firstprivate(evc,rseq,raux,raux2)
    do ik = 1, f%qe%nks
       rseq = 0d0
       !$omp critical (readio)
       ikg = ikg + 1
       ik0 = ikg
       read (luevc(ispin)) evc(1:f%qe%ngk(ikg))
       !$omp end critical (readio)
       rseq(f%qe%nl(f%qe%igk_k(1:f%qe%ngk(ik0),ik0))) = evc(1:f%qe%ngk(ik0))
       raux = reshape(rseq,shape(raux))
       call cfftnd(3,f%n,+1,raux)

       ! phase factor
       do k = 1, f%n(3)
          do j = 1, f%n(2)
             do i = 1, f%n(1)
                raux(i,j,k) = raux(i,j,k) * exp(tpi*img*(f%qe%kpt(1,ik0)*real(i-1,8)/real(f%n(1),8)+&
                   f%qe%kpt(2,ik0)*real(j-1,8)/real(f%n(2),8)+f%qe%kpt(3,ik0)*real(k-1,8)/real(f%n(3),8)))
             end do
          end do
       end do

       ! add the contribution from this k-point
       do ikk = 1, f%qe%nks
          xkpt = f%qe%kpt(:,ikk) - floor(f%qe%kpt(:,ikk))
          ik1 = mod(nint(xkpt(1) * f%qe%nk(1)),f%qe%nk(1))
          ik2 = mod(nint(xkpt(2) * f%qe%nk(2)),f%qe%nk(2))
          ik3 = mod(nint(xkpt(3) * f%qe%nk(3)),f%qe%nk(3))
          ilat = 1 + ik3 + f%qe%nk(3) * (ik2 + f%qe%nk(2) * ik1)
          raux2 = raux * exp(-tpi*img*(f%qe%kpt(1,ik0)*ik1+f%qe%kpt(2,ik0)*ik2+f%qe%kpt(3,ik0)*ik3))
          if (ilat < 1 .or. ilat > f%qe%nks) then
             !$omp critical (ioerror)
             write (uout,*) "kpoint number ", ikk
             write (uout,*) "kpoint coords ", f%qe%kpt(:,ikk)
             write (uout,*) "kpoint coords (main cell) ", xkpt
             write (uout,*) "ik ", ik1, ik2, ik3
             write (uout,*) "ilat ", ilat
             call ferror("get_qe_wnr","could not classify a k-point, non-uniform grid?",faterr)
             !$omp end critical (ioerror)
          end if
          !$omp critical (sum)
          fout(:,:,:,ilat) = fout(:,:,:,ilat) + raux2
          !$omp end critical (sum)
       end do
    end do
    !$omp end parallel do

    ! clean up
    deallocate(rseq,raux,raux2)
    luevc_ibnd(ispin) = luevc_ibnd(ispin) + 1
    if (allocated(evc)) deallocate(evc)

    ! normalize
    fout = fout / real(f%qe%nks,8) / sqrt(omega)

 end subroutine get_qe_wnr

  !> Build the Wannier functions for band ibnd and spin ispin from QEs
  !> Bloch coefficients, standalone version. Returns the Wannier
  !> function for lattice vector inr in cell grid fout. This version
  !> works on its own (c.f. get_qe_wnr) and is therefore slower. omega
  !> is the cell volume (used for normalization).
  module subroutine get_qe_wnr_standalone(f,omega,ibnd,ispin,inr,rotate,fout,ti)
    use tools_io, only: fopen_read, fopen_scratch, fclose, ferror, faterr
    use param, only: tpi, img
    class(grid3), intent(in) :: f
    real*8, intent(in) :: omega
    integer, intent(in) :: ibnd
    integer, intent(in) :: ispin
    integer, intent(in) :: inr(3)
    logical, intent(in) :: rotate
    complex*16, intent(out) :: fout(:,:,:)
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, k, is, ik, jbnd, luc, ireg
    complex*16, allocatable :: evcaux(:), evc(:,:), rseq(:), raux(:,:,:)

    ! some checks
    if (f%n(1) /= size(fout,1).or.f%n(2) /= size(fout,2).or.f%n(3) /= size(fout,3)) &
       call ferror("get_qe_wnr_standalone","inconsistent grid size",faterr)
    if (any(abs(f%qe%wk - f%qe%wk(1)) > 1d-5)) &
       call ferror("get_qe_wnr_standalone","wannier transformation only possible with uniform grids (no symmetry)",faterr)

    ! open the pwc file
    luc = fopen_read(f%qe%fpwc,form="unformatted",ti=ti)

    ! skip the header of the pwc file
    if (f%qe%gamma_only) then
       ireg = 18
    else
       ireg = 17
    end if

    ! read the psik coefficients from the file and perform the rotation; calculates evc
    allocate(evcaux(maxval(f%qe%ngk(1:f%qe%nks))))
    allocate(evc(maxval(f%qe%ngk(1:f%qe%nks)),f%qe%nks))
    evc = 0d0
    do i = 1, ireg
       read (luc)
    end do
    do is = 1, min(f%qe%nspin,ispin)
       do ik = 1, f%qe%nks
          do jbnd = 1, f%qe%nbnd
             read (luc) evcaux(1:f%qe%ngk(ik))
             if (f%iswan .and. (ibnd > f%qe%nbndw(is) .or. jbnd > f%qe%nbndw(is) .or. is < ispin)) cycle ! nbndw(:) only from wannier
             if (ibnd > f%qe%nbnd) cycle

             if (rotate) then
                evc(1:f%qe%ngk(ik),ik) = evc(1:f%qe%ngk(ik),ik) + f%qe%u(jbnd,ibnd,ik,is) * evcaux(1:f%qe%ngk(ik))
             elseif (ibnd == jbnd .and. is == ispin) then
                evc(1:f%qe%ngk(ik),ik) = evcaux(1:f%qe%ngk(ik))
             end if
          end do
       end do
    end do
    deallocate(evcaux)
    call fclose(luc)

    ! allocate auxiliary arrays
    allocate(rseq(f%n(1)*f%n(2)*f%n(3)))
    allocate(raux(f%n(1),f%n(2),f%n(3)))
    rseq = 0d0
    raux = 0d0
    fout = 0d0

    !$omp parallel do firstprivate(rseq,raux)
    do ik = 1, f%qe%nks
       ! calculate unk(r) as FT of the coefficients
       rseq = 0d0
       rseq(f%qe%nl(f%qe%igk_k(1:f%qe%ngk(ik),ik))) = evc(1:f%qe%ngk(ik),ik)
       raux = reshape(rseq,shape(raux))
       call cfftnd(3,f%n,+1,raux)

       ! multiply by the phase factor (e(ik*r))
       do k = 1, f%n(3)
          do j = 1, f%n(2)
             do i = 1, f%n(1)
                raux(i,j,k) = raux(i,j,k) * exp(tpi*img*(f%qe%kpt(1,ik)*real(i-1,8)/real(f%n(1),8)+&
                   f%qe%kpt(2,ik)*real(j-1,8)/real(f%n(2),8)+f%qe%kpt(3,ik)*real(k-1,8)/real(f%n(3),8)))
             end do
          end do
       end do

       ! the phase factor for this lattice vector (e(-ik*R))
       raux = raux * exp(-tpi*img*(f%qe%kpt(1,ik)*inr(1)+f%qe%kpt(2,ik)*inr(2)+f%qe%kpt(3,ik)*inr(3)))

       !$omp critical (sum)
       fout = fout + raux
       !$omp end critical (sum)
    end do
    !$omp end parallel do

    ! cleanup
    deallocate(rseq,raux)

    ! normalize
    fout = fout / real(f%qe%nks,8) / sqrt(omega)

  end subroutine get_qe_wnr_standalone

  !> Build the unk(r) functions on a real grid for band ibnd, k-point
  !> ik, and spin ispin from QEs Bloch coefficients, standalone
  !> version. Returns the unk(r) in cell grid fout. omega is the cell
  !> volume (used for normalization).
  module subroutine get_qe_psink_standalone(f,omega,ibnd,ik,ispin,usephase,inr,fout,ti)
    use tools_io, only: fopen_read, fopen_scratch, fclose, ferror, faterr
    use param, only: tpi, img
    class(grid3), intent(in) :: f
    real*8, intent(in) :: omega
    integer, intent(in) :: ibnd
    integer, intent(in) :: ik
    integer, intent(in) :: ispin
    logical :: usephase
    integer, intent(in) :: inr(3)
    complex*16, intent(out) :: fout(:,:,:)
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, k, is, ik_, jbnd, luc, ireg
    complex*16, allocatable :: evc(:), evcaux(:), rseq(:)

    ! some checks
    if (f%n(1) /= size(fout,1).or.f%n(2) /= size(fout,2).or.f%n(3) /= size(fout,3)) &
       call ferror("get_qe_psink_standalone","inconsistent grid size",faterr)

    ! open the pwc file
    luc = fopen_read(f%qe%fpwc,form="unformatted",ti=ti)

    ! skip the header of the pwc file
    if (f%qe%gamma_only) then
       ireg = 18
    else
       ireg = 17
    end if

    ! read the psik coefficients from the file; calculates evc
    allocate(evc(f%qe%ngk(ik)),evcaux(maxval(f%qe%ngk(1:f%qe%nks))))
    evc = 0d0
    do i = 1, ireg
       read (luc)
    end do
    main: do is = 1, f%qe%nspin
       do ik_ = 1, f%qe%nks
          do jbnd = 1, f%qe%nbnd
             read (luc) evcaux(1:f%qe%ngk(ik_))
             if (is == ispin .and. ik == ik_ .and. ibnd == jbnd) then
                evc(1:f%qe%ngk(ik_)) = evcaux(1:f%qe%ngk(ik_))
                exit main
             end if
          end do
       end do
    end do main
    deallocate(evcaux)
    call fclose(luc)

    ! calculate the unk(r) as FT of the Bloch coefficients
    allocate(rseq(f%n(1)*f%n(2)*f%n(3)))
    rseq = 0d0
    rseq(f%qe%nl(f%qe%igk_k(1:f%qe%ngk(ik),ik))) = evc(1:f%qe%ngk(ik))
    fout = reshape(rseq,shape(fout))
    call cfftnd(3,f%n,+1,fout)
    deallocate(rseq)

    ! the phase factor for this k-point (e(ik*r))
    if (usephase) then
       do k = 1, f%n(3)
          do j = 1, f%n(2)
             do i = 1, f%n(1)
                fout(i,j,k) = fout(i,j,k) * exp(tpi*img*(f%qe%kpt(1,ik)*(real(i-1,8)/real(f%n(1),8))+&
                   f%qe%kpt(2,ik)*(real(j-1,8)/real(f%n(2),8))+&
                   f%qe%kpt(3,ik)*(real(k-1,8)/real(f%n(3),8))))
             end do
          end do
       end do
    end if

    ! the phase factor for this lattice vector (e(-ik*R))
    fout = fout * exp(-tpi*img*(f%qe%kpt(1,ik)*inr(1)+f%qe%kpt(2,ik)*inr(2)+f%qe%kpt(3,ik)*inr(3)))

    ! normalize
    fout = fout / sqrt(omega)

  end subroutine get_qe_psink_standalone

  !xx! private procedures

  !> Interpolate by giving the value of the grid at the nearest point
  !> (or the would-be nearest, it is nearest only in orthogonal
  !> grids). f is the grid field to interpolate. xi is the point in
  !> crystallographic coordinates. y is the interpolated value at xi.
  !> This routine is thread-safe.
  subroutine grinterp_nearest(f,x0,y)
    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value

    real*8 :: x(3)
    integer :: idx(3)

    x = modulo(x0,1d0)
    idx = grid_near(f,x)
    y = f%f(idx(1),idx(2),idx(3))

  end subroutine grinterp_nearest

  !> Interpolate using a trilinear interpolant. f is the grid field to
  !> interpolate. xi is the point in crystallographic coordinates. y
  !> and yp are the value and gradient at xi.  This routine is
  !> thread-safe.
  subroutine grinterp_trilinear(f,x0,y,yp)
    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative

    integer :: idx(3), iidx(3), i, j, k
    real*8 :: ff(0:2,0:2,0:2), r(3), s(3), x(3)

    ! compute value at the cube vertices
    ff = 0d0
    x = modulo(x0,1d0)
    idx = grid_floor(f,x)
    do i = 0, 1
       do j = 0, 1
          do k = 0, 1
             iidx = modulo(idx+(/i,j,k/)-1,f%n)+1
             ff(i,j,k) = f%f(iidx(1),iidx(2),iidx(3))
          end do
       end do
    end do

    ! x and 1-x
    r = f%n * x - idx + 1
    s = 1d0 - r

    ! trilinear interpolation
    do i = 0, 1
       do j = 0, 1
          ff(i,j,2) = ff(i,j,0) * s(3) + ff(i,j,1) * r(3)
          ff(i,2,j) = ff(i,0,j) * s(2) + ff(i,1,j) * r(2)
          ff(2,i,j) = ff(0,i,j) * s(1) + ff(1,i,j) * r(1)
       end do
    end do
    do i = 0, 1
       ff(i,2,2) = ff(i,0,2) * s(2) + ff(i,1,2) * r(2)
       ff(2,i,2) = ff(2,i,0) * s(3) + ff(2,i,1) * r(3)
       ff(2,2,i) = ff(0,2,i) * s(1) + ff(1,2,i) * r(1)
    end do
    ff(2,2,2) = ff(0,2,2) * s(1) + ff(1,2,2) * r(1)
    y = ff(2,2,2)
    yp(1) = ff(1,2,2) - ff(0,2,2)
    yp(2) = ff(2,1,2) - ff(2,0,2)
    yp(3) = ff(2,2,1) - ff(2,2,0)

    ! from cell to crystallographic coordinates
    yp = yp * f%n

  end subroutine grinterp_trilinear

  !> Using the global cubic spline information (c2 array) with
  !> periodic boundary conditions, calculate the density and c2
  !> coeffs.  of the 6 points on the cube faces, and then average the
  !> spline prediction. This is a modified version of the abinit
  !> density interpolation subroutine. f is the grid field to
  !> interpolate. xi is the point in crystallographic coordinates. y,
  !> yp, and ypp are the value, gradient, and Hessian at xi.  This
  !> routine is thread-safe.
  subroutine grinterp_trispline(f,x0,y,yp,ypp)
    class(grid3), intent(inout), target :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    integer :: i, ii, jj, kk, ll, nn, indx(3), oii, onn, omm
    integer :: inii(4,3)
    real*8 :: bbb, ddu(2), hrh(2), hh(4,2), grd(4), lder(4)
    real*8 :: cof(2,3), ddstar(6), rhstar(6), sqder(6,4), sqvlr(6,4)
    real*8 :: pomsq(2,3), pom2sq(2,3)
    real*8,pointer :: ptddx(:,:,:),ptddy(:,:,:),ptddz(:,:,:),ptrho(:,:,:)
    real*8 :: xx(3), dix(3)

    real*8, parameter :: eps = 1d-12

    !$omp critical (checkalloc)
    if (.not.allocated(f%c2)) then
       call init_trispline(f)
    end if
    !$omp end critical (checkalloc)

    nullify(ptddx,ptddy,ptddz,ptrho)

    xx = modulo(x0,1d0)

    do i=1,3
       dix(i)=1d0/f%n(i)
    end do

    ! determine the index in the grid
    do ii=1,3
       indx(ii)=int(xx(ii)*f%n(ii))
       bbb=(xx(ii)-indx(ii)*dix(ii))*f%n(ii)
       if (indx(ii)==f%n(ii)) then
          indx(ii)=1
          xx(ii)=0.d0
       else
          indx(ii)=indx(ii)+1
       end if
       !  Explicit handling to avoid numeric problems
       if (bbb > 1.d0+eps ) then
          cof(1,ii)=0.d0
          cof(2,ii)=1.d0
       elseif (bbb < -eps ) then
          cof(1,ii)=1.d0
          cof(2,ii)=0.d0
       else
          cof(1,ii)=1.d0-bbb
          cof(2,ii)=bbb
       end if
    end do

    ! 3d interpolation of the valence density

    ! determination of the values of density and of its second derivative
    ! at the "star" = constructed at vv with primitive directions
    ! To interpolation the values at the faces of the grid cell are needed

    rhstar(:)=0.d0
    sqder(:,:)=0.d0
    sqvlr(:,:)=0.d0
    ddstar(:)=0.d0
    pomsq(:,:)=0.d0
    pom2sq(:,:)=0.d0

    oii=1; onn=1; omm=1
    if (indx(1)==f%n(1)) oii=1-f%n(1)
    if (indx(2)==f%n(2)) onn=1-f%n(2)
    if (indx(3)==f%n(3)) omm=1-f%n(3)

    ! the values in the corners of the grid cell

    ptddx=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,1)
    ptddy=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,2)
    ptddz=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,3)
    ptrho=>f%f(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)

    ! the coefficients for spline interpolation of density and its derivation
    do ii=1,3
       do jj=1,2
          pomsq(jj,ii)=(cof(jj,ii)*cof(jj,ii)*cof(jj,ii)-cof(jj,ii))/6.d0*dix(ii)*dix(ii)
          pom2sq(jj,ii)=(3.d0*cof(jj,ii)*cof(jj,ii)-1.d0)/6.d0*dix(ii)
          if (jj==1) pom2sq(jj,ii)=-pom2sq(jj,ii)
       end do
    end do


    do ii=1,2
       do jj=1,2
          do kk=1,2
             ddstar(ii)=ddstar(ii)+cof(jj,2)*cof(kk,3)*ptddx(ii,jj,kk)
             ddstar(ii+2)=ddstar(ii+2)+cof(jj,3)*cof(kk,1)*ptddy(kk,ii,jj)
             ddstar(ii+4)=ddstar(ii+4)+cof(jj,1)*cof(kk,2)*ptddz(jj,kk,ii)
             sqder(ii,jj)=sqder(ii,jj)+cof(kk,2)*ptddz(ii,kk,jj)
             sqder(ii,jj+2)=sqder(ii,jj+2)+cof(kk,3)*ptddy(ii,jj,kk)
             sqder(ii+2,jj)=sqder(ii+2,jj)+cof(kk,3)*ptddx(jj,ii,kk)
             sqder(ii+2,jj+2)=sqder(ii+2,jj+2)+cof(kk,1)*ptddz(kk,ii,jj)
             sqder(ii+4,jj)=sqder(ii+4,jj)+cof(kk,1)*ptddy(kk,jj,ii)
             sqder(ii+4,jj+2)=sqder(ii+4,jj+2)+cof(kk,2)*ptddx(jj,kk,ii)
             sqvlr(ii,jj)=sqvlr(ii,jj)+cof(kk,2)*ptrho(ii,kk,jj)+pomsq(kk,2)*ptddy(ii,kk,jj)
             sqvlr(ii,jj+2)=sqvlr(ii,jj+2)+cof(kk,3)*ptrho(ii,jj,kk)+pomsq(kk,3)*ptddz(ii,jj,kk)
             sqvlr(ii+2,jj+2)=sqvlr(ii+2,jj+2)+cof(kk,1)*ptrho(kk,ii,jj)+pomsq(kk,1)*ptddx(kk,ii,jj)
          end do
       end do
    end do

    do ii=1,2
       do jj=1,2
          sqvlr(ii+2,jj)=sqvlr(jj,ii+2)
          sqvlr(ii+4,jj)=sqvlr(jj+2,ii+2)
          sqvlr(ii+4,jj+2)=sqvlr(jj,ii)
       end do
    end do

    do ii=1,2
       do jj=1,2
          rhstar(ii)=rhstar(ii)+cof(jj,3)*sqvlr(ii,jj)+pomsq(jj,3)*sqder(ii,jj)+&
             &    cof(jj,2)*sqvlr(ii,jj+2)+pomsq(jj,2)*sqder(ii,jj+2)
          rhstar(ii+2)=rhstar(ii+2)+cof(jj,1)*sqvlr(ii+2,jj)+pomsq(jj,1)*sqder(ii+2,jj)+&
             &    cof(jj,3)*sqvlr(ii+2,jj+2)+pomsq(jj,3)*sqder(ii+2,jj+2)
          rhstar(ii+4)=rhstar(ii+4)+cof(jj,2)*sqvlr(ii+4,jj)+pomsq(jj,2)*sqder(ii+4,jj)+&
             &    cof(jj,1)*sqvlr(ii+4,jj+2)+pomsq(jj,1)*sqder(ii+4,jj+2)
       end do
    end do
    rhstar(:)=rhstar(:)/2.d0

    y = 0d0
    yp = 0d0
    ypp = 0d0
    kk=1; nn=1
    do ii=1,5,2
       do jj=1,2
          nn=-nn
          y=y+cof(jj,kk)*rhstar(ii+jj-1)+pomsq(jj,kk)*ddstar(ii+jj-1)
          yp(kk)=yp(kk)+pom2sq(jj,kk)*ddstar(ii+jj-1)
          ypp(kk,kk)=ypp(kk,kk)+cof(jj,kk)*ddstar(ii+jj-1)
          yp(kk)=yp(kk)+nn*rhstar(ii+jj-1)/dix(kk)
       end do
       kk=kk+1
    end do
    y=y/3.d0

    ! Off-diagonal elements of the hessian

    ! for the speed reasons the polynomial interpolation
    ! for second derivation fields is used in this case
    ! but the last step is always done by spline interpolation.
    do ii=1,3
       do jj=-1,2
          inii(jj+2,ii)=indx(ii)+jj
          if (inii(jj+2,ii) < 1) inii(jj+2,ii)=inii(jj+2,ii)+f%n(ii)
          if (inii(jj+2,ii) > f%n(ii)) inii(jj+2,ii)=inii(jj+2,ii)-f%n(ii)
       end do
    end do

    ! Not very nice
    do ii=1,3
       select case (ii)
       case (1)
          do jj=1,4
             ddu(1)=cof(1,2)*f%c2(inii(jj,1),inii(2,2),inii(2,3),3)+cof(2,2)*f%c2(inii(jj,1),inii(3,2),inii(2,3),3)
             ddu(2)=cof(1,2)*f%c2(inii(jj,1),inii(2,2),inii(3,3),3)+cof(2,2)*f%c2(inii(jj,1),inii(3,2),inii(3,3),3)
             hrh(1)=cof(1,2)*f%f(inii(jj,1),inii(2,2),inii(2,3))+cof(2,2)*f%f(inii(jj,1),inii(3,2),inii(2,3))+&
                &      pomsq(1,2)*f%c2(inii(jj,1),inii(2,2),inii(2,3),2)+pomsq(2,2)*f%c2(inii(jj,1),inii(3,2),inii(2,3),2)
             hrh(2)=cof(1,2)*f%f(inii(jj,1),inii(2,2),inii(3,3))+cof(2,2)*f%f(inii(jj,1),inii(3,2),inii(3,3))+&
                &      pomsq(1,2)*f%c2(inii(jj,1),inii(2,2),inii(3,3),2)+pomsq(2,2)*f%c2(inii(jj,1),inii(3,2),inii(3,3),2)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)

             ddu(1)=cof(1,3)*f%c2(inii(jj,1),inii(2,2),inii(2,3),2)+cof(2,3)*f%c2(inii(jj,1),inii(2,2),inii(3,3),2)
             ddu(2)=cof(1,3)*f%c2(inii(jj,1),inii(3,2),inii(2,3),2)+cof(2,3)*f%c2(inii(jj,1),inii(3,2),inii(3,3),2)
             hrh(1)=cof(1,3)*f%f(inii(jj,1),inii(2,2),inii(2,3))+cof(2,3)*f%f(inii(jj,1),inii(2,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(jj,1),inii(2,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(jj,1),inii(2,2),inii(3,3),3)
             hrh(2)=cof(1,3)*f%f(inii(jj,1),inii(3,2),inii(2,3))+cof(2,3)*f%f(inii(jj,1),inii(3,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(jj,1),inii(3,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(jj,1),inii(3,2),inii(3,3),3)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)
          end do
       case (2)
          do jj=1,4
             ddu(1)=cof(1,3)*f%c2(inii(2,1),inii(jj,2),inii(2,3),1)+cof(2,3)*f%c2(inii(2,1),inii(jj,2),inii(3,3),1)
             ddu(2)=cof(1,3)*f%c2(inii(3,1),inii(jj,2),inii(2,3),1)+cof(2,3)*f%c2(inii(3,1),inii(jj,2),inii(3,3),1)
             hrh(1)=cof(1,3)*f%f(inii(2,1),inii(jj,2),inii(2,3))+cof(2,3)*f%f(inii(2,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(2,1),inii(jj,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(2,1),inii(jj,2),inii(3,3),3)
             hrh(2)=cof(1,3)*f%f(inii(3,1),inii(jj,2),inii(2,3))+cof(2,3)*f%f(inii(3,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(3,1),inii(jj,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(3,1),inii(jj,2),inii(3,3),3)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)

             ddu(1)=cof(1,1)*f%c2(inii(2,1),inii(jj,2),inii(2,3),3)+cof(2,1)*f%c2(inii(3,1),inii(jj,2),inii(2,3),3)
             ddu(2)=cof(1,1)*f%c2(inii(2,1),inii(jj,2),inii(3,3),3)+cof(2,1)*f%c2(inii(3,1),inii(jj,2),inii(3,3),3)
             hrh(1)=cof(1,1)*f%f(inii(2,1),inii(jj,2),inii(2,3))+cof(2,1)*f%f(inii(3,1),inii(jj,2),inii(2,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(jj,2),inii(2,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(jj,2),inii(2,3),1)
             hrh(2)=cof(1,1)*f%f(inii(2,1),inii(jj,2),inii(3,3))+cof(2,1)*f%f(inii(3,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(jj,2),inii(3,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(jj,2),inii(3,3),1)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)
          end do
       case (3)
          do jj=1,4
             ddu(1)=cof(1,1)*f%c2(inii(2,1),inii(2,2),inii(jj,3),2)+cof(2,1)*f%c2(inii(3,1),inii(2,2),inii(jj,3),2)
             ddu(2)=cof(1,1)*f%c2(inii(2,1),inii(3,2),inii(jj,3),2)+cof(2,1)*f%c2(inii(3,1),inii(3,2),inii(jj,3),2)
             hrh(1)=cof(1,1)*f%f(inii(2,1),inii(2,2),inii(jj,3))+cof(2,1)*f%f(inii(3,1),inii(2,2),inii(jj,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(2,2),inii(jj,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(2,2),inii(jj,3),1)
             hrh(2)=cof(1,1)*f%f(inii(2,1),inii(3,2),inii(jj,3))+cof(2,1)*f%f(inii(3,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(3,2),inii(jj,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(3,2),inii(jj,3),1)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)

             ddu(1)=cof(1,2)*f%c2(inii(2,1),inii(2,2),inii(jj,3),1)+cof(2,2)*f%c2(inii(2,1),inii(3,2),inii(jj,3),1)
             ddu(2)=cof(1,2)*f%c2(inii(3,1),inii(2,2),inii(jj,3),1)+cof(2,2)*f%c2(inii(3,1),inii(3,2),inii(jj,3),1)
             hrh(1)=cof(1,2)*f%f(inii(2,1),inii(2,2),inii(jj,3))+cof(2,2)*f%f(inii(2,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,2)*f%c2(inii(2,1),inii(2,2),inii(jj,3),2)+pomsq(2,2)*f%c2(inii(2,1),inii(3,2),inii(jj,3),2)
             hrh(2)=cof(1,2)*f%f(inii(3,1),inii(2,2),inii(jj,3))+cof(2,2)*f%f(inii(3,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,2)*f%c2(inii(3,1),inii(2,2),inii(jj,3),2)+pomsq(2,2)*f%c2(inii(3,1),inii(3,2),inii(jj,3),2)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)
          end do
       end select
       do jj=-2,1
          grd(jj+3)=(indx(ii)+jj)*dix(ii)
       end do

       !  write(6,'("hh: ",/,4F16.8,/,4F16.8)') ((hh(kk,jj),kk=1,4),jj=1,2)
       !  write(6,'("grad: ",3F16.8)') (grad(kk),kk=1,3)
       !  write(6,'("dix: ",3F16.8)') (dix(kk),kk=1,3)
       !  write(6,'("grd: ",4F16.8)') (grd(kk),kk=1,4)
       !  write(6,'("inii: ",4I4)') (inii(kk,ii),kk=1,4)

       do jj=1,2

          !   polynomial interpolation

          do kk=1,3
             do ll=4,kk+1,-1
                hh(ll,jj)=(hh(ll,jj)-hh(ll-1,jj))/(grd(ll)-grd(ll-1))
             end do
          end do
          lder(4)=hh(4,jj)
          do kk=3,1,-1
             lder(kk)=hh(kk,jj)+(xx(ii)-grd(kk))*lder(kk+1)
          end do
          do kk=1,2
             do ll=3,kk+1,-1
                lder(ll)=lder(ll)+(xx(ii)-grd(ll-kk))*lder(ll+1)
             end do
          end do
          nn=ii+jj
          if (nn > 3) nn=nn-3
          ypp(ii,nn)=ypp(ii,nn)+lder(2)
          ypp(nn,ii)=ypp(nn,ii)+lder(2)
       end do
    end do

    ! averaging of the mixed derivations obtained in different order
    do ii=1,3
       do jj=1,3
          if (ii /= jj) ypp(ii,jj)=ypp(ii,jj)/2.d0
       end do
    end do

    nullify(ptddx,ptddy,ptddz,ptrho)

  end subroutine grinterp_trispline

  !> Tricubic interpolation based on: Lekien and Marsden,
  !> Int. J. Numer. Meth. Engng, 63 (2005) 455-471.  This
  !> interpolation is C^1 and local (the interpolant uses information
  !> of the grid points close to the point).  This subroutine has been
  !> adapted from the likely code, by David Kirkby, University of
  !> California, Irvine (https://github.com/deepzot/likely). f is the
  !> grid field to interpolate. xi is the point in crystallographic
  !> coordinates. y, yp, and ypp are the value, gradient, and Hessian
  !> at xi.  This routine is thread-safe.
  subroutine grinterp_tricubic(f,xi,y,yp,ypp)
    class(grid3), intent(inout), target :: f !< Input grid
    real*8, intent(in) :: xi(3) !< Target point
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    integer :: idx(3), iidx(3), i, j, k, l
    real*8 :: g(-1:2,-1:2,-1:2), x(3)
    real*8 :: a(64), b(64)
    real*8, dimension(0:3) :: aa, bb, bbx, bbxx, aax, aay, aaxy
    real*8, dimension(0:3) :: aaxx, aayy

    ! initialize
    y = 0d0
    yp = 0d0
    ypp = 0d0

    ! fetch values at the cube vertices
    x = modulo(xi,1d0)
    idx = grid_floor(f,x)
    do i = -1, 2
       do j = -1, 2
          do k = -1, 2
             iidx = modulo(idx+(/i,j,k/)-1,f%n)+1
             g(i,j,k) = f%f(iidx(1),iidx(2),iidx(3))
          end do
       end do
    end do

    ! Fill the values in the b-vector (see Lekien and Marsden)
    ! f
    b(1) = g(0,0,0)
    b(2) = g(1,0,0)
    b(3) = g(0,1,0)
    b(4) = g(1,1,0)
    b(5) = g(0,0,1)
    b(6) = g(1,0,1)
    b(7) = g(0,1,1)
    b(8) = g(1,1,1)

    ! fx
    b(9)  = 0.5d0*(g(1,0,0)-g(-1,0,0))
    b(10) = 0.5d0*(g(2,0,0)-g(0,0,0))
    b(11) = 0.5d0*(g(1,1,0)-g(-1,1,0))
    b(12) = 0.5d0*(g(2,1,0)-g(0,1,0))
    b(13) = 0.5d0*(g(1,0,1)-g(-1,0,1))
    b(14) = 0.5d0*(g(2,0,1)-g(0,0,1))
    b(15) = 0.5d0*(g(1,1,1)-g(-1,1,1))
    b(16) = 0.5d0*(g(2,1,1)-g(0,1,1))

    ! fy
    b(17) = 0.5d0*(g(0,1,0)-g(0,-1,0))
    b(18) = 0.5d0*(g(1,1,0)-g(1,-1,0))
    b(19) = 0.5d0*(g(0,2,0)-g(0,0,0))
    b(20) = 0.5d0*(g(1,2,0)-g(1,0,0))
    b(21) = 0.5d0*(g(0,1,1)-g(0,-1,1))
    b(22) = 0.5d0*(g(1,1,1)-g(1,-1,1))
    b(23) = 0.5d0*(g(0,2,1)-g(0,0,1))
    b(24) = 0.5d0*(g(1,2,1)-g(1,0,1))

    ! fz
    b(25) = 0.5d0*(g(0,0,1)-g(0,0,-1))
    b(26) = 0.5d0*(g(1,0,1)-g(1,0,-1))
    b(27) = 0.5d0*(g(0,1,1)-g(0,1,-1))
    b(28) = 0.5d0*(g(1,1,1)-g(1,1,-1))
    b(29) = 0.5d0*(g(0,0,2)-g(0,0,0))
    b(30) = 0.5d0*(g(1,0,2)-g(1,0,0))
    b(31) = 0.5d0*(g(0,1,2)-g(0,1,0))
    b(32) = 0.5d0*(g(1,1,2)-g(1,1,0))

    ! fxy
    b(33) = 0.25d0*(g(1,1,0)-g(-1,1,0)-g(1,-1,0)+g(-1,-1,0))
    b(34) = 0.25d0*(g(2,1,0)-g(0,1,0)-g(2,-1,0)+g(0,-1,0))
    b(35) = 0.25d0*(g(1,2,0)-g(-1,2,0)-g(1,0,0)+g(-1,0,0))
    b(36) = 0.25d0*(g(2,2,0)-g(0,2,0)-g(2,0,0)+g(0,0,0))
    b(37) = 0.25d0*(g(1,1,1)-g(-1,1,1)-g(1,-1,1)+g(-1,-1,1))
    b(38) = 0.25d0*(g(2,1,1)-g(0,1,1)-g(2,-1,1)+g(0,-1,1))
    b(39) = 0.25d0*(g(1,2,1)-g(-1,2,1)-g(1,0,1)+g(-1,0,1))
    b(40) = 0.25d0*(g(2,2,1)-g(0,2,1)-g(2,0,1)+g(0,0,1))

    ! fxz
    b(41) = 0.25d0*(g(1,0,1)-g(-1,0,1)-g(1,0,-1)+g(-1,0,-1))
    b(42) = 0.25d0*(g(2,0,1)-g(0,0,1)-g(2,0,-1)+g(0,0,-1))
    b(43) = 0.25d0*(g(1,1,1)-g(-1,1,1)-g(1,1,-1)+g(-1,1,-1))
    b(44) = 0.25d0*(g(2,1,1)-g(0,1,1)-g(2,1,-1)+g(0,1,-1))
    b(45) = 0.25d0*(g(1,0,2)-g(-1,0,2)-g(1,0,0)+g(-1,0,0))
    b(46) = 0.25d0*(g(2,0,2)-g(0,0,2)-g(2,0,0)+g(0,0,0))
    b(47) = 0.25d0*(g(1,1,2)-g(-1,1,2)-g(1,1,0)+g(-1,1,0))
    b(48) = 0.25d0*(g(2,1,2)-g(0,1,2)-g(2,1,0)+g(0,1,0))

    ! fyz
    b(49) = 0.25d0*(g(0,1,1)-g(0,-1,1)-g(0,1,-1)+g(0,-1,-1))
    b(50) = 0.25d0*(g(1,1,1)-g(1,-1,1)-g(1,1,-1)+g(1,-1,-1))
    b(51) = 0.25d0*(g(0,2,1)-g(0,0,1)-g(0,2,-1)+g(0,0,-1))
    b(52) = 0.25d0*(g(1,2,1)-g(1,0,1)-g(1,2,-1)+g(1,0,-1))
    b(53) = 0.25d0*(g(0,1,2)-g(0,-1,2)-g(0,1,0)+g(0,-1,0))
    b(54) = 0.25d0*(g(1,1,2)-g(1,-1,2)-g(1,1,0)+g(1,-1,0))
    b(55) = 0.25d0*(g(0,2,2)-g(0,0,2)-g(0,2,0)+g(0,0,0))
    b(56) = 0.25d0*(g(1,2,2)-g(1,0,2)-g(1,2,0)+g(1,0,0))

    ! fxyz
    b(57) = 0.125d0*(g(1,1,1)-g(-1,1,1)-g(1,-1,1)+g(-1,-1,1)-g(1,1,-1)+g(-1,1,-1)+g(1,-1,-1)-g(-1,-1,-1))
    b(58) = 0.125d0*(g(2,1,1)-g(0,1,1)-g(2,-1,1)+g(0,-1,1)-g(2,1,-1)+g(0,1,-1)+g(2,-1,-1)-g(0,-1,-1))
    b(59) = 0.125d0*(g(1,2,1)-g(-1,2,1)-g(1,0,1)+g(-1,0,1)-g(1,2,-1)+g(-1,2,-1)+g(1,0,-1)-g(-1,0,-1))
    b(60) = 0.125d0*(g(2,2,1)-g(0,2,1)-g(2,0,1)+g(0,0,1)-g(2,2,-1)+g(0,2,-1)+g(2,0,-1)-g(0,0,-1))
    b(61) = 0.125d0*(g(1,1,2)-g(-1,1,2)-g(1,-1,2)+g(-1,-1,2)-g(1,1,0)+g(-1,1,0)+g(1,-1,0)-g(-1,-1,0))
    b(62) = 0.125d0*(g(2,1,2)-g(0,1,2)-g(2,-1,2)+g(0,-1,2)-g(2,1,0)+g(0,1,0)+g(2,-1,0)-g(0,-1,0))
    b(63) = 0.125d0*(g(1,2,2)-g(-1,2,2)-g(1,0,2)+g(-1,0,2)-g(1,2,0)+g(-1,2,0)+g(1,0,0)-g(-1,0,0))
    b(64) = 0.125d0*(g(2,2,2)-g(0,2,2)-g(2,0,2)+g(0,0,2)-g(2,2,0)+g(0,2,0)+g(2,0,0)-g(0,0,0))

    ! calculate the coefficient vector
    a = matmul(c,b)

    ! interpolation in integer coordinates
    x = (x * f%n - (idx-1))

    l = 1 ! packed coefficient vector index
    do k = 0, 3
       do j = 0, 3
          ! horner's rule on x
          bb(j) = a(l) + x(1) * (a(l+1) + x(1) * (a(l+2) + x(1) * a(l+3)))
          bbx(j) = a(l+1) + x(1) * (2d0 * a(l+2) +  x(1) * 3d0 * a(l+3))
          bbxx(j) = 2d0 * a(l+2) + 6d0 * x(1) * a(l+3)

          ! advance
          l = l + 4
       end do
       ! horner's rule on y
       aa(k) = bb(0) + x(2) * (bb(1) + x(2) * (bb(2) + x(2) * bb(3)))

       aax(k) = bbx(0) + x(2) * (bbx(1) + x(2) * (bbx(2) + x(2) * bbx(3)))
       aay(k) = bb(1) + x(2) * (2d0 * bb(2) + x(2) * 3d0 * bb(3))
       aaxy(k) = bbx(1) + x(2) * (2d0 * bbx(2) + x(2) * 3d0 * bbx(3))

       aaxx(k) = bbxx(0) + x(2) * (bbxx(1) + x(2) * (bbxx(2) + x(2) * bbxx(3)))
       aayy(k) = 2d0 * bb(2) + 6d0 * x(2) * bb(3)
    end do

    ! field value
    y = aa(0) + x(3) * (aa(1) + x(3) * (aa(2) + x(3) * aa(3)))

    ! gradient
    yp(1) = aax(0) + x(3) * (aax(1) + x(3) * (aax(2) + x(3) * aax(3)))
    yp(2) = aay(0) + x(3) * (aay(1) + x(3) * (aay(2) + x(3) * aay(3)))
    yp(3) = aa(1) + x(3) * (2d0 * aa(2) + x(3) * 3d0 * aa(3))

    ! hessian
    ypp(1,1) = aaxx(0) + x(3) * (aaxx(1) + x(3) * (aaxx(2) + x(3) * aaxx(3)))
    ypp(1,2) = aaxy(0) + x(3) * (aaxy(1) + x(3) * (aaxy(2) + x(3) * aaxy(3)))
    ypp(1,3) = aax(1) + x(3) * (2d0 * aax(2) + x(3) * 3d0 * aax(3))
    ypp(2,2) = aayy(0) + x(3) * (aayy(1) + x(3) * (aayy(2) + x(3) * aayy(3)))
    ypp(2,3) = aay(1) + x(3) * (2d0 * aay(2) + x(3) * 3d0 * aay(3))
    ypp(3,3) = 2d0 * aa(2) + 6d0 * x(3) * aa(3)

    ! transform back to fractional coordinates and fill the Hessian
    do i = 1, 3
       yp(i) = yp(i) * f%n(i)
       do j = i, 3
          ypp(i,j) = ypp(i,j) * f%n(i) * f%n(j)
          ypp(j,i) = ypp(i,j)
       end do
    end do

  end subroutine grinterp_tricubic

  !> Smoothrho interpolation. The grid f is interpolated at point xi
  !> (cryst. coords.). Returns the interpolated field (y), first (yp),
  !> and second (ypp) derivatives, all in Cartesian coordinates.
  module subroutine grinterp_smr(f,xi,y,yp,ypp)
    use tools_io, only: string
    use types, only: realloc
    use param, only: icrd_crys, VSMALL
    class(grid3), intent(inout), target :: f !< Input grid
    real*8, intent(in) :: xi(3) !< Target point
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    integer :: i, j, k
    real*8 :: x1(3), xh(3)
    integer :: i0(3), ih(3)
    real*8 :: d, dd, ff, fp, fpp
    real*8, allocatable :: flist(:), w(:), f0list(:)
    real*8 :: yb, ypb(3), yppb(3,3)
    real*8 :: ypro, yppro(3), ypppro(3,3)
    real*8 :: y1, yp1(3), ypp1(3,3)
    real*8 :: yprom, ybm
    integer, allocatable :: i0list(:,:), eid(:)
    real*8, allocatable :: dist(:), wei(:), weip(:), weipp(:)
    real*8, allocatable :: ys(:), ysp(:,:), yspp(:,:,:)
    real*8, allocatable :: q(:), qp(:,:), qpp(:,:,:)
    real*8 :: u, up(3), upp(3,3)
    real*8 :: x0(3), dfin, swei, sweip(3), sweipp(3,3)
    integer :: nat, ierr, lvec(3), ii0

    !$omp critical (checkalloc)
    if (.not.allocated(f%smr_ilist)) then
       call init_smr(f)
    end if
    !$omp end critical (checkalloc)

    ! reserve memory
    allocate(flist(f%smr_nlist+4),w(f%smr_nlist+4),f0list(f%smr_nlist+4))

    ! calculate the contributing grid nodes and weights
    x0 = (xi - floor(xi)) * f%n
    if (f%smr_fdmax >= 1d0) then
       dfin = f%smr_fdmax * f%dmax
       call f%env%list_near_atoms(x0,icrd_crys,.true.,nat,ierr,eid=eid,dist=dist,lvec=lvec,up2d=dfin)
       allocate(i0list(3,nat),wei(nat),weip(nat),weipp(nat))
       do i = 1, nat
          i0list(:,i) = nint(f%env%xr2x(f%env%at(eid(i))%x) + lvec)
          call weifun(dist(i),dfin,wei(i),weip(i),weipp(i))
       end do
    else
       call f%env%list_near_atoms(x0,icrd_crys,.true.,nat,ierr,eid=eid,dist=dist,lvec=lvec,up2n=1)
       allocate(i0list(3,nat),wei(nat),weip(nat),weipp(nat))
       i0list(:,1) = nint(f%env%xr2x(f%env%at(eid(1))%x) + lvec)
       wei(1) = 1d0
       weip(1) = 0d0
       weipp(1) = 0d0
    end if

    ! run the interpolations
    allocate(ys(nat),ysp(3,nat),yspp(3,3,nat))
    allocate(q(nat),qp(3,nat),qpp(3,3,nat))
    ys = 0d0
    ysp = 0d0
    yspp = 0d0
    do ii0 = 1, nat
       ! the grid point
       i0 = i0list(:,ii0)

       ! get the rho/rho0 values at the grid points
       do i = 1, f%smr_nlist
          ih = i0 + f%smr_ilist(:,i)
          ih = modulo(ih,f%n) + 1
          flist(i) = log(max(f%f(ih(1),ih(2),ih(3)),VSMALL)/f%smr_rho0(ih(1),ih(2),ih(3)))
       end do
       flist(f%smr_nlist+1:) = 0d0

       ! solve the system of equations
       w = matmul(f%smr_phiinv,flist)

       ! calculate the rho/rho0 interpolation
       x1 = matmul(f%x2cg,xi * f%n - i0)
       y1 = 0d0
       yp1 = 0d0
       ypp1 = 0d0
       do i = 1, f%smr_nlist
          xh = x1 - f%smr_xlist(:,i)
          dd = dot_product(xh,xh)
          d = max(sqrt(dd),1d-40)

          call smr_kernelfun(d,smr_kkern,ff,fp,fpp)
          y1 = y1 + w(i) * ff
          yp1 = yp1 + w(i) * fp * xh / d
          ypp1(1,1) = ypp1(1,1) + w(i) * (fpp * (xh(1)/d)**2 + fp * (1/d - xh(1)**2/d**3))
          ypp1(2,2) = ypp1(2,2) + w(i) * (fpp * (xh(2)/d)**2 + fp * (1/d - xh(2)**2/d**3))
          ypp1(3,3) = ypp1(3,3) + w(i) * (fpp * (xh(3)/d)**2 + fp * (1/d - xh(3)**2/d**3))
          ypp1(1,2) = ypp1(1,2) + w(i) * (fpp * (xh(1)*xh(2)/d**2) - fp * xh(1)*xh(2)/d**3)
          ypp1(1,3) = ypp1(1,3) + w(i) * (fpp * (xh(1)*xh(3)/d**2) - fp * xh(1)*xh(3)/d**3)
          ypp1(2,3) = ypp1(2,3) + w(i) * (fpp * (xh(2)*xh(3)/d**2) - fp * xh(2)*xh(3)/d**3)
          ypp1(2,1) = ypp1(1,2)
          ypp1(3,1) = ypp1(1,3)
          ypp1(3,2) = ypp1(2,3)
       end do
       y1 = y1 + w(f%smr_nlist+1) + w(f%smr_nlist+2) * x1(1) + w(f%smr_nlist+3) * x1(2) + w(f%smr_nlist+4) * x1(3)
       yp1(1) = yp1(1) + w(f%smr_nlist+2)
       yp1(2) = yp1(2) + w(f%smr_nlist+3)
       yp1(3) = yp1(3) + w(f%smr_nlist+4)

       ! the promolecular density & derivatives
       call f%atenv%promolecular(xi,icrd_crys,ypro,yppro,ypppro,2)

       ! unroll the smoothing function
       if (ypro < 1d-40 .or. y1 > 100d0) then
          yb = 0d0
          ypb = 0d0
          yppb = 0d0
       else
          yb = exp(y1) * ypro
          yprom = max(ypro,1d-40)
          ybm = max(yb,1d-40)

          ypb(1) = yb * (yp1(1) + yppro(1) / yprom)
          ypb(2) = yb * (yp1(2) + yppro(2) / yprom)
          ypb(3) = yb * (yp1(3) + yppro(3) / yprom)
          yppb = 0d0
          yppb(1,1) = yb * (ypp1(1,1) + ypb(1)*ypb(1) / (ybm*ybm) + ypppro(1,1) / yprom - yppro(1)*yppro(1)/(yprom*yprom))
          yppb(1,2) = yb * (ypp1(1,2) + ypb(1)*ypb(2) / (ybm*ybm) + ypppro(1,2) / yprom - yppro(1)*yppro(2)/(yprom*yprom))
          yppb(1,3) = yb * (ypp1(1,3) + ypb(1)*ypb(3) / (ybm*ybm) + ypppro(1,3) / yprom - yppro(1)*yppro(3)/(yprom*yprom))
          yppb(2,2) = yb * (ypp1(2,2) + ypb(2)*ypb(2) / (ybm*ybm) + ypppro(2,2) / yprom - yppro(2)*yppro(2)/(yprom*yprom))
          yppb(2,3) = yb * (ypp1(2,3) + ypb(2)*ypb(3) / (ybm*ybm) + ypppro(2,3) / yprom - yppro(2)*yppro(3)/(yprom*yprom))
          yppb(3,3) = yb * (ypp1(3,3) + ypb(3)*ypb(3) / (ybm*ybm) + ypppro(3,3) / yprom - yppro(3)*yppro(3)/(yprom*yprom))
          yppb(2,1) = yppb(1,2)
          yppb(3,1) = yppb(1,3)
          yppb(3,2) = yppb(2,3)
       end if

       ! save these values for later
       ys(ii0) = yb
       ysp(:,ii0) = ypb
       yspp(:,:,ii0) = yppb

       ! weights and derivatives
       d = max(dist(ii0),1d-80)
       dd = d * d
       q(ii0) = wei(ii0)
       qp(:,ii0) = weip(ii0) * x1 / d
       qpp(1,1,ii0) = weip(ii0) / d + x1(1)*x1(1) / dd * (weipp(ii0) - weip(ii0) / d)
       qpp(2,2,ii0) = weip(ii0) / d + x1(2)*x1(2) / dd * (weipp(ii0) - weip(ii0) / d)
       qpp(3,3,ii0) = weip(ii0) / d + x1(3)*x1(3) / dd * (weipp(ii0) - weip(ii0) / d)
       qpp(1,2,ii0) = x1(1)*x1(2) / dd * (weipp(ii0) - weip(ii0) / d)
       qpp(1,3,ii0) = x1(1)*x1(3) / dd * (weipp(ii0) - weip(ii0) / d)
       qpp(2,3,ii0) = x1(2)*x1(3) / dd * (weipp(ii0) - weip(ii0) / d)
       qpp(2,1,ii0) = qpp(1,2,ii0)
       qpp(3,1,ii0) = qpp(1,3,ii0)
       qpp(3,2,ii0) = qpp(2,3,ii0)
    end do

    ! sums of weights
    swei = sum(q)
    do i = 1, 3
       sweip(i) = sum(qp(i,:))
       do j = 1, 3
          sweipp(i,j) = sum(qpp(i,j,:))
       end do
    end do

    ! put everything together
    y = 0d0
    yp = 0d0
    ypp = 0d0
    do i = 1, nat
       u = wei(i) / swei
       up(1) = (qp(1,i) - u * sweip(1)) / swei
       up(2) = (qp(2,i) - u * sweip(2)) / swei
       up(3) = (qp(3,i) - u * sweip(3)) / swei

       do j = 1, 3
          do k = j, 3
             upp(j,k) = qpp(j,k,i) / swei - qp(j,i) * sweip(k) / swei**2 - up(k) * sweip(j) / swei - &
                u * (sweipp(j,k) / swei - sweip(j) * sweip(k) / swei**2)
             upp(k,j) = upp(j,k)
          end do
       end do

       y = y + u * ys(i)
       yp = yp + up * ys(i) + u * ysp(:,i)
       do j = 1, 3
          do k = j, 3
             ypp(j,k) = ypp(j,k) + upp(j,k) * ys(i) + up(j) * ysp(k,i) + up(k) * ysp(j,i) + u * yspp(j,k,i)
             ypp(k,j) = ypp(j,k)
          end do
       end do
    end do

  contains
    subroutine weifun(x,a,w,wp,wpp)
      real*8, intent(in) :: x, a
      real*8, intent(out) :: w, wp, wpp

      real*8 :: x7, x4, x3, x2, a3, denom1, denom2

      if (x > a - 1d-14) then ! good bc goes to zero pretty quickly
         w = 0d0
      else
         x2 = x * x
         x3 = x2 * x
         x4 = x3 * x
         x7 = x4 * x3
         a3 = a*a*a
         denom2 = 1d0 / (x3 - a3)
         denom1 = denom2 / a3
         w = exp(x3 * denom1)
         wp = -3d0 * x2 * denom1 * a3 * denom2 * w
         wpp = -3d0 * x2 * denom1 * a3 * denom2 * wp + (6d0 * x - 24d0 * x4 * denom2 + 18d0 * x7 * denom2*denom2) * denom1 * w
      end if

    end subroutine weifun

  end subroutine grinterp_smr

  !> Pseudo-nearest grid point of a point x (crystallographic) (only
  !> nearest in orthogonal grids). If shift, the first point in the
  !> grid has index 1 (for indexing arrays, default is .true.). If
  !> main, translate the point to the main cell (returns indices
  !> between 1 and n if shift is true, default is .true.).
  function grid_near(f,x,main,shift) result(res)
    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x(3) !< Target point (cryst. coords.)
    logical, intent(in), optional :: main !< optional translation
    logical, intent(in), optional :: shift !< optional shift
    integer :: res(3)

    logical :: shift0, main0

    shift0 = .true.
    if (present(shift)) shift0 = shift
    main0 = .true.
    if (present(main)) main0 = main

    res = nint(x * f%n)
    if (main0) res = modulo(res,f%n)
    if (shift0) res = res + 1

  end function grid_near

  !> Floor grid point of a point x in crystallographic coords. If
  !> shift, the first point in the grid has index 1 (for indexing
  !> arrays, default is .true.). If main, translate the point to the
  !> main cell (returns indices between 1 and n if shift is true,
  !> default is .true.).
  function grid_floor(f,x,main,shift) result(res)
    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x(3) !< Target point (cryst. coords.)
    logical, intent(in), optional :: main !< optional translation
    logical, intent(in), optional :: shift !< optional shift
    integer :: res(3)

    logical :: shift0, main0

    shift0 = .true.
    if (present(shift)) shift0 = shift
    main0 = .true.
    if (present(main)) main0 = main

    res = floor(x * f%n)
    if (main0) res = modulo(res,f%n)
    if (shift0) res = res + 1

  end function grid_floor

  !> ceiling grid point of a point x in crystallographic coords. If
  !> shift, the first point in the grid has index 1 (for indexing
  !> arrays, default is .true.). If main, translate the point to the
  !> main cell (returns indices between 1 and n if shift is true,
  !> default is .true.).
  function grid_ceiling(f,x,main,shift) result(res)
    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x(3) !< Target point (cryst. coords.)
    logical, intent(in), optional :: main !< optional translation
    logical, intent(in), optional :: shift !< optional shift
    integer :: res(3)

    logical :: shift0, main0

    shift0 = .true.
    if (present(shift)) shift0 = shift
    main0 = .true.
    if (present(main)) main0 = main

    res = ceiling(x * f%n)
    if (main0) res = modulo(res,f%n)
    if (shift0) res = res + 1

  end function grid_ceiling

  !> Nearest grid point of a point x (crystallographic) in Euclidean
  !> distance. If shift, the first point in the grid has index 1 (for
  !> indexing arrays, default is .true.). If main, translate the point
  !> to the main cell (returns indices between 1 and n if shift is
  !> true, default is .true.).
  function euclidean_near(f,x,main,shift) result(res)
    use param, only: icrd_crys
    class(grid3), intent(in) :: f !< Input grid
    real*8, intent(in) :: x(3) !< Target point (cryst. coords.)
    logical, intent(in), optional :: main !< optional translation
    logical, intent(in), optional :: shift !< optional shift
    integer :: res(3)

    logical :: shift0, main0
    real*8 :: x0(3)
    integer :: idum, lt(3)
    real*8 :: dd

    shift0 = .true.
    if (present(shift)) shift0 = shift
    main0 = .true.
    if (present(main)) main0 = main

    lt = floor(x)
    x0 = (x - lt) * f%n
    call f%env%nearest_atom(x0,icrd_crys,idum,dd,lvec=res)
    if (main0) then
       res = modulo(res,f%n)
    else
       res = res + lt * f%n
    end if
    if (shift0) res = res + 1

  end function euclidean_near

  !> Initialize the geometry variables and the environment for the
  !> calculation of distances from the x2c matrix of the crystal and
  !> the number of points in each direction. Sets the variables...
  subroutine init_geometry(f,x2c,n,env)
    use tools, only: wscell
    use tools_math, only: matinv
    class(grid3), intent(inout) :: f
    real*8, intent(in) :: x2c(3,3)
    integer, intent(in) :: n(3)
    type(environ), intent(in), target :: env

    integer :: i
    real*8 :: xx(3)
    real*8, parameter :: nmaxenv = 15d0

    ! number of point and crystal matrices
    f%n(:) = n
    f%x2c = x2c
    f%c2x = x2c
    call matinv(f%c2x,3)

    ! grid x2c matrix
    f%x2cg = x2c
    do i = 1, 3
       f%x2cg(:,i) = f%x2cg(:,i) / f%n(i)
    end do

    ! voronoi relevant vectors and areas
    call wscell(f%x2cg,.true.,nf=f%nvec,ineighx=f%vec,area=f%area)

    ! dmax of the grid
    f%dmax = 0d0
    do i = 1, f%nvec
       xx = matmul(f%x2cg,f%vec(:,i))
       f%dmax = max(f%dmax,norm2(xx))
    end do

    ! grid environment
    call f%env%build_lattice(f%x2cg,f%dmax*nmaxenv)
    f%c2xg = f%x2cg
    call matinv(f%c2xg,3)

    ! atom environment
    f%atenv => env

  end subroutine init_geometry

  !> Copy the geometry variables from grid g to grid f.
  subroutine copy_geometry(f,g)
    class(grid3), intent(inout) :: f
    class(grid3), intent(in) :: g

    f%n = g%n
    f%x2c = g%x2c
    f%c2x = g%c2x
    f%x2cg = g%x2cg
    f%dmax = g%dmax
    f%env = g%env
    f%c2xg = g%c2xg
    f%atenv => g%atenv

  end subroutine copy_geometry

  !> Initialize the grid for the trispline interpolation. This is a
  !> modified version of the corresponding subroutine in abinit.
  subroutine init_trispline(f)
    use tools_io, only: ferror, faterr
    class(grid3), intent(inout) :: f !< Input grid

    integer :: istat
    integer :: d, nmax, i, i1, i2
    real*8, allocatable :: l(:,:), fg(:)
    real*8 :: fprev, fnext, fuse, fone

    if (allocated(f%c2)) return

    allocate(f%c2(f%n(1),f%n(2),f%n(3),3),stat=istat)
    if (istat /= 0) &
       call ferror('init_trispline','Error allocating c2 grids',faterr)

    nmax = maxval(f%n)
    allocate(l(nmax,nmax))
    allocate(fg(nmax))

    ! cholesky decomposition of the matrix:
    ! A =
    !  ( 4 1 0 ... 0 1 )
    !  ( 1 4 1 ... 0 0 )
    !  ( 0 1 4 ... 0 0 )
    !      .......
    !  ( 0 0 ... 1 4 1 )
    !  ( 1 0 ... 0 1 4 )
    ! that is the coefficient matrix of the system that determines
    ! the x^2 coefficients of the cubic spline (array c2 / 2).
    ! L is a lower-triangular matrix, A = L * L^t

    ! direction x->y->z
    do d = 1, 3
       nmax = f%n(d)
       l = 0d0
       l(1,1) = 2d0
       l(2,1) = 1d0 / l(1,1)
       l(2,2) = sqrt(15d0) / 2d0
       l(3,2) = 1d0 / l(2,2)
       l(nmax,1) = 0.5d0
       l(nmax,2) = - 0.5d0 / sqrt(15d0)
       do i = 3, nmax-1
          l(i,i) = sqrt(4d0 - 1d0 / l(i-1,i-1)**2)
          l(i+1,i) = 1d0 / l(i,i)
          l(nmax,i) = - l(nmax,i-1) * l(i,i-1) / l(i,i)
       end do
       l(nmax,nmax-1) = (1d0 - l(nmax,nmax-2)*l(nmax-1,nmax-2)) /&
          l(nmax-1,nmax-1)
       l(nmax,nmax) = sqrt(4d0 - sum(l(nmax,1:nmax-1)**2))
       l = l / sqrt(6d0 * nmax**2)

       ! for each of the grid points in the plane, find the spline c2
       ! (c2*x^2) coefficients
       do i1 = 1,f%n(mod(d,3)+1)
          do i2 = 1,f%n(mod(d+1,3)+1)
             select case(d)
             case(1)
                fg(1:nmax) = f%f(:,i1,i2)
             case(2)
                fg(1:nmax) = f%f(i2,:,i1)
             case(3)
                fg(1:nmax) = f%f(i1,i2,:)
             end select

             ! constant terms for the system
             fone = fg(1)
             fprev = fg(nmax)
             do i = 1, nmax-1
                fnext = fg(i+1)
                fuse = fprev + fnext - fg(i) - fg(i)
                fprev = fg(i)
                fg(i) = fuse
             end do
             fg(nmax) = fprev + fone - fg(nmax) - fg(nmax)

             ! solve by direct substitution, L^t b = c
             fg(1) = fg(1) / l(1,1)
             do i = 2, nmax-1
                fg(i) = (fg(i) - l(i,i-1) * fg(i-1)) / l(i,i)
             end do
             fg(nmax) = (fg(nmax) - &
                dot_product(l(nmax,1:nmax-1),fg(1:nmax-1))) / l(nmax,nmax)

             ! again, L a = b
             fg(nmax) = fg(nmax) / l(nmax,nmax)
             fg(nmax-1) = (fg(nmax-1) - l(nmax,nmax-1) * fg(nmax)) / l(nmax-1,nmax-1)
             do i = nmax-2, 1, -1
                fg(i) = (fg(i) - l(i+1,i) * fg(i+1) - l(nmax,i)*fg(nmax)) / l(i,i)
             end do

             ! write down the second derivatives
             select case(d)
             case(1)
                f%c2(:,i1,i2,d) = fg(1:nmax)
             case(2)
                f%c2(i2,:,i1,d) = fg(1:nmax)
             case(3)
                f%c2(i1,i2,:,d) = fg(1:nmax)
             end select

          end do
       end do

    end do

    deallocate(l,fg)

  end subroutine init_trispline

  !> Initialize the grid for smoothrho interpolation.
  subroutine init_smr(f)
    use environmod, only: environ
    use tools_math, only: matinvsym
    use types, only: realloc
    use param, only: icrd_cart
    class(grid3), intent(inout) :: f !< Input grid

    real*8 :: d, dd, xh(3), ff, fp, fpp
    integer :: i, j, nn, ierr
    real*8, allocatable :: dlist(:)
    integer, allocatable :: eid(:)

    ! do not init if already initialized
    if (allocated(f%smr_ilist)) return

    ! prepare
    if (allocated(f%smr_ilist)) deallocate(f%smr_ilist)
    if (allocated(f%smr_xlist)) deallocate(f%smr_xlist)
    if (allocated(f%smr_phiinv)) deallocate(f%smr_phiinv)

    ! calculate the stencil
    call f%env%list_near_atoms((/0d0,0d0,0d0/),icrd_cart,.true.,nn,ierr,eid,dlist,up2n=f%smr_nenv)
    f%smr_nlist = nn
    allocate(f%smr_xlist(3,nn),f%smr_ilist(3,nn))
    do i = 1, nn
       f%smr_ilist(:,i) = nint(f%env%at(eid(i))%x)
       f%smr_xlist(:,i) = matmul(f%x2cg,real(f%smr_ilist(:,i),8))
    end do

    ! calculate the inverse phi matrix
    allocate(f%smr_phiinv(f%smr_nlist+4,f%smr_nlist+4))
    f%smr_phiinv = 0d0
    do i = 1, f%smr_nlist
       do j = i, f%smr_nlist
          xh = f%smr_xlist(:,i) - f%smr_xlist(:,j)
          dd = dot_product(xh,xh)
          d = sqrt(dd)
          call smr_kernelfun(d,smr_kkern,ff,fp,fpp)
          f%smr_phiinv(i,j) = ff
          f%smr_phiinv(j,i) = f%smr_phiinv(i,j)
       end do
       f%smr_phiinv(f%smr_nlist+1,i) = 1d0
       f%smr_phiinv(f%smr_nlist+2:f%smr_nlist+4,i) = f%smr_xlist(:,i)
       f%smr_phiinv(i,f%smr_nlist+1:) = f%smr_phiinv(f%smr_nlist+1:,i)
    end do
    call matinvsym(f%smr_phiinv,f%smr_nlist+4)

   end subroutine init_smr

   !> Kernel function for smoothrho interpolation. Input: distance
   !> and degree. Output: value, first, and second derivative.
   subroutine smr_kernelfun(r,k,f,fp,fpp)
     real*8, intent(in) :: r
     integer, intent(in) :: k
     real*8, intent(out) :: f
     real*8, intent(out) :: fp
     real*8, intent(out) :: fpp

     if (r < 1d-40) then
        f = 0d0
        fp = 0d0
        fpp = 0d0
     elseif (mod(k,2) == 0) then
        f = r**k * log(r)
        fp = r**(k-1) * (k*log(r) + 1)
        fpp = (k*(k-1)*log(r) + 2*k - 1) * r**(k-2)
     elseif (mod(k,2) == 1) then
        f = r**k
        fp = k*r**(k-1)
        fpp = k*(k-1)*r**(k-2)
     end if

   end subroutine smr_kernelfun

end submodule proc
