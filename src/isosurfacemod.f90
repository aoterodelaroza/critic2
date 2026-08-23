! Copyright (c) 2015-2026 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> Isosurface triangulation of a scalar field given on a 3D grid
!> (marching cubes).
module isosurfacemod
  implicit none

  private

  public :: marching_cubes

  interface
     module subroutine marching_cubes(n,f,x2c,c2x,isoval,periodic,nv,xv,nrm,nf,idx)
       integer, intent(in) :: n(3)
       real*8, intent(in) :: f(:,:,:)
       real*8, intent(in) :: x2c(3,3)
       real*8, intent(in) :: c2x(3,3)
       real*8, intent(in) :: isoval
       logical, intent(in) :: periodic
       integer, intent(out) :: nv
       real*8, allocatable, intent(inout) :: xv(:,:)
       real*8, allocatable, intent(inout) :: nrm(:,:)
       integer, intent(out) :: nf
       integer, allocatable, intent(inout) :: idx(:,:)
     end subroutine marching_cubes
  end interface

end module isosurfacemod
