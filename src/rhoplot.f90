! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
!

!> Basic plotting capabilities: contour diagrams, 1d, 2d and 3d representations.
module rhoplot
  implicit none

  private

  public :: rhoplot_point
  public :: rhoplot_line
  public :: rhoplot_cube
  public :: rhoplot_plane
  public :: rhoplot_grdvec
  
  interface
     module subroutine rhoplot_point(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_point
     module subroutine rhoplot_line(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_line
     module subroutine rhoplot_cube(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_cube
     module subroutine rhoplot_plane(line)
       character*(*), intent(in) :: line
     end subroutine rhoplot_plane
     module subroutine rhoplot_grdvec()
     end subroutine rhoplot_grdvec
  end interface

end module rhoplot

