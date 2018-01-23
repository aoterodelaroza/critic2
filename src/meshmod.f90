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

! meshmod: becke-style meshes for molecular integrals.
! The meshmod module has been adapted from postg: 
! Copyright (c) 2013 Alberto Otero de la Roza
! <aoterodelaroza@ucmerced.edu>, Felix Kannemann
! <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider
! <hs7@post.queensu.ca>, and Axel D. Becke <axel.becke@dal.ca>
module meshmod
  use fieldmod, only: field
  use crystalmod, only: crystal
  implicit none

  private

  !> Becke-style mesh for molecular/crystal integration
  type mesh
     integer :: n = 0 !< Number of mesh points
     integer :: type !< Type of mesh
     real*8, allocatable :: w(:) !< Mesh weights
     real*8, allocatable :: x(:,:) !< Cartesian coordinates of the mesh points
     real*8, allocatable :: f(:,:) !< Scalar field values on the mesh
   contains
     procedure :: end => endmesh
     procedure :: gen => genmesh
     procedure :: gen_becke => genmesh_becke
     procedure :: gen_franchini => genmesh_franchini
     procedure :: fill => fillmesh
  end type mesh
  public :: mesh
  
  interface
     module subroutine endmesh(m)
       class(mesh), intent(inout) :: m
     end subroutine endmesh
     module subroutine genmesh(m,c,type)
       class(mesh), intent(inout) :: m
       type(crystal), intent(inout) :: c
       integer, intent(in), optional :: type
     end subroutine genmesh
     module subroutine genmesh_becke(m,c)
       class(mesh), intent(inout) :: m
       type(crystal), intent(in) :: c
     end subroutine genmesh_becke
     module subroutine genmesh_franchini(m,c,lvl)
       class(mesh), intent(inout) :: m
       type(crystal), intent(in) :: c
       integer, intent(in) :: lvl
     end subroutine genmesh_franchini
     module subroutine fillmesh(m,ff,prop,periodic)
       class(mesh), intent(inout) :: m
       type(field), intent(inout) :: ff
       integer, intent(in) :: prop(:)
       logical, intent(in) :: periodic
     end subroutine fillmesh
  end interface

end module meshmod
