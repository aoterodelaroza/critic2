! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
module meshmod
  use fieldmod, only: field
  use crystalmod, only: crystal
  implicit none

  private

  !> Type of mesh
  integer, parameter, public :: mesh_type_becke = 1
  integer, parameter, public :: mesh_type_franchini = 2

  !> Size of mesh
  integer, parameter, public :: mesh_level_small = 1
  integer, parameter, public :: mesh_level_normal = 2
  integer, parameter, public :: mesh_level_good = 3
  integer, parameter, public :: mesh_level_vgood = 4
  integer, parameter, public :: mesh_level_amazing = 5

  !> Becke-style mesh for molecular/crystal integration
  type mesh
     integer :: n = 0 !< Number of mesh points
     integer :: type !< Type of mesh
     integer :: lvl !< Level of the mesh
     real*8, allocatable :: w(:) !< Mesh weights
     real*8, allocatable :: x(:,:) !< Cartesian coordinates of the mesh points
     real*8, allocatable :: f(:,:) !< Scalar field values on the mesh
   contains
     procedure :: end => endmesh
     procedure :: gen => genmesh
     procedure :: gen_becke => genmesh_becke
     procedure :: gen_franchini => genmesh_franchini
     procedure :: fill => fillmesh
     procedure :: report
  end type mesh
  public :: mesh

  interface
     module subroutine endmesh(m)
       class(mesh), intent(inout) :: m
     end subroutine endmesh
     module subroutine genmesh(m,c,type,lvl)
       class(mesh), intent(inout) :: m
       type(crystal), intent(inout) :: c
       integer, intent(in), optional :: type
       integer, intent(in), optional :: lvl
     end subroutine genmesh
     module subroutine genmesh_becke(m,c,lvl)
       class(mesh), intent(inout) :: m
       type(crystal), intent(in) :: c
       integer, intent(in) :: lvl
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
     module subroutine report(m)
       class(mesh), intent(inout) :: m
     end subroutine report
  end interface

end module meshmod
