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

  !> Type of mesh (partition function)
  integer, parameter, public :: mesh_type_becke = 1
  integer, parameter, public :: mesh_type_franchini = 2

  !> Size of mesh
  integer, parameter, public :: mesh_level_small = 1
  integer, parameter, public :: mesh_level_normal = 2
  integer, parameter, public :: mesh_level_good = 3
  integer, parameter, public :: mesh_level_vgood = 4
  integer, parameter, public :: mesh_level_amazing = 5

  !> Becke-style mesh for molecular/crystal integration. The mesh is
  !> the union of atomic grids (radial shells times Lebedev angular
  !> quadratures) multiplied by a partition function. The points are
  !> stored atom-major, shell-major: atom k owns the shells
  !> ishoff(k):ishoff(k+1)-1, and shell ish owns the points
  !> ipsh(ish):ipsh(ish+1)-1.
  type mesh
     integer :: n = 0 !< Number of mesh points
     integer :: type = mesh_type_becke !< Type of mesh
     integer :: lvl = mesh_level_good !< Level of the mesh
     real*8, allocatable :: w(:) !< Mesh weights (w = wa * wp)
     real*8, allocatable :: x(:,:) !< Cartesian coordinates of the mesh points
     real*8, allocatable :: f(:,:) !< Scalar field values on the mesh
     ! atomic-subgrid bookkeeping
     integer :: nat = 0 !< Number of atoms carrying an atomic grid
     integer, allocatable :: idat(:) !< (nat) index of the atom in c%atcel
     integer, allocatable :: zat(:) !< (nat) atomic number
     integer, allocatable :: zeffat(:) !< (nat) effective charge used for the grid choice
     logical, allocatable :: fromtable(:) !< (nat) grid from the tuned tables (true) or legacy scheme
     real*8, allocatable :: xat(:,:) !< (3,nat) center of the atomic grid
     integer, allocatable :: ishoff(:) !< (nat+1) first shell of each atom
     real*8, allocatable :: rsh(:) !< (nsh) radius of each shell
     real*8, allocatable :: wrsh(:) !< (nsh) radial weight of each shell (4*pi*r^2*dr)
     integer, allocatable :: ipsh(:) !< (nsh+1) first point of each shell
     real*8, allocatable :: wa(:) !< (n) atomic quadrature weight (radial x angular)
     real*8, allocatable :: wp(:) !< (n) partition weight
   contains
     procedure :: end => endmesh
     procedure :: gen => genmesh
     procedure :: fill => fillmesh
     procedure :: spherical_average
     procedure :: report
  end type mesh
  public :: mesh

  interface
     pure module subroutine endmesh(m)
       class(mesh), intent(inout) :: m
     end subroutine endmesh
     module subroutine genmesh(m,c,type,lvl,zpsp)
       class(mesh), intent(inout) :: m
       type(crystal), intent(inout) :: c
       integer, intent(in), optional :: type
       integer, intent(in), optional :: lvl
       integer, intent(in), optional :: zpsp(:)
     end subroutine genmesh
     module subroutine fillmesh(m,ff,prop,periodic)
       class(mesh), intent(inout) :: m
       type(field), intent(inout) :: ff
       integer, intent(in) :: prop(:)
       logical, intent(in) :: periodic
     end subroutine fillmesh
     module function spherical_average(m,iat,f) result(fr)
       class(mesh), intent(in) :: m
       integer, intent(in) :: iat
       real*8, intent(in) :: f(:)
       real*8, allocatable :: fr(:)
     end function spherical_average
     module subroutine report(m)
       class(mesh), intent(inout) :: m
     end subroutine report
  end interface

end module meshmod
