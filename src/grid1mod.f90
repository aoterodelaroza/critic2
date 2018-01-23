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

!> One-dimensional grid class
module grid1mod
  implicit none

  private

  !> Radial grid type.
  type grid1
     logical :: isinit = .false. !< Is initialized?
     real*8 :: a !< Logarithmic grid parameter ri = a * exp(b * (i-1))
     real*8 :: b !< Logarithmic grid parameter ri = a * exp(b * (i-1))
     real*8 :: rmax !< Max. grid distance
     real*8 :: rmax2 !< Squared max. grid distance
     integer :: ngrid !< Number of nodes
     real*8, allocatable :: r(:) !< Node positions
     real*8, allocatable :: f(:) !< Grid values, f = 4*pi*r^2*rho
     real*8, allocatable :: fp(:) !< First derivative of f
     real*8, allocatable :: fpp(:) !< Second derivative of f 
     integer :: z 
     integer :: qat
   contains
     procedure :: read_db !< Read a one-dimesional grid from the density tables
     procedure :: interp !< Interpolate value and derivatives from the grid
  end type grid1
  public :: grid1

  ! table of all-electron and core 1d grids
  type(grid1), target, allocatable, public :: agrid(:)
  type(grid1), target, allocatable, public :: cgrid(:,:)
  public :: grid1_register_core
  public :: grid1_register_ae
  public :: grid1_clean_grids
  
  interface
     module subroutine grid1_end(g)
       class(grid1), intent(inout) :: g
     end subroutine grid1_end
     module subroutine read_db(g,z,q)
       class(grid1), intent(inout) :: g
       integer, intent(in) :: z
       integer, intent(in) :: q
     end subroutine read_db
     module subroutine interp(g,r0,f,fp,fpp)
       class(grid1), intent(in) :: g
       real*8, intent(in) :: r0
       real*8, intent(out) :: f
       real*8, intent(out) :: fp
       real*8, intent(out) :: fpp
     end subroutine interp
     module subroutine grid1_register_core(iz,iq)
       integer, intent(in) :: iz, iq
     end subroutine grid1_register_core
     module subroutine grid1_register_ae(iz)
       integer, intent(in) :: iz
     end subroutine grid1_register_ae
     module subroutine grid1_clean_grids()
     end subroutine grid1_clean_grids
  end interface

end module grid1mod
