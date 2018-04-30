! Copyright (c) 2007-2018 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.

! Atomic environment class
module environmod
  use types, only: celatom, species
  implicit none

  private

  !> Atomic environment type.
  !> The first ncell atoms correspond to the atoms in the main
  !> cell (c%ncel). The main cell atoms are centered around the
  !> origin.
  type environ
     real*8 :: dmax0 !< The environment contains all atoms at a distance of dmax0 or more to every point in the main cell
     real*8 :: boxsize !< length of the region side (bohr)
     real*8 :: xmin(3), xmax(3) !< encompassing box (environment)
     real*8 :: xminc(3), xmaxc(3) !< encompassing box (main cell)
     real*8 :: x0(3) !< origin of the to-region transformation
     integer :: n = 0 !< Number of atoms in the environment
     integer :: ncell = 0 !< Number of atoms in the unit cell
     integer :: nregc(3) !< Number of regions that cover the unit cell
     integer :: nreg(3) !< Number of regions that cover the environment
     integer :: nmin(3), nmax(3) !< Minimum and maximum region id
     integer :: nregion !< number of regions
     integer, allocatable :: nrlo(:), nrup(:) !< lower and upper bound for list of atoms in a region (applied to imap)
     integer, allocatable :: imap(:) !< atoms ordered by region, c2i(imap(1->n)) is ordered
     type(celatom), allocatable :: at(:) !< Atoms 
   contains
     procedure :: end => environ_end !< Deallocate arrays and nullify variables
     procedure :: build_mol => environ_build_from_molecule !< Build an environment from molecule input
     procedure :: build_crys => environ_build_from_crystal !< Build an environment from molecule input
     procedure :: c2p !< Cartesian to region
     procedure :: p2i !< Region to integer region
     procedure :: c2i !< Cartesian to integer region
  end type environ
  public :: environ

  ! module procedure interfaces
  interface
     module subroutine environ_end(e)
       class(environ), intent(inout) :: e
     end subroutine environ_end
     module subroutine environ_build_from_molecule(e,n,at)
       class(environ), intent(inout) :: e
       integer, intent(in) :: n
       type(celatom), intent(in) :: at(n)
     end subroutine environ_build_from_molecule
     module subroutine environ_build_from_crystal(e,nspc,spc,n,at,m_xr2c,m_x2xr,dmax0)
       class(environ), intent(inout) :: e
       integer, intent(in) :: nspc
       type(species), intent(in) :: spc(nspc)
       integer, intent(in) :: n
       type(celatom), intent(in) :: at(n)
       real*8, intent(in) :: m_xr2c(3,3)
       real*8, intent(in) :: m_x2xr(3,3)
       real*8, intent(in), optional :: dmax0
     end subroutine environ_build_from_crystal
     pure module function c2p(e,xx) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       integer :: res(3)
     end function c2p
     pure module function p2i(e,xx) result(res)
       class(environ), intent(in) :: e
       integer, intent(in)  :: xx(3)
       integer :: res
     end function p2i
     pure module function c2i(e,xx) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       integer :: res
     end function c2i
  end interface

end module environmod
