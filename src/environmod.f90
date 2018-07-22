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
  use types, only: anyatom, celatom, species
  implicit none

  private

  !> --- Atomic environment type. ---
  !> An atomic environment is built as the collection of all atoms whose distance is e%dmax0 or less to any point in the main
  !> cell. In molecules, the main cell is x: 0 -> 1, with x in reduced crystallographic coordinates. In crystals, 
  !> it is x: -0.5 -> 0.5. I.e. to translate to the main cell, use:
  !>   x = x - floor(x) in molecules
  !>   x = x - nint(x) in crystals
  !> 
  !> The atomic environment is encased in a box (vertices xmin/xmax), which is then divided into regions to facilitate distance
  !> searches. There are e%nregion regions, e%nreg in each direction, with indices between e%nmin and e%nmax. The regions are cubes with
  !> side length equal to boxsize. In addition to Cartesian (c), crystallographic (x), and reduced crystallographic (rx), we define
  !> two additional sets of coordinates: 
  !>   region partition (p): three integer indices that give the region to which a point belongs
  !>   region index (i): the integer index (1..e%nregion) of the encompassing region
  !> The routine c2p, c2i, and p2i convert between Cartesian and these two sets of coordinates. The origin of the c2p transformation
  !> is e%x0.
  !>
  !> The environment has env%n atoms. The first env%ncell atoms correspond to the atoms in the main cell (c%ncel). For atom i in
  !> the environment, at(i) contains:
  !>   %x - coordinates in the reduced cell (-0.5 to 0.5 if i = 1..ncell in crystals, 0 to 1 in mols).
  !>   %r - Cartesian coordinates corresponding to %x
  !>   %is - species
  !>   %idx - id from the non-equivalent list
  !>   %cidx - id from the complete list
  !>   %lvec - atcel(%cidx)%x + %lenv = xr2x(%x)
  !>
  !> To perform atomic distance calculations, the atoms are ordered by their region index (i coordinate) in the array e%imap.
  !> If k runs from 1 to env%n, c2i(at(imap(k))%r) is in ascending order. If l is a region, k1 = nrlo(l) and k2 = nrup(l)
  !> give the slice of the c2i(at(imap(1:n))%r) array corresponding to region l.
  !> 
  !> A number of region offsets (e%nregs) is stored. If l is an offset index, e%iaddregs(l) contains a packed index
  !> for the region offset (the packing/unpacking operations are handled by the packoffset and unpackoffset routines).
  !> All points in the current region are at a distance of at least e%rcutregs(l) from all points in the region given
  !> by offset e%iaddregs(l).
  type environ
     logical :: ismolecule !< Is this a molecule or a crystal?
     integer :: nspc !< Number of species
     type(species), allocatable :: spc(:) !< Species
     real*8 :: dmax0 !< Environment contains all atoms within dmax0 of every point in the cell
     real*8 :: sphmax !< The reduced cell is circumscribed by a sphere of radius sphmax
     real*8 :: boxsize !< Length of the region side (bohr)
     real*8 :: xmin(3), xmax(3) !< Encompassing box (environment)
     real*8 :: xminc(3), xmaxc(3) !< Encompassing box (main cell)
     real*8 :: x0(3) !< Origin of the to-region transformation
     integer :: n = 0 !< Number of atoms in the environment
     integer :: ncell = 0 !< Number of atoms in the unit cell
     integer :: nregc(3) !< Number of regions that cover the unit cell
     integer :: nreg(3) !< Number of regions that cover the environment
     integer :: nmin(3), nmax(3) !< Minimum and maximum region id
     integer :: nregion !< Number of regions
     ! coordinate conversion matrices
     real*8 :: m_x2xr(3,3) !< cryst. -> reduced cryst.
     real*8 :: m_xr2x(3,3) !< reduced cryst. -> cryst.
     real*8 :: m_c2xr(3,3) !< cart. -> reduced cryst.
     real*8 :: m_xr2c(3,3) !< reduced cryst. -> cart.
     real*8 :: m_c2x(3,3) !< cart. -> cryst.
     real*8 :: m_x2c(3,3) !< cryst. -> cart.
     ! atom/region mappings 
     integer, allocatable :: imap(:) !< atoms ordered by region, c2i(at(imap(1->n))%r) is ordered
     integer, allocatable :: nrlo(:) !< nrlo(ireg) = i, at(imap(i)) is the first atom in region ireg
     integer, allocatable :: nrhi(:) !< nrlo(ireg) = i, at(imap(i)) is the last atom in region ireg
     ! offset for region search
     integer :: rs_imax !< all offsets from -imax to +imax are included
     integer :: rs_2imax1 !< 2*rs_imax + 1
     real*8 :: rs_dmax !< max interaction dist. in -imax->imax cube (rs_2imax1 * sqrt(3)/2 * boxsize)
     integer :: rs_nreg !< Number of region offsets (rs_2imax1**3)
     integer, allocatable :: rs_ioffset(:) !< packed offsets for the search
     real*8, allocatable :: rs_rcut(:) !< All atoms in the cell with offset i are at a distance of at least rs_cut(i) from the main cell
     ! atom structure
     type(anyatom), allocatable :: at(:) !< Atoms (first ncell in the main cell)
   contains
     procedure :: end => environ_end !< Deallocate arrays and nullify variables
     ! initialization routines
     procedure :: build_mol => environ_build_from_molecule !< Build an environment from molecule input
     procedure :: build_crys => environ_build_from_crystal !< Build an environment from molecule input
     ! coordinate transformations
     procedure :: xr2c !< Reduced crystallographic to Cartesian
     procedure :: c2xr !< Cartesian to reduced crystallographic
     procedure :: xr2x !< Reduced crystallographic to crystallographic
     procedure :: x2xr !< Crystallographic to reduced crystallographic
     procedure :: c2x  !< Cartesian to crystallographic
     procedure :: x2c  !< Crystallographic to Cartesian
     procedure :: y2z !< Any to any (c,x,xr)
     procedure :: y2z_center !< Any to any (c,x,xr) and move to the center of the environment
     procedure :: c2p !< Cartesian to region
     procedure :: p2i !< Region to integer region
     procedure :: c2i !< Cartesian to integer region
     ! calculation routines
     procedure :: nearest_atom !< Returns the ID of the atom nearest to a given point
     procedure :: list_near_atoms !< Returns a list of atoms nearest to a given point
     procedure :: promolecular !< Calculates the promolecular or core density at a point
     ! utility routines
     procedure :: report => environ_report !< Write a report to stdout about the environment.
  end type environ
  public :: environ

  ! module procedure interfaces
  interface
     module subroutine environ_end(e)
       class(environ), intent(inout) :: e
     end subroutine environ_end
     module subroutine environ_build_from_molecule(e,nspc,spc,n,at,m_xr2c,m_x2xr,m_x2c)
       class(environ), intent(inout) :: e
       integer, intent(in) :: nspc
       type(species), intent(in) :: spc(nspc)
       integer, intent(in) :: n
       type(celatom), intent(in) :: at(n)
       real*8, intent(in) :: m_xr2c(3,3)
       real*8, intent(in) :: m_x2xr(3,3)
       real*8, intent(in) :: m_x2c(3,3)
     end subroutine environ_build_from_molecule
     module subroutine environ_build_from_crystal(e,nspc,spc,n,at,m_xr2c,m_x2xr,m_x2c,dmax0)
       class(environ), intent(inout) :: e
       integer, intent(in) :: nspc
       type(species), intent(in) :: spc(nspc)
       integer, intent(in) :: n
       type(celatom), intent(in) :: at(n)
       real*8, intent(in) :: m_xr2c(3,3)
       real*8, intent(in) :: m_x2xr(3,3)
       real*8, intent(in) :: m_x2c(3,3)
       real*8, intent(in), optional :: dmax0
     end subroutine environ_build_from_crystal
     pure module function xr2c(e,xx) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function xr2c
     pure module function c2xr(e,xx) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function c2xr
     pure module function xr2x(e,xx) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function xr2x
     pure module function x2xr(e,xx) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function x2xr
     pure module function c2x(e,xx) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function c2x
     pure module function x2c(e,xx) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       real*8 :: res(3)
     end function x2c
     pure module function y2z(e,xx,icrd,ocrd) result(res)
       class(environ), intent(in) :: e
       real*8, intent(in)  :: xx(3)
       integer, intent(in) :: icrd, ocrd
       real*8 :: res(3)
     end function y2z
     pure module subroutine y2z_center(e,xx,icrd,ocrd,lvec)
       class(environ), intent(in) :: e
       real*8, intent(inout) :: xx(3)
       integer, intent(in) :: icrd
       integer, intent(in) :: ocrd
       integer, intent(out), optional :: lvec(3)
     end subroutine y2z_center
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
     module subroutine nearest_atom(e,xp,icrd,nid,dist,lvec,nid0,id0,nozero)
       class(environ), intent(in) :: e
       real*8, intent(in) :: xp(3)
       integer, intent(in) :: icrd
       integer, intent(out) :: nid
       real*8, intent(out) :: dist
       integer, intent(out), optional :: lvec(3)
       integer, intent(in), optional :: nid0
       integer, intent(in), optional :: id0
       logical, intent(in), optional :: nozero
     end subroutine nearest_atom
     module subroutine list_near_atoms(e,xp,icrd,sorted,nat,eid,dist,lvec,ishell0,up2d,up2dsp,up2sh,up2n,nid0,id0,nozero)
       use param, only: icrd_rcrys
       class(environ), intent(in) :: e
       real*8, intent(in) :: xp(3)
       integer, intent(in) :: icrd
       logical, intent(in) :: sorted
       integer, intent(out) :: nat
       integer, allocatable, intent(inout) :: eid(:)
       real*8, allocatable, intent(inout) :: dist(:)
       integer, intent(out) :: lvec(3)
       integer, allocatable, intent(inout), optional :: ishell0(:)
       real*8, intent(in), optional :: up2d
       real*8, intent(in), optional :: up2dsp(:)
       integer, intent(in), optional :: up2sh
       integer, intent(in), optional :: up2n
       integer, intent(in), optional :: nid0
       integer, intent(in), optional :: id0
       logical, intent(in), optional :: nozero
     end subroutine list_near_atoms
     module subroutine promolecular(e,x0,icrd,f,fp,fpp,nder,zpsp,fr)
       use fragmentmod, only: fragment
       class(environ), intent(in) :: e
       real*8, intent(in) :: x0(3) 
       integer, intent(in) :: icrd
       real*8, intent(out) :: f
       real*8, intent(out) :: fp(3)
       real*8, intent(out) :: fpp(3,3)
       integer, intent(in) :: nder
       integer, intent(in), optional :: zpsp(:) 
       type(fragment), intent(in), optional :: fr
     end subroutine promolecular
     module subroutine environ_report(e)
       class(environ), intent(in) :: e
     end subroutine environ_report
  end interface

end module environmod
