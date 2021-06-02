!> Parts of code taken from pi7, by VLC, EFM, AMP, MFA, MBV, MAB.
!> (c) Victor Lua~na Cabal, Universidad de Oviedo, 1987--

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

!> Interface to aiPI (pi7) densities.
module pi_private
  use environmod, only: environ
  use grid1mod, only: grid1
  implicit none

  private

  type piatom
     character*(20) :: piname = ""
     logical :: pi_used = .false.
     integer :: nsym
     integer, allocatable :: naos(:)
     integer, allocatable :: naaos(:)
     integer, allocatable :: nsto(:)
     integer, allocatable :: nasto(:)
     integer, allocatable :: nn(:)
     real*8, allocatable :: z(:)
     real*8, allocatable :: xnsto(:)
     real*8, allocatable :: c(:,:)
     real*8, allocatable :: nelec(:)
     type(grid1) :: pgrid
  end type piatom

  type piwfn
     type(piatom), allocatable :: bas(:)
     real*8 :: globalcutoff = 0d0
     real*8, allocatable :: spcutoff(:,:)
     logical :: isealloc = .false.
     type(environ), pointer :: e
   contains
     procedure :: end => pi_end
     procedure :: read => pi_read
     procedure :: rho2
  end type piwfn
  public :: piwfn

  interface
     module subroutine pi_end(f)
       class(piwfn), intent(inout) :: f
     end subroutine pi_end
     module subroutine pi_read(f,nfile,piat,file,env,errmsg)
       use param, only: mlen
       class(piwfn), intent(inout) :: f
       integer, intent(in) :: nfile
       character*10, intent(in) :: piat(:)
       character(len=mlen), intent(in) :: file(:)
       type(environ), intent(in), target :: env
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine pi_read
     module subroutine rho2(f,xpos,exact,rho,grad,h)
       class(piwfn), intent(in) :: f
       real*8, intent(in) :: xpos(3)
       logical, intent(in) :: exact
       real*8, intent(out) :: rho, grad(3), h(3,3)
     end subroutine rho2
  end interface

end module pi_private

