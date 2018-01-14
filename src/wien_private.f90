!> This module contains routines  adapted (with permission) from the WIEN2k
!> source by P. Blaha, K. Schwarz et al. 

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

!> Interface to WIEN2k densities.
module wien_private
  implicit none

  private

  type wienwfn
     ! field quantities
     logical :: cnorm = .true.
     integer, allocatable :: lm(:,:,:)
     integer, allocatable :: lmmax(:)
     real*8, allocatable :: slm(:,:,:)
     real*8, allocatable :: sk(:)
     real*8, allocatable :: ski(:)
     real*8, allocatable :: tauk(:)
     real*8, allocatable :: tauki(:)
     integer :: lastind
     integer :: lastniz
     ! structural
     logical :: ishlat
     integer :: nat
     real*8, allocatable  :: rotloc(:,:,:)
     real*8, allocatable :: rnot(:)
     real*8, allocatable :: rmt(:)
     real*8, allocatable :: dx(:)
     integer, allocatable :: jri(:)
     integer, allocatable :: multw(:)
     real*8 :: br1(3,3)
     real*8 :: br2(3,3)
     real*8 :: br3(3,3)
     logical :: ortho
     integer :: ndat
     integer, allocatable :: iatnr(:)
     real*8, allocatable :: pos(:,:)
     integer, allocatable :: iop(:)
     integer :: niord
     integer, allocatable :: iz(:,:,:)
     real*8, allocatable :: tau(:,:)
     integer :: npos 
     real*8  :: atp(3,343)
     integer :: nwav
     real*8, allocatable :: krec(:,:)
     integer, allocatable :: inst(:)
     logical :: cmpl
     real*8 :: a(3)
   contains
     procedure :: end => wien_end
     procedure :: rmt_atom
     procedure :: read_clmsum
     procedure :: rho2
     procedure :: tolap
  end type wienwfn
  public :: wienwfn
  
  integer, parameter  :: ncom = 122 !< maximum number of LM pairs
  integer, parameter  :: nrad = 781 !< maximum value of npt
  integer, parameter  :: nsym = 48 !< max. number of symmetry operations
  integer, parameter  :: lmax2 = 14 !< maximum l for LM lattice harmonics expansion of charte inside muffins.

  interface
     module subroutine wien_end(f)
       class(wienwfn), intent(inout) :: f
     end subroutine wien_end
     module function rmt_atom(f,x)
       class(wienwfn), intent(in) :: f
       real*8, intent(in) :: x(3)
       real*8 :: rmt_atom
     end function rmt_atom
     module subroutine wien_read_struct(f,file)
       class(wienwfn), intent(inout) :: f
       character*(*), intent(in) :: file
     end subroutine wien_read_struct
     module subroutine read_clmsum(f,file,file2)
       class(wienwfn), intent(inout) :: f
       character*(*), intent(in) :: file, file2
     end subroutine read_clmsum
     module subroutine readslm(f,lu)
       class(wienwfn), intent(inout) :: f
       integer, intent(in) :: lu !< Input logical unit
     end subroutine readslm
     module subroutine readk(f,lu)
       class(wienwfn), intent(inout) :: f
       integer, intent(in) :: lu !< Input logical unit
     end subroutine readk
     module subroutine gbass(rbas,gbas)
       real*8, dimension(3,3) :: rbas, gbas
     end subroutine gbass
     module subroutine rotdef(f)
       class(wienwfn), intent(inout) :: f
     end subroutine rotdef
     module subroutine gener(f)
       class(wienwfn), intent(inout) :: f
     end subroutine gener
     module subroutine sternb(f)
       class(wienwfn), intent(in) :: f
     end subroutine sternb
     module subroutine rotator(vt,iz,tau,a)
       real*8 :: vt(3), tau(3), a(3)
       integer :: iz(3,3)
     end subroutine rotator
     module subroutine rotato(vt,iz,tau)
       real*8 :: vt(3), tau(3)
       integer :: iz(3,3)
     end subroutine rotato
     module subroutine rotat(vt,rotloc)
       real*8 :: vt(3),rotloc(3,3)
     end subroutine rotat
     module subroutine reduc(v,atp,npos,pos,iat,rmt)
       integer :: npos, iat
       real*8 :: v(3),atp(3,*),pos(3,*), rmt
     end subroutine reduc
     module function vnorm(v)
       real*8 :: vnorm
       real*8 :: v(3)
     end function vnorm
     module subroutine charge(f,chg,grad,hess,ir,r,jatom,v,natnr)
       class(wienwfn), intent(in) :: f
       real*8, intent(out) :: chg
       real*8, dimension(3), intent(out) :: grad
       real*8, dimension(6), intent(out) :: hess
       integer, intent(in) :: ir
       real*8,intent(in) :: r
       integer, intent(in) :: jatom
       real*8, dimension(3), intent(in) :: v
       integer, intent(in) :: natnr
     end subroutine charge
     module subroutine radial(f,rlm,rho,rho1,rho2,r,ir,jatom)
       class(wienwfn), intent(in) :: f
       real*8, dimension(nrad), intent(in) :: rlm
       real*8, intent(out) :: rho, rho1, rho2
       real*8, intent(in) :: r
       integer, intent(in) :: ir, jatom
     end subroutine radial
     module subroutine rhoout(f,v,chg,grad,hess)
       class(wienwfn), intent(in) :: f
       real*8, dimension(3), intent(in) :: v
       real*8, intent(out) :: chg
       real*8, dimension(3), intent(out) :: grad
       real*8, dimension(6), intent(out) :: hess
     end subroutine rhoout
     module subroutine ylm(v,lmax,y)
       integer            lmax
       double precision   v(3)
       complex*16         y(*)
     end subroutine ylm
     module subroutine rho2(f,v0,rho,grad,h)
       class(wienwfn), intent(in) :: f
       real*8, dimension(3), intent(in) :: v0
       real*8, intent(out) :: rho, grad(3), h(3,3)
     end subroutine rho2
     module subroutine tolap(f)
       class(wienwfn), intent(inout) :: f
     end subroutine tolap
  end interface

end module wien_private

