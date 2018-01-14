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

!> Code for the FLUXPRINT environment: 3d plotting.
module flux
  implicit none

  private

  !xx! list of routines
  public  :: fluxprint
  private :: flx_initialize
  public  :: flx_end
  private :: flx_printpath
  private :: flx_point
  private :: flx_ncp
  private :: flx_bcp
  private :: flx_graph
  private :: flx_findthetagrid
  private :: ball_radius

  interface
     module subroutine fluxprint()
     end subroutine fluxprint
     module subroutine flx_initialize(nosym,noballs,nocell)
       logical, intent(in) :: nosym
       logical, intent(in), optional :: noballs, nocell
     end subroutine flx_initialize
     module subroutine flx_end(molmotif)
       logical, intent(in) :: molmotif
     end subroutine flx_end
     module subroutine flx_printpath(rgb0)
       integer, intent(in) :: rgb0(3)
     end subroutine flx_printpath
     module subroutine flx_symprintpath(x,flxsym,rgb)
       integer, intent(in) :: flxsym
       real*8, dimension(3), intent(in) :: x
       integer, intent(in) :: rgb(3)
     end subroutine flx_symprintpath
     module subroutine flx_point(x,iup,flxsym,rgb)
       integer, intent(in) :: iup
       real*8, dimension(3), intent(in) :: x
       integer, intent(in) :: flxsym
       integer, intent(in) :: rgb(3)
     end subroutine flx_point
     module subroutine flx_ncp(id,ntheta,nphi,flxsym,lvec,rgb)
       integer, intent(in) :: id, ntheta, nphi
       integer, intent(in) :: flxsym
       integer, intent(in), dimension(3), optional :: lvec
       integer, intent(in) :: rgb(3)
     end subroutine flx_ncp
     module subroutine flx_bcp(id,iup,npoints,flxsym,bcpmethod,lvec,rgb)
       integer, intent(in) :: id, iup, npoints
       integer, intent(in) :: flxsym
       character*3, intent(in) :: bcpmethod
       integer, intent(in), dimension(3), optional :: lvec
       integer, intent(in) :: rgb(3)
     end subroutine flx_bcp
     module subroutine flx_graph(flxsym,igraph,cpid,lvec,rgb)
       integer, intent(in) :: flxsym
       integer, intent(in) :: igraph
       integer, intent(in), optional :: cpid
       integer, intent(in), dimension(3), optional :: lvec
       integer, intent(in) :: rgb(3)
     end subroutine flx_graph
     module subroutine flx_findthetagrid(lx,ly,r0,R,n,thetavec)
       real*8, intent(in) :: lx, ly, r0, R
       integer, intent(in) :: n
       real*8, dimension(n), intent(out) :: thetavec
     end subroutine flx_findthetagrid
     module function ball_radius(i,cel)
       integer, intent(in) :: i
       logical, intent(in) :: cel
       real*8 :: ball_radius
     end function ball_radius
  end interface

end module flux
