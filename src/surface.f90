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

!> Surface and mini-surface user-defined types and tools to work with them.
module surface
  implicit none

  private

  !> A face on the surface
  type miniface
     integer, dimension(4) :: v !< Pointers to vertices
     integer :: nv !< Number of vertices
  end type miniface

  !> Surface type
  type minisurf
     integer :: isinit !< was initialized? 0 no, 1 only v, 2 all
     integer :: rgb(3) !< optionally used for global color
     real*8 :: n(3)  !< center on which the surface is based
     real*8, allocatable :: r(:) !< vertices, r
     real*8, allocatable :: th(:) !< vertices, theta
     real*8, allocatable :: ph(:) !< vertices, phi
     type(miniface), allocatable :: f(:) !< faces
     integer :: nv !< number of vertex
     integer :: nf !< number of faces
     integer :: mv !< maximum number of vertex
     integer :: mf !< maximum number of faces
   contains
     procedure :: init => minisurf_begin
     procedure :: end => minisurf_close
     procedure :: clean => minisurf_clean
     procedure :: spheresphere
     procedure :: spheretriang
     procedure :: spherecub
     procedure :: writeint
     procedure :: readint
     procedure :: gauleg_nodes
     procedure :: lebedev_nodes
  end type minisurf
  public :: minisurf

  interface
     module subroutine minisurf_begin(s,m,f)
       class(minisurf), intent(inout) :: s
       integer, intent(in) :: m, f
     end subroutine minisurf_begin
     module subroutine minisurf_close(s)
       class(minisurf), intent(inout) :: s
     end subroutine minisurf_close
     module subroutine minisurf_clean(s)
       class(minisurf), intent(inout) :: s
     end subroutine minisurf_clean
     module subroutine spheresphere(s,xnuc,ntheta,nphi)
       class(minisurf), intent(inout) :: s
       real*8, dimension(3), intent(in) :: xnuc
       integer, intent(in) :: ntheta, nphi
     end subroutine spheresphere
     module subroutine spheretriang(s,xnuc,level)
       class(minisurf), intent(inout) :: s
       real*8, dimension(3), intent(in) :: xnuc
       integer, intent(in) :: level
     end subroutine spheretriang
     module subroutine spherecub(s,xnuc,level)
       class(minisurf), intent(inout) :: s
       real*8, dimension(3), intent(in) :: xnuc
       integer, intent(in) :: level
     end subroutine spherecub
     module subroutine writeint (s,n1,n2,meth,offfile)
       class(minisurf), intent(inout) :: s
       integer, intent(in) :: n1, n2, meth
       character*(*), intent(in) :: offfile
     end subroutine writeint
     module subroutine readint (s,n1,n2,meth,offfile,ierr)
       class(minisurf), intent(inout) :: s
       integer, intent(in) :: n1, n2, meth
       character*(*), intent(in) :: offfile
       integer, intent(out) :: ierr
     end subroutine readint
     module subroutine gauleg_nodes(srf,ntheta,nphi)
       class(minisurf), intent(inout) :: srf
       integer, intent(in) :: ntheta
       integer, intent(in) :: nphi
     end subroutine gauleg_nodes
     module subroutine lebedev_nodes(srf,npts)
       class(minisurf), intent(inout) :: srf
       integer, intent(in) :: npts
     end subroutine lebedev_nodes
  end interface

end module surface

