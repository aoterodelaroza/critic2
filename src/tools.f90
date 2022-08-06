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

! Toolbox with miscellaneous tools.
module tools
  implicit none

  private

  !xx! sorting and uniqing
  public :: uniqc
  public :: qcksort
  public :: mergesort
  public :: tiny_atom_type
  public :: delaunay_reduction
  public :: wscell

  interface qcksort
     module procedure qcksort_r8_inplace
     module procedure qcksort_i4_inplace
     module procedure qcksort_r8
     module procedure qcksort_i4
  end interface qcksort
  interface mergesort
     module procedure mergesort_r8
     module procedure mergesort_i4
  end interface mergesort

  interface
     module subroutine uniqc(v, nini, nfin, eps)
       integer, intent(in) :: nini
       integer, intent(inout) :: nfin
       real*8, intent(inout), dimension(3,*) :: v
       real*8, intent(in) :: eps
     end subroutine uniqc
     module subroutine qcksort_r8(arr, iord, first, last)
       real*8, dimension(:), intent(in) :: arr
       integer, dimension(:), intent(inout) :: iord
       integer, intent(in) :: first
       integer, intent(in) :: last
     end subroutine qcksort_r8
     module subroutine qcksort_i4(iarr, iord, first, last)
       integer, dimension(:), intent(in) :: iarr
       integer, dimension(:), intent(inout) :: iord
       integer, intent(in) :: first
       integer, intent(in) :: last
     end subroutine qcksort_i4
     module subroutine qcksort_r8_inplace(arr)
       real*8, intent(inout) :: arr(:)
     end subroutine qcksort_r8_inplace
     module subroutine qcksort_i4_inplace(arr)
       integer, intent(inout) :: arr(:)
     end subroutine qcksort_i4_inplace
     module subroutine mergesort_r8(arr,iord,ini,n)
       real*8, dimension(:), intent(in) :: arr
       integer, dimension(:), intent(inout) :: iord
       integer, intent(in) :: ini, n
     end subroutine mergesort_r8
     module subroutine mergesort_i4(arr, iord, ini, n)
       integer, dimension(:), intent(in) :: arr
       integer, dimension(:), intent(inout) :: iord
       integer, intent(in) :: ini, n
     end subroutine mergesort_i4
     module function tiny_atom_type(iz,nn)
       integer, intent(in) :: iz, nn
       integer :: tiny_atom_type
     end function tiny_atom_type
     module subroutine delaunay_reduction(x2c,rmat,rbas)
       real*8, intent(in) :: x2c(3,3)
       real*8, intent(out) :: rmat(3,4)
       real*8, intent(out), optional :: rbas(3,3)
     end subroutine delaunay_reduction
     module subroutine wscell(m_x2c,doreduction,nf,nv,mnfv,iside,nside,x,ineighc,ineighx,area,&
        isortho,m_xr2x,m_x2xr,m_xr2c,m_c2xr,n2_xr2x,n2_x2xr,n2_xr2c,n2_c2xr,ineighxr,&
        isortho_del)
       real*8, intent(in) :: m_x2c(3,3)
       logical, intent(in) :: doreduction
       integer, intent(out), optional :: nf
       integer, intent(out), optional :: nv
       integer, intent(out), optional :: mnfv
       integer, allocatable, intent(inout), optional :: iside(:,:)
       integer, allocatable, intent(inout), optional :: nside(:)
       real*8, allocatable, intent(inout), optional :: x(:,:)
       real*8, allocatable, intent(inout), optional :: ineighc(:,:)
       integer, allocatable, intent(inout), optional :: ineighx(:,:)
       real*8, allocatable, intent(inout), optional :: area(:)
       logical, intent(out), optional :: isortho
       real*8, intent(out), optional :: m_xr2x(3,3)
       real*8, intent(out), optional :: m_x2xr(3,3)
       real*8, intent(out), optional :: m_xr2c(3,3)
       real*8, intent(out), optional :: m_c2xr(3,3)
       real*8, intent(out), optional :: n2_xr2x
       real*8, intent(out), optional :: n2_x2xr
       real*8, intent(out), optional :: n2_xr2c
       real*8, intent(out), optional :: n2_c2xr
       integer, allocatable, intent(inout), optional :: ineighxr(:,:)
       logical, intent(out), optional :: isortho_del
     end subroutine wscell
  end interface

end module tools
