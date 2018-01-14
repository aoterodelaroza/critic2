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
!
module nci
  implicit none

  private
  public :: nciplot
  private :: write_cube_header
  private :: write_cube_body
  private :: check_no_extra_word
  private :: read_fragment
  private :: readchk
  private :: writechk
  
  interface
     module subroutine nciplot()
     end subroutine nciplot
     module subroutine write_cube_header(lu,l1,l2,periodic,nfrag,frag,x0x,x1x,x0,x1,nstep,xmat)
       use fragmentmod, only: fragment
       integer, intent(in) :: lu
       character*(*), intent(in) :: l1, l2
       logical, intent(in) :: periodic
       integer, intent(in) :: nfrag
       type(fragment), intent(in) :: frag(nfrag)
       real*8, intent(in) :: x0x(3), x1x(3), x0(3), x1(3)
       integer, intent(in) :: nstep(3)
       real*8, intent(in) :: xmat(3,3)
     end subroutine write_cube_header
     module subroutine write_cube_body(lu,n,c)
       integer, intent(in) :: lu
       integer, intent(in) :: n(3)
       real*8, intent(in) :: c(0:n(3)-1,0:n(2)-1,0:n(1)-1)
     end subroutine write_cube_body
     module function check_no_extra_word(line,lp,routine)
       character(len=:), allocatable :: aux2
       character*(*), intent(in) :: line, routine
       integer, intent(inout) :: lp
       logical :: check_no_extra_word
     end function check_no_extra_word
     module function read_fragment(lu) result(fr)
       use fragmentmod, only: fragment
       type(fragment) :: fr
       integer, intent(in) :: lu
     end function read_fragment
     module function readchk(x0,xmat,nstep,crho,cgrad,rhoat,nfrag,rhofrag) result(lchk)
       logical :: lchk
       real*8, intent(in) :: x0(3), xmat(3,3)
       integer, intent(in) :: nstep(3)
       real*8, allocatable, dimension(:,:,:), intent(inout) :: crho, cgrad, rhoat
       integer, intent(in) :: nfrag
       real*8, allocatable, dimension(:,:,:,:), intent(inout) :: rhofrag
     end function readchk
     module subroutine writechk(x0,xmat,nstep,crho,cgrad,rhoat,nfrag,rhofrag)
       real*8, intent(in) :: x0(3), xmat(3,3)
       integer, intent(in) :: nstep(3)
       real*8, allocatable, dimension(:,:,:), intent(in) :: crho, cgrad, rhoat
       integer, intent(in) :: nfrag
       real*8, allocatable, dimension(:,:,:,:), intent(in) :: rhofrag
     end subroutine writechk
  end interface

end module nci
