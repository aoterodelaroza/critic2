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

! Evaluation of arithmetic expressions. This entire module is
! THREAD-SAFE (or should be, at least).
module arithmetic
  implicit none

  public :: eval
  public :: eval_grid
  public :: fields_in_eval
  public :: setvariable
  public :: isvariable
  public :: clearvariable
  public :: clearallvariables
  public :: listvariables
  public :: listlibxc

  private

  ! module procedure interfaces
  interface
     recursive module function eval(expr,hardfail,iok,x0,sptr,periodic)
       use iso_c_binding, only: c_ptr
       real*8 :: eval
       character(*), intent(in) :: expr
       logical, intent(in) :: hardfail
       logical, intent(out) :: iok
       real*8, intent(in), optional :: x0(3)
       type(c_ptr), intent(in), optional :: sptr
       logical, intent(in), optional :: periodic
     end function eval
     module subroutine eval_grid(n,expr,sptr,f,iok)
       use iso_c_binding, only: c_ptr
       integer, intent(in) :: n(3)
       character(*), intent(in) :: expr
       type(c_ptr), intent(in) :: sptr
       real*8, intent(out) :: f(n(1),n(2),n(3))
       logical, intent(out) :: iok
     end subroutine eval_grid
     module subroutine fields_in_eval(expr,n,idlist,sptr)
       use iso_c_binding, only: c_ptr
       use param, only: mlen
       character(*), intent(in) :: expr
       integer, intent(out) :: n
       character(len=mlen), allocatable, intent(inout) :: idlist(:)
       type(c_ptr), intent(in) :: sptr
     end subroutine fields_in_eval
     module subroutine setvariable(ikey,ival)
       character*(*), intent(in) :: ikey
       real*8, intent(in) :: ival
     end subroutine setvariable
     module function isvariable(ikey,ival)
       character*(*), intent(in) :: ikey
       real*8, intent(out) :: ival
       logical :: isvariable
     end function isvariable
     module subroutine clearvariable(ikey)
       character*(*), intent(in) :: ikey
     end subroutine clearvariable
     module subroutine clearallvariables()
     end subroutine clearallvariables
     module subroutine listvariables()
     end subroutine listvariables
     module subroutine listlibxc(doref,doname,doflags)
       logical, intent(in) :: doref
       logical, intent(in) :: doname
       logical, intent(in) :: doflags
     end subroutine listlibxc
  end interface

end module arithmetic
