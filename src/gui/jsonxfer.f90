! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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


! Transfer of values between variables and the keys of a JSON object,
! in either direction, so that one list of calls both writes an object
! and reads it back (used by the scene files).
module jsonxfer
  use iso_c_binding, only: c_float
  use json_module, only: json_core, json_value
  implicit none

  private

  ! largest number of elements taken from a file for a list of objects
  ! (xfer_list) or a matrix (x_r8m_alloc)
  integer, parameter, public :: xfer_maxlist = 100000

  !> Transfer between variables and a JSON object. With reading =
  !> .false., each x call adds the value under the key; with reading =
  !> .true., it takes the value of the key, if the key is there with a
  !> value of the right type and (arrays) the same number of elements
  !> as the variable, and leaves the variable as it was otherwise.
  type json_xfer
     type(json_core), pointer :: json => null() ! the json core of the tree
     logical :: reading = .false. ! read from the tree (else write to it)
   contains
     procedure :: x_r8, x_r8v, x_r8m, x_r4, x_r4v, x_r4m, x_i4, x_i4v, x_i4m, x_i8v
     procedure :: x_l, x_lv, x_lm, x_s, x_sv
     generic :: x => x_r8, x_r8v, x_r8m, x_r4, x_r4v, x_r4m, x_i4, x_i4v, x_i4m, x_i8v,&
        x_l, x_lv, x_lm, x_s, x_sv ! a value (scalar or fixed-size array; scalars: optional range)
     procedure :: xs => x_salloc ! an allocatable string
     procedure :: xa => x_r8m_alloc ! an allocatable real*8 matrix (number of columns from the file)
     procedure :: obj => xfer_obj ! a child object
     procedure :: list => xfer_list ! an array of child objects
     procedure :: item => xfer_item ! a child object of an array
  end type json_xfer
  public :: json_xfer

  interface
     module subroutine x_r8(io,p,key,v,lo,hi)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       real*8, intent(inout) :: v
       real*8, intent(in), optional :: lo, hi
     end subroutine x_r8
     module subroutine x_r8v(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       real*8, intent(inout) :: v(:)
     end subroutine x_r8v
     module subroutine x_r8m(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       real*8, intent(inout) :: v(:,:)
     end subroutine x_r8m
     module subroutine x_r4(io,p,key,v,lo,hi)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       real(c_float), intent(inout) :: v
       real*8, intent(in), optional :: lo, hi
     end subroutine x_r4
     module subroutine x_r4v(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       real(c_float), intent(inout) :: v(:)
     end subroutine x_r4v
     module subroutine x_r4m(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       real(c_float), intent(inout) :: v(:,:)
     end subroutine x_r4m
     module subroutine x_i4(io,p,key,v,lo,hi)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       integer, intent(inout) :: v
       integer, intent(in), optional :: lo, hi
     end subroutine x_i4
     module subroutine x_i4v(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       integer, intent(inout) :: v(:)
     end subroutine x_i4v
     module subroutine x_i4m(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       integer, intent(inout) :: v(:,:)
     end subroutine x_i4m
     module subroutine x_i8v(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       integer*8, intent(inout) :: v(:)
     end subroutine x_i8v
     module subroutine x_l(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       logical, intent(inout) :: v
     end subroutine x_l
     module subroutine x_lv(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       logical, intent(inout) :: v(:)
     end subroutine x_lv
     module subroutine x_lm(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       logical, intent(inout) :: v(:,:)
     end subroutine x_lm
     module subroutine x_s(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       character*(*), intent(inout) :: v
     end subroutine x_s
     module subroutine x_sv(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       character*(*), intent(inout) :: v(:)
     end subroutine x_sv
     module subroutine x_salloc(io,p,key,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       character(len=:), allocatable, intent(inout) :: v
     end subroutine x_salloc
     module subroutine x_r8m_alloc(io,p,key,n1,v)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       integer, intent(in) :: n1
       real*8, allocatable, intent(inout) :: v(:,:)
     end subroutine x_r8m_alloc
     module subroutine xfer_obj(io,p,key,q,found)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       type(json_value), pointer, intent(out) :: q
       logical, intent(out) :: found
     end subroutine xfer_obj
     module subroutine xfer_list(io,p,key,a,n)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: p
       character*(*), intent(in) :: key
       type(json_value), pointer, intent(out) :: a
       integer, intent(inout) :: n
     end subroutine xfer_list
     module subroutine xfer_item(io,a,i,q)
       class(json_xfer), intent(inout) :: io
       type(json_value), pointer, intent(in) :: a
       integer, intent(in) :: i
       type(json_value), pointer, intent(out) :: q
     end subroutine xfer_item
  end interface

end module jsonxfer
