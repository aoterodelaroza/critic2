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

!> Hash table implementation (type hash)
!> Waiting for unlimited polymorphism (class(*)) to be implemented to generalize this.
module hashmod
  implicit none

  !> Node of a linked list
  type node
     type(node), pointer :: next => NULL()
     character(len=:), allocatable :: key
     character*1, allocatable :: val(:)
   contains
     procedure :: free => free_node
     procedure :: get => get_node
     procedure :: put => put_node
     procedure :: iskey => iskey_node
     procedure :: delkey => delkey_node
     procedure :: assign => assign_node
     !! If I use this, ifort 12.1 freezes when compiling other modules
     !! (like qtree) that has nothing to do with this, so I'm sticking
     !! with node%assign() for now.
     ! generic :: assignment(=) => assign
  end type node

  !> Hash table
  type hash
     type(node), allocatable :: root(:)
     integer :: len = 0
   contains
     procedure :: init => init_hash
     procedure :: free => free_hash
     procedure, private :: put_hash_i_
     procedure, private :: put_hash_r_
     procedure, private :: put_hash_d_
     procedure, private :: put_hash_s_
     generic :: put => put_hash_i_, put_hash_r_, put_hash_d_, put_hash_s_
     procedure, private :: get_hash_i_
     procedure, private :: get_hash_r_
     procedure, private :: get_hash_d_
     procedure, private :: get_hash_s_
     generic :: get => get_hash_i_, get_hash_r_, get_hash_d_, get_hash_s_
     procedure :: type => type_hash
     procedure :: iskey => iskey_hash
     procedure :: delkey => delkey_hash
     procedure :: keys => keys_hash
     procedure :: getkey => getkey_hash
  end type hash

  private
  public :: hash
  private :: init_hash
  private :: free_hash
  private :: get_hash_i_
  private :: get_hash_r_
  private :: get_hash_s_
  private :: get_hash_d_
  private :: type_hash
  private :: put_hash_i_
  private :: put_hash_r_
  private :: put_hash_s_
  private :: put_hash_d_
  private :: hashfun
  private :: free_node
  private :: get_node
  private :: put_node
  private :: iskey_node
  private :: delkey_hash
  private :: delkey_node
  private :: assign_node
  
  interface
     module subroutine init_hash(self,lenh0)
       class(hash), intent(inout) :: self
       integer, intent(in), optional :: lenh0
     end subroutine init_hash
     module subroutine free_hash(self)
       class(hash), intent(inout) :: self
     end subroutine free_hash
     module function get_hash_i_(self,key,typ) result(val)
       class(hash), intent(in) :: self
       character*(*), intent(in) :: key
       integer, intent(in) :: typ
       integer :: val
     end function get_hash_i_
     module function get_hash_r_(self,key,typ) result(val)
       class(hash), intent(in) :: self
       character*(*), intent(in) :: key
       real*4, intent(in) :: typ
       real*4 :: val
     end function get_hash_r_
     module function get_hash_s_(self,key,typ) result(val)
       class(hash), intent(in) :: self
       character*(*), intent(in) :: key
       character*(*), intent(in) :: typ
       character(len=:), allocatable :: val
     end function get_hash_s_
     module function get_hash_d_(self,key,typ) result(val)
       class(hash), intent(in) :: self
       character*(*), intent(in) :: key
       real*8, intent(in) :: typ
       real*8 :: val
     end function get_hash_d_
     module function type_hash(self,key) result(val)
       class(hash), intent(inout) :: self
       character*(*), intent(in) :: key
       character*1, allocatable :: oval(:)
       character*2 :: val
     end function type_hash
     module function iskey_hash(self,key) result(iskey)
       class(hash), intent(in) :: self
       character*(*), intent(in) :: key
       logical :: iskey
     end function iskey_hash
     module subroutine put_hash_i_(self,key,val0)
       class(hash), intent(inout) :: self
       character*(*), intent(in) :: key
       integer, intent(in) :: val0
     end subroutine put_hash_i_
     module subroutine put_hash_r_(self,key,val0)
       class(hash), intent(inout) :: self
       character*(*), intent(in) :: key
       real*4, intent(in) :: val0
     end subroutine put_hash_r_
     module subroutine put_hash_d_(self,key,val0)
       class(hash), intent(inout) :: self
       character*(*), intent(in) :: key
       real*8, intent(in) :: val0
     end subroutine put_hash_d_
     module subroutine put_hash_s_(self,key,val0)
       class(hash), intent(inout) :: self
       character*(*), intent(in) :: key
       character*(*), intent(in) :: val0
     end subroutine put_hash_s_
     module function hashfun(str,lenh)
       character*(*), intent(in) :: str
       integer, intent(in) :: lenh
       integer :: hashfun
     end function hashfun
     module recursive subroutine free_node(self)
       class(node), intent(inout) :: self
     end subroutine free_node
     module recursive subroutine get_node(self,key,val)
       class(node), intent(in) :: self
       character*(*), intent(in) :: key
       character*1, allocatable, intent(inout) :: val(:)
     end subroutine get_node
     module recursive subroutine put_node(self,key,val)
       class(node), intent(inout) :: self
       character*(*), intent(in) :: key
       character*1, intent(in) :: val(:)
     end subroutine put_node
     module recursive function iskey_node(self,key) result(iskey)
       class(node), intent(in) :: self
       character*(*), intent(in) :: key
       logical :: iskey
     end function iskey_node
     module subroutine delkey_hash(self,key) 
       class(hash), intent(inout) :: self
       character*(*), intent(in) :: key
     end subroutine delkey_hash
     module recursive subroutine delkey_node(self,key)
       class(node), intent(inout) :: self
       character*(*), intent(in) :: key
     end subroutine delkey_node
     module subroutine assign_node(to, from)
       class(node), intent(inout) :: to
       type(node), intent(in) :: from
     end subroutine assign_node
     module function keys_hash(self) result (nkeys)
       class(hash), intent(in) :: self
       integer :: nkeys
     end function keys_hash
     module function getkey_hash(self,ikey) result(key)
       class(hash), intent(in) :: self
       integer, intent(in) :: ikey
       character(len=:), allocatable :: key
     end function getkey_hash
  end interface

end module hashmod
