! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (hashmod) proc
  implicit none

contains

  !> Initialize the hash table
  module subroutine init_hash(self,lenh0)
    class(hash), intent(inout) :: self
    integer, intent(in), optional :: lenh0

    integer :: lenh
    integer, parameter :: lenh_default = 40

    if (allocated(self%root)) deallocate(self%root)
    lenh = lenh_default
    if (present(lenh0)) lenh = lenh0
    allocate(self%root(lenh))
    self%len = lenh

  end subroutine init_hash

  !> Destroy the hash table
  module subroutine free_hash(self)
    class(hash), intent(inout) :: self

    integer :: i

    if (self%len == 0) return
    do i = 1, self%len
       call self%root(i)%free()
    end do
    deallocate(self%root)
    self%len = 0

  end subroutine free_hash

  !> Get a value from the hash table, return integer
  module function get_hash_i_(self,key,typ) result(val)
    class(hash), intent(in) :: self
    character*(*), intent(in) :: key
    integer, intent(in) :: typ
    integer :: val

    character*1, allocatable :: str(:)
    integer :: fh

    val = 0d0
    if (self%len == 0) return
    fh = hashfun(key,self%len)
    call self%root(fh)%get(key,str)
    val = transfer(str(3:),val)

  end function get_hash_i_

  !> Get a value from the hash table, return real*4
  module function get_hash_r_(self,key,typ) result(val)
    class(hash), intent(in) :: self
    character*(*), intent(in) :: key
    real*4, intent(in) :: typ
    real*4 :: val

    character*1, allocatable :: str(:)
    integer :: fh

    val = 0d0
    if (self%len == 0) return
    fh = hashfun(key,self%len)
    call self%root(fh)%get(key,str)
    val = transfer(str(3:),val)

  end function get_hash_r_

  !> Get a value from the hash table, return string
  module function get_hash_s_(self,key,typ) result(val)
    class(hash), intent(in) :: self
    character*(*), intent(in) :: key
    character*(*), intent(in) :: typ
    character(len=:), allocatable :: val

    character*1, allocatable :: str(:)
    integer :: fh

    if (self%len == 0) return
    fh = hashfun(key,self%len)
    call self%root(fh)%get(key,str)
    allocate(character(len=size(str)-2)::val)
    val = transfer(str(3:),val)

  end function get_hash_s_

  !> Get a value from the hash table, return real*8
  module function get_hash_d_(self,key,typ) result(val)
    class(hash), intent(in) :: self
    character*(*), intent(in) :: key
    real*8, intent(in) :: typ
    real*8 :: val

    character*1, allocatable :: str(:)
    integer :: fh

    val = 0d0
    if (self%len == 0) return
    fh = hashfun(key,self%len)
    call self%root(fh)%get(key,str)
    val = transfer(str(3:),val)

  end function get_hash_d_

  !> Get the data type for key from the hash table
  module function type_hash(self,key) result(val)
    class(hash), intent(inout) :: self
    character*(*), intent(in) :: key
    character*1, allocatable :: oval(:)
    character*2 :: val

    integer :: fh

    if (self%len == 0) return
    fh = hashfun(key,self%len)
    call self%root(fh)%get(key,oval)
    val = oval(1) // oval(2)

  end function type_hash

  !> Check whether the given key exists in the hash
  module function iskey_hash(self,key) result(iskey)
    class(hash), intent(in) :: self
    character*(*), intent(in) :: key
    logical :: iskey

    integer :: fh

    iskey = .false.
    if (self%len == 0) return
    fh = hashfun(key,self%len)
    iskey = self%root(fh)%iskey(key)

  end function iskey_hash

  !> Put a value in the hash table, integer
  module subroutine put_hash_i_(self,key,val0)
    class(hash), intent(inout) :: self
    character*(*), intent(in) :: key
    integer, intent(in) :: val0

    character*1 :: typ(1)
    character*1, allocatable :: val(:)
    integer :: isiz

    isiz = size(transfer(val0,typ))
    allocate(val(isiz+2))
    val(1:2) = (/"i","_"/)
    val(3:) = transfer(val0,val)

    if (self%len == 0) call init_hash(self)
    call self%root(hashfun(key,self%len))%put(key,val)

  end subroutine put_hash_i_

  !> Put a value in the hash table, real*4
  module subroutine put_hash_r_(self,key,val0)
    class(hash), intent(inout) :: self
    character*(*), intent(in) :: key
    real*4, intent(in) :: val0

    character*1 :: typ(1)
    character*1, allocatable :: val(:)
    integer :: isiz

    isiz = size(transfer(val0,typ))
    allocate(val(isiz+2))
    val(1:2) = (/"r","_"/)
    val(3:) = transfer(val0,val)

    if (self%len == 0) call init_hash(self)
    call self%root(hashfun(key,self%len))%put(key,val)

  end subroutine put_hash_r_

  !> Put a value in the hash table, real*8
  module subroutine put_hash_d_(self,key,val0)
    class(hash), intent(inout) :: self
    character*(*), intent(in) :: key
    real*8, intent(in) :: val0

    character*1 :: typ(1)
    character*1, allocatable :: val(:)
    integer :: isiz

    isiz = size(transfer(val0,typ))
    allocate(val(isiz+2))
    val(1:2) = (/"d","_"/)
    val(3:) = transfer(val0,val)

    if (self%len == 0) call init_hash(self)
    call self%root(hashfun(key,self%len))%put(key,val)

  end subroutine put_hash_d_

  !> Put a value in the hash table, string
  module subroutine put_hash_s_(self,key,val0)
    class(hash), intent(inout) :: self
    character*(*), intent(in) :: key
    character*(*), intent(in) :: val0

    character*1 :: typ(1)
    character*1, allocatable :: val(:)
    integer :: isiz

    isiz = size(transfer(val0,typ))
    allocate(val(isiz+2))
    val(1:2) = (/"s","_"/)
    val(3:) = transfer(val0,val)

    if (self%len == 0) call init_hash(self)
    call self%root(hashfun(key,self%len))%put(key,val)

  end subroutine put_hash_s_

  !> Hashing function
  module function hashfun(str,lenh)
    character*(*), intent(in) :: str
    integer, intent(in) :: lenh
    integer :: hashfun

    integer :: i
    integer*8 :: isum

    isum = 0
    do i = 1, len(str)
       isum = isum + ichar(str(i:i))
    end do
    hashfun = int(mod(isum,int(lenh,8)),kind(hashfun)) + 1

  end function hashfun

  !> Free a node
  recursive module subroutine free_node(self)
    class(node), intent(inout) :: self

    if (associated(self%next)) then
       call free_node(self%next)
       deallocate(self%next)
    end if
    self%next => NULL()
    if (allocated(self%key)) deallocate(self%key)
    if (allocated(self%val)) deallocate(self%val)

  end subroutine free_node

  !> Get a value from a node or its children
  recursive module subroutine get_node(self,key,val)
    class(node), intent(in) :: self
    character*(*), intent(in) :: key
    character*1, allocatable, intent(inout) :: val(:)

    if (allocated(self%key)) then
       if (self%key == key) then
          if (allocated(val)) deallocate(val)
          allocate(val(size(self%val)))
          val = self%val
          return
       endif
    endif
    if (.not.associated(self%next)) return
    call get_node(self%next,key,val)

  end subroutine get_node

  !> Put a value in a node or its children
  recursive module subroutine put_node(self,key,val)
    class(node), intent(inout) :: self
    character*(*), intent(in) :: key
    character*1, intent(in) :: val(:)

    if (allocated(self%key)) then
       if (self%key /= key) then
          if (.not.associated(self%next)) allocate(self%next)
          call put_node(self%next,key,val)
       else
          if (allocated(self%val)) deallocate(self%val)
          allocate(self%val(size(val)))
          self%val = val
       end if
    else
       self%key = key
       allocate(self%val(size(val)))
       self%val = val
    end if

  end subroutine put_node

  !> Check whether this node or any of its children is the key
  recursive module function iskey_node(self,key) result(iskey)
    class(node), intent(in) :: self
    character*(*), intent(in) :: key
    logical :: iskey

    iskey = .false.
    if (allocated(self%key)) then
       if (self%key == key) then
          iskey = .true.
          return
       endif
    endif
    if (.not.associated(self%next)) return
    iskey = iskey_node(self%next,key)

  end function iskey_node

  !> Delete a key/value pair from the hash. No action if the key
  !> does not exist.
  module subroutine delkey_hash(self,key) 
    class(hash), intent(inout) :: self
    character*(*), intent(in) :: key

    integer :: fh

    if (self%len == 0) return
    fh = hashfun(key,self%len)
    call delkey_node(self%root(fh),key)

  end subroutine delkey_hash

  !> Delete the node that matches the key (doesn't work on the first)
  recursive module subroutine delkey_node(self,key)
    class(node), intent(inout) :: self
    character*(*), intent(in) :: key

    if (allocated(self%key)) then
       if (self%key == key) then
          if (associated(self%next)) then
             call self%assign(self%next)
          else
             deallocate(self%key)
             if (allocated(self%val)) deallocate(self%val)
          endif
          return
       endif
    endif
    if (.not.associated(self%next)) return
    call delkey_node(self%next,key)

  end subroutine delkey_node

  !> Assignment between nodes
  module subroutine assign_node(to, from)
    class(node), intent(inout) :: to
    type(node), intent(in) :: from

    if (allocated(to%key)) deallocate(to%key)
    if (allocated(to%val)) deallocate(to%val)
    if (allocated(from%key)) to%key = from%key
    if (allocated(from%val)) then
       allocate(to%val(size(from%val)))
       to%val = from%val
    endif
    to%next => from%next

  end subroutine assign_node

  !> Return the number of keys in the hash
  module function keys_hash(self) result (nkeys)
    class(hash), intent(in) :: self
    integer :: nkeys

    integer :: i
    type(node), pointer :: cur

    nkeys = 0
    do i = 1, self%len
       if (allocated(self%root(i)%key)) nkeys = nkeys + 1
       cur => self%root(i)%next
       do while(associated(cur))
          if (allocated(cur%key)) nkeys = nkeys + 1
          cur => cur%next
       end do
    end do

  end function keys_hash

  !> Return the key corresponding to integer ikey (1...keys_hash())
  module function getkey_hash(self,ikey) result(key)
    class(hash), intent(in) :: self
    integer, intent(in) :: ikey
    character(len=:), allocatable :: key

    integer :: i, nkeys
    type(node), pointer :: cur

    key = ""
    nkeys = 0
    do i = 1, self%len
       ! not having straightforward allocatable arrays of pointers
       ! is a major pain in the butt
       if (allocated(self%root(i)%key)) then
          nkeys = nkeys + 1
          if (nkeys == ikey) then
             key = self%root(i)%key
             return
          end if
       endif
       cur => self%root(i)%next
       do while(associated(cur))
          if (allocated(cur%key)) then
             nkeys = nkeys + 1
             if (nkeys == ikey) then
                key = cur%key
                return
             end if
          endif
          cur => cur%next
       end do
    end do

  end function getkey_hash

end submodule proc
