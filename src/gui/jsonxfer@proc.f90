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


submodule (jsonxfer) proc
  implicit none

contains

  module subroutine x_r8(io,p,key,v,lo,hi)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    real*8, intent(inout) :: v
    real*8, intent(in), optional :: lo, hi

    real*8 :: aux
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (present(lo)) then
          if (aux < lo) return
       end if
       if (present(hi)) then
          if (aux > hi) return
       end if
       v = aux
    else
       call io%json%add(p,key,v)
    end if

  end subroutine x_r8

  module subroutine x_r8v(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    real*8, intent(inout) :: v(:)

    real*8, allocatable :: aux(:)
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (size(aux) == size(v)) v = aux
    else
       call io%json%add(p,key,v)
    end if

  end subroutine x_r8v

  module subroutine x_r8m(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    real*8, intent(inout) :: v(:,:)

    real*8, allocatable :: aux(:)
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (size(aux) == size(v)) v = reshape(aux,shape(v))
    else
       call io%json%add(p,key,reshape(v,(/size(v)/)))
    end if

  end subroutine x_r8m

  ! c_float values go through real*8, and are taken only if they fit
  module subroutine x_r4(io,p,key,v,lo,hi)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    real(c_float), intent(inout) :: v
    real*8, intent(in), optional :: lo, hi

    real*8 :: aux

    aux = real(v,8)
    call io%x_r8(p,key,aux,lo,hi)
    if (io%reading .and. abs(aux) <= huge(1._c_float)) v = real(aux,c_float)

  end subroutine x_r4

  module subroutine x_r4v(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    real(c_float), intent(inout) :: v(:)

    real*8 :: aux(size(v))

    aux = real(v,8)
    call io%x_r8v(p,key,aux)
    if (io%reading .and. all(abs(aux) <= huge(1._c_float))) v = real(aux,c_float)

  end subroutine x_r4v

  module subroutine x_r4m(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    real(c_float), intent(inout) :: v(:,:)

    real*8 :: aux(size(v,1),size(v,2))

    aux = real(v,8)
    call io%x_r8m(p,key,aux)
    if (io%reading .and. all(abs(aux) <= huge(1._c_float))) v = real(aux,c_float)

  end subroutine x_r4m

  module subroutine x_i4(io,p,key,v,lo,hi)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    integer, intent(inout) :: v
    integer, intent(in), optional :: lo, hi

    integer :: aux
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (present(lo)) then
          if (aux < lo) return
       end if
       if (present(hi)) then
          if (aux > hi) return
       end if
       v = aux
    else
       call io%json%add(p,key,v)
    end if

  end subroutine x_i4

  module subroutine x_i4v(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    integer, intent(inout) :: v(:)

    integer, allocatable :: aux(:)
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (size(aux) == size(v)) v = aux
    else
       call io%json%add(p,key,v)
    end if

  end subroutine x_i4v

  module subroutine x_i4m(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    integer, intent(inout) :: v(:,:)

    integer, allocatable :: aux(:)
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (size(aux) == size(v)) v = reshape(aux,shape(v))
    else
       call io%json%add(p,key,reshape(v,(/size(v)/)))
    end if

  end subroutine x_i4m

  ! integer*8 values go through reals (exact below 2**53)
  module subroutine x_i8v(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    integer*8, intent(inout) :: v(:)

    real*8 :: aux(size(v))

    aux = real(v,8)
    call io%x_r8v(p,key,aux)
    if (io%reading .and. all(abs(aux) < 2d0**53)) v = nint(aux,8)

  end subroutine x_i8v

  module subroutine x_l(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    logical, intent(inout) :: v

    logical :: aux, found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (found) v = aux
    else
       call io%json%add(p,key,v)
    end if

  end subroutine x_l

  module subroutine x_lv(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    logical, intent(inout) :: v(:)

    logical, allocatable :: aux(:)
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (size(aux) == size(v)) v = aux
    else
       call io%json%add(p,key,v)
    end if

  end subroutine x_lv

  module subroutine x_lm(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    logical, intent(inout) :: v(:,:)

    logical, allocatable :: aux(:)
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (size(aux) == size(v)) v = reshape(aux,shape(v))
    else
       call io%json%add(p,key,reshape(v,(/size(v)/)))
    end if

  end subroutine x_lm

  module subroutine x_s(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    character*(*), intent(inout) :: v

    character(len=:), allocatable :: aux
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (found) v = aux
    else
       call io%json%add(p,key,trim(v))
    end if

  end subroutine x_s

  module subroutine x_sv(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    character*(*), intent(inout) :: v(:)

    type(json_value), pointer :: a, c
    character(len=:), allocatable :: aux
    character(len=len(v)) :: vnew(size(v))
    integer :: i
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,a,found)
       if (.not.found) return
       if (io%json%count(a) /= size(v)) return
       do i = 1, size(v)
          call io%json%get_child(a,i,c)
          if (.not.associated(c)) return
          call io%json%get(c,aux)
          if (io%json%failed()) then
             call io%json%clear_exceptions()
             return
          end if
          vnew(i) = aux
       end do
       v = vnew
    else
       call io%json%add(p,key,v,trim_str=.true.)
    end if

  end subroutine x_sv

  module subroutine x_salloc(io,p,key,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    character(len=:), allocatable, intent(inout) :: v

    character(len=:), allocatable :: aux
    logical :: found

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (found) v = aux
    else
       if (allocated(v)) call io%json%add(p,key,v)
    end if

  end subroutine x_salloc

  !> An allocatable matrix with n1 rows; on reading, the number of
  !> columns comes from the file (at most xfer_maxlist).
  module subroutine x_r8m_alloc(io,p,key,n1,v)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    integer, intent(in) :: n1
    real*8, allocatable, intent(inout) :: v(:,:)

    real*8, allocatable :: aux(:)
    logical :: found
    integer :: n2

    if (io%reading) then
       call io%json%get(p,key,aux,found)
       if (.not.found) return
       if (mod(size(aux),n1) /= 0) return
       n2 = size(aux) / n1
       if (n2 > xfer_maxlist) return
       if (allocated(v)) deallocate(v)
       allocate(v(n1,n2))
       v = reshape(aux,(/n1,n2/))
    else
       if (allocated(v)) call io%json%add(p,key,reshape(v,(/size(v)/)))
    end if

  end subroutine x_r8m_alloc

  !> A child object of p: made on writing, looked for on reading
  !> (found is .false. and q null if it is not there).
  module subroutine xfer_obj(io,p,key,q,found)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    type(json_value), pointer, intent(out) :: q
    logical, intent(out) :: found

    nullify(q)
    if (io%reading) then
       call io%json%get(p,key,q,found)
       if (.not.found) nullify(q)
    else
       call io%json%create_object(q,key)
       call io%json%add(p,q)
       found = .true.
    end if

  end subroutine xfer_obj

  !> An array of child objects of p: made on writing (n unchanged),
  !> looked for on reading (n is its number of elements, at most
  !> xfer_maxlist, 0 if it is not there).
  module subroutine xfer_list(io,p,key,a,n)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: key
    type(json_value), pointer, intent(out) :: a
    integer, intent(inout) :: n

    logical :: found

    nullify(a)
    if (io%reading) then
       n = 0
       call io%json%get(p,key,a,found)
       if (found) then
          n = min(io%json%count(a),xfer_maxlist)
       else
          nullify(a)
       end if
    else
       call io%json%create_array(a,key)
       call io%json%add(p,a)
    end if

  end subroutine xfer_list

  !> Element i of an array of objects a (xfer_list): appended on
  !> writing, looked for on reading (null if it is not there).
  module subroutine xfer_item(io,a,i,q)
    class(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: a
    integer, intent(in) :: i
    type(json_value), pointer, intent(out) :: q

    nullify(q)
    if (io%reading) then
       call io%json%get_child(a,i,q)
    else
       call io%json%create_object(q,'')
       call io%json%add(a,q)
    end if

  end subroutine xfer_item

end submodule proc
