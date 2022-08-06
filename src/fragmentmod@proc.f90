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

! Molecular fragment class.
submodule (fragmentmod) proc
  implicit none

contains

  !> Initialize a fragment
  module subroutine fragment_init(fr)
    class(fragment), intent(inout) :: fr

    if (allocated(fr%at)) deallocate(fr%at)
    if (allocated(fr%spc)) deallocate(fr%spc)
    allocate(fr%at(1))
    fr%nat = 0
    allocate(fr%spc(1))
    fr%nspc = 0

  end subroutine fragment_init

  !> Merge two or more fragments, delete repeated atoms. If fr already
  !> has a fragment, then add to it if add = .true. (default: .true.).
  !> Assumes all fragments have the same atomic species.
  module subroutine merge_array(fr,fra,add)
    use tools_io, only: equal, ferror, faterr
    use types, only: realloc
    class(fragment), intent(inout) :: fr
    type(fragment), intent(in) :: fra(:)
    logical, intent(in), optional :: add

    real*8, parameter :: eps = 1d-10

    integer :: i, j, k, nat0, nat1
    real*8 :: x(3)
    logical :: found, add0, ok

    add0 = .true.
    if (present(add)) add0 = add

    nat0 = 0
    do i = 1, size(fra)
       nat0 = nat0 + fra(i)%nat
    end do
    if (.not.add) then
       if (allocated(fr%at)) deallocate(fr%at)
       allocate(fr%at(nat0))
       fr%nat = 0
       if (allocated(fr%spc)) deallocate(fr%spc)
       allocate(fr%spc(fra(1)%nspc))
       fr%nspc = fra(1)%nspc
       fr%spc = fra(1)%spc
    end if

    do i = 1, size(fra)
       ok = (fra(i)%nspc == fr%nspc)
       if (ok) then
          do j = 1, fr%nspc
             ok = ok .and. equal(fr%spc(j)%name,fra(i)%spc(j)%name) .and. (fr%spc(j)%z == fra(i)%spc(j)%z)
          end do
       end if
       if (.not.ok) &
          call ferror("merge_array","inconsistent atomic species",faterr)

       nat0 = fr%nat
       nat1 = fr%nat + fra(i)%nat
       if (nat1 > size(fr%at)) call realloc(fr%at,2*nat1)
       do j = 1, fra(i)%nat
          found = .false.
          do k = 1, nat0
             x = abs(fra(i)%at(j)%r - fr%at(k)%r)
             found = all(x < eps)
             if (found) exit
          end do
          if (.not.found) then
             nat0 = nat0 + 1
             fr%at(nat0) = fra(i)%at(j)
          end if
       end do
       fr%nat = nat0
    end do
    call realloc(fr%at,fr%nat)

  end subroutine merge_array

  !> Append a fragment to the current fragment, delete repeated atoms.
  module subroutine append(fr,fra)
    use tools_io, only: ferror, faterr, equal
    use types, only: realloc
    class(fragment), intent(inout) :: fr
    class(fragment), intent(in) :: fra

    real*8, parameter :: eps = 1d-10

    integer :: j, k, nat0, nat1
    real*8 :: x(3)
    logical :: found, ok

    if (.not.allocated(fr%at)) then
       allocate(fr%at(fra%nat))
       allocate(fr%spc(fra%nspc))
       fr%nspc = fra%nspc
       fr%spc = fra%spc
    else
       call realloc(fr%at,fr%nat+fra%nat)
    end if

    ok = (fra%nspc == fr%nspc)
    if (ok) then
       do j = 1, fr%nspc
          ok = ok .and. equal(fr%spc(j)%name,fra%spc(j)%name) .and. (fr%spc(j)%z == fra%spc(j)%z)
       end do
    end if
    if (.not.ok) &
       call ferror("append","inconsistent atomic species",faterr)

    nat0 = fr%nat
    nat1 = fr%nat + fra%nat
    if (nat1 > size(fr%at)) call realloc(fr%at,2*nat1)
    do j = 1, fra%nat
       found = .false.
       do k = 1, nat0
          x = abs(fra%at(j)%r - fr%at(k)%r)
          found = all(x < eps)
          if (found) exit
       end do
       if (.not.found) then
          nat0 = nat0 + 1
          fr%at(nat0) = fra%at(j)
       end if
    end do
    fr%nat = nat0
    call realloc(fr%at,fr%nat)

  end subroutine append

  !> Returns the center of mass (in Cartesian coordinates).  If
  !> weight0 is false, then all atoms have the same weight.
  module function cmass(fr,weight0) result (x)
    use param, only: atmass
    class(fragment), intent(in) :: fr
    logical, intent(in), optional :: weight0
    real*8 :: x(3)

    integer :: i
    logical :: weight
    real*8 :: sum

    weight = .true.
    if (present(weight0)) weight = weight0

    x = 0d0
    sum = 0d0
    if (weight) then
       do i = 1, fr%nat
          if (fr%spc(fr%at(i)%is)%z > 0) then
             x = x + atmass(fr%spc(fr%at(i)%is)%z) * fr%at(i)%r
             sum = sum + atmass(fr%spc(fr%at(i)%is)%z)
          end if
       end do
    else
       do i = 1, fr%nat
          x = x + fr%at(i)%r
          sum = sum + 1d0
       end do
    end if
    x = x / max(sum,1d-40)

  end function cmass

  !> write an xyz-style file from an array of atomic coordinates.
  module subroutine writexyz(fr,file,usenames)
    use tools_io, only: fopen_write, nameguess, fclose
    use param, only: bohrtoa
    class(fragment), intent(in) :: fr
    character*(*), intent(in) :: file
    logical, intent(in) :: usenames

    integer :: i, lu

    ! write it
    lu = fopen_write(file)
    write (lu,*) fr%nat
    write (lu,*)
    do i = 1, fr%nat
       if (fr%spc(fr%at(i)%is)%z >= 0) then
          if (usenames) then
             write (lu,'(A,3(F20.10,X))') trim(fr%spc(fr%at(i)%is)%name), fr%at(i)%r * bohrtoa
          else
             write (lu,'(A,3(F20.10,X))') trim(nameguess(fr%spc(fr%at(i)%is)%z,.true.)), fr%at(i)%r * bohrtoa
          end if
       end if
    end do
    call fclose(lu)

  end subroutine writexyz

  !> Write a cml file (molecule) from an array of atomic coordinates.
  module subroutine writecml(fr,file,r,luout)
    use tools_math, only: matinv
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: pi, bohrtoa
    class(fragment), intent(in) :: fr
    character*(*), intent(in) :: file
    real*8, intent(in), optional :: r(3,3)
    integer, intent(out), optional :: luout

    integer :: i, j, lu, iz
    real*8 :: g(3,3), aa(3), bb(3), x(3), ri(3,3)

    ! write it
    lu = fopen_write(file)
    write (lu,'("<molecule>")')

    ! crystal structure
    if (present(r)) then
       ri = r
       call matinv(ri,3)
       g = matmul(transpose(r),r)
       do i = 1, 3
          aa(i) = sqrt(g(i,i))
       end do
       bb(1) = acos(g(2,3) / aa(2) / aa(3)) * 180d0 / pi
       bb(2) = acos(g(1,3) / aa(1) / aa(3)) * 180d0 / pi
       bb(3) = acos(g(1,2) / aa(1) / aa(2)) * 180d0 / pi
       aa = aa * bohrtoa
       write (lu,'(" <crystal>")')
       write (lu,'("  <scalar title=""a"" units=""units:angstrom"">",A,"</scalar>")') string(aa(1),'f',decimal=8)
       write (lu,'("  <scalar title=""b"" units=""units:angstrom"">",A,"</scalar>")') string(aa(2),'f',decimal=8)
       write (lu,'("  <scalar title=""c"" units=""units:angstrom"">",A,"</scalar>")') string(aa(3),'f',decimal=8)
       write (lu,'("  <scalar title=""alpha"" units=""units:degree"">",A,"</scalar>")') string(bb(1),'f',decimal=4)
       write (lu,'("  <scalar title=""beta"" units=""units:degree"">",A,"</scalar>")') string(bb(2),'f',decimal=4)
       write (lu,'("  <scalar title=""gamma"" units=""units:degree"">",A,"</scalar>")') string(bb(3),'f',decimal=4)
       write (lu,'("  <symmetry spaceGroup=""P 1"">")')
       write (lu,'("   <transform3>1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1</transform3>")')
       write (lu,'("  </symmetry>")')
       write (lu,'(" </crystal>")')
    end if

    write (lu,'(" <atomArray>")')
    do i = 1, fr%nat
       iz = fr%spc(fr%at(i)%is)%z
       if (iz >= 0) then
          if (present(r)) then
             x = matmul(ri,fr%at(i)%r)
             write (lu,'("<atom id=""a",A,""" elementType=""",A,""" xFract=""",A,""" yFract=""",A,""" zFract=""",A,"""/>")') &
                string(i), trim(nameguess(iz,.true.)), (trim(string(x(j),'f',18,10)),j=1,3)
          else
             write (lu,'("<atom id=""a",A,""" elementType=""",A,""" x3=""",A,""" y3=""",A,""" z3=""",A,"""/>")') &
                string(i), trim(nameguess(iz,.true.)), &
                (trim(string(fr%at(i)%r(j) * bohrtoa,'f',18,10)),j=1,3)
          end if
       end if
    end do
    if (present(luout)) then
       luout = lu
    else
       write (lu,'(" </atomArray>")')
       write (lu,'("</molecule>")')
       call fclose(lu)
    end if

  end subroutine writecml

  !> write an Gaussian-style input file from an array of atomic coordinates.
  module subroutine writegjf(fr,file)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: bohrtoa
    class(fragment), intent(in) :: fr
    character*(*), intent(in) :: file

    character(len=:), allocatable :: aux
    integer :: i, lu, isum, iz

    aux = file

    ! write it
    lu = fopen_write(aux)
    write (lu,'("#p b3lyp sto-3g")')
    write (lu,'("")')
    write (lu,'("title")')
    write (lu,'("")')

    isum = 0
    do i = 1, fr%nat
       iz = fr%spc(fr%at(i)%is)%z
       if (iz > 0) isum = isum + iz
    end do
    write (lu,'("0 ",A)') string(2*modulo(isum,2)+1)
    do i = 1, fr%nat
       iz = fr%spc(fr%at(i)%is)%z
       if (iz > 0) then
          write (lu,*) nameguess(iz,.true.), fr%at(i)%r * bohrtoa
       end if
    end do
    write (lu,'("")')
    call fclose(lu)

  end subroutine writegjf

  !> Adapt the size of an allocatable 1D type(fragment) array
  module subroutine realloc_fragment(a,nnew)
    use tools_io, only: ferror, faterr
    type(fragment), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(fragment), allocatable :: temp(:)
    integer :: nold

    if (.not.allocated(a)) &
       call ferror('realloc_fragment','array not allocated',faterr)
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_fragment

end submodule proc
