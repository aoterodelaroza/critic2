! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

  !> Build a fragment from atomic data. nspc/spc give the atomic
  !> species. nat is the number of atoms, with coordinates xyz(3,nat)
  !> given in the system indicated by icrd (icrd_cart for Cartesian or
  !> icrd_crys for crystallographic); the complementary set of
  !> coordinates is computed using the crystallographic -> Cartesian
  !> matrix m_x2c. is(nat) are the species and idx(nat)/cidx(nat) the
  !> non-equivalent/complete atom-list indices. Optionally, lvec(3,nat)
  !> are the lattice vectors to each atom in the complete atom list
  !> (default: zero); discrete flags whether the fragment is discrete
  !> (default: .true.); and, for non-discrete fragments, nlvec/mollvec
  !> give the number and list of connecting lattice vectors.
  module subroutine fragment_build(fr,nspc,spc,nat,xyz,icrd,is,idx,cidx,m_x2c,&
     lvec,discrete,nlvec,mollvec)
    use tools_math, only: matinv
    use tools_io, only: ferror, faterr
    use param, only: icrd_cart, icrd_crys
    class(fragment), intent(inout) :: fr
    integer, intent(in) :: nspc
    type(species), intent(in) :: spc(nspc)
    integer, intent(in) :: nat
    real*8, intent(in) :: xyz(3,nat)
    integer, intent(in) :: icrd
    integer, intent(in) :: is(nat)
    integer, intent(in) :: idx(nat)
    integer, intent(in) :: cidx(nat)
    real*8, intent(in) :: m_x2c(3,3)
    integer, intent(in), optional :: lvec(3,nat)
    logical, intent(in), optional :: discrete
    integer, intent(in), optional :: nlvec
    integer, intent(in), optional :: mollvec(:,:)

    integer :: i
    real*8 :: m_c2x(3,3)

    ! crystallographic -> Cartesian matrix
    fr%m_x2c = m_x2c
    if (icrd == icrd_cart) then
       m_c2x = m_x2c
       call matinv(m_c2x,3)
    elseif (icrd /= icrd_crys) then
       call ferror("fragment_build","unsupported coordinate system",faterr)
    end if

    ! atomic species
    fr%nspc = nspc
    if (allocated(fr%spc)) deallocate(fr%spc)
    allocate(fr%spc(nspc))
    fr%spc = spc(1:nspc)

    ! atoms
    fr%nat = nat
    if (allocated(fr%at)) deallocate(fr%at)
    allocate(fr%at(nat))
    do i = 1, nat
       if (icrd == icrd_cart) then
          fr%at(i)%r = xyz(:,i)
          fr%at(i)%x = matmul(m_c2x,xyz(:,i))
       else
          fr%at(i)%x = xyz(:,i)
          fr%at(i)%r = matmul(m_x2c,xyz(:,i))
       end if
       fr%at(i)%is = is(i)
       fr%at(i)%idx = idx(i)
       fr%at(i)%cidx = cidx(i)
       if (present(lvec)) then
          fr%at(i)%lvec = lvec(:,i)
       else
          fr%at(i)%lvec = 0
       end if
    end do

    ! discreteness and connecting lattice vectors
    fr%discrete = .true.
    if (present(discrete)) fr%discrete = discrete
    if (allocated(fr%lvec)) deallocate(fr%lvec)
    fr%nlvec = 0
    if (.not.fr%discrete) then
       if (present(nlvec)) fr%nlvec = nlvec
       if (fr%nlvec > 0) then
          if (.not.present(mollvec)) &
             call ferror("fragment_build","non-discrete fragment without connecting lattice vectors",faterr)
          allocate(fr%lvec(3,fr%nlvec))
          fr%lvec = mollvec(:,1:fr%nlvec)
       end if
    end if

  end subroutine fragment_build

  !> Build a fragment containing a single atom of species spc, located
  !> at coordinates xyz given in the system indicated by icrd
  !> (icrd_cart or icrd_crys). idx/cidx are the non-equivalent/complete
  !> atom-list indices and m_x2c the crystallographic -> Cartesian
  !> matrix. The fragment has a single species, so the atom species
  !> index is 1 and its lattice vector is zero.
  module subroutine fragment_build_atom(fr,spc,xyz,icrd,idx,cidx,m_x2c)
    class(fragment), intent(inout) :: fr
    type(species), intent(in) :: spc
    real*8, intent(in) :: xyz(3)
    integer, intent(in) :: icrd
    integer, intent(in) :: idx
    integer, intent(in) :: cidx
    real*8, intent(in) :: m_x2c(3,3)

    call fr%build(1,(/spc/),1,reshape(xyz,(/3,1/)),icrd,(/1/),(/idx/),(/cidx/),m_x2c)

  end subroutine fragment_build_atom

  !> Translate the fragment by the lattice vector lvec (crystallographic
  !> coordinates). Uses the crystallographic -> Cartesian matrix (m_x2c)
  !> saved at build time.
  module subroutine fragment_translate(fr,lvec)
    class(fragment), intent(inout) :: fr
    integer, intent(in) :: lvec(3)

    integer :: j

    do j = 1, fr%nat
       fr%at(j)%x = fr%at(j)%x + lvec
       fr%at(j)%r = matmul(fr%m_x2c,fr%at(j)%x)
       fr%at(j)%lvec = fr%at(j)%lvec + lvec
    end do

  end subroutine fragment_translate

  !> Translate the fragment so that its center of mass lies in the main
  !> cell (crystallographic coordinates between 0 and 1). Uses the
  !> crystallographic -> Cartesian matrix (m_x2c) saved at build time.
  module subroutine fragment_translate_to_main_cell(fr)
    use tools_math, only: matinv
    class(fragment), intent(inout) :: fr

    integer :: j, newl(3)
    real*8 :: xcm(3), m_c2x(3,3)

    ! Cartesian -> crystallographic matrix
    m_c2x = fr%m_x2c
    call matinv(m_c2x,3)

    ! lattice translation that brings the center of mass into the main cell
    xcm = fr%cmass()
    newl = floor(matmul(m_c2x,xcm))

    ! apply it to all atoms
    do j = 1, fr%nat
       fr%at(j)%x = fr%at(j)%x - newl
       fr%at(j)%r = matmul(fr%m_x2c,fr%at(j)%x)
       fr%at(j)%lvec = fr%at(j)%lvec - newl
    end do

  end subroutine fragment_translate_to_main_cell

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

    if (.not.allocated(fr%at) .or. .not.allocated(fr%spc) .or. fr%nspc == 0 .or. fr%nat == 0) then
       if (allocated(fr%at)) deallocate(fr%at)
       if (allocated(fr%spc)) deallocate(fr%spc)
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
  module subroutine writexyz(fr,file,usenames,ti)
    use tools_io, only: fopen_write, nameguess, fclose
    use param, only: bohrtoa
    class(fragment), intent(in) :: fr
    character*(*), intent(in) :: file
    logical, intent(in) :: usenames
    type(thread_info), intent(in), optional :: ti

    integer :: i, lu

    ! write it
    lu = fopen_write(file,ti=ti)
    write (lu,*) fr%nat
    write (lu,*)
    do i = 1, fr%nat
       if (fr%spc(fr%at(i)%is)%z >= 0) then
          if (usenames) then
             write (lu,'(A,3(F20.10," "))') trim(fr%spc(fr%at(i)%is)%name), fr%at(i)%r * bohrtoa
          else
             write (lu,'(A,3(F20.10," "))') trim(nameguess(fr%spc(fr%at(i)%is)%z,.true.)), fr%at(i)%r * bohrtoa
          end if
       end if
    end do
    call fclose(lu)

  end subroutine writexyz

  !> Write a cml file (molecule) from an array of atomic coordinates.
  module subroutine writecml(fr,file,r,luout,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: pi, bohrtoa
    class(fragment), intent(in) :: fr
    character*(*), intent(in) :: file
    real*8, intent(in), optional :: r(3,3)
    integer, intent(out), optional :: luout
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, lu, iz
    real*8 :: g(3,3), aa(3), bb(3), x(3), ri(3,3)

    ! write it
    lu = fopen_write(file,ti=ti)
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
  module subroutine writegjf(fr,file,ti)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: bohrtoa
    class(fragment), intent(in) :: fr
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: aux
    integer :: i, lu, isum, iz

    aux = file

    ! write it
    lu = fopen_write(aux,ti=ti)
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
