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

! Routines for handling crystal symmetry and related calculations
submodule (crystalmod) symmetry
  implicit none

  !xx! private procedures
  ! subroutine typeop(rot,type,vec,order)
  ! function equiv_tetrah(c,x0,t1,t2,leqv,lrotm,eps)
  ! function perm3(p,r,t) result(res)

  ! symmetry operation symbols
  integer, parameter :: ident=0 !< identifier for sym. operations
  integer, parameter :: inv=1 !< identifier for sym. operations
  integer, parameter :: c2=2 !< identifier for sym. operations
  integer, parameter :: c3=3 !< identifier for sym. operations
  integer, parameter :: c4=4 !< identifier for sym. operations
  integer, parameter :: c6=5 !< identifier for sym. operations
  integer, parameter :: s3=6 !< identifier for sym. operations
  integer, parameter :: s4=7 !< identifier for sym. operations
  integer, parameter :: s6=8 !< identifier for sym. operations
  integer, parameter :: sigma=9 !< identifier for sym. operations

contains

  !> Determines the site symmetry for a point x0 in cryst. coords.
  !> Two points are the same if their distance is less than eps0.
  !> Returns the site symmetry group symbol (sitesymm), the
  !> number of operations in this group (leqv) and the rotation
  !> operations (lrotm)
  module function sitesymm(c,x0,eps0,leqv,lrotm)
    use tools_io, only: string
    use param, only: eye
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in) :: x0(3) !< Input point in cryst. coords.
    real*8, intent(in), optional :: eps0 !< Two points are different if distance is > eps
    character*3 :: sitesymm !< point group symbol
    integer, optional :: leqv !< Number of operations in the group
    real*8, optional :: lrotm(3,3,48) !< Point group operations

    integer :: i, m
    real*8 :: dumy(3), dist2, eps, vec(3)
    integer :: type
    logical :: ok
    integer :: highest, highests
    integer :: nnsym, order
    integer :: ordersym(c%neqv), masksym(0:9)

    real*8, parameter :: eps_default = 1d-2

    if (c%ismolecule .or. c%havesym == 0 .or. (c%neqv == 1 .and. c%ncv == 1)) then
       sitesymm='C1'
       if (present(leqv).and.present(lrotm)) then
          leqv = 1
          lrotm(:,:,1) = eye
       end if
       return
    end if

    if (present(eps0)) then
       eps = eps0
    else
       eps = eps_default
    end if

    ! Run over all proper symmetry elements of symmetry
    nnsym = 0
    masksym = 0
    ordersym = 0
    do i = 1, c%neqv
       ok = .false.
       do m = 1, c%ncv
          dumy = matmul(c%rotm(1:3,1:3,i),x0)
          dumy = - dumy + x0 - c%rotm(1:3,4,i) - c%cen(:,m)
          call c%shortest(dumy,dist2)
          ok = (dist2 < eps)
          if (ok) exit
       end do
       if (ok) then
          ! A symmetry operation at location x has been found.
          nnsym = nnsym+1
          call typeop(c%rotm(:,:,i),type,vec,order)
          ordersym(nnsym) = order
          masksym(type) = masksym(type) + 1
          if (present(leqv).and.present(lrotm)) then
             leqv = nnsym
             lrotm(:,:,leqv) = c%rotm(1:3,1:3,i)
          end if
       endif
    enddo

    ! calculate the point group
    sitesymm = ""
    if (masksym(c3) > 2) then
       ! cubic groups
       if (masksym(c4) /= 0) then
          if (masksym(inv) /= 0) then
             sitesymm='Oh'
          else
             sitesymm='O'
          endif
       elseif (masksym(s4).ne.0) then
          sitesymm= 'Td'
       elseif (masksym(inv).ne.0) then
          sitesymm = 'Th'
       else
          sitesymm = 'T'
       endif
    else
       !Compute highest order proper axis.
       highest=0
       highests=0
       if (masksym(c2) /= 0) highest=2
       if (masksym(c3) /= 0) highest=3
       if (masksym(c4) /= 0) highest=4
       if (masksym(c6) /= 0) highest=6
       if (masksym(s3) /= 0) highests=3
       if (masksym(s4) /= 0) highests=4
       if (masksym(s6) /= 0) highests=6
       if (highest == 0) then
          if (masksym(inv) /= 0) then
             sitesymm='i'
          elseif (masksym(sigma) /= 0) then
             sitesymm='Cs'
          else
             sitesymm='C1'
          endif
       elseif (masksym(c2) >= highest) then
          if (masksym(sigma) == 0) then
             sitesymm='D' // string(highest)
          elseif (masksym(inv) .eq. 1) then
             if (highest == 3) then
                sitesymm= 'D3d'
             else
                sitesymm= 'D' // string(highest) // 'h'
             endif
          else
             if (highest .eq. 3) then
                sitesymm= 'D3h'
             else
                sitesymm= 'D' // string(highest) // 'd'
             endif
          endif
       elseif (masksym(sigma) == 0) then
          if (highests /= 0) then
             sitesymm= 'S' // string(highests/2)
          else
             sitesymm= 'C' // string(highest)
          endif
       elseif (masksym(sigma) .lt. highest) then
          sitesymm= 'C' // string(highest) // 'h'
       else
          sitesymm= 'C' // string(highest) // 'v'
       endif
    endif

  end function sitesymm

  !> Obtain symmetry equivalent positions of xp0 and write them to
  !> vec. Write the multiplicity of the xp0 position to mmult. xp0 and
  !> vec are in crystallographic coordinates. irotm and icenv contain
  !> the index of the rotation matrix and centering vectors
  !> responsible for the transformation of xp into the corresponding
  !> vector in vec. eps is the minimum distance to consider two points
  !> equivalent (in bohr). vec, irotm, icenv, and eps0 are optional.
  module subroutine symeqv(c,xp0,mmult,vec,irotm,icenv,eps0)
    use types, only: realloc
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in) :: xp0(3) !< input position (crystallographic)
    integer, intent(out) :: mmult !< multiplicity
    real*8, allocatable, intent(inout), optional :: vec(:,:) !< sym-eq positions (crystallographic)
    integer, allocatable, intent(inout), optional :: irotm(:) !< index of the operation
    integer, allocatable, intent(inout), optional :: icenv(:) !< index of the cent. vector
    real*8, intent(in), optional :: eps0 !< Minimum distance to consider two vectors different (bohr)

    real*8 :: avec(3,c%neqv*c%ncv)
    integer :: arotm(c%neqv*c%ncv)
    integer :: acenv(c%neqv*c%ncv)
    integer :: i, j
    integer :: mrot, mrot0
    real*8 :: xp(3)
    real*8 :: eps

    real*8, parameter :: eps_default = 1d-2

    ! the eps
    if (present(eps0)) then
       eps = eps0
    else
       eps = eps_default
    end if

    !.Translate position to main cell
    xp = xp0 - floor(xp0)

    !.Run over symmetry operations and create the list of copies
    mrot0 = 0
    do i = 1, c%neqv
       do  j = 1, c%ncv
          mrot0 = mrot0 + 1
          avec(:,mrot0) = matmul(c%rotm(1:3,1:3,i),xp) + c%rotm(:,4,i) + c%cen(:,j)
          avec(:,mrot0) = avec(:,mrot0) - floor(avec(:,mrot0))
          arotm(mrot0) = i
          acenv(mrot0) = j
       enddo
    enddo

    ! calculate distances and (possibly) write vec
    if (present(vec)) call realloc(vec,3,mrot0)
    if (present(irotm)) call realloc(irotm,mrot0)
    if (present(icenv)) call realloc(icenv,mrot0)
    mrot = 0
    d: do i = 1, mrot0
       do j = 1,i-1
          if (c%are_lclose(avec(:,i),avec(:,j),eps)) cycle d
       end do
       mrot = mrot + 1
       if (present(vec)) vec(:,mrot) = avec(:,i)
       if (present(irotm)) irotm(mrot) = arotm(i)
       if (present(icenv)) icenv(mrot) = acenv(i)
    end do d

    mmult=mrot
    if(present(vec)) call realloc(vec,3,mmult)
    if(present(irotm)) call realloc(irotm,mmult)
    if(present(icenv)) call realloc(icenv,mmult)

  end subroutine symeqv

  !> Calculate the multiplicity of the point x0 (cryst. coord.)
  module function get_mult(c,x0) result (mult)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    integer :: mult

    call c%symeqv(x0,mult)

  end function get_mult

  !> Use the spg library to find information about the space group,
  !> encapsulated in the spg user-defined type (spg).
  !> Input: cell vectors (m_x2c), ncel, atcel(:), at(:).
  !> If usenneq, use nneq and at(:) instead of ncel and atcel(:).
  !> If error, return the error in errmsg. Otherwise, return
  !> a zero-length string.
  module subroutine spglib_wrap(c,spg,usenneq,errmsg,ti)
    use iso_c_binding, only: c_double
    use spglib, only: spg_get_dataset, spg_get_error_message
    use global, only: symprec
    use tools_io, only: string, equal
    use param, only: maxzat0
    use types, only: realloc
    class(crystal), intent(in) :: c
    type(SpglibDataset), intent(inout) :: spg
    logical, intent(in) :: usenneq
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    real(c_double) :: lattice(3,3)
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: typ(:)
    integer :: ntyp, nat
    integer :: i, iz(maxzat0)
    character(len=32) :: error

    ! get the dataset from spglib
    errmsg = ""
    lattice = transpose(c%m_x2c)
    iz = 0
    ntyp = 0
    if (usenneq) then
       nat = c%nneq
       ntyp = c%nspc
       allocate(x(3,c%nneq),typ(c%nneq))
       do i = 1, c%nneq
          x(:,i) = c%at(i)%x
          typ(i) = c%at(i)%is
       end do
    else
       nat = c%ncel
       ntyp = c%nspc
       allocate(x(3,c%ncel),typ(c%ncel))
       do i = 1, c%ncel
          x(:,i) = c%atcel(i)%x
          typ(i) = c%atcel(i)%is
       end do
    end if

    spg = spg_get_dataset(lattice,x,typ,nat,symprec)
    deallocate(x,typ)

    ! check error messages (only if this is not threaded, because
    ! this routine is not thread-safe)
    if (.not.present(ti)) then
       error = trim(spg_get_error_message(c%spg%spglib_error))
       if (.not.equal(error,"no error")) then
          errmsg = string(error)
          return
       end if
    end if

  end subroutine spglib_wrap

  !> Set the Wyckoff positions in crystal c from the information in spg.
  !> Assumes the atoms are in the same order. If spg is not present
  !> or not initialized, then fill with "?".
  module subroutine spgtowyc(c,spg)
    class(crystal), intent(inout) :: c
    type(SpglibDataset), intent(inout), optional :: spg

    integer :: i, idx, ilet
    logical :: doit

    character(len=26), parameter :: wycklet = "abcdefghijklmnopqrstuvwxyz"

    doit = .false.
    if (present(spg)) then
       if (spg%n_atoms > 0) doit = .true.
    end if

    ! fill wyckoff letters
    if (doit) then
       do i = 1, c%ncel
          idx = c%atcel(i)%idx
          if (c%at(idx)%wyc == "?" .and. c%spg%wyckoffs(i) < 26) then
             ilet = c%spg%wyckoffs(i)+1
             c%at(idx)%wyc = wycklet(ilet:ilet)
          end if
       end do
    else
       do i = 1, c%nneq
          c%at(i)%wyc = "?"
       end do
    end if

  end subroutine spgtowyc

  !> Calculate the crystal symmetry operations.
  !> Input: cell vectors (m_x2c), ncel, atcel(:), at(:)
  !> Output: neqv, rotm, spg, at()%mult, ncel, and atcel().
  !> If usenneq, use nneq and at(:) instead of ncel and atcel(:).
  !> If error, errmsg in output has length > 0 and contains the error
  !> message.
  module subroutine calcsym(c,usenneq,errmsg,ti)
    use iso_c_binding, only: c_double
    use spglib, only: spg_get_error_message
    use global, only: symprec
    use tools_io, only: string, equal
    use param, only: eyet, eye
    use types, only: realloc
    class(crystal), intent(inout) :: c
    logical, intent(in) :: usenneq
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer, allocatable :: iidx(:)
    integer :: i, j, k, iat, idx
    logical :: found
    real*8 :: rotm(3,3), x0(3)

    c%spgavail = .false.
    call c%spglib_wrap(c%spg,usenneq,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) return
    c%spgavail = .true.

    ! make a copy of nneq into ncel, if appropriate
    if (usenneq) then
       c%ncel = c%nneq
       if (allocated(c%atcel)) deallocate(c%atcel)
       allocate(c%atcel(c%ncel))
       do i = 1, c%ncel
          c%atcel(i)%x = c%at(i)%x
          c%atcel(i)%r = c%at(i)%r

          c%atcel(i)%rxc = c%x2xr(c%atcel(i)%x)
          c%atcel(i)%rxc = c%atcel(i)%rxc - floor(c%atcel(i)%rxc)
          c%atcel(i)%rxc = c%xr2c(c%atcel(i)%rxc)

          c%atcel(i)%is = c%at(i)%is
          c%atcel(i)%cidx = i
       end do
    end if

    ! re-write nneq list based on the information from spg
    allocate(iidx(c%ncel))
    iidx = 0
    c%nneq = 0
    do i = 1, c%spg%n_atoms
       idx = c%spg%equivalent_atoms(i) + 1
       if (iidx(idx) == 0) then
          c%nneq = c%nneq + 1
          iidx(idx) = c%nneq
          if (c%nneq > size(c%at,1)) call realloc(c%at,2*c%nneq)
          c%at(c%nneq)%x = c%atcel(idx)%x
          c%at(c%nneq)%r = c%atcel(idx)%r
          c%at(c%nneq)%is = c%atcel(idx)%is
          c%at(c%nneq)%mult = 1
          c%at(c%nneq)%wyc = "?"
       else
          c%at(iidx(idx))%mult = c%at(iidx(idx))%mult + 1
       end if
       c%atcel(i)%idx = iidx(idx)
       c%atcel(i)%cidx = i
    end do
    call realloc(c%at,c%nneq)

    ! unpack spglib's output into pure translations and symops
    c%neqv = 1
    c%rotm(:,:,1) = eyet
    c%ncv = 1
    if (.not.allocated(c%cen)) allocate(c%cen(3,1))
    c%cen = 0d0
    do i = 1, c%spg%n_operations
       rotm = transpose(c%spg%rotations(:,:,i))

       ! is this a pure translation?
       if (all(abs(rotm - eye) < symprec)) then
          found = .false.
          do j = 1, c%ncv
             if (all(abs(c%spg%translations(:,i) - c%cen(:,j)) < symprec)) then
                found = .true.
                exit
             end if
          end do
          if (.not.found) then
             c%ncv = c%ncv + 1
             if (c%ncv > size(c%cen,2)) &
                call realloc(c%cen,3,2*c%ncv)
             c%cen(:,c%ncv) = c%spg%translations(:,i)
          end if
          cycle
       end if

       ! Do I have this rotation already?
       found = .false.
       do j = 1, c%neqv
          if (all(abs(rotm - c%rotm(1:3,1:3,j)) < symprec)) then
             found = .true.
          end if
       end do
       if (.not.found) then
          c%neqv = c%neqv + 1
          c%rotm(1:3,1:3,c%neqv) = rotm
          c%rotm(1:3,4,c%neqv) = c%spg%translations(:,i)
       end if
    end do
    call realloc(c%cen,3,c%ncv)
    c%havesym = 1

    ! generate symmetry operation info for the complete atom list
    ! use the spg%equivalent_atoms information
    do k = 1, c%ncel
       idx = c%spg%equivalent_atoms(k) + 1
       iat = iidx(idx)

       found = .false.
       loopi: do i = 1, c%neqv
          do j = 1, c%ncv
             x0 = matmul(c%rotm(1:3,1:3,i),c%at(iat)%x) + c%rotm(:,4,i) + c%cen(:,j)
             if (c%are_lclose(x0,c%atcel(k)%x,2d0*symprec)) then
                c%atcel(k)%ir = i
                c%atcel(k)%ic = j
                c%atcel(k)%lvec = nint(c%atcel(k)%x - x0)
                found = .true.
                exit loopi
             end if
          end do
       end do loopi
       if (.not.found) then
          errmsg = "error building rotation and center information for atom list"
          return
       end if
    end do
    deallocate(iidx)

    ! write down the Wyckoff positions from the spg
    call c%spgtowyc(c%spg)

  end subroutine calcsym

  !> Clear the symmetry information from the crystal and rebuild the
  !> atom list. If cel2neq, copy the atom information from atcel(:) to
  !> at(:). If neq2cel, copy from at(:) to atcel(:). Note that if
  !> neq2cel is used, then the environment and asterisms become
  !> obsolete and need to be recalculated.
  module subroutine clearsym(c,cel2neq,neq2cel)
    use types, only: neqatom, realloc
    use param, only: eyet
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: cel2neq
    logical, intent(in), optional :: neq2cel

    type(neqatom), allocatable :: aux(:)
    integer :: idir, i

    idir = 0
    if (present(cel2neq)) then
       if (cel2neq) idir = 1
    end if
    if (present(neq2cel)) then
       if (neq2cel) idir = -1
    end if

    ! nullify the space group and all the symmetry info
    c%spgavail = .false.
    c%spg%n_atoms = 0
    c%havesym = 0
    c%neqv = 1
    c%rotm = 0d0
    c%rotm(:,:,1) = eyet
    c%ncv = 1
    if (allocated(c%cen)) deallocate(c%cen)
    allocate(c%cen(3,4))
    c%cen = 0d0
    if (allocated(c%at) .and. c%nneq > 0) then
       do i = 1, c%nneq
          c%at(i)%mult = 1
          c%at(i)%wyc = "?"
       end do
    end if

    if (idir == 1) then
       ! convert ncel to nneq
       if (c%nneq /= c%ncel) then
          allocate(aux(c%nneq))
          aux = c%at(1:c%nneq)
          c%nneq = c%ncel
          call realloc(c%at,c%ncel)
          do i = 1, c%ncel
             c%at(i) = aux(c%atcel(i)%idx)
             c%at(i)%x = c%atcel(i)%x
             c%at(i)%is = c%atcel(i)%is
             c%at(i)%mult = 1
             c%at(i)%wyc = "?"
          end do
          deallocate(aux)
       end if
    elseif (idir == -1) then
       ! convert nneq to ncel
       c%ncel = c%nneq
       if (allocated(c%atcel)) deallocate(c%atcel)
       allocate(c%atcel(c%ncel))
       do i = 1, c%ncel
          c%atcel(i)%x = c%at(i)%x
          c%atcel(i)%r = c%at(i)%r

          c%atcel(i)%rxc = c%x2xr(c%atcel(i)%x)
          c%atcel(i)%rxc = c%atcel(i)%rxc - floor(c%atcel(i)%rxc)
          c%atcel(i)%rxc = c%xr2c(c%atcel(i)%rxc)

          c%atcel(i)%idx = i
          c%atcel(i)%cidx = i
          c%atcel(i)%ir = 1
          c%atcel(i)%ic = 1
          c%atcel(i)%lvec = 0
          c%atcel(i)%is = c%at(i)%is
       end do
    end if

  end subroutine clearsym

  ! Check that the space group operations are consistent.
  module subroutine checkgroup(c)
    use tools_math, only: matinv
    use tools_io, only: uout, string
    class(crystal), intent(inout) :: c

    integer :: i, j, i1, j1, ifound(2)
    real*8 :: op(4,4), op1(4,4)
    real*8, parameter :: eye4(4,4) = reshape((/1d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0,1d0/),shape(eye4))

    write (uout,'("* Checking the consistency of the space group")')

    ! list of operations
    write (uout,'("+ List of operations")')
    do i = 1, c%neqv
       do  j = 1, c%ncv
          call assignop(i,j,op)
          write (uout,'("  Operation ",A,"/",A)') string(i), string(j)
          op(1:3,4) = c%x2c(op(1:3,4)) * 0.52917720859d0
          write (uout,'(3("    ",4(A," ")/),"    ",4(A," "))') &
             ((string(op(i1,j1),'f',length=9,decimal=6,justify=3), j1=1,4), i1=1,4)
       enddo
    end do

    ! identity
    call checkexists(eye4,ifound)
    if (all(ifound > 0)) then
       write (uout,'("+ Is there an identity operation? ... yes")')
    else
       write (uout,'("+ Is there an identity operation? ... no")')
       goto 999
    end if

    ! inverse
    write (uout,'("+ Inverse operations")')
    do i = 1, c%neqv
       do  j = 1, c%ncv
          call assignop(i,j,op)
          call matinv(op,4)
          call checkexists(op,ifound)
          if (all(ifound > 0)) then
             write (uout,'("  Op. ",A,"/",A," has inverse: ",A,"/",A)') string(i), string(j),&
                string(ifound(1)), string(ifound(2))
          else
             write (uout,'("  Op. ",A,"/",A," has inverse: not found!")') string(i), string(j)
             goto 999
          end if
       end do
    end do

    ! closure
    write (uout,'("+ Closure")')
    do i = 1, c%neqv
       do  j = 1, c%ncv
          call assignop(i,j,op)
          do i1 = 1, c%neqv
             do  j1 = 1, c%ncv
                call assignop(i1,j1,op1)
                op1 = matmul(op,op1)
                call checkexists(op1,ifound)
                if (all(ifound > 0)) then
                   write (uout,'("  Op. ",A,"/",A," times ",A,"/",A," is in the closure: ",A,"/",A)') &
                      string(i), string(j), string(i1), string(j1), string(ifound(1)), string(ifound(2))
                else
                   write (uout,'("  Op. ",A,"/",A," times ",A,"/",A," is not in the closure")') &
                      string(i), string(j), string(i1), string(j1)
                   goto 999
                end if
             end do
          end do
       end do
    end do

999 continue
    write (uout,*)

  contains
    subroutine assignop(ieqv,icv,op)
      real*8, intent(out) :: op(4,4)
      integer, intent(in) :: ieqv, icv

      op = 0d0
      op(1:3,1:3) = c%rotm(:,1:3,ieqv)
      op(1:3,4) = c%rotm(:,4,ieqv) + c%cen(:,icv)
      op(4,4) = 1d0

    end subroutine assignop
    subroutine checkexists(op,ifound)
      real*8, intent(in) :: op(4,4)
      integer, intent(out) :: ifound(2)

      real*8 :: op1(4,4), x(3)
      integer :: i, j

      ifound = 0
      do i = 1, c%neqv
         do  j = 1, c%ncv
            call assignop(i,j,op1)
            x = op(1:3,4) - op1(1:3,4)
            x = x - nint(x)
            if (all(abs(op(1:3,1:3) - op1(1:3,1:3)) < 1d-4) .and. all(abs(x) < 1d-4)) then
               ifound(1) = i
               ifound(2) = j
               return
            end if
         end do
      end do

    end subroutine checkexists
  end subroutine checkgroup

  !> Calculate the irreducible WS wedge around point xorigin (cryst coords)
  !> and partition it into tetrahedra.
  module subroutine getiws(c,xorigin,ntetrag,tetrag)
    use tools_math, only: mixed
    use types, only: realloc
    class(crystal), intent(in) :: c !< the crystal structure
    real*8, intent(in) :: xorigin(3)
    integer, intent(out), optional :: ntetrag
    real*8, allocatable, intent(inout), optional :: tetrag(:,:,:)

    integer :: leqv
    real*8 :: lrotm(3,3,48)
    logical, allocatable :: active(:)
    character*3 :: pg
    integer :: i, j, n, m
    real*8 :: x0(3), xp1(3), xp2(3), xp3(3), xoriginc(3)

    real*8, parameter :: eps_wspesca = 1d-5 !< Criterion for tetrahedra equivalence
    real*8, parameter :: ws_eps_vol = 1d-5 !< Reject tetrahedra smaller than this.

    ! origin in Cartesian
    xoriginc = c%x2c(xorigin)

    ! local symmetry group
    pg = c%sitesymm(xorigin,leqv=leqv,lrotm=lrotm)

    ! count tetrahedra
    ntetrag = 0
    do i = 1, c%ws_nf
       ntetrag = ntetrag + 2 * c%ws_nside(i)
    end do

    ! build all the tetrahedra
    m = 0
    if (allocated(tetrag)) deallocate(tetrag)
    allocate(tetrag(3,4,ntetrag))
    do i = 1, c%ws_nf
       n = c%ws_nside(i)
       ! calculate middle point of the face
       x0 = 0d0
       do j = 1, n
          x0 = x0 + c%ws_x(:,c%ws_iside(j,i))
       end do
       x0 = c%x2c(x0) / n

       do j = 1, n
          xp1 = c%x2c(c%ws_x(:,c%ws_iside(j,i)))
          xp2 = c%x2c(c%ws_x(:,c%ws_iside(mod(j,n)+1,i)))
          ! tetrah 1
          m = m + 1
          tetrag(:,1,m) = (/ 0d0, 0d0, 0d0 /) + xoriginc
          tetrag(:,2,m) = x0 + xoriginc
          tetrag(:,3,m) = xp1 + xoriginc
          tetrag(:,4,m) = 0.5d0 * (xp1 + xp2) + xoriginc
          ! tetrah 2
          m = m + 1
          tetrag(:,1,m) = (/ 0d0, 0d0, 0d0 /) + xoriginc
          tetrag(:,2,m) = x0 + xoriginc
          tetrag(:,3,m) = xp2 + xoriginc
          tetrag(:,4,m) = 0.5d0 * (xp1 + xp2) + xoriginc
       end do
    end do

    ! check volumes and convert to crystallographic
    allocate(active(ntetrag))
    active = .true.
    do i = 1, ntetrag
       ! calculate volume
       xp1 = tetrag(:,2,i) - tetrag(:,1,i)
       xp2 = tetrag(:,3,i) - tetrag(:,1,i)
       xp3 = tetrag(:,4,i) - tetrag(:,1,i)
       if (abs(mixed(xp1,xp2,xp3)) / 6d0 < ws_eps_vol) then
          active(i) = .false.
          cycle
       end if
       do j = 1, 4
          tetrag(:,j,i) = matmul(c%m_c2x,tetrag(:,j,i))
       end do
    end do

    ! reduce tetrahedra equivalent by symmetry
    do i = 1, ntetrag
       if (.not.active(i)) cycle
       do j = i+1, ntetrag
          if (equiv_tetrah(c,xorigin,tetrag(:,:,i),tetrag(:,:,j),leqv,lrotm,eps_wspesca)) then
             active(j) = .false.
          end if
       end do
    end do

    ! rebuild the tetrahedra list
    in: do i = 1, ntetrag
       if (active(i)) cycle
       do j = i+1, ntetrag
          if (active(j)) then
             tetrag(:,:,i) = tetrag(:,:,j)
             active(j) = .false.
             cycle in
          end if
       end do
       ntetrag = i - 1
       exit
    end do in
    call realloc(tetrag,3,4,ntetrag)

    deallocate(active)

  end subroutine getiws

  !xx! private procedures

  !> For the symmetry operation rot, compute the order and
  !> characteristic direction of the operation. The output is given in
  !> type (the type of operation: ident, c2, c3,...)  and vec (the
  !> direction of the axis or the normal to the plane, 0 for a point
  !> symmetry element). Used in the structure initialization.
  subroutine typeop(rot,type,vec,order)
    use tools_math, only: eig
    use tools_io, only: ferror, faterr
    use param, only: tpi, eye

    real*8, intent(in) :: rot(3,4) !< rotm operation
    integer, intent(out) :: type !< output type
    real*8, dimension(3), intent(out) :: vec !< output eigenvector
    integer, intent(out) :: order

    real*8, parameter :: eps = 1e-3

    integer, dimension(2,2) :: iord
    real*8 :: trace, tone
    real*8, dimension(3) :: eigen, eigeni
    real*8, dimension(3,3) :: mat, mat3
    integer :: nm, nones, nminusones
    integer :: i

    trace = rot(1,1) + rot(2,2) + rot(3,3)
    eigen = 0d0
    eigeni = 0d0
    mat = rot(1:3,1:3)
    iord = 0
    vec = 0d0

    !.Determine type of operation
    order = 1
    if (abs(trace+3) .lt. eps) then
       ! inversion
       order = 2
       type = inv
       return
    elseif (abs(trace-3) .lt.eps) then
       ! identity
       type = ident
       return
    else
       nm=3
       ! It is an axis or a plane. Diagonalize.
       call eig(mat,3,eigen,eigeni)

       ! Determine actual operation.
       nones = 0
       nminusones = 0
       do i = 1, 3
          if (abs(eigen(i)-1) < eps .and. abs(eigeni(i)) < eps) then
             nones = nones + 1
             iord(1,nones) = i
          elseif (abs(eigen(i)+1) < eps .and. abs(eigeni(i)) < eps) then
             nminusones = nminusones +1
             iord(2,nminusones) = i
          endif
       enddo

       ! What is it?
       if (nones .eq. 1) then
          !.A proper axis. Order?
          tone = min(max((trace-1) / 2,-1d0),1d0)
          order = nint(abs(tpi/(acos(tone))))
          do i = 1,3
             vec(i) = mat(i,iord(1,1))
          enddo
          if (order .eq. 2) then
             type = c2
          elseif (order .eq. 3) then
             type = c3
          elseif (order .eq. 4) then
             type = c4
          elseif (order .eq. 6) then
             type = c6
          else
             call ferror ('typeop','Axis unknown',faterr)
          endif
       elseif (nones .eq. 2) then
          !A plane, Find directions.
          order = 2
          type = sigma
          do i = 1,3
             vec(i) = mat(i,iord(2,1))
          enddo
       elseif (nminusones .eq. 1) then
          !.An improper axes. Order?
          tone = min(max((trace+1) / 2,-1d0),1d0)
          order = nint(abs(tpi/(acos(tone))))
          do i = 1,3
             vec(i) = mat(i,iord(2,1))
          enddo
          if (order .eq. 3) then
             type = s3
          elseif (order .eq. 4) then
             type = s4
          elseif (order .eq. 6) then
             mat3 = rot(1:3,1:3)
             mat3 = matmul(matmul(mat3,mat3),mat3)
             if (all(abs(mat3+eye) < 1d-5)) then
                type = s3
             else
                type = s6
             end if
          else
             call ferror ('typeop','Axis unknown',faterr)
          endif
       else
          call ferror ('typeop', 'Sym. Element unknown', faterr)
       endif
    endif
    vec = vec / norm2(vec)

  end subroutine typeop

  !> Private for wigner. Determines if two tetrahedra are equivalent.
  function equiv_tetrah(c,x0,t1,t2,leqv,lrotm,eps)
    logical :: equiv_tetrah
    type(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    real*8, dimension(3,0:3), intent(in) :: t1, t2
    integer, intent(in) :: leqv
    real*8, intent(in) :: lrotm(3,3,48), eps

    integer :: i, k, p
    real*8 :: r1(3,0:3), d1(3,0:3), dist2(0:3), xdum(3)

    equiv_tetrah = .false.

    do i = 1, leqv
       do k = 0, 3
          r1(:,k) = matmul(lrotm(:,1:3,i),t1(:,k)-x0) + x0
       end do

       do p = 1, 6
          d1 = perm3(p,r1,t2)
          do k = 0, 3
             xdum = d1(:,k)
             d1(:,k) = c%x2c(xdum)
             dist2(k) = norm2(d1(:,k))
          end do
          if (all(dist2 < eps)) then
             equiv_tetrah = .true.
             return
          end if
       end do
    end do

  end function equiv_tetrah

  !> Private for equiv_tetrah, wigner. 3! permutations.
  function perm3(p,r,t) result(res)
    integer, intent(in) :: p
    real*8, intent(in) :: r(3,0:3), t(3,0:3)
    real*8 :: res(3,0:3)

    select case(p)
    case(1)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,1)
       res(:,2) = r(:,2) - t(:,2)
       res(:,3) = r(:,3) - t(:,3)
    case(2)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,1)
       res(:,2) = r(:,2) - t(:,3)
       res(:,3) = r(:,3) - t(:,2)
    case(3)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,2)
       res(:,2) = r(:,2) - t(:,1)
       res(:,3) = r(:,3) - t(:,3)
    case(4)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,2)
       res(:,2) = r(:,2) - t(:,3)
       res(:,3) = r(:,3) - t(:,1)
    case(5)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,3)
       res(:,2) = r(:,2) - t(:,1)
       res(:,3) = r(:,3) - t(:,2)
    case(6)
       res(:,0) = r(:,0) - t(:,0)
       res(:,1) = r(:,1) - t(:,3)
       res(:,2) = r(:,2) - t(:,2)
       res(:,3) = r(:,3) - t(:,1)
    end select

  end function perm3

  !> Calculate the number of lattice vectors in each direction in
  !> order to be sure that the main cell is surrounded by a shell at
  !> least rmax thick. x2r is the cryst-to-car matrix for the
  !> lattice. The heuristic method has been adapted from gulp.
  module subroutine search_lattice(x2r,rmax,imax,jmax,kmax)
    use param, only: pi

    real*8, intent(in) :: x2r(3,3), rmax
    integer, intent(out) :: imax, jmax, kmax

    integer :: i, nadd
    real*8 :: acell(3), bcell(3), gmat(3,3)

    gmat = matmul(transpose(x2r),x2r)
    do i = 1, 3
       acell(i) = sqrt(gmat(i,i))
    end do
    bcell(1) = acos(gmat(2,3) / acell(2) / acell(3)) * 180d0 / pi
    bcell(2) = acos(gmat(1,3) / acell(1) / acell(3)) * 180d0 / pi
    bcell(3) = acos(gmat(1,2) / acell(1) / acell(2)) * 180d0 / pi

    ! determine cell list, trick adapted from gulp
    nadd = 0
    if (bcell(1)<30 .or. bcell(2)<30 .or. bcell(3)<30 .or. bcell(1)>150 .or. bcell(2)>150 .or. bcell(3)>150) then
       nadd = 3
    else if (bcell(1)<50 .or. bcell(2)<50 .or. bcell(3)<50 .or. bcell(1)>130 .or. bcell(2)>130 .or. bcell(3)>130) then
       nadd = 2
    else if (bcell(1)<70 .or. bcell(2)<70 .or. bcell(3)<70 .or. bcell(1)>110 .or. bcell(2)>110 .or. bcell(3)>110) then
       nadd = 1
    else
       nadd = 1
    end if
    imax = ceiling(rmax / acell(1)) + nadd
    jmax = ceiling(rmax / acell(2)) + nadd
    kmax = ceiling(rmax / acell(3)) + nadd

  end subroutine search_lattice

  !> Get the holohedry and the Laue class from the Hermann-Mauguin
  !> point group label. Adapted from spglib, takes spglib HM point
  !> group labels.
  module subroutine pointgroup_info(hmpg,schpg,holo,laue)
    use tools_io, only: equal
    character*(*), intent(in) :: hmpg
    character(len=3), intent(out) :: schpg
    integer, intent(out) :: holo
    integer, intent(out) :: laue

    if (equal(hmpg,"")) then
       schpg = ""
       holo = holo_unk
       laue = laue_unk
    elseif (equal(hmpg,"1")) then
       schpg = "C1"
       holo = holo_tric
       laue = laue_1
    elseif (equal(hmpg,"-1")) then
       schpg = "Ci"
       holo = holo_tric
       laue = laue_1
    elseif (equal(hmpg,"2")) then
       schpg = "C2"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"m")) then
       schpg = "Cs"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"2/m")) then
       schpg = "C2h"
       holo = holo_mono
       laue = laue_2m
    elseif (equal(hmpg,"222")) then
       schpg = "D2"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"mm2")) then
       schpg = "C2v"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"mmm")) then
       schpg = "D2h"
       holo = holo_ortho
       laue = laue_mmm
    elseif (equal(hmpg,"4")) then
       schpg = "C4"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"-4")) then
       schpg = "S4"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"4/m")) then
       schpg = "C4h"
       holo = holo_tetra
       laue = laue_4m
    elseif (equal(hmpg,"422")) then
       schpg = "D4"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"4mm")) then
       schpg = "C4v"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"-42m")) then
       schpg = "D2d"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"4/mmm")) then
       schpg = "D4h"
       holo = holo_tetra
       laue = laue_4mmm
    elseif (equal(hmpg,"3")) then
       schpg = "C3"
       holo = holo_trig
       laue = laue_3
    elseif (equal(hmpg,"-3")) then
       schpg = "C3i"
       holo = holo_trig
       laue = laue_3
    elseif (equal(hmpg,"32")) then
       schpg = "D3"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"3m")) then
       schpg = "C3v"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"-3m")) then
       schpg = "D3d"
       holo = holo_trig
       laue = laue_3m
    elseif (equal(hmpg,"6")) then
       schpg = "C6"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"-6")) then
       schpg = "C3h"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"6/m")) then
       schpg = "C6h"
       holo = holo_hex
       laue = laue_6m
    elseif (equal(hmpg,"622")) then
       schpg = "D6"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"6mm")) then
       schpg = "C6v"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"-6m2")) then
       schpg = "D3h"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"6/mmm")) then
       schpg = "D6h"
       holo = holo_hex
       laue = laue_6mmm
    elseif (equal(hmpg,"23")) then
       schpg = "T"
       holo = holo_cub
       laue = laue_m3
    elseif (equal(hmpg,"m-3")) then
       schpg = "Th"
       holo = holo_cub
       laue = laue_m3
    elseif (equal(hmpg,"432")) then
       schpg = "O"
       holo = holo_cub
       laue = laue_m3m
    elseif (equal(hmpg,"-43m")) then
       schpg = "Td"
       holo = holo_cub
       laue = laue_m3m
    elseif (equal(hmpg,"m-3m")) then
       schpg = "Oh"
       holo = holo_cub
       laue = laue_m3m
    end if

  end subroutine pointgroup_info

end submodule symmetry
