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
  ! subroutine symelem_enum(c,se,iop)
  ! subroutine symelem_replicate(c,ncell,border,se)
  ! subroutine symelem_in_box(skind,x0,dirf,b,ok,xmid,onfar)
  ! function symelem_dirf(c,skind,dirc)
  ! function glide_letter(tred)
  ! function symelem_findtype(se,skind,sorder,dirc,slabel)
  ! subroutine symelem_addtype(se,skind,sorder,dirc,slabel,sdirlabel,nop,iop)
  ! subroutine symelem_add(se,skind,sorder,dirc,slabel,sdirlabel,xc,ifrom,iop,itype)
  ! subroutine symelem_trim(se)
  ! function symelem_key(skind,dirc,xc)
  ! subroutine primitive_int(v,iv,ok)
  ! function is_lattice_vec(c,v)
  ! subroutine inplane_lattice_vecs(c,ihkl,hklr,nv,vlist)
  ! subroutine reduce_glide(c,nv,vlist,t,tred)
  ! function glide_rank(t)
  ! subroutine typeop(rot,type,vec,order)
  ! function equiv_tetrah(c,x0,t1,t2,leqv,lrotm,eps)
  ! function perm3(p,r,t) result(res)

  ! half-width of the box of lattice translations added to an operation when
  ! enumerating the distinct elements it generates in one cell. The element
  ! moves by a fraction of the translation, never by more than the translation
  ! itself, so reaching one cell away needs no more than a couple of cells of
  ! translation; 4 leaves a wide margin (verified against the 530 structures
  ! of tests/zz_source/cif/allspg).
  integer, parameter :: tlim_def = 4

  ! tolerance for the element clipping, in box coordinates, so a fraction of
  ! the box rather than a length
  real*8, parameter :: clip_eps = 1d-5

  ! most symmetry operations recorded per element type
  integer, parameter :: maxop_def = 16

  integer, parameter :: inplane_mbox = 3
  integer, parameter :: inplane_nbox = (2*inplane_mbox+1)**3

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

  !> Determine the rotations in the point group from the space group
  !> information. The number of rotations is leqv and the rotation
  !> matrices are rotm(3,3,leqv).
  module subroutine pointgroup(c,leqv,rotm)
    use param, only: eye
    class(crystal), intent(in) :: c
    integer, intent(out) :: leqv
    real*8, allocatable, intent(inout) :: rotm(:,:,:)

    integer :: i, j
    logical :: found

    real*8, parameter :: eps = 1d-5

    ! initialize and check
    leqv = 1
    if (allocated(rotm)) deallocate(rotm)
    if (.not.c%spgavail .or. c%neqv == 1) then
       allocate(rotm(3,3,1))
       rotm(:,:,1) = eye
       return
    else
       allocate(rotm(3,3,c%neqv))
       rotm(:,:,1) = eye
    end if

    do i = 1, c%neqv
       found = .false.
       do j = 1, leqv
          if (all(abs(c%rotm(1:3,1:3,i) - rotm(:,:,j)) < eps)) then
             found = .true.
             exit
          end if
       end do
       if (.not.found) then
          leqv = leqv + 1
          rotm(:,:,leqv) = c%rotm(1:3,1:3,i)
       end if
    end do

  end subroutine pointgroup

  !> Determines the site symmetry for a point x0 in cryst. coords.
  !> Two points are the same if their distance is less than eps0.
  !> Returns the site symmetry group symbol (sitesymm), the
  !> number of operations in this group (leqv) and the rotation
  !> operations (lrotm)
  module function sitesymm(c,x0,eps0,leqv,lrotm)
    use param, only: eye
    class(crystal), intent(in) :: c !< Input crystal
    real*8, intent(in) :: x0(3) !< Input point in cryst. coords.
    real*8, intent(in), optional :: eps0 !< Two points are different if distance is > eps
    character*3 :: sitesymm !< point group symbol
    integer, optional :: leqv !< Number of operations in the group
    real*8, optional :: lrotm(3,3,48) !< Point group operations

    integer :: i, m, nnsym
    real*8 :: dumy(3), dist2, eps
    real*8 :: lrot(3,3,c%neqv)
    logical :: ok

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

    ! collect the operations that leave x0 in place
    nnsym = 0
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
          nnsym = nnsym + 1
          lrot(:,:,nnsym) = c%rotm(1:3,1:3,i)
       endif
    enddo
    if (present(leqv).and.present(lrotm)) then
       leqv = nnsym
       lrotm(:,:,1:nnsym) = lrot(:,:,1:nnsym)
    end if

    ! and name the group they form
    sitesymm = pointgroup_symbol(nnsym,lrot(:,:,1:nnsym))

  end function sitesymm

  !> Schoenflies symbol of the point group formed by the nops rotations
  !> rotm(3,3,nops) (the set must be a group). The site symmetry of a point and the symmetry
  !> compatible with a supercell lattice are both named through this.
  module function pointgroup_symbol(nops,rotm) result(symb)
    use tools_io, only: string
    integer, intent(in) :: nops
    real*8, intent(in) :: rotm(3,3,nops)
    character*3 :: symb

    integer :: i, type, order, highest, highests
    integer :: masksym(0:9)
    real*8 :: vec(3), op(3,4)

    ! count the operations of each type (typeop only reads the rotation)
    masksym = 0
    op = 0d0
    do i = 1, nops
       op(:,1:3) = rotm(:,:,i)
       call typeop(op,type,vec,order)
       masksym(type) = masksym(type) + 1
    end do

    ! calculate the point group
    symb = ""
    if (masksym(c3) > 2) then
       ! cubic groups
       if (masksym(c4) /= 0) then
          if (masksym(inv) /= 0) then
             symb='Oh'
          else
             symb='O'
          endif
       elseif (masksym(s4).ne.0) then
          symb= 'Td'
       elseif (masksym(inv).ne.0) then
          symb = 'Th'
       else
          symb = 'T'
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
             symb='i'
          elseif (masksym(sigma) /= 0) then
             symb='Cs'
          else
             symb='C1'
          endif
       elseif (masksym(c2) >= highest) then
          if (masksym(sigma) == 0) then
             symb='D' // string(highest)
          elseif (masksym(inv) .eq. 1) then
             if (highest == 3) then
                symb= 'D3d'
             else
                symb= 'D' // string(highest) // 'h'
             endif
          else
             if (highest .eq. 3) then
                symb= 'D3h'
             else
                symb= 'D' // string(highest) // 'd'
             endif
          endif
       elseif (masksym(sigma) == 0) then
          if (highests /= 0) then
             symb = 'S' // string(highests)
          else
             symb= 'C' // string(highest)
          endif
       elseif (masksym(sigma) .lt. highest) then
          symb= 'C' // string(highest) // 'h'
       else
          symb= 'C' // string(highest) // 'v'
       endif
    endif

  end function pointgroup_symbol

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
  module subroutine spglib_wrap(c,spg,usenneq,errmsg,useidx,ti)
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
    logical, intent(in), optional :: useidx
    type(thread_info), intent(in), optional :: ti

    real(c_double) :: lattice(3,3)
    real(c_double), allocatable :: x(:,:)
    integer, allocatable :: typ(:)
    integer :: ntyp, nat
    integer :: i, iz(maxzat0)
    character(len=32) :: error
    logical :: idxc

    ! color the atoms by non-equivalent-atom index instead of species
    idxc = .false.
    if (present(useidx)) idxc = useidx

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
          if (idxc) then
             typ(i) = c%atcel(i)%idx
          else
             typ(i) = c%atcel(i)%is
          end if
       end do
    end if

    spg = spg_get_dataset(lattice,x,typ,nat,symprec)
    deallocate(x,typ)

    ! check error messages (only if this is not threaded, because
    ! this routine is not thread-safe)
    if (.not.present(ti)) then
       error = trim(spg_get_error_message(spg%spglib_error))
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
    character*10, allocatable :: name(:)
    real*8, allocatable :: occs(:)
    type(siteocc), allocatable :: mixs(:)

    c%spgavail = .false.
    call c%spglib_wrap(c%spg,usenneq,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) return
    if (c%spg%n_atoms == 0 .or. c%spg%n_operations == 0) then
       errmsg = "error detecting space group"
       return
    end if
    c%spgavail = .true.

    ! make a copy of nneq into ncel, if appropriate
    ! carry the mixed-site co-occupant lists through the rebuild, keyed by the
    ! cell atom they come from (they are not stored in the cell list)
    if (usenneq) then
       c%ncel = c%nneq
       allocate(name(c%ncel),occs(c%ncel),mixs(c%ncel))
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
          name(i) = c%at(i)%name
          occs(i) = c%at(i)%occ
          if (allocated(c%at(i)%mix)) mixs(i) = c%at(i)%mix
       end do
    else
       allocate(name(c%ncel),occs(c%ncel),mixs(c%ncel))
       do i = 1, c%ncel
          name(i) = c%at(c%atcel(i)%idx)%name
          occs(i) = c%at(c%atcel(i)%idx)%occ
          if (allocated(c%at(c%atcel(i)%idx)%mix)) mixs(i) = c%at(c%atcel(i)%idx)%mix
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
          c%at(c%nneq)%name = name(idx)
          c%at(c%nneq)%mult = 1
          c%at(c%nneq)%wyc = "?"
          c%at(c%nneq)%occ = occs(idx)
          ! restore (or clear stale) mixed-site co-occupants for this representative
          if (mixs(idx)%nocc >= 2) then
             if (.not.allocated(c%at(c%nneq)%mix)) allocate(c%at(c%nneq)%mix)
             c%at(c%nneq)%mix = mixs(idx)
          else if (allocated(c%at(c%nneq)%mix)) then
             deallocate(c%at(c%nneq)%mix)
          end if
       else
          c%at(iidx(idx))%mult = c%at(iidx(idx))%mult + 1
       end if
       c%atcel(i)%idx = iidx(idx)
       c%atcel(i)%cidx = i
    end do
    call realloc(c%at,c%nneq)
    deallocate(name,occs,mixs)
    call c%set_haveocc()

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

  !> Clear the symmetry information from the crystal and set it to P1.
  !> If cel2neq, also copy the complete cell list atcel(:) over to the
  !> non-equivalent atom list at(:) (ncel -> nneq).
  module subroutine clearsym(c,cel2neq)
    use types, only: neqatom, realloc
    use param, only: eyet
    class(crystal), intent(inout) :: c
    logical, intent(in), optional :: cel2neq

    type(neqatom), allocatable :: aux(:)
    integer :: i, j
    character*10, allocatable :: name(:)
    logical :: cel2neq_

    cel2neq_ = .false.
    if (present(cel2neq)) cel2neq_ = cel2neq

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

    ! optionally convert ncel to nneq
    if (cel2neq_ .and. c%nneq /= c%ncel) then
       allocate(name(c%ncel))
       do i = 1, c%ncel
          name(i) = c%at(c%atcel(i)%idx)%name
       end do

       allocate(aux(c%nneq))
       aux = c%at(1:c%nneq)
       c%nneq = c%ncel
       call realloc(c%at,c%ncel)
       ! every cell atom becomes its own non-equivalent atom
       do i = 1, c%ncel
          c%at(i) = aux(c%atcel(i)%idx)
          c%at(i)%x = c%atcel(i)%x
          c%at(i)%r = c%atcel(i)%r
          c%at(i)%is = c%atcel(i)%is
          c%at(i)%mult = 1
          c%at(i)%wyc = "?"
          c%at(i)%name = name(i)
          c%atcel(i)%idx = i
          c%atcel(i)%cidx = i
          c%atcel(i)%lvec = 0
          c%atcel(i)%ir = 1
          c%atcel(i)%ic = 1
       end do
       deallocate(aux,name)

       ! the molecular fragments carry their own copy of the mapping to the
       ! non-equivalent list
       if (allocated(c%mol) .and. c%nmol > 0) then
          do i = 1, c%nmol
             do j = 1, c%mol(i)%nat
                c%mol(i)%at(j)%idx = c%mol(i)%at(j)%cidx
             end do
          end do
          call c%calculate_molecular_equivalence()
       end if
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

  !> Find the largest subgroup that does not contain any of the
  !> operations flagged in del (the identity, operation 1, is never
  !> removable). Returns lkeep, a mask of the operations of the chosen
  !> subgroup. Ties (several subgroups of the same, maximal size
  !> avoiding the deleted operations) are broken by keeping the most
  !> non-deleted operations, then by lowest operation indices.
  module subroutine symop_subgroup(c,del,lkeep)
    class(crystal), intent(in) :: c
    logical, intent(in) :: del(:)
    logical, intent(out) :: lkeep(:)

    integer, parameter :: nsubmax = 4096
    integer :: i, j, k, n, nsub, ifound(2), bestpop, pop, ov, bestov
    integer, allocatable :: prod(:,:)
    logical, allocatable :: subs(:,:), bestset(:), tmp(:)
    real*8 :: op(4,4), op1(4,4)
    logical :: better

    ! initialize
    n = c%neqv
    if (n <= 1) return
    lkeep = .false.
    lkeep(1:n) = .true.

    ! calculate the multiplication table prod(i,j)
    allocate(prod(n,n))
    do i = 1, n
       op = 0d0
       op(1:3,1:3) = c%rotm(:,1:3,i)
       op(1:3,4) = c%rotm(:,4,i)
       op(4,4) = 1d0
       do j = 1, n
          op1 = 0d0
          op1(1:3,1:3) = c%rotm(:,1:3,j)
          op1(1:3,4) = c%rotm(:,4,j)
          op1(4,4) = 1d0
          call opindex(matmul(op,op1),ifound)
          prod(i,j) = ifound(1)
       end do
    end do
    if (any(prod == 0)) return

    ! enumerate the subgroups as logical sets. Seed with the trivial
    ! subgroup and each cyclic subgroup, then join pairs of subgroups
    ! until no new ones appear.
    allocate(subs(n,nsubmax),bestset(n),tmp(n))
    nsub = 0
    tmp = .false.
    tmp(1) = .true.
    call addset(tmp)
    do i = 2, n
       tmp = .false.
       tmp(1) = .true.
       tmp(i) = .true.
       call addset(closeset(tmp))
    end do
    k = 1
    do while (k <= nsub)
       do j = 1, nsub
          call addset(closeset(subs(:,k) .or. subs(:,j)))
          if (nsub >= nsubmax) exit
       end do
       k = k + 1
    end do

    ! pick the largest subgroup that contains no deleted operation; tie-break by
    ! the most kept (non-deleted) operations, then prefer including the lowest-
    ! index operation
    bestpop = -1
    bestov = -1
    bestset = .false.
    bestset(1) = .true.
    do k = 1, nsub
       if (any(subs(:,k) .and. del(1:n))) cycle
       pop = count(subs(:,k))
       ov = count(subs(:,k) .and. .not.del(1:n))
       better = .false.
       if (pop > bestpop) then
          better = .true.
       elseif (pop == bestpop .and. ov > bestov) then
          better = .true.
       elseif (pop == bestpop .and. ov == bestov) then
          do i = 1, n
             if (subs(i,k) .neqv. bestset(i)) then
                better = subs(i,k)
                exit
             end if
          end do
       end if
       if (better) then
          bestpop = pop
          bestov = ov
          bestset = subs(:,k)
       end if
    end do

    lkeep = .false.
    lkeep(1:n) = bestset(1:n)
    deallocate(prod,subs,bestset,tmp)

  contains
    ! Find the index of the 4x4 operation op
    subroutine opindex(op,ifound)
      real*8, intent(in) :: op(4,4)
      integer, intent(out) :: ifound(2)

      real*8 :: op1(4,4), x(3)
      integer :: ii, jj

      ifound = 0
      do ii = 1, c%neqv
         do jj = 1, c%ncv
            op1 = 0d0
            op1(1:3,1:3) = c%rotm(:,1:3,ii)
            op1(1:3,4) = c%rotm(:,4,ii) + c%cen(:,jj)
            op1(4,4) = 1d0
            x = op(1:3,4) - op1(1:3,4)
            x = x - nint(x)
            if (all(abs(op(1:3,1:3) - op1(1:3,1:3)) < 1d-4) .and. all(abs(x) < 1d-4)) then
               ifound(1) = ii
               ifound(2) = jj
               return
            end if
         end do
      end do

    end subroutine opindex

    ! close a set of operations under composition (result is a subgroup)
    function closeset(set0) result(set)
      logical, intent(in) :: set0(:)
      logical :: set(n)

      integer :: ii, jj, kk
      logical :: changed

      set = set0(1:n)
      changed = .true.
      do while (changed)
         changed = .false.
         do ii = 1, n
            if (.not.set(ii)) cycle
            do jj = 1, n
               if (.not.set(jj)) cycle
               kk = prod(ii,jj)
               if (.not.set(kk)) then
                  set(kk) = .true.
                  changed = .true.
               end if
            end do
         end do
      end do

    end function closeset

    ! add a subgroup set to the list if not already present
    subroutine addset(set)
      logical, intent(in) :: set(:)
      integer :: ii

      do ii = 1, nsub
         if (all(subs(1:n,ii) .eqv. set(1:n))) return
      end do
      if (nsub < nsubmax) then
         nsub = nsub + 1
         subs(1:n,nsub) = set(1:n)
      end if

    end subroutine addset
  end subroutine symop_subgroup

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

  !> Mark in mask the element types of the list se that the symmetry
  !> operations iop generate (see list_symelems for how the operations are
  !> indexed). Used to turn a selection of operations into a selection of
  !> the element types drawn for them. mask has one entry per type in se.
  module subroutine symelem_type_mask(se,iop,mask)
    type(symelem_list), intent(in) :: se
    integer, intent(in) :: iop(:)
    logical, intent(out) :: mask(:)

    integer :: i, j

    mask = .false.
    do i = 1, min(se%ntype,size(mask,1))
       do j = 1, se%nop(i)
          if (any(iop == se%iop(j,i))) then
             mask(i) = .true.
             exit
          end if
       end do
    end do

  end subroutine symelem_type_mask

  !> Classify the rotation part rmat of a symmetry operation (crystallographic
  !> coordinates) into cl: the kind of geometric element it defines, the order
  !> and Hermann-Mauguin symbol of the point operation, the axis or plane
  !> normal in several frames, the projector on the space it leaves fixed, and
  !> what the screw or glide part of the symbol is measured against. cl%kind is
  !> zero for the identity, which has no element.
  module subroutine symop_classify(c,rmat,cl)
    use tools_math, only: eig, det3
    use tools_io, only: string
    use param, only: pi, eye
    class(crystal), intent(in) :: c
    real*8, intent(in) :: rmat(3,3)
    type(symop_class), intent(inout) :: cl

    real*8, parameter :: eps = 1d-5

    integer :: j, k, ord, idx, ier, ivec(3)
    real*8 :: rmatv(3,3), wacc(3,3), wsum(3,3), eval(3), evali(3)
    real*8 :: trace, det, ang, ridx
    logical :: isplane, ok

    cl%kind = 0
    cl%order = 0
    cl%rotnum = 0
    cl%base = ""
    cl%axis = 0d0
    cl%dirc = 0d0
    cl%dirf = 0d0
    cl%dirlabel = ""
    cl%vsh = 0d0
    cl%nv = 0
    cl%pmat = eye

    trace = rmat(1,1) + rmat(2,2) + rmat(3,3)
    det = det3(rmat)
    cl%isproper = (det > 0d0)

    ! the identity has no element, and its projector is the identity
    if (abs(trace - 3d0) < eps) then
       cl%base = "1"
       return
    end if

    ! order of the rotation part and projector on the space it fixes:
    ! P = (1/k) sum_m W^m, with k the order of W
    wacc = eye
    wsum = 0d0
    ord = 0
    do k = 1, 6
       wsum = wsum + wacc
       wacc = matmul(rmat,wacc)
       ord = ord + 1
       if (all(abs(wacc - eye) < eps)) exit
    end do
    cl%pmat = wsum / real(ord,8)

    ! the inversion center has no direction
    if (abs(trace + 3d0) < eps) then
       cl%kind = symop_kind_point
       cl%base = "-1"
       return
    end if

    ! rotation angle -> order of the rotation; an improper operation -n has
    ! cos(ang) = -(trace+1)/2, so a mirror comes out as n = 2
    if (cl%isproper) then
       ang = 0.5d0*(trace-1d0)
    else
       ang = -0.5d0*(trace+1d0)
    end if
    ang = acos(max(min(ang,1d0),-1d0))
    if (abs(ang) < eps) then
       cl%rotnum = 1
    else
       cl%rotnum = nint((2d0*pi)/ang)
    end if
    isplane = (.not.cl%isproper .and. cl%rotnum == 2)

    ! base Hermann-Mauguin symbol of the point operation
    if (cl%isproper) then
       cl%base = string(cl%rotnum)
    elseif (isplane) then
       cl%base = "m"
    else
       cl%base = string(-cl%rotnum)
    end if

    ! the axis (rotations and rotoinversions) or the plane normal is the
    ! eigenvector of W with eigenvalue det; a failed diagonalization leaves
    ! the axis at zero and the operation without an element
    rmatv = rmat
    call eig(rmatv,3,eval,evali,ier)
    idx = 0
    if (ier == 0) then
       do j = 1, 3
          if (abs(evali(j)) < eps .and. abs(eval(j)-det) < eps) then
             idx = j
             exit
          end if
       end do
    end if
    if (idx == 0) return
    cl%axis = rmatv(:,idx)

    ! scale the axis by its smallest non-zero component, so that it comes out
    ! as the integer indices of the direction whenever they are commensurate
    ridx = 1d40
    do j = 1, 3
       if (abs(cl%axis(j)) > eps .and. abs(cl%axis(j)) < ridx) ridx = abs(cl%axis(j))
    end do
    if (ridx < 1d39) cl%axis = cl%axis / ridx

    cl%dirc = c%x2c(cl%axis)
    if (norm2(cl%dirc) < 1d-10) then
       cl%axis = 0d0
       return
    end if
    cl%dirc = cl%dirc / norm2(cl%dirc)

    if (isplane) then
       ! the plane normal in reciprocal coordinates gives both the Miller
       ! indices and the plane equation in crystallographic coordinates
       cl%kind = symop_kind_plane
       cl%dirf = c%rc2rx(cl%dirc)
       call primitive_int(cl%dirf,ivec,ok)
       if (ok) cl%dirlabel = "(" // string(ivec(1)) // string(ivec(2)) // string(ivec(3)) // ")"

       ! the lattice translations the glide vector is reduced by
       call inplane_lattice_vecs(c,ivec,cl%dirf,cl%nv,cl%vlist)
    else
       cl%kind = symop_kind_axis
       cl%order = cl%rotnum
       cl%dirf = cl%axis
       call primitive_int(cl%axis,ivec,ok)
       if (ok) then
          cl%dirlabel = "[" // string(ivec(1)) // string(ivec(2)) // string(ivec(3)) // "]"

          ! shortest lattice translation along a proper rotation axis: the unit
          ! the screw component of the symbol is measured in
          if (cl%isproper) then
             cl%vsh = real(ivec,8)
             do j = 2, 6
                if (is_lattice_vec(c,real(ivec,8)/real(j,8))) cl%vsh = real(ivec,8)/real(j,8)
             end do

          end if
       end if
    end if

  end subroutine symop_classify

  !> Hermann-Mauguin symbol of the symmetry operation classified in cl whose
  !> intrinsic (screw/glide) translation is tint = cl%pmat . w, with w the
  !> translation part of the operation: the symbol of the point operation plus
  !> the screw subscript of a rotation axis or the letter of a glide plane.
  !> tred returns the translation the symbol names (crystallographic
  !> coordinates), which is tint reduced by the lattice translations the
  !> element contains -- the glide vector of the plane or the screw
  !> translation of the axis, and zero for a mirror or a pure rotation.
  module function symop_symbol(c,cl,tint,tred) result(slabel)
    use tools_io, only: string
    class(crystal), intent(in) :: c
    type(symop_class), intent(in) :: cl
    real*8, intent(in) :: tint(3)
    real*8, intent(out), optional :: tred(3)
    character(len=symlen) :: slabel

    real*8, parameter :: eps = 1d-5

    integer :: ip
    real*8 :: tred_(3), fscrew

    slabel = cl%base
    tred_ = 0d0
    if (cl%kind == symop_kind_plane) then
       ! mirror or glide: the letter of the reduced glide vector
       call reduce_glide(c,cl%nv,cl%vlist,tint,tred_)
       if (norm2(c%x2c(tred_)) > eps) then
          slabel = glide_letter(tred_)
       else
          tred_ = 0d0
       end if
    elseif (cl%kind == symop_kind_axis .and. norm2(cl%vsh) > 0d0) then
       ! rotation or screw: the component along the axis, in units of the
       ! shortest lattice translation along it
       fscrew = dot_product(tint,cl%vsh) / dot_product(cl%vsh,cl%vsh)
       ip = modulo(nint(cl%rotnum*fscrew),cl%rotnum)
       if (ip /= 0) then
          slabel = trim(cl%base) // "_" // string(ip)
          tred_ = (real(ip,8) / real(cl%rotnum,8)) * cl%vsh
       end if
    end if
    if (present(tred)) tred = tred_

  end function symop_symbol

  !> Deallocate and reset a symmetry-element list.
  module subroutine symelem_list_end(se)
    class(symelem_list), intent(inout) :: se

    se%ntype = 0
    se%n = 0
    if (allocated(se%kind)) deallocate(se%kind)
    if (allocated(se%order)) deallocate(se%order)
    if (allocated(se%dir)) deallocate(se%dir)
    if (allocated(se%label)) deallocate(se%label)
    if (allocated(se%dirlabel)) deallocate(se%dirlabel)
    if (allocated(se%nop)) deallocate(se%nop)
    if (allocated(se%iop)) deallocate(se%iop)
    se%maxop = 0
    if (allocated(se%itype)) deallocate(se%itype)
    if (allocated(se%x)) deallocate(se%x)
    if (allocated(se%key)) deallocate(se%key)
    if (allocated(se%tint)) deallocate(se%tint)

  end subroutine symelem_list_end

  !> Return in se the geometric symmetry elements (mirror/glide
  !> planes, rotation/screw/rotoinversion axes, and inversion centers)
  !> of the structure. For crystals, an element is included if it
  !> intersects the box spanned by ncell unit cells ([0,ncell] in
  !> crystallographic coordinates); if border, the elements lying
  !> exactly on the far faces of that box are included as well. The
  !> type list (kind/order/direction/symbol) is always built from the
  !> one-cell region, so it does not change with ncell or border. For
  !> molecules, ncell and border are ignored and the point-group
  !> elements, all of which pass through the center of mass, are
  !> returned. If iop is given, return only the elements of the
  !> symmetry operations it lists, the type list is then built from
  !> those operations alone. If typesonly, return only the element
  !> types, and ncell and border are ignored.
  module subroutine list_symelems(c,ncell,border,se,iop,typesonly)
    use types, only: molsymop_plane, molsymop_rotation, molsymop_imp_rotation, molsymop_inversion
    class(crystal), intent(in) :: c
    integer, intent(in) :: ncell(3)
    logical, intent(in) :: border
    type(symelem_list), intent(inout) :: se
    integer, intent(in), optional :: iop(:)
    logical, intent(in), optional :: typesonly

    integer :: i, skind, sorder, ioptype
    real*8 :: raxx(3)
    logical :: dotypes

    call se%end()
    dotypes = .false.
    if (present(typesonly)) dotypes = typesonly

    if (c%ismolecule) then
       !! molecule: the point-group operations, all through the center of mass
       if (.not.c%pg%avail) return
       do i = 1, c%pg%nop
          if (present(iop)) then
             if (.not.any(iop == i)) cycle
          end if
          ioptype = c%pg%op(i)%type
          if (ioptype == molsymop_inversion) then
             ! the inversion center is the center of mass, and has no direction
             call symelem_add(se,symop_kind_point,0,(/0d0,0d0,0d0/),c%pg%op(i)%sym,"",c%pg%xcm,iop=i)
             cycle
          elseif (ioptype == molsymop_plane) then
             skind = symop_kind_plane
             sorder = 0
          elseif (ioptype == molsymop_rotation .or. ioptype == molsymop_imp_rotation) then
             skind = symop_kind_axis
             sorder = c%pg%op(i)%opn
          else
             cycle
          end if
          raxx = c%pg%op(i)%axis
          if (norm2(raxx) < 1d-10) cycle
          raxx = raxx / norm2(raxx)
          call symelem_add(se,skind,sorder,raxx,c%pg%op(i)%sym,"",c%pg%xcm,iop=i)
       end do
       call symelem_trim(se)
    else
       !! crystal: the one-cell region is a fundamental domain for the lattice,
       !! so the elements it contains are a complete set of representatives,
       !! and the type list they give does not depend on ncell or border
       call symelem_enum(c,se,iop)

       !! the elements in the requested region are their lattice translates
       if (.not.dotypes) call symelem_replicate(c,ncell,border,se)
    end if

    !! the caller only wants the element types: drop the instances the
    !! enumeration built on the way (the types are a by-product of it)
    if (dotypes) then
       se%n = 0
       if (allocated(se%itype)) deallocate(se%itype)
       if (allocated(se%x)) deallocate(se%x)
       if (allocated(se%key)) deallocate(se%key)
       call symelem_trim(se)
    end if

  end subroutine list_symelems

  !xx! private procedures

  !> For the symmetry operation rot, compute the order and
  !> characteristic direction of the operation. The output is given in
  !> type (the type of operation: ident, c2, c3,...)  and vec (the
  !> direction of the axis or the normal to the plane, 0 for a point
  !> symmetry element). Used in the structure initialization.
  subroutine typeop(rot,type,vec,order)
    use tools_math, only: eig
    use tools_io, only: ferror, faterr
    use param, only: tpi

    real*8, intent(in) :: rot(3,4) !< rotm operation
    integer, intent(out) :: type !< output type
    real*8, dimension(3), intent(out) :: vec !< output eigenvector
    integer, intent(out) :: order

    real*8, parameter :: eps = 1e-3

    integer, dimension(2,2) :: iord
    integer :: ier
    real*8 :: trace, tone
    real*8, dimension(3) :: eigen, eigeni
    real*8, dimension(3,3) :: mat
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
       call eig(mat,3,eigen,eigeni,ier)
       if (ier /= 0) &
          call ferror('typeop','Error in diagonalization',faterr)

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
          ! the order of the rotoinversion as a group element: a 6-bar
          ! axis (-C6) has order 3 (S3), a 4-bar (-C4) order 4 (S4) and
          ! a 3-bar (-C3) order 6 (S6)
          if (order .eq. 3) then
             type = s3
          elseif (order .eq. 4) then
             type = s4
          elseif (order .eq. 6) then
             type = s6
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
    use tools_math, only: cellpar_from_metric

    real*8, intent(in) :: x2r(3,3), rmax
    integer, intent(out) :: imax, jmax, kmax

    integer :: nadd
    real*8 :: acell(3), bcell(3), gmat(3,3)

    gmat = matmul(transpose(x2r),x2r)
    call cellpar_from_metric(gmat,acell,bcell)

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

  ! Calculate the symmetry operations and point group of the crystal
  ! (treated as a molecule, using the complete atom list). This is a thin
  ! wrapper around molsymmod's calc_point_group.
  module subroutine calcmolsym(c,pg,errmsg)
    use molsymmod, only: calc_point_group
    class(crystal), intent(inout) :: c
    type(point_group), intent(inout) :: pg
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, nat
    real*8, allocatable :: x(:,:)
    integer, allocatable :: z(:)

    nat = c%ncel
    allocate(x(3,nat),z(nat))
    do i = 1, nat
       x(:,i) = c%atcel(i)%r
       z(i) = c%spc(c%atcel(i)%is)%z
    end do
    call calc_point_group(nat,x,z,pg,errmsg)

  end subroutine calcmolsym

  !> Symmetrize the positions of the atoms in a molecule using the
  !> current molecular point group (c%pg). errmsg has length > 0 on
  !> error.
  module subroutine symmetrize_molecule(c,errmsg)
    use tools_math, only: matinv
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, k, iop, kmin, n, it, ier
    integer, allocatable :: icount(:)
    real*8, allocatable :: xc(:,:), xacc(:,:), xnew(:,:)
    real*8 :: y(3), d2, d2min, rot(3,3), rinv(3,3)

    integer, parameter :: maxit_orth = 20 ! max Newton iterations for orthogonalization
    real*8, parameter :: eps_orth = 1d-14 ! orthogonalization convergence criterion

    errmsg = ""

    ! consistency checks
    if (.not.c%ismolecule) then
       errmsg = "symmetrize_molecule can only be used with molecules"
       return
    end if
    if (.not.c%pg%avail) then
       errmsg = "the molecular point group is not available"
       return
    end if

    ! atomic positions relative to the point-group center
    n = c%ncel
    allocate(xc(3,n),xacc(3,n),icount(n))
    do i = 1, n
       xc(:,i) = c%atcel(i)%r - c%pg%xcm
    end do

    ! accumulate the image of every atom under every idealized operation
    ! onto the nearest atom of the same species
    xacc = 0d0
    icount = 0
    do iop = 1, c%pg%nop
       ! idealize the operation: Newton iteration for the orthogonal
       ! factor of the polar decomposition, R <- (R + R^-T) / 2
       rot = c%pg%op(iop)%m
       do it = 1, maxit_orth
          rinv = rot
          call matinv(rinv,3,ier)
          if (ier /= 0) then
             errmsg = "singular point-group operation matrix"
             return
          end if
          rinv = transpose(rinv)
          if (all(abs(rot - rinv) < eps_orth)) exit
          rot = 0.5d0 * (rot + rinv)
       end do

       do i = 1, n
          y = matmul(rot,xc(:,i))
          kmin = 0
          d2min = huge(1d0)
          do k = 1, n
             if (c%atcel(k)%is /= c%atcel(i)%is) cycle
             d2 = dot_product(y - xc(:,k),y - xc(:,k))
             if (d2 < d2min) then
                d2min = d2
                kmin = k
             end if
          end do
          xacc(:,kmin) = xacc(:,kmin) + y
          icount(kmin) = icount(kmin) + 1
       end do
    end do

    ! every operation permutes the atoms, so each atom must have
    ! received exactly one image per operation; otherwise the point
    ! group is inconsistent with the current geometry
    if (any(icount /= c%pg%nop)) then
       errmsg = "the point group operations are inconsistent with the molecular geometry"
       return
    end if

    ! average the images, apply the new positions, and rebuild the
    ! structure (update_positions alone is a display-only update)
    allocate(xnew(3,n))
    do i = 1, n
       xnew(:,i) = c%pg%xcm + xacc(:,i) / real(c%pg%nop,8)
    end do
    call c%update_positions(xnew)
    ! the symmetrization only nudges the positions: keep the bonding
    call c%rebuild_after_move(copybonding=.true.,errmsg=errmsg)

  end subroutine symmetrize_molecule

  !> Reduce the complete cell atom list (c%atcel) to the non-equivalent
  !> atom list (c%at) given the symmetry operations already set in the
  !> crystal (neqv, ncv, rotm, cen). On input c%atcel(1:c%ncel) must be
  !> populated and name(:) holds the per-cell-atom names (celatom has no name
  !> field). For non-trivial symmetry the block spatial index must be built
  !> (c%build_env) on input; the P1 case (neqv=ncv=1) takes a fast path that
  !> needs no index. On output nneq, the reduced at() list (with
  !> name/mult/wyc), and the complete-list symmetry mapping
  !> (atcel%idx/cidx/ir/ic/lvec) are filled. If error, errmsg has length > 0.
  module subroutine reduceatoms(c,name,errmsg,occ)
    use param, only: icrd_crys
    use global, only: symprec
    use types, only: realloc
    class(crystal), intent(inout) :: c
    character*10, intent(in) :: name(:)
    character(len=:), allocatable, intent(out) :: errmsg
    real*8, intent(in), optional :: occ(:)

    integer :: k, io, it, id
    real*8 :: x0(3)
    logical, allocatable :: done(:)

    errmsg = ""

    ! P1: each atom is its own representative (1:1 at <-> atcel), so
    ! no orbit walk or spatial index is needed.
    if (c%neqv == 1 .and. c%ncv == 1) then
       c%nneq = c%ncel
       call realloc(c%at,c%nneq)
       do k = 1, c%ncel
          c%at(k)%x = c%atcel(k)%x
          c%at(k)%r = c%atcel(k)%r
          c%at(k)%is = c%atcel(k)%is
          c%at(k)%name = name(k)
          c%at(k)%wyc = "?"
          c%at(k)%mult = 1
          if (present(occ)) then
             c%at(k)%occ = occ(k)
          else
             c%at(k)%occ = 1d0
          end if
          c%atcel(k)%idx = k
          c%atcel(k)%cidx = k
          c%atcel(k)%ir = 1
          c%atcel(k)%ic = 1
          c%atcel(k)%lvec = 0
       end do
       return
    end if

    ! Orbit reduction: build the non-equivalent atom list, the
    ! multiplicities, and the complete-list symmetry mapping in a single
    ! O(ncel*neqv*ncv) pass. Each atom is assigned to the first
    ! representative whose orbit reaches it, using O(1) identify_atom lookups.
    allocate(done(c%ncel))
    done = .false.
    c%nneq = 0
    do k = 1, c%ncel
       if (done(k)) cycle
       c%nneq = c%nneq + 1
       if (c%nneq > size(c%at,1)) call realloc(c%at,2*c%nneq)
       c%at(c%nneq)%x = c%atcel(k)%x
       c%at(c%nneq)%r = c%atcel(k)%r
       c%at(c%nneq)%is = c%atcel(k)%is
       c%at(c%nneq)%name = name(k)
       c%at(c%nneq)%wyc = "?"
       c%at(c%nneq)%mult = 0
       if (present(occ)) then
          c%at(c%nneq)%occ = occ(k)
       else
          c%at(c%nneq)%occ = 1d0
       end if
       do io = 1, c%neqv
          do it = 1, c%ncv
             x0 = matmul(c%rotm(1:3,1:3,io),c%atcel(k)%x) + c%rotm(:,4,io) + c%cen(:,it)
             id = c%identify_atom(x0,icrd_crys,distmax=symprec)
             if (id == 0) then
                errmsg = "incomplete orbit while building the complete atom list"
                return
             end if
             if (.not.done(id)) then
                done(id) = .true.
                c%atcel(id)%idx = c%nneq
                c%atcel(id)%cidx = id
                c%atcel(id)%ir = io
                c%atcel(id)%ic = it
                c%atcel(id)%lvec = nint(c%atcel(id)%x - x0)
                c%at(c%nneq)%mult = c%at(c%nneq)%mult + 1
             end if
          end do
       end do
    end do
    call realloc(c%at,c%nneq)
    deallocate(done)

  end subroutine reduceatoms


  !> Enumerate the symmetry elements of crystal c that intersect the box
  !> [0,ncell] in crystallographic coordinates (plus the elements on its far
  !> faces, if border) and return them in se.
  !>
  !> For an operation (W,w) of order k (of the rotation part W), the projector
  !> on the space left fixed by W is P = (1/k) sum_m W^m. The operation splits
  !> into the intrinsic (screw/glide) translation t = P.w, which gives the
  !> Hermann-Mauguin symbol, and the location part wl = w - t. The element is
  !> the set of fixed points of (W,wl), and it passes through
  !>    x0 = (W - I + P)^-1 . (-wl)
  !> Every lattice translation u added to w gives another element of
  !> the group, generally at a different location and with a different
  !> symbol (e.g. the mirror and n-glide planes that alternate along
  !> [110] in P-42_1m), so the enumeration runs over u as well.
  subroutine symelem_enum(c,se,iop)
    use tools_math, only: matinv
    use param, only: eye
    type(crystal), intent(in) :: c
    type(symelem_list), intent(inout) :: se
    integer, intent(in), optional :: iop(:)

    real*8, parameter :: eps = 1d-5

    integer :: ieqv, icv, iop1, j, ier
    integer :: tlim(3), it1, it2, it3
    real*8 :: nmat(3,3), wvec(3), tvec(3), tint(3), tred(3), wl(3), x0(3), xmid(3)
    character(len=symlen) :: slabel
    logical :: ok
    type(symop_class) :: cl
    type(elem_box) :: box

    call se%end()

    ! the elements are enumerated in one cell, which is a fundamental domain
    ! for the lattice; symelem_replicate then translates them over the cells
    call box%set((/0d0,0d0,0d0/),eye,ok)
    if (.not.ok) return

    do ieqv = 1, c%neqv
       ! the rotation part fixes the kind, the direction, and the projector,
       ! and is shared by every centering vector
       call symop_classify(c,c%rotm(:,1:3,ieqv),cl)
       if (cl%kind == 0) cycle

       ! (W - I + P) is invertible and inverts W - I on the complement of the
       ! fixed space, where the location part lives
       nmat = c%rotm(:,1:3,ieqv) - eye + cl%pmat
       call matinv(nmat,3,ier)
       if (ier /= 0) cycle

       ! a lattice translation along the element does not move it, so only the
       ! directions with a transverse component need to be scanned
       do j = 1, 3
          tvec = -cl%pmat(:,j)
          tvec(j) = tvec(j) + 1d0
          if (norm2(c%x2c(tvec)) < eps) then
             tlim(j) = 0
          else
             tlim(j) = tlim_def
          end if
       end do

       do icv = 1, c%ncv
          ! only the requested operations, if a list was given
          iop1 = (icv-1)*c%neqv + ieqv
          if (present(iop)) then
             if (.not.any(iop == iop1)) cycle
          end if
          wvec = c%rotm(:,4,ieqv) + c%cen(:,icv)

          do it1 = -tlim(1), tlim(1)
             do it2 = -tlim(2), tlim(2)
                do it3 = -tlim(3), tlim(3)
                   tvec = wvec + real((/it1,it2,it3/),8)
                   tint = cl%pmat(:,1)*tvec(1) + cl%pmat(:,2)*tvec(2) + cl%pmat(:,3)*tvec(3)
                   wl = tint - tvec
                   x0 = nmat(:,1)*wl(1) + nmat(:,2)*wl(2) + nmat(:,3)*wl(3)

                   ! the location part of the operation, and the point of the
                   ! element it puts there (the matmuls are written out: gfortran
                   ! emits a library call for matmul at the -Og the debug build
                   ! uses, and this is the innermost loop)

                   ! keep only the elements that reach the cell, and move the
                   ! point to the middle of the part of the element inside it
                   call symelem_in_box(cl%kind,x0,cl%dirf,box,ok,xmid)
                   if (.not.ok) cycle

                   ! the symbol of this element depends on the translation
                   slabel = symop_symbol(c,cl,tint,tred=tred)
                   call symelem_add(se,cl%kind,cl%order,cl%dirc,slabel,cl%dirlabel,&
                      c%x2c(xmid),iop=iop1,tintc=c%x2c(tred))
                end do
             end do
          end do
       end do
    end do

  end subroutine symelem_enum

  !> Replace the element instances in se (a complete set of representatives,
  !> one per lattice class) by all their lattice translates that intersect the
  !> box [0,ncell] in crystallographic coordinates, plus those on its far faces
  !> if border. The element types are left untouched.
  subroutine symelem_replicate(c,ncell,border,se)
    type(crystal), intent(in) :: c
    integer, intent(in) :: ncell(3)
    logical, intent(in) :: border
    type(symelem_list), intent(inout) :: se

    real*8, parameter :: eps = 1d-5

    type(symelem_list) :: sen
    integer :: i, j, it, it1, it2, it3, tlo(3), thi(3), skind, ifrom
    real*8 :: x0f(3), xf(3), xmid(3), dirc(3), dirf(3), ev(3), bv(3,3)
    logical :: ok, onfar
    type(elem_box) :: box

    if (se%n == 0) return

    ! the box the elements are kept in: the displayed cells
    bv = 0d0
    do i = 1, 3
       bv(i,i) = real(ncell(i),8)
    end do
    call box%set((/0d0,0d0,0d0/),bv,ok)
    if (.not.ok) return

    ! carry over the type list, so that symelem_add reuses it
    do i = 1, se%ntype
       call symelem_addtype(sen,se%kind(i),se%order(i),se%dir(:,i),se%label(i),se%dirlabel(i),&
          se%nop(i),se%iop(:,i))
    end do

    ! one type at a time: the geometry is the same for every element of a type,
    ! and the instances of a type come out contiguous, so the duplicate scan in
    ! symelem_add only has to look at the type being built
    do it = 1, se%ntype
       ifrom = sen%n + 1
       skind = se%kind(it)
       dirc = se%dir(:,it)
       dirf = symelem_dirf(c,skind,dirc)

       ! a lattice translation along the element does not move it: the
       ! representative it generates is the element itself. Every other
       ! translate that reaches the box is generated by a translation in
       ! [-1,ncell] (see the discussion in symelem_enum).
       do j = 1, 3
          ev = c%m_x2c(:,j)
          if (skind == symop_kind_plane) then
             ok = (abs(dirf(j)) < eps)
          elseif (skind == symop_kind_axis) then
             ok = (norm2(ev - dot_product(ev,dirc)*dirc) < eps)
          else
             ok = .false.
          end if
          if (ok) then
             tlo(j) = 0
             thi(j) = 0
          else
             tlo(j) = -1
             thi(j) = ncell(j)
          end if
       end do

       do i = 1, se%n
          if (se%itype(i) /= it) cycle
          x0f = c%c2x(se%x(:,i))
          do it1 = tlo(1), thi(1)
             do it2 = tlo(2), thi(2)
                do it3 = tlo(3), thi(3)
                   xf = x0f + real((/it1,it2,it3/),8)

                   ! keep only the elements that reach the box, and move the
                   ! point to the middle of the part of the element inside it
                   call symelem_in_box(skind,xf,dirf,box,ok,xmid,onfar)
                   if (.not.ok) cycle
                   if (.not.border .and. onfar) cycle
                   call symelem_add(sen,skind,se%order(it),dirc,se%label(it),se%dirlabel(it),&
                      c%x2c(xmid),ifrom=ifrom,itype=it,tintc=se%tint(:,i))
                end do
             end do
          end do
       end do
    end do
    call symelem_trim(sen)
    call se%end()
    se%ntype = sen%ntype
    se%maxop = sen%maxop
    se%n = sen%n
    call move_alloc(sen%kind,se%kind)
    call move_alloc(sen%order,se%order)
    call move_alloc(sen%dir,se%dir)
    call move_alloc(sen%label,se%label)
    call move_alloc(sen%dirlabel,se%dirlabel)
    call move_alloc(sen%nop,se%nop)
    call move_alloc(sen%iop,se%iop)
    call move_alloc(sen%tint,se%tint)
    call move_alloc(sen%itype,se%itype)
    call move_alloc(sen%x,se%x)
    call move_alloc(sen%key,se%key)

  end subroutine symelem_replicate

  !> Direction of a symmetry element of kind skind in the crystallographic
  !> frame: the reciprocal-space normal for a plane (so that the plane is
  !> dirf.(x-x0) = 0), the direction for an axis, and zero for an inversion
  !> center. dirc is the Cartesian unit direction or normal.
  function symelem_dirf(c,skind,dirc) result(dirf)
    type(crystal), intent(in) :: c
    integer, intent(in) :: skind
    real*8, intent(in) :: dirc(3)
    real*8 :: dirf(3)

    if (skind == symop_kind_plane) then
       dirf = c%rc2rx(dirc)
    elseif (skind == symop_kind_axis) then
       dirf = c%c2x(dirc)
    else
       dirf = 0d0
    end if

  end function symelem_dirf

  !> Hermann-Mauguin letter of a glide plane whose (reduced, non-zero) glide
  !> vector in crystallographic coordinates is tred.
  function glide_letter(tred) result(c1)
    real*8, intent(in) :: tred(3)
    character(len=1) :: c1

    integer :: nhalf

    if (any(abs(abs(tred)-0.25d0) < 1d-3)) then
       c1 = "d"
       return
    end if
    nhalf = count(abs(abs(tred)-0.5d0) < 1d-3)
    if (nhalf >= 2) then
       c1 = "n"
    elseif (nhalf == 0) then
       c1 = "g"
    elseif (abs(abs(tred(1))-0.5d0) < 1d-3) then
       c1 = "a"
    elseif (abs(abs(tred(2))-0.5d0) < 1d-3) then
       c1 = "b"
    else
       c1 = "c"
    end if

  end function glide_letter

  !> Index of the element type (skind,sorder,dirc,slabel) in the element list
  !> se, or zero if it is not there. Axes and planes are the same type if they
  !> are parallel, whichever way their direction vector points.
  function symelem_findtype(se,skind,sorder,dirc,slabel) result(it)
    type(symelem_list), intent(in) :: se
    integer, intent(in) :: skind, sorder
    real*8, intent(in) :: dirc(3)
    character*(*), intent(in) :: slabel
    integer :: it

    integer :: i

    it = 0
    do i = 1, se%ntype
       if (se%kind(i) /= skind) cycle
       if (se%order(i) /= sorder) cycle
       if (se%label(i) /= slabel) cycle
       if (skind /= symop_kind_point) then
          if (abs(abs(dot_product(se%dir(:,i),dirc)) - 1d0) > 1d-6) cycle
       end if
       it = i
       return
    end do

  end function symelem_findtype

  !> Add the element type (skind,sorder,dirc,slabel,sdirlabel) to the element
  !> list se. dirc is the Cartesian unit direction (axes) or normal (planes).
  subroutine symelem_addtype(se,skind,sorder,dirc,slabel,sdirlabel,nop,iop)
    use types, only: realloc
    type(symelem_list), intent(inout) :: se
    integer, intent(in) :: skind, sorder
    real*8, intent(in) :: dirc(3)
    character*(*), intent(in) :: slabel, sdirlabel
    integer, intent(in) :: nop
    integer, intent(in) :: iop(:)

    integer :: n

    se%ntype = se%ntype + 1
    if (.not.allocated(se%kind)) then
       se%maxop = maxop_def
       allocate(se%kind(10),se%order(10),se%dir(3,10),se%label(10),se%dirlabel(10))
       allocate(se%nop(10),se%iop(se%maxop,10))
    elseif (se%ntype > size(se%kind,1)) then
       n = 2*se%ntype
       call realloc(se%kind,n)
       call realloc(se%order,n)
       call realloc(se%dir,3,n)
       call realloc(se%label,n)
       call realloc(se%dirlabel,n)
       call realloc(se%nop,n)
       call realloc(se%iop,se%maxop,n)
    end if
    se%kind(se%ntype) = skind
    se%order(se%ntype) = sorder
    se%dir(:,se%ntype) = dirc
    se%label(se%ntype) = slabel
    se%dirlabel(se%ntype) = sdirlabel
    se%nop(se%ntype) = min(nop,se%maxop)
    se%iop(:,se%ntype) = 0
    if (se%nop(se%ntype) > 0) se%iop(1:se%nop(se%ntype),se%ntype) = iop(1:se%nop(se%ntype))

  end subroutine symelem_addtype

  !> Add one symmetry element of kind skind, rotation order sorder, Cartesian
  !> unit direction/normal dirc, symbol slabel, and direction label sdirlabel,
  !> passing through the Cartesian point xc, to the element list se. Creates
  !> the element type if it is not in the list yet. Elements already in the
  !> list (same type and same position transverse to the element) are skipped;
  !> ifrom is the first instance the duplicate scan has to look at.
  subroutine symelem_add(se,skind,sorder,dirc,slabel,sdirlabel,xc,ifrom,iop,itype,tintc)
    use types, only: realloc
    type(symelem_list), intent(inout) :: se
    integer, intent(in) :: skind, sorder
    real*8, intent(in) :: dirc(3), xc(3)
    character*(*), intent(in) :: slabel, sdirlabel
    integer, intent(in), optional :: ifrom
    integer, intent(in), optional :: iop
    integer, intent(in), optional :: itype
    real*8, intent(in), optional :: tintc(3)

    real*8, parameter :: epskey = 1d-3 ! bohr

    integer :: i, it, n, i0
    real*8 :: key(3)

    ! find the element type, or create it
    ! the caller may already know the type (the replication does)
    if (present(itype)) then
       it = itype
    else
       it = symelem_findtype(se,skind,sorder,dirc,slabel)
       if (it == 0) then
          call symelem_addtype(se,skind,sorder,dirc,slabel,sdirlabel,0,(/0/))
          it = se%ntype
       end if
    end if

    ! remember which symmetry operations generate this element type
    if (present(iop)) then
       if (.not.any(se%iop(1:se%nop(it),it) == iop) .and. se%nop(it) < se%maxop) then
          se%nop(it) = se%nop(it) + 1
          se%iop(se%nop(it),it) = iop
       end if
    end if

    ! discard the element if it is in the list already (the key of every
    ! instance is cached, so this scan does no geometry). A caller that adds
    ! the elements of one type at a time passes the index where that type
    ! starts, so the scan does not run over the types already done.
    i0 = 1
    if (present(ifrom)) i0 = ifrom
    key = symelem_key(skind,se%dir(:,it),xc)
    do i = i0, se%n
       if (se%itype(i) /= it) cycle
       if (all(abs(se%key(:,i) - key) < epskey)) return
    end do

    ! add it
    se%n = se%n + 1
    if (.not.allocated(se%itype)) then
       allocate(se%itype(100),se%x(3,100),se%key(3,100),se%tint(3,100))
    elseif (se%n > size(se%itype,1)) then
       n = 2*se%n
       call realloc(se%itype,n)
       call realloc(se%x,3,n)
       call realloc(se%key,3,n)
       call realloc(se%tint,3,n)
    end if
    se%itype(se%n) = it
    se%x(:,se%n) = xc
    se%tint(:,se%n) = 0d0
    if (present(tintc)) se%tint(:,se%n) = tintc
    se%key(:,se%n) = key

  end subroutine symelem_add

  !> Trim the arrays in symmetry-element list se to their actual size.
  subroutine symelem_trim(se)
    use types, only: realloc
    type(symelem_list), intent(inout) :: se

    if (allocated(se%kind)) then
       if (se%ntype /= size(se%kind,1)) then
          call realloc(se%kind,se%ntype)
          call realloc(se%order,se%ntype)
          call realloc(se%dir,3,se%ntype)
          call realloc(se%label,se%ntype)
          call realloc(se%dirlabel,se%ntype)
          call realloc(se%nop,se%ntype)
          call realloc(se%iop,se%maxop,se%ntype)
       end if
    end if
    if (allocated(se%itype)) then
       if (se%n /= size(se%itype,1)) then
          call realloc(se%itype,se%n)
          call realloc(se%x,3,se%n)
          call realloc(se%key,3,se%n)
          call realloc(se%tint,3,se%n)
       end if
    end if

  end subroutine symelem_trim

  !> Position of the Cartesian point xc transverse to a symmetry element of
  !> kind skind with Cartesian unit direction/normal dirc. Two elements of the
  !> same type are the same element if and only if this vector is the same.
  function symelem_key(skind,dirc,xc) result(key)
    integer, intent(in) :: skind
    real*8, intent(in) :: dirc(3), xc(3)
    real*8 :: key(3)

    if (skind == symop_kind_plane) then
       key = dot_product(xc,dirc) * dirc
    elseif (skind == symop_kind_axis) then
       key = xc - dot_product(xc,dirc) * dirc
    else
       key = xc
    end if

  end function symelem_key

  !> Scale vector v to the smallest collinear vector with integer components,
  !> returned in iv and with the first non-zero component positive. ok is
  !> false if v is zero or its components are not commensurate.
  subroutine primitive_int(v,iv,ok)
    use tools_math, only: gcd
    real*8, intent(in) :: v(3)
    integer, intent(out) :: iv(3)
    logical, intent(out) :: ok

    real*8, parameter :: eps = 1d-4

    integer :: i, m, ig
    real*8 :: vv(3), vmax

    iv = 0
    ok = .false.
    vmax = maxval(abs(v))
    if (vmax < 1d-10) return
    vv = v / vmax

    ! smallest multiplier that makes all the components integers
    do m = 1, 24
       if (all(abs(m*vv - nint(m*vv)) < eps)) then
          iv = nint(m*vv)
          ok = .true.
          exit
       end if
    end do
    if (.not.ok) return

    ! reduce by the greatest common divisor
    ig = gcd(abs(iv),3)
    if (ig > 1) iv = iv / ig

    ! first non-zero component positive
    do i = 1, 3
       if (iv(i) /= 0) then
          if (iv(i) < 0) iv = -iv
          exit
       end if
    end do

  end subroutine primitive_int

  !> True if v (crystallographic coordinates) is a lattice translation of
  !> crystal c, centering vectors included.
  function is_lattice_vec(c,v) result(ok)
    type(crystal), intent(in) :: c
    real*8, intent(in) :: v(3)
    logical :: ok

    integer :: i
    real*8 :: d(3)

    ok = .true.
    do i = 1, c%ncv
       d = v - c%cen(:,i)
       if (all(abs(d - nint(d)) < 1d-5)) return
    end do
    ok = .false.

  end function is_lattice_vec

  !> Return in vlist the nv short lattice translations of crystal c (centering
  !> vectors included) that are contained in the plane with primitive Miller
  !> indices ihkl and reciprocal-space normal hklr. Crystallographic
  !> coordinates. The search box has to reach past the Miller indices: the
  !> shortest lattice translation in a plane such as (2-10) is (1,2,0), so a
  !> box of +/-1 would find no in-plane translation at all and leave the glide
  !> vectors unreduced. A plane with indices (h,k,l) always has a basis of its
  !> lattice translations with components no larger than max(|h|,|k|,|l|), so
  !> the box is taken from the indices themselves and not guessed; that also
  !> covers the non-conventional cells a NEWCELL transformation can produce,
  !> where the indices are large.
  subroutine inplane_lattice_vecs(c,ihkl,hklr,nv,vlist)
    type(crystal), intent(in) :: c
    integer, intent(in) :: ihkl(3)
    real*8, intent(in) :: hklr(3)
    integer, intent(out) :: nv
    real*8, allocatable, intent(inout) :: vlist(:,:)

    integer :: i, j1, j2, j3, mb

    real*8 :: v(3), hnorm

    mb = max(maxval(abs(ihkl)),1)
    if (allocated(vlist)) then
       if (size(vlist,2) < (2*mb+1)**3 * c%ncv) deallocate(vlist)
    end if
    if (.not.allocated(vlist)) allocate(vlist(3,(2*mb+1)**3 * c%ncv))

    hnorm = max(norm2(hklr),1d-10)
    nv = 0
    do i = 1, c%ncv
       do j1 = -mb, mb
          do j2 = -mb, mb
             do j3 = -mb, mb
                v = c%cen(:,i) + real((/j1,j2,j3/),8)
                if (abs(dot_product(hklr,v)) > 1d-5 * hnorm * max(norm2(v),1d0)) cycle
                nv = nv + 1
                vlist(:,nv) = v
             end do
          end do
       end do
    end do

  end subroutine inplane_lattice_vecs

  !> Reduce the glide vector t of a plane modulo the nv lattice translations
  !> vlist contained in that plane (as given by inplane_lattice_vecs), and
  !> return the shortest representative in tred. Crystallographic coordinates.
  subroutine reduce_glide(c,nv,vlist,t,tred)
    type(crystal), intent(in) :: c
    integer, intent(in) :: nv
    real*8, intent(in) :: vlist(3,nv), t(3)
    real*8, intent(out) :: tred(3)

    integer :: i, iter, ir, irbest
    real*8 :: tt(3), tbest(3), dmin, d
    logical :: again

    ! greedy reduction: subtract the in-plane lattice translation that shortens
    ! the glide vector the most, and repeat while it keeps getting shorter
    tred = t
    dmin = norm2(c%x2c(t))
    do iter = 1, 20
       again = .false.
       do i = 1, nv
          tt = tred - vlist(:,i)
          d = norm2(c%x2c(tt))
          if (d < dmin - 1d-6) then
             dmin = d
             tred = tt
             again = .true.
          end if
       end do
       if (.not.again) exit
    end do

    ! Among the lattice-equivalent glide vectors, keep the one that carries the
    ! conventional Hermann-Mauguin letter: the shortest one need not be it. In
    ! R-3c, for instance, the c glide (0,0,1/2) has equivalents that are shorter
    ! in a hexagonal cell but that no letter describes.
    irbest = glide_rank(tred)
    tbest = tred
    do i = 1, nv
       tt = tred - vlist(:,i)
       ir = glide_rank(tt)
       d = norm2(c%x2c(tt))
       if (ir < irbest .or. (ir == irbest .and. d < dmin - 1d-6)) then
          irbest = ir
          dmin = d
          tbest = tt
       end if
    end do
    tred = tbest

  end subroutine reduce_glide

  !> Rank of the Hermann-Mauguin letter of the glide vector t: the lower, the
  !> more conventional. Used to choose among the lattice-equivalent glide
  !> vectors of one plane: a mirror first, then an axial glide, then n, then d,
  !> and last a glide no letter describes. A vector outside the fundamental
  !> range ranks last, as the letters do not apply to it.
  function glide_rank(t) result(ir)
    real*8, intent(in) :: t(3)
    integer :: ir

    ! the letters only describe a glide vector reduced to the fundamental
    ! range: (-3/2,0,-1/2) would be read as a c glide, but it is not one
    if (any(abs(t) > 0.5d0 + 1d-3)) then
       ir = 5
       return
    end if
    if (all(abs(t) < 1d-5)) then
       ir = 0
       return
    end if
    select case (glide_letter(t))
    case ("a","b","c")
       ir = 1
    case ("n")
       ir = 2
    case ("d")
       ir = 3
    case default
       ir = 4
    end select

  end function glide_rank

  !> Set the parallelepiped box b from its origin o and edge vectors v
  !> (columns), both in the caller's coordinate frame. ok is false if the edge
  !> vectors are linearly dependent, and the box is then left zeroed, so that a
  !> caller that ignores ok gets an obviously empty box and not a half-valid
  !> one.
  module subroutine elem_box_set(b,o,v,ok)
    use tools_math, only: matinv
    class(elem_box), intent(inout) :: b
    real*8, intent(in) :: o(3)
    real*8, intent(in) :: v(3,3)
    logical, intent(out) :: ok

    integer :: ier

    b%o = o
    b%v = v
    b%vinv = v
    call matinv(b%vinv,3,ier)
    ok = (ier == 0)
    if (.not.ok) then
       b%v = 0d0
       b%vinv = 0d0
    end if

  end subroutine elem_box_set

  !> Clip the point x0 to the box b: ok is true if it is inside, and onfar
  !> says it sits on a far face of the box. See clip_plane_box for the frame
  !> conventions.
  module subroutine clip_point_box(b,x0,ok,onfar)
    type(elem_box), intent(in) :: b
    real*8, intent(in) :: x0(3)
    logical, intent(out) :: ok
    logical, intent(out), optional :: onfar

    real*8 :: xs(3), xa(3)

    if (present(onfar)) onfar = .false.

    ! into the box frame, where the box is the unit cube
    xa = x0 - b%o
    xs = b%vinv(:,1)*xa(1) + b%vinv(:,2)*xa(2) + b%vinv(:,3)*xa(3)

    ok = all(xs >= -clip_eps) .and. all(xs <= 1d0+clip_eps)
    if (ok .and. present(onfar)) onfar = any(xs > 1d0 - clip_eps)

  end subroutine clip_point_box

  !> Clip the line through x0 with direction dir to the box b. ok is false if
  !> the line misses the box, or only grazes it at a point or along an edge, as
  !> that draws nothing; otherwise the part inside the box runs from
  !> x0 + s1 * dir to x0 + s2 * dir. onfar says that part lies on a far face of
  !> the box. See clip_plane_box for the frame conventions.
  module subroutine clip_line_box(b,x0,dir,ok,s1,s2,onfar)
    type(elem_box), intent(in) :: b
    real*8, intent(in) :: x0(3)
    real*8, intent(in) :: dir(3)
    logical, intent(out) :: ok
    real*8, intent(out) :: s1
    real*8, intent(out) :: s2
    logical, intent(out), optional :: onfar

    integer :: j
    real*8 :: xs(3), ds(3), xa(3), xb(3), sa, sb, t1, t2, umax

    ok = .false.
    s1 = 0d0
    s2 = 0d0
    if (present(onfar)) onfar = .false.

    ! into the box frame, where the box is the unit cube: a direction
    ! transforms like the point, with the inverse of v
    xa = x0 - b%o
    xs = b%vinv(:,1)*xa(1) + b%vinv(:,2)*xa(2) + b%vinv(:,3)*xa(3)
    ds = b%vinv(:,1)*dir(1) + b%vinv(:,2)*dir(2) + b%vinv(:,3)*dir(3)

    ! the cube is the intersection of three slabs, so clip against each
    umax = maxval(abs(ds))
    if (umax < 1d-10) return
    sa = -1d40
    sb = 1d40
    do j = 1, 3
       if (abs(ds(j)) > 1d-10 * umax) then
          t1 = -xs(j) / ds(j)
          t2 = (1d0 - xs(j)) / ds(j)
          sa = max(sa,min(t1,t2))
          sb = min(sb,max(t1,t2))
       else
          if (xs(j) < -clip_eps .or. xs(j) > 1d0+clip_eps) return
       end if
    end do
    if ((sb - sa) * umax < clip_eps) return

    ! the parameter is untouched by the change of frame, so sa and sb are the
    ! bounds in the caller's own parameterization
    ok = .true.
    s1 = sa
    s2 = sb
    if (present(onfar)) then
       xa = xs + sa * ds
       xb = xs + sb * ds
       onfar = any(min(xa,xb) > 1d0 - clip_eps)
    end if

  end subroutine clip_line_box

  !> Clip the plane through x0 with normal nrm to the box b. x0, nrm and b are
  !> in the same frame, which is the caller's to choose: crystallographic
  !> coordinates with the box spanning the displayed cells, or Cartesian with
  !> the box around the drawn scene. nrm is the covector of the plane in that
  !> frame, i.e. the plane is the set of x with nrm . (x - x0) = 0.
  !> ok is false if the plane misses the box, or meets it only at a point or
  !> along an edge, as that draws nothing. Otherwise cen is the centroid of the
  !> intersection polygon, npt its number of vertices and xpt those vertices in
  !> order around it, and onfar says the polygon lies on a far face of the box.
  module subroutine clip_plane_box(b,x0,nrm,ok,cen,npt,xpt,onfar)
    use tools, only: qcksort
    type(elem_box), intent(in) :: b
    real*8, intent(in) :: x0(3)
    real*8, intent(in) :: nrm(3)
    logical, intent(out) :: ok
    real*8, intent(out) :: cen(3)
    integer, intent(out), optional :: npt
    real*8, intent(out), optional :: xpt(3,elem_maxpt)
    logical, intent(out), optional :: onfar

    ! the eight corners of the unit cube, and its twelve edges as pairs of them
    real*8, parameter :: cor(3,8) = reshape((/&
       0d0,0d0,0d0, 1d0,0d0,0d0, 0d0,1d0,0d0, 1d0,1d0,0d0,&
       0d0,0d0,1d0, 1d0,0d0,1d0, 0d0,1d0,1d0, 1d0,1d0,1d0/),shape(cor))
    integer, parameter :: cedge(2,12) = reshape((/&
       1,2, 3,4, 5,6, 7,8,& ! along the first axis
       1,3, 2,4, 5,7, 6,8,& ! along the second
       1,5, 2,6, 3,7, 4,8/),shape(cedge)) ! along the third

    integer :: j, k, np, iperm(elem_maxpt)
    real*8 :: xs(3), ds(3), xa(3), xb(3), fa, fb, epsf
    real*8 :: pts(3,elem_maxpt), cs(3), e1(3), e2(3), vv(3)
    real*8 :: ang(elem_maxpt), dd, dmax

    ok = .false.
    cen = x0
    if (present(npt)) npt = 0
    if (present(xpt)) xpt = 0d0
    if (present(onfar)) onfar = .false.

    ! into the box frame, where the box is the unit cube. The point transforms
    ! with the inverse of v, but the normal is a covector and transforms with
    ! the transpose, so that nrm . (x - x0) keeps its value
    xa = x0 - b%o
    xs = b%vinv(:,1)*xa(1) + b%vinv(:,2)*xa(2) + b%vinv(:,3)*xa(3)
    do j = 1, 3
       ds(j) = dot_product(b%v(:,j),nrm)
    end do

    ! the polygon vertices: where the plane meets the twelve cube edges
    epsf = 1d-6 * max(norm2(ds),1d-10)
    np = 0
    do k = 1, 12
       xa = cor(:,cedge(1,k))
       xb = cor(:,cedge(2,k))
       fa = dot_product(ds,xa - xs)
       fb = dot_product(ds,xb - xs)
       if (abs(fa) < epsf) call addpt(xa)
       if (abs(fb) < epsf) call addpt(xb)
       if (abs(fa) >= epsf .and. abs(fb) >= epsf .and. fa*fb < 0d0) &
          call addpt(xa + fa/(fa-fb) * (xb - xa))
    end do
    if (np < 3) return

    ! centroid of the polygon
    cs = 0d0
    do j = 1, np
       cs = cs + pts(:,j)
    end do
    cs = cs / real(np,8)

    ! the polygon has to extend in two independent directions, or the plane
    ! only grazes the box. e1 and e2 come out of that test as an in-plane
    ! orthonormal pair, and order the vertices below
    dmax = 0d0
    do j = 1, np
       dd = norm2(pts(:,j) - cs)
       if (dd > dmax) then
          dmax = dd
          e1 = pts(:,j) - cs
       end if
    end do
    if (dmax < clip_eps) return
    e1 = e1 / dmax
    dmax = 0d0
    do j = 1, np
       vv = pts(:,j) - cs
       vv = vv - dot_product(vv,e1) * e1
       dd = norm2(vv)
       if (dd > dmax) then
          dmax = dd
          e2 = vv
       end if
    end do
    if (dmax < clip_eps) return
    e2 = e2 / dmax

    ok = .true.
    cen = b%o + b%v(:,1)*cs(1) + b%v(:,2)*cs(2) + b%v(:,3)*cs(3)
    if (present(onfar)) then
       do j = 1, 3
          if (all(pts(j,1:np) > 1d0 - clip_eps)) onfar = .true.
       end do
    end if
    if (present(npt)) npt = np
    if (present(xpt)) then
       ! order the vertices around the polygon by their angle in the (e1,e2)
       ! plane; the map back to the caller's frame is affine, so the cyclic
       ! order it gives is the cyclic order there too
       do j = 1, np
          vv = pts(:,j) - cs
          ang(j) = atan2(dot_product(vv,e2),dot_product(vv,e1))
          iperm(j) = j
       end do
       call qcksort(ang,iperm,1,np)
       do j = 1, np
          vv = pts(:,iperm(j))
          xpt(:,j) = b%o + b%v(:,1)*vv(1) + b%v(:,2)*vv(2) + b%v(:,3)*vv(3)
       end do
    end if

  contains
    !> Append the box-frame point x to the polygon, unless it is already there
    subroutine addpt(x)
      real*8, intent(in) :: x(3)

      integer :: i

      do i = 1, np
         if (all(abs(pts(:,i) - x) < clip_eps)) return
      end do
      if (np >= elem_maxpt) return
      np = np + 1
      pts(:,np) = x

    end subroutine addpt
  end subroutine clip_plane_box

  !> Clip the symmetry element of kind skind through the point x0 to the box b,
  !> dispatching on the kind. dirf is the plane covector or the axis direction
  !> in the frame of b. ok is false if the element does not meet the box in
  !> more than a point or an edge; otherwise xmid is a point of the element
  !> inside the box, so the element can be brought to any other cell by a
  !> translation of at most one cell, and onfar says its intersection with the
  !> box lies on a far face.
  subroutine symelem_in_box(skind,x0,dirf,b,ok,xmid,onfar)
    integer, intent(in) :: skind
    real*8, intent(in) :: x0(3)
    real*8, intent(in) :: dirf(3)
    type(elem_box), intent(in) :: b
    logical, intent(out) :: ok
    real*8, intent(out) :: xmid(3)
    logical, intent(out), optional :: onfar

    real*8 :: s1, s2

    xmid = x0
    if (skind == symop_kind_point) then
       call clip_point_box(b,x0,ok,onfar)
    elseif (skind == symop_kind_axis) then
       call clip_line_box(b,x0,dirf,ok,s1,s2,onfar)
       if (ok) xmid = x0 + 0.5d0 * (s1 + s2) * dirf
    elseif (skind == symop_kind_plane) then
       call clip_plane_box(b,x0,dirf,ok,xmid,onfar=onfar)
    else
       ok = .false.
       if (present(onfar)) onfar = .false.
    end if

  end subroutine symelem_in_box

end submodule symmetry
