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

! Routines for editing molecular and crystal structures
submodule (crystalmod) mols
  implicit none

contains

  !> Identify a fragment in the unit cell. Input: cartesian coords. Output:
  !> A fragment object. This routine is thread-safe.
  module function identify_fragment(c,nat,x0) result(fr)
    use types, only: realloc
    use param, only: icrd_cart
    class(crystal), intent(in) :: c
    integer, intent(in) :: nat
    real*8, intent(in) :: x0(3,nat)
    type(fragment) :: fr

    integer :: id, i, n

    fr%nat = nat
    allocate(fr%at(nat))
    fr%nspc = c%nspc
    allocate(fr%spc(c%nspc))
    fr%spc = c%spc

    n = 0
    do i = 1, nat
       id = c%identify_atom(x0(:,i),icrd_cart)
       if (id > 0) then
          n = n + 1
          fr%at(n)%r = x0(:,i)
          fr%at(n)%x = c%c2x(x0(:,i))
          fr%at(n)%cidx = id
          fr%at(n)%idx = c%atcel(id)%idx
          fr%at(n)%lvec = nint(fr%at(n)%x - c%atcel(id)%x)
          fr%at(n)%is = c%atcel(id)%is
       end if
    end do
    fr%nat = n
    call realloc(fr%at,n)

  end function identify_fragment

  !> Identify a fragment in the unit cell from an external xyz
  !> file. An instance of a fragment object is returned. If any of the
  !> atoms is not correctly identified, return 0.
  module function identify_fragment_from_xyz(c,file,errmsg,ti) result(fr)
    use tools_io, only: fopen_read, string, fclose
    use param, only: bohrtoa, icrd_cart, mlen
    use types, only: realloc

    class(crystal), intent(in) :: c
    character*(*) :: file
    type(thread_info), intent(in), optional :: ti
    character(len=:), allocatable, intent(out) :: errmsg
    type(fragment) :: fr

    integer :: lu, nat
    integer :: id, i
    real*8 :: x0(3)
    character(len=mlen) :: word

    errmsg = "Error reading xyz file"
    lu = fopen_read(file,ti=ti)
    if (lu < 0) goto 999
    read(lu,*,err=999,end=999) nat
    read(lu,*,err=999,end=999)

    fr%nat = nat
    allocate(fr%at(nat))
    fr%nspc = c%nspc
    allocate(fr%spc(fr%nspc))
    fr%spc = c%spc
    do i = 1, nat
       read(lu,*,err=999,end=999) word, x0
       x0 = x0 / bohrtoa - c%molx0
       id = c%identify_atom(x0,icrd_cart)
       if (id == 0) then
          fr%nat = 0
          deallocate(fr%at)
          fr%nspc = 0
          deallocate(fr%spc)
          errmsg = "Could not identify fragment from xyz file: " // trim(file)
          return
       endif
       fr%at(i)%r = x0
       fr%at(i)%x = c%c2x(x0)
       fr%at(i)%cidx = id
       fr%at(i)%idx = c%atcel(id)%idx
       fr%at(i)%lvec = nint(fr%at(i)%x - c%atcel(id)%x)
       fr%at(i)%is = c%atcel(id)%is
    end do
    call fclose(lu)
    call realloc(fr%at,fr%nat)

    errmsg = ""
    return
999 continue
    if (lu > 0) call fclose(lu)

  end function identify_fragment_from_xyz

  !> Using the calculated asterisms for each atom determine the
  !> molecules in the system by finding the connected components of
  !> the graph. This routine also determines whether the crystal is
  !> extended or molecular. Fills c%nmol and c%mol.
  module subroutine fill_molecular_fragments(c)
    use fragmentmod, only: realloc_fragment
    use tools_math, only: lattice_direction
    use tools_io, only: ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    class(crystal), intent(inout) :: c

    integer :: i, j, k, l, jid, newid, newl(3)
    integer :: nat, nlvec, lwork, info
    logical :: found, fdisc
    integer, allocatable :: id(:), lvec(:,:), iord(:), icidx(:)
    logical, allocatable :: ldone(:), used(:)
    real*8, allocatable :: rlvec(:,:), sigma(:), uvec(:,:), vvec(:,:), work(:)
    real*8 :: xcm(3)

    ! initialize
    c%nmol = 0
    if (allocated(c%mol)) deallocate(c%mol)
    if (allocated(c%idatcelmol)) deallocate(c%idatcelmol)
    if (c%ncel == 0) return

    ! checks andallocate
    if (.not.allocated(c%nstar)) &
       call ferror('fill_molecular_fragments','no asterisms found',faterr)
    allocate(used(c%ncel))
    used = .false.
    allocate(c%mol(1),id(10),lvec(3,10),ldone(10))

    ! run over atoms in the unit cell
    do i = 1, c%ncel
       if (used(i)) cycle

       ! increment the fragment counter
       c%nmol = c%nmol + 1
       if (c%nmol > size(c%mol)) then
          call realloc_fragment(c%mol,2*c%nmol)
       end if
       c%mol(c%nmol)%discrete = .true.
       c%mol(c%nmol)%nlvec = 0

       ! initialize the stack with atom i in the seed
       nat = 1
       id(1) = i
       lvec(:,1) = 0d0
       ldone(1) = .false.
       ! run the stack
       do while (.not.all(ldone(1:nat)))
          ! find the next atom that is not done
          do j = 1, nat
             if (.not.ldone(j)) exit
          end do
          ldone(j) = .true.
          jid = id(j)

          ! run over all neighbors of j
          do k = 1, c%nstar(jid)%ncon
             ! id for the atom and lattice vector
             newid = c%nstar(jid)%idcon(k)
             newl = c%nstar(jid)%lcon(:,k) + lvec(:,j)

             ! Is this atom in the fragment already? -> skip it. If it
             ! has a different lattice vector, mark the fragment as
             ! not discrete.
             found = .false.
             do l = 1, nat
                found = (newid == id(l))
                fdisc = all(newl == lvec(:,l))
                if (found) exit
             end do
             if (found) then
                if (.not.fdisc.and..not.c%ismolecule) then
                   if (.not.allocated(c%mol(c%nmol)%lvec)) allocate(c%mol(c%nmol)%lvec(3,1))

                   newl = newl - lvec(:,l)
                   found = .false.
                   do l = 1, c%mol(c%nmol)%nlvec
                      if (all(c%mol(c%nmol)%lvec(:,l) == newl) .or. all(c%mol(c%nmol)%lvec(:,l) == -newl)) then
                         found = .true.
                         exit
                      end if
                   end do
                   if (.not.found) then
                      c%mol(c%nmol)%nlvec = c%mol(c%nmol)%nlvec + 1
                      if (c%mol(c%nmol)%nlvec > size(c%mol(c%nmol)%lvec,2)) then
                         call realloc(c%mol(c%nmol)%lvec,3,2*c%mol(c%nmol)%nlvec)
                      end if
                      c%mol(c%nmol)%lvec(:,c%mol(c%nmol)%nlvec) = newl
                   end if

                   c%mol(c%nmol)%discrete = .false.
                end if
                cycle
             end if

             ! Have we used this atom already?
             if (used(newid)) cycle

             ! Add the atom to the stack and mark it as used.
             nat = nat + 1
             if (nat > size(ldone)) then
                call realloc(id,2*nat)
                call realloc(lvec,3,2*nat)
                call realloc(ldone,2*nat)
             end if
             id(nat) = newid
             lvec(:,nat) = newl
             used(newid) = .true.
             ldone(nat) = .false.
          end do
       end do

       ! add this fragment to the list
       used(i) = .true.
       allocate(c%mol(c%nmol)%spc(c%nspc))
       c%mol(c%nmol)%nspc = c%nspc
       c%mol(c%nmol)%spc = c%spc(1:c%nspc)
       allocate(c%mol(c%nmol)%at(nat))
       c%mol(c%nmol)%nat = nat
       do j = 1, nat
          ! if not a discrete fragment, do not translate
          if (.not.c%mol(c%nmol)%discrete) lvec(:,j) = 0
          c%mol(c%nmol)%at(j)%x = c%atcel(id(j))%x + lvec(:,j)
          c%mol(c%nmol)%at(j)%r = c%x2c(c%mol(c%nmol)%at(j)%x)
          c%mol(c%nmol)%at(j)%cidx = id(j)
          c%mol(c%nmol)%at(j)%idx = c%atcel(id(j))%idx
          c%mol(c%nmol)%at(j)%lvec = lvec(:,j)
          c%mol(c%nmol)%at(j)%is = c%atcel(id(j))%is
       end do
       if (c%mol(c%nmol)%nlvec > 0) &
          call realloc(c%mol(c%nmol)%lvec,3,c%mol(c%nmol)%nlvec)

       ! reorder the fragment in order of increasing cidx
       allocate(iord(nat),icidx(nat))
       do j = 1, nat
          iord(j) = j
          icidx(j) = c%mol(c%nmol)%at(j)%cidx
       end do
       call qcksort(icidx,iord,1,nat)
       c%mol(c%nmol)%at = c%mol(c%nmol)%at(iord)
       deallocate(iord,icidx)
    end do
    call realloc_fragment(c%mol,c%nmol)
    deallocate(used)

    ! translate all fragments to the main cell
    if (.not.c%ismolecule) then
       do i = 1, c%nmol
          xcm = c%mol(i)%cmass()
          newl = floor(c%c2x(xcm))
          do j = 1, c%mol(i)%nat
             c%mol(i)%at(j)%x = c%mol(i)%at(j)%x - newl
             c%mol(i)%at(j)%r = c%x2c(c%mol(i)%at(j)%x)
             c%mol(i)%at(j)%lvec = c%mol(i)%at(j)%lvec - newl
          end do
       end do
    end if

    ! accumulate all the lattice vectors and calculate the vacuum lattice vectors
    nlvec = 0
    c%nlvac = 0
    do i = 1, c%nmol
       nlvec = nlvec + c%mol(i)%nlvec
    end do
    if (nlvec > 0) then
       allocate(rlvec(3,nlvec))
       k = 0
       do i = 1, c%nmol
          do j = 1, c%mol(i)%nlvec
             k = k + 1
             rlvec(:,k) = c%mol(i)%lvec(:,j)
          end do
       end do

       lwork = 5*3 + max(nlvec,3)
       allocate(sigma(3),uvec(3,3),vvec(1,1),work(lwork))
       sigma = 0d0
       call dgesvd('A','N',3,k,rlvec,3,sigma,uvec,3,vvec,1,work,lwork,info)
       if (info /= 0) &
          call ferror("fill_molecular_fragments","dgesvd failed!",faterr)

       c%nlvac = count(abs(sigma) < 1d-12)
       c%lvac = 0d0
       c%lcon = 0d0
       if (c%nlvac == 2) then
          c%lvac(:,1) = lattice_direction(uvec(:,3),.true.)
          c%lvac(:,2) = lattice_direction(uvec(:,2),.true.)
          c%lcon(:,1) = lattice_direction(uvec(:,1),.true.)
       elseif (c%nlvac == 1) then
          c%lvac(:,1) = lattice_direction(uvec(:,3),.true.)
          c%lcon(:,1) = lattice_direction(uvec(:,2),.true.)
          c%lcon(:,2) = lattice_direction(uvec(:,1),.true.)
       end if
       deallocate(rlvec,sigma,uvec,vvec,work)
    end if

    ! fill the mapping between atoms and molecules
    if (allocated(c%idatcelmol)) deallocate(c%idatcelmol)
    allocate(c%idatcelmol(c%ncel))
    c%idatcelmol = 0
    do i = 1, c%nmol
       do j = 1, c%mol(i)%nat
          c%idatcelmol(c%mol(i)%at(j)%cidx) = i
       end do
    end do

  end subroutine fill_molecular_fragments

  !> Calculate whether the crystal is a molecular crystal. If it is,
  !> for each molecule, calculate whether it is contained as a whole
  !> in the asymmetric unit (idxmol=-1) or not (idxmol >= 0). If it
  !> is not, then the molecule can be symmetry-unique idxmol = 0 or
  !> equivalent to the molecule with index idxmol>0. Fills ismol3d
  !> and idxmol.
  module subroutine calculate_molecular_equivalence(c)
    use tools, only: qcksort
    class(crystal), intent(inout) :: c

    integer :: mmax, i, j, nat
    integer, allocatable :: midx(:,:)

    if (c%nmol == 0 .or. c%ncel == 0) return

    ! 3d molecular crystals calculations
    c%ismol3d = all(c%mol(1:c%nmol)%discrete).and..not.c%ismolecule
    if (c%ismol3d) then
       ! assign fractional (-1)/symmetry unique (0)/symmetry equivalent (>0)
       if (allocated(c%idxmol)) deallocate(c%idxmol)
       allocate(c%idxmol(c%nmol))
       c%idxmol = 0

       ! find the maximum number of atoms
       mmax = 0
       do i = 1, c%nmol
          mmax = max(mmax,c%mol(i)%nat)
       end do

       ! sort the atomic indices
       allocate(midx(mmax,c%nmol))
       do i = 1, c%nmol
          nat = c%mol(i)%nat
          do j = 1, nat
             midx(j,i) = c%mol(i)%at(j)%idx
          end do
          call qcksort(midx(1:nat,i))

          ! check if the molecule is fractional
          do j = 1, c%mol(i)%nat-1
             if (midx(j,i) == midx(j+1,i)) then
                c%idxmol(i) = -1
                exit
             end if
          end do
       end do

       ! assign idxmol based on the atomic sequence
       do i = 1, c%nmol
          nat = c%mol(i)%nat
          if (c%idxmol(i) < 0) cycle
          do j = 1, i-1
             if (c%idxmol(j) < 0) cycle
             if (all(midx(1:nat,i) == midx(1:nat,j))) then
                c%idxmol(i) = j
                exit
             end if
          end do
       end do
       deallocate(midx)
    end if

  end subroutine calculate_molecular_equivalence

  !> List atoms in a number of cells around the main cell (nx cells),
  !> possibly with border (doborder).
  module function listatoms_cells(c,nx,doborder) result(fr)
    use types, only: realloc
    class(crystal), intent(in) :: c
    integer, intent(in) :: nx(3)
    logical, intent(in) :: doborder
    type(fragment) :: fr

    real*8, parameter :: rthr = 0.01d0
    real*8, parameter :: rthr1 = 1-rthr

    integer :: ix, iy, iz, i
    logical :: if1

    fr%nspc = c%nspc
    allocate(fr%spc(c%nspc))
    fr%spc(1:fr%nspc) = c%spc(1:c%nspc)
    allocate(fr%at(1))
    fr%nat = 0

    ! All atoms in these cells
    fr%nat = 0
    do ix = 0,nx(1)-1
       do iy = 0,nx(2)-1
          do iz = 0,nx(3)-1
             do i = 1, c%ncel
                fr%nat = fr%nat + 1
                if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                fr%at(fr%nat)%cidx = i
                fr%at(fr%nat)%idx = c%atcel(i)%idx
                fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                fr%at(fr%nat)%is = c%atcel(i)%is
             end do
          end do
       end do
    end do

    ! border: pick atoms
    if (doborder) then
       do ix = -1,nx(1)
          do iy = -1,nx(2)
             do iz = -1,nx(3)
                if (ix > -1 .and. ix < nx(1) .and. iy > -1 .and. iy < nx(2) .and.&
                   iz > -1 .and. iz < nx(3)) cycle
                do i = 1, c%ncel
                   ! border
                   if1 = (ix == -1 .and. c%atcel(i)%x(1)<rthr1 .or. &
                      ix == nx(1) .and. c%atcel(i)%x(1)>rthr .or. &
                      iy == -1 .and. c%atcel(i)%x(2)<rthr1 .or. &
                      iy == nx(2) .and. c%atcel(i)%x(2)>rthr .or. &
                      iz == -1 .and. c%atcel(i)%x(3)<rthr1 .or. &
                      iz == nx(3) .and. c%atcel(i)%x(3)>rthr)
                   if (.not.if1) then
                      fr%nat = fr%nat + 1
                      if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                      fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                      fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                      fr%at(fr%nat)%cidx = i
                      fr%at(fr%nat)%idx = c%atcel(i)%idx
                      fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                      fr%at(fr%nat)%is = c%atcel(i)%is
                      cycle
                   end if
                end do
             end do
          end do
       end do
    end if
    call realloc(fr%at,fr%nat)

  end function listatoms_cells

  !> List atoms inside a sphere of radius rsph and center xsph
  !> (cryst.)  or a cube of side rcub and center xcub (cryst.). Return
  !> the list of atomic positions (Cartesian) in x, the atomic numbers
  !> in z and the number of atoms in nat.
  module function listatoms_sphcub(c,rsph,xsph,rcub,xcub) result(fr)
    use tools_io, only: ferror, faterr
    use types, only: realloc
    class(crystal), intent(in) :: c
    real*8, intent(in), optional :: rsph, xsph(3)
    real*8, intent(in), optional :: rcub, xcub(3)
    type(fragment) :: fr

    integer :: ix, iy, iz, i, nn
    real*8 :: x0(3), d
    logical :: doagain, dosph

    if (.not.(present(rsph).and.present(xsph)).and..not.(present(rcub).and.present(xcub))) &
       call ferror("listatoms_sphcub","Need sphere or cube input",faterr)
    dosph = present(rsph)

    allocate(fr%at(1))
    fr%nat = 0
    allocate(fr%spc(c%nspc))
    fr%nspc = c%nspc
    fr%spc = c%spc

    ! all atoms in a sphere
    doagain = .true.
    nn = -1
    do while(doagain)
       doagain = .false.
       nn = nn + 1
       do ix = -nn, nn
          do iy = -nn, nn
             do iz = -nn, nn
                if (abs(ix)/=nn .and. abs(iy)/=nn .and. abs(iz)/=nn) cycle
                do i = 1, c%ncel
                   if (dosph) then
                      x0 = c%x2c(c%atcel(i)%x + (/ix,iy,iz/) - xsph)
                      if (all(abs(x0) > rsph)) cycle
                      d = norm2(x0)
                      if (d >= rsph) cycle
                   else
                      x0 = c%x2c(c%atcel(i)%x + (/ix,iy,iz/) - xcub)
                      if (any(abs(x0) > rcub)) cycle
                   endif

                   ! add this atom
                   fr%nat = fr%nat + 1
                   if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
                   fr%at(fr%nat)%x = c%atcel(i)%x + (/ix,iy,iz/)
                   fr%at(fr%nat)%r = c%x2c(fr%at(fr%nat)%x)
                   fr%at(fr%nat)%cidx = i
                   fr%at(fr%nat)%idx = c%atcel(i)%idx
                   fr%at(fr%nat)%lvec = (/ix,iy,iz/)
                   fr%at(fr%nat)%is = c%atcel(i)%is
                   doagain = .true.
                end do
             end do
          end do
       end do
    end do
    call realloc(fr%at,fr%nat)

  end function listatoms_sphcub

  !> List all molecules resulting from completing the initial fragment
  !> fri by adding adjacent atoms that are covalently bonded. Return
  !> the number of fragments (nfrag), the fragments themselves (fr),
  !> and whether the fragments are discrete (not connected to copies
  !> of themselves in a different cell).
  module subroutine listmolecules(c,fri,nfrag,fr,isdiscrete)
    use fragmentmod, only: realloc_fragment
    use types, only: realloc
    class(crystal), intent(inout) :: c
    type(fragment), intent(in) :: fri
    integer, intent(out) :: nfrag
    type(fragment), intent(out), allocatable :: fr(:)
    logical, intent(out), allocatable :: isdiscrete(:)

    integer :: i, j, k, l, newid, newl(3), jid
    integer :: nat
    integer, allocatable :: id(:), lvec(:,:)
    logical, allocatable :: ldone(:)
    logical :: found, ldist
    integer :: nseed
    integer, allocatable :: idseed(:), lseed(:,:)
    logical, allocatable :: fseed(:)

    ! unwrap the input fragment
    nseed = fri%nat
    allocate(idseed(nseed),lseed(3,nseed),fseed(nseed))
    do j = 1, nseed
       idseed(j) = fri%at(j)%cidx
       lseed(:,j) = fri%at(j)%lvec
    end do
    fseed = .false.

    ! allocate stuff
    nfrag = 0
    allocate(fr(1),isdiscrete(1),id(10),lvec(3,10),ldone(10))
    isdiscrete = .true.

    do i = 1, nseed
       if (fseed(i)) cycle

       ! initialize the stack with atom i in the seed
       nat = 1
       id(1) = idseed(i)
       lvec(:,1) = lseed(:,i)
       ldone(1) = .false.
       ldist = .true.
       ! run the stack
       do while (.not.all(ldone(1:nat)))
          ! find the next atom that is not done
          do j = 1, nat
             if (.not.ldone(j)) exit
          end do
          ldone(j) = .true.
          jid = id(j)

          ! run over all neighbors of j
          do k = 1, c%nstar(jid)%ncon
             ! id for the atom and lattice vector
             newid = c%nstar(jid)%idcon(k)
             newl = c%nstar(jid)%lcon(:,k) + lvec(:,j)

             ! is this atom in the fragment already? -> skip it
             found = .false.
             do l = 1, nat
                found = (newid == id(l)) .and. all(newl == lvec(:,l))
                if (found) exit
             end do
             if (found) cycle

             ! is this atom in the fragment already with a different
             ! lattice vector?  -> add it to the list but not to the
             ! stack, and mark the fragment as non-discrete
             found = .false.
             do l = 1, nat
                found = (newid == id(l))
                if (found) exit
             end do
             nat = nat + 1
             if (nat > size(ldone)) then
                call realloc(id,2*nat)
                call realloc(lvec,3,2*nat)
                call realloc(ldone,2*nat)
             end if
             id(nat) = newid
             lvec(:,nat) = newl

             if (found) then
                ldone(nat) = .true.
                ldist = .false.
             else
                ! if it wasn't found, then add the atom to the stack
                ldone(nat) = .false.
             end if
          end do
       end do

       ! add this fragment to the list
       fseed(i) = .true.
       nfrag = nfrag + 1
       if (nfrag > size(fr)) then
          call realloc_fragment(fr,2*nfrag)
          call realloc(isdiscrete,2*nfrag)
       end if
       fr(nfrag)%nat = nat
       allocate(fr(nfrag)%at(nat))
       if (allocated(fr(nfrag)%spc)) deallocate(fr(nfrag)%spc)
       fr(nfrag)%nspc = c%nspc
       allocate(fr(nfrag)%spc(c%nspc))
       fr(nfrag)%spc = c%spc
       isdiscrete(nfrag) = ldist
       do j = 1, nat
          ! if not a discrete fragment, do not translate
          if (.not.ldist) lvec(:,j) = 0
          fr(nfrag)%at(j)%x = c%atcel(id(j))%x + lvec(:,j)
          fr(nfrag)%at(j)%r = c%x2c(fr(nfrag)%at(j)%x)
          fr(nfrag)%at(j)%cidx = id(j)
          fr(nfrag)%at(j)%idx = c%atcel(id(j))%idx
          fr(nfrag)%at(j)%lvec = lvec(:,j)
          fr(nfrag)%at(j)%is = c%atcel(id(j))%is
       end do

       ! run over all atoms in the new fragment and mark those atoms in the seed
       do j = 1, nat
          do k = 1, nseed
             if (fseed(k)) cycle
             if (id(j) == idseed(k) .and. all(lvec(:,j) == lseed(:,k))) &
                fseed(k) = .true.
          end do
       end do
    end do
    call realloc_fragment(fr,nfrag)
    call realloc(isdiscrete,nfrag)

  end subroutine listmolecules

end submodule mols
