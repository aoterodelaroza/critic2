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
    use param, only: icrd_cart
    class(crystal), intent(inout) :: c
    integer, intent(in) :: nat
    real*8, intent(in) :: x0(3,nat)
    type(fragment) :: fr

    integer :: id, i, n
    real*8, allocatable :: fx(:,:)
    integer, allocatable :: fis(:), fidx(:), fcidx(:)

    ! gather the atoms that can be identified in the cell
    allocate(fx(3,nat),fis(nat),fidx(nat),fcidx(nat))
    n = 0
    do i = 1, nat
       id = c%identify_atom(x0(:,i),icrd_cart)
       if (id > 0) then
          n = n + 1
          fx(:,n) = x0(:,i)
          fcidx(n) = id
          fidx(n) = c%atcel(id)%idx
          fis(n) = c%atcel(id)%is
       end if
    end do

    ! build the fragment and assign the lattice vectors to the cell atoms
    call fr%build(c%nspc,c%spc,n,fx(:,1:n),icrd_cart,fis(1:n),fidx(1:n),fcidx(1:n),c%m_x2c)
    do i = 1, fr%nat
       fr%at(i)%lvec = nint(fr%at(i)%x - c%atcel(fr%at(i)%cidx)%x)
    end do
    deallocate(fx,fis,fidx,fcidx)

  end function identify_fragment

  !> Identify a fragment in the unit cell from an external xyz
  !> file. An instance of a fragment object is returned. If any of the
  !> atoms is not correctly identified, return 0.
  module function identify_fragment_from_xyz(c,file,errmsg,ti) result(fr)
    use tools_io, only: fopen_read, string, fclose
    use param, only: bohrtoa, icrd_cart, mlen

    class(crystal), intent(inout) :: c
    character*(*) :: file
    type(thread_info), intent(in), optional :: ti
    character(len=:), allocatable, intent(out) :: errmsg
    type(fragment) :: fr

    integer :: lu, nat
    integer :: id, i
    real*8 :: x0(3)
    real*8, allocatable :: fx(:,:)
    integer, allocatable :: fis(:), fidx(:), fcidx(:)
    character(len=mlen) :: word

    ! leave an empty fragment in case of error
    fr%nat = 0
    fr%nspc = 0

    errmsg = "Error reading xyz file"
    lu = fopen_read(file,ti=ti)
    if (lu < 0) goto 999
    read(lu,*,err=999,end=999) nat
    read(lu,*,err=999,end=999)

    ! gather and identify all the atoms in the xyz file
    allocate(fx(3,nat),fis(nat),fidx(nat),fcidx(nat))
    do i = 1, nat
       read(lu,*,err=999,end=999) word, x0
       x0 = x0 / bohrtoa - c%molx0
       id = c%identify_atom(x0,icrd_cart)
       if (id == 0) then
          errmsg = "Could not identify fragment from xyz file: " // trim(file)
          call fclose(lu)
          return
       endif
       fx(:,i) = x0
       fcidx(i) = id
       fidx(i) = c%atcel(id)%idx
       fis(i) = c%atcel(id)%is
    end do
    call fclose(lu)

    ! build the fragment and assign the lattice vectors to the cell atoms
    call fr%build(c%nspc,c%spc,nat,fx,icrd_cart,fis,fidx,fcidx,c%m_x2c)
    do i = 1, fr%nat
       fr%at(i)%lvec = nint(fr%at(i)%x - c%atcel(fr%at(i)%cidx)%x)
    end do

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
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c

    integer :: i, j, k, nat
    integer :: nlvec, lwork, info
    integer, allocatable :: lvec(:,:), lmol(:,:), nlmol(:)
    integer, allocatable :: fis(:), fidx(:), fcidx(:), flvec(:,:)
    logical, allocatable :: isdiscrete(:)
    real*8, allocatable :: rlvec(:,:), sigma(:), uvec(:,:), vvec(:,:), work(:)
    real*8, allocatable :: fx(:,:)

    ! initialize
    c%nmol = 0
    if (allocated(c%mol)) deallocate(c%mol)
    if (allocated(c%idatcelmol)) deallocate(c%idatcelmol)
    if (c%ncel == 0) return

    ! checks and allocate
    if (.not.allocated(c%nstar)) &
       call ferror('fill_molecular_fragments','no asterisms found',faterr)
    allocate(c%idatcelmol(2,c%ncel),lvec(3,c%ncel),isdiscrete(20),lmol(3,20),nlmol(20))
    nlmol = 0
    c%idatcelmol = 0
    isdiscrete = .true.

    ! depth-first search for the connected components
    do i = 1, c%ncel
       if (c%idatcelmol(1,i) > 0) cycle

       c%nmol = c%nmol + 1
       if (c%nmol > size(isdiscrete,1)) then
          call realloc(nlmol,2*c%nmol)
          call realloc(lmol,3,2*c%nmol)
          call realloc(isdiscrete,2*c%nmol)
          isdiscrete(c%nmol:) = .true.
          nlmol(c%nmol:) = 0
       end if
       call explore_node(i,c%nmol,(/0,0,0/))
    end do

    ! build the fragment for each molecule by collecting its atoms from
    ! the complete atom list (work arrays sized to fit any molecule)
    allocate(c%mol(c%nmol))
    allocate(fx(3,c%ncel),fis(c%ncel),fidx(c%ncel),fcidx(c%ncel),flvec(3,c%ncel))
    do i = 1, c%nmol
       nat = 0
       do k = 1, c%ncel
          if (c%idatcelmol(1,k) /= i) cycle
          nat = nat + 1
          c%idatcelmol(2,k) = nat
          if (isdiscrete(i)) then
             fx(:,nat) = c%atcel(k)%x + lvec(:,k)
             flvec(:,nat) = lvec(:,k)
          else
             fx(:,nat) = c%atcel(k)%x
             flvec(:,nat) = 0
          end if
          fis(nat) = c%atcel(k)%is
          fidx(nat) = c%atcel(k)%idx
          fcidx(nat) = k
       end do
       call c%mol(i)%build(c%nspc,c%spc,nat,fx(:,1:nat),icrd_crys,fis(1:nat),&
          fidx(1:nat),fcidx(1:nat),c%m_x2c,flvec(:,1:nat),&
          isdiscrete(i),nlmol(i),lmol(:,1:nlmol(i)))
    end do
    deallocate(nlmol,lmol,isdiscrete,lvec,fx,fis,fidx,fcidx,flvec)

    ! translate all fragments to the main cell
    if (.not.c%ismolecule) then
       do i = 1, c%nmol
          call c%mol(i)%translate_to_main_cell()
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

  contains
    subroutine explore_node(istart,nmol,lvecstart)
      integer, intent(in) :: istart
      integer, intent(in) :: nmol
      integer, intent(in) :: lvecstart(3)

      integer :: k, newid, lnew(3), m, i, sp
      integer :: lveci(3)
      logical :: found
      ! explicit DFS stack (replaces recursion to avoid stack overflow
      ! in large molecules); stores node, neighbor cursor, lattice
      ! vector
      integer, allocatable :: stk_i(:), stk_k(:), stk_lvec(:,:)

      allocate(stk_i(64),stk_k(64),stk_lvec(3,64))
      sp = 1
      stk_i(1) = istart
      stk_k(1) = 0
      stk_lvec(:,1) = lvecstart
      c%idatcelmol(1,istart) = nmol
      lvec(:,istart) = lvecstart

      do while (sp > 0)
         i = stk_i(sp)
         lveci = stk_lvec(:,sp)
         stk_k(sp) = stk_k(sp) + 1
         k = stk_k(sp)
         if (k > c%nstar(i)%ncon) then
            sp = sp - 1    ! exhausted this node, pop
            cycle
         end if

         newid = c%nstar(i)%idcon(k)
         if (c%idatcelmol(1,newid) == 0) then
            ! discover newid, then descend into it (push)
            c%idatcelmol(1,newid) = nmol
            lvec(:,newid) = lveci + c%nstar(i)%lcon(:,k)
            sp = sp + 1
            if (sp > size(stk_i,1)) then
               call realloc(stk_i,2*sp)
               call realloc(stk_k,2*sp)
               call realloc(stk_lvec,3,2*sp)
            end if
            stk_i(sp) = newid
            stk_k(sp) = 0
            stk_lvec(:,sp) = lvec(:,newid)
         else
            lnew = abs(lveci + c%nstar(i)%lcon(:,k) - lvec(:,newid))
            if (any(lnew /= 0)) then
               isdiscrete(nmol) = .false.
               found = .false.
               do m = 1, nlmol(nmol)
                  if (all(lmol(:,m) == lnew)) then
                     found = .true.
                     exit
                  end if
               end do
               if (.not.found) then
                  nlmol(nmol) = nlmol(nmol) + 1
                  if (size(lmol,2) < nlmol(nmol)) call realloc(lmol,3,2*nlmol(nmol))
                  lmol(:,nlmol(nmol)) = lnew
               end if
            end if
         end if
      end do

      deallocate(stk_i,stk_k,stk_lvec)
    end subroutine explore_node

  end subroutine fill_molecular_fragments

  !> Calculate whether the crystal is a molecular crystal. If it is,
  !> for each molecule, calculate whether it is symmetry-unique
  !> (idxmol = 0) or equivalent to molecule with index idxmol > 0.
  !> Fills ismol3d and idxmol in c.
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
       end do

       ! assign idxmol based on the atomic sequence
       do i = 1, c%nmol
          nat = c%mol(i)%nat
          do j = 1, i-1
             if (all(midx(1:nat,i) == midx(1:nat,j))) then
                c%idxmol(i) = j
                exit
             end if
          end do
       end do
       deallocate(midx)
    end if

  end subroutine calculate_molecular_equivalence

  !> Calculate the periodicity of this crystal and set c%iperiod.
  module subroutine calculate_periodicity(c)
    class(crystal), intent(inout) :: c

    integer :: nvac

    if (c%ismolecule) then
       if (c%nmol > 1) then
          c%iperiod = iperiod_mol_cluster
       else
          c%iperiod = iperiod_mol_single
       end if
    else
       nvac = count(c%vaclength > iperiod_vacthr)
       if (nvac == 3) then
          c%iperiod = iperiod_0d
       elseif (nvac == 2) then
          c%iperiod = iperiod_1d
       elseif (nvac == 1) then
          c%iperiod = iperiod_2d
       else
          if (c%ismol3d .or. c%nlvac == 3) then
             c%iperiod = iperiod_3d_molecular
          elseif (c%nlvac == 2) then
             c%iperiod = iperiod_3d_chain
          elseif (c%nlvac == 1) then
             c%iperiod = iperiod_3d_layered
          else
             c%iperiod = iperiod_3d_crystal
          end if
       end if
    end if

  end subroutine calculate_periodicity

  !> List atoms in a number of cells around the main cell (nx cells),
  !> possibly with border (doborder).
  module function listatoms_cells(c,nx,doborder) result(fr)
    use param, only: icrd_crys
    class(crystal), intent(in) :: c
    integer, intent(in) :: nx(3)
    logical, intent(in) :: doborder
    type(fragment) :: fr

    real*8, parameter :: rthr = 0.01d0
    real*8, parameter :: rthr1 = 1-rthr

    integer :: ix, iy, iz, i, n, nmax
    logical :: if1
    real*8, allocatable :: fx(:,:)
    integer, allocatable :: fis(:), fidx(:), fcidx(:), flvec(:,:)

    ! work arrays large enough to fit all the atoms (with border)
    if (doborder) then
       nmax = (nx(1)+2)*(nx(2)+2)*(nx(3)+2)*c%ncel
    else
       nmax = nx(1)*nx(2)*nx(3)*c%ncel
    end if
    allocate(fx(3,nmax),fis(nmax),fidx(nmax),fcidx(nmax),flvec(3,nmax))

    ! All atoms in these cells
    n = 0
    do ix = 0,nx(1)-1
       do iy = 0,nx(2)-1
          do iz = 0,nx(3)-1
             do i = 1, c%ncel
                n = n + 1
                fx(:,n) = c%atcel(i)%x + (/ix,iy,iz/)
                fcidx(n) = i
                fidx(n) = c%atcel(i)%idx
                flvec(:,n) = (/ix,iy,iz/)
                fis(n) = c%atcel(i)%is
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
                      n = n + 1
                      fx(:,n) = c%atcel(i)%x + (/ix,iy,iz/)
                      fcidx(n) = i
                      fidx(n) = c%atcel(i)%idx
                      flvec(:,n) = (/ix,iy,iz/)
                      fis(n) = c%atcel(i)%is
                      cycle
                   end if
                end do
             end do
          end do
       end do
    end if

    call fr%build(c%nspc,c%spc,n,fx(:,1:n),icrd_crys,fis(1:n),fidx(1:n),fcidx(1:n),c%m_x2c,flvec(:,1:n))
    deallocate(fx,fis,fidx,fcidx,flvec)

  end function listatoms_cells

  !> List atoms inside a sphere of radius rsph and center xsph
  !> (cryst.)  or a cube of side rcub and center xcub (cryst.). Return
  !> the list of atomic positions (Cartesian) in x, the atomic numbers
  !> in z and the number of atoms in nat.
  module function listatoms_sphcub(c,rsph,xsph,rcub,xcub) result(fr)
    use tools_io, only: ferror, faterr
    use types, only: realloc
    use param, only: icrd_crys
    class(crystal), intent(in) :: c
    real*8, intent(in), optional :: rsph, xsph(3)
    real*8, intent(in), optional :: rcub, xcub(3)
    type(fragment) :: fr

    integer :: ix, iy, iz, i, nn, n
    real*8 :: x0(3), d
    logical :: doagain, dosph
    real*8, allocatable :: fx(:,:)
    integer, allocatable :: fis(:), fidx(:), fcidx(:), flvec(:,:)

    if (.not.(present(rsph).and.present(xsph)).and..not.(present(rcub).and.present(xcub))) &
       call ferror("listatoms_sphcub","Need sphere or cube input",faterr)
    dosph = present(rsph)

    n = 0
    allocate(fx(3,10),fis(10),fidx(10),fcidx(10),flvec(3,10))

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
                   n = n + 1
                   if (n > size(fis)) then
                      call realloc(fx,3,2*n)
                      call realloc(fis,2*n)
                      call realloc(fidx,2*n)
                      call realloc(fcidx,2*n)
                      call realloc(flvec,3,2*n)
                   end if
                   fx(:,n) = c%atcel(i)%x + (/ix,iy,iz/)
                   fcidx(n) = i
                   fidx(n) = c%atcel(i)%idx
                   flvec(:,n) = (/ix,iy,iz/)
                   fis(n) = c%atcel(i)%is
                   doagain = .true.
                end do
             end do
          end do
       end do
    end do

    call fr%build(c%nspc,c%spc,n,fx(:,1:n),icrd_crys,fis(1:n),fidx(1:n),fcidx(1:n),c%m_x2c,flvec(:,1:n))
    deallocate(fx,fis,fidx,fcidx,flvec)

  end function listatoms_sphcub

  !> List all molecules resulting from completing the initial fragment
  !> fri by adding adjacent atoms that are covalently bonded. Return
  !> the number of fragments (nfrag), the fragments themselves (fr),
  !> and whether the fragments are discrete (not connected to copies
  !> of themselves in a different cell).
  module subroutine listmolecules(c,fri,nfrag,fr,isdiscrete)
    use fragmentmod, only: realloc_fragment
    use types, only: realloc
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    type(fragment), intent(in) :: fri
    integer, intent(out) :: nfrag
    type(fragment), intent(out), allocatable :: fr(:)
    logical, intent(out), allocatable :: isdiscrete(:)

    integer :: i, j, k, l, newid, newl(3), jid
    integer :: nat
    integer, allocatable :: id(:), lvec(:,:)
    integer, allocatable :: fis(:), fidx(:), fcidx(:)
    real*8, allocatable :: fx(:,:)
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
       isdiscrete(nfrag) = ldist
       allocate(fx(3,nat),fis(nat),fidx(nat),fcidx(nat))
       do j = 1, nat
          ! if not a discrete fragment, do not translate
          if (.not.ldist) lvec(:,j) = 0
          fx(:,j) = c%atcel(id(j))%x + lvec(:,j)
          fcidx(j) = id(j)
          fidx(j) = c%atcel(id(j))%idx
          fis(j) = c%atcel(id(j))%is
       end do
       call fr(nfrag)%build(c%nspc,c%spc,nat,fx,icrd_crys,fis,fidx,fcidx,c%m_x2c,lvec(:,1:nat))
       deallocate(fx,fis,fidx,fcidx)

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
