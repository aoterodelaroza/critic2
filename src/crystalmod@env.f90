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

! Routines for atomic environments and distance calculations.
submodule (crystalmod) env
  implicit none

  !xx! private procedures
  ! subroutine make_block_shell(c,n,nb,idb,dmax)
  ! subroutine blocks_inscribed_sphere(c,x,rmax,i0,i1)

contains

  !> Fill the environment variables (nblock, iblock0, blockrmax) as
  !> well as the %inext field in the complete atom list
  !> atoms. Requires the reduced lattice vectors as well as the
  !> complete atom list.
  module subroutine build_env(c)
    use tools_math, only: det3, cross
    class(crystal), intent(inout) :: c

    integer :: i, ix(3)
    real*8 :: xred(3,3), xlen(3), xr(3)

    real*8, parameter :: lenmax = 5d0

    ! calculate the number of cells in each direction and the lattice vectors
    do i = 1, 3
       c%nblock(i) = ceiling(norm2(c%m_xr2c(:,i))/lenmax)
       xred(:,i) = c%m_xr2c(:,i) / c%nblock(i)
    end do

    ! block volume and block radius
    c%blockomega = det3(xred)
    c%blockrmax = 0.5d0 * c%blockomega / max(norm2(cross(xred(:,1),xred(:,2))),&
       norm2(cross(xred(:,1),xred(:,3))),norm2(cross(xred(:,2),xred(:,3))))
    do i = 1, 3
       xlen(i) = norm2(xred(:,i))
    end do
    c%blockcv(1) = norm2(cross(xred(:,2),xred(:,3)))
    c%blockcv(2) = norm2(cross(xred(:,1),xred(:,3)))
    c%blockcv(3) = norm2(cross(xred(:,1),xred(:,2)))

    ! assign atomic indices
    if (allocated(c%iblock0)) deallocate(c%iblock0)
    allocate(c%iblock0(0:c%nblock(1)-1,0:c%nblock(2)-1,0:c%nblock(3)-1))
    c%iblock0 = 0
    c%atcel(:)%inext = 0
    do i = 1, c%ncel
       xr = c%x2xr(c%atcel(i)%x)
       xr = xr - floor(xr)
       ix = min(max(floor(xr * c%nblock),0),c%nblock-1)
       c%atcel(i)%inext = c%iblock0(ix(1),ix(2),ix(3))
       c%iblock0(ix(1),ix(2),ix(3)) = i
    end do

  end subroutine build_env

  !> Given the point xp (in icrd coordinates), calculate the list of
  !> nearest atoms. If sorted is true, the list is sorted by distance
  !> and shell. One (and only one) of these criteria must be given:
  !> - up2d = list up to a distance up2d.
  !> - up2dsp = between a minimum and a maximum species-dependent distance.
  !> - up2dcidx = up to a distance to each atom in the complete list.
  !> - up2sh = up to a number of shells.
  !> - up2n = up to a number of atoms.
  !> If up2n, up2sh, or ishell0 are present, the atom list is sorted by distance on output
  !> regardless of the value of sorted.
  !>
  !> Optional input:
  !> - nid0 = consider only atoms with index nid0 from the nneq list.
  !> - id0 = consider only atoms with index id0 from the complete list.
  !> - iz0 = consider only atoms with atomic number iz0.
  !> - ispc0 = consider only species with species ID ispc0.
  !> - nozero = disregard zero-distance atoms.
  !>
  !> Output:
  !> - nat = the list contains nat atoms.
  !>
  !> Optional output:
  !> - eid(nid) = atom ID from the complete atom list.
  !> - dist(nid) = distance from xp.
  !> - lvec(3,nid) = lattice vectors for the atoms.
  !> - ishell0(nid) = shell ID.
  !>
  !> This routine is thread-safe.
  module subroutine list_near_atoms(c,xp,icrd,sorted,nat,eid,dist,lvec,ishell0,up2d,&
     up2dsp,up2dcidx,up2sh,up2n,nid0,id0,iz0,ispc0,nozero)
    use tools, only: mergesort
    use tools_io, only: ferror, faterr
    use tools_math, only: cross, det3
    use types, only: realloc
    use param, only: icrd_cart, icrd_crys
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    logical, intent(in) :: sorted
    integer, intent(out) :: nat
    integer, allocatable, intent(inout), optional :: eid(:)
    real*8, allocatable, intent(inout), optional :: dist(:)
    integer, allocatable, intent(inout), optional :: lvec(:,:)
    integer, allocatable, intent(inout), optional :: ishell0(:)
    real*8, intent(in), optional :: up2d
    real*8, intent(in), optional :: up2dsp(1:c%nspc,2)
    real*8, intent(in), optional :: up2dcidx(1:c%ncel)
    integer, intent(in), optional :: up2sh
    integer, intent(in), optional :: up2n
    integer, intent(in), optional :: nid0
    integer, intent(in), optional :: id0
    integer, intent(in), optional :: iz0
    integer, intent(in), optional :: ispc0
    logical, intent(in), optional :: nozero

    logical :: ok, nozero_, docycle, sorted_
    real*8 :: x(3), xorigc(3), dmax, dd, lvecx(3), xdif(3), dsqrt
    integer :: i, j, k, ix(3), i0(3), i1(3), i0prev(3), i1prev(3), idx, nn
    integer :: ib(3), nsafe, up2n_
    integer, allocatable :: at_id(:), at_lvec(:,:)
    real*8, allocatable :: at_dist(:), rshel(:)
    integer, allocatable :: iaux(:), iord(:)
    integer :: nshel
    real*8 :: up2d_2, rcur
    real*8,allocatable :: up2dsp_2(:,:), up2dcidx_2(:)

    real*8, parameter :: eps = 1d-20
    integer, parameter :: ixmol_cut = 10 ! cutoff for switching to no-block molecule dist search
    real*8, parameter :: shell_eps = 1d-5

    ! process input
    nozero_ = .false.
    sorted_ = sorted .or. present(up2n) .or. present(up2sh) .or. present(ishell0)
    dmax = huge(1d0)
    if (present(nozero)) nozero_ = nozero
    if (present(up2d)) then
       dmax = up2d
       up2d_2 = up2d * up2d
    elseif (present(up2dsp)) then
       dmax = maxval(up2dsp)
       up2dsp_2 = up2dsp * up2dsp
    elseif (present(up2dcidx)) then
       dmax = maxval(up2dcidx)
       up2dcidx_2 = up2dcidx * up2dcidx
    elseif (present(up2n)) then
       !
    elseif (present(up2sh)) then
       !
    else
       call ferror("list_near_atoms","must give one of up2d, up2dsp, up2dcidx, up2sh, or up2n",faterr)
    end if

    ! locate the block containing the point
    if (icrd == icrd_crys) then
       x = c%x2xr(xp)
    elseif (icrd == icrd_cart) then
       x = c%c2xr(xp)
    else
       x = xp
    end if
    xorigc = c%xr2c(x)
    ix = floor(x * c%nblock)

    ! allocate space for atoms
    nat = 0
    nn = 20
    if (present(up2n)) nn = min(up2n,20)
    allocate(at_id(nn),at_dist(nn),at_lvec(3,nn))
    at_id = 0
    at_dist = 0d0
    at_lvec = 0

    ! run the search
    if (present(up2sh).or.present(up2n)) then
       ! search a number of shells or a number of atoms (up2n, up2sh) !!!

       up2n_ = 0
       if (present(up2sh)) then
          ! initialize shell list
          nshel = 0
          allocate(rshel(up2sh))
       else
          up2n_ = up2n
          if (c%ismolecule) up2n_ = min(up2n,c%ncel)
       end if

       if (c%ismolecule .and. maxval(min(abs(ix), abs(ix-c%nblock))) > ixmol_cut) then
          ! this point is too far away from the molecule: calculate
          ! the distance to every atom in the molecule
          nat = c%ncel
          nsafe = nat
          dmax = huge(1d0)
          call realloc(at_id,nat)
          call realloc(at_dist,nat)
          call realloc(at_lvec,3,nat)
          do i = 1, c%ncel
             dsqrt = norm2(c%atcel(i)%r - xorigc)
             at_id(i) = i
             at_dist(i) = dsqrt
             at_lvec(:,i) = 0
             if (present(up2sh)) &
                call add_shell_to_output_list()
          end do
       else
          ! start with zero current distance
          rcur = 0d0
          i0prev = huge(1)
          i1prev = -huge(1)+1

          ! Run over shells of nearby blocks and find atoms until we
          ! have at least the given number of shells
          do while(.true.)
             ! get the new region
             rcur = rcur + 1d-5
             call blocks_inscribed_sphere(c,x,rcur,i0,i1)

             ! add the atoms in the selected blocks in the region not already explored
             do i = i0(1), i1(1)
                do j = i0(2), i1(2)
                   do k = i0(3), i1(3)
                      if (i>=i0prev(1) .and. i<=i1prev(1) .and. j>=i0prev(2) .and. j<=i1prev(2) .and.&
                          k>=i0prev(3) .and. k<=i1prev(3)) cycle
                      if (c%ismolecule) then
                         if (i < 0 .or. i >= c%nblock(1) .or. j < 0 .or. j >= c%nblock(2) .or.&
                            k < 0 .or. k >= c%nblock(3)) cycle
                         ib = (/i, j, k/)
                         lvecx = 0
                         xdif = -xorigc
                      else
                         ib = (/modulo(i,c%nblock(1)), modulo(j,c%nblock(2)), modulo(k,c%nblock(3))/)
                         lvecx = matmul(c%m_xr2c,real(((/i,j,k/) - ib) / c%nblock,8))
                         xdif = lvecx - xorigc
                      end if

                      ! run over atoms in this block
                      idx = c%iblock0(ib(1),ib(2),ib(3))
                      do while (idx /= 0)
                         ! apply filters
                         docycle = .false.
                         if (present(nid0)) &
                            docycle = (c%atcel(idx)%idx /= nid0)
                         if (present(ispc0)) &
                            docycle = (c%atcel(idx)%is /= ispc0)
                         if (present(id0)) &
                            docycle = (idx /= id0)
                         if (present(iz0)) &
                            docycle = (c%spc(c%atcel(idx)%is)%z /= iz0)

                         if (.not.docycle) then
                            ! calculate distance
                            if (c%ismolecule) then
                               dd = dot_product(c%atcel(idx)%r  - xorigc,c%atcel(idx)%r - xorigc)
                            else
                               dd = dot_product(c%atcel(idx)%rxc + xdif,c%atcel(idx)%rxc + xdif)
                            end if

                            ! check if we should add the atom to the list
                            ok = .true.
                            if (nozero_ .and. dd < eps) then
                               ok = .false.
                               if (c%ismolecule) up2n_ = min(up2n,c%ncel-1)
                            end if
                            if (ok) then
                               dsqrt = sqrt(dd)
                               call add_atom_to_output_list()
                               if (present(up2sh)) &
                                  call add_shell_to_output_list()
                            end if
                         end if

                         ! next atom
                         idx = c%atcel(idx)%inext
                      end do ! while idx
                   end do
                end do
             end do

             if (present(up2sh)) then
                ! exit if we have enough shells
                nsafe = count(rshel(1:nshel) < rcur)
                if (nsafe >= up2sh) exit
             else
                ! check if we have enough atoms
                nsafe = count(at_dist(1:nat) < rcur)
                if (nsafe >= up2n_) exit
             end if

             ! exit if we have examined all atoms in the molecule
             if (c%ismolecule .and. nsafe == c%ncel) then
                up2n_ = min(up2n_,nat)
                exit
             end if

             ! next shell
             i0prev = i0
             i1prev = i1
          end do ! while loop
       end if

       if (present(up2n) .and. up2n_ == 1) then
          ! bypass for efficiency: handle the case when up2n = 1 (nearest atom)
          idx = minloc(at_dist(1:nat),1)
          nat = 1
          at_dist(1) = at_dist(idx)
          at_id(1) = at_id(idx)
          at_lvec(:,1) = at_lvec(:,idx)
          sorted_ = .false.
       else
          if (present(up2sh)) then
             ! re-calculate the dmax based on the shell radii
             allocate(iord(nshel))
             do i = 1, nshel
                iord(i) = i
             end do
             call mergesort(rshel,iord,1,nshel)
             if (c%ismolecule .and. nshel < up2sh) then
                dmax = rshel(iord(nshel)) + shell_eps
             else
                dmax = rshel(iord(up2sh)) + shell_eps
             end if
             deallocate(iord)
          end if

          ! filter out the unneeded atoms
          nsafe = count(at_dist(1:nat) <= rcur)
          allocate(iaux(nsafe))
          nsafe = 0
          do i = 1, nat
             if (at_dist(i) <= rcur) then
                nsafe = nsafe + 1
                iaux(nsafe) = i
             end if
          end do
          nat = nsafe
          at_id(1:nsafe) = at_id(iaux)
          at_dist(1:nsafe) = at_dist(iaux)
          at_lvec(:,1:nsafe) = at_lvec(:,iaux)
          call realloc(at_id,nsafe)
          call realloc(at_dist,nsafe)
          call realloc(at_lvec,3,nsafe)
          deallocate(iaux)
       end if
    else
       ! search atoms up to a given distance (up2d*) !!!

       ! calculate the region of blocks in each direction that fits a
       ! sphere with radius dmax
       call blocks_inscribed_sphere(c,x,dmax,i0,i1)

       ! run over the atoms and compile the list
       do i = i0(1), i1(1)
          do j = i0(2), i1(2)
             do k = i0(3), i1(3)
                if (c%ismolecule) then
                   if (i < 0 .or. i >= c%nblock(1) .or. j < 0 .or. j >= c%nblock(2) .or.&
                      k < 0 .or. k >= c%nblock(3)) cycle
                   ib = (/i, j, k/)
                   lvecx = 0
                   xdif = -xorigc
                else
                   ib = (/modulo(i,c%nblock(1)), modulo(j,c%nblock(2)), modulo(k,c%nblock(3))/)
                   lvecx = matmul(c%m_xr2c,real(((/i,j,k/) - ib) / c%nblock,8))
                   xdif = lvecx - xorigc
                end if

                ! run over atoms in this block
                idx = c%iblock0(ib(1),ib(2),ib(3))
                do while (idx /= 0)
                   ! apply filters
                   docycle = .false.
                   if (present(nid0)) &
                      docycle = (c%atcel(idx)%idx /= nid0)
                   if (present(ispc0)) &
                      docycle = (c%atcel(idx)%is /= ispc0)
                   if (present(id0)) &
                      docycle = (idx /= id0)
                   if (present(iz0)) &
                      docycle = (c%spc(c%atcel(idx)%is)%z /= iz0)

                   if (.not.docycle) then
                      ! calculate distance
                      if (c%ismolecule) then
                         x = c%atcel(idx)%r
                         dd = dot_product(x - xorigc,x - xorigc)
                      else
                         dd = dot_product(c%atcel(idx)%rxc + xdif,c%atcel(idx)%rxc + xdif)
                      end if

                      ! check if we should add the atom to the list
                      ok = .true.
                      if (nozero_ .and. dd < eps) then
                         ok = .false.
                      elseif (present(up2d)) then
                         ok = (dd <= up2d_2)
                      elseif (present(up2dsp)) then
                         ok = (dd >= up2dsp_2(c%atcel(idx)%is,1) .and. dd <= up2dsp_2(c%atcel(idx)%is,2))
                      elseif (present(up2dcidx)) then
                         ok = (dd <= up2dcidx_2(idx))
                      end if
                      if (ok) then
                         dsqrt = sqrt(dd)
                         call add_atom_to_output_list()
                      end if
                   end if

                   ! next atom
                   idx = c%atcel(idx)%inext
                end do ! while idx
             end do ! k = i0(3), i1(3)
          end do ! j = i0(2), i1(2)
       end do ! i = i0(1), i1(1)
    end if

    ! sort if necessary
    if (sorted_ .and. nat > 1) then
       ! permutation vector
       allocate(iord(nat))
       do i = 1, nat
          iord(i) = i
       end do

       ! sort by distance
       call mergesort(at_dist,iord,1,nat)

       ! reorder and clean up
       at_id = at_id(iord)
       at_dist = at_dist(iord)
       at_lvec = at_lvec(:,iord)
       deallocate(iord)
    end if

    ! assign shells if requested
    if (present(ishell0)) then
       if (allocated(ishell0)) deallocate(ishell0)
       allocate(ishell0(nat))
       dd = -1d0
       k = 0
       do i = 1, nat
          if (abs(at_dist(i) - dd) > shell_eps) then
             k = k + 1
             dd = at_dist(i)
          end if
          ishell0(i) = k
       end do
    end if

    ! reduce the list if up2n (always sorted)
    if (present(up2n)) nat = up2n_

    ! assign optional outputs
    if (present(eid)) eid = at_id(1:nat)
    if (present(dist)) dist = at_dist(1:nat)
    if (present(lvec)) lvec = at_lvec(:,1:nat)

  contains
    ! add current atom (idx) to the output list
    subroutine add_atom_to_output_list()
      nat = nat + 1
      if (nat > size(at_id,1)) then
         call realloc(at_id,2*nat)
         call realloc(at_dist,2*nat)
         call realloc(at_lvec,3,2*nat)
      end if
      at_id(nat) = idx
      at_dist(nat) = dsqrt
      if (c%ismolecule) then
         at_lvec(:,nat) = 0
      else
         at_lvec(:,nat) = nint(matmul(c%m_c2x,c%atcel(idx)%rxc + lvecx) - c%atcel(idx)%x)
      end if
    end subroutine add_atom_to_output_list

    subroutine add_shell_to_output_list()
      use types, only: realloc

      integer :: i

      do i = 1, nshel
         if (abs(rshel(i) - dsqrt) < shell_eps) return
      end do
      nshel = nshel + 1
      if (nshel > size(rshel,1)) call realloc(rshel,2*nshel)
      rshel(nshel) = dsqrt

    end subroutine add_shell_to_output_list

  end subroutine list_near_atoms

  !> Given the point xp (in icrd coordinates), calculate the
  !> nearest atom up to a distance distmax (or up to infinity, if no
  !> distmax is given). The nearest atom has ID nid from the complete
  !> list (atcel) and is at a distance dist, or nid=0 and dist=distmax
  !> if the search did not produce any atoms up to distmax. On output,
  !> the optional argument lvec contains the lattice vector to the
  !> nearest atom (i.e. its position is atcel(nid)%x + lvec).
  !> Optional input:
  !> - nid0 = consider only atoms with index nid0 from the nneq list.
  !> - id0 = consider only atoms with index id0 from the complete list.
  !> - iz0 = consider only atoms with atomic number iz0.
  !> - ispc0 = consider only species with species ID ispc0.
  !> - nozero = disregard zero-distance atoms.
  !>
  !> This routine is thread-safe.
  module subroutine nearest_atom(c,xp,icrd,nid,dist,distmax,lvec,nid0,id0,iz0,ispc0,nozero)
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    integer, intent(out) :: nid
    real*8, intent(out) :: dist
    real*8, intent(in), optional :: distmax
    integer, intent(out), optional :: lvec(3)
    integer, intent(in), optional :: nid0
    integer, intent(in), optional :: id0
    integer, intent(in), optional :: iz0
    integer, intent(in), optional :: ispc0
    logical, intent(in), optional :: nozero

    integer :: nat
    integer, allocatable :: eid(:), lvec_(:,:)
    real*8, allocatable :: dist_(:)
    logical :: dodist

    ! whether to use the distance (cannot be too long)
    dodist = present(distmax)
    if (dodist) dodist = (distmax < 1d0)

    ! get just one atom or all sorted atoms up to distmax
    if (dodist) then
       call c%list_near_atoms(xp,icrd,.true.,nat,eid=eid,dist=dist_,lvec=lvec_,up2d=distmax,&
          nid0=nid0,id0=id0,iz0=iz0,ispc0=ispc0,nozero=nozero)
    else
       call c%list_near_atoms(xp,icrd,.false.,nat,eid=eid,dist=dist_,lvec=lvec_,up2n=1,&
          nid0=nid0,id0=id0,iz0=iz0,ispc0=ispc0,nozero=nozero)
    end if

    ! if no atoms in output, return
    nid = 0
    dist = huge(1d0)
    if (nat == 0) return

    ! if distance is higher than distmax, return
    if (present(distmax).and..not.dodist) then
       if (dist_(1) > distmax) return
    end if

    ! write and finish
    nid = eid(1)
    dist = dist_(1)
    if (present(lvec)) lvec = lvec_(:,1)

  end subroutine nearest_atom

  !> Returns half of the nearest neighbor distance for non-equivalent
  !> atom ineq.
  module function get_rnn2(c,ineq)
    use param, only: icrd_cart
    class(crystal), intent(inout) :: c
    integer, intent(in) :: ineq
    real*8 :: get_rnn2

    integer :: iaux

    ! special case: only 1 atom in the molecule
    if (c%ismolecule .and. c%ncel <= 1) then
       get_rnn2 = 0d0
       return
    end if

    ! normal case
    call c%nearest_atom(c%at(ineq)%r,icrd_cart,iaux,get_rnn2,nozero=.true.)
    get_rnn2 = 0.5d0 * get_rnn2

  end function get_rnn2

  !> Given point x0 (with icrd input coordinates), if x0 corresponds
  !> to an atomic position (to within distmax or atomeps if distmax is
  !> not given), return the complete-list ID of the atom. Otherwise,
  !> return 0. Optionally, return the lattice vector translation
  !> (lvec) and the distance (dist) to the closest atom. This routine
  !> is thread-safe.
  module function identify_atom(c,x0,icrd,lvec,dist,distmax)
    use global, only: atomeps
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: icrd
    integer, intent(out), optional :: lvec(3)
    real*8, intent(out), optional :: dist
    real*8, intent(in), optional :: distmax
    integer :: identify_atom

    real*8 :: distmax_, dist0
    integer :: lvec0(3)

    identify_atom = 0
    distmax_ = atomeps
    if (present(distmax)) distmax_ = distmax

    call c%nearest_atom(x0,icrd,identify_atom,dist0,distmax=distmax_,lvec=lvec0)
    if (present(lvec)) lvec = lvec0
    if (present(dist)) dist = dist0

  end function identify_atom

  !> Calculate the core (if zpsp is present) or promolecular densities
  !> at a point x0 (coord format given by icrd) using atomic radial
  !> grids up to a number of derivatives nder (max: 2). Returns the
  !> density (f), gradient (fp, nder >= 1), and Hessian (fpp, nder >=
  !> 2). If a fragment (fr) is given, then only the atoms in it
  !> contribute. This routine is thread-safe.
  module subroutine promolecular_atom(c,x0,icrd,f,fp,fpp,nder,zpsp,fr)
    use grid1mod, only: cgrid, agrid, grid1
    use global, only: cutrad
    use fragmentmod, only: fragment
    use tools_io, only: ferror, faterr
    use param, only: icrd_crys, icrd_rcrys, maxzat
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: icrd
    real*8, intent(out) :: f
    real*8, intent(out) :: fp(3)
    real*8, intent(out) :: fpp(3,3)
    integer, intent(in) :: nder
    integer, intent(in), optional :: zpsp(:)
    type(fragment), intent(in), optional :: fr

    integer :: i, j, k, ii, iz
    real*8 :: xc(3), xx(3), rlvec(3), r, rinv1, rinv1rp
    real*8 :: rho, rhop, rhopp, rfac, radd, rmax
    logical :: iscore
    type(grid1), pointer :: g
    integer :: nat
    integer, allocatable :: nid(:), lvec(:,:)
    real*8, allocatable :: dist(:), rcutmax(:,:)
    logical, allocatable :: isinfr(:)

    ! initialize
    f = 0d0
    fp = 0d0
    fpp = 0d0
    iscore = present(zpsp)
    if (iscore) then
       if (all(zpsp <= 0)) return
    end if
    if (iscore.and..not.allocated(cgrid)) then
       call ferror("promolecular","cgrid not allocated",faterr)
    elseif (.not.iscore.and..not.allocated(agrid)) then
       call ferror("promolecular","agrid not allocated",faterr)
    end if

    ! convert to Cartesian
    if (icrd == icrd_crys) then
       xc = c%x2c(x0)
    elseif (icrd == icrd_rcrys) then
       xc = c%xr2c(x0)
    else
       xc = x0
    end if

    ! calculate the cutoffs
    allocate(rcutmax(c%nspc,2))
    rcutmax = 0d0
    rmax = 0d0
    do i = 1, c%nspc
       iz = c%spc(i)%z
       if (iz == 0 .or. iz > maxzat) cycle
       if (iscore) then
          if (zpsp(i) <= 0) cycle
          g => cgrid(iz,zpsp(i))
       else
          g => agrid(iz)
       end if
       if (g%isinit) then
          rcutmax(i,2) = min(cutrad(iz),g%rmax)
       else
          rcutmax(i,2) = cutrad(iz)
       end if
       rmax = max(rcutmax(i,2),rmax)
    end do

    ! compute the list of atoms that contribute to the point
    call c%list_near_atoms(x0,icrd,.false.,nat,nid,dist,lvec,up2dsp=rcutmax)
    deallocate(rcutmax)
    if (nat == 0) return

    ! if fragment is provided, use only the atoms in it (by their cell index)
    if (present(fr)) then
       allocate(isinfr(c%ncel))
       isinfr = .false.
       do i = 1, fr%nat
          isinfr(fr%at(i)%cidx) = .true.
       end do
    end if

    ! Do the density and derivatives sum
    do ii = 1, nat
       i = nid(ii)
       if (present(fr)) then
          if (.not.isinfr(i)) cycle
       end if

       iz = c%spc(c%atcel(i)%is)%z
       if (iz == 0 .or. iz > maxzat) cycle
       r = dist(ii)
       if (r > cutrad(iz)) cycle

       if (iscore) then
          if (zpsp(c%atcel(i)%is) <= 0) cycle
          g => cgrid(iz,zpsp(c%atcel(i)%is))
       else
          g => agrid(iz)
       end if
       if (.not.g%isinit) cycle

       r = max(max(r,g%r(1)),1d-14)
       call g%interp(r,rho,rhop,rhopp)
       rho = max(rho,0d0)
       f = f + rho

       if (nder < 1) cycle
       rlvec = real(lvec(:,ii),8)
       rlvec = c%x2c(rlvec)
       xx = xc - (c%atcel(i)%r + rlvec)
       rinv1 = 1d0 / r
       rinv1rp = rinv1 * rhop
       fp = fp + xx * rinv1rp

       if (nder < 2) cycle
       rfac = (rhopp - rinv1rp) * rinv1 * rinv1
       do j = 1, 3
          fpp(j,j) = fpp(j,j) + rinv1rp + rfac * xx(j) * xx(j)
          do k = 1, j-1
             radd = rfac * xx(j) * xx(k)
             fpp(j,k) = fpp(j,k) + radd
             fpp(k,j) = fpp(k,j) + radd
          end do
       end do
    end do

    if (allocated(isinfr)) deallocate(isinfr)

  end subroutine promolecular_atom

  !> Find the covalent connectivity and return the bonds in the
  !> c%nstar array. Two atoms i and j are connected if:
  !> a) Their distance is less than bondfac * (atmrad(i) + atmrad(j))
  !> b) The distance is at most dbond (0.2 ang) higher than the minimum distance.
  !> If i and j are metals, use b. If i and j are non-metals, use a.
  !> If one is metal and the other non-metal, apply both.
  module subroutine find_asterisms(c,nstar,atmrad,bondfac,rij,bonddelta)
    use global, only: bonddelta_def
    use param, only: maxzat0, ismetal
    use tools_io, only: string, ferror, faterr
    use tools, only: qcksort
    use types, only: realloc
    use param, only: icrd_crys, atmcov, maxzat
    class(crystal), intent(inout) :: c
    type(neighstar), allocatable, intent(inout) :: nstar(:)
    real*8, intent(in), optional :: atmrad(0:maxzat0)
    real*8, intent(in), optional :: bondfac
    real*8, intent(in), optional :: rij(:,:,:)
    real*8, intent(in), optional :: bonddelta

    integer :: i, j, iz, jz, jid
    real*8 :: ri, rj, dmax, dd
    logical :: bonded, ism, jsm
    real*8, allocatable :: rij2(:,:,:)
    real*8 :: dbond

    type atenv_type
       integer :: nat
       integer, allocatable :: eid(:)
       integer, allocatable :: lvec(:,:)
       real*8, allocatable :: dist(:)
    end type atenv_type
    type(atenv_type), allocatable :: atenv(:)

    ! return if there are no atoms
    if (c%ncel == 0) return

    ! initialize
    if (present(bonddelta)) then
       dbond = bonddelta
    else
       dbond = bonddelta_def
    end if

    ! allocate the asterism arrays
    if (allocated(nstar)) deallocate(nstar)
    allocate(nstar(c%ncel))
    do i = 1, c%ncel
       nstar(i)%ncon = 0
       if (allocated(nstar(i)%idcon)) deallocate(nstar(i)%idcon)
       if (allocated(nstar(i)%lcon)) deallocate(nstar(i)%lcon)
       if (allocated(nstar(i)%ordcon)) deallocate(nstar(i)%ordcon)
       if (allocated(nstar(i)%aromdir)) deallocate(nstar(i)%aromdir)
       allocate(nstar(i)%idcon(20))
       allocate(nstar(i)%lcon(3,20))
       allocate(nstar(i)%ordcon(20))
       allocate(nstar(i)%aromdir(3,20))
    end do

    ! pre-calculate the distance matrix
    allocate(rij2(c%nspc,2,c%nspc))
    rij2 = 0d0
    if (present(atmrad).and.present(bondfac)) then
       ! from atmrad and bondfac
       do i = 1, c%nspc
          if (c%spc(i)%z <= 0) cycle
          if (c%spc(i)%z > maxzat) cycle ! CPs do not form bonds
          ri = atmrad(c%spc(i)%z)
          do j = 1, c%nspc
             if (c%spc(j)%z <= 0) cycle
             if (c%spc(j)%z > maxzat) cycle ! CPs do not form bonds
             rj = atmrad(c%spc(j)%z)

             rij2(j,2,i) = (ri+rj) * bondfac
             rij2(j,1,i) = 0d0
             rij2(j,2,i) = rij2(j,2,i) * rij2(j,2,i)
             rij2(j,1,i) = rij2(j,1,i) * rij2(j,1,i)
          end do
       end do
    elseif (present(rij)) then
       ! from rij
       if (any(shape(rij) /= shape(rij2))) &
          call ferror('find_asterisms','error in rij shape',faterr)
       rij2 = rij * rij
    else
       call ferror('find_asterisms','error in input arguments',faterr)
    end if
    rij2 = sqrt(rij2)

    ! calculate the covalent radius
    dmax = 0d0
    do i = 1, c%nspc
       if (c%spc(i)%z <= 0) cycle
       if (c%spc(i)%z > maxzat) cycle
       dmax = max(atmcov(c%spc(i)%z),dmax)
    end do

    ! calculate the atomic environments
    allocate(atenv(c%ncel))
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       if (iz > maxzat) cycle
       call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.true.,atenv(i)%nat,atenv(i)%eid,&
          atenv(i)%dist,atenv(i)%lvec,up2d=(atmcov(iz) + dmax)*2d0,nozero=.true.)
    end do

    ! process the environments and generate the bonds
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       if (iz > maxzat) cycle
       ism = ismetal(iz)

       do j = 1, atenv(i)%nat
          jid = atenv(i)%eid(j)
          jz = c%spc(c%atcel(jid)%is)%z
          if (jz > maxzat) cycle
          jsm = ismetal(jz)
          bonded = .false.

          dd = atenv(i)%dist(j)
          if (ism .and. jsm) then
             if (atenv(i)%nat > 0 .and. atenv(jid)%nat > 0) then
                bonded = (dd < atenv(i)%dist(1) + dbond) .or. (dd < atenv(jid)%dist(1) + dbond)
             end if
          elseif (ism) then
             if (atenv(i)%nat > 0) then
                bonded = (dd < atenv(i)%dist(1) + dbond) .or.&
                   (dd > rij2(c%atcel(jid)%is,1,c%atcel(i)%is) .and. dd < rij2(c%atcel(jid)%is,2,c%atcel(i)%is))
             end if
          elseif (jsm) then
             if (atenv(jid)%nat > 0) then
                bonded = (dd < atenv(jid)%dist(1) + dbond) .or.&
                   (dd > rij2(c%atcel(jid)%is,1,c%atcel(i)%is) .and. dd < rij2(c%atcel(jid)%is,2,c%atcel(i)%is))
             end if
          else
             bonded = (dd > rij2(c%atcel(jid)%is,1,c%atcel(i)%is) .and. dd < rij2(c%atcel(jid)%is,2,c%atcel(i)%is))
          end if

          ! add the bonded atom
          if (bonded) then
             nstar(i)%ncon = nstar(i)%ncon + 1
             if (nstar(i)%ncon > size(nstar(i)%idcon,1)) then
                call realloc(nstar(i)%idcon,2*nstar(i)%ncon)
                call realloc(nstar(i)%lcon,3,2*nstar(i)%ncon)
                call realloc(nstar(i)%ordcon,2*nstar(i)%ncon)
                call realloc(nstar(i)%aromdir,3,2*nstar(i)%ncon)
             end if
             nstar(i)%idcon(nstar(i)%ncon) = jid
             nstar(i)%lcon(:,nstar(i)%ncon) = atenv(i)%lvec(:,j)
             nstar(i)%ordcon(nstar(i)%ncon) = 1
             nstar(i)%aromdir(:,nstar(i)%ncon) = 0d0
          end if
       end do

       ! reallocate
       call realloc(nstar(i)%idcon,nstar(i)%ncon)
       call realloc(nstar(i)%lcon,3,nstar(i)%ncon)
       call realloc(nstar(i)%ordcon,nstar(i)%ncon)
       call realloc(nstar(i)%aromdir,3,nstar(i)%ncon)
    end do

    ! perceive bond orders for the organic subgraph
    call perceive_bond_orders(c,nstar)

  end subroutine find_asterisms

  !> Given the point xp (in icrd coordinates), calculate the list of
  !> nearest lattice points. If sorted is true, the list is sorted by
  !> distance. One (and only one) of these criteria must be
  !> given:
  !> - up2d = list up to a distance up2d.
  !> - up2n = up to a number of atoms.
  !> If up2n is present, the list is sorted by distance on output
  !> regardless of the value of sorted.
  !>
  !> Optional input:
  !> - ndiv = divide the parent lattice vectors by ndiv; useful for grids
  !> - nozero = disregard zero-distance lattice points.
  !>
  !> Output:
  !> - nat = the list contains nat lattice points.
  !>
  !> Optional output:
  !> - dist(nid) = distance from xp.
  !> - lvec(3,nid) = lattice vector coordinates (crystallographic).
  !>
  !> This routine is thread-safe.
  module subroutine list_near_lattice_points(c,xp,icrd,sorted,nat,dist,lvec,x2c,ndiv,&
     up2d,up2n,nozero)
    use types, only: realloc
    use tools, only: mergesort, wscell
    use tools_math, only: det3, cross
    use tools_io, only: ferror, faterr
    use param, only: icrd_cart, icrd_crys
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    logical, intent(in) :: sorted
    integer, intent(out) :: nat
    real*8, allocatable, intent(inout), optional :: dist(:)
    integer, allocatable, intent(inout), optional :: lvec(:,:)
    real*8, intent(in), optional :: x2c(3,3)
    integer, intent(in), optional :: ndiv(3)
    real*8, intent(in), optional :: up2d
    integer, intent(in), optional :: up2n
    logical, intent(in), optional :: nozero

    logical :: nozero_
    logical :: sorted_, ok
    integer :: nb, nshellb, nsafe, idx
    integer :: nx(3), nn, i0(3), i1(3), i, j, k, lvecx(3)
    real*8 :: up2d_2, x(3), lvecini(3), blockcv(3), blockomega, xorigc(3)
    real*8 :: xdif(3), dd, dmax, x0(3)
    integer, allocatable :: idb(:,:), iaux(:)
    integer, allocatable :: at_lvec(:,:), iord(:)
    real*8, allocatable :: at_dist(:)
    real*8 :: x2c_(3,3), x2xr(3,3), xr2x(3,3), xr2c(3,3), c2xr(3,3)

    real*8, parameter :: eps = 1d-20

    ! process input
    nozero_ = .false.
    sorted_ = sorted .or. present(up2n)
    if (present(nozero)) nozero_ = nozero
    if (present(up2d)) then
       up2d_2 = up2d * up2d
    elseif (present(up2n)) then
       !
    else
       call ferror("list_near_lattice_points","must give one of up2d or up2n",faterr)
    end if

    if (present(x2c)) then
       if (present(ndiv)) then
          do i = 1, 3
             x2c_(:,i) = x2c(:,i) / real(ndiv(i),8)
          end do
       else
          x2c_ = x2c
       end if
       x0 = xp * ndiv
       call wscell(x2c_,.true.,m_xr2x=xr2x,m_xr2c=xr2c,m_x2xr=x2xr,m_c2xr=c2xr)
    elseif (present(ndiv)) then
       do i = 1, 3
          x2c_(:,i) = c%m_x2c(:,i) / real(ndiv(i),8)
       end do
       x0 = xp * ndiv
       call wscell(x2c_,.true.,m_xr2x=xr2x,m_xr2c=xr2c,m_x2xr=x2xr,m_c2xr=c2xr)
    else
       x0 = xp
       xr2c = c%m_xr2c
       c2xr = c%m_c2xr
       xr2x = c%m_xr2x
       x2xr = c%m_x2xr
    end if

    ! translate the point to the main reduced cell
    if (icrd == icrd_crys) then
       x = matmul(x2xr,x0)
    elseif (icrd == icrd_cart) then
       x = matmul(c2xr,x0)
    else
       x = x0
    end if
    lvecini = floor(x)
    x = x - lvecini
    xorigc = matmul(xr2c,x)

    ! allocate space for lattice points
    nat = 0
    nn = 20
    if (present(up2n)) nn = min(up2n,20)
    allocate(at_dist(nn),at_lvec(3,nn))
    at_dist = 0d0
    at_lvec = 0

    ! run the search
    if (present(up2n)) then
       ! search a number of shells or a number of points (up2n) !!!

       ! Run over shells of nearby blocks and find atoms until we
       ! have at least the given number of shells
       nshellb = 0
       do while(.true.)
          call make_block_shell(c,nshellb,nb,idb,dmax)

          do i = 1, nb
             xdif = matmul(xr2c,real(idb(:,i),8)) - xorigc
             dd = dot_product(xdif,xdif)

             ! check if we should add the atom to the list
             ok = .true.
             if (nozero_ .and. dd < eps) ok = .false.
             if (ok) then
                nat = nat + 1
                if (nat > size(at_dist,1)) then
                   call realloc(at_dist,2*nat)
                   call realloc(at_lvec,3,2*nat)
                end if
                at_dist(nat) = sqrt(dd)
                at_lvec(:,nat) = nint(matmul(xr2x,lvecini + idb(:,i)))
             end if
          end do

          ! check if we have enough atoms
          nsafe = count(at_dist(1:nat) <= dmax)
          if (nsafe >= up2n) exit

          ! next shell
          nshellb = nshellb + 1
       end do ! while loop

       if (up2n == 1) then
          ! bypass for efficiency: handle the case when up2n = 1 (nearest point)
          idx = minloc(at_dist(1:nat),1)
          nat = 1
          at_dist(1) = at_dist(idx)
          at_lvec(:,1) = at_lvec(:,idx)
          sorted_ = .false.
       else
          ! filter out the unneeded atoms
          nsafe = count(at_dist(1:nat) <= dmax)
          allocate(iaux(nsafe))
          nsafe = 0
          do i = 1, nat
             if (at_dist(i) <= dmax) then
                nsafe = nsafe + 1
                iaux(nsafe) = i
             end if
          end do
          nat = nsafe
          at_dist(1:nsafe) = at_dist(iaux)
          at_lvec(:,1:nsafe) = at_lvec(:,iaux)
          call realloc(at_dist,nsafe)
          call realloc(at_lvec,3,nsafe)
          deallocate(iaux)
       end if
    else
       ! search atoms up to a given distance (up2d) !!!

       ! calculate the number of blocks in each direction required for satifying
       ! that the largest sphere in the super-block has radius > dmax.
       ! r = Vblock / 2 / max(cv(3)/n3,cv(2)/n2,cv(1)/n1)
       blockcv(1) = norm2(cross(xr2c(:,2),xr2c(:,3)))
       blockcv(2) = norm2(cross(xr2c(:,1),xr2c(:,3)))
       blockcv(3) = norm2(cross(xr2c(:,1),xr2c(:,2)))
       blockomega = det3(xr2c)
       nx = ceiling(blockcv / (0.5d0 * blockomega / max(up2d,1d-40)))

       ! define the search space
       do i = 1, 3
          if (mod(nx(i),2) == 0) then
             i0(i) = - nx(i)/2
             i1(i) = + nx(i)/2
          else
             i0(i) = - (nx(i)/2+1)
             i1(i) = + (nx(i)/2+1)
          end if
       end do

       ! run over all lattice points
       do i = i0(1), i1(1)
          do j = i0(2), i1(2)
             do k = i0(3), i1(3)
                lvecx = (/i,j,k/)
                xdif = matmul(xr2c,real(lvecx,8)) - xorigc
                dd = dot_product(xdif,xdif)

                ! check if we should add the lattice point to the list
                ok = .true.
                if (nozero_ .and. dd < eps) then
                   ok = .false.
                elseif (present(up2d)) then
                   ok = (dd <= up2d_2)
                end if
                if (ok) then
                   nat = nat + 1
                   if (nat > size(at_dist,1)) then
                      call realloc(at_dist,2*nat)
                      call realloc(at_lvec,3,2*nat)
                   end if
                   at_dist(nat) = sqrt(dd)
                   at_lvec(:,nat) = nint(matmul(xr2x,lvecini + lvecx))
                end if
             end do ! k = i0(3), i1(3)
          end do ! j = i0(2), i1(2)
       end do ! i = i0(1), i1(1)
    end if

    ! sort if necessary
    if (sorted_ .and. nat > 1) then
       ! permutation vector
       allocate(iord(nat))
       do i = 1, nat
          iord(i) = i
       end do

       ! sort by distance
       call mergesort(at_dist,iord,1,nat)

       ! reorder and clean up
       at_dist = at_dist(iord)
       at_lvec = at_lvec(:,iord)
       deallocate(iord)
    end if

    ! reduce the list if up2n (always sorted)
    if (present(up2n)) nat = up2n

    ! output
    if (present(dist)) dist = at_dist(1:nat)
    if (present(lvec)) lvec = at_lvec(:,1:nat)

  end subroutine list_near_lattice_points

  !> Given the point xp (in icrd coordinates), calculate the nearest
  !> lattice point. The lattice point (in cryst. coords.) is returne
  !> in lvec and the distance in dist. If nozero, disregard
  !> zero-distance lattice points.  If ndiv is given, divide the
  !> parent lattice vectors by ndiv; useful for grids. This routine
  !> is thread-safe.
  module subroutine nearest_lattice_point(c,xp,icrd,dist,lvec,ndiv,nozero)
    class(crystal), intent(inout) :: c
    real*8, intent(in) :: xp(3)
    integer, intent(in) :: icrd
    real*8, intent(out) :: dist
    integer, intent(out), optional :: lvec(3)
    integer, intent(in), optional :: ndiv(3)
    logical, intent(in), optional :: nozero

    integer :: nat
    integer, allocatable :: lvec_(:,:)
    real*8, allocatable :: dist_(:)

    ! get just one atom
    call c%list_near_lattice_points(xp,icrd,.false.,nat,dist=dist_,&
       lvec=lvec_,ndiv=ndiv,up2n=1,nozero=nozero)

    ! if no atoms in output, return
    if (present(lvec)) lvec = 0
    dist = huge(1d0)
    if (nat == 0) return

    ! write and finish
    dist = dist_(1)
    if (present(lvec)) lvec = lvec_(:,1)

  end subroutine nearest_lattice_point

  ! Given atoms with complete cell IDs i and j and corresponding
  ! lattice vectors il and jl, return .true. if they are in the same
  ! molecule.
  module function in_same_molecule(c,i,il,j,jl)
    class(crystal), intent(inout) :: c
    integer, intent(in) :: i, j
    integer, intent(in) :: il(3), jl(3)
    logical :: in_same_molecule

    integer :: imol, jmol

    imol = c%idatcelmol(1,i)
    jmol = c%idatcelmol(1,j)
    if (imol /= jmol) then
       in_same_molecule = .false.
    elseif (.not.c%mol(imol)%discrete) then
       in_same_molecule = .true.
    else
       in_same_molecule = all(il - c%mol(imol)%at(c%idatcelmol(2,i))%lvec == &
          jl - c%mol(jmol)%at(c%idatcelmol(2,j))%lvec)
    end if

  end function in_same_molecule

  !xx! private procedures

  ! Find the indices for the nth shell of blocks. Sets the number of indices (nb),
  ! the indices themselves (idb(3,nb)), and the rmax for this shell (dmax).
  subroutine make_block_shell(c,n,nb,idb,dmax)
    type(crystal), intent(in) :: c
    integer, intent(in) :: n
    integer, intent(out) :: nb
    integer, allocatable, intent(inout) :: idb(:,:)
    real*8, intent(out) :: dmax

    integer :: i, j, ic

    ! special case: n = 0
    if (n == 0) then
       nb = 1
       if (allocated(idb)) deallocate(idb)
       allocate(idb(3,nb))
       idb = 0
       dmax = 0d0
       return
    end if

    ! allocate space for (2n+1)^3 - (2n-1)^3 blocks
    nb = 2 + 24 * n * n
    if (allocated(idb)) deallocate(idb)
    allocate(idb(3,nb))

    ! list them
    ic = 0
    do i = -n+1, n-1
       do j = -n+1, n-1
          ic = ic + 1
          idb(:,ic) = (/n,i,j/)
          ic = ic + 1
          idb(:,ic) = (/-n,i,j/)

          ic = ic + 1
          idb(:,ic) = (/i,n,j/)
          ic = ic + 1
          idb(:,ic) = (/i,-n,j/)

          ic = ic + 1
          idb(:,ic) = (/i,j,n/)
          ic = ic + 1
          idb(:,ic) = (/i,j,-n/)
       end do
       ic = ic + 1
       idb(:,ic) = (/n,n,i/)
       ic = ic + 1
       idb(:,ic) = (/n,-n,i/)
       ic = ic + 1
       idb(:,ic) = (/-n,n,i/)
       ic = ic + 1
       idb(:,ic) = (/-n,-n,i/)

       ic = ic + 1
       idb(:,ic) = (/n,i,n/)
       ic = ic + 1
       idb(:,ic) = (/n,i,-n/)
       ic = ic + 1
       idb(:,ic) = (/-n,i,n/)
       ic = ic + 1
       idb(:,ic) = (/-n,i,-n/)

       ic = ic + 1
       idb(:,ic) = (/i,n,n/)
       ic = ic + 1
       idb(:,ic) = (/i,n,-n/)
       ic = ic + 1
       idb(:,ic) = (/i,-n,n/)
       ic = ic + 1
       idb(:,ic) = (/i,-n,-n/)
    end do

    ic = ic + 1
    idb(:,ic) = (/n,n,n/)
    ic = ic + 1
    idb(:,ic) = (/n,n,-n/)
    ic = ic + 1
    idb(:,ic) = (/n,-n,n/)
    ic = ic + 1
    idb(:,ic) = (/n,-n,-n/)
    ic = ic + 1
    idb(:,ic) = (/-n,n,n/)
    ic = ic + 1
    idb(:,ic) = (/-n,n,-n/)
    ic = ic + 1
    idb(:,ic) = (/-n,-n,n/)
    ic = ic + 1
    idb(:,ic) = (/-n,-n,-n/)

    ! calculate the radius of the largest sphere contained in the blocked
    ! region (dmax)
    dmax = 0.5d0 * c%blockomega / maxval(c%blockcv / real(n,8))

  end subroutine make_block_shell

  !> Given a position in reduced crystallographic coordinates x,
  !> calculate the region of blocks (between i0 and i1) required to
  !> inscribe a sphere of radius rmax centered on x. On output, rmax
  !> is increased to the radius of the largest sphere centered on x
  !> inscribed in i0 -> i1.
  !>
  !> Calculation of distance between point inside the block and its six
  !> sides, assuming the corners of the block have (0,0,0) and (1,1,1)
  !> coordinates:
  !>   d1l = |(x-x0) * (a x b)| / ||a x b|| = |x3|   * V / ||a x b||
  !>   d1u = |(x-x0) * (a x b)| / ||a x b|| = |x3-1| * V / ||a x b||
  !>   d2l = |(x-x0) * (b x c)| / ||b x c|| = |x1|   * V / ||b x c||
  !>   d2u = |(x-x0) * (b x c)| / ||b x c|| = |x1-1| * V / ||b x c||
  !>   d3l = |(x-x0) * (c x a)| / ||c x a|| = |x2|   * V / ||c x a||
  !>   d3u = |(x-x0) * (c x a)| / ||c x a|| = |x2-1| * V / ||c x a||
  !> which results in:
  !>  dist = (xn(i) - i0(i)) * V / cv(i)
  !>  dist = (i1(i) + 1 - xn(i)) * V / cv(i)
  subroutine blocks_inscribed_sphere(c,x,rmax,i0,i1)
    type(crystal), intent(in) :: c
    real*8, intent(in) :: x(3)
    real*8, intent(inout) :: rmax
    integer :: i0(3), i1(3)

    real*8 :: xn(3), fac(3)

    ! block coordinates of the  point
    xn = x * c%nblock

    ! calculate the region
    fac = rmax * c%blockcv / c%blockomega
    i0 = floor(xn - fac)
    i1 = ceiling(xn - 1 + fac)

    ! calculate the new rmax
    rmax = minval((xn - i0) * c%blockomega / c%blockcv)
    rmax = min(rmax,minval((i1 + 1 - xn) * c%blockomega / c%blockcv))

  end subroutine blocks_inscribed_sphere

  !> Perceive integer bond orders from the molecular geometry and
  !> store them in nstar(:)%ordcon. Uses a rule-based algorithm
  !> heavily adapted from:
  !>   Zhang et al., J. Cheminform. 4 (2012) 26, doi:10.1186/1758-2946-4-26
  !> The initial connectivity is derived from critic2's own rules.
  subroutine perceive_bond_orders(c,nstar)
    use types, only: neighstar
    use tools, only: qcksort
    use tools_math, only: cross
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    type(neighstar), intent(inout) :: nstar(:)

    integer :: i, k, ia, ib, nb, mb, maxdeg, za, zb, ncand
    integer, allocatable :: bia(:), bib(:), blv(:,:), bord(:)
    logical, allocatable :: bfix(:)
    real*8, allocatable :: bdist(:)
    integer, allocatable :: deg(:), cap(:), iopen(:), nincid(:), incid(:,:)
    integer, allocatable :: cand(:), iord(:)
    real*8, allocatable :: ckey(:)
    real*8 :: r1, r2, r3, pyr
    integer :: o
    logical :: islin
    logical, allocatable :: arom(:)  ! per-atom aromaticity flag
    logical, allocatable :: barom(:) ! per-bond aromatic-ring flag
    real*8, allocatable :: bdir(:,:) ! per-bond outward (ring-exterior) direction
    ! maximum-weight matching workspace (shared with the contained solver)
    integer, allocatable :: ca(:), cb(:), cbond(:), wopen(:)
    logical, allocatable :: cursel(:), bestsel(:)
    real*8, allocatable :: cw(:), remw(:)
    real*8 :: bestscore
    integer :: nodecount

    ! recursive double-bond search
    integer, parameter :: maxnode = 2000000 ! maximum node count
    integer, parameter :: ncand_max = 300 ! maximum number of candidates (otherwise, greedy)

    if (c%ncel == 0) return

    ! per-atom capacity (max valence) of the supported elements; <0 = out of
    ! scope (metals, unsupported elements) -> their bonds are left single
    ! supported elements: H, C, N, O, F, P, S, Cl, Br, I
    allocate(cap(c%ncel),deg(c%ncel))
    cap = -1
    deg = 0
    do ia = 1, c%ncel
       cap(ia) = zhang_maxval(c%spc(c%atcel(ia)%is)%z)
    end do

    ! count neighbors (degree) and total bonds; return if no bonds
    nb = 0
    do ia = 1, c%ncel
       if (cap(ia) < 0) cycle
       do k = 1, nstar(ia)%ncon
          ib = nstar(ia)%idcon(k)
          if (cap(ib) < 0) cycle
          if (ib == ia) cycle ! skip self-image bonds (left single, see below)
          deg(ia) = deg(ia) + 1
          if (ia < ib) nb = nb + 1
       end do
    end do
    if (nb == 0) return
    maxdeg = maxval(deg)

    ! build the unique bond list and the per-atom incidence list
    allocate(bia(nb),bib(nb),blv(3,nb),bdist(nb),bord(nb),bfix(nb))
    allocate(nincid(c%ncel),incid(maxdeg,c%ncel))
    bord = 1
    bfix = .false.
    nincid = 0
    mb = 0
    do ia = 1, c%ncel
       if (cap(ia) < 0) cycle
       do k = 1, nstar(ia)%ncon
          ib = nstar(ia)%idcon(k)
          if (cap(ib) < 0) cycle
          if (ib <= ia) cycle ! count each bond once; skip self-image bonds
          mb = mb + 1
          bia(mb) = ia
          bib(mb) = ib
          blv(:,mb) = nstar(ia)%lcon(:,k)
          bdist(mb) = norm2(c%x2c(c%atcel(ib)%x + real(blv(:,mb),8) - c%atcel(ia)%x))
          nincid(ia) = nincid(ia) + 1
          incid(nincid(ia),ia) = mb
          nincid(ib) = nincid(ib) + 1
          incid(nincid(ib),ib) = mb
       end do
    end do

    ! Open valence available for upgrading bonds beyond single (cap - degree)
    ! No over-connected correction: bonds are controlled by the user
    allocate(iopen(c%ncel))
    iopen = 0
    do ia = 1, c%ncel
       if (cap(ia) < 0) cycle
       iopen(ia) = max(cap(ia) - deg(ia),0)
    end do

    ! Hard rule 3.3 (reinterpretation): pyramidal carbons cannot form
    ! multiple bonds. The planarity criterion is different from the
    ! article.
    do ia = 1, c%ncel
       if (cap(ia) < 0) cycle
       if (deg(ia) == 3 .and. c%spc(c%atcel(ia)%is)%z == 6) then
          pyr = pyramidality(ia)
          if (pyr > 0.4d0) iopen(ia) = 0
       end if
    end do

    ! Hard rule 2.3, triple bonds (reintrepretation): detected by a
    ! combination of length and linearity, rather than the paper's
    ! heuristic.
    do i = 1, nb
       if (bfix(i)) cycle
       ia = bia(i)
       ib = bib(i)
       if (iopen(ia) < 2 .or. iopen(ib) < 2) cycle
       if (deg(ia) > 2 .or. deg(ib) > 2) cycle
       za = c%spc(c%atcel(ia)%is)%z
       zb = c%spc(c%atcel(ib)%is)%z
       call zhang_refbond(za,zb,r1,r2,r3)
       if (r3 < 0d0) cycle
       if (bdist(i) >= r3 + 0.1d0/bohrtoa) cycle
       islin = .true.
       if (deg(ia) == 2) islin = islin .and. (atom_angle(ia) > 155d0)
       if (deg(ib) == 2) islin = islin .and. (atom_angle(ib) > 155d0)
       if (.not.islin) cycle
       bord(i) = 3
       bfix(i) = .true.
       iopen(ia) = iopen(ia) - 2
       iopen(ib) = iopen(ib) - 2
    end do

    ! Hard rule: oxyacid groups (N/P/S/Cl/Br/I with terminal O; Zhang rules
    ! 3.1/3.2/4.2). Assign localized X=O double bonds (nitro/nitrate, phosphate,
    ! sulfate, perchlorate/chlorate/iodate, ...)
    do ia = 1, c%ncel
       za = -1
       if (cap(ia) >= 0) za = c%spc(c%atcel(ia)%is)%z
       if (za /= 7 .and. za /= 15 .and. za /= 16 .and. za /= 17 .and. za /= 35 .and. za /= 53) cycle
       call assign_oxyacid(ia,za)
    end do

    ! length rule (reinterpreted): a candidate bond is unfixed, both
    ! ends are still open, and it is shorter than the single bond
    ! reference. Write down the bonds that are candidates for
    ! promotion to double in cand(.). ckey(.) is the bond distance
    ! minus the reference single-bond distance (>0).
    allocate(cand(nb),ckey(nb))
    ncand = 0
    do i = 1, nb
       if (bfix(i)) cycle
       ia = bia(i)
       ib = bib(i)
       if (iopen(ia) < 1 .or. iopen(ib) < 1) cycle
       za = c%spc(c%atcel(ia)%is)%z
       zb = c%spc(c%atcel(ib)%is)%z
       call zhang_refbond(za,zb,r1,r2,r3)
       if (r2 < 0d0) cycle
       if (bdist(i) >= r1) cycle
       ncand = ncand + 1
       cand(ncand) = i
       ckey(ncand) = bdist(i) - r1 ! <=0; ascending = shortest (most double) first
    end do

    if (ncand > 0) then
       ! sort candidates by distance, from most to least double-like
       allocate(iord(ncand))
       do k = 1, ncand
          iord(k) = k
       end do
       call qcksort(ckey(1:ncand),iord,1,ncand)

       ! Write down information for the candidates: which bond
       ! (cbond), endpoints (ca, cb), and weight (r - single_ref > 0).
       ! The first candidate has the shortest (most different) bond
       ! relative to the reference.
       ! remw is an accumulated deviation wrt reference:
       !   remw(k) = sum_{i = k}^ncand (ref_i - dist_i) > 0
       allocate(ca(ncand),cb(ncand),cbond(ncand),cw(ncand),remw(ncand+1))
       do k = 1, ncand
          cbond(k) = cand(iord(k))
          ca(k) = bia(cbond(k))
          cb(k) = bib(cbond(k))
          cw(k) = -ckey(iord(k)) ! = r1 - dist > 0, larger for shorter bonds
       end do
       remw(ncand+1) = 0d0
       do k = ncand, 1, -1
          remw(k) = remw(k+1) + cw(k)
       end do

       ! Greedy baseline: assign all double bonds that are currently
       ! possible in order from higest score (most different from the
       ! single reference) to lowest. This is also the fallback if
       ! the search is truncated.
       ! wopen = working iopen
       ! cursel = currently selected candidates
       ! bestsel = currently best set of selected candidates in the search
       allocate(cursel(ncand),bestsel(ncand),wopen(c%ncel))
       wopen = iopen
       cursel = .false.
       bestscore = 0d0
       do k = 1, ncand
          if (wopen(ca(k)) >= 1 .and. wopen(cb(k)) >= 1) then
             cursel(k) = .true.
             wopen(ca(k)) = wopen(ca(k)) - 1
             wopen(cb(k)) = wopen(cb(k)) - 1
             bestscore = bestscore + cw(k)
          end if
       end do
       bestsel = cursel

       ! branch-and-bound improvement over the greedy seed (bounded in both
       ! recursion depth and node count; the greedy result is kept otherwise)
       if (ncand <= ncand_max) then
          nodecount = 0
          wopen = iopen
          cursel = .false.
          call solve(1,0d0)
       end if

       ! apply the best matching found
       do k = 1, ncand
          if (bestsel(k)) then
             i = cbond(k)
             bord(i) = 2
             bfix(i) = .true.
             iopen(ca(k)) = iopen(ca(k)) - 1
             iopen(cb(k)) = iopen(cb(k)) - 1
          end if
       end do
       deallocate(iord,ca,cb,cbond,cw,remw,cursel,bestsel,wopen)
    end if

    ! flag aromatic atoms and bonds: those in a planar, conjugated sp2 ring of
    ! size 5-7. Geometric + Kekule criterion (not a Huckel 4n+2 electron count).
    allocate(arom(c%ncel),barom(nb),bdir(3,nb))
    arom = .false.
    barom = .false.
    bdir = 0d0
    call flag_aromatic(arom,barom,bdir)

    ! write the perceived orders back into both reciprocal nstar entries (using
    ! -1 = aromatic/1.5 for bonds inside an aromatic ring), and store the
    ! per-atom aromaticity flag
    do i = 1, nb
       o = bord(i)
       if (barom(i)) o = -1
       call set_bond_data(bia(i),bib(i),blv(:,i),o,bdir(:,i))
       call set_bond_data(bib(i),bia(i),-blv(:,i),o,bdir(:,i))
    end do
    do ia = 1, c%ncel
       nstar(ia)%isaromatic = arom(ia)
    end do

    deallocate(bia,bib,blv,bdist,bord,bfix,barom,bdir)
    deallocate(cap,deg,iopen,nincid,incid,cand,ckey,arom)

  contains

    !> True if the bond (ia -> ib at lattice vector lv) should be counted from
    !> ia's side, i.e. exactly once over the doubly-stored neighbor lists.
    !> Outward Cartesian vector from atom ja along bond ibnd.
    function out_vector(ja,ibnd) result(v)
      integer, intent(in) :: ja, ibnd
      real*8 :: v(3)
      if (bia(ibnd) == ja) then
         v = c%x2c(c%atcel(bib(ibnd))%x + real(blv(:,ibnd),8) - c%atcel(ja)%x)
      else
         v = c%x2c(c%atcel(bia(ibnd))%x - real(blv(:,ibnd),8) - c%atcel(ja)%x)
      end if
    end function out_vector

    !> Largest bond angle (degrees) at a degree-2 atom (its two incident bonds).
    function atom_angle(ja) result(ang)
      integer, intent(in) :: ja
      real*8 :: ang, v1(3), v2(3), n1, n2, ca
      ang = 0d0
      if (nincid(ja) < 2) return
      v1 = out_vector(ja,incid(1,ja))
      v2 = out_vector(ja,incid(2,ja))
      n1 = norm2(v1)
      n2 = norm2(v2)
      if (n1 < 1d-10 .or. n2 < 1d-10) return
      ca = dot_product(v1,v2) / (n1*n2)
      ca = max(min(ca,1d0),-1d0)
      ang = acos(ca) * 180d0 / acos(-1d0)
    end function atom_angle

    !> Pyramidality of a degree-3 atom: |triple product| of its three unit
    !> bond vectors. 0 for a planar (sp2) center, ~0.77 for ideal sp3.
    function pyramidality(ja) result(pyr)
      integer, intent(in) :: ja
      real*8 :: pyr, u1(3), u2(3), u3(3), nn
      pyr = 0d0
      if (nincid(ja) /= 3) return
      u1 = out_vector(ja,incid(1,ja))
      u2 = out_vector(ja,incid(2,ja))
      u3 = out_vector(ja,incid(3,ja))
      nn = norm2(u1)
      if (nn < 1d-10) return
      u1 = u1 / nn
      nn = norm2(u2)
      if (nn < 1d-10) return
      u2 = u2 / nn
      nn = norm2(u3)
      if (nn < 1d-10) return
      u3 = u3 / nn
      pyr = abs(dot_product(u1,cross(u2,u3)))
    end function pyramidality

    !> Assign localized X=O double bonds for an oxyacid center ia (X=N,P,S)
    !> with terminal (degree-1) oxygen neighbors.
    subroutine assign_oxyacid(ia,za)
      integer, intent(in) :: ia, za
      integer :: j, jb, no, kkk, ndbl
      integer, allocatable :: obnd(:), ord2(:)
      real*8, allocatable :: odist(:)

      ! List unfixed terminal-O bonds at this center.
      ! no = terminal unfixed O atoms, obnd(.) = bond ID, odist(.) = distance
      no = 0
      allocate(obnd(nincid(ia)),odist(nincid(ia)))
      do j = 1, nincid(ia)
         jb = incid(j,ia)
         if (bfix(jb)) cycle
         if (bia(jb) == ia) then
            if (c%spc(c%atcel(bib(jb))%is)%z /= 8) cycle
            if (deg(bib(jb)) /= 1) cycle
         else
            if (c%spc(c%atcel(bia(jb))%is)%z /= 8) cycle
            if (deg(bia(jb)) /= 1) cycle
         end if
         no = no + 1
         obnd(no) = jb
         odist(no) = bdist(jb)
      end do
      if (no == 0) then
         deallocate(obnd,odist)
         return
      end if

      ! ndbl: number of X=O double bonds to assign.
      if (za == 7) then
         ndbl = 1 ! nitro / nitrate: one localized N=O
      elseif (za == 17 .or. za == 35 .or. za == 53) then
         ! perhalate series: one X-O stays single (formal charge), rest double
         ! -> perchlorate 3, chlorate 2, chlorite 1, hypochlorite 0
         ndbl = min(no - 1,max(cap(ia) - deg(ia),0))
      else
         ndbl = max(cap(ia) - deg(ia),0) ! phosphate, sulfate
      end if
      ndbl = min(ndbl,no,max(iopen(ia),merge(1,0,za==7)))
      if (ndbl <= 0) return

      ! assign the double bonds to the shortest terminal-O bonds
      allocate(ord2(no))
      do j = 1, no
         ord2(j) = j
      end do
      call qcksort(odist(1:no),ord2,1,no)
      do kkk = 1, ndbl
         jb = obnd(ord2(kkk))
         bord(jb) = 2
         bfix(jb) = .true.
         if (bia(jb) == ia) then
            iopen(bib(jb)) = max(iopen(bib(jb))-1,0)
         else
            iopen(bia(jb)) = max(iopen(bia(jb))-1,0)
         end if
      end do
      iopen(ia) = max(iopen(ia)-ndbl,0)

    end subroutine assign_oxyacid

    !> Set the bond order (ord) and aromatic outward direction (vec) for the
    !> entry of atom ja pointing to jb at lattice vector lv.
    subroutine set_bond_data(ja,jb,lv,ord,vec)
      integer, intent(in) :: ja, jb, lv(3), ord
      real*8, intent(in) :: vec(3)
      integer :: kkk
      do kkk = 1, nstar(ja)%ncon
         if (nstar(ja)%idcon(kkk) == jb .and. all(nstar(ja)%lcon(:,kkk) == lv)) then
            nstar(ja)%ordcon(kkk) = ord
            nstar(ja)%aromdir(:,kkk) = vec
            return
         end if
      end do
    end subroutine set_bond_data

    !> Branch-and-bound search for the maximum-weight set of double
    !> bonds (candidate index idx onwards) that respects the valences
    !> in wopen. Updates bestscore/bestsel with the best assignment
    !> found.
    recursive subroutine solve(idx,score)
      integer, intent(in) :: idx
      real*8, intent(in) :: score
      integer :: pa, pb

      ! stop if the search is taking too long
      nodecount = nodecount + 1
      if (nodecount > maxnode) return

      ! stop if we are out of candidates; write it down if it is an improvement
      if (idx > ncand) then
         if (score > bestscore + 1d-12) then
            bestscore = score
            bestsel = cursel
         end if
         return
      end if

      ! return if even taking all remaining candidates cannot beat the best
      if (score + remw(idx) <= bestscore + 1d-12) return

      ! branch 1: make this bond a double (if capacity allows)
      pa = ca(idx)
      pb = cb(idx)
      if (wopen(pa) >= 1 .and. wopen(pb) >= 1) then
         wopen(pa) = wopen(pa) - 1
         wopen(pb) = wopen(pb) - 1
         cursel(idx) = .true.
         call solve(idx+1,score+cw(idx))
         wopen(pa) = wopen(pa) + 1
         wopen(pb) = wopen(pb) + 1
      end if

      ! branch 2: leave this bond single
      cursel(idx) = .false.
      call solve(idx+1,score)

    end subroutine solve

    !> Flag atoms and bonds that belong to a planar, conjugated sp2 ring of
    !> size 5-7. For each bond, the smallest ring through it is found by a
    !> periodic-aware BFS (the ring must close in real space, net lattice vector
    !> zero), then tested: every ring atom must be planar and sp2/conjugated.
    !> Sets arom(a) for the ring atoms, barom(e) for the ring bonds, and
    !> bdir(:,e) the outward (ring-exterior) cartesian direction for each ring
    !> bond (used to draw the dashed sub-bond on the inside of the ring).
    subroutine flag_aromatic(arom,barom,bdir)
      logical, intent(inout) :: arom(:)
      logical, intent(inout) :: barom(:)
      real*8, intent(inout) :: bdir(:,:)

      integer, parameter :: maxring = 7    ! largest ring size tested
      integer, parameter :: maxq = 50000   ! BFS state cap (safety guard)
      integer :: ib0, u, v, e, j, p, q, head, nq, m, kk, egoal, ne, kn, eb
      integer :: lb(3), nl(3)
      integer, allocatable :: qat(:), qlv(:,:), qpar(:), qdep(:), qedge(:)
      integer, allocatable :: ratom(:), rlv(:,:), redge(:)
      real*8 :: cen(3), rpos(3,maxring), mid(3), ax(3), out(3), an
      logical :: found, seen

      ! Up to maxq atoms in the queue
      ! qat = atom identifier
      ! qlv = accumulated lattice vector qlv
      ! qpar = parent index (root has qpar = 0)
      ! qdep = depth, number of bonds traversed
      ! qedge = bond index used to reach this state (root: undefined)
      allocate(qat(maxq),qlv(3,maxq),qpar(maxq),qdep(maxq),qedge(maxq))
      allocate(ratom(maxring),rlv(3,maxring),redge(maxring))

      ! run over bonds
      do ib0 = 1, nb
         u = bia(ib0)
         v = bib(ib0)
         lb = blv(:,ib0)

         ! BFS from (u,0) seeking (v,lb) over bonds /= ib0, ring size <= maxring
         nq = 1
         qat(1) = u
         qlv(:,1) = 0
         qpar(1) = 0
         qdep(1) = 0
         head = 0
         found = .false.

         ! run until the queue is exhausted or a ring is found
         do while (head < nq .and. .not.found)
            ! next item in the queue
            head = head + 1

            ! do not exceed max ring depth
            if (qdep(head) >= maxring-1) cycle

            ! consider all bonds incident to p. Skip the starting bond ib0 only
            ! at the root (u,0): that excludes the specific edge instance being
            ! closed, but still allows the same bond *index* at other periodic
            ! images (a periodic ring, e.g. a graphene hexagon, traverses the
            ! same translational bond twice).
            p = qat(head)
            do j = 1, nincid(p)
               e = incid(j,p)
               if (e == ib0 .and. qdep(head) == 0) cycle

               ! consider the atom at the other end, and calculate accum. lattice vector
               if (bia(e) == p) then
                  q = bib(e)
                  nl = qlv(:,head) + blv(:,e)
               else
                  q = bia(e)
                  nl = qlv(:,head) - blv(:,e)
               end if

               ! ring closes back to the other end of ib0? -> found! exit
               if (q == v .and. all(nl == lb)) then
                  found = .true.
                  egoal = e
                  exit
               end if

               ! enqueue (q,nl) if not already visited
               seen = .false.
               do kk = 1, nq
                  if (qat(kk) == q .and. all(qlv(:,kk) == nl)) then
                     seen = .true.
                     exit
                  end if
               end do
               if (.not.seen) then
                  if (nq >= maxq) cycle ! give up enqueuing (safety)
                  nq = nq + 1
                  qat(nq) = q
                  qlv(:,nq) = nl
                  qpar(nq) = head
                  qdep(nq) = qdep(head) + 1
                  qedge(nq) = e
               end if
            end do
         end do
         if (.not.found) cycle

         ! found: reconstruct the ring atoms, parent chain (p..u) reversed, then v
         ! only if 5-7 ring
         m = qdep(head) + 2
         if (m < 5 .or. m > maxring) cycle
         kk = m - 1
         ne = 0
         p = head
         do while (p > 0)
            ratom(kk) = qat(p)
            rlv(:,kk) = qlv(:,p)
            kk = kk - 1
            if (qpar(p) > 0) then ! edge from parent into p (root has none)
               ne = ne + 1
               redge(ne) = qedge(p)
            end if
            p = qpar(p)
         end do
         ratom(m) = v
         rlv(:,m) = lb
         ne = ne + 1
         redge(ne) = egoal ! closing edge p -> v
         ne = ne + 1
         redge(ne) = ib0   ! the bond we started from (v -> u)

         if (ring_is_aromatic(ratom,rlv,m)) then
            ! ring atom cartesian positions and centroid
            cen = 0d0
            do kk = 1, m
               rpos(:,kk) = c%x2c(c%atcel(ratom(kk))%x + real(rlv(:,kk),8))
               cen = cen + rpos(:,kk)
            end do
            cen = cen / real(m,8)

            do kk = 1, m
               arom(ratom(kk)) = .true.
            end do
            do kk = 1, ne
               barom(redge(kk)) = .true.
            end do

            ! outward (ring-exterior) in-plane direction for each ring bond
            ! (consecutive ring atoms, cyclic); keep the first one assigned so
            ! a bond shared by two fused rings is not overwritten
            do kk = 1, m
               kn = kk + 1
               if (kn > m) kn = 1
               mid = 0.5d0 * (rpos(:,kk) + rpos(:,kn))
               ax = rpos(:,kn) - rpos(:,kk)
               an = norm2(ax)
               if (an < 1d-10) cycle
               ax = ax / an
               out = mid - cen
               out = out - dot_product(out,ax) * ax ! perpendicular to the bond
               an = norm2(out)
               if (an < 1d-10) cycle
               out = out / an
               eb = bond_between(ratom(kk),ratom(kn),rlv(:,kn)-rlv(:,kk))
               if (eb > 0) then
                  if (norm2(bdir(:,eb)) < 1d-10) bdir(:,eb) = out
               end if
            end do
         end if
      end do

      deallocate(qat,qlv,qpar,qdep,qedge,ratom,rlv,redge)
    end subroutine flag_aromatic

    !> Return the local bond-list index connecting cell atoms a and b with
    !> relative lattice vector dl (b at lattice vector dl from a), or 0 if none.
    function bond_between(a,b,dl) result(eb)
      integer, intent(in) :: a, b, dl(3)
      integer :: eb, j, e
      eb = 0
      do j = 1, nincid(a)
         e = incid(j,a)
         if (bia(e) == a .and. bib(e) == b .and. all(blv(:,e) == dl)) then
            eb = e
            return
         elseif (bib(e) == a .and. bia(e) == b .and. all(blv(:,e) == -dl)) then
            eb = e
            return
         end if
      end do
    end function bond_between

    !> Test whether the ring with m atoms ratom(1:m) at lvecs rlv(:,1:m) is
    !> aromatic: planar (all atoms near the mean ring plane) and every atom sp2
    !> and contributing to the pi system (a double bond, or a lone-pair donor
    !> heteroatom: pyrrole-type N deg 3, or furan/thiophene O/S deg 2).
    function ring_is_aromatic(ratom,rlv,m) result(arom1)
      use tools_math, only: plane_from_points
      use param, only: bohrtoa
      integer, intent(in) :: m, ratom(m), rlv(3,m)
      logical :: arom1

      real*8, parameter :: dplane = 0.20d0 / bohrtoa ! max out-of-plane deviation (ang)
      real*8 :: pos(3,m), cen(3), nor(3), dev
      integer :: k, j, a, z
      logical :: hasdbl

      arom1 = .false.

      ! Cartesian positions of the ring atom images
      do k = 1, m
         pos(:,k) = c%x2c(c%atcel(ratom(k))%x + real(rlv(:,k),8))
      end do

      ! centroid, best-fit ring-plane normal, and planarity test: the maximum
      ! out-of-plane deviation must be within dplane
      call plane_from_points(pos,m,cen,nor,dev)
      if (norm2(nor) < 1d-10) return ! degenerate (collinear/coincident)
      if (dev > dplane) return       ! not planar

      ! conjugation: every ring atom must be sp2 / contribute to the pi system
      do k = 1, m
         a = ratom(k)
         hasdbl = .false.
         do j = 1, nincid(a)
            if (bord(incid(j,a)) == 2) then
               hasdbl = .true.
               exit
            end if
         end do
         if (hasdbl) cycle
         z = c%spc(c%atcel(a)%is)%z
         if (z == 7 .and. deg(a) == 3) cycle             ! pyrrole-type N
         if ((z == 8 .or. z == 16) .and. deg(a) == 2) cycle ! furan O / thiophene S
         return ! this atom is not part of a closed pi system -> not aromatic
      end do

      arom1 = .true.
    end function ring_is_aromatic

  end subroutine perceive_bond_orders

  !> Maximum valence (capacity) for the supported organic elements; returns
  !> a negative value for unsupported elements (metals, etc.).
  function zhang_maxval(z) result(mv)
    integer, intent(in) :: z
    integer :: mv
    select case (z)
    case (1)  ! H
       mv = 1
    case (6)  ! C
       mv = 4
    case (7)  ! N
       mv = 3
    case (8)  ! O
       mv = 2
    case (9)  ! F
       mv = 1
    case (15) ! P
       mv = 5
    case (16) ! S
       mv = 6
    case (17) ! Cl
       mv = 7
    case (35) ! Br
       mv = 7
    case (53) ! I
       mv = 7
    case default
       mv = -1
    end select
    ! Note: Cl/Br/I are given their hypervalent (oxyacid) capacity 7 rather than
    ! the halide valence 1. This is safe because there are no halogen
    ! multiple-bond references in zhang_refbond, so the triple/matching stages
    ! never form a halogen multiple bond; the only consumer of this capacity is
    ! the acid model (assign_oxyacid), so an ordinary halide (C-Cl, ...) still
    ! stays single. F (mv=1) is never hypervalent.
    ! The capacities here are standard textbook valences, NOT from Zhang et al.;
    ! the paper only specifies max connections = 4 for C/N/P/S and refers to a
    ! "maximum valence" (Eq 7) without tabulating values.
  end function zhang_maxval

  !> Reference single (r1), double (r2) and triple (r3) bond lengths (bohr) for
  !> an element pair. A negative value means no reference for that order.
  !> Sources (values in angstrom in the code, converted to bohr at the end):
  !>  - The six single/double pairs C-C, C-N, C-O, C-S, N-N, N-O are from
  !>    Zhang et al., J. Cheminform. 2012, 4:26, Table 3.
  !>  - Triple bonds and the C-P, N-S pairs are sums of Pyykko & Atsumi additive
  !>    covalent radii (single/double: Chem. Eur. J. 2009, 15, 186 and 12770;
  !>    triple: Pyykko, Riedel, Patzschke, Chem. Eur. J. 2005, 11, 3511), in pm:
  !>    C 75/67/60, N 71/60/54, O 63/57/53, P 111/102/94, S 103/94/95.
  !>  - O-O, O-S, O-P use observed values (additive radii deviate for these
  !>    polar/pi bonds): O2 1.208 / H2O2 1.475; sulfone S=O ~1.44, S-O ~1.57;
  !>    phosphate P-O ~1.62, P=O ~1.50.
  subroutine zhang_refbond(z1,z2,r1,r2,r3)
    use param, only: bohrtoa
    integer, intent(in) :: z1, z2
    real*8, intent(out) :: r1, r2, r3
    integer :: za, zb

    za = min(z1,z2)
    zb = max(z1,z2)
    r1 = -1d0
    r2 = -1d0
    r3 = -1d0
    if (za == 6 .and. zb == 6) then       ! C-C
       r1 = 1.49d0
       r2 = 1.33d0
       r3 = 1.20d0 ! Pyykko (0.60+0.60)
    elseif (za == 6 .and. zb == 7) then   ! C-N
       r1 = 1.43d0
       r2 = 1.29d0
       r3 = 1.14d0 ! Pyykko (0.60+0.54)
    elseif (za == 6 .and. zb == 8) then   ! C-O
       r1 = 1.41d0
       r2 = 1.28d0
       r3 = 1.13d0 ! Pyykko (0.60+0.53)
    elseif (za == 6 .and. zb == 16) then  ! C-S
       r1 = 1.80d0
       r2 = 1.72d0
    elseif (za == 6 .and. zb == 15) then  ! C-P
       r1 = 1.86d0 ! Pyykko (0.75+1.11)
       r2 = 1.69d0 ! Pyykko (0.67+1.02)
    elseif (za == 7 .and. zb == 7) then   ! N-N
       r1 = 1.38d0
       r2 = 1.29d0
       r3 = 1.08d0 ! Pyykko (0.54+0.54)
    elseif (za == 7 .and. zb == 8) then   ! N-O
       r1 = 1.39d0
       r2 = 1.24d0
    elseif (za == 8 .and. zb == 8) then   ! O-O
       r1 = 1.48d0 ! obs, H2O2 peroxide
       r2 = 1.21d0 ! obs, O2
    elseif (za == 8 .and. zb == 16) then  ! O-S
       r1 = 1.57d0 ! obs, S-O single
       r2 = 1.44d0 ! obs, sulfone S=O
    elseif (za == 8 .and. zb == 15) then  ! O-P
       r1 = 1.62d0 ! obs, phosphate P-O
       r2 = 1.50d0 ! obs, P=O
    elseif (za == 7 .and. zb == 16) then  ! N-S
       r1 = 1.74d0 ! Pyykko (0.71+1.03)
       r2 = 1.54d0 ! Pyykko (0.60+0.94)
    end if
    if (r1 > 0d0) r1 = r1 / bohrtoa
    if (r2 > 0d0) r2 = r2 / bohrtoa
    if (r3 > 0d0) r3 = r3 / bohrtoa
  end subroutine zhang_refbond

end submodule env
