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
    integer :: i, j, k, ix(3), i0(3), i1(3), idx, nn
    integer :: ib(3), ithis(3), nsafe, up2n_, nseen
    integer, allocatable :: at_id(:), at_lvec(:,:)
    real*8, allocatable :: at_dist(:), rshel(:)
    integer :: nb, nshellb
    integer, allocatable :: idb(:,:), iaux(:), iord(:)
    integer :: nshel
    real*8 :: up2d_2
    real*8,allocatable :: up2dsp_2(:,:), up2dcidx_2(:)

    real*8, parameter :: eps = 1d-20
    integer, parameter :: ixmol_cut = 10 ! cutoff for switching to no-block molecule dist search
    real*8, parameter :: shell_eps = 1d-5

    ! process input
    nozero_ = .false.
    sorted_ = sorted .or. present(up2n) .or. present(up2sh) .or. present(ishell0)
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
          ! Run over shells of nearby blocks and find atoms until we
          ! have at least the given number of shells
          nseen = 0
          nshellb = 0
          do while(.true.)
             call make_block_shell(c,nshellb,nb,idb,dmax)

             do i = 1, nb
                ithis = idb(:,i) + ix

                if (c%ismolecule) then
                   if (ithis(1) < 0 .or. ithis(1) >= c%nblock(1) .or.&
                      ithis(2) < 0 .or. ithis(2) >= c%nblock(2) .or.&
                      ithis(3) < 0 .or. ithis(3) >= c%nblock(3)) cycle
                   ib = ithis
                   lvecx = 0
                   xdif = -xorigc
                else
                   ib = (/modulo(ithis(1),c%nblock(1)), modulo(ithis(2),c%nblock(2)), modulo(ithis(3),c%nblock(3))/)
                   lvecx = real((ithis - ib) / c%nblock,8)
                   lvecx = matmul(c%m_xr2c,lvecx)
                   xdif = lvecx - xorigc
                end if

                ! run over atoms in this block
                idx = c%iblock0(ib(1),ib(2),ib(3))
                do while (idx /= 0)
                   ! increment the number of "seen" atoms
                   nseen = nseen + 1

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
                         if (c%ismolecule) up2n_ = min(up2n,c%ncel-1)
                      end if
                      if (ok) then
                         call add_atom_to_output_list()
                         dsqrt = at_dist(nat)
                         if (present(up2sh)) &
                            call add_shell_to_output_list()
                      end if
                   end if

                   ! next atom
                   idx = c%atcel(idx)%inext
                end do ! while idx
             end do

             if (present(up2sh)) then
                ! exit if we have enough shells
                nsafe = count(rshel(1:nshel) <= dmax)
                if (nsafe >= up2sh) exit
             else
                ! check if we have enough atoms
                nsafe = count(at_dist(1:nat) <= dmax)
                if (nsafe >= up2n_) exit
             end if

             ! exit if we have examined all atoms in the molecule
             if (c%ismolecule .and. nsafe == c%ncel) then
                up2n_ = min(up2n_,nat)
                exit
             end if

             ! next shell
             nshellb = nshellb + 1
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
                      if (ok) call add_atom_to_output_list()
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
      at_dist(nat) = sqrt(dd)
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

    ! get just one atom or all sorted atoms up to distmax
    if (present(distmax)) then
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
  !> c%nstar array. Two atoms i and j are bonded if their distance is
  !> less than bondfac * (atmrad(i) + atmrad(j))
  module subroutine find_asterisms(c,nstar,atmrad,bondfac,rij)
    use tools_io, only: string, ferror, faterr
    use types, only: realloc, celatom
    use param, only: maxzat0
    class(crystal), intent(inout) :: c
    type(neighstar), allocatable, intent(inout) :: nstar(:)
    real*8, intent(in), optional :: atmrad(0:maxzat0)
    real*8, intent(in), optional :: bondfac
    real*8, intent(in), optional :: rij(:,:,:)

    integer :: i, j, nx(3), i0shift(3), i1shift(3)
    real*8 :: ri, rj, dmax, dd, xi(3), xj(3), xdelta(3)
    real*8, allocatable :: rij2(:,:,:)
    integer :: iblock_stride, nconi, nconj
    integer :: ix, iy, iz, jx, jy, jz, iid, jid, iidx(3), jidx(3)
    integer :: iat, jat, iaux(3), lvecx(3)
    logical :: sameblock

    ! return if there are no atoms
    if (c%ncel == 0) return

    ! allocate the asterism arrays
    if (allocated(nstar)) deallocate(nstar)
    allocate(nstar(c%ncel))
    do i = 1, c%ncel
       nstar(i)%ncon = 0
       if (allocated(nstar(i)%idcon)) deallocate(nstar(i)%idcon)
       if (allocated(nstar(i)%lcon)) deallocate(nstar(i)%lcon)
       allocate(nstar(i)%idcon(20))
       allocate(nstar(i)%lcon(3,20))
    end do

    ! pre-calculate the distance matrix
    allocate(rij2(c%nspc,2,c%nspc))
    rij2 = 0d0
    dmax = 0d0
    if (present(atmrad).and.present(bondfac)) then
       ! from atmrad and bondfac
       do i = 1, c%nspc
          if (c%spc(i)%z <= 0) cycle
          ri = atmrad(c%spc(i)%z)
          do j = 1, c%nspc
             if (c%spc(j)%z <= 0) cycle
             rj = atmrad(c%spc(j)%z)

             rij2(j,2,i) = (ri+rj) * bondfac
             rij2(j,1,i) = 0d0
             rij2(j,2,i) = rij2(j,2,i) * rij2(j,2,i)
             rij2(j,1,i) = rij2(j,1,i) * rij2(j,1,i)
             dmax = max(dmax,(ri+rj) * bondfac)
          end do
       end do
    elseif (present(rij)) then
       ! from rij
       if (any(shape(rij) /= shape(rij2))) &
          call ferror('find_asterisms','error in rij shape',faterr)
       rij2 = rij * rij
       dmax = maxval(rij)
    else
       call ferror('find_asterisms','error in input arguments',faterr)
    end if

    ! calculate the number of blocks in each direction required for satifying
    ! that the largest sphere in the super-block has radius > dmax.
    ! r = Vblock / 2 / max(cv(3)/n3,cv(2)/n2,cv(1)/n1)
    nx = ceiling(c%blockcv / (0.5d0 * c%blockomega / max(dmax,1d-40)))

    ! define the search space
    do i = 1, 3
       if (mod(nx(i),2) == 0) then
          i0shift(i) = - nx(i)/2
          i1shift(i) = nx(i)/2
       else
          i0shift(i) = - (nx(i)/2+1)
          i1shift(i) = + (nx(i)/2+1)
       end if
    end do

    iblock_stride = maxval(c%nblock) + maxval(abs(i1shift)) + maxval(abs(i0shift)) + 10

    do ix = 0, c%nblock(1)-1
       do iy = 0, c%nblock(2)-1
          do iz = 0, c%nblock(3)-1
             if (c%iblock0(ix,iy,iz) == 0) cycle
             iid = (ix * iblock_stride + iy) * iblock_stride + iz
             iidx = (/ix,iy,iz/)

             do jx = ix + i0shift(1), ix + i1shift(1)
                do jy = iy + i0shift(2), iy + i1shift(2)
                   do jz = iz + i0shift(2), iz + i1shift(3)
                      ! skip if this block pair has been considered already, or outside the region
                      ! skip if there are no atoms
                      if (jx >= 0 .and. jx < c%nblock(1) .and.&
                          jy >= 0 .and. jy < c%nblock(2) .and.&
                          jz >= 0 .and. jz < c%nblock(3)) then
                         jid = (jx * iblock_stride + jy) * iblock_stride + jz
                         if (jid < iid) cycle
                         jidx = (/jx,jy,jz/)
                         if (c%iblock0(jidx(1),jidx(2),jidx(3)) == 0) cycle
                         xdelta = 0d0
                      elseif (c%ismolecule) then
                         cycle
                      else
                         iaux = (/jx,jy,jz/)
                         jidx = (/modulo(jx,c%nblock(1)), modulo(jy,c%nblock(2)), modulo(jz,c%nblock(3))/)
                         jid = (jidx(1) * iblock_stride + jidx(2)) * iblock_stride + jidx(3)
                         if (jid < iid) cycle
                         xdelta = real((iaux - jidx) / c%nblock,8)
                         xdelta = matmul(c%m_xr2c,xdelta)
                      end if
                      if (c%iblock0(jidx(1),jidx(2),jidx(3)) == 0) cycle
                      sameblock = (ix == jx) .and. (iy == jy) .and. (iz == jz)

                      jat = c%iblock0(jidx(1),jidx(2),jidx(3))
                      do while (jat /= 0)
                         xj = c%atcel(jat)%rxc + xdelta

                         if (sameblock) then
                            iat = c%atcel(jat)%inext
                         else
                            iat = c%iblock0(iidx(1),iidx(2),iidx(3))
                         end if
                         do while (iat /= 0)
                            xi = c%atcel(iat)%rxc

                            dd = dot_product(xi-xj,xi-xj)
                            if (dd > rij2(c%atcel(jat)%is,1,c%atcel(iat)%is) .and. &
                               dd < rij2(c%atcel(jat)%is,2,c%atcel(iat)%is)) then
                               ! add this atom
                               nstar(iat)%ncon = nstar(iat)%ncon + 1
                               nconi = nstar(iat)%ncon
                               if (nconi > size(nstar(iat)%idcon,1)) then
                                  call realloc(nstar(iat)%idcon,2*nconi)
                                  call realloc(nstar(iat)%lcon,3,2*nconi)
                               end if
                               nstar(iat)%idcon(nconi) = jat

                               nstar(jat)%ncon = nstar(jat)%ncon + 1
                               nconj = nstar(jat)%ncon
                               if (nconj > size(nstar(jat)%idcon,1)) then
                                  call realloc(nstar(jat)%idcon,2*nconj)
                                  call realloc(nstar(jat)%lcon,3,2*nconj)
                               end if
                               nstar(jat)%idcon(nconj) = iat
                               if (c%ismolecule) then
                                  nstar(iat)%lcon(:,nconi) = 0
                                  nstar(jat)%lcon(:,nconj) = 0
                               else
                                  lvecx = nint(matmul(c%m_c2x,xj-xi) + c%atcel(iat)%x - c%atcel(jat)%x)
                                  nstar(iat)%lcon(:,nconi) = lvecx
                                  nstar(jat)%lcon(:,nconj) = -lvecx
                               end if
                            end if

                            iat = c%atcel(iat)%inext
                         end do
                         jat = c%atcel(jat)%inext
                      end do
                   end do
                end do
             end do

          end do
       end do
    end do

    ! reallocate
    do i = 1, c%ncel
       call realloc(nstar(i)%idcon,nstar(i)%ncon)
       call realloc(nstar(i)%lcon,3,nstar(i)%ncon)
    end do

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
  module subroutine list_near_lattice_points(c,xp,icrd,sorted,nat,dist,lvec,ndiv,&
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
    real*8 :: x2c(3,3), x2xr(3,3), xr2x(3,3), xr2c(3,3), c2xr(3,3)

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

    if (present(ndiv)) then
       do i = 1, 3
          x2c(:,i) = c%m_x2c(:,i) / real(ndiv(i),8)
       end do
       x0 = xp * ndiv
       call wscell(x2c,.true.,m_xr2x=xr2x,m_xr2c=xr2c,m_x2xr=x2xr,m_c2xr=c2xr)
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
  !> inscribe a sphere of radius rmax.
  !>
  !> Calculation of distance between point inside the block and its six
  !> sides, assuming the corners of the block have (0,0,0) and (1,1,1)
  !> coordinates:
  !> d1l = |(x-x0) * (a x b)| / ||a x b|| = |x3|   * V / ||a x b||
  !> d1u = |(x-x0) * (a x b)| / ||a x b|| = |x3-1| * V / ||a x b||
  !> d2l = |(x-x0) * (b x c)| / ||b x c|| = |x1|   * V / ||b x c||
  !> d2u = |(x-x0) * (b x c)| / ||b x c|| = |x1-1| * V / ||b x c||
  !> d3l = |(x-x0) * (c x a)| / ||c x a|| = |x2|   * V / ||c x a||
  !> d3u = |(x-x0) * (c x a)| / ||c x a|| = |x2-1| * V / ||c x a||
  subroutine blocks_inscribed_sphere(c,x,rmax,i0,i1)
    type(crystal), intent(in) :: c
    real*8, intent(in) :: x(3)
    real*8, intent(in) :: rmax
    integer :: i0(3), i1(3)

    integer :: ix(3), i
    real*8 :: dist, xn(3)
    logical :: again

    ! find block coordinates and the home block
    xn = x * c%nblock
    ix = floor(xn)
    i0 = ix
    i1 = ix

    ! expand the blocks in each direction until a sphere with radius rmax
    ! fits in the region
    again = .true.
    do while (again)
       again = .false.
       do i = 1, 3
          dist = (xn(i) - i0(i)) * c%blockomega / c%blockcv(i)
          if (dist < rmax) then
             i0(i) = i0(i) - 1
             again = .true.
          end if
          dist = (i1(i) + 1 - xn(i)) * c%blockomega / c%blockcv(i)
          if (dist < rmax) then
             i1(i) = i1(i) + 1
             again = .true.
          end if
       end do
    end do

  end subroutine blocks_inscribed_sphere

end submodule env
