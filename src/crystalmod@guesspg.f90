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

! Old (pre-spglib) in-house crystal symmetry guesser. Used instead of
! spglib when the SYM OLD keyword is given before a structure is read.
! Based on the space-group operations guessing algorithm by teVelde
! (described in his PhD thesis).
submodule (crystalmod) guesspg
  implicit none

  ! private state for the old (pre-spglib) symmetry guesser (guessspg)
  real*8, parameter :: pusheps = 1d-14 !< Reduce atom coords. to (-pusheps,1-pusheps]

contains

  !> Calculate the crystal symmetry operations using the old (pre-spglib)
  !> guesser. Mirrors calcsym (usenneq=.true.) but uses guessspg instead of
  !> spglib: it does not provide the space-group number/symbol, Hall symbol,
  !> Wyckoff letters, or conventional-cell standardization (spgavail stays
  !> .false.). On input c%nneq/c%at hold the as-read atom basis. On output
  !> neqv, rotm, ncv, cen, at()%mult, the reduced at() list, and ncel/atcel
  !> are filled. If error, errmsg has length > 0.
  module subroutine calcsym_old(c,errmsg)
    use param, only: icrd_crys
    use global, only: symprec
    use types, only: realloc
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, k, io, it, id
    real*8 :: x0(3)
    character*10, allocatable :: name(:)
    logical, allocatable :: done(:)

    errmsg = ""
    c%spgavail = .false.

    if (c%nneq == 0) then
       errmsg = "no atoms to determine the symmetry"
       return
    end if

    ! Materialize the complete cell list (a copy of the as-read atoms) and
    ! build the block spatial index, so that locating the atom at a given
    ! position becomes an O(1) operation (c%identify_atom) instead of an
    ! O(N) scan. All the guesser routines below use this index.
    c%ncel = c%nneq
    allocate(name(c%ncel))
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
       c%atcel(i)%idx = i
       name(i) = c%at(i)%name
    end do
    call c%build_env()

    ! find the symmetry operations (neqv, rotm, ncv, cen) from the geometry
    c%havesym = 0
    call c%guessspg(2)

    ! Orbit reduction: build the non-equivalent atom list, the
    ! multiplicities, and the complete-list symmetry mapping in a single
    ! O(ncel*neqv*ncv) pass. Each atom is assigned to the first
    ! representative whose orbit reaches it (same convention as the
    ! previous reduceatoms+symeqv), using O(1) lookups.
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
    deallocate(done,name)

    c%havesym = 1

    ! no spglib dataset available: clear the Wyckoff letters
    call c%spgtowyc()

  end subroutine calcsym_old

  !xx! guessspg - crystal symmetry guessing module (old, pre-spglib)
  !> Guesses the symmetry operations from the geometry of the unit
  !> cell and the positions of the atoms in it. In: cell vectors
  !> (m_x2c). Inout: nneq, at(:). Out: neqv, ncv, cen, rotm. If level
  !> = 0, use no symmetry. If level = 1, find only the centering
  !> vectors. Level = 2, full symmetry.
  module subroutine guessspg(c,level)
    use param, only: eyet
    class(crystal), intent(inout) :: c
    integer, intent(in) :: level

    real*8 :: rmat(3,3), rot(3,3,48)
    integer :: isref, iref

    ! no symmetry
    if (level == 0 .and. c%havesym == 0) then
       c%neqv = 1
       c%rotm(:,:,1) = eyet
       c%ncv = 1
       if (allocated(c%cen)) deallocate(c%cen)
       allocate(c%cen(3,1))
       c%cen = 0d0
       c%havesym = 0
       return
    end if

    ! check if we already have this level
    if (c%havesym >= level) return

    ! write down the new level
    c%havesym = level

    ! Reference species (fewest atoms) used as the anchor for the
    ! centering search and the translation search, to minimize the number
    ! of candidate vectors that have to be tested.
    call refspecies(c,isref,iref)

    ! Find the centering group by transforming to a primitive cell
    call cenbuild(c,isref,iref)

    ! Find the non-centering operations?
    if (level > 1) then
       ! Find the rotation parts of the operations
       rmat = transpose(c%m_x2c)
       call lattpg(rmat,c%neqv,rot)
       c%rotm(1:3,1:3,1:c%neqv) = rot(:,:,1:c%neqv)

       ! the first operation is the identity, with zero translation
       ! (lattpg only sets the rotation parts, and filltrans skips op 1)
       c%rotm(:,:,1) = eyet

       ! Calculate the translation vectors
       call filltrans(c,isref,iref)
    end if

  end subroutine guessspg

  !> Choose the reference species (the one with the fewest atoms in the
  !> complete cell list) and the index of its first atom. Used to anchor
  !> the centering and translation searches with as few candidate vectors
  !> as possible. Requires the complete atom list (atcel/ncel).
  subroutine refspecies(c,isref,iref)
    type(crystal), intent(inout) :: c
    integer, intent(out) :: isref !< reference species
    integer, intent(out) :: iref !< first atom (atcel index) of that species

    integer :: i, nmin
    integer, allocatable :: nsp(:)

    allocate(nsp(c%nspc))
    nsp = 0
    do i = 1, c%ncel
       nsp(c%atcel(i)%is) = nsp(c%atcel(i)%is) + 1
    end do
    isref = 1
    nmin = huge(1)
    do i = 1, c%nspc
       if (nsp(i) > 0 .and. nsp(i) < nmin) then
          nmin = nsp(i)
          isref = i
       end if
    end do
    deallocate(nsp)

    iref = 1
    do i = 1, c%ncel
       if (c%atcel(i)%is == isref) then
          iref = i
          exit
       end if
    end do

  end subroutine refspecies

  !> Find all centering vectors in the crystal and write them to c%cen
  !> (c%ncv of them). A centering vector maps the reference atom (iref,
  !> species isref) onto another atom of the same species, so the
  !> candidates are the difference vectors from iref to every other atom
  !> of the reference species; each is validated with iscelltr. The
  !> validated candidates already form the complete (closed) centering
  !> group, so no group-closure step is needed. Cost: O(m*ncel) with O(1)
  !> lookups, m = number of reference-species atoms (true positives run a
  !> full O(ncel) check, false candidates exit early).
  subroutine cenbuild(c,isref,iref)
    use global, only: symprec
    use types, only: realloc
    type(crystal), intent(inout) :: c
    integer, intent(in) :: isref, iref

    real*8  :: tr(3)
    integer :: j, k
    logical :: rep

    ! the trivial centering vector
    if (allocated(c%cen)) deallocate(c%cen)
    allocate(c%cen(3,4))
    c%ncv = 1
    c%cen = 0d0

    do k = 1, c%ncel
       if (c%atcel(k)%is /= isref) cycle
       if (k == iref) cycle
       tr = c%atcel(k)%x - c%atcel(iref)%x
       call reduce(tr)

       ! skip if already in the list (ncv is small)
       rep = .false.
       do j = 1, c%ncv
          if (c%are_lclose(c%cen(:,j),tr,symprec)) then
             rep = .true.
             exit
          end if
       end do
       if (rep) cycle

       ! a true centering vector maps every atom onto a same-species atom
       if (iscelltr(c,tr)) then
          c%ncv = c%ncv + 1
          if (c%ncv > size(c%cen,2)) call realloc(c%cen,3,2*c%ncv)
          c%cen(:,c%ncv) = tr
       end if
    end do
    call realloc(c%cen,3,c%ncv)

  end subroutine cenbuild

  !> Given a vector in crystallographic coordinates, reduce
  !> it to the main cell (xi \in (-pusheps,1-pusheps]).
  subroutine reduce(tr)
    real*8, intent(inout) :: tr(3) !< Input/output vector to reduce

    integer :: i

    do i = 1, 3
       if (tr(i) .gt. -pusheps) then
          tr(i) = tr(i) - int(tr(i)+pusheps)
       else
          tr(i) = tr(i) - int(tr(i)+pusheps) + 1d0
       end if
    end do

  end subroutine reduce

  !> Returns .true. if the translation vector tr is a symmetry operation
  !> by itself, i.e. it maps every atom in the complete cell list onto an
  !> atom of the same species. Uses the O(1) atom locator (identify_atom),
  !> so the test is O(ncel) and exits early on the first atom with no
  !> image. This is part of the spg operations guessing algorithm by
  !> teVelde, described in his PhD thesis.
  function iscelltr(c,tr)
    use param, only: icrd_crys
    use global, only: symprec
    logical :: iscelltr
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: tr(3) !< Cell translation vector to check

    integer :: i, id

    iscelltr = .false.
    do i = 1, c%ncel
       id = c%identify_atom(c%atcel(i)%x + tr,icrd_crys,distmax=symprec)
       if (id == 0) return
       if (c%atcel(id)%is /= c%atcel(i)%is) return
    end do
    iscelltr = .true.

  end function iscelltr

  !> Given the rotational parts of the space group operators
  !> (rotm(1:3,1:3,i)), fill their translational part (rotm(:,4,i)).
  !> The translation of a valid operation maps the reference atom (iref,
  !> species isref) onto a same-species atom, so the candidate
  !> translations are obtained from iref's image; each is validated with
  !> goodop. Operations with no valid translation are discarded. This is
  !> part of the spg operations guessing algorithm by teVelde, described
  !> in his PhD thesis.
  subroutine filltrans(c,isref,iref)
    use global, only: symprec
    type(crystal), intent(inout) :: c
    integer, intent(in) :: isref, iref

    integer :: op, j, k, n
    real*8  :: xnew(3)
    logical :: saveop(48), doj, dok

    ! skip identity
    saveop(1) = .true.
    do op = 2, c%neqv
       saveop(op) = .false.
       ! the transformed reference atom
       xnew = matmul(c%rotm(1:3,1:3,op),c%atcel(iref)%x)
       doj = .true.
       j = 0
       do while (doj .and. j .lt. c%ncel)
          j = j + 1
          if (c%atcel(j)%is .eq. isref) then
             !.propose a translation vector (reference atom -> atom j)
             c%rotm(:,4,op) = c%atcel(j)%x - xnew - floor(c%atcel(j)%x - xnew)
             !.if it is truly an operation, stop searching
             if (goodop(c,op)) then
                saveop(op) = .true.
                doj = .false.
             end if
          end if
       end do
    end do

    ! rewrite the list of operations
    n = 0
    do op = 1, c%neqv
       if (saveop(op)) then
          n = n + 1
          c%rotm(:,:,n) = c%rotm(:,:,op)
          call reduce(c%rotm(:,4,n))
          dok = .true.
          k = 1
          do while (dok .and. k<c%ncv)
             k = k + 1
             if (c%are_lclose(c%rotm(:,4,n),c%cen(:,k),symprec)) then
                c%rotm(:,4,n) = 0d0
                dok = .false.
             end if
          end do
       end if
    end do
    c%neqv = n

  end subroutine filltrans

  !> Check if a symmetry operation (rotm(:,:,op)) is consistent with the
  !> atomic positions, i.e. it maps every atom in the complete cell list
  !> onto a same-species atom. Uses the O(1) atom locator
  !> (identify_atom), so the test is O(ncel) and exits early. This is
  !> part of the spg operations guessing algorithm by teVelde, described
  !> in his PhD thesis.
  function goodop(c,op)
    use param, only: icrd_crys
    use global, only: symprec
    logical :: goodop
    type(crystal), intent(inout) :: c
    integer, intent(in) :: op !< Operation identifier

    integer :: i, id
    real*8  :: xnew(3)

    goodop = .false.
    do i = 1, c%ncel
       xnew = matmul(c%rotm(1:3,1:3,op),c%atcel(i)%x) + c%rotm(:,4,op)
       id = c%identify_atom(xnew,icrd_crys,distmax=symprec)
       if (id == 0) return
       if (c%atcel(id)%is /= c%atcel(i)%is) return
    end do
    goodop = .true.

  end function goodop

  !> Determine the rotational parts of the symmetry operations from
  !> the lattice point group, using the lattice metric tensor. A
  !> crystallographic rotation is an integer matrix R (in the
  !> crystallographic basis) with |det R| = 1 that preserves the
  !> metric tensor G = rmat.rmat^T to within the symmetry
  !> tolerance. The columns of R are the images of the basis vectors,
  !> so each column is an integer lattice vector with the same length
  !> as the corresponding basis vector. rmat (rows = lattice vectors)
  !> is the input; nn (number of operations) and rot (rotation
  !> matrices, crystallographic, identity first) are the output.
  subroutine lattpg(rmat,nn,rot)
    use tools_math, only: det3, matinv
    use global, only: symprec
    use types, only: realloc
    use param, only: eye
    real*8, intent(in) :: rmat(3,3)
    integer, intent(out) :: nn
    real*8, intent(out) :: rot(3,3,48)

    real*8 :: gmat(3,3), ginv(3,3), aal(3), maxl, gp(3,3), rop(3,3), d
    integer :: i, nbox(3), ia, ib, ic, ncand, n1, n2, n3, j1, j2, j3, ier
    integer, allocatable :: cand(:,:), idx1(:), idx2(:), idx3(:)
    real*8, allocatable :: clen(:)

    gmat = matmul(rmat,transpose(rmat))
    do i = 1, 3
       aal(i) = sqrt(gmat(i,i))
    end do
    maxl = maxval(aal)

    ! Box of integer lattice vectors that can have Cartesian length <=
    ! maxl+symprec. The i-th integer index of a vector v is n_i = b_i.v
    ! (b_i = reciprocal lattice vector), so |n_i| <= |v| |b_i| = |v|
    ! sqrt(ginv_ii). This bound is correct for any (incl. sheared,
    ! non-reduced) cell, unlike maxl/aal(i).
    ginv = gmat
    call matinv(ginv,3,ier)
    do i = 1, 3
       nbox(i) = int((maxl+symprec)*sqrt(ginv(i,i))) + 1
    end do

    ! collect candidate columns: integer lattice vectors with Cartesian
    ! length <= maxl + symprec (the columns of any rotation are images of
    ! the basis vectors, so no longer than the longest basis vector)
    ncand = 0
    allocate(cand(3,32),clen(32))
    do ia = -nbox(1), nbox(1)
       do ib = -nbox(2), nbox(2)
          do ic = -nbox(3), nbox(3)
             if (ia == 0 .and. ib == 0 .and. ic == 0) cycle
             d = dot_product(real((/ia,ib,ic/),8),matmul(gmat,real((/ia,ib,ic/),8))) ! length^2
             if (d <= (maxl+symprec)**2) then
                ncand = ncand + 1
                if (ncand > size(clen)) then
                   call realloc(cand,3,2*ncand)
                   call realloc(clen,2*ncand)
                end if
                cand(:,ncand) = (/ia,ib,ic/)
                clen(ncand) = sqrt(max(d,0d0))
             end if
          end do
       end do
    end do

    ! candidate columns matching each basis vector's length
    allocate(idx1(ncand),idx2(ncand),idx3(ncand))
    n1 = 0; n2 = 0; n3 = 0
    do i = 1, ncand
       if (abs(clen(i) - aal(1)) <= symprec) then; n1 = n1+1; idx1(n1) = i; end if
       if (abs(clen(i) - aal(2)) <= symprec) then; n2 = n2+1; idx2(n2) = i; end if
       if (abs(clen(i) - aal(3)) <= symprec) then; n3 = n3+1; idx3(n3) = i; end if
    end do

    ! the identity is always present; place it first (guessspg relies on
    ! op 1 = identity)
    nn = 1
    rot(:,:,1) = eye

    ! enumerate triples of candidate columns (one matching each basis
    ! vector's length) forming an integer matrix R with |det| = 1 that
    ! preserves the metric (R^T G R = G)
    do j1 = 1, n1
       rop(:,1) = real(cand(:,idx1(j1)),8)
       do j2 = 1, n2
          rop(:,2) = real(cand(:,idx2(j2)),8)
          do j3 = 1, n3
             rop(:,3) = real(cand(:,idx3(j3)),8)
             if (abs(abs(det3(rop)) - 1d0) > 1d-6) cycle

             ! skip the identity (already placed first)
             if (all(abs(rop - eye) < 1d-6)) cycle

             ! the rotated metric R^T G R must match G
             gp = matmul(transpose(rop),matmul(gmat,rop))
             if (.not.metric_ok(gp,gmat,symprec)) cycle

             nn = nn + 1
             if (nn > 48) then
                ! a crystallographic point group has at most 48 operations;
                ! guard against overflow (should not happen for a lattice)
                nn = 48
                return
             end if
             rot(:,:,nn) = rop
          end do
       end do
    end do

  contains
    !> Port of spglib's is_identity_metric: are two metric tensors
    !> equal within a Cartesian distance tolerance eps? Lengths are
    !> compared directly; angle differences are recast as
    !> displacements (length * sin(dtheta)) and compared to eps.
    function metric_ok(g1,g2,eps) result(ok)
      real*8, intent(in) :: g1(3,3), g2(3,3), eps
      logical :: ok

      integer :: ii, jj, kk
      integer, parameter :: es(2,3) = reshape((/1,2, 1,3, 2,3/),(/2,3/))
      real*8 :: l1(3), l2(3), cos1, cos2, x, sind2, lave2

      ok = .false.
      do ii = 1, 3
         l1(ii) = sqrt(g1(ii,ii))
         l2(ii) = sqrt(g2(ii,ii))
         if (abs(l1(ii) - l2(ii)) > eps) return
      end do
      do ii = 1, 3
         jj = es(1,ii)
         kk = es(2,ii)
         cos1 = g2(jj,kk) / (l2(jj)*l2(kk))
         cos2 = g1(jj,kk) / (l1(jj)*l1(kk))
         x = cos1*cos2 + sqrt(max(1d0-cos1*cos1,0d0)) * sqrt(max(1d0-cos2*cos2,0d0))
         sind2 = 1d0 - x*x
         lave2 = ((l2(jj)+l1(jj)) * (l2(kk)+l1(kk))) / 4d0
         if (sind2 > 1d-12) then
            if (sind2*lave2 > eps*eps) return
         end if
      end do
      ok = .true.

    end function metric_ok
  end subroutine lattpg

end submodule guesspg
