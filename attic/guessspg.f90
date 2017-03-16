  ! private to guessspg
  integer, parameter :: maxch=1000 !< max. possible unit cell translations
  real*8, allocatable  :: disctr(:,:) !< possible unit cell translations
  real*8, parameter :: pusheps = 1d-14 !< Reduce atom coords. to (-pusheps,1-pusheps]
  integer :: nch !< number of checked translations

  ! private for guessspg crystal symmetry guessing module
  private :: cenbuild
  private :: cenclosure
  private :: cenreduce
  private :: filltrans
  private :: goodop
  private :: reduce
  private :: iscelltr
  private :: isrepeated

  !xx! guessspg - crystal symmetry guessing module
  !> Guesses the symmetry operations from the geometry of the unit
  !> cell and the positions of the atoms in it. Transform the atom
  !> list into the non-equivalent atom list.  In: cell parameters (aa,
  !> bb) Inout: nneq, at(:)%z, at(:)%x, at(:)%name Out: lcent, neqv,
  !> rotm. If level = 0, use no symmetry. If level = 1, find only
  !> the centering vectors. Level = 2, full symmetry.
  subroutine guessspg(c,level)
    use tools_math, only: crys2car_from_cellpar, crys2car_from_cellpar
    use param, only: eyet

    class(crystal), intent(inout) :: c
    integer, intent(in) :: level

    real*8 :: rmat(3,3)

    ! no symmetry
    if (level == 0 .and. c%havesym == 0) then
       c%neqv = 1
       c%rotm(:,:,1) = eyet
       c%ncv = 1
       c%cen = 0d0
       c%lcent = 0
       c%havesym = 0
       return
    end if

    ! check if we already have this level
    if (c%havesym >= level) return

    ! write down the new level
    c%havesym = level

    ! Find the centering group by transforming to a primitive cell
    call cenbuild(c)

    ! Find the non-centering operations?
    if (level > 1) then
       ! Find the rotation parts of the operations
       rmat = transpose(crys2car_from_cellpar(c%aa,c%bb))
       call lattpg(rmat,c%ncv,c%cen,c%neqv,c%rotm(1:3,1:3,:))

       ! Calculate the translation vectors
       call filltrans(c)
    end if

  end subroutine guessspg

  !> Find all centering vectors in the crystal. Uses c%nneq and c%at,
  !> and writes c%cen and c%ncv.
  subroutine cenbuild (c)
    use global, only: atomeps
    use types, only: realloc
    type(crystal), intent(inout) :: c

    real*8  :: tr(3)
    integer :: i, j, k
    logical :: checked

    ! add the trivial centering vector
    if (allocated(c%cen)) deallocate(c%cen)
    allocate(c%cen(3,4))
    c%ncv = 1
    c%cen = 0d0

    if (allocated(disctr)) deallocate(disctr)
    allocate(disctr(3,maxch))
    disctr = 0d0
    ! check all possible translations
    nch = 0
    do i = 1, c%nneq
       do j = i+1, c%nneq
          if (c%at(i)%z .eq. c%at(j)%z) then
             tr = c%at(i)%x - c%at(j)%x
             call reduce(tr)

             checked = .false.
             k = 0
             do while (.not.checked .and. k .lt. nch)
                k = k + 1
                if (c%are_close(disctr(:,k),tr,atomeps)) then
                   checked = .true.
                end if
             end do

             ! if it is a true translation
             if (.not. checked) then
                nch = nch + 1
                if (nch > size(disctr,2)) &
                   call realloc(disctr,3,2*size(disctr,2))
                disctr(:,nch) = tr
                !
                if (.not.isrepeated(c,tr)) then
                   if (iscelltr(c,tr)) then
                      c%ncv = c%ncv + 1
                      if (c%ncv > size(c%cen,2)) call realloc(c%cen,3,2*c%ncv)
                      c%cen(:,c%ncv) = tr
                      call cenclosure(c)
                   end if
                end if
             endif

             ! also its inverse may be a translation
             tr(1) = -tr(1)
             tr(2) = -tr(2)
             tr(3) = -tr(3)
             call reduce(tr)

             if (.not. checked) then
                nch = nch + 1
                if (nch > size(disctr,2)) &
                   call realloc(disctr,3,2*size(disctr,2))
                disctr(:,nch) = tr

                if (.not.isrepeated(c,tr)) then
                   if (iscelltr(c,tr)) then
                      c%ncv = c%ncv + 1
                      if (c%ncv > size(c%cen,2)) call realloc(c%cen,3,2*c%ncv)
                      c%cen(:,c%ncv) = tr
                      call cenclosure(c)
                   end if
                endif
             end if
          end if
       end do
    end do
    deallocate(disctr)
    call realloc(c%cen,3,c%ncv)

  end subroutine cenbuild

  !> Determine the full cetering group using the
  !> (possibly repeated) centering vectors.
  subroutine cenclosure(c)
    use types, only: realloc
    type(crystal), intent(inout) :: c

    integer :: fnc
    integer :: i, j, k
    logical :: doagain

    if (c%ncv == 0) return

    ! purge the equivalent vectors
    call cenreduce(c,c%ncv,c%cen)

    ! generate all possible combinations until consistency
    fnc = -1
    doagain = .true.
    do while(doagain)
       fnc = c%ncv
       ! abelian group -> only over pairs (perhaps same-pairs)
       do i = 1, c%ncv
          do j = i, c%ncv
             fnc = fnc + 1
             if (fnc > size(c%cen,2)) call realloc(c%cen,3,2*fnc)
             c%cen(:,fnc) = c%cen(:,i) + c%cen(:,j)
             call reduce(c%cen(:,fnc))

             do k = 1, nch
                if (abs(disctr(1,k) - c%cen(1,fnc)) .lt. 1d-5 .and. &
                   abs(disctr(2,k) - c%cen(2,fnc)) .lt. 1d-5 .and. &
                   abs(disctr(3,k) - c%cen(3,fnc)) .lt. 1d-5) then
                   fnc = fnc - 1
                   ! cycle the j-loop
                   goto 1
                end if
             end do

             nch = nch + 1
             if (nch > size(disctr,2)) &
                call realloc(disctr,3,2*size(disctr,2))
             disctr(:,nch) = c%cen(:,fnc)

1            continue
          end do
       end do
       call cenreduce(c,fnc,c%cen)
       if (fnc .gt. c%ncv) then
          c%ncv = fnc
          doagain = .true.
       else
          doagain = .false.
       end if
    end do

  end subroutine cenclosure

  !> Purges the repeated vectors in a list of centering
  !> vectors
  subroutine cenreduce(c,nc,cv)
    use global, only: atomeps
    type(crystal), intent(inout) :: c
    integer, intent(inout) :: nc !< Number of vectors in the list
    real*8, intent(inout) :: cv(3,nc) !< Array of vectors

    integer :: i, j, nnc
    logical :: found
    real*8  :: v1(3), v2(3)

    nnc = 0
    do i = 1, nc
       found = .false.
       j = 0
       do while(.not.found .and. j.lt.nnc)
          j = j + 1
          v1 = cv(:,i)
          v2 = cv(:,j)
          if (c%are_lclose(v1,v2,atomeps)) then
             found = .true.
          end if
       end do
       if (.not.found) then
          nnc = nnc + 1
          cv(:,nnc) = cv(:,i)
          call reduce(cv(:,nnc))
       end if
    end do
    nc = nnc

  end subroutine cenreduce

  !> Given a vector in crystallographic coordinates, reduce
  !> it to the main cell (xi \in (-pusheps,1-pusheps]).
  subroutine reduce(tr)

    real*8, intent(inout) :: tr(3) !< Inpout/output vector to reduce

    integer :: i

    do i = 1, 3
       if (tr(i) .gt. -pusheps) then
          tr(i) = tr(i) - int(tr(i)+pusheps)
       else
          tr(i) = tr(i) - int(tr(i)+pusheps) + 1d0
       end if
    end do

  end subroutine reduce

  !> Calculates if a given translation vector (tr)
  !> is a symmetry operation by itself in the current cell setup.
  !>
  !> This routine is part of the spg operations guessing algorithm by
  !> teVelde, described in his PhD thesis.
  function iscelltr(c,tr)
    use global, only: atomeps
    logical :: iscelltr
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: tr(3) !< Cell translation vector to check

    integer :: i, j
    real*8  :: v1(3), v2(3)

    iscelltr = .true.
    do i = 1, c%nneq
       iscelltr = .false.
       j = 0
       do while (.not.iscelltr .and. j < c%nneq)
          j = j + 1
          if (i /= j) then
             if (c%at(i)%z == c%at(j)%z) then
                v1 = c%at(i)%x + tr
                v2 = c%at(j)%x
                if (c%are_lclose(v1,v2,atomeps)) then
                   iscelltr = .true.
                end if
             end if
          end if
       end do
       if (.not.iscelltr) return
    end do

  end function iscelltr

  !> For a given trasnlation vector, determines if it
  !> is contained in the ncv / cen
  function isrepeated(c,tr)
    use global, only: atomeps
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: tr(3) !< Vector to check
    logical :: isrepeated

    integer :: i

    isrepeated = .false.
    do i = 1, c%ncv
       if (c%are_lclose(c%cen(:,i),tr,atomeps)) then
          isrepeated = .true.
          return
       end if
    end do

  end function isrepeated

  !> Given {abc}red, {xyz}neq, ncv, cen, neqv and the
  !> rotational parts of the space group operators (rotm(1:3,1:3,i)),
  !> fill their translational part (rotm(:,4,i)).
  !>
  !> This routine is part of the spg operations guessing algorithm by
  !> teVelde, described in his PhD thesis.
  subroutine filltrans(c)
    use global, only: atomeps
    type(crystal), intent(inout) :: c

    integer :: op, i, j, k, n
    real*8  :: xnew(3)
    logical :: saveop(48), doi, doj, dok

    ! skip identity
    saveop(1) = .true.
    do op = 2, c%neqv
       saveop(op) = .false.
       doi = .true.
       i = 0
       do while (doi .and. i .lt. c%nneq)
          i = i + 1
          ! the transformed i atom
          xnew = matmul(c%rotm(1:3,1:3,op),c%at(i)%x)
          doj = .true.
          j = 0
          do while (doj .and. j .lt. c%nneq)
             j = j + 1
             if (c%at(i)%z .eq. c%at(j)%z) then
                !.propose a translation vector
                c%rotm(:,4,op) = c%at(j)%x - xnew - floor(c%at(j)%x - xnew)
                !.if it is truly an operation exit the i- and j- loops
                if (goodop(c,op)) then
                   saveop(op) = .true.
                   doj = .false.
                end if
             end if
          end do
          if (saveop(op)) doi = .false.
       end do
    end do

    ! rewrite the list of operations
    n = 0
    do op = 1, c%neqv
       if (saveop(op)) then
          n = n + 1
          do k = 1, 4
             do j = 1, 3
                c%rotm(j,k,n) = c%rotm(j,k,op)
             end do
          end do
          call reduce(c%rotm(:,4,n))
          dok = .true.
          k = 1
          do while (dok .and. k<c%ncv)
             k = k + 1
             if (c%are_lclose(c%rotm(:,4,n),c%cen(:,k),atomeps)) then
                c%rotm(1,4,n) = 0d0
                c%rotm(2,4,n) = 0d0
                c%rotm(3,4,n) = 0d0
                dok = .false.
             end if
          end do
       end if
    end do
    c%neqv = n

  end subroutine filltrans

  !> Check if a symmetry operation (rotm) is consistent with
  !> the atomic positions.
  !>
  !> This routine is part of the spg operations guessing algorithm by
  !> teVelde, described in his PhD thesis.
  function goodop(c,op)
    use global, only: atomeps
    logical :: goodop
    type(crystal), intent(inout) :: c
    integer, intent(in) :: op !< Operation identifier

    integer :: i, j
    real*8  :: xnew(3), v1(3), v2(3)
    logical :: found

    goodop = .true.
    do i = 1, c%nneq
       xnew = matmul(c%rotm(1:3,1:3,op),c%at(i)%x) + c%rotm(:,4,op)
       found = .false.
       do j = 1, c%nneq
          if (c%at(i)%z .eq. c%at(j)%z) then
             v1 = xnew
             v2 = c%at(j)%x

             found = (c%are_lclose(v1,v2,atomeps))
             if (found) exit
          end if
       end do
       if (.not.found) then
          goodop = .false.
          return
       end if
    end do

  end function goodop

