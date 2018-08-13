  !> Calculates the neighbor environment of a given point x0 (cryst.
  !> coords) up to shell shmax. Return the number of neighbors for
  !> each shell in nneig, the non-equivalent atom index in wat,
  !> and the distance in dist. If the argument xenv is present,
  !> return the position of a representative atom from each shell.
  module subroutine pointshell(c,x0,shmax,nneig,wat,dist,xenv)
    use global, only: atomeps
    use types, only: realloc
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    integer, intent(in) :: shmax
    integer, intent(out) :: nneig(shmax)
    integer, intent(out) :: wat(shmax)
    real*8, intent(out) :: dist(shmax)
    real*8, intent(inout), allocatable, optional :: xenv(:,:,:)

    integer :: j, l, n
    real*8 :: aux(shmax), d2, x0c(3)
    integer :: iaux(shmax)
    real*8, allocatable :: aux2(:,:,:)

    x0c = c%x2c(x0 - floor(x0))
    dist = 1d30
    if (present(xenv)) then
       if (allocated(xenv)) deallocate(xenv)
       allocate(xenv(3,5,shmax),aux2(3,5,shmax))
       xenv = 0d0
    endif
    nneig = 0
    wat = 0
    do j = 1, c%env%n
       d2 = norm2(c%env%at(j)%r-x0c)
       if (d2 < atomeps) cycle
       do l = 1, shmax
          if (abs(d2 - dist(l)) < atomeps) then
             nneig(l) = nneig(l) + 1
             if (present(xenv)) then
                if (nneig(l) > size(xenv,2)) then
                   n = size(xenv)
                   call realloc(xenv,3,2*n,shmax)
                   call realloc(aux2,3,2*n,shmax)
                   xenv(:,n+1:,:) = 0d0
                endif
                xenv(:,nneig(l),l) = c%env%at(j)%x
             endif
             exit
          else if (d2 < dist(l)) then
             aux(l:shmax-1) = dist(l:shmax-1)
             dist(l+1:shmax) = aux(l:shmax-1)
             iaux(l:shmax-1) = nneig(l:shmax-1)
             nneig(l+1:shmax) = iaux(l:shmax-1)
             iaux(l:shmax-1) = wat(l:shmax-1)
             wat(l+1:shmax) = iaux(l:shmax-1)
             if (present(xenv)) then
                aux2(:,:,l:shmax-1) = xenv(:,:,l:shmax-1)
                xenv(:,:,l+1:shmax) = aux2(:,:,l:shmax-1)
             endif

             dist(l) = d2
             nneig(l) = 1
             wat(l) = c%env%at(j)%idx
             if (present(xenv)) then
                xenv(:,1,l) = c%env%at(j)%x
             endif
             exit
          end if
       end do
    end do
    if (allocated(aux2)) deallocate(aux2)

  end subroutine pointshell

