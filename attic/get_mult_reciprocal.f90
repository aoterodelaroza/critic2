  !> Calculate the multiplicity of the point x0 in reciprocal space
  !> (fractional coordinates). 
  module function get_mult_reciprocal(c,x0) result (mult)
    class(crystal), intent(in) :: c
    real*8, intent(in) :: x0(3)
    integer :: mult

    real*8 :: xp(3), x(3), d(3)
    real*8 :: avec(3,48)
    integer :: i, j, nvec
    logical :: found

    mult = 0

    !.Translate position to main cell
    xp = x0 - floor(x0)

    !.Run over symmetry operations and create the list of copies
    nvec = 0
    do i = 1, c%neqvg
       x = matmul(c%rotg(:,:,i),xp)
       found = .false.
       do j = 1, nvec
          d = x - avec(:,j)
          d = abs(d - nint(d))
          if (all(d < 1d-5)) then
             found = .true.
             exit
          end if
       end do
       if (.not.found) then
          nvec = nvec + 1
          avec(:,nvec) = x
       end if
    enddo
    mult = nvec

  end function get_mult_reciprocal


