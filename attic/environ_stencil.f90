  module function make_stencil(e,rcut) result(st)
    class(environ), intent(in) :: e
    real*8, intent(in) :: rcut
    type(stencil) :: st

    real*8 :: rsphmin
    integer :: i1, i2, i3, ix(3), ibase

    st%rcut = rcut
    st%rsph = rcut + e%boxsize / sqrt(2d0)

    rsphmin = 1.5d0 * e%boxsize

    if (st%rsph <= rsphmin) then
       st%rsph = rsphmin
       st%rcut = rsphmin - e%boxsize / sqrt(2d0)

       allocate(st%iadd(27))
       st%nreg = 1
       st%iadd(1) = 0
       ix = 0
       ibase = e%p2i(ix)

       do i1 = -1, 1
          ix(1) = i1
          do i2 = -1, 1
             ix(2) = i2
             do i3 = -1, 1
                ix(3) = i3
                if (all(ix == 0)) cycle
                st%nreg = st%nreg + 1
                st%iadd(st%nreg) = e%p2i(ix) - ibase
             end do
          end do
       end do
    else
       write (*,*) "bleh2 not implemented in make_stencil"
    end if

  end function make_stencil

