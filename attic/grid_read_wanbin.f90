  !> Read a collection of Wannier functions in binary format
  subroutine grid_read_wanbin(file,f,omega)
    use tools_io, only: fopen_read, getline_raw, lgetword, equal, ferror, faterr, &
       fclose
    use types, only: field, realloc

    character*(*), intent(in) :: file !< Input file
    type(field), intent(inout) :: f
    real*8, intent(in) :: omega

    integer :: luc
    integer :: n(3), nwan(3), ns(3), nbnd, nspin
    integer :: ispin, ibnd
    integer :: i, j, k
    real*8 :: fspin

    ! open file for reading
    luc = fopen_read(file,form="unformatted")

    ! number of grid points and wannier functions
    read (luc) n, nwan, nbnd, nspin

    if (nspin == 1) then
       fspin = 2d0
    else
       fspin = 1d0
       call ferror('grid_read_wanbin','nspin = 2 not implemented yet',faterr)
    end if

    ! allocate space for the wannier functions in memory
    f%n = n
    f%nwan = nwan
    ns = n * nwan
    if (allocated(f%f)) deallocate(f%f)
    allocate(f%f(f%n(1),f%n(2),f%n(3)))

    f%f = 0d0
    do ispin = 1, nspin
       do ibnd = 1, nbnd
          read (luc) f%fwan(:,:,:,ibnd,ispin)
          do i = 1, nwan(1)
             do j = 1, nwan(2)
                do k = 1, nwan(3)
                   f%f = f%f + real(f%fwan((i-1)*n(1)+1:i*n(1),(j-1)*n(2)+1:j*n(2),(k-1)*n(3)+1:k*n(3),ibnd,ispin),8)**2 + &
                      aimag(f%fwan((i-1)*n(1)+1:i*n(1),(j-1)*n(2)+1:j*n(2),(k-1)*n(3)+1:k*n(3),ibnd,ispin))**2
                end do
             end do
          end do
       end do
    end do
    f%f = f%f * fspin / omega

    ! wrap up
    call fclose(luc)
    f%init = .true.
    f%mode = mode_default

  end subroutine grid_read_wanbin

