! Copyright (c) 2015 Alberto Otero de la Roza <alberto@fluor.quimica.uniovi.es>,
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

! tools for grids
module grid_tools
  implicit none

  private

  public :: grid_from_array3
  public :: grid_read_cube
  public :: grid_read_siesta
  public :: grid_read_abinit
  public :: grid_read_vasp
  public :: grid_read_qub
  public :: grid_read_xsf
  public :: grid_read_elk
  public :: grinterp
  public :: init_trispline
  public :: grid_gradrho
  public :: grid_rhoat
  public :: grid_laplacian
  public :: grid_hxx

  integer, parameter, public :: mode_nearest = 1 !< interpolation mode: nearest grid node
  integer, parameter, public :: mode_trilinear = 2 !< interpolation mode: trilinear
  integer, parameter, public :: mode_trispline = 3 !< interpolation mode: trispline
  integer, parameter, public :: mode_default = mode_trispline

contains

  !> Build a grid field from a three-dimensional array
  subroutine grid_from_array3(g,f)
    use tools_io
    use types

    real*8, intent(in) :: g(:,:,:)
    type(field), intent(out) :: f

    integer :: istat, n(3)

    f%init = .true.
    f%mode = mode_default
    n = ubound(g) - lbound(g) + 1
    f%n(:) = n
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('grid_from_array3','Error allocating grid',faterr)
    f%f = g

  end subroutine grid_from_array3

  !> Read a grid in gaussian CUBE format
  subroutine grid_read_cube(file,f,verbose)
    use tools_io
    use types

    character*(*), intent(in) :: file !< Input file
    type(field), intent(out) :: f
    logical, intent(in) :: verbose

    integer :: luc
    integer :: nat
    integer :: istat, n(3), i, j, k

    if (verbose) then
       write (uout,'("* GRID input, from a CUBE file")')
       write (uout,'("  In file: ",A)') trim(file)
    end if

    luc = fopen_read(file)

    read (luc,*) 
    read (luc,*) 
    read (luc,*,iostat=istat) nat
    if (istat /= 0) &
       call ferror('grid_read_cube','Error reading nat',faterr,file)
    do i = 1, 3
       read (luc,*,iostat=istat) n(i)
       if (istat /= 0) &
          call ferror('grid_read_cube','Error reading nx, ny, nz',faterr,file)
    end do
    do i = 1, nat
       read (luc,*)
    end do

    f%init = .true.
    f%mode = mode_default
    f%n(:) = n
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('grid_read_cube','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),k=1,n(3)),j=1,n(2)),i=1,n(1))
    if (istat /= 0) &
       call ferror('grid_read_cube','Error reading grid',faterr,file)

    if (verbose) then
       write (uout,'("  Grid dimensions : ",3(A,2X))') (string(f%n(j)),j=1,3)
       write (uout,'("  First elements... ",3(A,2X))') (string(f%f(1,1,j),'e',decimal=12),j=1,3)
       write (uout,'("  Last elements... ",3(A,2X))') (string(f%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
       write (uout,'("  Sum of elements... ",A)') string(sum(f%f(:,:,:)),'e',decimal=12)
       write (uout,'("  Sum of squares of elements... ",A)') string(sum(f%f(:,:,:)**2),'e',decimal=12)
       write (uout,'("  GRID input successful ")') 
       write (uout,*)
    end if

    call fclose(luc)

  end subroutine grid_read_cube

  !> Read a grid in siesta RHO format
  subroutine grid_read_siesta(file,f,verbose)
    use tools_io
    use types

    character*(*), intent(in) :: file !< Input file
    type(field), intent(out) :: f
    logical, intent(in) :: verbose

    integer :: luc, nspin, istat, j
    integer :: i, iy, iz, n(3)
    real*8 :: r(3,3)
    real*4, allocatable :: g(:)

    if (verbose) then
       write (uout,'("* GRID input, from a SIESTA grid file")')
       write (uout,'("  In file: ",A)') trim(file)
    end if

    ! open file
    luc = fopen_read(file,'unformatted')

    ! initialize
    f%init = .true.
    f%mode = mode_default

    ! assume unformatted
    read (luc) r
    read (luc) n, nspin
    f%n = n

    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('grid_read_siesta','Error allocating grid',faterr,file)
    allocate(g(n(1)),stat=istat)
    if (istat /= 0) &
       call ferror('grid_read_siesta','Error allocating auxiliary grid',faterr,file)
    f%f = 0d0
    do i = 1, nspin
       do iz = 1, n(3)
          do iy = 1, n(2)
             read (luc) g
             f%f(:,iy,iz) = f%f(:,iy,iz) + g
          end do
       end do
    end do
    deallocate(g)

    if (verbose) then
       write (uout,'("  Grid dimensions : ",3(A,2X))') (string(f%n(j)),j=1,3)
       write (uout,'("  First elements... ",3(A,2X))') (string(f%f(1,1,j),'e',decimal=12),j=1,3)
       write (uout,'("  Last elements... ",3(A,2X))') (string(f%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
       write (uout,'("  Sum of elements... ",A)') string(sum(f%f(:,:,:)),'e',decimal=12)
       write (uout,'("  Sum of squares of elements... ",A)') string(sum(f%f(:,:,:)**2),'e',decimal=12)
       write (uout,'("  GRID input successful ")') 
       write (uout,*)
    end if

    call fclose(luc)

  end subroutine grid_read_siesta

  !> Read a grid in abinit format
  subroutine grid_read_abinit(file,f,verbose)
    use types
    use tools_io
    use abinit_private

    character*(*), intent(in) :: file !< Input file
    type(field), intent(out) :: f
    logical, intent(in) :: verbose

    integer :: luc
    integer :: fform0, istat, n(3), j
    type(hdr_type) :: hdr
    real*8, allocatable :: g(:,:,:)

    if (verbose) then
       write (uout,'("* GRID input, from an abinit-style file")')
       write (uout,'("  In file: ",A)') trim(file)
    end if

    luc = fopen_read(file,'unformatted')

    ! read the header
    call hdr_io(fform0,hdr,1,luc)

    f%init = .true.
    f%mode = mode_default
    f%n(:) = hdr%ngfft(:)
    n = f%n
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('grid_read_abinit','Error allocating grid',faterr,file)
    allocate(g(n(1),n(2),n(3)))
    read(luc,iostat=istat) g
    f%f = g
    deallocate(g)
    if (istat /= 0) &
       call ferror('grid_read_abinit','Error reading grid',faterr,file)

    if (verbose) then
       write (uout,'("  Grid dimensions : ",3(A,2X))') (string(f%n(j)),j=1,3)
       write (uout,'("  First elements... ",3(A,2X))') (string(f%f(1,1,j),'e',decimal=12),j=1,3)
       write (uout,'("  Last elements... ",3(A,2X))') (string(f%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
       write (uout,'("  Sum of elements... ",A)') string(sum(f%f(:,:,:)),'e',decimal=12)
       write (uout,'("  Sum of squares of elements... ",A)') string(sum(f%f(:,:,:)**2),'e',decimal=12)
       write (uout,'("  GRID input successful ")') 
       write (uout,*)
    end if

    call fclose(luc)

  end subroutine grid_read_abinit

  !> Read a grid in VASP format
  subroutine grid_read_vasp(file,f,omega,verbose)
    use types
    use tools_io

    character*(*), intent(in) :: file !< Input file
    type(field), intent(out) :: f
    real*8, intent(in) :: omega
    logical, intent(in) :: verbose

    integer :: luc
    integer :: istat, n(3), i, j, k
    character(len=:), allocatable :: line
    logical :: ok

    if (verbose) then
       write (uout,'("* GRID input, from a VASP-style file")')
       write (uout,'("  In file: ",A)') trim(file)
    end if

    luc = fopen_read(file)

    do while(.true.)
       ok = getline_raw(luc,line,.true.)
       if (len(trim(line)) == 0) exit
    end do

    read (luc,*,iostat=istat) n
    if (istat /= 0) &
       call ferror('grid_read_vasp','Error reading nx, ny, nz',faterr,file)

    f%init = .true.
    f%mode = mode_default
    f%n(:) = n
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('grid_read_vasp','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('grid_read_vasp','Error reading grid',faterr,file)
    f%f(:,:,:) = f%f(:,:,:) / omega

    if (verbose) then
       write (uout,'("  Grid dimensions : ",3(A,2X))') (string(f%n(j)),j=1,3)
       write (uout,'("  First elements... ",3(A,2X))') (string(f%f(1,1,j),'e',decimal=12),j=1,3)
       write (uout,'("  Last elements... ",3(A,2X))') (string(f%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
       write (uout,'("  Sum of elements... ",A)') string(sum(f%f(:,:,:)),'e',decimal=12)
       write (uout,'("  Sum of squares of elements... ",A)') string(sum(f%f(:,:,:)**2),'e',decimal=12)
       write (uout,'("  GRID input successful ")') 
       write (uout,*)
    end if

    call fclose(luc)

  end subroutine grid_read_vasp

  !> Read a grid in aimpac qub format
  subroutine grid_read_qub(file,f,verbose)
    use tools_io
    use types

    character*(*), intent(in) :: file !< Input file
    type(field), intent(out) :: f
    logical, intent(in) :: verbose

    integer :: luc
    integer :: istat, n(3), i, j, k

    if (verbose) then
       write (uout,'("* GRID input, from a AIMPAC QUB file")')
       write (uout,'("  In file: ",A)') trim(file)
    end if

    luc = fopen_read(file)

    read (luc,*,iostat=istat) f%n
    if (istat /= 0) &
       call ferror('grid_read_qub','Error reading n1, n2, n3',faterr,file)

    f%init = .true.
    f%mode = mode_default
    n = f%n(:)
    allocate(f%f(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('grid_read_qub','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((f%f(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('grid_read_qub','Error reading grid',faterr,file)

    if (verbose) then
       write (uout,'("  Grid dimensions : ",3(A,2X))') (string(f%n(j)),j=1,3)
       write (uout,'("  First elements... ",3(A,2X))') (string(f%f(1,1,j),'e',decimal=12),j=1,3)
       write (uout,'("  Last elements... ",3(A,2X))') (string(f%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
       write (uout,'("  Sum of elements... ",A)') string(sum(f%f(:,:,:)),'e',decimal=12)
       write (uout,'("  Sum of squares of elements... ",A)') string(sum(f%f(:,:,:)**2),'e',decimal=12)
       write (uout,'("  GRID input successful ")') 
       write (uout,*)
    end if

    call fclose(luc)

  end subroutine grid_read_qub

  !> Read a grid in xcrysden xsf format -- only first 3d grid in first 3d block
  subroutine grid_read_xsf(file,f,verbose,nwan,nin,omega)
    use tools_io
    use types

    character*(*), intent(in) :: file !< Input file
    type(field), intent(inout) :: f
    logical, intent(in) :: verbose
    integer, intent(in), optional :: nwan
    integer, intent(in), optional :: nin(3)
    real*8, intent(in), optional :: omega

    integer :: luc
    integer :: istat, n(3), lp, i, j, k, ix, jx, kx
    character(len=:), allocatable :: line, word
    logical :: found, ok, ok2, iswan
    real*8, dimension(3) :: x0, x1, x2, x3
    real*8 :: pmat(3,3)

    real*8, allocatable :: ggloc(:,:,:), otherloc(:,:,:)

    iswan = (present(nwan).and.present(nin).and.present(omega))

    ! header
    if (verbose) then
       if (.not.iswan) then
          write (uout,'("* GRID input, from a xcrysden XSF file")')
          write (uout,'("  In file: ",A)') trim(file)
       else if (nwan == 1) then
          write (uout,'("* Multiple GRID input, from xcrysden XSF file generated by wannier90")')
       end if
    end if

    ! open file for reading
    luc = fopen_read(file)

    ! position at the beginning of the first grid, ignore the rest
    found = .false.
    ok = .false.
    do while (getline_raw(luc,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,'primvec'))then
          ok = .true.
          read(luc,*,iostat=istat) pmat
          if (istat /= 0) &
             call ferror('grid_read_xsf','Error PRIMVEC',faterr,file)
       else if (equal(word,'begin_block_datagrid_3d'))then
          found = .true.
          exit
       end if
    end do
    if (.not.found) call ferror('grid_read_xsf','BEGIN_BLOCK_DATAGRID_3D not found',faterr,file)
    if (.not.ok) call ferror('grid_read_xsf','PRIMVEC not found',faterr,file)

    found = .false.
    do while (getline_raw(luc,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word(1:min(17,len(word))),'begin_datagrid_3d') .or. equal(word(1:min(11,len(word))),'datagrid_3d')) then
          found = .true.
          exit
       end if
    end do
    if (.not.found) call ferror('grid_read_xsf','BEGIN_DATAGRID_3D... not found',faterr,file)

    ! grid dimension
    read (luc,*,iostat=istat) n 
    if (istat /= 0) &
       call ferror('grid_read_xsf','Error reading n1, n2, n3',faterr,file)
    if (iswan) then
       if (any(nint(real(n,8)/nin) /= (n/nin))) &
          call ferror('grid_read_xsf','Wannier supercell not congruent with grid size',faterr,file)
       if (nwan > 1) then
          if (any(f%n * nin /= n)) &
             call ferror('grid_read_xsf','Wannier grid sizes not equal',faterr,file)
       else
          f%n = n / nin
       end if
    else
       f%n = n - 1
       iswan = .false.
    endif

    ! origin and edge vectors
    read (luc,*,iostat=istat) x0, x1, x2, x3

    f%init = .true.
    f%mode = mode_default
    allocate(ggloc(n(1),n(2),n(3)),stat=istat)
    if (istat /= 0) &
       call ferror('grid_read_xsf','Error allocating grid',faterr,file)
    if (istat /= 0) &
       call ferror('grid_read_xsf','Error allocating grid',faterr,file)
    read(luc,*,iostat=istat) (((ggloc(i,j,k),i=1,n(1)),j=1,n(2)),k=1,n(3))
    if (istat /= 0) &
       call ferror('grid_read_xsf','Error reading grid',faterr,file)

    allocate(f%f(f%n(1),f%n(2),f%n(3)),stat=istat)
    if (.not.iswan) then
       f%f = ggloc(1:n(1)-1,1:n(2)-1,1:n(3)-1)
       deallocate(ggloc)
    else
       ! re-order the grid from wannier90
       allocate(otherloc(n(1),n(2),n(3)))
       do i = 1, n(1)
          ix = mod(i,n(1))+1
          do j = 1, n(2)
             jx = mod(j,n(2))+1
             do k = 1, n(3)
                kx = mod(k,n(3))+1
                otherloc(i,j,k) = ggloc(ix,jx,kx)
             end do
          end do
       end do
       deallocate(ggloc)
       call move_alloc(otherloc,ggloc)

       ! save the wannier function
       if (.not.allocated(f%fwan)) allocate(f%fwan(n(1),n(2),n(3),1))
       if (nwan > size(f%fwan,4)) call realloc(f%fwan,n(1),n(2),n(3),nwan)
       f%fwan(:,:,:,nwan) = ggloc

       ! convert to density contribution
       ggloc = ggloc**2 * 2d0 / omega

       ! add to the density
       if (nwan == 1) then
          f%nwan = nin
          f%f = 0d0
       end if
       do i = 1, nin(1)
          do j = 1, nin(2)
             do k = 1, nin(3)
                f%f = f%f + ggloc((i-1)*f%n(1)+1:i*f%n(1),(j-1)*f%n(2)+1:j*f%n(2),(k-1)*f%n(3)+1:k*f%n(3))
             end do
          end do
       end do
       deallocate(ggloc)
    endif
    n = f%n

    if (verbose) then
       if (.not.iswan) then
          write (uout,'("  Grid dimensions : ",3(A,2X))') (string(f%n(j)),j=1,3)
          write (uout,'("  First elements... ",3(A,2X))') (string(f%f(1,1,j),'e',decimal=12),j=1,3)
          write (uout,'("  Last elements... ",3(A,2X))') (string(f%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
          write (uout,'("  Sum of elements... ",A)') string(sum(f%f(:,:,:)),'e',decimal=12)
          write (uout,'("  Sum of squares of elements... ",A)') string(sum(f%f(:,:,:)**2),'e',decimal=12)
          write (uout,'("  GRID input successful ")') 
          write (uout,*)
       else
          write (uout,'("  Wannier supercell-grid read from file ",A)') file
          write (uout,'("  Band index ",A)') string(nwan)
          write (uout,'("  Grid dimensions... ",3(A,X))') (string(f%nwan(j)),j=1,3)
          write (uout,'("  Accum. First elements... ",3(A,2X))') (string(f%f(1,1,j),'e',decimal=12),j=1,3)
          write (uout,'("  Accum. Last elements... ",3(A,2X))') (string(f%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
          write (uout,'("  Accum. Sum of elements... ",A)') string(sum(f%f(:,:,:)),'e',decimal=12)
       end if
    end if

    call fclose(luc)

  end subroutine grid_read_xsf

  !> Read a grid in elk format -- only first 3d grid in first 3d block
  subroutine grid_read_elk(file,f,verbose)
    use tools_io
    use types

    character*(*), intent(in) :: file !< Input file
    type(field), intent(out) :: f
    logical, intent(in) :: verbose

    integer :: luc, ios
    integer :: n(3), i, j, k
    real*8 :: dum(3)

    ! header
    if (verbose) then
       write (uout,'("* GRID input, from a ELK file")')
       write (uout,'("  In file: ",A)') trim(file)
    end if

    ! open file for reading
    luc = fopen_read(file)

    ! grid dimension
    read (luc,*,iostat=ios) n
    if (ios /= 0) &
       call ferror('grid_read_elk','Error reading n1, n2, n3',faterr,file)

    f%n = n - 1
    allocate(f%f(n(1)-1,n(2)-1,n(3)-1),stat=ios)
    if (ios /= 0) &
       call ferror('grid_read_elk','Error allocating grid',faterr,file)
    do k = 1, n(3)-1
       do j = 1, n(2)-1
          do i = 1, n(1)-1
             read (luc,*) dum, f%f(i,j,k)
          end do
       end do
    end do

    n = f%n
    f%init = .true.
    f%mode = mode_default

    if (verbose) then
       write (uout,'("  Grid dimensions : ",3(A,2X))') (string(f%n(j)),j=1,3)
       write (uout,'("  First elements... ",3(A,2X))') (string(f%f(1,1,j),'e',decimal=12),j=1,3)
       write (uout,'("  Last elements... ",3(A,2X))') (string(f%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
       write (uout,'("  Sum of elements... ",A)') string(sum(f%f(:,:,:)),'e',decimal=12)
       write (uout,'("  Sum of squares of elements... ",A)') string(sum(f%f(:,:,:)**2),'e',decimal=12)
       write (uout,'("  GRID input successful ")') 
       write (uout,*)
    end if

    call fclose(luc)

  end subroutine grid_read_elk

  !> Interpolate the function value, first and second derivative at
  !> point x0 (crystallographic coords.) using the grid g.
  subroutine grinterp(f,x0,y,yp,ypp) 
    use types

    type(field), intent(inout) :: f !< Grid to interpolate
    real*8, intent(in) :: x0(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    y = 0d0
    yp = 0d0
    ypp = 0d0

    if (f%mode == mode_nearest) then
       call grinterp_nearest(f,x0,y)
       yp = 0d0
       ypp = 0d0
    else if (f%mode == mode_trilinear) then
       call grinterp_trilinear(f,x0,y,yp)
       ypp = 0d0
    else  
       call grinterp_trispline(f,x0,y,yp,ypp)
    end if

    !    ! transform to derivatives wrt cryst coordinates
    !    yp = matmul(transpose(f%x2c),yp)
    !    ypp = matmul(matmul(transpose(f%x2c),ypp),f%x2c)

  end subroutine grinterp

  !> Interpolate from the nearest point (or the would-be nearest, it is
  !> nearest only in orthogonal grids). 
  subroutine grinterp_nearest(f,x0,y)
    use types

    type(field), intent(in) :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value

    real*8 :: x(3)
    integer :: idx(3)

    x = modulo(x0,1d0)
    idx = grid_near(f,x)
    y = f%f(idx(1),idx(2),idx(3))

  end subroutine grinterp_nearest

  !> Interpolate using a trilinear interpolation.
  subroutine grinterp_trilinear(f,x0,y,yp)
    use types

    type(field), intent(in) :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point (cryst. coords.)
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative

    integer :: idx(3), iidx(3), i, j, k
    real*8 :: ff(0:2,0:2,0:2), r(3), s(3), x(3)

    ! compute value at the cube vertices
    x = modulo(x0,1d0)
    idx = grid_floor(f,x)
    do i = 0, 1
       do j = 0, 1
          do k = 0, 1
             iidx = modulo(idx+(/i,j,k/)-1,f%n)+1
             ff(i,j,k) = f%f(iidx(1),iidx(2),iidx(3))
          end do
       end do
    end do

    ! x and 1-x
    r = f%n * x - idx + 1
    s = 1d0 - r

    ! trilinear interpolation
    do i = 0, 1
       do j = 0, 1
          ff(i,j,2) = ff(i,j,0) * s(3) + ff(i,j,1) * r(3)
          ff(2,i,j) = ff(0,i,j) * s(1) + ff(1,i,j) * r(1)
       end do
       ff(i,2,2) = ff(i,0,2) * s(2) + ff(i,1,2) * r(2)
       ff(2,i,2) = ff(0,i,2) * s(1) + ff(1,i,2) * r(1)
       ff(2,2,i) = ff(2,0,i) * s(2) + ff(2,1,i) * r(2)
    end do
    ff(2,2,2) = ff(0,2,2) * s(1) + ff(1,2,2) * r(1)
    y = ff(2,2,2)
    yp(1) = ff(1,2,2) - ff(0,2,2)
    yp(2) = ff(2,1,2) - ff(2,0,2)
    yp(3) = ff(2,2,1) - ff(2,2,0)
    ! yp = matmul((f%n-1)*yp,f%c2x)

  end subroutine grinterp_trilinear

  !> Abinit trispline interpolation. Using the global cubic spline information (c2
  !> array) with periodic boundary conditions, calculate the density and c2 coeffs.
  !> of the 6 points on the cube faces, and then average the spline prediction. This 
  !> is a modified version of the abinit density interpolation subroutine.
  subroutine grinterp_trispline(f,x0,y,yp,ypp)
    use types
    use tools_io

    type(field), intent(inout), target :: f !< Input grid
    real*8, intent(in) :: x0(3) !< Target point
    real*8, intent(out) :: y !< Interpolated value
    real*8, intent(out) :: yp(3) !< First derivative
    real*8, intent(out) :: ypp(3,3) !< Second derivative

    integer :: i, ii, jj, kk, ll, nn, indx(3), oii, onn, omm
    integer :: inii(4,3)
    real*8 :: bbb, ddu(2), hrh(2), hh(4,2), grd(4), lder(4)
    real*8 :: cof(2,3), ddstar(6), rhstar(6), sqder(6,4), sqvlr(6,4)
    real*8 :: pomsq(2,3), pom2sq(2,3)
    real*8,pointer :: ptddx(:,:,:),ptddy(:,:,:),ptddz(:,:,:),ptrho(:,:,:)
    real*8 :: xx(3), dix(3)

    real*8, parameter :: eps = 1d-12

    !$omp critical (checkalloc)
    if (.not.allocated(f%c2)) then
       call init_trispline(f)
    end if
    !$omp end critical (checkalloc)

    nullify(ptddx,ptddy,ptddz,ptrho)

    xx = modulo(x0,1d0)

    do i=1,3
       dix(i)=1d0/f%n(i)
    end do

    ! determine the index in the grid
    do ii=1,3
       indx(ii)=int(xx(ii)*f%n(ii))
       bbb=(xx(ii)-indx(ii)*dix(ii))*f%n(ii)
       if (indx(ii)==f%n(ii)) then
          indx(ii)=1
          xx(ii)=0.d0
       else
          indx(ii)=indx(ii)+1
       end if
       !  Explicit handling to avoid numeric problems
       if (bbb > 1.d0+eps ) then
          cof(1,ii)=0.d0
          cof(2,ii)=1.d0
       elseif (bbb < -eps ) then
          cof(1,ii)=1.d0
          cof(2,ii)=0.d0
       else
          cof(1,ii)=1.d0-bbb
          cof(2,ii)=bbb
       end if
    end do

    ! 3d interpolation of the valence density

    ! determination of the values of density and of its second derivative
    ! at the "star" = constructed at vv with primitive directions
    ! To interpolation the values at the faces of the grid cell are needed

    rhstar(:)=0.d0
    sqder(:,:)=0.d0
    sqvlr(:,:)=0.d0
    ddstar(:)=0.d0
    pomsq(:,:)=0.d0
    pom2sq(:,:)=0.d0

    oii=1; onn=1; omm=1
    if (indx(1)==f%n(1)) oii=1-f%n(1)
    if (indx(2)==f%n(2)) onn=1-f%n(2)
    if (indx(3)==f%n(3)) omm=1-f%n(3)

    ! the values in the corners of the grid cell

    ptddx=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,1)
    ptddy=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,2)
    ptddz=>f%c2(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm,3)
    ptrho=>f%f(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)

    ! the coefficients for spline interpolation of density and its derivation
    do ii=1,3
       do jj=1,2
          pomsq(jj,ii)=(cof(jj,ii)*cof(jj,ii)*cof(jj,ii)-cof(jj,ii))/6.d0*dix(ii)*dix(ii)
          pom2sq(jj,ii)=(3.d0*cof(jj,ii)*cof(jj,ii)-1.d0)/6.d0*dix(ii)
          if (jj==1) pom2sq(jj,ii)=-pom2sq(jj,ii)
       end do
    end do


    do ii=1,2
       do jj=1,2
          do kk=1,2
             ddstar(ii)=ddstar(ii)+cof(jj,2)*cof(kk,3)*ptddx(ii,jj,kk)
             ddstar(ii+2)=ddstar(ii+2)+cof(jj,3)*cof(kk,1)*ptddy(kk,ii,jj)
             ddstar(ii+4)=ddstar(ii+4)+cof(jj,1)*cof(kk,2)*ptddz(jj,kk,ii)
             sqder(ii,jj)=sqder(ii,jj)+cof(kk,2)*ptddz(ii,kk,jj)
             sqder(ii,jj+2)=sqder(ii,jj+2)+cof(kk,3)*ptddy(ii,jj,kk)
             sqder(ii+2,jj)=sqder(ii+2,jj)+cof(kk,3)*ptddx(jj,ii,kk)
             sqder(ii+2,jj+2)=sqder(ii+2,jj+2)+cof(kk,1)*ptddz(kk,ii,jj)
             sqder(ii+4,jj)=sqder(ii+4,jj)+cof(kk,1)*ptddy(kk,jj,ii)
             sqder(ii+4,jj+2)=sqder(ii+4,jj+2)+cof(kk,2)*ptddx(jj,kk,ii)
             sqvlr(ii,jj)=sqvlr(ii,jj)+cof(kk,2)*ptrho(ii,kk,jj)+pomsq(kk,2)*ptddy(ii,kk,jj)
             sqvlr(ii,jj+2)=sqvlr(ii,jj+2)+cof(kk,3)*ptrho(ii,jj,kk)+pomsq(kk,3)*ptddz(ii,jj,kk)
             sqvlr(ii+2,jj+2)=sqvlr(ii+2,jj+2)+cof(kk,1)*ptrho(kk,ii,jj)+pomsq(kk,1)*ptddx(kk,ii,jj)
          end do
       end do
    end do

    do ii=1,2
       do jj=1,2
          sqvlr(ii+2,jj)=sqvlr(jj,ii+2)
          sqvlr(ii+4,jj)=sqvlr(jj+2,ii+2)
          sqvlr(ii+4,jj+2)=sqvlr(jj,ii)
       end do
    end do

    do ii=1,2
       do jj=1,2
          rhstar(ii)=rhstar(ii)+cof(jj,3)*sqvlr(ii,jj)+pomsq(jj,3)*sqder(ii,jj)+&
             &    cof(jj,2)*sqvlr(ii,jj+2)+pomsq(jj,2)*sqder(ii,jj+2)
          rhstar(ii+2)=rhstar(ii+2)+cof(jj,1)*sqvlr(ii+2,jj)+pomsq(jj,1)*sqder(ii+2,jj)+&
             &    cof(jj,3)*sqvlr(ii+2,jj+2)+pomsq(jj,3)*sqder(ii+2,jj+2)
          rhstar(ii+4)=rhstar(ii+4)+cof(jj,2)*sqvlr(ii+4,jj)+pomsq(jj,2)*sqder(ii+4,jj)+&
             &    cof(jj,1)*sqvlr(ii+4,jj+2)+pomsq(jj,1)*sqder(ii+4,jj+2)
       end do
    end do
    rhstar(:)=rhstar(:)/2.d0

    y = 0d0
    yp = 0d0
    ypp = 0d0
    kk=1; nn=1
    do ii=1,5,2
       do jj=1,2
          nn=-nn
          y=y+cof(jj,kk)*rhstar(ii+jj-1)+pomsq(jj,kk)*ddstar(ii+jj-1)
          yp(kk)=yp(kk)+pom2sq(jj,kk)*ddstar(ii+jj-1)
          ypp(kk,kk)=ypp(kk,kk)+cof(jj,kk)*ddstar(ii+jj-1)
          yp(kk)=yp(kk)+nn*rhstar(ii+jj-1)/dix(kk)
       end do
       kk=kk+1
    end do
    y=y/3.d0

    ! Off-diagonal elements of the hessian

    ! for the speed reasons the polynomial interpolation
    ! for second derivation fields is used in this case
    ! but the last step is always done by spline interpolation.


    do ii=1,3
       do jj=-1,2
          inii(jj+2,ii)=indx(ii)+jj
          if (inii(jj+2,ii) < 1) inii(jj+2,ii)=inii(jj+2,ii)+f%n(ii)
          if (inii(jj+2,ii) > f%n(ii)) inii(jj+2,ii)=inii(jj+2,ii)-f%n(ii)
       end do
    end do

    ! Not very nice

    do ii=1,3
       select case (ii)
       case (1)
          do jj=1,4
             ddu(1)=cof(1,2)*f%c2(inii(jj,1),inii(2,2),inii(2,3),3)+cof(2,2)*f%c2(inii(jj,1),inii(3,2),inii(2,3),3)
             ddu(2)=cof(1,2)*f%c2(inii(jj,1),inii(2,2),inii(3,3),3)+cof(2,2)*f%c2(inii(jj,1),inii(3,2),inii(3,3),3)
             hrh(1)=cof(1,2)*f%f(inii(jj,1),inii(2,2),inii(2,3))+cof(2,2)*f%f(inii(jj,1),inii(3,2),inii(2,3))+&
                &      pomsq(1,2)*f%c2(inii(jj,1),inii(2,2),inii(2,3),2)+pomsq(2,2)*f%c2(inii(jj,1),inii(3,2),inii(2,3),2)
             hrh(2)=cof(1,2)*f%f(inii(jj,1),inii(2,2),inii(3,3))+cof(2,2)*f%f(inii(jj,1),inii(3,2),inii(3,3))+&
                &      pomsq(1,2)*f%c2(inii(jj,1),inii(2,2),inii(3,3),2)+pomsq(2,2)*f%c2(inii(jj,1),inii(3,2),inii(3,3),2)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)

             ddu(1)=cof(1,3)*f%c2(inii(jj,1),inii(2,2),inii(2,3),2)+cof(2,3)*f%c2(inii(jj,1),inii(2,2),inii(3,3),2)
             ddu(2)=cof(1,3)*f%c2(inii(jj,1),inii(3,2),inii(2,3),2)+cof(2,3)*f%c2(inii(jj,1),inii(3,2),inii(3,3),2)
             hrh(1)=cof(1,3)*f%f(inii(jj,1),inii(2,2),inii(2,3))+cof(2,3)*f%f(inii(jj,1),inii(2,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(jj,1),inii(2,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(jj,1),inii(2,2),inii(3,3),3)
             hrh(2)=cof(1,3)*f%f(inii(jj,1),inii(3,2),inii(2,3))+cof(2,3)*f%f(inii(jj,1),inii(3,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(jj,1),inii(3,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(jj,1),inii(3,2),inii(3,3),3)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)
          end do
       case (2)
          do jj=1,4
             ddu(1)=cof(1,3)*f%c2(inii(2,1),inii(jj,2),inii(2,3),1)+cof(2,3)*f%c2(inii(2,1),inii(jj,2),inii(3,3),1)
             ddu(2)=cof(1,3)*f%c2(inii(3,1),inii(jj,2),inii(2,3),1)+cof(2,3)*f%c2(inii(3,1),inii(jj,2),inii(3,3),1)
             hrh(1)=cof(1,3)*f%f(inii(2,1),inii(jj,2),inii(2,3))+cof(2,3)*f%f(inii(2,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(2,1),inii(jj,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(2,1),inii(jj,2),inii(3,3),3)
             hrh(2)=cof(1,3)*f%f(inii(3,1),inii(jj,2),inii(2,3))+cof(2,3)*f%f(inii(3,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,3)*f%c2(inii(3,1),inii(jj,2),inii(2,3),3)+pomsq(2,3)*f%c2(inii(3,1),inii(jj,2),inii(3,3),3)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)

             ddu(1)=cof(1,1)*f%c2(inii(2,1),inii(jj,2),inii(2,3),3)+cof(2,1)*f%c2(inii(3,1),inii(jj,2),inii(2,3),3)
             ddu(2)=cof(1,1)*f%c2(inii(2,1),inii(jj,2),inii(3,3),3)+cof(2,1)*f%c2(inii(3,1),inii(jj,2),inii(3,3),3)
             hrh(1)=cof(1,1)*f%f(inii(2,1),inii(jj,2),inii(2,3))+cof(2,1)*f%f(inii(3,1),inii(jj,2),inii(2,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(jj,2),inii(2,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(jj,2),inii(2,3),1)
             hrh(2)=cof(1,1)*f%f(inii(2,1),inii(jj,2),inii(3,3))+cof(2,1)*f%f(inii(3,1),inii(jj,2),inii(3,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(jj,2),inii(3,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(jj,2),inii(3,3),1)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)
          end do
       case (3)
          do jj=1,4
             ddu(1)=cof(1,1)*f%c2(inii(2,1),inii(2,2),inii(jj,3),2)+cof(2,1)*f%c2(inii(3,1),inii(2,2),inii(jj,3),2)
             ddu(2)=cof(1,1)*f%c2(inii(2,1),inii(3,2),inii(jj,3),2)+cof(2,1)*f%c2(inii(3,1),inii(3,2),inii(jj,3),2)
             hrh(1)=cof(1,1)*f%f(inii(2,1),inii(2,2),inii(jj,3))+cof(2,1)*f%f(inii(3,1),inii(2,2),inii(jj,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(2,2),inii(jj,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(2,2),inii(jj,3),1)
             hrh(2)=cof(1,1)*f%f(inii(2,1),inii(3,2),inii(jj,3))+cof(2,1)*f%f(inii(3,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,1)*f%c2(inii(2,1),inii(3,2),inii(jj,3),1)+pomsq(2,1)*f%c2(inii(3,1),inii(3,2),inii(jj,3),1)
             hh(jj,2)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)

             ddu(1)=cof(1,2)*f%c2(inii(2,1),inii(2,2),inii(jj,3),1)+cof(2,2)*f%c2(inii(2,1),inii(3,2),inii(jj,3),1)
             ddu(2)=cof(1,2)*f%c2(inii(3,1),inii(2,2),inii(jj,3),1)+cof(2,2)*f%c2(inii(3,1),inii(3,2),inii(jj,3),1)
             hrh(1)=cof(1,2)*f%f(inii(2,1),inii(2,2),inii(jj,3))+cof(2,2)*f%f(inii(2,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,2)*f%c2(inii(2,1),inii(2,2),inii(jj,3),2)+pomsq(2,2)*f%c2(inii(2,1),inii(3,2),inii(jj,3),2)
             hrh(2)=cof(1,2)*f%f(inii(3,1),inii(2,2),inii(jj,3))+cof(2,2)*f%f(inii(3,1),inii(3,2),inii(jj,3))+&
                &      pomsq(1,2)*f%c2(inii(3,1),inii(2,2),inii(jj,3),2)+pomsq(2,2)*f%c2(inii(3,1),inii(3,2),inii(jj,3),2)
             hh(jj,1)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)
          end do
       end select
       do jj=-2,1
          grd(jj+3)=(indx(ii)+jj)*dix(ii)
       end do

       !  write(6,'("hh: ",/,4F16.8,/,4F16.8)') ((hh(kk,jj),kk=1,4),jj=1,2)
       !  write(6,'("grad: ",3F16.8)') (grad(kk),kk=1,3)
       !  write(6,'("dix: ",3F16.8)') (dix(kk),kk=1,3)
       !  write(6,'("grd: ",4F16.8)') (grd(kk),kk=1,4)
       !  write(6,'("inii: ",4I4)') (inii(kk,ii),kk=1,4)

       do jj=1,2

          !   polynomial interpolation

          do kk=1,3
             do ll=4,kk+1,-1
                hh(ll,jj)=(hh(ll,jj)-hh(ll-1,jj))/(grd(ll)-grd(ll-1))
             end do
          end do
          lder(4)=hh(4,jj)
          do kk=3,1,-1
             lder(kk)=hh(kk,jj)+(xx(ii)-grd(kk))*lder(kk+1)
          end do
          do kk=1,2
             do ll=3,kk+1,-1
                lder(ll)=lder(ll)+(xx(ii)-grd(ll-kk))*lder(ll+1)
             end do
          end do
          nn=ii+jj
          if (nn > 3) nn=nn-3
          ypp(ii,nn)=ypp(ii,nn)+lder(2)
          ypp(nn,ii)=ypp(nn,ii)+lder(2)
       end do
    end do

    ! averaging of the mixed derivations obtained in different order
    do ii=1,3
       do jj=1,3
          if (ii /= jj) ypp(ii,jj)=ypp(ii,jj)/2.d0
       end do
    end do

    nullify(ptddx,ptddy,ptddz,ptrho)

    ! back to cartesian
    ! atrans = transpose(f%c2x)
    ! yp = matmul(atrans,yp)
    ! ypp = matmul(matmul(atrans,ypp),transpose(atrans))

  end subroutine grinterp_trispline

  !> Pseudo-nearest grid point of a x (crystallographic) (only nearest in 
  !> orthogonal grids).
  function grid_near(f,x)
    use types

    type(field), intent(in) :: f !< Input grid
    real*8, intent(in) :: x(3) !< Target point (cryst. coords.)
    integer :: grid_near(3) 

    grid_near = modulo(nint(x * f%n),f%n)+1

  end function grid_near

  !> Floor grid point of a x (crystallographic)
  function grid_floor(f,x)
    use types

    type(field), intent(in) :: f !< Input grid
    real*8, intent(in) :: x(3) !< Target point (cryst. coords.)
    integer :: grid_floor(3)

    grid_floor = modulo(floor(x * f%n),f%n)+1

  end function grid_floor

  !> Initialize the grid for the abinit trispline interpolation.
  subroutine init_trispline(f)
    use tools_io
    use types

    type(field), intent(inout) :: f !< Input grid

    integer :: istat
    integer :: d, nmax, i, i1, i2
    real*8, allocatable :: l(:,:), fg(:)
    real*8 :: fprev, fnext, fuse, fone

    if (allocated(f%c2)) return

    allocate(f%c2(f%n(1),f%n(2),f%n(3),3),stat=istat)
    if (istat /= 0) &
       call ferror('init_trispline','Error allocating c2 grids',faterr)

    nmax = maxval(f%n)
    allocate(l(nmax,nmax))
    allocate(fg(nmax))

    ! cholesky decomposition of the matrix: 
    ! A = 
    !  ( 4 1 0 ... 0 1 )
    !  ( 1 4 1 ... 0 0 )
    !  ( 0 1 4 ... 0 0 )
    !      .......
    !  ( 0 0 ... 1 4 1 )
    !  ( 1 0 ... 0 1 4 )
    ! that is the coefficient matrix of the system that determines
    ! the x^2 coefficients of the cubic spline (array c2 / 2).
    ! L is a lower-triangular matrix, A = L * L^t

    ! direction x->y->z
    do d = 1, 3
       nmax = f%n(d)
       l = 0d0
       l(1,1) = 2d0
       l(2,1) = 1d0 / l(1,1)
       l(2,2) = sqrt(15d0) / 2d0
       l(3,2) = 1d0 / l(2,2)
       l(nmax,1) = 0.5d0
       l(nmax,2) = - 0.5d0 / sqrt(15d0)
       do i = 3, nmax-1
          l(i,i) = sqrt(4d0 - 1d0 / l(i-1,i-1)**2)
          l(i+1,i) = 1d0 / l(i,i)
          l(nmax,i) = - l(nmax,i-1) * l(i,i-1) / l(i,i)
       end do
       l(nmax,nmax-1) = (1d0 - l(nmax,nmax-2)*l(nmax-1,nmax-2)) /&
          l(nmax-1,nmax-1)
       l(nmax,nmax) = sqrt(4d0 - sum(l(nmax,1:nmax-1)**2))
       l = l / sqrt(6d0 * nmax**2)

       ! for each of the grid points in the plane, find the spline c2
       ! (c2*x^2) coefficients
       do i1 = 1,f%n(mod(d,3)+1)
          do i2 = 1,f%n(mod(d+1,3)+1)
             select case(d)
             case(1)
                fg(1:nmax) = f%f(:,i1,i2)
             case(2)
                fg(1:nmax) = f%f(i2,:,i1)
             case(3)
                fg(1:nmax) = f%f(i1,i2,:)
             end select

             ! constant terms for the system
             fone = fg(1)
             fprev = fg(nmax)
             do i = 1, nmax-1
                fnext = fg(i+1)
                fuse = fprev + fnext - fg(i) - fg(i)
                fprev = fg(i)
                fg(i) = fuse
             end do
             fg(nmax) = fprev + fone - fg(nmax) - fg(nmax)

             ! solve by direct substitution, L^t b = c
             fg(1) = fg(1) / l(1,1)
             do i = 2, nmax-1
                fg(i) = (fg(i) - l(i,i-1) * fg(i-1)) / l(i,i)
             end do
             fg(nmax) = (fg(nmax) - &
                dot_product(l(nmax,1:nmax-1),fg(1:nmax-1))) / l(nmax,nmax)

             ! again, L a = b
             fg(nmax) = fg(nmax) / l(nmax,nmax)
             fg(nmax-1) = (fg(nmax-1) - l(nmax,nmax-1) * fg(nmax)) / l(nmax-1,nmax-1)
             do i = nmax-2, 1, -1
                fg(i) = (fg(i) - l(i+1,i) * fg(i+1) - l(nmax,i)*fg(nmax)) / l(i,i)
             end do

             ! write down the second derivatives
             select case(d)
             case(1)
                f%c2(:,i1,i2,d) = fg(1:nmax)
             case(2)
                f%c2(i2,:,i1,d) = fg(1:nmax)
             case(3)
                f%c2(i1,i2,:,d) = fg(1:nmax)
             end select

          end do
       end do

    end do

    deallocate(l,fg)

  end subroutine init_trispline

  !> Given the electron density in the isrho slot, calculate the laplacian 
  !> in islap.
  subroutine grid_laplacian(frho,flap)
    use tools_io
    use tools_math
    use param
    use types

    type(field), intent(in) :: frho
    type(field), intent(out) :: flap

    integer :: n(3), i1, i2, i3

    complex(8) :: zaux(frho%n(1)*frho%n(2)*frho%n(3))
    real*8 :: bvec(3,3), vol
    integer :: ig, ntot
    integer :: j1, j2, j3
    integer, allocatable :: ivg(:,:), igfft(:)
    real*8, allocatable :: vgc(:,:)

    if (.not.frho%init) &
       call ferror('grid_laplacian','no density grid',faterr)

    ! allocate slot
    n = frho%n    
    flap = frho
    flap%usecore = .false.

    ! calculate the c2 coefficients again in the first call
    if (allocated(flap%c2)) deallocate(flap%c2)

    ntot = n(1) * n(2) * n(3)

    ! bvec
    bvec(:,1) = cross(frho%x2c(:,3),frho%x2c(:,2))
    bvec(:,2) = cross(frho%x2c(:,1),frho%x2c(:,3))
    bvec(:,3) = cross(frho%x2c(:,2),frho%x2c(:,1))
    vol = abs(det(frho%x2c))
    bvec = 2d0 * pi / vol * bvec

    allocate(ivg(3,ntot))
    ig = 0
    do i1 = n(1)/2-n(1)+1, n(1)/2
       do i2 = n(2)/2-n(2)+1, n(2)/2
          do i3 = n(3)/2-n(3)+1, n(3)/2
             ig = ig + 1
             ivg(1,ig)=i1
             ivg(2,ig)=i2
             ivg(3,ig)=i3
          end do
       end do
    end do
    allocate(igfft(ntot),vgc(3,ntot))
    do ig = 1, ntot
       i1=ivg(1,ig)
       i2=ivg(2,ig)
       i3=ivg(3,ig)
       if (i1.ge.0) then
          j1=i1
       else
          j1=n(1)+i1
       end if
       if (i2.ge.0) then
          j2=i2
       else
          j2=n(2)+i2
       end if
       if (i3.ge.0) then
          j3=i3
       else
          j3=n(3)+i3
       end if
       igfft(ig)=j3*n(2)*n(1)+j2*n(1)+j1+1
       vgc(:,ig)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
    end do

    zaux = 0d0
    zaux = reshape(frho%f,shape(zaux))
    call cfftnd(3,n,-1,zaux)

    do ig = 1, ntot
       zaux(igfft(ig)) = -dot_product(vgc(:,ig),vgc(:,ig)) * zaux(igfft(ig))
    end do

    call cfftnd(3,n,+1,zaux)
    flap%f = real(reshape(real(zaux,8),shape(flap%f)),8)

    deallocate(igfft,vgc,ivg)

  end subroutine grid_laplacian

  !> Calculate the gradient norm of a scalar field.
  subroutine grid_gradrho(frho,fgrho)
    use tools_io
    use tools_math
    use types
    use param

    type(field), intent(in) :: frho
    type(field), intent(out) :: fgrho

    integer :: n(3), i, i1, i2, i3, nr1, nr2, iq1, iq2, n12
    complex(8) :: zaux(frho%n(1)*frho%n(2)*frho%n(3))
    real*8 :: bvec(3,3), vol
    real*8 :: vgc
    integer :: ig, ntot
    integer :: j1, j2, j3, igfft

    if (.not.frho%init) &
       call ferror('grid_gradgrho','no density grid',faterr)

    ! allocate slot
    n = frho%n    
    fgrho = frho
    fgrho%usecore = .false.

    ! calculate the c2 coefficients again in the first call
    if (allocated(fgrho%c2)) deallocate(fgrho%c2)

    ntot = n(1) * n(2) * n(3)

    ! bvec
    bvec(:,1) = cross(frho%x2c(:,3),frho%x2c(:,2))
    bvec(:,2) = cross(frho%x2c(:,1),frho%x2c(:,3))
    bvec(:,3) = cross(frho%x2c(:,2),frho%x2c(:,1))
    vol = abs(det(frho%x2c))
    bvec = 2d0 * pi / vol * bvec

    n12 = n(1)*n(2)
    fgrho%f = 0d0
    do i = 1, 3
       zaux = reshape(frho%f,shape(zaux))
       call cfftnd(3,n,-1,zaux)

       do ig = 1, ntot
          iq1 = (ig-1) / (n(2)*n(3))
          i1 = iq1 - n(1)/2 + 1
          nr1 = (ig-1) - n(2)*n(3) * iq1
          iq2 = nr1 / n(3)
          i2 = iq2 - n(2)/2 + 1
          nr2 = nr1 - n(3) * iq2
          i3 = nr2 - n(3)/2 + 1
          if (i1.ge.0) then
             j1=i1
          else
             j1=n(1)+i1
          end if
          if (i2.ge.0) then
             j2=i2
          else
             j2=n(2)+i2
          end if
          if (i3.ge.0) then
             j3=i3
          else
             j3=n(3)+i3
          end if
          igfft=j3*n(2)*n(1)+j2*n(1)+j1+1
          vgc = dble(i1)*bvec(i,1)+dble(i2)*bvec(i,2)+dble(i3)*bvec(i,3)
          zaux(igfft) = vgc * cmplx(-aimag(zaux(igfft)),dble(zaux(igfft)),8)
       end do

       call cfftnd(3,n,+1,zaux)
       fgrho%f = fgrho%f + (real(reshape(real(zaux,8),shape(fgrho%f)),8))**2
    end do
    fgrho%f = sqrt(fgrho%f)

  end subroutine grid_gradrho

  !> Given the electron density in the isrho slot, calculate the
  !> promolecular density in isrhoat. If a fragment is passed, then
  !> only the atoms in it contribute. The itype can be 1: Hirshfeld
  !> weight, 2: rho_promolecular 3: rho_core. frho serves as a
  !> template; only the frho%n is used except if itype == 1.
  subroutine grid_rhoat(frho,frhoat,itype,fr)
    use grd_atomic
    use tools_io
    use tools_math
    use types
    use param

    type(field), intent(in) :: frho
    type(field), intent(inout) :: frhoat
    integer, intent(in) :: itype ! 1: hirsh weight, 2: rho_promolecular 3: rho_core
    type(fragment), intent(in), optional :: fr

    integer :: n(3), i, j, k
    real*8 :: x(3), xdelta(3,3), rdum1(3), rdum2(3,3)
    real*8 :: rhoat, rhocore, rfin

    n = frho%n    
    frhoat = frho
    if (.not.allocated(frhoat%f)) allocate(frhoat%f(n(1),n(2),n(3)))
    if (allocated(frhoat%c2)) deallocate(frhoat%c2)
    frhoat%usecore = .false.

    do i = 1, 3
       xdelta(:,i) = 0d0
       xdelta(i,i) = 1d0 / real(n(i),8)
    end do

    !$omp parallel do private(x,rhoat,rhocore,rdum1,rdum2,rfin,k,j,i)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)

             if (itype == 1) then
                call grda_promolecular(x,rhoat,rdum1,rdum2,0,.false.,fr)
                call grda_promolecular(x,rhocore,rdum1,rdum2,0,.true.,fr)
                rfin = (frho%f(i,j,k)+rhocore) / max(rhoat,1d-14)
             else if (itype == 2) then
                call grda_promolecular(x,rhoat,rdum1,rdum2,0,.false.,fr)
                rfin = rhoat
             else
                call grda_promolecular(x,rhocore,rdum1,rdum2,0,.true.,fr)
                rfin = rhocore
             end if
             !$omp critical(write)
             frhoat%f(i,j,k) = rfin
             !$omp end critical(write)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine grid_rhoat

  !> Given the electron density in the isrho slot, calculate the sign
  !> of the second eigenvalue of the hessian.
  subroutine grid_hxx(frho,fxx,ix)
    use tools_io
    use tools_math
    use param
    use types

    type(field), intent(in) :: frho
    type(field), intent(out) :: fxx
    integer, intent(in) :: ix

    integer :: n(3), i1, i2, i3

    complex(8) :: zaux(frho%n(1)*frho%n(2)*frho%n(3))
    real*8 :: bvec(3,3), vol
    integer :: ig, ntot
    integer :: j1, j2, j3
    integer, allocatable :: ivg(:,:), igfft(:)
    real*8, allocatable :: vgc(:,:)

    if (.not.frho%init) &
       call ferror('grid_hxx','no density grid',faterr)

    ! allocate slot
    n = frho%n    
    fxx = frho
    fxx%usecore = .false.

    ! calculate the c2 coefficients again in the first call
    if (allocated(fxx%c2)) deallocate(fxx%c2)

    ntot = n(1) * n(2) * n(3)

    ! bvec
    bvec(:,1) = cross(frho%x2c(:,3),frho%x2c(:,2))
    bvec(:,2) = cross(frho%x2c(:,1),frho%x2c(:,3))
    bvec(:,3) = cross(frho%x2c(:,2),frho%x2c(:,1))
    vol = abs(det(frho%x2c))
    bvec = 2d0 * pi / vol * bvec

    allocate(ivg(3,ntot))
    ig = 0
    do i1 = n(1)/2-n(1)+1, n(1)/2
       do i2 = n(2)/2-n(2)+1, n(2)/2
          do i3 = n(3)/2-n(3)+1, n(3)/2
             ig = ig + 1
             ivg(1,ig)=i1
             ivg(2,ig)=i2
             ivg(3,ig)=i3
          end do
       end do
    end do
    allocate(igfft(ntot),vgc(3,ntot))
    do ig = 1, ntot
       i1=ivg(1,ig)
       i2=ivg(2,ig)
       i3=ivg(3,ig)
       if (i1.ge.0) then
          j1=i1
       else
          j1=n(1)+i1
       end if
       if (i2.ge.0) then
          j2=i2
       else
          j2=n(2)+i2
       end if
       if (i3.ge.0) then
          j3=i3
       else
          j3=n(3)+i3
       end if
       igfft(ig)=j3*n(2)*n(1)+j2*n(1)+j1+1
       vgc(:,ig)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
    end do

    zaux = 0d0
    zaux = reshape(frho%f,shape(zaux))
    call cfftnd(3,n,-1,zaux)

    do ig = 1, ntot
       zaux(igfft(ig)) = -vgc(ix,ig) * vgc(ix,ig) * zaux(igfft(ig))
    end do

    call cfftnd(3,n,+1,zaux)
    fxx%f = real(reshape(real(zaux,8),shape(fxx%f)),8)

    deallocate(igfft,vgc,ivg)

  end subroutine grid_hxx

end module grid_tools
