! based on topologically correct marching cubes (see weber
! et al.'s Exploring Scalar Fields Using Critical Isovalues,
! http://dx.doi.org/10.1109/VISUAL.2002.1183772
! Visualization, 2002, 171-178
  ! xxxx
  subroutine grid_cplinear(ff,nmax,xmax,nmin,xmin,nsad,xsad)
    use tools_io
    use types
    type(field), intent(in) :: ff
    integer, intent(out) :: nmax, nmin, nsad
    real*8, allocatable :: xmax(:,:), xmin(:,:), xsad(:,:)
    
    integer, parameter :: m0 = 10

    integer :: i, j, k, l, idx(3), didx(3)
    real*8 :: g0, gx(-1:1,-1:1,-1:1)
    real*8 :: a, b, c, d, e, f, g, h, k0, k1, k2, disc
    real*8 :: x(3)

    integer, parameter :: istar(3,20) = reshape((/&
       -1,  0,  0,&
        0, -1,  0,&
        0,  0, -1,&
        0,  0,  0,&
        1,  0,  0,&
        0,  1,  0,&
        0,  0,  1,&
        1,  1,  0,&
        1,  0,  1,&
        0,  1,  1,&
        1,  1,  1,&
        1,  0, -1,&
        0,  1, -1,&
        1,  1, -1,&
        1, -1,  0,&
        0, -1,  1,&
        1, -1,  1,&
       -1,  1,  0,&
       -1,  0,  1,&
       -1,  1,  1&
       /),shape(istar))

    if (.not.ff%init) &
       call ferror('grid_cplinear','grid not initialized',faterr)
       
    if (allocated(xmax)) deallocate(xmax)
    if (allocated(xmin)) deallocate(xmin)
    if (allocated(xsad)) deallocate(xsad)
    allocate(xmax(3,m0),xmin(3,m0),xsad(3,m0))
    nmax = 0
    nmin = 0
    nsad = 0
    gx = 0d0
    do i = 1, ff%n(1)
       idx(1) = i
       do j = 1, ff%n(2)
          idx(2) = j
          do k = 1, ff%n(3)
             idx(3) = k

             ! collect the necessary values from the grid
             do l = 1, size(istar,2)
                didx = modulo(idx - 1 + istar(:,l),ff%n)+1
                gx(istar(1,l),istar(2,l),istar(3,l)) = ff%f(didx(1),didx(2),didx(3))
             end do

             ! check if this is a maximum
             if (gx(0,0,0) > gx(1,0,0)  .and. gx(0,0,0) > gx(0,1,0)  .and. gx(0,0,0) > gx(0,0,1) .and. &
                 gx(0,0,0) > gx(-1,0,0) .and. gx(0,0,0) > gx(0,-1,0) .and. gx(0,0,0) > gx(0,0,-1)) then
                nmax = nmax + 1
                if (nmax > size(xmax,2)) call realloc(xmax,3,2*nmax)
                xmax(:,nmax) = real(idx-1,8) / ff%n
             end if

             ! check if this is a minimum
             if (gx(0,0,0) < gx(1,0,0)  .and. gx(0,0,0) < gx(0,1,0)  .and. gx(0,0,0) < gx(0,0,1) .and. &
                 gx(0,0,0) < gx(-1,0,0) .and. gx(0,0,0) < gx(0,-1,0) .and. gx(0,0,0) < gx(0,0,-1)) then
                nmin = nmin + 1
                if (nmin > size(xmin,2)) call realloc(xmin,3,2*nmin)
                xmin(:,nmin) = real(idx-1,8) / ff%n
             end if

             ! trilinear:
             ! T = a + e*x + c*y + b*z + g*x*y + f*x*z + d*y*z + h*x*y*z

             ! check if there is a body saddle point
             a = gx(0,0,0)
             b = gx(0,0,1)-a
             c = gx(0,1,0)-a
             e = gx(1,0,0)-a
             g = gx(1,1,0)-a-e-c
             f = gx(1,0,1)-a-e-b
             d = gx(0,1,1)-a-c-b
             h = gx(1,1,1)-a-b-c-d-e-f-g

             k0 = c * f - b * g
             k1 = d * f - b * h
             k2 = d * g - c * h
             if (k1 < 0) cycle

             disc = g*g - h/k1 * (e*k2 + g*k0)
             if (disc < 0) cycle

             ! first root
             x(3) = (-g+sqrt(disc))/h
             x(2) = (k0 + k1 * x(3)) / k2
             x(1) = -(c + d * x(3)) / (g + h * x(3))
             if (all(x < 0d0) .and. all(x > 1d0)) then
                nsad = nsad + 1
                if (nsad > size(xsad,2)) call realloc(xsad,3,2*nsad)
                xsad(:,nsad) = (idx-1+x) / ff%n 
             end if

             ! second root
             x(3) = (-g-sqrt(disc))/h
             x(2) = (k0 + k1 * x(3)) / k2
             x(1) = -(c + d * x(3)) / (g + h * x(3))
             if (all(x < 0d0) .and. all(x > 1d0)) then
                nsad = nsad + 1
                if (nsad > size(xsad,2)) call realloc(xsad,3,2*nsad)
                xsad(:,nsad) = (idx-1+x) / ff%n 
             end if

             ! face saddle, xy
             x = facesad(3,gx(0,0,0),gx(1,0,0),gx(0,1,0),gx(1,1,0),&
                           gx(0,0,-1),gx(1,0,-1),gx(0,1,-1),gx(1,1,-1),&
                           gx(0,0,1),gx(1,0,1),gx(0,1,1),gx(1,1,1))
             if (x(1) >= 0d0 .and. x(1) <= 1d0 .and. x(2) >= 0d0 .and. x(2) <= 1d0) then
                nsad = nsad + 1
                if (nsad > size(xsad,2)) call realloc(xsad,3,2*nsad)
                xsad(:,nsad) = (idx-1+x) / ff%n 
             end if

             ! face saddle, xz
             x = facesad(2,gx(0,0,0),gx(1,0,0),gx(0,0,1),gx(1,0,1),&
                           gx(0,-1,0),gx(1,-1,0),gx(0,-1,1),gx(1,-1,1),&
                           gx(0,1,0),gx(1,1,0),gx(0,1,1),gx(1,1,1))
             if (x(1) >= 0d0 .and. x(1) <= 1d0 .and. x(3) >= 0d0 .and. x(3) <= 1d0) then
                nsad = nsad + 1
                if (nsad > size(xsad,2)) call realloc(xsad,3,2*nsad)
                xsad(:,nsad) = (idx-1+x) / ff%n 
             end if

             ! face saddle, yz
             x = facesad(1,gx(0,0,0),gx(0,1,0),gx(0,0,1),gx(0,1,1),&
                           gx(-1,0,0),gx(-1,1,0),gx(-1,0,1),gx(-1,1,1),&
                           gx(1,0,0),gx(1,1,0),gx(1,0,1),gx(1,1,1))
             if (x(2) >= 0d0 .and. x(2) <= 1d0 .and. x(3) >= 0d0 .and. x(3) <= 1d0) then
                nsad = nsad + 1
                if (nsad > size(xsad,2)) call realloc(xsad,3,2*nsad)
                xsad(:,nsad) = (idx-1+x) / ff%n 
             end if
          end do
       end do
    end do

  end subroutine grid_cplinear

  ! bilinear:
  ! T = a*(1-x)*(1-y) + b*x*(1-y) + c*(1-x)*y + d*x*y
  function facesad(ix,a,b,c,d,am1,bm1,cm1,dm1,ap1,bp1,cp1,dp1) result(x)
    integer, intent(in) :: ix
    real*8, intent(in) :: am1,bm1,cm1,dm1,a,b,c,d,ap1,bp1,cp1,dp1
    real*8 :: x(3)
    
    real*8 :: s1, s2, y(2)
    
    ! critical point of the bilinear interpolant
    y(1) = 1d0 / (1d0 + (b-d)/(c-a))
    y(2) = 1d0 / (1d0 + (c-d)/(b-a))
    x = 0d0
    if (ix == 1) then
       x(2:3) = y
    elseif (ix == 3) then
       x(1:2) = y
    else
       x(1) = y(1)
       x(3) = y(2)
    end if

    ! sign condition - both have to have the same sign to be a 
    ! critical point of the trilinear interpolant
    s1 = c*(ap1-a) + a*(cp1-c) - d*(bp1-b) - b*(dp1-d)
    s2 = c*(am1-a) + a*(cm1-c) - d*(bm1-b) - b*(dm1-d)
    if (s1 * s2 < 0d0) x = -1d0

  end function facesad

  !> xxxxx
  subroutine autogrid(line)
    use fields
    use struct_basic
    use grid_tools
    use global
    use tools_io

    character*(*), intent(in) :: line
    integer :: nmax, nmin, nsad
    real*8, allocatable :: xmax(:,:), xmin(:,:), xsad(:,:)
    real*8 :: x0(3), cpeps0
    integer :: i
    character(len=:), allocatable :: discexpr

    ! initialize
    if (.not.quiet) then
       call tictac("Start AUTOGRID")
       write (uout,*)
    end if
    discexpr = ""

    ! save the cp eps
    cpeps0 = CP_eps_cp
    CP_eps_cp = 2d0 * maxval(cr%aa / f(refden)%n)

    ! Check that the reference field makes sense
    if (f(refden)%type /= type_grid) then
       call ferror("autogrid","AUTOGRID can be used with grids only",faterr,syntax=.true.)
       return
    end if

    call grid_cplinear(f(refden),nmax,xmax,nmin,xmin,nsad,xsad)

    do i = 1, nmax
       x0 = cr%x2c(xmax(:,i))
       call addcp(x0,discexpr,-3)
    end do
    
    do i = 1, nmin
       x0 = cr%x2c(xmin(:,i))
       call addcp(x0,discexpr,3)
    end do

    do i = 1, nsad
       x0 = cr%x2c(xsad(:,i))
       call addcp(x0,discexpr)
    end do

    ! Short report of non-equivalent cp-list
    call cp_short_report()
    call critshell(1)

    ! put back the cp eps
    CP_eps_cp = cpeps0

    if (.not.quiet) then
       call tictac("End AUTOGRID")
       write (uout,*)
    end if

  end subroutine autogrid

     ! autogrid
     elseif (equal(word,'autogrid')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before autogrid',faterr,line,syntax=.true.)
           cycle
        end if
        call autogrid(subline)

