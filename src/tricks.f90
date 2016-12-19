! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> Dirty hacks for dirty jobs. Abandon hope all ye who enter here.
module tricks
  implicit none

  private
  public :: trick
  private :: trick_grid_sphere
  private :: trick_stephens_nnm_channel
  private :: trick_cell_integral

contains

  subroutine trick(line0)
    character*(*), intent(in) :: line0

    ! call trick_recalculate_xdm()
    ! call trick_grid_sphere()
    ! call trick_stephens_nnm_channel(line0)
    call trick_cell_integral()
    
  end subroutine trick

  !> Test for the calculation of energies using a sphere plus grid approach.
  subroutine trick_grid_sphere()
    use rhoplot, only: rhoplot_cube
    use global, only: refden
    use fields, only: f, fused, type_grid, grd
    use integration, only: int_reorder_gridout
    use yt, only: yt_integrate, yt_weights
    use struct_basic, only: cr
    use grid_tools, only: grid_from_array3
    use tools_math, only: select_lebedev, gauleg
    use tools_io, only: string
    use types, only: scalar_value, realloc
    use param, only: fh, pi

    integer, parameter :: nleb = 770, nr = 50

    integer :: i, j, k, l
    integer :: ix, iy, iz
    real*8 :: rsph(cr%ncel)
    real*8 :: xleb(nleb), yleb(nleb), zleb(nleb), wleb(nleb)
    real*8 :: rp(nr), rw(nr), x(3), x0(3), unit(3)
    real*8 :: ff, fs, rsumf, rsums, ssumf, ssums
    real*8 :: d2, atp(cr%ncel), ntot
    type(scalar_value) :: res, res2
    ! for yt
    integer :: nbasin, luw
    real*8, allocatable :: xcoord(:,:)
    integer, allocatable :: idg(:,:,:), icp(:)
    real*8, allocatable :: w(:,:,:)

    ntot = f(2)%n(1)*f(2)%n(2)*f(2)%n(3)

    !! model input:
    ! molecule benzene.wfx 5
    ! load benzene.wfx
    ! load rho.cube
    ! load lap.cube
    ! trick

    ! integrate with yt
    write (*,*) "Doing YT on field 2"
    call yt_integrate(cr,f(2)%f,"",.true.,2d0,nbasin,xcoord,idg,luw)
    call int_reorder_gridout(cr,f(2),nbasin,xcoord,idg,.true.,2d0,luw,icp)
    ! call bader_integrate(cr,f(2),nbasin,xcoord,idg)

    ! defining the weight fields
    call realloc(fused,25)
    call realloc(f,25)
    allocate(w(f(2)%n(1),f(2)%n(2),f(2)%n(3)))
    do i = 1, cr%ncel
       call yt_weights(luw=luw,idb=i,w=w)
       call grid_from_array3(w,f(10+i))
       fused(10+i) = .true.
       call fh%put(string(10+i),10+i)
       f(10+i)%type = type_grid
    end do
    close(luw)
    deallocate(w)

    call rhoplot_cube("grid field 11 file bleh.cube")

    ! determine the radius for each atom
    write (*,*) "Calculating the radii"
    do i = 1, cr%ncel
       rsph(i) = cr%at(cr%atcel(i)%idx)%rnn2 * 0.99d0
    end do

    ! field 4 will contain the smooth grid
    f(4) = f(3)
    fused(4) = .true.
    call fh%put("4",4)

    ! smooth out the regions inside the spheres
    write (*,*) "Calculating the smooth grid"
    !$omp parallel do private(x0,x,d2) schedule(guided)
    do ix = 1, f(4)%n(1)
       x0(1) = real(ix-1,8) / f(4)%n(1)
       do iy = 1, f(4)%n(2)
          x0(2) = real(iy-1,8) / f(4)%n(2)
          do iz = 1, f(4)%n(3)
             x0(3) = real(iz-1,8) / f(4)%n(3)
             do i = 1, cr%ncel
                x = cr%atcel(i)%x - x0
                call cr%shortest(x,d2)
                if (d2 <= rsph(i)*rsph(i)) then
                   !$omp critical (smooth)
                   f(4)%f(ix,iy,iz) = f(4)%f(ix,iy,iz) * sin(0.5d0*pi*sqrt(d2)/rsph(i))**6
                   !$omp end critical (smooth)
                   exit
                end if
             end do
          end do
       end do
    end do
    !$omp end parallel do

    ! integrate the basins
    atp = 0d0
    call select_lebedev(nleb,xleb,yleb,zleb,wleb)
    do i = 1, cr%ncel
       write (*,*) "integrating atom ", i
       ! integrate the smooth field
       atp(i) = sum(f(10+i)%f * f(4)%f) * cr%omega / ntot

       ! calculate the difference between the smooth field and the actual field
       ! inside the basin. Contributions from all the spheres.
       do j = 1, cr%ncel
          ssumf = 0d0
          ssums = 0d0
          !$omp parallel do private(unit,rp,rw,rsumf,rsums,x,fs,ff,res,res2) schedule(guided)
          do k = 1, nleb
             unit = (/xleb(k),yleb(k),zleb(k)/)
             call gauleg(0d0,rsph(j),rp,rw,nr)
             rsumf = 0d0
             rsums = 0d0
             do l = 1, nr
                x = cr%atcel(j)%r + rp(l) * unit
                call grd(f(refden),x,2,res)
                call grd(f(10+i),x,2,res2)
                ff = res%del2f * res2%f
                fs = ff * sin(0.5d0*pi*rp(l)/rsph(j))**6
                rsumf = rsumf + rp(l)**2 * rw(l) * ff
                rsums = rsums + rp(l)**2 * rw(l) * fs
             end do
             !$omp critical (acum_ang)
             ssumf = ssumf + rsumf * wleb(k)
             ssums = ssums + rsums * wleb(k)
             !$omp end critical (acum_ang)
          end do
          !$omp end parallel do
          write (*,*) "sphere ", j, ssumf, ssums
          atp(i) = atp(i) + ssumf - ssums
       end do
       write (*,*) "atomic laplacian: ", atp(i)
    end do

    write (*,*) "all ", sum(atp)

  end subroutine trick_grid_sphere

  subroutine trick_stephens_nnm_channel(line)
    use struct_basic, only: cr
    use fields, only: f, grd0
    use graphics, only: graphics_open, graphics_ball, graphics_stick, graphics_close
    use global, only: eval_next
    use tools_math, only: norm, cross
    use tools_io, only: faterr, ferror, uout, string
    use param, only: bohrtoa, pi
    character*(*), intent(in) :: line
    
    ! parameters for input
    integer, parameter :: nang = 30
    integer, parameter :: nline = 30
    real*8, parameter :: rho0_def = 0.01d0
    real*8, parameter :: zeronorm = 1d-10 ! zero norm 
    real*8, parameter :: dmax = 10d0 ! maximum distance
    ! bracketing and bisection
    real*8, parameter :: bstep0 = 0.05d0
    real*8, parameter :: bfactor = 1.2d0
    real*8, parameter :: blimit = 1d-5

    integer :: luobj, lumtl, lp, i, j
    logical :: ok
    real*8, dimension(3) :: x0, x1, x, xp1, xp2, a1, a2, xd, xa, xb
    real*8 :: ang, rho, rhoi, step, dist, rhoa, rhob
    real*8 :: xlimit(3,nang,nline), mindist, maxdist
    real*8 :: rho0

    ! read input
    lp = 1
    ok = eval_next(x0(1),line,lp)
    ok = ok .and. eval_next(x0(2),line,lp)
    ok = ok .and. eval_next(x0(3),line,lp)
    ok = ok .and. eval_next(x1(1),line,lp)
    ok = ok .and. eval_next(x1(2),line,lp)
    ok = ok .and. eval_next(x1(3),line,lp)
    if (.not.ok) &
       call ferror('trick_stephens_nnm_channel','syntax: trick x0 y0 z0 x1 y1 z1',faterr)
    ok = eval_next(rho0,line,lp)
    if (.not.ok) rho0 = rho0_def

    ! convert to cartesian
    x0 = cr%x2c(x0)
    x1 = cr%x2c(x1)

    ! calculate the promolecular density profile
    ok = .true.
    write (uout,'("+ Promolecular density along the line")')
    write (uout,'("# Distance(ang)  rho(promol) (a.u.)")')
    do i = 1, nline
       x = x0 + real(i-1,8) / real(nline-1,8) * (x1-x0)
       rho = grd0(f(0),x)
       write (uout,'(4X,99(A,X))') string(norm(x-x0)*bohrtoa,'f',14,6), string(rho,'e',12,6)
       ok = ok .and. rho < rho0
    end do
    if (.not.ok) &
       call ferror('trick_stephens_nnm_channel','promolecular density along the line higher than threshold',faterr)

    ! vector in the line direction, normalized
    a1 = x1-x0
    a1 = a1 / norm(a1)

    write (uout,'("# Id  Distance(ang)  Avg. diam(ang)  Circle area(ang^2)")')
    do i = 1, nline
       ! position and density at this point along the line
       x = x0 + real(i-1,8) / real(nline-1,8) * (x1-x0)
       rhoi = grd0(f(0),x)

       ! find a vector normal to x-x0
       ok = .false.
       do j = 1, 3 
          a2 = 0d0
          a2(j) = 1d0
          xp1 = cross(a1,a2)
          if (norm(xp1) > zeronorm) then
             ok = .true.
             exit 
          end if
       end do
       if (.not.ok) &
          call ferror('trick_stephens_nnm_channel','could not find a vector orthogonal to the line',faterr)
       xp1 = xp1 / norm(xp1)
       
       ! find the other normal vector
       xp2 = cross(a1,xp1)
       if (norm(xp2) < zeronorm) &
          call ferror('trick_stephens_nnm_channel','could not find the second orthogonal vector',faterr)
       xp2 = xp2 / norm(xp2)
       
       ! run over the angles
       mindist = 1d30
       maxdist = 0d0
       do j = 1, nang
          ang = 2d0 * pi * real(j-1,8) / real(nang,8)

          ! direction along this angle
          xd = cos(ang) * xp1 + sin(ang) * xp2

          ! bracket
          step = bstep0 / bfactor
          xb = x
          rhob = rhoi
          do while (rhob <= rho0)
             step = step * bfactor
             xa = xb
             rhoa = rhob
             xb = x + xd * step
             rhob = grd0(f(0),xb)
             if (step > dmax) then
                xb = x + xd * dmax
                exit
             endif
          end do

          ! bisection
          dist = norm(xb - xa)
          do while (dist >= blimit)
             xd = 0.5d0 * (xb + xa)
             rho = grd0(f(0),xb)
             if (rho > rho0) then
                rhob = rho
                xb = xd
             else
                rhoa = rho
                xa = xd
             endif
             dist = norm(xb - xa)
          end do
          dist = norm(xd-x)
          maxdist = max(maxdist,dist)
          mindist = min(mindist,dist)
          xlimit(:,j,i) = xd
       end do

       ! Calculate the average radius and output
       dist = 0d0
       do j = 1, nang
          dist = dist + norm(xlimit(:,j,i) - x)
       end do
       dist = dist / nang
       write (uout,'(2X,99(A,X))') string(i), string(norm(x-x0)*bohrtoa,'f',14,6,4), &
          string(2d0*dist*bohrtoa,'f',14,6,4), string(pi*dist*dist*bohrtoa*bohrtoa,'f',14,6,4)
    end do

    ! write the obj file
    call graphics_open("obj","tube.obj",luobj,lumtl)
    do i = 1, nline
       x = x0 + real(i-1,8) / real(nline-1,8) * (x1-x0)
       call graphics_ball("obj",luobj,x,(/0,0,255/),0.2d0)
    end do
    do i = 1, nang
       call graphics_stick("obj",luobj,x0,xlimit(:,i,1),(/0,0,0/),0.1d0)
       call graphics_stick("obj",luobj,x1,xlimit(:,i,nline),(/0,0,0/),0.1d0)
       do j = 1, nline-1
          call graphics_stick("obj",luobj,xlimit(:,i,j),xlimit(:,i,j+1),(/0,0,0/),0.05d0)
       end do
    end do
    do i = 1, nline
       do j = 1, nang-1
          call graphics_stick("obj",luobj,xlimit(:,j,i),xlimit(:,j+1,i),(/0,0,0/),0.05d0)
       end do
       call graphics_stick("obj",luobj,xlimit(:,1,i),xlimit(:,nang,i),(/0,0,0/),0.05d0)
    end do
    call graphics_close("obj",luobj,lumtl)

  end subroutine trick_stephens_nnm_channel

  !> Calculate the cell integral of the reference field using 
  !> Franchini et al.'s Becke-style mesh
  subroutine trick_cell_integral()
    use fields, only: f
    use struct_basic, only: cr
    use meshmod, only: genmesh, fillmesh
    use global, only: refden
    use types, only: molmesh
    use tools_io, only: uout, string
    use param, only: im_rho

    type(molmesh) :: m
    integer :: prop(1), id(1)

    write (uout,'("* Trick: cell integral")')
    m = genmesh(cr)
    write (uout,'("  mesh size      ",A)') string(m%n)
    id(1) = 1
    prop(1) = im_rho
    call fillmesh(m,f(refden),id,prop,.not.cr%ismolecule)
    write (uout,'("cell integral ",A)') string(sum(m%f(:,1) * m%w),'f',12,6)

  end subroutine trick_cell_integral

end module tricks
