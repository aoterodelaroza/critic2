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

!> Dirty hacks for dirty jobs. Abandon hope all ye who enter here.
module tricks
  implicit none

  private
  public :: trick
  private :: trick_recalculate_xdm
  private :: trick_grid_sphere

contains

  subroutine trick(line0)
    use tools_io
    use struct_basic
    character*(*), intent(in) :: line0

    ! call trick_recalculate_xdm()
    call trick_grid_sphere()

  end subroutine trick

  !> Recalculate the XDM dispersion energy from the information in the 
  !> espresso output from which the crystal structure was read.
  subroutine trick_recalculate_xdm()
    use tools_io
    use struct_basic

    integer :: i, j, ii, jj, nn, i3
    integer :: lu, idx, idx1, idx2
    character(len=:), allocatable :: line, str
    logical :: ok
    real*8 :: a1, a2, maxc6, rmax, rmax2, x(3), ri, ri2, ri3
    real*8 :: cn0, rvdwx, rix, exx, eat, etotal
    real*8 :: c6(cr%ncel,cr%ncel), c8(cr%ncel,cr%ncel), c10(cr%ncel,cr%ncel)
    real*8 :: rc(cr%ncel,cr%ncel), rvdw(cr%ncel,cr%ncel)
    
    real*8, parameter :: ecut = 1d-11

    write (uout,'("* Trick: sum the XDM dispersion energy")')
    write (uout,'("+ Reading coefficients from the file: ",A)') string(cr%file)

    lu = fopen_read(string(cr%file))
    call cr%build_env(.false.,150d0)

    main: do while (getline(lu,line))
       ! read the parameters
       if (trim(line) == "* XDM dispersion") then
          ok = getline(lu,line)
          idx = index(line,"=")
          str = trim(line(idx+1:))
          read(str,*) a1
          ok = getline(lu,line)
          ok = getline(lu,line)
          idx = index(line,"=")
          str = trim(line(idx+1:))
          read(str,*) a2
       end if
       ! read the dispersion coefficients and the radii
       if (trim(line) == "+ Dispersion coefficients") then
          do i = 1, cr%ncel
             do j = 1, i
                ok = getline(lu,line)
                read(line,*) idx1, idx2, c6(i,j), c8(i,j), c10(i,j), rc(i,j), rvdw(i,j)
                if (idx1 /= i .or. idx2 /= j) then
                   write (*,*) idx1, i
                   write (*,*) idx2, j
                   call ferror("trick","read indices do not match",faterr)
                end if
                c6(j,i) = c6(i,j)
                c8(j,i) = c8(i,j)
                c10(j,i) = c10(i,j)
                rc(j,i) = rc(i,j)
                rvdw(j,i) = rvdw(i,j)
             end do
          end do
       end if
    end do main
    call fclose(lu)

    ! do stuff
    do i = 1, cr%ncel
       if (cr%atcel(i)%x(3) < 0.25d0) then
          c6(i,:) = 0d0
          c6(:,i) = 0d0
          c8(i,:) = 0d0
          c8(:,i) = 0d0
          c10(i,:) = 0d0
          c10(:,i) = 0d0
       endif
    end do

    ! calculate the energy
    maxc6 = maxval(c6)
    ! set the atomic environment for the sum
    rmax = (maxc6/ecut)**(1d0/6d0)
    rmax2 = rmax*rmax
    etotal = 0d0
    do i = 1, cr%ncel
       eat = 0d0
       do jj = 1, cr%nenv
          j = cr%atenv(jj)%cidx
          x = cr%atenv(jj)%r - cr%atcel(i)%r
          ri2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
          if (ri2 < 1d-15 .or. ri2>rmax2) cycle
          ri = sqrt(ri2)
          ri3 = ri2 * ri
    
          do nn = 6, 10, 2
             ! order R^nn
             i3 = nn / 2 - 1
    
             ! check the ii and the jj, also how the env is generated
             if (nn == 6) then
                cn0 = c6(i,j)
             elseif (nn == 8) then
                cn0 = c8(i,j)
             elseif (nn == 10) then
                cn0 = c10(i,j)
             end if
    
             ! energy contribution
             rvdwx = rvdw(i,j)**nn
             rix = ri**nn
             exx = cn0 / (rvdwx + rix)
             eat = eat + exx
          end do
       end do ! jj
       etotal = etotal + eat 
    end do ! ii
    etotal= -0.5d0 * etotal
    
    write (uout,'("  Evdw = ",A," Hartree, ",A," Ry")') &
       string(etotal,'e',decimal=10), string(etotal*2,'e',decimal=10)
    
    write (uout,*)

  end subroutine trick_recalculate_xdm

  !> Test for the calculation of energies using a sphere plus grid approach.
  subroutine trick_grid_sphere()
    use bisect
    use rhoplot
    use global
    use varbas
    use fields
    use bader
    use integration
    use yt
    use struct_basic
    use grid_tools
    use arithmetic
    use tools_math
    use tools_io
    use types
    use param

    integer, parameter :: nleb = 770, nr = 50
    integer, parameter :: nx = 100

    character(len=:), allocatable :: saux
    integer :: i, j, k, l, idx, idf, ids
    integer :: ix, iy, iz
    real*8 :: rsph(cr%ncel), xsphf(cr%ncel), xsphs(cr%ncel)
    real*8 :: xleb(nleb), yleb(nleb), zleb(nleb), wleb(nleb)
    real*8 :: rp(nr), rw(nr), lprop, x(3), x0(3), unit(3)
    real*8 :: ff, fs, rsumf, rsums, ssumf, ssums
    real*8 :: gsums, gsumf, d2, atp(cr%ncel), ntot
    type(scalar_value) :: res, res2
    ! for yt
    integer :: nbasin, luw
    real*8, allocatable :: xcoord(:,:)
    integer, allocatable :: idg(:,:,:), icp(:)
    real(kind=gk), allocatable :: w(:,:,:)

    ntot = f(2)%n(1)*f(2)%n(2)*f(2)%n(3)

    !! model input:
    ! molecule benzene.wfx 5
    ! load benzene.wfx
    ! load rho.cube
    ! load lap.cube
    ! trick

    ! integrate with yt
    write (*,*) "Doing YT on field 2"
    call yt_integrate(cr,f(2),nbasin,xcoord,idg,luw)
    call int_reorder_gridout(cr,f(2),nbasin,xcoord,idg,2d0,luw,icp)
    ! call bader_integrate(cr,f(2),nbasin,xcoord,idg)

    ! defining the weight fields
    call realloc(fused,25)
    call realloc(f,25)
    allocate(w(f(2)%n(1),f(2)%n(2),f(2)%n(3)))
    do i = 1, cr%ncel
       call yt_weights(luw,i,w)
       call grid_from_array3(w,f(10+i))
       fused(10+i) = .true.
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

end module tricks
