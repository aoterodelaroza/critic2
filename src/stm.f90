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

! This module is based on the code kindly provided by Enrico Benassi,
! Scuola Normale Superiore di Pisa, Italy <ebenassi3@gmail.com>

! Calculate the STM map of the surface slab in the input using the
! Tersoff-Hamann approximation and either a constant height or a constant
! current. The reference field is assumed to be the DOS at the Fermi level
! (although the density will also give reasonable results). 
module stm
  implicit none

  private

  public :: stm_driver
  private :: detect_vacuum
  private :: stm_bisect
  private :: stm_bisect_grid

contains

  
  subroutine stm_driver(line)
    use fields, only: f, type_grid, grd0, fields_typestring
    use crystalmod, only: cr
    use global, only: refden, eval_next, fileroot
    use tools_io, only: ferror, faterr, uout, lgetword, equal, string,&
       fopen_write, fclose
    use param, only: bohrtoa

    character*(*), intent(inout) :: line

    character(len=:), allocatable :: word, fname
    logical :: ok, iscur, doline, isgrid
    integer :: fid, lp, i, j, i1, i2, idx
    integer :: ix, nx(2), nn, ip1, ip2, np1, np2
    real*8 :: rtop, rtop0, rcur, rhei, rtmp, x0(2), x1(2), xx(3), xr(3), axlen
    real*8 :: faux, xmin, xmax, ymin, ymax, dist
    real*8, allocatable :: ff(:,:), fl(:)
    character*1, parameter :: laxis(3) = (/"a","b","c"/)

    if (cr%ismolecule) then
       call ferror("stm_driver","STM can not be used with molecules",faterr,syntax=.true.)
       return
    end if

    ! header
    write (uout,'("* Scanning tunneling microscopy plots ")')

    ! defaults
    iscur = .true.
    rcur = 1d-3
    rhei = -1d0
    rtop0 = -1d0
    nx = 1
    doline = .false.
    x0 = 0d0
    x1 = 1d0
    nn = 201
    if (f(refden)%type == type_grid) then
       np1 = -1
       np2 = -1
    else
       np1 = 51
       np2 = 51
    endif

    ! parse the input and check sanity
    lp = 1
    do while(.true.)
       word = lgetword(line,lp)
       if (equal(word,"current")) then
          iscur = .true.
          ok = eval_next(rtmp,line,lp)
          if (ok) rcur = rtmp
       elseif (equal(word,"height")) then
          iscur = .false.
          ok = eval_next(rtmp,line,lp)
          if (ok) rhei = rtmp
       elseif (equal(word,"top")) then
          ok = eval_next(rtop0,line,lp)
          if (.not.ok) then
             call ferror('stm_driver','Wrong syntax in STM/TOP',faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,"npts")) then
          ok = eval_next(np1,line,lp)
          ok = ok .and. eval_next(np2,line,lp)
          if (.not.ok) then
             call ferror('stm_driver','Wrong syntax in STM/NPTS',faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,"cells").or.equal(word,"cell")) then
          ok = eval_next(nx(1),line,lp)
          ok = ok .and. eval_next(nx(2),line,lp)
          if (.not.ok) then
             call ferror('stm_driver','Wrong syntax in STM/CELLS',faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,"line")) then
          doline = .true.
          ok = eval_next(x0(1),line,lp)
          ok = ok .and. eval_next(x0(2),line,lp)
          ok = ok .and. eval_next(x1(1),line,lp)
          ok = ok .and. eval_next(x1(2),line,lp)
          ok = ok .and. eval_next(nn,line,lp)
          if (.not.ok) then
             call ferror('stm_driver','Wrong syntax in STM/LINE',faterr,syntax=.true.)
             return
          end if
       else
          exit
       end if
    end do

    ! detect the vacuum direction, assign the two in-plane directions
    call detect_vacuum(ix,rtop)
    ip1 = 0
    do i = 1, 3
       if (i /= ix) then
          if (ip1 == 0) then
             ip1 = i
          else
             ip2 = i
          endif
       endif
    end do
    axlen = cr%aa(ix)
    if (rhei < 0d0) rhei = rtop
    if (rtop0 < 0d0) rtop0 = rtop
       
    if (.not.doline) then
       ! Calculate the STM signal in a plane
       ! Allocate the 2d plane
       if (np1 == -1 .and. np2 == -1) then
          np1 = f(refden)%n(ip1)
          np2 = f(refden)%n(ip2)
          isgrid = .true.
       else
          isgrid = .false.
       endif
       allocate(ff(np1,np2))

       ! Calculate the values on the plane
       if (iscur) then
          write (uout,'("  Mode: constant current, ",A," density.")') &
             trim(string(rcur,'e',9,3))
          write (uout,'("  In-plane points: ",A,", ",A)') &
             string(np1), string(np2)
          !$omp parallel do private (xx,faux) schedule(dynamic)
          do i = 1, np1
             do j = 1, np2
                xx(ip1) = real(i-1,8)/np1
                xx(ip2) = real(j-1,8)/np2
                xx(ix) = rtop + 0.2d0 / axlen
                if (isgrid) then
                   idx = nint(xx(ix) * f(refden)%n(ix))
                   faux = (stm_bisect_grid(i,j,idx,ip1,ip2,ix,rcur) - rtop0 * axlen) * bohrtoa
                else
                   faux = (stm_bisect(xx,ix,rcur) - rtop0 * axlen) * bohrtoa
                endif
                !$omp critical (ffwrite)
                ff(i,j) = faux
                !$omp end critical (ffwrite)
             end do
          end do
          !$omp end parallel do
       else
          ! constant height
          write (uout,'("  Mode: constant height, ",A,"-axis coord. = ",A)') &
             laxis(ix), trim(string(rhei,'f',8,4))
          write (uout,'("  In-plane points: ",A,", ",A)') &
             string(np1), string(np2)
          !$omp parallel do private (xx,faux) schedule(dynamic)
          do i = 1, np1
             do j = 1, np2
                xx(ip1) = real(i-1,8)/np1
                xx(ip2) = real(j-1,8)/np2
                xx(ix) = rhei
                xx = cr%x2c(xx)
                faux = max(grd0(f(refden),xx),0d0)
                !$omp critical (ffwrite)
                ff(i,j) = faux
                !$omp end critical (ffwrite)
             end do
          end do
          !$omp end parallel do
       endif

       ! Write the data to a file 
       fname = trim(fileroot) // "_stm.dat"
       write (uout,'("+ Name of the STM data file: ",A)') string(fname)
       fid = fopen_write(fname)
       write (fid,'("## STM data generated by critic2.")')
       write (fid,'("## Using field number ",A," of type ",A)') string(refden), fields_typestring(f(refden)%type)
       if (iscur) then
          write (fid,'("## Mode: constant current, ",A," density.")') &
             trim(string(rcur,'e',9,3))
       else
          write (fid,'("## Mode: constant height, ",A,"-axis coord. = ",A)') &
             laxis(ix), trim(string(rhei,'f',8,4))
       endif
       write (fid,'("## In-plane points: ",A,", ",A)') string(np1), string(np2)
       write (fid,'("## All cartesian units are angstrom. fstm is a.u. for constant-height/ang for constant-current.")')
       write (fid,'("## Columns: xcryst, ycryst, xcart, ycart, fstm")')
       xmin = 1d40
       xmax = -1d40
       ymin = 1d40
       ymax = -1d40
       do i1 = 1, nx(1)
          do i = 1, np1
             do i2 = 1, nx(2)
                do j = 1, np2
                   xx(ip1) = i1-1 + real(i-1,8)/np1
                   xx(ip2) = i2-1 + real(j-1,8)/np2
                   xx(ix) = rtop0
                   xr = cr%x2c(xx) * bohrtoa
                   xmin = min(xmin,xr(ip1))
                   xmax = max(xmax,xr(ip1))
                   ymin = min(ymin,xr(ip2))
                   ymax = max(ymax,xr(ip2))
                   write (fid,'(5(A,X))') string(xx(ip1),'e',15,7,4), string(xx(ip2),'e',15,7,4), &
                      string(xr(ip1),'e',15,7,4), string(xr(ip2),'e',15,7,4), &
                      string(ff(i,j),'e',20,12,6)
                end do
             end do
             write (fid,*)
          end do
       end do
       call fclose(fid)

       ! Write the gnuplot script
       fname = trim(fileroot) // "_stm.gnu"
       write (uout,'("+ Name of the STM gnuplot script: ",A)') string(fname)
       fid = fopen_write(fname)
       write (fid,'("set encoding iso_8859_1")')
       write (fid,'("set terminal postscript eps color enhanced ''Helvetica, 20''")')
       write (fid,'("set output ''",A,"''")') trim(fileroot) // "_stm.eps"
       write (fid,*)
       ! write (fid,'("set palette model XYZ functions gray**0.35, gray**0.5, gray**0.8")')
       write (fid,'("set palette rgb 34,35,36")')
       write (fid,'("set style line 1 lt 1 lw 1 lc rgb ''#000000''")')
       write (fid,'("set style line 2 lt 1 lw 1 lc rgb ''#000000''")')
       write (fid,*)
       if (iscur) then
          write (fid,'("set title ''Constant current plot, ",A," density")') &
             trim(string(rcur,'e',9,3))
       else
          write (fid,'("set title ''Constant height plot, ",A,"-axis coord. = ",A,"''")') &
             laxis(ix), trim(string(rhei,'f',8,4))
       endif
       write (fid,*)
       write (fid,'("set xrange [",A,":",A,"]")') trim(string(xmin,'e',12,6)), trim(string(xmax,'e',12,6))
       write (fid,'("set yrange [",A,":",A,"]")') trim(string(ymin,'e',12,6)), trim(string(ymax,'e',12,6))
       write (fid,'("set xtics nomirror out")')
       write (fid,'("set ytics nomirror out")')
       write (fid,'("set xlabel ''x (\305)''")')
       write (fid,'("set ylabel ''y (\305)''")')
       write (fid,*)
       write (fid,'("unset key")')
       write (fid,'("set size ratio -1")')
       write (fid,*)
       write (fid,'("# set pm3d at b map interpolate 5,5")')
       write (fid,'("set pm3d at b map")')
       write (fid,*)
       if (iscur) then
          write (fid,'("set cbtics")')
       else
          write (fid,'("set cbtics format ''%.2f''")')
       end if
       if (iscur) then
          write (fid,'("set cblabel ''Distance to the surface (\305)''")')
       else
          write (fid,'("set cblabel ''Electron density (a.u.)''")')
       endif
       write (fid,*)
       write (fid,'("splot ''",A,"'' u 3:4:5 w pm3d notitle ls 1")') &
          trim(fileroot) // "_stm.dat"
       call fclose(fid)

       ! wrap up
       deallocate(ff)
    else
       allocate(fl(nn))
       ! Calculate the STM signal along a line
       !$omp parallel do private (xx,faux) schedule(dynamic)
       do i = 1, nn
          xx(ip1) = x0(1) + (x1(1)-x0(1)) * real(i-1,8) / (nn-1)
          xx(ip2) = x0(2) + (x1(2)-x0(2)) * real(i-1,8) / (nn-1)
          if (iscur) then
             ! constant-current 
             xx(ix) = rtop + 0.2d0 / axlen
             faux = (stm_bisect(xx,ix,rcur) - rtop0 * axlen) * bohrtoa
          else
             xx(ix) = rhei
             xx = cr%x2c(xx)
             faux = max(grd0(f(refden),xx),0d0)
          endif
          !$omp critical (flwrite)
          fl(i) = faux
          !$omp end critical (flwrite)
       end do
       !$omp end parallel do

       ! Write the line to a file 
       fname = trim(fileroot) // "_stm_line.dat"
       write (uout,'("+ Name of the STM data file (line): ",A)') string(fname)
       fid = fopen_write(fname)
       write (fid,'("## STM data generated by critic2.")')
       write (fid,'("## Using field number ",A," of type ",A)') string(refden), fields_typestring(f(refden)%type)
       if (iscur) then
          write (fid,'("## Mode: constant current, ",A," density.")') &
             trim(string(rcur,'e',9,3))
       else
          write (fid,'("## Mode: constant height, ",A,"-axis coord. = ",A)') &
             laxis(ix), trim(string(rhei,'f',8,4))
       endif
       write (fid,'("## Initial point (cryst): ",A,", ",A)') trim(string(x0(1),'e',12,5)), trim(string(x0(2),'e',12,5))
       write (fid,'("## Final point (cryst): ",A,", ",A)') trim(string(x1(1),'e',12,5)), trim(string(x1(2),'e',12,5))
       write (fid,'("## Number of  points: ",A)') trim(string(nn))
       write (fid,'("## All cartesian units are angstrom. fstm is a.u. for constant-height/ang for constant-current.")')
       write (fid,'("## Columns: xcryst, ycryst, xcart, ycart, distance, fstm")')
       do i = 1, nn
          xx(ip1) = x0(1) + (x1(1)-x0(1)) * real(i-1,8) / (nn-1)
          xx(ip2) = x0(2) + (x1(2)-x0(2)) * real(i-1,8) / (nn-1)
          xx(ix) = rtop0
          xr = cr%x2c(xx) * bohrtoa
          dist = sqrt(xr(ip2)**2+xr(ip1)**2)
          write (fid,'(6(A,X))') string(xx(ip1),'e',15,7,4), string(xx(ip2),'e',15,7,4), &
             string(xr(ip1),'e',15,7,4), string(xr(ip2),'e',15,7,4), &
             string(dist,'e',15,7,4), string(fl(i),'e',20,12,6)
       end do
       call fclose(fid)
       deallocate(fl)
    endif

    ! final blank line
    write (uout,*)

  end subroutine stm_driver

  subroutine detect_vacuum(ix,rtop)
    use crystalmod, only: cr
    use tools_io, only: faterr, ferror
    integer, intent(out) :: ix
    real*8, intent(out) :: rtop
    
    integer :: i, j, jj, idx, nempty
    integer, allocatable :: ibin(:)
    real*8 :: rcur
    logical :: iszero

    integer, parameter :: nbin = 100
    integer, parameter :: nbin_empty = 25

    allocate(ibin(nbin))

    ! Apply a 0.01-binning in each of the crystallographic directions
    do i = 3, 1, -1
       ! Bin all the atoms according to their position along that
       ! crystallographic axis
       ibin = 0
       do j = 1, cr%ncel
          idx = nint(cr%atcel(j)%x(i) / (1d0/nbin))
          idx = modulo(idx,nbin)+1
          ibin(idx) = ibin(idx) + 1
       end do

       ! detect if there is a 25%-length contiguous segment of ibin with
       ! zero population
       iszero = .false.
       nempty = 0
       do jj = 1, nbin+nbin_empty
          j = modulo(jj-1,nbin)+1
          if (ibin(j) > 0) then
             iszero = .false.
             nempty = 0
          else
             if (.not.iszero) then
                iszero = .true.
                rcur = real(j-1,8) / nbin
             end if
             nempty = nempty + 1
             if (nempty >= nbin_empty) then
                ix = i
                rtop = rcur
                if (allocated(ibin)) deallocate(ibin)
                return
             endif
          endif
       end do
    end do
    if (allocated(ibin)) deallocate(ibin)

    call ferror('detect_vacuum','Could not detect vacuum for STM calculation',faterr)

  end subroutine detect_vacuum

  !> Starting at point x in cryst coords, bracket and bisect in the
  !> ix direction until the rho0 value of the reference field is found
  !> The reference field is assumed to be locally monotonic along ix.
  !> The output z is in cartesian.
  function stm_bisect(x,ix,rho0) result(z)
    use fields, only: f, grd0
    use crystalmod, only: cr
    use global, only: refden
    use tools_io, only: faterr, ferror
    
    real*8, intent(in) :: x(3)
    integer, intent(in) :: ix
    real*8, intent(in) :: rho0
    real*8 :: z

    real*8 :: x0(3), ff, xbra(2), fbra(2), axlen, xx(3)
    integer :: isign

    real*8, parameter :: step = 0.3d0

    axlen = cr%aa(ix)

    ! determine the bracketing direction
    x0 = cr%x2c(x)
    ff = grd0(f(refden),x0)
    xbra(1) = x0(ix)
    fbra(1) = ff
    if (ff > rho0) then
       isign = 1
    else
       isign = -1
    endif

    ! bracket 
    xx = x0
    do while(abs(xx(ix)-x0(ix)) < 0.5d0*axlen)
       xx(ix) = xx(ix) + isign * step
       ff = grd0(f(refden),xx)
       if ((isign == 1 .and. ff <= rho0) .or. (isign == -1 .and. ff > rho0)) exit
    end do
    if (abs(xx(ix)-x0(ix)) >= 0.5d0*axlen) &
       call ferror('stm_bisect','Could not bracket target density',faterr)
    if (xx(3) < xbra(1)) then
       xbra(2) = xbra(1)
       fbra(2) = fbra(1)
       xbra(1) = xx(3)
       fbra(1) = ff
    else
       xbra(2) = xx(3)
       fbra(2) = ff
    endif
    fbra = fbra - rho0

    ! bisection
    do while (abs(xbra(1)-xbra(2)) > 1d-4)
       xx(3) = 0.5d0 * (xbra(1) + xbra(2))
       ff = grd0(f(refden),xx) - rho0
       if (ff * fbra(1) > 0d0) then
          xbra(1) = xx(3)
          fbra(1) = ff
       else
          xbra(2) = xx(3)
          fbra(2) = ff
       endif
    end do
    z = 0.5d0*(xbra(1)+xbra(2))

  endfunction stm_bisect

  !> Find the grid point in the ix direction whose density is closest
  !> to rho0, then do a linear interpolation to find the fractional
  !> coordinate corresponding to that point. ip1 and ip2 are the indices
  !> for the other two crystallographic axes. The ix-line is given 
  !> by grid points i1 (ip1-axis) and i2 (ip2-axis) and the search starts
  !> at iz. The output z is in cartesian. The reference field is
  !> assumed to be locally monotonic along ix. 
  function stm_bisect_grid(i1,i2,iz,ip1,ip2,ix,rho0) result(z)
    use fields, only: f, type_grid
    use crystalmod, only: cr
    use global, only: refden
    use tools_io, only: ferror, faterr
    
    integer, intent(in) :: i1, i2, iz, ip1, ip2, ix
    real*8, intent(in) :: rho0
    real*8 :: z

    integer :: isign, i, k, nx, ibra(2)
    real*8 :: axlen, f0, f1, fbra(2)
    logical :: found

    if (f(refden)%type /= type_grid) &
       call ferror('stm_bisect_grid','not a grid',faterr)
    axlen = cr%aa(ix)
    nx = f(refden)%n(ix)

    ! determine the stepping direction
    f0 = feval_grid(i1,i2,iz,ip1,ip2,ix)
    if (f0 > rho0) then
       isign = 1
    else
       isign = -1
    endif
    
    ! find the interval
    found = .false.
    k = iz
    do i = 1, nx/2
       k = modulo(k+isign-1,nx)+1
       f1 = feval_grid(i1,i2,k,ip1,ip2,ix)
       if ((isign == 1 .and. f1 <= rho0).or.(isign == -1 .and. f1 > rho0)) then
          if (isign == 1) then
             ibra(1) = k-1
             fbra(1) = feval_grid(i1,i2,modulo(k-1-1,nx)+1,ip1,ip2,ix)
             ibra(2) = k
             fbra(2) = f1
          else
             ibra(1) = k
             fbra(1) = f1
             ibra(2) = k+1
             fbra(2) = feval_grid(i1,i2,modulo(k+1-1,nx)+1,ip1,ip2,ix)
          endif
          found =.true.
          exit
       endif
    end do
    if (.not.found) call ferror('stm_bisect_grid','Could not find density bracket',faterr)
    
    ! interpolate the point in grid coordinates (1 to nx)
    z = ibra(1) + (rho0-fbra(1))/(fbra(2)-fbra(1)) * (ibra(2)-ibra(1))

    ! transform to crystallographic
    z = (z-1) / nx
    z = z - floor(z)

    ! transform to cartesian
    z = z * axlen

  contains
    function feval_grid(i,j,k,n1,n2,n3) result(ff)
      integer, intent(in) :: i, j, k, n1, n2, n3
      real*8 :: ff

      if (n1 == 1 .and. n2 == 2 .and. n3 == 3) then
         ff = f(refden)%f(i,j,k)
      elseif (n1 == 1 .and. n2 == 3 .and. n3 == 2) then
         ff = f(refden)%f(i,k,j)
      elseif (n1 == 2 .and. n2 == 1 .and. n3 == 3) then
         ff = f(refden)%f(j,i,k)
      elseif (n1 == 2 .and. n2 == 3 .and. n3 == 1) then
         ff = f(refden)%f(k,i,j)
      elseif (n1 == 3 .and. n2 == 1 .and. n3 == 2) then
         ff = f(refden)%f(j,k,i)
      elseif (n1 == 3 .and. n2 == 2 .and. n3 == 1) then
         ff = f(refden)%f(k,j,i)
      else
         call ferror('feval_grid','Unknown grid direction indices',faterr)
      endif
    endfunction feval_grid
  endfunction stm_bisect_grid

end module stm
