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
  ! private :: trick_grid_sphere
  ! private :: trick_stephens_nnm_channel
  ! private :: trick_cell_integral
  private :: trick_uspex_unpack

contains

  subroutine trick(line0)
    use tools_io, only: lgetword, equal, ferror, faterr
    character*(*), intent(in) :: line0
    character(len=:), allocatable :: word

    integer :: lp

    lp = 1
    word = lgetword(line0,lp)

    if (equal(word,'uspex_unpack')) then
       call trick_uspex_unpack(line0(lp:))
    else
       call ferror('trick','Unknown keyword: ' // trim(word),faterr,line0,syntax=.true.)
       return
    end if

    ! call trick_grid_sphere()
    ! call trick_stephens_nnm_channel(line0)
    ! call trick_cell_integral()
    ! call trick_check_shortest()
    ! call trick_test_environment()

  end subroutine trick

  ! !> Test for the calculation of energies using a sphere plus grid approach.
  ! subroutine trick_grid_sphere()
  !   use rhoplot, only: rhoplot_cube
  !   use integration, only: int_reorder_gridout
  !   use yt, only: yt_integrate, yt_weights
  !   use systemmod, only: sy
  !   use fieldmod, only: type_grid
  !   use grid3mod, only: grid3
  !   use tools_math, only: select_lebedev, gauleg
  !   use tools_io, only: string
  !   use types, only: scalar_value_noalloc, realloc
  !   use param, only: pi

  !   integer, parameter :: nleb = 770, nr = 50

  !   integer :: i, j, k, l
  !   integer :: ix, iy, iz
  !   real*8 :: rsph(cr%ncel)
  !   real*8 :: xleb(nleb), yleb(nleb), zleb(nleb), wleb(nleb)
  !   real*8 :: rp(nr), rw(nr), x(3), x0(3), unit(3)
  !   real*8 :: ff, fs, rsumf, rsums, ssumf, ssums
  !   real*8 :: d2, atp(cr%ncel), ntot
  !   type(scalar_value_noalloc) :: res, res2
  !   ! for yt
  !   integer :: nbasin, luw
  !   real*8, allocatable :: xcoord(:,:)
  !   integer, allocatable :: idg(:,:,:), icp(:)
  !   real*8, allocatable :: w(:,:,:)

  !   call sy%c%checkflags(.false.,nn0=.true.)

  !   ntot = f(2)%n(1)*f(2)%n(2)*f(2)%n(3)

  !   !! model input:
  !   ! molecule benzene.wfx 5
  !   ! load benzene.wfx
  !   ! load rho.cube
  !   ! load lap.cube
  !   ! trick

  !   ! integrate with yt
  !   write (*,*) "Doing YT on field 2"
  !   call yt_integrate(cr,f(2)%f,"",.true.,2d0,nbasin,xcoord,idg,luw)
  !   call int_reorder_gridout(cr,f(2),nbasin,xcoord,idg,.true.,2d0,luw,icp)
  !   ! call bader_integrate(cr,f(2),nbasin,xcoord,idg)

  !   ! defining the weight fields
  !   call realloc(fused,25)
  !   call realloc(f,25)
  !   allocate(w(f(2)%n(1),f(2)%n(2),f(2)%n(3)))
  !   do i = 1, cr%ncel
  !      call yt_weights(luw=luw,idb=i,w=w)
  !      call grid_from_array3(w,f(10+i))
  !      fused(10+i) = .true.
  !      f(10+i)%type = type_grid
  !   end do
  !   close(luw)
  !   deallocate(w)

  !   call rhoplot_cube("grid field 11 file bleh.cube")

  !   ! determine the radius for each atom
  !   write (*,*) "Calculating the radii"
  !   do i = 1, cr%ncel
  !      rsph(i) = cr%at(cr%atcel(i)%idx)%rnn2 * 0.99d0
  !   end do

  !   ! field 4 will contain the smooth grid
  !   f(4) = f(3)
  !   fused(4) = .true.

  !   ! smooth out the regions inside the spheres
  !   write (*,*) "Calculating the smooth grid"
  !   !$omp parallel do private(x0,x,d2)
  !   do ix = 1, f(4)%n(1)
  !      x0(1) = real(ix-1,8) / f(4)%n(1)
  !      do iy = 1, f(4)%n(2)
  !         x0(2) = real(iy-1,8) / f(4)%n(2)
  !         do iz = 1, f(4)%n(3)
  !            x0(3) = real(iz-1,8) / f(4)%n(3)
  !            do i = 1, cr%ncel
  !               x = cr%atcel(i)%x - x0
  !               call cr%shortest(x,d2)
  !               if (d2 <= rsph(i)) then
  !                  !$omp critical (smooth)
  !                  f(4)%f(ix,iy,iz) = f(4)%f(ix,iy,iz) * sin(0.5d0*pi*d2/rsph(i))**6
  !                  !$omp end critical (smooth)
  !                  exit
  !               end if
  !            end do
  !         end do
  !      end do
  !   end do
  !   !$omp end parallel do

  !   ! integrate the basins
  !   atp = 0d0
  !   call select_lebedev(nleb,xleb,yleb,zleb,wleb)
  !   do i = 1, cr%ncel
  !      write (*,*) "integrating atom ", i
  !      ! integrate the smooth field
  !      atp(i) = sum(f(10+i)%f * f(4)%f) * cr%omega / ntot

  !      ! calculate the difference between the smooth field and the actual field
  !      ! inside the basin. Contributions from all the spheres.
  !      do j = 1, cr%ncel
  !         ssumf = 0d0
  !         ssums = 0d0
  !         !$omp parallel do private(unit,rp,rw,rsumf,rsums,x,fs,ff,res,res2)
  !         do k = 1, nleb
  !            unit = (/xleb(k),yleb(k),zleb(k)/)
  !            call gauleg(0d0,rsph(j),rp,rw,nr)
  !            rsumf = 0d0
  !            rsums = 0d0
  !            do l = 1, nr
  !               x = cr%atcel(j)%r + rp(l) * unit
  !               call grd(f(refden),x,2,res0_noalloc=res)
  !               call grd(f(10+i),x,2,res0_noalloc=res2)
  !               ff = res%del2f * res2%f
  !               fs = ff * sin(0.5d0*pi*rp(l)/rsph(j))**6
  !               rsumf = rsumf + rp(l)**2 * rw(l) * ff
  !               rsums = rsums + rp(l)**2 * rw(l) * fs
  !            end do
  !            !$omp critical (acum_ang)
  !            ssumf = ssumf + rsumf * wleb(k)
  !            ssums = ssums + rsums * wleb(k)
  !            !$omp end critical (acum_ang)
  !         end do
  !         !$omp end parallel do
  !         write (*,*) "sphere ", j, ssumf, ssums
  !         atp(i) = atp(i) + ssumf - ssums
  !      end do
  !      write (*,*) "atomic laplacian: ", atp(i)
  !   end do

  !   write (*,*) "all ", sum(atp)

  ! end subroutine trick_grid_sphere

  ! subroutine trick_stephens_nnm_channel(line)
  !   use crystalmod, only: cr
  !   use fields, only: f, grd0
  !   use graphics, only: graphics_open, graphics_ball, graphics_stick, graphics_close
  !   use global, only: eval_next
  !   use tools_math, only: norm, cross
  !   use tools_io, only: faterr, ferror, uout, string
  !   use param, only: bohrtoa, pi
  !   character*(*), intent(in) :: line

  !   ! parameters for input
  !   integer, parameter :: nang = 30
  !   integer, parameter :: nline = 30
  !   real*8, parameter :: rho0_def = 0.01d0
  !   real*8, parameter :: zeronorm = 1d-10 ! zero norm
  !   real*8, parameter :: dmax = 10d0 ! maximum distance
  !   ! bracketing and bisection
  !   real*8, parameter :: bstep0 = 0.05d0
  !   real*8, parameter :: bfactor = 1.2d0
  !   real*8, parameter :: blimit = 1d-5

  !   integer :: luobj, lumtl, lp, i, j
  !   logical :: ok
  !   real*8, dimension(3) :: x0, x1, x, xp1, xp2, a1, a2, xd, xa, xb
  !   real*8 :: ang, rho, rhoi, step, dist, rhoa, rhob
  !   real*8 :: xlimit(3,nang,nline), mindist, maxdist
  !   real*8 :: rho0

  !   ! read input
  !   lp = 1
  !   ok = eval_next(x0(1),line,lp)
  !   ok = ok .and. eval_next(x0(2),line,lp)
  !   ok = ok .and. eval_next(x0(3),line,lp)
  !   ok = ok .and. eval_next(x1(1),line,lp)
  !   ok = ok .and. eval_next(x1(2),line,lp)
  !   ok = ok .and. eval_next(x1(3),line,lp)
  !   if (.not.ok) &
  !      call ferror('trick_stephens_nnm_channel','syntax: trick x0 y0 z0 x1 y1 z1',faterr)
  !   ok = eval_next(rho0,line,lp)
  !   if (.not.ok) rho0 = rho0_def

  !   ! convert to cartesian
  !   x0 = cr%x2c(x0)
  !   x1 = cr%x2c(x1)

  !   ! calculate the promolecular density profile
  !   ok = .true.
  !   write (uout,'("+ Promolecular density along the line")')
  !   write (uout,'("# Distance(ang)  rho(promol) (a.u.)")')
  !   do i = 1, nline
  !      x = x0 + real(i-1,8) / real(nline-1,8) * (x1-x0)
  !      rho = grd0(f(0),x)
  !      write (uout,'(4X,99(A,X))') string(norm(x-x0)*bohrtoa,'f',14,6), string(rho,'e',12,6)
  !      ok = ok .and. rho < rho0
  !   end do
  !   if (.not.ok) &
  !      call ferror('trick_stephens_nnm_channel','promolecular density along the line higher than threshold',faterr)

  !   ! vector in the line direction, normalized
  !   a1 = x1-x0
  !   a1 = a1 / norm(a1)

  !   write (uout,'("# Id  Distance(ang)  Avg. diam(ang)  Circle area(ang^2)")')
  !   do i = 1, nline
  !      ! position and density at this point along the line
  !      x = x0 + real(i-1,8) / real(nline-1,8) * (x1-x0)
  !      rhoi = grd0(f(0),x)

  !      ! find a vector normal to x-x0
  !      ok = .false.
  !      do j = 1, 3
  !         a2 = 0d0
  !         a2(j) = 1d0
  !         xp1 = cross(a1,a2)
  !         if (norm(xp1) > zeronorm) then
  !            ok = .true.
  !            exit
  !         end if
  !      end do
  !      if (.not.ok) &
  !         call ferror('trick_stephens_nnm_channel','could not find a vector orthogonal to the line',faterr)
  !      xp1 = xp1 / norm(xp1)

  !      ! find the other normal vector
  !      xp2 = cross(a1,xp1)
  !      if (norm(xp2) < zeronorm) &
  !         call ferror('trick_stephens_nnm_channel','could not find the second orthogonal vector',faterr)
  !      xp2 = xp2 / norm(xp2)

  !      ! run over the angles
  !      mindist = 1d30
  !      maxdist = 0d0
  !      do j = 1, nang
  !         ang = 2d0 * pi * real(j-1,8) / real(nang,8)

  !         ! direction along this angle
  !         xd = cos(ang) * xp1 + sin(ang) * xp2

  !         ! bracket
  !         step = bstep0 / bfactor
  !         xb = x
  !         rhob = rhoi
  !         do while (rhob <= rho0)
  !            step = step * bfactor
  !            xa = xb
  !            rhoa = rhob
  !            xb = x + xd * step
  !            rhob = grd0(f(0),xb)
  !            if (step > dmax) then
  !               xb = x + xd * dmax
  !               exit
  !            endif
  !         end do

  !         ! bisection
  !         dist = norm(xb - xa)
  !         do while (dist >= blimit)
  !            xd = 0.5d0 * (xb + xa)
  !            rho = grd0(f(0),xb)
  !            if (rho > rho0) then
  !               rhob = rho
  !               xb = xd
  !            else
  !               rhoa = rho
  !               xa = xd
  !            endif
  !            dist = norm(xb - xa)
  !         end do
  !         dist = norm(xd-x)
  !         maxdist = max(maxdist,dist)
  !         mindist = min(mindist,dist)
  !         xlimit(:,j,i) = xd
  !      end do

  !      ! Calculate the average radius and output
  !      dist = 0d0
  !      do j = 1, nang
  !         dist = dist + norm(xlimit(:,j,i) - x)
  !      end do
  !      dist = dist / nang
  !      write (uout,'(2X,99(A,X))') string(i), string(norm(x-x0)*bohrtoa,'f',14,6,4), &
  !         string(2d0*dist*bohrtoa,'f',14,6,4), string(pi*dist*dist*bohrtoa*bohrtoa,'f',14,6,4)
  !   end do

  !   ! write the obj file
  !   call graphics_open("obj","tube.obj",luobj,lumtl)
  !   do i = 1, nline
  !      x = x0 + real(i-1,8) / real(nline-1,8) * (x1-x0)
  !      call graphics_ball("obj",luobj,x,(/0,0,255/),0.2d0)
  !   end do
  !   do i = 1, nang
  !      call graphics_stick("obj",luobj,x0,xlimit(:,i,1),(/0,0,0/),0.1d0)
  !      call graphics_stick("obj",luobj,x1,xlimit(:,i,nline),(/0,0,0/),0.1d0)
  !      do j = 1, nline-1
  !         call graphics_stick("obj",luobj,xlimit(:,i,j),xlimit(:,i,j+1),(/0,0,0/),0.05d0)
  !      end do
  !   end do
  !   do i = 1, nline
  !      do j = 1, nang-1
  !         call graphics_stick("obj",luobj,xlimit(:,j,i),xlimit(:,j+1,i),(/0,0,0/),0.05d0)
  !      end do
  !      call graphics_stick("obj",luobj,xlimit(:,1,i),xlimit(:,nang,i),(/0,0,0/),0.05d0)
  !   end do
  !   call graphics_close("obj",luobj,lumtl)

  ! end subroutine trick_stephens_nnm_channel

  ! !> Calculate the cell integral of the reference field using
  ! !> Franchini et al.'s Becke-style mesh
  ! subroutine trick_cell_integral()
  !   use fields, only: f
  !   use crystalmod, only: cr
  !   use meshmod, only: genmesh, fillmesh
  !   use global, only: refden
  !   use types, only: molmesh
  !   use tools_io, only: uout, string
  !   use param, only: im_rho

  !   type(molmesh) :: m
  !   integer :: prop(1), id(1)

  !   write (uout,'("* Trick: cell integral")')
  !   m = genmesh(cr)
  !   write (uout,'("  mesh size      ",A)') string(m%n)
  !   id(1) = 1
  !   prop(1) = im_rho
  !   call fillmesh(m,f(refden),id,prop,.not.cr%ismolecule)
  !   write (uout,'("cell integral ",A)') string(sum(m%f(:,1) * m%w),'f',12,6)

  ! end subroutine trick_cell_integral

  ! !> Check that the shortest vector is inside the calculated WS cell.
  ! subroutine trick_check_shortest()
  !   use systemmod, only: sy
  !   use tools_io, only: uout, string, ioj_right

  !   integer :: i, j, k, l
  !   real*8 :: x(3), xs(3), xin(3), dd
  !   integer :: nvws, nnok
  !   real*8, allocatable :: xvws(:,:)
  !   logical :: inside

  !   integer, parameter :: npts = 10
  !   real*8, parameter :: eps = 1d-10

  !   write (uout,*) "* Checking the implementation of the shortest routine"
  !   write (uout,*)

  !   write (uout,'("Voronoi-relevant vectors in the WS cell (cryst. coords.): ",A)') string(nvws)
  !   nvws = sy%c%ws_nf
  !   do i = 1, nvws
  !      write (uout,'(A,":",3(X,A))') string(i), (string(sy%c%ws_ineighx(j,i)),j=1,3)
  !   end do
  !   write (uout,*)

  !   allocate(xvws(3,nvws))
  !   write (uout,'("Half-Voronoi-relevant vectors in the WS cell (Cartesian coords.): ",A)') string(nvws)
  !   do i = 1, nvws
  !      x = 0.5d0 * sy%c%ws_ineighc(:,i)
  !      xvws(:,i) = x
  !      write (uout,'(A,":",3(X,A))') string(i), (string(x(j),'f',10,5,ioj_right),j=1,3)
  !   end do
  !   write (uout,*)

  !   write (uout,'("Testing points")')
  !   nnok = 0
  !   do i = 1, npts
  !      do j = 1, npts
  !         do k = 1, npts
  !            xs = (/real(i-1,8),real(j-1,8),real(k-1,8)/) / real(npts,8)
  !            xin = xs
  !            call sy%c%shortest(xs,dd)
  !            x = sy%c%c2x(xs)

  !            inside = .true.
  !            do l = 1, nvws
  !               inside = inside .and. (dot_product(xvws(:,l),xs-xvws(:,l)) <= eps)
  !            end do
  !            if (.not.inside) then
  !               nnok = nnok + 1
  !               write (uout,'("Orig:",3(X,A)," | Tran:",3(X,A)," | ",A)') (string(xin(l),'f',10,5,ioj_right),l=1,3), &
  !                  (string(x(l),'f',10,5,ioj_right),l=1,3), string(inside)
  !            end if

  !            write (uout,'("Orig:",3(X,A)," | L:",3(X,A))') (string(xin(l),'f',10,5,ioj_right),l=1,3), &
  !               (string(x(l)-xin(l),'f',10,5,ioj_right),l=1,3)
  !         end do
  !      end do
  !   end do
  !   if (nnok == 0) then
  !      write (uout,'("All points OK.")')
  !   else
  !      write (uout,'("There were ",A," points outside the WS cell")') string(nnok)
  !   end if
  !   write(uout,*)

  !   deallocate(xvws)

  ! end subroutine trick_check_shortest

  ! subroutine trick_test_environment()
  !   use systemmod, only: sy
  !   use global, only: atomeps
  !   use tools_io, only: uout, string, tictac
  !   use param, only: icrd_cart, icrd_crys, icrd_rcrys, atmcov, bohrtoa

  !   integer :: i, j
  !   real*8 :: xx(3), x(3), dist1, dist2
  !   real*8 :: f1, fp1(3), fpp1(3,3)
  !   real*8 :: f2, fp2(3), fpp2(3,3)
  !   integer :: nid1, nid2, lvec1(3), lvec2(3), nid3, nid4
  !   integer :: nneig1(20), wat1(20), ierr, ierr1, ierr2
  !   real*8 :: dd1(20)
  !   integer :: nat, lveca(3)
  !   integer, allocatable :: nida(:), ishella(:)
  !   real*8, allocatable :: dista(:)
  !   logical :: ok1, ok2, ok3

  !   associate(env => sy%c%env, cr => sy%c)

  !     write (uout,'("* Testing the environment....")')

  !     ! ! Test the nearest_atom routine
  !     ! do i = 1, 100
  !     !    call random_number(x)
  !     !    x = x * 10d0 - 5d0
  !     !    call cr%nearest_atom(x,nid1,dist1,lvec1)
  !     !    call env%nearest_atom(x,icrd_crys,nid2,dist2,lvec2)
  !     !    write (*,*) "point ", i
  !     !    write (*,*) "x = ", x
  !     !    ! write (*,*) "nid1 = ", nid1
  !     !    ! write (*,*) "nid2 = ", nid2
  !     !    write (*,*) "nid = ", abs(nid1-nid2)
  !     !    ! write (*,*) "dist1 = ", dist1
  !     !    ! write (*,*) "dist2 = ", dist2
  !     !    write (*,*) "dist = ", abs(dist1-dist2)
  !     !    ! write (*,*) "lvec1 = ", lvec1
  !     !    ! write (*,*) "lvec2 = ", lvec2
  !     !    write (*,*) "lvecdif = ", abs(lvec1 - lvec2)
  !     ! end do

  !     ! ! Test the list_near_atoms routine
  !     ! do i = 1, 1
  !     !    call random_number(x)
  !     !    x = x * 10d0 - 5d0

  !     !    ! slow test
  !     !    write (*,*) "point ", i, x
  !     !    call cr%pointshell(x,10,nneig1,wat1,dd1)
  !     !    write (*,*) "pointshell environment"
  !     !    do j = 1, 10
  !     !       write (*,*) j, nneig1(j), wat1(j), dd1(j)
  !     !    end do

  !     !    call env%list_near_atoms(x,icrd_crys,.true.,nat,ierr,nida,dista,lveca,ishella,up2n=10)
  !     !    write (*,*) "list_near_atoms environment (n), nat = ", nat, "(",ierr,")"
  !     !    do j = 1, nat
  !     !       write (*,*) j, nida(j), dista(j), ishella(j)
  !     !    end do

  !     !    call env%list_near_atoms(x,icrd_crys,.true.,nat,ierr,nida,dista,lveca,ishella,up2n=10000)
  !     !    write (*,*) "list_near_atoms environment (n), nat = ", nat, "(",ierr,")"
  !     !    do j = 1, nat
  !     !       write (*,*) j, nida(j), dista(j), ishella(j)
  !     !    end do

  !     !    call env%list_near_atoms(x,icrd_crys,.true.,nat,ierr,nida,dista,lveca,ishella,up2sh=10)
  !     !    write (*,*) "list_near_atoms environment (sh), nat = ", nat, "(",ierr,")"
  !     !    do j = 1, nat
  !     !       write (*,*) j, nida(j), dista(j), ishella(j)
  !     !    end do

  !     !    call env%list_near_atoms(x,icrd_crys,.true.,nat,ierr,nida,dista,lveca,ishella,up2sh=10000)
  !     !    write (*,*) "list_near_atoms environment (sh), nat = ", nat, "(",ierr,")"
  !     !    do j = 1, nat
  !     !       write (*,*) j, nida(j), dista(j), ishella(j)
  !     !    end do

  !     !    call env%list_near_atoms(x,icrd_crys,.true.,nat,ierr,nida,dista,lveca,ishella,up2d=5d0)
  !     !    write (*,*) "list_near_atoms environment (dist), nat = ", nat, "(",ierr,")"
  !     !    do j = 1, nat
  !     !       write (*,*) j, nida(j), dista(j), ishella(j)
  !     !    end do

  !     !    call env%list_near_atoms(x,icrd_crys,.true.,nat,ierr,nida,dista,lveca,ishella,up2d=3000d0)
  !     !    write (*,*) "list_near_atoms environment (dist), nat = ", nat, "(",ierr,")"
  !     !    do j = 1, nat
  !     !       write (*,*) j, nida(j), dista(j), ishella(j)
  !     !    end do
  !     ! end do

  !     ! ! Test the list_near_atoms routine in a high-symmetry environment
  !     ! x = 0.5d0
  !     ! ! write (*,*) "point ", i, x
  !     ! ! call cr%pointshell(x,20,nneig1,wat1,dd1)
  !     ! ! write (*,*) "pointshell environment"
  !     ! ! do j = 1, 20
  !     ! !    write (*,*) j, nneig1(j), wat1(j), dd1(j)
  !     ! ! end do
  !     ! ! x = 100d0
  !     ! ! x = 0.6d0
  !     ! call env%list_near_atoms(x,icrd_crys,.true.,nat,ierr,nida,dista,lveca,ishella,up2d=20d0)
  !     ! write (*,*) "list_near_atoms environment, nat = ", nat, " ierr = ", ierr
  !     ! do j = 1, nat
  !     !    write (*,*) j, env%at(nida(j))%idx, nida(j), dista(j), ishella(j)
  !     ! end do

  !     ! ! Test the promolecular routine
  !     ! do i = 1, 100
  !     !    call random_number(x)
  !     !    x = x * 10d0 - 5d0
  !     !    x = cr%x2c(x)
  !     !    call cr%promolecular(x,f1,fp1,fpp1,2)
  !     !    call env%promolecular(x,icrd_cart,f2,fp2,fpp2,2)
  !     !    ! write (*,*) "point ", i
  !     !    ! write (*,*) "x = ", x
  !     !    write (*,*) "f ", abs(f1-f2)
  !     !    ! write (*,*) "fp1 ", fp1
  !     !    ! write (*,*) "fp2 ", fp2
  !     !    ! write (*,*) "fpp1 ", fpp1(1,:)
  !     !    ! write (*,*) "fpp2 ", fpp2(1,:)
  !     !    ! write (*,*) "fpp1 ", fpp1(2,:)
  !     !    ! write (*,*) "fpp2 ", fpp2(2,:)
  !     !    ! write (*,*) "fpp1 ", fpp1(3,:)
  !     !    ! write (*,*) "fpp2 ", fpp2(3,:)
  !     ! end do

  !     ! ! Test the promolecular routine timing
  !     ! call tictac("1")
  !     ! do i = 1, 100
  !     !    call random_number(x)
  !     !    x = x * 10d0 - 5d0
  !     !    x = cr%x2c(x)
  !     !    call cr%promolecular(x,f1,fp1,fpp1,2,periodic=.true.)
  !     ! end do
  !     ! call tictac("2")
  !     ! do i = 1, 100
  !     !    call random_number(x)
  !     !    x = x * 10d0 - 5d0
  !     !    x = cr%x2c(x)
  !     !    call env%promolecular(x,icrd_cart,f2,fp2,fpp2,2)
  !     ! end do
  !     ! call tictac("3")

  !     ! ! Test the promolecular routine: going away from the system
  !     ! xx = 0.5d0
  !     ! do i = 1, 100
  !     !    x = 0.5d0 + 2000d0 * real(i-1,8) / real(1000-1,8)
  !     !    call cr%promolecular(x,f1,fp1,fpp1,2,periodic=.true.)
  !     !    call env%promolecular(x,icrd_cart,f2,fp2,fpp2,2)
  !     !    write (*,*) i, abs(f1-f2)
  !     ! end do

  !     ! ! Test the identify_atom routine
  !     ! do i = 1, 1000
  !     !    ! test that it gives the same result as with the crystal version
  !     !    call random_number(x)
  !     !    call random_number(xx)
  !     !    x = cr%atcel(floor(xx(2)*cr%ncel+1))%r + x / norm2(x) * (xx(1)/1000d0)

  !     !    nid1 = cr%identify_atom(x,.true.)
  !     !    nid2 = env%identify_atom(x,icrd_cart,.true.)
  !     !    nid3 = cr%identify_atom(x,.false.)
  !     !    nid4 = env%identify_atom(x,icrd_cart,.false.)
  !     !    write (*,*) i, abs(nid1-nid2), abs(nid3-nid4), x

  !     !    ! ! test for not finding the atom
  !     !    ! call random_number(x)
  !     !    ! call random_number(xx)
  !     !    ! x = cr%atcel(floor(xx(2)*cr%ncel+1))%r + x / norm2(x) * (xx(1)/30d0)
  !     !    ! nid1 = env%identify_atom(x,icrd_cart,.true.)
  !     !    ! nid2 = env%identify_atom(x,icrd_cart,.false.)
  !     !    ! write (*,*) i, nid1, nid2
  !     ! end do

  !     ! ! find asterisms
  !     ! write (*,*) "old asterisms:"
  !     ! do i = 1, cr%ncel
  !     !    write (*,*) "atom ", i, " is bonded to ", cr%nstar(i)%ncon
  !     !    do j = 1, cr%nstar(i)%ncon
  !     !       write (*,*) "-", cr%nstar(i)%idcon(j), " (", cr%nstar(i)%lcon(:,j), ")"
  !     !    end do
  !     ! end do
  !     ! write (*,*) "new asterisms:"
  !     ! call env%find_asterisms_covalent(cr%nstar)
  !     ! do i = 1, cr%ncel
  !     !    write (*,*) "atom ", i, " is bonded to ", cr%nstar(i)%ncon
  !     !    do j = 1, cr%nstar(i)%ncon
  !     !       write (*,*) "-", cr%nstar(i)%idcon(j), " (", cr%nstar(i)%lcon(:,j), ")"
  !     !    end do
  !     ! end do

  !     ! ! find asterisms, timing info
  !     ! call env%report()
  !     ! call tictac("1")
  !     ! call env%find_asterisms_covalent(cr%nstar)
  !     ! call tictac("2")

  !     do i = 0, 100
  !        ! module subroutine nearest_atom_short(e,x0,icrd,distmax,eid,dist,ierr,cidx0,idx0,nozero)
  !        xx = cr%x2c(cr%atcel(1)%x + (cr%atcel(2)%x-cr%atcel(1)%x) * real(i,8) / 100d0 + (/4,-8,12/))
  !        write (*,*) "input: ", xx
  !        call env%nearest_atom_long(xx,icrd_cart,env%dmax0-1d-10,nid1,lvec1,dist1,ierr1)
  !        call env%nearest_atom_short(xx,icrd_cart,0.5d0*env%boxsize-1d-10,nid2,lvec2,dist2,ierr2)
  !        if (ierr1 == 0 .and. ierr2 == 0) then
  !           write (*,*) i/100d0, ierr1, ierr2, nid1, nid2, lvec1, lvec2, dist1, dist2, norm2(cr%x2c(cr%atcel(nid1)%x + lvec1) - xx)
  !        else if (ierr1 == 0) then
  !           write (*,*) i/100d0, ierr1, nid1, lvec1, dist1, norm2(cr%x2c(cr%atcel(nid1)%x + lvec1) - xx)
  !        else
  !           write (*,*) i/100d0, ierr1, ierr2
  !        end if
  !     end do

  !   end associate

  ! end subroutine trick_test_environment

  ! Unpack uspex structures. Syntax:
  !   TRICK USPEX_UNPACK individuals.s file.POSCAR [template.xyz] [MAXDE maxde.r] [NONEG]
  subroutine trick_uspex_unpack(line0)
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed, read_seeds_from_file
    use global, only: fileroot, rborder_def, eval_next
    use tools_io, only: lgetword, getword, ferror, faterr, getline_raw, fopen_read,&
       fclose, zatguess, isinteger, string, uout, equal
    use tools_math, only: crosscorr_triangle, umeyama_graph_matching, rmsd_walker
    use tools, only: mergesort, qcksort
    use types, only: realloc
    use param, only: bohrtoa, isformat_xyz
    character*(*), intent(in) :: line0

    character*1024 :: bleh
    character(len=:), allocatable :: fileind, fileposcar, filexyz, fileout, line, errmsg, word
    character(len=:), allocatable :: linespc, linensp, msg
    integer :: lp, lu, idxb, idum, nat, is
    integer :: ii, jj, i, j, k, iz, nst
    integer, allocatable :: idst(:), iord(:), nsp(:)
    real*8, allocatable :: vst(:), hst(:)
    logical, allocatable :: active(:)
    type(crystalseed), allocatable :: seed(:)
    real*8 :: rdum, rprim(3,3), adv, adh, minh
    logical :: ok, firstpass
    character*1 :: let
    type(crystal) :: ci, cj, templt
    real*8, allocatable :: t(:), ihi(:), ihj(:), th2p(:), ip(:)
    integer, allocatable :: hvecp(:,:)
    real*8 :: tini, tend, nor, xnormi, xnormj, diffij
    integer :: ndigit
    integer :: ns, n
    type(crystal), allocatable :: mol(:)
    type(crystalseed) :: molseed, xseed
    logical :: usetemplate, noneg
    integer, allocatable :: isperm(:,:), cidxorig(:,:), saveperm(:)
    real*8, allocatable :: dref(:,:), ddg(:,:), ddh(:,:)
    real*8, allocatable :: x1(:,:), x2(:,:), mrot(:,:,:)
    real*8 :: xcm1(3), xcm2(3), xnew(3), maxde

    real*8, parameter :: vdiff_thr = 0.1d0 ! volume difference threshold, ang^3
    real*8, parameter :: ediff_thr = 0.01d0 ! energy difference threshold, eV
    real*8, parameter :: powdiff_thr = 0.07d0 ! powdiff threshold

    real*8, parameter :: th2ini = 5d0
    integer, parameter :: npts = 10001
    real*8, parameter :: lambda0 = 1.5406d0
    real*8, parameter :: fpol0 = 0d0
    real*8, parameter :: sigma = 0.05d0
    real*8, parameter :: xend = 50d0
    real*8, parameter :: h = (xend-th2ini) / real(npts-1,8)

    ! read the file names
    lp = 1
    fileind = getword(line0,lp)
    fileposcar = getword(line0,lp)
    filexyz = getword(line0,lp)
    if ((len_trim(fileind) == 0) .or. (len_trim(fileposcar) == 0)) &
       call ferror('trick_uspex_unpack','Error: individuals or POSCAR file not found',faterr)

    ! read the options
    noneg = .false.
    maxde = 1d40
    do while (.true.)
       word = lgetword(line0,lp)
       if (equal(word,'maxde')) then
          ok = eval_next(maxde,line0,lp)
          if (.not.ok) then
             call ferror('trick_uspex_unpack','Invalid MAXDE keyword',faterr,line0,syntax=.true.)
             return
          end if
       elseif (equal(word,'noneg')) then
          noneg = .true.
       elseif (len_trim(word) > 0) then
          call ferror('trick_uspex_unpack','Unknown extra keyword',faterr,line0,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! read the template
    usetemplate = (len_trim(filexyz) /= 0)
    if (usetemplate) then
       call molseed%read_mol(filexyz,isformat_xyz,rborder_def,.false.,errmsg)
       if (len_trim(errmsg) > 0) then
          errmsg = "error reading xyz file"
          goto 999
       end if
       call templt%struct_new(molseed,.true.)
       allocate(saveperm(templt%nneq))
    end if

    ! prepare some space
    nst = 0
    allocate(idst(100),vst(100),hst(100))

    ! read the individuals file
    lu = fopen_read(fileind)
    do while (getline_raw(lu,line))
       idxb = index(line,"]")
       if (idxb > 0) then
          ! add a new structure
          nst = nst + 1
          if (nst > size(idst,1)) then
             call realloc(idst,2*nst)
             call realloc(vst,2*nst)
             call realloc(hst,2*nst)
          end if

          ! read the id
          read (line,*) idum, idst(nst)

          ! read the volume and the enthalpy
          line = line(idxb+1:)
          read (line,*) hst(nst), vst(nst)
       end if
    end do
    call fclose(lu)
    call realloc(idst,nst)
    call realloc(vst,nst)
    call realloc(hst,nst)

    ! calculate the number of digits for output
    ndigit = ceiling(log10(nst+0.1d0))

    ! sort first by volume, then by enthalpy
    allocate(iord(nst))
    do i = 1, nst
       iord(i) = i
    end do
    call mergesort(hst,iord,1,nst)
    call mergesort(vst,iord,1,nst)

    ! Read all the seeds; USPEX format is expected, so this is not a
    ! full-fledged POSCAR reader.
    errmsg = "error reading POSCAR file"
    allocate(seed(nst))
    lu = fopen_read(fileposcar)
    do i = 1, nst
       read (lu,*,err=999) bleh
       read (lu,*,err=999) rdum
       if (abs(rdum - 1d0) > 1d-14) then
          errmsg = "scale in POSCAR file different from 1"
          goto 999
       end if

       ! read the cell vectors
       do j = 1, 3
          read (lu,*) rprim(1:3,j)
       end do
       rprim = rprim / bohrtoa
       seed(i)%m_x2c = rprim
       seed(i)%useabr = 2

       ! read the atomic species
       ok = getline_raw(lu,line,.false.)
       if (i == 1) then
          if (.not.ok) goto 999
          lp = 1
          seed(i)%nspc = 0
          allocate(seed(i)%spc(2))
          do while (.true.)
             word = getword(line,lp)
             if (len_trim(word) == 0) exit
             iz = zatguess(word)
             if (iz <= 0) then
                errmsg = "unknown atomic symbol in POSCAR"
                goto 999
             end if

             seed(i)%nspc = seed(i)%nspc + 1
             if (seed(i)%nspc > size(seed(i)%spc,1)) &
                call realloc(seed(i)%spc,2*seed(i)%nspc)
             seed(i)%spc(seed(i)%nspc)%name = word
             seed(i)%spc(seed(i)%nspc)%z = iz
          end do
          call realloc(seed(i)%spc,seed(i)%nspc)
          linespc = line
       elseif (line /= linespc) then
          errmsg = "unexpected change in atomic species line in POSCAR"
          goto 999
       else
          seed(i)%nspc = seed(1)%nspc
          seed(i)%spc = seed(1)%spc
       end if

       ! read the number of atoms
       ok = getline_raw(lu,line,.true.)
       if (.not.ok) goto 999
       if (i == 1) then
          ! read the number of atoms for each species
          allocate(nsp(seed(i)%nspc))
          lp = 1
          do j = 1, seed(i)%nspc
             ok = isinteger(nsp(j),line,lp)
             if (.not.ok) then
                errmsg = "error reading number of atomic species in POSCAR"
                goto 999
             end if
          end do

          ! assign nat and species mapping
          seed(i)%nat = 0
          allocate(seed(i)%is(sum(nsp)))
          do j = 1, seed(i)%nspc
             do k = 1, nsp(j)
                seed(i)%nat = seed(i)%nat + 1
                seed(i)%is(seed(i)%nat) = j
             end do
          end do

          ! deallocate and save the line
          deallocate(nsp)
          linensp = line
       elseif (line /= linensp) then
          errmsg = "unexpected change in number of atomic species line in POSCAR"
          goto 999
       else
          seed(i)%nat = seed(1)%nat
          seed(i)%is = seed(1)%is
       end if

       ! check the "direct" keyword
       read(lu,*,err=999) let
       if (let /= 'd' .and. let /= 'D') then
          errmsg = "Direct keyword expected but not found in POSCAR"
          goto 999
       end if

       ! read atomic coordinates
       allocate(seed(i)%x(3,seed(i)%nat))
       do j = 1, seed(i)%nat
          read(lu,*,err=999) seed(i)%x(:,j)
       enddo

       ! symmetry
       seed(i)%havesym = 0
       seed(i)%findsym = -1
       seed(i)%checkrepeats = .false.

       ! rest of the seed(i) information
       seed(i)%isused = .true.
       seed(i)%ismolecule = .false.
       seed(i)%cubic = .false.
       seed(i)%border = 0d0
       seed(i)%havex0 = .false.
       seed(i)%molx0 = 0d0
       seed(i)%file = fileposcar
       seed(i)%name = ""
    end do
    call fclose(lu)

    ! Prepare for the list of structures
    write (uout,'("* List of structures")')
    write (uout,'("# Volume difference threshold (ang) = ",A)') string(vdiff_thr,'e',decimal=5)
    write (uout,'("# Enthalpy difference threshold (eV) = ",A)') string(ediff_thr,'e',decimal=5)
    write (uout,'("# POWDIFF threshold = ",A)') string(powdiff_thr,'e',decimal=5)
    write (uout,'("# A structure is pruned if the three conditions hold concurrently.")')
    write (uout,'("id    fate     file or reason for pruning")')
    allocate(active(nst))
    active = .true.

    ! pre-pruning
    minh = minval(hst,.not.noneg.or.(hst > 0))
    do ii = 1, nst
       i = iord(ii)
       if (noneg .and. hst(i) < 0d0) then
          active(i) = .false.
          write (uout,'(A," pruned   because: NONEG and energy is ",A)') &
             string(i,ndigit), string(hst(i),'f',decimal=3)
       else
          adh = abs(hst(i)-minh)
          if (adh > maxde) then
             active(i) = .false.
             write (uout,'(A," pruned   because: MAXDE and energy is ",A,", ",A," above the minimum (",A,")")') &
                string(i,ndigit), string(hst(i),'f',decimal=3),&
                string(adh,'f',decimal=3), string(minh,'f',decimal=3)
          end if
       end if
    end do

    ! prune
    firstpass = .true.
    main: do ii = 1, nst
       i = iord(ii)
       if (.not.active(i)) cycle
       call ci%struct_new(seed(i),.true.)

       ! powder for structure i
       call ci%powder(th2ini,xend,.false.,npts,lambda0,fpol0,sigma,t,ihi,th2p,ip,hvecp)
       tini = ihi(1)*ihi(1)
       tend = ihi(npts)*ihi(npts)
       nor = (2d0 * sum(ihi(2:npts-1)*ihi(2:npts-1)) + tini + tend) * (xend - th2ini) / 2d0 / real(npts-1,8)
       ihi = ihi / sqrt(nor)
       xnormi = sqrt(abs(crosscorr_triangle(h,ihi,ihi,1d0)))

       ! begin applying the template
       if (usetemplate) then
          ! check the fragments
          nat = templt%nneq
          if (any(ci%mol(:)%nat /= nat)) then
             write (uout,'(A," pruned   because: molecules blew up")') string(i,ndigit)
             active(i) = .false.
             cycle main
          end if

          ! read the fragments from ci into structures
          ns = ci%nmol
          if (allocated(mol)) deallocate(mol)
          allocate(mol(ns))
          do j = 1, ns
             call molseed%from_fragment(ci%mol(j),.true.)
             molseed%border = 0d0
             call mol(j)%struct_new(molseed,.true.)
          end do

          ! allocate the use and permutation arrays
          if (allocated(isperm)) deallocate(isperm)
          allocate(isperm(templt%nneq,ns))

          ! make the mapping between structure+atom and the original cidx
          ! sort because the from_fragment routine also sorts
          if (allocated(cidxorig)) deallocate(cidxorig)
          allocate(cidxorig(nat,ns))
          do is = 1, ns
             do j = 1, nat
                cidxorig(j,is) = ci%mol(is)%at(j)%cidx
             end do
             call qcksort(cidxorig(:,is))
          end do

          ! allocate space and get reference distance matrix
          allocate(dref(nat,nat),ddg(nat,nat),ddh(nat,nat))
          call templt%distmatrix(dref,conn=.true.)

          ! calculate the permutations
          do is = 1, ns
             ddg = dref
             call mol(is)%distmatrix(ddh,conn=.true.)
             call umeyama_graph_matching(nat,ddg,ddh,isperm(:,is))
             if (any(templt%spc(templt%at(1:nat)%is)%z /= mol(is)%spc(mol(is)%at(isperm(1:nat,is))%is)%z)) then
                write (uout,'(A," pruned   because: atomic numbers in permutation did not match")') &
                   string(i,ndigit)
                active(i) = .false.
                cycle main
             end if
          end do
          deallocate(dref,ddg,ddh)

          ! Get the rmsd from walker. If moveatoms, also save the rotation.
          if (allocated(mrot)) deallocate(mrot)
          allocate(x1(3,nat),x2(3,nat),mrot(3,3,ns))
          do j = 1, nat
             x1(:,j) = templt%at(j)%r
          end do
          do is = 1, ns
             do j = 1, nat
                x2(:,j) = mol(is)%at(isperm(j,is))%r
             end do
             rdum = rmsd_walker(x1,x2,mrot(:,:,is))
          end do
          deallocate(x1,x2,mol)

          ! write the final structure
          xcm1 = templt%mol(1)%cmass(.false.)
          call ci%makeseed(xseed,.false.)
          n = 0
          do is = 1, ns
             xcm2 = ci%mol(is)%cmass(.false.)
             do j = 1, nat
                n = n + 1

                xnew = templt%mol(1)%at(j)%r - xcm1
                xnew = matmul(mrot(:,:,is),xnew)
                xnew = xnew + xcm2
                xnew = ci%c2x(xnew)

                xseed%x(:,n) = xnew
                xseed%is(n) = ci%atcel(cidxorig(isperm(j,is),is))%is
             end do
          end do
          call ci%struct_new(xseed,.true.)

          ! check the permutations are always the same
          if (firstpass) then
             saveperm = isperm(:,1)
             firstpass = .false.
          end if
          do is = 1, ns
             if (any(isperm(:,is) /= saveperm)) then
                write (uout,'(A," pruned   because: permutation does not match previous structures")') &
                   string(i,ndigit)
                active(i) = .false.
                cycle main
             end if
          end do
          ! do is = 1, ns
          !    msg = ""
          !    do j = 1, nat
          !       msg = msg // " " // string(isperm(j,is))
          !    end do
          !    write (uout,'(X,A," success (perm=",A,")")') string(is,2), string(msg)
          ! end do

          deallocate(isperm,cidxorig,mrot)
       end if

       ! remove duplicates from the list of structures
       do jj = ii+1, nst
          j = iord(jj)
          if (.not.active(j)) cycle

          adv = abs(vst(i)-vst(j))
          adh = abs(hst(i)-hst(j))
          if (adv > vdiff_thr) then
             ! sorted by volume, so it is pointless to continue
             exit
          else if (adh < ediff_thr) then
             ! make structure j
             call cj%struct_new(seed(j),.true.)

             ! powder for structure j
             call cj%powder(th2ini,xend,.false.,npts,lambda0,fpol0,sigma,t,ihj,th2p,ip,hvecp)
             tini = ihj(1)*ihj(1)
             tend = ihj(npts)*ihj(npts)
             nor = (2d0 * sum(ihj(2:npts-1)*ihj(2:npts-1)) + tini + tend) * (xend - th2ini) / 2d0 / real(npts-1,8)
             ihj = ihj / sqrt(nor)
             xnormj = sqrt(abs(crosscorr_triangle(h,ihj,ihj,1d0)))

             ! compare structures i and j
             diffij = max(1d0 - crosscorr_triangle(h,ihi,ihj,1d0) / xnormi / xnormj,0d0)
             if (diffij < powdiff_thr) then
                write (uout,'(A," pruned   because: same as ",A," (deltaV=",A,",deltaH=",A,",pow=",A,")")') &
                   string(j,ndigit), string(i), string(adv,'e',decimal=4), string(adh,'e',decimal=4),&
                   string(diffij,'e',decimal=4)
                active(j) = .false.
             end if
          end if
       end do

       ! this structure is active, so write it
       fileout = fileroot // "-" // string(i,length=ndigit,pad0=.true.) // ".res"
       ci%file = "EA" // string(i) // " V= " // string(vst(i),'f',decimal=3) // " H= " // string(hst(i),'f',decimal=3)
       call ci%wholemols()
       call ci%write_res(fileout,-1)
       write (uout,'(A," written  ",A," (V=",A,",H=",A,")")') string(i,ndigit), string(fileout), &
          string(vst(i),'f',decimal=3), string(hst(i),'f',decimal=3)
    end do main
    write (uout,'("+ Structures (pruned/written/total): ",A,"/",A,"/",A)') string(nst-count(active)), &
       string(count(active)), string(nst)
    write (uout,*)

    return
999 continue

    call ferror('trick_uspex_unpack',errmsg,faterr,line,syntax=.true.)
    return

  end subroutine trick_uspex_unpack

end module tricks
