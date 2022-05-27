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
  private :: trick_reduce
  private :: trick_makecif_ccdc
  private :: trick_bfgs
  private :: trick_compare_deformed

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
    else if (equal(word,'reduce')) then
       call trick_reduce(line0(lp:))
    else if (equal(word,'makecif')) then
       call trick_makecif_ccdc(line0(lp:))
    else if (equal(word,'bfgs')) then
       call trick_bfgs(line0(lp:))
    else if (equal(word,'compare')) then
       call trick_compare_deformed(line0(lp:))
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
  !   TRICK USPEX_UNPACK individuals.s file.POSCAR
  !         MOLECULES mol1.xyz mol2.xyz ...
  !         PATTERN i1.i i2.i ...
  !         [MOVEATOMS] [MAXDE maxde.r] [NONEG] [OUTPUTPOSCAR]
  subroutine trick_uspex_unpack(line0)
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed, read_seeds_from_file
    use global, only: fileroot, rborder_def, eval_next
    use tools_io, only: lgetword, getword, ferror, faterr, getline_raw, fopen_read,&
       fclose, zatguess, isinteger, string, uout, equal, lower
    use tools_math, only: crosscorr_triangle, rmsd_walker
    use tools, only: mergesort, qcksort
    use types, only: realloc
    use param, only: bohrtoa, isformat_xyz
    character*(*), intent(in) :: line0

    character*1024 :: sdum
    character(len=:), allocatable :: fileind, fileposcar, fileout, line
    character(len=:), allocatable :: linespc, linensp, errmsg, word, lword
    integer :: lp, lu, idxb, idum, nat, is
    integer :: ii, jj, i, j, k, iz, nst, ntemplate, npattern, pp
    integer, allocatable :: idst(:), iord(:), nsp(:), pattern(:)
    real*8, allocatable :: vst(:), hst(:)
    logical, allocatable :: active(:)
    type(crystalseed), allocatable :: seed(:)
    real*8 :: rdum, rprim(3,3), adv, adh, minh
    logical :: ok, firstpass, outputposcar
    character*1 :: let
    type(crystal) :: ci, cj
    real*8, allocatable :: t(:), ihi(:), ihj(:), th2p(:), ip(:)
    integer, allocatable :: hvecp(:,:)
    real*8 :: tini, tend, nor, xnormi, xnormj, diffij
    integer :: ndigit
    integer :: ns, n
    type(crystal), allocatable :: mol(:)
    type(crystalseed) :: molseed, xseed
    type(crystal), allocatable :: xtemplate(:), xaux(:)
    logical :: noneg, inmolecules, inpattern, moveatoms
    integer, allocatable :: iperm(:)
    real*8, allocatable :: x1(:,:), x2(:,:)
    real*8 :: xcm1(3), xcm2(3), xnew(3), maxde, mrot(3,3)
    integer :: ibins(4), nbins(4)

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
    if ((len_trim(fileind) == 0) .or. (len_trim(fileposcar) == 0)) &
       call ferror('trick_uspex_unpack','Error: individuals or POSCAR file not found',faterr)

    ! read the options
    ntemplate = 0
    allocate(xtemplate(1))
    npattern = 0
    allocate(pattern(10))
    inmolecules = .false.
    inpattern = .false.
    noneg = .false.
    maxde = 1d40
    moveatoms = .false.
    outputposcar = .false.
    do while (.true.)
       word = getword(line0,lp)
       lword = lower(word)
       if (equal(lword,'maxde')) then
          ok = eval_next(maxde,line0,lp)
          if (.not.ok) then
             call ferror('trick_uspex_unpack','Invalid MAXDE keyword',faterr,line0,syntax=.true.)
             return
          end if
          inmolecules = .false.
          inpattern = .false.
       elseif (equal(lword,'moveatoms')) then
          moveatoms = .true.
       elseif (equal(lword,'noneg')) then
          noneg = .true.
          inmolecules = .false.
          inpattern = .false.
       elseif (equal(lword,'molecules')) then
          inmolecules = .true.
          inpattern = .false.
       elseif (equal(lword,'pattern')) then
          inpattern = .true.
          inmolecules = .false.
       elseif (equal(lword,'outputposcar')) then
          outputposcar = .true.
       elseif (len_trim(word) > 0) then
          if (inmolecules) then
             ! read a molecule and prepare the structure in xtemplate
             ntemplate = ntemplate + 1
             if (ntemplate > size(xtemplate,1)) then
                if (allocated(xaux)) deallocate(xaux)
                allocate(xaux(2*ntemplate))
                xaux(1:size(xtemplate,1)) = xtemplate
                call move_alloc(xaux,xtemplate)
             end if
             call molseed%read_mol(word,isformat_xyz,rborder_def,.false.,errmsg)
             if (len_trim(errmsg) > 0) then
                errmsg = "error reading xyz file"
                goto 999
             end if
             call xtemplate(ntemplate)%struct_new(molseed,.true.)
          elseif (inpattern) then
             npattern = npattern + 1
             if (npattern > size(pattern,1)) call realloc(pattern,2*npattern)
             ok = isinteger(pattern(npattern),word)
             if (.not.ok) then
                errmsg = "error reading pattern"
                goto 999
             end if
          else
             call ferror('trick_uspex_unpack','Unknown extra keyword',faterr,line0,syntax=.true.)
             return
          end if
       else
          exit
       end if
    end do
    call realloc(pattern,npattern)

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

    ! prepare the permutation for reordering the POSCAR files
    if (npattern == 0 .or. ntemplate == 0) then
       errmsg = "pattern or molecular templates not found"
       goto 999
    end if
    nat = 0
    nbins = 0
    do i = 1, npattern
       pp = pattern(i)
       nbins(1) = nbins(1) + count(xtemplate(pp)%spc(xtemplate(pp)%atcel(1:xtemplate(pp)%ncel)%is)%z == 1)
       nbins(2) = nbins(2) + count(xtemplate(pp)%spc(xtemplate(pp)%atcel(1:xtemplate(pp)%ncel)%is)%z == 6)
       nbins(3) = nbins(3) + count(xtemplate(pp)%spc(xtemplate(pp)%atcel(1:xtemplate(pp)%ncel)%is)%z == 8)
       nbins(4) = nbins(4) + count(xtemplate(pp)%spc(xtemplate(pp)%atcel(1:xtemplate(pp)%ncel)%is)%z == 7)
       nat = nat + xtemplate(pp)%ncel
    end do
    if (sum(nbins) /= nat) then
       errmsg = "inconsistent sum of atoms, maybe atoms other than H, C, O, N?"
       goto 999
    end if
    allocate(iperm(nat))
    ibins = 0
    nat = 0
    do i = 1, npattern
       pp = pattern(i)
       do j = 1, xtemplate(pp)%ncel
          nat = nat + 1
          if (xtemplate(pp)%spc(xtemplate(pp)%atcel(j)%is)%z == 1) then
             ibins(1) = ibins(1) + 1
             iperm(ibins(1)) = nat
          else if (xtemplate(pp)%spc(xtemplate(pp)%atcel(j)%is)%z == 6) then
             ibins(2) = ibins(2) + 1
             iperm(nbins(1)+ibins(2)) = nat
          else if (xtemplate(pp)%spc(xtemplate(pp)%atcel(j)%is)%z == 8) then
             ibins(3) = ibins(3) + 1
             iperm(nbins(1)+nbins(2)+ibins(3)) = nat
          else if (xtemplate(pp)%spc(xtemplate(pp)%atcel(j)%is)%z == 7) then
             ibins(4) = ibins(4) + 1
             iperm(nbins(1)+nbins(2)+nbins(3)+ibins(4)) = nat
          else
          end if
       end do
    end do

    ! Read all the seeds; USPEX format is expected, so this is not a
    ! full-fledged POSCAR reader. Reorder the POSCAR atoms as we go.
    errmsg = "error reading POSCAR file"
    allocate(seed(nst))
    lu = fopen_read(fileposcar)
    do i = 1, nst
       read (lu,*,err=999) sdum
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
          if (any(nsp /= nbins)) then
             errmsg = "inconsistent atom types in POSCAR and molecular template pattern"
             goto 999
          end if

          ! assign nat and species mapping
          seed(i)%nat = 0
          allocate(seed(i)%is(sum(nsp)))
          do j = 1, seed(i)%nspc
             do k = 1, nsp(j)
                seed(i)%nat = seed(i)%nat + 1
                seed(i)%is(iperm(seed(i)%nat)) = j
             end do
          end do
          if (seed(i)%nat /= size(iperm,1)) then
             errmsg = "number of atoms in POSCAR inconsistent with template pattern"
             goto 999
          end if

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
          read(lu,*,err=999) seed(i)%x(:,iperm(j))
       enddo

       ! symmetry
       seed(i)%havesym = 0
       seed(i)%findsym = 0
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

       !! begin applying the template
       ! read the fragments from ci into structures
       ns = ci%nmol
       if (ns /= npattern) then
          write (uout,'(A," pruned   because: incorrect number of fragments")') string(i,ndigit)
          active(i) = .false.
          cycle main
       end if
       if (allocated(mol)) deallocate(mol)
       allocate(mol(ns))
       nat = 0
       do j = 1, ns
          do k = 1, ci%mol(j)%nat
             nat = nat + 1
             if (nat /= ci%mol(j)%at(k)%cidx) then
                write (uout,'(A," pruned   because: incorrect permutation (molecules blew up?)")') string(i,ndigit)
                active(i) = .false.
                cycle main
             end if
          end do
       end do

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

       ! move the atoms to conform with the template, if necessary
       if (moveatoms) then
          call ci%makeseed(xseed,.false.)
          xseed%findsym = 0
          n = 0
          do is = 1, ns
             pp = pattern(is)
             nat = xtemplate(pp)%ncel
             allocate(x1(3,nat),x2(3,nat))

             xcm1 = xtemplate(pp)%mol(1)%cmass(.false.)
             xcm2 = ci%mol(is)%cmass(.false.)
             do j = 1, nat
                x1(:,j) = xtemplate(pp)%at(j)%r - xcm1
                x2(:,j) = ci%mol(is)%at(j)%r - xcm2
             end do
             rdum = rmsd_walker(x1,x2,mrot)
             do j = 1, nat
                n = n + 1
                xnew = xtemplate(pp)%at(j)%r - xcm1
                xnew = matmul(mrot,xnew)
                xnew = xnew + xcm2
                xnew = ci%c2x(xnew)
                xseed%x(:,n) = xnew
             end do
             deallocate(x1,x2)
          end do
          call ci%struct_new(xseed,.true.)
       end if

       ! run wholemols
       call ci%wholemols()

       ! this structure is active, so write it
       if (outputposcar) then
          fileout = fileroot // "-pruned.POSCAR"
          ci%file = "EA" // string(i) // " V= " // string(vst(i),'f',decimal=3) // " H= " // string(hst(i),'f',decimal=3)
          call ci%write_vasp(fileout,.false.,.not.firstpass)
          firstpass = .false.
       else
          fileout = fileroot // "-" // string(i,length=ndigit,pad0=.true.) // ".res"
          ci%file = "EA" // string(i) // " V= " // string(vst(i),'f',decimal=3) // " H= " // string(hst(i),'f',decimal=3)
          call ci%write_res(fileout,-1)
       end if

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

  ! Reduce a list of structures by pruning the identical crystals. Syntax:
  !   TRICK REDUCE list.s
  ! The list.s file contains rows of the form: structure.s energy.r
  ! where structure.s is the file name containing the structure and
  ! energy.r is the energy in kcal/mol per molecule.
  subroutine trick_reduce(line0)
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use tools_io, only: getword, ferror, faterr, getline, fopen_read,&
       fclose, string, uout, isreal, ioj_center
    use tools_math, only: crosscorr_triangle
    use tools, only: mergesort
    use types, only: realloc
    use param, only: bohrtoa
    character*(*), intent(in) :: line0

    character(len=:), allocatable :: line, errmsg, filelist, file
    integer :: lu, lp, i, ii, j, jj
    integer :: nst
    integer, allocatable :: iord(:)
    real*8, allocatable :: vst(:), hst(:)
    type(crystalseed) :: seed
    type(crystal), allocatable :: st(:), staux(:)
    logical :: ok
    logical, allocatable :: active(:)

    real*8 :: adv, adh
    real*8, allocatable :: t(:), ihi(:), ihj(:), th2p(:), ip(:)
    integer, allocatable :: hvecp(:,:)
    real*8 :: tini, tend, nor, xnormi, xnormj, diffij
    integer :: ndigit, maxlen

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

    ! initialize
    errmsg = ""

    ! read the file name
    filelist = line0
    if ((len_trim(filelist) == 0)) &
       call ferror('trick_reduce','Error: list file not found',faterr)

    ! read the list file
    nst = 0
    maxlen = 0
    allocate(vst(10),hst(10),st(10))
    lu = fopen_read(filelist)
    do while (getline(lu,line))
       ! skip blank lines
       if (len_trim(line) == 0) cycle

       ! read this line's contents
       nst = nst + 1
       if (nst > size(vst,1)) then
          call realloc(vst,2*nst)
          call realloc(hst,2*nst)
          allocate(staux(2*nst))
          staux(1:size(st,1)) = st
          call move_alloc(staux,st)
       end if
       lp = 1
       file = getword(line,lp)
       ok = isreal(hst(nst),line,lp)
       if (.not.ok .or. len_trim(file) == 0) then
          errmsg = "Error reading list file"
          goto 999
       end if

       ! read and build the crystal structure
       call seed%read_any_file(file,0,errmsg)
       seed%havesym = 0
       seed%findsym = 0
       if (len_trim(errmsg) > 0) goto 999
       call st(nst)%struct_new(seed,.true.)
       vst(nst) = st(nst)%omega/real(st(nst)%nmol,8)
       maxlen = max(maxlen,len_trim(st(nst)%file))
    end do
    ! final realloc
    call fclose(lu)
    call realloc(vst,nst)
    call realloc(hst,nst)
    allocate(staux(nst))
    staux = st(1:nst)
    call move_alloc(staux,st)

    ! calculate the number of digits for output
    ndigit = ceiling(log10(nst+0.1d0))

    ! sort first by volume, then by enthalpy
    allocate(iord(nst))
    do i = 1, nst
       iord(i) = i
    end do
    call mergesort(hst,iord,1,nst)
    call mergesort(vst,iord,1,nst)

    ! header for the list of structures
    write (uout,'("* List of structures")')
    write (uout,'("# Volume difference threshold (ang) = ",A)') string(vdiff_thr,'e',decimal=5)
    write (uout,'("# Enthalpy difference threshold (eV) = ",A)') string(ediff_thr,'e',decimal=5)
    write (uout,'("# POWDIFF threshold = ",A)') string(powdiff_thr,'e',decimal=5)
    write (uout,'("# A structure is pruned if the three conditions hold concurrently.")')
    write (uout,'("id    fate     file or reason for pruning")')
    write (uout,'(A,"  decision    reason")') string("name",maxlen,ioj_center)

    ! prune
    allocate(active(nst))
    active = .true.
    main: do ii = 1, nst
       i = iord(ii)
       if (.not.active(i)) cycle

       ! powder for structure i
       call st(i)%powder(th2ini,xend,.false.,npts,lambda0,fpol0,sigma,t,ihi,th2p,ip,hvecp)
       tini = ihi(1)*ihi(1)
       tend = ihi(npts)*ihi(npts)
       nor = (2d0 * sum(ihi(2:npts-1)*ihi(2:npts-1)) + tini + tend) * (xend - th2ini) / 2d0 / real(npts-1,8)
       ihi = ihi / sqrt(nor)
       xnormi = sqrt(abs(crosscorr_triangle(h,ihi,ihi,1d0)))

       ! remove duplicates from the list of structures
       do jj = ii+1, nst
          j = iord(jj)
          if (.not.active(j)) cycle

          adv = abs(vst(i)-vst(j)) * bohrtoa**3    ! to ang^3
          adh = abs(hst(i)-hst(j)) * 0.043364104d0 ! to eV
          if (adv > vdiff_thr) then
             ! sorted by volume, so it is pointless to continue
             exit
          else if (adh < ediff_thr) then
             ! powder for structure j
             call st(j)%powder(th2ini,xend,.false.,npts,lambda0,fpol0,sigma,t,ihj,th2p,ip,hvecp)
             tini = ihj(1)*ihj(1)
             tend = ihj(npts)*ihj(npts)
             nor = (2d0 * sum(ihj(2:npts-1)*ihj(2:npts-1)) + tini + tend) * (xend - th2ini) / 2d0 / real(npts-1,8)
             ihj = ihj / sqrt(nor)
             xnormj = sqrt(abs(crosscorr_triangle(h,ihj,ihj,1d0)))

             ! compare structures i and j
             diffij = max(1d0 - crosscorr_triangle(h,ihi,ihj,1d0) / xnormi / xnormj,0d0)
             if (diffij < powdiff_thr) then
                write (uout,'(A,"  PRUNED  because: same as ",A," (deltaV=",A,",deltaH=",A,",pow=",A,")")') &
                   string(st(j)%file,maxlen), string(st(i)%file), string(adv,'e',decimal=4),&
                   string(adh,'e',decimal=4), string(diffij,'e',decimal=4)
                active(j) = .false.
             end if
          end if
       end do

       ! this structure is active, so write it
       write (uout,'(A,"   KEPT   with: (V=",A,",H=",A,")")') string(st(i)%file,maxlen), &
          string(vst(i),'f',decimal=3), string(hst(i),'f',decimal=3)
    end do main
    write (uout,'("+ Structures (pruned/written/total): ",A,"/",A,"/",A)') string(nst-count(active)), &
       string(count(active)), string(nst)
    write (uout,*)

    return
999 continue

    call ferror('trick_reduce',errmsg,faterr,line,syntax=.true.)
    return

  end subroutine trick_reduce

  ! Write a cif file for the current structure using the CCDC
  ! blind test instructions. Use:
  !  TRICK MAKECIF file.s [irank.i]
  ! where irank is the rank for the title of the data_ block. If
  ! no rank is present, use just "data_".
  subroutine trick_makecif_ccdc(line0)
    use systemmod, only: sy
    use tools_io, only: ferror, faterr, fopen_write, fclose, string,&
       nameguess, equal, isinteger, getword
    use param, only: maxzat, mlen, bohrtoa
    character*(*), intent(in) :: line0

    logical :: usesym, ok
    integer :: i, j, idx, lu, natmol, is, irank
    integer, allocatable :: atc(:,:)
    character(len=:), allocatable :: str, hmpg, file
    character*2 :: sym
    character*12 :: holo
    real*8 :: matdum(3,3)
    character(len=mlen), allocatable :: strfin(:)
    integer :: datvalues(8), lp

    ! Hill order for chemical formula. First C, then H, then all the other
    ! elements in alphabetical order.
    integer, parameter :: hillord(maxzat) = &
       (/6,1,89,47,13,95,18,33,85,79,5,56,4,107,83,97,35,20,48,58,98,17,&
       96,112,27,24,55,29,105,110,66,68,99,63,9,26,114,100,87,31,64,32,2,72,80,&
       67,108,53,49,77,19,36,57,3,103,71,116,115,101,12,25,42,109,7,11,41,60,10,&
       113,28,102,93,8,118,76,15,91,82,46,61,84,59,78,94,88,37,75,104,111,45,86,&
       44,16,51,21,34,106,14,62,50,38,73,65,43,52,90,22,81,69,117,92,23,74,54,&
       39,70,30,40/)

    ! consistency checks
    if (.not.associated(sy)) &
       call ferror('trick_makecif_ccdc','system not defined',faterr)
    if (.not.associated(sy%c)) &
       call ferror('trick_makecif_ccdc','crystal structure not defined',faterr)
    if (.not.sy%c%isinit) &
       call ferror('trick_makecif_ccdc','crystal structure not initialized',faterr)
    if (sy%c%ismolecule) &
       call ferror('trick_makecif_ccdc','structure is a molecule',faterr)
    if (.not.sy%c%ismol3d) &
       call ferror('trick_makecif_ccdc','structure is not a molecular crystal',faterr)

    ! transform to the standard cell
    matdum = sy%c%cell_standard(.false.,.false.,.true.)

    ! use symmetry?
    usesym = sy%c%spgavail

    ! open output file
    lp = 1
    file = getword(line0,lp)
    ok = isinteger(irank,line0,lp)
    if (.not.ok) then
       irank = -1
    end if

    lu = fopen_write(file)

    ! header
    if (irank > 0) then
       write (lu,'("data_",A)') string(irank)
    else
       write (lu,'("data_")')
    end if

    ! timestamp
    call date_and_time(values=datvalues)
    write (lu,'("_audit_creation_date ",A,"-",A,"-",A)') &
       (string(datvalues(i)),i=1,3)

    ! chemical formula sum
    allocate(atc(maxzat,sy%c%nmol))
    atc = 0
    do i = 1, sy%c%nmol
       do j = 1, sy%c%mol(i)%nat
          atc(sy%c%mol(i)%spc(sy%c%mol(i)%at(j)%is)%z,i) = atc(sy%c%mol(i)%spc(sy%c%mol(i)%at(j)%is)%z,i) + 1
       end do
       if (i > 2) then
          if (any(atc(:,i) - atc(:,1) /= 0)) then
             call ferror('trick_makecif_ccdc','inconsistent molecular fragments',faterr)
          end if
       end if
    end do
    str = ""
    do i = 1, maxzat
       idx = hillord(i)
       if (atc(idx,1) > 0) then
          sym = nameguess(idx,.true.)
          if (atc(idx,1) > 1) then
             str = str // trim(sym) // string(atc(idx,1)) // " "
          else
             str = str // trim(sym) //  " "
          end if
       end if
    end do
    natmol = sum(atc(:,1))
    deallocate(atc)
    str = trim(str)
    write (lu,'("_chemical_formula_sum ''",A,"''")') str

    ! cell dimensions
    write (lu,'("_cell_length_a ",F20.10)') sy%c%aa(1)*bohrtoa
    write (lu,'("_cell_length_b ",F20.10)') sy%c%aa(2)*bohrtoa
    write (lu,'("_cell_length_c ",F20.10)') sy%c%aa(3)*bohrtoa
    write (lu,'("_cell_angle_alpha ",F14.4)') sy%c%bb(1)
    write (lu,'("_cell_angle_beta ",F14.4)') sy%c%bb(2)
    write (lu,'("_cell_angle_gamma ",F14.4)') sy%c%bb(3)
    write (lu,'("_cell_volume ",F20.6)') sy%c%omega * bohrtoa**3

    ! number of molecules
    if (abs(real(sy%c%ncel,8)/real(natmol,8) - sy%c%ncel/natmol) > 1d-10) &
       call ferror('trick_makecif_ccdc','non-integer Z',faterr)
    write (lu,'("_cell_formula_units_Z ",A)') string(sy%c%ncel/natmol)

    ! get the strings for the symmetry operations
    if (usesym) then
       allocate(strfin(sy%c%neqv*sy%c%ncv))
       call sy%c%struct_report_symxyz(strfin)
       do i = 1, sy%c%neqv*sy%c%ncv
          if (index(strfin(i),"not found") > 0) then
             usesym = .false.
             exit
          end if
       end do
    end if

    ! crystal system
    hmpg = sy%c%spg%pointgroup_symbol
    if (equal(hmpg,"")) then
       holo = "unknown     "
    elseif (equal(hmpg,"1").or.equal(hmpg,"-1")) then
       holo = "triclinic   "
    elseif (equal(hmpg,"2").or.equal(hmpg,"m").or.equal(hmpg,"2/m")) then
       holo = "monoclinic  "
    elseif (equal(hmpg,"222").or.equal(hmpg,"mm2").or.equal(hmpg,"mmm")) then
       holo = "orthorhombic"
    elseif (equal(hmpg,"4").or.equal(hmpg,"-4").or.equal(hmpg,"4/m").or.&
       equal(hmpg,"422").or.equal(hmpg,"4mm").or.equal(hmpg,"-42m").or.&
       equal(hmpg,"4/mmm")) then
       holo = "tetragonal  "
    elseif (equal(hmpg,"3").or.equal(hmpg,"-3").or.equal(hmpg,"32").or.&
       equal(hmpg,"3m").or.equal(hmpg,"-3m")) then
       holo = "trigonal    "
    elseif (equal(hmpg,"6").or.equal(hmpg,"-6").or.equal(hmpg,"6/m").or.&
       equal(hmpg,"622").or.equal(hmpg,"6mm").or.equal(hmpg,"-6m2").or.&
       equal(hmpg,"6/mmm")) then
       holo = "hexagonal   "
    elseif (equal(hmpg,"23").or.equal(hmpg,"m-3").or.equal(hmpg,"432").or.&
       equal(hmpg,"-43m").or.equal(hmpg,"m-3m")) then
       holo = "cubic       "
    end if

    ! write the symmetry, if applicable
    if (usesym) then
       write (lu,'("_space_group_crystal_system ",A)') string(holo)
       write (lu,'("_space_group_name_H-M_alt ''",A,"''")') string(sy%c%spg%international_symbol)
       write (lu,'("_space_group_IT_number ",A)') string(sy%c%spg%spacegroup_number)

       write (lu,'("loop_")')
       write (lu,'("_symmetry_equiv_pos_site_id")')
       write (lu,'("_symmetry_equiv_pos_as_xyz")')
       do i = 1, sy%c%neqv*sy%c%ncv
          write (lu,'(" ",A," ''",A,"''")') string(i), string(strfin(i))
       end do
    else
       write (lu,'("_space_group_crystal_system triclinic")')
       write (lu,'("_space_group_name_H-M_alt ''P 1''")')
       write (lu,'("_space_group_IT_number 1")')

       write (lu,'("loop_")')
       write (lu,'("_symmetry_equiv_pos_site_id")')
       write (lu,'("_symmetry_equiv_pos_as_xyz")')
       write (lu,'(" 1 ''x,y,z''")')
    end if

    write (lu,'("loop_")')
    write (lu,'("_atom_site_label")')
    write (lu,'("_atom_site_type_symbol")')
    write (lu,'("_atom_site_fract_x")')
    write (lu,'("_atom_site_fract_y")')
    write (lu,'("_atom_site_fract_z")')
    if (usesym) then
       do i = 1, sy%c%nneq
          is = sy%c%at(i)%is
          write (lu,'(A5,X,A3,X,3(F20.14,X))') sy%c%spc(is)%name, nameguess(sy%c%spc(is)%z,.true.), sy%c%at(i)%x
       end do
    else
       do i = 1, sy%c%ncel
          is = sy%c%atcel(i)%is
          write (lu,'(A5,X,A3,X,3(F20.14,X))') sy%c%spc(is)%name, nameguess(sy%c%spc(is)%z,.true.), sy%c%atcel(i)%x
       end do
    end if
    call fclose(lu)

  end subroutine trick_makecif_ccdc

  ! Adapted from QE. Copyright (C) 2003-2010 Quantum ESPRESSO group.
  ! See the bfgs_module.f90 in the QE distribution for more details
  ! re implementation and authorship.
  !   trick bfgs file.geom file.bfgs file.calc [tightqe|looseqe|aims2]
  ! - file.geom contains the geometry and is updated unless the BFGS has converged.
  ! - file.bfgs contains the BFGS history
  ! - file.calc must contain:
  !   1 line, the energy in Ry
  !   3 lines, the stress tensor in Ry/bohr^3
  !   nat lines, the atomic forces in Ry/bohr
  ! - looseqe : the QE default convergence criteria (default)
  ! - tightqe : 10x tigther force and energy criteria
  ! - aims2 : FHIaims, 0.01 eV/ang
  subroutine trick_bfgs(line0)
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use tools_io, only: getword, ferror, faterr, fopen_read, fclose, fopen_write, uout,&
       string, ioj_left, ioj_center, lower, equal
    use tools_math, only: matinv
    character*(*), intent(in) :: line0

    type(crystalseed) :: seed
    type(crystal) :: c
    character(len=:), allocatable :: filegeom, filebfgs, filecalc
    character(len=:), allocatable :: errmsg, word, lword
    integer :: lp, lu
    integer :: i, j, k
    logical :: lwolfe, ok, conv_bfgs, energy_wolfe_condition, gradient_wolfe_condition
    integer :: n, nat
    real*8 :: energy, energysave, st(3,3), dE0s, den
    real*8 :: h(3,3), hinv(3,3), omega, g(3,3), ginv(3,3), fcell(3,3)
    integer :: scf_iter, bfgs_iter, tr_min_hit
    real*8 :: energy_p, nr_step_length_old, nr_step_length, trust_radius_old, trust_radius
    real*8 :: energy_error, grad_error, cell_error
    real*8, allocatable :: force(:,:)
    real*8, allocatable :: pos(:), grad(:), pos_old(:), grad_old(:), inv_hess(:,:)
    real*8, allocatable :: pos_p(:), grad_p(:), step(:), step_old(:), pos_best(:)
    real*8, allocatable :: hinv_block(:,:), metric(:,:)
    character*10 :: msgsum
    real*8 :: energy_thr, grad_thr, cell_thr, trust_radius_min
    ! bfgs history
    integer :: nstep
    real*8, allocatable :: step_tr(:), step_ene(:), step_de(:), step_df(:), step_dc(:)
    integer, allocatable :: step_bfgs_it(:), step_scf_it(:), step_idec(:)
    integer :: idec

    real*8, parameter :: kbarau = 147105.08d0
    real*8, parameter :: trust_radius_ini = 0.5d0

    real*8, parameter :: energy_thr_looseqe = 1d-4
    real*8, parameter :: grad_thr_looseqe = 1d-3
    real*8, parameter :: cell_thr_looseqe = 0.5d0 / kbarau
    real*8, parameter :: trust_radius_min_looseqe = 1d-4

    real*8, parameter :: energy_thr_tightqe = 1d-5
    real*8, parameter :: grad_thr_tightqe = 1d-4
    real*8, parameter :: cell_thr_tightqe = 0.5d0 / kbarau
    real*8, parameter :: trust_radius_min_tightqe = 1d-5

    real*8, parameter :: energy_thr_aims2 = 3.9d-5
    real*8, parameter :: grad_thr_aims2 = 3.9d-4
    real*8, parameter :: cell_thr_aims2 = 0.5d0 / kbarau
    real*8, parameter :: trust_radius_min_aims2 = 1d-4

    real*8, parameter :: w_1 = 0.01d0
    real*8, parameter :: w_2 = 0.5d0
    real*8, parameter :: trust_radius_max = 0.8d0

    ! header and initialization
    write (uout,'("* BFGS optimization (vc-relax), algorithm from QE")')
    energy_thr = energy_thr_looseqe
    grad_thr = grad_thr_looseqe
    cell_thr = cell_thr_looseqe
    trust_radius_min = trust_radius_min_looseqe

    ! parse input
    lp = 1
    filegeom = getword(line0,lp)
    filebfgs = getword(line0,lp)
    filecalc = getword(line0,lp)
    do while (.true.)
       word = getword(line0,lp)
       lword = lower(word)
       if (equal(lword,'looseqe')) then
          energy_thr = energy_thr_looseqe
          grad_thr = grad_thr_looseqe
          cell_thr = cell_thr_looseqe
          trust_radius_min = trust_radius_min_looseqe
       elseif (equal(lword,'tightqe')) then
          energy_thr = energy_thr_tightqe
          grad_thr = grad_thr_tightqe
          cell_thr = cell_thr_tightqe
          trust_radius_min = trust_radius_min_tightqe
       elseif (equal(lword,'aims2')) then
          energy_thr = energy_thr_aims2
          grad_thr = grad_thr_aims2
          cell_thr = cell_thr_aims2
          trust_radius_min = trust_radius_min_aims2
       elseif (len_trim(word) > 0) then
          call ferror('trick_bfgs','Unknown extra keyword',faterr,line0,syntax=.true.)
       else
          exit
       end if
    end do

    ! read crystal structure
    write (uout,'("+ Reading the structure from: ",A)') trim(filegeom)
    call seed%read_any_file(filegeom,-1,errmsg)
    if (len_trim(errmsg) > 0) then
       call ferror('trick_bfgs','error reading geometry file: ' // filegeom,faterr)
    end if
    call c%struct_new(seed,.true.)

    ! initialize
    lwolfe=.false.
    nat = c%ncel
    n = 3 * nat + 9

    ! allocate work space
    allocate(pos(n))
    allocate(grad(n))
    allocate(grad_old(n))
    allocate(pos_old(n))
    allocate(inv_hess(n,n))
    allocate(pos_p(n))
    allocate(grad_p(n))
    allocate(step(n))
    allocate(step_old(n) )
    allocate(pos_best(n) )
    allocate(hinv_block(n-9,n-9))
    allocate(metric(n,n))

    ! h, hinv, volume
    h = c%m_x2c
    hinv = c%m_c2x
    omega = c%omega
    hinv_block = 0.d0
    forall (k=0:nat-1,i=1:3,j=1:3) hinv_block(i+3*k,j+3*k) = hinv(i,j)

    ! generate metric to work with scaled ionic coordinates
    g = c%gtensor
    ginv = c%grtensor
    metric = 0.d0
    forall (k=0:nat-1,i=1:3,j=1:3) metric(i+3*k,j+3*k) = g(i,j)
    forall (k=nat:nat+2,i=1:3,j=1:3) metric(i+3*k,j+3*k) = 0.04d0 * omega * ginv(i,j)

    ! read the calc file
    write (uout,'("+ Reading the calculated values from: ",A)') trim(filecalc)
    lu = fopen_read(filecalc)
    allocate(force(3,nat))
    read (lu,*) energy, st, force
    call fclose(lu)
    energysave = energy

    ! calculate the cell force
    do j=1,3
       do i=1,3
          fcell(i,j) = -hinv(j,1)*st(i,1) - hinv(j,2)*st(i,2) - hinv(j,3)*st(i,3)
       end do
    end do
    fcell = omega * fcell

    ! convert forces to Cartesian gradient
    do i = 1, nat
       force(:,i) = -matmul(force(:,i),h)
    end do

    ! generate bfgs vectors for the degrees of freedom and their gradients
    pos = 0.0d0
    forall (k=0:nat-1,i=1:3) pos(3*k+i) = c%atcel(k+1)%x(i)
    forall (i=1:3,j=1:3) pos(n-9+j+3*(i-1)) = h(i,j)
    grad = 0.0
    forall (k=0:nat-1,i=1:3) grad(3*k+i) = force(i,k+1)
    forall (i=1:3,j=1:3) grad(n-9+j+3*(i-1)) = fcell(i,j)

    ! read the bfgs file, if it exists
    inquire(file=filebfgs,exist=ok)
    if (ok) then
       write (uout,'("+ Reading the BFGS file: ",A)') trim(filebfgs)
       lu =  fopen_read(filebfgs,"unformatted")
       read (lu) pos_p
       read (lu) grad_p
       read (lu) scf_iter
       read (lu) bfgs_iter
       read (lu) energy_p
       read (lu) pos_old
       read (lu) grad_old
       read (lu) inv_hess
       read (lu) tr_min_hit
       read (lu) nr_step_length
       read (lu) nstep
       allocate(step_tr(nstep+1),step_ene(nstep+1),step_de(nstep+1),step_df(nstep+1),step_dc(nstep+1))
       allocate(step_bfgs_it(nstep+1),step_scf_it(nstep+1),step_idec(nstep+1))
       read (lu) step_tr(1:nstep)
       read (lu) step_ene(1:nstep)
       read (lu) step_de(1:nstep)
       read (lu) step_df(1:nstep)
       read (lu) step_dc(1:nstep)
       read (lu) step_bfgs_it(1:nstep)
       read (lu) step_scf_it(1:nstep)
       read (lu) step_idec(1:nstep)
       call fclose(lu)

       step_old = pos - pos_p
       trust_radius_old = scnorm(step_old)
       step_old = step_old / trust_radius_old
    else
       ! initialize the inv_hess to the inverse of the metric
       write (uout,'("+ The BFGS file ",A," does not exist, initializing")') trim(filebfgs)
       inv_hess = metric
       call matinv(inv_hess,n)
       pos_p = 0d0
       grad_p = 0d0
       scf_iter = 0
       bfgs_iter = 0
       energy_p = energy
       step_old = 0d0
       nr_step_length = 0d0
       trust_radius_old = trust_radius_ini
       tr_min_hit = 0

       ! initialize the BFGS history (for reporting)
       nstep = 0
       allocate(step_tr(1),step_ene(1),step_de(1),step_df(1),step_dc(1))
       allocate(step_bfgs_it(1),step_scf_it(1),step_idec(1))
    end if
    scf_iter = scf_iter + 1

    ! check convergence
    energy_error = abs(energy_p - energy)
    grad_error = maxval(abs(matmul(transpose(hinv_block),grad(1:n-9))))
    cell_error = maxval(abs(matmul(transpose(reshape(grad(n-8:n), (/3,3/))),transpose(h)))) / omega
    conv_bfgs = energy_error < energy_thr
    conv_bfgs = conv_bfgs .and. (grad_error < grad_thr)
    conv_bfgs = conv_bfgs .and. (cell_error < cell_thr)
    IF (.not.conv_bfgs .and. (tr_min_hit > 1)) then
       write (uout,'("history already reset at previous step: stopping")')
       goto 99
    end IF
    IF (conv_bfgs) goto 99

    ! some output
    write (uout,'("  Number of scf cycles = ",A)') string(scf_iter)
    write (uout,'("  Number of bfgs steps = ",A)') string(bfgs_iter)
    if (scf_iter > 1) &
       write (uout,'("  Old energy = ",A)') string(energy_p,'f',decimal=10)
    write (uout,'("  New energy = ",A)') string(energy,'f',decimal=10)

    ! wolfe conditions
    energy_wolfe_condition = (energy-energy_p) < w_1 * (dot_product(grad_p,step_old)) * trust_radius_old
    if (.not.energy_wolfe_condition .and. scf_iter > 1) then
       ! failed energy condition, keep looking
       write (uout,'("+ Step not accepted, line search continues")')
       step = step_old
       dE0s = dot_product(grad_p,step) * trust_radius_old
       den = energy - energy_p - dE0s

       ! estimate new trust radius by interpolation
       trust_radius = - 0.5d0*dE0s*trust_radius_old / den
       write (uout,'("  Old trust radius = ",A)') string(trust_radius_old,'f',decimal=10)
       write (uout,'("  New trust radius = ",A)') string(trust_radius,'f',decimal=10)

       ! values from the last successful bfgs step are restored
       pos  = pos_p
       energy  = energy_p
       grad = grad_p

       if (trust_radius < trust_radius_min) then
          idec = 2 ! l-s,low tr
          ! history reset (at most twice in a row)
          write (uout,'("+ trust_radius < trust_radius_min, resetting BFGS history")')

          ! ... if tr_min_hit=1 the history has already been reset at the
          ! ... previous step : something is going wrong
          if ( tr_min_hit == 1 ) then
             write (uout,'("!! BFGS history already reset at previous step: stopping")')
             tr_min_hit = 2
          else
             tr_min_hit = 1
          end if
          inv_hess = metric
          call matinv(inv_hess,n)
          step = -matmul(inv_hess,grad)

          ! normalize step but remember its length
          nr_step_length = scnorm(step)
          step = step / nr_step_length
          trust_radius = min(trust_radius_ini, nr_step_length)
       else
          idec = 1 ! linesearch
          tr_min_hit = 0
       end if
    else
       ! a brave new step
       idec = 0 ! new step
       write (uout,'("+ Step accepted")')
       bfgs_iter = bfgs_iter + 1
       if (bfgs_iter > 1) then
          nr_step_length_old = nr_step_length

          ! check wolfe conditions
          gradient_wolfe_condition = abs(dot_product(grad,step_old)) < -w_2 * dot_product(grad_p,step_old)
          lwolfe = energy_wolfe_condition .and. gradient_wolfe_condition

          ! update inverse hessian
          call update_inverse_hessian(pos,grad,n)
       end if

       ! newton-raphson
       step = -matmul(inv_hess,grad)

       ! check uphill step condition
       if (dot_product(grad,step) > 0d0) then
          write (uout,'("!! Uphill step: resetting BFGS history")')
          inv_hess = metric
          call matinv(inv_hess,n)
          step = -matmul(inv_hess,grad)
       end if

       ! normalize
       nr_step_length = scnorm(step)
       step = step / nr_step_length

       ! new trust radius
       if (bfgs_iter == 1) then
          trust_radius = min(trust_radius_ini, nr_step_length)
          tr_min_hit = 0
       else
          call compute_trust_radius(lwolfe,energy,grad,n)
       end if

       write (uout,'("  Old trust radius = ",A)') string(trust_radius_old,'f',decimal=10)
       write (uout,'("  New trust radius = ",A)') string(trust_radius,'f',decimal=10)
    end if

    ! check step length
    if (nr_step_length < 1d-16) &
       call ferror("trick_bfgs","step length too low",faterr)

    ! write to the BFGS history
    nstep = nstep + 1
    step_tr(nstep) = trust_radius
    if (idec == 0) then
       step_ene(nstep) = energy
    else
       step_ene(nstep) = energysave
    endif
    step_de(nstep) = energy_error
    step_df(nstep) = grad_error
    step_dc(nstep) = cell_error
    step_bfgs_it(nstep) = bfgs_iter
    step_scf_it(nstep) = scf_iter
    step_idec(nstep) = idec

    ! write the bfgs file
    write (uout,'("+ Writing the BFGS file: ",A)') trim(filebfgs)
    lu =  fopen_write(filebfgs,"unformatted",errstop=.true.)
    write (lu) pos
    write (lu) grad
    write (lu) scf_iter
    write (lu) bfgs_iter
    write (lu) energy
    write (lu) pos_old
    write (lu) grad_old
    write (lu) inv_hess
    write (lu) tr_min_hit
    write (lu) nr_step_length
    write (lu) nstep
    write (lu) step_tr
    write (lu) step_ene
    write (lu) step_de
    write (lu) step_df
    write (lu) step_dc
    write (lu) step_bfgs_it
    write (lu) step_scf_it
    write (lu) step_idec
    call fclose(lu)

    ! update positions
    pos = pos + trust_radius * step
    forall(i=1:3,j=1:3) h(i,j) = pos(n-9+j+3*(i-1))

    ! write the new crystal structure
    write (uout,'("+ Over-writing the geometry file: ",A)') trim(filegeom)
    call c%makeseed(seed,.false.)
    seed%m_x2c = h
    forall (k=0:nat-1,i=1:3) seed%x(i,k+1) = pos(3*k+i)
    call c%struct_new(seed,.true.)
    call c%write_simple_driver(filegeom)

99  continue ! skip here if converged

    ! bfgs history
    write (uout,*)
    write (uout,'("+ BFGS history")')
    write (uout,'("# Energies in Ry, gradients in Ry/bohr, stress in kbar")')
    write (uout,'("#i bfgs scf step_type trust_radius    energy(Ry)       delta-E     C?   delta-grad  C?   delta-cell  C?")')
    do i = 1, nstep
       if (step_idec(i) == 0) then
          msgsum = " new step "
       elseif (step_idec(i) == 1) then
          msgsum = "linesearch"
       elseif (step_idec(i) == 2) then
          msgsum = "l-s,low tr"
       end if

       write (uout,'(99(A,X))') string(i,3,ioj_left), string(step_bfgs_it(i),3,ioj_left), string(step_scf_it(i),3,ioj_left), &
          msgsum, string(step_tr(i),'f',12,8), string(step_ene(i),'f',17,10),&
          string(step_de(i),'e',12,5), string(string(step_de(i) < energy_thr),3,ioj_center),&
          string(step_df(i),'e',12,5), string(string(step_df(i) < grad_thr),3,ioj_center),&
          string(step_dc(i)*kbarau,'e',12,5), string(string(step_dc(i) < cell_thr),3,ioj_center)
    end do
    if (conv_bfgs) then
       ! last step
       write (uout,'(99(A,X))') string(nstep+1,3,ioj_left), string(bfgs_iter,3,ioj_left),&
          string(scf_iter,3,ioj_left), &
          msgsum, string(trust_radius,'f',12,5), string(energy,'f',17,10),&
          string(energy_error,'e',12,5), string(string(energy_error < energy_thr),3,ioj_center),&
          string(grad_error,'e',12,5), string(string(grad_error < grad_thr),3,ioj_center),&
          string(cell_error*kbarau,'e',12,5), string(string(cell_error < cell_thr),3,ioj_center)
    end if

    ! summary
    write (uout,*)
    write (uout,'("+ Convergence summary")')
    write (uout,'("  Trust radius (min/current/max) = ",3(A,X))') string(trust_radius_min,'f',12,5),&
         string(trust_radius,'f',12,5), string(trust_radius_max,'f',12,5)
    write (uout,'("  Energy (Ry) = ",A," (thr=",A,")")') string(energy_error,'e',decimal=5), &
       string(energy_thr,'e',decimal=5)
    write (uout,'("  Gradient (Ry/bohr) = ",A," (thr=",A,")")') string(grad_error,'e',decimal=5), &
       string(grad_thr,'e',decimal=5)
    write (uout,'("  Stress (Ry/bohr^3) = ",A," (thr=",A,")")') string(cell_error*kbarau,'e',decimal=5), &
       string(cell_thr*kbarau,'e',decimal=5)
    if (conv_bfgs) then
       write (uout,'("+ BFGS is CONVERGED")')
    else
       write (uout,'("+ BFGS not finished")')
    end if

    ! wrap up
    write (uout,*)

    ! crash if minimum step reached
    if (.not.conv_bfgs .and. (tr_min_hit > 1)) &
       call ferror('trick_bfgs',"history already reset at previous step: stopping",faterr)

  contains

    function scnorm(vect)
      real*8 :: scnorm
      real*8, intent(in) :: vect(:)
      real*8 :: ss
      integer :: i,k,l,n

      scnorm = 0d0
      n = size(vect) / 3
      do i=1,n
         ss = 0d0
         do k=1,3
            do l=1,3
               ss = ss + vect(k+(i-1)*3)*metric(k+(i-1)*3,l+(i-1)*3)*vect(l+(i-1)*3)
            end do
         end do
         scnorm = max(scnorm, sqrt(ss))
      end do

    end function scnorm

    subroutine update_inverse_hessian(pos,grad,n)
      real*8, intent(in) :: pos(:)
      real*8, intent(in) :: grad(:)
      integer, intent(in) :: n
      !
      integer :: info
      real*8, allocatable :: y(:), s(:)
      real*8, allocatable :: Hy(:), yH(:)
      real*8 :: sdoty, sBs, Theta
      real*8, allocatable :: B(:,:)

      allocate(y(n), s(n), Hy(n), yH(n))
      s = pos - pos_p
      y = grad - grad_p
      sdoty = dot_product(s,y)

      if (abs(sdoty) < 1d-16) then
         write (*,*) "unexpected behavior in update_inverse_hessian, resetting bfgs history"
         inv_hess = metric
         call matinv(inv_hess,n)
         return
      else
         ! Conventional Curvature Trap here
         ! See section 18.2 (p538-539 ) of Nocedal and Wright "Numerical
         ! Optimization"for instance
         ! LDM Addition, April 2011
         !
         ! While with the Wolfe conditions the Hessian in most cases
         ! remains positive definite, if one is far from the minimum
         ! and/or "bonds" are being made/broken the curvature condition
         !        Hy = s ; or s = By
         ! cannot be satisfied if s.y < 0. In addition, if s.y is small
         ! compared to s.B.s too greedy a step is taken.
         !
         ! The trap below is conventional and "OK", and has been around
         ! for ~ 30 years but, unfortunately, is rarely mentioned in
         ! introductory texts and hence often neglected.
         !
         ! First, solve for inv_hess*t = s ; i.e. t = B*s
         ! Use yH as workspace here

         allocate(B(n,n))
         B = inv_hess
         yH = s
         call dposv('U',n,1,B,n,yH,n,info)
         ! Info .ne. 0 should be trapped ...
         if (info/=0) write (uout,*) "info = ", info, "for hessian"
         deallocate(B)
         ! Calculate s.B.s
         sBs = dot_product(s,yH)

         ! Now the trap itself
         if (sdoty < 0.20D0*sBs) then
            ! Conventional damping
            Theta = 0.8D0*sBs/(sBs-sdoty)
            write (uout,*) "warning: bfgs curvature condition failed, theta= ", theta
            y = Theta*y + (1.D0 - Theta)*yH
         endif
      end if
      Hy = matmul(inv_hess,y)
      yH = matmul(y,inv_hess)

      ! BFGS update
      inv_hess = inv_hess + 1.d0 / sdoty * ((1d0 + dot_product(y,Hy) / sdoty) * matrix(s,s) - &
         (matrix(s,yH) + matrix(Hy,s)))

    end subroutine update_inverse_hessian

    subroutine compute_trust_radius(lwolfe,energy,grad,n)
      logical, intent(in) :: lwolfe
      real*8, intent(in) :: energy
      real*8, intent(in) :: grad(:)
      integer, intent(in) :: n

      real*8 :: a
      logical :: ltest

      ltest = (energy - energy_p) < w_1 * dot_product(grad_p,step_old) * trust_radius_old

      ! The instruction below replaces the original instruction:
      !    ltest = ltest .AND. ( nr_step_length_old > trust_radius_old )
      ! which gives a random result if trust_radius was set equal to
      ! nr_step_length at previous step. I am not sure what the best
      ! action should be in that case, though (PG)
      !
      ltest = ltest .and. ( nr_step_length_old > trust_radius_old + 1d-8 )

      if (ltest) then
         a = 1.5d0
      else
         a = 1.1d0
      end if
      if (lwolfe) a = 2d0 * a

      trust_radius = min(trust_radius_max, a*trust_radius_old, nr_step_length)

      if (trust_radius < trust_radius_min) then
         ! ... the history is reset
         !
         ! ... if tr_min_hit the history has already been reset at the
         ! ... previous step : something is going wrong
         !
         if (tr_min_hit == 1) then
            write (uout,*) 'history already reset at previous step: stopping'
            tr_min_hit = 2
         else
            tr_min_hit = 1
         end if
         write (uout,*) "small trust radius: resetting bfgs history"
         inv_hess = metric
         call matinv(inv_hess,n)
         step = -matmul(inv_hess,grad)
         nr_step_length = scnorm(step)
         step = step / nr_step_length
         trust_radius = min(trust_radius_min, nr_step_length )
      else
         tr_min_hit = 0
      end if

    end subroutine compute_trust_radius

    function matrix( vec1, vec2 )
      real*8, intent(in) :: vec1(:), vec2(:)
      real*8 :: matrix(SIZE(vec1),SIZE(vec2))
      integer :: dim1, dim2

      dim1 = size( vec1 )
      dim2 = size( vec2 )
      matrix = 0d0
      call dger( dim1, dim2, 1d0, vec1, 1, vec2, 1, matrix, dim1)

    end function matrix

  end subroutine trick_bfgs

  ! Compare structures, allowing deformation of one crystal into the other such
  ! as may be cause by temperature or pressure effects without a phase change.
  !   TRICK COMPARE struct1 struct2 [THR thr.r] [WRITE] [NOH]
  ! If thr.r, then stop if the POWDIFF is below thr.r. If WRITE, write stucture 2
  ! to a file. If NOH, strip the hydrogens from the structure before comparing.
  subroutine trick_compare_deformed(line0)
    use iso_c_binding, only: c_double
    use spglib, only: spg_delaunay_reduce, spg_standardize_cell
    use environmod, only: environ
    use global, only: iunitname0, dunit0, iunit, fileroot
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use tools_math, only: matinv, m_c2x_from_cellpar, det3, crosscorr_triangle, &
       m_x2c_from_cellpar
    use tools_io, only: getword, faterr, ferror, uout, string, ioj_left, ioj_right,&
       isreal, equal, lgetword
    use param, only: pi, icrd_crys, eye
    character*(*), intent(in) :: line0

    type(crystalseed) :: seed
    type(environ) :: e
    integer :: lp, ierr, i, j
    character(len=:), allocatable :: file1, file2, errmsg, abc, word
    type(crystal) :: c1, c2, c2del
    real*8 :: xd2(3,3), cd2(3,3), dmax0, xx(3)
    real*8 :: aa2(3), bb2(3), cc2(3), dd
    real*8, allocatable :: dist(:)
    integer, allocatable :: eid(:), irange(:,:)
    integer :: nat, n1, n2, n3, i1, i2, i3
    real*8, allocatable :: iha1(:), iha2(:)
    real*8, allocatable :: t(:), th2p(:), ip(:)
    integer, allocatable :: hvecp(:,:)
    real*8 :: tini, tend, nor, diff, xnorm1, xnorm2, h, mindiff
    real*8 :: x0std1(3,3), x0std2(3,3), x0del1(3,3), x0del2(3,3), xd2min(3,3)
    logical :: ok, dowrite, noh
    real*8 :: powdiff_thr

    real*8, parameter :: th2ini = 5d0
    real*8, parameter :: th2end = 50d0
    integer, parameter :: npts = 1001
    real*8, parameter :: lambda0 = 1.5406d0
    real*8, parameter :: fpol0 = 0d0
    real*8, parameter :: sigma0 = 0.05d0

    real*8, parameter :: max_elong = 0.3d0 ! at most 30% elongation of cell lengths
    real*8, parameter :: max_ang = 20d0    ! at most 20 degrees change in angle
    character*1, parameter :: lvecname(3) = (/"a","b","c"/)

    ! header and initalization
    write (uout,'("* COMPARE, allowing for deformed cells")')

    ! read the input files
    lp = 1
    file1 = getword(line0,lp)
    if (len_trim(file1) == 0) &
       call ferror('trick_compare_deformed','Missing first structure file',faterr)
    file2 = getword(line0,lp)
    if (len_trim(file2) == 0) &
       call ferror('trick_compare_deformed','Missing second structure file',faterr)

    powdiff_thr = -1d0
    dowrite = .false.
    noh = .false.
    do while (.true.)
       word = lgetword(line0,lp)
       if (equal(word,'thr')) then
          ok = isreal(powdiff_thr,line0,lp)
          if (.not.ok) call ferror('trick_compare_deformed','Wrong THR',faterr)
       elseif (equal(word,'write')) then
          dowrite = .true.
       elseif (equal(word,'noh')) then
          noh = .true.
       elseif (len_trim(word) > 0) then
          if (.not.ok) call ferror('trick_compare_deformed','Unknown keyword',faterr)
       else
          exit
       end if
    end do

    ! read the structures, force symmetry recalculation
    write (uout,'("+ Reading the structure from: ",A)') trim(file1)
    call seed%read_any_file(file1,-1,errmsg)
    if (len_trim(errmsg) > 0) &
       call ferror('trick_compare_deformed','error reading geometry file: ' // file1,faterr)
    if (noh) call seed%strip_hydrogens()
    call c1%struct_new(seed,.true.)
    call c1%calcsym(.false.,errmsg)
    if (len_trim(errmsg) > 0) &
       call ferror('trick_compare_deformed','error recalculating symmetry: ' // file1,faterr)

    write (uout,'("+ Reading the structure from: ",A)') trim(file2)
    call seed%read_any_file(file2,-1,errmsg)
    if (len_trim(errmsg) > 0) &
       call ferror('trick_compare_deformed','error reading geometry file: ' // file2,faterr)
    if (noh) call seed%strip_hydrogens()
    call c2%struct_new(seed,.true.)
    call c2%calcsym(.false.,errmsg)
    if (len_trim(errmsg) > 0) &
       call ferror('trick_compare_deformed','error recalculating symmetry: ' // file1,faterr)

    ! check both are molecular crystals
    if (c1%ismolecule) &
       call ferror('trick_compare_defomred','structure 1 is a molecule',faterr)
    if (c2%ismolecule) &
       call ferror('trick_compare_defomred','structure 2 is a molecule',faterr)

    ! get the Delaunay cell of both cells
    x0std1 = c1%cell_standard(.true.,.false.,.false.)
    !x0del1 = c1%cell_delaunay()
    x0del1 = c1%cell_niggli()
    x0std2 = c2%cell_standard(.true.,.false.,.false.)
    !x0del2 = c2%cell_delaunay()
    x0del2 = c2%cell_niggli()
    if (all(abs(x0std1) < 1d-5)) x0std1 = eye
    if (all(abs(x0del1) < 1d-5)) x0del1 = eye
    if (all(abs(x0std2) < 1d-5)) x0std2 = eye
    if (all(abs(x0del2) < 1d-5)) x0del2 = eye

    ! some output for the first structure
    write (uout,'("+ Delaunay lattice vectors (",A,")")') iunitname0(iunit)
    write (uout,'("# Structure 1: ")')
    do i = 1, 3
       write (uout,'(4X,A,": ",3(A,X)," length = ",A)') lvecname(i),&
          (string(c1%m_x2c(j,i)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3),&
          string(c1%aa(i)*dunit0(iunit),'f',length=16,decimal=10,justify=5)
    end do
    write (uout,'("  Angles (deg): ",3(A,X))') (string(c1%bb(i),'f',length=8,decimal=3),i=1,3)
    write (uout,'("# Structure 2: ")')
    do i = 1, 3
       c2%aa(i) = norm2(c2%m_x2c(:,i))
       write (uout,'(4X,A,": ",3(A,X)," length = ",A)') lvecname(i),&
          (string(c2%m_x2c(j,i)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3),&
          string(c2%aa(i)*dunit0(iunit),'f',length=16,decimal=10,justify=5)
    end do
    write (uout,'("  Angles (deg): ",3(A,X))') (string(c2%bb(i),'f',length=8,decimal=3),i=1,3)
    write (uout,*)

    ! build the lattice vector environment for the second crystal
    dmax0 = 0d0
    do i = 1, 3
       dmax0 = max(dmax0,norm2(c1%m_x2c(:,i)))
    end do
    dmax0 = dmax0 * (1d0 + max_elong)
    call e%build_lattice(c2%m_x2c,dmax0*2d0)
    call e%list_near_atoms((/0d0,0d0,0d0/),icrd_crys,.true.,nat,ierr,eid=eid,dist=dist,up2d=dmax0*1.25d0,nozero=.true.)

    ! output for the second structure and determine integer ranges
    n1 = 0
    n2 = 0
    n3 = 0
    allocate(irange(nat,3))
    write (uout,'("+ Candidate lattice vectors for structure 2 (referred to the Delaunay basis): ")')
    write (uout,'("#Id        x        y        z       length   used-by")')
    do i = 1, nat
       xx = e%xr2x(e%at(eid(i))%x)

       abc = ""
       if (abs(dist(i) / c1%aa(1) - 1d0) < max_elong) then
          abc = trim(abc) // "1 "
          n1 = n1 + 1
          irange(n1,1) = i
       end if
       if (abs(dist(i) / c1%aa(2) - 1d0) < max_elong) then
          abc = trim(abc) // "2 "
          n2 = n2 + 1
          irange(n2,2) = i
       end if
       if (abs(dist(i) / c1%aa(3) - 1d0) < max_elong) then
          abc = trim(abc) // "3 "
          n3 = n3 + 1
          irange(n3,3) = i
       end if
       if (len_trim(abc) > 0) abc = "(" // trim(abc) // ")"

       write (uout,'(2X,6(A,X))') string(i,3,ioj_left), (string(xx(j),'f',8,2,ioj_right),j=1,3), &
          string(dist(i),'f',12,6,ioj_right), abc

       if (all((dist(i) / c1%aa - 1d0) > max_elong)) exit
    end do
    write (uout,*)

    ! calculate the powder of structure 1
    h =  (th2end-th2ini) / real(npts-1,8)
    call c1%powder(th2ini,th2end,.false.,npts,lambda0,fpol0,sigma0,t,iha1,th2p,ip,hvecp)
    tini = iha1(1)**2
    tend = iha1(npts)**2
    nor = (2d0 * sum(iha1(2:npts-1)**2) + tini + tend) * (th2end - th2ini) / 2d0 / real(npts-1,8)
    iha1 = iha1 / sqrt(nor)
    xnorm1 = crosscorr_triangle(h,iha1,iha1,1d0)
    xnorm1 = sqrt(abs(xnorm1))

    ! calculate the powder of structure 2 and prepare
    c2del = c2
    c2del%aa = c1%aa
    c2del%bb = c1%bb
    c2del%m_x2c = m_x2c_from_cellpar(c2del%aa,c2del%bb)
    c2del%grtensor = matmul(transpose(c2del%m_x2c),c2del%m_x2c)
    call matinv(c2del%grtensor,3)
    do i = 1, 3
       c2del%ar(i) = sqrt(c2del%grtensor(i,i))
    end do
    call c2del%powder(th2ini,th2end,.false.,npts,lambda0,fpol0,sigma0,t,iha2,th2p,ip,hvecp)
    tini = iha2(1)**2
    tend = iha2(npts)**2
    nor = (2d0 * sum(iha2(2:npts-1)**2) + tini + tend) * (th2end - th2ini) / 2d0 / real(npts-1,8)
    iha2 = iha2 / sqrt(nor)
    xnorm2 = crosscorr_triangle(h,iha2,iha2,1d0)
    xnorm2 = sqrt(abs(xnorm2))

    ! calculate baseline powdiff
    mindiff = max(1d0 - crosscorr_triangle(h,iha1,iha2,1d0) / xnorm1 / xnorm2,0d0)
    xd2min = eye

    ! run over all permutations
    write (uout,'("+ Structural comparison of candidate structures")')
    write (uout,'("# Reference structure is 1.")')
    write (uout,'("# Structure 2 takes lattice vectors (a,b,c) from the list above.")')
    write (uout,'("# max-dlen = maximum difference in cell lengths (bohr).")')
    write (uout,'("# max-dang = maximum difference in angles (degree).")')
    write (uout,'("#a  b  c max-dlen max-dang  powdiff")')
    write (uout,'("+ INITIAL POWDIFF = ",A)') string(mindiff,'f',12,9)
    if (mindiff < powdiff_thr) goto 999
    do i1 = 1, n1
       cd2(:,1) = e%xr2c(e%at(eid(irange(i1,1)))%x)
       aa2(1) = norm2(cd2(:,1))
       do i2 = 1, n2
          if (irange(i1,1) == irange(i2,2)) cycle
          cd2(:,2) = e%xr2c(e%at(eid(irange(i2,2)))%x)
          aa2(2) = norm2(cd2(:,2))
          do i3 = 1, n3
             if (irange(i1,1) == irange(i3,3) .or. irange(i2,2) == irange(i3,3)) cycle
             cd2(:,3) = e%xr2c(e%at(eid(irange(i3,3)))%x)
             aa2(3) = norm2(cd2(:,3))

             ! check collinear
             cc2(1) = dot_product(cd2(:,2),cd2(:,3)) / aa2(2) / aa2(3)
             cc2(2) = dot_product(cd2(:,1),cd2(:,3)) / aa2(1) / aa2(3)
             cc2(3) = dot_product(cd2(:,1),cd2(:,2)) / aa2(1) / aa2(2)
             if (any(abs(cc2) > 0.9999d0)) cycle

             ! check angle conditions
             bb2(1) = acos(cc2(1)) * 180d0 / pi
             bb2(2) = acos(cc2(2)) * 180d0 / pi
             bb2(3) = acos(cc2(3)) * 180d0 / pi
             if (any(abs(c1%bb - bb2) > max_ang)) cycle
             xd2 = matmul(c2%m_c2x,cd2)

             ! check determinant
             dd = det3(xd2)
             if (abs(dd) < 1d-5) cycle
             if (dd < 0d0) xd2 = -xd2

             ! make the new crystal
             c2del = c2
             call c2del%newcell(xd2,noenv=.true.)
             c2del%aa = c1%aa
             c2del%bb = c1%bb
             c2del%m_x2c = m_x2c_from_cellpar(c2del%aa,c2del%bb)
             c2del%grtensor = matmul(transpose(c2del%m_x2c),c2del%m_x2c)
             call matinv(c2del%grtensor,3)
             do i = 1, 3
                c2del%ar(i) = sqrt(c2del%grtensor(i,i))
             end do

             ! calculate the powder of structure 2
             call c2del%powder(th2ini,th2end,.false.,npts,lambda0,fpol0,sigma0,t,iha2,th2p,ip,hvecp)
             tini = iha2(1)**2
             tend = iha2(npts)**2
             nor = (2d0 * sum(iha2(2:npts-1)**2) + tini + tend) * (th2end - th2ini) / 2d0 / real(npts-1,8)
             iha2 = iha2 / sqrt(nor)
             xnorm2 = crosscorr_triangle(h,iha2,iha2,1d0)
             xnorm2 = sqrt(abs(xnorm2))

             ! calculate the powdiff
             diff = max(1d0 - crosscorr_triangle(h,iha1,iha2,1d0) / xnorm1 / xnorm2,0d0)
             if (diff < mindiff) then
                mindiff = diff
                xd2min = xd2
             end if

             ! write output
             write (uout,'(99(A,X))') string(irange(i1,1),2,ioj_right), string(irange(i2,2),2,ioj_right),&
                string(irange(i3,3),2,ioj_right),&
                string(maxval(abs(c1%aa-aa2)),'f',8,4,ioj_right), string(maxval(abs(c1%bb-bb2)),'f',8,3,ioj_right),&
                string(diff,'f',10,7)
             if (mindiff < powdiff_thr) goto 999
          end do
       end do
    end do

999 continue
    if (mindiff < powdiff_thr) &
       write (uout,'("--- Last POWDIFF satisfies the threshold requirement for matching structures, skipping...")')
    write (uout,'("+ FINAL POWDIFF = ",A)') string(mindiff,'f',12,9)
    write (uout,*)

    if (dowrite) then
       ! make the minimum-powdiff structure 2
       c2del = c2
       call c2del%newcell(xd2min)
       call c2del%makeseed(seed,.false.)
       seed%useabr = 1
       seed%aa = c1%aa
       seed%bb = c1%bb
       seed%findsym = 0
       call c2del%struct_new(seed,.true.)

       ! write both to a res file
       file1 = fileroot // "_structure_1.res"
       write (uout,'("+ Structure 1 written to file: ",A)') trim(file1)
       call c1%write_res(file1,-1)

       file2 = fileroot // "_structure_2.res"
       write (uout,'("+ Structure 2 written to file: ",A)') trim(file2)
       call c2del%write_res(file2,-1)
    end if

  end subroutine trick_compare_deformed

end module tricks
