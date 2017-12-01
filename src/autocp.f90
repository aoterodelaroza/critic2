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

!> Tools for automatic or manual critical point determination.
module autocp
  implicit none

  private 

  public :: autocritic
  public :: cpreport
  private :: critshell
  private :: writechk
  private :: readchk
  private :: atomic_connect_report
  private :: cp_short_report
  private :: barycentric
  private :: barycentric_divide
  private :: seed_from_simplex
  private :: cp_long_report
  private :: cp_vlong_report
  private :: graph_short_report
  private :: makegraph
  private :: scale_ws

  ! private for barycentric, initialized at the beginning of auto
  integer :: nstack
  integer, allocatable :: barstack(:,:,:), depstack(:)
  integer, parameter :: bardenom(2:4) = (/ 2, 6, 12 /)

  ! number of discarded degenerate critical points
  integer :: ndegenr

  ! clip the CPs and the seeds
  integer :: iclip = 0
  real*8 :: x0clip(3), x1clip(3), rclip

  ! variables that control the CP search
  integer :: dograph = 1 !< attempt build the topological graph after CP search. 

contains

  !> Automatic search for all the critical points in the crystal unit cell.
  !> Uses the IWS, with barycentric subdivision.
  subroutine autocritic(line)
    use systemmod, only: sy
    use fieldmod, only: type_grid
    use graphics, only: grhandle
    use surface, only: minisurf
    use global, only: quiet, cp_hdegen, eval_next, dunit0, iunit, iunitname0, fileroot
    use tools, only: uniqc
    use tools_math, only: norm
    use tools_io, only: uout, ferror, faterr, lgetword, equal, isexpression_or_word,&
       string, warning, tictac
    use types, only: realloc
    use param, only: pi
    character*(*), intent(in) :: line

    integer, parameter :: styp_ws = 1      ! recursive subdivision of the WS cell
    integer, parameter :: styp_pair = 2    ! pairs
    integer, parameter :: styp_triplet = 3 ! triplets
    integer, parameter :: styp_line = 4    ! line
    integer, parameter :: styp_sphere = 5  ! sphere
    integer, parameter :: styp_oh = 6  ! octahedron subdivision
    integer, parameter :: styp_point = 7   ! point
    type seed_
       integer :: typ        ! type of seeding strategy
       integer :: depth = 1  ! WS recursive subdivision level
       real*8 :: x0(3) = 0d0 ! origin of the search (crystallographic)
       real*8 :: rad = -1d0  ! radius of the search
       real*8 :: dist = 15d0 ! max. distance between atom pairs
       real*8 :: x1(3) = 0d0 ! end point for LINE
       integer :: npts = 1   ! number of points
       integer :: nr = 0     ! number of radial points
       integer :: ntheta = 0 ! number of theta (polar) points
       integer :: nphi = 0   ! number of phi (azimuthal) points
       integer :: nseed = 0  ! number of seeds generated
    end type seed_
    real*8, parameter :: seed_eps = 1d-5
    integer, parameter :: maxpointstr(0:7) =  (/ 6, 18,  66, 258, 1026, 4098, 16386, 66003  /)
    integer, parameter :: maxfacestr(0:7) =   (/ 8, 32, 128, 512, 2048, 8192, 32768, 131072 /)

    real*8  :: iniv(4,3), xdum(3)
    integer :: nt
    logical :: ok, dryrun, dochk
    character(len=:), allocatable :: word, str, discexpr
    real*8 :: gfnormeps
    integer :: ntetrag
    real*8, allocatable :: tetrag(:,:,:)
    logical :: cpdebug
    integer :: nseed 
    integer :: i, j, k, i1, i2, i3, lp, lpo, n0
    integer :: m, nf, ntheta, nphi, nr
    type(seed_), allocatable :: seed(:), saux(:)
    logical :: firstseed, hadx1
    integer :: nn, ilag, nphiact, ier, nss
    real*8 :: x0(3), x1(3), dist, r, theta, phi, delta_phi, delta_theta, x(3)
    real*8, allocatable :: xseed(:,:)
    logical, allocatable :: keep(:)
    type(minisurf) :: srf
    real*8 :: cpeps
    real*8 :: nuceps
    real*8 :: nucepsh
    type(grhandle) :: gr

    if (.not.quiet) then
       call tictac("Start AUTO")
       write (uout,*)
    end if

    ! defaults
    dochk = .false.
    gfnormeps = 1d-12
    cpdebug = .false.
    dryrun = .false.
    if (.not.sy%c%ismolecule) then
       nseed = 1
       allocate(seed(1))
       seed(1)%typ = styp_ws
    else
       nseed = 1
       allocate(seed(1))
       seed(1)%typ = styp_pair
    endif
    firstseed = .true.
    hadx1 = .false.
    iclip = 0
    CP_hdegen = 1d-8
    cpeps = 1d-2
    if (sy%f(sy%iref)%type == type_grid) then
       nuceps = 2d0 * maxval(sy%c%aa / sy%f(sy%iref)%grid%n)
       nucepsh = 2d0 * maxval(sy%c%aa / sy%f(sy%iref)%grid%n)
    else
       nuceps = 1d-1
       nucepsh = 2d-1
    end if
    dograph = 1

    ! parse the input
    lp = 1
    discexpr = ""
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'dry')) then
          dryrun = .true.
       elseif (equal(word,'verbose')) then
          cpdebug = .true.
       elseif (equal(word,'gradeps')) then
          ok = eval_next(gfnormeps,line,lp)
          if (.not.ok) then
             call ferror('autocritic','bad AUTO/GRADEPS syntax',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'cpeps')) then
          ok = eval_next(cpeps,line,lp)
          if (.not.ok) then
             call ferror('autocritic','bad AUTO/CPEPS syntax',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'epsdegen')) then
          ok = eval_next(CP_hdegen,line,lp)
          if (.not.ok) then
             call ferror('autocritic','bad AUTO/CPEPS syntax',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'nuceps')) then
          ok = eval_next(nuceps,line,lp)
          if (.not.ok) then
             call ferror('autocritic','bad AUTO/NUCEPS syntax',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'nucepsh')) then
          ok = eval_next(nucepsh,line,lp)
          if (.not.ok) then
             call ferror('autocritic','bad AUTO/NUCEPS syntax',faterr,line,syntax=.true.)
             return
          end if
       else if (equal(word,'chk')) then
          dochk = .true.
       else if (equal(word,'clip')) then
          word = lgetword(line,lp)
          if (equal(word,'cube')) then
             iclip = 1
             ok = eval_next(x0clip(1),line,lp)
             ok = ok .and. eval_next(x0clip(2),line,lp)
             ok = ok .and. eval_next(x0clip(3),line,lp)
             ok = ok .and. eval_next(x1clip(1),line,lp)
             ok = ok .and. eval_next(x1clip(2),line,lp)
             ok = ok .and. eval_next(x1clip(3),line,lp)
             if (.not. ok) then
                call ferror('critic','Wrong AUTO/CLIP/CUBE coordinates.',faterr,line,syntax=.true.)
                return
             end if
             if (sy%c%ismolecule) then
                x0clip = sy%c%c2x(x0clip / dunit0(iunit) - sy%c%molx0)
                x1clip = sy%c%c2x(x1clip / dunit0(iunit) - sy%c%molx0)
             endif
          elseif (equal(word,'sphere')) then
             iclip = 2
             ok = eval_next(x0clip(1),line,lp)
             ok = ok .and. eval_next(x0clip(2),line,lp)
             ok = ok .and. eval_next(x0clip(3),line,lp)
             ok = ok .and. eval_next(rclip,line,lp)
             if (.not. ok) then
                call ferror('critic','Wrong AUTO/CLIP/SPHERE values.',faterr,line,syntax=.true.)
                return
             end if
             if (sy%c%ismolecule) x0clip = sy%c%c2x(x0clip / dunit0(iunit) - sy%c%molx0)
             rclip = rclip / dunit0(iunit)
          else
             call ferror('critic','Wrong AUTO/CLIP option.',faterr,line,syntax=.true.)
             return
          end if

       else if (equal(word,'seed')) then
          if (firstseed) then
             nseed = 0
             firstseed = .false.
          end if
          nseed = nseed + 1
          if (nseed > size(seed)) then
             allocate(saux(nseed*2))
             saux(1:nseed-1) = seed(1:nseed-1)
             call move_alloc(saux,seed)
          end if
          word = lgetword(line,lp)
          if (equal(word,'ws')) then
             seed(nseed)%typ = styp_ws
          elseif (equal(word,'pair')) then
             seed(nseed)%typ = styp_pair
          elseif (equal(word,'triplet')) then
             seed(nseed)%typ = styp_triplet
          elseif (equal(word,'line')) then
             seed(nseed)%typ = styp_line
          elseif (equal(word,'sphere')) then
             seed(nseed)%typ = styp_sphere
          elseif (equal(word,'oh')) then
             seed(nseed)%typ = styp_oh
          elseif (equal(word,'point')) then
             seed(nseed)%typ = styp_point
          else
             call ferror('autocritic','Unknown keyword in AUTO/SEED',faterr,line,syntax=.true.)
             return
          endif
          ! convert the default positions to the center of the cell if this is a molecule
          if (sy%c%ismolecule) then
             seed(nseed)%x0 = sy%c%c2x(sy%c%x2c(seed(nseed)%x0) - sy%c%molx0)
             seed(nseed)%x1 = sy%c%c2x(sy%c%x2c(seed(nseed)%x1) - sy%c%molx0)
          endif
          do while (.true.)
             lpo = lp
             word = lgetword(line,lp)
             if (equal(word,'depth')) then
                ok = eval_next(seed(nseed)%depth,line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/DEPTH',faterr,line,syntax=.true.)
                   return
                end if
             elseif (equal(word,'x0')) then
                ok = eval_next(seed(nseed)%x0(1),line,lp)
                ok = ok.and.eval_next(seed(nseed)%x0(2),line,lp)
                ok = ok.and.eval_next(seed(nseed)%x0(3),line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/X0',faterr,line,syntax=.true.)
                   return
                end if
                if (sy%c%ismolecule) seed(nseed)%x0 = sy%c%c2x(seed(nseed)%x0 / dunit0(iunit) - sy%c%molx0)
             elseif (equal(word,'x1')) then
                ok = eval_next(seed(nseed)%x1(1),line,lp)
                ok = ok.and.eval_next(seed(nseed)%x1(2),line,lp)
                ok = ok.and.eval_next(seed(nseed)%x1(3),line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/X1',faterr,line,syntax=.true.)
                   return
                end if
                if (sy%c%ismolecule) seed(nseed)%x1 = sy%c%c2x(seed(nseed)%x1 / dunit0(iunit) - sy%c%molx0)
                hadx1 = .true.
             elseif (equal(word,'npts')) then
                ok = eval_next(seed(nseed)%npts,line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/NPTS',faterr,line,syntax=.true.)
                   return
                end if
             elseif (equal(word,'ntheta')) then
                ok = eval_next(seed(nseed)%ntheta,line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/NTHETA',faterr,line,syntax=.true.)
                   return
                end if
             elseif (equal(word,'nphi')) then
                ok = eval_next(seed(nseed)%nphi,line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/NPHI',faterr,line,syntax=.true.)
                   return
                end if
             elseif (equal(word,'nr')) then
                ok = eval_next(seed(nseed)%nr,line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/NR',faterr,line,syntax=.true.)
                   return
                end if
             elseif (equal(word,'dist')) then
                ok = eval_next(seed(nseed)%dist,line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/NPTS',faterr,line,syntax=.true.)
                   return
                end if
                seed(nseed)%dist = seed(nseed)%dist / dunit0(iunit)
             elseif (equal(word,'radius')) then
                ok = eval_next(seed(nseed)%rad,line,lp)
                if (.not.ok) then
                   call ferror('autocritic','Wrong AUTO/SEED/NPTS',faterr,line,syntax=.true.)
                   return
                end if
                seed(nseed)%rad = seed(nseed)%rad / dunit0(iunit)
             else
                lp = lpo
                exit
             endif
          end do
       elseif (equal(word,"discard")) then
          ok = isexpression_or_word(discexpr,line,lp)
          if (.not. ok) then
             call ferror("autocritic","wrong DISCARD keyword",faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror('autocritic','Unknown keyword in AUTO',faterr,line,syntax=.true.)
          return
       else
          exit
       endif
    end do

    ! build the seed list
    nn = 0
    allocate(xseed(3,10))
    do i = 1, nseed
       n0 = nn
       seed(i)%nseed = 0
       if (seed(i)%typ == styp_ws) then
          ! Determine the WS cell
          call sy%c%wigner(seed(i)%x0,ntetrag=ntetrag,tetrag=tetrag)
          if (seed(i)%rad > 0) call scale_ws(seed(i)%rad,seed(i)%x0,ntetrag,tetrag)
          ! Recursively subdivide each tetrahedron
          do nt = 1, ntetrag     
             do j = 1, 4
                xdum = tetrag(j,:,nt)
                iniv(j,:) = sy%c%x2c(xdum)
             end do
             call barycentric(iniv,seed(i)%depth,nn,xseed)
          enddo
          if (allocated(tetrag)) deallocate(tetrag)

       elseif (seed(i)%typ == styp_pair) then
          ! between all pairs of atoms
          do i1 = 1, sy%c%ncel
             do i2 = 1, sy%c%ncel
                if (i1 == i2) cycle
                dist = norm(sy%c%atcel(i1)%r - sy%c%atcel(i2)%r)
                if (dist > seed(i)%dist) cycle
                do k = 1, seed(i)%npts
                   nn = nn + 1
                   if (nn > size(xseed,2)) call realloc(xseed,3,2*nn)
                   xseed(:,nn) = sy%c%atcel(i1)%x + real(k,8) / real(seed(i)%npts+1,8) * &
                      (sy%c%atcel(i2)%x - sy%c%atcel(i1)%x)
                end do
             end do
          end do
       elseif (seed(i)%typ == styp_triplet) then
          ! between all atom triplets
          do i1 = 1, sy%c%ncel
             do i2 = 1, sy%c%ncel
                if (i1 == i2) cycle
                dist = norm(sy%c%atcel(i1)%r - sy%c%atcel(i2)%r)
                if (dist > seed(i)%dist) cycle
                do i3 = 1, sy%c%ncel
                   if (i1 == i3 .or. i2 == i3) cycle
                   dist = norm(sy%c%atcel(i1)%r - sy%c%atcel(i3)%r)
                   if (dist > seed(i)%dist) cycle
                   dist = norm(sy%c%atcel(i2)%r - sy%c%atcel(i3)%r)
                   if (dist > seed(i)%dist) cycle

                   nn = nn + 1
                   if (nn > size(xseed,2)) call realloc(xseed,3,2*nn)
                   xseed(:,nn) = 1d0/3d0 * (sy%c%atcel(i1)%x + sy%c%atcel(i2)%x + sy%c%atcel(i3)%x)
                end do
             end do
          end do
       elseif (seed(i)%typ == styp_line) then
          if (.not.hadx1) then
             call ferror('autocritic','SEED/LINE requires X1',faterr,syntax=.true.)
             return
          end if
          ! add a line of points
          do j = 1, seed(i)%npts
             nn = nn + 1
             if (nn > size(xseed,2)) call realloc(xseed,3,2*nn)
             xseed(:,nn) = seed(i)%x0 + real(j-1,8) / real(seed(i)%npts-1,8) * &
                (seed(i)%x1 - seed(i)%x0)
          end do
       elseif (seed(i)%typ == styp_sphere) then
          if (seed(i)%nr == 0.or.seed(i)%ntheta == 0.or.seed(i)%nphi == 0) then
             call ferror('autocritic','SEED/SPHERE requires NR, NPHI, and NTHETA',faterr,syntax=.true.)
             return
          end if
          if (seed(i)%rad < 0) then
             call ferror('autocritic','SEED/SPHERE requires a radius',faterr,syntax=.true.)
             return
          end if
          ! inside a sphere
          nr = seed(i)%nr
          nphi = seed(i)%nphi
          ntheta = seed(i)%ntheta
          x1 = sy%c%x2c(seed(i)%x0)
          delta_theta = pi/2d0/ntheta
          theta = delta_theta
          nphiact = nphi
          do i1 = 1, ntheta
             delta_phi = 2d0*pi/nphiact
             phi = 0d0
             do i2 = 1, nphiact
                do i3 = 1, nr
                   r = seed(i)%rad * real(i3,8) / real(nr,8)

                   nn = nn + 1
                   if (nn > size(xseed,2)) call realloc(xseed,3,2*nn)
                   xseed(:,nn) = x1 + r * &
                      (/ sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) /)
                   xseed(:,nn) = sy%c%c2x(xseed(:,nn))

                   nn = nn + 1
                   if (nn > size(xseed,2)) call realloc(xseed,3,2*nn)
                   xseed(:,nn) = x1 + r * &
                      (/ sin(pi-theta)*cos(phi), sin(pi-theta)*sin(phi), cos(pi-theta) /)
                   xseed(:,nn) = sy%c%c2x(xseed(:,nn))
                end do
                phi=phi+delta_phi
             end do
             theta = theta + delta_theta
             nphiact = nphiact + nphiact
          end do
       elseif (seed(i)%typ == styp_oh) then
          if (seed(i)%nr == 0) then
             call ferror('autocritic','SEED/OH requires NR',faterr,syntax=.true.)
             return
          end if
          if (seed(i)%rad < 0) then
             call ferror('autocritic','SEED/OH requires a radius',faterr,syntax=.true.)
             return
          end if
          if (seed(i)%depth > 7) then
             call ferror('autocritic','The depth in SEED/OH can not be > 7',faterr,syntax=.true.)
             return
          end if
          ! recursive subdivision of an octahedron
          nr = seed(i)%nr
          m = maxpointstr(seed(i)%depth) + 1
          nf = maxfacestr(seed(i)%depth) + 1
          x0 = sy%c%x2c(seed(i)%x0)
          call srf%init(m,nf)
          call srf%clean()
          call srf%spheretriang(x0,seed(i)%depth)

          ! assign the nodes
          do j = 1, srf%nv
             do k = 1, nr
                r = seed(i)%rad * real(k,8) / real(nr,8)

                nn = nn + 1
                if (nn > size(xseed,2)) call realloc(xseed,3,2*nn)
                x0 = srf%n + r * &
                   (/ srf%r(j) * sin(srf%th(j)) * cos(srf%ph(j)),&
                   srf%r(j) * sin(srf%th(j)) * sin(srf%ph(j)),&
                   srf%r(j) * cos(srf%th(j)) /)
                xseed(:,nn) = sy%c%c2x(x0)
             end do
          end do

          ! clean up
          call srf%end()
       elseif (seed(i)%typ == styp_point) then
          ! add a point
          nn = nn + 1
          if (nn > size(xseed,2)) call realloc(xseed,3,2*nn)
          xseed(:,nn) = seed(i)%x0
       endif
       seed(i)%nseed = nn - n0
    end do
    call realloc(xseed,3,nn)

    ! write the header to the output 
    write (uout,'("* Automatic determination of CPs")')
    write (uout,'("  Discard new CPs if another CP was found at a distance less than: ",A,X,A)') &
       string(cpeps*dunit0(iunit),'e',decimal=3), iunitname0(iunit)
    write (uout,'("  Discard new CPs if a nucleus was found at a distance less than: ",A,X,A)') &
       string(nuceps*dunit0(iunit),'e',decimal=3), iunitname0(iunit)
    write (uout,'("  Discard new CPs if a hydrogen was found at a distance less than: ",A,X,A)') &
       string(nucepsh*dunit0(iunit),'e',decimal=3), iunitname0(iunit)
    write (uout,'("  CPs are degenerate if any Hessian element abs value is less than: ",A)') &
       string(CP_hdegen,'e',decimal=3)
    if (len_trim(discexpr) > 0) &
       write (uout,'("  Discard CP expression: ",A)') trim(discexpr)
    write (uout,'("  Discard CPs if grad(f) is above: ",A)') string(gfnormeps,'e',decimal=3)
    write (uout,'("+ List of seeding actions")')
    write (uout,'("  Id nseed     Type           Description")')
    do i = 1, nseed
       x0 = seed(i)%x0
       x1 = seed(i)%x1
       r = seed(i)%rad * dunit0(iunit)
       dist = seed(i)%dist * dunit0(iunit)
       if (sy%c%ismolecule) then
          x0 = (sy%c%x2c(x0) + sy%c%molx0) * dunit0(iunit)
          x1 = (sy%c%x2c(x1) + sy%c%molx0) * dunit0(iunit)
       endif
       str = "  " // string(i,2)
       str = str // string(seed(i)%nseed,7)
       if (seed(i)%typ == styp_ws) then
          str = str // " WS recursive  "
          str = str // "  depth=" // string(seed(i)%depth)
          str = trim(str) // ", x0=" // string(x0(1),'f',7,4) // " " // &
             string(x0(2),'f',7,4) // " " // string(x0(3),'f',7,4)
          if (r > 0d0) then
             str = trim(str) // ", radius=" // trim(string(r,'f',10,4))
          else
             str = trim(str) // ", no radius "
          endif
       elseif (seed(i)%typ == styp_pair) then
          str = str // " Atom pairs    "
          str = str // "  dist=" // trim(string(dist,'f',10,4))
          str = str // ", npts=" // string(seed(i)%npts)
       elseif (seed(i)%typ == styp_triplet) then
          str = str // " Atom triplets "
          str = str // "  dist=" // trim(string(dist,'f',10,4))
       elseif (seed(i)%typ == styp_line) then
          str = str // " Line          "
          str = str // "  x0=" // string(x0(1),'f',7,4) // " " // &
             string(x0(2),'f',7,4) // " " // string(x0(3),'f',7,4)
          str = trim(str) // ", x1=" // string(x1(1),'f',7,4) // " " // &
             string(x1(2),'f',7,4) // " " // string(x1(3),'f',7,4)
          str = str // ", npts=" // string(seed(i)%npts)
       elseif (seed(i)%typ == styp_sphere) then
          str = str // " Sphere        "
          str = str // "  x0=" // string(x0(1),'f',7,4) // " " // &
             string(x0(2),'f',7,4) // " " // string(x0(3),'f',7,4)
          str = trim(str) // ", radius=" // trim(string(r,'f',10,4))
          str = str // ", ntheta=" // string(seed(i)%ntheta)
          str = str // ", nphi=" // string(seed(i)%nphi)
          str = str // ", nr=" // string(seed(i)%nr)
       elseif (seed(i)%typ == styp_oh) then
          str = str // " Oh recursive  "
          str = str // "  depth=" // string(seed(i)%depth)
          str = trim(str) // ", x0=" // string(x0(1),'f',7,4) // " " // &
             string(x0(2),'f',7,4) // " " // string(x0(3),'f',7,4)
          str = trim(str) // ", radius=" // trim(string(r,'f',10,4))
          str = str // ", nr=" // string(seed(i)%nr)
       elseif (seed(i)%typ == styp_point) then
          str = str // " Point         "
          str = str // "  x0=" // string(x0(1),'f',7,4) // " " // &
             string(x0(2),'f',7,4) // " " // string(x0(3),'f',7,4)
       endif
       write (uout,'(A)') str
    end do
    write (uout,'("+ Number of seeds before pruning and clipping: ",A)') string(nn)

    ! move all the seeds to the main cell
    do i = 1, nn
       xseed(:,i) = xseed(:,i) - floor(xseed(:,i))
    end do

    ! eliminate the seeds outside the molecular cell
    if (sy%c%ismolecule) then
       allocate(keep(nn))
       keep = .true.
       do i = 1, nn
          x0 = xseed(:,i)
          if (x0(1) < sy%c%molborder(1) .or. x0(1) > (1d0-sy%c%molborder(1)) .or.&
             x0(2) < sy%c%molborder(2) .or. x0(2) > (1d0-sy%c%molborder(2)) .or.&
             x0(3) < sy%c%molborder(3) .or. x0(3) > (1d0-sy%c%molborder(3))) then
             keep(i) = .false.
          endif
       end do
       ilag = 0
       do i = 1, nn
          if (keep(i)) then
             ilag = ilag + 1
             xseed(:,ilag) = xseed(:,i)
          endif
       end do
       nn = ilag
       call realloc(xseed,3,nn)
    endif

    ! clip the cube or the sphere
    if (iclip > 0) then
       if (.not.allocated(keep)) then
          allocate(keep(nn))
          keep = .true.
       end if
       if (iclip == 1) then
          ! cube
          do i = 1, nn
             if (xseed(1,i) < x0clip(1).or.xseed(1,i) > x1clip(1).or.&
                xseed(2,i) < x0clip(2).or.xseed(2,i) > x1clip(2).or.&
                xseed(3,i) < x0clip(3).or.xseed(3,i) > x1clip(3)) &
                keep(i) = .false.
          end do
       else
          ! sphere
          do i = 1, nn
             if (.not.sy%c%are_lclose(xseed(:,i),x0clip,rclip)) &
                keep(i) = .false.
          end do
       end if
       ilag = 0
       do i = 1, nn
          if (keep(i)) then
             ilag = ilag + 1
             xseed(:,ilag) = xseed(:,i)
          endif
       end do
       nn = ilag
       call realloc(xseed,3,nn)
       deallocate(keep)
    end if

    ! transform to cartesian
    do i = 1, nn
       xseed(:,i) = sy%c%x2c(xseed(:,i))
    end do

    ! uniq the list
    call uniqc(xseed,1,nn,seed_eps)

    ! this is the final list of seeds
    write (uout,'("+ Number of seeds: ",A)') string(nn)

    ! write the seeds to an obj file
    if (cpdebug) then
       str = trim(fileroot) // "_seeds.obj" 
       write (uout,'("+ Writing seeds to file: ",A)') str
       call sy%c%write_3dmodel(str,"obj",(/1,1,1/),.true.,.false.,.false.,&
          .true.,.true.,-1d0,(/0d0,0d0,0d0/),-1d0,(/0d0,0d0,0d0/),gr)
       do i = 1, nn
          call gr%ball(xseed(:,i) + sy%c%molx0,(/100,100,255/),0.3d0)
       end do
       call gr%close()
    endif

    ! Initialize the CP search
    if (.not.allocated(sy%f(sy%iref)%cp)) call sy%f(sy%iref)%init_cplist

    ! Read cps from external file
    if (dochk) call readchk()

    ! If not a dry run, do the search
    if (.not.dryrun) then
       write (uout,'("+ Searching for CPs")')
       if (cpdebug) &
          write (uout,'("  CP localization progress:")')
       nss = max(nn / 25,1)
       ndegenr = 0
       !$omp parallel do private(ier,x0) schedule(dynamic)
       do i = 1, nn
          if (cpdebug) then
             !$omp critical (progress)
             if (mod(i,nss) == 1) then
                write (uout,'("  [",A,"/",A,"]")') string(i), string(nn)
             end if
             !$omp end critical (progress)
          end if
          x0 = xseed(:,i)
          call sy%f(sy%iref)%newton(x0,gfnormeps,ier)
          if (ier <= 0) then
             ! Check if it's inside the sphere
             ok = .true.
             if (iclip > 0) then
                x = sy%c%c2x(x0)
                if (iclip == 1) then
                   ! cube
                   ok = .not.(x(1) < x0clip(1).or.x(1) > x1clip(1).or.&
                      x(2) < x0clip(2).or.x(2) > x1clip(2).or.&
                      x(3) < x0clip(3).or.x(3) > x1clip(3))
                else
                   ! sphere
                   x = x - x0clip
                   x = sy%c%x2c(x)
                   ok = .not.(norm(x) >= rclip) 
                end if
             end if

             if (ok) then
                !$omp critical (addcp)
                call sy%addcp(sy%iref,x0,discexpr,cpeps,nuceps,nucepsh)
                !$omp end critical (addcp)
             end if
          end if
       end do
       !$omp end parallel do

       if (ndegenr > 0) then
          call ferror('autocritic',string(ndegenr) // " degenerate critical points discarded.",warning)
          write (uout,'("This is normal if your system has a vacuum (i.e. a molecule or a slab)")')
          write (uout,'("Consider: ")')
          write (uout,'("  - changing the threshold for degenerate critical points (EPSDEGEN)")')
          write (uout,'("  - checking the sanity of the scalar field.")')
          write (uout,'("  - using the DISCARD keyword to delete vacuum critical points.")')
          if (sy%c%ismolecule) &
             write (uout,'("  - reducing the size of the molecular cell")')
          write (uout,*)
       end if

       ! Sort the cp list, using the value of the reference field
       write (uout,'("+ Sorting the CPs")')
       call sy%f(sy%iref)%sortcps(cpeps)
    end if
    write (uout,*)

    ! Short report of non-equivalent cp-list
    call cp_short_report()

    ! Graph
    if (dograph > 0) then
       if (.not.dryrun) then
          call makegraph()
       end if
       call graph_short_report()
    end if

    ! Long report of all the CPs in the unit cell
    call cp_long_report()

    ! Attractor connectivity matrix
    if (dograph > 0) call atomic_connect_report()

    !.Particular information for every critical point found
    call cp_vlong_report()

    ! Write CPs to an external file
    if (dochk) call writechk()

    ! reallocate to the new CP list and clean up
    call realloc(sy%f(sy%iref)%cp,sy%f(sy%iref)%ncp)
    call realloc(sy%f(sy%iref)%cpcel,sy%f(sy%iref)%ncpcel)
    iclip = 0

    if (.not.quiet) then
       call tictac("End AUTO")
       write (uout,*)
    end if

  end subroutine autocritic

  !> Report the results of the critical point search.
  subroutine cpreport(line)
    use systemmod, only: sy, system
    use struct_drivers, only: struct_write
    use crystalseedmod, only: crystalseed
    use global, only: eval_next, prunedist, gcpchange
    use tools_io, only: lgetword, equal, getword, ferror, faterr, nameguess
    use types, only: realloc, gpathp
    character*(*), intent(in) :: line

    integer :: lp, n, lp2
    integer :: i, iup, nstep, ier
    character(len=:), allocatable :: word
    character(len=:), allocatable :: line2, aux
    logical :: ok
    logical :: agraph
    type(crystalseed) :: seed
    type(system) :: syaux
    real*8 :: x(3), plen
    type(gpathp), allocatable :: xpath(:)

    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'short')) then
          call cp_short_report()
       elseif (equal(word,'long')) then
          if (.not.sy%c%ismolecule) then
             call cp_long_report()
          else
             call cp_short_report()
          end if
       elseif (equal(word,'verylong')) then
          call cp_vlong_report()
       elseif (equal(word,'shells')) then
          ok = eval_next(n,line,lp)
          if (.not.ok) n = 10
          call critshell(n)
       elseif (len_trim(word) > 0) then
          ! parse for special keywords and build the line for write
          lp2 = 1
          line2 = ""
          agraph = .false.
          do while (.true.)
             word = getword(line,lp2)
             if (equal(word,'graph')) then
                agraph = .true.
             elseif (len_trim(word) > 0) then
                aux = line2 // " " // trim(word)
                line2 = aux
             else
                exit
             end if
          end do

          ! build the crystal structure containing the crystal points
          seed%isused = .true.
          seed%file = sy%c%file
          seed%nat = sy%f(sy%iref)%ncpcel
          allocate(seed%x(3,sy%f(sy%iref)%ncpcel))
          allocate(seed%z(sy%f(sy%iref)%ncpcel))
          allocate(seed%name(sy%f(sy%iref)%ncpcel))
          do i = 1, sy%f(sy%iref)%ncpcel
             seed%x(:,i) = sy%f(sy%iref)%cpcel(i)%x
             if (i <= sy%c%ncel) then
                seed%z(i) = sy%c%at(sy%c%atcel(i)%idx)%z
                seed%name(i) = sy%c%at(sy%c%atcel(i)%idx)%name
             else
                if (sy%f(sy%iref)%cp(sy%f(sy%iref)%cpcel(i)%idx)%typ == -3) then
                   seed%z(i) = 119
                elseif (sy%f(sy%iref)%cp(sy%f(sy%iref)%cpcel(i)%idx)%typ == -1) then
                   seed%z(i) = 120
                elseif (sy%f(sy%iref)%cp(sy%f(sy%iref)%cpcel(i)%idx)%typ == 1) then
                   seed%z(i) = 121
                elseif (sy%f(sy%iref)%cp(sy%f(sy%iref)%cpcel(i)%idx)%typ == 3) then
                   seed%z(i) = 122
                else
                   call ferror("cpreport","unclassified critical point",faterr)
                end if
                seed%name(i) = nameguess(seed%z(i))
             end if
          end do
          seed%usezname = 3
          seed%useabr = 2
          seed%crys2car = sy%c%crys2car
          seed%havesym = 0
          seed%findsym = 0

          ! calculate gradient paths and add them to the seed
          if (agraph) then
             if (allocated(xpath)) deallocate(xpath)
             allocate(xpath(1))
             !$omp parallel do private(iup,x,nstep,ier,plen) firstprivate(xpath) schedule(dynamic)
             do i = sy%c%ncel+1, sy%f(sy%iref)%ncpcel

                if (sy%f(sy%iref)%typnuc == -3 .and. sy%f(sy%iref)%cp(sy%f(sy%iref)%cpcel(i)%idx)%typ == -1) then
                   iup = 1
                else if (sy%f(sy%iref)%typnuc == 3 .and. sy%f(sy%iref)%cp(sy%f(sy%iref)%cpcel(i)%idx)%typ == 1) then
                   iup = -1
                else
                   iup = 0
                end if

                if (iup /= 0) then
                   x = sy%f(sy%iref)%cpcel(i)%r + 0.5d0 * gcpchange * sy%f(sy%iref)%cpcel(i)%brvec
                   call sy%f(sy%iref)%gradient(x,iup,nstep,ier,.false.,plen,xpath,prunedist,pathini=sy%f(sy%iref)%cpcel(i)%r)
                   !$omp critical (add)
                   call addpath(size(xpath,1),xpath)
                   !$omp end critical (add)
                   x = sy%f(sy%iref)%cpcel(i)%r - 0.5d0 * gcpchange * sy%f(sy%iref)%cpcel(i)%brvec
                   call sy%f(sy%iref)%gradient(x,iup,nstep,ier,.false.,plen,xpath,prunedist,pathini=sy%f(sy%iref)%cpcel(i)%r)
                   !$omp critical (add)
                   call addpath(size(xpath,1),xpath)
                   !$omp end critical (add)
                end if
             end do
             !$omp end parallel do
          end if

          ! molecule information
          seed%ismolecule = sy%c%ismolecule
          seed%cubic = .false.
          seed%border = 0d0
          seed%havex0 = .true.
          seed%molx0 = sy%c%molx0

          ! write the structure to the external file
          call syaux%init()
          call syaux%c%struct_new(seed,.true.)
          call struct_write(syaux,line2)
          call syaux%end()
          return
       else
          exit
       endif
    end do

  contains
    subroutine addpath(nstep,xpath)
      use types, only: realloc
      integer, intent(in) :: nstep
      type(gpathp), intent(in) :: xpath(nstep)
      integer :: i, n

      call realloc(seed%x,3,seed%nat+nstep)
      call realloc(seed%z,seed%nat+nstep)
      call realloc(seed%name,seed%nat+nstep)
      n = seed%nat
      do i = 1, nstep
         n = n + 1
         seed%x(:,n) = xpath(i)%x
         seed%z(n) = 123
         seed%name(n) = "Xz"
      end do
      seed%nat = n

    end subroutine addpath
  end subroutine cpreport

  !> Calculates the neighbor environment of each non-equivalent CP.
  subroutine critshell(shmax)
    use systemmod, only: sy
    use global, only: iunit, iunitname0, dunit0
    use tools_io, only: uout, string, ioj_center
    integer, intent(in) :: shmax

    integer :: i, j, k, l, ln, lvec(3), ilo, ihi
    real*8 :: x0(3), xq(3)
    real*8 :: dist2(sy%f(sy%iref)%ncp,shmax), d2, dmin
    integer :: nneig(sy%f(sy%iref)%ncp,shmax), imin
    integer :: wcp(sy%f(sy%iref)%ncp,shmax)

    character*1, parameter :: namecrit(0:4) = &
       (/'n','b','r','c','?'/)

    dist2 = 1d30
    dmin = 1d30
    nneig = 0
    wcp = 0
    if (.not.sy%c%ismolecule) then
       ilo = 0
       ihi = 26
    else
       ilo = (1 * 3 + 1) * 3 + 1
       ihi = ilo
    end if
    do i = 1, sy%f(sy%iref)%ncp
       x0 = sy%c%x2c(sy%f(sy%iref)%cp(i)%x)
       do j = 1, sy%f(sy%iref)%ncpcel
          do k = ilo, ihi
             ln = k
             lvec(1) = mod(ln,3) - 1
             ln = ln / 3
             lvec(2) = mod(ln,3) - 1
             ln = ln / 3
             lvec(3) = mod(ln,3) - 1
             xq = sy%c%x2c(sy%f(sy%iref)%cpcel(j)%x + lvec)
             d2 = dot_product(xq-x0,xq-x0)
             if (d2 < 1d-12) cycle
             do l = 1, shmax
                if (abs(d2 - dist2(i,l)) < 1d-12) then
                   nneig(i,l) = nneig(i,l) + 1
                   exit
                else if (d2 < dist2(i,l)) then 
                   dist2(i,l+1:shmax) = dist2(i,l:shmax-1)
                   nneig(i,l+1:shmax) = nneig(i,l:shmax-1)
                   wcp(i,l+1:shmax) = wcp(i,l:shmax-1)
                   dist2(i,l) = d2
                   nneig(i,l) = 1
                   wcp(i,l) = sy%f(sy%iref)%cpcel(j)%idx
                   exit
                end if
             end do
          end do
       end do
       if (dist2(i,1) < dmin) then
          dmin = dist2(i,1)
          imin = i
       end if
    end do

    write (uout,'("* Environments of the critical points")')
    write (uout,'("# ncp typ neig   dist(",A,")  ncp typ")') string(iunitname0(iunit))
    do i = 1, sy%f(sy%iref)%ncp
       do j = 1, shmax
          if (j == 1) then
             write (uout,'(6(A,X))') &
                string(i,length=6,justify=ioj_center), namecrit(sy%f(sy%iref)%cp(i)%typind),&
                string(nneig(i,j),length=5,justify=ioj_center), &
                string(sqrt(dist2(i,j))*dunit0(iunit),'f',length=12,decimal=8,justify=3), &
                string(wcp(i,j),length=4,justify=ioj_center), &
                namecrit(sy%f(sy%iref)%cp(wcp(i,j))%typind)
          else
             if (wcp(i,j) /= 0) then
                write (uout,'(6X,"...",6(A,X))') &
                   string(nneig(i,j),length=5,justify=ioj_center), &
                   string(sqrt(dist2(i,j))*dunit0(iunit),'f',length=12,decimal=8,justify=3), &
                   string(wcp(i,j),length=4,justify=ioj_center), &
                   namecrit(sy%f(sy%iref)%cp(wcp(i,j))%typind)
             end if
          end if
       end do
    end do
    write (uout,*)
    write (uout,'("* Minimum CP distance is ",A,X,A," between CP# ",A," and ",A)') &
       string(sqrt(dmin)*dunit0(iunit),'f',length=12,decimal=7,justify=ioj_center), &
       string(iunitname0(iunit)), string(imin), string(wcp(imin,1))
    write (uout,*)

  end subroutine critshell

  !> Write the CP information to the checkpoint file.
  subroutine writechk()
    use systemmod, only: sy
    use global, only: fileroot
    use tools_io, only: uout, fopen_write, fclose, string
    character(len=:), allocatable :: cpfile
    integer :: lucp, i

    cpfile = trim(fileroot) // ".chk_cps" 

    write (uout,'("* Writing CP file : ",A)') string(cpfile)
    write (uout,*)

    lucp = fopen_write(cpfile,"unformatted")
    write (lucp) sy%f(sy%iref)%ncp, sy%f(sy%iref)%ncpcel
    write (lucp) (sy%f(sy%iref)%cp(i),i=1,sy%f(sy%iref)%ncp)
    write (lucp) (sy%f(sy%iref)%cpcel(i),i=1,sy%f(sy%iref)%ncpcel)
    call fclose(lucp)

  end subroutine writechk

  !> Read the CP information from the checkpoint file.
  subroutine readchk()
    use systemmod, only: sy
    use global, only: fileroot
    use tools_io, only: uout, fopen_write, string, fclose
    use types, only: realloc
    integer :: lucp, i

    character(len=:), allocatable :: cpfile
    logical :: existcpfile

    cpfile = trim(fileroot) // ".chk_cps" 

    inquire(file=cpfile,exist=existcpfile)
    if (existcpfile) then
       write (uout,'("* Reading checkpoint file : ",A)') string(cpfile)
       write (uout,*) 

       lucp = fopen_write(cpfile,"unformatted")
       read (lucp) sy%f(sy%iref)%ncp, sy%f(sy%iref)%ncpcel

       if (.not.allocated(sy%f(sy%iref)%cp)) then
          allocate(sy%f(sy%iref)%cp(sy%f(sy%iref)%ncp))
       else if (size(sy%f(sy%iref)%cp) < sy%f(sy%iref)%ncp) then
          call realloc(sy%f(sy%iref)%cp,sy%f(sy%iref)%ncp)
       end if
       read (lucp) (sy%f(sy%iref)%cp(i),i=1,sy%f(sy%iref)%ncp)

       if (.not.allocated(sy%f(sy%iref)%cpcel)) then
          allocate(sy%f(sy%iref)%cpcel(sy%f(sy%iref)%ncpcel))
       else if (size(sy%f(sy%iref)%cpcel) < sy%f(sy%iref)%ncpcel) then
          call realloc(sy%f(sy%iref)%cpcel,sy%f(sy%iref)%ncpcel)
       end if
       read (lucp) (sy%f(sy%iref)%cpcel(i),i=1,sy%f(sy%iref)%ncpcel)

       call fclose(lucp)

    end if

  end subroutine readchk

  subroutine atomic_connect_report()
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_center
    integer :: i, j, k, i1, i2, m1, m2
    integer :: ind
    integer :: connectm(0:sy%f(sy%iref)%ncp,0:sy%f(sy%iref)%ncp), ncon, icon(0:sy%f(sy%iref)%ncp)
    character(len=1024) :: auxsmall(0:sy%f(sy%iref)%ncp)
    character*(1) ::    smallnamecrit(0:4)

    data   smallnamecrit   /'n','b','r','c','?'/

    connectm = 0
    do i = 1, sy%f(sy%iref)%ncp
       if (sy%f(sy%iref)%cp(i)%typ /= sign(1,sy%f(sy%iref)%typnuc)) cycle
       i1 = sy%f(sy%iref)%cp(i)%ipath(1)
       i2 = sy%f(sy%iref)%cp(i)%ipath(2)
       if (i1 <= 0 .or. i2 <= 0) cycle
       m1 = sy%f(sy%iref)%cp(i1)%mult
       m2 = sy%f(sy%iref)%cp(i2)%mult

       if (i1 == i2) then
          connectm(i1,i2) = connectm(i1,i2) + 2 * sy%f(sy%iref)%cp(i)%mult / m1
       else
          connectm(i1,i2) = connectm(i1,i2) + sy%f(sy%iref)%cp(i)%mult / m1
          connectm(i2,i1) = connectm(i2,i1) + sy%f(sy%iref)%cp(i)%mult / m2
       end if
    end do

    ! determine the any-nonzero rows / columns
    ncon = -1
    do i = 1, sy%f(sy%iref)%ncp
       if (any(connectm(:,i) /= 0) .or. any(connectm(i,:) /= 0)) then
          ncon = ncon + 1
          icon(ncon) = i
       end if
    end do

    if (ncon >= 0) then
       write (uout,'("* Attractor connectivity matrix")')
       write (uout,'("  to be read: each cp in a row i connects to mij cps in column j")')
       do j = 0, ncon
          ind = (sy%f(sy%iref)%cp(icon(j))%typ + 3) / 2
          write (auxsmall(j),'(a1,"(",A,")")') smallnamecrit(ind), string(icon(j))
       end do
       do i = 0, ncon / 6
          write (uout,'(13x,6(A))') (string(auxsmall(j),6,ioj_center), j = 6*i,min(6*(i+1)-1,ncon))
          write (uout,'(13x,6(A))') (string(sy%f(sy%iref)%cp(icon(j))%name,6,ioj_center), j = 6*i,min(6*(i+1)-1,ncon))
          do j = 0, ncon
             write (uout,'(A,x,A,x,6(A))') string(auxsmall(j),6,ioj_center), &
                string(sy%f(sy%iref)%cp(icon(j))%name,6,ioj_center), &
                (string(connectm(icon(j),icon(k)),6,ioj_center), k = 6*i, min(6*(i+1)-1,ncon))
          end do
       end do
       write (uout,*)
    end if

  end subroutine atomic_connect_report

  ! write the short report about the non-equivalent cps
  subroutine cp_short_report()
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_left, ioj_center, ioj_right
    use types, only: scalar_value
    use global, only: iunit, iunitname0, dunit0
    integer :: i, j
    integer :: numclass(0:3), multclass(0:3)
    character*(8) :: namecrit(0:6)
    integer :: ind

    data namecrit /'nucleus','bond','ring','cage','degener','nnattr','eff.nuc'/

    write (uout,'("* Critical point list, final report (non-equivalent cps)")')

    ! determine the topological class
    numclass = 0
    multclass = 0
    do i = 1, sy%f(sy%iref)%ncp
       numclass(sy%f(sy%iref)%cp(i)%typind) = numclass(sy%f(sy%iref)%cp(i)%typind) + 1
       multclass(sy%f(sy%iref)%cp(i)%typind) = multclass(sy%f(sy%iref)%cp(i)%typind) + sy%f(sy%iref)%cp(i)%mult
    end do
    write (uout,'("  Topological class (n|b|r|c): ",4(A,"(",A,") "))') &
       (string(numclass(i)),string(multclass(i)),i=0,3)
    if (sy%c%ismolecule) then
       write (uout,'("  Poincare-Hopf sum: ",A)') string(multclass(0)-multclass(1)+multclass(2)-multclass(3))
       write (uout,'("# ncp   type   CPname              position (",A,")               name            f             |grad|           lap")') &
          string(iunitname0(iunit))
    else
       write (uout,'("  Morse sum: ",A)') string(multclass(0)-multclass(1)+multclass(2)-multclass(3))
       write (uout,'("# ncp   pg  type   CPname         position (cryst. coords.)       mult  name            f             |grad|           lap")')
    endif

    ! report
    do i = 1, sy%f(sy%iref)%ncp
       ind = sy%f(sy%iref)%cp(i)%typind
       if (sy%f(sy%iref)%cp(i)%isdeg) then
          ind = 4
       else if (sy%f(sy%iref)%cp(i)%isnnm) then
          ind = 5
       end if

       if (.not.sy%c%ismolecule) then
          write (uout,'(2X,A,1x,A," (3,",A,") ",A,1X,3(A,1X),A,1x,A,3(1X,A))') &
             string(i,length=4,justify=ioj_left), string(sy%f(sy%iref)%cp(i)%pg,length=3,justify=ioj_center),&
             string(sy%f(sy%iref)%cp(i)%typ,length=2), string(namecrit(ind),length=8,justify=ioj_center),&
             (string(sy%f(sy%iref)%cp(i)%x(j),'f',length=12,decimal=8,justify=3),j=1,3), &
             string(sy%f(sy%iref)%cp(i)%mult,length=3,justify=ioj_right), &
             string(sy%f(sy%iref)%cp(i)%name,length=10,justify=ioj_center),&
             string(sy%f(sy%iref)%cp(i)%s%f,'e',decimal=8,length=15,justify=4),&
             string(sy%f(sy%iref)%cp(i)%s%gfmod,'e',decimal=8,length=15,justify=4),&
             string(sy%f(sy%iref)%cp(i)%s%del2f,'e',decimal=8,length=15,justify=4)
       else
          write (uout,'(2X,A," (3,",A,") ",A,1X,3(A,1X),A,3(1X,A))') &
             string(i,length=4,justify=ioj_left),&
             string(sy%f(sy%iref)%cp(i)%typ,length=2), string(namecrit(ind),length=8,justify=ioj_center),&
             (string((sy%f(sy%iref)%cp(i)%r(j)+sy%c%molx0(j))*dunit0(iunit),'f',length=12,decimal=8,justify=3),j=1,3), &
             string(sy%f(sy%iref)%cp(i)%name,length=10,justify=ioj_center),&
             string(sy%f(sy%iref)%cp(i)%s%f,'e',decimal=8,length=15,justify=4),&
             string(sy%f(sy%iref)%cp(i)%s%gfmod,'e',decimal=8,length=15,justify=4),&
             string(sy%f(sy%iref)%cp(i)%s%del2f,'e',decimal=8,length=15,justify=4)
       endif
    enddo
    write (uout,*)

  end subroutine cp_short_report

  !> Subdivides the initial simplex tetrahedron by iniv. uses the
  !> common variables nstack, barstack and depstack to store the list
  !> of simplex where a search is launched. Accumulates the seeds in
  !> xseed (number of seeds nn).
  subroutine barycentric(iniv,depthmax,nn,xseed)
    real*8, dimension(4,3), intent(in) :: iniv
    integer, intent(in) :: depthmax
    integer, intent(inout) :: nn
    real*8, intent(inout), allocatable :: xseed(:,:)

    real*8 :: simp(4,3)
    integer :: i, j, d
    integer :: isimp(4,3), dim
    integer :: sizestack
    real*8 :: denom
    real*8 :: edge(6,3), zedge(6,3)
    real*8 :: face(4,2,3), zface(4,3)
    real*8 :: body(3,3), zbody(3)

    ! build the edges, faces and body
    edge(1,:) = iniv(2,:) - iniv(1,:)
    edge(2,:) = iniv(3,:) - iniv(1,:)
    edge(3,:) = iniv(4,:) - iniv(1,:)
    edge(4,:) = iniv(3,:) - iniv(2,:)
    edge(5,:) = iniv(4,:) - iniv(2,:)
    edge(6,:) = iniv(4,:) - iniv(3,:)
    zedge(1,:) = iniv(1,:)
    zedge(2,:) = iniv(1,:)
    zedge(3,:) = iniv(1,:)
    zedge(4,:) = iniv(2,:)
    zedge(5,:) = iniv(2,:)
    zedge(6,:) = iniv(3,:)
    ! faces
    face(1,1,:) = iniv(2,:) - iniv(1,:)
    face(1,2,:) = iniv(3,:) - iniv(1,:)
    face(2,1,:) = iniv(2,:) - iniv(1,:)
    face(2,2,:) = iniv(4,:) - iniv(1,:)
    face(3,1,:) = iniv(3,:) - iniv(1,:)
    face(3,2,:) = iniv(4,:) - iniv(1,:)
    face(4,1,:) = iniv(3,:) - iniv(2,:)
    face(4,2,:) = iniv(4,:) - iniv(2,:)
    zface(1,:) =  iniv(1,:)
    zface(2,:) =  iniv(1,:)
    zface(3,:) =  iniv(1,:)
    zface(4,:) =  iniv(2,:)
    ! body
    body(1,:) = iniv(2,:) - iniv(1,:)
    body(2,:) = iniv(3,:) - iniv(1,:)
    body(3,:) = iniv(4,:) - iniv(1,:)
    zbody(:) = iniv(1,:)

    ! allocate and nullify the stack of simplex
    d = depthmax+1
    sizestack = (2**d-1) + (6**d-1) / 5 + (24**d-1) / 23
    allocate(barstack(sizestack,4,3))
    allocate(depstack(sizestack))
    barstack = 0 
    nstack = 0

    ! dim = 1 (vertex of the tetrahedron)
    simp = 0d0
    do i = 1, 4
       simp(1,:) = iniv(i,1:3)
       call seed_from_simplex(simp,1,nn,xseed)
    end do

    ! barycentric subdivision
    isimp(1,:) = (/0,0,0/)
    isimp(2,:) = (/1,0,0/)
    isimp(3,:) = (/0,1,0/)
    isimp(4,:) = (/0,0,1/)
    do dim = 2, 4
       call barycentric_divide(dim,isimp(1:dim,1:dim-1),depthmax)
    end do

    ! search the simplex in the stack
    do i = 1, nstack
       dim = count ((/any(barstack(i,:,1) > 0),any(barstack(i,:,2) > 0),any(barstack(i,:,3) > 0)/)) + 1
       denom = real(bardenom(dim)**(depthmax-depstack(i)),8)
       if (dim == 2) then
          ! all 6 edges
          do j = 1, 6
             simp(1,:) = zedge(j,:) + barstack(i,1,1) * edge(j,:) / denom
             simp(2,:) = zedge(j,:) + barstack(i,2,1) * edge(j,:) / denom
             call seed_from_simplex(simp,dim,nn,xseed)
          end do
       else if (dim == 3) then
          ! all 4 faces
          do j = 1, 4
             simp(1,:) = zface(j,:) + (barstack(i,1,1)*face(j,1,:) + barstack(i,1,2)*face(j,2,:)) / denom
             simp(2,:) = zface(j,:) + (barstack(i,2,1)*face(j,1,:) + barstack(i,2,2)*face(j,2,:)) / denom
             simp(3,:) = zface(j,:) + (barstack(i,3,1)*face(j,1,:) + barstack(i,3,2)*face(j,2,:)) / denom
             call seed_from_simplex(simp,dim,nn,xseed)
          end do
       else if (dim == 4) then
          ! the body
          simp = matmul(barstack(i,:,:),body) / denom
          do j = 1, 4
             simp(j,:) = simp(j,:) + zbody
          end do
          call seed_from_simplex(simp,dim,nn,xseed)
       end if
    end do

    ! deallocate
    deallocate(barstack,depstack)

  end subroutine barycentric

  !> Recursive subdivision of the simplex simp up to some depth. uses
  !> convex coordinates and integer arithmetic.
  recursive subroutine barycentric_divide(n,simp,depth)
    use tools_io, only: ferror, faterr
    integer, intent(in) :: n
    integer, intent(in) :: simp(:,:)
    integer, intent(in) :: depth

    integer :: newsimp(n,n-1)

    if (size(simp,1) /= n .or. size(simp,2) /= n-1) &
       call ferror("barycentric_divide","inconsistent size of simp",faterr)

    ! add this simplex to the stack
    nstack = nstack + 1
    barstack(nstack,1:n,1:n-1) = simp(:,:)
    depstack(nstack) = depth

    ! subdivide?
    if (depth == 0) return

    select case (n)
    case (2)
       ! a line
       newsimp(1,1) = simp(1,1)*2
       newsimp(2,1) = simp(1,1)+simp(2,1)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,1) = simp(1,1)+simp(2,1)
       newsimp(2,1) = 2*simp(2,1)
       call barycentric_divide(n,newsimp,depth-1)
    case (3)
       ! a triangle
       newsimp(1,:) = 6*simp(1,:)
       newsimp(2,:) = 3*simp(1,:) + 3*simp(2,:)
       newsimp(3,:) = 2*simp(1,:) + 2*simp(2,:) + 2*simp(3,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) = 6*simp(1,:)
       newsimp(2,:) = 3*simp(1,:) + 3*simp(3,:)
       newsimp(3,:) = 2*simp(1,:) + 2*simp(3,:) + 2*simp(2,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) = 6*simp(2,:)
       newsimp(2,:) = 3*simp(2,:) + 3*simp(1,:)
       newsimp(3,:) = 2*simp(2,:) + 2*simp(1,:) + 2*simp(3,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) = 6*simp(2,:)
       newsimp(2,:) = 3*simp(2,:) + 3*simp(3,:)
       newsimp(3,:) = 2*simp(2,:) + 2*simp(3,:) + 2*simp(1,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) = 6*simp(3,:)
       newsimp(2,:) = 3*simp(3,:) + 3*simp(1,:)
       newsimp(3,:) = 2*simp(3,:) + 2*simp(1,:) + 2*simp(2,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) = 6*simp(3,:)
       newsimp(2,:) = 3*simp(3,:) + 3*simp(2,:)
       newsimp(3,:) = 2*simp(3,:) + 2*simp(2,:) + 2*simp(1,:)
       call barycentric_divide(n,newsimp,depth-1)
    case (4)
       ! a tetrahedron
       newsimp(1,:) =12*simp(1,:)
       newsimp(2,:) = 6*simp(1,:) + 6*simp(2,:)
       newsimp(3,:) = 4*simp(1,:) + 4*simp(2,:) + 4*simp(3,:)
       newsimp(4,:) = 3*simp(1,:) + 3*simp(2,:) + 3*simp(3,:) + 3*simp(4,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(1,:)
       newsimp(2,:) = 6*simp(1,:) + 6*simp(2,:)
       newsimp(3,:) = 4*simp(1,:) + 4*simp(2,:) + 4*simp(4,:)
       newsimp(4,:) = 3*simp(1,:) + 3*simp(2,:) + 3*simp(4,:) + 3*simp(3,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(1,:)
       newsimp(2,:) = 6*simp(1,:) + 6*simp(3,:)
       newsimp(3,:) = 4*simp(1,:) + 4*simp(3,:) + 4*simp(2,:)
       newsimp(4,:) = 3*simp(1,:) + 3*simp(3,:) + 3*simp(2,:) + 3*simp(4,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(1,:)
       newsimp(2,:) = 6*simp(1,:) + 6*simp(3,:)
       newsimp(3,:) = 4*simp(1,:) + 4*simp(3,:) + 4*simp(4,:)
       newsimp(4,:) = 3*simp(1,:) + 3*simp(3,:) + 3*simp(4,:) + 3*simp(2,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(1,:)
       newsimp(2,:) = 6*simp(1,:) + 6*simp(4,:)
       newsimp(3,:) = 4*simp(1,:) + 4*simp(4,:) + 4*simp(2,:)
       newsimp(4,:) = 3*simp(1,:) + 3*simp(4,:) + 3*simp(2,:) + 3*simp(3,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(1,:)
       newsimp(2,:) = 6*simp(1,:) + 6*simp(4,:)
       newsimp(3,:) = 4*simp(1,:) + 4*simp(4,:) + 4*simp(3,:)
       newsimp(4,:) = 3*simp(1,:) + 3*simp(4,:) + 3*simp(3,:) + 3*simp(2,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(2,:)
       newsimp(2,:) = 6*simp(2,:) + 6*simp(1,:)
       newsimp(3,:) = 4*simp(2,:) + 4*simp(1,:) + 4*simp(3,:)
       newsimp(4,:) = 3*simp(2,:) + 3*simp(1,:) + 3*simp(3,:) + 3*simp(4,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(2,:)
       newsimp(2,:) = 6*simp(2,:) + 6*simp(1,:)
       newsimp(3,:) = 4*simp(2,:) + 4*simp(1,:) + 4*simp(4,:)
       newsimp(4,:) = 3*simp(2,:) + 3*simp(1,:) + 3*simp(4,:) + 3*simp(3,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(2,:)
       newsimp(2,:) = 6*simp(2,:) + 6*simp(3,:)
       newsimp(3,:) = 4*simp(2,:) + 4*simp(3,:) + 4*simp(1,:)
       newsimp(4,:) = 3*simp(2,:) + 3*simp(3,:) + 3*simp(1,:) + 3*simp(4,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(2,:)
       newsimp(2,:) = 6*simp(2,:) + 6*simp(3,:)
       newsimp(3,:) = 4*simp(2,:) + 4*simp(3,:) + 4*simp(4,:)
       newsimp(4,:) = 3*simp(2,:) + 3*simp(3,:) + 3*simp(4,:) + 3*simp(1,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(2,:)
       newsimp(2,:) = 6*simp(2,:) + 6*simp(4,:)
       newsimp(3,:) = 4*simp(2,:) + 4*simp(4,:) + 4*simp(1,:)
       newsimp(4,:) = 3*simp(2,:) + 3*simp(4,:) + 3*simp(1,:) + 3*simp(3,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(2,:)
       newsimp(2,:) = 6*simp(2,:) + 6*simp(4,:)
       newsimp(3,:) = 4*simp(2,:) + 4*simp(4,:) + 4*simp(3,:)
       newsimp(4,:) = 3*simp(2,:) + 3*simp(4,:) + 3*simp(3,:) + 3*simp(1,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(3,:)
       newsimp(2,:) = 6*simp(3,:) + 6*simp(1,:)
       newsimp(3,:) = 4*simp(3,:) + 4*simp(1,:) + 4*simp(2,:)
       newsimp(4,:) = 3*simp(3,:) + 3*simp(1,:) + 3*simp(2,:) + 3*simp(4,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(3,:)
       newsimp(2,:) = 6*simp(3,:) + 6*simp(1,:)
       newsimp(3,:) = 4*simp(3,:) + 4*simp(1,:) + 4*simp(4,:)
       newsimp(4,:) = 3*simp(3,:) + 3*simp(1,:) + 3*simp(4,:) + 3*simp(2,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(3,:)
       newsimp(2,:) = 6*simp(3,:) + 6*simp(2,:)
       newsimp(3,:) = 4*simp(3,:) + 4*simp(2,:) + 4*simp(1,:)
       newsimp(4,:) = 3*simp(3,:) + 3*simp(2,:) + 3*simp(1,:) + 3*simp(4,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(3,:)
       newsimp(2,:) = 6*simp(3,:) + 6*simp(2,:)
       newsimp(3,:) = 4*simp(3,:) + 4*simp(2,:) + 4*simp(4,:)
       newsimp(4,:) = 3*simp(3,:) + 3*simp(2,:) + 3*simp(4,:) + 3*simp(1,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(3,:)
       newsimp(2,:) = 6*simp(3,:) + 6*simp(4,:)
       newsimp(3,:) = 4*simp(3,:) + 4*simp(4,:) + 4*simp(1,:)
       newsimp(4,:) = 3*simp(3,:) + 3*simp(4,:) + 3*simp(1,:) + 3*simp(2,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(3,:)
       newsimp(2,:) = 6*simp(3,:) + 6*simp(4,:)
       newsimp(3,:) = 4*simp(3,:) + 4*simp(4,:) + 4*simp(2,:)
       newsimp(4,:) = 3*simp(3,:) + 3*simp(4,:) + 3*simp(2,:) + 3*simp(1,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(4,:)
       newsimp(2,:) = 6*simp(4,:) + 6*simp(1,:)
       newsimp(3,:) = 4*simp(4,:) + 4*simp(1,:) + 4*simp(2,:)
       newsimp(4,:) = 3*simp(4,:) + 3*simp(1,:) + 3*simp(2,:) + 3*simp(3,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(4,:)
       newsimp(2,:) = 6*simp(4,:) + 6*simp(1,:)
       newsimp(3,:) = 4*simp(4,:) + 4*simp(1,:) + 4*simp(3,:)
       newsimp(4,:) = 3*simp(4,:) + 3*simp(1,:) + 3*simp(3,:) + 3*simp(2,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(4,:)
       newsimp(2,:) = 6*simp(4,:) + 6*simp(2,:)
       newsimp(3,:) = 4*simp(4,:) + 4*simp(2,:) + 4*simp(1,:)
       newsimp(4,:) = 3*simp(4,:) + 3*simp(2,:) + 3*simp(1,:) + 3*simp(3,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(4,:)
       newsimp(2,:) = 6*simp(4,:) + 6*simp(2,:)
       newsimp(3,:) = 4*simp(4,:) + 4*simp(2,:) + 4*simp(3,:)
       newsimp(4,:) = 3*simp(4,:) + 3*simp(2,:) + 3*simp(3,:) + 3*simp(1,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(4,:)
       newsimp(2,:) = 6*simp(4,:) + 6*simp(3,:)
       newsimp(3,:) = 4*simp(4,:) + 4*simp(3,:) + 4*simp(1,:)
       newsimp(4,:) = 3*simp(4,:) + 3*simp(3,:) + 3*simp(1,:) + 3*simp(2,:)
       call barycentric_divide(n,newsimp,depth-1)
       newsimp(1,:) =12*simp(4,:)
       newsimp(2,:) = 6*simp(4,:) + 6*simp(3,:)
       newsimp(3,:) = 4*simp(4,:) + 4*simp(3,:) + 4*simp(2,:)
       newsimp(4,:) = 3*simp(4,:) + 3*simp(3,:) + 3*simp(2,:) + 3*simp(1,:)
       call barycentric_divide(n,newsimp,depth-1)
    end select

  end subroutine barycentric_divide

  !> Using simplex simp(:,:), launch a critical point search. dim is
  !> the dimension of the search. depth, current depth of the
  !> barycentric subdivision. id, point number count. mid,
  !> max. number of seeds. depth, id and mid are only used in the
  !> output to stdout. this routine directly writes the cp to the cp
  !> list if one is found.
  subroutine seed_from_simplex(simp, dim, nn, xseed)
    use systemmod, only: sy
    use types, only: realloc
    real*8, intent(in) :: simp(4,3)
    integer, intent(in) :: dim
    integer, intent(inout) :: nn
    real*8, intent(inout), allocatable :: xseed(:,:)

    integer :: i
    real*8 :: xc(3), xx(3), xc0(3)

    if (dim == 1) then
       xc = simp(1,:)
       xx = sy%c%c2x(xc)
    else
       xc = 0d0
       do i = 1, dim
          xc = xc + simp(i,:) / real(dim,8)
       end do
       xc0 = xc
       xx = sy%c%c2x(xc)
    endif
    nn = nn + 1
    if (nn > size(xseed,2)) call realloc(xseed,3,2*nn)
    xseed(:,nn) = xx

  end subroutine seed_from_simplex

  subroutine cp_long_report()
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_left, ioj_center, ioj_right
    integer :: i, j
    character :: neqlab
    character*(1) :: smallnamecrit(0:3)

    data smallnamecrit   /'n','b','r','c'/

    if (sy%c%ismolecule) return

    ! Symmetry information 
    write (uout,'("* Complete CP list")')
    write (uout,'("# (x symbols are the non-equivalent representative atoms)")')
    write (uout,'("#  cp   ncp  typ      position (cryst. coords.)        op.    (lvec+cvec)")')
    do i = 1, sy%f(sy%iref)%ncpcel
       if (sy%f(sy%iref)%cpcel(i)%ir == 1 .and. sy%f(sy%iref)%cpcel(i)%ic == 1) then
          neqlab = "x"
       else
          neqlab = " "
       end if
       write (uout,'(11(A,X))') &
          string(neqlab,length=1), string(i,length=6,justify=ioj_left),&
          string(sy%f(sy%iref)%cpcel(i)%idx,length=4,justify=ioj_left),&
          string(smallnamecrit(sy%f(sy%iref)%cpcel(i)%typind),length=1),&
          (string(sy%f(sy%iref)%cpcel(i)%x(j),'f',length=12,decimal=8,justify=3),j=1,3), &
          string(sy%f(sy%iref)%cpcel(i)%ir,length=3,justify=ioj_center),&
          (string(sy%f(sy%iref)%cpcel(i)%lvec(j) + sy%c%cen(j,sy%f(sy%iref)%cpcel(i)%ic),'f',length=5,decimal=1,justify=3),j=1,3)
    end do
    write (uout,*)

    ! graph information 
    write (uout,'("* Complete CP list, bcp and rcp connectivity table")')
    write (uout,'("# (cp(end)+lvec connected to bcp/rcp)")')
    write (uout,'("#cp  ncp   typ        position (cryst. coords.)            end1 (lvec)      end2 (lvec)")')
    do i = 1, sy%f(sy%iref)%ncpcel
       if (sy%f(sy%iref)%cpcel(i)%typ == -1 .or. sy%f(sy%iref)%cpcel(i)%typ == 1 .and. dograph > 0) then
          write (uout,'(7(A,X),"(",3(A,X),") ",A," (",3(A,X),")")') &
             string(i,length=6,justify=ioj_left),&
             string(sy%f(sy%iref)%cpcel(i)%idx,length=4,justify=ioj_left),&
             string(smallnamecrit(sy%f(sy%iref)%cpcel(i)%typind),length=1),&
             (string(sy%f(sy%iref)%cpcel(i)%x(j),'f',length=13,decimal=8,justify=4),j=1,3), &
             string(sy%f(sy%iref)%cpcel(i)%ipath(1),length=4,justify=ioj_center),&
             (string(sy%f(sy%iref)%cpcel(i)%ilvec(j,1),length=2,justify=ioj_right),j=1,3),&
             string(sy%f(sy%iref)%cpcel(i)%ipath(2),length=4,justify=ioj_center),&
             (string(sy%f(sy%iref)%cpcel(i)%ilvec(j,2),length=2,justify=ioj_right),j=1,3)
       else
          write (uout,'(20(A,X))') &
             string(i,length=6,justify=ioj_left),&
             string(sy%f(sy%iref)%cpcel(i)%idx,length=4,justify=ioj_left),&
             string(smallnamecrit(sy%f(sy%iref)%cpcel(i)%typind),length=1),&
             (string(sy%f(sy%iref)%cpcel(i)%x(j),'f',length=13,decimal=8,justify=4),j=1,3)
       end if
    end do
    write (uout,*)

  end subroutine cp_long_report

  subroutine cp_vlong_report()
    use systemmod, only: sy
    use tools_io, only: uout, string
    use types, only: scalar_value
    use global, only: iunitname0, iunit
    use param, only: bohrtoa
    integer :: i, j
    real*8 :: minden, maxbden, fness, xx(3)
    type(scalar_value) :: res

    write (uout,'("* Additional properties at the critical points")')
    minden = 1d30
    maxbden = 1d-30
    do i = 1, sy%f(sy%iref)%ncp
       write (uout,'("+ Critical point no. ",A)') string(i)
       if (.not.sy%c%ismolecule) then
          write (uout,'("  Crystallographic coordinates: ",3(A,X))') &
             (string(sy%f(sy%iref)%cp(i)%x(j),'f',decimal=10),j=1,3)
          write (uout,'("  Cartesian coordinates (",A,"): ",3(A,X))') &
             iunitname0(iunit), (string(sy%f(sy%iref)%cp(i)%r(j),'f',decimal=10),j=1,3)
       else
          xx = (sy%f(sy%iref)%cp(i)%r + sy%c%molx0) * bohrtoa
          write (uout,'("  Coordinates (",A,"): ",3(A,X))') &
             iunitname0(iunit), (string(xx(j),'f',decimal=10),j=1,3)
       end if
       call sy%propty(sy%iref,sy%f(sy%iref)%cp(i)%x,res,.true.,.true.)
       if (res%f < minden) minden = res%f
       if (res%f > maxbden .and. sy%f(sy%iref)%cp(i)%typ == -1) maxbden = res%f
    enddo
    if (maxbden > 1d-12) then
       fness = minden / maxbden
    else
       fness = 0d0
    end if
    if (.not.sy%c%ismolecule) &
       write (uout,'("+ Flatness (rho_min / rho_b,max): ",A)') string(fness,'f',decimal=6)
    write (uout,*)

  end subroutine cp_vlong_report

  !> Write to the stdout information about the graph.
  subroutine graph_short_report()
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_center
    use global, only: iunit, iunitname0, dunit0
    integer :: i, j, k, maxlen
    character*(20) :: nam(2)
    logical :: isbcp

    ! bonds
    do k = 1, 2
       isbcp = (k==1)
       if (isbcp) then
          write (uout,'("* Analysis of system bonds")') 
          write (uout,'("# ncp is the bond from the non-equivalent CP list.")') 
          write (uout,'("# End-1 and End-2 are the bond path terminator from the non-equivalent")') 
          write (uout,'("#   CP list (index in parentheses).")') 
          write (uout,'("# r1 and r2 are the geometric distances between bond and terminators.")') 
          write (uout,'("# r1-B-r2 is the geometric angle between bond and terminators (angle).")') 
          write (uout,'("# p1 and p2 are the bond path lengths.")') 
          write (uout,'("# ncp  End-1    End-2   r1(",A,") r2(",A,")  r1/r2  r1-B-r2  p1(",A,") p2(",A,")")') &
             string(iunitname0(iunit)), string(iunitname0(iunit)), string(iunitname0(iunit)),&
             string(iunitname0(iunit))

       else
          write (uout,'("* Analysis of system rings")') 
          write (uout,'("# ncp is the ring from the non-equivalent CP list.")') 
          write (uout,'("# End-1 and End-2 are the ring path terminator from the non-equivalent")') 
          write (uout,'("#   CP list (index in parentheses).")') 
          write (uout,'("# r1 and r2 are the geometric distances between ring and terminators.")') 
          write (uout,'("# r1-B-r2 is the geometric angle between bond and terminators (angle).")') 
          write (uout,'("# p1 and p2 are the ring path lengths.")') 
          write (uout,'("# ncp  End-1    End-2   r1(",A,") r2(",A,")  r1/r2  r1-B-r2  p1(",A,") p2(",A,")")') &
             string(iunitname0(iunit)), string(iunitname0(iunit)), string(iunitname0(iunit)),&
             string(iunitname0(iunit))
       end if
       do i = 1, sy%f(sy%iref)%ncp
          if (sy%f(sy%iref)%cp(i)%typ /= -1 .and. isbcp .or. sy%f(sy%iref)%cp(i)%typ /= 1 .and..not.isbcp) cycle

          do j = 1, 2
             if (sy%f(sy%iref)%cp(i)%ipath(j) == 0) then
                nam(j) = ' ?? '
             elseif (sy%f(sy%iref)%cp(i)%ipath(j) == -1) then
                nam(j) = 'Inf.'
             else
                nam(j) = string(sy%f(sy%iref)%cp(sy%f(sy%iref)%cp(i)%ipath(j))%name) // &
                   " (" // string(sy%f(sy%iref)%cp(i)%ipath(j)) // ")"
             end if
          end do
          maxlen = max(8,len_trim(nam(1)),len_trim(nam(2)))
          write (uout,'(99(A,X))') string(i,length=5,justify=ioj_center),&
             string(nam(1),length=maxlen,justify=ioj_center), string(nam(2),length=maxlen,justify=ioj_center),&
             (string(sy%f(sy%iref)%cp(i)%brdist(j)*dunit0(iunit),'f',length=8,decimal=4,justify=3),j=1,2), &
             string(sy%f(sy%iref)%cp(i)%brdist(1)/sy%f(sy%iref)%cp(i)%brdist(2),'f',length=8,decimal=4,justify=3), &
             string(sy%f(sy%iref)%cp(i)%brang,'f',length=7,decimal=2,justify=3),&
             (string(sy%f(sy%iref)%cp(i)%brpathlen(j)*dunit0(iunit),'f',length=8,decimal=4,justify=3),j=1,2)
       enddo
       write (uout,*)
    end do

  end subroutine graph_short_report

  !> Attempt to build the complete graph for the system. First, start
  !> an upwards gradient path from the bcps to determine the bond
  !> paths. Then calculate the ring paths with a pair of downwards
  !> gradient paths. If dograph is >1, determine the stable and
  !> unstable cps on each 2d manifold associated to every
  !> non-equivalent bcp and rcp.
  subroutine makegraph()
    use systemmod, only: sy
    use tools_math, only: norm, eig
    use types, only: scalar_value
    use param, only: pi
    integer :: i, j, k
    integer :: nstep
    real*8 :: xdif(3,2), v(3)
    real*8, dimension(3,3) :: evec
    real*8, dimension(3) :: reval
    real*8 :: dist, xdtemp(3,2), xx(3), plen(2)
    integer :: wcp
    integer :: ier, idir
    logical :: isbcp
    type(scalar_value) :: res
    real*8, allocatable :: xdis(:,:,:), xplen(:,:)

    real*8, parameter :: change = 1d-2

    associate(f => sy%f(sy%iref), cr => sy%c)

      allocate(xdis(3,2,f%ncp),xplen(2,f%ncp))
      xdis = 0d0
      xplen = 0d0

      ! run over known non-equivalent cps  
      !$omp parallel do private(res,isbcp,evec,reval,idir,xdtemp,nstep,ier,xx,plen) schedule(dynamic)
      do i = 1, f%ncp
         ! BCP/RCP paths
         if (abs(f%cp(i)%typ) == 1) then
            isbcp = (f%cp(i)%typ == -1)

            ! diagonalize hessian at the bcp/rcp, calculate starting points
            ! the third component of the hessian is up/down direction
            call f%grd(f%cp(i)%r,2,res)
            evec = res%hf
            call eig(evec,reval)
            if (isbcp) then
               xx = evec(:,3)
            else
               xx = evec(:,1)
            end if
            xx = xx / norm(xx)
            xdtemp(:,1) = f%cp(i)%r + change * xx
            xdtemp(:,2) = f%cp(i)%r - change * xx

            ! follow up/down both directions
            if (isbcp) then
               idir = 1
            else
               idir = -1
            end if
            do j = 1, 2
               call f%gradient(xdtemp(:,j),idir,nstep,ier,.true.,plen(j)) 
            end do
            !$omp critical (xdis1)
            f%cp(i)%brvec = xx 
            xdis(:,:,i) = xdtemp
            xplen(:,i) = plen
            !$omp end critical (xdis1)
         else
            !$omp critical (xdis2)
            f%cp(i)%brvec = 0d0
            !$omp end critical (xdis2)
         end if
      end do
      !$omp end parallel do

      ! Fill the eigenvectors
      do i = 1, f%ncpcel
         if (abs(f%cp(f%cpcel(i)%idx)%typ) == 1) then
            isbcp = (f%cp(f%cpcel(i)%idx)%typ == -1)
            call f%grd(f%cpcel(i)%r,2,res)
            evec = res%hf
            call eig(evec,reval)
            if (isbcp) then
               xx = evec(:,3)
            else
               xx = evec(:,1)
            end if
            f%cpcel(i)%brvec = xx / norm(xx)
         else
            f%cpcel(i)%brvec = 0d0
         end if
      end do

      ! run over known non-equivalent cps  
      do i = 1, f%ncp
         ! BCP paths
         if (abs(f%cp(i)%typ) == 1) then
            do j = 1, 2
               ! save difference vector and distance
               xdif(:,j) = xdis(:,j,i) - f%cp(i)%r
               f%cp(i)%brdist(j) = sqrt(dot_product(xdif(:,j),xdif(:,j)))
               f%cp(i)%brpathlen(j) = xplen(j,i)

               ! check the identity of the cp
               xdis(:,j,i) = cr%c2x(xdis(:,j,i))
               call f%nearest_cp(xdis(:,j,i),wcp,dist,type=f%cp(i)%typ*3)
               if (wcp /= 0) then
                  f%cp(i)%ipath(j) = f%cpcel(wcp)%idx
               else
                  if (cr%ismolecule .and.( &
                     xdis(1,j,i) < cr%molborder(1) .or. xdis(1,j,i) > (1d0-cr%molborder(1)) .or.&
                     xdis(2,j,i) < cr%molborder(2) .or. xdis(2,j,i) > (1d0-cr%molborder(2)) .or.&
                     xdis(3,j,i) < cr%molborder(3) .or. xdis(3,j,i) > (1d0-cr%molborder(3)))) then
                     f%cp(i)%ipath(j) = -1
                  else
                     f%cp(i)%ilvec(:,j) = 0
                     do k = 1, f%ncpcel
                        if (f%cpcel(k)%idx /= i) cycle
                        f%cpcel(k)%ipath(j) = 0
                        f%cpcel(k)%ilvec(:,j) = 0
                     end do
                  endif
                  cycle
               end if

               ! generate equivalent bond/ring paths
               do k = 1, f%ncpcel
                  if (f%cpcel(k)%idx /= i) cycle
                  v = matmul(cr%rotm(1:3,1:3,f%cpcel(k)%ir),xdis(:,j,i)) + &
                     cr%rotm(:,4,f%cpcel(k)%ir) + &
                     cr%cen(:,f%cpcel(k)%ic) + f%cpcel(k)%lvec
                  call f%nearest_cp(v,wcp,dist,type=f%cp(i)%typ*3)
                  f%cpcel(k)%ipath(j) = wcp
                  f%cpcel(k)%ilvec(:,j) = nint(v - f%cpcel(wcp)%x)
               end do
            end do
            f%cp(i)%brang = dot_product(xdif(:,1),xdif(:,2)) / &
               (f%cp(i)%brdist(1)+1d-12) / (f%cp(i)%brdist(2)+1d-12)
            if (abs(f%cp(i)%brang) > 1d0) &
               f%cp(i)%brang = f%cp(i)%brang / abs(f%cp(i)%brang+1d-12)
            f%cp(i)%brang = acos(f%cp(i)%brang) * 180d0 / pi
         end if
      end do

      deallocate(xdis)

    end associate
  end subroutine makegraph

  !> Scale the Wigner-Seitz cell to make it fit inside a sphere of radius
  !> rad. 
  subroutine scale_ws(rad,wso,ntetrag,tetrag)
    use systemmod, only: sy
    use tools_math, only: norm
    real*8, intent(in) :: rad
    real*8, intent(in) :: wso(3)
    integer, intent(inout) :: ntetrag
    real*8, intent(inout) :: tetrag(:,:,:)

    integer :: i, j
    real*8 :: x(3), maxn

    do i = 1, ntetrag
       ! calculate the max norm for this tetrag
       maxn = -1d30
       do j = 1, 4
          x = tetrag(j,:,i) - wso
          x = sy%c%x2c(x)
          maxn = max(norm(x),maxn)
       end do
       ! scale to rad
       do j = 1, 4
          x = tetrag(j,:,i) - wso
          x = x / maxn * rad
          tetrag(j,:,i) = x + wso
       end do
    end do

  end subroutine scale_ws

end module autocp

