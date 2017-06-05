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

!> Qtree, gradient path tracing.
module qtree_gpaths
  implicit none

  private
  public :: gradient_qtree
  public :: gradient_color
  public :: gradient_full

contains

  !> Trace a GP, projecting on grid points if the trajectory comes
  !> close enough. xp in cartesian, f, gf and gfmod must be known
  !> output xp is garbage
  subroutine gradient_qtree (xp,base_t,res,iidx,cpout,ier,usestack,trm,fgr,lapgr)
    use qtree_basic, only: qtreeidx, qtreei, qtreer, nnuc, maxl, minlen,&
       ode_type, ode_b, ode_o, ode_a, ngrd_term, ode_b2, safety, qinv, q1inv,&
       periodic, r_betagp, r_betaint, ode_fsal, nder, savefgr, savelapgr, neargp
    use varbas, only: nearest_cp, cpcel
    use fields, only: f, GRD
    use global, only: color_allocate, qtreefac, stepsize, refden, ode_abserr,&
       killext, mpstep
    use struct_basic, only: cr
    use tools_io, only: ferror, faterr, uout
    use types, only: scalar_value

    real*8, intent(inout) :: xp(3)
    integer, intent(inout) :: base_t
    type(scalar_value), intent(inout) :: res
    integer(qtreeidx), intent(in) :: iidx
    integer, intent(out) :: cpout
    integer, intent(out) :: ier
    logical, intent(in) :: usestack
    integer(qtreei), intent(inout) :: trm(:,:)
    real(qtreer), intent(inout) :: fgr(:,:), lapgr(:,:)

    integer :: nn
    real*8 :: dy(3), grds(3,10), ystp(3), xp2(3)
    real*8 :: xerr(3), derr, xtemp(3)

    integer, parameter :: mstack = 5000
    real*8, parameter :: eps = 1d-12
    integer, parameter :: mstep = 40000

    real*8 :: maxdist
    integer :: nid, l2
    real*8 :: dist, h0, xpr(3)
    integer :: i, j, nstep
    integer :: lrot
    integer :: bstack(mstack), nstack
    integer(qtreeidx) :: istack(mstack), idx
    logical :: gridp, doproj
    integer :: base_to

    if (color_allocate == 0) then
       base_to = 1
    else
       base_to = base_t
    end if

    ! if starting at a cp, mark it as unassigned
    if (res%gfmod < eps) then
       cpout = nnuc+1
       return
    end if

    doproj = (color_allocate == 1)
    l2 = 2**maxl
    cpout = 0
    ier = 0
    maxdist = min(minlen / 2**maxl / qtreefac,0.1d0)
    h0 = stepsize
    gridp = .false.
    lrot = 1

    nstack = 1
    istack(1) = iidx
    idx = iidx
    bstack(1) = base_to

    do nstep = 1, mstep

       ! step
       if (ode_type == 0) then
          ! fixed-step
          grds = 0d0
          grds(:,1) = res%gf / (res%gfmod+1d-80)
          ystp = ode_b(1) * grds(:,1)
          do i = 2, ode_o
             dy = 0d0
             do j = 1, i-1
                dy = dy + ode_a(i,j) * grds(:,j)
             end do
             xtemp = xp + h0 * dy
             call grd(f(refden),xtemp,1,res0=res)
             grds(:,i) = res%gf / (res%gfmod+1d-80)
             ystp = ystp + ode_b(i) * grds(:,i)
          end do
          xp = xp + h0 * ystp
          ngrd_term = ngrd_term + ode_o - 1
       else
          ! embedded methods
          grds(:,1) = res%gf / (res%gfmod+1d-80)
          do nn = 1, 40
             ! fixed-step
             ystp = ode_b(1) * grds(:,1)
             xerr = (ode_b(1) - ode_b2(1)) * grds(:,1)
             do i = 2, ode_o
                dy = 0d0
                do j = 1, i-1
                   dy = dy + ode_a(j,i) * grds(:,j)
                end do
                xtemp = xp + h0 * dy
                call grd(f(refden),xtemp,1,res0=res)
                grds(:,i) = res%gf / (res%gfmod+1d-80)
                ystp = ystp + ode_b(i) * grds(:,i)
                xerr = xerr + (ode_b(i) - ode_b2(i)) * grds(:,i)
             end do
             xp2 = xp + h0 * ystp
             ngrd_term = ngrd_term + ode_o - 1
             xerr = h0 * xerr
             derr = sqrt(dot_product(xerr,xerr)) / ode_abserr

             ! adapt the step, gsl version
             if (derr > 1.1d0) then
                derr = min(safety / derr**qinv,0.2d0)
                h0 = h0 * derr
             else if (derr < 0.5d0) then
                derr = max(min(safety / derr**q1inv,5d0), 1d0)
                ! do not limit step-size
                !                h0 = min(h0 * derr, stepsize)
                h0 = h0 * derr
                xp = xp2
                exit 
             else
                xp = xp2
                exit
             end if
          end do
          if (nn == 40) then
             call ferror('gradient_qtree','could not find good step',faterr)
          end if
       end if

       ! get nearest grid point
       if (doproj) then
          xpr = xp
          call neargp(xpr,base_t,lrot,idx,dist)
          base_to = base_t
       end if
       if (base_t == 0) then
          if (killext .and. .not.periodic) then
             ! Check if the point is known
             cpout = nnuc+2
             if (usestack) then
                do i = 2, nstack
                   trm(istack(i),bstack(i)) = int(cpout,1)
                end do
             end if
             ier = 0
             return
          else if (periodic) then
             call ferror('grdient_qtree','did not locate tetrah (tetraeps?)',faterr)
          else
             doproj = .false.
          end if
       end if

       ! if the nearest point is different from the last one, project
       ! to the grid. Do not project in the first mpstep steps or if the
       ! distance is < maxdist.
       gridp = .false.
       if (doproj .and. idx /= istack(nstack) .and. nstep > mpstep .and. dist < maxdist) then
          if (trm(idx,base_to) /= nnuc+1) then
             if (trm(idx,base_to) /= 0) then
                ! Check if the point is known
                cpout = abs(trm(idx,base_to))
                if (usestack) then
                   if (gridp) then
                      do i = 2, nstack-1
                         trm(istack(i),bstack(i)) = int(cpout,1)
                      end do
                      trm(istack(nstack),bstack(nstack)) = trm(idx,base_to)
                   else
                      do i = 2, nstack
                         trm(istack(i),bstack(i)) = int(cpout,1)
                      end do
                   end if
                end if
                ier = 0
                return
             end if

             ! project and add to the stack
             xp = xpr
             nstack = nstack + 1
             if (nstack > mstack) then
                !$omp critical (IO)
                write (uout,*) 'Too many stack points'
                write (uout,*) 'Try with a better ODE method or a smaller stepsize to avoid bouncing'
                !$omp end critical (IO)
                call ferror('grdient_qtree','too many stack points',faterr)
             end if
             istack(nstack) = idx
             bstack(nstack) = base_to

             gridp = .true.
          end if
       end if

       ! tasks in crystallographic
       xp = cr%c2x(xp)

       ! get nearest nuclear CP
       call nearest_cp(xp,nid,dist,type=f(refden)%typnuc)

       ! is it inside a beta-sphere?
       if (dist <= r_betagp(cpcel(nid)%idx)) then
          xp = cpcel(nid)%x
          cpout = cpcel(nid)%idx
          if (usestack) then
             if (gridp) then
                do i = 2, nstack-1
                   trm(istack(i),bstack(i)) = int(cpout,1)
                end do
                if (dist <= r_betaint(cpcel(nid)%idx)) then
                   trm(istack(nstack),bstack(nstack)) = int(-cpout,1)
                else
                   trm(istack(nstack),bstack(nstack)) = int(cpout,1)
                end if
             else
                do i = 2, nstack
                   trm(istack(i),bstack(i)) = int(cpout,1)
                end do
             end if
          end if
          ier = 0
          return
       end if

       ! tasks in cartesian
       xp = cr%x2c(xp)

       ! accept the new point and recalculate f, if applicable
       if (gridp .or. .not.ode_fsal) then
          call grd(f(refden),xp,nder,res0=res)
          ngrd_term = ngrd_term + 1
          if (gridp) then
             if (savefgr) fgr(idx,base_to) = res%fval
             if (savelapgr) lapgr(idx,base_to) = -res%del2fval
          end if
       end if

    end do

    ! ! if (debug > 0) then
    !    xp = cr%c2x(xp)
    !    call nearest_cp(xp,nid,dist,type=f(refden)%typnuc)
    !    !$omp critical (IO)
    !    write (uout,'("+ Unknown term., marked as unknown : ",3(F20.12,2X))') xp
    !    write (uout,'("  Closest nuc.: ",I4," at distance ",F20.12)') nid, dist
    !    !$omp end critical (IO)
    ! ! end if
    cpout = nnuc+3

    ! write (uout,'("* Last point: ",3(F20.12,2X))') xp
    ! write (uout,'("* Closest nuc.: ",I4," at distance ",F20.12)') nid, dist
    ! call ferror('grdient_qtree','too many steps in gradient',faterr)

  end subroutine gradient_qtree

  !> Trace a GP, checking adjacent grid point colors. xp in cartesian,
  !> f, gf and gfmod must be known output xp is garbage.
  subroutine gradient_color (xp,base_t,rver,res,cpout,ier,trm)
    use qtree_basic, only: qtreei, qtreeidx, nnuc, maxl, ode_type, ode_b,&
       ode_o, ode_a, ngrd_term, ode_b2, safety, periodic,&
       qinv, q1inv, lrotm, cindex, ode_fsal, r_betagp, crys2convex, locate_tetrah
    use varbas, only: nearest_cp, cpcel
    use fields, only: f, grd
    use global, only: color_allocate, stepsize, refden, ode_abserr,&
       ws_origin, killext
    use struct_basic, only: cr
    use tools_io, only: ferror, faterr
    use types, only: scalar_value

    real*8, intent(inout) :: xp(3)
    integer, intent(inout) :: base_t
    real*8, intent(inout) :: rver(3)
    type(scalar_value), intent(inout) :: res
    integer, intent(out) :: cpout
    integer, intent(out) :: ier
    integer(qtreei), intent(inout) :: trm(:,:)

    real*8, parameter :: eps = 1d-12
    integer, parameter :: mstep = 4000

    integer :: neightrm
    integer :: nid, l2, lrot
    integer(qtreeidx) :: idx
    real*8 :: dist, h0, xnew(3)
    integer :: i, j, nstep
    integer :: nbase, neigh(3,8), ntrm
    logical :: agree
    integer :: nn
    real*8 :: dy(3), grds(3,10), ystp(3), xp2(3)
    real*8 :: xerr(3), derr, xtemp(3)
    logical :: docolor
    integer :: base_to

    if (color_allocate == 0) then
       base_to = 1
    else
       base_to = base_t
    end if

    ! if starting at a cp, mark it as unassigned
    if (res%gfmod < eps) then
       cpout = nnuc+1
       return
    end if

    docolor = .true.
    l2 = 2**maxl
    cpout = 0
    ier = 0
    h0 = stepsize
    lrot = 1

    do nstep = 1, mstep

       ! step
       if (ode_type == 0) then
          ! fixed-step
          grds = 0d0
          grds(:,1) = res%gf / (res%gfmod+1d-80)
          ystp = ode_b(1) * grds(:,1)
          do i = 2, ode_o
             dy = 0d0
             do j = 1, i-1
                dy = dy + ode_a(i,j) * grds(:,j)
             end do
             xtemp = xp + h0 * dy
             call grd(f(refden),xtemp,1,res0=res)
             grds(:,i) = res%gf / (res%gfmod+1d-80)
             ystp = ystp + ode_b(i) * grds(:,i)
          end do
          xp = xp + h0 * ystp
          ngrd_term = ngrd_term + ode_o - 1
       else
          ! embedded methods
          grds(:,1) = res%gf / (res%gfmod+1d-80)
          do nn = 1, 40
             ! fixed-step
             ystp = ode_b(1) * grds(:,1)
             xerr = (ode_b(1) - ode_b2(1)) * grds(:,1)
             do i = 2, ode_o
                dy = 0d0
                do j = 1, i-1
                   dy = dy + ode_a(j,i) * grds(:,j)
                end do
                xtemp = xp + h0 * dy
                call grd(f(refden),xtemp,1,res0=res)
                grds(:,i) = res%gf / (res%gfmod+1d-80)
                ystp = ystp + ode_b(i) * grds(:,i)
                xerr = xerr + (ode_b(i) - ode_b2(i)) * grds(:,i)
             end do
             xp2 = xp + h0 * ystp
             ngrd_term = ngrd_term + ode_o - 1
             xerr = h0 * xerr
             derr = sqrt(dot_product(xerr,xerr)) / ode_abserr

             ! adapt the step, gsl version
             if (derr > 1.1d0) then
                derr = min(safety / derr**qinv,0.2d0)
                h0 = h0 * derr
             else if (derr < 0.5d0) then
                derr = max(min(safety / derr**q1inv,5d0), 1d0)
                ! do not limit step-size
                ! h0 = min(h0 * derr, stepsize)
                h0 = h0 * derr
                xp = xp2
                exit 
             else
                xp = xp2
                exit
             end if
          end do
          if (nn == 40) then
             call ferror('gradient_color','could not find good step',faterr)
          end if
       end if

       ! tasks in crystallographic
       xp = cr%c2x(xp)

       ! transform to WS and find base tetrahedron and rotation matrix.
       if (docolor) then
          xnew = matmul(lrotm(:,1:3,lrot),xp-ws_origin) + ws_origin
          call crys2convex(xnew,nbase,rver)
          if (nbase /= base_t) then
             call locate_tetrah(xp,base_t,rver,lrot)
             if (base_t == 0) then
                if (periodic) then
                   call ferror('gradient_color','did not find tetrah',faterr)
                else
                   if (killext) then
                      cpout = nnuc+2
                      ier = 0
                      return
                   else
                      docolor = .false.
                      goto 1 ! step the search for neighbors
                   end if
                end if
             end if
             if (color_allocate == 0) then
                docolor = .false.
                goto 1 ! step the search for neighbors
             else
                base_to = base_t
             end if
          end if

          ! Check neighbors
          rver = rver * l2
          where(rver < 0d0) rver = 0d0
          where(rver > l2) rver = real(l2,8)
          neigh(:,1) = (/floor(rver(1)),   floor(rver(2)),   floor(rver(3))  /)
          neigh(:,2) = (/ceiling(rver(1)), floor(rver(2)),   floor(rver(3))  /)
          neigh(:,3) = (/floor(rver(1)),   ceiling(rver(2)), floor(rver(3))  /)
          neigh(:,4) = (/ceiling(rver(1)), ceiling(rver(2)), floor(rver(3))  /)
          neigh(:,5) = (/floor(rver(1)),   floor(rver(2)),   ceiling(rver(3))/)
          neigh(:,6) = (/ceiling(rver(1)), floor(rver(2)),   ceiling(rver(3))/)
          neigh(:,7) = (/floor(rver(1)),   ceiling(rver(2)), ceiling(rver(3))/)
          neigh(:,8) = (/ceiling(rver(1)), ceiling(rver(2)), ceiling(rver(3))/)

          ntrm = 0
          neightrm = 0
          agree = .true.
          do i = 1, 8
             if (any(neigh(:,i) < 0) .or. sum(neigh(:,i)) > l2) cycle
             idx = cindex(neigh(:,i),maxl)
             if (trm(idx,base_to) == 0) then
                agree = .false.
                exit
             end if
             ntrm = ntrm + 1
             if (ntrm == 1) then
                neightrm = abs(trm(idx,base_to))
             else
                agree = agree .and. (abs(trm(idx,base_to))==neightrm)
             end if
             if (.not.agree) exit
          end do
          if (agree) then
             cpout = neightrm
             ier = 0
             return
          end if
       end if
1      continue

       ! get nearest nuclear CP
       call nearest_cp(xp,nid,dist,type=f(refden)%typnuc)

       ! is it inside a beta-sphere?
       if (dist <= r_betagp(cpcel(nid)%idx)) then
          cpout = cpcel(nid)%idx
          ier = 0
          return
       end if

       ! tasks in cartesian
       xp = cr%x2c(xp)

       ! grd at point and next step
       if (.not.ode_fsal) then
          call grd(f(refden),xp,1,res0=res)
          ngrd_term = ngrd_term + 1
       end if

    end do

    ! ! if (debug > 0) then
    !    xp = cr%c2x(xp)
    !    call nearest_cp(xp,nid,dist,type=f(refden)%typnuc)
    !    !$omp critical (IO)
    !    write (uout,'("+ Unknown term., marked as unknown : ",3(F20.12,2X))') xp
    !    write (uout,'("  Closest nuc.: ",I4," at distance ",F20.12)') nid, dist
    !    !$omp end critical (IO)
    ! ! end if
    cpout = nnuc+3

    ! write (uout,'("* Last point: ",3(F20.12,2X))') xp
    ! write (uout,'("* Closest nuc.: ",I4," at distance ",F20.12)') nid, dist
    ! call ferror('gradient_color','too many steps in gradient',faterr)

  end subroutine gradient_color

  !> Trace a GP, without approximations. xp in cartesian,
  !> f, gf and gfmod must be known output xp is garbage.
  subroutine gradient_full (xp,base_t,rver,res,cpout,ier)
    use qtree_basic, only: nnuc, ode_o, ode_a, ode_type, ode_b, ngrd_term,&
       ode_b2, safety, qinv, q1inv, periodic, lrotm, r_betagp, ode_fsal,&
       crys2convex, locate_tetrah
    use varbas, only: nearest_cp, cpcel
    use fields, only: f, grd
    use global, only: stepsize, refden, ode_abserr, killext, ws_origin
    use struct_basic, only: cr
    use tools_io, only: ferror, faterr
    use types, only: scalar_value

    real*8, intent(inout) :: xp(3)
    integer, intent(inout) :: base_t
    real*8, intent(inout) :: rver(3)
    type(scalar_value), intent(inout) :: res
    integer, intent(out) :: cpout
    integer, intent(out) :: ier

    real*8, parameter :: eps = 1d-12
    integer, parameter :: mstep = 4000

    integer :: nid
    real*8 :: dist, h0
    integer :: i, j, nstep
    integer :: nn
    real*8 :: dy(3), grds(3,10), ystp(3), xp2(3)
    real*8 :: xerr(3), derr, xtemp(3), xnew(3)
    integer :: nbase, lrot

    ! if starting at a cp, mark it as unassigned
    if (res%gfmod < eps) then
       cpout = nnuc+1
       return
    end if

    cpout = 0
    h0 = stepsize
    lrot = 1

    do nstep = 1, mstep
       ! step
       if (ode_type == 0) then
          ! fixed-step
          grds = 0d0
          grds(:,1) = res%gf / (res%gfmod+1d-80)
          ystp = ode_b(1) * grds(:,1)
          do i = 2, ode_o
             dy = 0d0
             do j = 1, i-1
                dy = dy + ode_a(i,j) * grds(:,j)
             end do
             xtemp = xp + h0 * dy
             call grd(f(refden),xtemp,1,res0=res)
             grds(:,i) = res%gf / (res%gfmod+1d-80)
             ystp = ystp + ode_b(i) * grds(:,i)
          end do
          xp = xp + h0 * ystp
          ngrd_term = ngrd_term + ode_o - 1
       else
          ! embedded methods
          grds(:,1) = res%gf / (res%gfmod+1d-80)
          do nn = 1, 40
             ! fixed-step
             ystp = ode_b(1) * grds(:,1)
             xerr = (ode_b(1) - ode_b2(1)) * grds(:,1)
             do i = 2, ode_o
                dy = 0d0
                do j = 1, i-1
                   dy = dy + ode_a(j,i) * grds(:,j)
                end do
                xtemp = xp + h0 * dy
                call grd(f(refden),xtemp,1,res0=res)
                grds(:,i) = res%gf / (res%gfmod+1d-80)
                ystp = ystp + ode_b(i) * grds(:,i)
                xerr = xerr + (ode_b(i) - ode_b2(i)) * grds(:,i)
             end do
             xp2 = xp + h0 * ystp
             ngrd_term = ngrd_term + ode_o - 1
             xerr = h0 * xerr
             derr = sqrt(dot_product(xerr,xerr)) / ode_abserr

             ! adapt the step, gsl version
             if (derr > 1.1d0) then
                derr = min(safety / derr**qinv,0.2d0)
                h0 = h0 * derr
             else if (derr < 0.5d0) then
                derr = max(min(safety / derr**q1inv,5d0), 1d0)
                ! do not limit step-size
                h0 = h0 * derr
                xp = xp2
                exit 
             else
                xp = xp2
                exit
             end if
          end do
          if (nn == 40) then
             call ferror('gradient_full','could not find good step',faterr)
          end if
       end if

       ! tasks in crystallographic
       xp = cr%c2x(xp)

       ! check if it has wandered out of the integration region
       if (killext .and. .not.periodic) then
          call crys2convex(xp,nbase,rver)
          if (nbase /= base_t) then
             xnew = matmul(lrotm(:,1:3,lrot),xp-ws_origin) + ws_origin
             call locate_tetrah(xp,base_t,rver,lrot)
             if (base_t == 0) then
                cpout = nnuc+2
                ier = 0
                return
             end if
          end if
       end if

       ! get nearest nuclear CP
       call nearest_cp(xp,nid,dist,type=f(refden)%typnuc)

       ! is it inside a beta-sphere?
       if (dist <= r_betagp(cpcel(nid)%idx)) then
          cpout = cpcel(nid)%idx
          ier = 0
          return
       end if

       ! tasks in cartesian
       xp = cr%x2c(xp)

       ! grd at point and next step
       if (.not.ode_fsal) then
          call grd(f(refden),xp,1,res0=res)
          ngrd_term = ngrd_term + 1
       end if

    end do

    cpout = nnuc+3

  end subroutine gradient_full

end module qtree_gpaths
