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

!> Numerical algorithms for the navigation of the electron density and related 
!> scalar fields. In particular, integration of the trajectories of the 
!> gradient and determination of the position of the CPs. Some additional
!> routines like the analysis of the atomic shells are also provided
module navigation
  implicit none

  private

  public :: gradient
  public :: newton
  private :: adaptive_stepper
  private :: stepper_euler1
  private :: stepper_heun
  private :: stepper_bs
  private :: stepper_rkck
  private :: stepper_dp
  public :: prunepath

contains

  !> Do a Newton-Raphson search at point r (cartesian). A CP is found
  !> when the gradient is less than gfnormeps. ier is the exit code:
  !> 0 (success), 1 (singular Hessian), and 2 (too many iterations).
  subroutine newton (r,gfnormeps,ier)
    use fields
    use varbas
    use global
    use struct_basic
    use tools_io
    use tools_math
    use types

    real*8, dimension(3), intent(inout) :: r
    integer, intent(out) :: ier
    real*8, intent(in) :: gfnormeps

    real*8 :: r1(3), xx(3), er
    integer :: iw(3)
    integer :: it
    logical :: ok
    type(scalar_value) :: res

    integer, parameter :: maxit = 200
  
    ! Transform to main unit cell
    r = cr%c2x(r)
    r = cr%x2c(r - floor(r))

    do it = 1, maxit
       ! Evaluate and stop criterion
       call grd(f(refden),r,2,res)
       if (res%gfmod < gfnormeps) then
          ier = 0
          return
       end if

       ! Invert h matrix and do a Newton-Raphson step (H^{-1}*grad).
       if (abs(detsym(res%hf)) < 1d-30) then
          ier = 1
          return
       end if
       call dgeco(res%hf,3,3,iw,er,r1)
       call dgedi(res%hf,3,3,iw,xx,r1,1)
       r = r - matmul(res%hf,res%gf)
    end do

    ! Too many iterations
    ier = 2

  end subroutine newton
  
  !> Generalized gradient tracing routine. The gp integration starts
  !> at xpoint (Cartesian) with step step (step < 0 for a slow, i.e.,
  !> 1d-3, start). iup = 1 if the gp is traced up the density, -1 if
  !> down. mstep = max. number of steps. nstep = actual number of
  !> steps (output). ier = 0 (correct), 1 (short step), 2 (too many
  !> iterations), 3 (outside molcell in molecules). extinf = .true.
  !> fills the ax (position, cryst), arho (density), agrad (gradient),
  !> ah (hessian) arrays describing the path. up2r (optional), trace
  !> the gp only up to a certain distance reffered to xref
  !> (cartesian). up2rho, up to a density. up2beta = .true.  , stop
  !> when reaching a beta-sphere. upflag = true (output) if ended
  !> through one of up2r or up2rho.
  subroutine gradient (fid, xpoint, iup, nstep, mstep, ier, extinf, &
    ax, arho, agrad, ah, up2r, xref, up2rho, up2beta, upflag)
    use fields
    use varbas
    use global
    use struct_basic
    use tools_math, only: eigns
    use tools_io
    use types
    use param

    type(field), intent(inout) :: fid
    real*8, dimension(3), intent(inout) :: xpoint
    integer, intent(in) :: iup
    integer, intent(out) :: nstep
    integer, intent(in) :: mstep
    integer, intent(out) :: ier
    integer, intent(in) :: extinf
    real*8, intent(out), dimension(mstep,3), optional :: ax
    real*8, intent(out), dimension(mstep), optional :: arho
    real*8, intent(out), dimension(mstep,3), optional :: agrad
    real*8, intent(out), dimension(mstep,3,3), optional :: ah
    real*8, intent(in), optional :: up2r, up2rho
    logical, intent(in), optional :: up2beta
    real*8, intent(in), dimension(3), optional :: xref
    logical, intent(out), optional :: upflag

    real*8, parameter :: minstep = 1d-7
    integer, parameter :: mhist = 5

    integer :: i, j
    real*8 :: t, h0, hini
    real*8 :: dx(3), scalhist(mhist)
    real*8 :: xlast(3), xlast2(3), len, xold(3), xini(3)
    real*8 :: sphrad, xnuc(3), cprad, xcp(3), dist2
    integer :: idnuc, idcp, nhist
    integer :: lvec(3)
    logical :: ok
    type(scalar_value) :: res

    ! initialization
    ier = 0
    t = 0d0
    h0 = abs(NAV_step) * iup
    hini = h0
    xini = cr%c2x(xpoint)
    scalhist = 1d0
    nhist = 0
    xlast2 = xpoint
    xlast = xpoint

    if (present(up2r)) then
       if (present(xref)) then
          xlast2 = xref
          xlast = xref
       else
          call ferror ('gradient','up2r but no reference', faterr)
       end if
       len = 0d0
    end if
    if (present(upflag)) then
       upflag = .false.
    end if

    ! initialize spherical coordinates navigation
    xold = 10d0 ! dummy!

    ! properties at point
    call grd(fid,xpoint,2,res)

    do nstep = 1, mstep
       ! tasks in crystallographic
       xpoint = cr%c2x(xpoint)

       ! get nearest nucleus
       idnuc = 0
       call cr%nearest_atom(xpoint,idnuc,sphrad,lvec)
       xnuc = cr%x2c(cr%atcel(idnuc)%x - lvec)
       idnuc = cr%atcel(idnuc)%idx

       ! get nearest -3 CP (idncp) and +3 CP (idccp), skip hydrogens
       if ((fid%typnuc==-3 .and. iup==1 .or. fid%typnuc==3 .and. iup==-1) .and.&
          cr%at(idnuc)%z /= 1) then
          idcp = idnuc
          cprad = sphrad
          xcp = xnuc
       else
          idcp = 0
          cprad = 1d15
          xcp = 0d0
       end if
       if (ncpcel > 0) then
          cprad = cprad * cprad
          do i = 1, ncpcel
             if (.not.(cpcel(i)%typ==-3 .and. iup==1 .or. cpcel(i)%typ==3 .and. iup==-1)) cycle  ! only cages if down, only ncps if up
             dx = xpoint - cpcel(i)%x
             call cr%shortest(dx,dist2)
             if (dist2 < cprad) then
                cprad = dist2
                idcp = cpcel(i)%idx
                xcp = cr%x2c(cpcel(i)%x + nint(xpoint-cpcel(i)%x-dx))
             end if
          end do
          cprad = sqrt(cprad)
       end if

       ! is it a nuclear position? 
       ok = .false.
       ! beta-sphere if up2beta is activated
       if (present(up2beta)) then
          if (up2beta .and. idcp/=0) then
             ok = ok .or. (cprad <= cp(idcp)%rbeta)
          else
             ok = ok .or. (cprad <= Rbetadef)
          end if
       else
          ok = ok .or. (cprad <= Rbetadef)
       end if
       if (ok) then
          xpoint = xcp
          if (extinf > 0) then
             ax(nstep,:) = cr%c2x(xpoint)
             if (extinf > 1) then
                arho(nstep) = res%f
                agrad(nstep,:) =  (/ 0d0, 0d0, 0d0 /)
                forall (i=1:3, j=1:3) ah(nstep,i,j) = 0d0
                forall (i=1:3) ah(nstep,i,i) = 1d0
             end if
          end if
          ier = 0
          return
       end if

       ! save info
       if (extinf > 0) then
          ax(nstep,:) = xpoint
          if (extinf > 1) then
             arho(nstep) = res%f
             agrad(nstep,:) = res%gf
             ah(nstep,:,:) = res%hf
          end if
       end if

       ! is it outside the molcell?
       if (cr%ismolecule .and. iup==-1) then
          if (xpoint(1) < cr%molborder(1) .or. xpoint(1) > (1d0-cr%molborder(1)) .or.&
             xpoint(2) < cr%molborder(2) .or. xpoint(2) > (1d0-cr%molborder(2)) .or.&
             xpoint(3) < cr%molborder(3) .or. xpoint(3) > (1d0-cr%molborder(3))) then
             ier = 3
             return
          end if
       endif

       ! tasks in cartesian
       xpoint = cr%x2c(xpoint)

       ! is it a cp?
       ok = (res%gfmod < NAV_gradeps)
       if (ok) then
          ier = 0
          return
       end if

       ! up 2 rho?
       if (present(up2rho)) then
          if (iup == 1 .and. res%f > up2rho .or. iup == -1 .and. res%f < up2rho) then
             if (present(upflag)) then
                upflag = .true.
             end if
             ier = 0
             return
          end if
       end if

       ! up 2 r?
       if (present(up2r)) then
          len = len + sqrt(dot_product(xpoint-xlast,xpoint-xlast))
          if (len > up2r) then
             if (present(upflag)) then
                upflag = .true.
             end if
             ier = 0
             return
          end if
       end if

       ! take step
       xlast2 = xlast
       xlast = xpoint
       ok = adaptive_stepper(fid,xpoint,h0,hini,NAV_gradeps,res)

       ! add to the trajectory angle history, terminate the gradient if the 
       ! trajectory bounces around mhist times
       nhist = mod(nhist,mhist) + 1
       scalhist(nhist) = dot_product(xlast-xlast2,xpoint-xlast)
       ok = ok .and. .not.(all(scalhist < 0d0))

       if (.not.ok .or. abs(h0) < minstep) then
          ier = 1
          return
       end if
    end do

    nstep = mstep
    ier = 2

  end subroutine gradient

  !> Integration using adaptive_stepper step, old scheme. Grow the step if
  !> the angle of two successive steps is almost 180, shrink if it is
  !> less than 90 degrees.
  function adaptive_stepper(fid,xpoint,h0,maxstep,eps,res)
    use fields
    use global
    use struct_basic
    use tools_math
    use param
    use types

    logical :: adaptive_stepper
    type(field), intent(inout) :: fid
    real*8, intent(inout) :: xpoint(3)
    real*8, intent(inout) :: h0
    real*8, intent(in) :: maxstep, eps
    type(scalar_value), intent(inout) :: res

    integer :: ier, iup
    real*8 :: grdt(3), ogrdtemp(3), ogrdt(3)
    real*8 :: xtemp(3), escalar, xerrv(3)
    real*8 :: nerr
    logical :: ok, first

    real*8, parameter :: h0break = 1.d-10

    adaptive_stepper = .true.
    ier = 1
    if (h0 > 0) then
       iup = 1
    else
       iup = -1
    end if

    grdt = res%gf / (res%gfmod + VSMALL)
    ogrdt = grdt

    first = .true.
    do while (ier /= 0)
       ! new point
       if (NAV_stepper == NAV_stepper_euler) then
          call stepper_euler1(xpoint,grdt,h0,xtemp)
       else if (NAV_stepper == NAV_stepper_heun) then
          call stepper_heun(fid,xpoint,grdt,h0,xtemp,xerrv,res)
       else if (NAV_stepper == NAV_stepper_rkck) then
          call stepper_rkck(fid,xpoint,grdt,h0,xtemp,xerrv,res)
       else if (NAV_stepper == NAV_stepper_dp) then
          call stepper_dp(fid,xpoint,grdt,h0,xtemp,xerrv,res)
       else if (NAV_stepper == NAV_stepper_bs) then
          call stepper_bs(fid,xpoint,grdt,h0,xtemp,xerrv,res)
       end if

       ! FSAL for BS stepper
       if (NAV_stepper /= NAV_stepper_bs) then
          call grd(fid,xtemp,2,res)
       end if

       ! poor man's adaptive step size in Euler
       if (NAV_stepper == NAV_stepper_euler .or. NAV_stepper == NAV_stepper_heun) then
          ! angle with next step
          escalar = dot_product(ogrdt,res%gf / (res%gfmod+VSMALL))

          ! gradient eps in cartesian
          ok = (res%gfmod < 0.99d0*eps)

          ! Check if they differ in > 90 deg.
          if (escalar < 0.d0.and..not.ok) then
             if (abs(h0) >= h0break) then
                h0 = 0.5d0 * h0
                ier = 1
             else
                adaptive_stepper = .false.
                return
             end if
          else
             ! Accept point. If angle is favorable, take longer steps
             if (escalar > 0.9 .and. first) &
                h0 = dsign(min(abs(maxstep), abs(1.6d0*h0)),maxstep)
             ier = 0
             xpoint = xtemp
          end if
       else
          ! use the error estimate
          nerr = norm(xerrv)
          if (nerr < NAV_maxerr) then
             ! accept point
             ier = 0
             xpoint = xtemp
             ! if this is the first time through, and the norm is very small, propose a longer step
             if (first .and. nerr < NAV_maxerr/10d0) &
                h0 = dsign(min(abs(maxstep), abs(1.6d0*h0)),maxstep)
          else
             ! propose a new shorter step using the error estimate
             h0 = 0.9d0 * h0 * NAV_maxerr / nerr
             if (abs(h0) < VSMALL) then
                adaptive_stepper = .false.
                return
             end if
          endif
       end if
       first = .false.
    enddo

  end function adaptive_stepper

  !> Euler stepper.
  subroutine stepper_euler1(xpoint,grdt,h0,xout)
    
    real*8, intent(in) :: xpoint(3), h0, grdt(3)
    real*8, intent(out) :: xout(3)
  
    xout = xpoint + h0 * grdt
  
  end subroutine stepper_euler1

  !> Heun stepper.
  subroutine stepper_heun(fid,xpoint,grdt,h0,xout,xerr,res)
    use types
    use fields
    use global
    use param
    
    type(field), intent(inout) :: fid
    real*8, intent(in) :: xpoint(3), h0, grdt(3)
    real*8, intent(out) :: xout(3), xerr(3)
    type(scalar_value), intent(inout) :: res
    
    real*8 :: ak2(3)

    xerr = xpoint + h0 * grdt
    call grd(fid,xerr,2,res)
    ak2 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + 0.5d0 * h0 * (ak2 + grdt)
    xerr = xout - xerr
  
  end subroutine stepper_heun

  !> Bogacki-Shampine embedded 2(3) method, fsal
  subroutine stepper_bs(fid,xpoint,grdt,h0,xout,xerr,res)
    use types
    use fields
    use global
    use param
    
    type(field), intent(inout) :: fid
    real*8, intent(in) :: xpoint(3), h0, grdt(3)
    real*8, intent(out) :: xout(3), xerr(3)
    type(scalar_value), intent(inout) :: res

    real*8, dimension(3) :: ak1, ak2, ak3, ak4

    ak1 = grdt

    xout = xpoint + h0 * (0.5d0*ak1)
    call grd(fid,xout,2,res)
    ak2 = res%gf / (res%gfmod+VSMALL)

    xout = xpoint + h0 * (0.75d0*ak2)
    call grd(fid,xout,2,res)
    ak3 = res%gf / (res%gfmod+VSMALL)

    xout = xpoint + h0 * (0.75d0*ak2)
    call grd(fid,xout,2,res)
    ak3 = res%gf / (res%gfmod+VSMALL)

    xout = xpoint + h0 * (2d0/9d0*ak1 + 1d0/3d0*ak2 + 4d0/9d0*ak3)
    call grd(fid,xout,2,res)
    ak4 = res%gf / (res%gfmod+VSMALL)

    xerr = xpoint + h0 * (7d0/24d0*ak1 + 1d0/4d0*ak2 + 1d0/3d0*ak3 + 1d0/8d0*ak4) - xout

  end subroutine stepper_bs

  !> Runge-Kutta-Cash-Karp embedded 4(5)-order, local extrapolation.
  subroutine stepper_rkck(fid,xpoint,grdt,h0,xout,xerr,res)
    use fields
    use global
    use types
    use param
    
    type(field), intent(inout) :: fid
    real*8, intent(in) :: xpoint(3), grdt(3), h0
    real*8, intent(out) :: xout(3), xerr(3)
    type(scalar_value), intent(inout) :: res

    real*8, parameter :: B21=.2d0, &
         B31=3.d0/40.d0, B32=9.d0/40.d0,&
         B41=.3d0, B42=-.9d0, B43=1.2d0,&
         B51=-11.d0/54.d0, B52=2.5d0, B53=-70.d0/27.d0, B54=35.d0/27.d0,&
         B61=1631.d0/55296.d0,B62=175.d0/512.d0, B63=575.d0/13824.d0, B64=44275.d0/110592.d0, B65=253.d0/4096.d0,&
         C1=37.d0/378.d0, C3=250.d0/621.d0, C4=125.d0/594.d0, C6=512.d0/1771.d0,&
         DC1=C1-2825.d0/27648.d0, DC3=C3-18575.d0/48384.d0, DC4=C4-13525.d0/55296.d0, DC5=-277.d0/14336.d0, DC6=C6-.25d0
    real*8, dimension(3) :: ak2, ak3, ak4, ak5, ak6
    
    xout = xpoint + h0*B21*grdt

    call grd(fid,xout,2,res)
    ak2 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(B31*grdt+B32*ak2)

    call grd(fid,xout,2,res)
    ak3 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(B41*grdt+B42*ak2+B43*ak3)

    call grd(fid,xout,2,res)
    ak4 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(B51*grdt+B52*ak2+B53*ak3+B54*ak4)

    call grd(fid,xout,2,res)
    ak5 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(B61*grdt+B62*ak2+B63*ak3+B64*ak4+B65*ak5)

    call grd(fid,xout,2,res)
    ak6 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(C1*grdt+C3*ak3+C4*ak4+C6*ak6)
    xerr = h0*(DC1*grdt+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)

  end subroutine stepper_rkck

  !> Doermand-Prince embedded 4(5)-order, local extrapolation.
  subroutine stepper_dp(fid,xpoint,grdt,h0,xout,xerr,res)
    use fields
    use global
    use types
    use param
    
    type(field), intent(inout) :: fid
    real*8, intent(in) :: xpoint(3), grdt(3), h0
    real*8, intent(out) :: xout(3), xerr(3)
    type(scalar_value), intent(inout) :: res

    real*8, parameter :: dp_a(7,7) = reshape((/&
       0.0d0,  0d0,0d0,0d0,0d0,0d0,0d0,&
       1d0/5d0,         0.0d0,0d0,0d0,0d0,0d0,0d0,&
       3d0/40d0,        9d0/40d0,       0.0d0,0d0,0d0,0d0,0d0,&
       44d0/45d0,      -56d0/15d0,      32d0/9d0,        0d0,0d0,0d0,0d0,&
       19372d0/6561d0, -25360d0/2187d0, 64448d0/6561d0, -212d0/729d0,  0d0,0d0,0d0,&
       9017d0/3168d0,  -355d0/33d0,     46732d0/5247d0,  49d0/176d0,  -5103d0/18656d0, 0d0,0d0,&
       35d0/384d0,      0d0,            500d0/1113d0,    125d0/192d0, -2187d0/6784d0,  11d0/84d0,      0d0&
       /),shape(dp_a))
    real*8, parameter :: dp_b2(7) = (/5179d0/57600d0, 0d0, 7571d0/16695d0, 393d0/640d0,&
       -92097d0/339200d0, 187d0/2100d0, 1d0/40d0/)
    real*8, parameter :: dp_b(7) = (/ 35d0/384d0, 0d0, 500d0/1113d0, 125d0/192d0, &
       -2187d0/6784d0, 11d0/84d0, 0d0 /)
    real*8, parameter :: dp_c(7) = dp_b2 - dp_b
    real*8, dimension(3) :: ak2, ak3, ak4, ak5, ak6, ak7

    xout = xpoint + h0*dp_a(2,1)*grdt

    call grd(fid,xout,2,res)
    ak2 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_a(3,1)*grdt+dp_a(3,2)*ak2)

    call grd(fid,xout,2,res)
    ak3 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_a(4,1)*grdt+dp_a(4,2)*ak2+dp_a(4,3)*ak3)

    call grd(fid,xout,2,res)
    ak4 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_a(5,1)*grdt+dp_a(5,2)*ak2+dp_a(5,3)*ak3+dp_a(5,4)*ak4)

    call grd(fid,xout,2,res)
    ak5 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_a(6,1)*grdt+dp_a(6,2)*ak2+dp_a(6,3)*ak3+dp_a(6,4)*ak4+dp_a(6,5)*ak5)

    call grd(fid,xout,2,res)
    ak6 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_b(1)*grdt+dp_b(2)*ak2+dp_b(3)*ak3+dp_b(4)*ak4+dp_b(5)*ak5+dp_b(6)*ak6)

    call grd(fid,xout,2,res)
    ak7 = res%gf / (res%gfmod+VSMALL)
    xerr = h0*(dp_c(1)*grdt+dp_c(2)*ak2+dp_c(3)*ak3+dp_c(4)*ak4+dp_c(5)*ak5+dp_c(6)*ak6+dp_c(7)*ak7)
    xout = xout + xerr

  end subroutine stepper_dp

  !> Prune a gradient path. A gradient path is given by the n points
  !> in x (fractional coordinates) referred to crystal structure c.
  !> In output, the number of points in the path is reduced so that
  !> the distances between adjacent points are at least fprune. The
  !> reduced number of points is returned in n, and the fractional
  !> coordinates in x(:,1:n).
  subroutine prunepath(c,n,x,fprune)
    use struct_basic
    
    type(crystal), intent(in) :: c
    integer, intent(inout) :: n
    real*8, intent(inout) :: x(n,3)
    real*8, intent(in) :: fprune

    integer :: i, nn
    real*8 :: x0(3), d

    ! prune the path
    x0 = x(1,:)
    nn = 1
    do i = 1, n
       d = c%distance(x(i,:),x0)
       if (d > fprune) then
          nn = nn + 1
          x(nn,:) = x(i,:)
          x0 = x(i,:)
       end if
    end do
    n = nn

  end subroutine prunepath

end module navigation
