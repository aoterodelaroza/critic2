! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (tools_math) proc
  implicit none

contains

  !> Calculate the cross-correlation between two functions with
  !> triangle weight. For similarity measures between spectra, as
  !> proposed here:
  !>   de Gelder et al., J. Comput. Chem., 22 (2001) 273.
  module function crosscorr_triangle(h,f,g,l) result(dfg)
    use tools_io, only: ferror, faterr
    real*8, intent(in) :: h, l
    real*8, intent(in) :: f(:), g(:)
    real*8 :: dfg

    integer :: n, m, i
    real*8 :: w
    
    ! check that the lengths are equal
    n = size(f)
    if (size(g) /= n) call ferror('crosscorr_triangle','inconsistent spectra',faterr)

    ! calculate maximum displacement
    m = floor(l / h)
    if (m <= 0 .or. m >= n) &
       call ferror('crosscorr_triangle','incorrect triangle slope',faterr)

    ! calculate the cross-correlation integral
    dfg = 0d0
    do i = 0, m
       ! triangle weight
       w = max(1d0 - real(i,8) * h / l,0d0)

       ! cfg(r) = int (f(x) * g(x+r))
       dfg = dfg + sum(f(1:n-i) * g(i+1:n)) * w

       ! cfg(-r) = int (f(x) * g(x-r))
       if (i == 0) cycle
       dfg = dfg + sum(g(1:n-i) * f(i+1:n)) * w

    end do
    dfg = dfg * h**2

  end function crosscorr_triangle

  !> Gives a crystallographic to cartesian conversion matrix from
  !> the cell parameters using the Cholesky decomposition of the
  !> metric tensor. Input angles in degrees.
  module function m_x2c_from_cellpar(aal,bbl) result(mat)
    use tools_io, only: ferror, faterr

    real*8, intent(in) :: aal(3),bbl(3)
    real*8 :: mat(3,3)

    real*8 :: deg2rad, calphl, cbetal, cgamml, HH1, HH2, HH3
    real*8 :: fac

    deg2rad = atan(1d0) / 45d0
    calphl = cos(bbl(1) * deg2rad)
    cbetal = cos(bbl(2) * deg2rad)
    cgamml = cos(bbl(3) * deg2rad)
    HH1 = sqrt(max(1d0-cgamml*cgamml,0d0))
    HH2 = (calphl-cbetal*cgamml) / HH1
    fac = 1d0 - cbetal*cbetal - HH2*HH2
    if (fac < 0d0) &
       call ferror("m_x2c_from_cellpar","invalid cell",faterr)
    HH3 = sqrt(max(1d0 - cbetal*cbetal - HH2*HH2,0d0))

    mat = 0d0
    mat(1,1) = aal(1)
    mat(1,2) = aal(2) * cgamml
    mat(1,3) = aal(3) * cbetal
    mat(2,2) = aal(2) * HH1
    mat(2,3) = aal(3) * HH2
    mat(3,3) = aal(3) * HH3

  end function m_x2c_from_cellpar

  !> Gives a cartesian to crystallographic conversion matrix from
  !> the cell parameters using the Cholesky decomposition of the
  !> metric tensor.
  module function m_c2x_from_cellpar(aal,bbl) result(mat)
    real*8, intent(in) :: aal(3),bbl(3)
    real*8 :: mat(3,3)

    real*8 :: deg2rad, calphl, cbetal, cgamml, HH1, HH2, HH3

    deg2rad = atan(1d0) / 45d0
    calphl = cos(bbl(1) * deg2rad)
    cbetal = cos(bbl(2) * deg2rad)
    cgamml = cos(bbl(3) * deg2rad)
    HH1 = sqrt(1d0-cgamml*cgamml)
    HH2 = (calphl-cbetal*cgamml) / sqrt(1d0-cgamml*cgamml)
    HH3 = sqrt(1d0 - cbetal*cbetal - HH2*HH2)

    mat = 0d0
    mat(1,1) = 1d0 / aal(1)
    mat(1,2) = - cgamml / (aal(1) * HH1)
    mat(1,3) = (HH2 * cgamml - HH1 * cbetal) / (aal(1) * HH1 * HH3)
    mat(2,2) = 1d0 / (aal(2) * HH1)
    mat(2,3) = - HH2 / (aal(2) * HH1 * HH3)
    mat(3,3) = 1d0 / (aal(3) * HH3)

  end function m_c2x_from_cellpar

  !> Factorial of an integer (returns real)
  module function factorial(n) result(f)
    use param, only: mfact, fact
    use tools_io, only: ferror, faterr

    real*8 :: f
    integer, intent(in) :: n

    if (n < 0 .or. n > mfact) then
       call ferror('factorial','Invalid argument',faterr)
    else
       f = fact(n)
    end if

  end function factorial

  !> Generates real regular solid harmonics from the spherical
  !> harmonics (Y_lm).  The index rlm(j) with j=l(l+1)+m+1. For a
  !> given (l,m) pair, if m>0, then rlm(j) is the cosine real regular
  !> solid harmonic:
  !>    C_lm = 1/sqrt(2) * ((-1)^m * Rlm + Rl-m)
  !> and if m<0 (abs(m) in the equation), rlm(j) is the sine real
  !> regular solid harmonic:
  !>    S_lm = 1/sqrt(2) * (-(-1)^m * i * Rlm + i * Rl-m)
  !> For m = 0, C_lm = R_lm and S_lm = 0. The Rlm are the regular
  !> solid harmonics calculated from the Ylm using:
  !>    R_lm = sqrt(4pi/(2l+1)) * r^l * Y_lm
  !>
  !> The order in the output is:
  !>  rrlm(:) = (/ C00, C11, C10, S11, C22, C21, C20, S21, S22,... /)
  !> and the first few are:
  !> C00 = 1 ; C11 = x ; C10 = z ; S11 = y
  !> C22 = sqrt(3)/2*(x^2-y^2) ; C21 = sqrt(3)*x*z 
  !> C20 = 1/2*(3*z^2-r^2)
  !> S21 = sqrt(3)*y*z ; S22 = sqrt(3)*x*y
  !>
  module subroutine genrlm_real(lmax,r,tp,rrlm)
    use param, only: pi, img

    integer, intent(in) :: lmax !< maximum angular momentum
    real*8, intent(in) :: r !< distance to origin
    real*8, intent(in) :: tp(2) !< (theta,phi) angles
    real*8, intent(out) :: rrlm((lmax+1)*(lmax+1)) !< array of real solid harmonics, (lmax+1)**2

    complex*16 :: rlm((lmax+1)*(lmax+1))
    integer :: l, imin, imax, m, ip, im, iphas
    real*8, parameter :: sh = 1d0/sqrt(2d0)

    call genylm(lmax,tp,rlm)
    do l = 0, lmax
       ! solid harmonics
       imin = l*l+1
       imax = (l+1)*(l+1)
       rlm(imin:imax) = rlm(imin:imax) * sqrt(4d0*pi/(2*l+1)) * r**l
       
       ! m = 0
       ip = imin+l
       rrlm(ip) = real(rlm(ip),8)

       ! m > 0
       do m = 1, l
          ip = imin + l + m 
          im = imin + l - m 
          iphas = (-1)**m
          rrlm(im) = sh * real(iphas*rlm(ip)+rlm(im),8)
          rrlm(ip) = sh * real(-iphas*img*rlm(ip)+img*rlm(im),8)
       end do
    end do
    
  end subroutine genrlm_real

  !> Generates a sequence of spherical harmonics, including the Condon-Shortley
  !> phase, evaluated at angles (theta,phi) for 0<l<lmax. The values
  !> are returned in a packed array ylm indexed with j=l(l+1)+m+1. The
  !> algorithm of Masters and Richards-Dinger is used, Geophys. J. Int.
  !> 135, 307 (1998). This routine is numerically stable and accurate to
  !> near machine precision for l<=50.
  module subroutine genylm(lmax,tp,ylm)

    integer, intent(in) :: lmax !< maximum angular momentum
    real*8, intent(in) :: tp(2) !< (theta,phi) coordinates
    complex*16, intent(out) :: ylm(*) !< array of spherical harmonics length=(lmax+1)**2

    ! local variables
    integer l,m,lm1,lm2
    real*8, parameter :: fourpi=12.566370614359172954d0
    real*8 sn,cs,dx,sum,t1
    ! automatic arrays
    real*8 x(0:lmax)
    complex*16 z(lmax)
    if ((lmax.lt.0).or.(lmax.gt.50)) then
       write(*,*)
       write(*,'("Error(genylm): lmax out of range : ",I8)') lmax
       write(*,*)
       stop
    end if
    ylm(1)=0.28209479177387814347d0
    if (lmax.eq.0) return
    sn=sin(tp(1))
    cs=cos(tp(1))
    ! phase factors exp(i*m*phi)
    do m=1,lmax
       t1=dble(m)*tp(2)
       z(m)=cmplx(cos(t1),sin(t1),8)
    end do
    do l=1,lmax
       if (mod(l,2).eq.0) then
          x(l)=1.d0
       else
          x(l)=-1.d0
       end if
       ! recursion loop
       dx=0.d0
       do m=l,1,-1
          t1=sqrt(dble((l+m)*(l-m+1)))
          x(m-1)=-(sn*dx+dble(2*m)*cs*x(m))/t1
          dx=sn*x(m)*t1
       end do
       ! rescale values and multiply with phase factors
       t1=sn
       sum=0.d0
       do m=1,l
          x(m)=t1*x(m)
          sum=sum+x(m)**2
          t1=t1*sn
       end do
       sum=2.d0*sum+x(0)**2
       t1=sqrt(dble(2*l+1)/(fourpi*sum))
       lm1=l*(l+1)+1
       lm2=lm1
       ylm(lm1)=t1*x(0)
       do m=1,l
          lm1=lm1+1
          lm2=lm2-1
          ylm(lm1)=t1*x(m)*z(m)
          ylm(lm2)=conjg(ylm(lm1))
          if (mod(m,2).ne.0) ylm(lm2)=-ylm(lm2)
       end do
    end do
    return
  end subroutine genylm

  !> Convert the input vector v to spherical coordinates. 
  !> v = (r*sin(th)*cos(ph),r*sin(th)*sin(ph),r*cos(th))
  module subroutine tosphere(v,r,tp)
    use param, only: pi
    real*8, intent(in) :: v(3) !< input vector
    real*8, intent(out) :: r !< distance
    real*8, intent(out) :: tp(2) !< (theta,phi)

    real*8, parameter :: eps=1.d-14
    real*8 :: t1

    r=sqrt(v(1)**2+v(2)**2+v(3)**2)
    tp = 0d0
    if (r > eps) then
       t1=v(3)/r
       if (t1 >= 1.d0) then
          tp(1)=0.d0
       else if (t1 <= -1.d0) then
          tp(1)=pi
       else
          tp(1)=acos(t1)
       end if
       if ((abs(v(1)) > eps).or.(abs(v(2)) > eps)) then
          tp(2)=atan2(v(2),v(1))
       end if
    end if

  end subroutine tosphere

  !> Calculate the derivatives of f(r)* Ylm. yl contains the list of
  !> spherical harmonics in standard (genylm) order. r is the distance
  !> to the origin, l and m are the angular momentum numbers, c, cp,
  !> and cpp are the value of f(r) and its first and second
  !> derivative. On output, grad and hess are the gradient and the
  !> hessian of f(r)*Ylm.
  module subroutine ylmderiv(yl,r,l,m,c,cp,cpp,grad,hess)
    use param, only: half

    complex*16, dimension(:), intent(in) :: yl
    real*8, intent(in) :: r
    integer, intent(in) :: l, m
    real*8, intent(in) :: c, cp, cpp
    complex*16, intent(out) :: grad(3), hess(6)

    !.Obtains derivatives of f(r) Ylm. See formulas used in notes and
    ! varshalovic.
    complex*16 :: dm, d0, dp
    complex*16 :: d00, dp0, dm0, dmp, dmm, dpp
    real*8 :: lpm,lmm,lpmm1,lpmm2,lpmm3,lpmp1,lpmp2,lpmp3
    real*8 :: lpmp4,lmmm1,lmmm2,lmmm3,lmmp1,lmmp2,lmmp3,lmmp4
    integer :: mm1, mm2, mp1, mp2, lp1, lp2, lm1, lm2
    real*8  :: sqrt2, r1, r2, fcoef1, fcoef2, d1, d2, d3, c1, c2, fden1, fden2
    real*8  :: dcoef1, dcoef2, dcoef3, dden1, dden2, dden3

    integer :: elem
    logical :: good

    ! inline functions ----
    elem(l,m)=l*(l+1)+m+1
    good(l,m)=l.ge.abs(m) .and. l.ge.0
    ! ---------------------

    sqrt2=1d0/sqrt(2d0)
    r1=1d0/r
    r2=r1*r1
    !.first derivatives
    fcoef1=(l+1)*r1*c+cp
    fcoef2=   -l*r1*c+cp
    fden1=0d0
    dden1=0d0
    if (l.gt.0) fden1=sqrt((2*l+1d0)*(2*l-1d0))
    fden2=sqrt((2*l+1d0)*(2*l+3d0))

    !.second derivatives
    dcoef1=(l*l-1)*r2*c + (2*l+1)*r1*cp +  cpp
    dcoef2=l*(l+1)*r2*c -     2d0*r1*cp -  cpp
    dcoef3=l*(l+2)*r2*c - (2*l+1)*r1*cp +  cpp

    if (l.ge.2) dden1=sqrt((2*l-3d0)*(2*l-1d0)**2*(2*l+1d0))
    dden2=(2*l-1) * (2*l+3)
    dden3=sqrt((2*l+1d0)*(2*l+3d0)**2 *(2*l+5d0))

    c1=0d0
    d1=0d0

    if (l.gt.0) c1=fcoef1/fden1
    c2=fcoef2/fden2

    if (l.ge.2) d1=dcoef1/dden1
    d2=dcoef2/dden2
    d3=dcoef3/dden3

    lm1=l-1
    lp1=l+1
    mm1=m-1
    mp1=m+1
    lm2=l-2
    lp2=l+2
    mm2=m-2
    mp2=m+2

    lpm=dble(l+m)
    lmm=dble(l-m)
    lpmm1=lpm-1d0
    lpmm2=lpm-2d0
    lpmm3=lpm-3d0
    lpmp1=lpm+1d0
    lpmp2=lpm+2d0
    lpmp3=lpm+3d0
    lpmp4=lpm+4d0
    lmmm1=lmm-1d0
    lmmm2=lmm-2d0
    lmmm3=lmm-3d0
    lmmp1=lmm+1d0
    lmmp2=lmm+2d0
    lmmp3=lmm+3d0
    lmmp4=lmm+4d0

    !.first derivatives
    ! follow varshalovich
    dm=(0d0,0d0)
    if (good(lm1,mm1)) dm=dm-sqrt2*sqrt(lpmm1*lpm  )*c1*yl(elem(lm1,mm1))
    if (good(lp1,mm1)) dm=dm+sqrt2*sqrt(lmmp1*lmmp2)*c2*yl(elem(lp1,mm1))

    d0=(0d0,0d0)
    if (good(lm1,m)) d0=d0+      sqrt(lpm*  lmm  )*c1*yl(elem(lm1,m  ))
    if (good(lp1,m)) d0=d0+      sqrt(lmmp1*lpmp1)*c2*yl(elem(lp1,m  ))

    dp=(0d0,0d0)
    if (good(lm1,mp1)) dp=dp-sqrt2*sqrt(lmmm1*lmm  )*c1*yl(elem(lm1,mp1))
    if (good(lp1,mp1)) dp=dp+sqrt2*sqrt(lpmp1*lpmp2)*c2*yl(elem(lp1,mp1))

    grad(1)=sqrt2*(dm-dp)
    grad(2)=(0d0,1d0)*sqrt2*(dm+dp)
    grad(3)=d0
    !.End first derivatives

    !.Second derivatives
    dm0=(0d0,0d0)
    if (good(lm2,mm1)) dm0=dm0-sqrt2*sqrt(lmm*lpm*lpmm1*lpmm2)*d1*yl(elem(lm2,mm1))
    if (good(l  ,mm1)) dm0=dm0+sqrt2*(2*m-1)*sqrt(lpm*lmmp1)*  d2*yl(elem(l  ,mm1))
    if (good(lp2,mm1)) dm0=dm0+sqrt2*sqrt(lmmp1*lpmp1*lmmp2*lmmp3)*d3*yl(elem(lp2,mm1))

    dp0=(0d0,0d0)
    if (good(lm2,mp1)) dp0=dp0-sqrt2*sqrt(lpm*lmm*lmmm1*lmmm2)*d1*yl(elem(lm2,mp1))
    if (good(l  ,mp1)) dp0=dp0-sqrt2*(2*m+1)*sqrt(lmm*lpmp1)*  d2*yl(elem(l  ,mp1))
    if (good(lp2,mp1)) dp0=dp0+sqrt2*sqrt(lpmp1*lmmp1*lpmp2*lpmp3)*d3*yl(elem(lp2,mp1))

    d00=(0d0,0d0)
    if (good(lm2,m)) d00=d00+sqrt(lmm*lmmm1*lpm*lpmm1)*d1*yl(elem(lm2,m))
    if (good(l,m)) d00=d00-(2*l*l+2*l-2*m*m-1)*d2*yl(elem(l,m))
    if (good(lp2,m)) d00=d00+sqrt(lmmp1*lmmp2*lpmp1*lpmp2)*d3*yl(elem(lp2,m))

    dmm=(0d0,0d0)
    if (good(lm2,mm2)) dmm=dmm+half*sqrt(lpm*lpmm1*lpmm2*lpmm3)*d1*yl(elem(lm2,mm2))
    if (good(l,mm2)) dmm=dmm+sqrt(lmmp1*lmmp2*lpm*lpmm1)*d2*yl(elem(l,mm2))
    if (good(lp2,mm2)) dmm=dmm+half*sqrt(lmmp1*lmmp2*lmmp3*lmmp4)*d3*yl(elem(lp2,mm2))

    dmp=(0d0,0d0)
    if (good(lm2,m)) dmp=dmp+half*sqrt(lpm*lpmm1*lmm*lmmm1)*d1*yl(elem(lm2,m))
    if (good(l,m)) dmp=dmp+(l*l+l+m*m-1)*d2*yl(elem(l,m))
    if (good(lp2,m)) dmp=dmp+half*sqrt(lmmp1*lmmp2*lpmp1*lpmp2)*d3*yl(elem(lp2,m))

    dpp=(0d0,0d0)
    if (good(lm2,mp2)) dpp=dpp+half*sqrt(lmm*lmmm1*lmmm2*lmmm3)*d1*yl(elem(lm2,mp2))
    if (good(l,mp2)) dpp=dpp+sqrt(lpmp1*lpmp2*lmm*lmmm1)*d2*yl(elem(l,mp2))
    if (good(lp2,mp2)) dpp=dpp+half*sqrt(lpmp1*lpmp2*lpmp3*lpmp4)*d3*yl(elem(lp2,mp2))

    hess(1)=          half*(dmm-dmp-dmp+dpp)
    hess(2)=(0d0,1d0)*half*(dmm        -dpp)
    hess(4)=         -half*(dmm+dmp+dmp+dpp)
    hess(6)=d00
    hess(3)=          sqrt2*(dm0-dp0)
    hess(5)=(0d0,1d0)*sqrt2*(dm0+dp0)

  end subroutine ylmderiv

  !> Given the value of a function (rlm) on an exponential grid defined
  !> by r_i = a * exp((i-1)*b), calculate the interpolated value (rho)
  !> and the first (rho1) and second (rho2) derivative at the
  !> distance r0. Does not apply to grid1_interp.
  module subroutine radial_derivs (rlm,rho,rho1,rho2,r0,a,b)
    real*8, dimension(:), intent(in) :: rlm
    real*8, intent(out) :: rho, rho1, rho2
    real*8, intent(in) :: r0
    real*8, intent(in) :: a, b
    
    ! radial grid derivation formulas
    integer, parameter :: noef(6,3) = reshape((/&
       0,  1,  2,  3,  4,  5,&
       -2, -1,  0,  1,  2,  3,&
       -5, -4, -3, -2, -1,  0/),shape(noef)) !< Node offsets for 6-point derivation formulas.
    real*8, parameter :: coef1(6,3) = reshape((/&
       -274,  600, -600,  400,  -150,  24,&
       6,  -60,  -40,  120,   -30,   4,&
       -24, 150, -400,  600,  -600, 274 /),shape(coef1)) !< Coefficients for first derivative.
    real*8, parameter :: coef2(6,3) = reshape((/&
       225, -770,  1070,  -780,   305,   -50,&
       -5,   80,  -150,    80,    -5,     0,&
       -50,  305,  -780 ,  1070,  -770,   225/),shape(coef2)) !< Coefficients for second derivative.
    real*8, parameter :: fac1=1d0/120d0 !< Prefactor for first derivative.
    real*8, parameter :: fac2=2d0/120d0 !< Prefactor for second derivative.

    integer :: ir, temp_ir
    integer :: nr
    real*8 :: r, rn, rrlm(4,0:2)
    real*8, dimension(4,4) :: x1dr12
    integer :: i, j
    integer :: ii, jj, ic
    real*8, dimension(4) :: r1, dr1
    real*8 :: prod

    nr = size(rlm)
    rn = a * exp(real(nr-1,8) * b)
    r = max(r0,a)
    if (r >= rn) then
       rho = 0d0
       rho1 = 0d0
       rho2 = 0d0
       return
    end if
    ir = min(max(floor(log(r / a) / b + 1),1),nr)

    ! careful with grid limits.
    if (ir <= 2) then
       temp_ir = 2
    else if (ir >= (nr - 2)) then
       temp_ir = nr - 2
    else
       temp_ir = ir
    end if

    rrlm(:,0) = rlm(temp_ir-1:temp_ir+2)
    x1dr12 = 0d0
    do i = 1, 4
       ii = temp_ir - 2 + i
       if (ii <= 2) then
          ic = 1
       else if ( ii >= (nr-2)) then
          ic = 3
       else
          ic = 2
       end if

       rrlm(i,1:2) = 0d0
       do j = 1, 6
          jj = ii + noef(j,ic)
          rrlm(i,1) = rrlm(i,1) + coef1(j,ic) * rlm(jj)
          rrlm(i,2) = rrlm(i,2) + coef2(j,ic) * rlm(jj)
       end do
       rrlm(i,1) = rrlm(i,1) * fac1
       rrlm(i,2) = rrlm(i,2) * fac2

       ! calculate factors and distances
       r1(i) = a*exp((ii-1)*b)
       dr1(i) = r - r1(i)
       do j = 1, i-1
          x1dr12(i,j) = 1.d0 / (r1(i) - r1(j))
          x1dr12(j,i) = -x1dr12(i,j)
       end do
    end do

    ! interpolate, lagrange 3rd order, 4 nodes
    rho = 0.d0
    rho1 = 0.d0
    rho2 = 0.d0
    do i = 1, 4
       prod = 1.d0
       do j = 1 ,4
          if (i == j) then
             cycle
          end if
          prod = prod * dr1(j) * x1dr12(i,j)
       end do
       rho = rho + rrlm(i,0) * prod
       rho1 = rho1 + rrlm(i,1) * prod
       rho2 = rho2 + rrlm(i,2) * prod
    end do

    rho2 = rho2 / (b * r)**2 - rho1 / b / r**2
    rho1 = rho1 / b / r

  end subroutine radial_derivs

  !> Compute x**i assuming that the value is:
  !>   (a) x**0 = 1D0  no matter what x is.
  !>   (b) 0D0**I = 0D0  for any I<>0.
  !>   (c) X**I  otherwise.
  module function ep(x,i)
    integer, intent(in) :: i !< Exponent
    real*8, intent(in) ::  x !< Base
    real*8 :: ep

    if (i.eq.0) then
       ep=1d0
    else if (x.eq.0d0) then
       ep=0d0
    else
       ep=x**i
    endif
    return
  end function ep

  !> Compute the greatest common divisor of an array of num integers (n). 
  module function gcdn(n,num)
    integer, intent(in) :: n(num)
    integer, intent(in) :: num
    integer :: gcdn

    integer :: i

    gcdn = 0
    if (num <= 0) return
    if (num == 1) then
       gcdn = n(1)
       return
    end if

    gcdn = gcd2(n(1),n(2))
    do i = 3, num
       gcdn = gcd2(gcdn,n(i))
    end do

  end function gcdn
  
  !> Compute greatest common divisor of two integers.
  module function gcd2(m,n)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer :: gcd2

    integer :: mhi, mlo, q

    gcd2 = 1

    mhi = max(m,n)
    mlo = min(m,n)
    do while (.true.)
       q = mod(mhi,mlo)
       if (q == 0) exit
       mhi = mlo
       mlo = q
    end do
    gcd2 = mlo

  end function gcd2
  
  !> Compute the least common multiple of an array of num integers (n). 
  module function lcmn(n,num)
    integer, intent(in) :: n(num)
    integer, intent(in) :: num
    integer :: lcmn

    integer :: i

    lcmn = 0
    if (num <= 0) return
    if (num == 1) then
       lcmn = n(1)
       return
    end if

    lcmn = lcm2(n(1),n(2))
    do i = 3, num
       lcmn = lcm2(lcmn,n(i))
    end do

  end function lcmn

  !> Compute the least common multiple of two integers.
  module function lcm2(m,n)
    integer, intent(in) :: m, n
    integer :: lcm2

    lcm2 = m * (n / gcd2(m,n))

  end function lcm2

  !> Compute the greatest common divisor of an array of num integers (n). 
  !> integer*8 version.
  module function gcdn_i8(n,num)
    integer*8, intent(in) :: n(num)
    integer, intent(in) :: num
    integer*8 :: gcdn_i8

    integer :: i

    gcdn_i8 = 0
    if (num <= 0) return
    if (num == 1) then
       gcdn_i8 = n(1)
       return
    end if

    gcdn_i8 = gcd2_i8(n(1),n(2))
    do i = 3, num
       gcdn_i8 = gcd2_i8(gcdn_i8,n(i))
    end do

  end function gcdn_i8
  
  !> Compute greatest common divisor of two integers.  integer*8
  !> version.
  module function gcd2_i8(m,n)
    integer*8, intent(in) :: m
    integer*8, intent(in) :: n
    integer*8 :: gcd2_i8

    integer*8 :: mhi, mlo, q

    gcd2_i8 = 1

    mhi = max(m,n)
    mlo = min(m,n)
    do while (.true.)
       q = mod(mhi,mlo)
       if (q == 0) exit
       mhi = mlo
       mlo = q
    end do
    gcd2_i8 = mlo

  end function gcd2_i8
  
  !> Compute the least common multiple of an array of num integers
  !> (n).  integer*8 version.
  module function lcmn_i8(n,num)
    integer*8, intent(in) :: n(num)
    integer, intent(in) :: num
    integer*8 :: lcmn_i8

    integer :: i

    lcmn_i8 = 0
    if (num <= 0) return
    if (num == 1) then
       lcmn_i8 = n(1)
       return
    end if

    lcmn_i8 = lcm2_i8(n(1),n(2))
    do i = 3, num
       lcmn_i8 = lcm2_i8(lcmn_i8,n(i))
    end do

  end function lcmn_i8

  !> Compute the least common multiple of two integers.  integer*8
  !> version.
  module function lcm2_i8(m,n)
    integer*8, intent(in) :: m, n

    lcm2_i8 = m * (n / gcd2_i8(m,n))

  end function lcm2_i8

  !> Find a rational number q/r that approximates x0 to within eps (r > 0).
  module subroutine rational_approx(x0,q,r,eps)
    real*8, intent(in) :: x0
    integer*8, intent(out):: q, r
    real*8, intent(in) :: eps

    real*8 :: x, diff, qr
    integer*8 :: nx, a, b, c, d

    ! convert to 0->1
    nx = floor(x0)
    x = x0 - nx

    ! set a and b and check if x0 is an integer
    a = 0
    b = 1
    diff = abs(real(a,8)/real(b,8) - x)
    if (diff < eps) then
       q = a
       r = b
       goto 999
    end if

    ! set c and d and check if x0 is an integer
    c = 1
    d = 1
    diff = abs(real(c,8)/real(d,8) - x)
    if (diff < eps) then
       q = c
       r = d
       goto 999
    end if

    ! bisect the fraction: a/b < (a+c)/(b+d) < c/d
    do while (.true.)
       q = a+c
       r = b+d

       qr = real(q,8) / real(r,8)
       diff = abs(qr - x)
       if (diff < eps) exit

       if (x > qr) then
          a = q
          b = r
       else
          c = q
          d = r
       end if
    end do

    ! recover the original value
999 continue
    q = r * nx + q

    return

  end subroutine rational_approx

  !> Find the lattice vector corresponding to the input real vector
  !> x0. If allowr is .true. yy can be either in x0 or -x0 direction,
  !> and the routine chooses the one with fewer minus signs.
  module function lattice_direction(x0,allowr) result(yy)
    real*8, intent(in) :: x0(3)
    logical, intent(in) :: allowr
    integer :: yy(3)

    integer*8 :: q(3), r(3), dl
    integer :: i, miny, nminus, nplus

    do i = 1, 3
       call rational_approx(x0(i),q(i),r(i),1d-5)
    end do
    dl = lcm(r,3)

    do i = 1, 3
       yy(i) = int(q(i) * (dl / r(i)),4)
    end do

    miny = maxval(abs(yy))
    do i = 1, 3
       if (yy(i) /= 0) then
          miny = min(miny,abs(yy(i)))
       end if
    end do
    yy = yy / miny

    if (allowr) then
       nminus = count(yy < 0)
       nplus = count(yy > 0)
       if (nminus > nplus) yy = - yy
    end if

  end function lattice_direction

  !> Find the eigenvectors and eigenvalues of a real nxn symmetric matrix.
  !> The eigenvectors are stored column-wise in mat.
  module subroutine eig(mat,eval)
    use tools_io, only: ferror, faterr
    ! diagonalize a real symmetric matrix
    real*8, intent(inout) :: mat(:,:) !< Input matrix and output eigenvectors (column-wise)
    real*8, intent(out), optional :: eval(:) !< Output eigenvalues

    real*8 :: work(size(mat,1))
    integer :: ifail
    real*8 :: loceval(size(mat,1))

    work = 0d0
    if (present(eval)) then
       call trace(mat,eval,work,size(mat,1),ifail)
    else
       call trace(mat,loceval,work,size(mat,1),ifail)
    end if

    if (ifail /= 0) call ferror('eig','Error in diagonalization',faterr)

  end subroutine eig

  !> Find the eigenvectors and eigenvalues of a real nxn matrix.
  !> The eigenvectors are stored column-wise in mat.
  module subroutine eigns(mat,eval,evali)
    use tools_io, only: ferror, faterr
    real*8, intent(inout) :: mat(:,:) !< Input matrix and output eigenvectors (column-wise)
    real*8, intent(out), optional :: eval(:) !< Real part of the eigenvalues
    real*8, intent(out), optional :: evali(:) !< Imaginary part of the eigenvalues

    integer :: ifail
    real*8, dimension(size(mat,1)) :: loceval, localevali, work1, work2, work3
    real*8, dimension(size(mat,1),size(mat,1)) :: mati, matt, matt2

    mati = 0d0

    if (present(eval)) then
       call cg(size(mat,1),size(mat,1),mat,mati,eval,evali,1,matt,matt2,work1,work2,work3,ifail)
    else
       call cg(size(mat,1),size(mat,1),mat,mati,loceval,localevali,1,matt,matt2,work1,work2,work3,ifail)
    end if
    mat = matt

    if (ifail /= 0) call ferror('eigns','Error in diagonalization',faterr)
  end subroutine eigns

  !> Given a point x0 (cartesian), calculate the rank and signature of
  !> the hessian of f. Adapted to the determination of CP type.  If
  !> the Hessian elements are less than eps (in absolute value), it is
  !> considered zero.
  module subroutine rsindex(mat,ehess,r,s,eps)
    real*8, intent(inout)  :: mat(3,3)
    real*8, intent(out) :: ehess(3)
    integer, intent(out) :: r, s
    real*8, intent(in) :: eps

    integer :: nhplus, nhminus, i

    call eig(mat,ehess)

    nhplus = 0
    nhminus = 0
    do i = 1, 3
       if (ehess(i) > eps) nhplus = nhplus+1
       if (ehess(i) < -eps) nhminus = nhminus+1
    end do
    r = nhplus+nhminus
    s = nhplus-nhminus

  end subroutine rsindex

  !> Cross product of two 3-vectors
  module function cross(v1,v2) result (vx)
    real*8, intent(in) :: v1(3) !< First vector
    real*8, intent(in) :: v2(3) !< Second vector
    real*8 :: vx(3)

    vx(1) = + v1(2)*v2(3) - v1(3)*v2(2)
    vx(2) = - v1(1)*v2(3) + v1(3)*v2(1)
    vx(3) = + v1(1)*v2(2) - v1(2)*v2(1)

  end function cross

  !> Mixed product of three 3-vectors
  module function mixed(v1,v2,v3)
    real*8, intent(in) :: v1(3) !< First vector
    real*8, intent(in) :: v2(3) !< Second vector
    real*8, intent(in) :: v3(3) !< Third vector
    real*8 :: mixed

    mixed = v1(1)*v2(2)*v3(3) - v1(1)*v2(3)*v3(2) - &
            v1(2)*v2(1)*v3(3) + v1(2)*v2(3)*v3(1) + &
            v1(3)*v2(1)*v3(2) - v1(3)*v2(2)*v3(1)

  end function mixed

  !> Norm-2 of a 3x3 matrix
  module function mnorm2(a)
    real*8, intent(in) :: a(3,3)
    real*8 :: mnorm2
    
    real*8 :: b(3,3), eval(3)

    b = matmul(transpose(a),a)
    call eig(b,eval)
    mnorm2 = sqrt(maxval(eval))
    
  end function mnorm2

  !> Determinant of a real 3x3 symmetric matrix
  module function detsym(m)
    real*8, intent(in) :: m(3,3) !< Input matrix
    real*8 :: detsym

    detsym = m(1,1) * m(2,2) * m(3,3) + &
       2d0 * m(1,2) * m(2,3) * m(1,3) - &
       m(1,3) * m(1,3) * m(2,2) - &
       m(2,3) * m(2,3) * m(1,1) - &
       m(1,2) * m(1,2) * m(3,3) 

  end function detsym

  !> Determinant of a real 3x3 matrix
  module function det(m)
    real*8, intent(in) :: m(3,3) !< Input matrix
    real*8 :: det

    det = m(1,1) * (m(2,2) * m(3,3) - m(2,3) * m(3,2)) + &
          m(1,2) * (m(2,3) * m(3,1) - m(2,1) * m(3,3)) + &
          m(1,3) * (m(2,1) * m(3,2) - m(2,2) * m(3,1))

  end function det

  !> Inverse of a 3x3 matrix
  module function matinv(m) result(mo)
    real*8, intent(in) :: m(3,3) !< Input matrix
    real*8 :: mo(3,3) !< Inverse of input

    mo(1,1) =   m(2,2)*m(3,3) - m(2,3)*m(3,2)
    mo(1,2) = - m(2,1)*m(3,3) + m(2,3)*m(3,1)
    mo(1,3) =   m(2,1)*m(3,2) - m(2,2)*m(3,1)
    mo(2,1) = - m(1,2)*m(3,3) + m(1,3)*m(3,2)
    mo(2,2) =   m(1,1)*m(3,3) - m(1,3)*m(3,1)
    mo(2,3) = - m(1,1)*m(3,2) + m(1,2)*m(3,1)
    mo(3,1) =   m(1,2)*m(2,3) - m(1,3)*m(2,2)
    mo(3,2) = - m(1,1)*m(2,3) + m(1,3)*m(2,1)
    mo(3,3) =   m(1,1)*m(2,2) - m(1,2)*m(2,1)
    mo = transpose(mo) / det(m)

  end function matinv

  !> Scale or extend the plane defined by points x0, x1, x2 (Cartesian)
  !> with scale factors sxi, syi (default: 1) and extend factors
  !> zx0i, zx1i in the x direction (default: 0) and zy0i, zy1i in the
  !> y direction (default: 0).
  module subroutine plane_scale_extend(x0,x1,x2,sxi,syi,zx0i,zx1i,zy0i,zy1i)
    use tools_io, only: ferror, faterr
    use param, only: VSMALL
    real*8, intent(inout) :: x0(3), x1(3), x2(3)
    real*8, intent(in), optional :: sxi, syi, zx0i, zx1i, zy0i, zy1i
    
    real*8 :: sx, sy, zx0, zx1, zy0, zy1
    real*8 :: ax(3), ay(3), ax0(3), ax1(3), ay0(3), ay1(3), dx, dy

    sx = 1d0
    sy = 1d0
    zx0 = 0d0
    zx1 = 0d0
    zy0 = 0d0
    zy1 = 0d0
    if (present(sxi)) sx = sxi
    if (present(syi)) sy = syi
    if (present(zx0i)) zx0 = zx0i
    if (present(zx1i)) zx1 = zx1i
    if (present(zy0i)) zy0 = zy0i
    if (present(zy1i)) zy1 = zy1i
    
    ax = (x1-x0)
    ay = (x2-x0)
    dx = norm2(ax)
    dy = norm2(ay)
    if (dx < VSMALL .or. dy < VSMALL) &
       call ferror('plane_scale_extend','zero-area plane',faterr)

    ax0 = -(ax/dx) * (0.5d0 * (sx-1d0) * dx + zx0)
    ax1 =  (ax/dx) * (0.5d0 * (sx-1d0) * dx + zx1)
    ay0 = -(ay/dy) * (0.5d0 * (sy-1d0) * dy + zy0)
    ay1 =  (ay/dy) * (0.5d0 * (sy-1d0) * dy + zy1)

    x0 = x0 + ax0 + ay0
    x1 = x1 + ax1 + ay0
    x2 = x2 + ax0 + ay1

  end subroutine plane_scale_extend

  !> This routine selects particular contour values (ziso(1:niso))
  !> based on the given selection scheme (niso_type). lin0 and lin1
  !> are the limits of the linear range fmax and fmin are the maximum
  !> and minimum values of the field. 
  module subroutine assign_ziso(niso_type,niso,ziso,lin0,lin1,fmax,fmin)
    use param, only: pi
    integer, intent(in) :: niso_type
    integer, intent(inout) :: niso
    real*8, allocatable, intent(inout) :: ziso(:)
    real*8, intent(in) :: lin0, lin1, fmax, fmin

    real*8 :: ffmin, ffmax, delta
    integer :: i
    integer :: nhalf
    real*8, parameter :: eps = 1d-12

    if (niso_type == niso_manual) then
       ! niso/ziso already given in input
    elseif (niso_type == niso_lin) then
       ! linear range
       if (allocated(ziso)) deallocate(ziso)
       allocate(ziso(niso))
       do i = 1, niso
          ziso(i) = lin0 + real(i-1,8) * (lin1 - lin0) / real(niso-1,8)
       end do
    elseif (niso_type == niso_log) then
       if (allocated(ziso)) deallocate(ziso)

       ! logarithmic range
       if (fmin < -eps) then
          ! negative contours
          ffmin = log(eps)
          ffmax = log(-fmin)
          nhalf = max(niso / 2,2)
          niso = 2 * nhalf + 1
          allocate(ziso(niso))
          delta = (ffmax-ffmin) / (nhalf-1)
          do i = 1, nhalf
             ziso(i) = -exp(ffmin+(nhalf-i)*delta)
          enddo

          ! zero
          ziso(nhalf+1) = 0d0

          ! positive contours
          ffmin = log(eps)
          ffmax = log(fmax)
          delta = (ffmax-ffmin) / (nhalf-1)
          do i = 1, nhalf
             ziso(nhalf+1+i) = exp(ffmin+(i-1)*delta)
          enddo
       else
          ! positive contours
          allocate(ziso(niso))
          ffmin = log(max(fmin,eps))
          ffmax = log(abs(fmax))
          delta = (ffmax-ffmin)/(niso-1)
          do i = 1, niso
             ziso(i)=exp(ffmin+(i-1)*delta)
          enddo
       endif
    elseif (niso_type == niso_atan) then
       ! arctangent mapping
       if (allocated(ziso)) deallocate(ziso)
       allocate(ziso(niso))
       ffmin = 2d0/pi*atan(fmin)
       ffmax = 2d0/pi*atan(fmax)
       delta = (ffmax-ffmin)/(niso-1)
       do i = 1, niso
          ziso(i) = tan(pi*(ffmin+(i-1)*delta)/2d0)
       enddo
    elseif (niso_type == niso_bader) then
       ! bader
       if (fmin < -eps) then
          niso = 40
       else
          niso = 20
       end if
       if (allocated(ziso)) deallocate(ziso)
       allocate(ziso(niso))
       ziso(1) = 1.d-3
       ziso(2) = 2.d-3
       ziso(3) = 4.d-3
       ziso(4) = 8.d-3
       ziso(5) = 1.d-2
       ziso(6) = 2.d-2
       ziso(7) = 4.d-2
       ziso(8) = 8.d-2
       ziso(9) = 1.d-1
       ziso(10)= 2.d-1
       ziso(11)= 4.d-1
       ziso(12)= 8.d-1
       ziso(13)= 1.d0 
       ziso(14)= 2.d0 
       ziso(15)= 4.d0 
       ziso(16)= 8.d0 
       ziso(17)= 1.d1 
       ziso(18)= 2.d1 
       ziso(19)= 4.d1 
       ziso(20)= 8.d1 
       if (fmin < -eps) then
          ziso(21:40) = ziso(1:20)
          ziso(1) = -8.d1
          ziso(2) = -4.d1
          ziso(3) = -2.d1
          ziso(4) = -1.d1
          ziso(5) = -8.d0
          ziso(6) = -4.d0
          ziso(7) = -2.d0
          ziso(8) = -1.d0
          ziso(9) = -8.d-1
          ziso(10)= -4.d-1
          ziso(11)= -2.d-1
          ziso(12)= -1.d-1
          ziso(13)= -8.d-2
          ziso(14)= -4.d-2
          ziso(15)= -2.d-2
          ziso(16)= -1.d-2
          ziso(17)= -8.d-3
          ziso(18)= -4.d-3
          ziso(19)= -2.d-3
          ziso(20)= -1.d-3
       end if
    end if

  end subroutine assign_ziso

  !> Generate the combination of n objects taken p at a time that
  !> corresponds to a given index in the lexicographic order. Buckles
  !> and Lybanon's algorithm (toms/515, B.P. Buckles and M.  Lybanon,
  !> ACM TOMS 3 (1977) 180-182).
  module subroutine comb(n, p, l, c)
    integer :: n, p, l, c(p)

    ! THIS SUBROUTINE FINDS THE COMBINATION SET OF N THINGS
    ! TAKEN P AT A TIME FOR A GIVEN LEXICOGRAPHICAL INDEX.
    ! N - NUMBER OF THINGS IN THE SET
    ! P - NUMBER OF THINGS IN EACH COMBINATION
    ! L - LEXICOGRAPHICAL INDEX OF COMBINATION SOUGHT
    ! C - OUTPUT ARRAY CONTAINING THE COMBINATION SET
    ! THE FOLLOWING RELATIONSHIPS MUST EXIST AMONG THE INPUT
    ! VARIABLES.  L MUST BE GREATER THAN OR EQUAL TO 1 AND LESS
    ! THAN OR EQUAL TO THE MAXIMUM LEXICOGRAPHICAL INDEX.
    ! P MUST BE LESS THAN OR EQUAL TO N AND GREATER THAN ZERO.
    INTEGER K, R, P1, I

    ! SPECIAL CASE CODE IF P = 1
    if (p <= 1) then
       C(1) = L
       return
    endif

    ! INITIALIZE LOWER BOUND INDEX AT ZERO
    K = 0
    ! LOOP TO SELECT ELEMENTS IN ASCENDING ORDER
    P1 = P - 1
    C(1) = 0
    do I=1,P1
       ! SET LOWER BOUND ELEMENT NUMBER FOR NEXT ELEMENT VALUE
       if (I.NE.1) C(I) = C(I-1)
       ! LOOP TO CHECK VALIDITY OF EACH ELEMENT VALUE
       do while (k < l)
          C(I) = C(I) + 1
          R = nchoosek(N-C(I),P-I)
          K = K + R
       end do
       K = K - R
    end do
    C(P) = C(P1) + L - K

  end subroutine comb

  !> Returns m choose n. Same source as comb.
  module function nchoosek(M, N)
    integer :: m, n
    integer :: nchoosek

    ! ACM ALGORITHM 160 TRANSLATED TO FORTRAN.  CALCULATES THE
    ! NUMBER OF COMBINATIONS OF M THINGS TAKEN N AT A TIME.
    INTEGER P, I, N1, R
    N1 = N
    P = M - N1
    IF (N1 < P) then
       P = N1
       N1 = M - P
    end IF
    R = N1 + 1
    IF (P.EQ.0) R = 1
    IF (P >= 2) then
       do I=2,P
          R = (R*(N1+I))/I
       end DO
    end IF
    nchoosek = R
  end function nchoosek


  !> Rotate and translate the set of points x1 (3,n) to match the set
  !> of points x2 (3,n), and calculate the rmsd of the resulting
  !> positions. This routine uses Walker's algorithm based on
  !> quaternion algebra. (Walker et al., CVGIP-Imag. Understan. 54
  !> (1991) 358.)
  module function rmsd_walker(x1o,x2o) result(rmsd)
    use tools_io, only: ferror, faterr
    real*8, intent(in) :: x1o(:,:), x2o(:,:)
    real*8 :: rmsd
    
    integer :: n, i, idx
    real*8, allocatable :: x1(:,:), x2(:,:)
    real*8 :: xcm1(3), xcm2(3), c1(4,4), c2, c3(4,4), xex(4)
    real*8 :: w(4,4), q(4,4), a(4,4), eval(4), evali(4)

    n = size(x1o,2)
    if (size(x1o,1) /= 3 .or. size(x2o,1) /= 3 .or. size(x2o,2) /= n) &
       call ferror("rmsd_walker","Inconsistent number of points",faterr)
    allocate(x1(3,n),x2(3,n))
    
    xcm1 = sum(x1o,2) / n
    xcm2 = sum(x2o,2) / n
    do i = 1, n
       x1(:,i) = x1o(:,i) - xcm1
       x2(:,i) = x2o(:,i) - xcm2
    end do

    ! The c1, c2, and c3 quantities (eqs. 35-37)
    c1 = 0d0
    c2 = 0.5d0 * real(n,8)
    c3 = 0d0
    xex = 0d0
    do i = 1, n
       xex(1:3) = x1(:,i)
       w = wmat(xex)
       xex(1:3) = x2(:,i)
       q = qmat(xex)
       c1 = c1 - matmul(transpose(q),w)
       c3 = c3 + (w - q)
    end do

    ! the a matrix (eq. 47), diagonalization, and rotation matrix
    a = matmul(transpose(c3),c3) * c2 - c1;
    call eigns(a,eval,evali)
    idx = maxloc(eval,1)
    xex = a(:,idx)
    xex = xex / norm2(xex)
    w = wmat(xex)
    q = qmat(xex)
    q = matmul(transpose(wmat(xex)),qmat(xex))

    ! rmsd and clean up
    x1 = matmul(q(1:3,1:3),x1)
    rmsd = sqrt(sum((x1 - x2)**2) / n)
    deallocate(x1,x2)

  contains
    ! the W and Q functions (eqs. 15 and 16)
    function wmat(x) result(w)
      real*8, intent(in) :: x(4)
      real*8 :: w(4,4)
      w(:,1) = (/  x(4), -x(3),  x(2), -x(1)/)
      w(:,2) = (/  x(3),  x(4), -x(1), -x(2)/)
      w(:,3) = (/ -x(2),  x(1),  x(4), -x(3)/)
      w(:,4) = (/  x(1),  x(2),  x(3),  x(4)/)
    end function wmat
    function qmat(x) result(q)
      real*8, intent(in) :: x(4)
      real*8 :: q(4,4)
      q(:,1) = (/  x(4),  x(3), -x(2), -x(1)/)
      q(:,2) = (/ -x(3),  x(4),  x(1), -x(2)/)
      q(:,3) = (/  x(2), -x(1),  x(4), -x(3)/)
      q(:,4) = (/  x(1),  x(2),  x(3),  x(4)/)
    end function qmat
  end function rmsd_walker

  !> Find the Gauss-Legendre nodes and weights for an interval.
  module subroutine gauleg(x1,x2,x,w,n)
    use param, only: pi
    real*8, intent(in) :: x1 !< Left limit of the interval
    real*8, intent(in) :: x2 !< Right limit of the interval
    real*8, dimension(n), intent(out) :: x !< Position of the nodes
    real*8, dimension(n), intent(out) :: w !< Weights
    integer, intent(in) :: n !< Number of nodes

    real*8, parameter :: eps = 3.0d-16

    integer :: m
    real*8 :: xm, xl
    real*8 :: z, p1, p2, p3, pp, z1
    integer :: i, j

    !.The roots are symmetric in the interval, so we only have to find
    ! half of them.

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)

    ! Loop over the desired roots.
    do i=1,m
       z=cos(pi*(i-.25d0)/(n+.5d0))
1      continue
       p1=1.d0
       p2=0.d0

       ! Loop up the recurrence relation to get the Legendre polynomial
       ! evaluated at z.
       do j=1,n
          p3=p2
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
       end do

       !.p1 is now the desired Legendre polynomial. We next compute pp,
       ! derivative , by a standard relation involving p2, the polyn-
       ! omial of one lower order.
       pp=n*(z*p1-p2)/(z*z-1.d0)
       z1=z

       !.Newton's method.
       !
       z=z1-p1/pp
       if(abs(z-z1).gt.eps)go to 1

       !.Scale the root to the desired interval.
       x(i)=xm-xl*z

       !.and put in its symmetric counterpart.
       x(n+1-i)=xm+xl*z

       !.compute the weight.
       w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)

       !.and its symmetric counterpart.
       w(n+1-i)=w(i)
    end do

  end subroutine gauleg

end submodule proc
