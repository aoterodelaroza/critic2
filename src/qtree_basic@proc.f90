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

submodule (qtree_basic) proc
  implicit none

  ! Butcher tableaus of ODE solvers
  ! c_1 |
  ! c_2 | a_21
  ! c_3 | a_31  a_32
  ! ... | ...   ...
  ! c_o | a_o1  a_o2  ...  a_oo-1
  !------------------------------
  !     | b_1   b_2        b_o
  !
  ! k_i = f(t_n+c_i*h, y_n+ a_i1 * h * k_1+...+ a_ii-1 * h * k_i-1)
  ! y_n+1 = y_n + h * sum_i=1..o (  b_i * k_i )
  !
  ! The step is calculated with the b (not b2) formula
  !
  !! Euler method 
  integer, target :: euler_o = 1
  real*8, target :: euler_c(1) = (/0d0/)
  real*8, target :: euler_a(1,1) = reshape((/1d0/),shape(euler_a))
  real*8, target :: euler_b(1) = (/1d0/)
  integer, target :: euler_type = 0
  !! Heun's method
  integer, target :: heun_o = 2
  real*8, target :: heun_c(2) = (/0d0, 1d0/)
  real*8, target :: heun_a(2,2) = reshape((/0d0, 1d0,&
     1d0, 0d0/),shape(heun_a))
  real*8, target :: heun_b(2) = (/0.5d0, 0.5d0/)
  integer, target :: heun_type = 0
  !! Kutta's method
  integer, target :: kutta_o = 3
  real*8, target :: kutta_c(3) = (/0d0, 0.5d0, 1d0/)
  real*8, target :: kutta_a(3,3) = reshape((/0.0d0, 0.5d0, -1.0d0,&
     0.5d0, 0.0d0,  2.0d0,&
     -1.0d0, 2.0d0,  0.0d0/),shape(kutta_a))
  real*8, target :: kutta_b(3) = (/1d0/6d0, 2d0/3d0, 1d0/6d0/)
  integer, target :: kutta_type = 0
  !! RK4 method
  integer, target :: rk4_o = 4
  real*8, target :: rk4_c(4) = (/0d0, 0.5d0, 0.5d0, 1d0/)
  real*8, target :: rk4_a(4,4) = reshape((/0.0d0, 0.5d0, 0.0d0, 0.0d0,&
     0.5d0, 0.0d0, 0.5d0, 0.0d0,&
     0.0d0, 0.5d0, 0.0d0, 1.0d0,&
     0.0d0, 0.0d0, 1.0d0, 0.0d0/),shape(rk4_a))
  real*8, target :: rk4_b(4) = (/1d0/6d0, 1d0/3d0, 1d0/3d0, 1d0/6d0/)
  integer, target :: rk4_type = 0
  !! Heun-Euler embedded method
  !! 1st order, 2nd order error est., 2 evals
  integer, target :: heul_o = 2
  integer, target :: heul_o2 = 1
  real*8, target :: heul_c(2) = (/0d0, 1d0/)
  real*8, target :: heul_a(2,2) = reshape((/0d0, 0d0,&
     1d0, 0d0/),shape(heul_a))
  real*8, target :: heul_b2(2) = (/0.5d0, 0.5d0/)
  real*8, target :: heul_b(2) = (/1.0d0, 0.0d0/) 
  integer, target :: heul_type = 1
  !! Bogacki-Shampine embedded method 
  !! 3th order, 5th order error est., 3(4) evals
  !! Local extrapolation, fsal
  integer, target :: bs_o = 4
  integer, target :: bs_o2 = 2
  real*8, target :: bs_c(4) = (/0d0, 0.5d0, 0.75d0, 1d0/)
  real*8, target :: bs_a(4,4) = reshape(&
     (/0.0d0,0d0,0d0,0d0,&
     1d0/2d0,        0.0d0,0d0,0d0,&
     0.0d0,          3d0/4d0,     0.0d0,0d0,&
     2d0/9d0,        1d0/3d0,     4d0/9d0,       0.0d0&
     /),shape(bs_a))
  real*8, target :: bs_b(4) = (/2d0/9d0, 1d0/3d0, 4d0/9d0, 0d0/)
  real*8, target :: bs_b2(4) = (/7d0/24d0, 1d0/4d0, 1d0/3d0, 1d0/8d0/)
  integer, target :: bs_type = 1
  !! Runge-Kutta-Cash-Karp embedded method 
  !! 4th order, 5th order error est., 6 evals
  integer, target :: rkck_o = 6
  integer, target :: rkck_o2 = 4
  real*8, target :: rkck_c(6) = (/0d0, 1d0/5d0, 3d0/10d0, 3d0/5d0, 1d0, 7d0/8d0/)
  real*8, target :: rkck_a(6,6) = reshape(&
     (/0.0d0,0d0,0d0,0d0,0d0,0d0,&
     1d0/5d0,        0.0d0,0d0,0d0,0d0,0d0,&
     3d0/40d0,       9d0/40d0,    0.0d0,0d0,0d0,0d0,&
     3d0/10d0,      -9d0/10d0,    6d0/5d0,       0.0d0,0d0,0d0,&
     -11d0/54d0,      5d0/2d0,    -70d0/27d0,     35d0/27d0,        0.0d0,0d0,&
     1631d0/55296d0, 175d0/512d0, 575d0/13824d0, 44275d0/110592d0, 253d0/4096d0, 0.0d0&
     /),shape(rkck_a))
  real*8, target :: rkck_b(6) = (/2825d0/27648d0, 0d0, 18575d0/48384d0, 13525d0/55296d0,&
     277d0/14336d0, 1d0/4d0 /)
  real*8, target :: rkck_b2(6) = (/ 37d0/378d0, 0d0, 250d0/621d0, 125d0/594d0, 0d0, 512d0/1771d0 /)
  integer, target :: rkck_type = 1
  !! Dormand-Prince embedded method 
  !! 4th order, 5th order error est., 6(7) evals
  !! Local extrapolation, fsal
  integer, target :: dp_o = 7
  integer, target :: dp_o2 = 4
  real*8, target :: dp_c(7) = (/0d0, 1d0/5d0, 3d0/10d0, 4d0/5d0, 8d0/9d0, 1d0, 1d0/)
  real*8, target :: dp_a(7,7) = reshape(&
     (/0.0d0,  0d0,0d0,0d0,0d0,0d0,0d0,&
     1d0/5d0,         0.0d0,0d0,0d0,0d0,0d0,0d0,&
     3d0/40d0,        9d0/40d0,       0.0d0,0d0,0d0,0d0,0d0,&
     44d0/45d0,      -56d0/15d0,      32d0/9d0,        0d0,0d0,0d0,0d0,&
     19372d0/6561d0, -25360d0/2187d0, 64448d0/6561d0, -212d0/729d0,  0d0,0d0,0d0,&
     9017d0/3168d0,  -355d0/33d0,     46732d0/5247d0,  49d0/176d0,  -5103d0/18656d0, 0d0,0d0,&
     35d0/384d0,      0d0,            500d0/1113d0,    125d0/192d0, -2187d0/6784d0,  11d0/84d0,      0d0&
     /),shape(dp_a))
  real*8, target :: dp_b2(7) = (/5179d0/57600d0, 0d0, 7571d0/16695d0, 393d0/640d0,&
     -92097d0/339200d0, 187d0/2100d0, 1d0/40d0/)
  real*8, target :: dp_b(7) = (/ 35d0/384d0, 0d0, 500d0/1113d0, 125d0/192d0, &
     -2187d0/6784d0, 11d0/84d0, 0d0 /)
  integer, target :: dp_type = 1

contains

  !> Initialize qtree, at level lvl and pre-split level plvl. This routine is called once
  !> at the beginning of the integration, or twice if the beta-sphere sizes are determined
  !> dynamically.
  module subroutine qtree_initialize(lvl,plvl,acum_atprop,trm,fgr,lapgr,vgr,verbose)
    use systemmod, only: sy
    use global, only: minl, prop_mode, integ_scheme, integ_mode, keastnum,&
       qtree_ode_mode, color_allocate, plot_mode, docontacts, ws_use, ws_origin,&
       ws_scale
    use tools_math, only: mixed, cross
    use tools_io, only: ferror, faterr, uout, warning
    use param, only: eye
    integer, intent(in) :: lvl, plvl
    logical, intent(in) :: verbose
    integer(qtreei), allocatable, intent(out) :: trm(:,:)
    real(qtreer), allocatable, intent(out) :: fgr(:,:), lapgr(:,:), vgr(:)
    real*8, allocatable, intent(out) :: acum_atprop(:,:)

    integer :: i, j, ntetrag
    real*8 :: xp1(3), xp2(3), xp3(3), xx(3), dist, iw(3), r1(3), er
    real*8 :: vtotal, sumi
    integer :: l2
    integer(qtreeidx) :: siz
    character*3 :: pg
    real*8, allocatable :: tetrag(:,:,:)

    ! count non-equivalent maxima (nuclear and non-nuclear)
    ! the list is ordered with maxima first.
    nnuc = 0
    do i = 1, sy%f(sy%iref)%ncp
       if (sy%f(sy%iref)%cp(i)%typ == sy%f(sy%iref)%typnuc) nnuc = nnuc + 1
    end do

    ! set the max. level with user input
    maxl = lvl
    if (maxl <= minl) call ferror('qtree_sphere_radii','maxl <= minl',faterr)
    l2 = 2**maxl

    ! interpret integration scheme
    if (prop_mode == 0) integ_scheme = 0
    if (integ_scheme == 0) then
       integ_mode = 0
       prop_mode = 0
    else if (integ_scheme == 1) then
       integ_mode = -1
       integ_mode(maxl) = 11
    else if (integ_scheme == 2) then
       integ_mode = -1
       integ_mode(maxl) = 1
    else if (integ_scheme == 3) then
       integ_mode = 1
    else if (integ_scheme == 4) then
       integ_mode = keastnum
    else if (integ_scheme == 5) then
       integ_mode = 12
    else if (integ_scheme == 6) then
       do i = minl+1,min(maxl,6)
          integ_mode(i) = 12
       end do
       do i = 7,maxl
          integ_mode(i) = -1
          if (i == maxl) integ_mode(i) = 11
       end do
    else if (integ_scheme == 7) then
       do i = minl+1,min(maxl,6)
          integ_mode(i) = 12
       end do
       do i = 7,maxl
          integ_mode(i) = -1
          if (i == maxl) integ_mode(i) = 1
       end do
    end if
    intcorner_deferred = all((integ_mode == -1) .or. (integ_mode == 11))

    ! interpret the ode-mode
    call map_ode_pointers(qtree_ode_mode)

    ! set properties to calculate and derivatives required
    savefgr = (prop_mode == 1 .or. prop_mode == 2)
    savelapgr = (prop_mode == 2)
    if (savelapgr) then
       nder = 2
    else
       nder = 1
    end if

    ! default for color allocation
    if (color_allocate < 0 .or. color_allocate > 1) then
       if (maxl > 8) then
          color_allocate = 0
       else
          color_allocate = 1
       end if
    end if
    if (plot_mode > 0 .and. color_allocate == 0) then
       write (uout,'("* COLOR_ALLOCATE has been set to 0.")')
       write (uout,'("* PLOT mode is set to 0.")') 
       call ferror('qtree','Non zero PLOT_MODE is not compatible with a single color array',warning)
       write (uout,*)
       plot_mode = 0
    end if
    if (docontacts .and. color_allocate == 0) then
       write (uout,'("* COLOR_ALLOCATE has been set to 0.")')
       write (uout,'("* The contacts are NOT going to be used.")') 
       call ferror('qtree','DOCONTACTS is not compatible with a single color array',warning)
       write (uout,*)
       docontacts = .false.
    end if

    ! Determine local point group
    if (ws_use) then
       ! Build initial tetrahedron list
       pg = sy%c%sitesymm(ws_origin,leqv=leqv,lrotm=lrotm)
       call sy%c%wigner(ws_origin,ntetrag=ntetrag,tetrag=tetrag) 
    else
       ws_scale = -1d0
       ws_origin = 0d0
       leqv = 1
       lrotm(:,:,1) = eye
       call sy%c%pmwigner(ntetrag=ntetrag,tetrag=tetrag)
    end if
    periodic = .true.
       
    ! Pre-split level
    if (verbose .and. plvl > 0) &
       write (uout,'("* Splitting the initial tetrahedra, ",I2," times."/)') plvl
    call presplit_ws(plvl,ntetrag,tetrag)

    ! Allocate
    nt_orig = ntetrag
    maxlen = -1d0
    minlen = 1d30
    allocate(torig(3,nt_orig),tvec(3,3,nt_orig),tvol(nt_orig))
    allocate(borig(3,nt_orig),bvec(3,3,nt_orig))
    allocate(cmat(3,3,nt_orig),dmat(3,3,nt_orig))

    ! Tetrahedra initialization
    do i = 1, nt_orig
       ! vertex in crystallographic coordinates
       torig(:,i) = tetrag(1,:,i) 
       tvec(:,1,i) = tetrag(2,:,i) - torig(:,i) 
       tvec(:,2,i) = tetrag(3,:,i) - torig(:,i) 
       tvec(:,3,i) = tetrag(4,:,i) - torig(:,i) 
       if (ws_scale > 0d0) then
          tvec(:,1:3,i) = tvec(:,1:3,i) / ws_scale
          periodic = .false.
       end if

       ! vertex in cartesian, maxlen, minlen, volume, handedness
       xp1 = sy%c%x2c(tvec(:,1,i))
       xp2 = sy%c%x2c(tvec(:,2,i))
       xp3 = sy%c%x2c(tvec(:,3,i))
       dist = dot_product(xp1,xp1)
       maxlen = max(maxlen,dist)
       minlen = min(minlen,dist)
       dist = dot_product(xp2,xp2)
       maxlen = max(maxlen,dist)
       minlen = min(minlen,dist)
       dist = dot_product(xp3,xp3)
       maxlen = max(maxlen,dist)
       minlen = min(minlen,dist)
       dist = dot_product(xp3-xp1,xp3-xp1)
       maxlen = max(maxlen,dist)
       minlen = min(minlen,dist)
       dist = dot_product(xp3-xp2,xp3-xp2)
       maxlen = max(maxlen,dist)
       minlen = min(minlen,dist)
       dist = dot_product(xp2-xp1,xp2-xp1)
       maxlen = max(maxlen,dist)
       minlen = min(minlen,dist)

       tvol(i) = mixed(xp1,xp2,xp3) / 6d0
       if (tvol(i) < 0d0) then
          xx = tvec(:,2,i)
          tvec(:,2,i) = tvec(:,3,i)
          tvec(:,3,i) = xx
          xx = xp2
          xp2 = xp3
          xp3 = xx
          tvol(i) = - tvol(i)
       end if

       ! toconvex matrix
       cmat(1,:,i) = cross(xp2,xp3) / (6d0 * tvol(i))
       cmat(2,:,i) = cross(xp3,xp1) / (6d0 * tvol(i))
       cmat(3,:,i) = cross(xp1,xp2) / (6d0 * tvol(i))

       dmat(:,:,i) = cmat(:,:,i)
       call dgeco(dmat(:,:,i),3,3,iw,er,r1)
       call dgedi(dmat(:,:,i),3,3,iw,xx,r1,1)

       ! fill cartesian / 2**maxl vectors
       borig(:,i) = sy%c%x2c(torig(:,i))
       bvec(:,1,i) = xp1 
       bvec(:,2,i) = xp2
       bvec(:,3,i) = xp3
       bvec(:,:,i) = bvec(:,:,i) / l2
    end do
    deallocate(tetrag)

    ! total volume and maximum and minimum lengths
    maxlen = sqrt(maxlen)
    minlen = sqrt(minlen)
    if (ws_scale > 0d0) then
       vtotal = sy%c%omega / ws_scale**3 
    else
       vtotal = sy%c%omega 
    end if
    
    if (verbose) then
       write (uout,'("* Cell list of tetrahedra (cryst. coord.) contains : ",I3)') ntetrag
       sumi = 0d0
       do i = 1, ntetrag
          write (uout,'("+ Tetrahedron : ",I3," with points at ")') i
          do j = 1, 4
             write (uout,'(4X,I3,3(1X,E20.13))') j, tetrag(j,1,i), tetrag(j,2,i), tetrag(j,3,i)
          end do
          write (uout,'(4X," Volume : ",1p,E20.12)') tvol(i)
          sumi = sumi + tvol(i)
       end do
       write (uout,'("+ Sum of tetrahedra volumes : ",1p,E20.12)') sumi
       write (uout,'("+ Cell volume               : ",1p,E20.12)') sy%c%omega / sy%c%ncv
       write (uout,*)
       write (uout,'("* END of cell construction.")')
       write (uout,*)
    end if

    ! eps for crys2convex
    crys2convex_eps = tetra_eps_warn * maxlen / 2**maxl
    crys2convex_eps1 = 1d0 + crys2convex_eps

    ! check symmetry
    if (ws_use) then
       call qtree_checksymmetry()
    end if

    ! tetrahedra contacts
    if (docontacts) then
       allocate(tcontact(nt_orig,4))
       tcontact = tcontact_void
       call find_tetrah_contacts()
    end if

    ! allocate atomic property arrays
    allocate(atprop(nnuc+3,sy%npropi))
    allocate(acum_atprop(nnuc+3,sy%npropi))
    atprop = 0d0
    acum_atprop = 0d0

    ! allocate big arrays
    siz = 2**maxl + 1
    siz = siz*(siz+1)*(siz+2)/6
    if (color_allocate == 1) then
       allocate(trm(siz,nt_orig))
    else
       allocate(trm(siz,1))
    end if
    if (savefgr) then
       if (color_allocate == 1) then
          allocate(fgr(siz,nt_orig))
       else
          allocate(fgr(siz,1))
       end if
    end if
    if (savelapgr) then
       if (color_allocate == 1) then
          allocate(lapgr(siz,nt_orig))
       else
          allocate(lapgr(siz,1))
       end if
    end if
    if (intcorner_deferred) then
       allocate(vgr(siz))
    end if
    trm = 0
    if (savefgr) fgr = -1d0
    if (savelapgr) lapgr = 0d0
    if (allocated(vgr)) vgr = 0d0

  end subroutine qtree_initialize

  !> Check that the symmetry of the cell and the tetrahedra are consistent.
  module subroutine qtree_checksymmetry()
    use systemmod, only: sy
    use global, only: ws_scale

    real*8 :: sumi
    integer :: i
    real*8 :: vtotal

    if (ws_scale > 0d0) then
       vtotal = sy%c%omega / ws_scale**3 
    else
       vtotal = sy%c%omega 
    end if

    sumi = 0d0
    do i = 1, nt_orig
       sumi = sumi + tvol(i)
    end do
    leqvf = vtotal / sumi

  end subroutine qtree_checksymmetry

  !> Clean up allocated memory.
  module subroutine qtree_cleanup()

    if (allocated(atprop)) deallocate(atprop)
    if (allocated(torig)) deallocate(torig)
    if (allocated(borig)) deallocate(borig)
    if (allocated(tvec)) deallocate(tvec) 
    if (allocated(tvol)) deallocate(tvol)
    if (allocated(tcontact)) deallocate(tcontact)
    if (allocated(bvec)) deallocate(bvec) 
    if (allocated(cmat)) deallocate(cmat)
    if (allocated(dmat)) deallocate(dmat)

  end subroutine qtree_cleanup

  !> Find the beta-spheres that fulfill Rodriguez et al. criterion:
  !> the gradient at a collection of points is at most 45 degrees from
  !> the radial direction.
  module subroutine find_beta_rodriguez(nuc,rbeta)
    use systemmod, only: sy
    use surface, only: minisurf
    use param, only: pi
    use types, only: scalar_value
    integer, intent(in) :: nuc
    real*8, intent(inout) :: rbeta

    real*8, parameter :: shrink = 0.9d0
    real*8, parameter :: angdev = 45d0

    type(minisurf) :: srf
    real*8 :: xnuc(3), unit(3), x(3), angle
    integer :: i
    logical :: accept
    type(scalar_value) :: res

    xnuc = sy%f(sy%iref)%cp(nuc)%x
    angle = cos(angdev * pi / 180d0)

    call srf%init(100,100)
    call srf%clean()
    call srf%spherecub(xnuc,2)
    
    accept = .false.
    do while(.not.accept)
       accept = .true.
       do i = 1, srf%nv
          unit = (/ sin(srf%th(i)) * cos(srf%ph(i)),&
             sin(srf%th(i)) * sin(srf%ph(i)),&
             cos(srf%th(i)) /)
          x = xnuc + unit * rbeta
          call sy%f(sy%iref)%grd(x,1,res)
          accept = accept .and. (dot_product(-res%gf/res%gfmod,unit) >= angle) 
          if (.not.accept) then
             rbeta = rbeta * shrink
             exit
          end if
       end do
    end do

    call srf%end()

  end subroutine find_beta_rodriguez

  !> Remap the pointers to the given one-step ODE solver.
  module subroutine map_ode_pointers(odem)
    use global, only: ode_abserr
    use tools_io, only: ferror, faterr

    integer, intent(in) :: odem

    ! pointers to ode solver
    if (odem == 1) then
       ode_o => euler_o
       ode_a => euler_a
       ode_b => euler_b
       ode_type => euler_type
       ode_fsal = .false.
    else if (odem == 2) then
       ode_o => heun_o
       ode_a => heun_a
       ode_b => heun_b
       ode_type => heun_type
       ode_fsal = .false.
    else if (odem == 3) then
       ode_o => kutta_o
       ode_a => kutta_a
       ode_b => kutta_b
       ode_type => kutta_type
       ode_fsal = .false.
    else if (odem == 4) then
       ode_o => rk4_o
       ode_a => rk4_a
       ode_b => rk4_b
       ode_type => rk4_type
       ode_fsal = .false.
    else if (odem == 5) then
       ode_o => heul_o
       ode_o2 => heul_o2
       ode_a => heul_a
       ode_b => heul_b
       ode_b2 => heul_b2
       ode_type => heul_type
       ode_fsal = .false.
    else if (odem == 6) then
       ode_o => bs_o
       ode_o2 => bs_o2
       ode_a => bs_a
       ode_b => bs_b
       ode_b2 => bs_b2
       ode_type => bs_type
       ode_fsal = .true.
    else if (odem == 7) then
       ode_o => rkck_o
       ode_o2 => rkck_o2
       ode_a => rkck_a
       ode_b => rkck_b
       ode_b2 => rkck_b2
       ode_type => rkck_type
       ode_fsal = .false.
    else if (odem == 8) then
       ode_o => dp_o
       ode_o2 => dp_o2
       ode_a => dp_a
       ode_b => dp_b
       ode_b2 => dp_b2
       ode_type => dp_type
       ode_fsal = .true.
    else
       call ferror('map_ode_pointers','unknown odem',faterr)
    end if
    if (ode_type == 1) then
       qinv = 1.d0 / real(ode_o2,8)
       q1inv = 1.d0 / real(ode_o2+1,8)
    end if
    ! apply defaults to ode_abserr
    if (ode_abserr < 0d0) then
       if (odem == 5) then
          ode_abserr = 1d-3
       else if (odem == 6) then
          ode_abserr = 1d-4
       else if (odem == 7) then
          ode_abserr = 1d-4
       else 
          ode_abserr = 1d-4
       end if
    end if

  end subroutine map_ode_pointers

  !> Given the scaled convex coordinates of a grid point, return its
  !> index in the 1D color or properties vectors.
  module function cindex(i,lvl)

    integer, intent(in) :: i(3)
    integer, intent(in) :: lvl
    integer(qtreeidx) :: cindex

    integer(qtreeidx) :: n1, n2, n3

    n1 = (i(1) + i(2) + i(3)) * 2**(maxl-lvl) - 1
    n2 = (i(2) + i(3)) * 2**(maxl-lvl) - 1
    n3 = i(3) * 2**(maxl-lvl) - 1
    
    cindex = 3 + n1 * (12+(n1+1)*(10+2*n1)) / 12 + (n2+1)*(n2+2)/2 + n3 
    ! (n1+1)*(n1+2)*(n1+3)/6 + (n2+1)*(n2+2)/2 + (n3+1) + 1

  end function cindex

  !> Transforms from crystallographic to convex coordinates. If the
  !> point is not inside any known IWST, base_t = 0 is returned.
  module subroutine crys2convex(x,base_t,rver)
    use systemmod, only: sy
    real*8, intent(in) :: x(3)
    integer, intent(out) :: base_t
    real*8, intent(out) :: rver(3)

    integer :: i
    real*8 :: xp(3)

    xp = matmul(sy%c%crys2car(1:3,1:3),x(1:3))
    base_t = 0
    do i = 1, nt_orig
       rver = matmul(cmat(:,:,i),xp - borig(:,i))
       !if (all(rver >= -crys2convex_eps) .and. (1d0-sum(rver) > -crys2convex_eps)) then
       if (all(rver >= -crys2convex_eps) .and. (sum(rver) < crys2convex_eps1)) then
          base_t = i
          return
       end if
    end do
    rver = 0d0

  end subroutine crys2convex

  !> Given a point in crystallographic coordinates, determine the IWST to
  !> which it belongs and the associated convex coordinates. Also,
  ! return the poi
  module subroutine locate_tetrah(x,base_t,rver,lrot)
    use global, only: ws_origin
    real*8, intent(inout) :: x(3)
    integer, intent(out) :: base_t
    real*8, intent(out) :: rver(3)
    integer, intent(out) :: lrot

    integer :: i, h1, h2 ,h3 
    real*8 :: xi(3), xic(3), xic0(3)
    real*8 :: eps

    nlocate = nlocate + 1

1   continue

    xic0 = x - nint(x)
    do h1 = -1, 1
       do h2 = -1, 1
          do h3 = -1, 1
             xic = xic0 + (/h1,h2,h3/)
             do i = 1, leqv
                xi = matmul(lrotm(:,1:3,i),xic-ws_origin) + ws_origin
                call crys2convex(xi,base_t,rver)
                if (base_t > 0) then
                   x = xic
                   lrot = i
                   return
                end if
             end do
          end do
       end do
    end do
    if (base_t == 0 .and. eps < tetra_eps_strict) then
       nlocate_sloppy = nlocate_sloppy + 1
       eps = tetra_eps_strict
       goto 1
    end if

  end subroutine locate_tetrah

  !> Given a point xp in cartesian coordinates, this routine
  !> detrmines the nearest grid point and returns its tetrahedron and
  !> index, the rotation matrix from the site-symmetry point group of
  !> the origin to use on xp and the distance
  module subroutine neargp(xp,base_t,lrot,idx,dist)
    use systemmod, only: sy
    use global, only: ws_origin
    use tools_io, only: ferror, faterr
    real*8, intent(inout) :: xp(3)
    integer, intent(inout) :: base_t
    integer, intent(inout) :: lrot
    integer(qtreeidx), intent(out) :: idx
    real*8, intent(out) :: dist

    real*8 :: xnew(3), rver(3), difver(3), ndist
    integer :: nbase, neigh(3,8), l2
    integer :: i, nn

    l2 = 2**maxl

    ! tasks in crystallographic
    xp = sy%c%c2x(xp)

    ! transform to WS and find base tetrahedron and rotation matrix.
    xnew = matmul(lrotm(:,1:3,lrot),xp-ws_origin) + ws_origin
    call crys2convex(xnew,nbase,rver)
    if (nbase /= base_t) then
       call locate_tetrah(xp,base_t,rver,lrot)
       if (base_t == 0) then
          if (periodic) then
             call ferror('neargp','locate_tetrah did not find tetrah (tetraeps?)',faterr)
          else
             return
          end if
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

    ! Calculate nearest neighbor and distance
    dist = 1d30
    nn = 0
    do i = 1, 8
       if (any(neigh(:,i) < 0) .or. sum(neigh(:,i)) > l2) cycle
       difver = (rver - real(neigh(:,i),8)) / l2
       xnew = borig(:,base_t) + matmul(dmat(:,:,base_t),difver)
       ndist = dot_product(xnew,xnew)
       if (ndist < dist) then
          nn = i
          dist = ndist
       end if
    end do
    rver = real(neigh(:,nn),8) / l2
    xp = borig(:,base_t) + matmul(dmat(:,:,base_t),rver)
    idx = cindex(neigh(:,nn),maxl)
    dist = sqrt(dist)

  end subroutine neargp

  !> Find the contacts between tetrahedron faces
  module subroutine find_tetrah_contacts()
    use systemmod, only: sy
    use tools_io, only: uout, ferror, faterr
    integer :: t, f, p, op, c, i, t2, f2, invp, invc, invop
    real*8 :: xface(3,3,4,nt_orig) ! xface(vcoords,vertex,face,tetrah)
    real*8 :: aux(3,3), aux2(3,3)
    logical :: cont
    integer :: ldiff(3)

    ! fill the faces
    do t = 1, nt_orig
       ! 0 1 2
       xface(:,1,1,t) = torig(:,t)
       xface(:,2,1,t) = torig(:,t) + tvec(:,1,t)
       xface(:,3,1,t) = torig(:,t) + tvec(:,2,t)

       ! 0 1 3
       xface(:,1,2,t) = torig(:,t)
       xface(:,2,2,t) = torig(:,t) + tvec(:,1,t)
       xface(:,3,2,t) = torig(:,t) + tvec(:,3,t)

       ! 0 2 3
       xface(:,1,3,t) = torig(:,t)
       xface(:,2,3,t) = torig(:,t) + tvec(:,2,t)
       xface(:,3,3,t) = torig(:,t) + tvec(:,3,t)

       ! 1 2 3
       xface(:,1,4,t) = torig(:,t) + tvec(:,1,t)
       xface(:,2,4,t) = torig(:,t) + tvec(:,2,t)
       xface(:,3,4,t) = torig(:,t) + tvec(:,3,t)
    end do

    ! run over faces and check for equivalence
    ! full point symmetry -> if 3 vertex are equivalent, the
    !  plane is equivalent so trms can be copied safely.
    do t = 1, nt_orig
       do f = 1, 4
          ! skip known faces
          if (tcontact(t,f) /= tcontact_void) cycle
          
          cont = .false.
          t1: do op = 1, leqv
             c = 1
             do p = 1, 6
                ! transformed face
                aux = matmul(lrotm(1:3,1:3,op),xface(:,:,f,t))
                do i = 1, 3
                   aux(:,i) = aux(:,i) + sy%c%cen(:,c)
                end do
                aux = matmul(aux,perm3(:,:,p))

                ! compare with the other faces in the same tetrah.
                do f2 = f, 4
                   aux2 = aux - xface(:,:,f2,t)
                   ldiff = nint(aux2(:,1))
                   ! careful with identity
                   if (f2 == f .and. op == 1 .and. c == 1) cycle
                   do i = 1, 3
                      aux2(:,i) = abs(aux2(:,i) - ldiff)
                   end do
                   if (all(aux2 < eps_tetrah_contact)) then
                      cont = .true.
                      tcontact(t,f) = (op-1)+leqv*((c-1)+1*((p-1)+6*((f2-1)+4*(t-1))))+1
                      call inverse_operation(p,op,c,invp,invop,invc)
                      tcontact(t,f2) = -((invop-1)+leqv*((invc-1)+1*((invp-1)+6*((f-1)+4*(t-1))))+1)
                      exit t1
                   end if
                end do
                ! compare with other tetrah.
                do t2 = t+1, nt_orig
                   do f2 = 1, 4
                      aux2 = aux - xface(:,:,f2,t2)
                      ldiff = nint(aux2(:,1))
                      do i = 1, 3
                         aux2(:,i) = abs(aux2(:,i) - ldiff)
                      end do
                      if (all(aux2 < eps_tetrah_contact)) then
                         cont = .true.
                         tcontact(t,f) = (op-1)+leqv*((c-1)+1*((p-1)+6*((f2-1)+4*(t2-1))))+1
                         call inverse_operation(p,op,c,invp,invop,invc)
                         tcontact(t2,f2) = -((invop-1)+leqv*((invc-1)+1*((invp-1)+6*((f-1)+4*(t-1))))+1)
                         exit t1
                      end if
                   end do
                end do
             end do
          end do t1
          if (.not.cont .and. periodic) then
             write (uout,'(" t = ",I3," f = ",I3)') t, f
             call ferror('find_tetrah_contacts','failed to determine a contact, use nocontacts',faterr)
          end if
       end do
    end do

  end subroutine find_tetrah_contacts

  !> Calculate the inverse of a tetrahedron face transformation, given
  !> by a rotation (op), a centering translation (c) and a permutation of the 
  !> vertex (p).
  module subroutine inverse_operation(p,op,c,invp,invop,invc)
    use systemmod, only: sy
    use tools_io, only: ferror, faterr
    use param, only: eye
    integer, intent(in) :: p, op, c
    integer, intent(out) :: invp, invop, invc

    logical :: found
    real*8 :: xaux(3)

    ! inverse permutation
    if (p == 2) then
       invp = 3
    else if (p == 3) then
       invp = 2
    else
       invp = p
    end if
    ! inverse rotation
    found = .false.
    do invop = 1, leqv
       if (all(abs(matmul(lrotm(1:3,1:3,op),lrotm(1:3,1:3,invop)) - eye) < eps_tetrah_contact)) then
          found = .true.
          exit
       end if
    end do
    if (.not.found) call ferror('inverse_operation','could not find matrix inverse operation',faterr)
    ! inverse translation
    found = .false.
    do invc = 1, sy%c%ncv
       xaux = sy%c%cen(:,c) - sy%c%cen(:,invc)
       xaux = xaux - nint(xaux)
       if (all(abs(xaux) < eps_tetrah_contact)) then
          found = .true.
          exit
       end if
    end do
    if (.not.found) call ferror('inverse_operation','could not find translation inverse',faterr)

  end subroutine inverse_operation

  !> Calculate the maximum and minimum lengths of the grid 
  !> at the highest subdivision level.
  module subroutine get_tlengths(minlen,maxlen)

    real*8, intent(out) :: minlen, maxlen

    integer :: sn
    integer :: s_iv(3,4,8*maxl)
    integer :: s_l(8*maxl)

    integer :: iv(3,4), l, tt, k, l2
    integer :: dif(3)
    real*8 :: d2, xx(3)

    l2 = 2**maxl

    sn = 1
    s_iv(:,1,1) = (/ 0, 0, 0 /)
    s_iv(:,2,1) = (/ 1, 0, 0 /)
    s_iv(:,3,1) = (/ 0, 1, 0 /)
    s_iv(:,4,1) = (/ 0, 0, 1 /)
    s_l(1) = 0
    minlen = 1d20
    maxlen = 0d0
    do while(sn > 0) 
       iv = s_iv(:,:,sn)
       l = s_l(sn)
       sn = sn - 1

       if (l < maxl) then
          ! 1, 1-2, 1-3, 1-4
          sn = sn + 1
          s_l(sn) = l+1
          s_iv(:,1,sn) = iv(:,1) * 2
          s_iv(:,2,sn) = iv(:,1) + iv(:,2)
          s_iv(:,3,sn) = iv(:,1) + iv(:,3)
          s_iv(:,4,sn) = iv(:,1) + iv(:,4)

          ! 2, 1-2, 2-3, 2-4
          sn = sn + 1
          s_l(sn) = l+1
          s_iv(:,1,sn) = iv(:,2) * 2
          s_iv(:,2,sn) = iv(:,1) + iv(:,2)
          s_iv(:,3,sn) = iv(:,2) + iv(:,3)
          s_iv(:,4,sn) = iv(:,2) + iv(:,4)

          ! 3, 1-3, 2-3, 3-4
          sn = sn + 1
          s_l(sn) = l+1
          s_iv(:,1,sn) = iv(:,3) * 2
          s_iv(:,2,sn) = iv(:,1) + iv(:,3)
          s_iv(:,3,sn) = iv(:,2) + iv(:,3)
          s_iv(:,4,sn) = iv(:,3) + iv(:,4)

          ! 4, 1-4, 2-4, 3-4
          sn = sn + 1
          s_l(sn) = l+1
          s_iv(:,1,sn) = iv(:,4) * 2
          s_iv(:,2,sn) = iv(:,1) + iv(:,4)
          s_iv(:,3,sn) = iv(:,2) + iv(:,4)
          s_iv(:,4,sn) = iv(:,3) + iv(:,4)

          ! 2-3, 1-2, 1-3, 1-4
          sn = sn + 1
          s_l(sn) = l+1
          s_iv(:,1,sn) = iv(:,2) + iv(:,3)
          s_iv(:,2,sn) = iv(:,1) + iv(:,2)
          s_iv(:,3,sn) = iv(:,1) + iv(:,3)
          s_iv(:,4,sn) = iv(:,1) + iv(:,4)

          ! 1-4, 1-2, 2-3, 2-4
          sn = sn + 1
          s_l(sn) = l+1
          s_iv(:,1,sn) = iv(:,1) + iv(:,4)
          s_iv(:,2,sn) = iv(:,1) + iv(:,2)
          s_iv(:,3,sn) = iv(:,2) + iv(:,3)
          s_iv(:,4,sn) = iv(:,2) + iv(:,4)

          ! 1-4, 1-3, 2-3, 3-4
          sn = sn + 1
          s_l(sn) = l+1
          s_iv(:,1,sn) = iv(:,1) + iv(:,4)
          s_iv(:,2,sn) = iv(:,1) + iv(:,3)
          s_iv(:,3,sn) = iv(:,2) + iv(:,3)
          s_iv(:,4,sn) = iv(:,3) + iv(:,4)

          ! 2-3, 1-4, 2-4, 3-4
          sn = sn + 1
          s_l(sn) = l+1
          s_iv(:,1,sn) = iv(:,2) + iv(:,3)
          s_iv(:,2,sn) = iv(:,1) + iv(:,4)
          s_iv(:,3,sn) = iv(:,2) + iv(:,4)
          s_iv(:,4,sn) = iv(:,3) + iv(:,4)
       else
          dif = iv(:,2) - iv(:,1)
          do tt = 1, nt_orig
             xx = 0d0
             do k = 1, 3
                xx = xx + dif(k) * bvec(:,k,tt)
             end do
             d2 = dot_product(xx,xx)
             minlen = min(d2,minlen)
             maxlen = max(d2,maxlen)
          end do
       end if
    end do
    minlen = sqrt(minlen)
    maxlen = sqrt(maxlen) 

  end subroutine get_tlengths

  !> Pre-split the initial tetrahedra list p times.
  module subroutine presplit_ws(plvl,ntetrag,tetrag)
    integer, intent(in) :: plvl
    integer, intent(inout) :: ntetrag
    real*8, intent(inout), allocatable :: tetrag(:,:,:)

    integer :: i, l, n
    real*8 :: tsave(0:3,3)

    do l = 1, plvl
       n = ntetrag
       do i = 1, n

          tsave = tetrag(:,1:3,i)

          ! 1, 1-2, 1-3, 1-4
          tetrag(1,:,i) = tsave(0,:) 
          tetrag(2,:,i) = tsave(0,:) + tsave(1,:)
          tetrag(3,:,i) = tsave(0,:) + tsave(2,:)
          tetrag(4,:,i) = tsave(0,:) + tsave(3,:)
          tetrag(2:4,:,i) = 0.5d0 * tetrag(2:4,:,i) 

          ! 2, 1-2, 2-3, 2-4
          ntetrag = ntetrag + 1
          tetrag(1,:,ntetrag) = tsave(1,:) 
          tetrag(2,:,ntetrag) = tsave(0,:) + tsave(1,:)
          tetrag(3,:,ntetrag) = tsave(1,:) + tsave(2,:)
          tetrag(4,:,ntetrag) = tsave(1,:) + tsave(3,:)
          tetrag(2:4,:,ntetrag) = 0.5d0 * tetrag(2:4,:,ntetrag) 

          ! 3, 1-3, 2-3, 3-4
          ntetrag = ntetrag + 1
          tetrag(1,:,ntetrag) = tsave(2,:) 
          tetrag(2,:,ntetrag) = tsave(0,:) + tsave(2,:)
          tetrag(3,:,ntetrag) = tsave(1,:) + tsave(2,:)
          tetrag(4,:,ntetrag) = tsave(2,:) + tsave(3,:)
          tetrag(2:4,:,ntetrag) = 0.5d0 * tetrag(2:4,:,ntetrag) 

          ! 4, 1-4, 2-4, 3-4
          ntetrag = ntetrag + 1
          tetrag(1,:,ntetrag) = tsave(3,:) 
          tetrag(2,:,ntetrag) = tsave(0,:) + tsave(3,:)
          tetrag(3,:,ntetrag) = tsave(1,:) + tsave(3,:)
          tetrag(4,:,ntetrag) = tsave(2,:) + tsave(3,:)
          tetrag(2:4,:,ntetrag) = 0.5d0 * tetrag(2:4,:,ntetrag) 

          ! 2-3, 1-2, 1-3, 1-4
          ntetrag = ntetrag + 1
          tetrag(1,:,ntetrag) = tsave(1,:) + tsave(2,:) 
          tetrag(2,:,ntetrag) = tsave(0,:) + tsave(1,:)
          tetrag(3,:,ntetrag) = tsave(0,:) + tsave(2,:)
          tetrag(4,:,ntetrag) = tsave(0,:) + tsave(3,:)
          tetrag(:,:,ntetrag) = 0.5d0 * tetrag(:,:,ntetrag) 

          ! 1-4, 1-2, 2-3, 2-4
          ntetrag = ntetrag + 1
          tetrag(1,:,ntetrag) = tsave(0,:) + tsave(3,:) 
          tetrag(2,:,ntetrag) = tsave(0,:) + tsave(1,:)
          tetrag(3,:,ntetrag) = tsave(1,:) + tsave(2,:)
          tetrag(4,:,ntetrag) = tsave(1,:) + tsave(3,:)
          tetrag(:,:,ntetrag) = 0.5d0 * tetrag(:,:,ntetrag) 

          ! 1-4, 1-3, 2-3, 3-4
          ntetrag = ntetrag + 1
          tetrag(1,:,ntetrag) = tsave(0,:) + tsave(3,:) 
          tetrag(2,:,ntetrag) = tsave(0,:) + tsave(2,:)
          tetrag(3,:,ntetrag) = tsave(1,:) + tsave(2,:)
          tetrag(4,:,ntetrag) = tsave(2,:) + tsave(3,:)
          tetrag(:,:,ntetrag) = 0.5d0 * tetrag(:,:,ntetrag) 

          ! 2-3, 1-4, 2-4, 3-4
          ntetrag = ntetrag + 1
          tetrag(1,:,ntetrag) = tsave(1,:) + tsave(2,:) 
          tetrag(2,:,ntetrag) = tsave(0,:) + tsave(3,:)
          tetrag(3,:,ntetrag) = tsave(1,:) + tsave(3,:)
          tetrag(4,:,ntetrag) = tsave(2,:) + tsave(3,:)
          tetrag(:,:,ntetrag) = 0.5d0 * tetrag(:,:,ntetrag) 
       end do
    end do

  end subroutine presplit_ws

end submodule proc
