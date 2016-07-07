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

!> Qtree, work with tetrahedra.
module qtree_tetrawork
  implicit none

  private
  public :: tetrah_subdivide
  public :: term_rec
  public :: tetrah_paint
  public :: integ_inner_keast
  public :: integ_inner_cubpack
  public :: integ_border_keast
  public :: integ_border_cubpack
  public :: integ_corner
  public :: integ_corner_deferred
  public :: integ_corner_sum
  public :: paint_inside_spheres
  private :: cubpack_f

contains

  !> Stack-based recursive subdivision of one IWST with in-line
  !> integration.
  subroutine tetrah_subdivide(base_t,iiv,il,acum_atprop,trm,fgr,lapgr,vgr)
    use qtree_basic
    use global
    use struct_basic
    use tools_io
    
    integer, intent(in) :: base_t
    integer, intent(in) :: iiv(3,4)
    integer, intent(in) :: il
    integer(qtreei), intent(inout) :: trm(:,:)
    real(qtreer), intent(inout) :: fgr(:,:), lapgr(:,:), vgr(:)
    real*8, intent(inout) :: acum_atprop(:,:)
    
    integer :: i, j, unk, ats(4)
    integer :: ts(4)
    integer(qtreeidx) :: idxx(4)
    integer :: base_to

    integer :: sn
    integer :: s_iv(3,4,8*maxl)
    integer :: s_l(8*maxl)

    integer :: iv(3,4)
    integer :: l, lrest, imin, imax
    real*8 :: xp(3,4)
    
    sn = 1
    s_iv(:,:,1) = iiv
    s_l(1) = il
    if (color_allocate == 0) then
       base_to = 1
    else
       base_to = base_t
    end if

    do while (sn > 0)
       iv = s_iv(:,:,sn)
       l = s_l(sn)
       lrest = 2**(maxl-l)
       sn = sn - 1

       if (plot_mode > 0 .and. plotsticks) then
          do i = 1, 4
             xp(:,i) = borig(:,base_t)
             do j = 1, 3
                xp(:,i) = xp(:,i) + bvec(:,j,base_t) * iv(j,i) * lrest
             end do
             xp(:,i) = cr%c2x(xp)
          end do
          do i = 1, 4
             do j = i+1, 4
                write (lustick(l),'(2(3(F8.5,X),X))') xp(:,i), xp(:,j)
             end do
          end do
       end if

       do i = 1, 4
          ts(i) = term_rec(base_t,iv(:,i),l,trm,fgr,lapgr)
          idxx(i) = cindex(iv(:,i),l)
          trm(idxx(i),base_to) = int(ts(i),1)
       end  do
       
       ats = abs(ts)
       if (all(ats(2:4) == ats(1)) .and. l > minl) then
          if (any(ts > 0) .and. l < maxl) then
             call tetrah_paint(base_t,iv,l,ats(1),trm)
          end if
          if (all(ts > 0)) then
             if (integ_mode(l) == 0 .or. ts(1) > nnuc) then
                ! only sum volume
                acum_atprop(ats(1),1) = acum_atprop(ats(1),1) + tvol(base_t) / 8**l
                cycle
             else if (integ_mode(l) >= 1 .and. integ_mode(l) <= 10) then
                ! keast 
                call integ_inner_keast(base_t,iv,l,ts(1),integ_mode(l),acum_atprop)
                cycle
             else if (integ_mode(l) == 11) then
                ! corner
                if (intcorner_deferred) then
                   call integ_corner_deferred(base_t,iv,l,ts,vgr)
                else
                   call integ_corner(base_t,iv,l,ts,acum_atprop,fgr,lapgr)
                end if
                cycle
             else if (integ_mode(l) == 12) then
                ! cubpack
                call integ_inner_cubpack(base_t,iv,l,ts(1),acum_atprop)
                cycle
             end if
          else if (all(ts < 0) .or. ts(1) > nnuc) then
             ! only sum volume
             acum_atprop(ats(1),1) = acum_atprop(ats(1),1) + tvol(base_t) / 8**l
             cycle
          end if
       end if

       if (l >= maxl) then
          ! check that no beta-sphere makes contact with an undecided tetrahedron
          if (any(ts < 0) .and. periodic .and. checkbeta) then
             imin = 99999
             imax = 0
             do i = 1, 4
                if (ts(i) == 0) cycle
                imin = min(abs(ts(i)),imin)
                imax = max(abs(ts(i)),imax)
             end do
             if (imin /= imax) then
                !$omp critical (IO)
                write (uout,*) " An undecided tetrahedron is overlapping with a beta-sphere. "
                write (uout,*) " Make beta-spheres smaller for atoms with negative index: "
                write (uout,*) " terms: ", ts
                write (uout,*) " Offending vertex: "
                write (uout,'("Tetrah: ",I3/,4(10X,"(",3(I5,X),")"/))') base_t, iv * lrest
                !$omp end critical (IO)
                call ferror('qtree_integration','beta-sphere leaks out of the basin',faterr)
             end if
          end if

          if (integ_mode(l) == 0) then
             ! only sum volume, factor out errors and in-cp grid points
             unk = count((ts == nnuc+1) .or. (ts == nnuc+3))
             do i = 1, 4
                if (ts(i) == 0) call ferror('tetrah_subdivide','zero term in border tetrahedron',faterr)
                if (ts(i) == nnuc+1 .or. ts(i) == nnuc+3) cycle
                acum_atprop(abs(ts(1)),1) = acum_atprop(abs(ts(1)),1) + tvol(base_t) / 8**l / (4d0-unk)
             end do
          else if (integ_mode(l) >= 1 .and. integ_mode(l) <= 10) then
             ! keast fixed rule
             call integ_border_keast(base_t,iv,l,ts,integ_mode(l),acum_atprop)
          else if (integ_mode(l) == 11) then
             ! corner
             if (intcorner_deferred) then
                call integ_corner_deferred(base_t,iv,l,ts,vgr)
             else
                call integ_corner(base_t,iv,l,ts,acum_atprop,fgr,lapgr)
             end if
          else if (integ_mode(l) == 12) then
             call integ_border_cubpack(base_t,iv,l,ts,acum_atprop)
          else
             call ferror("tetrah_subdivide","unknown integ_mode at maxl",faterr)
          end if
          cycle
       end if

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
    end do

  end subroutine tetrah_subdivide

  !> Determines the color of a given grid point, if it is not known. 
  function term_rec(base_t,iver,l,trm,fgr,lapgr)
    use qtree_basic
    use qtree_gpaths
    use varbas
    use fields
    use global
    use struct_basic
    use tools_io
    use types

    integer, intent(in) :: base_t, iver(3), l
    integer(qtreei), intent(inout) :: trm(:,:)
    real(qtreer), intent(inout) :: fgr(:,:), lapgr(:,:)
    integer :: term_rec

    integer :: j, ier
    integer :: nbase, termi(2), nini
    integer(qtreeidx) :: idx
    real*8 :: xp(3), lden, lrest, rver(3)
    real*8 :: xaux(3), xcrys(3)
    real*8 :: faux, gaux(3), gmodaux, raux(3)
    integer :: nid, ode_save, base_to
    real*8 :: dist
    type(scalar_value) :: res

    if (color_allocate == 0) then
       base_to = 1
    else
       base_to = base_t
    end if

    idx = cindex(iver,l)
    if (trm(idx,base_to) /= 0) then
       term_rec = trm(idx,base_to)
       return
    end if

    ! calculate point in cryst. coords
    lden = 2**l
    lrest = 2**(maxl-l)
    rver = real(iver,8) / lden
    xp = torig(:,base_t)
    do j = 1, 3
       xp = xp + tvec(:,j,base_t) * rver(j)
    end do
    xcrys = xp

    ! inside a beta-sphere?
    call nearest_cp(xp,nid,dist,type=f(refden)%typnuc)
    if (dist <= r_betagp(cpcel(nid)%idx)) then
       if (dist <= r_betaint(cpcel(nid)%idx)) then
          term_rec = -cpcel(nid)%idx
       else
          term_rec = cpcel(nid)%idx
       end if
       return
    end if

    ! value at point
    xp =cr%x2c(xp)
    call grd(f(refden),xp,nder,res)
    ngrd_term = ngrd_term + 1
    if (savefgr) fgr(idx,base_to) = res%fval
    if (savelapgr) lapgr(idx,base_to) = -res%del2fval

    nterm = nterm + 1
    nbase = base_t
    if (gradient_mode == 1) then
       call gradient_full(xp,nbase,rver,res,term_rec,ier)
    else if (gradient_mode == 2) then
       call gradient_color(xp,nbase,rver,res,term_rec,ier,trm)
    else if (gradient_mode == 3) then
       call gradient_qtree(xp,nbase,res,idx,term_rec,ier,.true.,trm,fgr,lapgr)
    else if (gradient_mode == -1) then
       xaux = xp
       faux = res%fval
       gaux = res%gf
       gmodaux = res%gfmod
       ode_save = qtree_ode_mode
       ! full with selected qtree_ode_mode
       gradient_mode = 1
       qtree_ode_mode = ode_save
       call map_ode_pointers(qtree_ode_mode)
       nini = ngrd_term
       xp = xaux
       res%fval = faux
       res%gf = gaux
       res%gfmod = gmodaux
       raux = rver
       nbase = base_t
       call gradient_full(xp,nbase,raux,res,termi(1),ier)
       ngrd1 = ngrd1 + ngrd_term - nini
       ! dormand-prince 4-5, embedded
       gradient_mode = 1
       qtree_ode_mode = 8
       call map_ode_pointers(qtree_ode_mode)
       nini = ngrd_term
       xp = xaux
       res%fval = faux
       res%gf = gaux
       res%gfmod = gmodaux
       raux = rver
       nbase = base_t
       call gradient_full(xp,nbase,raux,res,termi(2),ier)
       ngrd2 = ngrd2 + ngrd_term - nini

       write (uout,'(I3,X,3(F10.6,X),3(F12.8,X),I4,X,I4)') &
          base_t, rver, xcrys, termi(1), termi(2)
       if (termi(1) /= termi(2)) then
          ndiff = ndiff + 1
          write (ludif,'(2X,"ball ",3(F10.6,X),"type ",I2)') &
             xcrys(1), xcrys(2), xcrys(3), 1
       end if
       gradient_mode = -1
       qtree_ode_mode = ode_save
       call map_ode_pointers(qtree_ode_mode)
       term_rec = termi(2)
    else if (gradient_mode == -2) then
       xaux = xp
       faux = res%fval
       gaux = res%gf
       gmodaux = res%gfmod
       ode_save = qtree_ode_mode
       ! color with selected qtree_ode_mode
       gradient_mode = 2
       qtree_ode_mode = ode_save
       call map_ode_pointers(qtree_ode_mode)
       nini = ngrd_term
       xp = xaux
       res%fval = faux
       res%gf = gaux
       res%gfmod = gmodaux
       raux = rver
       nbase = base_t
       call gradient_color(xp,nbase,raux,res,termi(1),ier,trm)
       ngrd1 = ngrd1 + ngrd_term - nini
       ! dormand-prince 4-5, embedded
       gradient_mode = 1
       qtree_ode_mode = 8
       call map_ode_pointers(qtree_ode_mode)
       nini = ngrd_term
       xp = xaux
       res%fval = faux
       res%gf = gaux
       res%gfmod = gmodaux
       raux = rver
       nbase = base_t
       call gradient_full(xp,nbase,raux,res,termi(2),ier)
       ngrd2 = ngrd2 + ngrd_term - nini

       write (uout,'(I3,X,3(F10.6,X),3(F12.8,X),I4,X,I4)') &
          base_t, rver, xcrys, termi(1), termi(2)
       if (termi(1) /= termi(2)) then
          ndiff = ndiff + 1
          write (ludif,'(2X,"ball ",3(F10.6,X),"type ",I2)') &
             xcrys(1), xcrys(2), xcrys(3), 1
       end if
       gradient_mode = -2
       qtree_ode_mode = ode_save
       call map_ode_pointers(qtree_ode_mode)
       term_rec = termi(2)
    else if (gradient_mode == -3) then
       xaux = xp
       faux = res%fval
       gaux = res%gf
       gmodaux = res%gfmod
       ode_save = qtree_ode_mode
       ! qtree with selected qtree_ode_mode
       gradient_mode = 3
       qtree_ode_mode = ode_save
       call map_ode_pointers(qtree_ode_mode)
       nini = ngrd_term
       xp = xaux
       res%fval = faux
       res%gf = gaux
       res%gfmod = gmodaux
       raux = rver
       nbase = base_t
       call gradient_qtree(xp,nbase,res,idx,termi(1),ier,.false.,trm,fgr,lapgr)
       ngrd1 = ngrd1 + ngrd_term - nini
       ! dormand-prince 4-5, embedded
       gradient_mode = 1
       qtree_ode_mode = 8
       call map_ode_pointers(qtree_ode_mode)
       nini = ngrd_term
       xp = xaux
       res%fval = faux
       res%gf = gaux
       res%gfmod = gmodaux
       raux = rver
       nbase = base_t
       call gradient_full(xp,nbase,raux,res,termi(2),ier)
       ngrd2 = ngrd2 + ngrd_term - nini

       write (uout,'(I3,X,3(F10.6,X),3(F12.8,X),I4,X,I4)') &
          base_t, rver, xcrys, termi(1), termi(2)
       if (termi(1) /= termi(2)) then
          ndiff = ndiff + 1
          write (ludif,'(2X,"ball ",3(F10.6,X),"type ",I2)') &
             xcrys(1), xcrys(2), xcrys(3), 1
       end if
       gradient_mode = -3
       qtree_ode_mode = ode_save
       call map_ode_pointers(qtree_ode_mode)
       term_rec = termi(2)
    else
       call ferror("term_rec","wrong gradient_mode",faterr)
    end if

  end function term_rec

  !> Paint the interior of a tetrahedral region
  subroutine tetrah_paint(base_t,iv,l,color,trm)
    use qtree_basic
    use global

    integer, intent(in) :: base_t
    integer, intent(in) :: iv(3,4)
    integer, intent(in) :: l
    integer, intent(in) :: color
    integer(qtreei), intent(inout) :: trm(:,:)

    real*8, parameter :: eps = 1d-10

    integer :: i, j, k, o, n
    integer :: lrest, l2, vin(3)
    real*8 :: p(3,3), cond(3), acond(3)
    real*8 :: er, r1(3), xx(3), iw(3)
    real*8 :: c1, c2, rmin(4), rmax(4), q(3,4), dot
    integer :: count
    logical :: ok
    integer :: hmin, hmax, kmin, kmax, lmin, lmax, base_to
    integer(qtreeidx) :: idx

    if (color_allocate == 0) then
       base_to = 1
    else
       base_to = base_t
    end if
    lrest = 2**(maxl-l)
    l2 = 2**maxl
    
    hmin = minval(iv(1,:)) * lrest
    hmax = maxval(iv(1,:)) * lrest
    kmin = minval(iv(2,:)) * lrest
    kmax = maxval(iv(2,:)) * lrest
    lmin = minval(iv(3,:)) * lrest
    lmax = maxval(iv(3,:)) * lrest
    ! determine the conditions imposed by the faces on the interior
    ! points
    count = 0
    do i = 1, 2
       do j = i+1, 3
          do k = j+1, 4
             if (i == 2) then
                o = 1
             else if (j == 3) then
                o = 2
             else if (k == 4) then
                o = 3
             else
                o = 4
             end if
             count = count + 1
             if (iv(1,i) == 0 .and. iv(1,j) == 0 .and. iv(1,k) == 0) then
                cond = (/ 1d0, 0d0, 0d0 /)
             else if (iv(2,i) == 0 .and. iv(2,j) == 0 .and. iv(2,k) == 0) then
                cond = (/ 0d0, 1d0, 0d0 /)
             else if (iv(3,i) == 0 .and. iv(3,j) == 0 .and. iv(3,k) == 0) then
                cond = (/ 0d0, 0d0, 1d0 /)
             else
                p(1,:) = iv(:,i)
                p(2,:) = iv(:,j)
                p(3,:) = iv(:,k)
                call dgeco(p,3,3,iw,er,r1)
                call dgedi(p,3,3,iw,xx,r1,1)
                cond = (/1d0, 1d0, 1d0/)
                cond = matmul(p,cond)
                acond = abs(cond)
                where (acond<1d-12) acond = 1d30
                cond = cond / minval(acond)
             end if
             c1 = cond(1)*iv(1,i) + cond(2)*iv(2,i) + cond(3)*iv(3,i)
             c2 = cond(1)*iv(1,o) + cond(2)*iv(2,o) + cond(3)*iv(3,o)
             rmin(count) = min(c1,c2)
             rmax(count) = max(c1,c2)
             q(:,count) = cond
          end do
       end do
    end do
    rmin = rmin * lrest - eps
    rmax = rmax * lrest + eps

    do i = hmin, hmax
       do j = kmin, kmax
          do k = lmin, lmax
             if (i+j+k > l2) cycle
             ok = .true.
             do n = 1, 4
                dot = q(1,n)*i+q(2,n)*j+q(3,n)*k
                if (rmin(n) > dot .or. rmax(n) < dot) then
                   ok = .false.
                   exit
                end if
             end do
             if (ok) then
                vin = (/i,j,k/)
                idx = cindex(vin,maxl)
                ! do not overwrite in-sphere points
                if (trm(idx,base_to) >= 0) then
                   trm(idx,base_to) = int(color,1)
                end if
             end if
          end do
       end do
    end do

  end subroutine tetrah_paint

  !> Keast integration for a tetrahedron completely contained inside
  !> the region represented by the grid.
  subroutine integ_inner_keast(base_t,iv,l,color,klvl,acum_atprop)
    use qtree_basic
    use keast
    use fields
    use global

    integer, intent(in) :: base_t
    integer, intent(in) :: iv(3,4)
    integer, intent(in) :: l
    integer, intent(in) :: color
    integer, intent(in) :: klvl
    real*8, intent(inout) :: acum_atprop(:,:)

    integer :: i, j, lrest, l8
    real*8 :: xp(3)
    real*8 :: lprop(Nprops)
    real*8 :: ccrd(4)

    ! volume
    l8 = 8**l
    acum_atprop(abs(color),1) = acum_atprop(abs(color),1) + tvol(base_t) / l8
    if (color < 0) then
       return
    end if

    lrest = 2**(maxl-l)
    do i = 1, korder(klvl)
       ccrd(1:3) = kxyz(:,i,klvl)
       ccrd(4) = 1 - sum(ccrd(1:3))
       xp = borig(:,base_t)
       do j = 1, 3
          xp = xp + bvec(:,j,base_t) * dot_product(ccrd,iv(j,:)*lrest)
       end do
       call grdall(xp,lprop)
       ngrd_int = ngrd_int + 1
       acum_atprop(color,2:Nprops) = acum_atprop(color,2:Nprops) + tvol(base_t) / l8 * lprop(2:Nprops) * kw(i,klvl)
    end do

  end subroutine integ_inner_keast

  !> Cubpack integration for a tetrahedron completely contained inside
  !> the region represented by the grid.
  subroutine integ_inner_cubpack(base_t,iv,l,color,acum_atprop)
    use qtree_basic
    use CUI
    use varbas
    use fields
    use global
    use tools_io

    integer, intent(in) :: base_t
    integer, intent(in) :: iv(3,4)
    integer, intent(in) :: l
    integer, intent(in) :: color
    real*8, intent(inout) :: acum_atprop(:,:)

    integer :: i, j, lrest, l8
    real*8 :: vert(3,4,1)
    real*8 :: lprop(2:Nprops), cub_abserr(2:Nprops)
    integer :: ier, neval

    ! volume
    l8 = 8**l
    acum_atprop(abs(color),1) = acum_atprop(abs(color),1) + tvol(base_t) / 8**l
    if (color < 0) then
       return
    end if

    lrest = 2**(maxl-l)
    do i = 1, 4
       vert(:,i,1) = borig(:,base_t)
       do j = 1, 3
          vert(:,i,1) = vert(:,i,1) + bvec(:,j,base_t) * iv(j,i) * lrest
       end do
    end do

    ier = 1
    lprop = 0d0
    cub_abserr = 0d0
    call cubatr(3,Nprops-1,cubpack_f,1,vert,(/1/),lprop,cub_abserr,&
       IFAIL=ier, EpsAbs=cub_abs, EpsRel=cub_rel, MaxPts=cub_mpts, Neval=neval)
    ngrd_int = ngrd_int + neval
    if (ier == 1) then
       !$omp critical (IO)
       write (uout,*) "too many evaluations in cubpack"
       !$omp end critical (IO)
    else if (ier > 1) then
       call ferror('integ_inner_cubpack','severe error in cubpack',faterr)
    end if
    acum_atprop(color,2:Nprops) = acum_atprop(color,2:Nprops) + lprop

  end subroutine integ_inner_cubpack

  !> Keast integration for a tetrahedron on a (external or
  !> internal) interatomic surface
  subroutine integ_border_keast(base_t,iv,l,ts,klvl,acum_atprop)
    use qtree_basic
    use varbas
    use fields
    use global
    use keast
    use struct_basic
    use tools_math, only: mixed
    use tools_io

    integer, intent(in) :: base_t
    integer, intent(in) :: iv(3,4)
    integer, intent(in) :: l
    integer, intent(in) :: ts(4)
    integer, intent(in) :: klvl
    real*8, intent(inout) :: acum_atprop(:,:)

    integer :: i, j, k, l8, unk, its
    real*8 :: xp(3), temp(3)
    real*8 :: lprop(Nprops), lprop1(Nprops)
    integer :: nid
    real*8 :: ccrd(4), dist, xnuc(3), r2
    real*8 :: xvec(3,4), xint(3,4), xdot(4,4), dcoef, den, lvol
    real*8 :: vvec(3,3)
    integer :: ins, in1, out1, out2, lvec(3)
    logical :: in(4)

    l8 = 8**l
    nid = 0
    in1 = 0
    out1 = 0
    out2 = 0
    if (any(ts < 0)) then
       ! a sphere-interior boundary. all the abs(ts) are the same 
       its = abs(ts(1))
       acum_atprop(its,1) = acum_atprop(its,1) + tvol(base_t) / l8

       ! find the nucleus, use the barycenter.
       xp = borig(:,base_t)
       do i = 1, 3
          xp = xp + bvec(:,i,base_t) * 0.25d0 * sum(iv(i,:))
       end do
       xp = cr%c2x(xp)
       r2 = 1d30
       do i=1,ncpcel
          if (cpcel(i)%idx /= its) cycle
          temp = cpcel(i)%x - xp
          call cr%shortest(temp,dist)
          if (dist < r2) then
             nid = i
             r2 = dist
             lvec = nint(cpcel(nid)%x - xp - temp)
          end if
       end do
       xnuc = cr%x2c(cpcel(nid)%x - lvec)
       if (abs(sqrt(r2)-r_betaint(its)) > 2d0*maxlen) then
          call ferror('integ_border_keast','unknown xnuc for the tetrahedron',faterr)
       end if
       r2 = r_betaint(its) * r_betaint(its)

       ! build the tetrahedron vertex and inner product matrix
       ins = 0
       in = .false.
       do i = 1, 4
          xvec(:,i) = borig(:,base_t)
          do j = 1, 3
             xvec(:,i) = xvec(:,i) + bvec(:,j,base_t) * iv(j,i)
          end do
          do j = 1, i
             xdot(i,j) = dot_product(xvec(:,i)-xnuc,xvec(:,j)-xnuc)
             xdot(j,i) = xdot(i,j)
          end do
          if (xdot(i,i) < r2) then 
             ins = ins + 1
             in(i) = .true.
          end if
       end do

       if (ins == 1) then
          ! one vertex in the sphere -> 2 integrations.

          ! determine the index of the interior vertex
          do i = 1, 4
             if (in(i)) then
                in1 = i
                exit
             end if
          end do

          ! calculate edge intersections with the sphere
          ! and the small tetrahedron
          k = 0
          do i = 1, 4
             if (in1 == i) cycle
             k = k + 1
             den = xdot(i,i) + xdot(in1,in1) - 2d0*xdot(i,in1)
             dcoef = (-xdot(i,in1)+xdot(in1,in1) + &
                sqrt(xdot(i,in1)**2 - xdot(i,i)*xdot(in1,in1) + r2*den)) / den
             xvec(:,i) = (1-dcoef)*xvec(:,in1) + dcoef * xvec(:,i)
             vvec(:,k) = xvec(:,i) - xvec(:,in1)
             xdot(i,i) = dot_product(xvec(:,i)-xnuc,xvec(:,i)-xnuc)
          end do
          lvol = abs(mixed(vvec(:,1),vvec(:,2),vvec(:,3)) / 6d0)

          ! integrate 
          do i = 1, korder(klvl)
             ccrd(1:3) = kxyz(:,i,klvl)
             ccrd(4) = 1 - sum(ccrd(1:3))
             lprop = 0d0

             ! add the complete tetrahedron
             xp = borig(:,base_t)
             do j = 1, 3
                xp = xp + bvec(:,j,base_t) * dot_product(ccrd,iv(j,:))
             end do
             call grdall(xp,lprop1)
             ngrd_int = ngrd_int + 1
             lprop(2:Nprops) = lprop1(2:Nprops) * tvol(base_t) / l8

             ! substract the tetrahedron inside the sphere
             xp = 0d0
             do j = 1, 4
                xp = xp + xvec(:,j) * ccrd(j)
             end do
             call grdall(xp,lprop1)
             ngrd_int = ngrd_int + 1
             lprop1(2:Nprops) = lprop1(2:Nprops) * lvol
             lprop(2:Nprops) = lprop(2:Nprops) - lprop1(2:Nprops)

             acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) + lprop(2:Nprops) * kw(i,klvl)
          end do

       else if (ins == 3) then
          ! one vertex out of the sphere -> 1 integration.

          ! determine the index of the exterior vertex
          do i = 1, 4
             if (.not.in(i)) then
                out1 = i
                exit
             end if
          end do

          ! calculate edge intersections with the sphere
          ! and the small tetrahedron
          k = 0
          do i = 1, 4
             if (out1 == i) cycle
             k = k + 1
             den = xdot(i,i) + xdot(out1,out1) - 2d0*xdot(i,out1)
             dcoef = (-xdot(i,out1)+xdot(out1,out1) - &
                sqrt(xdot(i,out1)**2 - xdot(i,i)*xdot(out1,out1) + r2*den)) / den
             xvec(:,i) = (1-dcoef)*xvec(:,out1) + dcoef * xvec(:,i)
             vvec(:,k) = xvec(:,i) - xvec(:,out1)
          end do
          lvol = abs(mixed(vvec(:,1),vvec(:,2),vvec(:,3)) / 6d0)

          ! integrate 
          do i = 1, korder(klvl)
             ccrd(1:3) = kxyz(:,i,klvl)
             ccrd(4) = 1 - sum(ccrd(1:3))
             xp = 0d0
             do j = 1, 4
                xp = xp + xvec(:,j) * ccrd(j)
             end do
             call grdall(xp,lprop)
             ngrd_int = ngrd_int + 1
             acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) + lprop(2:Nprops) * lvol * kw(i,klvl)
          end do

       else if (ins == 2) then
          ! two vertex out of the sphere -> 1 integration.

          ! determine the index of the exterior vertex
          out1 = 0
          do i = 1, 4
             if (.not.in(i)) then
                if (out1 == 0) then
                   out1 = i
                else
                   out2 = i
                   exit
                end if
             end if
          end do

          ! calculate edge intersections of the vertex out1
          ! xint contains the four points, 1 and 2 associated to out1, 3 and 4 to out2
          ! 2 and 3 are opposite vertex (the interior vertex are different)
          k = 0
          do i = 1, 4
             if (out1 == i .or. out2 == i) cycle
             k = k + 1
             den = xdot(i,i) + xdot(out1,out1) - 2d0*xdot(i,out1)
             dcoef = (-xdot(i,out1)+xdot(out1,out1) - &
                sqrt(xdot(i,out1)**2 - xdot(i,i)*xdot(out1,out1) + r2*den)) / den
             xint(:,k) = (1-dcoef)*xvec(:,out1) + dcoef*xvec(:,i)
          end do
          do i = 1, 4
             if (out1 == i .or. out2 == i) cycle
             k = k + 1
             den = xdot(i,i) + xdot(out2,out2) - 2d0*xdot(i,out2)
             dcoef = (-xdot(i,out2)+xdot(out2,out2) - &
                sqrt(xdot(i,out2)**2 - xdot(i,i)*xdot(out2,out2) + r2*den)) / den
             xint(:,k) = (1-dcoef)*xvec(:,out2) + dcoef*xvec(:,i)
          end do

          ! integrate 
          do i = 1, korder(klvl)
             ccrd(1:3) = kxyz(:,i,klvl)
             ccrd(4) = 1 - sum(ccrd(1:3))
             lprop = 0d0

             ! add out1-1-2-3
             xp = 0d0
             xp = xp + xvec(:,out1) * ccrd(1)
             xp = xp + xint(:,1) * ccrd(2)
             xp = xp + xint(:,2) * ccrd(3)
             xp = xp + xint(:,3) * ccrd(4)
             lvol = abs(mixed(xint(:,1)-xvec(:,out1),&
                xint(:,2)-xvec(:,out1),&
                xint(:,3)-xvec(:,out1)) / 6d0)
             call grdall(xp,lprop1)
             ngrd_int = ngrd_int + 1
             lprop1(2:Nprops) = lprop1(2:Nprops) * lvol
             lprop(2:Nprops) = lprop(2:Nprops) + lprop1(2:Nprops)
             ! add out2-3-4-2
             xp = 0d0
             xp = xp + xvec(:,out2) * ccrd(1)
             xp = xp + xint(:,3) * ccrd(2)
             xp = xp + xint(:,4) * ccrd(3)
             xp = xp + xint(:,2) * ccrd(4)
             lvol = abs(mixed(xint(:,3)-xvec(:,out2),&
                xint(:,4)-xvec(:,out2),&
                xint(:,2)-xvec(:,out2)) / 6d0)
             call grdall(xp,lprop1)
             ngrd_int = ngrd_int + 1
             lprop1(2:Nprops) = lprop1(2:Nprops) * lvol
             lprop(2:Nprops) = lprop(2:Nprops) + lprop1(2:Nprops)
             ! add out1-out2-2-3
             xp = 0d0
             xp = xp + xvec(:,out1) * ccrd(1)
             xp = xp + xvec(:,out2) * ccrd(2)
             xp = xp + xint(:,2) * ccrd(3)
             xp = xp + xint(:,3) * ccrd(4)
             lvol = abs(mixed(xvec(:,out2)-xvec(:,out1),&
                xint(:,2)-xvec(:,out1),&
                xint(:,3)-xvec(:,out1)) / 6d0)
             call grdall(xp,lprop1)
             ngrd_int = ngrd_int + 1
             lprop1(2:Nprops) = lprop1(2:Nprops) * lvol
             lprop(2:Nprops) = lprop(2:Nprops) + lprop1(2:Nprops)

             acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) + lprop(2:Nprops) * kw(i,klvl)
          end do

       else
          !$omp critical (IO)
          write (uout,'("ins : ",I2)') ins
          write (uout,'("in_vector : ",4(L2,X))') in
          write (uout,'("trm_vector : ",4(I2,X))') ts
          write (uout,'("r : ",F14.9)') sqrt(r2)
          write (uout,'("d1 : ",F14.9)') sqrt(xdot(1,1))
          write (uout,'("d2 : ",F14.9)') sqrt(xdot(2,2))
          write (uout,'("d3 : ",F14.9)') sqrt(xdot(3,3))
          write (uout,'("d4 : ",F14.9)') sqrt(xdot(4,4))
          !$omp end critical (IO)
          call ferror('integ_border_keast','ins /= 1,2,3',faterr)
       end if

    else
       ! integrate the complete tetrahedron
       lprop = 0d0
       lprop(1) = tvol(base_t) / l8
       do i = 1, korder(klvl)
          ccrd(1:3) = kxyz(:,i,klvl)
          ccrd(4) = 1 - sum(ccrd(1:3))
          xp = borig(:,base_t)
          do j = 1, 3
             xp = xp + bvec(:,j,base_t) * dot_product(ccrd,iv(j,:))
          end do
          call grdall(xp,lprop1)
          ngrd_int = ngrd_int + 1
          lprop(2:Nprops) = lprop(2:Nprops) + lprop1(2:Nprops) * tvol(base_t) / l8 * kw(i,klvl)
       end do

       ! assign equally for each vertex sharing the tetrahedron
       unk = count((ts == nnuc+1) .or. (ts == nnuc+3))
       lprop = lprop / (4d0-real(unk,8))
       do i = 1, 4
          if (ts(i) <= 0) then
             call ferror('integ_border_keast','<=0 term in border tetrahedron',faterr)
          end if
          if (ts(i) == nnuc+1 .or. ts(i) == nnuc+3) cycle
          acum_atprop(ts(i),:) = acum_atprop(ts(i),:) + lprop(:)
       end do

    end if

  end subroutine integ_border_keast

  !> Cubpack integration for a tetrahedron on a (external or
  !> internal) interatomic surface
  subroutine integ_border_cubpack(base_t,iv,l,ts,acum_atprop)
    use qtree_basic
    use CUI
    use varbas
    use fields
    use global
    use struct_basic
    use tools_math, only: mixed
    use tools_io

    integer, intent(in) :: base_t
    integer, intent(in) :: iv(3,4)
    integer, intent(in) :: l
    integer, intent(in) :: ts(4)
    real*8, intent(inout) :: acum_atprop(:,:)

    integer :: i, j, k, l8, unk, its
    real*8 :: xp(3), temp(3)
    real*8 :: lprop(2:Nprops), cub_abserr(2:Nprops)
    integer :: ier, nid
    real*8 :: dist, xnuc(3), r2
    real*8 :: xvec(3,4), xint(3,4), xdot(4,4), dcoef, den, lvol
    real*8 :: vvec(3,3), vert(3,4,1)
    integer :: ins, in1, out1, out2, lvec(3), neval
    logical :: in(4)

    l8 = 8**l
    nid = 0
    in1 = 0
    out1 = 0
    out2 = 0
    if (any(ts < 0)) then
       ! a sphere-interior boundary. All the abs(ts) are the same
       its = abs(ts(1))
       acum_atprop(its,1) = acum_atprop(its,1) + tvol(base_t) / l8

       ! find the nucleus, use the barycenter.
       xp = borig(:,base_t)
       do i = 1, 3
          xp = xp + bvec(:,i,base_t) * 0.25d0 * sum(iv(i,:))
       end do
       xp = cr%c2x(xp)
       r2 = 1d30
       do i=1,ncpcel
          if (cpcel(i)%idx /= its) cycle
          temp = cpcel(i)%x - xp
          call cr%shortest(temp,dist)
          if (dist < r2) then
             nid = i
             r2 = dist
             lvec = nint(cpcel(nid)%x - xp - temp)
          end if
       end do
       xnuc = cr%x2c(cpcel(nid)%x - lvec)
       if (abs(sqrt(r2)-r_betaint(its)) > 2d0*maxlen) then
          call ferror('integ_border_keast','unknown xnuc for the tetrahedron',faterr)
       end if
       r2 = r_betaint(its) * r_betaint(its)

       ! build the tetrahedron vertex and inner product matrix
       ins = 0
       in = .false.
       do i = 1, 4
          xvec(:,i) = borig(:,base_t)
          do j = 1, 3
             xvec(:,i) = xvec(:,i) + bvec(:,j,base_t) * iv(j,i)
          end do
          do j = 1, i
             xdot(i,j) = dot_product(xvec(:,i)-xnuc,xvec(:,j)-xnuc)
             xdot(j,i) = xdot(i,j)
          end do
          if (xdot(i,i) < r2) then 
             ins = ins + 1
             in(i) = .true.
          end if
       end do

       if (ins == 1) then
          ! one vertex in the sphere -> 2 integrations.

          ! determine the index of the interior vertex
          do i = 1, 4
             if (in(i)) then
                in1 = i
                exit
             end if
          end do

          ! calculate edge intersections with the sphere
          ! and the small tetrahedron
          k = 0
          do i = 1, 4
             if (in1 == i) cycle
             k = k + 1
             den = xdot(i,i) + xdot(in1,in1) - 2d0*xdot(i,in1)
             dcoef = (-xdot(i,in1)+xdot(in1,in1) + &
                sqrt(xdot(i,in1)**2 - xdot(i,i)*xdot(in1,in1) + r2*den)) / den
             xvec(:,i) = (1-dcoef)*xvec(:,in1) + dcoef * xvec(:,i)
             vvec(:,k) = xvec(:,i) - xvec(:,in1)
             xdot(i,i) = dot_product(xvec(:,i)-xnuc,xvec(:,i)-xnuc)
          end do
          lvol = abs(mixed(vvec(:,1),vvec(:,2),vvec(:,3)) / 6d0)

          ! add the complete tetrahedron
          do i = 1, 4
             vert(:,i,1) = borig(:,base_t)
             do j = 1, 3
                vert(:,i,1) = vert(:,i,1) + bvec(:,j,base_t) * iv(j,i)
             end do
          end do
          ier = 1
          lprop = 0d0
          cub_abserr = 0d0
          call cubatr(3,Nprops-1,cubpack_f,1,vert,(/1/),lprop,cub_abserr,&
             IFAIL=ier, EpsAbs=cub_abs, EpsRel=cub_rel, MaxPts=cub_mpts, Neval=neval)
          ngrd_int = ngrd_int + neval
          if (ier == 1) then
             !$omp critical (IO)
             write (uout,*) "too many evaluations in cubpack"
             !$omp end critical (IO)
          else if (ier > 1) then
             call ferror('integ_border_cubpack','severe error in cubpack',faterr)
          end if
          acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) + lprop

          ! substract the tetrahedron inside the sphere
          vert(:,:,1) = xvec
          ier = 1
          lprop = 0d0
          cub_abserr = 0d0
          call cubatr(3,Nprops-1,cubpack_f,1,vert,(/1/),lprop,cub_abserr,&
             IFAIL=ier, EpsAbs=cub_abs, EpsRel=cub_rel, MaxPts=cub_mpts, Neval=neval)
          ngrd_int = ngrd_int + neval
          if (ier == 1) then
             !$omp critical (IO)
             write (uout,*) "too many evaluations in cubpack"
             !$omp end critical (IO)
          else if (ier > 1) then
             call ferror('integ_border_cubpack','severe error in cubpack',faterr)
          end if
          acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) - lprop

       else if (ins == 3) then
          ! one vertex out of the sphere -> 1 integration.

          ! determine the index of the exterior vertex
          do i = 1, 4
             if (.not.in(i)) then
                out1 = i
                exit
             end if
          end do

          ! calculate edge intersections with the sphere
          ! and the small tetrahedron
          k = 0
          do i = 1, 4
             if (out1 == i) cycle
             k = k + 1
             den = xdot(i,i) + xdot(out1,out1) - 2d0*xdot(i,out1)
             dcoef = (-xdot(i,out1)+xdot(out1,out1) - &
                sqrt(xdot(i,out1)**2 - xdot(i,i)*xdot(out1,out1) + r2*den)) / den
             xvec(:,i) = (1-dcoef)*xvec(:,out1) + dcoef * xvec(:,i)
             vvec(:,k) = xvec(:,i) - xvec(:,out1)
          end do
          lvol = abs(mixed(vvec(:,1),vvec(:,2),vvec(:,3)) / 6d0)

          ! integrate 
          vert(:,:,1) = xvec
          ier = 1
          lprop = 0d0
          cub_abserr = 0d0
          call cubatr(3,Nprops-1,cubpack_f,1,vert,(/1/),lprop,cub_abserr,&
             IFAIL=ier, EpsAbs=cub_abs, EpsRel=cub_rel, MaxPts=cub_mpts, Neval=neval)
          ngrd_int = ngrd_int + neval
          if (ier == 1) then
             !$omp critical (IO)
             write (uout,*) "too many evaluations in cubpack"
             !$omp end critical (IO)
          else if (ier > 1) then
             call ferror('integ_border_cubpack','severe error in cubpack',faterr)
          end if
          acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) + lprop

       else if (ins == 2) then
          ! two vertex out of the sphere -> 1 integration.

          ! determine the index of the exterior vertex
          out1 = 0
          do i = 1, 4
             if (.not.in(i)) then
                if (out1 == 0) then
                   out1 = i
                else
                   out2 = i
                   exit
                end if
             end if
          end do

          ! calculate edge intersections of the vertex out1
          ! xint contains the four points, 1 and 2 associated to out1, 3 and 4 to out2
          ! 2 and 3 are opposite vertex (the interior vertex are different)
          k = 0
          do i = 1, 4
             if (out1 == i .or. out2 == i) cycle
             k = k + 1
             den = xdot(i,i) + xdot(out1,out1) - 2d0*xdot(i,out1)
             dcoef = (-xdot(i,out1)+xdot(out1,out1) - &
                sqrt(xdot(i,out1)**2 - xdot(i,i)*xdot(out1,out1) + r2*den)) / den
             xint(:,k) = (1-dcoef)*xvec(:,out1) + dcoef*xvec(:,i)
          end do
          do i = 1, 4
             if (out1 == i .or. out2 == i) cycle
             k = k + 1
             den = xdot(i,i) + xdot(out2,out2) - 2d0*xdot(i,out2)
             dcoef = (-xdot(i,out2)+xdot(out2,out2) - &
                sqrt(xdot(i,out2)**2 - xdot(i,i)*xdot(out2,out2) + r2*den)) / den
             xint(:,k) = (1-dcoef)*xvec(:,out2) + dcoef*xvec(:,i)
          end do

          ! integrate 
          ! add out1-1-2-3
          vert(:,1,1) = xvec(:,out1)
          vert(:,2,1) = xint(:,1)
          vert(:,3,1) = xint(:,2)
          vert(:,4,1) = xint(:,3)
          ier = 1
          lprop = 0d0
          cub_abserr = 0d0
          call cubatr(3,Nprops-1,cubpack_f,1,vert,(/1/),lprop,cub_abserr,&
             IFAIL=ier, EpsAbs=cub_abs, EpsRel=cub_rel, MaxPts=cub_mpts, Neval=neval)
          ngrd_int = ngrd_int + neval
          if (ier == 1) then
             !$omp critical (IO)
             write (uout,*) "too many evaluations in cubpack"
             !$omp end critical (IO)
          else if (ier > 1) then
             call ferror('integ_border_cubpack','severe error in cubpack',faterr)
          end if
          acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) + lprop
          ! add out2-3-4-2
          vert(:,1,1) = xvec(:,out2)
          vert(:,2,1) = xint(:,3)
          vert(:,3,1) = xint(:,4)
          vert(:,4,1) = xint(:,2)
          ier = 1
          lprop = 0d0
          cub_abserr = 0d0
          call cubatr(3,Nprops-1,cubpack_f,1,vert,(/1/),lprop,cub_abserr,&
             IFAIL=ier, EpsAbs=cub_abs, EpsRel=cub_rel, MaxPts=cub_mpts, Neval=neval)
          ngrd_int = ngrd_int + neval
          if (ier == 1) then
             !$omp critical (IO)
             write (uout,*) "too many evaluations in cubpack"
             !$omp end critical (IO)
          else if (ier > 1) then
             call ferror('integ_border_cubpack','severe error in cubpack',faterr)
          end if

          acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) + lprop
          ! add out1-out2-2-3
          vert(:,1,1) = xvec(:,out1)
          vert(:,2,1) = xvec(:,out2)
          vert(:,3,1) = xint(:,2)
          vert(:,4,1) = xint(:,3)
          ier = 1
          lprop = 0d0
          cub_abserr = 0d0
          call cubatr(3,Nprops-1,cubpack_f,1,vert,(/1/),lprop,cub_abserr,&
             IFAIL=ier, EpsAbs=cub_abs, EpsRel=cub_rel, MaxPts=cub_mpts, Neval=neval)
          ngrd_int = ngrd_int + neval
          if (ier == 1) then
             !$omp critical (IO)
             write (uout,*) "too many evaluations in cubpack"
             !$omp end critical (IO)
          else if (ier > 1) then
             call ferror('integ_border_cubpack','severe error in cubpack',faterr)
          end if
          acum_atprop(its,2:Nprops) = acum_atprop(its,2:Nprops) + lprop

       else
          !$omp critical (IO)
          write (uout,'("ins : ",I2)') ins
          write (uout,'("in_vector : ",4(L2,X))') in
          write (uout,'("trm_vector : ",4(I2,X))') ts
          write (uout,'("r : ",F14.9)') sqrt(r2)
          write (uout,'("d1 : ",F14.9)') sqrt(xdot(1,1))
          write (uout,'("d2 : ",F14.9)') sqrt(xdot(2,2))
          write (uout,'("d3 : ",F14.9)') sqrt(xdot(3,3))
          write (uout,'("d4 : ",F14.9)') sqrt(xdot(4,4))
          !$omp end critical (IO)
          call ferror('integ_border_cubpack','ins /= 1,2,3',faterr)
       end if

    else

       ! add the complete tetrahedron
       do i = 1, 4
          vert(:,i,1) = borig(:,base_t)
          do j = 1, 3
             vert(:,i,1) = vert(:,i,1) + bvec(:,j,base_t) * iv(j,i)
          end do
       end do
       ier = 1
       lprop = 0d0
       cub_abserr = 0d0
       call cubatr(3,Nprops-1,cubpack_f,1,vert,(/1/),lprop,cub_abserr,&
          IFAIL=ier, EpsAbs=cub_abs, EpsRel=cub_rel, MaxPts=cub_mpts, Neval=neval)
       ngrd_int = ngrd_int + neval
       if (ier == 1) then
          !$omp critical (IO)
          write (uout,*) "too many evaluations in cubpack"
          !$omp end critical (IO)
       else if (ier > 1) then
          call ferror('integ_border_cubpack','severe error in cubpack',faterr)
       end if

       ! assign equally for each vertex sharing the tetrahedron
       unk = count((ts == nnuc+1) .or. (ts == nnuc+3))
       lprop = lprop / (4d0-real(unk,8))
       do i = 1, 4
          if (ts(i) <= 0) then
             call ferror('integ_border_cubpack','<=0 term in border tetrahedron',faterr)
          end if
          if (ts(i) == nnuc+1 .or. ts(i) == nnuc+3) cycle
          acum_atprop(ts(i),1) = acum_atprop(ts(i),1) + tvol(base_t) / l8 / (4d0-real(unk,8))
          acum_atprop(ts(i),2:Nprops) = acum_atprop(ts(i),2:Nprops) + lprop
       end do

    end if

  end subroutine integ_border_cubpack

  !> Vertex integration for any tetrahedron.
  subroutine integ_corner(base_t,iv,l,ts,acum_atprop,fgr,lapgr)
    use qtree_basic
    use varbas
    use fields
    use global
    use tools_io
    use types

    integer, intent(in) :: base_t
    integer, intent(in) :: iv(3,4)
    integer, intent(in) :: l
    integer, intent(in) :: ts(4)
    real(qtreer), intent(inout) :: fgr(:,:), lapgr(:,:)
    real*8, intent(inout) :: acum_atprop(:,:)

    integer :: unk
    integer :: i, j, lrest, base_to
    real*8 :: lprop(Nprops), vfac, xp(3), lf(4), llap(4)
    integer(qtreeidx) :: idx(4)
    type(scalar_value) :: res

    if (color_allocate == 0) then
       base_to = 1
    else
       base_to = base_t
    end if
    lrest = 2**(maxl-l)
    lprop = 0d0
    do i = 1, 4
       idx(i) = cindex(iv(:,i),l)
       if ((prop_mode == 1 .or. prop_mode == 2)) then
          lf(i) = -1d0
          llap(i) = 0d0
          if (savefgr) then
             lf(i) = fgr(idx(i),base_to)
          end if
          if (savelapgr) then
             llap(i) = lapgr(idx(i),base_to)
          end if

          if (ts(i) > 0 .and. lf(i) < 0d0) then
             xp = borig(:,base_t)
             do j = 1, 3
                xp = xp + bvec(:,j,base_t) * iv(j,i) * lrest
             end do
             call grd(f(refden),xp,nder,res)
             lf(i) = res%fval
             if (savefgr) then
                fgr(idx(i),base_to) = res%fval
             end if
             if (savelapgr) then
                llap(i) = -res%del2fval
                lapgr(idx(i),base_to) = -res%del2fval
             end if
          end if
       end if
    end do

    unk = count((ts == nnuc+1) .or. (ts == nnuc+3))
    if (unk == 4) return
    vfac = tvol(base_t) / 8**l / (4d0-real(unk,8))

    ! sum volume
    do i = 1, 4
       if (ts(i) /= 0) then
          acum_atprop(abs(ts(i)),1) = acum_atprop(abs(ts(i)),1) + vfac
       end if
    end do

    if (prop_mode == 1) then
       do i = 1, 4
          if (ts(i) > 0) then
             acum_atprop(ts(i),2) = acum_atprop(ts(i),2) + vfac * lf(i)
          end if
       end do
    else if (prop_mode == 2) then
       do i = 1, 4
          if (ts(i) > 0) then
             acum_atprop(ts(i),2) = acum_atprop(ts(i),2) + vfac * lf(i)
             acum_atprop(ts(i),3) = acum_atprop(ts(i),3) + vfac * llap(i)
          end if
       end do
    else if (prop_mode == 3) then
       do i = 1, 4
          if (ts(i) > 0) then
             xp = borig(:,base_t)
             do j = 1, 3
                xp = xp + bvec(:,j,base_t) * iv(j,i) * lrest
             end do
             call grdall(xp,lprop)
             ngrd_int = ngrd_int + 1
             acum_atprop(ts(i),2:Nprops) = acum_atprop(ts(i),2:Nprops) + vfac * lprop(2:Nprops)
          end if
       end do
    else
       call ferror("integ_corner","wrong integ option",faterr)
    end if

  end subroutine integ_corner

  !> Vertex integration for any tetrahedron. This routine merely
  !> writes down the volume associated to each corner. The sum
  !> is calculated after the tetrahedron painting. Used if the flag
  !> intcorner_deferred is active (all int modes equal to -1 or 11).
  subroutine integ_corner_deferred(base_t,iv,l,ts,vgr)
    use qtree_basic
    use tools_io

    integer, intent(in) :: base_t
    integer, intent(in) :: iv(3,4)
    integer, intent(in) :: l
    integer, intent(in) :: ts(4)
    real(qtreer), intent(inout) :: vgr(:)

    integer :: unk, i
    integer(qtreeidx) :: idx
    real*8 :: vfac

    unk = count((ts == nnuc+1) .or. (ts == nnuc+3))
    if (unk == 4) return
    vfac = tvol(base_t) / 8**l / (4d0-real(unk,8))

    do i = 1, 4
       idx = cindex(iv(:,i),l)
       vgr(idx) = vgr(idx) + vfac
    end do

  end subroutine integ_corner_deferred

  !> Sum the contribution of each grid to the atomic properties. Executed
  !> after the color assignment.
  subroutine integ_corner_sum(base_t,trm,vgr,acum_atprop)
    use qtree_basic
    use varbas
    use fields
    use global
    use tools_io
    use types
    
    integer, intent(in) :: base_t
    integer(qtreei), intent(inout) :: trm(:,:)
    real(qtreer), intent(inout) :: vgr(:)
    real*8, intent(inout) :: acum_atprop(:,:)
    
    integer :: i, j, k, l2, vin(3)
    real*8 :: xx(3), lprop(Nprops), vfac
    integer :: ts, base_to
    integer(qtreeidx) :: idx
    type(scalar_value) :: res

    real*8, parameter :: vfac_low = 1d-15

    if (color_allocate == 0) then
       base_to = 1
    else
       base_to = base_t
    end if
    l2 = 2**maxl
    lprop = 0d0

    do i = 0, l2
       do j = 0, l2-i
          do k = 0, l2-i-j
             vin = (/i,j,k/)
             idx = cindex(vin,maxl)
             ts = trm(idx,base_to)
             vfac = vgr(idx)
             if (vfac < vfac_low) cycle

             if (ts /= 0) then
                acum_atprop(abs(ts),1) = acum_atprop(abs(ts),1) + vfac
             end if

             if (ts > 0) then
                xx = borig(:,base_t)
                xx = xx + bvec(:,1,base_t) * real(i,8)
                xx = xx + bvec(:,2,base_t) * real(j,8)
                xx = xx + bvec(:,3,base_t) * real(k,8)
                if (prop_mode == 1) then
                   call grd(f(refden),xx,0,res)
                   ngrd_int = ngrd_int + 1
                   acum_atprop(ts,2) = acum_atprop(ts,2) + vfac * res%fval
                else if (prop_mode == 2) then
                   call grd(f(refden),xx,2,res)
                   acum_atprop(ts,2) = acum_atprop(ts,2) + vfac * res%fval
                   acum_atprop(ts,3) = acum_atprop(ts,3) + vfac * res%del2fval
                   ngrd_int = ngrd_int + 1
                else if (prop_mode == 3) then
                   call grdall(xx,lprop)
                   ngrd_int = ngrd_int + 1
                   acum_atprop(ts,2:Nprops) = acum_atprop(ts,2:Nprops) + vfac * lprop(2:Nprops)
                else
                   call ferror("integ_corner_sum","wrong integ option",faterr)
                end if
             end if
          end do
       end do
    end do

  end subroutine integ_corner_sum

  !> Paint grid points from one IWST that are inside a beta-sphere. 
  subroutine paint_inside_spheres(tt,tto,trm)
    use qtree_basic
    use fields
    use global
    use varbas

    integer, intent(in) :: tt, tto
    integer(qtreei), intent(inout) :: trm(:,:)

    integer :: i, j, k, l2
    real*8 :: xx(3)
    integer :: nid, vin(3)
    real*8 :: dist
    integer(qtreeidx) :: idx

    l2 = 2**maxl

    do i = 0, l2
       do j = 0, l2-i
          do k = 0, l2-i-j
             xx = torig(:,tt)
             xx = xx + tvec(:,1,tt) * real(i,8) / l2
             xx = xx + tvec(:,2,tt) * real(j,8) / l2
             xx = xx + tvec(:,3,tt) * real(k,8) / l2
             call nearest_cp(xx,nid,dist,type=f(refden)%typnuc)
             if (dist <= r_betagp(cpcel(nid)%idx)) then
                vin = (/i,j,k/)
                idx = cindex(vin,maxl)
                if (dist <= r_betaint(cpcel(nid)%idx)) then
                   trm(idx,tto) = int(-cpcel(nid)%idx,1)
                else
                   trm(idx,tto) = int(cpcel(nid)%idx,1)
                end if
             end if
          end do
       end do
    end do
    
  end subroutine paint_inside_spheres

  !> Wrapper function for cubpack.
  function cubpack_f(numfun,x) result(value)
    use qtree_basic
    use fields
    use global
    USE Precision_Model
    integer, intent(in) :: numfun
    real(kind=stnd), dimension(:), intent(in) :: x
    real(kind=stnd), dimension(numfun) :: value

    real*8 :: lprop(Nprops)

    call grdall(x,lprop)
    value = lprop(2:Nprops)
    
  end function cubpack_f

end module qtree_tetrawork
