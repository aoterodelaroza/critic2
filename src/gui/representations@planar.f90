! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Planar shapes (reptype_planar): the geometry of the shapes drawn on
! the screen, their editing handles and hit tests, and their
! tessellation into the triangles of the draw list. Everything is in
! the NDC of the render buffer.
submodule (representations) planar
  implicit none

  ! ellipse outlines: segment count bounds, and the target NDC length of
  ! a segment (about two pixels in a 1024-pixel render buffer)
  integer, parameter :: nseg_ellipse_min = 32
  integer, parameter :: nseg_ellipse_max = 512
  real*8, parameter :: seglen_ellipse = 0.004d0

  ! curves: minimum number of segments (the target length is that of the ellipses)
  integer, parameter :: nseg_curve_min = 8

  ! dashed outlines: dash and gap lengths; dotted outlines: distance
  ! between the dot centers (the dots are as wide as the outline). All in
  ! outline widths, before they are stretched to fit the path.
  real*8, parameter :: dash_on = 4d0
  real*8, parameter :: dash_off = 2.5d0
  real*8, parameter :: dot_step = 2d0

  ! segments in a round stroke join or cap
  integer, parameter :: nseg_disk = 16

  ! longest miter join, in half-widths; sharper corners are rounded
  real*8, parameter :: miter_limit = 4d0

  ! distance from the top edge to the rotation handle of an ellipse or rectangle
  real*8, parameter :: rothandle_off = 0.08d0

  ! depth (NDC) of the first planar shape, and the step to the next one
  ! (each later shape a bit closer, so it covers the ones before it)
  real*8, parameter :: zfirst = 1d0 - 2d-6
  real*8, parameter :: zstep = 4d-6

  ! signs of the corners (handles 1-4) and edge midpoints (handles 5-8)
  ! of an ellipse or rectangle, in its own axes
  real*8, parameter :: hsign(2,8) = reshape((/&
     -1d0,-1d0,  1d0,-1d0,  1d0, 1d0, -1d0, 1d0,&
      0d0,-1d0,  1d0, 0d0,  0d0, 1d0, -1d0, 0d0/),shape(hsign))

contains

  !> Whether the shape is closed (it has an inside that can be filled).
  module function planar_isclosed(sh) result(ok)
    type(planar_shape), intent(in) :: sh
    logical :: ok

    ok = (sh%kind == planarkind_ellipse .or. sh%kind == planarkind_rect .or.&
       sh%kind == planarkind_polygon)

  end function planar_isclosed

  !> Whether the shape is given by its list of points (sh%x), rather than
  !> by a center, half-sizes and an angle.
  module function planar_haspoints(sh) result(ok)
    type(planar_shape), intent(in) :: sh
    logical :: ok

    ok = (sh%kind == planarkind_polygon .or. sh%kind == planarkind_polyline .or.&
       sh%kind == planarkind_arrow .or. sh%kind == planarkind_curve .or.&
       sh%kind == planarkind_freehand)

  end function planar_haspoints

  !> The center line of the shape sh: the n points in x(2,n), joined in
  !> order (and the last back to the first if the shape is closed). An
  !> ellipse is sampled finely enough to look smooth.
  module subroutine planar_path(sh,n,x)
    use param, only: tpi
    type(planar_shape), intent(in) :: sh
    integer, intent(out) :: n
    real*8, allocatable, intent(inout) :: x(:,:)

    real*8 :: u(2), v(2), r, t, c(2)
    integer :: i

    n = 0
    if (sh%kind == planarkind_ellipse .or. sh%kind == planarkind_rect) then
       u = (/cos(sh%ang),sin(sh%ang)/)
       v = (/-u(2),u(1)/)
       if (sh%kind == planarkind_ellipse) then
          r = max(abs(sh%hs(1)),abs(sh%hs(2)))
          n = min(max(ceiling(tpi * r / seglen_ellipse),nseg_ellipse_min),nseg_ellipse_max)
          if (allocated(x)) deallocate(x)
          allocate(x(2,n))
          do i = 1, n
             t = tpi * real(i-1,8) / real(n,8)
             x(:,i) = sh%xc + sh%hs(1) * cos(t) * u + sh%hs(2) * sin(t) * v
          end do
       else
          n = 4
          if (allocated(x)) deallocate(x)
          allocate(x(2,n))
          do i = 1, 4
             x(:,i) = sh%xc + hsign(1,i) * sh%hs(1) * u + hsign(2,i) * sh%hs(2) * v
          end do
       end if
    elseif (sh%kind == planarkind_curve .and. sh%npt == 3 .and. allocated(sh%x)) then
       ! the quadratic Bezier through the middle point: its control point
       ! is c, and it passes through x(:,2) at t = 1/2
       c = 2d0 * sh%x(:,2) - 0.5d0 * (sh%x(:,1) + sh%x(:,3))
       r = norm2(c - sh%x(:,1)) + norm2(sh%x(:,3) - c)
       n = min(max(ceiling(r / seglen_ellipse),nseg_curve_min),nseg_ellipse_max) + 1
       if (allocated(x)) deallocate(x)
       allocate(x(2,n))
       do i = 1, n
          t = real(i-1,8) / real(n-1,8)
          x(:,i) = (1d0-t)**2 * sh%x(:,1) + 2d0 * t * (1d0-t) * c + t**2 * sh%x(:,3)
       end do
    elseif (sh%npt > 0 .and. allocated(sh%x)) then
       n = sh%npt
       if (allocated(x)) deallocate(x)
       allocate(x(2,n))
       x = sh%x(:,1:n)
    end if

  end subroutine planar_path

  !> Place the middle point of the curve sh at the default bend: off the
  !> middle of the chord between its ends, to the left of it.
  module subroutine planar_curve_default_bend(sh)
    type(planar_shape), intent(inout) :: sh

    real*8 :: d(2)

    if (sh%npt < 3 .or. .not.allocated(sh%x)) return
    d = sh%x(:,3) - sh%x(:,1)
    sh%x(:,2) = 0.5d0 * (sh%x(:,1) + sh%x(:,3)) + planar_curve_bend_def * (/-d(2),d(1)/)

  end subroutine planar_curve_default_bend

  !> The editing handles of the shape sh: nh positions in xh(2,nh). An
  !> ellipse or rectangle has its four corners (1-4), the midpoints of
  !> its four edges (5-8), and a rotation handle (9) beyond the top
  !> edge. A shape given by points has one handle per point, except a
  !> freehand line, which has none (it can only be moved).
  module subroutine planar_handles(sh,nh,xh)
    type(planar_shape), intent(in) :: sh
    integer, intent(out) :: nh
    real*8, allocatable, intent(inout) :: xh(:,:)

    real*8 :: u(2), v(2)
    integer :: i

    nh = 0
    if (sh%kind == planarkind_ellipse .or. sh%kind == planarkind_rect) then
       u = (/cos(sh%ang),sin(sh%ang)/)
       v = (/-u(2),u(1)/)
       nh = 9
       if (allocated(xh)) deallocate(xh)
       allocate(xh(2,nh))
       do i = 1, 8
          xh(:,i) = sh%xc + hsign(1,i) * sh%hs(1) * u + hsign(2,i) * sh%hs(2) * v
       end do
       xh(:,9) = sh%xc + (abs(sh%hs(2)) + rothandle_off) * v
    elseif (sh%kind /= planarkind_freehand .and. sh%npt > 0 .and. allocated(sh%x)) then
       nh = sh%npt
       if (allocated(xh)) deallocate(xh)
       allocate(xh(2,nh))
       xh = sh%x(:,1:nh)
    end if

  end subroutine planar_handles

  !> The rotation handle of the shape sh (0 if it has none).
  module function planar_rotation_handle(sh) result(ih)
    type(planar_shape), intent(in) :: sh
    integer :: ih

    ih = 0
    if (sh%kind == planarkind_ellipse .or. sh%kind == planarkind_rect) ih = 9

  end function planar_rotation_handle

  !> Drag handle ih (see planar_handles) of the shape to position x. sh0
  !> is the shape as it was when the drag started, and the result goes to
  !> sh. A corner moves with the opposite corner fixed, and an edge
  !> midpoint with the opposite edge fixed. If constrain: a corner keeps
  !> the two sizes equal (circle, square), an edge scales both, the
  !> rotation snaps to 15 degrees, and a point snaps the segment to its
  !> neighbor to multiples of 45 degrees.
  module subroutine planar_drag_handle(sh,sh0,ih,x,constrain)
    use param, only: pi
    type(planar_shape), intent(inout) :: sh
    type(planar_shape), intent(in) :: sh0
    integer, intent(in) :: ih
    real*8, intent(in) :: x(2)
    logical, intent(in) :: constrain

    real*8 :: u(2), v(2), f(2), d(2), dl(2), a(2), d1, m, ang, ref(2)
    integer :: k

    real*8, parameter :: snaprot = 15d0 * pi / 180d0
    real*8, parameter :: snappt = 45d0 * pi / 180d0

    if (sh0%kind == planarkind_ellipse .or. sh0%kind == planarkind_rect) then
       u = (/cos(sh0%ang),sin(sh0%ang)/)
       v = (/-u(2),u(1)/)
       if (ih >= 1 .and. ih <= 4) then
          ! corner: the opposite corner stays
          f = sh0%xc - hsign(1,ih) * sh0%hs(1) * u - hsign(2,ih) * sh0%hs(2) * v
          d = x - f
          dl = (/dot_product(d,u),dot_product(d,v)/)
          if (constrain) then
             m = max(abs(dl(1)),abs(dl(2)))
             dl = sign(m,dl)
          end if
          sh%hs = 0.5d0 * abs(dl)
          sh%xc = f + 0.5d0 * (dl(1) * u + dl(2) * v)
       elseif (ih >= 5 .and. ih <= 8) then
          ! edge midpoint: the opposite edge stays
          if (hsign(1,ih) /= 0d0) then
             k = 1
             a = hsign(1,ih) * u
          else
             k = 2
             a = hsign(2,ih) * v
          end if
          f = sh0%xc - sh0%hs(k) * a
          d1 = dot_product(x - f,a)
          sh%hs(k) = 0.5d0 * abs(d1)
          if (constrain) sh%hs(3-k) = sh%hs(k)
          sh%xc = f + 0.5d0 * d1 * a
       elseif (ih == 9) then
          ! rotation: the handle sits on the +v axis
          d = x - sh0%xc
          if (norm2(d) < 1d-10) return
          ang = atan2(d(2),d(1)) - 0.5d0 * pi
          if (constrain) ang = snaprot * anint(ang / snaprot)
          sh%ang = modulo(ang + pi,2d0*pi) - pi
       end if
    elseif (sh0%kind == planarkind_curve .and. ih == 2 .and. sh0%npt == 3) then
       ! the middle of a curve; constrained, it stays on the perpendicular
       ! bisector of the chord (a symmetric bend)
       sh%x(:,2) = x
       d = sh0%x(:,3) - sh0%x(:,1)
       if (constrain .and. norm2(d) > 1d-10) then
          f = 0.5d0 * (sh0%x(:,1) + sh0%x(:,3))
          a = (/-d(2),d(1)/) / norm2(d)
          sh%x(:,2) = f + dot_product(x - f,a) * a
       end if
    elseif (sh0%kind == planarkind_curve .and. sh0%npt == 3) then
       ! an end of a curve; constrained, the chord snaps to 45 degrees
       sh%x(:,ih) = x
       if (constrain) then
          ref = sh0%x(:,4-ih)
          d = x - ref
          if (norm2(d) > 1d-10) then
             ang = snappt * anint(atan2(d(2),d(1)) / snappt)
             sh%x(:,ih) = ref + norm2(d) * (/cos(ang),sin(ang)/)
          end if
       end if
    elseif (ih >= 1 .and. ih <= sh0%npt) then
       sh%x(:,ih) = x
       if (constrain .and. sh0%npt >= 2) then
          ! the segment from the previous point (the next, for the first)
          if (ih > 1) then
             ref = sh0%x(:,ih-1)
          else
             ref = sh0%x(:,2)
          end if
          d = x - ref
          if (norm2(d) > 1d-10) then
             ang = snappt * anint(atan2(d(2),d(1)) / snappt)
             sh%x(:,ih) = ref + norm2(d) * (/cos(ang),sin(ang)/)
          end if
       end if
    end if

  end subroutine planar_drag_handle

  !> Translate the shape sh by d.
  module subroutine planar_move(sh,d)
    type(planar_shape), intent(inout) :: sh
    real*8, intent(in) :: d(2)

    integer :: i

    sh%xc = sh%xc + d
    if (allocated(sh%x)) then
       do i = 1, sh%npt
          sh%x(:,i) = sh%x(:,i) + d
       end do
    end if

  end subroutine planar_move

  !> Whether position x hits the shape sh. If onstroke, test whether x
  !> is on the outline, up to a distance tol beyond its edge; otherwise,
  !> whether it is inside a closed shape (filled or not).
  module function planar_hit(sh,x,tol,onstroke) result(ok)
    type(planar_shape), intent(in) :: sh
    real*8, intent(in) :: x(2)
    real*8, intent(in) :: tol
    logical, intent(in) :: onstroke
    logical :: ok

    integer :: n, i, j, nseg
    real*8, allocatable :: xp(:,:)
    real*8 :: dmax
    logical :: closed

    ok = .false.
    call planar_path(sh,n,xp)
    if (n == 0) return
    closed = planar_isclosed(sh)

    if (onstroke) then
       dmax = tol
       if (sh%stroke) dmax = dmax + 0.5d0 * sh%width
       if (n == 1) then
          ok = (norm2(x - xp(:,1)) <= dmax)
          return
       end if
       nseg = n - 1
       if (closed) nseg = n
       do i = 1, nseg
          j = modulo(i,n) + 1
          if (segdist(x,xp(:,i),xp(:,j)) <= dmax) then
             ok = .true.
             return
          end if
       end do
    elseif (closed .and. n >= 3) then
       ok = inside_polygon(x,n,xp)
    end if

  end function planar_hit

  !> Simplify the polyline of n points x(2,n) (Ramer-Douglas-Peucker):
  !> drop the points closer than tol to the line through their kept
  !> neighbors. The end points are always kept; n is updated and the
  !> kept points are packed at the start of x.
  module subroutine planar_simplify(n,x,tol)
    integer, intent(inout) :: n
    real*8, intent(inout) :: x(:,:)
    real*8, intent(in) :: tol

    logical, allocatable :: keep(:)
    integer, allocatable :: stk(:,:)
    integer :: ns, i1, i2, k, imax, m
    real*8 :: d, dmax

    if (n <= 2) return
    allocate(keep(n),stk(2,2*n))
    keep = .false.
    keep(1) = .true.
    keep(n) = .true.
    ns = 1
    stk(:,1) = (/1,n/)
    do while (ns > 0)
       i1 = stk(1,ns)
       i2 = stk(2,ns)
       ns = ns - 1
       imax = 0
       dmax = -1d0
       do k = i1+1, i2-1
          d = segdist(x(:,k),x(:,i1),x(:,i2))
          if (d > dmax) then
             dmax = d
             imax = k
          end if
       end do
       if (imax > 0 .and. dmax > tol) then
          keep(imax) = .true.
          stk(:,ns+1) = (/i1,imax/)
          stk(:,ns+2) = (/imax,i2/)
          ns = ns + 2
       end if
    end do

    m = 0
    do k = 1, n
       if (keep(k)) then
          m = m + 1
          x(:,m) = x(:,k)
       end if
    end do
    n = m

  end subroutine planar_simplify

  !> A new, shown shape with the style of the selected shape in p (the
  !> defaults if none is selected), so consecutive shapes look alike.
  !> Its geometry is empty (no points, zero size), for the caller to set.
  module function planar_template(p) result(sh)
    type(rep_planar), intent(in) :: p
    type(planar_shape) :: sh

    if (p%isel >= 1 .and. p%isel <= p%nshape) then
       sh = p%shape(p%isel)
    else
       sh = planar_shape()
    end if
    sh%shown = .true.
    sh%xc = 0d0
    sh%hs = 0d0
    sh%ang = 0d0
    sh%npt = 0
    if (allocated(sh%x)) deallocate(sh%x)

  end function planar_template

  !> Append the shape sh to the list in p and select it.
  module subroutine planar_append(p,sh)
    type(rep_planar), intent(inout) :: p
    type(planar_shape), intent(in) :: sh

    type(planar_shape), allocatable :: aux(:)

    allocate(aux(p%nshape+1))
    if (p%nshape > 0) aux(1:p%nshape) = p%shape(1:p%nshape)
    aux(p%nshape+1) = sh
    call move_alloc(aux,p%shape)
    p%nshape = p%nshape + 1
    p%isel = p%nshape

  end subroutine planar_append

  !> Remove shape idel from the list in p. The selection stays on the
  !> shape it was on, or is cleared if that was the one removed.
  module subroutine planar_delete(p,idel)
    type(rep_planar), intent(inout) :: p
    integer, intent(in) :: idel

    integer :: i

    if (idel < 1 .or. idel > p%nshape) return
    do i = idel, p%nshape-1
       p%shape(i) = p%shape(i+1)
    end do
    p%nshape = p%nshape - 1
    if (p%isel == idel) then
       p%isel = 0
    elseif (p%isel > idel) then
       p%isel = p%isel - 1
    end if

  end subroutine planar_delete

  !> Tessellate the shape sh into triangles and append them to the
  !> planar draw lists of obj (on top of everything or behind the
  !> scene). The outline is drawn with miter joins (rounded where the
  !> corner is too sharp) and round caps; the arrowheads are triangles.
  !> Every triangle of the outline is at the same depth, so its overlaps
  !> are drawn once; the fill goes first, just behind it, so that a
  !> translucent outline blends over the fill.
  module subroutine planar_tessellate(sh,obj)
    type(planar_shape), intent(in) :: sh
    type(scene_objects), intent(inout) :: obj

    integer :: n, i
    real*8, allocatable :: xp(:,:)
    real*8 :: z, hw
    real(c_float) :: rgba(4)
    logical :: closed

    call planar_path(sh,n,xp)
    if (n == 0) return
    closed = planar_isclosed(sh)
    if ((sh%kind == planarkind_arrow .or. sh%kind == planarkind_curve) .and. n < 2) return

    ! the depth of this shape
    obj%nflatshape = obj%nflatshape + 1
    z = max(zfirst - zstep * (obj%nflatshape - 1),-0.999d0)
    hw = 0.5d0 * sh%width

    ! the fill, at the depth of the shape (the outline goes half a step closer)
    if (closed .and. sh%fill .and. sh%fillalpha > 0._c_float .and. n >= 3) then
       rgba = (/sh%fillrgb,sh%fillalpha/)
       if (sh%kind == planarkind_polygon) then
          call fill_polygon(n,xp)
       else
          ! ellipse and rectangle: convex
          do i = 2, n-1
             call emit_tri(xp(:,1),xp(:,i),xp(:,i+1))
          end do
       end if
    end if
    z = z - 0.5d0 * zstep

    ! the outline
    if (sh%stroke .and. hw > 0d0 .and. sh%alpha > 0._c_float) then
       rgba = (/sh%rgb,sh%alpha/)
       if (sh%kind == planarkind_arrow .or. sh%kind == planarkind_curve) then
          call stroke_with_heads(n,xp)
       else
          call stroke_styled(n,xp,closed,.true.,.true.)
       end if
    end if


  contains
    !> Stroke the open path of m points x (an arrow or a curve) with
    !> arrowheads at the ends sh%heads says. The line stops at the base
    !> of each head, measured along the path, so a head sits tangent to a
    !> curve; the heads shrink if the path is shorter than both together.
    subroutine stroke_with_heads(m,x)
      integer, intent(in) :: m
      real*8, intent(in) :: x(2,m)

      real*8 :: seg(m-1), ltot, hl(2), hwh, scal, q1(2), q2(2), tip(2), base(2), u(2)
      real*8, allocatable :: y(:,:)
      integer :: k, k1, k2, ny
      logical :: hashead(2)

      do k = 1, m-1
         seg(k) = norm2(x(:,k+1) - x(:,k))
      end do
      ltot = sum(seg)
      if (ltot < 1d-12) return
      hashead(1) = (sh%heads == planarheads_start .or. sh%heads == planarheads_both)
      hashead(2) = (sh%heads == planarheads_end .or. sh%heads == planarheads_both)
      hl = merge(sh%headl * sh%width,0d0,hashead)
      scal = 1d0
      if (sum(hl) > ltot) scal = ltot / sum(hl)
      hl = hl * scal
      hwh = 0.5d0 * sh%headw * sh%width * scal

      ! the line, between the bases of the heads
      call path_point(m,x,seg,hl(1),q1,k1)
      call path_point(m,x,seg,ltot-hl(2),q2,k2)
      if (ltot - sum(hl) > 1d-12) then
         allocate(y(2,k2-k1+2))
         ny = 1
         y(:,1) = q1
         do k = k1+1, k2
            ny = ny + 1
            y(:,ny) = x(:,k)
         end do
         ny = ny + 1
         y(:,ny) = q2
         call stroke_styled(ny,y,.false.,.not.hashead(1),.not.hashead(2))
      end if

      ! the heads: a triangle from the base of each to its tip
      do k = 1, 2
         if (.not.hashead(k)) cycle
         if (k == 1) then
            tip = x(:,1)
            base = q1
         else
            tip = x(:,m)
            base = q2
         end if
         u = tip - base
         if (norm2(u) < 1d-12) cycle
         u = (/-u(2),u(1)/) / norm2(u)
         call emit_tri(tip,base + hwh * u,base - hwh * u)
      end do

    end subroutine stroke_with_heads

    !> Stroke the path of m points x (closed if closed) in the outline
    !> style of the shape: solid (stroke_path, with round caps cap1/cap2
    !> at the ends of an open path), dashed (butt-ended pieces of the
    !> path) or dotted (round dots along it). The pattern is stretched to
    !> fit the path: a whole number of periods around a closed path, and
    !> a dash or a dot at both ends of an open one.
    subroutine stroke_styled(m,x,closed,cap1,cap2)
      integer, intent(in) :: m
      real*8, intent(in) :: x(2,m)
      logical, intent(in) :: closed, cap1, cap2

      real*8, allocatable :: xo(:,:), seg(:), y(:,:)
      real*8 :: ltot, on, per, s0, s1, q0(2), q1(2)
      integer :: mo, k, i, nper, k0, k1, ny

      if (sh%dash == planardash_solid .or. m < 2) then
         call stroke_path(m,x,closed,cap1,cap2)
         return
      end if

      ! the path as an open polyline (a closed one ends where it starts)
      mo = m
      if (closed) mo = m + 1
      allocate(xo(2,mo),seg(mo-1))
      xo(:,1:m) = x
      if (closed) xo(:,mo) = x(:,1)
      do k = 1, mo-1
         seg(k) = norm2(xo(:,k+1) - xo(:,k))
      end do
      ltot = sum(seg)
      if (ltot < 1d-12) return

      if (sh%dash == planardash_dotted) then
         per = dot_step * sh%width
         nper = max(nint(ltot / per),1)
         per = ltot / nper
         do i = 0, merge(nper-1,nper,closed)
            call path_point(mo,xo,seg,i*per,q0,k0)
            call emit_disk(q0,hw)
         end do
      else
         on = dash_on * sh%width
         per = on + dash_off * sh%width
         if (closed) then
            nper = max(nint(ltot / per),1)
            on = on * ltot / (nper * per)
            per = ltot / nper
            nper = nper - 1
         elseif (ltot <= on) then
            nper = 0
            on = ltot
         else
            nper = max(nint((ltot - on) / per),1)
            on = on * ltot / (nper * per + on)
            per = (ltot - on) / nper
         end if
         allocate(y(2,mo+2))
         do i = 0, nper
            s0 = i * per
            s1 = min(s0 + on,ltot)
            call path_point(mo,xo,seg,s0,q0,k0)
            call path_point(mo,xo,seg,s1,q1,k1)
            ny = 1
            y(:,1) = q0
            do k = k0+1, k1
               ny = ny + 1
               y(:,ny) = xo(:,k)
            end do
            ny = ny + 1
            y(:,ny) = q1
            call stroke_path(ny,y,.false.,.false.,.false.)
         end do
      end if

    end subroutine stroke_styled

    !> Stroke the path of m points x (closed if closed) with half-width
    !> hw. cap1/cap2 = round caps at the first/last point of an open path.
    subroutine stroke_path(m,x,closed,cap1,cap2)
      integer, intent(in) :: m
      real*8, intent(in) :: x(2,m)
      logical, intent(in) :: closed, cap1, cap2

      real*8, allocatable :: y(:,:), nrm(:,:), mv(:,:)
      logical, allocatable :: okm(:)
      real*8 :: d(2), dl, mm(2), lm, ch, la(2), ra(2), lb(2), rb(2)
      integer :: i, j, k, ny, nseg, kin, kout

      ! drop repeated points (and the closing point of a closed path)
      allocate(y(2,m))
      ny = 0
      do i = 1, m
         if (ny > 0) then
            if (norm2(x(:,i) - y(:,ny)) < 1d-10) cycle
         end if
         ny = ny + 1
         y(:,ny) = x(:,i)
      end do
      if (closed .and. ny > 1) then
         if (norm2(y(:,ny) - y(:,1)) < 1d-10) ny = ny - 1
      end if
      if (ny == 0) return
      if (ny == 1) then
         if (cap1 .or. cap2) call emit_disk(y(:,1),hw)
         return
      end if

      ! segment normals
      nseg = ny - 1
      if (closed .and. ny >= 3) then
         nseg = ny
      end if
      allocate(nrm(2,nseg))
      do k = 1, nseg
         j = modulo(k,ny) + 1
         d = y(:,j) - y(:,k)
         dl = norm2(d)
         nrm(:,k) = (/-d(2),d(1)/) / dl
      end do

      ! miter offsets at the vertices with two segments, where not too long
      allocate(mv(2,ny),okm(ny))
      okm = .false.
      mv = 0d0
      do i = 1, ny
         kout = i
         kin = i - 1
         if (kin == 0) kin = merge(nseg,0,nseg == ny)
         if (kout > nseg .or. kin == 0) cycle
         mm = nrm(:,kin) + nrm(:,kout)
         lm = norm2(mm)
         if (lm < 1d-10) cycle
         mm = mm / lm
         ch = dot_product(mm,nrm(:,kout))
         if (ch * miter_limit > 1d0) then
            mv(:,i) = mm * hw / ch
            okm(i) = .true.
         end if
      end do

      ! one quad per segment
      do k = 1, nseg
         j = modulo(k,ny) + 1
         if (okm(k)) then
            la = y(:,k) + mv(:,k)
            ra = y(:,k) - mv(:,k)
         else
            la = y(:,k) + hw * nrm(:,k)
            ra = y(:,k) - hw * nrm(:,k)
         end if
         if (okm(j)) then
            lb = y(:,j) + mv(:,j)
            rb = y(:,j) - mv(:,j)
         else
            lb = y(:,j) + hw * nrm(:,k)
            rb = y(:,j) - hw * nrm(:,k)
         end if
         call emit_tri(la,ra,rb)
         call emit_tri(la,rb,lb)
      end do

      ! round joins where there is no miter, and round caps
      do i = 1, ny
         if (okm(i)) cycle
         if (nseg < ny .and. (i == 1 .or. i == ny)) then
            ! an end of an open path
            if ((i == 1 .and. cap1) .or. (i == ny .and. cap2)) call emit_disk(y(:,i),hw)
         else
            call emit_disk(y(:,i),hw)
         end if
      end do

    end subroutine stroke_path

    !> Fill the (possibly concave) polygon of m points x by ear
    !> clipping. If no ear is found (a self-intersecting polygon), the
    !> rest is filled as a fan.
    subroutine fill_polygon(m,x)
      integer, intent(in) :: m
      real*8, intent(in) :: x(2,m)

      integer, allocatable :: idx(:)
      integer :: nr, i, ia, ib, ic, j
      real*8 :: area, sgn
      logical :: found, isear

      allocate(idx(m))
      do i = 1, m
         idx(i) = i
      end do
      nr = m

      ! orientation
      area = 0d0
      do i = 1, m
         j = modulo(i,m) + 1
         area = area + x(1,i) * x(2,j) - x(1,j) * x(2,i)
      end do
      sgn = sign(1d0,area)

      do while (nr > 3)
         found = .false.
         do i = 1, nr
            ia = idx(modulo(i-2,nr)+1)
            ib = idx(i)
            ic = idx(modulo(i,nr)+1)
            ! convex corner
            if (sgn * cross2(x(:,ib)-x(:,ia),x(:,ic)-x(:,ib)) <= 1d-14) cycle
            ! no other vertex inside
            isear = .true.
            do j = 1, nr
               if (idx(j) == ia .or. idx(j) == ib .or. idx(j) == ic) cycle
               if (in_triangle(x(:,idx(j)),x(:,ia),x(:,ib),x(:,ic))) then
                  isear = .false.
                  exit
               end if
            end do
            if (.not.isear) cycle
            call emit_tri(x(:,ia),x(:,ib),x(:,ic))
            idx(i:nr-1) = idx(i+1:nr)
            nr = nr - 1
            found = .true.
            exit
         end do
         if (.not.found) exit
      end do
      do i = 2, nr-1
         call emit_tri(x(:,idx(1)),x(:,idx(i)),x(:,idx(i+1)))
      end do

    end subroutine fill_polygon

    !> Emit a disk of radius r around x, as a triangle fan.
    subroutine emit_disk(x,r)
      use param, only: tpi
      real*8, intent(in) :: x(2), r

      integer :: i
      real*8 :: t1, t2

      do i = 1, nseg_disk
         t1 = tpi * real(i-1,8) / real(nseg_disk,8)
         t2 = tpi * real(i,8) / real(nseg_disk,8)
         call emit_tri(x,x + r * (/cos(t1),sin(t1)/),x + r * (/cos(t2),sin(t2)/))
      end do

    end subroutine emit_disk

    !> Emit the triangle (a,b,c), with the current color, at the depth of
    !> this shape, to the draw list of its layer.
    subroutine emit_tri(a,b,c)
      real*8, intent(in) :: a(2), b(2), c(2)

      if (sh%infront) then
         call emit_vertices(obj%flatfront,obj%nflatfront,a,b,c)
      else
         call emit_vertices(obj%flatback,obj%nflatback,a,b,c)
      end if

    end subroutine emit_tri

    !> Append the vertices of the triangle (a,b,c) to the vertex list lst
    !> with count nv (allocated by scene_objects_reset).
    subroutine emit_vertices(lst,nv,a,b,c)
      use shapes, only: flat_vert_nf
      use types, only: realloc
      real(c_float), allocatable, intent(inout) :: lst(:,:)
      integer, intent(inout) :: nv
      real*8, intent(in) :: a(2), b(2), c(2)

      integer :: k

      if (nv + 3 > size(lst,2)) call realloc(lst,int(flat_vert_nf),2*(nv+3))
      lst(1:2,nv+1) = real(a,c_float)
      lst(1:2,nv+2) = real(b,c_float)
      lst(1:2,nv+3) = real(c,c_float)
      do k = nv+1, nv+3
         lst(3,k) = real(z,c_float)
         lst(4:7,k) = rgba
      end do
      nv = nv + 3

    end subroutine emit_vertices

  end subroutine planar_tessellate

  !xx! private procedures

  !> The point q at arc length s along the path of m points x, with
  !> segment lengths seg; k is the segment it falls on (between x(:,k)
  !> and x(:,k+1)).
  subroutine path_point(m,x,seg,s,q,k)
    integer, intent(in) :: m
    real*8, intent(in) :: x(2,m), seg(m-1), s
    real*8, intent(out) :: q(2)
    integer, intent(out) :: k

    real*8 :: acc, t

    acc = 0d0
    do k = 1, m-1
       if (acc + seg(k) >= s .or. k == m-1) then
          t = 0d0
          if (seg(k) > 1d-14) t = min(max((s - acc) / seg(k),0d0),1d0)
          q = x(:,k) + t * (x(:,k+1) - x(:,k))
          return
       end if
       acc = acc + seg(k)
    end do

  end subroutine path_point

  !> Distance from point p to the segment ab.
  function segdist(p,a,b) result(d)
    real*8, intent(in) :: p(2), a(2), b(2)
    real*8 :: d

    real*8 :: ab(2), t, l2

    ab = b - a
    l2 = dot_product(ab,ab)
    if (l2 < 1d-24) then
       d = norm2(p - a)
    else
       t = min(max(dot_product(p - a,ab) / l2,0d0),1d0)
       d = norm2(p - a - t * ab)
    end if

  end function segdist

  !> z component of the cross product of the 2D vectors a and b.
  function cross2(a,b) result(c)
    real*8, intent(in) :: a(2), b(2)
    real*8 :: c

    c = a(1) * b(2) - a(2) * b(1)

  end function cross2

  !> Whether the point p is inside (or on the edge of) the triangle abc,
  !> of either orientation.
  function in_triangle(p,a,b,c) result(ok)
    real*8, intent(in) :: p(2), a(2), b(2), c(2)
    logical :: ok

    real*8 :: d1, d2, d3

    d1 = cross2(b - a,p - a)
    d2 = cross2(c - b,p - b)
    d3 = cross2(a - c,p - c)
    ok = .not.((d1 < 0d0 .or. d2 < 0d0 .or. d3 < 0d0) .and. (d1 > 0d0 .or. d2 > 0d0 .or. d3 > 0d0))

  end function in_triangle

  !> Whether the point p is inside the polygon of n points x (even-odd rule).
  function inside_polygon(p,n,x) result(ok)
    real*8, intent(in) :: p(2)
    integer, intent(in) :: n
    real*8, intent(in) :: x(2,n)
    logical :: ok

    integer :: i, j

    ok = .false.
    j = n
    do i = 1, n
       if ((x(2,i) > p(2)) .neqv. (x(2,j) > p(2))) then
          if (p(1) < (x(1,j) - x(1,i)) * (p(2) - x(2,i)) / (x(2,j) - x(2,i)) + x(1,i)) &
             ok = .not.ok
       end if
       j = i
    end do

  end function inside_polygon

end submodule planar
