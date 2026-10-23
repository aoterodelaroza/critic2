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

! Maps of a field on a plane (the PLANE keyword): the definition shared
! by the plane map window and the plane map objects, and the drawing of
! the objects (reptype_slice).
submodule (representations) plane
  implicit none

  ! the number of colors of the map
  integer, parameter :: nlut = 256
  integer, parameter :: nkeyreg = 13 ! entries of plane_def_key for the region and nx, ny; the quantity follows

contains

  !> The chemical function of quantity nti (its index in arith_names),
  !> or zero if nti is a property of the field. Only the functions that
  !> take just the field (elf(1)) are quantities.
  module function quant_chem(nti) result(j)
    use arithmetic, only: arith_names, arithkind_chem
    integer, intent(in) :: nti
    integer :: j

    j = nti - nquant
    if (j < 1 .or. j > size(arith_names,1)) then
       j = 0
    elseif (arith_names(j)%kind /= arithkind_chem .or. arith_names(j)%args /= "(id)") then
       j = 0
    end if

  end function quant_chem

  !> How the quantity of a form is evaluated by the rhoplot kernels:
  !> the arithmetic expression expr (if isexpr), or quantity nti of
  !> field ifield. Returns the field id and its property ntik, or (id
  !> = -1) the expression qexpr: expr, or that of the chemical
  !> function.
  module subroutine quant_eval(isexpr,expr,ifield,nti,id,qexpr,ntik)
    logical, intent(in) :: isexpr
    character(len=*), intent(in) :: expr
    integer, intent(in) :: ifield, nti
    integer, intent(out) :: id, ntik
    character(len=:), allocatable, intent(out) :: qexpr

    id = -1
    ntik = 0
    if (isexpr) then
       qexpr = trim(expr)
    elseif (quant_chem(nti) > 0) then
       qexpr = quant_label(.false.,"",ifield,nti)
    else
       id = ifield
       ntik = max(min(nti,nquant-1),0)
       qexpr = ""
    end if

  end subroutine quant_eval

  !> The label of the quantity of a form: the arithmetic expression
  !> expr (if isexpr), or quantity nti of field ifield ($1:lap, or
  !> elf(1)).
  module function quant_label(isexpr,expr,ifield,nti) result(str)
    use arithmetic, only: arith_names
    use tools_io, only: string
    logical, intent(in) :: isexpr
    character(len=*), intent(in) :: expr
    integer, intent(in) :: ifield, nti
    character(len=:), allocatable :: str

    integer :: j

    j = quant_chem(nti)
    if (isexpr) then
       str = trim(expr)
    elseif (j > 0) then
       str = trim(arith_names(j)%name) // "(" // string(ifield) // ")"
    else
       str = "$" // string(ifield) // trim(quant_short(max(min(nti,nquant-1),0)))
    end if

  end function quant_label

  !> The map of a new system isys. The quantity is the value of the
  !> reference field. The region: for crystals, the parallelogram of
  !> a cell face; for molecules, the plane of the first three atoms (if
  !> they are not on a line) on the rectangle that covers the whole
  !> molecule with some margin, or a face of the molecular box. The
  !> vacuum directions of a crystal (axes with a long empty gap) are
  !> treated as in molecules: a slab gets the face of its two periodic
  !> axes through the first atom; a chain, its periodic lattice vector
  !> and the perpendicular direction along which the atoms spread the
  !> most, covering them with the margin; a molecule in a box, the
  !> plane of the molecule.
  module subroutine plane_def_set_defaults(d,isys)
    use systems, only: sys
    use crystalmod, only: iperiod_vacthr
    use tools_math, only: cross
    class(plane_def), intent(inout) :: d
    integer, intent(in) :: isys

    real*8, parameter :: mmol = 2d0 / bohrtoa ! margin around the molecule (bohr)
    real*8 :: c(3), e1(3), e2(3), f1(3), f2(3), umin, umax, vmin, vmax, p(3), xf(3)
    real*8 :: sxx, syy, sxy, theta
    real*8, allocatable :: xat(:,:), q(:,:)
    integer :: i, nvac, jper(3), ncel
    logical :: done, ismol, isvac(3)

    d%isexpr = .false.
    d%nti = 0
    call d%fref%set(sys(isys),sys(isys)%iref)
    ismol = sys(isys)%c%ismolecule
    ncel = sys(isys)%c%ncel

    ! the vacuum axes, and the periodic ones
    isvac = .false.
    if (.not.ismol) isvac = (sys(isys)%c%vaclength > iperiod_vacthr)
    nvac = count(isvac)
    jper = pack((/1,2,3/),.not.isvac,(/0,0,0/))

    ! the atoms (cell-frame Cartesian), on the same side of the vacuum
    ! in the vacuum axes
    allocate(xat(3,ncel))
    do i = 1, ncel
       xf = sys(isys)%c%atcel(i)%x
       xf = xf - floor(xf)
       where (isvac .and. xf >= sys(isys)%c%vactop - 1d-10) xf = xf - 1d0
       xat(:,i) = sys(isys)%c%x2c(xf)
    end do

    ! the region: the parallelogram of a face for crystals and slabs,
    ! a rectangle for molecules and chains
    d%region = merge(0,1,ismol .or. nvac >= 2)

    done = .false.
    if (nvac == 1) then
       ! slab: the face of the periodic axes through the first atom
       call face(jper(1),jper(2))
       done = .true.
    elseif (nvac == 2 .and. ncel >= 1) then
       ! chain: the periodic lattice vector (from 0 to 1), and the
       ! perpendicular direction along which the atoms spread the most
       ! (the principal axis of their perpendicular positions; a vacuum
       ! lattice vector if all are on the chain), covering them
       e1 = sys(isys)%c%m_x2c(:,jper(1)) / sys(isys)%c%aa(jper(1))
       f1 = sys(isys)%c%m_x2c(:,findloc(isvac,.true.,dim=1))
       f1 = f1 - dot_product(f1,e1) * e1
       f1 = f1 / norm2(f1)
       f2 = cross(e1,f1)
       allocate(q(2,ncel))
       q(1,:) = matmul(f1,xat)
       q(2,:) = matmul(f2,xat)
       q(1,:) = q(1,:) - sum(q(1,:)) / ncel
       q(2,:) = q(2,:) - sum(q(2,:)) / ncel
       sxx = dot_product(q(1,:),q(1,:))
       syy = dot_product(q(2,:),q(2,:))
       sxy = dot_product(q(1,:),q(2,:))
       ! (any direction if the spread is isotropic, or there is none)
       if (abs(sxy) + abs(sxx - syy) > 1d-8) then
          theta = 0.5d0 * atan2(2d0 * sxy,sxx - syy)
          e2 = cos(theta) * f1 + sin(theta) * f2
       else
          e2 = f1
       end if
       call extent(e2,vmin,vmax)
       p = xat(:,1) + (vmin - mmol) * e2
       ! the side along the chain starts at the origin of the cell
       xf = sys(isys)%c%c2x(p)
       xf(jper(1)) = 0d0
       p = sys(isys)%c%x2c(xf)
       d%xp(:,1) = form(p)
       d%xp(:,2) = form(p + sys(isys)%c%m_x2c(:,jper(1)))
       d%xp(:,3) = form(p + (vmax - vmin + 2*mmol) * e2)
       done = .true.
    elseif ((ismol .or. nvac == 3) .and. ncel >= 3) then
       ! molecule: the plane of the first three atoms, covering all of them
       e1 = xat(:,2) - xat(:,1)
       if (norm2(e1) > 1d-4) then
          e1 = e1 / norm2(e1)
          c = xat(:,3) - xat(:,1)
          c = c - dot_product(c,e1) * e1
          if (norm2(c) > 1d-4) then
             e2 = c / norm2(c)
             call extent(e1,umin,umax)
             call extent(e2,vmin,vmax)
             p = xat(:,1) + (umin - mmol) * e1 + (vmin - mmol) * e2
             d%xp(:,1) = form(p)
             d%xp(:,2) = form(p + (umax - umin + 2*mmol) * e1)
             d%xp(:,3) = form(p + (vmax - vmin + 2*mmol) * e2)
             done = .true.
          end if
       end if
    end if

    if (.not.done) then
       if (ismol) then
          d%xp(:,1) = form(sys(isys)%c%x2c((/0d0,0d0,0.5d0/)))
          d%xp(:,2) = form(sys(isys)%c%x2c((/1d0,0d0,0.5d0/)))
          d%xp(:,3) = form(sys(isys)%c%x2c((/0d0,1d0,0.5d0/)))
       else
          call face(1,2)
       end if
    end if
    ! no negative zeros in the forms
    where (abs(d%xp) < 1d-12) d%xp = 0d0

  contains
    !> The point xc (cell-frame Cartesian bohr) in the coordinates of
    !> the forms.
    function form(xc) result(x)
      real*8, intent(in) :: xc(3)
      real*8 :: x(3)
      if (ismol) then
         x = (xc + sys(isys)%c%molx0) * bohrtoa
      else
         x = sys(isys)%c%c2x(xc)
      end if
    end function form

    !> The face of the cell spanned by axes j1 and j2 (crystals), at
    !> the coordinates of the first atom along the vacuum axes.
    subroutine face(j1,j2)
      integer, intent(in) :: j1, j2
      d%xp(:,1) = 0d0
      if (ncel > 0) then
         xf = sys(isys)%c%c2x(xat(:,1))
         d%xp(:,1) = merge(xf,0d0,isvac)
         d%xp(j1,1) = 0d0
         d%xp(j2,1) = 0d0
      end if
      d%xp(:,2) = d%xp(:,1)
      d%xp(j1,2) = d%xp(j1,2) + 1d0
      d%xp(:,3) = d%xp(:,1)
      d%xp(j2,3) = d%xp(j2,3) + 1d0
    end subroutine face

    !> The extent of the atoms along the unit vector e, from the first
    !> atom.
    subroutine extent(e,emin,emax)
      real*8, intent(in) :: e(3)
      real*8, intent(out) :: emin, emax
      real*8 :: pe(ncel)
      pe = matmul(e,xat) - dot_product(e,xat(:,1))
      emin = minval(pe)
      emax = maxval(pe)
    end subroutine extent
  end subroutine plane_def_set_defaults

  !> The grid of the region of the map in system isys (cell-frame
  !> Cartesian bohr): the origin x0 and the steps uu, vv (the points are
  !> x0 + (ix-1)*uu + (iy-1)*vv), the first point a0, and the
  !> orthonormal frame of the plane (e1 along the first side, e2
  !> perpendicular to it toward the third point). Returns a non-empty
  !> errmsg if the points do not define a region.
  module subroutine plane_def_grid(d,isys,x0,uu,vv,a0,e1,e2,errmsg)
    use systems, only: sys
    use tools_math, only: plane_scale_extend
    class(plane_def), intent(in) :: d
    integer, intent(in) :: isys
    real*8, intent(out) :: x0(3), uu(3), vv(3), a0(3), e1(3), e2(3)
    character(len=:), allocatable, intent(out) :: errmsg

    real*8 :: x1(3), x2(3), d1(3), c(3), cperp(3), l1, h, m

    errmsg = ""
    x0 = 0d0
    uu = 0d0
    vv = 0d0
    e1 = 0d0
    e2 = 0d0
    a0 = cart(d%xp(:,1))
    x1 = cart(d%xp(:,2))
    x2 = cart(d%xp(:,3))
    d1 = x1 - a0
    l1 = norm2(d1)
    c = x2 - a0
    if (l1 < 1d-8) then
       errmsg = "The origin and the X end are the same point"
       return
    end if
    e1 = d1 / l1
    cperp = c - dot_product(c,e1) * e1
    h = norm2(cperp)
    if (h < 1d-8) then
       errmsg = "The three points are on a line: they do not define a plane"
       return
    end if
    e2 = cperp / h

    ! the second side perpendicular to the first, or to the third
    ! point; with the margin on every side
    if (d%region == 0) x2 = a0 + cperp
    m = d%margin / bohrtoa
    x0 = a0
    call plane_scale_extend(x0,x1,x2,1d0,1d0,m,m,m,m)
    ! a negative margin of half a side or more would turn the region over
    if (dot_product(x1 - x0,e1) <= 1d-8 .or. dot_product(x2 - x0,e2) <= 1d-8) then
       errmsg = "The negative margin is larger than the region"
       return
    end if
    uu = (x1 - x0) / real(d%nx-1,8)
    vv = (x2 - x0) / real(d%ny-1,8)

  contains
    !> The point x of the forms in cell-frame Cartesian bohr.
    function cart(x) result(xc)
      real*8, intent(in) :: x(3)
      real*8 :: xc(3)
      if (sys(isys)%c%ismolecule) then
         xc = x / bohrtoa - sys(isys)%c%molx0
      else
         xc = sys(isys)%c%x2c(x)
      end if
    end function cart
  end subroutine plane_def_grid

  !> If autonpts, the number of points along each side of the region
  !> of the map in system isys (nx, ny): the cells of the grid are
  !> (nearly) square, with sides of plane_auto_spacing, or less if the
  !> longer side would have fewer than plane_auto_nmin points, or more
  !> if there would be more than plane_auto_ntot points in all, more
  !> than plane_npts_max along a side, or the map would take longer
  !> than plane_auto_tmax (the cost of the quantity is probed on a few
  !> points when it changes).
  module subroutine plane_def_autopoints(d,isys)
    use systems, only: sysc, field_cost
    use rhoplot, only: rhoplot_nti_request
    use param, only: bohrtoa
    class(plane_def), intent(inout) :: d
    integer, intent(in) :: isys

    integer, parameter :: nprobe = 9 ! points along each side of the probe of the cost

    real*8 :: x0(3), uu(3), vv(3), a0(3), e1(3), e2(3), l1, l2, h, t, np, nmax
    real*8, allocatable :: ff(:,:)
    integer*8, allocatable :: key(:)
    integer*8 :: c0, c1, crate
    integer :: idq, ntiq
    logical :: stale
    character(len=:), allocatable :: errmsg, qexpr

    if (.not.d%autonpts) return
    call d%grid(isys,x0,uu,vv,a0,e1,e2,errmsg)
    if (len_trim(errmsg) > 0) return

    ! the sides (bohr), and the spacing
    l1 = norm2(uu) * (d%nx - 1)
    l2 = norm2(vv) * (d%ny - 1)
    if (max(l1,l2) <= 0d0) return
    h = min(plane_auto_spacing / bohrtoa,max(l1,l2) / (plane_auto_nmin - 1))

    ! the cost of the quantity (wall time per point): the cost the
    ! field has measured for a property of a field (not measured during
    ! a dynamics run, which resets the fields on every frame: the value
    ! in hand is kept until there is a new one); for an expression,
    ! measured on a coarse grid of the region when the quantity changes
    call quant_eval(d%isexpr,d%expr,d%fref%id,d%nti,idq,qexpr,ntiq)
    if (idq >= 0) then
       t = field_cost(isys,idq,.not.sysc(isys)%md_run,rhoplot_nti_request(ntiq))
       if (t > 0d0) d%tpoint = t
       ! the probe in hand, if any, was for another quantity
       if (allocated(d%tpointkey)) deallocate(d%tpointkey)
       stale = .false.
    else
       call d%key(isys,key)
       stale = .not.allocated(d%tpointkey) .or. .not.allocated(d%tpointexpr)
       if (.not.stale) stale = (size(d%tpointkey) /= size(key) - nkeyreg)
       if (.not.stale) stale = any(d%tpointkey /= key(nkeyreg+1:)) .or. (d%tpointexpr /= trim(d%expr))
    end if
    if (stale) then
       d%tpointkey = key(nkeyreg+1:)
       d%tpointexpr = trim(d%expr)
       d%tpoint = 0d0
       call system_clock(c0,crate)
       call d%values(isys,x0,uu * (d%nx - 1) / (nprobe - 1),vv * (d%ny - 1) / (nprobe - 1),ff,errmsg,&
          npts=nprobe)
       call system_clock(c1)
       if (allocated(ff)) d%tpoint = real(c1 - c0,8) / real(crate,8) / nprobe**2
    end if

    ! at most plane_auto_ntot points, and fewer if the map would take
    ! longer than plane_auto_tmax
    nmax = plane_auto_ntot
    if (d%tpoint > 0d0) nmax = min(nmax,plane_auto_tmax / d%tpoint)
    np = (l1 / h + 1d0) * (l2 / h + 1d0)
    if (np > nmax) h = h * sqrt(np / nmax)
    h = max(h,max(l1,l2) / (plane_npts_max - 1))

    ! the sides come from the grid of the current nx and ny, with a
    ! rounding error that depends on them: keep them unless they are
    ! more than half a point away, or a side near k + 1/2 points could
    ! alternate between two values (and the map be calculated again in
    ! every frame)
    t = l1 / h + 1d0
    if (abs(t - d%nx) > 0.5d0 + 1d-6) d%nx = int(max(min(nint(t),plane_npts_max),plane_npts_min),c_int)
    t = l2 / h + 1d0
    if (abs(t - d%ny) > 0.5d0 + 1d-6) d%ny = int(max(min(nint(t),plane_npts_max),plane_npts_min),c_int)

  end subroutine plane_def_autopoints

  !> Everything the values of the map in system isys depend on: the
  !> region and the number of points (the first nkeyreg entries), then
  !> the quantity, the geometry, and the fields (all of them, for an
  !> expression; which enters through the size of key). The
  !> expression itself is not in it.
  module subroutine plane_def_key(d,isys,key)
    use systems, only: sys, sysc
    class(plane_def), intent(in) :: d
    integer, intent(in) :: isys
    integer*8, allocatable, intent(inout) :: key(:)

    integer :: j, nf

    nf = 0
    if (d%isexpr) nf = sys(isys)%nf + 1
    if (allocated(key)) deallocate(key)
    allocate(key(19+nf))
    key(1:9) = transfer(d%xp,key(1:9))
    key(10) = d%region
    key(11) = transfer(d%margin,key(11))
    key(12) = d%nx
    key(13) = d%ny
    key(14) = merge(1_8,0_8,d%isexpr)
    key(15) = d%fref%uid
    key(16) = d%fref%id
    key(17) = d%nti
    key(18) = sys(isys)%field_cuid(d%fref%id)
    key(19) = transfer(sysc(isys)%timelastchange_geometry,key(19))
    do j = 0, nf-1
       key(20+j) = sys(isys)%field_cuid(j)
    end do

  end subroutine plane_def_key

  !> Calculate the values of the map in system isys on the grid x0,
  !> uu, vv (see plane_def_grid) in ff (nx,ny). The expression is
  !> checked again (the fields it uses may have come or gone). If they
  !> cannot be calculated, ff is not allocated and errmsg says why
  !> (empty for an empty or wrong expression: exprerr says it). If
  !> npts is present, the grid has npts points along each side
  !> instead of nx and ny (the probe of the cost).
  module subroutine plane_def_values(d,isys,x0,uu,vv,ff,errmsg,npts)
    use systems, only: sys
    use rhoplot, only: rhoplot_plane_values
    class(plane_def), intent(inout) :: d
    integer, intent(in) :: isys
    real*8, intent(in) :: x0(3), uu(3), vv(3)
    real*8, allocatable, intent(inout) :: ff(:,:)
    character(len=:), allocatable, intent(inout) :: errmsg
    integer, intent(in), optional :: npts

    integer :: id, nx, ny, idq, ntiq
    logical :: ok
    character(len=:), allocatable :: lerr, qexpr

    errmsg = ""
    if (allocated(ff)) deallocate(ff)
    nx = d%nx
    ny = d%ny
    if (present(npts)) then
       nx = npts
       ny = npts
    end if
    if (d%isexpr) then
       id = -1
       d%exprerr = ""
       if (len_trim(d%expr) > 0) then
          call sys(isys)%check_expression(trim(d%expr),lerr)
          d%exprerr = lerr
       end if
       ok = (len_trim(d%expr) > 0 .and. len_trim(d%exprerr) == 0)
    else
       id = d%fref%id
       ok = sys(isys)%goodfield(id)
       if (.not.ok) errmsg = "The field is not available in this system"
    end if
    if (.not.ok) return

    allocate(ff(nx,ny))
    call quant_eval(d%isexpr,d%expr,id,d%nti,idq,qexpr,ntiq)
    call rhoplot_plane_values(sys(isys),x0,uu,vv,nx,ny,idq,qexpr,ntiq,ff,errmsg)
    if (len_trim(errmsg) > 0) deallocate(ff)

  end subroutine plane_def_values

  !> The contour levels of the map with values ff: niso of them, in
  !> ziso (none for an empty list).
  module subroutine plane_def_levels(d,ff,niso,ziso)
    use tools_math, only: assign_ziso, niso_manual, niso_bader
    use tools_io, only: isreal
    use types, only: realloc
    class(plane_def), intent(in) :: d
    real*8, intent(in) :: ff(:,:)
    integer, intent(out) :: niso
    real*8, allocatable, intent(inout) :: ziso(:)

    integer :: i, lp
    real*8 :: fmin, fmax, rval
    character(len=:), allocatable :: llist

    niso = 0
    if (d%ctype == niso_manual) then
       ! a list of values, separated by spaces or commas
       if (allocated(ziso)) deallocate(ziso)
       allocate(ziso(10))
       llist = d%llist
       do i = 1, len(llist)
          if (llist(i:i) == ",") llist(i:i) = " "
       end do
       lp = 1
       do while (isreal(rval,llist,lp))
          niso = niso + 1
          if (niso > size(ziso)) call realloc(ziso,2*niso)
          ziso(niso) = rval
       end do
       call realloc(ziso,niso)
    else
       if (d%lrange .and. d%ctype /= niso_bader) then
          fmin = d%lmin
          fmax = d%lmax
       else
          fmin = minval(ff)
          fmax = maxval(ff)
       end if
       ! a constant map has no contours
       if (fmax <= fmin) return
       niso = max(d%nlev,2)
       call assign_ziso(d%ctype,niso,ziso,fmin,fmax)
    end if

  end subroutine plane_def_levels

  !> The contour lines of the map with values ff on the grid with
  !> steps uu, vv: the niso levels (ziso) and the nline lines (see
  !> rhoplot_contour_lines), in the frame of the region.
  module subroutine plane_def_contours(d,ff,uu,vv,niso,ziso,nline,lptr,xy,liso)
    use rhoplot, only: rhoplot_contour_lines
    class(plane_def), intent(in) :: d
    real*8, intent(in) :: ff(:,:), uu(3), vv(3)
    integer, intent(out) :: niso, nline
    real*8, allocatable, intent(inout) :: ziso(:), xy(:,:)
    integer, allocatable, intent(inout) :: lptr(:), liso(:)

    integer :: nx, ny
    real*8 :: r01, r02, cosa

    nline = 0
    call d%levels(ff,niso,ziso)
    if (niso == 0) return
    nx = size(ff,1)
    ny = size(ff,2)
    r01 = norm2(uu) * (nx-1)
    r02 = norm2(vv) * (ny-1)
    cosa = dot_product(uu,vv) / (norm2(uu) * norm2(vv))
    call rhoplot_contour_lines(ff,nx,ny,r01,r02,cosa,niso,ziso,nline,lptr,xy,liso)

  end subroutine plane_def_contours

  !> What the contour lines of the map depend on, other than the values
  !> and the list of levels (llist): whether they are drawn, and the
  !> levels.
  module subroutine plane_def_contkey(d,key)
    class(plane_def), intent(in) :: d
    integer*8, allocatable, intent(inout) :: key(:)

    if (allocated(key)) deallocate(key)
    allocate(key(6))
    key(1) = merge(1_8,0_8,d%showcont)
    key(2) = d%ctype
    key(3) = d%nlev
    key(4) = merge(1_8,0_8,d%lrange)
    key(5) = transfer(d%lmin,key(5))
    key(6) = transfer(d%lmax,key(6))

  end subroutine plane_def_contkey

  !> The colors of the values ff of the map, in rgb (3,nx*ny; the
  !> first index of ff runs fastest), and the range of the values on the
  !> scale of the colors (gmin, gmax).
  module subroutine plane_def_colors(d,ff,gmin,gmax,rgb)
    use ieee_arithmetic, only: ieee_is_finite
    class(plane_def), intent(in) :: d
    real*8, intent(in) :: ff(:,:)
    real*8, intent(out) :: gmin, gmax
    real(c_float), allocatable, intent(inout) :: rgb(:,:)

    integer :: nx, ny, ix, iy, k, il
    real*8 :: gv
    real*8, allocatable :: g(:,:)
    real(c_float) :: lut(3,nlut)

    ! the values on the scale of the colors, and their range
    nx = size(ff,1)
    ny = size(ff,2)
    allocate(g(nx,ny))
    gmin = huge(1d0)
    gmax = -huge(1d0)
    do iy = 1, ny
       do ix = 1, nx
          g(ix,iy) = d%cvalue(ff(ix,iy))
          if (.not.ieee_is_finite(g(ix,iy))) cycle
          gmin = min(gmin,g(ix,iy))
          gmax = max(gmax,g(ix,iy))
       end do
    end do
    if (d%crange) then
       gmin = d%cvalue(d%cmin)
       gmax = d%cvalue(d%cmax)
    end if
    if (gmin > gmax) then
       gmin = 0d0
       gmax = 1d0
    end if
    if (gmax <= gmin) gmax = gmin + max(abs(gmin),1d0) * 1d-6

    ! the colors
    call iw_colormap_lut(d%icmap,lut)
    if (allocated(rgb)) deallocate(rgb)
    allocate(rgb(3,nx*ny))
    k = 0
    do iy = 1, ny
       do ix = 1, nx
          k = k + 1
          gv = (g(ix,iy) - gmin) / (gmax - gmin)
          if (.not.ieee_is_finite(gv)) gv = 0d0
          il = nint(max(min(gv,1d0),0d0) * (nlut-1)) + 1
          rgb(:,k) = lut(:,il)
       end do
    end do

  end subroutine plane_def_colors

  !> What the colors of the map depend on, other than the values.
  module subroutine plane_def_colorkey(d,key)
    class(plane_def), intent(in) :: d
    integer*8, allocatable, intent(inout) :: key(:)

    if (allocated(key)) deallocate(key)
    allocate(key(5))
    key(1) = d%icmap
    key(2) = d%cscale
    key(3) = merge(1_8,0_8,d%crange)
    key(4) = transfer(d%cmin,key(4))
    key(5) = transfer(d%cmax,key(5))

  end subroutine plane_def_colorkey

  !> The value f on the scale of the colors of the map.
  module function plane_def_cvalue(d,f)
    class(plane_def), intent(in) :: d
    real*8, intent(in) :: f
    real*8 :: plane_def_cvalue

    plane_def_cvalue = scale_value(f,d%cscale)

  end function plane_def_cvalue

  !> Add the plane map object r to the draw lists obj: the colored
  !> mesh of the region and the contour lines (flat cylinders on both
  !> sides of the map). The mesh is flat, or a relief of the values
  !> (r%slice%rl). The values, the mesh, its
  !> colors, its shape, and the lines are cached, and calculated again
  !> only if what they depend on changed. The values are not calculated
  !> while the initialization threads run, or if they are given
  !> (nocalc).
  module subroutine slice_draw(r,obj)
    use systems, only: sys, are_threads_running
    use tools_math, only: cross
    use ieee_arithmetic, only: ieee_is_finite
    class(representation), intent(inout) :: r
    type(scene_objects), intent(inout) :: obj

    integer*8, allocatable :: key(:)
    integer :: isys, nx, ny, ix, iy, k, i, j, niso, is, nab
    logical :: stale
    real*8 :: a0(3), nrm(3), a(2,2), det, ainv(2,2), gmin, gmax, dd(2)
    real*8 :: ab(2,3), pp(3,3), nn(3,3), off
    type(dl_cylinder) :: dcyl

    isys = r%id
    if (sys(isys)%c%ncel == 0) return
    if (.not.allocated(r%slice%c)) allocate(r%slice%c)

    associate(d => r%slice%d, c => r%slice%c)
      ! the values (and with them, the mesh and the lines go stale)
      if (.not.r%slice%nocalc .and. .not.are_threads_running()) call d%autopoints(isys)
      call d%key(isys,key)
      stale = .not.allocated(c%key) .or. .not.allocated(c%keyexpr)
      if (.not.stale) stale = (size(c%key) /= size(key))
      if (.not.stale) stale = any(c%key /= key) .or. (c%keyexpr /= trim(d%expr))
      ! the values given by a window are never calculated here, nor any
      ! while the initialization threads run (the old ones are drawn)
      if (r%slice%nocalc .or. are_threads_running()) stale = .false.
      if (stale) then
         call move_alloc(key,c%key)
         c%keyexpr = trim(d%expr)
         if (allocated(c%ff)) deallocate(c%ff)
         if (allocated(c%keycont)) deallocate(c%keycont)
         if (allocated(c%keycolor)) deallocate(c%keycolor)
         c%mesh = dl_mesh()
         call d%grid(isys,c%x0,c%uu,c%vv,a0,c%e1,c%e2,c%errmsg)
         if (len_trim(c%errmsg) == 0) call d%values(isys,c%x0,c%uu,c%vv,c%ff,c%errmsg)
      end if
      if (.not.allocated(c%ff)) return
      nx = size(c%ff,1)
      ny = size(c%ff,2)

      ! the normal of the map, and the transformation from the frame
      ! of the region to the coordinates of the grid (0 at its origin)
      nrm = cross(c%e1,c%e2)
      a(1,1) = dot_product(c%uu,c%e1)
      a(1,2) = dot_product(c%vv,c%e1)
      a(2,1) = dot_product(c%uu,c%e2)
      a(2,2) = dot_product(c%vv,c%e2)
      det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
      ainv = reshape((/a(2,2), -a(2,1), -a(1,2), a(1,1)/),(/2,2/)) / det

      ! the mesh: one vertex per point (the first index runs fastest),
      ! two triangles per cell of the grid; the positions and the
      ! normals are its shape (below)
      if (c%mesh%nv == 0) then
         c%mesh%nv = nx * ny
         c%mesh%nf = 2 * (nx-1) * (ny-1)
         allocate(c%mesh%x(3,c%mesh%nv),c%mesh%nrm(3,c%mesh%nv),c%mesh%idx(3,c%mesh%nf))
         k = 0
         do iy = 1, ny-1
            do ix = 1, nx-1
               i = (iy-1) * nx + ix - 1 ! 0-based index of the vertex (ix,iy)
               c%mesh%idx(:,k+1) = int((/i, i+1, i+nx+1/),c_int)
               c%mesh%idx(:,k+2) = int((/i, i+nx+1, i+nx/),c_int)
               k = k + 2
            end do
         end do
         if (allocated(c%keyshape)) deallocate(c%keyshape)
      end if

      ! its colors
      call d%colorkey(key)
      stale = .not.allocated(c%keycolor) .or. .not.allocated(c%mesh%rgbv)
      if (.not.stale) stale = any(c%keycolor /= key)
      if (stale) then
         call move_alloc(key,c%keycolor)
         call d%colors(c%ff,gmin,gmax,c%mesh%rgbv)
      end if
      c%mesh%alpha = r%slice%alpha

      ! its shape: flat, or the relief
      associate(rl => r%slice%rl)
        if (rl%on .and. rl%height > 0d0) then
           key = [1_8,transfer(rl%height,1_8),int(rl%scale,8),merge(1_8,0_8,rl%frange),&
              transfer(rl%fmin,1_8),transfer(rl%fmax,1_8),merge(1_8,0_8,rl%flip)]
        else
           key = [0_8]
        end if
        stale = .not.allocated(c%keyshape)
        if (.not.stale) stale = (size(c%keyshape) /= size(key))
        if (.not.stale) stale = any(c%keyshape /= key)
        if (stale) then
           call move_alloc(key,c%keyshape)
           call slice_shape(c,rl,c%keyshape(1) == 1_8)
        end if
      end associate
      call dl_append(obj%msh,obj%nmsh,c%mesh)

      ! the contour lines
      if (.not.d%showcont) return
      call d%contkey(key)
      stale = .not.allocated(c%keycont) .or. .not.allocated(c%keylist)
      if (.not.stale) stale = any(c%keycont /= key) .or. (c%keylist /= trim(d%llist))
      if (stale) then
         call move_alloc(key,c%keycont)
         c%keylist = trim(d%llist)
         call d%contours(c%ff,c%uu,c%vv,niso,c%ziso,c%nline,c%lptr,c%xy,c%liso)
      end if

      ! one flat cylinder per segment, on each side of the map (a flat
      ! cylinder is drawn with half its radius), off the mesh along its
      ! normal; the coordinates of the lines are in the frame of the
      ! region, from the origin of the grid. On a relief, a segment that
      ! crosses the diagonal of its cell is split there, so that it
      ! follows the two triangles of the cell.
      dcyl%r = real(r%slice%crad,c_float)
      do i = 1, c%nline
         if (c%ziso(c%liso(i)) > 0d0) then
            dcyl%rgb = d%rgbpos
         else
            dcyl%rgb = d%rgbneg
         end if
         do j = c%lptr(i), c%lptr(i+1)-2
            ab(:,1) = matmul(ainv,c%xy(:,j))
            ab(:,2) = matmul(ainv,c%xy(:,j+1))
            ix = min(max(floor(0.5d0 * (ab(1,1) + ab(1,2))),0),nx-2)
            iy = min(max(floor(0.5d0 * (ab(2,1) + ab(2,2))),0),ny-2)
            dd = (ab(1,1:2) - ix) - (ab(2,1:2) - iy)
            nab = 2
            if (c%isrelief .and. dd(1) * dd(2) < 0d0) then
               nab = 3
               ab(:,3) = ab(:,2)
               ab(:,2) = ab(:,1) + dd(1) / (dd(1) - dd(2)) * (ab(:,3) - ab(:,1))
            end if
            do k = 1, nab
               call surface_point(c,ab(:,k),ix,iy,pp(:,k),nn(:,k))
            end do
            do is = -1, 1, 2
               off = real(is,8) * 0.5d0 * r%slice%crad
               do k = 1, nab-1
                  dcyl%x1 = real(pp(:,k) + off * nn(:,k),c_float)
                  dcyl%x2 = real(pp(:,k+1) + off * nn(:,k+1),c_float)
                  call dl_append(obj%cylflat,obj%ncylflat,dcyl)
               end do
            end do
         end do
      end do
    end associate

  contains

    !> The shape of the mesh of the map c: the positions and normals of
    !> its vertices. Flat, or (relief) the relief rl: each point rises
    !> along the normal with its value on the scale of the relief, from
    !> zero at the bottom of its range (fixed, in either order, or that
    !> of the values on the map) to the height at the top (clipped to
    !> it; zero where the value is not finite). Against the normal if
    !> the relief is flipped.
    subroutine slice_shape(c,rl,relief)
      type(slice_cache), intent(inout) :: c
      type(relief_def), intent(in) :: rl
      logical, intent(in) :: relief

      integer :: nx, ny, ix, iy, k, i0, i1, j0, j1
      real*8 :: nv(3), hmin, hmax
      real*8, allocatable :: xx(:,:,:), hh(:,:)
      logical, allocatable :: ok(:,:)

      nx = size(c%ff,1)
      ny = size(c%ff,2)
      c%isrelief = relief

      ! the heights; the range of the map is that of its finite values
      ! (on the log scale, of those not zero)
      allocate(hh(nx,ny))
      hh = 0d0
      if (relief) then
         do iy = 1, ny
            do ix = 1, nx
               hh(ix,iy) = scale_value(c%ff(ix,iy),rl%scale)
            end do
         end do
         if (rl%frange) then
            hmin = min(scale_value(rl%fmin,rl%scale),scale_value(rl%fmax,rl%scale))
            hmax = max(scale_value(rl%fmin,rl%scale),scale_value(rl%fmax,rl%scale))
         else
            ok = ieee_is_finite(c%ff) .and. ieee_is_finite(hh)
            if (rl%scale == 1) ok = ok .and. abs(c%ff) > 1d-30
            hmin = minval(hh,mask=ok)
            hmax = maxval(hh,mask=ok)
         end if
         if (hmax > hmin) then
            hh = (hh - hmin) / (hmax - hmin)
            where (ieee_is_finite(c%ff) .and. ieee_is_finite(hh))
               hh = merge(-rl%height,rl%height,rl%flip) * max(min(hh,1d0),0d0)
            elsewhere
               hh = 0d0
            end where
         else
            hh = 0d0
         end if
      end if

      ! the positions
      allocate(xx(3,nx,ny))
      do iy = 1, ny
         do ix = 1, nx
            xx(:,ix,iy) = c%x0 + real(ix-1,8) * c%uu + real(iy-1,8) * c%vv + hh(ix,iy) * nrm
         end do
      end do

      ! the vertices; the normals are that of the plane if flat, else
      ! from the differences between the neighbors along each side
      k = 0
      do iy = 1, ny
         do ix = 1, nx
            k = k + 1
            c%mesh%x(:,k) = real(xx(:,ix,iy),c_float)
            nv = nrm
            if (relief) then
               i0 = max(ix-1,1)
               i1 = min(ix+1,nx)
               j0 = max(iy-1,1)
               j1 = min(iy+1,ny)
               nv = cross(xx(:,i1,iy) - xx(:,i0,iy),xx(:,ix,j1) - xx(:,ix,j0))
               if (norm2(nv) > 1d-12) then
                  nv = nv / norm2(nv)
               else
                  nv = nrm
               end if
            end if
            c%mesh%nrm(:,k) = real(nv,c_float)
         end do
      end do

    end subroutine slice_shape

    !> The point p of the mesh of the map c at the grid coordinates ab
    !> (0 at the origin of the grid) and the normal of the mesh there
    !> (nl). On a relief, the point is in the triangle of the cell
    !> (ix,iy) (0-based) that holds it, and both are interpolated
    !> linearly from the vertices of that triangle.
    subroutine surface_point(c,ab,ix,iy,p,nl)
      type(slice_cache), intent(in) :: c
      real*8, intent(in) :: ab(2)
      integer, intent(in) :: ix, iy
      real*8, intent(out) :: p(3), nl(3)

      integer :: k, iv(3)
      real*8 :: t, s, w(3)

      if (.not.c%isrelief) then
         p = c%x0 + ab(1) * c%uu + ab(2) * c%vv
         nl = nrm
         return
      end if

      ! the triangle and the weights of its vertices
      t = min(max(ab(1) - ix,0d0),1d0)
      s = min(max(ab(2) - iy,0d0),1d0)
      k = iy * nx + ix + 1 ! the vertex (ix+1,iy+1)
      if (t >= s) then
         iv = (/k, k+1, k+nx+1/)
         w = (/1d0-t, t-s, s/)
      else
         iv = (/k, k+nx+1, k+nx/)
         w = (/1d0-s, t, s-t/)
      end if

      ! the point and the normal
      p = w(1) * c%mesh%x(:,iv(1)) + w(2) * c%mesh%x(:,iv(2)) + w(3) * c%mesh%x(:,iv(3))
      nl = w(1) * c%mesh%nrm(:,iv(1)) + w(2) * c%mesh%nrm(:,iv(2)) + w(3) * c%mesh%nrm(:,iv(3))
      if (norm2(nl) > 1d-12) then
         nl = nl / norm2(nl)
      else
         nl = nrm
      end if

    end subroutine surface_point

  end subroutine slice_draw

  !> The value f on a scale of a plane map (0 = linear, 1 = log|f|,
  !> 2 = atan).
  function scale_value(f,scale)
    real*8, intent(in) :: f
    integer, intent(in) :: scale
    real*8 :: scale_value

    if (scale == 1) then
       scale_value = log10(max(abs(f),1d-30))
    elseif (scale == 2) then
       scale_value = 2d0 / pi * atan(f)
    else
       scale_value = f
    end if

  end function scale_value

end submodule plane
