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

! Isosurface representations (reptype_isosurface): field binding, grid
! regions and sizes, cost estimates, and the sampled-field bookkeeping.
submodule (representations) iso
  implicit none

contains

  !> Default isovalue for field ifield of system isys. The policy, in
  !> order:
  !>
  !> - a promolecular density is a density by construction: in a
  !>   molecule use the conventional contour (iso_isoval_dens). A grid
  !>   field is not tested for being a density by integrating it: only
  !>   a smooth pseudo-density integrates to the electron count on its
  !>   grid (an all-electron one is tens of percent off, since a uniform
  !>   grid samples the nuclear cusps badly), and the valence charge is
  !>   only known if the user set the pseudopotential charges.
  !> - a grid that looks like an all-electron density
  !>   (field_guess_allelectron) uses iso_isoval_ae.
  !> - a spike-dominated grid (large max|f| over mean|f|: all-electron
  !>   or molecular densities, orbitals, laplacians) uses the level that
  !>   encloses iso_qcharge_def of the integral of |f|. Unlike the mean
  !>   it is not inflated by the core cusps, and unlike a plain value
  !>   quantile it does not follow the amount of vacuum in the box.
  !> - anything flat or bounded (ELF and the like, smooth valence
  !>   densities) keeps the old rule: twice the average of a
  !>   non-negative field, twice its rms if it is signed.
  !>
  !> The result is clamped to the range of the field. Non-grid fields
  !> other than the promolecular density have no cheap statistics, and
  !> fall back to a fixed value.
  module function iso_default_isovalue(isys,ifield) result(isoval)
    use systems, only: sys, sys_init, ok_system
    use fieldmod, only: type_grid, type_promol, type_promol_frag
    integer, intent(in) :: isys
    integer, intent(in) :: ifield
    real*8 :: isoval

    real*8 :: fmin, fmax, fmean, frms, qlevel, amean, rnuc, rval

    isoval = iso_isoval_def
    if (.not.ok_system(isys,sys_init)) return
    if (.not.sys(isys)%goodfield(ifield)) return

    ! the promolecular density is known to be a density
    if (sys(isys)%f(ifield)%type == type_promol .or.&
       sys(isys)%f(ifield)%type == type_promol_frag) then
       if (sys(isys)%c%ismolecule) isoval = iso_isoval_dens
       return
    end if

    ! everything else needs grid statistics
    if (sys(isys)%f(ifield)%type /= type_grid) return
    if (.not.sys(isys)%f(ifield)%grid%isinit) return
    call sys(isys)%f(ifield)%grid%stats(fmin=fmin,fmax=fmax,fmean=fmean,famean=amean,&
       frms=frms,qlevel=qlevel,qfrac=iso_qcharge_def)
    if (fmax <= fmin) return

    if (sys(isys)%f(ifield)%guess_allelectron(rnuc,rval)) then
       isoval = iso_isoval_ae
    elseif (max(abs(fmin),abs(fmax)) > iso_spikeratio * amean) then
       ! qlevel is a level of |f|: use it on the side the field lives on
       if (fmax > 0d0) then
          isoval = qlevel
       else
          isoval = -qlevel
       end if
    elseif (fmin >= 0d0) then
       isoval = 2d0 * fmean
    else
       isoval = 2d0 * frms
    end if
    if (isoval <= fmin .or. isoval >= fmax) isoval = 0.5d0 * (fmin + fmax)

  end function iso_default_isovalue

  !> Edge lengths (angstrom) of the box an isosurface of system isys
  !> samples: the given region box, else the valid window of the grid
  !> domain when the target field is a grid (a partial grid's box does
  !> not span the cell, and only part of it can feed the interpolation;
  !> matches iso_sample_domain), else the unit cell. Zero if the system
  !> is not available.
  module function iso_box_lengths(isys,ifield,box) result(alen)
    use systems, only: sys, sys_init, ok_system
    integer, intent(in) :: isys
    integer, intent(in), optional :: ifield
    real*8, intent(in), optional :: box(3,0:3)
    real*8 :: alen(3)

    integer :: i
    real*8 :: xmat(3,3), x0c(3), flo(3), fhi(3)

    alen = 0d0
    if (.not.ok_system(isys,sys_init)) return

    xmat = sys(isys)%c%m_x2c
    flo = 0d0
    fhi = 1d0
    if (present(box)) then
       xmat = box(:,1:3)
    elseif (present(ifield)) then
       if (iso_isgridfield(isys,ifield)) &
          call sys(isys)%f(ifield)%grid%get_domain(xmat,x0c,flo=flo,fhi=fhi)
    end if
    do i = 1, 3
       alen(i) = norm2(xmat(:,i)) * max(fhi(i)-flo(i),0d0) * bohrtoa
    end do

  end function iso_box_lengths

  !> Sampling-grid dimensions for an isosurface of system isys at
  !> coarseness level ilevel: 0 = native grid (returns all-zero),
  !> 1..iso_nlevel = named levels (points per angstrom over the sampled
  !> box), iso_level_custom = a custom grid, given either as the
  !> ncustom dimensions or, if ptsang is present, as a resolution in
  !> points per angstrom. On a named level ptsang instead overrides the
  !> level's tabulated resolution, for a caller working to a budget of
  !> its own. Never returns all-zero for a non-native level (that value
  !> is the native-grid sentinel).
  module function iso_grid_size(isys,ilevel,ncustom,capped,ifield,box,ptsang) result(n)
    use systems, only: sys_init, ok_system
    integer, intent(in) :: isys
    integer, intent(in) :: ilevel
    integer, intent(in), optional :: ncustom(3)
    logical, intent(out), optional :: capped
    integer, intent(in), optional :: ifield
    real*8, intent(in), optional :: box(3,0:3)
    real*8, intent(in), optional :: ptsang
    integer :: n(3)

    integer :: i, it, nmin, nmax, nc(3)
    real*8 :: pa, fac, alen(3)

    n = 0
    if (present(capped)) capped = .false.
    if (ilevel <= 0) return
    if (.not.ok_system(isys,sys_init)) return

    if (ilevel > iso_nlevel .and. .not.present(ptsang)) then
       ! custom, points mode: the dimensions come straight from the
       ! user, and asking for a custom level without them (and without a
       ! resolution) is a caller error -- answer with the same "no grid"
       ! zero a non-positive level gives, rather than silently returning
       ! a useless minimum-sized one
       if (.not.present(ncustom)) return
       nmin = iso_npts_custom_min
       nc = ncustom
       n = min(max(nc,nmin),iso_npts_custom_max)
    else
       ! points per angstrom along the axes of the sampled box: the
       ! level's tabulated value, or the one the caller gives. A custom
       ! resolution is taken literally, down to the floor a custom grid
       ! obeys; a named level keeps the coarser floor of the presets
       if (ilevel > iso_nlevel) then
          nmin = iso_npts_custom_min
          nmax = iso_npts_custom_max
          pa = ptsang
       else
          nmin = iso_npts_axmin
          nmax = huge(1)
          pa = iso_level_ptsang(ilevel)
          if (present(ptsang)) pa = ptsang
       end if
       alen = iso_box_lengths(isys,ifield=ifield,box=box)
       do i = 1, 3
          n(i) = min(max(nint(alen(i) * pa),nmin),nmax)
       end do
    end if

    ! total-points cap: re-scale until honest, since the per-axis floor
    ! can push the total back up in very anisotropic cells
    do it = 1, 10
       if (product(real(n,8)) <= real(iso_maxpts_total,8)) exit
       fac = (real(iso_maxpts_total,8) / product(real(n,8)))**(1d0/3d0)
       n = max(nint(n * fac),nmin)
       if (present(capped)) capped = .true.
    end do

  end function iso_grid_size

  !> Points per angstrom to sample a box of edge lengths alen
  !> (angstrom) within the iso_autogrid_secs budget, given costest,
  !> the measured cost of one sample point in seconds.  The budget
  !> only ever coarsens: the answer is pamax, the resolution the
  !> caller would use anyway, dropped to whatever the budget allows
  !> when that would take too long.
  module function iso_auto_ptsang(costest,alen,pamax) result(ppa)
    real*8, intent(in) :: costest
    real*8, intent(in) :: alen(3)
    real*8, intent(in) :: pamax
    real*8 :: ppa

    real*8 :: pa

    ppa = pamax
    if (costest <= 0d0) return
    pa = iso_ptsang_from_npts(iso_autogrid_secs / costest,alen)
    if (pa <= 0d0) return
    ppa = min(pa,pamax)

  end function iso_auto_ptsang

  !> The resolution (points per angstrom) that a total of ntot sample
  !> points amounts to over a box of edge lengths alen (angstrom).
  module function iso_ptsang_from_npts(ntot,alen) result(ppa)
    real*8, intent(in) :: ntot
    real*8, intent(in) :: alen(3)
    real*8 :: ppa

    real*8 :: vol

    ppa = 0d0
    vol = alen(1) * alen(2) * alen(3)
    if (vol <= 0d0 .or. ntot <= 0d0) return
    ppa = (ntot / vol)**(1d0/3d0)

  end function iso_ptsang_from_npts

  !> Name of the named coarseness level ilevel, as it appears in
  !> iso_level_optstr. That string is the one place the names live; this
  !> unpacks one of them so callers building their own option list do not
  !> have to know it is null-separated.
  module function iso_level_label(ilevel) result(str)
    integer, intent(in) :: ilevel
    character(len=:), allocatable :: str

    integer :: i, i0, i1

    str = ""
    if (ilevel < 1 .or. ilevel > iso_nlevel) return
    i0 = 1
    do i = 1, ilevel
       i1 = index(iso_level_optstr(i0:),c_null_char) + i0 - 1
       if (i == ilevel) str = iso_level_optstr(i0:i1-1)
       i0 = i1 + 1
    end do

  end function iso_level_label

  !> The field-evaluation request that samples the molecular orbital
  !> selected on this isosurface. Meaningful only when imosel is nonzero;
  !> the sampling loop, the cost benchmark and the MO window all have to
  !> ask for exactly the same thing or they describe different work.
  module function iso_mo_request(r) result(request)
    use types, only: field_evaluation_avail, fieldeval_category_mo
    class(rep_isosurface), intent(in) :: r
    type(field_evaluation_avail) :: request

    call request%clear()
    request%avail(fieldeval_category_mo) = .true.
    request%moini = r%imosel
    request%moend = r%imoidx

  end function iso_mo_request

  !> Convert the staged region inputs (mode iregion, origin/corner or
  !> center in column 0 of x; far corner, edge endpoints, or
  !> half-lengths in columns 1-3, per mode) of system isys into the
  !> cell-frame Cartesian box (origin in column 0, edge vectors in
  !> columns 1-3). ok is false if the box is degenerate (near-zero
  !> volume). Cartesian inputs are in user-frame angstrom (the
  !> molecular frame for molecules). The frac mode mirrors the CUBE
  !> keyword's x0/x1 convention (rhoplot_cube).
  module subroutine iso_region_to_box(isys,iregion,x,box,ok)
    use systems, only: sys, sys_init, ok_system
    use tools_math, only: det3
    integer, intent(in) :: isys
    integer, intent(in) :: iregion
    real*8, intent(in) :: x(3,0:3)
    real*8, intent(out) :: box(3,0:3)
    logical, intent(out), optional :: ok

    integer :: i
    real*8 :: u0(3), u1(3)

    real*8, parameter :: voleps = 1d-8 ! minimum region volume relative to the edge lengths

    box = 0d0
    if (present(ok)) ok = .false.
    if (.not.ok_system(isys,sys_init)) return
    associate(c => sys(isys)%c)
      if (iregion == iso_region_frac) then
         ! cell-aligned box between fractional points x(:,0) and x(:,1)
         box(:,0) = c%x2c(x(:,0))
         do i = 1, 3
            box(:,i) = c%m_x2c(:,i) * (x(i,1) - x(i,0))
         end do
      elseif (iregion == iso_region_ortho .or. iregion == iso_region_simplebox .or.&
         iregion == iso_region_cube) then
         ! axis-aligned Cartesian box, between corners x(:,0) and x(:,1)
         ! (ortho) or centered on x(:,0) with half-lengths x(:,1)
         ! (simple box and cube; the cube's single half-length is
         ! replicated across x(:,1) by its writers)
         if (iregion == iso_region_ortho) then
            u0 = x(:,0) / bohrtoa
            u1 = x(:,1) / bohrtoa
         else
            u0 = (x(:,0) - x(:,1)) / bohrtoa
            u1 = (x(:,0) + x(:,1)) / bohrtoa
         end if
         if (c%ismolecule) then
            u0 = u0 - c%molx0
            u1 = u1 - c%molx0
         end if
         box(:,0) = u0
         do i = 1, 3
            box(i,i) = u1(i) - u0(i)
         end do
      elseif (iregion == iso_region_bbox) then
         ! bounding box of the atoms, grown by the buffer x(1,1) in every
         ! direction. Derived from the atoms on every call, so it follows
         ! geometry edits; only the buffer is user state.
         if (c%ncel > 0) then
            u0 = c%atcel(1)%r
            u1 = c%atcel(1)%r
            do i = 2, c%ncel
               u0 = min(u0,c%atcel(i)%r)
               u1 = max(u1,c%atcel(i)%r)
            end do
            u0 = u0 - max(x(1,1),0d0) / bohrtoa
            u1 = u1 + max(x(1,1),0d0) / bohrtoa
            box(:,0) = u0
            do i = 1, 3
               box(i,i) = u1(i) - u0(i)
            end do
         end if
      elseif (iregion == iso_region_parallel) then
         ! origin x(:,0) plus the endpoints of the three edges (the
         ! molecular frame shift cancels in the edge vectors)
         box(:,0) = x(:,0) / bohrtoa
         if (c%ismolecule) box(:,0) = box(:,0) - c%molx0
         do i = 1, 3
            box(:,i) = (x(:,i) - x(:,0)) / bohrtoa
         end do
      else
         ! whole cell
         box(:,1:3) = c%m_x2c
      end if
    end associate
    ! degenerate when the volume is tiny relative to the edge lengths
    ! (also catches zero-length edges)
    if (present(ok)) ok = (abs(det3(box(:,1:3))) > &
       voleps * norm2(box(:,1)) * norm2(box(:,2)) * norm2(box(:,3)))

  end subroutine iso_region_to_box

  !> Seed the staged region inputs of mode iregion for system isys with
  !> a sensible starting box: the whole-cell equivalent, except for the
  !> molecule-centered modes (simple box and cube), which start centered
  !> on the molecular centroid and covering the atoms plus a margin.
  module subroutine iso_region_seed(isys,iregion,x)
    use systems, only: sys, sys_init, ok_system
    integer, intent(in) :: isys
    integer, intent(in) :: iregion
    real*8, intent(out) :: x(3,0:3)

    integer :: i
    real*8 :: xmin(3), xmax(3), xsh(3), xc(3), hl(3)

    real*8, parameter :: molmargin = 3d0 ! half-length margin beyond the atoms (ang)

    x = 0d0
    if (.not.ok_system(isys,sys_init)) return
    associate(c => sys(isys)%c)
      xsh = 0d0
      if (c%ismolecule) xsh = c%molx0
      if (iregion == iso_region_frac) then
         x(:,1) = 1d0
      elseif (iregion == iso_region_ortho) then
         ! Cartesian bounding box of the cell, in user-frame angstrom: the
         ! min/max of each component over the corners is the sum of the
         ! negative/positive entries of that row of the lattice matrix
         do i = 1, 3
            xmin(i) = sum(min(c%m_x2c(i,:),0d0))
            xmax(i) = sum(max(c%m_x2c(i,:),0d0))
         end do
         x(:,0) = (xmin + xsh) * bohrtoa
         x(:,1) = (xmax + xsh) * bohrtoa
      elseif (iregion == iso_region_parallel) then
         x(:,0) = xsh * bohrtoa
         do i = 1, 3
            x(:,i) = (c%m_x2c(:,i) + xsh) * bohrtoa
         end do
      elseif (iregion == iso_region_bbox) then
         ! only the buffer is stored; the box itself comes from the atoms
         x(:,1) = iso_bbox_buffer_def
      elseif (iregion == iso_region_simplebox .or. iregion == iso_region_cube) then
         ! origin at the molecular centroid; half-lengths covering the
         ! atoms plus a margin
         xc = 0d0
         do i = 1, c%ncel
            xc = xc + c%atcel(i)%r
         end do
         if (c%ncel > 0) xc = xc / real(c%ncel,8)
         hl = 0d0
         do i = 1, c%ncel
            hl = max(hl,abs(c%atcel(i)%r - xc))
         end do
         x(:,0) = (xc + xsh) * bohrtoa
         x(:,1) = hl * bohrtoa + molmargin
         ! the cube's single half-length is stored replicated in x(:,1)
         if (iregion == iso_region_cube) x(:,1) = maxval(x(:,1))
      end if
    end associate

  end subroutine iso_region_seed

  !> Store the cell-frame Cartesian point xc (bohr) into a staged
  !> region coordinate row x of mode iregion, converting to the mode's
  !> units (fractional for cell fractions, user-frame angstrom
  !> otherwise): the inverse of iso_region_to_box's point convention,
  !> kept next to it so the units of rgn_x live in this module only.
  module subroutine iso_region_point_from_cart(isys,iregion,xc,x)
    use systems, only: sys, sys_init, ok_system
    integer, intent(in) :: isys
    integer, intent(in) :: iregion
    real*8, intent(in) :: xc(3)
    real*8, intent(out) :: x(3)

    x = 0d0
    if (.not.ok_system(isys,sys_init)) return
    associate(c => sys(isys)%c)
      if (iregion == iso_region_frac) then
         x = c%c2x(xc)
      elseif (c%ismolecule) then
         x = (xc + c%molx0) * bohrtoa
      else
         x = xc * bohrtoa
      end if
    end associate

  end subroutine iso_region_point_from_cart

  !> Measure the field-evaluation cost (wall seconds per sample point)
  !> of sampling field ifield of system isys over the staged options
  !> (region mode iregion, coordinates x, n points per axis): evaluate
  !> the field at random points of the box the build would sample --
  !> in batches of increasing size, with the same OpenMP parallelism as
  !> the real sampling loop, until the benchmark has run long enough to
  !> be meaningful. The caller multiplies by its point count to show a
  !> total. Returns a negative value if the estimate cannot be made
  !> (bad system/field or degenerate region).
  module function iso_estimate_cost(isys,ifield,iregion,x,n,request) result(secs)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_init, ok_system
    use types, only: field_evaluation_avail, scalar_value
    integer, intent(in) :: isys
    integer, intent(in) :: ifield
    integer, intent(in) :: iregion
    real*8, intent(in) :: x(3,0:3)
    integer, intent(in) :: n(3)
    type(field_evaluation_avail), intent(in), optional :: request
    real*8 :: secs

    integer :: i, nb, ntot
    logical :: ok, per0, pereval, useres
    real*8 :: xmat(3,3), cmat(3,3), x0c(3), xp(3), rdum, t
    real*8, allocatable :: xr(:,:)
    type(scalar_value) :: res

    real*8, parameter :: timetarget = 0.1d0 ! benchmark until a batch takes this long
    integer, parameter :: nbatch0 = 32 ! initial batch size (a slow field stops after one batch)
    integer, parameter :: nptsmax = 65536 ! total benchmark points cap (fast fields)

    secs = -1d0
    if (.not.ok_system(isys,sys_init)) return
    if (.not.sys(isys)%goodfield(ifield)) return
    ! benchmark over the exact box the build would sample, with the same
    ! evaluation periodicity (iso_sample_domain owns the domain policy)
    call iso_sample_domain(isys,ifield,iregion,x,n,xmat,x0c,cmat,per0,pereval,ok)
    if (.not.ok) return

    ! batches of growing size: a slow field stops after the first one, a
    ! fast one grows until the timing is meaningful. The rate comes from
    ! the last batch alone, so the earlier ones (dominated by OpenMP
    ! startup) only calibrate the batch size.
    ! a caller that samples something other than the plain field -- a
    ! single molecular orbital, say -- passes the same request the real
    ! loop uses, or the benchmark would time the wrong quantity
    useres = present(request)
    ntot = 0
    nb = nbatch0
    do while (.true.)
       allocate(xr(3,nb))
       call random_number(xr)
       t = glfwGetTime()
       !$omp parallel do private(xp,res,rdum) schedule(static)
       do i = 1, nb
          xp = x0c + matmul(xmat,xr(:,i))
          if (useres) then
             call sys(isys)%f(ifield)%grd(xp,request,res,periodic=pereval)
          else
             rdum = sys(isys)%f(ifield)%grd0(xp,periodic=pereval)
          end if
       end do
       !$omp end parallel do
       t = glfwGetTime() - t
       ntot = ntot + nb
       deallocate(xr)
       if (t >= timetarget .or. ntot >= nptsmax) exit
       nb = min(4 * nb,nptsmax - ntot)
    end do
    ! seconds per sample point; the caller scales by its point count
    secs = t / real(nb,8)

  end function iso_estimate_cost

  !> Commit a staged sampling grid (n points per axis, region mode
  !> iregion with its staged coordinates x) as the applied state of
  !> isosurface iso. The applied region is stored as user coordinates
  !> (not as a Cartesian box) so it tracks later cell/molecule edits.
  !> The counterpart of grid_isapplied; keep the two in sync.
  module subroutine iso_apply_grid(iso,n,iregion,x)
    use interfaces_glfw, only: glfwGetTime
    class(rep_isosurface), intent(inout) :: iso
    integer, intent(in) :: n(3)
    integer, intent(in) :: iregion
    real*8, intent(in) :: x(3,0:3)

    iso%nptsxyz = n
    iso%iregion_ap = iregion
    iso%rgn_x_ap = x
    ! the single writer of the applied state, so one stamp here is all the
    ! renderer needs to know its samples are for an older grid or region
    iso%timelastapply_grid = glfwGetTime()

  end subroutine iso_apply_grid

  !> Stage and apply the sampling grid of isosurface iso (system isys)
  !> so that a newly created isosurface draws by itself, without the
  !> user having to commit one with the editor's Calculate grid
  !> button.
  module subroutine iso_autogrid(iso,isys)
    use systems, only: sys, sys_init, ok_system
    class(rep_isosurface), intent(inout) :: iso
    integer, intent(in) :: isys

    integer :: ilev0, n(3)
    real*8 :: box(3,0:3), pamax, pa
    logical :: okbox

    if (.not.ok_system(isys,sys_init)) return
    if (.not.sys(isys)%goodfield(iso%ifield)) return
    if (iso_isgridfield(isys,iso%ifield)) return

    ! the box the sampling covers; only a region mode needs one
    box = 0d0
    if (iso%iregion /= iso_region_cell) then
       call iso_region_to_box(isys,iso%iregion,iso%rgn_x,box,okbox)
       if (.not.okbox) return
    end if

    ! the level the isosurface was created with is the ceiling: this
    ! only coarsens a grid that would take too long, it does not
    ! override the preference with a finer one
    ilev0 = max(min(iso%ilevel,iso_nlevel),1)
    n = iso%staged_grid_size(isys,ilev0,box)
    if (all(n == 0)) return

    ! seconds per sample point over that grid
    iso%costest = iso%measure_cost(isys,n)
    if (iso%costest < 0d0) return

    ! the resolution the budget affords. If it is the level's own, that
    ! level is what gets applied -- a named level says more in the
    ! editor than the same grid spelled out as a custom one
    pamax = iso_level_ptsang(ilev0)
    pa = iso_auto_ptsang(iso%costest,iso%staged_box_lengths(isys,box),pamax)
    if (pa < pamax) then
       iso%ilevel = iso_level_custom
       iso%icustom = iso_custom_ptsang
       iso%ptsangcustom = min(max(pa,iso_ptsang_min),iso_ptsang_max)
       n = iso%staged_grid_size(isys,iso_level_custom,box)
    else
       iso%ilevel = ilev0
    end if
    if (all(n == 0)) return
    call iso%apply_grid(n,iso%iregion,iso%rgn_x)

  end subroutine iso_autogrid

  !> Sampling-grid dimensions of coarseness level ilevel over the
  !> staged options of isosurface iso. box is the staged region box,
  !> read only when the region is not the whole cell. A custom level is
  !> given by whichever of the two custom representations iso%icustom
  !> selects, or by icustom if the caller overrides it. capped is
  !> iso_grid_size's.
  module function iso_staged_grid_size(iso,isys,ilevel,box,capped,icustom) result(n)
    class(rep_isosurface), intent(in) :: iso
    integer, intent(in) :: isys
    integer, intent(in) :: ilevel
    real*8, intent(in) :: box(3,0:3)
    logical, intent(out), optional :: capped
    integer, intent(in), optional :: icustom
    integer :: n(3)

    integer :: ic

    ic = iso%icustom
    if (present(icustom)) ic = icustom

    if (ilevel > iso_nlevel .and. ic == iso_custom_ptsang) then
       if (iso%iregion == iso_region_cell) then
          n = iso_grid_size(isys,ilevel,capped=capped,ptsang=iso%ptsangcustom,ifield=iso%ifield)
       else
          n = iso_grid_size(isys,ilevel,capped=capped,ptsang=iso%ptsangcustom,box=box)
       end if
    elseif (iso%iregion == iso_region_cell) then
       n = iso_grid_size(isys,ilevel,iso%nptscustom,capped,ifield=iso%ifield)
    else
       n = iso_grid_size(isys,ilevel,iso%nptscustom,capped,box=box)
    end if

  end function iso_staged_grid_size

  !> Edge lengths (angstrom) of the box the staged options of
  !> isosurface iso sample. box is the staged region box, read only
  !> when the region is not the whole cell.
  module function iso_staged_box_lengths(iso,isys,box) result(alen)
    class(rep_isosurface), intent(in) :: iso
    integer, intent(in) :: isys
    real*8, intent(in) :: box(3,0:3)
    real*8 :: alen(3)

    if (iso%iregion == iso_region_cell) then
       alen = iso_box_lengths(isys,ifield=iso%ifield)
    else
       alen = iso_box_lengths(isys,box=box)
    end if

  end function iso_staged_box_lengths

  !> Benchmark the staged sampling of isosurface iso over an n(1) x
  !> n(2) x n(3) grid and return the seconds per sample point (negative
  !> when it cannot be measured). Prices what the sampling loop
  !> actually evaluates: an isosurface of a single molecular orbital is
  !> sampled with an MO request, not the plain field, and costs several
  !> times less per point, so benchmarking the density instead would
  !> buy the wrong grid.
  module function iso_measure_cost(iso,isys,n) result(secs)
    use types, only: field_evaluation_avail
    class(rep_isosurface), intent(in) :: iso
    integer, intent(in) :: isys
    integer, intent(in) :: n(3)
    real*8 :: secs

    type(field_evaluation_avail) :: request

    if (iso%imosel /= 0) then
       request = iso%mo_request()
       secs = iso_estimate_cost(isys,iso%ifield,iso%iregion,iso%rgn_x,n,request)
    else
       secs = iso_estimate_cost(isys,iso%ifield,iso%iregion,iso%rgn_x,n)
    end if

  end function iso_measure_cost

  !> Return true if the staged sampling grid (n, iregion, x) is already
  !> the applied state of isosurface iso. The region coordinates are
  !> irrelevant in whole-cell mode.
  module function iso_grid_isapplied(iso,n,iregion,x) result(isap)
    class(rep_isosurface), intent(in) :: iso
    integer, intent(in) :: n(3)
    integer, intent(in) :: iregion
    real*8, intent(in) :: x(3,0:3)
    logical :: isap

    isap = all(iso%nptsxyz == n) .and. (iso%iregion_ap == iregion)
    if (isap .and. iregion /= iso_region_cell) isap = all(iso%rgn_x_ap == x)

  end function iso_grid_isapplied

  !> The box sampled by an isosurface over field ifield of system isys
  !> with region mode iregion (coordinates x) on an n(1) x n(2) x n(3)
  !> grid: the region box if one is set, else the valid window of the
  !> grid domain for grid fields (a partial grid's box does not span
  !> the cell, and its faces cannot feed the interpolation stencil),
  !> else the unit cell. Returns the box origin x0c, its edge matrix
  !> xmat and inverse cmat, the mesh periodicity per0 (whether the
  !> triangulation wraps; a region never does), and the
  !> field-evaluation periodicity pereval (a property of the system:
  !> always periodic in crystals -- grd aborts otherwise -- so a region
  !> can span several cells; never periodic in molecules). Region and
  !> partial-grid boxes are scaled per axis by n/(n-1) so the n samples
  !> span them inclusively and the surface reaches the box faces; the
  !> periodic whole-cell paths keep the endpoint-exclusive convention.
  !> ok is false for a degenerate box, an empty interpolation window,
  !> or fewer than 2 points on an axis that needs the inclusive scaling.
  module subroutine iso_sample_domain(isys,ifield,iregion,x,n,xmat,x0c,cmat,per0,pereval,ok)
    use systems, only: sys
    use tools_math, only: matinv
    integer, intent(in) :: isys
    integer, intent(in) :: ifield
    integer, intent(in) :: iregion
    real*8, intent(in) :: x(3,0:3)
    integer, intent(in) :: n(3)
    real*8, intent(out) :: xmat(3,3)
    real*8, intent(out) :: x0c(3)
    real*8, intent(out) :: cmat(3,3)
    logical, intent(out) :: per0
    logical, intent(out) :: pereval
    logical, intent(out) :: ok

    integer :: i, ier
    real*8 :: fac(3), box(3,0:3), flo(3), fhi(3)

    ok = .true.
    fac = 1d0
    associate(c => sys(isys)%c)
      pereval = .not.c%ismolecule
      if (iregion /= iso_region_cell) then
         call iso_region_to_box(isys,iregion,x,box,ok)
         x0c = box(:,0)
         xmat = box(:,1:3)
         cmat = xmat
         if (ok) then
            call matinv(cmat,3,ier)
            ok = (ier == 0)
         end if
         per0 = .false.
         fac = real(n,8) / real(max(n-1,1),8)
      elseif (iso_isgridfield(isys,ifield)) then
         associate(g => sys(isys)%f(ifield)%grid)
           call g%get_domain(xmat,x0c,cmat,flo=flo,fhi=fhi)
           per0 = .not.c%ismolecule .and. .not.g%partial
           if (g%partial) then
              ! restrict the sampling to the window of the domain box
              ! where interpolation is valid (get_domain owns that rule)
              ! and sample it endpoint-inclusively, like a region box
              ok = all(fhi > flo)
              x0c = x0c + matmul(xmat,flo)
              fac = (fhi - flo) * real(n,8) / real(max(n-1,1),8)
           end if
         end associate
      else
         xmat = c%m_x2c
         cmat = c%m_c2x
         x0c = 0d0
         per0 = .not.c%ismolecule
      end if
    end associate

    ! apply the axis scaling; the inclusive paths need at least 2 points
    ! per axis to place samples on both faces
    if (any(fac /= 1d0)) ok = ok .and. all(n >= 2)
    if (ok) then
       do i = 1, 3
          xmat(:,i) = xmat(:,i) * fac(i)
          cmat(i,:) = cmat(i,:) / fac(i)
       end do
    end if

  end subroutine iso_sample_domain

  !> The box sampled by the applied state of isosurface iso in system
  !> isys (derived fresh from the applied user coordinates, so cell
  !> edits are tracked). Thin wrapper over iso_sample_domain, which
  !> owns the domain policy.
  module subroutine iso_sampled_box(iso,isys,xmat,x0c,cmat,per0,pereval,ok)
    class(rep_isosurface), intent(in) :: iso
    integer, intent(in) :: isys
    real*8, intent(out) :: xmat(3,3)
    real*8, intent(out) :: x0c(3)
    real*8, intent(out) :: cmat(3,3)
    logical, intent(out) :: per0
    logical, intent(out) :: pereval
    logical, intent(out) :: ok

    call iso_sample_domain(isys,iso%ifield,iso%iregion_ap,iso%rgn_x_ap,iso%nptsxyz,&
       xmat,x0c,cmat,per0,pereval,ok)

  end subroutine iso_sampled_box

  !> Return true if field ifield of system isys is an initialized grid
  !> field (the predicate behind the "native grid" coarseness level).
  module function iso_isgridfield(isys,ifield) result(isg)
    use systems, only: sys
    use fieldmod, only: type_grid
    integer, intent(in) :: isys
    integer, intent(in) :: ifield
    logical :: isg

    isg = sys(isys)%goodfield(ifield,type=type_grid)
    if (isg) isg = sys(isys)%f(ifield)%grid%isinit

  end function iso_isgridfield

  !> Select field ifield of system isys for isosurface iso: seed the
  !> default isovalue and reset the grid coarseness (native for grid
  !> fields, the preference level otherwise) together with its applied
  !> dimensions. All paths that change the isosurface field go through
  !> here, so the level/field invariant lives in one place.
  module subroutine iso_set_field(iso,isys,ifield)
    use systems, only: sys, sys_init, ok_system
    class(rep_isosurface), intent(inout) :: iso
    integer, intent(in) :: isys
    integer, intent(in) :: ifield

    logical :: isgrid

    if (.not.ok_system(isys,sys_init)) return
    iso%ifield = max(ifield,0)
    ! an MO selection belongs to the previous field
    iso%imosel = 0
    iso%imoidx = 0
    ! the isovalues of the previous field mean nothing for this one: keep
    ! a single isosurface (with its color) at the new default level
    if (iso%niso < 1) call iso%add_iso()
    iso%niso = 1
    iso%slot(1)%isoval = iso_default_isovalue(isys,iso%ifield)
    iso%slot(1)%built = .false.
    ! the colors mapped on the previous field's surface describe values that
    ! are gone. Dropping them is the whole statement: it is what the editor
    ! reads to know it has no legend to draw, and the rest of the map state
    ! (the values, the stamps) is rebuilt by color_slots on the resample this
    ! field change forces. A new quantity does get its colormap and its range
    ! suggested again
    if (allocated(iso%slot(1)%mesh%rgbv)) deallocate(iso%slot(1)%mesh%rgbv)
    iso%slot(1)%shown = .true. ! the surviving isosurface is the only one there is
    iso%slot(1)%icmap_auto = .true.
    iso%slot(1)%maprange_auto = .true.
    iso%slot(1)%maperr = ""
    isgrid = iso_isgridfield(isys,iso%ifield)
    if (isgrid) then
       iso%ilevel = 0
    else
       iso%ilevel = iso_defaultlevel
    end if
    ! a field change resets the region to the default for this system
    ! kind, and drops the cost estimate and the stale out-of-domain
    ! warning. A grid field comes with its own data domain, so the whole
    ! cell shows it at its native resolution right away; a field that
    ! has to be sampled gets, in a molecule (whose cell is an artifact),
    ! the box around the atoms instead of that cell
    if (isgrid .or. .not.sys(isys)%c%ismolecule) then
       iso%iregion = iso_region_cell
    else
       iso%iregion = iso_region_bbox
    end if
    call iso_region_seed(isys,iso%iregion,iso%rgn_x)
    iso%costest = -1d0
    iso%outdomain = .false.
    ! the histogram and the range describe the previous field's data
    iso%nhist = 0
    iso%hist_xscale = 0 ! the new field suggests its own axis scale
    iso%frange = (/1d0,-1d0/)
    ! apply an all-zero grid: for a grid field this is the native grid
    ! (cheap, built immediately); for any other field it means "not
    ! generated yet" -- sampling can be expensive, so nothing is built
    ! until the user commits a grid with the editor's Calculate grid
    ! button
    call iso%apply_grid((/0,0,0/),iso%iregion,iso%rgn_x)

  end subroutine iso_set_field

  !> Add an isosurface to isosurface object iso. Without isoval, the new
  !> level is chosen from the ones already there: the negative
  !> counterpart of the last isovalue if the data reach that far and no
  !> negative level is shown yet (the +/- pair of an orbital or a
  !> deformation density), otherwise half of it; the first isosurface of
  !> all takes the last-resort default. Without rgb, the color is the
  !> first one of the palette that no isosurface is using, so a new
  !> isosurface is visibly distinct from the others without the user
  !> having to pick.
  module subroutine iso_add_iso(iso,isoval,rgb)
    class(rep_isosurface), intent(inout) :: iso
    real*8, intent(in), optional :: isoval
    real(c_float), intent(in), optional :: rgb(3)

    integer :: i, j
    real*8 :: val
    logical :: used, haverange
    type(iso_slot), allocatable :: aux(:)

    ! grow the list
    if (.not.allocated(iso%slot)) allocate(iso%slot(4))
    if (iso%niso >= size(iso%slot)) then
       allocate(aux(2*size(iso%slot)))
       aux(1:iso%niso) = iso%slot(1:iso%niso)
       call move_alloc(aux,iso%slot)
    end if
    iso%niso = iso%niso + 1

    ! the level
    haverange = (iso%frange(1) <= iso%frange(2))
    if (present(isoval)) then
       val = isoval
    elseif (iso%niso == 1) then
       val = iso_isoval_def
    else
       val = iso%slot(iso%niso-1)%isoval
       if (haverange .and. all(iso%slot(1:iso%niso-1)%isoval > 0d0) .and. -val > iso%frange(1)) then
          val = -val
       else
          val = 0.5d0 * val
       end if
       ! a derived level outside the data shows nothing: fall back to
       ! halfway between the data limit and the outermost level there
       ! already, which is inside the data and distinct from all of them.
       ! (Not the middle of the range: value distributions are strongly
       ! skewed -- the reason the default level is a charge quantile --
       ! so the midpoint of a density is deep in the nuclear cusps.)
       if (haverange) then
          if (val <= iso%frange(1)) &
             val = 0.5d0 * (iso%frange(1) + minval(iso%slot(1:iso%niso-1)%isoval))
          if (val >= iso%frange(2)) &
             val = 0.5d0 * (iso%frange(2) + maxval(iso%slot(1:iso%niso-1)%isoval))
       end if
    end if

    ! the color: first palette entry no other isosurface is using
    iso%slot(iso%niso) = iso_slot(isoval=val)
    if (present(rgb)) then
       iso%slot(iso%niso)%rgb = rgb
    else
       do i = 1, iso_npalette
          used = .false.
          do j = 1, iso%niso-1
             used = used .or. all(abs(iso%slot(j)%rgb - iso_rgb_palette(:,i)) < iso_rgb_tol)
          end do
          if (.not.used) exit
       end do
       if (i > iso_npalette) i = 1 ! all of them are in use: start over
       iso%slot(iso%niso)%rgb = iso_rgb_palette(:,i)
    end if

  end subroutine iso_add_iso

  !> Remove isosurface i from isosurface object iso. The isosurfaces
  !> above it move down whole (mesh and build stamp included), so they
  !> are not rebuilt. An isosurface object with no isosurfaces left is
  !> valid but draws nothing, so the caller decides whether to keep the
  !> last one (the editor does).
  module subroutine iso_del_iso(iso,i)
    class(rep_isosurface), intent(inout) :: iso
    integer, intent(in) :: i

    integer :: j

    if (i < 1 .or. i > iso%niso) return
    do j = i, iso%niso-1
       iso%slot(j) = iso%slot(j+1)
    end do
    iso%niso = iso%niso - 1
    ! the vacated entry still holds a copy of the triangulation that moved
    ! down; drop it instead of keeping a whole mesh alive unreachable
    iso%slot(iso%niso+1) = iso_slot()

  end subroutine iso_del_iso

  !> Recompute the value range and the plottable histogram of
  !> isosurface object iso from the data ff. Stamped with the samples,
  !> so it describes whatever the meshes were last built from.
  module subroutine iso_stamp_histogram(iso,ff)
    use grid3mod, only: field_stats
    class(rep_isosurface), intent(inout) :: iso
    real*8, intent(in) :: ff(:,:,:)

    integer :: i, is, idef
    real*8 :: hy(iso_nhist,hscale_num), dh, e1, e2

    call field_stats(ff,fmin=iso%frange(1),fmax=iso%frange(2),hist=hy,&
       hrange=iso%hist_range,hhave=iso%hist_have,hcumq=iso%hist_q,&
       hcumv=iso%hist_v,hcumrange=iso%hist_cumrange,hdef=idef)
    ! take the suggested scale while the user has not chosen one, and
    ! whenever the chosen one is not available for this data
    if (iso%hist_xscale < 1 .or. iso%hist_xscale > hscale_num) then
       iso%hist_xscale = idef
    elseif (.not.iso%hist_have(iso%hist_xscale)) then
       iso%hist_xscale = idef
    end if
    ! staircase per scale, in field units (the axis applies the transform)
    do is = 1, hscale_num
       if (.not.iso%hist_have(is)) cycle
       dh = (iso%hist_range(2,is) - iso%hist_range(1,is)) / real(iso_nhist,8)
       do i = 1, iso_nhist
          e1 = iso%hist_range(1,is) + real(i-1,8)*dh
          e2 = iso%hist_range(1,is) + real(i,8)*dh
          if (is == hscale_log) then
             e1 = 10d0**e1
             e2 = 10d0**e2
          elseif (is == hscale_asinh) then
             e1 = 2d0 * sinh(0.5d0*e1)
             e2 = 2d0 * sinh(0.5d0*e2)
          end if
          iso%hist_x(2*i-1,is) = real(e1,c_double)
          iso%hist_x(2*i,is) = real(e2,c_double)
          iso%hist_y(2*i-1,is) = real(max(hy(i,is),0.5d0),c_double)
          iso%hist_y(2*i,is) = real(max(hy(i,is),0.5d0),c_double)
       end do
    end do
    ! a constant (or empty) field has no distribution to show
    if (any(iso%hist_have)) then
       iso%nhist = 2*iso_nhist
    else
       iso%nhist = 0
    end if

  end subroutine iso_stamp_histogram

  !> Stamp the keys of isosurface object iso in system isys that say
  !> its field samples are current: the field they were taken from, the
  !> MO they show, the generation of the system's field set, and the
  !> time. The single writer of the sample state, read back by the
  !> staleness test in add_isosurface_meshes.
  module subroutine iso_stamp_built(iso,isys)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys
    class(rep_isosurface), intent(inout) :: iso
    integer, intent(in) :: isys

    iso%ifield_built = iso%ifield
    iso%imosel_built = iso%imosel
    iso%imoidx_built = iso%imoidx
    iso%fieldgen_built = sys(isys)%fieldgen
    iso%time_built = glfwGetTime()

  end subroutine iso_stamp_built

  !> Sample the field of isosurface iso (system isys), or its MO imoidx
  !> if imosel is nonzero, on a grid of n points over the region iregion
  !> with coordinates x (the applied state or a staged one: iso is not
  !> modified, so a blocking job can sample first and apply the grid
  !> only if it was not cancelled), into ff. outdomain: some samples
  !> fell where the field cannot be evaluated (they are zeros). ok is
  !> false if the box is degenerate. The loop can be cancelled (ff is
  !> then incomplete) and reports its progress.
  module subroutine iso_sample(iso,isys,n,iregion,x,imosel,imoidx,ff,outdomain,ok)
    use systems, only: sys
    use global, only: abort_poll, progress_start, progress_step
    use types, only: scalar_value, field_evaluation_avail, fieldeval_category_mo
    class(rep_isosurface), intent(in) :: iso
    integer, intent(in) :: isys
    integer, intent(in) :: n(3)
    integer, intent(in) :: iregion
    real*8, intent(in) :: x(3,0:3)
    integer, intent(in) :: imosel
    integer, intent(in) :: imoidx
    real*8, allocatable, intent(inout) :: ff(:,:,:)
    logical, intent(out) :: outdomain
    logical, intent(out) :: ok

    integer :: j, k, l
    logical :: per0, pereval, lval, linvalid
    real*8 :: xp(3), xmat(3,3), cmat(3,3), x0c(3)
    type(scalar_value) :: res
    type(field_evaluation_avail) :: request

    outdomain = .false.
    call iso_sample_domain(isys,iso%ifield,iregion,x,n,xmat,x0c,cmat,per0,pereval,ok)
    if (.not.ok) return
    if (allocated(ff)) deallocate(ff)
    allocate(ff(n(1),n(2),n(3)))

    ! samples where the field cannot be evaluated -- outside its domain
    ! (e.g. a region beyond the cell of a grid field) or too close to a
    ! partial grid's edge for the interpolation stencil -- come back as
    ! zeros; remember that it happened so the editor can warn. An MO
    ! selection samples an orbital of the field instead of the field
    ! itself; an MO-only request makes grd skip the density work (the
    ! request of iso_mo_request)
    call request%clear()
    request%avail(fieldeval_category_mo) = .true.
    request%moini = imosel
    request%moend = imoidx
    linvalid = .false.
    call progress_start(n(3)*n(2),"lines")
    !$omp parallel do private(xp,res,lval) schedule(dynamic) collapse(2) reduction(.or.:linvalid)
    do l = 1, n(3)
       do k = 1, n(2)
          if (abort_poll()) cycle
          do j = 1, n(1)
             xp = x0c + matmul(xmat,(/real(j-1,8)/n(1),real(k-1,8)/n(2),real(l-1,8)/n(3)/))
             if (imosel /= 0) then
                call sys(isys)%f(iso%ifield)%grd(xp,request,res,periodic=pereval)
                lval = res%satisfied
                if (lval) then
                   ff(j,k,l) = res%fspc
                else
                   ff(j,k,l) = 0d0
                end if
             else
                ff(j,k,l) = sys(isys)%f(iso%ifield)%grd0(xp,periodic=pereval,valid=lval)
             end if
             linvalid = linvalid .or. .not.lval
          end do
          call progress_step()
       end do
    end do
    !$omp end parallel do
    outdomain = linvalid

  end subroutine iso_sample

  !> Install ff as the field samples of isosurface object iso in system
  !> isys, as if the renderer had just taken them: refresh the
  !> histogram, stamp the keys that say the samples are current for the
  !> selected field, MO, and applied grid, and drop the cached
  !> triangulations (the data changed even though the levels did not).
  !> The counterpart of the sampling loop in add_draw_elements, for a
  !> producer that samples or keeps its own grids -- the blocking jobs
  !> of the editor and the MO window, and the MO window's cache.
  module subroutine iso_set_samples(iso,isys,ff,outdomain)
    class(rep_isosurface), intent(inout) :: iso
    integer, intent(in) :: isys
    real*8, intent(in) :: ff(:,:,:)
    logical, intent(in) :: outdomain

    integer :: i

    iso%ff = ff
    iso%outdomain = outdomain
    call iso%stamp_histogram(ff)

    ! the samples describe this field, this MO, and this applied grid
    call iso%stamp_built(isys)

    ! the levels did not move but the data under them did (slot may be
    ! unallocated, so no array section here)
    do i = 1, iso%niso
       iso%slot(i)%built = .false.
    end do

  end subroutine iso_set_samples

  !> Return true if the applied state of isosurface iso in system isys
  !> describes a generated isosurface: a sampled grid has been
  !> committed, or the field is a native grid shown over the whole cell
  !> (built directly from its data). False means nothing is drawn until
  !> the user commits a grid with the editor's Calculate grid button.
  !> The single decoder of the all-zero nptsxyz sentinel; the renderer
  !> and the editor warning both call it.
  module function iso_isgenerated(iso,isys) result(gen)
    class(rep_isosurface), intent(in) :: iso
    integer, intent(in) :: isys
    logical :: gen

    gen = any(iso%nptsxyz /= 0) .or. &
       (iso_isgridfield(isys,iso%ifield) .and. iso%iregion_ap == iso_region_cell)

  end function iso_isgenerated

end submodule iso
