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

! Scene object and GL rendering utilities
submodule (representations) proc
  use param, only: bohrtoa
  implicit none

  ! extension of unit cell in the vacuum direction
  real*8, parameter :: vacextension = 2d0 / bohrtoa

contains

  !> Initialize a representation for system isys with ID irep.
  !> itype is the representation type, flavor is the representation
  !> flavor, and icount is the count array of the calling scene.
  module subroutine representation_init(r,isys,irep,itype,flavor,icount)
    use systems, only: sys_ready, ok_system
    use tools_io, only: string
    class(representation), intent(inout) :: r
    integer, intent(in) :: isys
    integer, intent(in) :: irep
    integer, intent(in) :: itype
    integer, intent(in) :: flavor
    integer, intent(inout) :: icount(0:repflavor_NUM)

    ! check the system is sane
    if (.not.ok_system(isys,sys_ready)) return

    ! save the basic information in the represntation from the arguments
    r%id = isys
    r%idrep = irep
    r%type = itype
    r%flavor = flavor

    ! default display options
    r%atoms%display = .true.
    r%bonds%display = .true.
    r%labels%display = .false.
    r%poly%display = .false.

    ! if type was given, mark as initialized and set name
    if (itype == reptype_atoms) then
       r%isinit = .true.
       r%shown = .true.
       if (flavor == repflavor_atoms_ballandstick) then
          r%name = "Ball and Stick"
       elseif (flavor == repflavor_atoms_sticks) then
          r%name = "Bonds"
          r%atoms%display = .false.
       elseif (flavor == repflavor_atoms_licorice) then
          r%name = "Licorice"
       elseif (flavor == repflavor_atoms_vdwcontacts) then
          r%name = "VdW Contacts"
          r%atoms%display = .false.
       elseif (flavor == repflavor_atoms_hbonds) then
          r%name = "Hydrogen Bonds"
          r%atoms%display = .false.
       elseif (flavor == repflavor_atoms_criticalpoints) then
          r%name = "Critical Points"
          r%bonds%display = .false.
          r%labels%display = .true.
       elseif (flavor == repflavor_atoms_gradientpaths) then
          r%name = "Gradient Paths"
          r%bonds%display = .false.
          r%labels%display = .false.
       elseif (flavor == repflavor_atoms_polyhedra) then
          r%name = "Polyhedra"
          r%bonds%display = .false.
          r%poly%display = .true.
       end if
    elseif (itype == reptype_unitcell) then
       r%isinit = .true.
       r%shown = .true.
       r%name = "Unit Cell"
    elseif (itype == reptype_axes) then
       r%isinit = .true.
       r%shown = .true.
       r%name = "Axes"
    elseif (itype == reptype_rotaxis) then
       r%isinit = .true.
       r%shown = .true.
       r%name = "Rotation axis"
    elseif (itype == reptype_symelem) then
       r%isinit = .true.
       r%shown = .true.
       r%name = "Symmetry elements"
    else
       r%isinit = .false.
       r%shown = .false.
       r%name = ""
    end if

    ! increment type counter and set name
    icount(flavor) = icount(flavor) + 1
    if (icount(flavor) > 1) then
       r%name = trim(r%name) // "/" // string(icount(flavor))
    end if

    ! increment global counter
    icount(0) = icount(0) + 1
    r%iord = icount(0)

    ! set all default values
    call r%set_defaults(0)

  end subroutine representation_init

  !> Reset a symmetry-element style.
  module subroutine symelem_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    class(symelem_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: nop
    logical, allocatable :: shownold(:)

    ! reset the time and remember the previous selection
    d%timelastreset = glfwGetTime()
    if (allocated(d%shown)) call move_alloc(d%shown,shownold)
    d%nop = 0
    d%isinit = .false.
    if (allocated(d%kind)) deallocate(d%kind)
    if (allocated(d%dir)) deallocate(d%dir)
    if (allocated(d%order)) deallocate(d%order)
    if (allocated(d%label)) deallocate(d%label)

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! recompute the operation snapshot
    call sys(r%id)%c%list_symops(nop,d%kind,d%dir,d%order,d%label)
    d%nop = nop

    ! restore the previous visibility if the operation count is unchanged,
    ! otherwise show all operations
    allocate(d%shown(nop))
    if (allocated(shownold)) then
       if (size(shownold,1) == nop) then
          d%shown = shownold
       else
          d%shown = .true.
       end if
    else
       d%shown = .true.
    end if
    d%isinit = .true.

  end subroutine symelem_style_reset

  !> Deallocate all arrays and end the symmetry-element style.
  module subroutine symelem_style_end(d)
    class(symelem_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    d%nop = 0
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%kind)) deallocate(d%kind)
    if (allocated(d%dir)) deallocate(d%dir)
    if (allocated(d%order)) deallocate(d%order)
    if (allocated(d%label)) deallocate(d%label)

  end subroutine symelem_style_end

  !> Set all values to default for the representation. Set a subset of
  !> defaults if itype = 0 (all), 1 (atom), 2 (bonds), 3 (labels),
  !> 4 (mol), 5 (unit cell), 6 (cartesian axes), 7 (rotation axes),
  !> 8 (coordination polyhedra), 9 (symmetry elements).
  module subroutine representation_set_defaults(r,itype)
    use systems, only: sys, sys_ready, ok_system
    use global, only: bondfactor_def, bonddelta_def
    use gui_main, only: ColorAtomBorder_def, ColorBond_def, ColorBondBorder_def,&
       ColorLabel_def, ColorRotaxis_def, ColorAxes_def, ColorVdwContacts_def,&
       ColorHbonds_def, ColorHbondStrong_def, ColorHbondModerate_def, ColorHbondWeak_def
    use param, only: atmcov0, atmvdw0
    class(representation), intent(inout) :: r
    integer, intent(in) :: itype

    integer :: isys, imol, iat

    ! check the system is sane
    isys = r%id
    if (.not.ok_system(isys,sys_ready)) return

    !! initialize an empty representation
    if (itype == 0) then
       ! selection group
       r%sel%pertype = 1
       if (sys(isys)%c%ismolecule) then
          r%sel%ncell = 0
       else
          r%sel%ncell = 1
       end if
       r%sel%origin = 0d0
       r%sel%tshift = 0d0
       r%sel%filter = ""
       r%sel%errfilter = ""
       if (sys(isys)%c%ismolecule) then
          r%sel%border = .false.
          r%sel%onemotif = .false.
       else
          r%sel%border = .true.
          ! show connected molecules by default if there is more than one
          ! fragment, or if a non-discrete fragment carries reconnection
          ! lattice vectors (dangling pieces split by the cell boundary)
          r%sel%onemotif = (sys(isys)%c%nmol > 1)
          if (.not.r%sel%onemotif) then
             moldangler: do imol = 1, sys(isys)%c%nmol
                if (sys(isys)%c%mol(imol)%discrete) cycle
                do iat = 1, sys(isys)%c%mol(imol)%nat
                   if (any(sys(isys)%c%mol(imol)%at(iat)%lvec /= 0)) then
                      r%sel%onemotif = .true.
                      exit moldangler
                   end if
                end do
             end do moldangler
          end if
       end if
    end if

    !--> atoms
    if (itype == 0 .or. itype == 1) then
       r%atoms%radii_type = 0
       r%atoms%radii_scale = atomcovradscale_def
       r%atoms%radii_value = atomconstantrad_def
       r%atoms%color_type = 0
       r%atoms%border_size = atomborder_def
       r%atoms%border_rgb = ColorAtomBorder_def
       if (r%flavor == repflavor_atoms_licorice) then
          r%atoms%radii_type = 2
          r%atoms%radii_value = atomrad_licorice_def
       elseif (r%flavor == repflavor_atoms_criticalpoints) then
          r%atoms%radii_type = 2
          r%atoms%radii_value = atomrad_criticalpoints_def
          r%atoms%border_size = atomborder_criticalpoints_def
       elseif (r%flavor == repflavor_atoms_gradientpaths) then
          r%atoms%radii_type = 2
          r%atoms%radii_value = atomrad_gradientpaths_def
          r%atoms%border_size = atomborder_gradientpaths_def
       end if
    end if

    !--> bonds
    if (itype == 0 .or. itype == 2) then
       r%bonds%atmrad = atmcov0
       r%bonds%bfactor = bondfactor_def
       r%bonds%bdelta = bonddelta_def
       r%bonds%color_style = 0
       r%bonds%border_size = bondborder_def
       r%bonds%rad = bondrad_def
       r%bonds%border_rgb = ColorBondBorder_def
       r%bonds%rgb = ColorBond_def
       r%bonds%order = 4 ! calculated (value from ordcon)
       r%bonds%imol = 0
       r%bonds%bothends = .true.
       r%bonds%hbond_classify = .false.
       if (r%flavor == repflavor_atoms_sticks) then
          r%bonds%color_style = 1
          r%bonds%border_size = bondborder_stickflav_def
       elseif (r%flavor == repflavor_atoms_licorice) then
          r%bonds%color_style = 1
          r%bonds%rad = bondrad_licorice_def
       elseif (r%flavor == repflavor_atoms_vdwcontacts) then
          ! van der waals contacts: dashed, intermolecular-only bonds using
          ! the sum of the van der Waals radii as the distance cutoff
          r%bonds%atmrad = atmvdw0
          r%bonds%bfactor = bondfactor_vdwcontacts_def
          r%bonds%order = 0 ! dashed
          r%bonds%imol = 2 ! intermolecular only
          r%bonds%bothends = .false.
          r%bonds%rad = bondrad_vdwcontacts_def
          r%bonds%border_size = 0d0
          r%bonds%rgb = ColorVdwContacts_def
       elseif (r%flavor == repflavor_atoms_hbonds) then
          ! hydrogen bonds: dashed, intermolecular-only contacts.
          ! Jeffrey-Steiner strength classification (distance + D-H...A angle) is applied at render time.
          r%bonds%atmrad = atmvdw0
          r%bonds%bfactor = bondfactor_hbonds_def
          r%bonds%order = 0 ! dashed
          r%bonds%imol = 2 ! intermolecular only
          r%bonds%bothends = .false.
          r%bonds%rad = bondrad_vdwcontacts_def
          r%bonds%border_size = 0d0
          r%bonds%rgb = ColorHbonds_def
          r%bonds%hbond_classify = .true.
          r%bonds%hbond_rgb(:,1) = ColorHbondStrong_def
          r%bonds%hbond_rgb(:,2) = ColorHbondModerate_def
          r%bonds%hbond_rgb(:,3) = ColorHbondWeak_def
          r%bonds%hbond_dist = hbond_dist_def
          r%bonds%hbond_ang = hbond_ang_def
       end if
    end if

    !--> labels
    if (itype == 0 .or. itype == 3) then
       r%labels%type = 0
       r%labels%scale = 0.5d0
       r%labels%rgb = ColorLabel_def
       r%labels%const_size = .false.
       r%labels%offset = (/0d0,0d0,0d0/)
       if (r%flavor == repflavor_atoms_criticalpoints) then
          r%labels%type = 2 ! cell atom
          r%labels%scale = 0.3d0
          r%labels%offset = (/0d0,0.25d0,0d0/)
       end if
    end if

    ! unit cell
    if (itype == 0 .or. itype == 5) then
       r%uc%inner = .true.
       r%uc%coloraxes = .true.
       r%uc%vaccutsticks = .true.
       r%uc%radius = uc_radius_def
       r%uc%radiusinner = uc_radiusinner_def
       r%uc%rgb = 0._c_float
       r%uc%innersteplen = uc_innersteplen_def
       r%uc%innerstipple = .true.
    end if

    ! cartesian axes
    if (itype == 0 .or. itype == 6) then
       r%axes%kind = 0 ! cartesian
       r%axes%rot = eye ! no extra orientation by default
       r%axes%placement = 1
       r%axes%origin = 0d0
       if (sys(isys)%c%ismolecule) then
          r%axes%coordtype = 1 ! cartesian (angstrom)
       else
          r%axes%coordtype = 0 ! crystallographic
       end if
       r%axes%winpos = axes_winpos_def
       r%axes%length = axes_length_def
       r%axes%radius = axes_radius_def
       r%axes%conelength = axes_conelength_def
       r%axes%coneradius = axes_coneradius_def
       r%axes%rgb = ColorAxes_def ! x = red, y = green, z = blue
       r%axes%showlabels = .true.
       r%axes%labelscale = axes_labelscale_def
       r%axes%labelconstsize = .false.
       r%axes%labeldistance = axes_labeldistance_def
       r%axes%labeloffset = 0d0
       r%axes%scalewithzoom = .false.
       r%axes%scale = 1d0
       r%axes%scale_auto = (r%axes%placement == 1)
       r%axes%labelrgb = 0._c_float
       r%axes%labelstr(1) = "x"
       r%axes%labelstr(2) = "y"
       r%axes%labelstr(3) = "z"
    end if

    ! rotation axis
    if (itype == 0 .or. itype == 7) then
       r%rotaxis%origin = 0d0
       r%rotaxis%dir = (/0d0,0d0,1d0/)
       r%rotaxis%length = 0d0
       r%rotaxis%radius = rotaxis_radius_def
       r%rotaxis%rgb = ColorRotaxis_def ! black
    end if

    ! coordination polyhedra
    if (itype == 0 .or. itype == 8) then
       r%poly%alpha = 0.5d0
       r%poly%usecentercolor = .true.
       r%poly%rgb = 0._c_float
       r%poly%edge_rad = 0.05d0
       r%poly%edge_rgb = 0._c_float
       r%poly%usecentercolor_edge = .true.
       r%poly%coplanar_eps = 0.1d0
       r%poly%showcorners = .true.
    end if

    ! symmetry elements
    if (itype == 0 .or. itype == 9) then
       r%symelem%origin_transient = 0d0
       if (sys(isys)%c%ismolecule) then
          r%symelem%coordtype = 2 ! cartesian (bohr)
          if (sys(isys)%c%pg%avail) then
             r%symelem%origin = sys(isys)%c%pg%xcm + sys(isys)%c%molx0
          else
             r%symelem%origin = 0d0
          end if
       else
          r%symelem%coordtype = 0 ! crystallographic
          r%symelem%origin = 0d0
       end if
       r%symelem%usecustomrgb = .false.
       r%symelem%rgb = symelem_rgb_def
       if (r%symelem%style%isinit) r%symelem%style%shown = .true.
    end if

    ! initialize the styles
    call r%reset_all_styles(itype)

  end subroutine representation_set_defaults

  !> Terminate a representation
  module subroutine representation_end(r)
    class(representation), intent(inout) :: r

    r%isinit = .false.
    r%shown = .false.
    r%type = reptype_none
    r%flavor = repflavor_unknown
    r%id = 0
    r%idrep = 0
    r%iord = 0
    r%name = ""
    r%sel%filter = ""
    r%sel%errfilter = ""

    call r%atoms%style%end()
    call r%bonds%style%end()
    call r%labels%style%end()
    call r%mols%style%end()
    call r%poly%style%end()
    call r%symelem%style%end()

  end subroutine representation_end

  !> Update the representation styles to respond to changes in the the
  !> associated system.
  module subroutine update_styles(r)
    use systems, only: sys_ready, ok_system, sysc
    class(representation), intent(inout) :: r

    logical :: doreset

    ! consistency checks
    if (.not.r%isinit .or. r%id == 0) return
    if (.not.ok_system(r%id,sys_ready)) return

    if (r%type == reptype_atoms) then
       ! check if we need to reset the representation styles
       ! atoms
       doreset = .not.r%atoms%style%isinit
       doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%atoms%style%timelastreset)
       if (doreset) call r%atoms%style%reset(r)

       ! bonds: if the geometry changed
       doreset = .not.r%bonds%style%isinit
       doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%bonds%style%timelastreset)
       if (doreset) call r%bonds%style%reset(r)

       ! bonds: if the system has been rebonded and this representation tracks
       ! the bonds in the system (%use_sys_nstar), recalculate the bond style
       doreset = r%bonds%style%use_sys_nstar .and. (sysc(r%id)%timelastchange_rebond > r%bonds%style%timelastreset)
       if (doreset) call r%bonds%style%copy_neighstars_from_system(r%id)

       ! molecules: if the geometry or the bonds changed
       doreset = .not.r%mols%style%isinit
       doreset = doreset .or. (sysc(r%id)%timelastchange_rebond > r%mols%style%timelastreset)
       if (doreset) call r%mols%style%reset(r)

       ! labels: if the geometry changed
       doreset = .not.r%labels%style%isinit
       doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%labels%style%timelastreset)
       if (doreset) call r%labels%style%reset(r)

       ! coordination polyhedra: if the geometry changed
       doreset = .not.r%poly%style%isinit
       doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%poly%style%timelastreset)
       if (doreset) call r%poly%style%reset(r)

    elseif (r%type == reptype_symelem) then
       ! symmetry elements: if the geometry changed
       doreset = .not.r%symelem%style%isinit
       doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%symelem%style%timelastreset)
       if (doreset) call r%symelem%style%reset(r)
    end if

  end subroutine update_styles

  !> Add the spheres, cylinder, etc. to the draw lists. Use nc number
  !> of cells and the data from representation r. If doanim, use qpt
  !> iqpt and frequency ifreq to animate the representation.
  module subroutine add_draw_elements(r,nc,obj,doanim,iqpt,ifreq)
    use systems, only: sys, sysc
    use systemmod, only: system
    use crystalmod, only: crystal, iperiod_vacthr, symop_kind_plane
    use arithmetic, only: pretokenize, token
    use gui_main, only: ColorAxes_def
    use tools_io, only: string
    use tools_math, only: cross, plane_from_points
    use types, only: realloc
    use tools, only: mergesort
    use param, only: tpi, img, atmass, icrd_crys, pi
    class(representation), intent(inout) :: r
    integer, intent(in) :: nc(3)
    type(scene_objects), intent(inout) :: obj
    logical, intent(in) :: doanim
    integer, intent(in) :: iqpt, ifreq

    logical, allocatable :: lshown(:,:,:,:)
    logical :: havefilter, step, isedge(3), usetshift, doanim_, dobonds, isvac(3)
    logical :: isvacdir, docycle, dovac(3)
    integer :: n(3), i, j, k, imol, lvec(3), id, idaux, n0(3), n1(3)
    integer :: i1, i2, i3, ix(3), idl
    integer :: ib, ineigh, ixn(3), ix1(3), ix2(3), nstep, vacshift(3), iord
    integer :: nimg, nres, nbond, mb, mbb
    real(c_float) :: rgb(3)
    real*8 :: rad1, rad2, dd, f1, f2, axsc
    integer, allocatable :: hbcat(:) ! per-bond H-bond class cache for the current atom
    real(c_float) :: bondrgb(3)
    type(crystal), pointer :: c ! the system's crystal structure (sys(r%id)%c)
    type(system), pointer :: syptr
    type(token), allocatable :: toklist(:) ! pre-tokenized filter expression
    complex*16, allocatable :: vibbase(:,:) ! per-atom vibration phasors (3,ncel)
    real*8 :: xx(3), xc(3), x0(3), x1(3), x2(3), res, uoriginc(3), xpolyc(3)
    real*8 :: ucini(3), ucend(3)
    complex*16 :: xdelta0(3), xdelta1(3), xdelta2(3)
    type(dl_sphere) :: dsph
    type(dl_cylinder) :: dcyl
    type(dl_cylinder_giz) :: dcylgiz
    type(dl_string) :: dstr
    type(dl_string_giz) :: dstrgiz
    logical :: fixed
    character(len=:), allocatable :: errmsg
    real*8, allocatable :: xvpoly(:,:), up2dsp(:,:)
    complex*16, allocatable :: dvpoly(:,:)
    integer, allocatable :: eidp(:), lvecp(:,:)
    real(c_float) :: rgbface(3), rgbedge(3)
    integer :: natp, idpoly, kp
    logical :: dopoly, corneractive
    integer, allocatable :: cornlist(:,:)
    integer :: ncorn, ica, idc, imolc

    interface
       subroutine runqhull_basintriangulate_step1(n,x0,xvert,nf,ctx,ier) bind(c)
         use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
         integer(c_int), value :: n
         real(c_double) :: x0(3)
         real(c_double) :: xvert(3,n)
         integer(c_int) :: nf
         type(c_ptr) :: ctx
         integer(c_int) :: ier
       end subroutine runqhull_basintriangulate_step1
       subroutine runqhull_basintriangulate_step2(nf,iface,ctx) bind(c)
         use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
         integer(c_int), value :: nf
         integer(c_int) :: iface(3,nf)
         type(c_ptr), value :: ctx
       end subroutine runqhull_basintriangulate_step2
    end interface

    real*8, parameter :: rthr = 0.01d0
    real*8, parameter :: rthr1 = 1-rthr
    integer, parameter :: uc(3,2,12) = reshape((/&
       0,0,0,  1,0,0,&
       0,0,0,  0,1,0,&
       0,0,0,  0,0,1,&
       1,1,0,  1,1,1,&
       1,0,1,  1,1,1,&
       0,1,1,  1,1,1,&
       1,0,0,  1,1,0,&
       0,1,0,  1,1,0,&
       0,1,0,  0,1,1,&
       0,0,1,  0,1,1,&
       0,0,1,  1,0,1,&
       1,0,0,  1,0,1/),shape(uc))
    integer, parameter :: ucdir(12) = (/1, 2, 3, 3, 2, 1, 2, 1, 3, 2, 1, 3/)

    ! system has been initialized: ensured by scene_build_lists, which
    ! calls this routine.

    ! return if not initialized
    if (.not.r%isinit) return
    if (.not.r%shown) return

    ! the system's crystal structure (also visible in the contained routines)
    c => sys(r%id)%c

    ! the draw lists have been reset/allocated by scene_build_lists
    ! (scene_objects%reset); dl_append grows them as needed
    doanim_ = doanim
    if (doanim_) doanim_ = doanim_ .and. (iqpt > 0 .and. ifreq > 0 .and. c%vib%hasvibs)

    ! precompute the per-atom vibration phasors: displacement of atom iat at
    ! lattice translation L is vibbase(:,iat) * exp(i 2 pi q.L) (see vibdelta)
    if (doanim_) then
       allocate(vibbase(3,c%ncel))
       do i = 1, c%ncel
          vibbase(:,i) = c%vib%vec(:,i,ifreq,iqpt) * &
             exp(img * tpi * dot_product(c%atcel(i)%x,c%vib%qpt(:,iqpt))) / &
             sqrt(atmass(c%spc(c%atcel(i)%is)%z))
       end do
    end if

    if (r%type == reptype_atoms) then
       !!! atoms and bonds representation !!!

       !! first, the atoms
       ! do we have a filter? If so, tokenize it once here; the evaluation for
       ! each atom image below reuses the token list (skips the string parsing)
       havefilter = (len_trim(r%sel%filter) > 0) .and. (len_trim(r%sel%errfilter) == 0)
       if (havefilter) then
          syptr => sys(r%id)
          errmsg = ""
          call pretokenize(r%sel%filter,toklist,errmsg,c_loc(syptr))
          if (len_trim(errmsg) > 0) then
             havefilter = .false.
             r%sel%errfilter = errmsg
          end if
       end if
       usetshift = any(abs(r%sel%tshift) > 1d-5)

       ! calculate the periodicity
       n = 1
       if (r%sel%pertype == 1) then
          n = nc
       elseif (r%sel%pertype == 2) then
          n = r%sel%ncell
       end if

       ! origin shift
       if (c%ismolecule) then
          uoriginc = r%sel%origin / bohrtoa
       else
          uoriginc = c%x2c(r%sel%origin)
       end if

       ! whether we will force the polyhedra corner atoms to be drawn (only when
       ! the representation displays atoms, so lshown tracks visible spheres)
       corneractive = r%poly%display .and. r%poly%showcorners .and. r%atoms%display .and.&
          r%poly%style%isinit

       ! whether we'll be doing bonds, allocate array to check whether
       ! an atom has been drawn (also used to deduplicate forced corner atoms)
       dobonds = r%bonds%display .and. r%bonds%style%isinit
       if (dobonds .or. corneractive) then
          ! bound the lattice translations reachable by the loops below so the
          ! atom and bond sites can index lshown directly (no bounds checks):
          ! the base image range is [-1,n] (border), plus the vacuum shift (1),
          ! the tshift rounding (the adjustment is -floor(x-tshift), bounded by
          ! ceiling(|tshift|)+1 because the atomic coordinates x are in [0,1)),
          ! the molecule lattice vectors (onemotif), and the neighbor-star
          ! connectivity vectors (bonds). Polyhedra corner images are not
          ! bounded by this; that path keeps using check_lshown, which grows
          ! the array as needed.
          mb = 1
          if (usetshift) mb = mb + maxval(ceiling(abs(r%sel%tshift))) + 1
          if (r%sel%onemotif) then
             mbb = 0
             do imol = 1, c%nmol
                do k = 1, c%mol(imol)%nat
                   mbb = max(mbb,maxval(abs(c%mol(imol)%at(k)%lvec)))
                end do
             end do
             mb = mb + mbb
          end if
          if (dobonds) then
             mbb = 0
             do i = 1, c%ncel
                if (r%bonds%style%nstar(i)%ncon > 0) &
                   mbb = max(mbb,maxval(abs(r%bonds%style%nstar(i)%lcon(:,1:r%bonds%style%nstar(i)%ncon))))
             end do
             mb = mb + mbb
          end if
          allocate(lshown(c%ncel,-1-mb:n(1)+mb,-1-mb:n(2)+mb,-1-mb:n(3)+mb))
          lshown = .false.
       end if

       ! presize the draw lists from the known atom/bond counts (the 2x growth
       ! in dl_append absorbs border/vacuum extras and forced polyhedra corners)
       nimg = product(n)
       nres = c%ncel * nimg
       if (r%atoms%display .or. dobonds) call obj%reserve(nsph = obj%nsph + nres)
       if (r%labels%display) call obj%reserve(nstring = obj%nstring + nres)
       if (dobonds) then
          nbond = sum(r%bonds%style%nstar(1:c%ncel)%ncon) / 2
          if (r%bonds%color_style /= 0) nbond = 2*nbond
          call obj%reserve(ncyl = obj%ncyl + nbond*nimg)
       end if
       ncorn = 0
       if (corneractive) allocate(cornlist(4,100))

       ! coordination polyhedra: per-species distance-window scratch
       if (r%poly%display .and. r%poly%style%isinit) &
          allocate(up2dsp(c%nspc,2))

       ! whether there is vacuum in any direction
       dovac = (c%vaclength > iperiod_vacthr)
       if (any(dovac)) then
          ucini = c%vactop - 1d0 - vacextension / c%aa
          ucend = c%vacbot + vacextension / c%aa
       end if

       ! run over atoms, either directly or per-molecule
       i = 0
       imol = 0
       do while(.true.)
          if (r%sel%onemotif) then
             ! this is a new molecule if there are no molecules or this is the last atom
             ! in the previous one
             step = (imol == 0)
             if (.not.step) step = (k == c%mol(imol)%nat)
             if (step) then
                imol = imol + 1
                k = 0
             end if

             ! we are finished if we have all molecules
             if (imol > c%nmol) exit

             ! Add the new atom, translated by the molecule lattice vector
             k = k + 1
             i = c%mol(imol)%at(k)%cidx
             lvec = c%mol(imol)%at(k)%lvec
          else
             ! next atom in the complete list, exit if done
             i = i + 1
             if (i > c%ncel) exit
             lvec = 0
             imol = c%idatcelmol(1,i)
          end if
          ! i is current atom from the complete atom list
          ! imol is the corresponding molecule

          ! skip hidden atoms
          id = sysc(r%id)%attype_celatom_to_id(r%atoms%style%type,i)
          if (.not.r%atoms%style%shown(id)) cycle

          ! skip hidden molecules
          if (.not.r%mols%style%shown(imol)) cycle

          ! calculate the border
          xx = c%atcel(i)%x
          n0 = 0
          n1 = n-1
          if (r%sel%border.and..not.r%sel%onemotif) then
             do j = 1, 3
                ! not in a vacuum direction
                if (.not.dovac(j)) then
                   if (xx(j) < rthr) then
                      n1(j) = n(j)
                   elseif (xx(j) > rthr1) then
                      n0(j) = -1
                   end if
                end if
             end do
          end if

          ! calculate the vacuum shift
          vacshift = 0
          if (any(dovac)) then
             xx = c%atcel(i)%x + lvec
             do j = 1, 3
                if (dovac(j)) then
                   if (xx(j) < ucini(j)) then
                      vacshift(j) = 1
                   elseif (xx(j) > ucend(j)) then
                      vacshift(j) = -1
                   end if
                end if
             end do
          end if

          ! draw the spheres and cylinders
          rgb = r%atoms%style%rgb(:,id) * r%mols%style%tint_rgb(:,imol)
          rad1 = r%atoms%style%rad(id) * r%mols%style%scale_rad(imol)

          ! coordination polyhedron: if this atom is a shown center, find its
          ! corner atoms once (translation-invariant; reused for every image)
          dopoly = .false.
          if (r%poly%display .and. r%poly%style%isinit) then
             idpoly = sysc(r%id)%attype_celatom_to_id(r%poly%style%type,i)
             if (r%poly%style%shown(idpoly) .and. r%poly%style%dmax(idpoly) > 0d0 .and.&
                any(r%poly%style%corner(:,idpoly))) then
                do j = 1, c%nspc
                   if (r%poly%style%corner(j,idpoly)) then
                      up2dsp(j,1) = r%poly%style%dmin(idpoly)
                      up2dsp(j,2) = r%poly%style%dmax(idpoly)
                   else
                      up2dsp(j,:) = 0d0
                   end if
                end do
                call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.true.,&
                   natp,eid=eidp,lvec=lvecp,up2dsp=up2dsp,nozero=.true.)
                dopoly = (natp >= 3)
                if (dopoly) then
                   ! grow-only corner buffers (reused across center atoms)
                   if (allocated(xvpoly)) then
                      if (size(xvpoly,2) < natp) deallocate(xvpoly)
                   end if
                   if (.not.allocated(xvpoly)) allocate(xvpoly(3,natp))
                   if (allocated(dvpoly)) then
                      if (size(dvpoly,2) < natp) deallocate(dvpoly)
                   end if
                   if (.not.allocated(dvpoly)) allocate(dvpoly(3,natp))
                   if (r%poly%usecentercolor) then
                      rgbface = rgb
                   else
                      rgbface = r%poly%rgb
                   end if
                   if (r%poly%usecentercolor_edge) then
                      rgbedge = rgb
                   else
                      rgbedge = r%poly%edge_rgb
                   end if
                end if
             end if
          end if

          ! lazy per-bond H-bond classification cache for this atom: the
          ! classification is translation-invariant, so it is computed once per
          ! unique bond (hbond_class) and reused for every periodic image
          if (dobonds .and. r%bonds%hbond_classify) then
             if (allocated(hbcat)) then
                if (size(hbcat,1) < r%bonds%style%nstar(i)%ncon) deallocate(hbcat)
             end if
             if (.not.allocated(hbcat)) allocate(hbcat(max(r%bonds%style%nstar(i)%ncon,1)))
             hbcat(1:r%bonds%style%nstar(i)%ncon) = -1
          end if

          do i1 = n0(1), n1(1)
             do i2 = n0(2), n1(2)
                do i3 = n0(3), n1(3)
                   ix = (/i1,i2,i3/) + lvec + vacshift
                   if (usetshift) then
                      xx = c%atcel(i)%x - r%sel%tshift
                      ix = ix + nint(xx - floor(xx) + r%sel%tshift - c%atcel(i)%x)
                   end if

                   xx = c%atcel(i)%x + ix
                   xc = c%x2c(xx)

                   ! apply the filter
                   if (havefilter) then
                      res = sys(r%id)%eval(r%sel%filter,errmsg,xc,toklist)
                      if (len_trim(errmsg) == 0) then
                         if (res == 0d0) cycle
                      else
                         havefilter = .false.
                         r%sel%errfilter = errmsg
                      end if
                   end if

                   ! draw the coordination polyhedron for this center image
                   ! (corners are the search result translated to this image)
                   if (dopoly) then
                      do kp = 1, natp
                         xpolyc = c%atcel(eidp(kp))%x + lvecp(:,kp) + ix
                         xvpoly(:,kp) = c%x2c(xpolyc) + uoriginc
                         dvpoly(:,kp) = vibdelta(eidp(kp),lvecp(:,kp)+ix) ! per-corner vibration delta
                      end do
                      call build_polyhedron(xvpoly(:,1:natp),dvpoly(:,1:natp),natp,rgbface,rgbedge,&
                         r%poly%alpha,r%poly%edge_rad,r%poly%coplanar_eps)

                      ! collect this polyhedron's corner atom images to force
                      ! them visible after the main loop (deduplicated there)
                      if (corneractive) then
                         do kp = 1, natp
                            ncorn = ncorn + 1
                            if (ncorn > size(cornlist,2)) call realloc(cornlist,4,2*ncorn)
                            cornlist(1,ncorn) = eidp(kp)
                            cornlist(2:4,ncorn) = lvecp(:,kp) + ix
                         end do
                      end if
                   end if

                   ! animation delta of this (center) atom
                   xdelta1 = vibdelta(i,ix)

                   ! draw the atom. If the atoms are hidden but the bonds shown,
                   ! add an invisible (ghost) pick target at the atom site so
                   ! hidden atoms can still be picked, measured, and box-selected;
                   ! it carries the real idx and is rendered only into the pick
                   ! buffer (skipped in the visible sphere pass).
                   if (r%atoms%display .or. dobonds) then
                      dsph%x = real(xc + uoriginc,c_float)
                      dsph%rgb = rgb
                      dsph%idx(1) = i
                      dsph%idx(2:4) = ix
                      dsph%xdelta = cmplx(xdelta1,kind=c_float_complex)
                      if (r%atoms%display) then
                         dsph%r = real(rad1,c_float)
                         dsph%border = real(r%atoms%border_size,c_float)
                         dsph%rgbborder = r%atoms%border_rgb
                         dsph%ghost = .false.
                      else
                         dsph%r = real(2d0 * r%bonds%rad,c_float) ! generous click radius
                         dsph%border = 0._c_float
                         dsph%rgbborder = rgb
                         dsph%ghost = .true.
                      end if
                      call dl_append(obj%sph,obj%nsph,dsph)
                   end if

                   ! mark this atom image as drawn (for bonds and for
                   ! deduplicating forced polyhedra corner atoms); ix is within
                   ! the lshown bounds by construction (see the mb margin above)
                   if (allocated(lshown)) lshown(i,ix(1),ix(2),ix(3)) = .true.

                   ! bonds
                   if (dobonds) then
                      do ib = 1, r%bonds%style%nstar(i)%ncon
                         ineigh = r%bonds%style%nstar(i)%idcon(ib)
                         if (.not.r%bonds%style%shown(c%atcel(ineigh)%is,c%atcel(i)%is)) cycle
                         ixn = ix + r%bonds%style%nstar(i)%lcon(:,ib)

                         if (r%bonds%imol == 1) then ! intramol
                            if (.not.c%in_same_molecule(i,ix,ineigh,ixn)) cycle
                         elseif (r%bonds%imol == 2) then ! intermol
                            if (c%in_same_molecule(i,ix,ineigh,ixn)) cycle
                         end if

                         if (r%bonds%bothends) then
                            ! skip if the atom has been represented already
                            ! (draws once, and only if both atoms are present)
                            if (.not.lshown(ineigh,ixn(1),ixn(2),ixn(3))) cycle
                         else
                            ! skip if the atom has not been represented already
                            ! (draws once, only one of the atoms need be present)
                            if (lshown(ineigh,ixn(1),ixn(2),ixn(3))) cycle
                         end if

                         ! per-bond order: in "Calculated" mode (4) use the order
                         ! stored in the connectivity (ordcon); otherwise use the fixed order
                         if (r%bonds%order == 4) then
                            iord = r%bonds%style%nstar(i)%ordcon(ib)
                         else
                            iord = r%bonds%order
                         end if

                         ! bond endpoints (Cartesian, bohr)
                         x1 = xc + uoriginc
                         x2 = c%atcel(ineigh)%x + ixn
                         x2 = c%x2c(x2) + uoriginc

                         ! bond color
                         bondrgb = r%bonds%rgb

                         ! Jeffrey-Steiner hydrogen-bond strength classification:
                         ! color each H...A contact by its strength (computed once
                         ! per unique bond in hbond_class, cached across images)
                         if (r%bonds%hbond_classify) then
                            if (hbcat(ib) < 0) hbcat(ib) = hbond_class(i,ib)
                            if (hbcat(ib) == 0) cycle ! not an H-bond
                            bondrgb = r%bonds%hbond_rgb(:,hbcat(ib))
                         end if

                         ! animation delta of the other end
                         xdelta2 = vibdelta(ineigh,ixn)

                         ! fields shared by all cylinders of this bond
                         dcyl%r = real(r%bonds%rad,c_float)
                         dcyl%order = iord
                         dcyl%border = real(r%bonds%border_size,c_float)
                         dcyl%rgbborder = r%bonds%border_rgb
                         dcyl%arvec = real(r%bonds%style%nstar(i)%aromdir(:,ib),c_float)

                         if (r%bonds%color_style == 0 .or. r%bonds%hbond_classify) then
                            ! single cylinder with the bond color
                            dcyl%x1 = real(x1,c_float)
                            dcyl%x1delta = cmplx(xdelta1,kind=c_float_complex)
                            dcyl%x2 = real(x2,c_float)
                            dcyl%x2delta = cmplx(xdelta2,kind=c_float_complex)
                            dcyl%rgb = bondrgb
                            call dl_append(obj%cyl,obj%ncyl,dcyl)
                         else
                            ! two half-cylinders, each colored like its end atom;
                            ! the split point balances the two atomic radii
                            idaux = sysc(r%id)%attype_celatom_to_id(r%atoms%style%type,ineigh)
                            rad2 = r%atoms%style%rad(idaux) * r%mols%style%scale_rad(c%idatcelmol(1,ineigh))
                            dd = norm2(x2 - x1)
                            f1 = min(max((0.5d0 + 0.5d0 * (rad2 - rad1) / dd),0._c_float),1._c_float)
                            f2 = 1._c_float - f1
                            x0 = f1 * x1 + f2 * x2
                            xdelta0 = f1 * xdelta1 + f2 * xdelta2

                            dcyl%x1 = real(x1,c_float)
                            dcyl%x1delta = cmplx(xdelta1,kind=c_float_complex)
                            dcyl%x2 = real(x0,c_float)
                            dcyl%x2delta = cmplx(xdelta0,kind=c_float_complex)
                            dcyl%rgb = rgb
                            call dl_append(obj%cyl,obj%ncyl,dcyl)

                            dcyl%x1 = real(x0,c_float)
                            dcyl%x1delta = cmplx(xdelta0,kind=c_float_complex)
                            dcyl%x2 = real(x2,c_float)
                            dcyl%x2delta = cmplx(xdelta2,kind=c_float_complex)
                            dcyl%rgb = r%atoms%style%rgb(:,idaux) * &
                               r%mols%style%tint_rgb(:,c%idatcelmol(1,ineigh))
                            call dl_append(obj%cyl,obj%ncyl,dcyl)
                         end if
                      end do ! ncon
                   end if

                   if (r%labels%display) then
                      select case(r%labels%type)
                      case (0,5,6)
                         idl = c%atcel(i)%is
                      case (2,3)
                         idl = i
                      case (1,4,8)
                         idl = c%atcel(i)%idx
                      case (7)
                         idl = c%idatcelmol(1,i)
                      end select

                      ! labels
                      if (r%labels%style%shown(idl)) then
                         dstr%x = real(xc + uoriginc,c_float)
                         dstr%xdelta = cmplx(xdelta1,kind=c_float_complex)
                         dstr%r = real(rad1,c_float)
                         dstr%rgb = r%labels%rgb
                         if (r%labels%const_size) then
                            dstr%scale = real(r%labels%scale,c_float)
                         else
                            dstr%scale = real(-r%labels%scale,c_float)
                         end if
                         dstr%offset = real(r%labels%offset,c_float)
                         dstr%str = trim(r%labels%style%str(idl))
                         if (r%labels%type == 3) then
                            ! add the lattice vectors
                            dstr%str = dstr%str // "[" //&
                               string(ix(1)) // "," // string(ix(2)) // "," //string(ix(3)) // "]"
                         end if
                         call dl_append(obj%string,obj%nstring,dstr)
                      end if ! label display conditions
                   end if ! label_display
                end do ! i3
             end do ! i2
          end do ! i1
       end do ! loop over complete atom list (i)

       ! draw the polyhedra corner atoms that the selection did not already draw
       ! (so every drawn polyhedron shows its corner atoms)
       if (corneractive) then
          do ica = 1, ncorn
             ix = cornlist(2:4,ica)
             call check_lshown(cornlist(1,ica),ix(1),ix(2),ix(3))
             if (lshown(cornlist(1,ica),ix(1),ix(2),ix(3))) cycle ! already drawn
             lshown(cornlist(1,ica),ix(1),ix(2),ix(3)) = .true.

             ! style and position of the corner atom
             idc = sysc(r%id)%attype_celatom_to_id(r%atoms%style%type,cornlist(1,ica))
             imolc = c%idatcelmol(1,cornlist(1,ica))
             rgb = r%atoms%style%rgb(:,idc) * r%mols%style%tint_rgb(:,imolc)
             rad1 = r%atoms%style%rad(idc) * r%mols%style%scale_rad(imolc)
             xx = c%atcel(cornlist(1,ica))%x + ix
             xc = c%x2c(xx)

             ! animation delta of the corner atom (so it moves with the polyhedron)
             xdelta1 = vibdelta(cornlist(1,ica),ix)

             dsph%x = real(xc + uoriginc,c_float)
             dsph%r = real(rad1,c_float)
             dsph%rgb = rgb
             dsph%idx(1) = cornlist(1,ica)
             dsph%idx(2:4) = ix
             dsph%xdelta = cmplx(xdelta1,kind=c_float_complex)
             dsph%border = real(r%atoms%border_size,c_float)
             dsph%rgbborder = r%atoms%border_rgb
             dsph%ghost = .false.
             call dl_append(obj%sph,obj%nsph,dsph)
          end do
       end if
       if (allocated(cornlist)) deallocate(cornlist)
       if (allocated(up2dsp)) deallocate(up2dsp)
    elseif (r%type == reptype_unitcell) then
       !!! unit cell representation !!!

       ! number of cells
       n = 1
       if (r%sel%pertype == 1) then
          n = nc
       elseif (r%sel%pertype == 2) then
          n = r%sel%ncell
       end if

       ! vacuum directions: we have a vacuum and only one cell in that direction
       isvac = .false.
       if (r%uc%vaccutsticks) then
          do i = 1, 3
             if (c%vaclength(i) > iperiod_vacthr .and. n(i) == 1) then
                isvac(i) = .true.
             end if
          end do
       end if

       ! external cell
       do i = 1, 12
          ix1 = uc(:,1,i) * n
          ix2 = uc(:,2,i) * n
          call process_vacuum_uc_sticks(ix1,ix2,x1,x2,docycle)
          if (docycle) cycle

          ! add the sticks
          dcyl%x1 = real(x1,c_float)
          dcyl%x2 = real(x2,c_float)
          dcyl%r = real(r%uc%radius,c_float)
          if (r%uc%coloraxes.and.i>=1.and.i<=3) then
             dcyl%rgb = ColorAxes_def(:,i)
          else
             dcyl%rgb = r%uc%rgb
          end if
          call dl_append(obj%cylflat,obj%ncylflat,dcyl)
       end do

       ! draw inner cylinders
       if (r%uc%inner) then
          do i1 = 0, n(1)-1
             do i2 = 0, n(2)-1
                do i3 = 0, n(3)-1
                   do i = 1, 12
                      ix1 = (/i1,i2,i3/) + uc(:,1,i)
                      ix2 = (/i1,i2,i3/) + uc(:,2,i)

                      ! skip outer cylinders
                      isedge = (ix1 == 0 .and. ix2 == 0) .or. (ix1 == n .and. ix2 == n)
                      isedge(ucdir(i)) = .true.
                      if (all(isedge)) cycle

                      ! process vacuum
                      call process_vacuum_uc_sticks(ix1,ix2,x1,x2,docycle)
                      if (docycle) cycle

                      ! stippled lines for the inner lines
                      dcyl%r = real(r%uc%radiusinner,c_float)
                      dcyl%rgb = r%uc%rgb
                      if (r%uc%innerstipple) then
                         nstep = ceiling(norm2(x2 - x1) / r%uc%innersteplen)
                         do j = 1, nstep
                            dcyl%x1 = real(x1 + real(2*j-1,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            dcyl%x2 = real(x1 + real(2*j,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            call dl_append(obj%cylflat,obj%ncylflat,dcyl)
                         end do
                      else
                         dcyl%x1 = real(x1 ,c_float)
                         dcyl%x2 = real(x2 ,c_float)
                         call dl_append(obj%cylflat,obj%ncylflat,dcyl)
                      end if
                   end do
                end do
             end do
          end do
       end if
    elseif (r%type == reptype_axes) then
       !!! cartesian axes representation !!!

       ! placement: at the cartesian origin (+shift), or anchored to a
       ! fixed window position. In the latter case the geometry is built
       ! around the local origin and positioned at render time; the window
       ! position is stamped onto each gizmo draw item below.
       fixed = (r%axes%placement == 1)
       if (fixed) then
          uoriginc = 0d0
          ! stamp the window placement onto the gizmo templates so every
          ! appended gizmo item carries its own position/scale flag
          dcylgiz%winpos = real(r%axes%winpos,c_float)
          dcylgiz%scalewithzoom = r%axes%scalewithzoom
          dstrgiz%winpos = real(r%axes%winpos,c_float)
          dstrgiz%scalewithzoom = r%axes%scalewithzoom
       else
          ! origin in the requested coordinate system, converted to
          ! cartesian (bohr)
          if (r%axes%coordtype == 2) then
             uoriginc = r%axes%origin ! cartesian (bohr)
          elseif (r%axes%coordtype == 0 .and. .not.c%ismolecule) then
             uoriginc = c%x2c(r%axes%origin) ! crystallographic
          else
             uoriginc = r%axes%origin / bohrtoa ! cartesian (angstrom)
          end if
          ! for molecules, cartesian coordinates are referred to the molecular center
          if (c%ismolecule) uoriginc = uoriginc - c%molx0
       end if

       ! global scale factor applied to the whole gizmo (arrows and labels)
       axsc = r%axes%scale

       ! arrowhead geometry (head length capped so it never exceeds the
       ! total axis length)
       rad1 = min(r%axes%conelength, r%axes%length) * axsc ! head length
       rad2 = r%axes%coneradius * axsc ! head radius

       do k = 1, 3
          ! unit direction for this axis: cartesian (lab-frame) or along
          ! the crystallographic lattice vector
          if (r%axes%kind == 1 .and. .not.c%ismolecule) then
             x0 = c%m_x2c(:,k)
             x0 = x0 / norm2(x0)
          else
             x0 = 0d0
             x0(k) = 1d0
          end if
          ! reorient the axis directions (identity unless a frame was requested)
          x0 = matmul(r%axes%rot,x0)

          ! shaft (round, lit cylinder)
          x1 = uoriginc
          x2 = uoriginc + max(r%axes%length * axsc - rad1,0d0) * x0
          dcyl%x1 = real(x1,c_float)
          dcyl%x2 = real(x2,c_float)
          dcyl%x1delta = cmplx(0d0,0d0,kind=c_float_complex)
          dcyl%x2delta = cmplx(0d0,0d0,kind=c_float_complex)
          dcyl%r = real(r%axes%radius * axsc,c_float)
          dcyl%rgb = r%axes%rgb(:,k)
          dcyl%order = 1
          dcyl%border = 0._c_float
          dcyl%rgbborder = 0._c_float
          if (fixed) then
             dcylgiz%dl_cylinder = dcyl
             call dl_append(obj%cylgiz,obj%ncylgiz,dcylgiz)
          else
             call dl_append(obj%cyl,obj%ncyl,dcyl)
          end if

          ! arrowhead (cone) from the shaft end to the axis tip
          dcyl%x1 = real(x2,c_float)
          dcyl%x2 = real(uoriginc + r%axes%length * axsc * x0,c_float)
          dcyl%r = real(rad2,c_float)
          if (fixed) then
             dcylgiz%dl_cylinder = dcyl
             call dl_append(obj%conegiz,obj%nconegiz,dcylgiz)
          else
             call dl_append(obj%cone,obj%ncone,dcyl)
          end if

          ! label: along the axis from the arrowhead tip by the shared
          ! distance, plus the per-axis cartesian offset
          if (r%axes%showlabels) then
             dstr%x = real(uoriginc + (r%axes%length + r%axes%labeldistance) * axsc * x0 + &
                r%axes%labeloffset(:,k) * axsc,c_float)
             dstr%xdelta = cmplx(0d0,0d0,kind=c_float_complex)
             dstr%r = real(rad2,c_float)
             dstr%rgb = r%axes%labelrgb
             if (r%axes%labelconstsize) then
                dstr%scale = real(r%axes%labelscale * axsc,c_float)
             else
                dstr%scale = real(-r%axes%labelscale * axsc,c_float)
             end if
             dstr%offset = 0._c_float
             dstr%str = trim(r%axes%labelstr(k))
             if (fixed) then
                dstrgiz%dl_string = dstr
                call dl_append(obj%stringgiz,obj%nstringgiz,dstrgiz)
             else
                call dl_append(obj%string,obj%nstring,dstr)
             end if
          end if
       end do
    elseif (r%type == reptype_rotaxis) then
       !!! rotation-axis representation (single black cylinder through the origin) !!!

       ! origin in cartesian (bohr); for molecules referred to the molecular center
       uoriginc = r%rotaxis%origin
       if (c%ismolecule) uoriginc = uoriginc - c%molx0

       ! double-ended cylinder (the axis line): origin +/- length*dir
       x0 = r%rotaxis%dir / max(norm2(r%rotaxis%dir),1d-10)
       dcyl%x1 = real(uoriginc - r%rotaxis%length * x0,c_float)
       dcyl%x2 = real(uoriginc + r%rotaxis%length * x0,c_float)
       dcyl%x1delta = cmplx(0d0,0d0,kind=c_float_complex)
       dcyl%x2delta = cmplx(0d0,0d0,kind=c_float_complex)
       dcyl%r = real(r%rotaxis%radius,c_float)
       dcyl%rgb = r%rotaxis%rgb
       dcyl%order = 1
       dcyl%border = 0._c_float
       dcyl%rgbborder = 0._c_float
       call dl_append(obj%cyl,obj%ncyl,dcyl)
    elseif (r%type == reptype_symelem) then
       !!! symmetry element(s) (plane/axis) !!!
       if (r%symelem%style%isinit) then
          ! persistent set
          if (r%symelem%coordtype == 2) then
             uoriginc = r%symelem%origin ! cartesian (bohr)
          elseif (r%symelem%coordtype == 0 .and. .not.c%ismolecule) then
             uoriginc = c%x2c(r%symelem%origin) ! crystallographic
          else
             uoriginc = r%symelem%origin / bohrtoa ! cartesian (angstrom)
          end if
          if (c%ismolecule) uoriginc = uoriginc - c%molx0
          do i1 = 1, r%symelem%style%nop
             if (.not.r%symelem%style%shown(i1)) cycle
             if (r%symelem%style%kind(i1) == 0) cycle
             call draw_symmetry_element(r%symelem%style%kind(i1),r%symelem%style%dir(:,i1),&
                r%symelem%style%order(i1),uoriginc,r%symelem%usecustomrgb,r%symelem%rgb)
          end do
       else
          ! transient single element (hover/selection preview)
          uoriginc = r%symelem%origin_transient
          if (c%ismolecule) uoriginc = uoriginc - c%molx0
          if (r%symelem%kind /= 0) &
             call draw_symmetry_element(r%symelem%kind,r%symelem%dir,r%symelem%order,&
                uoriginc,.false.,r%symelem%rgb)
       end if
    end if ! reptype
  contains

    !> Jeffrey-Steiner hydrogen-bond strength classification of bond ib of cell
    !> atom i. Returns 0 if the contact is not an H-bond (skip it), or the
    !> strength class 1..3 (strong/moderate/weak). The geometry is
    !> translation-invariant, so this runs once per unique bond (zero-cell
    !> image) and the result is cached across periodic images.
    !> G. A. Jeffrey, An Introduction to Hydrogen Bonding, Oxford University Press, 1997
    !> T. Steiner, Angew. Chem. Intl. Ed. 41 (2002) 48, doi:10.1002/1521-3773(20020104)41:1<48::AID-ANIE48>3.0.CO;2-U
    function hbond_class(i,ib) result(cat)
      integer, intent(in) :: i, ib
      integer :: cat

      integer :: ineigh, hbzi, hbzn, hbhcel, hblv(3), hbacel, hbalv(3)
      integer :: hbdon, hbdc, hbac, hbib
      real*8 :: xhb_h(3), xhb_a(3), xhb_d(3), hbang, dd

      cat = 0
      ineigh = r%bonds%style%nstar(i)%idcon(ib)

      ! identify the H end (the other end is the acceptor); work in fractional
      ! coordinates with atom i in the zero cell
      hbzi = c%spc(c%atcel(i)%is)%z
      hbzn = c%spc(c%atcel(ineigh)%is)%z
      if (hbzi == 1) then
         hbhcel = i
         hblv = 0
         hbacel = ineigh
         hbalv = r%bonds%style%nstar(i)%lcon(:,ib)
      elseif (hbzn == 1) then
         hbhcel = ineigh
         hblv = r%bonds%style%nstar(i)%lcon(:,ib)
         hbacel = i
         hbalv = 0
      else
         return ! neither end is hydrogen: not an H-bond
      end if
      xhb_h = c%atcel(hbhcel)%x + hblv
      xhb_a = c%atcel(hbacel)%x + hbalv

      ! find the donor atom D in D-H...A
      hbdon = 0
      if (allocated(c%nstar)) then
         do hbib = 1, c%nstar(hbhcel)%ncon
            if (c%spc(c%atcel(&
               c%nstar(hbhcel)%idcon(hbib))%is)%z > 1) then
               hbdon = hbib
               exit
            end if
         end do
      end if
      if (hbdon == 0) return ! no donor atom: skip

      ! H...A distance class
      dd = c%distance(xhb_a,xhb_h)
      if (dd < r%bonds%hbond_dist(1)) then
         hbdc = 1
      elseif (dd < r%bonds%hbond_dist(2)) then
         hbdc = 2
      else
         hbdc = 3
      end if

      ! D-H...A angle class (the weaker of the two wins)
      xhb_d = c%atcel(c%nstar(hbhcel)%idcon(hbdon))%x +&
         hblv + c%nstar(hbhcel)%lcon(:,hbdon)
      hbang = c%angle(xhb_d,xhb_h,xhb_a) * 180d0 / pi
      if (hbang < hbond_angmin_def) return ! too bent to be a H-bond
      if (hbang >= r%bonds%hbond_ang(2)) then
         hbac = 1
      elseif (hbang >= r%bonds%hbond_ang(1)) then
         hbac = 2
      else
         hbac = 3
      end if
      cat = max(hbdc,hbac)

    end function hbond_class

    !> Append an opaque thin edge cylinder between cartesian points p and q
    !> (uses the host dcyl template for radius/color).
    subroutine append_edge(p,q)
      real*8, intent(in) :: p(3), q(3)

      dcyl%x1 = real(p,c_float)
      dcyl%x2 = real(q,c_float)
      call dl_append(obj%cyl,obj%ncyl,dcyl)

    end subroutine append_edge

    !> Draw one symmetry element of kind skind (symop_kind_plane/axis), with
    !> unit direction/normal sdir, rotation order sorder, passing through the
    !> cartesian-bohr point uoriginc. If usecustom, everything is drawn in
    !> customrgb; otherwise planes use the default color and axes are colored by
    !> rotation order. Uses the host r%symelem%size/r%symelem%cen for sizing.
    subroutine draw_symmetry_element(skind,sdir,sorder,uoriginc,usecustom,customrgb)
      integer, intent(in) :: skind, sorder
      real*8, intent(in) :: sdir(3), uoriginc(3)
      logical, intent(in) :: usecustom
      real(c_float), intent(in) :: customrgb(3)

      real*8 :: lx0(3), lxx(3), lx1(3), lx2(3), lxc(3), le1v(3), le2v(3), lres
      real(c_float) :: rgbel(3)
      integer :: j1, j2, j3, m1(3)
      type(dl_plane) :: dpl

      ! unit direction: the plane normal or the axis direction (cartesian)
      lx0 = sdir / max(norm2(sdir),1d-10)

      ! opaque thin-cylinder template (plane frame edges, axis shafts)
      dcyl%x1delta = cmplx(0d0,0d0,kind=c_float_complex)
      dcyl%x2delta = cmplx(0d0,0d0,kind=c_float_complex)
      dcyl%r = real(symelem_frame_radius,c_float)
      dcyl%alpha = 1._c_float
      dcyl%order = 1
      dcyl%border = 0._c_float
      dcyl%rgbborder = 0._c_float

      if (skind == symop_kind_plane) then
         ! mirror/glide plane: translucent fill + opaque border frame
         rgbel = symelem_rgb_def
         if (usecustom) rgbel = customrgb
         dcyl%rgb = rgbel

         ! in-plane orthonormal basis perpendicular to the plane normal lx0
         if (abs(lx0(1)) < 0.9d0) then
            lxx = (/1d0,0d0,0d0/)
         else
            lxx = (/0d0,1d0,0d0/)
         end if
         lx1 = cross(lx0,lxx)
         lx1 = lx1 / norm2(lx1)
         lx2 = cross(lx0,lx1) ! unit, (lx0,lx1,lx2) orthonormal

         ! rectangle center = projection of the system center onto the plane
         lxc = r%symelem%cen - dot_product(r%symelem%cen - uoriginc,lx0) * lx0
         lres = symelem_margin * r%symelem%size
         le1v = lres * lx1
         le2v = lres * lx2

         ! translucent fill
         dpl%x = real(lxc,c_float)
         dpl%e1 = real(le1v,c_float)
         dpl%e2 = real(le2v,c_float)
         dpl%rgb = rgbel
         dpl%alpha = symelem_alpha
         call dl_append(obj%plane,obj%nplane,dpl)

         ! opaque border frame (4 edge cylinders)
         call append_edge(lxc - le1v - le2v, lxc + le1v - le2v)
         call append_edge(lxc + le1v - le2v, lxc + le1v + le2v)
         call append_edge(lxc + le1v + le2v, lxc - le1v + le2v)
         call append_edge(lxc - le1v + le2v, lxc - le1v - le2v)
      else
         ! rotation/rotoinversion axis: a thick opaque shaft (colored by the
         ! rotation order, unless a custom color is set) through every visible
         ! lattice point (crystals) or the molecular center (molecules)
         rgbel = symelem_rgb_def
         if (usecustom) then
            rgbel = customrgb
         elseif (sorder >= lbound(symelem_rgb_order,2) .and. sorder <= ubound(symelem_rgb_order,2)) then
            if (any(symelem_rgb_order(:,sorder) /= 0._c_float)) rgbel = symelem_rgb_order(:,sorder)
         end if
         dcyl%rgb = rgbel
         dcyl%r = real(symelem_axis_radius,c_float)

         lres = symelem_margin * r%symelem%size
         m1 = 0
         if (.not.c%ismolecule) m1 = nc
         do j1 = 0, m1(1)
            do j2 = 0, m1(2)
               do j3 = 0, m1(3)
                  lxc = uoriginc + c%x2c(real((/j1,j2,j3/),8))
                  lxc = lxc + dot_product(r%symelem%cen - lxc,lx0) * lx0 ! foot of the center on the axis
                  dcyl%x1 = real(lxc - lres * lx0,c_float)
                  dcyl%x2 = real(lxc + lres * lx0,c_float)
                  call dl_append(obj%cyl,obj%ncyl,dcyl)
               end do
            end do
         end do
      end if

    end subroutine draw_symmetry_element

    !> Build a coordination polyhedron from nvv vertex positions xv
    !> (cartesian, bohr). The vertex centroid is used as the interior
    !> reference point. Adds translucent triangular faces (color rgbf,
    !> opacity alphaf) and opaque edge cylinders (color rgbe, radius
    !> rade) to the draw lists. If the vertices are coplanar to within
    !> eps, a filled polygon is drawn instead of a 3D convex hull.
    subroutine build_polyhedron(xv,dv,nvv,rgbf,rgbe,alphaf,rade,eps)
      use iso_c_binding, only: c_ptr, c_int, c_float_complex
      integer, intent(in) :: nvv
      real*8, intent(in) :: xv(3,nvv)
      complex*16, intent(in) :: dv(3,nvv) ! per-vertex vibration deltas
      real(c_float), intent(in) :: rgbf(3), rgbe(3)
      real*8, intent(in) :: alphaf, rade, eps

      real*8, parameter :: edge_coplanar_cos = 0.9986d0 ! ~3 degrees

      integer :: a, b, k, kk, ntri, nedge, ip, iq, e
      integer, allocatable :: itri(:,:), edgei(:,:), iord_(:)
      real*8, allocatable :: edgenrm(:,:), ang(:)
      logical, allocatable :: edgekeep(:)
      real*8 :: cen0(3), nrm(3), dev, e1u(3), e2u(3), cc, nrm_a(3)
      complex*16 :: dcen(3)
      type(c_ptr) :: ctx
      integer(c_int) :: ier
      integer :: nf
      type(dl_cylinder) :: dedge

      ! centroid (an interior reference point), best-fit plane unit normal, and
      ! the maximum out-of-plane deviation
      call plane_from_points(xv,nvv,cen0,nrm,dev)

      if (dev < eps .and. norm2(nrm) > 1d-10) then
         ! planar polygon: order vertices by angle about the normal and
         ! fan-triangulate from the centroid
         e1u = xv(:,1) - cen0
         e1u = e1u - dot_product(e1u,nrm)*nrm
         if (norm2(e1u) < 1d-10) return
         e1u = e1u / norm2(e1u)
         e2u = cross(nrm,e1u)

         allocate(ang(nvv),iord_(nvv))
         do k = 1, nvv
            ang(k) = atan2(dot_product(xv(:,k)-cen0,e2u),dot_product(xv(:,k)-cen0,e1u))
            iord_(k) = k
         end do
         call mergesort(ang,iord_,1,nvv)

         ! the fan apex sits at the centroid; animate it by the mean vertex delta
         dcen = 0d0
         do k = 1, nvv
            dcen = dcen + dv(:,k)
         end do
         dcen = dcen / nvv

         dedge%r = real(rade,c_float)
         dedge%rgb = rgbe
         do k = 1, nvv
            kk = mod(k,nvv) + 1
            call append_triangle(cen0,xv(:,iord_(k)),xv(:,iord_(kk)),&
               dcen,dv(:,iord_(k)),dv(:,iord_(kk)),rgbf,alphaf)
            dedge%x1 = real(xv(:,iord_(k)),c_float)
            dedge%x2 = real(xv(:,iord_(kk)),c_float)
            dedge%x1delta = cmplx(dv(:,iord_(k)),kind=c_float_complex)
            dedge%x2delta = cmplx(dv(:,iord_(kk)),kind=c_float_complex)
            call dl_append(obj%cylflat,obj%ncylflat,dedge)
         end do
         return
      end if

      ! 3D convex hull of the vertices, seen from the interior centroid
      call runqhull_basintriangulate_step1(nvv,cen0,xv,nf,ctx,ier)
      if (ier /= 0 .or. nf <= 0) return
      ntri = nf
      allocate(itri(3,ntri))
      call runqhull_basintriangulate_step2(ntri,itri,ctx)

      ! faces
      do a = 1, ntri
         call append_triangle(xv(:,itri(1,a)),xv(:,itri(2,a)),xv(:,itri(3,a)),&
            dv(:,itri(1,a)),dv(:,itri(2,a)),dv(:,itri(3,a)),rgbf,alphaf)
      end do

      ! edges: collect unique triangle edges and the adjacent-face normals;
      ! suppress edges shared by two near-coplanar triangles (triangulation
      ! diagonals across a flat polyhedron face)
      allocate(edgei(2,3*ntri),edgenrm(3,3*ntri),edgekeep(3*ntri))
      nedge = 0
      do a = 1, ntri
         nrm_a = cross(xv(:,itri(2,a))-xv(:,itri(1,a)),xv(:,itri(3,a))-xv(:,itri(1,a)))
         if (norm2(nrm_a) > 1d-10) nrm_a = nrm_a / norm2(nrm_a)
         do k = 1, 3
            ip = itri(k,a)
            iq = itri(mod(k,3)+1,a)
            if (ip > iq) then
               b = ip
               ip = iq
               iq = b
            end if
            e = 0
            do b = 1, nedge
               if (edgei(1,b) == ip .and. edgei(2,b) == iq) then
                  e = b
                  exit
               end if
            end do
            if (e == 0) then
               nedge = nedge + 1
               edgei(:,nedge) = (/ip,iq/)
               edgenrm(:,nedge) = nrm_a
               edgekeep(nedge) = .true.
            else
               cc = abs(dot_product(edgenrm(:,e),nrm_a))
               edgekeep(e) = (cc < edge_coplanar_cos)
            end if
         end do
      end do
      dedge%r = real(rade,c_float)
      dedge%rgb = rgbe
      do e = 1, nedge
         if (.not.edgekeep(e)) cycle
         dedge%x1 = real(xv(:,edgei(1,e)),c_float)
         dedge%x2 = real(xv(:,edgei(2,e)),c_float)
         dedge%x1delta = cmplx(dv(:,edgei(1,e)),kind=c_float_complex)
         dedge%x2delta = cmplx(dv(:,edgei(2,e)),kind=c_float_complex)
         call dl_append(obj%cylflat,obj%ncylflat,dedge)
      end do

    end subroutine build_polyhedron

    !> Append a translucent triangular face (vertices p1,p2,p3, cartesian bohr)
    !> with per-vertex vibration deltas d1,d2,d3.
    subroutine append_triangle(p1,p2,p3,d1,d2,d3,rgb_,alpha_)
      use iso_c_binding, only: c_float_complex
      real*8, intent(in) :: p1(3), p2(3), p3(3)
      complex*16, intent(in) :: d1(3), d2(3), d3(3)
      real(c_float), intent(in) :: rgb_(3)
      real*8, intent(in) :: alpha_
      type(dl_triangle) :: dtri

      dtri%x1 = real(p1,c_float)
      dtri%x2 = real(p2,c_float)
      dtri%x3 = real(p3,c_float)
      dtri%rgb = rgb_
      dtri%alpha = real(alpha_,c_float)
      dtri%x1delta = cmplx(d1,kind=c_float_complex)
      dtri%x2delta = cmplx(d2,kind=c_float_complex)
      dtri%x3delta = cmplx(d3,kind=c_float_complex)
      call dl_append(obj%triangle,obj%ntriangle,dtri)

    end subroutine append_triangle

    !> Vibration-animation displacement of the periodic image of cell atom iat
    !> at lattice translation ix; zero if the scene is not animating. Uses the
    !> phasors precomputed in vibbase (mass, mode vector and atom-position
    !> phase); only the lattice-translation phase depends on the image.
    function vibdelta(iat,ix) result(dv)
      integer, intent(in) :: iat
      integer, intent(in) :: ix(3)
      complex*16 :: dv(3)

      dv = 0d0
      if (.not.doanim_) return
      dv = vibbase(:,iat) * exp(img * tpi * dot_product(real(ix,8),c%vib%qpt(:,iqpt)))

    end function vibdelta

    subroutine process_vacuum_uc_sticks(ix1,ix2,x1,x2,docycle)
      integer, intent(in) :: ix1(3), ix2(3)
      real*8, intent(out) :: x1(3), x2(3)
      logical, intent(out) :: docycle

      real*8 :: ucini(3), ucend(3)

      ! initialize
      docycle = .false.

      ! check whether this stick goes in the vacuum direction
      isvacdir = any((ix1 /= ix2) .and. isvac .and..not.all(isvac))

      ucini = real(ix1,8)
      ucend = real(ix2,8)
      if (isvacdir) then
         ! there is vacuum and this stick goes in the vacuum direction: cut it
         do j = 1, 3
            if (isvac(j)) then
               ucini(j) = c%vactop(j) - 1d0 - vacextension / c%aa(j)
               ucend(j) = c%vacbot(j) + vacextension / c%aa(j)
            end if
         end do
      elseif (any(isvac).and..not.all(isvac)) then
         ! there is vacuum and this stick does not go in the vacuum direction:
         ! remove it if "uc" has a 1 in that coordinate; shift it if it has 0 coordinate
         do j = 1, 3
            if (isvac(j)) then
               if (ix1(j) /= 0) then
                  docycle = .true.
                  return
               else
                  ucini(j) = 0.5d0 * (c%vactop(j) - 1d0 + c%vacbot(j))
                  ucend(j) = 0.5d0 * (c%vactop(j) - 1d0 + c%vacbot(j))
               end if
            end if
         end do
      end if

      ! stick ends
      x1 = ucini + r%sel%origin
      x1 = c%x2c(x1)
      x2 = ucend + r%sel%origin
      x2 = c%x2c(x2)

    end subroutine process_vacuum_uc_sticks

    subroutine check_lshown(i,i1,i2,i3)
      integer, intent(in) :: i, i1, i2, i3

      integer :: l, l1, l2, l3, u, u1, u2, u3
      logical, allocatable :: lshown_aux(:,:,:,:)

      if (i < lbound(lshown,1) .or. i > ubound(lshown,1) .or.&
         i1 < lbound(lshown,2) .or. i1 > ubound(lshown,2) .or.&
         i2 < lbound(lshown,3) .or. i2 > ubound(lshown,3) .or.&
         i3 < lbound(lshown,4) .or. i3 > ubound(lshown,4)) then
         l = min(i,lbound(lshown,1))
         u = max(i,ubound(lshown,1))
         l1 = min(i1,lbound(lshown,2))
         u1 = max(i1,ubound(lshown,2))
         l2 = min(i2,lbound(lshown,3))
         u2 = max(i2,ubound(lshown,3))
         l3 = min(i3,lbound(lshown,4))
         u3 = max(i3,ubound(lshown,4))
         allocate(lshown_aux(l:u,l1:u1,l2:u2,l3:u3))
         lshown_aux = .false.
         lshown_aux(lbound(lshown,1):ubound(lshown,1),lbound(lshown,2):ubound(lshown,2),&
            lbound(lshown,3):ubound(lshown,3),lbound(lshown,4):ubound(lshown,4)) = &
            lshown
         call move_alloc(lshown_aux,lshown)
      end if

    end subroutine check_lshown

  end subroutine add_draw_elements

  !> Reset all styles. Reset if itype = 0 (all), 1 (atom), 2 (bonds),
  !> 3 (labels), 4 (mol).
  module subroutine reset_all_styles(r,itype)
    class(representation), intent(inout) :: r
    integer, intent(in) :: itype

    if (itype == 0 .or. itype == 1) &
       call r%atoms%style%reset(r)
    if (itype == 0 .or. itype == 2) &
       call r%bonds%style%reset(r)
    if (itype == 0 .or. itype == 3) &
       call r%labels%style%reset(r)
    if (itype == 0 .or. itype == 4) &
       call r%mols%style%reset(r)
    if (itype == 0 .or. itype == 8) &
       call r%poly%style%reset(r)

  end subroutine reset_all_styles

  !> Reset atom style to the parameters and the contents of the system
  !> point at by representation r. Uses d%type to fill the arrays.
  module subroutine atom_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sysc, sys_ready, ok_system, atlisttype_species
    use gui_main, only: ColorElement
    use param, only: atmcov, atmvdw, jmlcol, jmlcol2, maxzat, maxzat0
    class(atom_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, ispc, iz

    ! if not initialized, set type
    if (.not.d%isinit) d%type = atlisttype_species

    ! set the atom style to zero
    d%ntype = 0
    d%isinit = .false.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%rgb)) deallocate(d%rgb)
    if (allocated(d%rad)) deallocate(d%rad)

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! fill data
    d%ntype = sysc(r%id)%attype_number(d%type)
    allocate(d%shown(d%ntype),d%rgb(3,d%ntype),d%rad(d%ntype))
    d%shown = .true.
    do i = 1, d%ntype
       ispc = sysc(r%id)%attype_species(d%type,i)
       iz = sys(r%id)%c%spc(ispc)%z

       ! color scheme or GUI colors
       if (r%atoms%color_type == 0) then
          d%rgb(:,i) = ColorElement(:,iz)
       elseif (r%atoms%color_type == 1) then
          d%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
       else
          d%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
       end if

       ! scale covalent, vdw, or absolute value
       if (r%atoms%radii_type == 0) then
          d%rad(i) = r%atoms%radii_scale * atmcov(iz)
       elseif (r%atoms%radii_type == 1) then
          d%rad(i) = r%atoms%radii_scale * atmvdw(iz)
       else
          d%rad(i) = r%atoms%radii_value
       endif

       if (r%flavor == repflavor_atoms_criticalpoints) then
          ! show only the critical point atoms
          if (iz <= maxzat .or. iz == maxzat0) d%shown(i) = .false.
       elseif (r%flavor == repflavor_atoms_gradientpaths) then
          ! show only the gradient paths
          if (iz /= maxzat0) d%shown(i) = .false.
       else
          ! do not shown the critical point atoms
          if (iz > maxzat) d%shown(i) = .false.
       end if
    end do
    d%isinit = .true.

  end subroutine atom_style_reset

  !> Reset colors in an atom style to defaults.
  module subroutine atom_style_reset_colors(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sysc, sys_ready, ok_system
    use gui_main, only: ColorElement
    class(atom_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, ispc, iz

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    do i = 1, d%ntype
       ispc = sysc(r%id)%attype_species(d%type,i)
       iz = sys(r%id)%c%spc(ispc)%z
       d%rgb(:,i) = ColorElement(:,iz)
    end do

  end subroutine atom_style_reset_colors

  !> Deallocate all arrays and end the atom syle.
  module subroutine atom_style_end(d)
    class(atom_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%rgb)) deallocate(d%rgb)
    if (allocated(d%rad)) deallocate(d%rad)

  end subroutine atom_style_end

  !> Reset molecule style with default values. Use the information in
  !> representation r, or leave it empty if system is uninitalized.
  module subroutine mol_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    class(mol_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i

    ! set the atom style to zero
    d%ntype = 0
    d%isinit = .false.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%tint_rgb)) deallocate(d%tint_rgb)
    if (allocated(d%scale_rad)) deallocate(d%scale_rad)

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! fill
    d%ntype = sys(r%id)%c%nmol
    allocate(d%shown(d%ntype),d%tint_rgb(3,d%ntype))
    allocate(d%scale_rad(d%ntype))
    do i = 1, sys(r%id)%c%nmol
       d%tint_rgb(:,i) = 1._c_float
       d%scale_rad(i) = 1d0
    end do
    d%shown = .true.
    d%isinit = .true.

  end subroutine mol_style_reset

  !> Deallocate all arrays and end the mol syle.
  module subroutine mol_style_end(d)
    class(mol_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%tint_rgb)) deallocate(d%tint_rgb)
    if (allocated(d%scale_rad)) deallocate(d%scale_rad)

  end subroutine mol_style_end

  !> Generate the neighbor stars from the data in the rij table using
  !> the geometry in system isys.
  module subroutine generate_neighstars(d,r)
    use systems, only: sys, sys_ready, ok_system
    class(bond_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    ! check all the info is available
    if (.not.d%isinit) return
    if (.not.ok_system(r%id,sys_ready)) return

    ! generate the new neighbor star using this representation's bonding
    ! parameters: per-species covalent radii and bond factor for non-metal
    ! bonds, bond delta for metal bonds (same criteria as the geometry window)
    call sys(r%id)%c%find_asterisms(d%nstar,atmrad=r%bonds%atmrad,bondfac=r%bonds%bfactor,&
       bonddelta=r%bonds%bdelta)

  end subroutine generate_neighstars

  !> Copy the neighbor stars from the given system.
  module subroutine copy_neighstars_from_system(d,isys)
    use systems, only: sys, sys_ready, ok_system
    class(bond_geom_style), intent(inout) :: d
    integer, intent(in) :: isys

    ! check all the info is available
    if (.not.d%isinit) return
    if (.not.ok_system(isys,sys_ready)) return

    ! copy the nstar
    if (allocated(sys(isys)%c%nstar)) &
       d%nstar = sys(isys)%c%nstar

  end subroutine copy_neighstars_from_system

  !> Reset bond style with default values, according to the given
  !> representation.
  module subroutine bond_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    class(bond_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, j, iz

    ! clear the bond style
    d%isinit = .false.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%nstar)) deallocate(d%nstar)
    d%use_sys_nstar = .true.

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return
    d%isinit = .true.

    ! fill temp options
    allocate(d%shown(sys(r%id)%c%nspc,sys(r%id)%c%nspc))
    d%shown = .true.

    ! fill data according to flavor
    if (r%flavor == repflavor_atoms_vdwcontacts) then
       ! van der waals contacts
       d%use_sys_nstar = .false.
       do i = 1, sys(r%id)%c%nspc
          if (sys(r%id)%c%spc(i)%z == 1) then
             d%shown(i,:) = .false.
             d%shown(:,i) = .false.
          end if
       end do
       call d%generate_neighstars(r)
    elseif (r%flavor == repflavor_atoms_hbonds) then
       ! hydrogen bonds
       d%use_sys_nstar = .false.
       d%shown = .false.
       do i = 1, sys(r%id)%c%nspc
          if (sys(r%id)%c%spc(i)%z /= 1) cycle
          do j = 1, sys(r%id)%c%nspc
             iz = sys(r%id)%c%spc(j)%z
             if (iz == 7 .or. iz == 8 .or. iz == 9 .or. iz == 16 .or. iz == 17) then
                d%shown(i,j) = .true.
                d%shown(j,i) = .true.
             end if
          end do
       end do
       call d%generate_neighstars(r)
    else
       ! other flavors are default (track system bonds)
       d%use_sys_nstar = .true.
       call d%copy_neighstars_from_system(r%id)
    end if

  end subroutine bond_style_reset

  !> Deallocate all arrays and end the bond syle.
  module subroutine bond_style_end(d)
    class(bond_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    d%use_sys_nstar = .true.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%nstar)) deallocate(d%nstar)

  end subroutine bond_style_end

  !> Reset label style with default values. Use the information in
  !> representation r.
  module subroutine label_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use tools_io, only: nameguess, string
    class(label_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i

    ! not initialized
    d%isinit = .false.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%str)) deallocate(d%str)

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! fill according to the style
    d%isinit = .true.
    select case(r%labels%type)
    case (0,5,6)
       d%ntype = sys(r%id)%c%nspc
    case (2,3)
       d%ntype = sys(r%id)%c%ncel
    case (1,4,8)
       d%ntype = sys(r%id)%c%nneq
    case (7)
       d%ntype = sys(r%id)%c%nmol
    end select
    allocate(d%shown(d%ntype))
    allocate(d%str(d%ntype))

    ! fill shown, exclude hydrogens
    d%shown = .true.
    do i = 1, d%ntype
       select case(r%labels%type)
       case (0,5,6)
          if (sys(r%id)%c%spc(i)%z == 1) d%shown(i) = .false.
       case (2,3)
          if (sys(r%id)%c%spc(sys(r%id)%c%atcel(i)%is)%z == 1) d%shown(i) = .false.
       case (1,4,8)
          if (sys(r%id)%c%spc(sys(r%id)%c%at(i)%is)%z == 1) d%shown(i) = .false.
       end select
    end do

    ! fill text
    do i = 1, d%ntype
       if (r%labels%type == 0) then ! 0 = atomic symbol
          d%str(i) = trim(nameguess(sys(r%id)%c%spc(i)%z,.true.))
       elseif (r%labels%type == 1) then ! 1 = atom name
          d%str(i) = trim(sys(r%id)%c%at(i)%name)
       elseif (r%labels%type == 6) then ! 6 = Z
          d%str(i) = string(sys(r%id)%c%spc(i)%z)
       elseif (r%labels%type == 8) then ! 8 = wyckoff
          d%str(i) = string(sys(r%id)%c%at(i)%mult) // string(sys(r%id)%c%at(i)%wyc)
       else
          d%str(i) = string(i)
       end if
    end do

  end subroutine label_style_reset

  !> Deallocate all arrays and end the label syle.
  module subroutine label_style_end(d)
    class(label_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%str)) deallocate(d%str)

  end subroutine label_style_end

  !> Reset the coordination-polyhedra style to defaults from the
  !> system pointed at by representation r.
  module subroutine coordpoly_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sysc, sys_ready, ok_system, atlisttype_species
    use param, only: atmcov0, atmeneg, maxzat
    use global, only: bondfactor_def
    class(coordpoly_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    ! typical-anion species used as default corners: N, O, F, S, Cl, Br, I
    integer, parameter :: zanion(7) = (/7,8,9,16,17,35,53/)

    integer :: i, j, ispc, iz, jz, nspc, navg, nvalid
    real*8 :: dd, avgeneg
    logical :: useavg
    logical, allocatable :: spccenter(:), spccorner(:)

    ! if not initialized, set type
    if (.not.d%isinit) d%type = atlisttype_species

    ! reset the style to zero
    d%ntype = 0
    d%isinit = .false.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%corner)) deallocate(d%corner)
    if (allocated(d%dmin)) deallocate(d%dmin)
    if (allocated(d%dmax)) deallocate(d%dmax)

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! count the number of species that can act as centers or corners
    nspc = sys(r%id)%c%nspc
    allocate(spccenter(nspc),spccorner(nspc))
    spccenter = .false.
    spccorner = .false.
    nvalid = 0
    do j = 1, nspc
       jz = sys(r%id)%c%spc(j)%z
       if (jz > 0 .and. jz <= maxzat) nvalid = nvalid + 1
    end do

    if (nvalid == 1) then
       ! single species: it acts as both center and corner
       do j = 1, nspc
          jz = sys(r%id)%c%spc(j)%z
          if (jz > 0 .and. jz <= maxzat) then
             spccenter(j) = .true.
             spccorner(j) = .true.
          end if
       end do
    else
       ! If typical species that act as anions are present, use them
       ! as corners. Otherwise, use the electronegativity average to
       ! determine centers and corners.
       useavg = .true.
       do j = 1, nspc
          jz = sys(r%id)%c%spc(j)%z
          if (jz > 0 .and. jz <= maxzat .and. any(zanion == jz)) then
             useavg = .false.
             exit
          end if
       end do

       if (.not.useavg) then
          ! corners are the anions, centers are everyone else
          do j = 1, nspc
             jz = sys(r%id)%c%spc(j)%z
             if (jz > 0 .and. jz <= maxzat) then
                spccorner(j) = any(zanion == jz)
                spccenter(j) = .not.spccorner(j)
             end if
          end do
          ! all present species are anions: fall back to average
          if (.not.any(spccenter)) then
             useavg = .true.
             spccenter = .false.
             spccorner = .false.
          end if
       end if

       if (useavg) then
          ! average electronegativity over species with a defined value
          avgeneg = 0d0
          navg = 0
          do j = 1, nspc
             jz = sys(r%id)%c%spc(j)%z
             if (jz > 0 .and. jz <= maxzat .and. atmeneg(jz) > 0d0) then
                avgeneg = avgeneg + atmeneg(jz)
                navg = navg + 1
             end if
          end do
          if (navg > 0) avgeneg = avgeneg / navg
          do j = 1, nspc
             jz = sys(r%id)%c%spc(j)%z
             if (jz > 0 .and. jz <= maxzat .and. atmeneg(jz) > 0d0) then
                spccenter(j) = atmeneg(jz) < avgeneg
                spccorner(j) = atmeneg(jz) > avgeneg
             end if
          end do
       end if
    end if

    ! fill the style, use the rcov sum times bondfactor as the distance cutoff
    d%ntype = sysc(r%id)%attype_number(d%type)
    allocate(d%shown(d%ntype),d%corner(nspc,d%ntype),d%dmin(d%ntype),d%dmax(d%ntype))
    d%dmin = 0d0
    d%dmax = 0d0
    do i = 1, d%ntype
       ispc = sysc(r%id)%attype_species(d%type,i)
       iz = sys(r%id)%c%spc(ispc)%z
       d%shown(i) = spccenter(ispc)
       do j = 1, nspc
          d%corner(j,i) = spccorner(j)
          if (d%corner(j,i)) then
             jz = sys(r%id)%c%spc(j)%z
             dd = (atmcov0(iz) + atmcov0(jz)) * bondfactor_def
             d%dmax(i) = max(d%dmax(i),dd)
          end if
       end do
    end do
    d%isinit = .true.

  end subroutine coordpoly_style_reset

  !> Deallocate all arrays and end the coordination-polyhedra style.
  module subroutine coordpoly_style_end(d)
    class(coordpoly_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    d%ntype = 0
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%corner)) deallocate(d%corner)
    if (allocated(d%dmin)) deallocate(d%dmin)
    if (allocated(d%dmax)) deallocate(d%dmax)

  end subroutine coordpoly_style_end

end submodule proc
