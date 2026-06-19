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
    r%atoms_display = .true.
    r%bonds_display = .true.
    r%labels_display = .false.
    r%poly_display = .false.

    ! if type was given, mark as initialized and set name
    if (itype == reptype_atoms) then
       r%isinit = .true.
       r%shown = .true.
       if (flavor == repflavor_atoms_ballandstick) then
          r%name = "Ball and Stick"
       elseif (flavor == repflavor_atoms_sticks) then
          r%name = "Bonds"
          r%atoms_display = .false.
       elseif (flavor == repflavor_atoms_licorice) then
          r%name = "Licorice"
       elseif (flavor == repflavor_atoms_vdwcontacts) then
          r%name = "VdW Contacts"
          r%atoms_display = .false.
       elseif (flavor == repflavor_atoms_hbonds) then
          r%name = "Hydrogen Bonds"
          r%atoms_display = .false.
       elseif (flavor == repflavor_atoms_criticalpoints) then
          r%name = "Critical Points"
          r%bonds_display = .false.
          r%labels_display = .true.
       elseif (flavor == repflavor_atoms_gradientpaths) then
          r%name = "Gradient Paths"
          r%bonds_display = .false.
          r%labels_display = .false.
       elseif (flavor == repflavor_atoms_polyhedra) then
          r%name = "Polyhedra"
          r%bonds_display = .false.
          r%poly_display = .true.
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
       r%name = "Symmetry element"
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

  !> Set all values to default for the representation. Set a subset of
  !> defaults if itype = 0 (all), 1 (atom), 2 (bonds), 3 (labels),
  !> 4 (mol), 5 (unit cell), 6 (cartesian axes), 7 (rotation axes),
  !> 8 (coordination polyhedra).
  module subroutine representation_set_defaults(r,itype)
    use systems, only: sys, sys_ready, ok_system
    use global, only: bondfactor_def, bonddelta_def
    use gui_main, only: ColorAtomBorder_def, ColorBond_def, ColorBondBorder_def,&
       ColorLabel_def, ColorRotaxis_def, ColorAxes_def
    use param, only: atmcov0
    class(representation), intent(inout) :: r
    integer, intent(in) :: itype

    integer :: isys

    ! check the system is sane
    isys = r%id
    if (.not.ok_system(isys,sys_ready)) return

    !! initialize an empty representation
    if (itype == 0) then
       ! global parameters
       r%pertype = 1
       if (sys(isys)%c%ismolecule) then
          r%ncell = 0
       else
          r%ncell = 1
       end if
       r%origin = 0d0
       r%tshift = 0d0
       ! atoms, bonds, labels
       r%filter = ""
       r%errfilter = ""
       if (sys(isys)%c%ismolecule) then
          r%border = .false.
          r%onemotif = .false.
       else
          r%border = .true.
          r%onemotif = (sys(isys)%c%nmol > 1)
       end if
    end if

    !--> atoms
    if (itype == 0 .or. itype == 1) then
       r%atom_radii_type = 0
       r%atom_radii_scale = atomcovradscale_def
       r%atom_radii_value = atomconstantrad_def
       r%atom_color_type = 0
       r%atom_border_size = atomborder_def
       r%atom_border_rgb = ColorAtomBorder_def
       if (r%flavor == repflavor_atoms_licorice) then
          r%atom_radii_type = 2
          r%atom_radii_value = atomrad_licorice_def
       elseif (r%flavor == repflavor_atoms_criticalpoints) then
          r%atom_radii_type = 2
          r%atom_radii_value = atomrad_criticalpoints_def
          r%atom_border_size = atomborder_criticalpoints_def
       elseif (r%flavor == repflavor_atoms_gradientpaths) then
          r%atom_radii_type = 2
          r%atom_radii_value = atomrad_gradientpaths_def
          r%atom_border_size = atomborder_gradientpaths_def
       end if
    end if

    !--> bonds
    if (itype == 0 .or. itype == 2) then
       r%bond_atmrad = atmcov0
       r%bond_bfactor = bondfactor_def
       r%bond_bdelta = bonddelta_def
       r%bond_color_style = 0
       r%bond_border_size = bondborder_def
       r%bond_rad = bondrad_def
       r%bond_border_rgb = ColorBondBorder_def
       r%bond_rgb = ColorBond_def
       r%bond_order = 4 ! calculated (value from ordcon)
       r%bond_imol = 0
       r%bond_bothends = .true.
       if (r%flavor == repflavor_atoms_sticks) then
          r%bond_color_style = 1
          r%bond_border_size = bondborder_stickflav_def
       elseif (r%flavor == repflavor_atoms_licorice) then
          r%bond_color_style = 1
          r%bond_rad = bondrad_licorice_def
       end if
    end if

    !--> labels
    if (itype == 0 .or. itype == 3) then
       r%label_type = 0
       r%label_scale = 0.5d0
       r%label_rgb = ColorLabel_def
       r%label_const_size = .false.
       r%label_offset = (/0d0,0d0,0d0/)
       if (r%flavor == repflavor_atoms_criticalpoints) then
          r%label_type = 2 ! cell atom
          r%label_scale = 0.3d0
          r%label_offset = (/0d0,0.25d0,0d0/)
       end if
    end if

    ! unit cell
    if (itype == 0 .or. itype == 5) then
       r%uc_inner = .true.
       r%uc_coloraxes = .true.
       r%uc_vaccutsticks = .true.
       r%uc_radius = uc_radius_def
       r%uc_radiusinner = uc_radiusinner_def
       r%uc_rgb = 0._c_float
       r%uc_innersteplen = uc_innersteplen_def
       r%uc_innerstipple = .true.
    end if

    ! cartesian axes
    if (itype == 0 .or. itype == 6) then
       r%axes_kind = 0 ! cartesian
       r%axes_rot = eye ! no extra orientation by default
       r%axes_placement = 1
       if (sys(isys)%c%ismolecule) then
          r%axes_coordtype = 1 ! cartesian (angstrom)
       else
          r%axes_coordtype = 0 ! crystallographic
       end if
       r%axes_winpos = axes_winpos_def
       r%axes_length = axes_length_def
       r%axes_radius = axes_radius_def
       r%axes_conelength = axes_conelength_def
       r%axes_coneradius = axes_coneradius_def
       r%axes_rgb = ColorAxes_def ! x = red, y = green, z = blue
       r%axes_showlabels = .true.
       r%axes_labelscale = axes_labelscale_def
       r%axes_labelconstsize = .false.
       r%axes_labeldistance = axes_labeldistance_def
       r%axes_labeloffset = 0d0
       r%axes_scalewithzoom = .false.
       r%axes_scale = 1d0
       r%axes_scale_auto = (r%axes_placement == 1)
       r%axes_labelrgb = 0._c_float
       r%axes_labelstr(1) = "x"
       r%axes_labelstr(2) = "y"
       r%axes_labelstr(3) = "z"
    end if

    ! rotation axis
    if (itype == 0 .or. itype == 7) then
       r%rotaxis_dir = (/0d0,0d0,1d0/)
       r%rotaxis_length = 0d0
       r%rotaxis_radius = rotaxis_radius_def
       r%rotaxis_rgb = ColorRotaxis_def ! black
    end if

    ! coordination polyhedra appearance (part of the atoms representation; the
    ! center/corner/distance geometry lives in coordpoly_style, reset below)
    if (itype == 0 .or. itype == 1) then
       r%poly_alpha = 0.5d0
       r%poly_usecentercolor = .true.
       r%poly_rgb = 0._c_float
       r%poly_edge_rad = 0.05d0
       r%poly_edge_rgb = 0._c_float
       r%poly_usecentercolor_edge = .true.
       r%poly_coplanar_eps = 0.1d0
       r%poly_showcorners = .true.
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
    r%filter = ""
    r%errfilter = ""

    call r%atom_style%end()
    call r%bond_style%end()
    call r%label_style%end()
    call r%mol_style%end()
    call r%coordpoly_style%end()

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

    ! check if we need to reset the representation styles
    ! atoms
    doreset = .not.r%atom_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%atom_style%timelastreset)
    if (doreset) call r%atom_style%reset(r)

    ! bonds: if the geometry changed
    doreset = .not.r%bond_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%bond_style%timelastreset)
    if (doreset) call r%bond_style%reset(r)

    ! bonds: if the system has been rebonded and this representation tracks
    ! the bonds in the system (%use_sys_nstar), recalculate the bond style
    doreset = r%bond_style%use_sys_nstar .and. (sysc(r%id)%timelastchange_rebond > r%bond_style%timelastreset)
    if (doreset) call r%bond_style%copy_neighstars_from_system(r%id)

    ! molecules: if the geometry or the bonds changed
    doreset = .not.r%mol_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_rebond > r%mol_style%timelastreset)
    if (doreset) call r%mol_style%reset(r)

    ! labels: if the geometry changed
    doreset = .not.r%label_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%label_style%timelastreset)
    if (doreset) call r%label_style%reset(r)

    ! coordination polyhedra: if the geometry changed
    doreset = .not.r%coordpoly_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%coordpoly_style%timelastreset)
    if (doreset) call r%coordpoly_style%reset(r)

  end subroutine update_styles

  !> Add the spheres, cylinder, etc. to the draw lists. Use nc number
  !> of cells and the data from representation r. If doanim, use qpt
  !> iqpt and frequency ifreq to animate the representation.
  module subroutine add_draw_elements(r,nc,obj,doanim,iqpt,ifreq)
    use systems, only: sys, sysc
    use crystalmod, only: iperiod_vacthr
    use gui_main, only: ColorAxes_def
    use tools_io, only: string, nameguess
    use tools_math, only: cross, plane_from_points
    use types, only: realloc
    use tools, only: mergesort
    use param, only: tpi, img, atmass, icrd_crys
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
    real(c_float) :: rgb(3)
    real*8 :: rad1, rad2, dd, f1, f2, axsc
    real*8 :: xx(3), xc(3), x0(3), x1(3), x2(3), res, uoriginc(3), xpolyc(3)
    real*8 :: ucini(3), ucend(3), e1v(3), e2v(3)
    real(c_float) :: rgbax(3)
    complex*16 :: xdelta0(3), xdelta1(3), xdelta2(3)
    type(dl_sphere), allocatable :: auxsph(:)
    type(dl_cylinder), allocatable :: auxcyl(:)
    type(dl_string), allocatable :: auxstr(:)
    type(dl_plane), allocatable :: auxplane(:)
    type(dl_cylinder) :: dcyl
    type(dl_string) :: dstr
    logical :: fixed
    character(len=:), allocatable :: errmsg
    ! coordination polyhedra
    real*8, allocatable :: xvpoly(:,:), up2dsp(:,:)
    complex*16, allocatable :: dvpoly(:,:)
    integer, allocatable :: eidp(:), lvecp(:,:)
    real(c_float) :: rgbface(3), rgbedge(3)
    integer :: natp, idpoly, kp
    logical :: dopoly, corneractive
    ! corner atoms forced visible (poly_showcorners)
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

    ! initialize the drawlists if not done already
    if (.not.allocated(obj%sph)) then
       allocate(obj%sph(100))
       obj%nsph = 0
    end if
    if (.not.allocated(obj%cyl)) then
       allocate(obj%cyl(100))
       obj%ncyl = 0
    end if
    if (.not.allocated(obj%cylflat)) then
       allocate(obj%cylflat(100))
       obj%ncylflat = 0
    end if
    if (.not.allocated(obj%cone)) then
       allocate(obj%cone(100))
       obj%ncone = 0
    end if
    if (.not.allocated(obj%cylgiz)) then
       allocate(obj%cylgiz(10))
       obj%ncylgiz = 0
    end if
    if (.not.allocated(obj%conegiz)) then
       allocate(obj%conegiz(10))
       obj%nconegiz = 0
    end if
    if (.not.allocated(obj%stringgiz)) then
       allocate(obj%stringgiz(10))
       obj%nstringgiz = 0
    end if
    if (.not.allocated(obj%string)) then
       allocate(obj%string(100))
       obj%nstring = 0
    end if
    doanim_ = doanim
    if (doanim_) doanim_ = doanim_ .and. (iqpt > 0 .and. ifreq > 0 .and. sys(r%id)%c%vib%hasvibs)

    if (r%type == reptype_atoms) then
       !!! atoms and bonds representation !!!

       !! first, the atoms
       ! do we have a filter?
       havefilter = (len_trim(r%filter) > 0) .and. (len_trim(r%errfilter) == 0)
       usetshift = any(abs(r%tshift) > 1d-5)

       ! calculate the periodicity
       n = 1
       if (r%pertype == 1) then
          n = nc
       elseif (r%pertype == 2) then
          n = r%ncell
       end if

       ! origin shift
       if (sys(r%id)%c%ismolecule) then
          uoriginc = r%origin / bohrtoa
       else
          uoriginc = sys(r%id)%c%x2c(r%origin)
       end if

       ! whether we will force the polyhedra corner atoms to be drawn (only when
       ! the representation displays atoms, so lshown tracks visible spheres)
       corneractive = r%poly_display .and. r%poly_showcorners .and. r%atoms_display .and.&
          r%coordpoly_style%isinit

       ! whether we'll be doing bonds, allocate array to check whether
       ! an atom has been drawn (also used to deduplicate forced corner atoms)
       dobonds = r%bonds_display .and. r%bond_style%isinit
       if (dobonds .or. corneractive) then
          allocate(lshown(sys(r%id)%c%ncel,-1:n(1)+1,-1:n(2)+1,-1:n(3)+1))
          lshown = .false.
       end if
       ncorn = 0
       if (corneractive) allocate(cornlist(4,100))

       ! coordination polyhedra: per-species distance-window scratch
       if (r%poly_display .and. r%coordpoly_style%isinit) &
          allocate(up2dsp(sys(r%id)%c%nspc,2))

       ! whether there is vacuum in any direction
       dovac = (sys(r%id)%c%vaclength > iperiod_vacthr)
       if (any(dovac)) then
          ucini = sys(r%id)%c%vactop - 1d0 - vacextension / sys(r%id)%c%aa
          ucend = sys(r%id)%c%vacbot + vacextension / sys(r%id)%c%aa
       end if

       ! run over atoms, either directly or per-molecule
       i = 0
       imol = 0
       do while(.true.)
          if (r%onemotif) then
             ! this is a new molecule if there are no molecules or this is the last atom
             ! in the previous one
             step = (imol == 0)
             if (.not.step) step = (k == sys(r%id)%c%mol(imol)%nat)
             if (step) then
                imol = imol + 1
                k = 0
             end if

             ! we are finished if we have all molecules
             if (imol > sys(r%id)%c%nmol) exit

             ! Add the new atom; only translate if the fragment is discrete
             k = k + 1
             i = sys(r%id)%c%mol(imol)%at(k)%cidx
             if (sys(r%id)%c%mol(imol)%discrete) then
                lvec = sys(r%id)%c%mol(imol)%at(k)%lvec
             else
                lvec = 0
             end if
          else
             ! next atom in the complete list, exit if done
             i = i + 1
             if (i > sys(r%id)%c%ncel) exit
             lvec = 0
             imol = sys(r%id)%c%idatcelmol(1,i)
          end if
          ! i is current atom from the complete atom list
          ! imol is the corresponding molecule

          ! skip hidden atoms
          id = sysc(r%id)%attype_celatom_to_id(r%atom_style%type,i)
          if (.not.r%atom_style%shown(id)) cycle

          ! skip hidden molecules
          if (.not.r%mol_style%shown(imol)) cycle

          ! calculate the border
          xx = sys(r%id)%c%atcel(i)%x
          n0 = 0
          n1 = n-1
          if (r%border.and..not.r%onemotif) then
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
             xx = sys(r%id)%c%atcel(i)%x
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
          rgb = r%atom_style%rgb(:,id) * r%mol_style%tint_rgb(:,imol)
          rad1 = r%atom_style%rad(id) * r%mol_style%scale_rad(imol)

          ! coordination polyhedron: if this atom is a shown center, find its
          ! corner atoms once (translation-invariant; reused for every image)
          dopoly = .false.
          if (r%poly_display .and. r%coordpoly_style%isinit) then
             idpoly = sysc(r%id)%attype_celatom_to_id(r%coordpoly_style%type,i)
             if (r%coordpoly_style%shown(idpoly) .and. r%coordpoly_style%dmax(idpoly) > 0d0 .and.&
                any(r%coordpoly_style%corner(:,idpoly))) then
                do j = 1, sys(r%id)%c%nspc
                   if (r%coordpoly_style%corner(j,idpoly)) then
                      up2dsp(j,1) = r%coordpoly_style%dmin(idpoly)
                      up2dsp(j,2) = r%coordpoly_style%dmax(idpoly)
                   else
                      up2dsp(j,:) = 0d0
                   end if
                end do
                call sys(r%id)%c%list_near_atoms(sys(r%id)%c%atcel(i)%x,icrd_crys,.true.,&
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
                   if (r%poly_usecentercolor) then
                      rgbface = rgb
                   else
                      rgbface = r%poly_rgb
                   end if
                   if (r%poly_usecentercolor_edge) then
                      rgbedge = rgb
                   else
                      rgbedge = r%poly_edge_rgb
                   end if
                end if
             end if
          end if

          do i1 = n0(1), n1(1)
             do i2 = n0(2), n1(2)
                do i3 = n0(3), n1(3)
                   ix = (/i1,i2,i3/) + lvec + vacshift
                   if (usetshift) then
                      xx = sys(r%id)%c%atcel(i)%x - r%tshift
                      ix = ix + nint(xx - floor(xx) + r%tshift - sys(r%id)%c%atcel(i)%x)
                   end if

                   xx = sys(r%id)%c%atcel(i)%x + ix
                   xc = sys(r%id)%c%x2c(xx)

                   ! apply the filter
                   if (havefilter) then
                      res = sys(r%id)%eval(r%filter,errmsg,xc)
                      if (len_trim(errmsg) == 0) then
                         if (res == 0d0) cycle
                      else
                         havefilter = .false.
                         r%errfilter = errmsg
                      end if
                   end if

                   ! draw the coordination polyhedron for this center image
                   ! (corners are the search result translated to this image)
                   if (dopoly) then
                      do kp = 1, natp
                         xpolyc = sys(r%id)%c%atcel(eidp(kp))%x + lvecp(:,kp) + ix
                         xvpoly(:,kp) = sys(r%id)%c%x2c(xpolyc) + uoriginc
                         dvpoly(:,kp) = vibdelta(eidp(kp),xpolyc) ! per-corner vibration delta
                      end do
                      call build_polyhedron(xvpoly(:,1:natp),dvpoly(:,1:natp),natp,rgbface,rgbedge,&
                         r%poly_alpha,r%poly_edge_rad,r%poly_coplanar_eps)

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
                   xdelta1 = vibdelta(i,xx)

                   ! draw the atom, reallocate if necessary
                   if (r%atoms_display) then
                      obj%nsph = obj%nsph + 1
                      if (obj%nsph > size(obj%sph,1)) then
                         allocate(auxsph(2*obj%nsph))
                         auxsph(1:size(obj%sph,1)) = obj%sph
                         call move_alloc(auxsph,obj%sph)
                      end if

                      ! write down the sphere
                      obj%sph(obj%nsph)%x = real(xc + uoriginc,c_float)
                      obj%sph(obj%nsph)%r = real(rad1,c_float)
                      obj%sph(obj%nsph)%rgb = rgb
                      obj%sph(obj%nsph)%idx(1) = i
                      obj%sph(obj%nsph)%idx(2:4) = ix
                      obj%sph(obj%nsph)%xdelta = cmplx(xdelta1,kind=c_float_complex)
                      obj%sph(obj%nsph)%border = real(r%atom_border_size,c_float)
                      obj%sph(obj%nsph)%rgbborder = r%atom_border_rgb
                   else if (dobonds) then
                      ! atoms hidden but bonds shown: add an invisible pick target at
                      ! the atom site so hidden atoms can still be picked, measured,
                      ! and box-selected. Carries the real idx; rendered only into the
                      ! pick buffer (skipped in the visible sphere pass).
                      obj%nsph = obj%nsph + 1
                      if (obj%nsph > size(obj%sph,1)) then
                         allocate(auxsph(2*obj%nsph))
                         auxsph(1:size(obj%sph,1)) = obj%sph
                         call move_alloc(auxsph,obj%sph)
                      end if

                      ! write down the ghost sphere
                      obj%sph(obj%nsph)%x = real(xc + uoriginc,c_float)
                      obj%sph(obj%nsph)%r = real(2d0 * r%bond_rad,c_float) ! generous click radius
                      obj%sph(obj%nsph)%rgb = rgb
                      obj%sph(obj%nsph)%idx(1) = i
                      obj%sph(obj%nsph)%idx(2:4) = ix
                      obj%sph(obj%nsph)%xdelta = cmplx(xdelta1,kind=c_float_complex)
                      obj%sph(obj%nsph)%border = 0._c_float
                      obj%sph(obj%nsph)%rgbborder = rgb
                      obj%sph(obj%nsph)%ghost = .true.
                   end if

                   ! mark this atom image as drawn (for bonds and for
                   ! deduplicating forced polyhedra corner atoms)
                   if (allocated(lshown)) then
                      call check_lshown(i,ix(1),ix(2),ix(3))
                      lshown(i,ix(1),ix(2),ix(3)) = .true.
                   end if

                   ! bonds
                   if (dobonds) then
                      ! bonds
                      do ib = 1, r%bond_style%nstar(i)%ncon
                         ineigh = r%bond_style%nstar(i)%idcon(ib)
                         if (.not.r%bond_style%shown(sys(r%id)%c%atcel(ineigh)%is,sys(r%id)%c%atcel(i)%is)) cycle
                         ixn = ix + r%bond_style%nstar(i)%lcon(:,ib)

                         if (r%bond_imol == 1) then ! intramol
                            if (.not.sys(r%id)%c%in_same_molecule(i,ix,ineigh,ixn)) cycle
                         elseif (r%bond_imol == 2) then ! intermol
                            if (sys(r%id)%c%in_same_molecule(i,ix,ineigh,ixn)) cycle
                         end if

                         call check_lshown(ineigh,ixn(1),ixn(2),ixn(3))
                         if (r%bond_bothends) then
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
                         if (r%bond_order == 4) then
                            iord = r%bond_style%nstar(i)%ordcon(ib)
                         else
                            iord = r%bond_order
                         end if

                         ! draw the bond, reallocate if necessary
                         if (r%bond_color_style == 0) then
                            obj%ncyl = obj%ncyl + 1
                         else
                            obj%ncyl = obj%ncyl + 2
                         end if
                         if (obj%ncyl > size(obj%cyl,1)) then
                            allocate(auxcyl(2*obj%ncyl))
                            auxcyl(1:size(obj%cyl,1)) = obj%cyl
                            call move_alloc(auxcyl,obj%cyl)
                         end if

                         ! animation delta of the other end
                         xdelta2 = vibdelta(ineigh,sys(r%id)%c%atcel(ineigh)%x + ixn)

                         x1 = xc + uoriginc
                         x2 = sys(r%id)%c%atcel(ineigh)%x + ixn
                         x2 = sys(r%id)%c%x2c(x2) + uoriginc
                         if (r%bond_color_style == 0) then
                            obj%cyl(obj%ncyl)%x1 = real(x1,c_float)
                            obj%cyl(obj%ncyl)%x1delta = cmplx(xdelta1,kind=c_float_complex)
                            obj%cyl(obj%ncyl)%x2 = real(x2,c_float)
                            obj%cyl(obj%ncyl)%x2delta = cmplx(xdelta2,kind=c_float_complex)
                            obj%cyl(obj%ncyl)%r = real(r%bond_rad,c_float)
                            obj%cyl(obj%ncyl)%rgb = r%bond_rgb
                            obj%cyl(obj%ncyl)%order = iord
                            obj%cyl(obj%ncyl)%border = real(r%bond_border_size,c_float)
                            obj%cyl(obj%ncyl)%rgbborder = r%bond_border_rgb
                            obj%cyl(obj%ncyl)%arvec = real(r%bond_style%nstar(i)%aromdir(:,ib),c_float)
                         else
                            idaux = sysc(r%id)%attype_celatom_to_id(r%atom_style%type,ineigh)
                            rad2 = r%atom_style%rad(idaux) * r%mol_style%scale_rad(sys(r%id)%c%idatcelmol(1,ineigh))
                            dd = norm2(x2 - x1)
                            f1 = min(max((0.5d0 + 0.5d0 * (rad2 - rad1) / dd),0._c_float),1._c_float)
                            f2 = 1._c_float - f1
                            x0 = f1 * x1 + f2 * x2
                            xdelta0 = f1 * xdelta1 + f2 * xdelta2

                            ! add the two cylinders to the list
                            obj%cyl(obj%ncyl-1)%x1 = real(x1,c_float)
                            obj%cyl(obj%ncyl-1)%x1delta = cmplx(xdelta1,kind=c_float_complex)
                            obj%cyl(obj%ncyl-1)%x2 = real(x0,c_float)
                            obj%cyl(obj%ncyl-1)%x2delta = cmplx(xdelta0,kind=c_float_complex)
                            obj%cyl(obj%ncyl-1)%r = real(r%bond_rad,c_float)
                            obj%cyl(obj%ncyl-1)%rgb = rgb
                            obj%cyl(obj%ncyl-1)%order = iord
                            obj%cyl(obj%ncyl-1)%border = real(r%bond_border_size,c_float)
                            obj%cyl(obj%ncyl-1)%rgbborder = r%bond_border_rgb
                            obj%cyl(obj%ncyl-1)%arvec = real(r%bond_style%nstar(i)%aromdir(:,ib),c_float)

                            obj%cyl(obj%ncyl)%x1 = real(x0,c_float)
                            obj%cyl(obj%ncyl)%x1delta = cmplx(xdelta0,kind=c_float_complex)
                            obj%cyl(obj%ncyl)%x2 = real(x2,c_float)
                            obj%cyl(obj%ncyl)%x2delta = cmplx(xdelta2,kind=c_float_complex)
                            obj%cyl(obj%ncyl)%r = real(r%bond_rad,c_float)
                            obj%cyl(obj%ncyl)%rgb = r%atom_style%rgb(:,idaux) * &
                               r%mol_style%tint_rgb(:,sys(r%id)%c%idatcelmol(1,ineigh))
                            obj%cyl(obj%ncyl)%order = iord
                            obj%cyl(obj%ncyl)%border = real(r%bond_border_size,c_float)
                            obj%cyl(obj%ncyl)%rgbborder = r%bond_border_rgb
                            obj%cyl(obj%ncyl)%arvec = real(r%bond_style%nstar(i)%aromdir(:,ib),c_float)
                         end if
                      end do ! ncon
                   end if

                   if (r%labels_display) then
                      select case(r%label_type)
                      case (0,5,6)
                         idl = sys(r%id)%c%atcel(i)%is
                      case (2,3)
                         idl = i
                      case (1,4,8)
                         idl = sys(r%id)%c%atcel(i)%idx
                      case (7)
                         idl = sys(r%id)%c%idatcelmol(1,i)
                      end select

                      ! labels
                      if (r%label_style%shown(idl)) then
                         obj%nstring = obj%nstring + 1
                         if (obj%nstring > size(obj%string,1)) then
                            allocate(auxstr(2*obj%nstring))
                            auxstr(1:size(obj%string,1)) = obj%string
                            call move_alloc(auxstr,obj%string)
                         end if

                         obj%string(obj%nstring)%x = real(xc + uoriginc,c_float)
                         obj%string(obj%nstring)%xdelta = cmplx(xdelta1,kind=c_float_complex)
                         obj%string(obj%nstring)%r = real(rad1,c_float)
                         obj%string(obj%nstring)%rgb = r%label_rgb
                         if (r%label_const_size) then
                            obj%string(obj%nstring)%scale = real(r%label_scale,c_float)
                         else
                            obj%string(obj%nstring)%scale = real(-r%label_scale,c_float)
                         end if
                         obj%string(obj%nstring)%offset = real(r%label_offset,c_float)
                         obj%string(obj%nstring)%str = trim(r%label_style%str(idl))
                         if (r%label_type == 3) then
                            ! add the lattice vectors
                            obj%string(obj%nstring)%str = obj%string(obj%nstring)%str // "[" //&
                               string(ix(1)) // "," // string(ix(2)) // "," //string(ix(3)) // "]"
                         end if
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
             idc = sysc(r%id)%attype_celatom_to_id(r%atom_style%type,cornlist(1,ica))
             imolc = sys(r%id)%c%idatcelmol(1,cornlist(1,ica))
             rgb = r%atom_style%rgb(:,idc) * r%mol_style%tint_rgb(:,imolc)
             rad1 = r%atom_style%rad(idc) * r%mol_style%scale_rad(imolc)
             xx = sys(r%id)%c%atcel(cornlist(1,ica))%x + ix
             xc = sys(r%id)%c%x2c(xx)

             ! animation delta of the corner atom (so it moves with the polyhedron)
             xdelta1 = vibdelta(cornlist(1,ica),xx)

             obj%nsph = obj%nsph + 1
             if (obj%nsph > size(obj%sph,1)) then
                allocate(auxsph(2*obj%nsph))
                auxsph(1:size(obj%sph,1)) = obj%sph
                call move_alloc(auxsph,obj%sph)
             end if
             obj%sph(obj%nsph)%x = real(xc + uoriginc,c_float)
             obj%sph(obj%nsph)%r = real(rad1,c_float)
             obj%sph(obj%nsph)%rgb = rgb
             obj%sph(obj%nsph)%idx(1) = cornlist(1,ica)
             obj%sph(obj%nsph)%idx(2:4) = ix
             obj%sph(obj%nsph)%xdelta = cmplx(xdelta1,kind=c_float_complex)
             obj%sph(obj%nsph)%border = real(r%atom_border_size,c_float)
             obj%sph(obj%nsph)%rgbborder = r%atom_border_rgb
          end do
       end if
       if (allocated(cornlist)) deallocate(cornlist)
       if (allocated(up2dsp)) deallocate(up2dsp)
    elseif (r%type == reptype_unitcell) then
       !!! unit cell representation !!!

       ! number of cells
       n = 1
       if (r%pertype == 1) then
          n = nc
       elseif (r%pertype == 2) then
          n = r%ncell
       end if

       ! vacuum directions: we have a vacuum and only one cell in that direction
       isvac = .false.
       if (r%uc_vaccutsticks) then
          do i = 1, 3
             if (sys(r%id)%c%vaclength(i) > iperiod_vacthr .and. n(i) == 1) then
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
          call increase_ncylflat()
          obj%cylflat(obj%ncylflat)%x1 = real(x1,c_float)
          obj%cylflat(obj%ncylflat)%x2 = real(x2,c_float)
          obj%cylflat(obj%ncylflat)%r = real(r%uc_radius,c_float)
          if (r%uc_coloraxes.and.i>=1.and.i<=3) then
             obj%cylflat(obj%ncylflat)%rgb = ColorAxes_def(:,i)
          else
             obj%cylflat(obj%ncylflat)%rgb = r%uc_rgb
          end if
       end do

       ! draw inner cylinders
       if (r%uc_inner) then
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
                      if (r%uc_innerstipple) then
                         nstep = ceiling(norm2(x2 - x1) / r%uc_innersteplen)
                         do j = 1, nstep
                            call increase_ncylflat()
                            obj%cylflat(obj%ncylflat)%x1 = &
                               real(x1 + real(2*j-1,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            obj%cylflat(obj%ncylflat)%x2 = &
                               real(x1 + real(2*j,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            obj%cylflat(obj%ncylflat)%r = real(r%uc_radiusinner,c_float)
                            obj%cylflat(obj%ncylflat)%rgb = r%uc_rgb
                         end do
                      else
                         call increase_ncylflat()
                         obj%cylflat(obj%ncylflat)%x1 = real(x1 ,c_float)
                         obj%cylflat(obj%ncylflat)%x2 = real(x2 ,c_float)
                         obj%cylflat(obj%ncylflat)%r = real(r%uc_radiusinner,c_float)
                         obj%cylflat(obj%ncylflat)%rgb = r%uc_rgb
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
       ! around the local origin and positioned at render time; record
       ! the requested window position in the gizmo draw list.
       fixed = (r%axes_placement == 1)
       if (fixed) then
          uoriginc = 0d0
          obj%gizwinpos = r%axes_winpos
          obj%gizscalewithzoom = r%axes_scalewithzoom
       else
          ! origin in the requested coordinate system, converted to
          ! cartesian (bohr)
          if (r%axes_coordtype == 2) then
             uoriginc = r%origin ! cartesian (bohr)
          elseif (r%axes_coordtype == 0 .and. .not.sys(r%id)%c%ismolecule) then
             uoriginc = sys(r%id)%c%x2c(r%origin) ! crystallographic
          else
             uoriginc = r%origin / bohrtoa ! cartesian (angstrom)
          end if
          ! for molecules, cartesian coordinates are referred to the molecular center
          if (sys(r%id)%c%ismolecule) uoriginc = uoriginc - sys(r%id)%c%molx0
       end if

       ! global scale factor applied to the whole gizmo (arrows and labels)
       axsc = r%axes_scale

       ! arrowhead geometry (head length capped so it never exceeds the
       ! total axis length)
       rad1 = min(r%axes_conelength, r%axes_length) * axsc ! head length
       rad2 = r%axes_coneradius * axsc ! head radius

       do k = 1, 3
          ! unit direction for this axis: cartesian (lab-frame) or along
          ! the crystallographic lattice vector
          if (r%axes_kind == 1 .and. .not.sys(r%id)%c%ismolecule) then
             x0 = sys(r%id)%c%m_x2c(:,k)
             x0 = x0 / norm2(x0)
          else
             x0 = 0d0
             x0(k) = 1d0
          end if
          ! reorient the axis directions (identity unless a frame was requested)
          x0 = matmul(r%axes_rot,x0)

          ! shaft (round, lit cylinder)
          x1 = uoriginc
          x2 = uoriginc + max(r%axes_length * axsc - rad1,0d0) * x0
          dcyl%x1 = real(x1,c_float)
          dcyl%x2 = real(x2,c_float)
          dcyl%x1delta = cmplx(0d0,0d0,kind=c_float_complex)
          dcyl%x2delta = cmplx(0d0,0d0,kind=c_float_complex)
          dcyl%r = real(r%axes_radius * axsc,c_float)
          dcyl%rgb = r%axes_rgb(:,k)
          dcyl%order = 1
          dcyl%border = 0._c_float
          dcyl%rgbborder = 0._c_float
          if (fixed) then
             call append_cyl(obj%cylgiz,obj%ncylgiz,dcyl)
          else
             call append_cyl(obj%cyl,obj%ncyl,dcyl)
          end if

          ! arrowhead (cone) from the shaft end to the axis tip
          dcyl%x1 = real(x2,c_float)
          dcyl%x2 = real(uoriginc + r%axes_length * axsc * x0,c_float)
          dcyl%r = real(rad2,c_float)
          if (fixed) then
             call append_cyl(obj%conegiz,obj%nconegiz,dcyl)
          else
             call append_cyl(obj%cone,obj%ncone,dcyl)
          end if

          ! label: along the axis from the arrowhead tip by the shared
          ! distance, plus the per-axis cartesian offset
          if (r%axes_showlabels) then
             dstr%x = real(uoriginc + (r%axes_length + r%axes_labeldistance) * axsc * x0 + &
                r%axes_labeloffset(:,k) * axsc,c_float)
             dstr%xdelta = cmplx(0d0,0d0,kind=c_float_complex)
             dstr%r = real(rad2,c_float)
             dstr%rgb = r%axes_labelrgb
             if (r%axes_labelconstsize) then
                dstr%scale = real(r%axes_labelscale * axsc,c_float)
             else
                dstr%scale = real(-r%axes_labelscale * axsc,c_float)
             end if
             dstr%offset = 0._c_float
             dstr%str = trim(r%axes_labelstr(k))
             if (fixed) then
                call append_str(obj%stringgiz,obj%nstringgiz,dstr)
             else
                call append_str(obj%string,obj%nstring,dstr)
             end if
          end if
       end do
    elseif (r%type == reptype_rotaxis) then
       !!! rotation-axis representation (single black cylinder through the origin) !!!

       ! origin in cartesian (bohr); for molecules referred to the molecular center
       uoriginc = r%origin
       if (sys(r%id)%c%ismolecule) uoriginc = uoriginc - sys(r%id)%c%molx0

       ! double-ended cylinder (the axis line): origin +/- length*dir
       x0 = r%rotaxis_dir / max(norm2(r%rotaxis_dir),1d-10)
       dcyl%x1 = real(uoriginc - r%rotaxis_length * x0,c_float)
       dcyl%x2 = real(uoriginc + r%rotaxis_length * x0,c_float)
       dcyl%x1delta = cmplx(0d0,0d0,kind=c_float_complex)
       dcyl%x2delta = cmplx(0d0,0d0,kind=c_float_complex)
       dcyl%r = real(r%rotaxis_radius,c_float)
       dcyl%rgb = r%rotaxis_rgb
       dcyl%order = 1
       dcyl%border = 0._c_float
       dcyl%rgbborder = 0._c_float
       call append_cyl(obj%cyl,obj%ncyl,dcyl)
    elseif (r%type == reptype_symelem) then
       !!! symmetry element (plane/axis) !!!

       ! unit direction: the plane normal or the axis direction (cartesian)
       x0 = r%symelem_dir / max(norm2(r%symelem_dir),1d-10)

       ! reference point in cartesian (bohr); for molecules referred to the
       ! molecular center
       uoriginc = r%origin
       if (sys(r%id)%c%ismolecule) uoriginc = uoriginc - sys(r%id)%c%molx0

       ! opaque thin-cylinder template (plane frame edges, axis shafts)
       dcyl%x1delta = cmplx(0d0,0d0,kind=c_float_complex)
       dcyl%x2delta = cmplx(0d0,0d0,kind=c_float_complex)
       dcyl%r = real(symelem_frame_radius,c_float)
       dcyl%rgb = r%symelem_rgb
       dcyl%alpha = 1._c_float
       dcyl%order = 1
       dcyl%border = 0._c_float
       dcyl%rgbborder = 0._c_float

       if (r%symelem_kind == symelem_kind_plane) then
          ! mirror/glide plane: translucent fill + opaque border frame

          ! in-plane orthonormal basis perpendicular to the plane normal x0
          if (abs(x0(1)) < 0.9d0) then
             xx = (/1d0,0d0,0d0/)
          else
             xx = (/0d0,1d0,0d0/)
          end if
          x1 = cross(x0,xx)
          x1 = x1 / norm2(x1)
          x2 = cross(x0,x1) ! unit, (x0,x1,x2) orthonormal

          ! rectangle center = projection of the system center onto the plane
          xc = r%symelem_cen - dot_product(r%symelem_cen - uoriginc,x0) * x0
          res = symelem_margin * r%symelem_size
          e1v = res * x1
          e2v = res * x2

          ! translucent fill
          obj%nplane = obj%nplane + 1
          if (obj%nplane > size(obj%plane,1)) then
             allocate(auxplane(2*obj%nplane))
             auxplane(1:size(obj%plane,1)) = obj%plane
             call move_alloc(auxplane,obj%plane)
          end if
          obj%plane(obj%nplane)%x = real(xc,c_float)
          obj%plane(obj%nplane)%e1 = real(e1v,c_float)
          obj%plane(obj%nplane)%e2 = real(e2v,c_float)
          obj%plane(obj%nplane)%rgb = r%symelem_rgb
          obj%plane(obj%nplane)%alpha = symelem_alpha

          ! opaque border frame (4 edge cylinders)
          call append_edge(xc - e1v - e2v, xc + e1v - e2v)
          call append_edge(xc + e1v - e2v, xc + e1v + e2v)
          call append_edge(xc + e1v + e2v, xc - e1v + e2v)
          call append_edge(xc - e1v + e2v, xc - e1v - e2v)
       else
          ! rotation/rotoinversion axis: a thick opaque shaft (colored by the
          ! rotation order) through every visible lattice point (crystals) or
          ! the molecular center (molecules)

          ! axis color from the rotation order
          rgbax = symelem_rgb_def
          if (r%symelem_order >= lbound(symelem_rgb_order,2) .and. &
              r%symelem_order <= ubound(symelem_rgb_order,2)) then
             if (any(symelem_rgb_order(:,r%symelem_order) /= 0._c_float)) &
                rgbax = symelem_rgb_order(:,r%symelem_order)
          end if
          dcyl%rgb = rgbax
          dcyl%r = real(symelem_axis_radius,c_float)

          ! each shaft spans the system bounding sphere plus a margin on each
          ! side, centered on the projection of the system center onto the axis
          ! line
          res = symelem_margin * r%symelem_size
          n0 = 0
          n1 = 0
          if (.not.sys(r%id)%c%ismolecule) n1 = nc
          do i1 = n0(1), n1(1)
             do i2 = n0(2), n1(2)
                do i3 = n0(3), n1(3)
                   xc = uoriginc + sys(r%id)%c%x2c(real((/i1,i2,i3/),8))
                   xc = xc + dot_product(r%symelem_cen - xc,x0) * x0 ! foot of the center on the axis
                   dcyl%x1 = real(xc - res * x0,c_float)
                   dcyl%x2 = real(xc + res * x0,c_float)
                   call append_cyl(obj%cyl,obj%ncyl,dcyl)
                end do
             end do
          end do
       end if
    end if ! reptype
  contains
    subroutine increase_ncone()

      obj%ncone = obj%ncone + 1
      if (obj%ncone > size(obj%cone,1)) then
         allocate(auxcyl(2*obj%ncone))
         auxcyl(1:size(obj%cone,1)) = obj%cone
         call move_alloc(auxcyl,obj%cone)
      end if

    end subroutine increase_ncone

    !> Append a cylinder record to an allocatable cylinder list,
    !> reallocating if necessary.
    subroutine append_cyl(lst,n,it)
      type(dl_cylinder), allocatable, intent(inout) :: lst(:)
      integer, intent(inout) :: n
      type(dl_cylinder), intent(in) :: it
      type(dl_cylinder), allocatable :: aux(:)

      n = n + 1
      if (n > size(lst,1)) then
         allocate(aux(2*n))
         aux(1:size(lst,1)) = lst
         call move_alloc(aux,lst)
      end if
      lst(n) = it

    end subroutine append_cyl

    !> Append an opaque thin edge cylinder between cartesian points p and q
    !> (uses the host dcyl template for radius/color).
    subroutine append_edge(p,q)
      real*8, intent(in) :: p(3), q(3)

      dcyl%x1 = real(p,c_float)
      dcyl%x2 = real(q,c_float)
      call append_cyl(obj%cyl,obj%ncyl,dcyl)

    end subroutine append_edge

    !> Build a coordination polyhedron from nv vertex positions xv (cartesian,
    !> bohr). The vertex centroid is used as the interior reference point (it is
    !> always inside the convex hull, unlike the cation for one-sided
    !> coordinations). Adds translucent triangular faces (color rgbf, opacity
    !> alphaf) and opaque edge cylinders (color rgbe, radius rade) to the draw
    !> lists. If the vertices are coplanar to within eps, a filled polygon is
    !> drawn instead of a 3D convex hull.
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

         do k = 1, nvv
            kk = mod(k,nvv) + 1
            call append_triangle(cen0,xv(:,iord_(k)),xv(:,iord_(kk)),&
               dcen,dv(:,iord_(k)),dv(:,iord_(kk)),rgbf,alphaf)
            call increase_ncylflat()
            obj%cylflat(obj%ncylflat)%x1 = real(xv(:,iord_(k)),c_float)
            obj%cylflat(obj%ncylflat)%x2 = real(xv(:,iord_(kk)),c_float)
            obj%cylflat(obj%ncylflat)%r = real(rade,c_float)
            obj%cylflat(obj%ncylflat)%rgb = rgbe
            obj%cylflat(obj%ncylflat)%x1delta = cmplx(dv(:,iord_(k)),kind=c_float_complex)
            obj%cylflat(obj%ncylflat)%x2delta = cmplx(dv(:,iord_(kk)),kind=c_float_complex)
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
      do e = 1, nedge
         if (.not.edgekeep(e)) cycle
         call increase_ncylflat()
         obj%cylflat(obj%ncylflat)%x1 = real(xv(:,edgei(1,e)),c_float)
         obj%cylflat(obj%ncylflat)%x2 = real(xv(:,edgei(2,e)),c_float)
         obj%cylflat(obj%ncylflat)%r = real(rade,c_float)
         obj%cylflat(obj%ncylflat)%rgb = rgbe
         obj%cylflat(obj%ncylflat)%x1delta = cmplx(dv(:,edgei(1,e)),kind=c_float_complex)
         obj%cylflat(obj%ncylflat)%x2delta = cmplx(dv(:,edgei(2,e)),kind=c_float_complex)
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
      type(dl_triangle), allocatable :: auxtri(:)

      obj%ntriangle = obj%ntriangle + 1
      if (obj%ntriangle > size(obj%triangle,1)) then
         allocate(auxtri(2*obj%ntriangle))
         auxtri(1:size(obj%triangle,1)) = obj%triangle
         call move_alloc(auxtri,obj%triangle)
      end if
      obj%triangle(obj%ntriangle)%x1 = real(p1,c_float)
      obj%triangle(obj%ntriangle)%x2 = real(p2,c_float)
      obj%triangle(obj%ntriangle)%x3 = real(p3,c_float)
      obj%triangle(obj%ntriangle)%rgb = rgb_
      obj%triangle(obj%ntriangle)%alpha = real(alpha_,c_float)
      obj%triangle(obj%ntriangle)%x1delta = cmplx(d1,kind=c_float_complex)
      obj%triangle(obj%ntriangle)%x2delta = cmplx(d2,kind=c_float_complex)
      obj%triangle(obj%ntriangle)%x3delta = cmplx(d3,kind=c_float_complex)

    end subroutine append_triangle

    !> Vibration-animation displacement of cell atom iat whose (animated)
    !> fractional position is xfrac; zero if the scene is not animating.
    function vibdelta(iat,xfrac) result(dv)
      integer, intent(in) :: iat
      real*8, intent(in) :: xfrac(3)
      complex*16 :: dv(3)

      real*8 :: vmass, vph

      dv = 0d0
      if (.not.doanim_) return
      vmass = atmass(sys(r%id)%c%spc(sys(r%id)%c%atcel(iat)%is)%z)
      vph = tpi * dot_product(xfrac,sys(r%id)%c%vib%qpt(:,iqpt))
      dv = sys(r%id)%c%vib%vec(:,iat,ifreq,iqpt) * exp(img * vph) / sqrt(vmass)

    end function vibdelta

    !> Append a string record to an allocatable string list,
    !> reallocating if necessary.
    subroutine append_str(lst,n,it)
      type(dl_string), allocatable, intent(inout) :: lst(:)
      integer, intent(inout) :: n
      type(dl_string), intent(in) :: it
      type(dl_string), allocatable :: aux(:)

      n = n + 1
      if (n > size(lst,1)) then
         allocate(aux(2*n))
         aux(1:size(lst,1)) = lst
         call move_alloc(aux,lst)
      end if
      lst(n) = it

    end subroutine append_str

    subroutine increase_ncylflat()

      obj%ncylflat = obj%ncylflat + 1
      if (obj%ncylflat > size(obj%cylflat,1)) then
         allocate(auxcyl(2*obj%ncylflat))
         auxcyl(1:size(obj%cylflat,1)) = obj%cylflat
         call move_alloc(auxcyl,obj%cylflat)
      end if

    end subroutine increase_ncylflat

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
               ucini(j) = sys(r%id)%c%vactop(j) - 1d0 - vacextension / sys(r%id)%c%aa(j)
               ucend(j) = sys(r%id)%c%vacbot(j) + vacextension / sys(r%id)%c%aa(j)
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
                  ucini(j) = 0.5d0 * (sys(r%id)%c%vactop(j) - 1d0 + sys(r%id)%c%vacbot(j))
                  ucend(j) = 0.5d0 * (sys(r%id)%c%vactop(j) - 1d0 + sys(r%id)%c%vacbot(j))
               end if
            end if
         end do
      end if

      ! stick ends
      x1 = ucini + r%origin
      x1 = sys(r%id)%c%x2c(x1)
      x2 = ucend + r%origin
      x2 = sys(r%id)%c%x2c(x2)

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
       call r%atom_style%reset(r)
    if (itype == 0 .or. itype == 2) &
       call r%bond_style%reset(r)
    if (itype == 0 .or. itype == 3) &
       call r%label_style%reset(r)
    if (itype == 0 .or. itype == 4) &
       call r%mol_style%reset(r)
    if (itype == 0 .or. itype == 1) &
       call r%coordpoly_style%reset(r)

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
       if (r%atom_color_type == 0) then
          d%rgb(:,i) = ColorElement(:,iz)
       elseif (r%atom_color_type == 1) then
          d%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
       else
          d%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
       end if

       ! scale covalent, vdw, or absolute value
       if (r%atom_radii_type == 0) then
          d%rad(i) = r%atom_radii_scale * atmcov(iz)
       elseif (r%atom_radii_type == 1) then
          d%rad(i) = r%atom_radii_scale * atmvdw(iz)
       else
          d%rad(i) = r%atom_radii_value
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
    call sys(r%id)%c%find_asterisms(d%nstar,atmrad=r%bond_atmrad,bondfac=r%bond_bfactor,&
       bonddelta=r%bond_bdelta)

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
    select case(r%label_type)
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
       select case(r%label_type)
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
       if (r%label_type == 0) then ! 0 = atomic symbol
          d%str(i) = trim(nameguess(sys(r%id)%c%spc(i)%z,.true.))
       elseif (r%label_type == 1) then ! 1 = atom name
          d%str(i) = trim(sys(r%id)%c%at(i)%name)
       elseif (r%label_type == 6) then ! 6 = Z
          d%str(i) = string(sys(r%id)%c%spc(i)%z)
       elseif (r%label_type == 8) then ! 8 = wyckoff
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
