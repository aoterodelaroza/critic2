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
  !> itype is the representation type, style is the scene style
  !> flavor is the representation flavor, and icount is the
  !> count array of the calling scene.
  module subroutine representation_init(r,isys,irep,itype,style,flavor,icount)
    use systems, only: sys_ready, ok_system
    use tools_io, only: string
    class(representation), intent(inout) :: r
    integer, intent(in) :: isys
    integer, intent(in) :: irep
    integer, intent(in) :: itype
    integer, intent(in) :: style
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
          r%name = "VdW contacts"
          r%atoms_display = .false.
       elseif (flavor == repflavor_atoms_hbonds) then
          r%name = "Hydrogen bonds"
          r%atoms_display = .false.
       end if
    elseif (itype == reptype_unitcell) then
       r%isinit = .true.
       r%shown = .true.
       r%name = "Unit Cell"
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
    call r%set_defaults(style,0)

  end subroutine representation_init

  !> Set all values to default for the representation with the given
  !> scene style. Set a subset of defaults if itype = 0 (all),
  !> 1 (atom), 2 (bonds), 3 (labels), 4 (mol), 5 (unit cell).
  module subroutine representation_set_defaults(r,style,itype)
    use scenes, only: style_phong
    use systems, only: sys, sys_ready, ok_system
    use global, only: bondfactor
    class(representation), intent(inout) :: r
    integer, intent(in) :: style
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
       r%atom_radii_value = bondrad_def
       r%atom_color_type = 0
       r%atom_border_size = atomborder_def
       r%atom_border_rgb = (/0._c_float,0._c_float,0._c_float/)
       if (r%flavor == repflavor_atoms_licorice) then
          r%atom_radii_type = 2
          r%atom_radii_value = atomrad_licorice_def
       end if
    end if

    !--> bonds
    if (itype == 0 .or. itype == 2) then
       r%bond_distancetype = 0
       r%bond_dmin = 0._c_float
       r%bond_dmax = 0._c_float
       r%bond_bfmin = 0._c_float
       r%bond_bfmax = real(bondfactor,c_float)
       r%bond_radtype = (/0_c_int,0_c_int/)
       r%bond_color_style = 0
       r%bond_border_size = bondborder_def
       r%bond_rad = bondrad_def
       r%bond_border_rgb = (/0._c_float,0._c_float,0._c_float/)
       r%bond_rgb = (/0._c_float,0._c_float,0._c_float/)
       r%bond_order = 1
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
       r%label_scale = 0.5_c_float
       r%label_rgb = (/0._c_float,0._c_float,0._c_float/)
       r%label_const_size = .false.
       r%label_offset = (/0._c_float,0._c_float,0._c_float/)
    end if

    !--> molecules

    ! unit cell
    if (itype == 0 .or. itype == 5) then
       r%uc_inner = .true.
       r%uc_coloraxes = .true.
       r%uc_vaccutsticks = .true.
       r%uc_radius = uc_radius_def
       r%uc_radiusinner = uc_radiusinner_def
       if (style == style_phong) then
          r%uc_rgb = 1._c_float
       else
          r%uc_rgb = 0._c_float
       end if
       r%uc_innersteplen = uc_innersteplen_def
       r%uc_innerstipple = .true.
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

    ! labels
    doreset = .not.r%label_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%label_style%timelastreset)
    if (doreset) call r%label_style%reset(r)

  end subroutine update_styles

  !> Add the spheres, cylinder, etc. to the draw lists. Use nc number
  !> of cells and the data from representation r. If doanim, use qpt
  !> iqpt and frequency ifreq to animate the representation.
  module subroutine add_draw_elements(r,nc,obj,doanim,iqpt,ifreq)
    use systems, only: sys
    use crystalmod, only: iperiod_vacthr
    use tools_io, only: string, nameguess
    use param, only: tpi, img, atmass
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
    integer :: ib, ineigh, ixn(3), ix1(3), ix2(3), nstep, vacshift(3)
    real(c_float) :: rgb(3)
    real*8 :: rad1, rad2, dd, f1, f2
    real*8 :: xx(3), xc(3), x0(3), x1(3), x2(3), res, uoriginc(3), phase, mass
    real*8 :: ucini(3), ucend(3)
    complex*16 :: xdelta0(3), xdelta1(3), xdelta2(3)
    type(dl_sphere), allocatable :: auxsph(:)
    type(dl_cylinder), allocatable :: auxcyl(:)
    type(dl_string), allocatable :: auxstr(:)
    character(len=:), allocatable :: errmsg

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

       ! whether we'll be doing bonds, allocate array to check whether
       ! an atoms has been drawn
       dobonds = r%bonds_display .and. r%bond_style%isinit
       if (dobonds) then
          allocate(lshown(sys(r%id)%c%ncel,-1:n(1)+1,-1:n(2)+1,-1:n(3)+1))
          lshown = .false.
       end if

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
          if (r%atom_style%type == 0) then ! species
             id = sys(r%id)%c%atcel(i)%is
          elseif (r%atom_style%type == 1) then ! nneq
             id = sys(r%id)%c%atcel(i)%idx
          else ! ncel
             id = i
          end if
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

                   ! calculate the animation delta
                   xdelta1 = 0d0
                   if (doanim_) then
                      mass = atmass(sys(r%id)%c%spc(sys(r%id)%c%atcel(i)%is)%z)
                      phase = tpi * dot_product(xx,sys(r%id)%c%vib%qpt(:,iqpt))
                      xdelta1 = sys(r%id)%c%vib%vec(:,i,ifreq,iqpt) * exp(img * phase) / sqrt(mass)
                   end if

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
                      obj%sph(obj%nsph)%border = r%atom_border_size
                      obj%sph(obj%nsph)%rgbborder = r%atom_border_rgb
                   end if

                   ! bonds
                   if (dobonds) then
                      ! mark this atom as drawn
                      call check_lshown(i,ix(1),ix(2),ix(3))
                      lshown(i,ix(1),ix(2),ix(3)) = .true.

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

                         ! calculate the animation delta of the other end
                         xdelta2 = 0d0
                         if (doanim_) then
                            mass = atmass(sys(r%id)%c%spc(sys(r%id)%c%atcel(ineigh)%is)%z)
                            xx = sys(r%id)%c%atcel(ineigh)%x + ixn
                            phase = tpi * dot_product(xx,sys(r%id)%c%vib%qpt(:,iqpt))
                            xdelta2 = sys(r%id)%c%vib%vec(:,ineigh,ifreq,iqpt) * exp(img * phase) / sqrt(mass)
                         end if

                         x1 = xc + uoriginc
                         x2 = sys(r%id)%c%atcel(ineigh)%x + ixn
                         x2 = sys(r%id)%c%x2c(x2) + uoriginc
                         if (r%bond_color_style == 0) then
                            obj%cyl(obj%ncyl)%x1 = real(x1,c_float)
                            obj%cyl(obj%ncyl)%x1delta = cmplx(xdelta1,kind=c_float_complex)
                            obj%cyl(obj%ncyl)%x2 = real(x2,c_float)
                            obj%cyl(obj%ncyl)%x2delta = cmplx(xdelta2,kind=c_float_complex)
                            obj%cyl(obj%ncyl)%r = r%bond_rad
                            obj%cyl(obj%ncyl)%rgb = r%bond_rgb
                            obj%cyl(obj%ncyl)%order = r%bond_order
                            obj%cyl(obj%ncyl)%border = r%bond_border_size
                            obj%cyl(obj%ncyl)%rgbborder = r%bond_border_rgb
                         else
                            ! calculate the midpoint, taking into account the atomic radii
                            if (r%atom_style%type == 0) then ! species
                               idaux = sys(r%id)%c%atcel(ineigh)%is
                            elseif (r%atom_style%type == 1) then ! nneq
                               idaux = sys(r%id)%c%atcel(ineigh)%idx
                            else ! ncel
                               idaux = ineigh
                            end if
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
                            obj%cyl(obj%ncyl-1)%r = r%bond_rad
                            obj%cyl(obj%ncyl-1)%rgb = rgb
                            obj%cyl(obj%ncyl-1)%order = r%bond_order
                            obj%cyl(obj%ncyl-1)%border = r%bond_border_size
                            obj%cyl(obj%ncyl-1)%rgbborder = r%bond_border_rgb

                            obj%cyl(obj%ncyl)%x1 = real(x0,c_float)
                            obj%cyl(obj%ncyl)%x1delta = cmplx(xdelta0,kind=c_float_complex)
                            obj%cyl(obj%ncyl)%x2 = real(x2,c_float)
                            obj%cyl(obj%ncyl)%x2delta = cmplx(xdelta2,kind=c_float_complex)
                            obj%cyl(obj%ncyl)%r = r%bond_rad
                            obj%cyl(obj%ncyl)%rgb = r%atom_style%rgb(:,idaux) * &
                               r%mol_style%tint_rgb(:,sys(r%id)%c%idatcelmol(1,ineigh))
                            obj%cyl(obj%ncyl)%order = r%bond_order
                            obj%cyl(obj%ncyl)%border = r%bond_border_size
                            obj%cyl(obj%ncyl)%rgbborder = r%bond_border_rgb
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
                            obj%string(obj%nstring)%scale = r%label_scale
                         else
                            obj%string(obj%nstring)%scale = -r%label_scale
                         end if
                         obj%string(obj%nstring)%offset = r%label_offset
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
          obj%cylflat(obj%ncylflat)%r = r%uc_radius
          if (r%uc_coloraxes.and.i==1) then
             obj%cylflat(obj%ncylflat)%rgb = (/1._c_float,0._c_float,0._c_float/)
          elseif (r%uc_coloraxes.and.i==2) then
             obj%cylflat(obj%ncylflat)%rgb = (/0._c_float,1._c_float,0._c_float/)
          elseif (r%uc_coloraxes.and.i==3) then
             obj%cylflat(obj%ncylflat)%rgb = (/0._c_float,0._c_float,1._c_float/)
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
                            obj%cylflat(obj%ncylflat)%r = r%uc_radiusinner
                            obj%cylflat(obj%ncylflat)%rgb = r%uc_rgb
                         end do
                      else
                         call increase_ncylflat()
                         obj%cylflat(obj%ncylflat)%x1 = real(x1 ,c_float)
                         obj%cylflat(obj%ncylflat)%x2 = real(x2 ,c_float)
                         obj%cylflat(obj%ncylflat)%r = r%uc_radiusinner
                         obj%cylflat(obj%ncylflat)%rgb = r%uc_rgb
                      end if
                   end do
                end do
             end do
          end do
       end if
    end if ! reptype
  contains
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

  end subroutine reset_all_styles

  !> Reset atom style to the parameters and the contents of the system
  !> point at by representation r. Uses d%type to fill the arrays.
  module subroutine atom_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use gui_main, only: ColorElement
    use param, only: atmcov, atmvdw, jmlcol, jmlcol2
    class(atom_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, ispc, iz

    ! if not initialized, set type
    if (.not.d%isinit) d%type = 0

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

    ! fill according to the style
    if (d%type == 0) then ! species
       d%ntype = sys(r%id)%c%nspc
       allocate(d%shown(d%ntype),d%rgb(3,d%ntype))
       allocate(d%rad(d%ntype))
       do i = 1, d%ntype
          iz = sys(r%id)%c%spc(i)%z
          if (r%atom_color_type == 0) then
             d%rgb(:,i) = ColorElement(:,iz)
          elseif (r%atom_color_type == 1) then
             d%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
          else
             d%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
          end if
          if (r%atom_radii_type == 0) then
             d%rad(i) = r%atom_radii_scale * real(atmcov(iz),c_float)
          elseif (r%atom_radii_type == 1) then
             d%rad(i) = r%atom_radii_scale * real(atmvdw(iz),c_float)
          else
             d%rad(i) = r%atom_radii_value
          endif
       end do
    elseif (d%type == 1) then ! nneq
       d%ntype = sys(r%id)%c%nneq
       allocate(d%shown(d%ntype),d%rgb(3,d%ntype))
       allocate(d%rad(d%ntype))
       do i = 1, sys(r%id)%c%nneq
          ispc = sys(r%id)%c%at(i)%is
          iz = sys(r%id)%c%spc(ispc)%z
          if (r%atom_color_type == 0) then
             d%rgb(:,i) = ColorElement(:,iz)
          elseif (r%atom_color_type == 1) then
             d%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
          else
             d%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
          end if
          if (r%atom_radii_type == 0) then
             d%rad(i) = r%atom_radii_scale * real(atmcov(iz),c_float)
          elseif (r%atom_radii_type == 1) then
             d%rad(i) = r%atom_radii_scale * real(atmvdw(iz),c_float)
          else
             d%rad(i) = r%atom_radii_value
          endif
       end do
    else ! ncel
       d%ntype = sys(r%id)%c%ncel
       allocate(d%shown(d%ntype),d%rgb(3,d%ntype))
       allocate(d%rad(d%ntype))
       do i = 1, sys(r%id)%c%ncel
          ispc = sys(r%id)%c%atcel(i)%is
          iz = sys(r%id)%c%spc(ispc)%z
          if (r%atom_color_type == 0) then
             d%rgb(:,i) = ColorElement(:,iz)
          elseif (r%atom_color_type == 1) then
             d%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
          else
             d%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
          end if
          if (r%atom_radii_type == 0) then
             d%rad(i) = r%atom_radii_scale * real(atmcov(iz),c_float)
          elseif (r%atom_radii_type == 1) then
             d%rad(i) = r%atom_radii_scale * real(atmvdw(iz),c_float)
          else
             d%rad(i) = r%atom_radii_value
          endif
       end do
    end if
    d%shown = .true.
    d%isinit = .true.

  end subroutine atom_style_reset

  !> Reset colors in an atom style to defaults.
  module subroutine atom_style_reset_colors(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use gui_main, only: ColorElement
    class(atom_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, ispc, iz

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! fill according to the style
    if (d%type == 0) then ! species
       do i = 1, d%ntype
          iz = sys(r%id)%c%spc(i)%z
          d%rgb(:,i) = ColorElement(:,iz)
       end do
    elseif (d%type == 1) then ! nneq
       do i = 1, sys(r%id)%c%nneq
          ispc = sys(r%id)%c%at(i)%is
          iz = sys(r%id)%c%spc(ispc)%z
          d%rgb(:,i) = ColorElement(:,iz)
       end do
    else ! ncel
       do i = 1, sys(r%id)%c%ncel
          ispc = sys(r%id)%c%atcel(i)%is
          iz = sys(r%id)%c%spc(ispc)%z
          d%rgb(:,i) = ColorElement(:,iz)
       end do
    end if

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
       d%scale_rad(i) = 1._c_float
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
    use param, only: atmcov, atmvdw
    class(bond_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, j
    real*8 :: r1cov, r1vdw, r2cov, r2vdw
    real*8, allocatable :: rij_t(:,:,:)

    ! check all the info is available
    if (.not.d%isinit) return
    if (.not.ok_system(r%id,sys_ready)) return

    ! allocate temporary data for rij table
    allocate(rij_t(sys(r%id)%c%nspc,2,sys(r%id)%c%nspc))

    ! fill table data
    do i = 1, sys(r%id)%c%nspc
       r1cov = atmcov(sys(r%id)%c%spc(i)%z)
       r1vdw = atmvdw(sys(r%id)%c%spc(i)%z)
       do j = i, sys(r%id)%c%nspc
          r2cov = atmcov(sys(r%id)%c%spc(j)%z)
          r2vdw = atmvdw(sys(r%id)%c%spc(j)%z)
          if (r%bond_distancetype == 0) then
             if (r%bond_radtype(1) == 0) then
                rij_t(i,1,j) = r%bond_bfmin * (r1cov + r2cov)
             else
                rij_t(i,1,j) = r%bond_bfmin * (r1vdw + r2vdw)
             end if
             if (r%bond_radtype(2) == 0) then
                rij_t(i,2,j) = r%bond_bfmax * (r1cov + r2cov)
             else
                rij_t(i,2,j) = r%bond_bfmax * (r1vdw + r2vdw)
             end if
          else
             rij_t(i,1,j) = r%bond_dmin / bohrtoa
             rij_t(i,2,j) = r%bond_dmax / bohrtoa
          end if
          rij_t(j,:,i) = rij_t(i,:,j)
       end do
    end do

    ! generate the new neighbor star
    call sys(r%id)%c%find_asterisms(d%nstar,rij=rij_t)

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

end submodule proc
