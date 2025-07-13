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
  implicit none


contains

  !> Initialize a representation for system isys with ID irep.
  !> itype is the representation type, style is the scene style
  !> flavor is the representation flavor, and icount is the
  !> count array of the calling scene.
  module subroutine representation_init(r,isys,irep,itype,style,flavor,icount)
    use scenes, only: style_phong
    use systems, only: sys, sys_ready, ok_system
    use tools_io, only: string
    class(representation), intent(inout), target :: r
    integer, intent(in) :: isys
    integer, intent(in) :: irep
    integer, intent(in) :: itype
    integer, intent(in) :: style
    integer, intent(in) :: flavor
    integer, intent(inout) :: icount(0:repflavor_NUM)

    ! check the system is sane
    if (.not.ok_system(isys,sys_ready)) return

    ! common settings
    r%isinit = .false.
    r%shown = .false.
    r%type = reptype_none
    r%flavor = repflavor_unknown
    r%id = isys
    r%idrep = irep
    r%name = ""
    r%filter = ""
    r%errfilter = ""
    r%pertype = 1
    r%ncell = 1
    r%border = .true.
    r%onemotif = .false.
    r%atoms_display = .true.
    r%bonds_display = .true.
    r%labels_display = .false.
    r%atom_radii_reset_type = 0
    r%atom_radii_reset_scale = 0.7_c_float
    r%atom_color_reset_type = 0
    r%uc_radius = 0.15_c_float
    r%uc_radiusinner = 0.15_c_float
    r%uc_innersteplen = 2d0
    r%uc_innerstipple = .true.
    r%uc_inner = .true.
    r%uc_coloraxes = .true.
    r%origin = 0._c_float
    r%tshift = 0._c_float

    ! style-dependent settings
    if (style == style_phong) then
       r%uc_rgb = 1._c_float
    else
       r%uc_rgb = 0._c_float
    end if

    ! type-dependent settings
    if (itype == reptype_atoms) then
       r%name = "Ball and Stick"
       r%isinit = .true.
       r%shown = .true.
       r%type = reptype_atoms
       if (sys(isys)%c%ismolecule) then
          r%ncell = 0
          r%border = .false.
          r%onemotif = .false.
       else
          r%border = .true.
          r%onemotif = (sys(isys)%c%nmol > 1)
          r%ncell = 1
       end if
    elseif (itype == reptype_unitcell) then
       r%isinit = .true.
       r%shown = .true.
       r%type = reptype_unitcell
       r%name = "Unit Cell"
    end if
    r%flavor = flavor

    ! increment type counter and set name
    icount(flavor) = icount(flavor) + 1
    if (icount(flavor) > 1) then
       r%name = trim(r%name) // "/" // string(icount(flavor))
    end if

    ! increment global counter
    icount(0) = icount(0) + 1
    r%iord = icount(0)

    ! apply flavors, global options
    if (flavor == repflavor_atoms_vdwcontacts) then
       r%name = "VdW contacts"
       r%atoms_display = .false.
       r%bonds_display = .true.
       r%labels_display = .false.
    elseif (flavor == repflavor_atoms_hbonds) then
       r%name = "Hydrogen bonds"
       r%atoms_display = .false.
       r%bonds_display = .true.
       r%labels_display = .false.
    end if

    ! initialize the styles
    call r%reset_all_styles()

  end subroutine representation_init

  !> Reset the representation to the default values.
  module subroutine representation_reset(r)
    use systems, only: sys, sys_ready, ok_system
    use tools_io, only: string
    class(representation), intent(inout), target :: r

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! common settings
    r%isinit = .false.
    r%filter = ""
    r%errfilter = ""
    r%pertype = 1
    r%ncell = 1
    r%border = .true.
    r%onemotif = .false.
    r%atoms_display = .true.
    r%bonds_display = .true.
    r%labels_display = .false.
    r%atom_radii_reset_type = 0
    r%atom_radii_reset_scale = 0.7_c_float
    r%atom_color_reset_type = 0
    r%uc_radius = 0.15_c_float
    r%uc_radiusinner = 0.15_c_float
    r%uc_innersteplen = 2d0
    r%uc_innerstipple = .true.
    r%uc_inner = .true.
    r%uc_coloraxes = .true.
    r%origin = 0._c_float
    r%tshift = 0._c_float
    r%uc_rgb = 0._c_float

    ! type-dependent settings
    if (r%type == reptype_atoms) then
       r%isinit = .true.
       if (sys(r%id)%c%ismolecule) then
          r%ncell = 0
          r%border = .false.
          r%onemotif = .false.
       else
          r%border = .true.
          r%onemotif = (sys(r%id)%c%nmol > 1)
          r%ncell = 1
       end if
    elseif (r%type == reptype_unitcell) then
       r%isinit = .true.
    end if

    ! apply flavors, global options
    if (r%flavor == repflavor_atoms_vdwcontacts) then
       r%atoms_display = .false.
       r%bonds_display = .true.
       r%labels_display = .false.
    elseif (r%flavor == repflavor_atoms_hbonds) then
       r%atoms_display = .false.
       r%bonds_display = .true.
       r%labels_display = .false.
    end if

    ! initialize the styles
    call r%reset_all_styles()

  end subroutine representation_reset

  !> Terminate a representation
  module subroutine representation_end(r)
    class(representation), intent(inout), target :: r

    r%name = ""
    r%filter = ""
    r%errfilter = ""
    r%isinit = .false.
    r%shown = .false.
    r%type = reptype_none
    r%flavor = repflavor_unknown
    r%id = 0
    r%idrep = 0
    r%iord = 0

  end subroutine representation_end

  !> Update the representation to respond to a change in the number
  !> of atoms or molecules in the associated system.
  module subroutine update_structure(r)
    use systems, only: sys_ready, ok_system, sysc
    class(representation), intent(inout), target :: r

    logical :: doreset

    ! consistency checks
    if (.not.r%isinit .or. r%id == 0) return
    if (.not.ok_system(r%id,sys_ready)) return

    ! check if we need to reset the representation styles
    ! atoms
    doreset = .not.r%atom_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%atom_style%timelastreset)
    if (doreset) call r%atom_style%reset(r%id)

    ! bonds: if the geometry changed
    doreset = doreset .or. .not.r%bond_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%bond_style%timelastreset)
    if (doreset) call r%bond_style%reset(r%id,r%flavor)

    ! bonds: if the system has been rebonded and this representation tracks
    ! the bonds in the system (%isdef), recalculate the bond style
    doreset = r%bond_style%isdef .and. (sysc(r%id)%timelastchange_rebond > r%bond_style%timelastreset)
    if (doreset) call r%bond_style%copy_neighstars_from_system(r%id)

    ! molecules: if the geometry or the bonds changed
    doreset = .not.r%mol_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_rebond > r%mol_style%timelastreset)
    if (doreset) call r%mol_style%reset(r%id)

    ! labels
    doreset = .not.r%label_style%isinit
    doreset = doreset .or. (sysc(r%id)%timelastchange_geometry > r%label_style%timelastreset)
    if (doreset) call r%label_style%reset(r%id)

  end subroutine update_structure

  !> Add the spheres, cylinder, etc. to the draw lists. Use nc number
  !> of cells and the data from representation r. If doanim, use qpt
  !> iqpt and frequency ifreq to animate the representation.
  module subroutine add_draw_elements(r,nc,nsph,drawlist_sph,ncyl,drawlist_cyl,&
     ncylflat,drawlist_cylflat,nstring,drawlist_string,doanim,iqpt,ifreq)
    use systems, only: sys
    use tools_io, only: string, nameguess
    use param, only: bohrtoa, tpi, img, atmass
    class(representation), intent(inout), target :: r
    integer, intent(in) :: nc(3)
    integer, intent(inout) :: nsph
    type(dl_sphere), intent(inout), allocatable :: drawlist_sph(:)
    integer, intent(inout) :: ncyl
    type(dl_cylinder), intent(inout), allocatable :: drawlist_cyl(:)
    integer, intent(inout) :: ncylflat
    type(dl_cylinder), intent(inout), allocatable :: drawlist_cylflat(:)
    integer, intent(inout) :: nstring
    type(dl_string), intent(inout), allocatable :: drawlist_string(:)
    logical, intent(in) :: doanim
    integer, intent(in) :: iqpt, ifreq

    logical, allocatable :: lshown(:,:,:,:)
    logical :: havefilter, step, isedge(3), usetshift, doanim_, dobonds
    integer :: n(3), i, j, k, imol, lvec(3), id, idaux, n0(3), n1(3)
    integer :: i1, i2, i3, ix(3), idl
    integer :: ib, ineigh, ixn(3), ix1(3), ix2(3), nstep
    real(c_float) :: rgb(3)
    real*8 :: rad1, rad2, dd, f1, f2
    real*8 :: xx(3), xc(3), x0(3), x1(3), x2(3), res, uoriginc(3), phase, mass
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
    if (.not.allocated(drawlist_sph)) then
       allocate(drawlist_sph(100))
       nsph = 0
    end if
    if (.not.allocated(drawlist_cyl)) then
       allocate(drawlist_cyl(100))
       ncyl = 0
    end if
    if (.not.allocated(drawlist_cylflat)) then
       allocate(drawlist_cylflat(100))
       ncylflat = 0
    end if
    if (.not.allocated(drawlist_string)) then
       allocate(drawlist_string(100))
       nstring = 0
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
          uoriginc = sys(r%id)%c%x2c(real(r%origin,8))
       end if

       ! whether we'll be doing bonds, allocate array to check whether
       ! an atoms has been drawn
       dobonds = r%bonds_display .and. r%bond_style%isinit
       if (dobonds) then
          allocate(lshown(sys(r%id)%c%ncel,-1:n(1)+1,-1:n(2)+1,-1:n(3)+1))
          lshown = .false.
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
          n0 = 0
          n1 = n-1
          if (r%border.and..not.r%onemotif) then
             xx = sys(r%id)%c%atcel(i)%x
             do j = 1, 3
                if (xx(j) < rthr) then
                   n1(j) = n(j)
                elseif (xx(j) > rthr1) then
                   n0(j) = -1
                end if
             end do
          end if

          ! draw the spheres and cylinders
          rgb = r%atom_style%rgb(:,id) * r%mol_style%tint_rgb(:,imol)
          rad1 = r%atom_style%rad(id) * r%mol_style%scale_rad(imol)
          do i1 = n0(1), n1(1)
             do i2 = n0(2), n1(2)
                do i3 = n0(3), n1(3)
                   ix = (/i1,i2,i3/) + lvec
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
                      nsph = nsph + 1
                      if (nsph > size(drawlist_sph,1)) then
                         allocate(auxsph(2*nsph))
                         auxsph(1:size(drawlist_sph,1)) = drawlist_sph
                         call move_alloc(auxsph,drawlist_sph)
                      end if

                      ! write down the sphere
                      drawlist_sph(nsph)%x = real(xc + uoriginc,c_float)
                      drawlist_sph(nsph)%r = real(rad1,c_float)
                      drawlist_sph(nsph)%rgb = rgb
                      drawlist_sph(nsph)%idx(1) = i
                      drawlist_sph(nsph)%idx(2:4) = ix
                      drawlist_sph(nsph)%xdelta = cmplx(xdelta1,kind=c_float_complex)
                      drawlist_sph(nsph)%border = r%atom_style%border_size
                      drawlist_sph(nsph)%rgbborder = r%atom_style%rgbborder
                   end if

                   ! bonds
                   if (dobonds) then
                      ! mark this atom as drawn
                      call check_lshown(i,ix(1),ix(2),ix(3))
                      lshown(i,ix(1),ix(2),ix(3)) = .true.

                      ! bonds
                      do ib = 1, r%bond_style%nstar(i)%ncon
                         ineigh = r%bond_style%nstar(i)%idcon(ib)
                         if (.not.r%bond_style%shown_g(sys(r%id)%c%atcel(ineigh)%is,sys(r%id)%c%atcel(i)%is)) cycle
                         ixn = ix + r%bond_style%nstar(i)%lcon(:,ib)

                         if (r%bond_style%imol_g == 1) then ! intramol
                            if (.not.sys(r%id)%c%in_same_molecule(i,ix,ineigh,ixn)) cycle
                         elseif (r%bond_style%imol_g == 2) then ! intermol
                            if (sys(r%id)%c%in_same_molecule(i,ix,ineigh,ixn)) cycle
                         end if

                         call check_lshown(ineigh,ixn(1),ixn(2),ixn(3))
                         if (r%bond_style%bothends_g) then
                            ! skip if the atom has been represented already
                            ! (draws once, and only if both atoms are present)
                            if (.not.lshown(ineigh,ixn(1),ixn(2),ixn(3))) cycle
                         else
                            ! skip if the atom has not been represented already
                            ! (draws once, only one of the atoms need be present)
                            if (lshown(ineigh,ixn(1),ixn(2),ixn(3))) cycle
                         end if

                         ! draw the bond, reallocate if necessary
                         if (r%bond_style%style_g == 0) then
                            ncyl = ncyl + 1
                         else
                            ncyl = ncyl + 2
                         end if
                         if (ncyl > size(drawlist_cyl,1)) then
                            allocate(auxcyl(2*ncyl))
                            auxcyl(1:size(drawlist_cyl,1)) = drawlist_cyl
                            call move_alloc(auxcyl,drawlist_cyl)
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
                         if (r%bond_style%style_g == 0) then
                            drawlist_cyl(ncyl)%x1 = real(x1,c_float)
                            drawlist_cyl(ncyl)%x1delta = cmplx(xdelta1,kind=c_float_complex)
                            drawlist_cyl(ncyl)%x2 = real(x2,c_float)
                            drawlist_cyl(ncyl)%x2delta = cmplx(xdelta2,kind=c_float_complex)
                            drawlist_cyl(ncyl)%r = r%bond_style%rad_g
                            drawlist_cyl(ncyl)%rgb = r%bond_style%rgb_g
                            drawlist_cyl(ncyl)%order = r%bond_style%order_g
                            drawlist_cyl(ncyl)%border = r%bond_style%border_g
                            drawlist_cyl(ncyl)%rgbborder = r%bond_style%rgbborder_g
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
                            drawlist_cyl(ncyl-1)%x1 = real(x1,c_float)
                            drawlist_cyl(ncyl-1)%x1delta = cmplx(xdelta1,kind=c_float_complex)
                            drawlist_cyl(ncyl-1)%x2 = real(x0,c_float)
                            drawlist_cyl(ncyl-1)%x2delta = cmplx(xdelta0,kind=c_float_complex)
                            drawlist_cyl(ncyl-1)%r = r%bond_style%rad_g
                            drawlist_cyl(ncyl-1)%rgb = rgb
                            drawlist_cyl(ncyl-1)%order = r%bond_style%order_g
                            drawlist_cyl(ncyl-1)%border = r%bond_style%border_g
                            drawlist_cyl(ncyl-1)%rgbborder = r%bond_style%rgbborder_g

                            drawlist_cyl(ncyl)%x1 = real(x0,c_float)
                            drawlist_cyl(ncyl)%x1delta = cmplx(xdelta0,kind=c_float_complex)
                            drawlist_cyl(ncyl)%x2 = real(x2,c_float)
                            drawlist_cyl(ncyl)%x2delta = cmplx(xdelta2,kind=c_float_complex)
                            drawlist_cyl(ncyl)%r = r%bond_style%rad_g
                            drawlist_cyl(ncyl)%rgb = r%atom_style%rgb(:,idaux) * &
                               r%mol_style%tint_rgb(:,sys(r%id)%c%idatcelmol(1,ineigh))
                            drawlist_cyl(ncyl)%order = r%bond_style%order_g
                            drawlist_cyl(ncyl)%border = r%bond_style%border_g
                            drawlist_cyl(ncyl)%rgbborder = r%bond_style%rgbborder_g
                         end if
                      end do ! ncon
                   end if

                   if (r%labels_display) then
                      select case(r%label_style%style)
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
                         nstring = nstring + 1
                         if (nstring > size(drawlist_string,1)) then
                            allocate(auxstr(2*nstring))
                            auxstr(1:size(drawlist_string,1)) = drawlist_string
                            call move_alloc(auxstr,drawlist_string)
                         end if

                         drawlist_string(nstring)%x = real(xc + uoriginc,c_float)
                         drawlist_string(nstring)%xdelta = cmplx(xdelta1,kind=c_float_complex)
                         drawlist_string(nstring)%r = real(rad1,c_float)
                         drawlist_string(nstring)%rgb = r%label_style%rgb
                         if (r%label_style%const_size) then
                            drawlist_string(nstring)%scale = r%label_style%scale
                         else
                            drawlist_string(nstring)%scale = -r%label_style%scale
                         end if
                         drawlist_string(nstring)%offset = r%label_style%offset
                         drawlist_string(nstring)%str = trim(r%label_style%str(idl))
                         if (r%label_style%style == 3) then
                            ! add the lattice vectors
                            drawlist_string(nstring)%str = drawlist_string(nstring)%str // "[" //&
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

       ! external cell
       do i = 1, 12
          x1 = real(uc(:,1,i) * n,8) + r%origin
          x1 = sys(r%id)%c%x2c(x1)
          x2 = real(uc(:,2,i) * n,8) + r%origin
          x2 = sys(r%id)%c%x2c(x2)

          call increase_ncylflat()
          drawlist_cylflat(ncylflat)%x1 = real(x1,c_float)
          drawlist_cylflat(ncylflat)%x2 = real(x2,c_float)
          drawlist_cylflat(ncylflat)%r = r%uc_radius
          if (r%uc_coloraxes.and.i==1) then
             drawlist_cylflat(ncylflat)%rgb = (/1._c_float,0._c_float,0._c_float/)
          elseif (r%uc_coloraxes.and.i==2) then
             drawlist_cylflat(ncylflat)%rgb = (/0._c_float,1._c_float,0._c_float/)
          elseif (r%uc_coloraxes.and.i==3) then
             drawlist_cylflat(ncylflat)%rgb = (/0._c_float,0._c_float,1._c_float/)
          else
             drawlist_cylflat(ncylflat)%rgb = r%uc_rgb
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

                      x1 = real(ix1,8) + r%origin
                      x1 = sys(r%id)%c%x2c(x1)
                      x2 = real(ix2,8) + r%origin
                      x2 = sys(r%id)%c%x2c(x2)

                      ! logical :: uc_innerstipple ! stippled lines for the inner lines
                      if (r%uc_innerstipple) then
                         nstep = ceiling(norm2(x2 - x1) / r%uc_innersteplen)
                         do j = 1, nstep
                            call increase_ncylflat()
                            drawlist_cylflat(ncylflat)%x1 = real(x1 + real(2*j-1,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            drawlist_cylflat(ncylflat)%x2 = real(x1 + real(2*j,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            drawlist_cylflat(ncylflat)%r = r%uc_radiusinner
                            drawlist_cylflat(ncylflat)%rgb = r%uc_rgb
                         end do
                      else
                         call increase_ncylflat()
                         drawlist_cylflat(ncylflat)%x1 = real(x1 ,c_float)
                         drawlist_cylflat(ncylflat)%x2 = real(x2 ,c_float)
                         drawlist_cylflat(ncylflat)%r = r%uc_radiusinner
                         drawlist_cylflat(ncylflat)%rgb = r%uc_rgb
                      end if
                   end do
                end do
             end do
          end do
       end if
    end if ! reptype
  contains
    subroutine increase_ncylflat()

      ncylflat = ncylflat + 1
      if (ncylflat > size(drawlist_cylflat,1)) then
         allocate(auxcyl(2*ncylflat))
         auxcyl(1:size(drawlist_cylflat,1)) = drawlist_cylflat
         call move_alloc(auxcyl,drawlist_cylflat)
      end if

    end subroutine increase_ncylflat

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

  !> Reset all styles.
  module subroutine reset_all_styles(r)
    class(representation), intent(inout), target :: r

    call r%atom_style%reset(r%id)
    call r%mol_style%reset(r%id)
    call r%bond_style%reset(r%id,r%flavor)
    call r%label_style%reset(r%id)

  end subroutine reset_all_styles

  !> Reset atom style to defaults consistent with system isys, or empty
  !> if system is not ready.
  module subroutine reset_atom_style(d,isys)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use gui_main, only: ColorElement
    use param, only: atmcov
    class(draw_style_atom), intent(inout), target :: d
    integer, intent(in) :: isys

    integer :: i, ispc, iz

    ! if not initialized, set type
    if (.not.d%isinit) d%type = 0

    ! set the atom style to zero
    d%ntype = 0
    d%isinit = .false.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%rgb)) deallocate(d%rgb)
    if (allocated(d%rad)) deallocate(d%rad)
    d%border_size = atomborder_def
    d%rgbborder = 0._c_float

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(isys,sys_ready)) return

    ! fill according to the style
    if (d%type == 0) then ! species
       d%ntype = sys(isys)%c%nspc
       allocate(d%shown(d%ntype),d%rgb(3,d%ntype))
       allocate(d%rad(d%ntype))
       do i = 1, d%ntype
          iz = sys(isys)%c%spc(i)%z
          d%rgb(:,i) = ColorElement(:,iz)
          d%rad(i) = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    elseif (d%type == 1) then ! nneq
       d%ntype = sys(isys)%c%nneq
       allocate(d%shown(d%ntype),d%rgb(3,d%ntype))
       allocate(d%rad(d%ntype))
       do i = 1, sys(isys)%c%nneq
          ispc = sys(isys)%c%at(i)%is
          iz = sys(isys)%c%spc(ispc)%z
          d%rgb(:,i) = ColorElement(:,iz)
          d%rad(i) = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    else ! ncel
       d%ntype = sys(isys)%c%ncel
       allocate(d%shown(d%ntype),d%rgb(3,d%ntype))
       allocate(d%rad(d%ntype))
       do i = 1, sys(isys)%c%ncel
          ispc = sys(isys)%c%atcel(i)%is
          iz = sys(isys)%c%spc(ispc)%z
          d%rgb(:,i) = ColorElement(:,iz)
          d%rad(i) = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    end if
    d%shown = .true.
    d%isinit = .true.

  end subroutine reset_atom_style

  !> Reset colors in an atom style to defaults.
  module subroutine reset_colors_atom_style(d,isys)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use gui_main, only: ColorElement
    class(draw_style_atom), intent(inout), target :: d
    integer, intent(in) :: isys

    integer :: i, ispc, iz

    ! check the system is sane
    if (.not.ok_system(isys,sys_ready)) return

    ! fill according to the style
    if (d%type == 0) then ! species
       do i = 1, d%ntype
          iz = sys(isys)%c%spc(i)%z
          d%rgb(:,i) = ColorElement(:,iz)
       end do
    elseif (d%type == 1) then ! nneq
       do i = 1, sys(isys)%c%nneq
          ispc = sys(isys)%c%at(i)%is
          iz = sys(isys)%c%spc(ispc)%z
          d%rgb(:,i) = ColorElement(:,iz)
       end do
    else ! ncel
       do i = 1, sys(isys)%c%ncel
          ispc = sys(isys)%c%atcel(i)%is
          iz = sys(isys)%c%spc(ispc)%z
          d%rgb(:,i) = ColorElement(:,iz)
       end do
    end if

  end subroutine reset_colors_atom_style

  !> Reset molecule style with default values. Use the information in
  !> system isys, or leave it empty if isys = 0.
  module subroutine reset_mol_style(d,isys)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    class(draw_style_molecule), intent(inout), target :: d
    integer, intent(in) :: isys

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
    if (.not.ok_system(isys,sys_ready)) return

    ! fill
    d%ntype = sys(isys)%c%nmol
    allocate(d%shown(d%ntype),d%tint_rgb(3,d%ntype))
    allocate(d%scale_rad(d%ntype))
    do i = 1, sys(isys)%c%nmol
       d%tint_rgb(:,i) = 1._c_float
       d%scale_rad(i) = 1._c_float
    end do
    d%shown = .true.
    d%isinit = .true.

  end subroutine reset_mol_style

  !> Generate the neighbor stars from the data in the rij table using
  !> the geometry in system isys.
  module subroutine generate_neighstars(d,isys)
    use systems, only: sys, sys_ready, ok_system
    use param, only: bohrtoa, atmcov, atmvdw
    class(draw_style_bond), intent(inout), target :: d
    integer, intent(in) :: isys

    integer :: i, j
    real*8 :: r1cov, r1vdw, r2cov, r2vdw
    real*8, allocatable :: rij_t(:,:,:)

    ! check all the info is available
    if (.not.d%isinit) return
    if (.not.ok_system(isys,sys_ready)) return

    ! allocate temporary data for rij table
    allocate(rij_t(sys(isys)%c%nspc,2,sys(isys)%c%nspc))

    ! fill table data
    do i = 1, sys(isys)%c%nspc
       r1cov = atmcov(sys(isys)%c%spc(i)%z)
       r1vdw = atmvdw(sys(isys)%c%spc(i)%z)
       do j = i, sys(isys)%c%nspc
          r2cov = atmcov(sys(isys)%c%spc(j)%z)
          r2vdw = atmvdw(sys(isys)%c%spc(j)%z)
          if (d%distancetype_g == 0) then
             if (d%radtype_g(1) == 0) then
                rij_t(i,1,j) = d%bfmin_g * (r1cov + r2cov)
             else
                rij_t(i,1,j) = d%bfmin_g * (r1vdw + r2vdw)
             end if
             if (d%radtype_g(2) == 0) then
                rij_t(i,2,j) = d%bfmax_g * (r1cov + r2cov)
             else
                rij_t(i,2,j) = d%bfmax_g * (r1vdw + r2vdw)
             end if
          else
             rij_t(i,1,j) = d%dmin_g / bohrtoa
             rij_t(i,2,j) = d%dmax_g / bohrtoa
          end if
          rij_t(j,:,i) = rij_t(i,:,j)
       end do
    end do

    ! generate the new neighbor star
    call sys(isys)%c%find_asterisms(d%nstar,rij=rij_t)

  end subroutine generate_neighstars

  !> Copy the neighbor stars from the given system.
  module subroutine copy_neighstars_from_system(d,isys)
    use systems, only: sys, sys_ready, ok_system
    class(draw_style_bond), intent(inout), target :: d
    integer, intent(in) :: isys

    ! check all the info is available
    if (.not.d%isinit) return
    if (.not.ok_system(isys,sys_ready)) return

    ! copy the nstar
    if (allocated(sys(isys)%c%nstar)) &
       d%nstar = sys(isys)%c%nstar

  end subroutine copy_neighstars_from_system

  !> Reset bond style with default values, according to the given
  !> flavor. Use the information in system isys, or leave it empty if
  !> isys = 0.
  module subroutine reset_bond_style(d,isys,flavor)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use global, only: bondfactor
    class(draw_style_bond), intent(inout), target :: d
    integer, intent(in) :: isys
    integer, intent(in), optional :: flavor

    integer :: i, j, iz
    integer :: flavor_

    ! optional arguments
    flavor_ = repflavor_unknown
    if (present(flavor)) flavor_ = flavor

    ! clear the bond style
    d%isinit = .false.
    if (allocated(d%shown_g)) deallocate(d%shown_g)
    if (allocated(d%nstar)) deallocate(d%nstar)
    d%isdef = .true.

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(isys,sys_ready)) return
    d%isinit = .true.

    ! fill temp options
    allocate(d%shown_g(sys(isys)%c%nspc,sys(isys)%c%nspc))
    d%distancetype_g = 0
    d%dmin_g = 0._c_float
    d%dmax_g = 0._c_float
    d%bfmin_g = 0._c_float
    d%bfmax_g = real(bondfactor,c_float)
    d%radtype_g = 0
    d%style_g = 0
    d%rad_g = bond_rad_def
    d%border_g = atomborder_def
    d%rgbborder_g = 0._c_float
    d%rgb_g = 0._c_float
    d%order_g = 1
    d%imol_g = 0
    d%bothends_g = .true.
    d%shown_g = .true.

    ! fill data according to flavor
    if (flavor_ == repflavor_atoms_vdwcontacts) then
       ! van der waals contacts
       d%isdef = .false.
       d%border_g = 0._c_float
       d%rgbborder_g = 0._c_float
       d%bothends_g = .false.
       d%distancetype_g = 0_c_int
       d%bfmin_g = 0._c_float
       d%bfmax_g = 1._c_float
       d%radtype_g(2) = 1_c_int
       d%imol_g = 2_c_int
       d%rgb_g = (/0.51_c_float,0.83_c_float,0.11_c_float/)
       do i = 1, sys(isys)%c%nspc
          if (sys(isys)%c%spc(i)%z == 1) then
             d%shown_g(i,:) = .false.
             d%shown_g(:,i) = .false.
          end if
       end do
       d%rad_g = 0.15_c_float
       d%order_g = 0
       call d%generate_neighstars(isys)
    elseif (flavor_ == repflavor_atoms_hbonds) then
       ! hydrogen bonds
       d%isdef = .false.
       d%border_g = 0._c_float
       d%rgbborder_g = 0._c_float
       d%bothends_g = .false.
       d%distancetype_g = 0_c_int
       d%bfmin_g = 1.2_c_float
       d%radtype_g(1) = 0_c_int
       d%bfmax_g = 1._c_float
       d%radtype_g(2) = 1_c_int
       d%imol_g = 2_c_int
       d%rgb_g = (/0.11_c_float,0.44_c_float,0.83_c_float/)
       d%shown_g = .false.
       do i = 1, sys(isys)%c%nspc
          if (sys(isys)%c%spc(i)%z /= 1) cycle
          do j = 1, sys(isys)%c%nspc
             iz = sys(isys)%c%spc(j)%z
             if (iz == 7 .or. iz == 8 .or. iz == 9 .or. iz == 16 .or. iz == 17) then
                d%shown_g(i,j) = .true.
                d%shown_g(j,i) = .true.
             end if
          end do
       end do
       d%rad_g = 0.15_c_float
       d%order_g = 0
       call d%generate_neighstars(isys)
    else
       ! default flavor
       d%isdef = .true.
       call d%copy_neighstars_from_system(isys)
    end if

  end subroutine reset_bond_style

  !> Reset label style with default values. Use the information in
  !> system isys, or leave it empty if isys = 0.
  module subroutine reset_label_style(d,isys)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use tools_io, only: nameguess, string
    class(draw_style_label), intent(inout), target :: d
    integer, intent(in) :: isys

    integer :: i

    ! if not initialized, set type
    if (.not.d%isinit) d%style = 0

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(isys,sys_ready)) return

    ! set the atom style to defaults
    d%isinit = .true.
    d%scale = 0.5_c_float
    d%rgb = 0._c_float
    d%offset = 0._c_float
    d%const_size = .false.

    ! fill according to the style
    select case(d%style)
    case (0,5,6)
       d%ntype = sys(isys)%c%nspc
    case (2,3)
       d%ntype = sys(isys)%c%ncel
    case (1,4,8)
       d%ntype = sys(isys)%c%nneq
    case (7)
       d%ntype = sys(isys)%c%nmol
    end select
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%str)) deallocate(d%str)
    allocate(d%shown(d%ntype))
    allocate(d%str(d%ntype))

    ! fill shown, exclude hydrogens
    d%shown = .true.
    do i = 1, d%ntype
       select case(d%style)
       case (0,5,6)
          if (sys(isys)%c%spc(i)%z == 1) d%shown(i) = .false.
       case (2,3)
          if (sys(isys)%c%spc(sys(isys)%c%atcel(i)%is)%z == 1) d%shown(i) = .false.
       case (1,4,8)
          if (sys(isys)%c%spc(sys(isys)%c%at(i)%is)%z == 1) d%shown(i) = .false.
       end select
    end do

    ! fill text
    do i = 1, d%ntype
       if (d%style == 0) then ! 0 = atomic symbol
          d%str(i) = trim(nameguess(sys(isys)%c%spc(i)%z,.true.))
       elseif (d%style == 1) then ! 1 = atom name
          d%str(i) = trim(sys(isys)%c%at(i)%name)
       elseif (d%style == 6) then ! 6 = Z
          d%str(i) = string(sys(isys)%c%spc(i)%z)
       elseif (d%style == 8) then ! 8 = wyckoff
          d%str(i) = string(sys(isys)%c%at(i)%mult) //&
             string(sys(isys)%c%at(i)%wyc)
       else
          d%str(i) = string(i)
       end if
    end do

  end subroutine reset_label_style

end submodule proc
