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


! Options and style rows of the objects, to and from the JSON objects
! of the scene files.
submodule (representations) json
  implicit none

contains

  !> Transfer the settings of object r to (io%reading = .false.) or
  !> from (.true.) the JSON object p: its options if rows is .false.,
  !> or its style rows (per species, atom, molecule,...) if .true.
  !> The type and flavor of the object are the caller's. On reading,
  !> the options are read into an object made with the defaults of its
  !> kind, and the rows after its styles were made for its system
  !> (update): rows that do not have the size the system gives are
  !> not taken. Transient state (caches, editor selections, highlights
  !> under the mouse) and the links between objects are not
  !> transferred. Measurements anchored to critical points carry the
  !> CPs by content (field, type, position), and find them again in the
  !> CP list of the field when it is read.
  module subroutine representation_xfer_json(r,io,p,rows)
    use jsonxfer, only: json_xfer
    use json_module, only: json_value
    use systems, only: sys, sysc, atlisttype_NUM, cp_anchor_add, cp_anchor_resolve
    use global, only: new_uid
    use utils, only: iw_ncmap
    use display, only: disp_maxncell
    use types, only: id_mo_id
    class(representation), intent(inout) :: r
    type(json_xfer), intent(inout) :: io
    type(json_value), pointer, intent(in) :: p
    logical, intent(in) :: rows

    type(json_value), pointer :: q, a, b
    integer :: i, n, labelmax
    logical :: found

    if (rows) then
       call xfer_rows()
    else
       call xfer_options()
    end if

  contains

    !> The options of the object.
    subroutine xfer_options()

      ! common
      call io%xs(p,'name',r%name)
      call io%x(p,'shown',r%shown)
      call io%x(p,'pertype',r%disp%pertype,0,2)
      call io%x(p,'ncell',r%disp%ncell)
      if (io%reading) r%disp%ncell = max(min(r%disp%ncell,disp_maxncell),1)

      ! atoms (the atom-based kinds, and the atoms at the ends of
      ! gradient paths)
      if (reptype_is_atombased(r%type) .or. r%type == reptype_gpaths) then
         call io%obj(p,'atoms',q,found)
         if (found) then
            call io%x(q,'radii_type',r%atoms%radii_type,0,2)
            call io%x(q,'radii_scale',r%atoms%radii_scale,lo=0d0)
            call io%x(q,'radii_value',r%atoms%radii_value,lo=0d0)
            call io%x(q,'color_type',r%atoms%color_type,0,2)
            call io%x(q,'border_size',r%atoms%border_size,lo=0d0)
            call io%x(q,'border_rgb',r%atoms%border_rgb)
            call io%x(q,'occ_sectors',r%atoms%occ_sectors)
            call io%x(q,'occ_empty_rgb',r%atoms%occ_empty_rgb)
            call io%x(q,'fitpoly',r%atoms%fitpoly)
            call io%x(q,'style_type',r%atoms%style%type,1,atlisttype_NUM)
         end if
      end if

      select case (r%type)
      case (reptype_bonds)
         call io%obj(p,'bonds',q,found)
         if (.not.found) return
         call io%x(q,'atmrad',r%bonds%atmrad)
         call io%x(q,'bfactor',r%bonds%bfactor)
         call io%x(q,'bdelta',r%bonds%bdelta)
         call io%x(q,'color_style',r%bonds%color_style,0,1)
         call io%x(q,'rad',r%bonds%rad,lo=0d0)
         call io%x(q,'border_size',r%bonds%border_size,lo=0d0)
         call io%x(q,'border_rgb',r%bonds%border_rgb)
         call io%x(q,'rgb',r%bonds%rgb)
         call io%x(q,'order',r%bonds%order,0,4)
         call io%x(q,'imol',r%bonds%imol,0,2)
         call io%x(q,'bothends',r%bonds%bothends)
         call io%x(q,'hbond_classify',r%bonds%hbond_classify)
         call io%x(q,'hbond_rgb',r%bonds%hbond_rgb)
         call io%x(q,'hbond_dist',r%bonds%hbond_dist)
         call io%x(q,'hbond_ang',r%bonds%hbond_ang)

      case (reptype_labels)
         call io%obj(p,'labels',q,found)
         if (.not.found) return
         ! label types: 8 (Wyckoff) is for crystals only
         labelmax = 8
         if (sys(r%id)%c%ismolecule) labelmax = 7
         call io%x(q,'type',r%labels%type,0,labelmax)
         call io%x(q,'scale',r%labels%scale,lo=0d0)
         call io%x(q,'scale_cp',r%labels%scale_cp,lo=0d0)
         call io%x(q,'rgb',r%labels%rgb)
         call io%x(q,'const_size',r%labels%const_size)
         call io%x(q,'offset',r%labels%offset)
         call xref(q,'field',r%labels%fref)

      case (reptype_polyhedra)
         call io%obj(p,'polyhedra',q,found)
         if (.not.found) return
         call io%x(q,'alpha',r%poly%alpha)
         call io%x(q,'usecentercolor',r%poly%usecentercolor)
         call io%x(q,'rgb',r%poly%rgb)
         call io%x(q,'edge_rad',r%poly%edge_rad,lo=0d0)
         call io%x(q,'edge_rgb',r%poly%edge_rgb)
         call io%x(q,'usecentercolor_edge',r%poly%usecentercolor_edge)
         call io%x(q,'coplanar_eps',r%poly%coplanar_eps)
         call io%x(q,'showcorners',r%poly%showcorners)
         call io%x(q,'ihighlight',r%poly%ihighlight,0)
         call io%x(q,'rgbhl',r%poly%rgbhl)
         call io%x(q,'style_type',r%poly%style%type,0,atlisttype_NUM)

      case (reptype_unitcell)
         call io%obj(p,'unitcell',q,found)
         if (.not.found) return
         call io%x(q,'inner',r%uc%inner)
         call io%x(q,'coloraxes',r%uc%coloraxes)
         call io%x(q,'vaccutsticks',r%uc%vaccutsticks)
         call io%x(q,'radius',r%uc%radius,lo=0d0)
         call io%x(q,'radiusinner',r%uc%radiusinner,lo=0d0)
         call io%x(q,'rgb',r%uc%rgb)
         call io%x(q,'innersteplen',r%uc%innersteplen,lo=1d-2)
         call io%x(q,'innerstipple',r%uc%innerstipple)
         call io%x(q,'origin',r%uc%origin)

      case (reptype_axes)
         call io%obj(p,'axes',q,found)
         if (.not.found) return
         call io%x(q,'kind',r%axes%kind,0,1)
         call io%x(q,'placement',r%axes%placement,axplace_scene,axplace_window)
         call io%x(q,'coordtype',r%axes%coordtype,0,2)
         call io%x(q,'origin',r%axes%origin)
         call io%x(q,'winpos',r%axes%winpos)
         call io%x(q,'length',r%axes%length,lo=0d0)
         call io%x(q,'radius',r%axes%radius,lo=0d0)
         call io%x(q,'conelength',r%axes%conelength,lo=0d0)
         call io%x(q,'coneradius',r%axes%coneradius,lo=0d0)
         call io%x(q,'rgb',r%axes%rgb)
         call io%x(q,'showlabels',r%axes%showlabels)
         call io%x(q,'labelscale',r%axes%labelscale,lo=0d0)
         call io%x(q,'labelconstsize',r%axes%labelconstsize)
         call io%x(q,'labelrgb',r%axes%labelrgb)
         call io%x(q,'labelstr',r%axes%labelstr)
         call io%x(q,'labeldistance',r%axes%labeldistance)
         call io%x(q,'labeloffset',r%axes%labeloffset)
         call io%x(q,'scale',r%axes%scale,lo=0d0)
         call io%x(q,'scale_auto',r%axes%scale_auto)
         call io%x(q,'scalewithzoom',r%axes%scalewithzoom)

      case (reptype_symelem)
         call io%obj(p,'symelem',q,found)
         if (.not.found) return
         call io%x(q,'origin',r%symelem%origin)
         call io%x(q,'coordtype',r%symelem%coordtype,0,2)
         call io%x(q,'usecustomrgb',r%symelem%usecustomrgb)
         call io%x(q,'rgb',r%symelem%rgb)

      case (reptype_shapes)
         n = r%shapes%nshape
         call io%list(p,'shapes',a,n)
         if (.not.associated(a)) return
         if (io%reading) then
            if (allocated(r%shapes%shape)) deallocate(r%shapes%shape)
            allocate(r%shapes%shape(n))
            r%shapes%nshape = n
            r%shapes%isel = 0
         end if
         do i = 1, n
            call io%item(a,i,b)
            if (.not.associated(b)) cycle
            associate(sh => r%shapes%shape(i))
              call io%x(b,'kind',sh%kind,1,shapekind_NUM)
              call io%x(b,'shown',sh%shown)
              call io%x(b,'x1',sh%x1)
              call io%x(b,'v',sh%v)
              call io%x(b,'rad',sh%rad,lo=0d0)
              call io%x(b,'headr',sh%headr)
              call io%x(b,'headl',sh%headl)
              call io%x(b,'rgb',sh%rgb)
              call io%x(b,'alpha',sh%alpha)
              call io%x(b,'rim',sh%rim)
            end associate
         end do

      case (reptype_planar)
         n = r%planar%nshape
         call io%list(p,'shapes',a,n)
         if (.not.associated(a)) return
         if (io%reading) then
            if (allocated(r%planar%shape)) deallocate(r%planar%shape)
            allocate(r%planar%shape(n))
            r%planar%nshape = n
            r%planar%isel = 0
         end if
         do i = 1, n
            call io%item(a,i,b)
            if (.not.associated(b)) cycle
            associate(sh => r%planar%shape(i))
              call io%x(b,'kind',sh%kind,1,planarkind_NUM)
              call io%x(b,'shown',sh%shown)
              call io%x(b,'xc',sh%xc)
              call io%x(b,'hs',sh%hs)
              call io%x(b,'ang',sh%ang)
              call io%xa(b,'x',2,sh%x)
              if (io%reading) then
                 sh%npt = 0
                 if (allocated(sh%x)) sh%npt = size(sh%x,2)
              end if
              call io%x(b,'stroke',sh%stroke)
              call io%x(b,'width',sh%width,lo=0d0)
              call io%x(b,'dash',sh%dash,planardash_solid,planardash_dotted)
              call io%x(b,'rgb',sh%rgb)
              call io%x(b,'alpha',sh%alpha)
              call io%x(b,'fillrgb',sh%fillrgb)
              call io%x(b,'fillalpha',sh%fillalpha)
              call io%x(b,'filltype',sh%filltype,planarfill_none,planarfill_crosshatched)
              call io%x(b,'hatchang',sh%hatchang)
              call io%x(b,'hatchsp',sh%hatchsp,lo=1d-3)
              call io%x(b,'hatchw',sh%hatchw,lo=0d0)
              call io%x(b,'heads',sh%heads,planarheads_none,planarheads_both)
              call io%x(b,'headl',sh%headl)
              call io%x(b,'headw',sh%headw)
              call io%x(b,'infront',sh%infront)
            end associate
         end do

      case (reptype_text)
         n = r%text%ntext
         call io%list(p,'texts',a,n)
         if (.not.associated(a)) return
         if (io%reading) then
            if (allocated(r%text%t)) deallocate(r%text%t)
            allocate(r%text%t(n))
            r%text%ntext = n
            r%text%isel = 0
         end if
         do i = 1, n
            if (io%reading) r%text%t(i)%str = ""
            call io%item(a,i,b)
            if (.not.associated(b)) cycle
            associate(t => r%text%t(i))
              call io%x(b,'shown',t%shown)
              call io%xs(b,'str',t%str)
              call io%x(b,'placement',t%placement,textpos_screen,textpos_bond)
              call io%x(b,'winpos',t%winpos)
              call io%x(b,'pos',t%pos)
              call io%x(b,'idx1',t%idx1)
              call io%x(b,'idx2',t%idx2)
              call io%x(b,'scale',t%scale,lo=0d0)
              call io%x(b,'rgb',t%rgb)
              call io%x(b,'scalewithzoom',t%scalewithzoom)
              call io%x(b,'offset',t%offset)
              call io%x(b,'depth',t%depth)
              call io%x(b,'infront',t%infront)
            end associate
         end do

      case (reptype_measure)
         n = r%measure%nitem
         call io%list(p,'items',a,n)
         if (.not.associated(a)) return
         if (io%reading) then
            if (allocated(r%measure%item)) deallocate(r%measure%item)
            allocate(r%measure%item(n))
            r%measure%nitem = n
            r%measure%isel = 0
         end if
         do i = 1, n
            call io%item(a,i,b)
            if (.not.associated(b)) cycle
            associate(m => r%measure%item(i))
              call io%x(b,'shown',m%shown)
              call io%x(b,'n',m%n,2,4)
              call io%x(b,'idx',m%idx)
              call io%x(b,'stamp',m%stamp)
              call io%x(b,'rgb',m%rgb)
              call io%x(b,'rad',m%rad,lo=0d0)
              call io%x(b,'sectorrad',m%sectorrad,lo=0d0)
              call io%x(b,'sectoralpha',m%sectoralpha)
              call io%x(b,'planealpha',m%planealpha)
              call io%x(b,'planeext',m%planeext)
              call io%x(b,'textscale',m%textscale,lo=0d0)
              call io%x(b,'ndec',m%ndec,0,10)
              call io%x(b,'dashed',m%dashed)
              call io%x(b,'dashlen',m%dashlen,lo=1d-3)
              call io%x(b,'orient',m%orient)
              call io%x(b,'scalesystem',m%scalesystem)
              call io%x(b,'offset',m%offset)
              call xfer_cpanchors(b,m)
            end associate
         end do
         ! drop the items whose number of atoms is missing or wrong, or
         ! whose critical points were not read
         if (io%reading) then
            r%measure%item = pack(r%measure%item(1:n),r%measure%item(1:n)%n >= 2)
            r%measure%nitem = size(r%measure%item)
         end if

      case (reptype_isosurface)
         call io%obj(p,'isosurface',q,found)
         if (.not.found) return
         call xref(q,'field',r%iso%fref)
         call io%x(q,'imosel',r%iso%imosel,id_mo_id,0)
         call io%x(q,'imoidx',r%iso%imoidx)
         call io%x(q,'ilevel',r%iso%ilevel,0,iso_level_custom)
         call io%x(q,'icustom',r%iso%icustom,iso_custom_npts,iso_custom_ptsang)
         call io%x(q,'nptscustom',r%iso%nptscustom)
         call io%x(q,'ptsangcustom',r%iso%ptsangcustom,iso_ptsang_min,iso_ptsang_max)
         call io%x(q,'iregion',r%iso%iregion,0,iso_region_NUM)
         call io%x(q,'rgn_x',r%iso%rgn_x)
         call io%x(q,'nptsxyz',r%iso%nptsxyz)
         call io%x(q,'iregion_ap',r%iso%iregion_ap,0,iso_region_NUM)
         call io%x(q,'rgn_x_ap',r%iso%rgn_x_ap)
         call io%x(q,'ihighlight',r%iso%ihighlight,0)
         call io%x(q,'rgbhl',r%iso%rgbhl)
         call io%x(q,'hist_xscale',r%iso%hist_xscale,0,hscale_num)
         call io%x(q,'hist_yscale',r%iso%hist_yscale,1,hscale_num)
         if (io%reading) then
            ! the applied grid: all zero (native), or at least the
            ! smallest custom grid along each axis and within the
            ! total-points cap
            r%iso%nptscustom = min(max(r%iso%nptscustom,iso_npts_custom_min),iso_npts_custom_max)
            if (any(r%iso%nptsxyz /= 0)) then
               if (any(r%iso%nptsxyz < iso_npts_custom_min) .or.&
                  product(real(r%iso%nptsxyz,8)) > iso_maxpts_total) r%iso%nptsxyz = 0
            end if
         end if

         n = r%iso%niso
         call io%list(q,'surfaces',a,n)
         if (.not.associated(a)) return
         if (io%reading) then
            if (allocated(r%iso%slot)) deallocate(r%iso%slot)
            allocate(r%iso%slot(n))
            r%iso%niso = n
            r%iso%isel = min(1,n)
         end if
         do i = 1, n
            if (io%reading) r%iso%slot(i)%uid = new_uid()
            call io%item(a,i,b)
            if (.not.associated(b)) cycle
            associate(sl => r%iso%slot(i))
              call io%x(b,'isoval',sl%isoval)
              call io%x(b,'shown',sl%shown)
              call io%x(b,'rgb',sl%rgb)
              call io%x(b,'alpha',sl%alpha)
              call io%x(b,'imap_mode',sl%imap_mode,iso_map_color,iso_map_expr)
              call xref(b,'mapfield',sl%mapref)
              call io%x(b,'mapexpr',sl%mapexpr)
              call io%x(b,'icmap',sl%icmap,1,iw_ncmap)
              call io%x(b,'icmap_auto',sl%icmap_auto)
              call io%x(b,'maprange',sl%maprange)
              call io%x(b,'maprange_auto',sl%maprange_auto)
            end associate
         end do

      case (reptype_cps)
         call io%obj(p,'cps',q,found)
         if (.not.found) return
         call xref(q,'field',r%cps%fref)
         call io%x(q,'show',r%cps%show)
         call io%x(q,'rgb',r%cps%rgb)
         call io%x(q,'rad',r%cps%rad)
         call io%x(q,'radscale',r%cps%radscale,lo=0d0)

      case (reptype_gpaths)
         call io%obj(p,'gpaths',q,found)
         if (.not.found) return
         call xref(q,'field',r%gpaths%fref)
         call io%x(q,'rgb',r%gpaths%rgb)
         call io%x(q,'rad',r%gpaths%rad,lo=0d0)
         call io%x(q,'style',r%gpaths%style,gpaths_style_tube,gpaths_style_beads)
         call io%x(q,'showends',r%gpaths%showends)
         call io%x(q,'fthr',r%gpaths%fthr)

      case (reptype_legend)
         call io%obj(p,'legend',q,found)
         if (.not.found) return
         call io%x(q,'corner',r%legend%corner,legcorner_topleft,legcorner_bottomright)
         call io%x(q,'scale',r%legend%scale,lo=0d0)
         call io%x(q,'textrgb',r%legend%textrgb)
         call io%x(q,'bgrgb',r%legend%bgrgb)
         call io%x(q,'bgalpha',r%legend%bgalpha)
         call io%x(q,'border',r%legend%border)
         call io%x(q,'borderrgb',r%legend%borderrgb)

      case (reptype_scalebar)
         call io%obj(p,'scalebar',q,found)
         if (.not.found) return
         call io%x(q,'value',r%scalebar%value,lo=0d0)
         call io%x(q,'unit',r%scalebar%unit,scaleunit_angstrom,scaleunit_pm)
         call io%x(q,'heads',r%scalebar%heads,planarheads_none,planarheads_both)
         call io%x(q,'corner',r%scalebar%corner,legcorner_topleft,scalepos_custom)
         call io%x(q,'winpos',r%scalebar%winpos)
         call io%x(q,'showtext',r%scalebar%showtext)
         call io%x(q,'customtext',r%scalebar%customtext)
         call io%x(q,'text',r%scalebar%text)
         call io%x(q,'scale',r%scalebar%scale,lo=0d0)
         call io%x(q,'width',r%scalebar%width,lo=0d0)
         call io%x(q,'rgb',r%scalebar%rgb)
         call io%x(q,'textrgb',r%scalebar%textrgb)

      case (reptype_colorbar)
         call io%obj(p,'colorbar',q,found)
         if (.not.found) return
         call io%x(q,'corner',r%colorbar%corner,legcorner_topleft,scalepos_custom)
         call io%x(q,'winpos',r%colorbar%winpos)
         call io%x(q,'vertical',r%colorbar%vertical)
         call io%x(q,'length',r%colorbar%length,lo=0d0)
         call io%x(q,'thick',r%colorbar%thick,lo=0d0)
         call io%x(q,'nticks',r%colorbar%nticks,2,100)
         call io%x(q,'scale',r%colorbar%scale,lo=0d0)
         call io%x(q,'textrgb',r%colorbar%textrgb)
         call io%x(q,'border',r%colorbar%border)
         call io%x(q,'borderrgb',r%colorbar%borderrgb)
         call io%x(q,'showtitle',r%colorbar%showtitle)
         call io%x(q,'customtitle',r%colorbar%customtitle)
         call io%x(q,'title',r%colorbar%title)
      end select

    end subroutine xfer_options

    !> The style rows of the object.
    subroutine xfer_rows()

      call io%obj(p,'rows',q,found)
      if (.not.found) return

      ! atom colors and radii, molecule tints and scales
      if (reptype_is_atombased(r%type) .or. r%type == reptype_gpaths) then
         if (allocated(r%atoms%style%rgb)) call io%x(q,'atom_rgb',r%atoms%style%rgb)
         if (allocated(r%atoms%style%rad)) call io%x(q,'atom_rad',r%atoms%style%rad)
      end if
      if (reptype_is_atombased(r%type)) then
         if (allocated(r%mols%style%tint_rgb)) call io%x(q,'mol_tint_rgb',r%mols%style%tint_rgb)
         if (allocated(r%mols%style%scale_rad)) call io%x(q,'mol_scale_rad',r%mols%style%scale_rad)
      end if

      select case (r%type)
      case (reptype_bonds)
         ! bonds of this object's own (not the system's) neighbor stars
         ! are made again with its bonding parameters
         call io%x(q,'use_sys_nstar',r%bonds%style%use_sys_nstar)
         if (allocated(r%bonds%style%shown)) call io%x(q,'bond_shown',r%bonds%style%shown)
         if (io%reading .and. r%bonds%style%isinit) call r%bonds%style%refresh_nstar(r)

      case (reptype_labels)
         if (allocated(r%labels%style%shown)) call io%x(q,'label_shown',r%labels%style%shown)
         if (allocated(r%labels%style%str)) call io%x(q,'label_str',r%labels%style%str)
         if (allocated(r%labels%style%cpshown)) call io%x(q,'label_cpshown',r%labels%style%cpshown)
         if (allocated(r%labels%style%cpstr)) call io%x(q,'label_cpstr',r%labels%style%cpstr)

      case (reptype_polyhedra)
         if (allocated(r%poly%style%shown)) call io%x(q,'poly_shown',r%poly%style%shown)
         if (allocated(r%poly%style%corner)) call io%x(q,'poly_corner',r%poly%style%corner)
         if (allocated(r%poly%style%dmin)) call io%x(q,'poly_dmin',r%poly%style%dmin)
         if (allocated(r%poly%style%dmax)) call io%x(q,'poly_dmax',r%poly%style%dmax)

      case (reptype_symelem)
         if (allocated(r%symelem%style%shown)) call io%x(q,'symelem_shown',r%symelem%style%shown)

      case (reptype_legend)
         if (allocated(r%legend%shown)) call io%x(q,'legend_shown',r%legend%shown)
         if (allocated(r%legend%label)) call io%x(q,'legend_label',r%legend%label)
      end select

    end subroutine xfer_rows

    !> A field reference, as the slot of the field (the scene files make
    !> the fields again in the same slots). On reading, a slot without
    !> a field leaves the reference as it was.
    subroutine xref(qq,key,fr)
      type(json_value), pointer, intent(in) :: qq
      character*(*), intent(in) :: key
      type(field_ref), intent(inout) :: fr

      integer :: id

      id = fr%id
      call io%x(qq,key,id)
      if (io%reading .and. id >= 0) then
         if (sys(r%id)%goodfield(id)) call fr%set(sys(r%id),id)
      end if

    end subroutine xref

    !> The critical points measurement m is anchored to (idx(1,j) < 0,
    !> entries of the system's anchor registry), by content: the field
    !> slot, the CP type, and its fractional position. On reading, the
    !> anchors are made again in the registry; an item with a CP anchor
    !> that cannot be read is dropped (n = 0).
    subroutine xfer_cpanchors(b,m)
      type(json_value), pointer, intent(in) :: b
      type(measurement_item), intent(inout) :: m

      type(json_value), pointer :: c, d
      integer :: j, k, nc, ifield, icel, typ, loff(3)
      real*8 :: x(3)
      logical :: ok

      if (m%n < 2 .or. m%n > 4) return
      if (all(m%idx(1,1:m%n) >= 0)) return
      nc = m%n
      call io%list(b,'cps',c,nc)
      ok = associated(c) .and. (nc == m%n)
      do j = 1, m%n
         if (.not.ok) exit
         ! one entry per anchor (empty for an atom)
         call io%item(c,j,d)
         if (m%idx(1,j) >= 0) cycle
         if (.not.associated(d)) then
            ok = .false.
            exit
         end if
         if (io%reading) then
            ifield = -1
            typ = huge(1)
            x = huge(1d0)
         else
            k = -m%idx(1,j)
            if (k > sysc(r%id)%ncpanc) cycle
            ! resolving the anchor finds its field again if it moved to
            ! another slot
            ! (an anchor whose field is gone is not saved: the item is
            ! dropped on reading, as it is not drawn now)
            if (.not.cp_anchor_resolve(r%id,m%idx(1,j),ifield,icel,loff)) cycle
            ifield = sysc(r%id)%cpanc(k)%ifield
            typ = sysc(r%id)%cpanc(k)%typ
            x = sysc(r%id)%cpanc(k)%x
         end if
         call io%x(d,'field',ifield,0,huge(1))
         call io%x(d,'typ',typ,-3,3)
         call io%x(d,'x',x)
         if (io%reading) then
            k = 0
            if (typ /= huge(1) .and. all(abs(x) <= 1d1)) k = cp_anchor_add(r%id,ifield,typ,x)
            if (k == 0) then
               ok = .false.
               exit
            end if
            ! the entry may be another image of the same CP
            m%idx(1,j) = -k
            m%idx(2:4,j) = m%idx(2:4,j) - nint(sysc(r%id)%cpanc(k)%x - x)
         end if
      end do
      if (io%reading .and. .not.ok) m%n = 0

    end subroutine xfer_cpanchors

  end subroutine representation_xfer_json

end submodule json
