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

! The plane map window (the PLANE keyword), the form of the maps it
! shares with the editor of the plane map objects, and that editor.
submodule (windows) plane
  use interfaces_cimgui
  use param, only: bohrtoa
  use tools_math, only: niso_lin, niso_log, niso_atan, niso_bader, niso_manual
  use representations, only: plane_npts_min, plane_npts_max, plane_nlev_max, plane_auto_spacing,&
     plane_auto_ntot, plane_auto_tmax
  implicit none

  ! transient markers in the view: the three points (origin, end of
  ! the first side, point on the second side), the region, and the
  ! point under the mouse in the map; opacity of the markers and of
  ! the region, and radius of the points (bohr), a fraction of the
  ! scene radius between the two limits
  real(c_float), parameter :: pt_rgb(3,3) = reshape((/&
     0.95_c_float,0.80_c_float,0.10_c_float,&
     0.95_c_float,0.35_c_float,0.20_c_float,&
     0.80_c_float,0.45_c_float,0.95_c_float/),(/3,3/))
  real(c_float), parameter :: region_rgb(3) = (/0.60_c_float,0.60_c_float,0.60_c_float/)
  real(c_float), parameter :: cur_rgb(3) = (/0.10_c_float,0.85_c_float,0.30_c_float/)
  real(c_float), parameter :: marker_alpha = 0.8_c_float
  real(c_float), parameter :: region_alpha = 0.12_c_float ! a flat box draws its face twice
  real*8, parameter :: marker_radmin = 0.20d0
  real*8, parameter :: marker_radmax = 0.50d0

  ! the names of the three points
  character(len=7), parameter :: pt_name(3) = (/"Origin ","X end  ","Y point"/)

  ! the kinds of contour levels, in the order of the combo
  integer, parameter :: ctype_of_combo(0:4) = (/niso_lin,niso_log,niso_atan,niso_bader,niso_manual/)

  ! most lattice translations tried for each candidate label
  integer, parameter :: maxtrans = 20000

contains

  !> Draw the plane map window: the values of a field of the system
  !> shown by the anchor view, one of its derivatives, or an
  !> arithmetic expression, on a region of a plane (defined by three
  !> points, typed or picked in the view), as PLANE calculates them,
  !> shown as a colored map with contour lines.
  module subroutine draw_plane(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8,&
       iw_calcwidth, iw_close_event, iw_setpos_bottomright, iw_calcheight
    use representations, only: rep_shape, shapekind_sphere, shapekind_box, bondrad_def
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged, ismol, ldum, xediting, refit, hover, fresh
    integer :: isys, iview, istat, i, ipick
    real*8 :: xc(3), rad, v(3,3), u, xhov(3)
    type(ImVec2) :: sz0
    real(c_float) :: htable
    type(rep_shape), allocatable :: shp(:)

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)
    ismol = .false.
    if (goodsys) ismol = sys(isys)%c%ismolecule
    u = 1d0
    if (ismol) u = bohrtoa

    ! initialize the state on the first pass
    if (w%firstpass) then
       call w%plane_free_texture()
       w%pl = plane_state()
       w%errmsg = ""
    end if

    ! a new system: the plane goes to a face of the cell (crystals),
    ! to the plane of the first three atoms (molecules), or in between
    ! for slabs and chains (see plane_def_set_defaults), and the
    ! quantity is the value of the reference field
    if (goodsys) then
       if (.not.w%pl%sys%is(isys) .or. syschanged) then
          call w%pl%sys%set(isys)
          call cancel_pick()
          w%errmsg = ""
          call w%pl%d%set_defaults(isys)
          w%pl%viewgen = -1
       end if
       ! the reference field, if the field is gone
       if (.not.w%pl%d%fref%ok(sys(isys))) call w%pl%d%fref%set(sys(isys),sys(isys)%iref)
    end if

    ! a pending pick in the view it was armed on; dropped, releasing the
    ! view, if the system is not ready or the window moved to another view
    if (w%pl%picking) then
       if (goodsys .and. iview == w%pl%pickview) then
          call view_pick_result(iview,w%id,isys,w%pl%pick,istat,xc)
          if (istat == ipick_point) then
             w%pl%d%xp(:,w%pl%pickend) = form_coords(isys,xc)
             w%pl%viewdirty = .true.
          end if
          if (istat /= ipick_pending) w%pl%picking = .false.
       else
          call cancel_pick()
       end if
    end if

    hover = .false.
    if (goodsys) then
       ! the system
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! the plane, the region, and the number of points
       xediting = .false.
       if (plane_def_form(w%pl%d,isys,ttshown,xediting,.true.,ipick,w%pl%picking,w%pl%showcval,w%pl%vscale,&
          w%pl%sec)) w%pl%viewdirty = .true.
       if (ipick > 0) then
          w%pl%picking = .true.
          w%pl%pickend = ipick
          w%pl%pickview = iview
          call w%pl%pick%arm()
          call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick " // trim(pt_name(ipick)) // ": an atom,&
             & a critical point, or a point in space",w%id,acceptempty=.true.)
       end if

       ! the labels, and the map in the view
       call draw_labels_option()
       if (sys(isys)%c%ncel == 0) &
          call iw_text("This system has no atoms: there are no field values to map",danger=.true.,wrap=.true.)

       ! the values on the plane and what is made from them, calculated
       ! again if anything they depend on changed
       call update_grid(refit,fresh)
       if (allocated(w%pl%ff)) then
          call update_contours()
          call update_texture()
          call update_marks()
       end if
       call sync_view(fresh)

       ! the map takes the rest of the window, leaving room for the
       ! note about the labels (if shown), the table of values (if
       ! open, it takes part of the space), and Close
       call igGetContentRegionAvail(sz0)
       sz0%y = sz0%y - iw_calcheight(merge(3,2,w%pl%showlabels .and. w%pl%marksmany),0,.true.)
       htable = 0._c_float
       if (w%pl%secdata%open .and. allocated(w%pl%ff)) then
          htable = max(0.35_c_float * sz0%y,iw_calcheight(4,0,.false.))
          sz0%y = sz0%y - htable - iw_calcheight(0,0,.true.)
       end if
       sz0%y = max(sz0%y,iw_calcheight(8,0,.false.))
       call draw_map(sz0,refit,hover,xhov)

       ! the table of values
       call draw_data_section()

       ! the three points, the region, and the point under the mouse in
       ! the anchor view, for as long as the window is open
       if (associated(win(iview)%sc)) then
          if (win(iview)%sc%isinit /= 0) then
             rad = min(max(0.01d0 * real(win(iview)%sc%scenerad,8),marker_radmin),marker_radmax)
             allocate(shp(3))
             do i = 1, 3
                shp(i) = rep_shape(kind=shapekind_sphere,x1=cart_to_abs(form_to_cart(isys,w%pl%d%xp(:,i))),&
                   rad=rad,rgb=pt_rgb(:,i),alpha=marker_alpha)
             end do
             ! the region (unless the map itself is in the view)
             if (allocated(w%pl%ff) .and. .not.(w%pl%showview .and. w%pl%viewgen >= 0)) then
                v = 0d0
                v(:,1) = real(size(w%pl%ff,1)-1,8) * w%pl%uu
                v(:,2) = real(size(w%pl%ff,2)-1,8) * w%pl%vv
                shp = [shp,rep_shape(kind=shapekind_box,x1=cart_to_abs(w%pl%x0),v=v,rad=0.15d0*rad,&
                   rgb=region_rgb,alpha=region_alpha)]
             end if
             if (hover) &
                shp = [shp,rep_shape(kind=shapekind_sphere,x1=cart_to_abs(xhov),&
                   rad=max(rad,1.6d0*bondrad_def),rgb=cur_rgb,alpha=1._c_float)]
             call win(iview)%sc%show_transient_shapes(w%uid,1,shp)
          end if
       end if

       if (allocated(w%pl%errmsg)) then
          if (len_trim(w%pl%errmsg) > 0) call iw_text(w%pl%errmsg,danger=.true.,wrap=.true.)
       end if
       if (allocated(w%errmsg)) then
          if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)
       end if
    end if

    ! bottom row: close
    call iw_setpos_bottomright(len("Close"),1)
    if (iw_button("Close")) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) then
       call cancel_pick()
       call w%end()
    end if

  contains

    !> The Display section: the labels of the atoms and critical
    !> points near the plane, and whether the map is shown in the
    !> view, and its opacity there. More: what is labeled, the color
    !> and size of the labels, their text, and the distance to the
    !> plane.
    subroutine draw_labels_option()
      use utils, only: iw_checkbox, iw_dragfloat_realc, iw_coloredit, iw_combo_simple, iw_separatortext

      call iw_separatortext("Display",w%pl%secdisp,hasmore=w%pl%showlabels)
      if (.not.w%pl%secdisp%open) return
      ldum = iw_checkbox("Labels##planeshowlabels",w%pl%showlabels,labelleft=.true.)
      call iw_tooltip("Label the atoms and the critical points of the reference field near the plane",ttshown)
      if (iw_checkbox("Show in view##planeshow",w%pl%showview,gap=.true.,labelleft=.true.)) w%pl%viewgen = -1
      call iw_tooltip("Draw this map in the view, following the changes made here, for as long as&
         & the window is open (Draw > Plane Map makes a map that stays)",ttshown)
      if (w%pl%showview) then
         call iw_text("Opacity",alignframe=.true.,gap=.true.)
         if (iw_dragfloat_realc("##planealpha",x1=w%pl%alpha,speed=0.01_c_float,min=0._c_float,&
            max=1._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)) &
            w%pl%viewdirty = .true.
         call iw_tooltip("Opacity of the map in the view (1 = opaque; the contour lines stay opaque)",ttshown)
      end if
      if (.not.w%pl%secdisp%more .or. .not.w%pl%showlabels) return

      ! more: what is labeled, and the color and size of the labels
      ldum = iw_checkbox("Critical points##planeshowcps",w%pl%showcps,labelleft=.true.)
      call iw_tooltip("Label also the critical points of the reference field other than the nuclei&
         & (n = non-nuclear maximum, b = bond, r = ring, c = cage)",ttshown)
      call iw_text("Color",alignframe=.true.,gap=.true.)
      ldum = iw_coloredit("##planergblabel",rgb=w%pl%rgblabel,sameline=.true.)
      call iw_tooltip("Color of the labels and their dots",ttshown)
      call iw_text("Size",alignframe=.true.,gap=.true.)
      ldum = iw_dragfloat_realc("##planelscale",x1=w%pl%lscale,speed=0.01_c_float,min=0.1_c_float,&
         max=10._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
      call iw_tooltip("Size of the labels and their dots, relative to the text of the interface",ttshown)

      ! more: the text, and the distance to the plane
      call iw_text("Text",alignframe=.true.)
      call iw_combo_simple("##planeltext","Species name"//c_null_char//"Atomic symbol"//c_null_char//&
         "Atom name"//c_null_char,w%pl%ltext,sameline=.true.)
      call iw_tooltip("Text in the labels of the atoms",ttshown)
      call iw_text("Distance",alignframe=.true.,gap=.true.)
      ldum = iw_dragfloat_real8("##planehmax",x1=w%pl%hmax,speed=0.01d0,min=0d0,decimal=3,sameline=.true.)
      call iw_tooltip("Largest distance to the plane of the labeled atoms and critical points (" //&
         lunit() // ")",ttshown)
      ldum = iw_checkbox("Show distance##planeshowdist",w%pl%showdist,gap=.true.,labelleft=.true.)
      call iw_tooltip("Add to each label the signed distance of the atom or critical point to the&
         & plane (" // lunit() // ", positive on the side of the normal; omitted if zero to two&
         & decimals)",ttshown)

    end subroutine draw_labels_option

    !> Show the map of the window in the anchor view, as a transient
    !> plane map (it goes away when the window stops showing it), armed
    !> in every frame. fresh = the values of the window are those of the
    !> form (else they are being edited, or wait for the threads, and
    !> the view keeps the last map). The view is given the values of the
    !> window (it never calculates them) when they are new, and the form
    !> when it changed; with no values (an error), it shows nothing.
    subroutine sync_view(fresh)
      use representations, only: slice_cache
      logical, intent(in) :: fresh

      integer :: itrep
      logical :: found

      if (.not.w%pl%showview .or. .not.associated(win(iview)%sc)) return
      if (win(iview)%sc%isinit == 0) return

      ! nothing to show: up to date but no values, or nothing given yet
      if (fresh .and. .not.allocated(w%pl%ff)) return
      if (.not.fresh .and. w%pl%viewgen < 0) return
      call win(iview)%sc%show_transient_slice(w%uid,1,itrep,found)
      if (itrep <= 0) return

      associate(r => win(iview)%sc%reptrans(itrep))
        r%slice%nocalc = .true.
        if (.not.found) then
           ! a new slot (or the old one was reaped): hidden until it has
           ! the values of the window
           w%pl%viewgen = -1
           r%shown = .false.
        end if
        if (.not.fresh) return
        if (w%pl%viewgen == w%pl%gen .and. .not.w%pl%viewdirty) return

        r%name = "Plane Map"
        r%shown = .true.
        r%slice%d = w%pl%d
        r%slice%alpha = w%pl%alpha
        ! the values of the window, if they are new to the view
        if (w%pl%viewgen /= w%pl%gen) &
           r%slice%c = slice_cache(ff=w%pl%ff,x0=w%pl%x0,uu=w%pl%uu,vv=w%pl%vv,&
              e1=w%pl%e1,e2=w%pl%e2,errmsg="",keyexpr=trim(w%pl%d%expr))
      end associate
      w%pl%viewgen = w%pl%gen
      w%pl%viewdirty = .false.
      win(iview)%sc%forcebuildlists = .true.

    end subroutine sync_view

    !> The map (size sz): the colors, the contour lines, and the labels,
    !> with the color bar on its right. If refit, fit the axes to the
    !> region. Returns in hover whether the mouse is on the region and
    !> the point under it (xhov, cell-frame Cartesian bohr).
    subroutine draw_map(sz,refit,hover,xhov)
      use utils, only: iw_colormap_id
      use gui_main, only: g
      type(ImVec2), intent(in) :: sz
      logical, intent(in) :: refit
      logical, intent(out) :: hover
      real*8, intent(out) :: xhov(3)

      character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
      integer :: j, nx, ny
      integer(c_int) :: flags, coltext
      real(c_double) :: xmin, xmax, ymin, ymax, crn(2,4), mu, mv
      real(c_float) :: px(4), py(4), wbar, fs
      real(c_float), allocatable :: box(:,:), ang(:)
      real*8, allocatable :: gap(:), sval(:)
      integer, allocatable :: ival(:), icand(:)
      integer :: nbox, nval, k
      integer(c_int) :: col(2)
      type(ImVec2) :: pmin, pmax, szp, pt, tsz
      type(c_ptr) :: dl

      hover = .false.
      xhov = 0d0
      if (.not.allocated(w%pl%ff)) return
      nx = size(w%pl%ff,1)
      ny = size(w%pl%ff,2)

      ! the corners of the region in plot coordinates
      crn(:,1) = to_plot(w%pl%x0)
      crn(:,2) = to_plot(w%pl%x0 + real(nx-1,8) * w%pl%uu)
      crn(:,3) = to_plot(w%pl%x0 + real(nx-1,8) * w%pl%uu + real(ny-1,8) * w%pl%vv)
      crn(:,4) = to_plot(w%pl%x0 + real(ny-1,8) * w%pl%vv)

      ! the plot leaves room for the color bar on its right
      wbar = iw_calcwidth(10,0)
      szp = sz
      szp%x = max(sz%x - wbar,iw_calcwidth(20,0))

      str1 = "##planemap" // c_null_char
      flags = ior(ImPlotFlags_NoTitle,ImPlotFlags_Equal)
      flags = ior(flags,ImPlotFlags_NoLegend)
      flags = ior(flags,ImPlotFlags_Crosshairs)
      if (ipBeginPlot(c_loc(str1),szp,flags)) then
         str2 = "u (" // lunit() // ")" // c_null_char
         str3 = "v (" // lunit() // ")" // c_null_char
         call ipSetupAxes(c_loc(str2),c_loc(str3),ImPlotAxisFlags_None,ImPlotAxisFlags_None)
         if (refit) then
            call ipSetupAxisLimits(ImAxis_X1,minval(crn(1,:)),maxval(crn(1,:)),ImPlotCond_Always)
            call ipSetupAxisLimits(ImAxis_Y1,minval(crn(2,:)),maxval(crn(2,:)),ImPlotCond_Always)
         end if

         ! the colors: the texture on the region, inside the plot
         call ipGetPlotCurrentLimits(xmin,xmax,ymin,ymax)
         call ipPlotToPixels(xmin,ymin,ImAxis_X1,ImAxis_Y1,px(1),py(1))
         call ipPlotToPixels(xmax,ymax,ImAxis_X1,ImAxis_Y1,px(2),py(2))
         pmin%x = min(px(1),px(2))
         pmin%y = min(py(1),py(2))
         pmax%x = max(px(1),px(2))
         pmax%y = max(py(1),py(2))
         dl = igGetWindowDrawList()
         call ImDrawList_PushClipRect(dl,pmin,pmax,.true._c_bool)
         if (w%pl%tex /= 0) then
            do j = 1, 4
               call ipPlotToPixels(crn(1,j),crn(2,j),ImAxis_X1,ImAxis_Y1,px(j),py(j))
            end do
            ! the texels are centered on the points of the grid
            call ImDrawList_AddImageQuad(dl,int(w%pl%tex,c_intptr_t),&
               ImVec2(px(1),py(1)),ImVec2(px(2),py(2)),ImVec2(px(3),py(3)),ImVec2(px(4),py(4)),&
               ImVec2(0.5_c_float/nx,0.5_c_float/ny),ImVec2(1._c_float-0.5_c_float/nx,0.5_c_float/ny),&
               ImVec2(1._c_float-0.5_c_float/nx,1._c_float-0.5_c_float/ny),&
               ImVec2(0.5_c_float/nx,1._c_float-0.5_c_float/ny),-1_c_int)
         end if
         call ImDrawList_PopClipRect(dl)

         ! the labels go where their dot and text are (ImGui scales text
         ! linearly in the font size, so the measurement at the
         ! interface size scales with it); then the values of the
         ! contour lines, where they do not cover a label or another
         ! value, along the line, which has a gap under each
         dl = igGetWindowDrawList()
         coltext = igGetColorU32_Vec4(ImVec4(w%pl%rgblabel(1),w%pl%rgblabel(2),w%pl%rgblabel(3),1._c_float))
         fs = w%pl%lscale * g%FontSize
         nbox = 0
         nval = 0
         allocate(box(4,w%pl%nmark+w%pl%ncval),ival(w%pl%ncval),icand(w%pl%ncval),sval(w%pl%ncval),&
            gap(w%pl%ncval),ang(w%pl%ncval))
         if (w%pl%showlabels) then
            do j = 1, w%pl%nmark
               call ipPlotToPixels(w%pl%xmark(1,j),w%pl%xmark(2,j),ImAxis_X1,ImAxis_Y1,px(1),py(1))
               str1 = trim(w%pl%lmark(j)) // c_null_char
               call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
               nbox = nbox + 1
               box(:,nbox) = (/px(1) - 3._c_float * w%pl%lscale,py(1) - (tsz%y + 1._c_float) * w%pl%lscale,&
                  px(1) + (4._c_float + tsz%x) * w%pl%lscale,py(1) + 3._c_float * w%pl%lscale/)
            end do
         end if
         if (w%pl%d%showcont) then
            if (w%pl%showcval .and. w%pl%ncval > 0) &
               call place_values(pmin,pmax,box,nbox,nval,ival,icand,sval,gap,ang)
            if (w%pl%ncpos > 0) call plot_contours(1,w%pl%d%rgbpos,w%pl%cpos,w%pl%ncpos,nval,ival,sval,gap)
            if (w%pl%ncneg > 0) call plot_contours(2,w%pl%d%rgbneg,w%pl%cneg,w%pl%ncneg,nval,ival,sval,gap)
         end if

         call ImDrawList_PushClipRect(dl,pmin,pmax,.true._c_bool)
         col(1) = igGetColorU32_Vec4(ImVec4(w%pl%d%rgbpos(1),w%pl%d%rgbpos(2),w%pl%d%rgbpos(3),1._c_float))
         col(2) = igGetColorU32_Vec4(ImVec4(w%pl%d%rgbneg(1),w%pl%d%rgbneg(2),w%pl%d%rgbneg(3),1._c_float))
         do j = 1, nval
            k = ival(j)
            call ipPlotToPixels(w%pl%xcval(1,icand(j),k),w%pl%xcval(2,icand(j),k),ImAxis_X1,ImAxis_Y1,px(1),py(1))
            str1 = fmt_level(w%pl%zcval(k)) // c_null_char
            call own_add_text_rotated(dl,g%Font,w%pl%vscale*g%FontSize,px(1),py(1),ang(j),col(w%pl%icval(1,k)),c_loc(str1))
         end do
         if (w%pl%showlabels) then
            do j = 1, w%pl%nmark
               call ipPlotToPixels(w%pl%xmark(1,j),w%pl%xmark(2,j),ImAxis_X1,ImAxis_Y1,px(1),py(1))
               call ImDrawList_AddCircleFilled(dl,ImVec2(px(1),py(1)),3._c_float*w%pl%lscale,coltext,8_c_int)
               str1 = trim(w%pl%lmark(j)) // c_null_char
               call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
               pt%x = px(1) + 4._c_float * w%pl%lscale
               pt%y = py(1) - (tsz%y + 1._c_float) * w%pl%lscale
               call ImDrawList_AddText_FontPtr(dl,g%Font,fs,pt,coltext,c_loc(str1),c_null_ptr,&
                  0._c_float,c_null_ptr)
            end do
         end if
         call ImDrawList_PopClipRect(dl)

         ! the point under the mouse
         if (ipIsPlotHovered()) then
            call ipGetPlotMousePos(mu,mv)
            call plot_point(mu,mv,hover,xhov)
         end if
         call ipEndPlot()
      end if

      ! the color bar
      call igSameLine(0._c_float,-1._c_float)
      str1 = "##planecbar" // c_null_char
      str2 = "%.3g" // c_null_char
      call ipColormapScale(c_loc(str1),w%pl%tmin,w%pl%tmax,ImVec2(wbar,szp%y),c_loc(str2),&
         iw_colormap_id(w%pl%d%icmap))
      if (w%pl%d%cscale == 1) then
         call iw_tooltip("Colors of the map: log10 of the absolute value",ttshown)
      elseif (w%pl%d%cscale == 2) then
         call iw_tooltip("Colors of the map: 2/pi atan of the value",ttshown)
      else
         call iw_tooltip("Colors of the map: the value",ttshown)
      end if

      if (w%pl%showlabels .and. w%pl%marksmany) &
         call iw_text("(The region is too large for labels)",disabled=.true.)

    end subroutine draw_map

    !> Where the values of the contour lines go, longest line first:
    !> at one of the places along each line long enough for them,
    !> inside the plot (pmin-pmax), where they do not cover a label or
    !> come close to a value already placed (the nbox boxes in box, in
    !> pixels: xmin, ymin, xmax, ymax). The place tried first moves
    !> along the line from one level to the next (by the golden ratio),
    !> so the values of nested lines do not pile up. Returns the nval
    !> values placed: their line in the list of values (ival), the
    !> place (icand) and its distance along the line (sval, plot
    !> units), half the gap they leave in the line (gap, plot units),
    !> and the angle of their text (ang, radians, clockwise on the
    !> screen).
    subroutine place_values(pmin,pmax,box,nbox,nval,ival,icand,sval,gap,ang)
      use param, only: pi
      type(ImVec2), intent(in) :: pmin, pmax
      real(c_float), intent(inout) :: box(:,:)
      integer, intent(inout) :: nbox
      integer, intent(out) :: nval
      integer, intent(out) :: ival(:), icand(:)
      real*8, intent(out) :: sval(:), gap(:)
      real(c_float), intent(out) :: ang(:)

      real*8, parameter :: golden = 0.6180339887498949d0

      character(kind=c_char,len=:), allocatable, target :: str1
      real(c_float) :: px, py, px1, py1, scl, bx(4), pad, sep, a, hx, hy
      type(ImVec2) :: tsz
      integer :: i, k, l, m, k0
      real*8 :: g, sc
      logical :: free

      ! pixels per plot unit (the axes have the same scale)
      call ipPlotToPixels(0d0,0d0,ImAxis_X1,ImAxis_Y1,px,py)
      call ipPlotToPixels(1d0,0d0,ImAxis_X1,ImAxis_Y1,px1,py1)
      scl = abs(px1 - px)
      pad = 3._c_float * w%pl%vscale
      nval = 0
      if (scl <= 0._c_float) return

      do i = 1, w%pl%ncval
         str1 = fmt_level(w%pl%zcval(i)) // c_null_char
         call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
         tsz%x = tsz%x * w%pl%vscale
         tsz%y = tsz%y * w%pl%vscale
         if (w%pl%lcval(i) * scl < 3._c_float * tsz%x) cycle
         g = (0.5d0 * tsz%x + pad) / scl
         sep = tsz%y

         ! the places, from the one nearest to the phase of this level,
         ! alternating to both sides
         k0 = floor(modulo(w%pl%icval(4,i) * golden,1d0) * plane_ncand)
         do m = 0, plane_ncand-1
            l = (m + 1) / 2
            if (mod(m,2) == 1) l = -l
            k = modulo(k0 + l,plane_ncand) + 1
            sc = (k - 0.5d0) / plane_ncand * w%pl%lcval(i)
            if (sc - g < 0d0 .or. sc + g > w%pl%lcval(i)) cycle
            call ipPlotToPixels(w%pl%xcval(1,k,i),w%pl%xcval(2,k,i),ImAxis_X1,ImAxis_Y1,px,py)

            ! the direction of the line on the screen (y goes down),
            ! turned so the text reads left to right (bottom to top if
            ! vertical)
            a = real(atan2(-w%pl%tcval(2,k,i),w%pl%tcval(1,k,i)),c_float)
            if (a >= 0.5_c_float * real(pi,c_float)) a = a - real(pi,c_float)
            if (a < -0.5_c_float * real(pi,c_float)) a = a + real(pi,c_float)

            ! the box around the turned text, all inside the plot
            hx = 0.5_c_float * (abs(cos(a)) * tsz%x + abs(sin(a)) * tsz%y) + pad
            hy = 0.5_c_float * (abs(sin(a)) * tsz%x + abs(cos(a)) * tsz%y) + pad
            bx = (/px - hx,py - hy,px + hx,py + hy/)
            if (bx(1) < pmin%x .or. bx(3) > pmax%x .or. bx(2) < pmin%y .or. bx(4) > pmax%y) cycle
            free = .true.
            do l = 1, nbox
               if (bx(1) < box(3,l) .and. bx(3) > box(1,l) .and. bx(2) < box(4,l) .and. bx(4) > box(2,l)) then
                  free = .false.
                  exit
               end if
            end do
            if (.not.free) cycle

            ! placed; the next values keep some distance from it
            nbox = nbox + 1
            box(:,nbox) = (/bx(1) - sep,bx(2) - sep,bx(3) + sep,bx(4) + sep/)
            nval = nval + 1
            ival(nval) = i
            icand(nval) = k
            sval(nval) = sc
            gap(nval) = g
            ang(nval) = a
            exit
         end do
      end do

    end subroutine place_values

    !> Plot the contour lines of array ia (1 = positive levels, in
    !> cpos; 2 = the rest, in cneg), the n rows of c, in color rgb, with
    !> a gap of half-length gap (plot units, along the line) at
    !> distance sval along the lines of the nval values placed (ival).
    subroutine plot_contours(ia,rgb,c,n,nval,ival,sval,gap)
      use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      integer, intent(in) :: ia
      real(c_float), intent(in) :: rgb(3)
      real(c_double), intent(in), target :: c(:,:)
      integer, intent(in) :: n, nval
      integer, intent(in) :: ival(:)
      real*8, intent(in) :: sval(:), gap(:)

      character(kind=c_char,len=:), allocatable, target :: str1
      real(c_double), allocatable, target :: d(:,:)
      integer, allocatable :: atrow(:)
      integer :: i, j, k, r, nd, n1, ncut, istate
      real*8 :: s, s1, ds, sa, sb, nan

      str1 = "##planec" // string(ia) // c_null_char
      call ipSetNextLineStyle(ImVec4(rgb(1),rgb(2),rgb(3),1._c_float),1._c_float)

      ! the values on lines of this array: none, plot the lines as they are
      ncut = 0
      allocate(atrow(n))
      atrow = 0
      do j = 1, nval
         k = ival(j)
         if (w%pl%icval(1,k) /= ia) cycle
         atrow(w%pl%icval(2,k)) = j
         ncut = ncut + 1
      end do
      if (ncut == 0) then
         call ipPlotLine(c_loc(str1),c_loc(c(1,1)),c_loc(c(1,2)),int(n,c_int),ImPlotLineFlags_None,0_c_int)
         return
      end if

      ! copy the lines, cutting them under their values
      nan = ieee_value(nan,ieee_quiet_nan)
      allocate(d(n+3*ncut,2))
      nd = 0
      r = 1
      do while (r <= n)
         if (atrow(r) == 0) then
            nd = nd + 1
            d(nd,:) = c(r,:)
            r = r + 1
            cycle
         end if
         j = atrow(r)
         k = ival(j)
         n1 = w%pl%icval(3,k)
         sa = sval(j) - gap(j)
         sb = sval(j) + gap(j)
         nd = nd + 1
         d(nd,:) = c(r,:)
         s = 0d0
         istate = 0
         do i = r+1, r+n1-1
            ds = norm2(c(i,:) - c(i-1,:))
            s1 = s + ds
            if (istate == 0 .and. s1 >= sa) then
               nd = nd + 1
               d(nd,:) = c(i-1,:) + (sa - s) / max(ds,1d-40) * (c(i,:) - c(i-1,:))
               nd = nd + 1
               d(nd,:) = nan
               istate = 1
            end if
            if (istate == 1 .and. s1 >= sb) then
               nd = nd + 1
               d(nd,:) = c(i-1,:) + (sb - s) / max(ds,1d-40) * (c(i,:) - c(i-1,:))
               istate = 2
            end if
            if (istate /= 1) then
               nd = nd + 1
               d(nd,:) = c(i,:)
            end if
            s = s1
         end do
         r = r + n1
      end do
      call ipPlotLine(c_loc(str1),c_loc(d(1,1)),c_loc(d(1,2)),int(nd,c_int),ImPlotLineFlags_None,0_c_int)

    end subroutine plot_contours

    !> A contour level written with three significant digits, as %.3g.
    function fmt_level(x) result(str)
      real*8, intent(in) :: x
      character(len=:), allocatable :: str

      character(len=40) :: aux, fmt
      integer :: e, idx

      if (x == 0d0) then
         str = "0"
         return
      end if
      e = floor(log10(abs(x)))
      if (e < -4 .or. e >= 3) then
         write (aux,'(ES0.2)') x
         idx = index(aux,"E")
         str = strip_zeros(aux(1:idx-1)) // "e" // trim(aux(idx+1:))
      else
         write (fmt,'("(F0.",I0,")")') max(2 - e,0)
         write (aux,fmt) x
         str = strip_zeros(trim(aux))
         if (str(1:1) == ".") then
            str = "0" // str
         elseif (str(1:min(2,len(str))) == "-.") then
            str = "-0" // str(2:)
         end if
      end if

    end function fmt_level

    !> Remove the trailing zeros of the decimals of a number, and the
    !> point if they are all zeros.
    function strip_zeros(s0) result(str)
      character(len=*), intent(in) :: s0
      character(len=:), allocatable :: str

      integer :: n

      str = trim(s0)
      if (index(str,".") == 0) return
      n = len(str)
      do while (n > 1 .and. str(n:n) == "0")
         n = n - 1
      end do
      if (str(n:n) == ".") n = n - 1
      str = str(1:n)

    end function strip_zeros

    !> The point of the region with plot coordinates (mu,mv): whether it
    !> is on the region (inside) and its position (x, cell-frame
    !> Cartesian bohr).
    subroutine plot_point(mu,mv,inside,x)
      real(c_double), intent(in) :: mu, mv
      logical, intent(out) :: inside
      real*8, intent(out) :: x(3)

      integer :: nx, ny
      real*8 :: d(2), m(2,2), det, s, t

      inside = .false.
      nx = size(w%pl%ff,1)
      ny = size(w%pl%ff,2)
      x = w%pl%a0 + (mu * w%pl%e1 + mv * w%pl%e2) / u

      ! the grid coordinates: x - x0 = s * (nx-1) * uu + t * (ny-1) * vv
      d(1) = dot_product(x - w%pl%x0,w%pl%e1)
      d(2) = dot_product(x - w%pl%x0,w%pl%e2)
      m(1,1) = dot_product(w%pl%uu,w%pl%e1) * (nx-1)
      m(2,1) = dot_product(w%pl%uu,w%pl%e2) * (nx-1)
      m(1,2) = dot_product(w%pl%vv,w%pl%e1) * (ny-1)
      m(2,2) = dot_product(w%pl%vv,w%pl%e2) * (ny-1)
      det = m(1,1) * m(2,2) - m(1,2) * m(2,1)
      if (abs(det) < 1d-20) return
      s = (d(1) * m(2,2) - d(2) * m(1,2)) / det
      t = (m(1,1) * d(2) - m(2,1) * d(1)) / det
      if (s < 0d0 .or. s > 1d0 .or. t < 0d0 .or. t > 1d0) return
      inside = .true.

    end subroutine plot_point

    !> The table of the values on the plane (opened and closed from
    !> its title): the position of each point (fractional for crystals, Å
    !> for molecules), its plot coordinates, and the value. The menu
    !> of the header (right click) writes it to the output, the
    !> console, or a file.
    subroutine draw_data_section()
      use utils, only: iw_table_column, iw_table_headers_row, iw_separatortext
      integer, parameter :: itable_plane = 1 ! the table, for table_write_*
      character(kind=c_char,len=:), allocatable, target :: strtab
      character(len=:), allocatable :: un
      integer(c_int) :: flags
      integer :: k, ix, iy, iwrite, nx, ny
      real*8 :: xp(3), xq(3), pc(2)
      type(ImVec2) :: szt
      type(c_ptr) :: clipper
      type(ImGuiListClipper), pointer :: clipper_f

      if (.not.allocated(w%pl%ff)) return
      call iw_separatortext("Values on the plane",w%pl%secdata)
      if (.not.w%pl%secdata%open) return

      nx = size(w%pl%ff,1)
      ny = size(w%pl%ff,2)
      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      flags = ior(flags,ImGuiTableFlags_ScrollY)
      flags = ior(flags,ImGuiTableFlags_ScrollX)
      strtab = "##planedatatable" // c_null_char
      szt%x = 0._c_float
      szt%y = htable
      un = lunit()
      if (igBeginTable(c_loc(strtab),6_c_int,flags,szt,0._c_float)) then
         if (ismol) then
            call iw_table_column("x (Å)",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("y (Å)",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("z (Å)",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         else
            call iw_table_column("x (frac)",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("y (frac)",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("z (frac)",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         end if
         call iw_table_column("u (" // un // ")",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("v (" // un // ")",id=4_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column(quantity_label(),id=5_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_headers_row(writemenu=iwrite,ttshown=ttshown,freezetop=.true.)

         ! the rows (the second index runs fastest, as in PLANE); if the
         ! table is being written as text, they are captured as they are
         ! drawn (all of them)
         clipper = ImGuiListClipper_ImGuiListClipper()
         call ImGuiListClipper_Begin(clipper,int(nx*ny,c_int),-1._c_float)
         call w%table_write_begin(itable_plane,iwrite,clipper)
         do while(ImGuiListClipper_Step(clipper))
            call c_f_pointer(clipper,clipper_f)
            do k = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
               ix = (k-1) / ny + 1
               iy = k - (ix-1) * ny
               call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
               xp = w%pl%x0 + real(ix-1,8) * w%pl%uu + real(iy-1,8) * w%pl%vv
               pc = to_plot(xp)
               where (abs(pc) < 1d-10) pc = 0d0
               if (ismol) then
                  xq = (xp + sys(isys)%c%molx0) * bohrtoa
                  where (abs(xq) < 1d-10) xq = 0d0
                  call cell(0,string(xq(1),'f',decimal=6))
                  call cell(1,string(xq(2),'f',decimal=6))
                  call cell(2,string(xq(3),'f',decimal=6))
               else
                  xq = sys(isys)%c%c2x(xp)
                  ! no negative zeros from the rounding of the plane
                  where (abs(xq) < 1d-10) xq = 0d0
                  call cell(0,string(xq(1),'f',decimal=7))
                  call cell(1,string(xq(2),'f',decimal=7))
                  call cell(2,string(xq(3),'f',decimal=7))
               end if
               call cell(3,string(pc(1),'f',decimal=6))
               call cell(4,string(pc(2),'f',decimal=6))
               call cell(5,string(w%pl%ff(ix,iy),'e',decimal=8))
            end do
         end do
         call ImGuiListClipper_End(clipper)
         call ImGuiListClipper_destroy(clipper)
         call w%table_write_end(itable_plane,table_title("Plane map",isys))
         call igEndTable()
      end if

    end subroutine draw_data_section

    !> Text str in column icol of the current row of a table.
    subroutine cell(icol,str)
      integer, intent(in) :: icol
      character(len=*), intent(in) :: str

      if (igTableSetColumnIndex(int(icol,c_int))) call iw_text(str)

    end subroutine cell

    !> The label of the quantity: $n with the property, or the
    !> expression.
    function quantity_label() result(str)
      character(len=:), allocatable :: str

      if (w%pl%d%isexpr) then
         str = trim(w%pl%d%expr)
      else
         str = "$" // string(w%pl%d%fref%id) // trim(quant_short(w%pl%d%nti))
      end if

    end function quantity_label

    !> Calculate the grid of the region and the values on it if the
    !> points, the region, the number of points, the quantity, the
    !> fields, or the geometry changed. Waits for the initialization
    !> threads to finish and for the inputs to be no longer edited.
    !> refit = the region changed (the plot fits the axes to it).
    !> fresh = the grid and the values (or their absence) are those of
    !> the form (not while it is edited or the threads run).
    subroutine update_grid(refit,fresh)
      use systems, only: are_threads_running
      logical, intent(out) :: refit, fresh

      integer*8, allocatable :: key(:)
      logical :: stale
      real*8 :: oldx0(3), olduu(3), oldvv(3), olda0(3)
      integer :: oldnx, oldny

      refit = .false.
      fresh = .false.
      if (xediting .or. are_threads_running()) return
      fresh = .true.

      ! everything the values depend on, and the system
      call w%pl%d%autopoints(isys)
      call w%pl%d%key(isys,key)
      key = [key,w%pl%sys%uid]
      stale = .not.allocated(w%pl%key) .or. .not.allocated(w%pl%keyexpr)
      if (.not.stale) stale = (size(w%pl%key) /= size(key))
      if (.not.stale) stale = any(w%pl%key /= key) .or. (w%pl%keyexpr /= trim(w%pl%d%expr))
      if (.not.stale) return
      call move_alloc(key,w%pl%key)
      w%pl%keyexpr = trim(w%pl%d%expr)

      ! the region before, to know whether the plot must fit it again
      oldx0 = w%pl%x0
      olda0 = w%pl%a0
      olduu = w%pl%uu
      oldvv = w%pl%vv
      oldnx = 0
      oldny = 0
      if (allocated(w%pl%ff)) then
         oldnx = size(w%pl%ff,1)
         oldny = size(w%pl%ff,2)
         deallocate(w%pl%ff)
      end if
      w%pl%errmsg = ""
      w%pl%gen = w%pl%gen + 1
      if (sys(isys)%c%ncel == 0) return

      ! the grid and the frame of the plot
      call w%pl%d%grid(isys,w%pl%x0,w%pl%uu,w%pl%vv,w%pl%a0,w%pl%e1,w%pl%e2,w%pl%errmsg)
      if (len_trim(w%pl%errmsg) > 0) return
      ! (relative to a0, the origin of the plot: moving the plane as
      ! a whole keeps the axes)
      refit = (oldnx == 0) .or. any(abs(oldx0 - olda0 - w%pl%x0 + w%pl%a0) > 1d-10) .or. &
         any(abs(olduu*(oldnx-1) - w%pl%uu*(w%pl%d%nx-1)) > 1d-10) .or.&
         any(abs(oldvv*(oldny-1) - w%pl%vv*(w%pl%d%ny-1)) > 1d-10)

      ! the values
      call w%pl%d%values(isys,w%pl%x0,w%pl%uu,w%pl%vv,w%pl%ff,w%pl%errmsg)

    end subroutine update_grid

    !> Calculate the contour lines again if the grid or the levels
    !> changed, in plot coordinates, the positive levels in cpos and
    !> the rest in cneg (NaN between two lines), and the places along
    !> each line where its value can go, longest line first.
    subroutine update_contours()
      use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      use tools, only: mergesort
      integer*8, allocatable :: key(:)
      integer :: i, ii, j, k, niso, nline, np
      integer, allocatable :: lptr(:), liso(:), iord(:), il(:,:)
      real*8, allocatable :: ziso(:), xy(:,:), ll(:)
      real*8 :: off(2), rval, nan, dl, acc, st
      logical :: pos

      call w%pl%d%contkey(key)
      key = [key,w%pl%gen]
      if (allocated(w%pl%keycont) .and. allocated(w%pl%keylist)) then
         if (size(w%pl%keycont) == size(key)) then
            if (all(w%pl%keycont == key) .and. w%pl%keylist == trim(w%pl%d%llist)) return
         end if
      end if
      call move_alloc(key,w%pl%keycont)
      w%pl%keylist = trim(w%pl%d%llist)
      w%pl%ncpos = 0
      w%pl%ncneg = 0
      w%pl%ncval = 0
      if (.not.w%pl%d%showcont) return

      ! the lines, in the frame of the region; to plot coordinates
      call w%pl%d%contours(w%pl%ff,w%pl%uu,w%pl%vv,niso,ziso,nline,lptr,xy,liso)
      if (nline == 0) return
      off = to_plot(w%pl%x0)
      nan = ieee_value(nan,ieee_quiet_nan)
      if (allocated(w%pl%cpos)) deallocate(w%pl%cpos)
      if (allocated(w%pl%cneg)) deallocate(w%pl%cneg)
      allocate(w%pl%cpos(size(xy,2)+nline,2),w%pl%cneg(size(xy,2)+nline,2))
      allocate(il(3,nline))
      do i = 1, nline
         pos = (ziso(liso(i)) > 0d0)
         np = merge(w%pl%ncpos,w%pl%ncneg,pos)
         il(:,i) = (/merge(1,2,pos),np+1,lptr(i+1)-lptr(i)/)
         do j = lptr(i), lptr(i+1)-1
            np = np + 1
            do k = 1, 2
               rval = xy(k,j) * u + off(k)
               if (pos) then
                  w%pl%cpos(np,k) = rval
               else
                  w%pl%cneg(np,k) = rval
               end if
            end do
         end do
         np = np + 1
         if (pos) then
            w%pl%cpos(np,:) = nan
            w%pl%ncpos = np
         else
            w%pl%cneg(np,:) = nan
            w%pl%ncneg = np
         end if
      end do

      ! the length of each line, longest first, and the places along
      ! it where its value can go, with the direction of the line there
      allocate(ll(nline),iord(nline))
      do i = 1, nline
         ll(i) = 0d0
         do j = lptr(i)+1, lptr(i+1)-1
            ll(i) = ll(i) + norm2(xy(:,j) - xy(:,j-1))
         end do
         iord(i) = i
      end do
      call mergesort(-ll,iord,1,nline)
      if (allocated(w%pl%xcval)) deallocate(w%pl%xcval)
      if (allocated(w%pl%tcval)) deallocate(w%pl%tcval)
      if (allocated(w%pl%lcval)) deallocate(w%pl%lcval)
      if (allocated(w%pl%zcval)) deallocate(w%pl%zcval)
      if (allocated(w%pl%icval)) deallocate(w%pl%icval)
      allocate(w%pl%xcval(2,plane_ncand,nline),w%pl%tcval(2,plane_ncand,nline),w%pl%lcval(nline),&
         w%pl%zcval(nline),w%pl%icval(4,nline))
      w%pl%ncval = nline
      do i = 1, nline
         ii = iord(i)
         do k = 1, plane_ncand
            w%pl%xcval(:,k,i) = xy(:,lptr(ii))
            w%pl%tcval(:,k,i) = (/1d0,0d0/)
         end do
         acc = 0d0
         k = 1
         do j = lptr(ii)+1, lptr(ii+1)-1
            dl = norm2(xy(:,j) - xy(:,j-1))
            do while (k <= plane_ncand .and. dl > 0d0)
               st = (k - 0.5d0) / plane_ncand * ll(ii)
               if (acc + dl < st) exit
               w%pl%tcval(:,k,i) = (xy(:,j) - xy(:,j-1)) / dl
               w%pl%xcval(:,k,i) = xy(:,j-1) + (st - acc) * w%pl%tcval(:,k,i)
               k = k + 1
            end do
            acc = acc + dl
         end do
         w%pl%xcval(:,:,i) = w%pl%xcval(:,:,i) * u + spread(off,2,plane_ncand)
         w%pl%lcval(i) = ll(ii) * u
         w%pl%zcval(i) = ziso(liso(ii))
         w%pl%icval(1:3,i) = il(:,ii)
         w%pl%icval(4,i) = liso(ii)
      end do

    end subroutine update_contours

    !> Make the texture of the map again if the grid or the colors
    !> changed, and keep the range of its colors in tmin, tmax.
    subroutine update_texture()
      use interfaces_opengl3
      integer*8, allocatable :: key(:)
      integer :: k, nx, ny
      real(c_float), allocatable :: rgb(:,:)
      character(kind=c_char), allocatable, target :: pix(:)

      call w%pl%d%colorkey(key)
      key = [key,w%pl%gen]
      if (allocated(w%pl%keytex) .and. w%pl%tex /= 0) then
         if (size(w%pl%keytex) == size(key)) then
            if (all(w%pl%keytex == key)) return
         end if
      end if
      call move_alloc(key,w%pl%keytex)

      ! the colors, as RGBA bytes (the first index of the grid runs
      ! fastest, along the width of the texture), and their range
      nx = size(w%pl%ff,1)
      ny = size(w%pl%ff,2)
      call w%pl%d%colors(w%pl%ff,w%pl%tmin,w%pl%tmax,rgb)
      allocate(pix(4*nx*ny))
      do k = 1, nx*ny
         pix(4*k-3) = char(nint(rgb(1,k) * 255))
         pix(4*k-2) = char(nint(rgb(2,k) * 255))
         pix(4*k-1) = char(nint(rgb(3,k) * 255))
         pix(4*k) = char(255)
      end do

      ! upload
      if (w%pl%tex == 0) call glGenTextures(1,c_loc(w%pl%tex))
      call glBindTexture(GL_TEXTURE_2D,w%pl%tex)
      call glPixelStorei(GL_UNPACK_ALIGNMENT,1)
      call glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,int(nx,c_int),int(ny,c_int),0,GL_RGBA,&
         GL_UNSIGNED_BYTE,c_loc(pix))
      call glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR)
      call glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR)
      call glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE)
      call glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE)
      call glBindTexture(GL_TEXTURE_2D,0)

    end subroutine update_texture

    !> Find the atoms, and the critical points of the reference field
    !> other than the nuclei, closer to the plane than hmax and inside
    !> the region (their lattice images, for crystals), and keep their
    !> plot coordinates and labels (with their signed distance to the
    !> plane, if showdist).
    subroutine update_marks()
      use tools_math, only: matinv, cross
      use types, only: realloc
      integer*8 :: key(8)
      integer :: k, l1, l2, l3, lo(3), hi(3), nc, ier, nx, ny
      real*8 :: am(3,3), nrm(3), xq(3), s(3), fc(3,4), hb, eps
      real*8, allocatable :: xf(:,:)
      character(len=12), allocatable :: lbl(:)
      character(len=:), allocatable :: sdist

      key(1) = w%pl%gen
      key(2) = merge(1_8,0_8,w%pl%showlabels)
      key(6) = merge(1_8,0_8,w%pl%showcps)
      key(7) = w%pl%ltext
      key(8) = merge(1_8,0_8,w%pl%showdist)
      key(3) = transfer(w%pl%hmax,key(3))
      key(4) = transfer(sysc(isys)%timelastchange_cplist,key(4))
      key(5) = sys(isys)%field_cuid(sys(isys)%iref)
      if (allocated(w%pl%keymark)) then
         if (all(w%pl%keymark == key)) return
      end if
      w%pl%keymark = key
      w%pl%nmark = 0
      w%pl%marksmany = .false.
      if (allocated(w%pl%xmark)) deallocate(w%pl%xmark)
      if (allocated(w%pl%lmark)) deallocate(w%pl%lmark)
      allocate(w%pl%xmark(2,10),w%pl%lmark(10))
      if (.not.w%pl%showlabels) return

      ! the candidates: the atoms, and the critical points of the
      ! reference field except the nuclei (crystallographic)
      call label_candidates(isys,nc,xf,lbl,cps=w%pl%showcps,short=.true.,atext=w%pl%ltext)

      ! the coordinates of a point in the frame of the region: x - x0 =
      ! s1 * side1 + s2 * side2 + s3 * normal
      nx = size(w%pl%ff,1)
      ny = size(w%pl%ff,2)
      nrm = cross(w%pl%e1,w%pl%e2)
      am(:,1) = w%pl%uu * (nx-1)
      am(:,2) = w%pl%vv * (ny-1)
      am(:,3) = nrm
      call matinv(am,3,ier)
      if (ier /= 0) return
      hb = w%pl%hmax / u
      eps = 1d-6

      ! the lattice translations that can bring a point of the cell to
      ! the region (crystals; only the box for molecules), if there are
      ! not too many
      lo = 0
      hi = 0
      if (.not.ismol) then
         fc(:,1) = sys(isys)%c%c2x(w%pl%x0)
         fc(:,2) = sys(isys)%c%c2x(w%pl%x0 + w%pl%uu * (nx-1))
         fc(:,3) = sys(isys)%c%c2x(w%pl%x0 + w%pl%vv * (ny-1))
         fc(:,4) = sys(isys)%c%c2x(w%pl%x0 + w%pl%uu * (nx-1) + w%pl%vv * (ny-1))
         lo = floor(minval(fc,2)) - 1
         hi = ceiling(maxval(fc,2)) + 1
         if (product(real(hi - lo + 1,8)) > maxtrans) then
            w%pl%marksmany = .true.
            return
         end if
      end if

      do k = 1, nc
         do l1 = lo(1), hi(1)
            do l2 = lo(2), hi(2)
               do l3 = lo(3), hi(3)
                  xq = sys(isys)%c%x2c(xf(:,k) + real((/l1,l2,l3/),8))
                  s = matmul(am,xq - w%pl%x0)
                  if (s(1) < -eps .or. s(1) > 1d0+eps .or. s(2) < -eps .or. s(2) > 1d0+eps) cycle
                  if (abs(s(3)) > hb + eps) cycle
                  w%pl%nmark = w%pl%nmark + 1
                  if (w%pl%nmark > size(w%pl%lmark)) then
                     call realloc(w%pl%xmark,2,2*w%pl%nmark)
                     call realloc(w%pl%lmark,2*w%pl%nmark)
                  end if
                  w%pl%xmark(:,w%pl%nmark) = to_plot(xq)
                  w%pl%lmark(w%pl%nmark) = lbl(k)
                  ! the signed distance to the plane, as C(+0.01), unless
                  ! it is zero to two decimals
                  if (w%pl%showdist) then
                     sdist = string(s(3) * u,'f',decimal=2)
                     if (sdist /= "0.00" .and. sdist /= "-0.00") then
                        if (s(3) > 0d0) sdist = "+" // sdist
                        w%pl%lmark(w%pl%nmark) = trim(lbl(k)) // "(" // sdist // ")"
                     end if
                  end if
               end do
            end do
         end do
      end do

    end subroutine update_marks

    !> The plot coordinates (Å for molecules, bohr for crystals) of the
    !> point x (cell-frame Cartesian bohr): along e1 and e2 from the
    !> origin the grid was calculated for (a0).
    function to_plot(x) result(pc)
      real*8, intent(in) :: x(3)
      real*8 :: pc(2)

      pc(1) = dot_product(x - w%pl%a0,w%pl%e1) * u
      pc(2) = dot_product(x - w%pl%a0,w%pl%e2) * u

    end function to_plot

    !> The unit of lengths in the window: Å for molecules, bohr for
    !> crystals.
    function lunit() result(str)
      character(len=:), allocatable :: str
      if (ismol) then
         str = "Å"
      else
         str = "bohr"
      end if
    end function lunit

    !> The point x (cell-frame Cartesian bohr) in the absolute frame
    !> (with molx0 for molecules), for the scene.
    function cart_to_abs(x) result(xa)
      real*8, intent(in) :: x(3)
      real*8 :: xa(3)
      xa = x
      if (ismol) xa = xa + sys(isys)%c%molx0
    end function cart_to_abs

    !> Forget a pending pick and release the view it was armed on.
    subroutine cancel_pick()
      if (.not.w%pl%picking) return
      if (w%pl%pickview >= 1 .and. w%pl%pickview <= nwin) &
         call win(w%pl%pickview)%viewmode_release_forced(w%id)
      w%pl%picking = .false.
    end subroutine cancel_pick

  end subroutine draw_plane

  !> The form of the map d of system isys, shared by the plane map
  !> window and the editor of the plane map objects: the three points,
  !> the region, the margin, the number of points, the buttons that
  !> move the plane, the quantity, the
  !> colors, and the contour lines. Returns whether anything changed.
  !> editing becomes true if a value is being dragged. If live, a drag
  !> is a change in every frame. If
  !> ipick is present, each point has a Pick button: ipick is the one
  !> pressed (0 = none), and pickoff disables them. Not live, a drag
  !> counts as a change only when it ends (the map of an object is
  !> calculated again on each change). If showval is present, the
  !> contour line has a Values checkbox for it (the labels of the
  !> levels in the window, which are not part of the map), and, if
  !> valscale is present too, a control for their size.
  module function plane_def_form(d,isys,ttshown,editing,live,ipick,pickoff,showval,valscale,sec) result(changed)
    use systems, only: sys
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_dragfloat_realc, iw_calcwidth,&
       iw_inputint, iw_combo_simple, iw_checkbox, iw_coloredit, iw_inputtext, iw_cmap_optstr,&
       iw_field_combo, iw_expression_helper, iw_exprctx_point, iw_separatortext
    use tools_io, only: string
    type(plane_def), intent(inout) :: d
    integer, intent(in) :: isys
    logical, intent(inout) :: ttshown
    logical, intent(inout) :: editing
    logical, intent(in) :: live
    integer, intent(out), optional :: ipick
    logical, intent(in), optional :: pickoff
    logical, intent(inout), optional :: showval
    real(c_float), intent(inout), optional :: valscale
    type(iw_section), intent(inout) :: sec(plane_nsec)
    logical :: changed

    integer :: i
    logical :: ismol, ch, pickoff_
    character(len=:), allocatable :: un

    changed = .false.
    if (present(ipick)) ipick = 0
    pickoff_ = .false.
    if (present(pickoff)) pickoff_ = pickoff
    ismol = sys(isys)%c%ismolecule
    if (ismol) then
       un = "Å"
    else
       un = "bohr"
    end if

    ! the plane, the region, and the number of points
    call iw_separatortext("Plane",sec(plane_sec_plane),hasmore=.true.)
    if (sec(plane_sec_plane)%open) then
       do i = 1, 3
          call draw_point(i)
       end do
       call draw_move()
       if (sec(plane_sec_plane)%more) call draw_region()
    end if

    ! the quantity (nothing more for an expression), the colors, and
    ! the contour lines (nothing more if they are not shown)
    call iw_separatortext("Property",sec(plane_sec_property),hasmore=.not.d%isexpr)
    if (sec(plane_sec_property)%open) call draw_quantity()
    call iw_separatortext("Colors",sec(plane_sec_colors),hasmore=.true.)
    if (sec(plane_sec_colors)%open) call draw_colors()
    call iw_separatortext("Contours",sec(plane_sec_contours),hasmore=d%showcont)
    if (sec(plane_sec_contours)%open) call draw_contours()

  contains
    !> A drag of a real value (x1 or the three of x3), live or not.
    subroutine drag(str,x1,x3,speed,decimal,sameline,min)
      character(len=*), intent(in) :: str
      real*8, intent(inout), optional :: x1, x3(3)
      real*8, intent(in) :: speed
      integer, intent(in) :: decimal
      logical, intent(in), optional :: sameline
      real*8, intent(in), optional :: min


      logical :: ch

      ch = iw_dragfloat_real8(str,x1=x1,x3=x3,speed=speed,decimal=decimal,sameline=sameline,min=min)
      editing = editing .or. igIsItemActive()
      ! not live: the change counts when the drag (or the typing) ends
      if (live) then
         changed = changed .or. ch
      else
         changed = changed .or. igIsItemDeactivated()
      end if

    end subroutine drag

    !> A fixed range of values: the checkbox flag and, if it is on,
    !> the drags of xmin and xmax, with the tooltips of the three.
    !> sfx = the ID suffix of the three widgets.
    subroutine draw_range(sfx,flag,xmin,xmax,ttflag,ttmin,ttmax)
      character(len=*), intent(in) :: sfx, ttflag, ttmin, ttmax
      logical, intent(inout) :: flag
      real*8, intent(inout) :: xmin, xmax

      changed = changed .or. iw_checkbox("Fixed range##plane" // sfx // "range",flag,labelleft=.true.)
      call iw_tooltip(ttflag,ttshown)
      if (.not.flag) return
      call iw_text("Min",alignframe=.true.,gap=.true.)
      call drag("##plane" // sfx // "min",x1=xmin,speed=0.01d0,decimal=4,sameline=.true.)
      call iw_tooltip(ttmin,ttshown)
      call iw_text("Max",alignframe=.true.,gap=.true.)
      call drag("##plane" // sfx // "max",x1=xmax,speed=0.01d0,decimal=4,sameline=.true.)
      call iw_tooltip(ttmax,ttshown)

    end subroutine draw_range

    !> One of the three points of the plane (ip = 1 origin, 2 end of
    !> the first side, 3 point on the second side): the label, the
    !> button to pick it in the view, and its coordinates.
    subroutine draw_point(ip)
      integer, intent(in) :: ip

      character(len=:), allocatable :: sfx

      sfx = "##planept" // string(ip)
      call iw_text(pt_name(ip),rgb=pt_rgb(:,ip),alignframe=.true.)
      if (ip == 1) then
         call iw_tooltip("The origin of the plane (shown in the view in this color)",ttshown)
      elseif (ip == 2) then
         call iw_tooltip("The end of the first side of the plane, which starts at the origin&
            & (shown in the view in this color)",ttshown)
      else
         call iw_tooltip("A point on the second side of the plane: its end, for the parallelogram;&
            & the plane through it, for the rectangle (shown in the view in this color)",ttshown)
      end if
      if (present(ipick)) then
         if (iw_button("Pick" // sfx // "pick",sameline=.true.,disabled=pickoff_)) ipick = ip
         call iw_tooltip("Pick this point in the view: an atom, a critical point, or a point&
            & in empty space",ttshown)
      end if
      if (ismol) then
         call drag("(Å)" // sfx,x3=d%xp(:,ip),speed=0.01d0,decimal=5,sameline=.true.)
      else
         call drag("(frac)" // sfx,x3=d%xp(:,ip),speed=0.001d0,decimal=6,sameline=.true.)
      end if
      call iw_tooltip("Coordinates of this point",ttshown)

    end subroutine draw_point

    !> The buttons that displace the whole plane by d%step along the
    !> Cartesian axes and, in crystals, along the lattice vectors.
    subroutine draw_move()
      character(len=1), parameter :: lbl(6) = (/"x","y","z","a","b","c"/)
      integer :: j, k
      real*8 :: dx(3), st
      character(len=1) :: sg

      call iw_text("Move",alignframe=.true.)
      call iw_tooltip("Displace the whole plane",ttshown)
      do j = 1, merge(3,6,ismol)
         do k = 1, 2
            sg = merge("+","-",k == 1)
            st = merge(d%step,-d%step,k == 1)
            if (iw_button(sg // lbl(j) // "##planemove",sameline=.true.)) then
               dx = 0d0
               if (j <= 3) then
                  dx(j) = st
                  ! crystals: fractional points, displacement in bohr
                  if (.not.ismol) dx = sys(isys)%c%c2x(dx)
               else
                  dx(j-3) = st / sys(isys)%c%aa(j-3)
               end if
               d%xp = d%xp + spread(dx,2,3)
               changed = .true.
            end if
            call iw_tooltip("Displace the plane along " // sg // lbl(j) // " (" //&
               merge("Cartesian axis","lattice vector",j <= 3) // ")",ttshown)
         end do
      end do
      call iw_text("Step",alignframe=.true.,gap=.true.)
      call igPushItemWidth(iw_calcwidth(7,1))
      ch = iw_dragfloat_real8("##planestep",x1=d%step,speed=0.01d0,decimal=3,sameline=.true.,min=0d0)
      call igPopItemWidth()
      call iw_tooltip("Length of each displacement (" // un // ")",ttshown)

    end subroutine draw_move

    !> More of the plane: the region (a rectangle or a parallelogram)
    !> and the margin around it, then the number of points along each
    !> side (shown, but not editable, if automatic).
    subroutine draw_region()

      call iw_text("Region",alignframe=.true.)
      call iw_combo_simple("##planeregion","Rectangle" // c_null_char // "Parallelogram" // c_null_char,&
         d%region,sameline=.true.,changed=ch)
      changed = changed .or. ch
      call iw_tooltip("Rectangle: the first side goes from the origin to the X end, and the&
         & second is perpendicular to it, as long as needed to reach the Y point. Parallelogram:&
         & the sides go from the origin to the X end and to the Y point, as in PLANE",ttshown)
      call iw_text("Margin",alignframe=.true.,gap=.true.)
      call drag("##planemargin",x1=d%margin,speed=0.01d0,decimal=3,sameline=.true.)
      call iw_tooltip("Length added to the region on every side (" // un //&
         "; negative removes it)",ttshown)

      call iw_text("Points",alignframe=.true.)
      if (d%autonpts) call igBeginDisabled(.true._c_bool)
      ch = iw_inputint("##planenx",d%nx,width=6,sameline=.true.,notlive=.true.)
      call iw_tooltip("Number of points along the first side (" // string(plane_npts_min) // " to " //&
         string(plane_npts_max) // ")",ttshown,whendisabled=.true.)
      call iw_text("x",alignframe=.true.,sameline=.true.)
      ch = iw_inputint("##planeny",d%ny,width=6,sameline=.true.,notlive=.true.) .or. ch
      call iw_tooltip("Number of points along the second side (" // string(plane_npts_min) // " to " //&
         string(plane_npts_max) // ")",ttshown,whendisabled=.true.)
      if (d%autonpts) then
         call igEndDisabled()
      else
         d%nx = max(min(d%nx,plane_npts_max),plane_npts_min)
         d%ny = max(min(d%ny,plane_npts_max),plane_npts_min)
         changed = changed .or. ch
      end if
      ch = iw_checkbox("Automatic##planeautonpts",d%autonpts,gap=.true.,labelleft=.true.)
      changed = changed .or. ch
      call iw_tooltip("Choose the number of points from the size of the region: square cells&
         & of " // string(plane_auto_spacing,'f',decimal=2) // " Å, coarser for very large&
         & regions (at most " // string(plane_auto_ntot) // " points) and for slow fields (so that&
         & a map takes at most about " // string(nint(plane_auto_tmax)) // " s)",ttshown)

    end subroutine draw_region

    !> The quantity on the map: a field, or an arithmetic expression
    !> (checked as it is typed). More: the property of the field.
    subroutine draw_quantity()
      integer :: ikind, ifield, iq
      logical :: echanged
      character(len=:), allocatable :: lerr

      call iw_text("Kind",alignframe=.true.)
      ikind = merge(1,0,d%isexpr)
      call iw_combo_simple("##planekind","Field" // c_null_char // "Expression" // c_null_char,ikind,&
         sameline=.true.)
      call iw_tooltip("A quantity of a field, or an arithmetic expression of the fields",ttshown)
      changed = changed .or. (d%isexpr .neqv. (ikind == 1))
      d%isexpr = (ikind == 1)
      if (.not.d%isexpr) then
         call iw_text("Field",alignframe=.true.,gap=.true.)
         ifield = d%fref%id
         if (iw_field_combo("##planefield",isys,ifield,width=iw_calcwidth(24,1),&
            nonestr="<field not available>",sameline=.true.)) then
            call d%fref%set(sys(isys),ifield)
            changed = .true.
         end if
         call iw_tooltip("Field on the map",ttshown)
         if (sec(plane_sec_property)%more) then
            call iw_text("Quantity",alignframe=.true.)
            iq = d%nti
            call iw_combo_simple("##planequant",quant_combo,iq,sameline=.true.)
            changed = changed .or. (iq /= d%nti)
            d%nti = max(min(iq,nquant-1),0)
            call iw_tooltip("Property of the field on the map",ttshown)
         end if
      else
         call iw_text("Expression",alignframe=.true.,gap=.true.)
         echanged = iw_inputtext("##planeexpr",bufsize=1023,textf=d%expr,width=36,sameline=.true.)
         call iw_tooltip("Arithmetic expression of the fields of the system (for instance,&
            & ""$1 - $2"" or ""$1:l""), evaluated at the points of the plane",ttshown)
         if (iw_expression_helper("##planehelp",isys,iw_exprctx_point,textf=d%expr,&
            ttshown=ttshown,ifield=sys(isys)%iref)) echanged = .true.
         if (echanged) then
            changed = .true.
            d%exprerr = ""
            if (len_trim(d%expr) > 0) then
               call sys(isys)%check_expression(trim(d%expr),lerr)
               d%exprerr = lerr
            end if
         end if
         if (len_trim(d%expr) == 0) then
            call iw_text("The expression is empty",danger=.true.,wrap=.true.)
         elseif (len_trim(d%exprerr) > 0) then
            call iw_text(trim(d%exprerr),danger=.true.,wrap=.true.)
         end if
      end if

    end subroutine draw_quantity

    !> The colors of the map. More: a fixed range of values for them.
    subroutine draw_colors()

      call iw_text("Colormap",alignframe=.true.)
      call iw_combo_simple("##planecmap",iw_cmap_optstr,d%icmap,startsatone=.true.,sameline=.true.,changed=ch)
      changed = changed .or. ch
      call iw_tooltip("Colormap of the map",ttshown)
      call iw_text("Scale",alignframe=.true.,gap=.true.)
      call iw_combo_simple("##planecscale","Linear" // c_null_char // "Log |f|" // c_null_char //&
         "Arctangent" // c_null_char,d%cscale,sameline=.true.,changed=ch)
      changed = changed .or. ch
      call iw_tooltip("How the values become colors: linearly, by the logarithm of their&
         & absolute value, or by their arctangent (as in PLANE COLORMAP)",ttshown)
      if (.not.sec(plane_sec_colors)%more) return

      ! more: the range
      call draw_range("c",d%crange,d%cmin,d%cmax,&
         "Colors for a fixed range of values (or the range of the values on the map)",&
         "The value of the first color","The value of the last color")

    end subroutine draw_colors

    !> The contour lines of the map: whether they are drawn, and the
    !> levels. More: a fixed range for the levels, the colors of the
    !> lines, and the values on them.
    subroutine draw_contours()
      integer :: ic

      changed = changed .or. iw_checkbox("Show##planeshowcont",d%showcont,labelleft=.true.)
      call iw_tooltip("Draw contour lines on the map",ttshown)
      if (.not.d%showcont) return

      ! the levels
      ic = 0
      do while (ic < 4 .and. ctype_of_combo(ic) /= d%ctype)
         ic = ic + 1
      end do
      call iw_text("Levels",alignframe=.true.,gap=.true.)
      call iw_combo_simple("##planectype","Linear" // c_null_char // "Logarithmic" // c_null_char //&
         "Arctangent" // c_null_char // "Bader" // c_null_char // "List" // c_null_char,ic,sameline=.true.)
      changed = changed .or. (d%ctype /= ctype_of_combo(ic))
      d%ctype = ctype_of_combo(ic)
      call iw_tooltip("The contour levels, as in PLANE CONTOUR: evenly spaced, logarithmic&
         & (positive and negative), evenly spaced in the arctangent, the Bader set (1, 2, 4,&
         & 8 x 10^n), or a list of values",ttshown)
      if (d%ctype == niso_manual) then
         call iw_text("List",alignframe=.true.,gap=.true.)
         changed = changed .or. iw_inputtext("##planellist",bufsize=1023,textf=d%llist,width=30,sameline=.true.)
         call iw_tooltip("The contour levels, separated by spaces or commas",ttshown)
      elseif (d%ctype /= niso_bader) then
         call iw_text("Number",alignframe=.true.,gap=.true.)
         ch = iw_inputint("##planenlev",d%nlev,width=5,sameline=.true.,notlive=.true.)
         d%nlev = max(min(d%nlev,plane_nlev_max),2_c_int)
         changed = changed .or. ch
         call iw_tooltip("Number of contour levels (at least 2)",ttshown)
      end if
      if (.not.sec(plane_sec_contours)%more) return

      ! more: the range of the levels
      if (d%ctype /= niso_manual .and. d%ctype /= niso_bader) &
         call draw_range("l",d%lrange,d%lmin,d%lmax,&
         "Contour levels in a fixed range of values (or the range of the values on the map)",&
         "The lowest contour level","The highest contour level")

      ! more: the colors of the lines, and the values on them
      call iw_text("Positive color",alignframe=.true.)
      changed = changed .or. iw_coloredit("##planergbpos",rgb=d%rgbpos,sameline=.true.)
      call iw_tooltip("Color of the contour lines of the positive levels",ttshown)
      call iw_text("Negative color",alignframe=.true.,gap=.true.)
      changed = changed .or. iw_coloredit("##planergbneg",rgb=d%rgbneg,sameline=.true.)
      call iw_tooltip("Color of the contour lines of the zero and negative levels",ttshown)
      if (present(showval)) then
         ch = iw_checkbox("Values##planeshowcval",showval,gap=.true.,labelleft=.true.)
         call iw_tooltip("Write the value of the level on the contour lines, in their color,&
            & where there is room",ttshown)
         if (showval .and. present(valscale)) then
            call iw_text("Size",alignframe=.true.,gap=.true.)
            ch = iw_dragfloat_realc("##planevscale",x1=valscale,speed=0.01_c_float,min=0.1_c_float,&
               max=10._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
            call iw_tooltip("Size of the values, relative to the text of the interface",ttshown)
         end if
      end if

    end subroutine draw_contours
  end function plane_def_form

  !> The editor of a plane map object (reptype_slice): the map, its
  !> opacity, the radius of the contour lines, and what is wrong with
  !> it. Returns whether anything changed.
  module function draw_editrep_slice(w,ttshown) result(changed)
    use utils, only: iw_text, iw_tooltip, iw_dragfloat_real8, iw_dragfloat_realc, iw_separatortext
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: editing

    editing = .false.
    changed = plane_def_form(w%rep%slice%d,w%rep%id,ttshown,editing,.false.,sec=w%editrep_slicesec(1:plane_nsec))

    ! how it is drawn
    call iw_separatortext("Style",w%editrep_slicesec(plane_nsec+1))
    if (.not.w%editrep_slicesec(plane_nsec+1)%open) return
    call iw_text("Opacity",alignframe=.true.)
    changed = changed .or. iw_dragfloat_realc("##sliceopacity",x1=w%rep%slice%alpha,speed=0.01_c_float,&
       min=0._c_float,max=1._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Opacity of the map (1 = opaque)",ttshown)
    call iw_text("Line radius",alignframe=.true.,gap=.true.)
    changed = changed .or. iw_dragfloat_real8("(Å)##slicecrad",x1=w%rep%slice%crad,speed=0.001d0,&
       min=0d0,scale=bohrtoa,decimal=3,sameline=.true.)
    call iw_tooltip("Radius of the contour lines",ttshown)

    ! what is wrong
    if (allocated(w%rep%slice%c)) then
       if (allocated(w%rep%slice%c%errmsg)) then
          if (len_trim(w%rep%slice%c%errmsg) > 0) &
             call iw_text(w%rep%slice%c%errmsg,danger=.true.,wrap=.true.)
       end if
    end if

  end function draw_editrep_slice

  !> Delete the texture of the map of the plane map window w.
  module subroutine plane_free_texture(w)
    use interfaces_opengl3
    class(window), intent(inout), target :: w

    if (w%pl%tex /= 0) call glDeleteTextures(1,c_loc(w%pl%tex))
    w%pl%tex = 0

  end subroutine plane_free_texture

end submodule plane
