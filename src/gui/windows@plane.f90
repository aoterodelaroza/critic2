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
  use representations, only: plane_npts_min, plane_npts_max, plane_nlev_max
  implicit none

  ! transient markers in the view: the three points (origin, end of
  ! the first side, point on the second side), the region, and the
  ! point under the mouse in the map; opacity of the markers and of
  ! the region, and radius of the points (bohr), a fraction of the
  ! scene radius between the two limits
  real(c_float), parameter :: pt_rgb(3,3) = reshape((/&
     0.95_c_float,0.80_c_float,0.10_c_float,&
     0.95_c_float,0.35_c_float,0.20_c_float,&
     0.20_c_float,0.55_c_float,0.95_c_float/),(/3,3/))
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
    real*8 :: xc(3), rad, v(3,3), u, xhov(3), fhov
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

    ! a new system: the plane goes to a face of the cell (crystals) or
    ! to the plane of the first three atoms (molecules), and the
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
       if (plane_def_form(w%pl%d,isys,ttshown,xediting,.true.,ipick,w%pl%picking)) &
          w%pl%viewdirty = .true.
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
       ! readout of the point under the mouse, the table of values (if
       ! open, it takes part of the space), and Close
       call igGetContentRegionAvail(sz0)
       sz0%y = sz0%y - iw_calcheight(3,0,.true.)
       htable = 0._c_float
       if (w%pl%showdata) then
          htable = max(0.35_c_float * sz0%y,iw_calcheight(4,0,.false.))
          sz0%y = sz0%y - htable - iw_calcheight(0,0,.true.)
       end if
       sz0%y = max(sz0%y,iw_calcheight(8,0,.false.))
       call draw_map(sz0,refit,hover,xhov,fhov)

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

    !> The labels of the atoms and critical points near the plane, and
    !> (in its own line) whether the map is shown in the view.
    subroutine draw_labels_option()
      use utils, only: iw_checkbox

      call iw_text("Labels",highlight=.true.,alignframe=.true.)
      ldum = iw_checkbox("Show##planeshowlabels",w%pl%showlabels,sameline=.true.)
      call iw_tooltip("Label the atoms and the critical points of the reference field near the plane",ttshown)
      if (w%pl%showlabels) then
         ldum = iw_dragfloat_real8("##planehmax",x1=w%pl%hmax,speed=0.01d0,min=0d0,decimal=3,sameline=.true.)
         call iw_tooltip("Largest distance to the plane of the labeled atoms and critical points (" //&
            lunit() // ")",ttshown)
         ldum = iw_checkbox("Show CPs##planeshowcps",w%pl%showcps,sameline=.true.)
         call iw_tooltip("Label also the critical points of the reference field other than the nuclei&
            & (n = non-nuclear maximum, b = bond, r = ring, c = cage)",ttshown)
      end if
      if (iw_checkbox("Show in View##planeshow",w%pl%showview)) w%pl%viewgen = -1
      call iw_tooltip("Draw this map in the view, following the changes made here, for as long as&
         & the window is open (Draw > Plane Map makes a map that stays)",ttshown)

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
    !> region. Returns in hover whether the mouse is on the region, the
    !> point under it (xhov, cell-frame Cartesian bohr), and the value
    !> there (fhov), which are also written below the map.
    subroutine draw_map(sz,refit,hover,xhov,fhov)
      use utils, only: iw_colormap_id
      type(ImVec2), intent(in) :: sz
      logical, intent(in) :: refit
      logical, intent(out) :: hover
      real*8, intent(out) :: xhov(3), fhov

      character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
      integer :: j, nx, ny
      integer(c_int) :: flags, coltext
      real(c_double) :: xmin, xmax, ymin, ymax, crn(2,4), mu, mv
      real(c_float) :: px(4), py(4), wbar
      type(ImVec2) :: pmin, pmax, szp, pt, tsz
      type(c_ptr) :: dl
      real*8 :: xyz(3)

      hover = .false.
      xhov = 0d0
      fhov = 0d0
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

         ! the contour lines
         if (w%pl%d%showcont) then
            if (w%pl%ncpos > 0) then
               str1 = "##planecpos" // c_null_char
               call ipSetNextLineStyle(ImVec4(w%pl%d%rgbpos(1),w%pl%d%rgbpos(2),w%pl%d%rgbpos(3),1._c_float),1._c_float)
               call ipPlotLine(c_loc(str1),c_loc(w%pl%cpos(1,1)),c_loc(w%pl%cpos(1,2)),int(w%pl%ncpos,c_int),&
                  ImPlotLineFlags_None,0_c_int)
            end if
            if (w%pl%ncneg > 0) then
               str1 = "##planecneg" // c_null_char
               call ipSetNextLineStyle(ImVec4(w%pl%d%rgbneg(1),w%pl%d%rgbneg(2),w%pl%d%rgbneg(3),1._c_float),1._c_float)
               call ipPlotLine(c_loc(str1),c_loc(w%pl%cneg(1,1)),c_loc(w%pl%cneg(1,2)),int(w%pl%ncneg,c_int),&
                  ImPlotLineFlags_None,0_c_int)
            end if
         end if

         ! the labels: a dot and the text
         if (w%pl%showlabels .and. w%pl%nmark > 0) then
            dl = igGetWindowDrawList()
            coltext = igGetColorU32_Col(ImGuiCol_Text,1._c_float)
            call ImDrawList_PushClipRect(dl,pmin,pmax,.true._c_bool)
            do j = 1, w%pl%nmark
               call ipPlotToPixels(w%pl%xmark(1,j),w%pl%xmark(2,j),ImAxis_X1,ImAxis_Y1,px(1),py(1))
               call ImDrawList_AddCircleFilled(dl,ImVec2(px(1),py(1)),3._c_float,coltext,8_c_int)
               str1 = trim(w%pl%lmark(j)) // c_null_char
               call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
               pt%x = px(1) + 4._c_float
               pt%y = py(1) - tsz%y - 1._c_float
               call ImDrawList_AddText_Vec2(dl,pt,coltext,c_loc(str1),c_null_ptr)
            end do
            call ImDrawList_PopClipRect(dl)
         end if

         ! the point under the mouse
         if (ipIsPlotHovered()) then
            call ipGetPlotMousePos(mu,mv)
            call plot_value(mu,mv,hover,xhov,fhov)
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

      ! the readout of the point under the mouse
      if (w%pl%showlabels .and. w%pl%marksmany) &
         call iw_text("(The region is too large for labels)",disabled=.true.)
      call iw_text("Cursor",highlight=.true.)
      if (hover) then
         if (ismol) then
            xyz = (xhov + sys(isys)%c%molx0) * bohrtoa
            call iw_text("(" // string(xyz(1),'f',decimal=4) // ", " // string(xyz(2),'f',decimal=4) //&
               ", " // string(xyz(3),'f',decimal=4) // ") Å",sameline=.true.)
         else
            xyz = sys(isys)%c%c2x(xhov)
            call iw_text("(" // string(xyz(1),'f',decimal=5) // ", " // string(xyz(2),'f',decimal=5) //&
               ", " // string(xyz(3),'f',decimal=5) // ") frac",sameline=.true.)
         end if
         call iw_text("value = " // string(fhov,'e',decimal=6),sameline=.true.)
      else
         call iw_text("move the mouse over the map to read the values there; it is also shown&
            & in the view",disabled=.true.,sameline=.true.)
      end if

    end subroutine draw_map

    !> The point of the region with plot coordinates (mu,mv): whether it
    !> is on the region (inside), its position (x, cell-frame Cartesian
    !> bohr), and the value there (f), interpolated from the four
    !> nearest points of the grid.
    subroutine plot_value(mu,mv,inside,x,f)
      real(c_double), intent(in) :: mu, mv
      logical, intent(out) :: inside
      real*8, intent(out) :: x(3), f

      integer :: nx, ny, i0, j0
      real*8 :: d(2), m(2,2), det, s, t, si, tj

      inside = .false.
      f = 0d0
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

      si = s * (nx-1)
      tj = t * (ny-1)
      i0 = min(int(si),nx-2) + 1
      j0 = min(int(tj),ny-2) + 1
      si = si - (i0-1)
      tj = tj - (j0-1)
      f = (1d0-si) * (1d0-tj) * w%pl%ff(i0,j0) + si * (1d0-tj) * w%pl%ff(i0+1,j0) +&
         (1d0-si) * tj * w%pl%ff(i0,j0+1) + si * tj * w%pl%ff(i0+1,j0+1)

    end subroutine plot_value

    !> The collapsible table of the values on the plane: the position
    !> of each point (fractional for crystals, Å for molecules), its
    !> plot coordinates, and the value. The menu of the header (right
    !> click) writes it to the output, the console, or a file.
    subroutine draw_data_section()
      use utils, only: iw_table_column, iw_table_headers_row
      integer, parameter :: itable_plane = 1 ! the table, for table_write_*
      character(kind=c_char,len=:), allocatable, target :: strsec, strtab
      character(len=:), allocatable :: un
      integer(c_int) :: flags
      integer :: k, ix, iy, iwrite, nx, ny
      real*8 :: xp(3), xq(3), pc(2)
      type(ImVec2) :: szt
      type(c_ptr) :: clipper
      type(ImGuiListClipper), pointer :: clipper_f

      w%pl%showdata = .false.
      if (.not.allocated(w%pl%ff)) return
      strsec = "Values on the plane##planedata" // c_null_char
      if (.not.igTreeNodeEx_Str(c_loc(strsec),ImGuiTreeNodeFlags_Framed)) return
      w%pl%showdata = .true.

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
      call igTreePop()

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
      real*8 :: oldx0(3), olduu(3), oldvv(3)
      integer :: oldnx, oldny

      refit = .false.
      fresh = .false.
      if (xediting .or. are_threads_running()) return
      fresh = .true.

      ! everything the values depend on, and the system
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
      refit = (oldnx == 0) .or. any(abs(oldx0 - w%pl%x0) > 1d-10) .or. &
         any(abs(olduu*(oldnx-1) - w%pl%uu*(w%pl%d%nx-1)) > 1d-10) .or.&
         any(abs(oldvv*(oldny-1) - w%pl%vv*(w%pl%d%ny-1)) > 1d-10)

      ! the values
      call w%pl%d%values(isys,w%pl%x0,w%pl%uu,w%pl%vv,w%pl%ff,w%pl%errmsg)

    end subroutine update_grid

    !> Calculate the contour lines again if the grid or the levels
    !> changed, in plot coordinates, the positive levels in cpos and
    !> the rest in cneg (NaN between two lines).
    subroutine update_contours()
      use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      integer*8, allocatable :: key(:)
      integer :: i, j, k, niso, nline, np
      integer, allocatable :: lptr(:), liso(:)
      real*8, allocatable :: ziso(:), xy(:,:)
      real*8 :: off(2), rval, nan
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
      if (.not.w%pl%d%showcont) return

      ! the lines, in the frame of the region; to plot coordinates
      call w%pl%d%contours(w%pl%ff,w%pl%uu,w%pl%vv,niso,ziso,nline,lptr,xy,liso)
      if (nline == 0) return
      off = to_plot(w%pl%x0)
      nan = ieee_value(nan,ieee_quiet_nan)
      if (allocated(w%pl%cpos)) deallocate(w%pl%cpos)
      if (allocated(w%pl%cneg)) deallocate(w%pl%cneg)
      allocate(w%pl%cpos(size(xy,2)+nline,2),w%pl%cneg(size(xy,2)+nline,2))
      do i = 1, nline
         pos = (ziso(liso(i)) > 0d0)
         np = merge(w%pl%ncpos,w%pl%ncneg,pos)
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
    !> plot coordinates and labels.
    subroutine update_marks()
      use tools_math, only: matinv, cross
      use types, only: realloc
      integer*8 :: key(6)
      integer :: k, l1, l2, l3, lo(3), hi(3), nc, ier, nx, ny
      real*8 :: am(3,3), nrm(3), xq(3), s(3), fc(3,4), hb, eps
      real*8, allocatable :: xf(:,:)
      character(len=12), allocatable :: lbl(:)

      key(1) = w%pl%gen
      key(2) = merge(1_8,0_8,w%pl%showlabels)
      key(6) = merge(1_8,0_8,w%pl%showcps)
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
      call label_candidates(isys,nc,xf,lbl,cps=w%pl%showcps,short=.true.)

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
  !> the region, the margin, the number of points, the quantity, the
  !> colors, and the contour lines. Returns whether anything changed.
  !> editing becomes true if a value is being dragged. If live, a drag
  !> is a change in every frame. If
  !> ipick is present, each point has a Pick button: ipick is the one
  !> pressed (0 = none), and pickoff disables them. Not live, a drag
  !> counts as a change only when it ends (the map of an object is
  !> calculated again on each change).
  module function plane_def_form(d,isys,ttshown,editing,live,ipick,pickoff) result(changed)
    use systems, only: sys
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_calcwidth,&
       iw_inputint, iw_combo_simple, iw_checkbox, iw_coloredit, iw_inputtext, iw_cmap_optstr,&
       iw_field_combo, iw_expression_helper, iw_exprctx_point
    use tools_io, only: string
    type(plane_def), intent(inout) :: d
    integer, intent(in) :: isys
    logical, intent(inout) :: ttshown
    logical, intent(inout) :: editing
    logical, intent(in) :: live
    integer, intent(out), optional :: ipick
    logical, intent(in), optional :: pickoff
    logical :: changed

    integer :: i, ival
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
    do i = 1, 3
       call draw_point(i)
    end do
    call draw_region()

    ! the quantity, the colors, and the contour lines
    call draw_quantity()
    call draw_colors()
    call draw_contours()

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

    !> The region of the plane (a rectangle or a parallelogram), the
    !> margin around it, and the number of points along each side.
    subroutine draw_region()

      call iw_text("Region",highlight=.true.,alignframe=.true.)
      call igSameLine(0._c_float,-1._c_float)
      call igPushItemWidth(iw_calcwidth(14,1))
      ival = d%region
      call iw_combo_simple("##planeregion","Rectangle" // c_null_char // "Parallelogram" // c_null_char,&
         d%region)
      changed = changed .or. (ival /= d%region)
      call igPopItemWidth()
      call iw_tooltip("Rectangle: the first side goes from the origin to the X end, and the&
         & second is perpendicular to it, as long as needed to reach the Y point. Parallelogram:&
         & the sides go from the origin to the X end and to the Y point, as in PLANE",ttshown)

      call iw_text("Margin",highlight=.true.,alignframe=.true.,sameline=.true.)
      call drag("##planemargin",x1=d%margin,speed=0.01d0,decimal=3,sameline=.true.)
      call iw_tooltip("Length added to the region on every side (" // un //&
         "; negative removes it)",ttshown)

      call iw_text("Points",highlight=.true.,alignframe=.true.,sameline=.true.)
      ch = iw_inputint("##planenx",d%nx,width=6,sameline=.true.,notlive=.true.)
      d%nx = max(min(d%nx,plane_npts_max),plane_npts_min)
      changed = changed .or. ch
      call iw_tooltip("Number of points along the first side (" // string(plane_npts_min) // " to " //&
         string(plane_npts_max) // ")",ttshown)
      ch = iw_inputint("##planeny",d%ny,width=6,sameline=.true.,notlive=.true.)
      d%ny = max(min(d%ny,plane_npts_max),plane_npts_min)
      changed = changed .or. ch
      call iw_tooltip("Number of points along the second side (" // string(plane_npts_min) // " to " //&
         string(plane_npts_max) // ")",ttshown)

    end subroutine draw_region

    !> The quantity on the map: a field and one of its properties, or
    !> an arithmetic expression (checked as it is typed).
    subroutine draw_quantity()
      integer :: ikind, ifield, iq
      logical :: echanged
      character(len=:), allocatable :: lerr

      call iw_text("Quantity",highlight=.true.,alignframe=.true.)
      ikind = merge(1,0,d%isexpr)
      call igSameLine(0._c_float,-1._c_float)
      call igPushItemWidth(iw_calcwidth(10,1))
      call iw_combo_simple("##planekind","Field" // c_null_char // "Expression" // c_null_char,ikind)
      call igPopItemWidth()
      call iw_tooltip("A quantity of a field, or an arithmetic expression of the fields",ttshown)
      changed = changed .or. (d%isexpr .neqv. (ikind == 1))
      d%isexpr = (ikind == 1)
      call igSameLine(0._c_float,-1._c_float)
      if (.not.d%isexpr) then
         ifield = d%fref%id
         if (iw_field_combo("##planefield",isys,ifield,width=iw_calcwidth(24,1),&
            nonestr="<field not available>")) then
            call d%fref%set(sys(isys),ifield)
            changed = .true.
         end if
         call iw_tooltip("Field on the map",ttshown)
         call igSameLine(0._c_float,-1._c_float)
         iq = d%nti
         call iw_combo_simple("##planequant",quant_combo,iq)
         changed = changed .or. (iq /= d%nti)
         d%nti = max(min(iq,nquant-1),0)
         call iw_tooltip("Quantity of the field on the map",ttshown)
      else
         echanged = iw_inputtext("##planeexpr",bufsize=1023,textf=d%expr,width=36)
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

    !> The colors of the map.
    subroutine draw_colors()

      call iw_text("Colors",highlight=.true.,alignframe=.true.)
      call igSameLine(0._c_float,-1._c_float)
      call igPushItemWidth(iw_calcwidth(10,1))
      ival = d%icmap
      call iw_combo_simple("##planecmap",iw_cmap_optstr,d%icmap,startsatone=.true.)
      changed = changed .or. (ival /= d%icmap)
      call igPopItemWidth()
      call iw_tooltip("Colormap of the map",ttshown)
      call igSameLine(0._c_float,-1._c_float)
      call igPushItemWidth(iw_calcwidth(10,1))
      ival = d%cscale
      call iw_combo_simple("##planecscale","Linear" // c_null_char // "Log |f|" // c_null_char //&
         "Arctangent" // c_null_char,d%cscale)
      changed = changed .or. (ival /= d%cscale)
      call igPopItemWidth()
      call iw_tooltip("How the values become colors: linearly, by the logarithm of their&
         & absolute value, or by their arctangent (as in PLANE COLORMAP)",ttshown)
      changed = changed .or. iw_checkbox("Range##planecrange",d%crange,sameline=.true.)
      call iw_tooltip("Colors for a fixed range of values (or the range of the values on the map)",ttshown)
      if (d%crange) then
         call drag("##planecmin",x1=d%cmin,speed=0.01d0,decimal=4,sameline=.true.)
         call iw_tooltip("The value of the first color",ttshown)
         call drag("##planecmax",x1=d%cmax,speed=0.01d0,decimal=4,sameline=.true.)
         call iw_tooltip("The value of the last color",ttshown)
      end if

    end subroutine draw_colors

    !> The contour lines of the map.
    subroutine draw_contours()
      integer :: ic

      call iw_text("Contours",highlight=.true.,alignframe=.true.)
      changed = changed .or. iw_checkbox("Show##planeshowcont",d%showcont,sameline=.true.)
      call iw_tooltip("Draw contour lines on the map",ttshown)
      if (.not.d%showcont) return
      ic = 0
      do while (ic < 4 .and. ctype_of_combo(ic) /= d%ctype)
         ic = ic + 1
      end do
      call igSameLine(0._c_float,-1._c_float)
      call igPushItemWidth(iw_calcwidth(12,1))
      call iw_combo_simple("##planectype","Linear" // c_null_char // "Logarithmic" // c_null_char //&
         "Arctangent" // c_null_char // "Bader" // c_null_char // "List" // c_null_char,ic)
      call igPopItemWidth()
      changed = changed .or. (d%ctype /= ctype_of_combo(ic))
      d%ctype = ctype_of_combo(ic)
      call iw_tooltip("The contour levels, as in PLANE CONTOUR: evenly spaced, logarithmic&
         & (positive and negative), evenly spaced in the arctangent, the Bader set (1, 2, 4,&
         & 8 x 10^n), or a list of values",ttshown)
      if (d%ctype == niso_manual) then
         changed = changed .or. iw_inputtext("##planellist",bufsize=1023,textf=d%llist,width=30,sameline=.true.)
         call iw_tooltip("The contour levels, separated by spaces or commas",ttshown)
      elseif (d%ctype /= niso_bader) then
         ch = iw_inputint("##planenlev",d%nlev,width=5,sameline=.true.,notlive=.true.)
         d%nlev = max(min(d%nlev,plane_nlev_max),2_c_int)
         changed = changed .or. ch
         call iw_tooltip("Number of contour levels (at least 2)",ttshown)
         changed = changed .or. iw_checkbox("Range##planelrange",d%lrange,sameline=.true.)
         call iw_tooltip("Contour levels in a fixed range of values (or the range of the&
            & values on the map)",ttshown)
         if (d%lrange) then
            call drag("##planelmin",x1=d%lmin,speed=0.01d0,decimal=4,sameline=.true.)
            call iw_tooltip("The lowest contour level",ttshown)
            call drag("##planelmax",x1=d%lmax,speed=0.01d0,decimal=4,sameline=.true.)
            call iw_tooltip("The highest contour level",ttshown)
         end if
      end if
      changed = changed .or. iw_coloredit("##planergbpos",rgb=d%rgbpos,sameline=.true.)
      call iw_tooltip("Color of the contour lines of the positive levels",ttshown)
      changed = changed .or. iw_coloredit("##planergbneg",rgb=d%rgbneg,sameline=.true.)
      call iw_tooltip("Color of the contour lines of the zero and negative levels",ttshown)

    end subroutine draw_contours
  end function plane_def_form

  !> The editor of a plane map object (reptype_slice): the map, its
  !> opacity, the radius of the contour lines, and what is wrong with
  !> it. Returns whether anything changed.
  module function draw_editrep_slice(w,ttshown) result(changed)
    use utils, only: iw_text, iw_tooltip, iw_dragfloat_real8, iw_dragfloat_realc
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: editing

    editing = .false.
    changed = plane_def_form(w%rep%slice%d,w%rep%id,ttshown,editing,.false.)

    ! how it is drawn
    call iw_text("Opacity",highlight=.true.,alignframe=.true.)
    changed = changed .or. iw_dragfloat_realc("##sliceopacity",x1=w%rep%slice%alpha,speed=0.01_c_float,&
       min=0._c_float,max=1._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
    call iw_tooltip("Opacity of the map (1 = opaque)",ttshown)
    call iw_text("Line radius",highlight=.true.,alignframe=.true.,sameline=.true.)
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
