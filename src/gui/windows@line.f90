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

! Routines for the line profile window (the LINE keyword).
submodule (windows) line
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! transient markers in the view: the endpoints (A, B), the segment,
  ! and the point of the cursor; opacity, and radius of the endpoints
  ! (bohr), a fraction of the scene radius between the two limits
  real(c_float), parameter :: enda_rgb(3) = (/0.95_c_float,0.35_c_float,0.20_c_float/)
  real(c_float), parameter :: endb_rgb(3) = (/0.20_c_float,0.55_c_float,0.95_c_float/)
  real(c_float), parameter :: seg_rgb(3) = (/0.85_c_float,0.85_c_float,0.85_c_float/)
  real(c_float), parameter :: cur_rgb(3) = (/1.00_c_float,0.80_c_float,0.10_c_float/)
  real(c_float), parameter :: marker_alpha = 0.8_c_float
  real*8, parameter :: marker_radmin = 0.20d0
  real*8, parameter :: marker_radmax = 0.50d0

  ! the colors of the curves, in the order they are added (Tableau 10)
  integer, parameter :: npalette = 8
  real(c_float), parameter :: palette(3,npalette) = reshape((/&
     0.122_c_float,0.467_c_float,0.706_c_float, 1.000_c_float,0.498_c_float,0.055_c_float,&
     0.173_c_float,0.627_c_float,0.173_c_float, 0.839_c_float,0.153_c_float,0.157_c_float,&
     0.580_c_float,0.404_c_float,0.741_c_float, 0.549_c_float,0.337_c_float,0.294_c_float,&
     0.890_c_float,0.467_c_float,0.761_c_float, 0.090_c_float,0.745_c_float,0.812_c_float/),&
     (/3,npalette/))

  ! an atom or a critical point is on the line if it is closer than
  ! this (bohr)
  real*8, parameter :: online_tol = 0.05d0

  ! the number of points
  integer(c_int), parameter :: npts_min = 2
  integer(c_int), parameter :: npts_max = 10000

contains

  !> Draw the line profile window: the values of fields of the system
  !> shown by the anchor view, their derivatives, or arithmetic
  !> expressions, along the segment between two points (typed or
  !> picked in the view), as LINE calculates them, plotted as curves.
  module subroutine draw_line(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8,&
       iw_calcwidth, iw_close_event, iw_setpos_bottomright, iw_calcheight,&
       iw_inputint, iw_combo_simple
    use representations, only: rep_shape, shapekind_sphere, shapekind_cylinder, bondrad_def
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged, ismol, ldum, xediting, refit
    integer :: isys, iview, istat
    real*8 :: xc(3), xa(3), xb(3), rad, dmax, v(3,3)
    type(ImVec2) :: sz0
    real(c_float) :: htable

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)
    ismol = .false.
    if (goodsys) ismol = sys(isys)%c%ismolecule

    ! initialize the state on the first pass
    if (w%firstpass) then
       w%ln = line_state()
       w%errmsg = ""
    end if

    ! a new system: the endpoints go to the first two atoms (or the
    ! diagonal of the cell) and there is one curve, the value of the
    ! reference field
    if (goodsys) then
       if (.not.w%ln%sys%is(isys) .or. syschanged) then
          call w%ln%sys%set(isys)
          call cancel_pick()
          w%errmsg = ""
          call default_endpoints()
          if (allocated(w%ln%curve)) deallocate(w%ln%curve)
          allocate(w%ln%curve(4))
          w%ln%ncurve = 0
          w%ln%ncolor = 0
          call add_curve()
          w%ln%dcur = -1d0 ! the middle of the line, when it is calculated
          w%ln%keysys = 0
       end if
       ! the reference for the curves whose field is gone
       call fix_fields()
    end if

    ! a pending pick in the view it was armed on; dropped, releasing the
    ! view, if the system is not ready or the window moved to another view
    if (w%ln%picking) then
       if (goodsys .and. iview == w%ln%pickview) then
          call view_pick_result(iview,w%id,isys,w%ln%pick,istat,xc)
          if (istat == ipick_point) then
             if (w%ln%pickend == 1) then
                w%ln%xa = form_coords(isys,xc)
             else
                w%ln%xb = form_coords(isys,xc)
             end if
          end if
          if (istat /= ipick_pending) w%ln%picking = .false.
       else
          call cancel_pick()
       end if
    end if

    if (goodsys) then
       ! the system
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! the endpoints and the number of points
       xediting = .false.
       call draw_endpoint(1,w%ln%xa)
       call draw_endpoint(2,w%ln%xb)
       if (iw_button("Swap##lineswap")) then
          xc = w%ln%xa
          w%ln%xa = w%ln%xb
          w%ln%xb = xc
       end if
       call iw_tooltip("Exchange the two endpoints of the line",ttshown)
       call iw_text("Points",highlight=.true.,alignframe=.true.,sameline=.true.)
       ldum = iw_inputint("##linenpts",w%ln%npts,width=7,sameline=.true.,notlive=.true.)
       w%ln%npts = max(min(w%ln%npts,npts_max),npts_min)
       call iw_tooltip("Number of points on the line, including the endpoints (" //&
          string(npts_min) // " to " // string(npts_max) // ")",ttshown)

       ! the curves
       call draw_curves_section()
       if (sys(isys)%c%ncel == 0) &
          call iw_text("This system has no atoms: there are no field values to plot",danger=.true.,wrap=.true.)

       ! the line and the values, calculated again if anything they
       ! depend on changed
       call update_line(refit)

       ! the plot takes the rest of the window, leaving room for the
       ! scale of the value axis, the cursor, the table of values (if
       ! open, it takes part of the space), and Close
       call igGetContentRegionAvail(sz0)
       sz0%x = -1._c_float
       sz0%y = sz0%y - iw_calcheight(5,0,.true.)
       htable = 0._c_float
       if (w%ln%showdata) then
          htable = max(0.4_c_float * sz0%y,iw_calcheight(4,0,.false.))
          sz0%y = sz0%y - htable - iw_calcheight(0,0,.true.)
       end if
       sz0%y = max(sz0%y,iw_calcheight(6,0,.false.))
       call draw_plot(sz0,refit)

       ! the table of values
       call draw_data_section()

       ! the endpoints, the segment, and the cursor in the anchor view,
       ! for as long as the window is open. The cursor is thicker than
       ! a bond, so it shows on a line between two bonded atoms
       if (associated(win(iview)%sc) .and. allocated(w%ln%d)) then
          if (win(iview)%sc%isinit /= 0) then
             rad = min(max(0.01d0 * real(win(iview)%sc%scenerad,8),marker_radmin),marker_radmax)
             xa = to_abs(w%ln%xa)
             xb = to_abs(w%ln%xb)
             v = 0d0
             v(:,1) = xb - xa
             dmax = w%ln%d(size(w%ln%d))
             xc = xa
             if (dmax > 0d0) xc = xa + (xb - xa) * w%ln%dcur / dmax
             call win(iview)%sc%show_transient_shapes(w%uid,1,(/&
                rep_shape(kind=shapekind_cylinder,x1=xa,v=v,rad=0.25d0*rad,rgb=seg_rgb,alpha=marker_alpha),&
                rep_shape(kind=shapekind_sphere,x1=xa,rad=rad,rgb=enda_rgb,alpha=marker_alpha),&
                rep_shape(kind=shapekind_sphere,x1=xb,rad=rad,rgb=endb_rgb,alpha=marker_alpha),&
                rep_shape(kind=shapekind_sphere,x1=xc,rad=max(rad,1.6d0*bondrad_def),rgb=cur_rgb,&
                alpha=1._c_float)/))
          end if
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

    !> One endpoint of the line (ie = 1 for A, 2 for B): the label, the
    !> button to pick it in the view, and its coordinates.
    subroutine draw_endpoint(ie,x)
      integer, intent(in) :: ie
      real*8, intent(inout) :: x(3)

      character(len=1) :: lbl
      real(c_float) :: rgb(3)

      lbl = merge("A","B",ie == 1)
      rgb = merge(enda_rgb,endb_rgb,ie == 1)
      call iw_text(lbl,rgb=rgb,alignframe=.true.)
      call iw_tooltip("Endpoint " // lbl // " of the line (shown in the view in this color)",ttshown)
      if (iw_button("Pick##linepick" // lbl,sameline=.true.,disabled=w%ln%picking)) then
         w%ln%picking = .true.
         w%ln%pickend = ie
         w%ln%pickview = iview
         call w%ln%pick%arm()
         call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick endpoint " // lbl // ": an atom,&
            & a critical point, or a point in space",w%id,acceptempty=.true.)
      end if
      call iw_tooltip("Pick endpoint " // lbl // " in the view: an atom, a critical point,&
         & or a point in empty space",ttshown)
      if (ismol) then
         ldum = iw_dragfloat_real8("(Å)##linex" // lbl,x3=x,speed=0.01d0,decimal=5,sameline=.true.)
      else
         ldum = iw_dragfloat_real8("(frac)##linex" // lbl,x3=x,speed=0.001d0,decimal=6,sameline=.true.)
      end if
      ! the calculation waits until the coordinates are no longer edited
      xediting = xediting .or. igIsItemActive()
      call iw_tooltip("Coordinates of endpoint " // lbl,ttshown)

    end subroutine draw_endpoint

    !> The table of the curves: whether each is shown, its color, and
    !> what it is (a field and a quantity of it, or an expression);
    !> the button that adds curves, the icons that remove them, and the
    !> errors.
    subroutine draw_curves_section()
      use gui_main, only: g
      use utils, only: iw_checkbox, iw_coloredit, iw_field_combo, iw_inputtext,&
         iw_expression_helper, iw_exprctx_point, iw_table_column, iw_close_button
      integer :: i, idel, ifield, ikind, iq
      logical :: echanged
      integer(c_int) :: flags
      character(kind=c_char,len=:), allocatable, target :: strtab
      character(len=:), allocatable :: lerr, sfx
      type(ImVec2) :: szt

      call iw_text("Curves",highlight=.true.)
      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      strtab = "##linecurves" // c_null_char
      szt%x = 0._c_float
      szt%y = 0._c_float
      idel = 0
      if (igBeginTable(c_loc(strtab),5,flags,szt,0._c_float)) then
         ! the width of the close icon, also when it is not there (one
         ! curve), so the table does not move when it appears
         call iw_table_column("",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed,&
            width=igGetTextLineHeight())
         call iw_table_column("Show",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Color",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Kind",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Definition",id=4_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call igTableHeadersRow()
         do i = 1, w%ln%ncurve
            associate(cv => w%ln%curve(i))
              sfx = "##linecurve" // string(i)
              call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)

              ! remove the curve (not the last one); the close icon is
              ! one line of text high: centered on the row frames
              if (w%ln%ncurve > 1 .and. igTableSetColumnIndex(0_c_int)) then
                 call igSetCursorPosY(igGetCursorPosY() + g%Style%FramePadding%y)
                 if (iw_close_button(sfx // "del")) idel = i
                 call iw_tooltip("Remove this curve",ttshown)
              end if

              ! shown, and the color
              if (igTableSetColumnIndex(1_c_int)) then
                 ldum = iw_checkbox(sfx // "show",cv%show)
                 call iw_tooltip("Show this curve in the plot",ttshown)
              end if
              if (igTableSetColumnIndex(2_c_int)) then
                 ldum = iw_coloredit(sfx // "rgb",rgb=cv%rgb,nolabel=.true.)
                 call iw_tooltip("Color of the curve",ttshown)
              end if

              ! a field or an expression
              if (igTableSetColumnIndex(3_c_int)) then
                 ikind = merge(1,0,cv%isexpr)
                 call igPushItemWidth(iw_calcwidth(10,1))
                 call iw_combo_simple(sfx // "kind","Field" // c_null_char // "Expression" // c_null_char,&
                    ikind)
                 call igPopItemWidth()
                 call iw_tooltip("A quantity of a field, or an arithmetic expression of the fields",ttshown)
                 cv%isexpr = (ikind == 1)
              end if

              ! the definition: the field and the quantity, or the
              ! expression (checked as it is typed)
              if (igTableSetColumnIndex(4_c_int)) then
                 if (.not.cv%isexpr) then
                    ifield = cv%fref%id
                    if (iw_field_combo(sfx // "field",isys,ifield,width=iw_calcwidth(24,1),&
                       nonestr="<field not available>")) call cv%fref%set(sys(isys),ifield)
                    call iw_tooltip("Field of the curve",ttshown)
                    call igSameLine(0._c_float,-1._c_float)
                    iq = cv%nti
                    call iw_combo_simple(sfx // "quant",quant_combo,iq)
                    cv%nti = max(min(iq,nquant-1),0)
                    call iw_tooltip("Quantity of the field plotted along the line",ttshown)
                 else
                    echanged = iw_inputtext(sfx // "expr",bufsize=1023,textf=cv%expr,width=36)
                    call iw_tooltip("Arithmetic expression of the fields of the system (for instance,&
                       & ""$1 - $2"" or ""$1:l""), evaluated at the points of the line",ttshown)
                    if (iw_expression_helper(sfx // "help",isys,iw_exprctx_point,textf=cv%expr,&
                       ttshown=ttshown,ifield=sys(isys)%iref)) echanged = .true.
                    if (echanged) then
                       cv%exprerr = ""
                       if (len_trim(cv%expr) > 0) then
                          call sys(isys)%check_expression(trim(cv%expr),lerr)
                          cv%exprerr = lerr
                       end if
                    end if
                 end if
              end if

            end associate
         end do
         call igEndTable()
      end if
      if (idel > 0) call remove_curve(idel)

      ! add a curve
      if (iw_button("Add Curve##lineadd")) call add_curve()
      call iw_tooltip("Add a curve: the value of the reference field, which can then be changed",ttshown)

      ! what is wrong with the curves
      do i = 1, w%ln%ncurve
         associate(cv => w%ln%curve(i))
           if (cv%isexpr) then
              if (len_trim(cv%expr) == 0) then
                 call iw_text("Curve " // string(i) // ": the expression is empty",danger=.true.,wrap=.true.)
              elseif (len_trim(cv%exprerr) > 0) then
                 call iw_text("Curve " // string(i) // ": " // trim(cv%exprerr),danger=.true.,wrap=.true.)
              end if
           end if
           if (allocated(cv%errmsg)) then
              if (len_trim(cv%errmsg) > 0) &
                 call iw_text("Curve " // string(i) // ": " // cv%errmsg,danger=.true.,wrap=.true.)
           end if
         end associate
      end do

    end subroutine draw_curves_section

    !> The plot of the curves (size sz), with the atoms and critical
    !> points on the line, the cursor, and the values of the curves at
    !> the cursor below it. If refit, fit the axes to the curves.
    subroutine draw_plot(sz,refit)
      type(ImVec2), intent(in) :: sz
      logical, intent(inout) :: refit

      character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
      integer :: i, j, np, ip
      integer(c_int) :: flags, yscale
      integer :: iscale
      type(ImVec2) :: tsz
      real(c_double) :: xmin, xmax, ymin, ymax, ylo, yhi, dmax
      real(c_float) :: px, py, pxa, pya, pxb, pyb
      type(ImVec2) :: pmin, pmax, pt
      type(c_ptr) :: dl
      integer(c_int) :: coldim
      real*8 :: t, val

      if (.not.allocated(w%ln%d)) return
      np = size(w%ln%d)
      dmax = w%ln%d(np)

      ! the value axis: linear, symmetric log, or log
      call iw_text("Value axis",highlight=.true.,alignframe=.true.)
      call igSameLine(0._c_float,-1._c_float)
      call igPushItemWidth(iw_calcwidth(16,1))
      iscale = w%ln%yscale
      call iw_combo_simple("##lineyscale","Linear" // c_null_char // "Symmetric log" // c_null_char //&
         "Log" // c_null_char,w%ln%yscale)
      call igPopItemWidth()
      if (w%ln%yscale /= iscale) refit = .true.
      call iw_tooltip("Scale of the value axis. The symmetric log scale is logarithmic away from&
         & zero and linear near it, for curves that change sign",ttshown)
      if (w%ln%yscale == 1) then
         yscale = ImPlotScale_SymLog
      elseif (w%ln%yscale == 2) then
         yscale = ImPlotScale_Log10
      else
         yscale = ImPlotScale_Linear
      end if

      str1 = "##lineplot" // c_null_char
      flags = ImPlotFlags_NoTitle
      if (ipBeginPlot(c_loc(str1),sz,flags)) then
         if (ismol) then
            str2 = "Distance from A (Å)" // c_null_char
         else
            str2 = "Distance from A (bohr)" // c_null_char
         end if
         str3 = "" // c_null_char
         call ipSetupAxes(c_loc(str2),c_loc(str3),ImPlotAxisFlags_None,ImPlotAxisFlags_None)
         call ipSetupAxisScale(ImAxis_Y1,yscale)

         ! fit the axes when the curves change; otherwise, keep the
         ! user's zoom
         if (refit) then
            call ipSetupAxisLimits(ImAxis_X1,0._c_double,max(dmax,1d-10),ImPlotCond_Always)
            if (curves_range(ylo,yhi)) &
               call ipSetupAxisLimits(ImAxis_Y1,ylo,yhi,ImPlotCond_Always)
         end if

         ! the curves
         do i = 1, w%ln%ncurve
            associate(cv => w%ln%curve(i))
              if (.not.cv%show .or. .not.allocated(cv%y)) cycle
              if (size(cv%y) /= np) cycle
              str1 = curve_label(i) // "##linecurve" // string(cv%uid) // c_null_char
              call ipSetNextLineStyle(ImVec4(cv%rgb(1),cv%rgb(2),cv%rgb(3),1._c_float),2._c_float)
              call ipPlotLine(c_loc(str1),c_loc(w%ln%d),c_loc(cv%y),int(np,c_int),&
                 ImPlotLineFlags_SkipNaN,0_c_int)
            end associate
         end do

         ! the atoms and critical points on the line: a vertical line
         ! and the label at the bottom of the plot (away from the
         ! legend), on the left of the line near the right edge
         call ipGetPlotCurrentLimits(xmin,xmax,ymin,ymax)
         call ipPlotToPixels(xmin,ymin,ImAxis_X1,ImAxis_Y1,pxa,pya)
         call ipPlotToPixels(xmax,ymax,ImAxis_X1,ImAxis_Y1,pxb,pyb)
         pmin%x = min(pxa,pxb)
         pmin%y = min(pya,pyb)
         pmax%x = max(pxa,pxb)
         pmax%y = max(pya,pyb)
         if (w%ln%nmark > 0 .and. pmax%x > pmin%x .and. pmax%y > pmin%y) then
            dl = igGetWindowDrawList()
            coldim = igGetColorU32_Col(ImGuiCol_Text,0.5_c_float)
            call ImDrawList_PushClipRect(dl,pmin,pmax,.true._c_bool)
            do j = 1, w%ln%nmark
               call ipPlotToPixels(real(w%ln%dmark(j),c_double),ymax,ImAxis_X1,ImAxis_Y1,px,py)
               call ImDrawList_AddLine(dl,ImVec2(px,pmin%y),ImVec2(px,pmax%y),coldim,1._c_float)
               str1 = trim(w%ln%lmark(j)) // c_null_char
               call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
               pt%x = px + 3._c_float
               if (pt%x + tsz%x > pmax%x) pt%x = px - 3._c_float - tsz%x
               pt%y = pmax%y - tsz%y - 2._c_float
               call ImDrawList_AddText_Vec2(dl,pt,coldim,c_loc(str1),c_null_ptr)
            end do
            call ImDrawList_PopClipRect(dl)
         end if

         ! the cursor, which moves the point in the view
         ldum = logical(ipDragLineX(0_c_int,w%ln%dcur,ImVec4(cur_rgb(1),cur_rgb(2),cur_rgb(3),&
            1._c_float),1.5_c_float,ImPlotDragToolFlags_None))
         w%ln%dcur = max(min(w%ln%dcur,dmax),0d0)
         call ipEndPlot()
      end if

      ! the values of the curves at the cursor, interpolated between
      ! the two nearest points
      call iw_text("Cursor",highlight=.true.)
      call iw_text("d = " // string(w%ln%dcur,'f',decimal=4),sameline=.true.)
      call iw_tooltip("Drag the vertical line in the plot to move the cursor, which is also&
         & shown in the view",ttshown)
      if (np >= 2 .and. dmax > 0d0) then
         t = max(min(w%ln%dcur / dmax,1d0),0d0) * real(np-1,8)
         ip = min(max(floor(t),0),np-2) + 1
         t = t - real(ip-1,8)
         do i = 1, w%ln%ncurve
            associate(cv => w%ln%curve(i))
              if (.not.cv%show .or. .not.allocated(cv%y)) cycle
              if (size(cv%y) /= np) cycle
              val = (1d0 - t) * cv%y(ip) + t * cv%y(ip+1)
              call iw_text("  " // curve_label(i) // " = " // string(val,'e',decimal=6),&
                 rgb=cv%rgb,sameline=.true.)
            end associate
         end do
      end if

    end subroutine draw_plot

    !> The collapsible table of the values along the line: the
    !> position and the distance of each point, and the value of
    !> each curve. The menu of the header (right click) chooses the
    !> curves in the table and writes it to the output, the console,
    !> or a file.
    subroutine draw_data_section()
      use utils, only: iw_table_column, iw_table_headers_row, iw_tblwrite_none
      integer, parameter :: itable_line = 1 ! the table, for table_write_*
      character(kind=c_char,len=:), allocatable, target :: strsec, strtab
      integer(c_int) :: flags
      integer :: i, k, iwrite, ncol, np, itoggle
      logical :: rclick
      real*8 :: xp(3), xa(3), xb(3)
      type(ImVec2) :: szt
      type(c_ptr) :: clipper
      type(ImGuiListClipper), pointer :: clipper_f

      w%ln%showdata = .false.
      if (.not.allocated(w%ln%d)) return
      strsec = "Values along the line##linedata" // c_null_char
      if (.not.igTreeNodeEx_Str(c_loc(strsec),ImGuiTreeNodeFlags_Framed)) return
      w%ln%showdata = .true.

      np = size(w%ln%d)
      ncol = 4
      do i = 1, w%ln%ncurve
         if (in_table(i)) ncol = ncol + 1
      end do
      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      flags = ior(flags,ImGuiTableFlags_ScrollY)
      flags = ior(flags,ImGuiTableFlags_ScrollX)
      strtab = "##linedatatable" // c_null_char
      szt%x = 0._c_float
      szt%y = htable
      if (igBeginTable(c_loc(strtab),int(ncol,c_int),flags,szt,0._c_float)) then
         if (ismol) then
            call iw_table_column("x (Å)",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("y (Å)",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("z (Å)",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("d (Å)",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         else
            call iw_table_column("x (frac)",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("y (frac)",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("z (frac)",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("d (bohr)",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         end if
         k = 3
         do i = 1, w%ln%ncurve
            if (.not.in_table(i)) cycle
            k = k + 1
            call iw_table_column(curve_label(i),id=int(k,c_int),flags=ImGuiTableColumnFlags_WidthFixed)
         end do
         call iw_table_headers_row(freezetop=.true.,rclicked=rclick)
         iwrite = iw_tblwrite_none
         call header_menu(rclick,iwrite,itoggle)

         ! the rows; if the table is being written as text, they are
         ! captured as they are drawn (all of them)
         ! the endpoints the values were calculated for
         xa = w%ln%keyxa
         xb = w%ln%keyxb
         clipper = ImGuiListClipper_ImGuiListClipper()
         call ImGuiListClipper_Begin(clipper,int(np,c_int),-1._c_float)
         call w%table_write_begin(itable_line,iwrite,clipper)
         do while(ImGuiListClipper_Step(clipper))
            call c_f_pointer(clipper,clipper_f)
            do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
               call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
               xp = xa + (xb - xa) * real(i-1,8) / real(np-1,8)
               if (ismol) then
                  call cell(0,string(xp(1),'f',decimal=6))
                  call cell(1,string(xp(2),'f',decimal=6))
                  call cell(2,string(xp(3),'f',decimal=6))
               else
                  call cell(0,string(xp(1),'f',decimal=7))
                  call cell(1,string(xp(2),'f',decimal=7))
                  call cell(2,string(xp(3),'f',decimal=7))
               end if
               call cell(3,string(w%ln%d(i),'f',decimal=6))
               call curve_cells(i)
            end do
         end do
         call ImGuiListClipper_End(clipper)
         call ImGuiListClipper_destroy(clipper)
         call w%table_write_end(itable_line,table_title("Line profile",isys))
         call igEndTable()
         ! a curve toggled in the menu: its column changes in the next
         ! frame, so the columns of this one match its cells
         if (itoggle > 0) w%ln%curve(itoggle)%intable = .not.w%ln%curve(itoggle)%intable
      end if
      call igTreePop()

    end subroutine draw_data_section

    !> The menu of the header of the table of values, opened by a
    !> right click (rclick): the curves in the table (itoggle = the
    !> curve clicked, 0 if none) and, last, the entries that ask (in
    !> iwrite, iw_tblwrite_*) for the table to be written as text.
    subroutine header_menu(rclick,iwrite,itoggle)
      use utils, only: iw_menuitem, iw_table_write_items
      logical, intent(in) :: rclick
      integer, intent(inout) :: iwrite
      integer, intent(out) :: itoggle

      integer :: i
      character(kind=c_char,len=:), allocatable, target :: strpop

      itoggle = 0
      strpop = "##linecolumnmenu" // c_null_char
      if (rclick) call igOpenPopup_Str(c_loc(strpop),ImGuiPopupFlags_None)
      if (igBeginPopup(c_loc(strpop),ImGuiWindowFlags_None)) then
         call igPushItemFlag(ImGuiItemFlags_SelectableDontClosePopup,.true._c_bool)
         do i = 1, w%ln%ncurve
            if (.not.has_values(i)) cycle
            if (iw_menuitem(curve_label(i) // "##lineintable" // string(w%ln%curve(i)%uid),&
               selected=w%ln%curve(i)%intable)) itoggle = i
            call iw_tooltip("Show the values of this curve in the table",ttshown)
         end do
         call igPopItemFlag()
         call igSeparator()
         call iw_table_write_items(iwrite,ttshown)
         call igEndPopup()
      end if

    end subroutine header_menu

    !> The values of the curves at point ip, in the columns of the
    !> table of values after the position and the distance.
    subroutine curve_cells(ip)
      integer, intent(in) :: ip
      integer :: i, k

      k = 3
      do i = 1, w%ln%ncurve
         if (.not.in_table(i)) cycle
         k = k + 1
         call cell(k,string(w%ln%curve(i)%y(ip),'e',decimal=8))
      end do

    end subroutine curve_cells

    !> Text str in column icol of the current row of a table.
    subroutine cell(icol,str)
      integer, intent(in) :: icol
      character(len=*), intent(in) :: str

      if (igTableSetColumnIndex(int(icol,c_int))) call iw_text(str)

    end subroutine cell

    !> Curve i has values for the current line.
    logical function has_values(i)
      integer, intent(in) :: i

      has_values = .false.
      if (.not.allocated(w%ln%d)) return
      if (.not.allocated(w%ln%curve(i)%y)) return
      has_values = (size(w%ln%curve(i)%y) == size(w%ln%d))

    end function has_values

    !> Curve i is a column of the table of values.
    logical function in_table(i)
      integer, intent(in) :: i

      in_table = has_values(i) .and. w%ln%curve(i)%intable

    end function in_table

    !> The label of curve i: $n with the quantity (:lap), or the
    !> expression.
    function curve_label(i) result(str)
      integer, intent(in) :: i
      character(len=:), allocatable :: str

      associate(cv => w%ln%curve(i))
        if (cv%isexpr) then
           str = trim(cv%expr)
        else
           str = "$" // string(cv%fref%id) // trim(quant_short(cv%nti))
        end if
      end associate

    end function curve_label

    !> The range of the values of the shown curves, with a margin, in
    !> ylo and yhi. Returns false if there are no values.
    logical function curves_range(ylo,yhi)
      use ieee_arithmetic, only: ieee_is_finite
      real(c_double), intent(out) :: ylo, yhi

      integer :: i, j
      real*8 :: dy

      ylo = huge(1d0)
      yhi = -huge(1d0)
      do i = 1, w%ln%ncurve
         if (.not.w%ln%curve(i)%show .or. .not.has_values(i)) cycle
         do j = 1, size(w%ln%curve(i)%y)
            if (.not.ieee_is_finite(w%ln%curve(i)%y(j))) cycle
            ! a log axis only shows the positive values
            if (w%ln%yscale == 2 .and. w%ln%curve(i)%y(j) <= 0d0) cycle
            ylo = min(ylo,w%ln%curve(i)%y(j))
            yhi = max(yhi,w%ln%curve(i)%y(j))
         end do
      end do
      curves_range = (yhi >= ylo)
      if (.not.curves_range) return
      if (w%ln%yscale == 2) then
         ylo = ylo / 2d0
         yhi = yhi * 2d0
      else
         dy = max(yhi - ylo,1d-10 * max(abs(ylo),abs(yhi)),1d-30)
         ylo = ylo - 0.05d0 * dy
         yhi = yhi + 0.05d0 * dy
      end if

    end function curves_range

    !> Calculate the line again if its endpoints, the number of points,
    !> or the system (geometry, critical points) changed, and the
    !> values of each curve if the line or anything the curve depends
    !> on changed. Waits for the initialization threads to finish and
    !> for the endpoints to be no longer edited. refit = some values
    !> changed (the plot fits the axes to them).
    subroutine update_line(refit)
      use systems, only: are_threads_running
      use rhoplot, only: rhoplot_line_values
      logical, intent(out) :: refit

      integer*8, allocatable :: key(:)
      integer*8 :: keyref
      integer :: i, j, np, id
      logical :: stale, ok
      real*8 :: x0(3), x1(3), u, dlen
      real*8, allocatable :: f(:)
      character(len=:), allocatable :: lerr

      refit = .false.
      if (xediting .or. are_threads_running()) return

      ! the line: the distances of the points and the markers (the
      ! critical points are those of the reference field)
      keyref = sys(isys)%field_cuid(sys(isys)%iref)
      stale = .not.allocated(w%ln%d) .or. w%ln%keysys /= w%ln%sys%uid .or.&
         any(w%ln%keyxa /= w%ln%xa) .or. any(w%ln%keyxb /= w%ln%xb) .or.&
         w%ln%keynpts /= w%ln%npts .or. w%ln%keytime /= sysc(isys)%timelastchange_geometry .or.&
         w%ln%keytimecp /= sysc(isys)%timelastchange_cplist .or. w%ln%keyref /= keyref
      if (stale) then
         np = w%ln%npts
         u = 1d0
         if (ismol) u = bohrtoa
         dlen = norm2(form_to_cart(isys,w%ln%xb) - form_to_cart(isys,w%ln%xa))
         if (allocated(w%ln%d)) deallocate(w%ln%d)
         allocate(w%ln%d(np))
         do i = 1, np
            w%ln%d(i) = dlen * u * real(i-1,8) / real(np-1,8)
         end do
         ! the cursor stays where it was if the line is long enough;
         ! otherwise, it goes to the middle
         if (w%ln%dcur < 0d0 .or. w%ln%dcur > dlen * u) w%ln%dcur = 0.5d0 * dlen * u
         call find_markers(dlen*u)
         w%ln%gen = w%ln%gen + 1
         w%ln%keysys = w%ln%sys%uid
         w%ln%keyxa = w%ln%xa
         w%ln%keyxb = w%ln%xb
         w%ln%keynpts = w%ln%npts
         w%ln%keytime = sysc(isys)%timelastchange_geometry
         w%ln%keytimecp = sysc(isys)%timelastchange_cplist
         w%ln%keyref = keyref
         refit = .true.
      end if

      ! the curves: a field curve depends on its field, and an
      ! expression on all of them (their number enters through the
      ! size of key). None in a system without atoms (an empty system):
      ! the evaluation of the fields needs the atomic environment
      if (sys(isys)%c%ncel == 0) return
      x0 = form_to_cart(isys,w%ln%keyxa)
      x1 = form_to_cart(isys,w%ln%keyxb)
      np = size(w%ln%d)
      do i = 1, w%ln%ncurve
         associate(cv => w%ln%curve(i))
           if (cv%isexpr) then
              allocate(key(2+sys(isys)%nf+1))
              key(1) = w%ln%gen
              key(2) = 1
              do j = 0, sys(isys)%nf
                 key(3+j) = sys(isys)%field_cuid(j)
              end do
           else
              allocate(key(6))
              key(1) = w%ln%gen
              key(2) = 0
              key(3) = cv%fref%uid
              key(4) = cv%fref%id
              key(5) = cv%nti
              key(6) = sys(isys)%field_cuid(cv%fref%id)
           end if
           stale = .not.allocated(cv%key) .or. .not.allocated(cv%keyexpr)
           if (.not.stale) stale = (size(cv%key) /= size(key))
           if (.not.stale) stale = any(cv%key /= key) .or. (cv%keyexpr /= trim(cv%expr))
           if (stale) then
              refit = .true.
              cv%errmsg = ""
              if (allocated(cv%y)) deallocate(cv%y)
              if (cv%isexpr) then
                 ! checked again: the fields it uses may have come or gone
                 id = -1
                 cv%exprerr = ""
                 if (len_trim(cv%expr) > 0) then
                    call sys(isys)%check_expression(trim(cv%expr),lerr)
                    cv%exprerr = lerr
                 end if
                 ok = (len_trim(cv%expr) > 0 .and. len_trim(cv%exprerr) == 0)
              else
                 id = cv%fref%id
                 ok = sys(isys)%goodfield(id)
                 if (.not.ok) cv%errmsg = "the field is not available in this system"
              end if
              if (ok) then
                 allocate(f(np),cv%y(np))
                 call rhoplot_line_values(sys(isys),x0,x1,np,id,trim(cv%expr),cv%nti,f,cv%y,cv%errmsg)
                 if (len_trim(cv%errmsg) > 0) deallocate(cv%y)
                 deallocate(f)
              end if
              call move_alloc(key,cv%key)
              cv%keyexpr = trim(cv%expr)
           end if
           if (allocated(key)) deallocate(key)
         end associate
      end do

    end subroutine update_line

    !> Find the atoms, and the critical points of the reference field
    !> other than the nuclei, on the line (of length dlen, in the units
    !> of the distances), and keep their distances from A and labels.
    subroutine find_markers(dlen)
      use types, only: realloc
      real*8, intent(in) :: dlen

      integer :: i, k, is, js, ns, nc, nhit, nfound, l(3), lfound(3,27)
      real*8 :: x0(3), x1(3), vl(3), vv, fa(3), fb(3), xs(3), xp(3), t, dd
      real*8, allocatable :: xf(:,:), dhit(:)
      integer, allocatable :: khit(:)
      character(len=12), allocatable :: lbl(:)
      logical :: seen

      w%ln%nmark = 0
      if (allocated(w%ln%dmark)) deallocate(w%ln%dmark)
      if (allocated(w%ln%lmark)) deallocate(w%ln%lmark)
      allocate(w%ln%dmark(0),w%ln%lmark(0))

      x0 = form_to_cart(isys,w%ln%xa)
      x1 = form_to_cart(isys,w%ln%xb)
      vl = x1 - x0
      vv = dot_product(vl,vl)
      if (vv < 1d-20) return

      ! the candidates: the atoms, and the critical points of the
      ! reference field except the nuclei (crystallographic)
      call label_candidates(isys,nc,xf,lbl)

      ! samples along the line, at most half a cell apart in each
      ! crystallographic direction: an image of a candidate on the
      ! line is in the cell of a sample minus the candidate's
      ! position, or the next one along each direction. Molecules
      ! have only the image in the box.
      fa = sys(isys)%c%c2x(x0)
      fb = sys(isys)%c%c2x(x1)
      ns = 1
      if (.not.ismol) ns = ceiling(2d0 * maxval(abs(fb - fa))) + 1

      allocate(dhit(10),khit(10))
      nhit = 0
      do k = 1, nc
         nfound = 0
         do is = 0, ns
            if (ismol) then
               if (is > 0) exit
               l = 0
            else
               xs = fa + (fb - fa) * real(is,8) / real(ns,8)
            end if
            do js = 0, merge(0,7,ismol)
               if (.not.ismol) l = floor(xs - xf(:,k)) + (/ibits(js,0,1),ibits(js,1,1),ibits(js,2,1)/)
               ! each image only once
               seen = .false.
               do i = 1, nfound
                  if (all(lfound(:,i) == l)) seen = .true.
               end do
               if (seen) cycle
               xp = sys(isys)%c%x2c(xf(:,k) + real(l,8))
               t = dot_product(xp - x0,vl) / vv
               if (t < -1d-8 .or. t > 1d0 + 1d-8) cycle
               dd = norm2(x0 + t * vl - xp)
               if (dd > online_tol) cycle
               if (nfound < size(lfound,2)) then
                  nfound = nfound + 1
                  lfound(:,nfound) = l
               end if
               nhit = nhit + 1
               if (nhit > size(dhit)) then
                  call realloc(dhit,2*nhit)
                  call realloc(khit,2*nhit)
               end if
               dhit(nhit) = max(min(t,1d0),0d0) * dlen
               khit(nhit) = k
            end do
         end do
      end do

      ! keep them
      deallocate(w%ln%dmark,w%ln%lmark)
      allocate(w%ln%dmark(nhit),w%ln%lmark(nhit))
      w%ln%nmark = nhit
      w%ln%dmark = dhit(1:nhit)
      w%ln%lmark = lbl(khit(1:nhit))

    end subroutine find_markers

    !> The endpoints of a new system: the first two atoms, or the
    !> diagonal of the cell (of the molecular box) if there are fewer.
    subroutine default_endpoints()

      if (sys(isys)%c%ncel >= 2) then
         w%ln%xa = form_coords(isys,sys(isys)%c%atcel(1)%r)
         w%ln%xb = form_coords(isys,sys(isys)%c%atcel(2)%r)
      else
         w%ln%xa = form_coords(isys,sys(isys)%c%x2c((/0d0,0d0,0d0/)))
         w%ln%xb = form_coords(isys,sys(isys)%c%x2c((/1d0,1d0,1d0/)))
      end if

    end subroutine default_endpoints

    !> Add a curve: the value of the reference field, with the next
    !> color of the palette.
    subroutine add_curve()
      type(line_curve), allocatable :: aux(:)

      if (w%ln%ncurve >= size(w%ln%curve)) then
         allocate(aux(2*size(w%ln%curve)))
         aux(1:w%ln%ncurve) = w%ln%curve(1:w%ln%ncurve)
         call move_alloc(aux,w%ln%curve)
      end if
      w%ln%ncurve = w%ln%ncurve + 1
      w%ln%curve(w%ln%ncurve) = line_curve()
      call w%ln%curve(w%ln%ncurve)%fref%set(sys(isys),sys(isys)%iref)
      w%ln%curve(w%ln%ncurve)%rgb = palette(:,modulo(w%ln%ncolor,npalette)+1)
      w%ln%ncolor = w%ln%ncolor + 1
      w%ln%curve(w%ln%ncurve)%uid = w%ln%ncolor

    end subroutine add_curve

    !> Remove curve i.
    subroutine remove_curve(i)
      integer, intent(in) :: i
      integer :: j

      do j = i, w%ln%ncurve-1
         w%ln%curve(j) = w%ln%curve(j+1)
      end do
      w%ln%curve(w%ln%ncurve) = line_curve()
      w%ln%ncurve = w%ln%ncurve - 1

    end subroutine remove_curve

    !> The curves whose field is gone go to the reference field.
    subroutine fix_fields()
      integer :: i

      do i = 1, w%ln%ncurve
         if (.not.w%ln%curve(i)%fref%ok(sys(isys))) &
            call w%ln%curve(i)%fref%set(sys(isys),sys(isys)%iref)
      end do

    end subroutine fix_fields

    !> The point xf of the form in the absolute frame (with molx0 for
    !> molecules), for the scene.
    function to_abs(xf) result(x)
      real*8, intent(in) :: xf(3)
      real*8 :: x(3)
      x = form_to_cart(isys,xf)
      if (ismol) x = x + sys(isys)%c%molx0
    end function to_abs

    !> Forget a pending pick and release the view it was armed on.
    subroutine cancel_pick()
      if (.not.w%ln%picking) return
      if (w%ln%pickview >= 1 .and. w%ln%pickview <= nwin) &
         call win(w%ln%pickview)%viewmode_release_forced(w%id)
      w%ln%picking = .false.
    end subroutine cancel_pick

  end subroutine draw_line

end submodule line
