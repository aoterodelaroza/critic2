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

! The editor of the plane map objects (the PLANE keyword): the form
! of the map and the plot of its values.
submodule (windows) plane
  use interfaces_cimgui
  use param, only: bohrtoa
  use tools_math, only: niso_lin, niso_log, niso_atan, niso_bader, niso_manual
  use representations, only: plane_npts_min, plane_npts_max, plane_nlev_max, plane_auto_spacing,&
     plane_auto_ntot, plane_auto_tmax, nquant, quant_chem, quant_label
  implicit none

  ! transient markers in the view: the three points (origin, end of
  ! the first side, point on the second side), the region (when the
  ! map itself is not in the view), and the point under the mouse in
  ! the map; opacity of the markers and of the region, and radius of
  ! the points (bohr), a fraction of the scene radius between the two
  ! limits
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

  ! the note of a map that is not updated (the Update map checkbox is off)
  character(len=*), parameter :: noupdate_msg = "The map is not updated: check Update map to&
     & calculate it again"

contains

  !> The form of the map d of system isys (the editor of the plane map
  !> objects): the three points, the region, the margin, the number of
  !> points, the buttons that move the plane, the quantity, the
  !> colors, and the contour lines. Returns whether anything changed;
  !> a drag counts as a change only when it ends (the map is
  !> calculated again on each change). Each point has a Pick button:
  !> ipick is the one pressed (0 = none), and pickoff disables them.
  !> The contour lines have a Values checkbox (showval; the labels of
  !> the levels in the plot of the map, which are not part of the map
  !> in the view) and a control for their size (valscale).
  module function plane_def_form(d,isys,ttshown,ipick,pickoff,showval,valscale,sec) result(changed)
    use systems, only: sys
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_dragfloat_realc, iw_calcwidth,&
       iw_inputint, iw_combo_simple, iw_checkbox, iw_coloredit, iw_inputtext, iw_cmap_optstr,&
       iw_separatortext
    use tools_io, only: string, lower
    type(plane_def), intent(inout) :: d
    integer, intent(in) :: isys
    logical, intent(inout) :: ttshown
    integer, intent(out) :: ipick
    logical, intent(in) :: pickoff
    logical, intent(inout) :: showval
    real(c_float), intent(inout) :: valscale
    type(iw_section), intent(inout) :: sec(plane_nsec)
    logical :: changed

    integer :: i
    logical :: ismol, ch, ch2, okaxes
    real*8 :: axes(3,0:2), a0(3), x1(3), x2(3), e2(3)
    character(len=:), allocatable :: errmsg, noaxes
    ! the axes of the plane for the Move and Extend buttons: 0 = its
    ! normal, 1 = u (the first side), 2 = v (the second side)
    character(len=1), parameter :: axlbl(0:2) = (/"n","u","v"/)
    character(len=*), parameter :: axdesc(0:2) = [character(len=19) :: "n (its normal)",&
       "u (its first side)","v (its second side)"]

    changed = .false.
    ipick = 0
    ismol = sys(isys)%c%ismolecule

    ! whether the map follows the form (on again: calculated if
    ! stale). A toggle is a change: the flag travels with the form
    ! between a window and the object it edits
    changed = iw_checkbox("Update map##planeupdate",d%update,labelleft=.true.)
    call iw_tooltip("Calculate the map again whenever the plane, the property, the colors, or the&
       & contours change. Uncheck it to make several changes without waiting for the map, then&
       & check it again to calculate it",ttshown)

    ! the plane, the region, and the number of points
    call iw_separatortext("Plane",sec(plane_sec_plane),hasmore=.true.)
    if (sec(plane_sec_plane)%open) then
       do i = 1, 3
          call draw_point(i)
       end do
       ! the axes (unit vectors, Cartesian), after the points may have
       ! changed; the buttons that use them are disabled if the points
       ! do not define a plane
       call d%corners(isys,a0,x1,x2,axes(:,1),e2,axes(:,0),errmsg)
       okaxes = (len(errmsg) == 0)
       if (okaxes) then
          axes(:,2) = (x2 - a0) / norm2(x2 - a0)
          noaxes = ""
       else
          noaxes = " (disabled: " // lower(errmsg(1:1)) // errmsg(2:) // ")"
       end if
       call draw_move()
       if (sec(plane_sec_plane)%more) call draw_region()
    end if

    ! the quantity, the colors, and the contour lines (nothing more if
    ! they are not shown)
    call iw_separatortext("Property",sec(plane_sec_property))
    if (sec(plane_sec_property)%open) call draw_quantity()
    call iw_separatortext("Colors",sec(plane_sec_colors),hasmore=.true.)
    if (sec(plane_sec_colors)%open) call draw_colors()
    call iw_separatortext("Contours",sec(plane_sec_contours),hasmore=d%showcont)
    if (sec(plane_sec_contours)%open) call draw_contours()

  contains
    !> A drag of a real value (x1 or the three of x3): a change when
    !> the drag (or the typing) ends.
    subroutine drag(str,x1,x3,speed,decimal,sameline,min,sig)
      character(len=*), intent(in) :: str
      real*8, intent(inout), optional :: x1, x3(3)
      real*8, intent(in) :: speed
      integer, intent(in), optional :: decimal
      logical, intent(in), optional :: sameline
      real*8, intent(in), optional :: min
      integer, intent(in), optional :: sig

      logical :: ch

      ch = iw_dragfloat_real8(str,x1=x1,x3=x3,speed=speed,decimal=decimal,sameline=sameline,min=min,&
         sig=sig)
      changed = changed .or. logical(igIsItemDeactivatedAfterEdit())

    end subroutine drag

    !> A fixed range of values: the checkbox flag and, if it is on,
    !> the drags of xmin and xmax, with the tooltips of the three.
    !> sfx = the ID suffix of the three widgets.
    subroutine draw_range(sfx,flag,xmin,xmax,ttflag,ttmin,ttmax)
      character(len=*), intent(in) :: sfx, ttflag, ttmin, ttmax
      logical, intent(inout) :: flag
      real*8, intent(inout) :: xmin, xmax

      real*8 :: smin, smax
      character(len=*), parameter :: tttype = " (ctrl+click to type a value, e.g. 1e-6)"

      ch = iw_checkbox("Fixed range##plane" // sfx // "range",flag,labelleft=.true.)
      changed = changed .or. ch
      call iw_tooltip(ttflag,ttshown)
      if (.not.flag) return
      ! general format, and a drag proportional to the value (the range
      ! may be very small, e.g. a density far from the nuclei); a zero
      ! starts at the scale of the other bound
      smin = abs(xmin)
      if (smin == 0d0) smin = abs(xmax)
      smax = abs(xmax)
      if (smax == 0d0) smax = abs(xmin)
      call iw_text("Min",alignframe=.true.,gap=.true.)
      call drag("##plane" // sfx // "min",x1=xmin,speed=0.01d0 * max(smin,1d-6),sig=4,sameline=.true.)
      call iw_tooltip(ttmin // tttype,ttshown)
      call iw_text("Max",alignframe=.true.,gap=.true.)
      call drag("##plane" // sfx // "max",x1=xmax,speed=0.01d0 * max(smax,1d-6),sig=4,sameline=.true.)
      call iw_tooltip(ttmax // tttype,ttshown)

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
         call iw_tooltip("The origin of the plane",ttshown)
      elseif (ip == 2) then
         call iw_tooltip("The end of the first side of the plane",ttshown)
      else
         call iw_tooltip("A point on the second side of the plane: its end, for the parallelogram;&
            & the plane through it, for the rectangle",ttshown)
      end if
      if (iw_button("Pick" // sfx // "pick",sameline=.true.,disabled=pickoff)) ipick = ip
      call iw_tooltip("Pick this point in the view: an atom, a critical point, or a point&
         & in empty space",ttshown)
      if (ismol) then
         call drag("(Å)" // sfx,x3=d%xp(:,ip),speed=0.01d0,decimal=5,sameline=.true.)
      else
         call drag("(frac)" // sfx,x3=d%xp(:,ip),speed=0.001d0,decimal=6,sameline=.true.)
      end if
      call iw_tooltip("Coordinates of this point",ttshown)

    end subroutine draw_point

    !> The buttons that displace the whole plane by d%step along the
    !> Cartesian axes and, in crystals, along the lattice vectors
    !> (first line), and along its normal and its two sides (second
    !> line), with the step.
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
                  dx(j) = st / bohrtoa
               else
                  ! along a lattice vector, as a unit vector
                  dx(j-3) = 1d0
                  dx = sys(isys)%c%x2c(dx) / sys(isys)%c%aa(j-3) * st / bohrtoa
               end if
               call shift(dx)
               changed = .true.
            end if
            call iw_tooltip("Displace the plane along " // sg // lbl(j) // " (" //&
               merge("Cartesian axis","lattice vector",j <= 3) // ")",ttshown)
         end do
      end do

      ! along the normal and the sides
      call iw_text("Normal",alignframe=.true.)
      call iw_tooltip("Displace the whole plane along its normal (the direction of the&
         & first side times the second)",ttshown)
      do j = 0, 2
         if (j == 1) then
            call iw_text("In plane",alignframe=.true.,gap=.true.)
            call iw_tooltip("Displace the whole plane along its sides: u is the first side and v the&
               & second",ttshown)
         end if
         do k = 1, 2
            sg = merge("+","-",k == 1)
            st = merge(d%step,-d%step,k == 1)
            if (iw_button(sg // axlbl(j) // "##planemove",sameline=.true.,disabled=.not.okaxes)) then
               call shift(st / bohrtoa * axes(:,j))
               changed = .true.
            end if
            call iw_tooltip("Displace the plane along " // sg // trim(axdesc(j)) // noaxes,ttshown,&
               whendisabled=.true.)
         end do
      end do
      call iw_text("Step (Å)",alignframe=.true.,gap=.true.)
      call igPushItemWidth(iw_calcwidth(7,1))
      ch = iw_dragfloat_real8("##planestep",x1=d%step,speed=0.01d0,decimal=3,sameline=.true.,min=0d0)
      call igPopItemWidth()
      call iw_tooltip("Length of each displacement and extension",ttshown)

    end subroutine draw_move

    !> Displace the points ip of the plane (all three if not given)
    !> by the Cartesian vector dxc (bohr).
    subroutine shift(dxc,ip)
      real*8, intent(in) :: dxc(3)
      integer, intent(in), optional :: ip

      integer :: i

      do i = 1, 3
         if (present(ip)) then
            if (i /= ip) cycle
         end if
         d%xp(:,i) = form_coords(isys,form_to_cart(isys,d%xp(:,i)) + dxc)
      end do

    end subroutine shift

    !> More of the plane: the region (a rectangle or a parallelogram)
    !> and the margin around it, then the number of points along each
    !> side (shown, but not editable, if automatic).
    subroutine draw_region()
      integer :: j, k
      real*8 :: dx(3)
      character(len=1) :: sg

      call iw_text("Region",alignframe=.true.)
      call iw_combo_simple("##planeregion","Rectangle" // c_null_char // "Parallelogram" // c_null_char,&
         d%region,sameline=.true.,changed=ch)
      changed = changed .or. ch
      call iw_tooltip("Rectangle: the first side goes from the origin to the X end, and the&
         & second is perpendicular to it, as long as needed to reach the Y point. Parallelogram:&
         & the sides go from the origin to the X end and to the Y point",ttshown)
      call iw_text("Margin (Å)",alignframe=.true.)
      call drag("##planemargin",x1=d%margin,speed=0.01d0,decimal=3,sameline=.true.)
      call iw_tooltip("Length added to the region on every side (negative removes it)",ttshown)

      ! extend the region by the step along its sides: + moves the end
      ! of the side, - moves the other two points back
      call iw_text("Extend",alignframe=.true.,gap=.true.)
      call iw_tooltip("Lengthen the region by the step on one side: u is the first side and v the&
         & second, + at the end of the side and - at the origin",ttshown)
      do j = 1, 2
         do k = 1, 2
            sg = merge("+","-",k == 1)
            if (iw_button(sg // axlbl(j) // "##planeextend",sameline=.true.,disabled=.not.okaxes)) then
               dx = d%step / bohrtoa * axes(:,j)
               if (k == 2) call shift(-dx)
               call shift(dx,j+1)
               changed = .true.
            end if
            call iw_tooltip("Extend the region by the step at the " // trim(merge("end   ","origin",k == 1)) //&
               " of its " // trim(merge("first ","second",j == 1)) // " side (" // sg // axlbl(j) // ")" //&
               noaxes,ttshown,whendisabled=.true.)
         end do
      end do

      call iw_text("Points",alignframe=.true.)
      if (d%autonpts) call igBeginDisabled(.true._c_bool)
      ch = iw_inputint("##planenx",d%nx,width=6,sameline=.true.,notlive=.true.)
      call iw_tooltip("Number of points along the first side (" // string(plane_npts_min) // " to " //&
         string(plane_npts_max) // ")",ttshown,whendisabled=.true.)
      call iw_text("x",alignframe=.true.,sameline=.true.)
      ch2 = iw_inputint("##planeny",d%ny,width=6,sameline=.true.,notlive=.true.)
      call iw_tooltip("Number of points along the second side (" // string(plane_npts_min) // " to " //&
         string(plane_npts_max) // ")",ttshown,whendisabled=.true.)
      if (d%autonpts) then
         call igEndDisabled()
      else
         d%nx = max(min(d%nx,plane_npts_max),plane_npts_min)
         d%ny = max(min(d%ny,plane_npts_max),plane_npts_min)
         changed = changed .or. ch .or. ch2
      end if
      ch = iw_checkbox("Automatic##planeautonpts",d%autonpts,gap=.true.,labelleft=.true.)
      changed = changed .or. ch
      call iw_tooltip("Choose the number of points from the size of the region: square cells&
         & of " // string(plane_auto_spacing,'f',decimal=2) // " Å, coarser for very large&
         & regions (at most " // string(plane_auto_ntot) // " points) and for slow fields (so that&
         & a map takes at most about " // string(nint(plane_auto_tmax)) // " s)",ttshown)

    end subroutine draw_region

    !> The quantity on the map: a quantity of a field, or an
    !> arithmetic expression (checked as it is typed).
    subroutine draw_quantity()

      ch = quantity_form(isys,d%isexpr,d%fref,d%nti,d%expr,d%exprerr,"plane","the map","the plane",ttshown)
      changed = changed .or. ch

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
         & absolute value, or by their arctangent",ttshown)
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

      ch = iw_checkbox("Show##planeshowcont",d%showcont,labelleft=.true.)
      changed = changed .or. ch
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
      call iw_tooltip("The contour levels: evenly spaced, logarithmic (positive and negative), evenly&
         & spaced in the arctangent, the Bader set (1, 2, 4, 8 x 10^n), or a list of values",ttshown)
      if (d%ctype == niso_manual) then
         call iw_text("List",alignframe=.true.,gap=.true.)
         ch = iw_inputtext("##planellist",bufsize=1023,textf=d%llist,width=30,sameline=.true.)
         changed = changed .or. ch
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
      ch = iw_coloredit("##planergbpos",rgb=d%rgbpos,sameline=.true.)
      changed = changed .or. ch
      call iw_tooltip("Color of the contour lines of the positive levels",ttshown)
      call iw_text("Negative color",alignframe=.true.,gap=.true.)
      ch = iw_coloredit("##planergbneg",rgb=d%rgbneg,sameline=.true.)
      changed = changed .or. ch
      call iw_tooltip("Color of the contour lines of the zero and negative levels",ttshown)
      ch = iw_checkbox("Values##planeshowcval",showval,gap=.true.,labelleft=.true.)
      call iw_tooltip("Write the value of the level on the contour lines, in their color,&
         & where there is room",ttshown)
      if (showval) then
         call iw_text("Size",alignframe=.true.,gap=.true.)
         ch = iw_dragfloat_realc("##planevscale",x1=valscale,speed=0.01_c_float,min=0.1_c_float,&
            max=10._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
         call iw_tooltip("Size of the values, relative to the text of the interface",ttshown)
      end if

    end subroutine draw_contours
  end function plane_def_form

  !> The editor of a plane map object (reptype_slice): the map (the
  !> three points, typed or picked in the view, the region, the number
  !> of points, the quantity, the colors, and the contour lines, as
  !> PLANE calculates them), how it is drawn in the view, the map
  !> itself as a plot with its color bar and the labels of the atoms
  !> and critical points near the plane, and the table of its values.
  !> The values are those of the object (slice_update_cache), never
  !> calculated here. Returns whether the object changed.
  module function draw_editrep_slice(w,ttshown) result(changed)
    use systems, only: sys, sysc, are_threads_running
    use gui_main, only: g
    use representations, only: slice_cache, plane_def, slice_update_cache, rep_shape,&
       shapekind_sphere, shapekind_box, bondrad_def
    use utils, only: iw_text, iw_tooltip, iw_dragfloat_real8, iw_dragfloat_realc, iw_separatortext,&
       iw_separatortext_height, iw_calcwidth, iw_calcheight, iw_checkbox, iw_coloredit, iw_combo_simple,&
       string_nozeros
    use tools_io, only: string
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch, ismol, hover, mapshown, refit, ldum
    integer :: isys, iview, istat, i, ipick, nline
    real*8 :: xc(3), rad, xhov(3), x0(3), uu(3), vv(3), a0(3), e1(3), e2(3), v(3,3)
    character(len=:), allocatable :: errmsg
    type(ImVec2) :: sz0
    real(c_float) :: htable, hmap, hres, frac
    type(rep_shape), allocatable :: shp(:)
    type(slice_cache), pointer :: c
    type(plane_def), pointer :: d

    integer, parameter :: map_minlines = 8 ! minimum height of the map (lines)

    changed = .false.
    isys = w%isys
    iview = w%anchor_view()
    ismol = sys(isys)%c%ismolecule
    if (.not.allocated(w%rep%slice%c)) allocate(w%rep%slice%c)
    c => w%rep%slice%c
    d => w%rep%slice%d

    ! a pending pick of one of the points, commanded to the anchor view
    if (w%editrep_pick_item > 0) then
       call view_pick_result(iview,w%id,isys,w%editrep_pick,istat,xc)
       if (istat == ipick_point) then
          d%xp(:,w%editrep_pick_item) = form_coords(isys,xc)
          changed = .true.
       end if
       if (istat /= ipick_pending) w%editrep_pick_item = 0
    end if

    ! the map: the plane, the region, the number of points, the
    ! quantity, the colors, and the contour lines; a Pick button arms
    ! a pick in the view, received above
    ch = plane_def_form(d,isys,ttshown,ipick,w%editrep_pick_item > 0,w%smap%showcval,w%smap%vscale,&
       w%editrep_slicesec(1:plane_nsec))
    changed = changed .or. ch
    if (ipick > 0) then
       w%editrep_pick_item = ipick
       call w%editrep_pick%arm()
       call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick " // trim(pt_name(ipick)) // ": an atom,&
          & a critical point, or a point in space",w%id,acceptempty=.true.)
    end if

    ! how it is drawn in the view
    call iw_separatortext("Style",w%editrep_slicesec(slicesec_style))
    if (w%editrep_slicesec(slicesec_style)%open) then
       call iw_text("Opacity",alignframe=.true.)
       ch = iw_dragfloat_realc("##sliceopacity",x1=w%rep%slice%alpha,speed=0.01_c_float,&
          min=0._c_float,max=1._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       changed = changed .or. ch
       call iw_tooltip("Opacity of the map (1 = opaque)",ttshown)
       call iw_text("Line radius",alignframe=.true.,gap=.true.)
       ch = iw_dragfloat_real8("(Å)##slicecrad",x1=w%rep%slice%crad,speed=0.001d0,&
          min=0d0,scale=bohrtoa,decimal=3,sameline=.true.)
       changed = changed .or. ch
       call iw_tooltip("Radius of the contour lines",ttshown)
       ch = relief_checkbox(w%rep%slice%rl,"slice",.false.,ttshown)
       changed = changed .or. ch
       ch = relief_form(w%rep%slice%rl,"slice",.true.,ttshown)
       changed = changed .or. ch
    end if

    ! the values of the object: the view calculates them when it draws
    ! the object; here if the object is hidden, while something in
    ! this window shows them
    if (.not.w%rep%shown .and. (w%editrep_slicesec(slicesec_map)%open .or.&
       w%editrep_slicesec(slicesec_values)%open)) call slice_update_cache(w%rep)
    ! a shown object whose values the view could not calculate while
    ! the initialization threads ran gets them once they finish
    if (w%rep%shown .and. c%gen == 0 .and. d%update .and. sys(isys)%c%ncel > 0 .and.&
       .not.are_threads_running()) win(iview)%sc%forcebuildlists = .true.
    mapshown = allocated(c%ff)
    if (sys(isys)%c%ncel == 0) &
       call iw_text("This system has no atoms: there are no field values to map",danger=.true.,wrap=.true.)

    ! the map, with the labels of the atoms and critical points near
    ! the plane (more: what is labeled and how)
    hover = .false.
    refit = .false.
    call iw_separatortext("Map",w%editrep_slicesec(slicesec_map),hasmore=w%smap%showlabels)
    if (w%editrep_slicesec(slicesec_map)%open) then
       call draw_labels_option()
       if (mapshown) then
          call check_grid(refit)
          call update_contours()
          call update_texture()
          call update_marks()
       end if
    end if

    ! the map takes the rest of the window, leaving room for what
    ! follows it: the note about the labels, the table of values (its
    ! title, and the table if it is open, which takes part of the
    ! space), the messages, and the frame of the plot. The window
    ! grows to fit its content:
    ! the stretch of the map and the table beyond their minimum
    ! heights is elastic, not content to fit
    call igGetContentRegionAvail(sz0)
    nline = 0
    if (w%editrep_slicesec(slicesec_map)%open .and. w%smap%showlabels .and. w%smap%marksmany) nline = nline + 1
    if (.not.d%update) nline = nline + 1
    if (allocated(c%errmsg)) then
       if (len_trim(c%errmsg) > 0) nline = nline + 1
    end if
    if (allocated(w%errmsg)) then
       if (len_trim(w%errmsg) > 0) nline = nline + 1
    end if
    hres = iw_calcheight(0,nline,.true.) + g%Style%ItemSpacing%y
    if (mapshown) hres = hres + iw_separatortext_height()
    htable = 0._c_float
    if (mapshown .and. w%editrep_slicesec(slicesec_values)%open) then
       frac = merge(0.35_c_float,1._c_float,w%editrep_slicesec(slicesec_map)%open)
       htable = max(frac * (sz0%y - hres),iw_calcheight(4,0,.false.))
       hres = hres + htable + iw_calcheight(0,0,.true.)
       w%heightslack = htable - iw_calcheight(4,0,.false.)
    end if
    if (w%editrep_slicesec(slicesec_map)%open) then
       hmap = iw_calcheight(map_minlines,0,.false.)
       sz0%y = max(sz0%y - hres,hmap)
       if (mapshown) then
          w%heightslack = w%heightslack + sz0%y - hmap
          call draw_map(sz0,refit,hover,xhov,.false.)
       elseif (are_threads_running()) then
          call iw_text("Waiting for the initialization threads to finish",disabled=.true.)
       end if
    end if

    ! the table of values
    call draw_data_section()

    ! the three points, the region (unless the map itself is in the
    ! view), and the point under the mouse in the anchor view, for as
    ! long as the window is open
    if (associated(win(iview)%sc)) then
       if (win(iview)%sc%isinit /= 0) then
          rad = min(max(0.01d0 * real(win(iview)%sc%scenerad,8),marker_radmin),marker_radmax)
          allocate(shp(3))
          do i = 1, 3
             shp(i) = rep_shape(kind=shapekind_sphere,x1=cart_to_abs(form_to_cart(isys,d%xp(:,i))),&
                rad=rad,rgb=pt_rgb(:,i),alpha=marker_alpha)
          end do
          if (.not.w%rep%shown .or. .not.mapshown) then
             call d%grid(isys,x0,uu,vv,a0,e1,e2,errmsg)
             if (len(errmsg) == 0) then
                v = 0d0
                v(:,1) = real(d%nx-1,8) * uu
                v(:,2) = real(d%ny-1,8) * vv
                shp = [shp,rep_shape(kind=shapekind_box,x1=cart_to_abs(x0),v=v,rad=0.15d0*rad,&
                   rgb=region_rgb,alpha=region_alpha)]
             end if
          end if
          if (hover) &
             shp = [shp,rep_shape(kind=shapekind_sphere,x1=cart_to_abs(xhov),&
             rad=max(rad,1.6d0*bondrad_def),rgb=cur_rgb,alpha=1._c_float)]
          call win(iview)%sc%show_transient_shapes(w%uid,1,shp)
       end if
    end if

    ! the map stays as it is until it is updated again, and what is wrong
    if (.not.d%update) call iw_text(noupdate_msg,wrap=.true.)
    if (allocated(c%errmsg)) then
       if (len_trim(c%errmsg) > 0) call iw_text(c%errmsg,danger=.true.,wrap=.true.)
    end if
    if (allocated(w%errmsg)) then
       if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)
    end if

    ! the image of the map, as drawn in the Map section this frame,
    ! to the file chosen in the dialog
    if (w%okfile_set .and. mapshown .and. w%editrep_slicesec(slicesec_map)%open) call export_map()
    w%okfile_set = .false.

  contains

    !> The labels of the atoms and critical points near the plane.
    !> More: what is labeled, the color and size of the labels, their
    !> text, and the distance to the plane.
    subroutine draw_labels_option()

      ldum = iw_checkbox("Labels##slicelabels",w%smap%showlabels,labelleft=.true.)
      call iw_tooltip("Label the atoms and the critical points of the reference field near the plane&
         & in the map below",ttshown)
      if (.not.w%editrep_slicesec(slicesec_map)%more .or. .not.w%smap%showlabels) return

      ! more: what is labeled, and the color and size of the labels
      ldum = iw_checkbox("Critical points##slicecps",w%smap%showcps,labelleft=.true.)
      call iw_tooltip("Label also the critical points of the reference field other than the nuclei&
         & (n = non-nuclear maximum, b = bond, r = ring, c = cage)",ttshown)
      call iw_text("Color",alignframe=.true.,gap=.true.)
      ldum = iw_coloredit("##slicergblabel",rgb=w%smap%rgblabel,sameline=.true.)
      call iw_tooltip("Color of the labels and their dots",ttshown)
      call iw_text("Size",alignframe=.true.,gap=.true.)
      ldum = iw_dragfloat_realc("##slicelscale",x1=w%smap%lscale,speed=0.01_c_float,min=0.1_c_float,&
         max=10._c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
      call iw_tooltip("Size of the labels and their dots, relative to the text of the interface",ttshown)

      ! more: the text, and the distance to the plane
      call iw_text("Text",alignframe=.true.)
      call iw_combo_simple("##sliceltext","Species name"//c_null_char//"Atomic symbol"//c_null_char//&
         "Atom name"//c_null_char,w%smap%ltext,sameline=.true.)
      call iw_tooltip("Text in the labels of the atoms",ttshown)
      call iw_text("Distance",alignframe=.true.,gap=.true.)
      ldum = iw_dragfloat_real8("##slicehmax",x1=w%smap%hmax,speed=0.01d0,min=0d0,decimal=3,sameline=.true.)
      call iw_tooltip("Largest distance to the plane of the labeled atoms and critical points (Å)",&
         ttshown)
      ldum = iw_checkbox("Show distance##sliceshowdist",w%smap%showdist,gap=.true.,labelleft=.true.)
      call iw_tooltip("Add to each label the signed distance of the atom or critical point to the&
         & plane (Å)",ttshown)

    end subroutine draw_labels_option

    !> Whether the plot must fit its axes to the region again: it is
    !> not the one last plotted (relative to a0, the origin of the
    !> plot: moving the plane as a whole keeps the axes). Keeps the
    !> region plotted.
    subroutine check_grid(refit)
      logical, intent(out) :: refit

      real*8 :: region(3,3)

      region(:,1) = c%x0 - c%a0
      region(:,2) = c%uu * (size(c%ff,1)-1)
      region(:,3) = c%vv * (size(c%ff,2)-1)
      refit = any(abs(region - w%smap%region) > 1d-10)
      if (refit) w%smap%region = region

    end subroutine check_grid

    !> The map (size sz): the colors, the contour lines, and the labels,
    !> with the color bar on its right. If refit, fit the axes to the
    !> region. Returns in hover whether the mouse is on the region and
    !> the point under it (xhov, cell-frame Cartesian bohr). If forprint,
    !> draw it for export_map: with the axes as last shown in the
    !> window and no mouse readouts.
    subroutine draw_map(sz,refit,hover,xhov,forprint)
      use utils, only: iw_colormap_id
      type(ImVec2), intent(in) :: sz
      logical, intent(in) :: refit
      logical, intent(in) :: forprint
      logical, intent(out) :: hover
      real*8, intent(out) :: xhov(3)

      character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
      integer :: j, nx, ny
      integer(c_int) :: flags, coltext
      real(c_double) :: xmin, xmax, ymin, ymax, crn(2,4), mu, mv
      real(c_float) :: px(4), py(4), wbar, fs, fv
      real(c_float), allocatable :: box(:,:), ang(:)
      real*8, allocatable :: gap(:), sval(:)
      integer, allocatable :: ival(:), icand(:)
      integer :: nbox, nval, k
      integer(c_int) :: col(2)
      type(ImVec2) :: pmin, pmax, szp, pt, tsz
      type(c_ptr) :: dl, font

      hover = .false.
      xhov = 0d0
      if (.not.allocated(c%ff)) return
      nx = size(c%ff,1)
      ny = size(c%ff,2)

      ! the corners of the region in plot coordinates
      crn(:,1) = to_plot(c%x0)
      crn(:,2) = to_plot(c%x0 + real(nx-1,8) * c%uu)
      crn(:,3) = to_plot(c%x0 + real(nx-1,8) * c%uu + real(ny-1,8) * c%vv)
      crn(:,4) = to_plot(c%x0 + real(ny-1,8) * c%vv)

      ! the plot leaves room for the color bar on its right
      wbar = iw_calcwidth(10,0)
      szp = sz
      szp%x = max(sz%x - wbar,iw_calcwidth(20,0))

      str1 = "##planemap" // c_null_char
      flags = ior(ImPlotFlags_NoTitle,ImPlotFlags_Equal)
      flags = ior(flags,ImPlotFlags_NoLegend)
      if (forprint) then
         flags = ior(flags,ImPlotFlags_NoInputs)
      else
         flags = ior(flags,ImPlotFlags_Crosshairs)
         w%smap%szmap = (/sz%x,sz%y/)
      end if
      if (ipBeginPlot(c_loc(str1),szp,flags)) then
         str2 = "u (Å)" // c_null_char
         str3 = "v (Å)" // c_null_char
         call ipSetupAxes(c_loc(str2),c_loc(str3),ImPlotAxisFlags_None,ImPlotAxisFlags_None)
         if (forprint) then
            call ipSetupAxisLimits(ImAxis_X1,w%smap%lims(1),w%smap%lims(2),ImPlotCond_Always)
            call ipSetupAxisLimits(ImAxis_Y1,w%smap%lims(3),w%smap%lims(4),ImPlotCond_Always)
         elseif (refit) then
            call ipSetupAxisLimits(ImAxis_X1,minval(crn(1,:)),maxval(crn(1,:)),ImPlotCond_Always)
            call ipSetupAxisLimits(ImAxis_Y1,minval(crn(2,:)),maxval(crn(2,:)),ImPlotCond_Always)
         end if

         ! the colors: the texture on the region, inside the plot
         call ipGetPlotCurrentLimits(xmin,xmax,ymin,ymax)
         if (.not.forprint) w%smap%lims = (/xmin,xmax,ymin,ymax/)
         call ipPlotToPixels(xmin,ymin,ImAxis_X1,ImAxis_Y1,px(1),py(1))
         call ipPlotToPixels(xmax,ymax,ImAxis_X1,ImAxis_Y1,px(2),py(2))
         pmin%x = min(px(1),px(2))
         pmin%y = min(py(1),py(2))
         pmax%x = max(px(1),px(2))
         pmax%y = max(py(1),py(2))
         dl = igGetWindowDrawList()
         call ImDrawList_PushClipRect(dl,pmin,pmax,.true._c_bool)
         if (w%smap%tex /= 0) then
            do j = 1, 4
               call ipPlotToPixels(crn(1,j),crn(2,j),ImAxis_X1,ImAxis_Y1,px(j),py(j))
            end do
            ! the texels are centered on the points of the grid
            call ImDrawList_AddImageQuad(dl,int(w%smap%tex,c_intptr_t),&
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
         coltext = igGetColorU32_Vec4(ImVec4(w%smap%rgblabel(1),w%smap%rgblabel(2),w%smap%rgblabel(3),1._c_float))
         font = igGetFont()
         fs = w%smap%lscale * igGetFontSize()
         fv = w%smap%vscale * igGetFontSize()
         nbox = 0
         nval = 0
         allocate(box(4,w%smap%nmark+w%smap%ncval),ival(w%smap%ncval),icand(w%smap%ncval),sval(w%smap%ncval),&
            gap(w%smap%ncval),ang(w%smap%ncval))
         if (w%smap%showlabels) then
            do j = 1, w%smap%nmark
               call ipPlotToPixels(w%smap%xmark(1,j),w%smap%xmark(2,j),ImAxis_X1,ImAxis_Y1,px(1),py(1))
               str1 = trim(w%smap%lmark(j)) // c_null_char
               call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
               nbox = nbox + 1
               box(:,nbox) = (/px(1) - 3._c_float * w%smap%lscale,py(1) - (tsz%y + 1._c_float) * w%smap%lscale,&
                  px(1) + (4._c_float + tsz%x) * w%smap%lscale,py(1) + 3._c_float * w%smap%lscale/)
            end do
         end if
         if (d%showcont) then
            if (w%smap%showcval .and. w%smap%ncval > 0) &
               call place_values(pmin,pmax,box,nbox,nval,ival,icand,sval,gap,ang)
            if (w%smap%ncpos > 0) call plot_contours(1,d%rgbpos,w%smap%cpos,w%smap%ncpos,nval,ival,sval,gap)
            if (w%smap%ncneg > 0) call plot_contours(2,d%rgbneg,w%smap%cneg,w%smap%ncneg,nval,ival,sval,gap)
         end if

         call ImDrawList_PushClipRect(dl,pmin,pmax,.true._c_bool)
         col(1) = igGetColorU32_Vec4(ImVec4(d%rgbpos(1),d%rgbpos(2),d%rgbpos(3),1._c_float))
         col(2) = igGetColorU32_Vec4(ImVec4(d%rgbneg(1),d%rgbneg(2),d%rgbneg(3),1._c_float))
         do j = 1, nval
            k = ival(j)
            call ipPlotToPixels(w%smap%xcval(1,icand(j),k),w%smap%xcval(2,icand(j),k),ImAxis_X1,ImAxis_Y1,px(1),py(1))
            str1 = fmt_level(w%smap%zcval(k)) // c_null_char
            call own_add_text_rotated(dl,font,fv,px(1),py(1),ang(j),col(w%smap%icval(1,k)),c_loc(str1))
         end do
         if (w%smap%showlabels) then
            do j = 1, w%smap%nmark
               call ipPlotToPixels(w%smap%xmark(1,j),w%smap%xmark(2,j),ImAxis_X1,ImAxis_Y1,px(1),py(1))
               call ImDrawList_AddCircleFilled(dl,ImVec2(px(1),py(1)),3._c_float*w%smap%lscale,coltext,8_c_int)
               str1 = trim(w%smap%lmark(j)) // c_null_char
               call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
               pt%x = px(1) + 4._c_float * w%smap%lscale
               pt%y = py(1) - (tsz%y + 1._c_float) * w%smap%lscale
               call ImDrawList_AddText_FontPtr(dl,font,fs,pt,coltext,c_loc(str1),c_null_ptr,&
                  0._c_float,c_null_ptr)
            end do
         end if
         call ImDrawList_PopClipRect(dl)

         ! the point under the mouse
         if (.not.forprint .and. ipIsPlotHovered()) then
            call ipGetPlotMousePos(mu,mv)
            call plot_point(mu,mv,hover,xhov)
         end if
         call ipEndPlot()
      end if

      ! the color bar
      call igSameLine(0._c_float,-1._c_float)
      str1 = "##planecbar" // c_null_char
      str2 = "%.3g" // c_null_char
      call ipColormapScale(c_loc(str1),w%smap%tmin,w%smap%tmax,ImVec2(wbar,szp%y),c_loc(str2),&
         iw_colormap_id(d%icmap))
      if (forprint) return
      if (d%cscale == 1) then
         call iw_tooltip("Colors of the map: log10 of the absolute value",ttshown)
      elseif (d%cscale == 2) then
         call iw_tooltip("Colors of the map: 2/pi atan of the value",ttshown)
      else
         call iw_tooltip("Colors of the map: the value",ttshown)
      end if

      if (w%smap%showlabels .and. w%smap%marksmany) &
         call iw_text("(The region is too large for labels)",disabled=.true.)

    end subroutine draw_map

    !> Save the map, as last drawn in the window, to the PNG file
    !> w%okfile, for a publication: drawn again in a separate ImGui
    !> context, black on white, and rendered into a framebuffer
    !> w%smap%imgscale times larger than the window.
    subroutine export_map()
      use interfaces_opengl3
      use interfaces_stb
      use gui_main, only: shader_version, font_ranges
      use tools_io, only: uout
      real(c_float), parameter :: bakemax = 128._c_float ! largest size of the font atlas glyphs (pixels)

      type(c_ptr) :: pc, dd
      character(kind=c_char,len=:), allocatable, target :: str
      character(kind=c_char), allocatable, target :: pix(:)
      integer(c_int), target :: fbo, rbo
      integer(c_int) :: wid, hei, idum, npix(2)
      real(c_float) :: wsz(2), xs
      real*8 :: xdum(3)
      logical :: ldum
      integer :: i

      ! the map in the print context, at its size in the window (two
      ! frames: the plot lays itself out in the first)
      call plane_image_size(w%smap%szmap,w%smap%imgscale,wsz,xs,npix)
      wid = npix(1)
      hei = npix(2)
      str = shader_version
      pc = ipPrintContextCreate(igGetFontSize(),xs,bakemax,font_dejavu_base85_ptr,c_loc(font_ranges),c_loc(str))
      do i = 1, 2
         call ipPrintFrameBegin(pc,wsz(1),wsz(2))
         call draw_map(ImVec2(w%smap%szmap(1),w%smap%szmap(2)),.false.,ldum,xdum,.true.)
         dd = ipPrintFrameEnd()
      end do

      ! render it to a framebuffer and read it back
      call glGenFramebuffers(1_c_int,c_loc(fbo))
      call glGenRenderbuffers(1_c_int,c_loc(rbo))
      call glBindRenderbuffer(GL_RENDERBUFFER,rbo)
      call glRenderbufferStorage(GL_RENDERBUFFER,GL_RGBA8,wid,hei)
      call glBindRenderbuffer(GL_RENDERBUFFER,0_c_int)
      call glBindFramebuffer(GL_FRAMEBUFFER,fbo)
      call glFramebufferRenderbuffer(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_RENDERBUFFER,rbo)
      if (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE) then
         call glViewport(0_c_int,0_c_int,wid,hei)
         call glClearColor(1._c_float,1._c_float,1._c_float,1._c_float)
         call glClear(GL_COLOR_BUFFER_BIT)
         call ImGui_ImplOpenGL3_RenderDrawData(dd)
         allocate(pix(4*wid*hei))
         call glReadPixels(0_c_int,0_c_int,wid,hei,GL_RGBA,GL_UNSIGNED_BYTE,c_loc(pix))
      end if
      call glBindFramebuffer(GL_FRAMEBUFFER,0_c_int)
      call glDeleteFramebuffers(1_c_int,c_loc(fbo))
      call glDeleteRenderbuffers(1_c_int,c_loc(rbo))
      call ipPrintContextDestroy(pc)

      ! write the file (the rows of the framebuffer go bottom to top)
      if (.not.allocated(pix)) then
         w%errmsg = "Could not create the framebuffer for the image (too large?)"
         return
      end if
      pix(4::4) = char(255,c_char)
      str = trim(w%okfile) // c_null_char
      call stbi_flip_vertically_on_write(1_c_int)
      idum = stbi_write_png(c_loc(str),wid,hei,4_c_int,c_loc(pix),4_c_int*wid)
      call stbi_flip_vertically_on_write(0_c_int)
      if (idum == 0) then
         w%errmsg = "Error writing the image file: " // trim(w%okfile)
      else
         w%errmsg = ""
         write (uout,'("* Plane map image (",A,"x",A," pixels) written to: ",A/)') string(wid), string(hei),&
            trim(w%okfile)
         call okfile_save_dir(w%okfile)
      end if

    end subroutine export_map

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
      pad = 3._c_float * w%smap%vscale
      nval = 0
      if (scl <= 0._c_float) return

      do i = 1, w%smap%ncval
         str1 = fmt_level(w%smap%zcval(i)) // c_null_char
         call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
         tsz%x = tsz%x * w%smap%vscale
         tsz%y = tsz%y * w%smap%vscale
         if (w%smap%lcval(i) * scl < 3._c_float * tsz%x) cycle
         g = (0.5d0 * tsz%x + pad) / scl
         sep = tsz%y

         ! the places, from the one nearest to the phase of this level,
         ! alternating to both sides
         k0 = floor(modulo(w%smap%icval(4,i) * golden,1d0) * plane_ncand)
         do m = 0, plane_ncand-1
            l = (m + 1) / 2
            if (mod(m,2) == 1) l = -l
            k = modulo(k0 + l,plane_ncand) + 1
            sc = (k - 0.5d0) / plane_ncand * w%smap%lcval(i)
            if (sc - g < 0d0 .or. sc + g > w%smap%lcval(i)) cycle
            call ipPlotToPixels(w%smap%xcval(1,k,i),w%smap%xcval(2,k,i),ImAxis_X1,ImAxis_Y1,px,py)

            ! the direction of the line on the screen (y goes down),
            ! turned so the text reads left to right (bottom to top if
            ! vertical)
            a = real(atan2(-w%smap%tcval(2,k,i),w%smap%tcval(1,k,i)),c_float)
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
    !> cpos; 2 = the rest, in cneg), the n rows of cl, in color rgb, with
    !> a gap of half-length gap (plot units, along the line) at
    !> distance sval along the lines of the nval values placed (ival).
    subroutine plot_contours(ia,rgb,cl,n,nval,ival,sval,gap)
      use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      integer, intent(in) :: ia
      real(c_float), intent(in) :: rgb(3)
      real(c_double), intent(in), target :: cl(:,:)
      integer, intent(in) :: n, nval
      integer, intent(in) :: ival(:)
      real*8, intent(in) :: sval(:), gap(:)

      character(kind=c_char,len=:), allocatable, target :: str1
      real(c_double), allocatable, target :: dd(:,:)
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
         if (w%smap%icval(1,k) /= ia) cycle
         atrow(w%smap%icval(2,k)) = j
         ncut = ncut + 1
      end do
      if (ncut == 0) then
         call ipPlotLine(c_loc(str1),c_loc(cl(1,1)),c_loc(cl(1,2)),int(n,c_int),ImPlotLineFlags_None,0_c_int)
         return
      end if

      ! copy the lines, cutting them under their values
      nan = ieee_value(nan,ieee_quiet_nan)
      allocate(dd(n+3*ncut,2))
      nd = 0
      r = 1
      do while (r <= n)
         if (atrow(r) == 0) then
            nd = nd + 1
            dd(nd,:) = cl(r,:)
            r = r + 1
            cycle
         end if
         j = atrow(r)
         k = ival(j)
         n1 = w%smap%icval(3,k)
         sa = sval(j) - gap(j)
         sb = sval(j) + gap(j)
         nd = nd + 1
         dd(nd,:) = cl(r,:)
         s = 0d0
         istate = 0
         do i = r+1, r+n1-1
            ds = norm2(cl(i,:) - cl(i-1,:))
            s1 = s + ds
            if (istate == 0 .and. s1 >= sa) then
               nd = nd + 1
               dd(nd,:) = cl(i-1,:) + (sa - s) / max(ds,1d-40) * (cl(i,:) - cl(i-1,:))
               nd = nd + 1
               dd(nd,:) = nan
               istate = 1
            end if
            if (istate == 1 .and. s1 >= sb) then
               nd = nd + 1
               dd(nd,:) = cl(i-1,:) + (sb - s) / max(ds,1d-40) * (cl(i,:) - cl(i-1,:))
               istate = 2
            end if
            if (istate /= 1) then
               nd = nd + 1
               dd(nd,:) = cl(i,:)
            end if
            s = s1
         end do
         r = r + n1
      end do
      call ipPlotLine(c_loc(str1),c_loc(dd(1,1)),c_loc(dd(1,2)),int(nd,c_int),ImPlotLineFlags_None,0_c_int)

    end subroutine plot_contours

    !> A contour level written with three significant digits, as %.3g.
    function fmt_level(x) result(str)
      real*8, intent(in) :: x
      character(len=:), allocatable :: str

      character(len=40) :: aux
      integer :: e, idx

      if (x == 0d0) then
         str = "0"
         return
      end if
      e = floor(log10(abs(x)))
      if (e < -4 .or. e >= 3) then
         write (aux,'(ES0.2)') x
         idx = index(aux,"E")
         str = string_nozeros(x / 10d0**e,2) // "e" // trim(aux(idx+1:))
      else
         str = string_nozeros(x,max(2 - e,0))
      end if

    end function fmt_level

    !> The point of the region with plot coordinates (mu,mv): whether it
    !> is on the region (inside) and its position (x, cell-frame
    !> Cartesian bohr).
    subroutine plot_point(mu,mv,inside,x)
      real(c_double), intent(in) :: mu, mv
      logical, intent(out) :: inside
      real*8, intent(out) :: x(3)

      integer :: nx, ny
      real*8 :: dd(2), m(2,2), det, s, t

      inside = .false.
      nx = size(c%ff,1)
      ny = size(c%ff,2)
      x = c%a0 + (mu * c%e1 + mv * c%e2) / bohrtoa

      ! the grid coordinates: x - x0 = s * (nx-1) * uu + t * (ny-1) * vv
      dd(1) = dot_product(x - c%x0,c%e1)
      dd(2) = dot_product(x - c%x0,c%e2)
      m(1,1) = dot_product(c%uu,c%e1) * (nx-1)
      m(2,1) = dot_product(c%uu,c%e2) * (nx-1)
      m(1,2) = dot_product(c%vv,c%e1) * (ny-1)
      m(2,2) = dot_product(c%vv,c%e2) * (ny-1)
      det = m(1,1) * m(2,2) - m(1,2) * m(2,1)
      if (abs(det) < 1d-20) return
      s = (dd(1) * m(2,2) - dd(2) * m(1,2)) / det
      t = (m(1,1) * dd(2) - m(2,1) * dd(1)) / det
      if (s < 0d0 .or. s > 1d0 .or. t < 0d0 .or. t > 1d0) return
      inside = .true.

    end subroutine plot_point

    !> The table of the values on the plane (opened and closed from
    !> its title): the position of each point (fractional for crystals, Å
    !> for molecules), its plot coordinates, and the value. The menu
    !> of the header (right click) writes it to the output, the
    !> console, or a file.
    subroutine draw_data_section()
      use utils, only: iw_table_column, iw_table_headers_row
      integer, parameter :: itable_plane = 1 ! the table, for table_write_*
      character(kind=c_char,len=:), allocatable, target :: strtab
      integer(c_int) :: flags
      integer :: k, ix, iy, iwrite, nx, ny
      real*8 :: xp(3), xq(3), pc(2)
      type(ImVec2) :: szt
      type(c_ptr) :: clipper
      type(ImGuiListClipper), pointer :: clipper_f

      if (.not.mapshown) return
      call iw_separatortext("Values on the plane",w%editrep_slicesec(slicesec_values))
      if (.not.w%editrep_slicesec(slicesec_values)%open) return

      nx = size(c%ff,1)
      ny = size(c%ff,2)
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
         call iw_table_column("u (Å)",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("v (Å)",id=4_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
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
               xp = c%x0 + real(ix-1,8) * c%uu + real(iy-1,8) * c%vv
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
               call cell(5,string(c%ff(ix,iy),'e',decimal=8))
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

      str = quant_label(d%isexpr,d%expr,d%fref%id,d%nti)

    end function quantity_label

    !> The contour lines of the object, in plot coordinates, the
    !> positive levels in cpos and the rest in cneg (NaN between two
    !> lines), and the places along each line where its value can go,
    !> longest line first; made again if the lines of the object are
    !> not those they were made from (gencont).
    subroutine update_contours()
      use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      use tools, only: mergesort
      integer :: i, ii, j, k, nline, np
      integer, allocatable :: iord(:), il(:,:)
      real*8, allocatable :: ll(:)
      real*8 :: off(2), rval, nan, dl, acc, st
      logical :: pos

      if (.not.d%showcont .or. c%gencont == 0) then
         w%smap%ncpos = 0
         w%smap%ncneg = 0
         w%smap%ncval = 0
         w%smap%gencont = 0
         return
      end if
      if (w%smap%gencont == c%gencont) return
      w%smap%gencont = c%gencont
      w%smap%ncpos = 0
      w%smap%ncneg = 0
      w%smap%ncval = 0
      nline = c%nline
      if (nline == 0) return

      ! the lines, in the frame of the region; to plot coordinates
      off = to_plot(c%x0)
      nan = ieee_value(nan,ieee_quiet_nan)
      if (allocated(w%smap%cpos)) deallocate(w%smap%cpos)
      if (allocated(w%smap%cneg)) deallocate(w%smap%cneg)
      allocate(w%smap%cpos(size(c%xy,2)+nline,2),w%smap%cneg(size(c%xy,2)+nline,2))
      allocate(il(3,nline))
      do i = 1, nline
         pos = (c%ziso(c%liso(i)) > 0d0)
         np = merge(w%smap%ncpos,w%smap%ncneg,pos)
         il(:,i) = (/merge(1,2,pos),np+1,c%lptr(i+1)-c%lptr(i)/)
         do j = c%lptr(i), c%lptr(i+1)-1
            np = np + 1
            do k = 1, 2
               rval = c%xy(k,j) * bohrtoa + off(k)
               if (pos) then
                  w%smap%cpos(np,k) = rval
               else
                  w%smap%cneg(np,k) = rval
               end if
            end do
         end do
         np = np + 1
         if (pos) then
            w%smap%cpos(np,:) = nan
            w%smap%ncpos = np
         else
            w%smap%cneg(np,:) = nan
            w%smap%ncneg = np
         end if
      end do

      ! the length of each line, longest first, and the places along
      ! it where its value can go, with the direction of the line there
      allocate(ll(nline),iord(nline))
      do i = 1, nline
         ll(i) = 0d0
         do j = c%lptr(i)+1, c%lptr(i+1)-1
            ll(i) = ll(i) + norm2(c%xy(:,j) - c%xy(:,j-1))
         end do
         iord(i) = i
      end do
      call mergesort(-ll,iord,1,nline)
      if (allocated(w%smap%xcval)) deallocate(w%smap%xcval)
      if (allocated(w%smap%tcval)) deallocate(w%smap%tcval)
      if (allocated(w%smap%lcval)) deallocate(w%smap%lcval)
      if (allocated(w%smap%zcval)) deallocate(w%smap%zcval)
      if (allocated(w%smap%icval)) deallocate(w%smap%icval)
      allocate(w%smap%xcval(2,plane_ncand,nline),w%smap%tcval(2,plane_ncand,nline),w%smap%lcval(nline),&
         w%smap%zcval(nline),w%smap%icval(4,nline))
      w%smap%ncval = nline
      do i = 1, nline
         ii = iord(i)
         do k = 1, plane_ncand
            w%smap%xcval(:,k,i) = c%xy(:,c%lptr(ii))
            w%smap%tcval(:,k,i) = (/1d0,0d0/)
         end do
         acc = 0d0
         k = 1
         do j = c%lptr(ii)+1, c%lptr(ii+1)-1
            dl = norm2(c%xy(:,j) - c%xy(:,j-1))
            do while (k <= plane_ncand .and. dl > 0d0)
               st = (k - 0.5d0) / plane_ncand * ll(ii)
               if (acc + dl < st) exit
               w%smap%tcval(:,k,i) = (c%xy(:,j) - c%xy(:,j-1)) / dl
               w%smap%xcval(:,k,i) = c%xy(:,j-1) + (st - acc) * w%smap%tcval(:,k,i)
               k = k + 1
            end do
            acc = acc + dl
         end do
         w%smap%xcval(:,:,i) = w%smap%xcval(:,:,i) * bohrtoa + spread(off,2,plane_ncand)
         w%smap%lcval(i) = ll(ii) * bohrtoa
         w%smap%zcval(i) = c%ziso(c%liso(ii))
         w%smap%icval(1:3,i) = il(:,ii)
         w%smap%icval(4,i) = c%liso(ii)
      end do

    end subroutine update_contours

    !> The texture of the map, from the colors of the object, made
    !> again if they are not those it was made from (gentex), and the
    !> range of the colors in tmin, tmax.
    subroutine update_texture()
      use interfaces_opengl3
      integer :: k, nx, ny
      character(kind=c_char), allocatable, target :: pix(:)

      if (c%gencolor == 0 .or. .not.allocated(c%mesh%rgbv)) return
      nx = size(c%ff,1)
      ny = size(c%ff,2)
      if (size(c%mesh%rgbv,2) /= nx*ny) return
      if (w%smap%gentex == c%gencolor .and. w%smap%tex /= 0) return
      w%smap%gentex = c%gencolor
      w%smap%tmin = c%cmin
      w%smap%tmax = c%cmax

      ! the colors, as RGBA bytes (the first index of the grid runs
      ! fastest, along the width of the texture)
      allocate(pix(4*nx*ny))
      do k = 1, nx*ny
         pix(4*k-3) = char(nint(c%mesh%rgbv(1,k) * 255))
         pix(4*k-2) = char(nint(c%mesh%rgbv(2,k) * 255))
         pix(4*k-1) = char(nint(c%mesh%rgbv(3,k) * 255))
         pix(4*k) = char(255)
      end do

      ! upload
      if (w%smap%tex == 0) call glGenTextures(1,c_loc(w%smap%tex))
      call glBindTexture(GL_TEXTURE_2D,w%smap%tex)
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

      key(1) = c%gen
      key(2) = merge(1_8,0_8,w%smap%showlabels)
      key(6) = merge(1_8,0_8,w%smap%showcps)
      key(7) = w%smap%ltext
      key(8) = merge(1_8,0_8,w%smap%showdist)
      key(3) = transfer(w%smap%hmax,key(3))
      key(4) = transfer(sysc(isys)%timelastchange_cplist,key(4))
      key(5) = sys(isys)%field_cuid(sys(isys)%iref)
      if (allocated(w%smap%keymark)) then
         if (all(w%smap%keymark == key)) return
      end if
      w%smap%keymark = key
      w%smap%nmark = 0
      w%smap%marksmany = .false.
      if (allocated(w%smap%xmark)) deallocate(w%smap%xmark)
      if (allocated(w%smap%lmark)) deallocate(w%smap%lmark)
      allocate(w%smap%xmark(2,10),w%smap%lmark(10))
      if (.not.w%smap%showlabels) return

      ! the candidates: the atoms, and the critical points of the
      ! reference field except the nuclei (crystallographic)
      call label_candidates(isys,nc,xf,lbl,cps=w%smap%showcps,short=.true.,atext=w%smap%ltext)

      ! the coordinates of a point in the frame of the region: x - x0 =
      ! s1 * side1 + s2 * side2 + s3 * normal
      nx = size(c%ff,1)
      ny = size(c%ff,2)
      nrm = cross(c%e1,c%e2)
      am(:,1) = c%uu * (nx-1)
      am(:,2) = c%vv * (ny-1)
      am(:,3) = nrm
      call matinv(am,3,ier)
      if (ier /= 0) return
      hb = w%smap%hmax / bohrtoa
      eps = 1d-6

      ! the lattice translations that can bring a point of the cell to
      ! the region (crystals; only the box for molecules), if there are
      ! not too many
      lo = 0
      hi = 0
      if (.not.ismol) then
         fc(:,1) = sys(isys)%c%c2x(c%x0)
         fc(:,2) = sys(isys)%c%c2x(c%x0 + c%uu * (nx-1))
         fc(:,3) = sys(isys)%c%c2x(c%x0 + c%vv * (ny-1))
         fc(:,4) = sys(isys)%c%c2x(c%x0 + c%uu * (nx-1) + c%vv * (ny-1))
         lo = floor(minval(fc,2)) - 1
         hi = ceiling(maxval(fc,2)) + 1
         if (product(real(hi - lo + 1,8)) > maxtrans) then
            w%smap%marksmany = .true.
            return
         end if
      end if

      do k = 1, nc
         do l1 = lo(1), hi(1)
            do l2 = lo(2), hi(2)
               do l3 = lo(3), hi(3)
                  xq = sys(isys)%c%x2c(xf(:,k) + real((/l1,l2,l3/),8))
                  s = matmul(am,xq - c%x0)
                  if (s(1) < -eps .or. s(1) > 1d0+eps .or. s(2) < -eps .or. s(2) > 1d0+eps) cycle
                  if (abs(s(3)) > hb + eps) cycle
                  w%smap%nmark = w%smap%nmark + 1
                  if (w%smap%nmark > size(w%smap%lmark)) then
                     call realloc(w%smap%xmark,2,2*w%smap%nmark)
                     call realloc(w%smap%lmark,2*w%smap%nmark)
                  end if
                  w%smap%xmark(:,w%smap%nmark) = to_plot(xq)
                  w%smap%lmark(w%smap%nmark) = lbl(k)
                  ! the signed distance to the plane, as C(+0.01), unless
                  ! it is zero to two decimals
                  if (w%smap%showdist) then
                     sdist = string(s(3) * bohrtoa,'f',decimal=2)
                     if (sdist /= "0.00" .and. sdist /= "-0.00") then
                        if (s(3) > 0d0) sdist = "+" // sdist
                        w%smap%lmark(w%smap%nmark) = trim(lbl(k)) // "(" // sdist // ")"
                     end if
                  end if
               end do
            end do
         end do
      end do

    end subroutine update_marks

    !> The plot coordinates (Å) of the point x (cell-frame Cartesian
    !> bohr): along e1 and e2 from the origin the grid was calculated
    !> for (a0).
    function to_plot(x) result(pc)
      real*8, intent(in) :: x(3)
      real*8 :: pc(2)

      pc(1) = dot_product(x - c%a0,c%e1) * bohrtoa
      pc(2) = dot_product(x - c%a0,c%e2) * bohrtoa

    end function to_plot

    !> The point x (cell-frame Cartesian bohr) in the absolute frame
    !> (with molx0 for molecules), for the scene.
    function cart_to_abs(x) result(xa)
      real*8, intent(in) :: x(3)
      real*8 :: xa(3)
      xa = x
      if (ismol) xa = xa + sys(isys)%c%molx0
    end function cart_to_abs

  end function draw_editrep_slice

  !> Delete the texture of the map of the editor of a plane map
  !> object (w).
  module subroutine plane_free_texture(w)
    use interfaces_opengl3
    class(window), intent(inout), target :: w

    if (w%smap%tex /= 0) call glDeleteTextures(1,c_loc(w%smap%tex))
    w%smap%tex = 0

  end subroutine plane_free_texture

  !> Size of the image exported from the editor of a plane map whose plot and
  !> color bar take szmap (pixels of the window), at the given scale:
  !> its size in the window, with the padding of the print context
  !> (wsz), the scale used (xs, from 1 to scale, so the image fits in
  !> maxpix), and its size in pixels (npix).
  module subroutine plane_image_size(szmap,scale,wsz,xs,npix)
    real(c_float), intent(in) :: szmap(2), scale
    real(c_float), intent(out) :: wsz(2), xs
    integer(c_int), intent(out) :: npix(2)

    real(c_float), parameter :: pad = 8._c_float ! window padding of the print context (ImGui default)
    real(c_float), parameter :: spc = 8._c_float ! item spacing of the print context (ImGui default)
    real(c_float), parameter :: maxpix = 8192._c_float ! largest side of the image (pixels)

    wsz = szmap + (/spc + 2 * pad,2 * pad/)
    xs = max(1._c_float,min(scale,maxpix / maxval(wsz)))
    npix = nint(xs * wsz)

  end subroutine plane_image_size

  !> The checkbox of the relief of a plane map in the view (rl, of a
  !> rep_slice), after the previous pair in the line if gap. id makes
  !> the widget ID unique. Returns whether it changed.
  function relief_checkbox(rl,id,gap,ttshown) result(changed)
    use utils, only: iw_tooltip, iw_checkbox
    use representations, only: relief_def
    type(relief_def), intent(inout) :: rl
    character(len=*), intent(in) :: id
    logical, intent(in) :: gap
    logical, intent(inout) :: ttshown
    logical :: changed

    changed = iw_checkbox("Relief##" // id // "relief",rl%on,gap=gap,labelleft=.true.)
    call iw_tooltip("Draw the map as a relief: each point rises from the plane with its value&
       & on the scale of the relief",ttshown)

  end function relief_checkbox

  !> The relief of a plane map in the view (rl, of a rep_slice), if
  !> it is drawn: its height, after the previous pair in the line if
  !> gap (else on a new line), its side of the plane, and its scale
  !> and range. id makes the
  !> widget IDs unique. Returns whether anything changed.
  function relief_form(rl,id,gap,ttshown) result(changed)
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_dragfloat_real8, iw_combo_simple
    use representations, only: relief_def
    type(relief_def), intent(inout) :: rl
    character(len=*), intent(in) :: id
    logical, intent(in) :: gap
    logical, intent(inout) :: ttshown
    logical :: changed

    logical :: ch
    real*8 :: smin, smax
    character(len=*), parameter :: tttype = " (ctrl+click to type a value, e.g. 1e-6)"

    changed = .false.
    if (.not.rl%on) return
    call iw_text("Height (Å)",alignframe=.true.,gap=gap)
    ch = iw_dragfloat_real8("##" // id // "rheight",x1=rl%height,speed=0.01d0,min=0d0,&
       scale=bohrtoa,decimal=2,sameline=.true.)
    changed = changed .or. ch
    call iw_tooltip("Height of the relief at the top of its range",ttshown)
    ch = iw_checkbox("Flip##" // id // "rflip",rl%flip,gap=.true.,labelleft=.true.)
    changed = changed .or. ch
    call iw_tooltip("Raise the relief on the other side of the plane",ttshown)

    ! the scale and the range
    call iw_text("Relief scale",alignframe=.true.)
    call iw_combo_simple("##" // id // "rscale","Linear" // c_null_char // "Log |f|" // c_null_char //&
       "Arctangent" // c_null_char,rl%scale,sameline=.true.,changed=ch)
    changed = changed .or. ch
    call iw_tooltip("How the values become heights: linearly, by the logarithm of their&
       & absolute value, or by their arctangent",ttshown)
    ch = iw_checkbox("Fixed range##" // id // "rfrange",rl%frange,gap=.true.,labelleft=.true.)
    changed = changed .or. ch
    call iw_tooltip("Heights for a fixed range of values (or the range of the values on the map)",ttshown)
    if (.not.rl%frange) return
    ! general format, and a drag proportional to the value; a zero
    ! starts at the scale of the other bound
    smin = abs(rl%fmin)
    if (smin == 0d0) smin = abs(rl%fmax)
    smax = abs(rl%fmax)
    if (smax == 0d0) smax = abs(rl%fmin)
    call iw_text("Min",alignframe=.true.)
    ch = iw_dragfloat_real8("##" // id // "rfmin",x1=rl%fmin,speed=0.01d0*max(smin,1d-6),sig=4,&
       sameline=.true.)
    changed = changed .or. ch
    call iw_tooltip("The value at the bottom of the relief, on the plane" // tttype,ttshown)
    call iw_text("Max",alignframe=.true.,gap=.true.)
    ch = iw_dragfloat_real8("##" // id // "rfmax",x1=rl%fmax,speed=0.01d0*max(smax,1d-6),sig=4,&
       sameline=.true.)
    changed = changed .or. ch
    call iw_tooltip("The value at the top of the relief" // tttype,ttshown)

  end function relief_form

  !> The quantity of a form (the plane map and the cube windows, the
  !> plane map object editor) in system isys: quantity nti of a field
  !> fref, or an arithmetic expression expr (if isexpr; checked as it
  !> is typed, into exprerr). sfx = the prefix of the widget IDs;
  !> place = where the quantity is (for the tooltips: "the map"), and
  !> pts = whose points the expression is evaluated at ("the
  !> plane"). Returns whether it changed.
  module function quantity_form(isys,isexpr,fref,nti,expr,exprerr,sfx,place,pts,ttshown,fieldw) result(changed)
    use systems, only: sys
    use utils, only: iw_text, iw_tooltip, iw_calcwidth, iw_combo_simple, iw_field_combo, iw_inputtext,&
       iw_expression_helper, iw_exprctx_point
    integer, intent(in) :: isys
    logical, intent(inout) :: isexpr
    type(field_ref), intent(inout) :: fref
    integer, intent(inout) :: nti
    character(len=*), intent(inout) :: expr, exprerr
    character(len=*), intent(in) :: sfx, place, pts
    logical, intent(inout) :: ttshown
    integer, intent(in), optional :: fieldw
    logical :: changed

    integer :: ikind, ifield, fw
    logical :: echanged
    character(len=:), allocatable :: lerr

    changed = .false.

    call iw_text("Kind",alignframe=.true.)
    ikind = merge(1,0,isexpr)
    call iw_combo_simple("##" // sfx // "kind","Field" // c_null_char // "Expression" // c_null_char,ikind,&
       sameline=.true.)
    call iw_tooltip("A quantity of a field, or an arithmetic expression of the fields",ttshown)
    changed = changed .or. (isexpr .neqv. (ikind == 1))
    isexpr = (ikind == 1)
    if (.not.isexpr) then
       ! the field and its quantity, in the same line as the kind
       ifield = fref%id
       ! the field combo is as wide as the window allows (fieldw characters)
       fw = 24
       if (present(fieldw)) fw = fieldw
       if (iw_field_combo("##" // sfx // "field",isys,ifield,width=iw_calcwidth(fw,1),&
          nonestr="<field not available>",sameline=.true.)) then
          call fref%set(sys(isys),ifield)
          changed = .true.
       end if
       call iw_tooltip("Field on " // place,ttshown)
       if (quantity_combo("##" // sfx // "quant",isys,fref%id,nti,ttshown,sameline=.true.)) changed = .true.
    else
       ! on its own line, with room for the Insert button
       call iw_text("Expression",alignframe=.true.)
       echanged = iw_inputtext("##" // sfx // "expr",bufsize=1023,textf=expr,width=36,sameline=.true.)
       call iw_tooltip("Arithmetic expression of the fields of the system (for instance,&
          & ""$1 - $2"" or ""$1:l""), evaluated at the points of " // pts,ttshown)
       if (iw_expression_helper("##" // sfx // "help",isys,iw_exprctx_point,textf=expr,&
          ttshown=ttshown,ifield=sys(isys)%iref)) echanged = .true.
       if (echanged) then
          changed = .true.
          exprerr = ""
          if (len_trim(expr) > 0) then
             call sys(isys)%check_expression(trim(expr),lerr)
             exprerr = lerr
          end if
       end if
       if (len_trim(expr) == 0) then
          call iw_text("The expression is empty",danger=.true.,wrap=.true.)
       elseif (len_trim(exprerr) > 0) then
          call iw_text(trim(exprerr),danger=.true.,wrap=.true.)
       end if
    end if

  end function quantity_form

  !> A combo for quantity nti of field ifield of system isys: its value
  !> and derivatives, and the chemical functions the field can
  !> evaluate (elf, lol,...), with their descriptions as tooltips. A
  !> chemical function the field cannot evaluate goes back to the
  !> value; while the field is not available, the choice is kept.
  !> strid = the widget ID. Returns whether nti changed.
  module function quantity_combo(strid,isys,ifield,nti,ttshown,sameline) result(changed)
    use systems, only: sys
    use arithmetic, only: arith_names
    use types, only: field_evaluation_avail
    use utils, only: iw_combo_simple, chemfun_avail
    character(len=*), intent(in) :: strid
    integer, intent(in) :: isys, ifield
    integer, intent(inout) :: nti
    logical, intent(inout) :: ttshown
    logical, intent(in), optional :: sameline
    logical :: changed

    integer :: j, n, ic, nti0
    integer :: map(0:nquant+size(arith_names,1)-1)
    logical :: good, ch
    character(len=:), allocatable :: str, tips, why
    type(field_evaluation_avail) :: av

    ! the value and the derivatives (nti = 0 to nquant-1), and their
    ! descriptions
    character(len=*,kind=c_char), parameter :: quant_combo = &
       "Value (f)" // c_null_char // "Gradient, x" // c_null_char // "Gradient, y" // c_null_char //&
       "Gradient, z" // c_null_char // "Gradient norm" // c_null_char // "Hessian, xx" // c_null_char //&
       "Hessian, xy" // c_null_char // "Hessian, xz" // c_null_char // "Hessian, yy" // c_null_char //&
       "Hessian, yz" // c_null_char // "Hessian, zz" // c_null_char // "Laplacian" // c_null_char
    character(len=*,kind=c_char), parameter :: quant_tips = &
       "Value of the field" // c_null_char // "x component of the gradient" // c_null_char //&
       "y component of the gradient" // c_null_char // "z component of the gradient" // c_null_char //&
       "Norm of the gradient" // c_null_char // "xx component of the Hessian" // c_null_char //&
       "xy component of the Hessian" // c_null_char // "xz component of the Hessian" // c_null_char //&
       "yy component of the Hessian" // c_null_char // "yz component of the Hessian" // c_null_char //&
       "zz component of the Hessian" // c_null_char // "Laplacian of the field" // c_null_char

    ! the entries: the properties, then the chemical functions
    nti0 = nti
    str = quant_combo
    tips = quant_tips
    do j = 0, nquant-1
       map(j) = j
    end do
    n = nquant
    good = sys(isys)%goodfield(ifield)
    if (good) then
       call sys(isys)%f(ifield)%eval_avail(av)
       do j = 1, size(arith_names,1)
          if (quant_chem(nquant+j) == 0) cycle
          if (.not.chemfun_avail(arith_names(j),sys(isys)%f(ifield),av,why)) cycle
          str = str // trim(arith_names(j)%name) // c_null_char
          tips = tips // trim(arith_names(j)%desc) // " (" // trim(arith_names(j)%name) // ")" // c_null_char
          map(n) = nquant + j
          n = n + 1
       end do
    end if

    ! the combo
    ic = 0
    do j = 0, n-1
       if (map(j) == nti) ic = j
    end do
    call iw_combo_simple(strid,str,ic,sameline=sameline,changed=ch,tooltips=tips,ttshown=ttshown)
    if (good .or. ch) nti = map(max(min(ic,n-1),0))
    changed = (nti /= nti0)

  end function quantity_combo

end submodule plane
