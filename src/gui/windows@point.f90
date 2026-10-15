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

! Routines for the point properties window (the POINT keyword).
submodule (windows) point
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! transient marker at the point: color, opacity, and radius (bohr),
  ! a fraction of the scene radius between the two limits
  real(c_float), parameter :: marker_rgb(3) = (/1.0_c_float,0.6_c_float,0.1_c_float/)
  real(c_float), parameter :: marker_alpha = 0.8_c_float
  real*8, parameter :: marker_radmin = 0.25d0
  real*8, parameter :: marker_radmax = 0.6d0

contains

  !> Draw the point properties window: the properties of a field of
  !> the system shown by the anchor view at a point, which is typed or
  !> picked in the view (an atom, a critical point, or empty space), as
  !> POINT calculates them.
  module subroutine draw_point(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_text, iw_button, iw_tooltip, iw_checkbox, iw_dragfloat_real8,&
       iw_field_combo, iw_inputtext, iw_calcwidth, iw_close_event, iw_setpos_bottomright,&
       iw_table_column, iw_table_headers_row, iw_cell_right
    use representations, only: rep_shape, shapekind_sphere
    use systemmod, only: propitem_header, propitem_row, propitem_text
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer, parameter :: itable_point = 1 ! the table, for table_write_*

    logical :: doquit, goodsys, syschanged, ismol, ldum
    integer :: i, isys, iview, istat, iwrite, ifield
    integer(c_int) :: tflags
    real*8 :: xc(3), rad
    type(ImVec2) :: sz0
    character(len=:,kind=c_char), allocatable, target :: str1

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)
    ismol = .false.
    if (goodsys) ismol = sys(isys)%c%ismolecule

    ! initialize the form on the first pass
    if (w%firstpass) then
       w%pt = point_state()
       w%pt%expr = ""
       w%errmsg = ""
    end if

    ! a new system: the point goes to the center of the cell (of the
    ! molecular box) and the field to the reference
    if (goodsys) then
       if (.not.w%pt%sys%is(isys) .or. syschanged) then
          call w%pt%sys%set(isys)
          call cancel_pick()
          w%pt%x = form_coords(isys,sys(isys)%c%x2c((/0.5d0,0.5d0,0.5d0/)))
          w%pt%fref = field_ref()
       end if
       ! the reference if no field is chosen, or if it is gone
       if (.not.w%pt%fref%ok(sys(isys))) call w%pt%fref%set(sys(isys),sys(isys)%iref)
    end if

    ! a pending pick in the view it was armed on; dropped, releasing the
    ! view, if the system is not ready or the window moved to another view
    if (w%pt%picking) then
       if (goodsys .and. iview == w%pt%pickview) then
          call view_pick_result(iview,w%id,isys,w%pt%pick,istat,xc)
          if (istat == ipick_point) w%pt%x = form_coords(isys,xc)
          if (istat /= ipick_pending) w%pt%picking = .false.
       else
          call cancel_pick()
       end if
    end if

    if (goodsys) then
       ! the system
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! the field and whether the others are reported too
       call iw_text("Field",highlight=.true.,alignframe=.true.)
       call igSameLine(0._c_float,-1._c_float)
       ifield = w%pt%fref%id
       if (iw_field_combo("##pointfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
          nonestr="<field not available>")) call w%pt%fref%set(sys(isys),ifield)
       call iw_tooltip("Field whose properties are calculated at the point",ttshown)
       ldum = iw_checkbox("All fields##pointall",w%pt%allfields,sameline=.true.)
       call iw_tooltip("Also give the value, gradient norm, and Laplacian of all other fields&
          & at the point (ALL)",ttshown)

       ! the point
       call iw_text("Point",highlight=.true.,alignframe=.true.)
       if (iw_button("Pick##pointpick",sameline=.true.,disabled=w%pt%picking)) then
          w%pt%picking = .true.
          w%pt%pickview = iview
          call w%pt%pick%arm()
          call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick an atom, a critical point,&
             & or a point in space",w%id,acceptempty=.true.)
       end if
       call iw_tooltip("Pick the point in the view: an atom, a critical point,&
          & or a point in empty space",ttshown)
       if (ismol) then
          ldum = iw_dragfloat_real8("(Å)##pointx",x3=w%pt%x,speed=0.01d0,decimal=5,sameline=.true.)
       else
          ldum = iw_dragfloat_real8("(frac)##pointx",x3=w%pt%x,speed=0.001d0,decimal=6,sameline=.true.)
       end if
       call iw_tooltip("Coordinates of the point",ttshown)

       ! an expression evaluated at the point
       call iw_text("Expression",highlight=.true.,alignframe=.true.)
       call igSameLine(0._c_float,-1._c_float)
       ldum = iw_inputtext("##pointexpr",bufsize=2047,texta=w%pt%expr,width=40,notlive=.true.)
       call iw_tooltip("Arithmetic expression to evaluate at the point, in addition to the&
          & properties of the field (e.g. $1+$2, or $1:l for the Laplacian of field 1)",ttshown)

       ! the point in the anchor view, for as long as the window is open
       if (associated(win(iview)%sc)) then
          if (win(iview)%sc%isinit /= 0) then
             rad = min(max(0.01d0 * real(win(iview)%sc%scenerad,8),marker_radmin),marker_radmax)
             call win(iview)%sc%show_transient_shapes(w%uid,1,(/rep_shape(kind=shapekind_sphere,&
                x1=to_abs(),rad=rad,rgb=marker_rgb,alpha=marker_alpha)/))
          end if
       end if

       ! the report, calculated again if anything it depends on changed
       call update_report()

       ! the report table: the position, then the properties
       call iw_text("Properties",highlight=.true.)
       tflags = ImGuiTableFlags_None
       tflags = ior(tflags,ImGuiTableFlags_NoSavedSettings)
       tflags = ior(tflags,ImGuiTableFlags_RowBg)
       tflags = ior(tflags,ImGuiTableFlags_Borders)
       tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
       str1 = "##pointtable" // c_null_char
       sz0%x = 0._c_float
       sz0%y = 0._c_float
       if (igBeginTable(c_loc(str1),4,tflags,sz0,0._c_float)) then
          call iw_table_column("Property",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Value",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_headers_row(writemenu=iwrite,ttshown=ttshown)
          call w%table_write_begin(itable_point,iwrite)

          ! the Cartesian coordinates in the absolute frame (molecules
          ! with their origin, as POINT reads and writes them)
          xc = to_abs()
          if (.not.ismol) &
             call row("Coordinates (frac)",vec3(w%pt%x,'f',7))
          call row("Coordinates (bohr)",vec3(xc,'f',7))
          call row("Coordinates (Å)",vec3(xc*bohrtoa,'f',7))
          if (sys(isys)%goodfield(w%pt%fref%id)) &
             call row("Field",string(w%pt%fref%id) // ": " // trim(sys(isys)%f(w%pt%fref%id)%name),&
                whole=.true.)
          do i = 1, w%pt%nitem
             associate(it => w%pt%item(i))
               select case (it%kind)
               case (propitem_header,propitem_text)
                  call row(it%label,"")
               case (propitem_row)
                  call row("",it%value)
               case default
                  call row(it%label,it%value)
               end select
             end associate
          end do
          if (allocated(w%pt%exprval)) then
             if (len(w%pt%exprval) > 0) &
                call row("Expression (" // trim(w%pt%expr) // ")",w%pt%exprval)
          end if

          call w%table_write_end(itable_point,table_title("Point properties",isys))
          call igEndTable()
       end if
       if (.not.sys(isys)%goodfield(w%pt%fref%id)) &
          call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       if (allocated(w%pt%errmsg)) then
          if (len_trim(w%pt%errmsg) > 0) call iw_text(w%pt%errmsg,danger=.true.,wrap=.true.)
       end if
    end if

    ! bottom row: the equivalent command and close
    call iw_setpos_bottomright(len("Copy Command") + len("Close"),2)
    if (iw_button("Copy Command",disabled=.not.goodsys)) then
       str1 = point_command() // c_null_char
       call igSetClipboardText(c_loc(str1))
    end if
    if (goodsys) call iw_tooltip("Copy to the clipboard the POINT command that calculates&
       & these properties:" // new_line('a') // point_command(),ttshown)
    if (iw_button("Close",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) then
       call cancel_pick()
       call w%end()
    end if

  contains

    !> A row of the report table: the label, then the value. A value
    !> with several numbers (separated by blanks) takes one column for
    !> each, right-aligned so that the rows of a tensor line up. If
    !> whole, the value is a single piece of text instead.
    subroutine row(label,value,whole)
      use tools_io, only: getword
      character(len=*), intent(in) :: label, value
      logical, intent(in), optional :: whole

      integer :: k, lp
      character(len=:), allocatable :: word

      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
      if (igTableSetColumnIndex(0_c_int)) call iw_text(label)
      if (present(whole)) then
         if (whole) then
            if (igTableSetColumnIndex(1_c_int)) call iw_text(value)
            return
         end if
      end if
      lp = 1
      do k = 1, 3
         word = getword(value,lp)
         if (len(word) == 0) exit
         if (igTableSetColumnIndex(int(k,c_int))) call iw_cell_right(word)
      end do

    end subroutine row

    !> Calculate the report at the point again if the inputs it was
    !> made from changed: the system, its geometry and fields (the
    !> other fields enter through ALL and the expression), the
    !> point properties, the options, and the point. Waits for the
    !> initialization threads to finish (the last report stays).
    subroutine update_report()
      use systems, only: are_threads_running
      use types, only: scalar_value
      integer*8, allocatable :: key(:)
      logical :: stale
      integer :: j
      real*8 :: x0(3), rdum
      type(scalar_value) :: res
      character(len=:), allocatable :: lerr, keyexpr

      ! the expression and the point properties, as text
      keyexpr = w%pt%expr
      do j = 1, sys(isys)%npropp
         keyexpr = keyexpr // c_null_char // trim(sys(isys)%propp(j)%name) // c_null_char //&
            trim(sys(isys)%propp(j)%expr) // c_null_char // string(sys(isys)%propp(j)%ispecial)
      end do

      ! the number of fields enters through the size of key
      allocate(key(5+sys(isys)%nf+1))
      key(1) = w%pt%sys%uid
      key(2) = w%pt%fref%uid
      key(3) = merge(1_8,0_8,w%pt%allfields)
      key(4) = sys(isys)%npropp
      key(5) = w%pt%fref%id
      do j = 0, sys(isys)%nf
         key(6+j) = sys(isys)%field_cuid(j)
      end do

      stale = .not.allocated(w%pt%key) .or. .not.allocated(w%pt%keyexpr)
      if (.not.stale) stale = (size(w%pt%key) /= size(key))
      if (.not.stale) stale = any(w%pt%key /= key) .or. any(w%pt%keyx /= w%pt%x) .or.&
         (w%pt%keytime /= sysc(isys)%timelastchange_geometry) .or. (w%pt%keyexpr /= keyexpr)
      if (.not.stale .or. are_threads_running()) return

      w%pt%nitem = 0
      w%pt%errmsg = ""
      w%pt%exprval = ""
      x0 = sys(isys)%c%c2x(form_to_cart(isys,w%pt%x))
      if (sys(isys)%goodfield(w%pt%fref%id)) &
         call sys(isys)%propty_items(w%pt%fref%id,x0,res,.false.,w%pt%allfields,&
            w%pt%item,w%pt%nitem,w%pt%errmsg)
      if (len_trim(w%pt%expr) > 0) then
         lerr = ""
         rdum = sys(isys)%eval(trim(w%pt%expr),lerr,form_to_cart(isys,w%pt%x))
         if (len_trim(lerr) > 0) then
            if (len_trim(w%pt%errmsg) > 0) w%pt%errmsg = w%pt%errmsg // new_line('a')
            w%pt%errmsg = w%pt%errmsg // "Error in the expression: " // trim(lerr)
         else
            w%pt%exprval = string(rdum,'e',decimal=9)
         end if
      end if

      call move_alloc(key,w%pt%key)
      w%pt%keyx = w%pt%x
      w%pt%keytime = sysc(isys)%timelastchange_geometry
      w%pt%keyexpr = keyexpr

    end subroutine update_report

    !> The POINT command(s) that calculate the report: the field (ALL)
    !> and the expression, if any. Molecular coordinates in the units
    !> POINT reads them in a script: Å, the default after MOLECULE,
    !> unless UNITS set them.
    function point_command() result(str)
      use global, only: iunit, iunit_isdef, iunit_bohr
      character(len=:), allocatable :: str

      character(len=:), allocatable :: sx

      if (ismol .and. .not.iunit_isdef .and. iunit == iunit_bohr) then
         sx = vec3(w%pt%x / bohrtoa,'f',10)
      else
         sx = vec3(w%pt%x,'f',10)
      end if
      str = "POINT " // sx // " FIELD " // string(w%pt%fref%id)
      if (w%pt%allfields) str = str // " ALL"
      if (len_trim(w%pt%expr) > 0) &
         str = str // new_line('a') // "POINT " // sx // ' FIELD "' // trim(w%pt%expr) // '"'

    end function point_command

    !> Three numbers separated by spaces.
    function vec3(x,fmt,ndec) result(str)
      real*8, intent(in) :: x(3)
      character*1, intent(in) :: fmt
      integer, intent(in) :: ndec
      character(len=:), allocatable :: str

      str = string(x(1),fmt,decimal=ndec) // " " // string(x(2),fmt,decimal=ndec) //&
         " " // string(x(3),fmt,decimal=ndec)

    end function vec3

    !> The point in the absolute frame (with molx0 for molecules), for
    !> the scene.
    function to_abs() result(x)
      real*8 :: x(3)
      x = form_to_cart(isys,w%pt%x)
      if (ismol) x = x + sys(isys)%c%molx0
    end function to_abs

    !> Forget a pending pick and release the view it was armed on.
    subroutine cancel_pick()
      if (.not.w%pt%picking) return
      if (w%pt%pickview >= 1 .and. w%pt%pickview <= nwin) &
         call win(w%pt%pickview)%viewmode_release_forced(w%id)
      w%pt%picking = .false.
    end subroutine cancel_pick

  end subroutine draw_point

end submodule point
