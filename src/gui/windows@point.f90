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

  ! name and definition (its tooltip) of the built-in properties
  ! (propid_*), for the list of properties calculated at the point
  character(len=*), parameter :: propname(propid_NUM) = (/&
     "Type                         ","Field value                  ",&
     "Valence field value          ","Gradient                     ",&
     "Gradient norm                ","Valence gradient norm        ",&
     "Kinetic energy density       ","Laplacian                    ",&
     "Valence Laplacian            ","Hessian                      ",&
     "Hessian eigenvalues          ","Ellipticity                  ",&
     "Spin field values            ","Spin kinetic energy densities",&
     "(point property)             ","Other fields                 "/)
  character(len=*), parameter :: propdef(propid_NUM) = (/&
     "rank and signature (r,s), or nucleus","f                                   ",&
     "fval                                ","grad f                              ",&
     "|grad f|                            ","|grad fval|                         ",&
     "G (tau)                             ","del2 f                              ",&
     "del2 fval                           ","d2f/dxi dxj                         ",&
     "l1, l2, l3                          ","l1/l2 - 1 (at bond-like points)     ",&
     "up and down: f, |grad f|, del2 f    ","up and down: G                      ",&
     "                                    ","f, |grad f|, del2 f of each (ALL)   "/)

contains

  !> Draw the point properties window: the properties of a field of
  !> the system shown by the anchor view at a point, which is typed or
  !> picked in the view (an atom, a critical point, or empty space), as
  !> POINT calculates them.
  module subroutine draw_point(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8,&
       iw_field_combo, iw_calcwidth, iw_close_event, iw_setpos_bottomright,&
       iw_table_column, iw_table_headers_row, iw_calcheight
    use representations, only: rep_shape, shapekind_sphere
    use systemmod, only: propitem_header, propitem_row, propitem_text, propid_none, propid_type,&
       propid_fval, propid_gfmodval, propid_lapval, propid_ellip, propid_pointprop, propid_field
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer, parameter :: itable_point = 1 ! the table, for table_write_*

    logical :: doquit, goodsys, syschanged, ismol, ldum, xediting
    integer :: i, isys, iview, istat, iwrite, ifield
    type(point_col) :: cc
    character(len=:), allocatable :: hdr
    real(c_float) :: wnew(3)
    integer(c_int) :: tflags
    real*8 :: xc(3), xrep(3), rad
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
       w%errmsg = ""
       ! hidden by default: the type and valence quantities, and the
       ! other fields
       w%pt%col(propid_type)%show = .false.
       w%pt%col(propid_fval)%show = .false.
       w%pt%col(propid_gfmodval)%show = .false.
       w%pt%col(propid_lapval)%show = .false.
       w%pt%col(propid_field)%show = .false.
       ! the ellipticity is of order one: fixed point
       w%pt%col(propid_ellip)%expo = .false.
       allocate(w%pt%ppcol(0))
    end if

    ! a new system: the point goes to the center of the cell (of the
    ! molecular box) and the field to the reference
    if (goodsys) then
       if (.not.w%pt%sys%is(isys) .or. syschanged) then
          call w%pt%sys%set(isys)
          call cancel_pick()
          w%errmsg = ""
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

       ! the field
       call iw_text("Field",highlight=.true.,alignframe=.true.)
       call igSameLine(0._c_float,-1._c_float)
       ifield = w%pt%fref%id
       if (iw_field_combo("##pointfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
          nonestr="<field not available>")) call w%pt%fref%set(sys(isys),ifield)
       call iw_tooltip("Field whose properties are calculated at the point",ttshown)

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
       ! the report waits until the coordinates are no longer being edited
       xediting = igIsItemActive()
       call iw_tooltip("Coordinates of the point",ttshown)

       ! which properties are calculated, and how they are written
       call draw_props_section()

       ! the point in the anchor view, for as long as the window is open
       if (associated(win(iview)%sc)) then
          if (win(iview)%sc%isinit /= 0) then
             rad = min(max(0.01d0 * real(win(iview)%sc%scenerad,8),marker_radmin),marker_radmax)
             call win(iview)%sc%show_transient_shapes(w%uid,1,(/rep_shape(kind=shapekind_sphere,&
                x1=to_abs(w%pt%x),rad=rad,rgb=marker_rgb,alpha=marker_alpha)/))
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
       tflags = ior(tflags,ImGuiTableFlags_ScrollY)
       str1 = "##pointtable" // c_null_char
       ! the rest of the window, as the results table of the critical
       ! points window, leaving room for a message and Close; but at
       ! least a few rows (the stretch beyond them is elastic)
       call igGetContentRegionAvail(sz0)
       sz0%x = 0._c_float
       sz0%y = sz0%y - iw_calcheight(3,0,.true.)
       sz0%y = max(sz0%y,iw_calcheight(7,0,.false.))
       w%heightslack = sz0%y - iw_calcheight(7,0,.false.)
       if (igBeginTable(c_loc(str1),4,tflags,sz0,0._c_float)) then
          call iw_table_column("Property",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Value",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_headers_row(writemenu=iwrite,ttshown=ttshown,freezetop=.true.)
          call w%table_write_begin(itable_point,iwrite)
          wnew = 0._c_float

          ! the point of the report, with the Cartesian coordinates in
          ! the absolute frame (molecules with their origin, as POINT
          ! reads and writes them)
          xrep = w%pt%x
          if (allocated(w%pt%key)) xrep = w%pt%keyx
          xc = to_abs(xrep)
          if (.not.ismol) &
             call row("Coordinates (frac)",vec3(xrep,'f',7))
          call row("Coordinates (bohr)",vec3(xc,'f',7))
          call row("Coordinates (Å)",vec3(xc*bohrtoa,'f',7))
          ! the label of a tensor goes on the row of its first line
          hdr = ""
          do i = 1, w%pt%nitem
             associate(it => w%pt%item(i))
               if (.not.item_col(it,cc)) cycle
               select case (it%kind)
               case (propitem_header)
                  hdr = it%label
               case (propitem_text)
                  call row(it%label,"")
               case (propitem_row)
                  call row(hdr,nums(it%x,cc))
                  hdr = ""
               case default
                  if (size(it%x,1) > 0) then
                     call row(it%label,nums(it%x,cc))
                  else
                     call row(it%label,it%value)
                  end if
               end select
             end associate
          end do

          w%pt%wcol = wnew
          call w%table_write_end(itable_point,table_title("Point properties",isys))
          call igEndTable()
       end if
       if (.not.sys(isys)%goodfield(w%pt%fref%id)) &
          call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       if (allocated(w%pt%errmsg)) then
          if (len_trim(w%pt%errmsg) > 0) call iw_text(w%pt%errmsg,danger=.true.,wrap=.true.)
       end if
       if (allocated(w%errmsg)) then
          if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)
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

    !> The collapsible list of the properties calculated at the point:
    !> the built-in ones of the field (those in the report), the point
    !> properties of the system, and the other fields, each with whether
    !> it is shown and the notation and decimal places of its numbers in
    !> the table; and the form to add point properties (as POINTPROP).
    subroutine draw_props_section()
      character(kind=c_char,len=:), allocatable, target :: strsec, strtab
      integer(c_int) :: flt
      integer :: jp, ip, nrow, idel
      logical :: other, del, ldum
      type(ImVec2) :: szt

      ! a framed header across the window, so it stands out
      strsec = "Properties calculated at the point##pointprops" // c_null_char
      if (.not.igTreeNodeEx_Str(c_loc(strsec),ImGuiTreeNodeFlags_Framed)) return

      ! whether there are other fields than the one of the window
      other = .false.
      do jp = 0, sys(isys)%nf
         if (jp /= w%pt%fref%id .and. sys(isys)%goodfield(jp)) other = .true.
      end do

      ! the table of properties: show, name, definition, notation and
      ! decimal places of the numbers in the table, delete
      nrow = sys(isys)%npropp
      if (other) nrow = nrow + 1
      do jp = 1, propid_NUM
         if (jp /= propid_pointprop .and. jp /= propid_field .and. in_report(jp)) nrow = nrow + 1
      end do
      flt = ImGuiTableFlags_None
      flt = ior(flt,ImGuiTableFlags_NoSavedSettings)
      flt = ior(flt,ImGuiTableFlags_RowBg)
      flt = ior(flt,ImGuiTableFlags_Borders)
      flt = ior(flt,ImGuiTableFlags_SizingFixedFit)
      flt = ior(flt,ImGuiTableFlags_ScrollY)
      strtab = "##tablepointprops" // c_null_char
      szt%x = 0._c_float
      szt%y = iw_calcheight(min(nrow,5)+1,0,.false.)
      idel = 0
      if (igBeginTable(c_loc(strtab),5_c_int,flt,szt,0._c_float)) then
         call iw_table_column("Show",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Property",id=1_c_int,flags=ImGuiTableColumnFlags_WidthStretch)
         call iw_table_column("Format",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Decimals",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("(delete)",id=4_c_int,flags=ior(ImGuiTableColumnFlags_WidthFixed,&
            ImGuiTableColumnFlags_NoHeaderLabel))
         call iw_table_headers_row(freezetop=.true.)

         ! the built-in properties in the report
         do jp = 1, propid_NUM
            if (jp == propid_pointprop .or. jp == propid_field) cycle
            if (.not.in_report(jp)) cycle
            call prop_row("b" // string(jp),trim(propname(jp)),trim(propdef(jp)),w%pt%col(jp),&
               numbers=(jp /= propid_type))
         end do

         ! the point properties of the system
         do ip = 1, sys(isys)%npropp
            associate(pp => sys(isys)%propp(ip))
              jp = ppcol_index(pp%name)
              if (pp%ispecial == 0) then
                 call prop_row("p" // string(ip),trim(pp%name) // " (" // trim(pp%expr) // ")",&
                    "Point property " // trim(pp%name) // ": " // trim(pp%expr),w%pt%ppcol(jp),del=del)
              else
                 call prop_row("p" // string(ip),trim(pp%name) // " (stress tensor)",&
                    "Point property " // trim(pp%name) // ": the Schrödinger stress tensor",&
                    w%pt%ppcol(jp),del=del)
              end if
              if (del) idel = ip
            end associate
         end do

         ! the other fields (ALL)
         if (other) call prop_row("f",trim(propname(propid_field)),trim(propdef(propid_field)),&
            w%pt%col(propid_field),calc=.true.)

         call igEndTable()
      end if
      if (idel > 0) then
         ! forget its settings, so a new property with its name starts afresh
         jp = ppcol_find(sys(isys)%propp(idel)%name)
         if (jp > 0) w%pt%ppcol = (/w%pt%ppcol(1:jp-1),w%pt%ppcol(jp+1:)/)
         call sys(isys)%delete_pointprop(idel)
      end if

      ! add a point property: an expression, or one of the keywords
      ! of POINTPROP, for the field of the window
      ldum = pointprop_add_form(isys,w%pt%fref%id,w%pt%ppkind,w%pt%ppname,w%pt%ppexpr,"point",&
         "the point",ttshown,w%errmsg)

      call igTreePop()

    end subroutine draw_props_section

    !> A row of the list of properties calculated at the point: label,
    !> its definition (the tooltip of the label), and settings c. Without numbers, no notation and
    !> decimals. If calc, Show also decides whether it is calculated
    !> (the others are always calculated, and only hidden). If del is
    !> present, the row has a delete button, and del is whether it was
    !> clicked.
    subroutine prop_row(suffix,label,def,c,numbers,calc,del)
      use gui_main, only: g
      use utils, only: iw_checkbox, iw_intstepper, iw_combo_simple, iw_close_button
      character(len=*), intent(in) :: suffix, label, def
      type(point_col), intent(inout) :: c
      logical, intent(in), optional :: numbers, calc
      logical, intent(out), optional :: del

      logical :: ldum, ch, numbers_, calc_
      integer :: ifmt
      character(len=:), allocatable :: suf

      numbers_ = .true.
      if (present(numbers)) numbers_ = numbers
      calc_ = .false.
      if (present(calc)) calc_ = calc
      suf = "##pointprop" // suffix

      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
      if (igTableSetColumnIndex(0_c_int)) then
         ldum = iw_checkbox(suf // "show",c%show)
         if (calc_) then
            call iw_tooltip("Calculate this property and show it in the table",ttshown)
         else
            call iw_tooltip("Show this property in the table",ttshown)
         end if
      end if
      if (igTableSetColumnIndex(1_c_int)) then
         call iw_text(label,alignframe=.true.)
         call iw_tooltip(def,ttshown)
      end if
      if (numbers_ .and. igTableSetColumnIndex(2_c_int)) then
         ifmt = merge(1,0,c%expo)
         call iw_combo_simple(suf // "fmt","Fixed" // c_null_char // "Exp" // c_null_char,&
            ifmt,changed=ch)
         if (ch) c%expo = (ifmt == 1)
         call iw_tooltip("Notation of the numbers of this property in the table: fixed&
            & point (0.0123) or exponential (1.23E-02)",ttshown)
      end if
      if (numbers_ .and. igTableSetColumnIndex(3_c_int)) then
         ldum = iw_intstepper(suf // "dec",c%ndec,minval=0_c_int,maxval=16_c_int,&
            ndigit=2,tooltip="Number of decimal places of this property in the table")
      end if
      if (present(del)) then
         del = .false.
         if (igTableSetColumnIndex(4_c_int)) then
            ! the close icon is one line of text high: centered on the row frames
            call igSetCursorPosY(igGetCursorPosY() + g%Style%FramePadding%y)
            del = iw_close_button(suf // "del")
            call iw_tooltip("Remove this point property from the system",ttshown)
         end if
      end if

    end subroutine prop_row

    !> Whether the report has an entry of the built-in property jp.
    function in_report(jp)
      integer, intent(in) :: jp
      logical :: in_report
      in_report = any(w%pt%item(1:w%pt%nitem)%iprop == jp)
    end function in_report

    !> The settings c of the property of report entry it (the defaults
    !> for a message); false if it is not shown: hidden, or a point
    !> property no longer in the system.
    function item_col(it,c) result(ok)
      type(propty_item), intent(in) :: it
      type(point_col), intent(out) :: c
      logical :: ok
      integer :: ic
      ok = .true.
      if (it%iprop == propid_none) return
      if (it%iprop == propid_pointprop) then
         ok = (it%idx >= 1 .and. it%idx <= sys(isys)%npropp)
         if (ok) then
            ! ppcol_index may grow ppcol, so not in the same statement
            ic = ppcol_index(sys(isys)%propp(it%idx)%name)
            c = w%pt%ppcol(ic)
         end if
      else
         c = w%pt%col(it%iprop)
      end if
      ok = ok .and. c%show
    end function item_col

    !> The settings of the point property with this name, created with
    !> the defaults if it has none yet.
    function ppcol_index(name) result(ic)
      character(len=*), intent(in) :: name
      integer :: ic
      type(point_col), allocatable :: aux(:)
      ic = ppcol_find(name)
      if (ic > 0) return
      allocate(aux(size(w%pt%ppcol,1)+1))
      aux(1:size(w%pt%ppcol,1)) = w%pt%ppcol
      aux(size(aux,1)) = point_col(name=name,ndec=9)
      call move_alloc(aux,w%pt%ppcol)
      ic = size(w%pt%ppcol,1)
    end function ppcol_index

    !> The settings of the point property with this name (0 if none). A
    !> loop: findloc on the name component would copy it every call.
    function ppcol_find(name) result(ic)
      character(len=*), intent(in) :: name
      integer :: ic
      do ic = 1, size(w%pt%ppcol,1)
         if (w%pt%ppcol(ic)%name == name) return
      end do
      ic = 0
    end function ppcol_find

    !> The numbers x written with the notation and decimal places of c,
    !> separated by spaces.
    function nums(x,c) result(str)
      real*8, intent(in) :: x(:)
      type(point_col), intent(in) :: c
      character(len=:), allocatable :: str
      integer :: k
      str = ""
      do k = 1, size(x,1)
         str = str // " " // string(x(k),merge('e','f',c%expo),decimal=int(c%ndec))
      end do
    end function nums

    !> A row of the report table: the label, then the value. A value
    !> with several numbers (separated by blanks) takes one column for
    !> each, right-aligned so that the rows of a tensor line up: on the
    !> widest number of the column in the last frame (w%pt%wcol), not on
    !> the column edge, so the column fits its numbers as they change.
    !> wnew collects the widths of this frame.
    subroutine row(label,value)
      use tools_io, only: getword
      use utils, only: iw_textwidth
      character(len=*), intent(in) :: label, value

      integer :: k, lp
      real(c_float) :: tw
      character(len=:), allocatable :: word

      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
      if (igTableSetColumnIndex(0_c_int)) call iw_text(label)
      lp = 1
      do k = 1, 3
         word = trim(getword(value,lp))
         if (len(word) == 0) exit
         if (igTableSetColumnIndex(int(k,c_int))) then
            tw = iw_textwidth(word)
            wnew(k) = max(wnew(k),tw)
            call igSetCursorPosX(igGetCursorPosX() + max(w%pt%wcol(k) - tw,0._c_float))
            call iw_text(word)
         end if
      end do

    end subroutine row

    !> Calculate the report at the point again if the inputs it was
    !> made from changed: the system, its geometry and fields (the
    !> other fields enter through ALL and the point properties), the
    !> point properties, the options, and the point. Waits for the
    !> initialization threads to finish and for the coordinates to be
    !> no longer edited (the last report stays).
    subroutine update_report()
      use systems, only: are_threads_running
      use types, only: scalar_value
      integer*8, allocatable :: key(:)
      logical :: stale
      integer :: j
      real*8 :: x0(3)
      type(scalar_value) :: res
      character(len=:), allocatable :: keyexpr

      ! the point properties, as text
      keyexpr = ""
      do j = 1, sys(isys)%npropp
         keyexpr = keyexpr // c_null_char // trim(sys(isys)%propp(j)%name) // c_null_char //&
            trim(sys(isys)%propp(j)%expr) // c_null_char // string(sys(isys)%propp(j)%ispecial)
      end do

      ! the number of fields enters through the size of key
      allocate(key(5+sys(isys)%nf+1))
      key(1) = w%pt%sys%uid
      key(2) = w%pt%fref%uid
      key(3) = merge(1_8,0_8,w%pt%col(propid_field)%show)
      key(4) = sys(isys)%npropp
      key(5) = w%pt%fref%id
      do j = 0, sys(isys)%nf
         key(6+j) = sys(isys)%field_cuid(j)
      end do

      stale = .not.allocated(w%pt%key) .or. .not.allocated(w%pt%keyexpr)
      if (.not.stale) stale = (size(w%pt%key) /= size(key))
      if (.not.stale) stale = any(w%pt%key /= key) .or. any(w%pt%keyx /= w%pt%x) .or.&
         (w%pt%keytime /= sysc(isys)%timelastchange_geometry) .or. (w%pt%keyexpr /= keyexpr)
      if (.not.stale .or. xediting .or. are_threads_running()) return

      w%pt%nitem = 0
      w%pt%errmsg = ""
      x0 = sys(isys)%c%c2x(form_to_cart(isys,w%pt%x))
      if (sys(isys)%goodfield(w%pt%fref%id)) &
         call sys(isys)%propty_items(w%pt%fref%id,x0,res,.false.,w%pt%col(propid_field)%show,&
            w%pt%item,w%pt%nitem,w%pt%errmsg)

      call move_alloc(key,w%pt%key)
      w%pt%keyx = w%pt%x
      w%pt%keytime = sysc(isys)%timelastchange_geometry
      w%pt%keyexpr = keyexpr

    end subroutine update_report

    !> The commands that calculate the report: the point properties of
    !> the system (POINTPROP), then POINT for the field (ALL). Molecular
    !> coordinates in the units POINT reads them in a script: Å, the
    !> default after MOLECULE, unless UNITS set them.
    function point_command() result(str)
      use global, only: iunit, iunit_isdef, iunit_bohr
      character(len=:), allocatable :: str

      character(len=:), allocatable :: sx
      integer :: j

      if (ismol .and. .not.iunit_isdef .and. iunit == iunit_bohr) then
         sx = vec3(w%pt%x / bohrtoa,'f',10)
      else
         sx = vec3(w%pt%x,'f',10)
      end if
      str = ""
      do j = 1, sys(isys)%npropp
         associate(pp => sys(isys)%propp(j))
           if (pp%ispecial == 0) then
              str = str // "POINTPROP " // trim(pp%name) // ' "' // trim(pp%expr) // '"' // new_line('a')
           else
              str = str // "POINTPROP stress" // new_line('a')
           end if
         end associate
      end do
      str = str // "POINT " // sx // " FIELD " // string(w%pt%fref%id)
      if (w%pt%col(propid_field)%show) str = str // " ALL"

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
      if (.not.w%pt%picking) return
      if (w%pt%pickview >= 1 .and. w%pt%pickview <= nwin) &
         call win(w%pt%pickview)%viewmode_release_forced(w%id)
      w%pt%picking = .false.
    end subroutine cancel_pick

  end subroutine draw_point

end submodule point
