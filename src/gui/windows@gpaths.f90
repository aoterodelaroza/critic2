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

! Trace gradient paths window (Tools > Trace Gradient Paths), the GUI
! counterpart of FLUXPRINT.
submodule (windows) gpaths
  use interfaces_cimgui
  use representations, only: gtrace_point, gtrace_sphere, gtrace_rings,&
     gtrace_axis_axial, gtrace_axis_plane, gpaths_trace
  use param, only: bohrtoa
  implicit none

  ! names of the seed kinds for the combo, in the order of the gtrace_*
  ! values, and their descriptions
  character(len=*), parameter :: kind_names = "Point" // c_null_char //&
     "Sphere around an NCP or CCP" // c_null_char //&
     "Axis and plane of a BCP or RCP" // c_null_char //&
     "Ring paths of all RCPs" // c_null_char
  character(len=*), parameter :: kind_tooltips = &
     "Paths up and/or down the gradient from a point (FLUXPRINT POINT)" // c_null_char //&
     "Paths leaving a nuclear or cage critical point in directions spread over a sphere&
     & (FLUXPRINT NCP, CCP): the skeleton of an atomic basin, or of a cage" // c_null_char //&
     "The two paths along the axis of a bond or ring critical point and the paths in the&
     & perpendicular plane (FLUXPRINT BCP, RCP): the bond paths and the lines of the interatomic&
     & surface of a BCP" // c_null_char //&
     "The two paths along the axis of every ring critical point, which go to the cage critical points&
     & (FLUXPRINT GRAPH 1)" // c_null_char

  ! paths of a BCP/RCP seed, in the order of the gtrace_axis_* values
  character(len=*), parameter :: axis_names = "Axis and plane" // c_null_char //&
     "Axis" // c_null_char // "Plane" // c_null_char
  character(len=*), parameter :: axis_tooltips = &
     "Both the paths along the axis and the paths in the perpendicular plane" // c_null_char //&
     "The two paths along the axis of the critical point (the bond paths of a BCP)" // c_null_char //&
     "The paths in the plane perpendicular to the axis (the interatomic surface of a BCP)" // c_null_char

  ! starting angles in the plane, in the order of the flux_method_* values
  character(len=*), parameter :: method_names = "Uniform" // c_null_char //&
     "Eigenvalue quotient" // c_null_char // "Dynamical" // c_null_char
  character(len=*), parameter :: method_tooltips = &
     "Directions spread uniformly around the critical point (BRAINDEAD)" // c_null_char //&
     "Uniform directions remapped with the quotient of the two Hessian eigenvalues in the plane,&
     & to sample the surface evenly (QUOTIENT)" // c_null_char //&
     "Directions from the linearized dynamical system around the critical point; traces 50 paths&
     & first to estimate the distance to the nearest cage (DYNAMICAL)" // c_null_char

  ! a seed CP picked in the view is the nearest one of the right type
  ! within this distance of the picked point (bohr)
  real*8, parameter :: pick_maxdist = 1d0 / bohrtoa

contains

  !> Draw the trace gradient paths window: the form (field, seed and
  !> its options, style), traced as a blocking job of the main loop
  !> (run_gpaths) into a gradient paths object of the view; the table
  !> of the groups of paths of that object; and the text export.
  module subroutine draw_gpaths(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use representations, only: field_has_cps
    use gui_main, only: ColorGpath
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_setpos_bottomright,&
       iw_field_combo, iw_calcwidth, iw_combo_simple, iw_dragfloat_real8, iw_inputint,&
       iw_radiobutton, iw_checkbox, iw_coloredit, iw_dragfloat_realc, iw_inputtext
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged, hascps, ok, ch, ldum, hov
    integer :: isys, iview, ifield, ihover, irep, iaux
    integer(c_int) :: iup
    type(field_ref) :: fref
    character(len=:), allocatable :: formerr

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)

    ! initialize the form on the first pass
    if (w%firstpass) then
       w%gp = gpaths_state()
       w%gp%rgb = ColorGpath
       w%errmsg = ""
    end if

    ! a new system: the reference field, no seed CP, the export file
    if (goodsys) then
       if (.not.w%gp%sys%is(isys) .or. syschanged) then
          call w%gp%sys%set(isys)
          w%gp%fref = field_ref()
          w%gp%icel = 0
          call cancel_pick(w)
          w%okfile = okfile_default(isys,"structure","paths.txt")
       end if
       if (.not.w%gp%fref%ok(sys(isys))) then
          call fref%set(sys(isys),sys(isys)%iref)
          if (.not.fref%same(w%gp%fref)) then
             w%gp%fref = fref
             w%gp%icel = 0
          end if
       end if
    end if

    ! a pending pick in the view it was armed on
    if (w%gp%picking) then
       if (goodsys .and. iview == w%gp%pickview) then
          call poll_pick(w,isys,iview)
       else
          call cancel_pick(w)
       end if
    end if

    ihover = 0
    if (goodsys) then
       ! the system and the field
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)
       call iw_text("Field",highlight=.true.,alignframe=.true.)
       call igSameLine(0._c_float,-1._c_float)
       ifield = w%gp%fref%id
       if (iw_field_combo("##gpfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
          nonestr="<field not available>")) then
          call w%gp%fref%set(sys(isys),ifield)
          w%gp%icel = 0
          call cancel_pick(w)
       end if
       call iw_tooltip("Field whose gradient paths are traced",ttshown)
       ifield = w%gp%fref%id
       ok = sys(isys)%goodfield(ifield)
       if (.not.ok) &
          call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
    end if

    if (goodsys .and. ok) then
       hascps = field_has_cps(isys,ifield)

       ! the kind of seed
       call iw_text("Seed",highlight=.true.)
       call igPushItemWidth(iw_calcwidth(32,1))
       call iw_combo_simple("##gpkind",kind_names,w%gp%kind,changed=ch,startsatone=.true.,&
          tooltips=kind_tooltips,ttshown=ttshown)
       call igPopItemWidth()
       call iw_tooltip("What the gradient paths start from",ttshown)
       if (ch) call cancel_pick(w)

       ! the seed: a point, or a CP of the field
       if (w%gp%kind == gtrace_point) then
          call iw_text("Position",alignframe=.true.)
          ldum = iw_dragfloat_real8("##gpx",x3=w%gp%x,speed=0.001d0,decimal=4,sameline=.true.)
          if (sys(isys)%c%ismolecule) then
             call iw_tooltip("Cartesian coordinates of the point (Å)",ttshown)
          else
             call iw_tooltip("Fractional coordinates of the point",ttshown)
          end if
          if (iw_button("Pick##gppickpoint",sameline=.true.,disabled=w%gp%picking)) call start_pick(w,iview)
          call iw_tooltip("Pick the point in the view: click an atom or a critical point to use its&
             & position, or empty space to use the clicked point (on the plane through the center&
             & of the scene)",ttshown)
          call iw_text("Direction",alignframe=.true.)
          iup = int(w%gp%iuppoint,c_int)
          if (iw_radiobutton("Up##gpup",int=iup,intval=1_c_int,sameline=.true.)) w%gp%iuppoint = 1
          call iw_tooltip("Trace the path up the gradient (towards the maxima)",ttshown)
          if (iw_radiobutton("Down##gpdown",int=iup,intval=-1_c_int,sameline=.true.)) w%gp%iuppoint = -1
          call iw_tooltip("Trace the path down the gradient (towards the minima)",ttshown)
          if (iw_radiobutton("Both##gpboth",int=iup,intval=0_c_int,sameline=.true.)) w%gp%iuppoint = 0
          call iw_tooltip("Trace both paths",ttshown)
       elseif (.not.hascps .and. w%gp%kind /= gtrace_sphere) then
          call iw_text("This field has no critical points other than the nuclei: search for them in&
             & the Critical Points window",disabled=.true.,wrap=.true.)
          if (iw_button("Open Critical Points...##gpopencp")) &
             iaux = stack_create_window(wintype_cp,.true.,idparent=w%parent(),orraise=-1)
          call iw_tooltip("Open the critical points window for this system",ttshown)
       elseif (w%gp%kind == gtrace_rings) then
          call iw_text("The paths along the axis of every ring critical point of the field",&
             disabled=.true.,wrap=.true.)
       else
          ! the CP: chosen in the table, or picked in the view
          call draw_cp_table(w,isys,ihover)
          if (iw_button("Pick##gppickcp",disabled=w%gp%picking)) call start_pick(w,iview)
          call iw_tooltip("Pick the critical point in the view: click it (or, for a nucleus, its atom)",&
             ttshown)
          if (w%gp%icel > 0) then
             call iw_text("Selected:",sameline=.true.)
             call iw_text(cp_name(isys,ifield,w%gp%icel),highlight=.true.,sameline=.true.)
          end if

          ! the options of the seed
          if (w%gp%kind == gtrace_sphere) then
             call iw_text("Directions",alignframe=.true.)
             ldum = iw_inputint("Azimuthal##gpntheta",w%gp%ntheta,sameline=.true.,width=7)
             call iw_tooltip("Number of azimuthal angles of the starting directions",ttshown)
             ldum = iw_inputint("Polar##gpnphi",w%gp%nphi,sameline=.true.,width=7)
             call iw_tooltip("Number of polar angles of the starting directions, besides the two&
                & poles",ttshown)
             w%gp%ntheta = max(w%gp%ntheta,1)
             w%gp%nphi = max(w%gp%nphi,0)
          else
             call iw_text("Paths",alignframe=.true.)
             call igSameLine(0._c_float,-1._c_float)
             call igPushItemWidth(iw_calcwidth(15,1))
             call iw_combo_simple("##gpaxis",axis_names,w%gp%iupaxis,tooltips=axis_tooltips,ttshown=ttshown)
             call igPopItemWidth()
             call iw_tooltip("Which paths of the critical point are traced",ttshown)
             if (w%gp%iupaxis /= gtrace_axis_axial) then
                ldum = iw_inputint("Number##gpnpts",w%gp%npts,sameline=.true.,width=7)
                call iw_tooltip("Number of paths in the plane (the Dynamical starting angles use&
                   & about this number)",ttshown)
                w%gp%npts = max(w%gp%npts,1)
                call iw_text("Angles",alignframe=.true.)
                call igSameLine(0._c_float,-1._c_float)
                call igPushItemWidth(iw_calcwidth(20,1))
                call iw_combo_simple("##gpmethod",method_names,w%gp%method,tooltips=method_tooltips,&
                   ttshown=ttshown)
                call igPopItemWidth()
                call iw_tooltip("How the starting directions in the plane are chosen",ttshown)
             end if
          end if
       end if

       ! symmetry and style
       ldum = iw_checkbox("All symmetry-equivalent copies##gpsym",w%gp%sym)
       call iw_tooltip("Also draw the paths of the seeds equivalent to this one by symmetry (traced&
          & once and carried by the symmetry operations); otherwise only those of this seed",ttshown)
       call iw_text("Color",alignframe=.true.)
       ldum = iw_coloredit("##gprgb",rgb=w%gp%rgb,sameline=.true.)
       call iw_tooltip("Color of the new paths",ttshown)
       ldum = iw_dragfloat_realc("Radius (Å)##gprad",x1=w%gp%rad,speed=0.002_c_float,&
          min=0.005_c_float,max=1._c_float,scale=real(bohrtoa,c_float),decimal=3,&
          sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Radius of the tubes or spheres of the new paths (Å)",ttshown)

       ! trace
       formerr = form_error(w,isys,ifield,hascps)
       if (iw_button("Trace##gptrace",disabled=(len(formerr) > 0 .or. w%gp%picking))) then
          if (w%request_block()) then
             w%gp%pending_view = iview
             w%errmsg = ""
          end if
       end if
       call iw_tooltip("Trace the gradient paths and add them to the gradient paths object of this&
          & field in the view (created if there is none)",ttshown)
       if (len(formerr) > 0) call iw_text(formerr,disabled=.true.,sameline=.true.)

       ! the groups of paths of the object of this field in the view
       irep = target_object(w,iview,.false.)
       if (irep > 0) then
          associate(r => win(iview)%sc%rep(irep))
            call iw_text("Gradient Paths",highlight=.true.)
            call iw_text("(object: " // trim(r%name) // ")",sameline=.true.)
            ch = iw_checkbox("Show the bond paths (AUTO)##gpshowbond",r%gpaths%showbond)
            call iw_tooltip("Also draw the bond paths that the critical point search found",ttshown)
            ch = draw_gpaths_trace_table(r%gpaths,"gpwin",w%id,ttshown,hov) .or. ch
            if (ch) call win(iview)%sc%undo_note("Edit " // trim(r%name))
            if (ch .or. hov) win(iview)%sc%forcebuildlists = .true.

            ! export
            if (r%gpaths%ntrace > 0) then
               call iw_text("Export",highlight=.true.)
               ldum = iw_inputtext("##gpexportfile",bufsize=1023,texta=w%okfile,width=36)
               call iw_tooltip("Text file the traced gradient paths are written to",ttshown)
               if (iw_button("Browse...##gpexportbrowse",sameline=.true.)) &
                  iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savegpathsfile,&
                  idparent=w%id,orraise=-1)
               call iw_tooltip("Choose the file with a file browser",ttshown)
               if (iw_button("Write##gpexportwrite",sameline=.true.,disabled=(len_trim(w%okfile) == 0))) then
                  call write_paths(isys,r,trim(w%okfile),w%errmsg)
                  if (len(w%errmsg) == 0) call okfile_save_dir(w%okfile)
               end if
               call iw_tooltip("Write the points of the traced gradient paths (Cartesian coordinates,&
                  & and fractional for crystals, and the field value) to the file",ttshown)
               call w%okfile_warn_overwrite()
            end if
          end associate
       end if
    end if

    ! the seed CP under the mouse in the table is highlighted in the view
    if (goodsys) then
       call cps_hover_set(w%gp%fref,iview,(/0,ihover/),w%gp%ihover,w%gp%hoverview,w%gp%hoverview_uid)
    else
       call w%clear_gpaths_hover()
    end if

    ! the error message from the last operation
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)

    ! close button
    call iw_setpos_bottomright(5,1)
    if (iw_button("Close")) doquit = .true.
    if (iw_close_event(w%focused())) doquit = .true.

    if (doquit) then
       call cancel_pick(w)
       call w%end()
    end if

  end subroutine draw_gpaths

  !> Draw the overlay of the pending trace of the trace gradient paths
  !> window.
  module subroutine block_gpaths(w)
    use systems, only: sys, sysc, sys_init
    use utils, only: iw_wait_overlay
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    character(len=:), allocatable :: info
    integer :: isys

    info = ""
    isys = w%gp%sys%get(sys_init)
    if (isys > 0) then
       info = "System: " // string(isys) // ": " // trim(sysc(isys)%seed%name)
       if (sys(isys)%goodfield(w%gp%fref%id)) &
          info = info // newline // "Field:  " // string(w%gp%fref%id) // ": " //&
          trim(sys(isys)%f(w%gp%fref%id)%name)
    end if
    call iw_wait_overlay("Tracing gradient paths...",info,.false.)

  end subroutine block_gpaths

  !> Run the pending trace of the trace gradient paths window, called
  !> by the main loop after the frame with the overlay: trace the
  !> paths of the form and add them, as a new group, to the gradient
  !> paths object of the field in the view the trace was requested
  !> in (created if there is none).
  module subroutine run_gpaths(w)
    use systems, only: sys, sys_init
    class(window), intent(inout), target :: w

    integer :: isys, iview, ifield, irep
    type(gpaths_trace) :: t
    type(rep_gpaths) :: gnew
    character(len=:), allocatable :: errmsg

    ! the system, the field, and the view
    isys = w%gp%sys%get(sys_init)
    iview = w%gp%pending_view
    ifield = w%gp%fref%id
    errmsg = "The system or the field is no longer available"
    if (isys <= 0) goto 999
    if (.not.w%gp%fref%ok(sys(isys))) goto 999
    if (iview < 1 .or. iview > nwin) goto 999
    if (.not.win(iview)%isopen .or. .not.associated(win(iview)%sc)) goto 999
    if (win(iview)%isys /= isys) goto 999
    errmsg = form_error(w,isys,ifield,.true.)
    if (len(errmsg) > 0) goto 999

    ! the recipe and style of the group
    t%kind = w%gp%kind
    if (t%kind == gtrace_point) then
       t%x = sys(isys)%c%c2x(form_to_cart(isys,w%gp%x))
       t%iup = w%gp%iuppoint
    elseif (t%kind /= gtrace_rings) then
       t%x = sys(isys)%f(ifield)%cpcel(w%gp%icel)%x
       t%iup = w%gp%iupaxis
    end if
    t%ntheta = w%gp%ntheta
    t%nphi = w%gp%nphi
    t%npts = w%gp%npts
    t%method = w%gp%method
    t%sym = w%gp%sym
    t%rgb = w%gp%rgb
    t%rad = w%gp%rad

    ! trace into the object of the field in the view, settled for the
    ! current field first (an update could drop stale groups); or, if
    ! there is none, into a new one, made only if the trace works
    irep = target_object(w,iview,.false.)
    if (irep > 0) then
       call win(iview)%sc%rep(irep)%update()
       call win(iview)%sc%rep(irep)%gpaths%run_trace(isys,t,errmsg)
    else
       gnew%fref = w%gp%fref
       call gnew%run_trace(isys,t,errmsg)
    end if
    if (len(errmsg) > 0) goto 999
    if (irep == 0) then
       irep = target_object(w,iview,.true.)
       associate(g => win(iview)%sc%rep(irep)%gpaths)
         call win(iview)%sc%rep(irep)%update()
         g%tfield = gnew%tfield
         g%tcuid = gnew%tcuid
         g%trev = gnew%trev
         g%ttime = gnew%ttime
       end associate
    end if
    call win(iview)%sc%rep(irep)%gpaths%add_trace(t)
    win(iview)%sc%rep(irep)%shown = .true.
    win(iview)%sc%forcebuildlists = .true.
    call win(iview)%sc%undo_note("Trace gradient paths",.false.)
    return

999 continue
    w%errmsg = errmsg

  end subroutine run_gpaths

  !> Stop highlighting the seed CP and the group of paths under the
  !> mouse in the tables of the trace gradient paths window.
  module subroutine clear_gpaths_hover(w)
    class(window), intent(inout), target :: w

    call cps_hover_set(w%gp%fref,0,(/0,0/),w%gp%ihover,w%gp%hoverview,w%gp%hoverview_uid)
    call gpaths_thover_release(w%id)

  end subroutine clear_gpaths_hover

  !> The table of the groups of traced gradient paths of object g
  !> (shared by the trace gradient paths window and the object editor;
  !> id makes the ImGui IDs unique): show, color, radius, description,
  !> and remove. The group under the mouse goes to g%thover, set by
  !> window owner (hoverchanged: the highlight changed); the table of
  !> another window does not clear it. Returns whether the groups
  !> changed.
  module function draw_gpaths_trace_table(g,id,owner,ttshown,hoverchanged) result(changed)
    use utils, only: iw_text, iw_tooltip, iw_checkbox, iw_coloredit, iw_dragfloat_realc,&
       iw_button, iw_highlight_selectable, iw_table_column, iw_table_headers_row, iw_calcheight
    use tools_io, only: string
    type(rep_gpaths), intent(inout) :: g
    character(len=*), intent(in) :: id
    integer, intent(in) :: owner
    logical, intent(inout) :: ttshown
    logical, intent(out) :: hoverchanged
    logical :: changed

    integer :: it, irem, ncol, ihover
    integer(c_int) :: flags
    real(c_float) :: xpos
    character(kind=c_char,len=:), allocatable, target :: str1, suffix
    type(ImVec2) :: sz

    changed = .false.
    hoverchanged = .false.
    ihover = 0
    if (g%ntrace == 0) then
       call iw_text("No traced gradient paths",disabled=.true.)
       call set_thover()
       return
    end if

    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1 = "##tabletraces" // id // c_null_char
    sz%x = 0._c_float
    sz%y = iw_calcheight(min(8,g%ntrace+1),0,.false.)
    irem = 0
    if (igBeginTable(c_loc(str1),5,flags,sz,0._c_float)) then
       ncol = -1
       call iw_table_column("Show",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Color",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Radius (Å)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_column("Paths",icol=ncol,flags=ImGuiTableColumnFlags_WidthStretch)
       call iw_table_column("",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
       call iw_table_headers_row(freezetop=.true.,autofit=.true.)

       do it = 1, g%ntrace
          associate(t => g%trace(it))
            suffix = "_" // id // string(it)
            call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)

            ! show, and the row selectable that highlights the group
            if (igTableSetColumnIndex(0_c_int)) then
               xpos = igGetCursorPosX()
               call iw_text("",alignframe=.true.)
               if (iw_highlight_selectable("##gptrrowsel" // suffix)) ihover = it
               call igSameLine(0._c_float,0._c_float)
               call igSetCursorPosX(xpos)
               if (iw_checkbox("##gptrshown" // suffix,t%shown)) changed = .true.
               call iw_tooltip("Show these gradient paths",ttshown)
            end if

            ! color and radius
            if (igTableSetColumnIndex(1_c_int)) then
               if (iw_coloredit("##gptrcolor" // suffix,rgb=t%rgb)) changed = .true.
               call iw_tooltip("Color of these gradient paths",ttshown)
            end if
            if (igTableSetColumnIndex(2_c_int)) then
               if (iw_dragfloat_realc("##gptrrad" // suffix,x1=t%rad,speed=0.002_c_float,&
                  min=0.005_c_float,max=1._c_float,scale=real(bohrtoa,c_float),decimal=3,&
                  flags=ImGuiSliderFlags_AlwaysClamp)) changed = .true.
               call iw_tooltip("Radius of these gradient paths in Å",ttshown)
            end if

            ! description
            if (igTableSetColumnIndex(3_c_int)) then
               if (t%traced .and. allocated(t%label)) then
                  call iw_text(t%label // merge(" (symmetry)","           ",t%sym))
               else
                  call iw_text("(not traced yet)",disabled=.true.)
               end if
            end if

            ! remove
            if (igTableSetColumnIndex(4_c_int)) then
               if (iw_button("Remove##gptrrem" // suffix,danger=.true.)) irem = it
               call iw_tooltip("Remove these gradient paths",ttshown)
            end if
          end associate
       end do
       call igEndTable()
    end if
    call set_thover()

    ! remove the group (after the table)
    if (irem > 0) then
       call g%remove_trace(irem)
       changed = .true.
    end if

    ! show or hide all
    if (iw_button("Show All##gptrshowall" // id)) then
       g%trace(1:g%ntrace)%shown = .true.
       changed = .true.
    end if
    call iw_tooltip("Show all the traced gradient paths",ttshown)
    if (iw_button("Hide All##gptrhideall" // id,sameline=.true.)) then
       g%trace(1:g%ntrace)%shown = .false.
       changed = .true.
    end if
    call iw_tooltip("Hide all the traced gradient paths",ttshown)
    if (iw_button("Remove All##gptrremall" // id,sameline=.true.,danger=.true.)) then
       call g%drop_traces()
       changed = .true.
    end if
    call iw_tooltip("Remove all the traced gradient paths",ttshown)

  contains
    ! the highlighted group: set by this table, or cleared if this
    ! table set it
    subroutine set_thover()
      if (ihover > 0) then
         hoverchanged = (g%thover /= ihover)
         g%thover = ihover
         g%thowner = owner
      elseif (g%thowner == owner .and. g%thover /= 0) then
         hoverchanged = .true.
         g%thover = 0
         g%thowner = 0
      end if
    end subroutine set_thover
  end function draw_gpaths_trace_table

  !> Stop highlighting the groups of traced gradient paths highlighted
  !> from the table of window owner, in the objects of all the views.
  module subroutine gpaths_thover_release(owner)
    use representations, only: reptype_gpaths
    integer, intent(in) :: owner

    integer :: iv, i

    do iv = 1, nwin
       if (.not.win(iv)%isopen .or. .not.associated(win(iv)%sc)) cycle
       associate(sc => win(iv)%sc)
         do i = 1, sc%nrep
            if (.not.sc%rep(i)%isinit .or. sc%rep(i)%type /= reptype_gpaths) cycle
            if (sc%rep(i)%gpaths%thowner /= owner) cycle
            sc%rep(i)%gpaths%thover = 0
            sc%rep(i)%gpaths%thowner = 0
            sc%forcebuildlists = .true.
         end do
       end associate
    end do

  end subroutine gpaths_thover_release

  !xx! private procedures

  !> The table of the cell CPs of the field that can seed paths of the
  !> chosen kind (NCPs and CCPs for a sphere, BCPs and RCPs for the
  !> axis and plane): click a row to choose the seed. Returns in ihover
  !> the cell CP under the mouse (0 = none).
  subroutine draw_cp_table(w,isys,ihover)
    use systems, only: sys
    use utils, only: iw_text, iw_highlight_selectable, iw_table_column, iw_table_headers_row,&
       iw_calcheight
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    integer, intent(out) :: ihover

    integer :: k, n, ncol
    integer, allocatable :: irow(:)
    integer(c_int) :: flags
    logical :: clicked
    real*8 :: x(3)
    real(c_float) :: xpos
    character(kind=c_char,len=:), allocatable, target :: str1, suffix
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    type(ImVec2) :: sz

    ihover = 0
    associate(f => sys(isys)%f(w%gp%fref%id), c => sys(isys)%c)
      ! the rows
      allocate(irow(f%ncpcel))
      n = 0
      do k = 1, f%ncpcel
         if (.not.cp_kind_ok(w%gp%kind,f%cpcel(k)%typ)) cycle
         n = n + 1
         irow(n) = k
      end do
      if (n == 0) then
         call iw_text("The field has no critical points of this type",disabled=.true.)
         return
      end if

      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      flags = ior(flags,ImGuiTableFlags_ScrollY)
      str1 = "##tablegpcps" // c_null_char
      sz%x = 0._c_float
      sz%y = iw_calcheight(min(8,n+1),0,.false.)
      if (igBeginTable(c_loc(str1),4,flags,sz,0._c_float)) then
         ncol = -1
         call iw_table_column("CP",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Type",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         if (c%ismolecule) then
            call iw_table_column("Position (Å)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         else
            call iw_table_column("Position (fractional)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         end if
         call iw_table_column("Field",icol=ncol,flags=ImGuiTableColumnFlags_WidthStretch)
         call iw_table_headers_row(freezetop=.true.,autofit=.true.)

         clipper = ImGuiListClipper_ImGuiListClipper()
         call ImGuiListClipper_Begin(clipper,n,-1._c_float)
         do while(ImGuiListClipper_Step(clipper))
            call c_f_pointer(clipper,clipper_f)
            do k = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
               call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
               suffix = "_" // string(irow(k))
               if (igTableSetColumnIndex(0_c_int)) then
                  ! the row selectable chooses the seed
                  xpos = igGetCursorPosX()
                  call iw_text("",alignframe=.true.)
                  if (iw_highlight_selectable("##gpcprow" // suffix,clicked=clicked,&
                     selected=(irow(k) == w%gp%icel))) ihover = irow(k)
                  if (clicked) w%gp%icel = irow(k)
                  call igSameLine(0._c_float,0._c_float)
                  call igSetCursorPosX(xpos)
                  call iw_text(cp_name(isys,w%gp%fref%id,irow(k)))
               end if
               if (igTableSetColumnIndex(1_c_int)) &
                  call iw_text(cp_typename(f%cpcel(irow(k))%typ,f%cp(f%cpcel(irow(k))%idx)%isnuc))
               if (igTableSetColumnIndex(2_c_int)) then
                  if (c%ismolecule) then
                     x = (c%x2c(f%cpcel(irow(k))%x) + c%molx0) * bohrtoa
                  else
                     x = f%cpcel(irow(k))%x
                  end if
                  call iw_text(string(x(1),'f',10,4) // " " // string(x(2),'f',10,4) // " " //&
                     string(x(3),'f',10,4))
               end if
               if (igTableSetColumnIndex(3_c_int)) &
                  call iw_text(string(f%cp(f%cpcel(irow(k))%idx)%s%f,'e',decimal=5))
            end do
         end do
         call ImGuiListClipper_End(clipper)
         call ImGuiListClipper_destroy(clipper)
         call igEndTable()
      end if
    end associate

  end subroutine draw_cp_table

  !> Whether a CP of type typ can seed paths of kind ikind.
  logical function cp_kind_ok(ikind,typ)
    integer, intent(in) :: ikind, typ

    if (ikind == gtrace_sphere) then
       cp_kind_ok = (typ == -3 .or. typ == 3)
    else
       cp_kind_ok = (typ == -1 .or. typ == 1)
    end if

  end function cp_kind_ok

  !> Name of the type of a CP of type typ (isnuc: a nucleus).
  function cp_typename(typ,isnuc) result(str)
    integer, intent(in) :: typ
    logical, intent(in) :: isnuc
    character(len=:), allocatable :: str

    select case (typ)
    case (-3)
       str = merge("nucleus","NCP    ",isnuc)
       str = trim(str)
    case (-1)
       str = "BCP"
    case (1)
       str = "RCP"
    case default
       str = "CCP"
    end select

  end function cp_typename

  !> Label of cell CP icel of field ifield (system isys): the name of its
  !> symmetry-unique CP and the cell index.
  function cp_name(isys,ifield,icel) result(str)
    use systems, only: sys
    use tools_io, only: string
    integer, intent(in) :: isys, ifield, icel
    character(len=:), allocatable :: str

    associate(f => sys(isys)%f(ifield))
      str = trim(f%cp(f%cpcel(icel)%idx)%name) // " (" // string(icel) // ")"
    end associate

  end function cp_name

  !> Why the form cannot be traced (empty if it can), for field ifield
  !> of system isys (hascps: the field has CPs other than the nuclei;
  !> a sphere can also start at a nucleus).
  function form_error(w,isys,ifield,hascps) result(err)
    use systems, only: sys
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, ifield
    logical, intent(in) :: hascps
    character(len=:), allocatable :: err

    err = ""
    if (w%gp%kind == gtrace_point) return
    if (.not.hascps .and. w%gp%kind /= gtrace_sphere) then
       err = "No critical points"
       return
    end if
    if (w%gp%kind == gtrace_rings) return
    if (w%gp%icel < 1 .or. w%gp%icel > sys(isys)%f(ifield)%ncpcel) then
       err = "Choose a critical point"
    elseif (.not.cp_kind_ok(w%gp%kind,sys(isys)%f(ifield)%cpcel(w%gp%icel)%typ)) then
       err = "Choose a critical point"
    end if

  end function form_error

  !> The gradient paths object of the field of the form in the scene of
  !> view iview: the first one; or, if there is none and create, a new
  !> one. Returns its index in the scene (0 = none).
  function target_object(w,iview,create) result(irep)
    use systems, only: sys
    use representations, only: reptype_gpaths, repflavor_gpaths
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview
    logical, intent(in) :: create
    integer :: irep

    integer :: i, isys

    irep = 0
    if (iview < 1 .or. iview > nwin) return
    if (.not.win(iview)%isopen .or. .not.associated(win(iview)%sc)) return
    associate(sc => win(iview)%sc)
      do i = 1, sc%nrep
         if (.not.sc%rep(i)%isinit .or. sc%rep(i)%type /= reptype_gpaths) cycle
         if (.not.sc%rep(i)%gpaths%fref%same(w%gp%fref)) cycle
         irep = i
         return
      end do
      if (.not.create) return
      isys = win(iview)%isys
      call sc%add_representation(reptype_gpaths,repflavor_gpaths,id=irep)
      call sc%rep(irep)%gpaths%fref%set(sys(isys),w%gp%fref%id)
    end associate

  end function target_object

  !> Start a pick of the seed in view iview.
  subroutine start_pick(w,iview)
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview

    if (iview < 1 .or. iview > nwin) return
    w%gp%picking = .true.
    w%gp%pickview = iview
    call w%gp%pick%arm()
    if (w%gp%kind == gtrace_point) then
       call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick an atom, a critical point, or a point&
          & to trace the gradient paths from",w%id,acceptempty=.true.)
    else
       call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the critical point to trace the&
          & gradient paths from",w%id,acceptempty=.true.)
    end if

  end subroutine start_pick

  !> Forget the pending pick of window w, releasing the view.
  subroutine cancel_pick(w)
    type(window), intent(inout), target :: w

    if (.not.w%gp%picking) return
    if (w%gp%pickview >= 1 .and. w%gp%pickview <= nwin) &
       call win(w%gp%pickview)%viewmode_release_forced(w%id)
    w%gp%picking = .false.

  end subroutine cancel_pick

  !> Handle the pending pick in view iview (system isys): the picked
  !> point is the point seed; or the seed is the clicked CP, or the
  !> nearest CP of the right type to the picked point (an atom).
  subroutine poll_pick(w,isys,iview)
    use systems, only: sys
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview

    integer :: istat, cpidx(5), k, kmin, it, it2
    real*8 :: xc(3), d, dmin
    integer, parameter :: typs(2,2) = reshape((/-3,3,-1,1/),shape(typs))

    call view_pick_result(iview,w%id,isys,w%gp%pick,istat,xc,cpidx)
    if (istat == ipick_point) then
       if (w%gp%kind == gtrace_point) then
          w%gp%x = form_coords(isys,xc)
       elseif (sys(isys)%goodfield(w%gp%fref%id)) then
          associate(f => sys(isys)%f(w%gp%fref%id))
            kmin = 0
            if (cpidx(1) > 0 .and. cpidx(5) == w%gp%fref%id .and. cpidx(1) <= f%ncpcel) then
               ! a critical point of the field was clicked
               if (cp_kind_ok(w%gp%kind,f%cpcel(cpidx(1))%typ)) kmin = cpidx(1)
            else
               ! the nearest CP of the right types to the picked point
               it = merge(1,2,w%gp%kind == gtrace_sphere)
               dmin = pick_maxdist
               do it2 = 1, 2
                  call f%nearest_cp(sys(isys)%c%c2x(xc),k,d,type=typs(it2,it))
                  if (k > 0 .and. d < dmin) then
                     dmin = d
                     kmin = k
                  end if
               end do
            end if
          end associate
          if (kmin > 0) then
             w%gp%icel = kmin
             w%errmsg = ""
          else
             w%errmsg = "There is no critical point of this type at the picked position"
          end if
       end if
    end if
    if (istat /= ipick_pending) w%gp%picking = .false.

  end subroutine poll_pick

  !> Write the points of the traced gradient paths of the gradient
  !> paths object r (system isys) to text file file: a block per path
  !> (each copy of each site of each group), with the Cartesian
  !> coordinates (Å; for molecules in the frame of the input), the
  !> fractional coordinates (crystals), and the field value. On error,
  !> errmsg is not empty.
  subroutine write_paths(isys,r,file,errmsg)
    use systems, only: sys, sysc
    use tools_io, only: fopen_write, fclose, string
    integer, intent(in) :: isys
    type(representation), intent(in) :: r
    character(len=*), intent(in) :: file
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: lu, it, is, ic, j, k, n
    real*8 :: xf(3), xc(3)
    character(len=:), allocatable :: str

    errmsg = ""
    lu = fopen_write(file,errstop=.false.)
    if (lu < 0) then
       errmsg = "Could not open file: " // file
       return
    end if

    associate(c => sys(isys)%c, g => r%gpaths)
      write (lu,'("# Gradient paths traced with critic2 (Trace Gradient Paths)")')
      write (lu,'("# System: ",A)') trim(sysc(isys)%seed%name)
      if (sys(isys)%goodfield(g%fref%id)) &
         write (lu,'("# Field: ",A,": ",A)') string(g%fref%id), trim(sys(isys)%f(g%fref%id)%name)
      if (c%ismolecule) then
         write (lu,'("# Columns: x y z (Cartesian, Å) field")')
      else
         write (lu,'("# Columns: x y z (Cartesian, Å) xf yf zf (fractional) field")')
      end if
      do it = 1, g%ntrace
         associate(t => g%trace(it))
           if (.not.t%traced) cycle
           str = ""
           if (allocated(t%label)) str = t%label
           write (lu,'("#")')
           write (lu,'("# Group ",A,": ",A)') string(it), str
           do is = 1, t%nsite
              associate(s => t%site(is))
                do ic = 1, s%ncopy
                   do j = 1, s%npath
                      n = s%path(j)%n
                      write (lu,'("# Path ",A," of site ",A,", copy ",A,": ",A," points")') &
                         string(j), string(is), string(ic), string(n)
                      do k = 1, n
                         xf = matmul(s%cop(:,1:3,ic),s%path(j)%x(:,k)) + s%cop(:,4,ic)
                         xc = c%x2c(xf)
                         if (c%ismolecule) then
                            xc = (xc + c%molx0) * bohrtoa
                            write (lu,'(3F16.8,1X,ES18.10)') xc, s%path(j)%f(k)
                         else
                            xc = xc * bohrtoa
                            write (lu,'(3F16.8,3F14.8,1X,ES18.10)') xc, xf, s%path(j)%f(k)
                         end if
                      end do
                      write (lu,*)
                   end do
                end do
              end associate
           end do
         end associate
      end do
    end associate
    call fclose(lu)

  end subroutine write_paths

end submodule gpaths
