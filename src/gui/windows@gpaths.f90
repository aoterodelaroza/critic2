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
  use representations, only: gtrace_sphere, gtrace_axis, gtrace_axis_both, gtrace_axis_axial,&
     gtrace_axis_plane, gpaths_trace
  use param, only: bohrtoa
  implicit none

  ! starting angles in the plane, in the order of the flux_method_* values
  character(len=*), parameter :: method_names = "Uniform" // c_null_char //&
     "Quotient" // c_null_char // "Dynamical" // c_null_char
  character(len=*), parameter :: method_tooltips = &
     "Directions spread uniformly around the critical point (BRAINDEAD)" // c_null_char //&
     "Uniform directions remapped with the quotient of the two Hessian eigenvalues in the plane,&
     & to sample the surface evenly (QUOTIENT)" // c_null_char //&
     "Directions from the linearized dynamical system around the critical point; traces 50 paths&
     & first to estimate the distance to the nearest cage (DYNAMICAL)" // c_null_char

  ! a seed CP picked as an atom in the view is the CP at the position
  ! of the atom (the nucleus), within this distance (bohr)
  real*8, parameter :: pick_maxdist = 1d-3

contains

  !> Draw the trace gradient paths window: the form (field, the table
  !> of the seed CPs, the BCP/RCP and NCP/CCP options, style), traced as a blocking job of the main loop
  !> (run_gpaths) into a gradient paths object of the view; the table
  !> of the groups of paths of that object; and the text export.
  module subroutine draw_gpaths(w)
    use systems, only: sys, sysc, sys_init, ok_system, atlisttype_ncel_frac, atlisttype_nneq
    use representations, only: field_has_cps
    use gui_main, only: ColorGpath, ColorHighlightScene
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_setpos_bottomright,&
       iw_field_combo, iw_calcwidth, iw_combo_simple, iw_inputint, iw_checkbox, iw_coloredit,&
       iw_dragfloat_realc, iw_inputtext
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged, ok, ch, ldum, hov
    integer :: isys, iview, ifield, ihover(2), irep, iaux, itype
    integer, allocatable :: irow(:)
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

    ! a new system: the reference field, the export file
    if (goodsys) then
       if (.not.w%gp%sys%is(isys) .or. syschanged) then
          call w%gp%sys%set(isys)
          w%gp%fref = field_ref()
          call cancel_pick(w)
          w%okfile = okfile_default(isys,"structure","paths.txt")
       end if
       if (.not.w%gp%fref%ok(sys(isys))) then
          call fref%set(sys(isys),sys(isys)%iref)
          if (.not.fref%same(w%gp%fref)) w%gp%fref = fref
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
          call cancel_pick(w)
       end if
       call iw_tooltip("Field whose gradient paths are traced",ttshown)
       ifield = w%gp%fref%id
       ok = sys(isys)%goodfield(ifield)
       if (.not.ok) &
          call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
    end if

    if (goodsys .and. ok) then
       ! the seed CPs: the table of the CPs (the nuclei only, if there
       ! was no search), and the selection
       call sync_sel(w,isys)
       if (.not.sys(isys)%c%ismolecule) then
          if (iw_checkbox("Generate symmetry-equivalent paths##gpsym",w%gp%sym)) call sym_changed(w,isys)
          call iw_tooltip("Trace the paths of the selected symmetry-unique critical points once and&
             & carry them to all the symmetry-equivalent copies in the cell; otherwise, the table lists&
             & every critical point in the cell and only the paths of the selected ones are drawn",ttshown)
       end if
       if (.not.field_has_cps(isys,ifield)) then
          call iw_text("This field has no critical points other than the nuclei: search for them in&
             & the Critical Points window",disabled=.true.,wrap=.true.)
          if (iw_button("Calculate Critical Points##gpopencp")) &
             iaux = stack_create_window(wintype_cp,.true.,idparent=w%parent(),orraise=-1)
          call iw_tooltip("Open the critical points window for this system",ttshown)
       end if
       call table_rows(w,isys,irow)
       call draw_cp_table(w,isys,iview,irow,ihover)
       if (iw_button("All##gpselall")) w%gp%sel(irow) = .true.
       call iw_tooltip("Select all the critical points",ttshown)
       if (iw_button("None##gpselnone",sameline=.true.)) w%gp%sel = .false.
       call iw_tooltip("Unselect all the critical points",ttshown)
       if (iw_button("Toggle##gpseltoggle",sameline=.true.)) w%gp%sel(irow) = .not.w%gp%sel(irow)
       call iw_tooltip("Select the critical points that are not selected, and vice versa",ttshown)
       if (iw_button("Pick##gppickcp",sameline=.true.,disabled=w%gp%picking)) call start_pick(w,iview)
       call iw_tooltip("Pick a critical point in the view and add it to the selection: click it&
          & (or, for a nucleus, its atom)",ttshown)
       call iw_text("(" // string(count(w%gp%sel)) // " selected)",sameline=.true.)

       ! the options of the BCP/RCP seeds
       call iw_text("Bond/Ring",alignframe=.true.)
       call iw_tooltip("The paths traced from the selected bond and ring critical points",ttshown)
       ldum = iw_checkbox("Bond path##gpbondpath",w%gp%bondpath,sameline=.true.)
       call iw_tooltip("Trace the two paths along the axis of the critical point (the bond&
          & paths of a BCP, the ring paths of an RCP)",ttshown)
       ldum = iw_checkbox("Plane##gpplane",w%gp%plane,sameline=.true.)
       call iw_tooltip("Trace the paths in the plane perpendicular to the axis (the&
          & interatomic surface of a BCP, the ring surface of an RCP)",ttshown)
       if (w%gp%plane) then
          ldum = iw_inputint("Number##gpnpts",w%gp%npts,sameline=.true.,width=3)
          call iw_tooltip("Number of paths in the plane (the Dynamical starting angles use&
             & about this number)",ttshown)
          w%gp%npts = max(w%gp%npts,1)
          call iw_text("Angles",sameline=.true.)
          call igSameLine(0._c_float,-1._c_float)
          call igPushItemWidth(iw_calcwidth(10,1))
          call iw_combo_simple("##gpmethod",method_names,w%gp%method,tooltips=method_tooltips,&
             ttshown=ttshown)
          call igPopItemWidth()
          call iw_tooltip("How the starting directions in the plane are chosen",ttshown)
       end if

       ! the options of the NCP/CCP seeds
       call iw_text("Nuc./Cage",alignframe=.true.)
       call iw_tooltip("The paths traced from the selected nuclei and nuclear and cage critical&
          & points, in directions spread over a sphere",ttshown)
       ldum = iw_inputint("Azimuthal##gpntheta",w%gp%ntheta,sameline=.true.,width=3)
       call iw_tooltip("Number of azimuthal angles of the starting directions",ttshown)
       ldum = iw_inputint("Polar##gpnphi",w%gp%nphi,sameline=.true.,width=3)
       call iw_tooltip("Number of polar angles of the starting directions, besides the two&
          & poles",ttshown)
       w%gp%ntheta = max(w%gp%ntheta,1)
       w%gp%nphi = max(w%gp%nphi,0)

       ! style
       call iw_text("Color",alignframe=.true.)
       ldum = iw_coloredit("##gprgb",rgb=w%gp%rgb,sameline=.true.)
       call iw_tooltip("Color of the new paths",ttshown)
       ldum = iw_dragfloat_realc("Radius (Å)##gprad",x1=w%gp%rad,speed=0.002_c_float,&
          min=0.005_c_float,max=1._c_float,scale=real(bohrtoa,c_float),decimal=3,&
          sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Radius of the tubes or spheres of the new paths (Å)",ttshown)

       ! trace
       formerr = form_error(w,isys)
       if (iw_button("Trace##gptrace",disabled=(len(formerr) > 0 .or. w%gp%picking))) then
          if (w%request_block()) then
             w%gp%pending_view = iview
             w%errmsg = ""
          end if
       end if
       call iw_tooltip("Trace the gradient paths of the selected critical points and add them to the&
          & gradient paths object of this field in the view (created if there is none), one group&
          & per critical point (per symmetry-unique one if all the copies are drawn)",ttshown)
       if (len(formerr) > 0) call iw_text(formerr,disabled=.true.,sameline=.true.)

       ! the groups of paths of the object of this field in the view
       irep = target_object(w,iview,.false.)
       if (irep > 0) then
          associate(r => win(iview)%sc%rep(irep))
            call iw_text("Gradient Paths",highlight=.true.)
            call iw_text("(object: " // trim(r%name) // ")",sameline=.true.)
            ch = iw_checkbox("Show the bond paths##gpshowbond",r%gpaths%showbond)
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

    ! the seed CP under the mouse in the table is highlighted in the
    ! view: (symmetry-unique CP,0) in symmetry mode, (0,cell CP)
    ! otherwise; a nucleus, as its atom
    if (goodsys) then
       iaux = 0
       if (ihover(1) > 0 .and. ihover(1) <= sys(isys)%c%nneq) then
          iaux = 1
          itype = atlisttype_nneq
       elseif (ihover(2) > 0 .and. ihover(2) <= sys(isys)%c%ncel) then
          iaux = 2
          itype = atlisttype_ncel_frac
       end if
       if (iaux > 0) then
          call sysc(isys)%highlight_atoms(.true.,ihover(iaux:iaux),itype,&
             reshape(ColorHighlightScene,(/4,1/)))
          ihover = 0
       end if
       call cps_hover_set(w%gp%fref,iview,ihover,w%gp%ihover,w%gp%hoverview,w%gp%hoverview_uid)
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
  !> paths of each selected CP (of each symmetry-unique one, if all the
  !> copies are drawn) and add them, as a new group, to the gradient
  !> paths object of the field in the view the trace was requested
  !> in (created if there is none).
  module subroutine run_gpaths(w)
    use systems, only: sys, sys_init
    class(window), intent(inout), target :: w

    integer :: isys, iview, ifield, irep, k, ntraced
    type(gpaths_trace) :: t, t0
    type(rep_gpaths) :: gnew
    character(len=:), allocatable :: errmsg, err1

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
    errmsg = form_error(w,isys)
    if (len(errmsg) > 0) goto 999

    ! the object of the field in the view, settled for the current
    ! field first (an update could drop stale groups)
    irep = target_object(w,iview,.false.)
    if (irep > 0) call win(iview)%sc%rep(irep)%update()

    ! the options and style shared by all the groups (neither Bond
    ! path nor Plane only reaches here with no BCP/RCP selected)
    if (.not.w%gp%plane) then
       t0%iup = gtrace_axis_axial
    elseif (.not.w%gp%bondpath) then
       t0%iup = gtrace_axis_plane
    else
       t0%iup = gtrace_axis_both
    end if
    t0%ntheta = w%gp%ntheta
    t0%nphi = w%gp%nphi
    t0%npts = w%gp%npts
    t0%method = w%gp%method
    t0%sym = sym_mode(w,isys)
    t0%rgb = w%gp%rgb
    t0%rad = w%gp%rad

    ! one group per selected CP (in symmetry mode, only the first copy
    ! of each symmetry-unique CP can be selected); on error, go on with
    ! the next one and report the first error
    ntraced = 0
    errmsg = ""
    associate(f => sys(isys)%f(ifield))
      do k = 1, f%ncpcel
         if (.not.w%gp%sel(k)) cycle
         if (f%cpcel(k)%typ == 0) cycle

         ! the recipe: the seed
         t = t0
         if (abs(f%cpcel(k)%typ) == 1) then
            t%kind = gtrace_axis
         else
            t%kind = gtrace_sphere
         end if
         t%x = f%cpcel(k)%x

         ! trace into the object; or, if there is none, into a new one,
         ! made only if the trace works
         if (irep > 0) then
            call win(iview)%sc%rep(irep)%gpaths%run_trace(isys,t,err1)
         else
            gnew%fref = w%gp%fref
            call gnew%run_trace(isys,t,err1)
         end if
         if (len(err1) > 0) then
            if (len(errmsg) == 0) errmsg = err1
            cycle
         end if
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
         ntraced = ntraced + 1
      end do
    end associate

    if (ntraced > 0) then
       win(iview)%sc%rep(irep)%shown = .true.
       win(iview)%sc%forcebuildlists = .true.
       call win(iview)%sc%undo_note("Trace gradient paths",.false.)
    end if
    if (len(errmsg) > 0) goto 999
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

  !> The table of the CPs of the field, nuclei included, in the colors
  !> of the view iview (the badges): the cell CPs in irow (the
  !> symmetry-unique CPs, as their first cell copy, in symmetry mode).
  !> A click on a row selects or unselects the CP. Returns in ihover
  !> the CP under the mouse: (symmetry-unique CP, 0) in symmetry mode,
  !> (0, cell CP) otherwise, (0,0) if none.
  subroutine draw_cp_table(w,isys,iview,irow,ihover)
    use systems, only: sys, atlisttype_ncel_frac, atlisttype_nneq
    use utils, only: iw_text, iw_highlight_selectable, iw_table_column, iw_table_headers_row,&
       iw_calcheight
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    integer, intent(in) :: irow(:)
    integer, intent(out) :: ihover(2)

    integer :: kr, k, iu, j, ncol, nrow
    integer(c_int) :: flags
    logical :: clicked, symm, hov
    real*8 :: x(3)
    real(c_float) :: xpos
    character(kind=c_char,len=:), allocatable, target :: str1, suffix, lbl
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    type(ImVec2) :: sz

    ihover = 0
    nrow = size(irow,1)
    if (nrow == 0) return
    symm = sym_mode(w,isys)
    associate(f => sys(isys)%f(w%gp%fref%id), c => sys(isys)%c)
      ! the field at the nuclei is calculated on first use
      call f%init_cplist_deferred()

      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      flags = ior(flags,ImGuiTableFlags_ScrollY)
      str1 = "##tablegpcps" // c_null_char
      sz%x = 0._c_float
      sz%y = iw_calcheight(min(10,nrow+1),0,.false.)
      if (igBeginTable(c_loc(str1),3,flags,sz,0._c_float)) then
         ncol = -1
         call iw_table_column("CP",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         if (c%ismolecule) then
            call iw_table_column("Position (Å)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         else
            call iw_table_column("Position (fractional)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         end if
         call iw_table_column("Field",icol=ncol,flags=ImGuiTableColumnFlags_WidthStretch)
         call iw_table_headers_row(freezetop=.true.,autofit=.true.)

         clipper = ImGuiListClipper_ImGuiListClipper()
         call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)
         do while(ImGuiListClipper_Step(clipper))
            call c_f_pointer(clipper,clipper_f)
            do kr = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
               k = irow(kr)
               iu = f%cpcel(k)%idx
               call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
               suffix = "_" // string(k)
               if (igTableSetColumnIndex(0_c_int)) then
                  ! the row selectable selects the seed
                  xpos = igGetCursorPosX()
                  call iw_text("",alignframe=.true.)
                  hov = iw_highlight_selectable("##gpcprow" // suffix,clicked=clicked,selected=w%gp%sel(k))
                  call igSameLine(0._c_float,0._c_float)
                  call igSetCursorPosX(xpos)

                  ! the badge (an atom for the nuclei), which takes the
                  ! mouse from the selectable: it acts as the row
                  lbl = trim(f%cp(iu)%name)
                  if (.not.symm) lbl = lbl // " " // string(k)
                  lbl = lbl // "##gpcpbadge" // suffix
                  if (symm .and. iu <= c%nneq) then
                     call atom_badge(iview,isys,iu,atlisttype_nneq,lbl)
                  elseif (k <= c%ncel) then
                     call atom_badge(iview,isys,k,atlisttype_ncel_frac,lbl)
                  else
                     call cp_badge(iview,isys,w%gp%fref%id,iu,lbl)
                  end if
                  if (igIsItemHovered(ImGuiHoveredFlags_None)) then
                     hov = .true.
                     if (igIsMouseClicked(ImGuiMouseButton_Left,.false._c_bool)) clicked = .true.
                  end if
                  if (hov) then
                     if (symm) then
                        ihover = (/iu,0/)
                     else
                        ihover = (/0,k/)
                     end if
                  end if

                  ! a click toggles the row and anchors a range; a
                  ! shift-click selects the rows from the anchor to this one
                  if (clicked) then
                     if (igIsKeyDown(ImGuiKey_ModShift) .and. w%lastselected >= 1 .and.&
                        w%lastselected <= nrow) then
                        do j = min(w%lastselected,kr), max(w%lastselected,kr)
                           w%gp%sel(irow(j)) = .true.
                        end do
                     else
                        w%gp%sel(k) = .not.w%gp%sel(k)
                        w%lastselected = kr
                     end if
                  end if
               end if
               if (igTableSetColumnIndex(1_c_int)) then
                  if (symm) then
                     x = f%cp(iu)%x
                  else
                     x = f%cpcel(k)%x
                  end if
                  if (c%ismolecule) x = (c%x2c(x) + c%molx0) * bohrtoa
                  call iw_text(string(x(1),'f',10,4) // " " // string(x(2),'f',10,4) // " " //&
                     string(x(3),'f',10,4))
               end if
               if (igTableSetColumnIndex(2_c_int)) &
                  call iw_text(string(f%cp(iu)%s%f,'e',decimal=5))
            end do
         end do
         call ImGuiListClipper_End(clipper)
         call ImGuiListClipper_destroy(clipper)
         call igEndTable()
      end if
    end associate

  end subroutine draw_cp_table

  !> Whether the form of window w (system isys) traces the paths of
  !> the symmetry-unique CPs and carries them to all their copies:
  !> crystals only.
  logical function sym_mode(w,isys)
    use systems, only: sys
    type(window), intent(in) :: w
    integer, intent(in) :: isys

    sym_mode = w%gp%sym .and. .not.sys(isys)%c%ismolecule

  end function sym_mode

  !> The cell CPs listed in the table of window w (system isys), in
  !> irow: the first cell copy of each symmetry-unique CP in symmetry
  !> mode, all of them otherwise.
  subroutine table_rows(w,isys,irow)
    use systems, only: sys
    type(window), intent(in) :: w
    integer, intent(in) :: isys
    integer, allocatable, intent(out) :: irow(:)

    integer :: k, n
    logical, allocatable :: done(:)

    associate(f => sys(isys)%f(w%gp%fref%id))
      if (.not.sym_mode(w,isys)) then
         irow = (/(k, k = 1, f%ncpcel)/)
         return
      end if
      allocate(irow(f%ncpcel),done(f%ncp))
      done = .false.
      n = 0
      do k = 1, f%ncpcel
         if (done(f%cpcel(k)%idx)) cycle
         done(f%cpcel(k)%idx) = .true.
         n = n + 1
         irow(n) = k
      end do
      irow = irow(1:n)
    end associate

  end subroutine table_rows

  !> The symmetry mode of the form of window w (system isys) changed:
  !> entering it, a symmetry-unique CP is selected if any of its copies
  !> was (the selection moves to its first copy, the row of the
  !> table). Clears the shift-click anchor, a row of the table.
  subroutine sym_changed(w,isys)
    use systems, only: sys
    type(window), intent(inout) :: w
    integer, intent(in) :: isys

    integer :: k
    integer, allocatable :: irow(:)
    logical, allocatable :: usel(:)

    w%lastselected = 0
    if (.not.sym_mode(w,isys)) return
    associate(f => sys(isys)%f(w%gp%fref%id))
      allocate(usel(f%ncp))
      usel = .false.
      do k = 1, f%ncpcel
         if (w%gp%sel(k)) usel(f%cpcel(k)%idx) = .true.
      end do
      call table_rows(w,isys,irow)
      w%gp%sel = .false.
      w%gp%sel(irow) = usel(f%cpcel(irow)%idx)
    end associate

  end subroutine sym_changed

  !> Keep the selection of seed CPs of window w in step with the field
  !> of the form and its CP list (system isys): unselect all if either
  !> changed.
  subroutine sync_sel(w,isys)
    use systems, only: sys, sysc
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    integer :: n

    n = sys(isys)%f(w%gp%fref%id)%ncpcel
    if (allocated(w%gp%sel)) then
       if (size(w%gp%sel,1) == n .and. w%gp%seltime == sysc(isys)%timelastchange_cplist .and.&
          w%gp%selfref%same(w%gp%fref)) return
       deallocate(w%gp%sel)
    end if
    allocate(w%gp%sel(n))
    w%gp%sel = .false.
    w%gp%seltime = sysc(isys)%timelastchange_cplist
    w%gp%selfref = w%gp%fref
    w%lastselected = 0

  end subroutine sync_sel

  !> Why the form cannot be traced (empty if it can), for system isys:
  !> no CPs selected, the selection is not for the current field and
  !> CP list, or a BCP/RCP is selected with neither Bond path nor
  !> Plane checked.
  function form_error(w,isys) result(err)
    use systems, only: sys, sysc
    type(window), intent(in), target :: w
    integer, intent(in) :: isys
    character(len=:), allocatable :: err

    err = "Select the critical points"
    if (.not.allocated(w%gp%sel)) return
    if (size(w%gp%sel,1) /= sys(isys)%f(w%gp%fref%id)%ncpcel) return
    if (w%gp%seltime /= sysc(isys)%timelastchange_cplist) return
    if (.not.w%gp%selfref%same(w%gp%fref)) return
    if (.not.any(w%gp%sel)) return
    err = ""
    if (.not.w%gp%bondpath .and. .not.w%gp%plane) then
       associate(f => sys(isys)%f(w%gp%fref%id))
         if (any(w%gp%sel .and. abs(f%cpcel(1:f%ncpcel)%typ) == 1)) err = "Check Bond path or Plane"
       end associate
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

  !> Start a pick of a seed CP in view iview.
  subroutine start_pick(w,iview)
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview

    if (iview < 1 .or. iview > nwin) return
    w%gp%picking = .true.
    w%gp%pickview = iview
    call w%gp%pick%arm()
    call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick a critical point to trace the gradient&
       & paths from",w%id,acceptempty=.true.)

  end subroutine start_pick

  !> Forget the pending pick of window w, releasing the view.
  subroutine cancel_pick(w)
    type(window), intent(inout), target :: w

    if (.not.w%gp%picking) return
    if (w%gp%pickview >= 1 .and. w%gp%pickview <= nwin) &
       call win(w%gp%pickview)%viewmode_release_forced(w%id)
    w%gp%picking = .false.

  end subroutine cancel_pick

  !> Handle the pending pick in view iview (system isys): the clicked
  !> CP of the field, or the CP at the clicked atom (its nucleus), is
  !> added to the selection.
  subroutine poll_pick(w,isys,iview)
    use systems, only: sys
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview

    integer :: istat, cpidx(5), k
    real*8 :: xc(3), d

    call view_pick_result(iview,w%id,isys,w%gp%pick,istat,xc,cpidx)
    if (istat == ipick_point .and. sys(isys)%goodfield(w%gp%fref%id)) then
       call sync_sel(w,isys)
       associate(f => sys(isys)%f(w%gp%fref%id))
         if (cpidx(1) > 0 .and. cpidx(5) == w%gp%fref%id .and. cpidx(1) <= f%ncpcel) then
            ! a critical point of the field was clicked
            k = cpidx(1)
         else
            ! the CP at the picked point (an atom)
            call f%nearest_cp(sys(isys)%c%c2x(xc),k,d)
            if (k > 0 .and. d > pick_maxdist) k = 0
         end if
       end associate
       if (k > 0) then
          if (sym_mode(w,isys)) then
             associate(f => sys(isys)%f(w%gp%fref%id))
               k = findloc(f%cpcel(1:f%ncpcel)%idx,f%cpcel(k)%idx,1)
             end associate
          end if
          w%gp%sel(k) = .true.
          w%lastselected = 0
          w%errmsg = ""
       else
          w%errmsg = "Pick a critical point of this field, or an atom"
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
